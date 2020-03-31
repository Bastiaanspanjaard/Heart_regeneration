# Description ####
# Iterative tree building based on the comparison of degrees and detection rates
# in the scar graph.
# The test we use is an extension of the naive "the-first-scar-should-have-the-
# highest-degree" scheme one would use if there are no dropouts. The extension
# is that we now look for the scar that would have had the highest degree if
# there wouldn't have been dropouts. Once this scar has been determined, it is
# removed from the dataset, in which we again determine the first scar, and so
# on.
# To find the first scar, we simply try out each scar. Under the assumption that
# a scar is the first created scar, we can calculate its detection rate as the
# amount of times we see it coincide with another scar, divided by the total
# amount of times we see cells with other scars. Given this detection rate and
# the amount of times every other scar is observed, we calculate the expected
# number of connections the first scar should have and compare it to the
# observed number of connections. The underlying distribution for the number
# of connections is the Poisson binomial poisson distribution, so we use this
# to calculate the p-value of the observed value under the assumption of the
# expected value. After calculating this p-value for all scars, we select the
# scar with the highest p-value as the first scar.

# Parameters ####
sample_name <- "Hr27"
library_names <- c(sample_name, paste(sample_name, "a", sep = ""), paste(sample_name, "v", sep = ""), paste(sample_name, "b", sep = ""))
# Fraction of doublets expected.
doublet.rate <- 0.1 # Default is 0.1, set to 0 to turn off.
# The minimum detection rate for a scar to be considered as top scar.
min.detection.rate <- 0.05 # Default value is 0.05
# Minimum cell number ratio between branches.
branch.size.ratio <- 0.125 # Default 0.125, set to 0 to turn off
# Maximum scar probability to include scar in tree building
max.scar.p <- 0.01
# Maximum number of embryos a scar can be present in to include in tree building
max.larvae <- 1

parameters <-
  data.frame(Doublet.rate = doublet.rate,
             Min.detection.rate = min.detection.rate,
             Branch.size.ratio = branch.size.ratio,
             Max.scar.p = max.scar.p,
             Max.larvae = max.larvae)

# For testing purposes: how many scars to include in tree building (takes the
# most frequent scars, set to NA to include all)
number.scars <- NA

# Dependencies ####
source("~/Documents/Projects/TOMO_scar/Scripts/linnaeus-scripts/scar_helper_functions.R")

# Load data ####
print("Loading data")
# mRNA larvae
tsne.coord.in <- read.csv("./Data/final_metadata.csv", stringsAsFactors = F)
colnames(tsne.coord.in)[1] <- "Cell"
# Cell type decisions as follows: all T-cell types and all macrophage types become T-cells/macrophages. 
# Split endocardium into 1 and 2, A and V but merge endocardium frzb.
tsne.coord.in$Cell.type <- tsne.coord.in$lineage.ident
tsne.coord.in$Cell.type[grepl("T-cell", tsne.coord.in$Cell.type)] <- "T-cells"
tsne.coord.in$Cell.type[grepl("Macrophage", tsne.coord.in$Cell.type)] <- "Macrophages"
tsne.coord.in$Cell.type[grepl("Endocardium", tsne.coord.in$Cell.type)] <- 
  tsne.coord.in$final.zoom[grepl("Endocardium", tsne.coord.in$Cell.type)]
tsne.coord.in$Cell.type[grepl("Endocardium frzb", tsne.coord.in$Cell.type)] <- "Endocardium (frzb)"
tsne.coord <- tsne.coord.in[tsne.coord.in$orig.ident %in% library_names, c("Cell", "Cell.type")]
# tsne.coord <- tsne.coord.in[tsne.coord.in$orig.ident %in% "Hr10", c("Cell", "Cell.type")]
N <- nrow(tsne.coord)

# Cell type colors
# Color decisions as follows: use color in setFibro and setCM, and if that does not exist, use color in set1.
# For perivascular cells, use the color from the fibroblast-subset.
# Duplicate colors for endocardium 1 and 2.
celltypes <- data.frame(table(tsne.coord.in$Cell.type))
colnames(celltypes)[1] <- c("Cell.type")

celltype_colors <- data.frame(lapply(readRDS("./Data/col.table.rds")[, ], as.character), stringsAsFactors = F)
celltype_colors <- celltype_colors[celltype_colors$set1 != "Perivascular cells" | is.na(celltype_colors$set1), ]
celltype_colors$Cell.type <- celltype_colors$setCM
celltype_colors$Cell.type[is.na(celltype_colors$Cell.type)] <- celltype_colors$setFibro[is.na(celltype_colors$Cell.type)]
# For the remaining cell types, take the color according to set1 if the color has not been taken yet.
celltype_colors$Cell.type[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)] <- 
  celltype_colors$set1[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)]


# celltype_colors$Cell.type[is.na(celltype_colors$Cell.type)] <- celltype_colors$set1[is.na(celltype_colors$Cell.type)]
endo_add <- celltype_colors[celltype_colors$Cell.type %in% c("Endocardium (V)", "Endocardium (A)"), ]
endo_add$Cell.type[endo_add$Cell.type == "Endocardium (V)"] <- 
  "Endocardium 2 (V)"
endo_add$Cell.type[endo_add$Cell.type == "Endocardium (A)"] <- 
  "Endocardium 2 (A)"
celltype_colors$Cell.type[celltype_colors$Cell.type == "Endocardium (V)"] <- 
  "Endocardium 1 (V)"
celltype_colors$Cell.type[celltype_colors$Cell.type == "Endocardium (A)"] <- 
  "Endocardium 1 (A)"
celltype_colors <- rbind(celltype_colors, endo_add)
celltype_colors <- merge(celltype_colors[, c("name", "color", "Cell.type")], celltypes, all.y = T)
# Endocardium frzb A/V still not on here.


# Scars
scar.input <- 
  read.csv(paste("./Data/scars/", sample_name, "_scars_compared.csv", sep = ""),
           stringsAsFactors = F)

scar.input <- merge(scar.input[, c("Cell", "Scar", "Presence", "p")],
                    tsne.coord)

if(!("Cell.type" %in% names(scar.input))){
  scar.input$Cell.type <- "Type.O.Negative"
}
if("p" %in% names(scar.input)){
  cells.in.tree.pre.f <- scar.input[scar.input$p <= max.scar.p, 
                              c("Cell", "Scar", "Cell.type")]
  if("Presence" %in% names(scar.input)){
    cells.in.tree.pre.f <- 
      scar.input[scar.input$p <= max.scar.p & scar.input$Presence <= max.larvae, 
                 c("Cell", "Scar", "Cell.type")]
  }
}else{
  cells.in.tree.pre.f <- scar.input
}
cells.in.tree.pre.f <- cells.in.tree.pre.f[!duplicated(cells.in.tree.pre.f), ]

scar.freqs <- data.frame(table(cells.in.tree.pre.f$Scar))
colnames(scar.freqs)[1] <- "Scar"
scar.freqs <- scar.freqs[order(-scar.freqs$Freq), ]
set.seed(1)
if(is.na(number.scars)){
  include.scars <- scar.freqs$Scar[scar.freqs$Freq > 1]
}else{
  include.scars <- scar.freqs$Scar[1:number.scars]
}
cells.in.tree.pre.f <- cells.in.tree.pre.f[cells.in.tree.pre.f$Scar %in% include.scars, ]


# Start tree building ####
# Try to build a tree and detect any weak scars until a tree has been built
# without finding weak scars.
weak.scars <- character()
repeat{
  print(paste("Tree building without scars", paste(weak.scars, collapse = ", ")))
  weak.scar.found <- F
  # Remove weak scars
  cells.in.tree <- 
    cells.in.tree.pre.f[!(cells.in.tree.pre.f$Scar %in% weak.scars), ]
  
  
  # Filter out low-frequency scar connections ####
  print("Filtering doublets")

  # Count how often every scar-scar connection is seen
  scar.connections <- connections.for.graph(cells.in.tree)
  only.once.connections <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
  # only.once.connections <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
  colnames(only.once.connections) <- c("Scar.A", "Scar.B")
  only.once.connections <- 
    merge(only.once.connections, 
          scar.connections)
  
  # NEW
  only.once.connections <- only.once.connections[only.once.connections$x_AB > 0, ]
  
  # Calculate how many doublets we'd expect given a general doublet rate, and 
  # calculate the doublet rate for every scar connection.
  only.once.connections$AB.doublets <- 
    2 * doublet.rate * only.once.connections$x_A * only.once.connections$x_B/N
  only.once.connections$AB.doublet.rate <- only.once.connections$AB.doublets/N
  
  # Calculate the probability that the connections we see can all be explained by
  # doublets.
  only.once.connections$Doublet.p <-
    apply(only.once.connections[, c("x_AB", "AB.doublet.rate")], 1,
          function(x){
            x_d <- as.integer(x[1])
            p_d <- as.numeric(x[2])
            binom.test(x_d, N, p_d, alternative = "greater")$p.value
          }
    )
  only.once.connections$Doublet.padj <-
    p.adjust(only.once.connections$Doublet.p, "fdr")
  
  ooc.cutoff <- 
    only.once.connections[only.once.connections$Doublet.padj > 0.01, ]
  
  # Identify incorrect connections
  incorrect.connections <- 
    only.once.connections[only.once.connections$Doublet.padj >= 0.01, ]
  
  # Identify and remove cells with incorrect connections
  # ic <- 1
  inc.cells <- character()
  if(nrow(incorrect.connections) > 0){
    for(ic in 1:nrow(incorrect.connections)){
      inc.scar.A <- incorrect.connections$Scar.A[ic]
      inc.scar.B <- incorrect.connections$Scar.B[ic]
      cells.scar.A <- cells.in.tree$Cell[cells.in.tree$Scar == inc.scar.A]
      cells.scar.B <- cells.in.tree$Cell[cells.in.tree$Scar == inc.scar.B]
      
      inc.cells <- unique(c(inc.cells, intersect(cells.scar.A, cells.scar.B)))
    }
  }
  cells.in.tree.f <- cells.in.tree[!(cells.in.tree$Cell %in% inc.cells), ]

  # Iteration start conditions ####
  scar.amount <- length(unique(cells.in.tree.f$Scar))
  it.tree.building <- vector("list", scar.amount)
  tree.summary <- 
    initialize.branches(cells.in.tree.f, scar.remove = "Root",
                        size.ratio = branch.size.ratio)

  # Iterative tree building ####
  scar.index <- 1
  while(scar.index <= scar.amount & !weak.scar.found){ 
    print(paste("Iterative tree building, identifying scar",
                scar.index, "of", scar.amount)) 
    
    # Select the topmost incomplete line in the tree summary to define the scar 
    # graph branch for which to determine the first scar created. Deduce the other 
    # scars removed for this branch and which component corresponds to the branch.
    top.incomplete.edge.index <- min(which(is.na(tree.summary$Node.2)))
    top.incomplete.edge <- tree.summary[top.incomplete.edge.index, ]
    
    # Determine which scars to remove
    starting.scar <- as.character(top.incomplete.edge$Node.1)
    scars.to.remove <- find.scars.to.remove(starting.scar, tree.summary)
    
    # Finding the correct component is not fully straightforward: it may be the
    # second connected component of the third connected component of the seventh
    # connected component, for example, when scars are removed starting from the 
    # top. In tree.summary we only record the component number of the last 
    # branching event, and we record that mostly to ensure we are covering all
    # components. While this means the 'component history' of any connected
    # component at any level can be reconstructed (simply follow the scars back up
    # and read out the component numbers), this is not necessary.
    # To select the correct component, we remove scars in two steps: we first 
    # remove all scars but the bottommost one (the [starting.scar]) and create the
    # graph with those scars removed. As a second step, we select the component
    # that includes the [starting.scar], remove that scar, and find the correct
    # component in the resulting disconnected graph.
    preceding.scars <- setdiff(scars.to.remove, starting.scar)
    current.cs <- cells.in.tree.f[!(cells.in.tree.f$Scar %in% preceding.scars), ]
    first.decomposition <- graph.and.decompose(current.cs)
    for(comp in 1:length(first.decomposition)){
      if(starting.scar %in% V(first.decomposition[[comp]])$name){
        break
      }
    }
    if(starting.scar == "Root"){
      second.decomposition <- 
        first.decomposition
    }else{
      second.decomposition <- 
        decompose(delete_vertices(first.decomposition[[comp]], starting.scar))
    }
    current.component <- top.incomplete.edge$Component
    current.graph <- second.decomposition[[current.component]]
    current.cs.component <- current.cs[current.cs$Scar %in% V(current.graph)$name, ]
    
    if(length(V(current.graph)) == 1){
      # Condition for last scar in branch - this will have a graph without any
      # connections, and with one vertex.
      scar.remove <- V(current.graph)$name
      it.tree.element <- list(Scar = scar.remove)
    }else{
      # Test congruence of degree and detection efficiency of all scars, assuming
      # they are the topmost scar of the current data.
      scar.lls <- create.degree.lls(current.cs.component, 
                                    current.graph)
      
      # Calculate the weighted average detection rate of all scars.
      average.det.rates <- 
        ddply(scar.lls, .(Scar),
              function(x) data.frame(Mean.p_A = weighted.mean(x$p_A, x$Total.other)))
      
      # Select only scars whose average detection rate is higher than cutoff
      scar.lls <- merge(scar.lls, average.det.rates)
      scar.lls <- scar.lls[order(-scar.lls$Degree.p,
                                 -scar.lls$Degree,
                                 -scar.lls$Mean.p_A), ]
      scar.lls.select <- scar.lls[scar.lls$Mean.p_A > min.detection.rate, ]

      if(nrow(scar.lls.select) > 0){
        scar.remove <- scar.lls.select$Scar[1]
        difficult.scars <-
          scar.lls$Scar[scar.lls$Degree.p >
                          scar.lls$Degree.p[scar.lls$Scar == scar.remove]]
      }else{
        # How to relate this to weak scars?
        print("No scar above minimum detection rate. Taking best scar under minimum detection rate")
        scar.remove <- scar.lls$Scar[1]
      }
      
      it.tree.element <- list(Scar = scar.remove,
                              LLS = scar.lls,
                              LLS.select = scar.lls.select)
    }
    it.tree.building[[scar.index]] <- it.tree.element
    
    # Add scar to tree.summary as Node.2
    tree.summary$Node.2[top.incomplete.edge.index] <- scar.remove
    
    ## Test if difficult scars would change the scar graph topology. Test for one
    # scar at a time, remove if it would indeed change the topology.
    remaining.cs <- current.cs.component[current.cs.component$Scar != scar.remove, ]

    if(nrow(remaining.cs) > 0){
      if(length(difficult.scars) > 0){
        print("Testing difficult scars")
      }
      graph.with <- graph.and.decompose(remaining.cs)
      for(difficult.scar in difficult.scars){
        # remove difficult scar from (sub)graph, see if it is still connected
        difficult.subgraph <-
          graph.with[[which(unlist(lapply(graph.with,
                                          function(x) difficult.scar %in% names(V(x)))))]]
        this.scar.weak <-
          !is.connected(delete.vertices(difficult.subgraph,
                        which(names(V(difficult.subgraph)) == difficult.scar)))
        if(this.scar.weak){
          weak.scars <- c(weak.scars, difficult.scar)
          weak.scar.found <- T
        }
      }
    }

    # Add newfound scar and components to tree.summary.
    remaining.cs <- current.cs.component[!(current.cs.component$Scar %in% scar.remove), ]
    if(nrow(remaining.cs) > 0){
      tree.summary.add <- 
        initialize.branches(remaining.cs, scar.remove = scar.remove, 
                            size.ratio = branch.size.ratio)
      tree.summary.add$Depth <- tree.summary$Depth[top.incomplete.edge.index] + 1
      # If current branch was not a main branch, set non-main flags for all
      # consecutive scars as well.
      if(!tree.summary$Main[top.incomplete.edge.index]){
        tree.summary.add$Main <- F
      }
      tree.summary <- rbind(tree.summary, tree.summary.add)
    }
    
    scar.index <- scar.index + 1
  }
  if(!weak.scar.found){break}
  print("One or more weak scars found - removing and restarting tree building")
}
rm(tree.summary.add, average.det.rates, current.cs, current.cs.component,
   current.graph, difficult.subgraph, first.decomposition, graph.with,
   it.tree.element, remaining.cs, scar.connections, scar.lls, scar.lls.select,
   second.decomposition, top.incomplete.edge)

# TEST: only consider "main" and "off-main" later
tree.summary <- tree.summary[, -6]

# Place cells in uncollapsed tree ####
print("Placing cells in uncollapsed tree")
# Name nodes - we need node names to easily determine whether doublet-flagged
# cells can be placed or not (all scars need to be in one branch)
# NEW
tree.summary.nodes <- tree.summary[tree.summary$Node.1 == "Root", ]
tree.summary.nodes$Node.uc <-
  paste("0", tree.summary.nodes$Component, sep = "_")
for(d in 1:max(tree.summary$Depth)){
  tree.summary.node.add <- tree.summary[tree.summary$Depth == d, ]
  tree.summary.node.add <-
    merge(tree.summary.node.add, tree.summary.nodes[, c("Node.2", "Node.uc")],
          by.x = "Node.1", by.y = "Node.2")
  tree.summary.node.add$Node.uc <-
    paste(tree.summary.node.add$Node.uc, tree.summary.node.add$Component, sep = "_")
  tree.summary.nodes <- rbind(tree.summary.nodes, tree.summary.node.add)
}
tree.summary <- tree.summary.nodes
rm(tree.summary.node.add, tree.summary.nodes)

# Place correct cells in the lowest possible place
correct.cell.placement.positions <- 
  merge(cells.in.tree.f, tree.summary[, c("Node.2", "Depth")],
        by.x = "Scar", by.y = "Node.2")
correct.cell.depths <-
  aggregate(correct.cell.placement.positions$Depth,
            by = list(Cell = correct.cell.placement.positions$Cell),
            max)
colnames(correct.cell.depths)[2] <- "Depth"                      
correct.cell.placement <- merge(correct.cell.placement.positions, correct.cell.depths)
correct.cell.placement <-
  merge(correct.cell.placement, tree.summary[, c("Node.2", "Node.uc")],
        by.x = "Scar", by.y = "Node.2")
rm(correct.cell.depths, correct.cell.placement.positions)

# Determine which doublet-flagged cells can be placed: 
# unplaceable cells will not have any scars in the actual tree; of the remainder,
#   incorrect cells will have conflicting scar placements; 
#   correct cells are the remaining cells.
if(length(inc.cells) > 0){
  cells.in.tree.flagged <- cells.in.tree[cells.in.tree$Cell %in% inc.cells, ]
  cells.in.tree.flagged <- 
    merge(cells.in.tree.flagged, tree.summary[, c("Node.2", "Depth")],
          by.x = "Scar", by.y = "Node.2", all.x = T)
  unplaceable.cells <- 
    unique(cells.in.tree.flagged$Cell[is.na(cells.in.tree.flagged$Depth)])
  placeable.cells <- cells.in.tree.flagged[!is.na(cells.in.tree.flagged$Depth), ]
  unplaceable.cells <- setdiff(unplaceable.cells, placeable.cells$Cell)
  placeable.cells <- merge(placeable.cells, tree.summary[, c("Node.2", "Node.uc")],
                           by.x = "Scar", by.y = "Node.2")
  correct.conflicting <- unique(placeable.cells[, c("Cell", "Cell.type")])
  correct.conflicting$Node.uc <- NA                               
  correct.conflicting$Node.uc <-
    sapply(correct.conflicting$Cell,
           function(x){
             this.cell <- placeable.cells[placeable.cells$Cell == x, ]
             if(max(table(this.cell$Depth)) > 1){
               # If a cell has more than one scar at the same depth, it belongs
               # to two separate branches and cannot be placed in the tree.
               return("-1")
             }else{
               # Test if a cell has scars that belong to separate branches.
               this.cell <- this.cell[order(this.cell$Depth), ]
               conflict <- F
               for(i in 2:nrow(this.cell)){
                 if(!grepl(this.cell$Node.uc[i-1], this.cell$Node.uc[i])){
                   conflict <- T
                   break
                 }
               }
               if(conflict){
                 return("-1")
               }else{
                 return(this.cell$Node.uc[nrow(this.cell)])
               }
             }
           }
    )
  really.conflicting <- correct.conflicting[correct.conflicting$Node.uc == "-1", ]
  actually.not.conflicting <- correct.conflicting[correct.conflicting$Node.uc != "-1", ]
  actually.not.conflicting <-
    merge(actually.not.conflicting, tree.summary[, c("Node.2", "Depth", "Node.uc")])
  colnames(actually.not.conflicting)[which(colnames(actually.not.conflicting) == "Node.2")] <-
    "Scar"
  
  # Place placeable doublet-flagged cells
  correct.cell.placement <- rbind(correct.cell.placement, actually.not.conflicting)
}else{
  really.conflicting <- cells.in.tree.f[0, ]
  unplaceable.cells <- character()
  actually.not.conflicting <- cells.in.tree.f[0, ]
}
# Note that all cells in correct.cell.placement can now be placed by the
# scar column, which is the lowest scar in the tree that that cell has.

# Calculate tree statistics
tree.statistics <- 
  data.frame(Cells = length(unique(cells.in.tree$Cell)),
             Placeable = nrow(correct.cell.placement),
             Suspected.doublets = nrow(really.conflicting),
             Recovered.suspected.doublets = nrow(actually.not.conflicting),
             Unplaceable = length(unplaceable.cells))

# Collapse tree ####
print("Collapsing tree")
# Remove 'singles' (i.e. scars in tree.summary$Node.1 that only occur once)
# from the tree.summary by collapsing them with their successors while changing
# the name to "[scar 1], [scar 2]". 
# After collapsing, make a dictionary that says which scar belongs to which
# node.
tree.summary.collapse <- tree.summary

# Loop over dataframe to find all singles.
index <- 1
while(index <= nrow(tree.summary.collapse)){
  non.singles <-
    unique(tree.summary.collapse$Node.1[duplicated(tree.summary.collapse$Node.1)])
  if(!(tree.summary.collapse$Node.1[index] %in% non.singles)){
    single.name <- tree.summary.collapse$Node.1[index]
    downstream.name <- tree.summary.collapse$Node.2[index]
    collapsed.name <- paste(single.name, downstream.name, sep = ",")
    tree.summary.collapse$Node.1[tree.summary.collapse$Node.1 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse$Node.2[tree.summary.collapse$Node.2 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse <- tree.summary.collapse[-index, ]
  }else{
    index <- index + 1
  }
  if(nrow(tree.summary.collapse) == 1){break}
}
rm(single.name, downstream.name, collapsed.name)

# New node names
tree.summary.collapse.nodes <- 
  tree.summary.collapse[tree.summary.collapse$Depth == min(tree.summary.collapse$Depth), ]
tree.summary.collapse.nodes$Node <-
  paste("0", tree.summary.collapse.nodes$Component, sep = "_")
for(d in (1 + min(tree.summary.collapse$Depth)):max(tree.summary.collapse$Depth)){
  tree.summary.collapse.node.add <- 
    tree.summary.collapse[tree.summary.collapse$Depth == d, ]
  if(nrow(tree.summary.collapse.node.add) > 0){
    tree.summary.collapse.node.add <- 
      merge(tree.summary.collapse.node.add, 
            tree.summary.collapse.nodes[, c("Node.2", "Node")],
            by.x = "Node.1", by.y = "Node.2")
    tree.summary.collapse.node.add$Node <-
      paste(tree.summary.collapse.node.add$Node, 
            tree.summary.collapse.node.add$Component, sep = "_")
    tree.summary.collapse.nodes <- 
      rbind(tree.summary.collapse.nodes, tree.summary.collapse.node.add)
  }
}
tree.summary.collapse <- tree.summary.collapse.nodes
rm(tree.summary.collapse.nodes, tree.summary.collapse.node.add, d)

# Dictionary of scars to new node names
scar.node.dictionary.1 <- data.frame(Scar = tree.summary.collapse$Node.2,
                                     Node = tree.summary.collapse$Node)
scar.node.dictionary.2 <- scar.node.dictionary.1[0, ]
scar.node.dictionary <- scar.node.dictionary.1[0, ]
for(i in 1:nrow(scar.node.dictionary.1)){
  if(grepl(",", scar.node.dictionary.1$Scar[i])){
    # If the node is collapsed, create a dictionary entry for every scar
    collapsed.scars <- unlist(strsplit(as.character(scar.node.dictionary.1$Scar[i]), ","))
    scar.node.dictionary.2 <- 
      rbind(scar.node.dictionary.2,
            data.frame(Scar = collapsed.scars,
                       Node = scar.node.dictionary.1$Node[i]))
  }else{
    scar.node.dictionary <- 
      rbind(scar.node.dictionary, scar.node.dictionary.1[i, ])
  }
}
scar.node.dictionary <- rbind(scar.node.dictionary, scar.node.dictionary.2)
rm(scar.node.dictionary.1, scar.node.dictionary.2, i, collapsed.scars)

# If any scars are collapsed onto the root node, add entries to the dictionary
# for those as well - the dictionary is used to attach cells to the correct
# nodes in the tree and if the root tree node also has one or more scars,
# cells must be added to that node as well.
root.name <-
  unique(tree.summary.collapse$Node.1[grepl("Root", 
                                            tree.summary.collapse$Node.1)])
if(root.name != "Root"){
  collapsed.scars <- unlist(strsplit(as.character(root.name), ","))[-1]
  root.dictionary <- data.frame(Scar = collapsed.scars,
                                Node = 0)
  scar.node.dictionary <- rbind(root.dictionary, scar.node.dictionary)
  rm(root.dictionary, collapsed.scars)
}
rm(root.name)

# Place cells in collapsed tree using dictionary ####
print("Placing cells in collapsed tree")
correct.cell.placement <- merge(correct.cell.placement, scar.node.dictionary)


# Calculate cumulative node counts with and without cell type stratification ####
# Note: while it's not strictly necessary at this stage to calculate the
# cumulative node counts with cell type stratification, this will be useful
# later on for enrichment calculations.
print("Calculating cumulative node sizes")
node.count <- 
  data.frame(table(correct.cell.placement$Node, correct.cell.placement$Cell.type))
colnames(node.count)[1:2] <- c("Node", "Cell.type")
node.count$Node <- as.character(node.count$Node)
node.count$Node.depth <- 
  sapply(node.count$Node, 
         function(x){
           y <- unlist(strsplit(x, "_"))
           return(length(y) - 1)
         }
  )

node.count.working <- node.count[, c("Node", "Cell.type", "Freq", "Node.depth")]
node.count.cumulative <- node.count[0, c("Node", "Cell.type", "Freq", "Node.depth")]

# Do for decreasing depth (start with maximal):
depth.order <- 
  sort(unique(node.count$Node.depth), 
        decreasing = T)[-length(unique(node.count$Node.depth))]
for(c.depth in depth.order){
  # Add nodes of current depth to separate (final output) dataframe.
  node.count.cumulative <- 
    rbind(node.count.working[node.count.working$Node.depth==c.depth, 
                             c("Node", "Cell.type", "Freq", "Node.depth")],
          node.count.cumulative)
  # Change nodes of current depth in working dataframe: depth-- and remove the
  # last node. Then aggregate frequencies.
  node.count.working$Node[node.count.working$Node.depth == c.depth] <- 
    sapply(node.count.working$Node[node.count.working$Node.depth == c.depth],
           function(x) {
             y <- unlist(strsplit(as.character(x),"_"))
             z <- y[-length(y)]
             return(paste(z, collapse = "_"))
           }
    )
  node.count.working$Node.depth[node.count.working$Node.depth == c.depth] <-
    c.depth - 1
  node.count.working <- 
    aggregate(node.count.working$Freq,
              by = list(Node = node.count.working$Node,
                        Cell.type = node.count.working$Cell.type,
                        Node.depth = node.count.working$Node.depth),
              sum)
  colnames(node.count.working)[4] <- "Freq"
}
node.count.cumulative <- 
  rbind(node.count.working[node.count.working$Node.depth == min(node.count$Node.depth), 
                           c("Node", "Cell.type", "Freq", "Node.depth")],
        node.count.cumulative)
rm(node.count, depth.order, node.count.working, c.depth)

node.count.cumulative.agg <-
  aggregate(node.count.cumulative$Freq,
            by = list(Node = node.count.cumulative$Node,
                      Node.depth = node.count.cumulative$Node.depth),
            sum)
colnames(node.count.cumulative.agg)[3] <- "Freq"

# Separate main and off-main nodes ####
print("Separating main and off-main")
# Do recursive: calculate what's main and off-main for a given dataframe of
# branches (these should be sister clades), determine main and off-main for
# children, return dataframe with main and off-main marked.
start.branches <- 
  node.count.cumulative.agg[node.count.cumulative.agg$Node.depth == 
                              min(node.count.cumulative.agg$Node.depth), ]

main.branches <- 
  determine.main(aggregates = node.count.cumulative.agg, 
                 current.branches = start.branches,
                 main.parent = T, branch.size.ratio = branch.size.ratio)
node.count.cumulative.agg <- main.branches
node.count.cumulative <-
  merge(node.count.cumulative, node.count.cumulative.agg[, c("Node", "Main")])

node.count.cumulative.agg.main <- 
  node.count.cumulative.agg[node.count.cumulative.agg$Main, ]
node.count.cumulative.main <- 
  node.count.cumulative[node.count.cumulative$Main, ]
tree.summary.collapse.main <- 
  merge(tree.summary.collapse, 
        node.count.cumulative.agg.main)[, c("Node", "Node.1", "Node.2")]

start.branches.2 <- 
  node.count.cumulative.agg[node.count.cumulative.agg$Node.depth == 
                              min(node.count.cumulative.agg$Node.depth), ]
tree.statistics$Main <- sum(start.branches.2$Freq[start.branches.2$Main])
tree.statistics$Off.main <- sum(start.branches.2$Freq[!start.branches.2$Main])
rm(main.branches, start.branches, start.branches.2)

# Move cells on off-main branches to the above on-main branches ####
# Note: this is a cell placement issue, the node counts do not change.

# Determine the lowest 'main' parent of non-main nodes. It is possible they
# do not exist - these non-main nodes are removed.
off.main.nodes <- 
  data.frame(Orig.node = 
               node.count.cumulative.agg[!node.count.cumulative.agg$Main, "Node"],
             stringsAsFactors = F)
off.main.nodes$Parent <-
  sapply(off.main.nodes$Orig.node,
         function(x) {
           y <- unlist(strsplit(x, "_"))
           return(paste(y[-length(y)], collapse = "_"))
         }
  )
repeat{
  off.main.nodes <- merge(off.main.nodes, 
                          node.count.cumulative.agg[, c("Node", "Main")],
                          by.x = "Parent", by.y = "Node")
  if(sum(off.main.nodes$Main) == nrow(off.main.nodes)){break}
  off.main.nodes$Parent[!off.main.nodes$Main] <-
    sapply(off.main.nodes$Parent[!off.main.nodes$Main],
           function(x) {
             y <- unlist(strsplit(x, "_"))
             return(paste(y[-length(y)], collapse = "_"))
           }
    )
  off.main.nodes <- off.main.nodes[, c("Parent", "Orig.node")]
}

# For cells currently placed in non-main nodes, determine to which main node
# they should be attached.
correct.cell.placement.main <- 
  merge(correct.cell.placement, node.count.cumulative.agg[, c("Node", "Main")])
correct.cell.placement.off <- 
  correct.cell.placement.main[!correct.cell.placement.main$Main, ]
correct.cell.placement.main <- 
  correct.cell.placement.main[correct.cell.placement.main$Main, ]
correct.cell.placement.off <-
  merge(off.main.nodes, correct.cell.placement.off[, -which(colnames(correct.cell.placement.off) == "Main")], 
        by.x = "Orig.node", by.y = "Node")
correct.cell.placement.off <-
  correct.cell.placement.off[, -which(colnames(correct.cell.placement.off) == "Orig.node")]
colnames(correct.cell.placement.off)[1] <- "Node"
correct.cell.placement.main <- rbind(correct.cell.placement.main, correct.cell.placement.off)

rm(correct.cell.placement.off, off.main.nodes)

# Make edgelists without and with cells ####
print("Making edgelists")
# Without
tree.plot <- tree.summary.collapse.main[, c("Node", "Node.2")]
colnames(tree.plot) <- c("Child", "Scar.acquisition")
tree.plot$Parent <-
  sapply(tree.plot$Child,
         function(x){
           y <- unlist(strsplit(x, "_"))
           z <- paste(y[-length(y)], collapse = "_")
         }
  )
root.scars <- 
  unique(tree.summary.collapse.main$Node.1[grepl("Root", tree.summary.collapse.main$Node.1)])
if(root.scars != "Root"){
  root.add <- data.frame(Parent = "Root",
                         Child = 0,
                         Scar.acquisition = root.scars,
                         stringsAsFactors = F)
  root.add$Scar.acquisition <- 
    sapply(root.add$Scar.acquisition,
           function(x) paste(unlist(strsplit(x, ","))[-1], collapse = ","))
  tree.plot <- rbind(root.add, tree.plot)
  rm(root.add)
}
tree.plot$Cell.type <- "NA"
tree.plot$fill <- "black"
tree.plot$size <- 1
rm(root.scars)

# With
if("Cell.type" %in% names(correct.cell.placement.main)){
  cells.add <- 
    correct.cell.placement.main[,
                           c("Node", "Cell", "Cell.type")]
}else{
  cells.add <- 
    correct.cell.placement.main[, 
                           c("Node", "Cell")]
  cell.add$Cell.type <- "Cell"
}
cells.add$Scar.acquisition <- ""
colnames(cells.add)[1:2] <- c("Parent", "Child")
cells.add$fill <- "lightgrey"
cells.add$size <- 0.5
cells.add$Child <- 
  sapply(cells.add$Child,
         function(x){
           x <- as.character(x)
           if(grepl(";", x)){
             return(paste(unlist(strsplit(x, ";")), collapse = "d"))
           }else{
             return(x)
           }
         }
  )

tree.plot.cells <- rbind(tree.plot, cells.add)

# Visualize trees ####
print("Visualization of full trees")
## With pie charts
tree.plot.cells <- 
  make.edgelist(tree.summary.collapse, node.count.cumulative.agg,
                correct.cell.placement, main.min = 5, off.main.min = 50)

tree.plot.cells.scar.blind <- tree.plot.cells
tree.plot.cells.scar.blind$Scar.acquisition <- ""
LINNAEUS.pie <- generate_tree(tree.plot.cells.scar.blind)
save(LINNAEUS.pie, file = paste("~/Documents/Projects/heart_Bo/Data/Trees/", sample_name, "_tree_pie.Robj", sep = ""))
# Without cells
LINNAEUS.pie.wg <-
  collapsibleTree(df = LINNAEUS.pie, root = LINNAEUS.pie$scar, pieNode = T,
                  pieSummary = T,collapsed = F,
                  width = 600, height = 500,
                  ctypes = celltype_colors$Cell.type,linkLength=50,
                  ct_colors = celltype_colors$color, angle = pi/2,
                  nodeSize_class = c(10, 20, 35), nodeSize_breaks = c(0, 50, 1000, 1e6))
print(LINNAEUS.pie.wg)
htmlwidgets::saveWidget(
  LINNAEUS.pie.wg,
  file = paste("~/Documents/Projects/heart_Bo/Images/Trees/tree_", sample_name, "_LINNAEUS_pie.html", sep = ""))

# Write parameters and statistics ####
parst.output <- data.frame(rbind(t(parameters), t(tree.statistics)))
write.csv(parst.output, paste("./Data/Trees/", sample_name, "_LINNAUS_par_stat.csv", sep = ""),
          quote = F)
