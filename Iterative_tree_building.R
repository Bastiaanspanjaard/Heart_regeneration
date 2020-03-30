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

# Dependencies ####
source("~/Documents/Projects/TOMO_scar/Scripts/linnaeus-scripts/scar_helper_functions.R")

# Parameters ####
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
tsne.coord <- tsne.coord.in[tsne.coord.in$orig.ident %in% "Hr10", c("Cell", "Cell.type")]
# tsne.coord.in <- read.csv("./Data/Larvae_data/Larvae_Seurat_batch_r_out_cells_2.csv")
# For Z1
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L1",
#                             c("Cell", "Cluster", "Cell.type")]
# For Z2
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library %in% c("L21", "L22"),
#                             c("Cell", "Cluster", "Cell.type")]
# For Z3
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L3",
#                             c("Cell", "Cluster", "Cell.type")]
# For Z4
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L4", c("Cell", "Cluster", "Cell.type")]
# For Z5
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library == "L5", c("Cell", "Cluster", "Cell.type")]
# mRNA adults
# tsne.coord.in.1 <- 
#   read.csv("./Data/Adult_data/Adults567_brain_cells.csv",
#            stringsAsFactors = F)
# tsne.coord.in.1$Cell.type <-
#   paste(tsne.coord.in.1$Cell.type, "brain")
# tsne.coord.in.2 <- 
#   read.csv("./Data/Adult_data/Adults567_heart_cells.csv",
#            stringsAsFactors = F)
# tsne.coord.in.2$Cell.type <-
#   paste(tsne.coord.in.2$Cell.type, "heart")
# tsne.coord.in.3 <- 
#   read.csv("./Data/Adult_data/Adults567_pancreas_cells.csv",
#            stringsAsFactors = F)
# tsne.coord.in.3$Cell.type <-
#   paste(tsne.coord.in.3$Cell.type, "pancreas")
# tsne.coord.in <- rbind(tsne.coord.in.1, tsne.coord.in.2, tsne.coord.in.3)
# Count total number of cells present even without scars
# For A5
# tsne.coord <- tsne.coord.in[tsne.coord.in$Library %in% c("B5", "H5", "P5"),
#                             c("Cell", "Cell.type")]
# For (simulated) tree B
# N <- 3000 #125 #
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

# larvae.colors <- read.csv("./Data/color_table_larvae.csv",
#                           stringsAsFactors = F)
# colnames(larvae.colors)[2] <- "Cell.type"
# adult.colors <- read.csv("./Data/color_table_adult-2.csv",
#                          stringsAsFactors = F)
# adult.colors$Cell.type <-
#   paste(adult.colors$Cell.type, adult.colors$color)

# Scars
scar.input <- 
  read.csv("./Data/scars/Hr10_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/Simulations/Tree_C2_100cellsout_detection03.csv")
  # read.csv("./Data/Simulations/Tree_B2_2000cellsout_d0_wweakint.csv")
  # read.csv("./Data/Simulations/Tree_B2_2000cellsout_d005_wweakint.csv")
  # read.csv("./Data/2017_10X_1/Z1_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L1", scar.input$Barcode, sep = "_")
  # read.csv("./Data/2017_10X_2/Z2_scars_compared.csv", stringsAsFactors = F)
  # read.csv("./Data/2017_10X_2/Z3_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L3", scar.input$Barcode, sep = "_")
# read.csv("./Data/2017_10X_10_CR/Z4_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L4", scar.input$Barcode, sep = "_")
# read.csv("./Data/2017_10X_10_CR/Z5_scars_compared.csv", stringsAsFactors = F)
# scar.input$Cell <- paste("L5", scar.input$Barcode, sep = "_")

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

# cells.in.tree <- cells.in.tree[!grepl(";", cells.in.tree$Cell), ]

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
  # dfilter.start <- Sys.time()
  
  # dataset.graph <- graph.and.decompose(cells.in.tree)
  # dataset.degrees <-
  #   data.frame(lapply(dataset.graph, function(x) degree(x, mode = "all", loops = F)))
  # colnames(dataset.degrees) <- "Degree"
  # dataset.degrees$Scar <- as.character(rownames(dataset.degrees))
  # dataset.counts <- data.frame(table(cells.in.tree$Scar))
  # colnames(dataset.counts) <- c("Scar", "Count")
  # dataset.dc <- merge(dataset.counts, dataset.degrees)
  # dataset.dc$scardeg.ratio <- dataset.dc$Count/dataset.dc$Degree
  # ggplot(dataset.dc) +
  #   geom_histogram(aes(x = scardeg.ratio))
  # dataset.dc$degscar.ratio <- dataset.dc$Degree/dataset.dc$Count
  # ggplot(dataset.dc) +
  #   geom_histogram(aes(x = log(degscar.ratio)))
  # dataset.dc$log.dsr <- log(dataset.dc$degscar.ratio)
  # 
  # 
  # cd.lm <- lm(Degree ~ Count + 0, dataset.dc)
  # dataset.dc$cd.lm.ratio <- cd.lm$residuals/dataset.dc$Degree
  # dataset.dc$cd.lm <- cd.lm$residuals
  # 
  # dataset.dc$Exp.d.max <- dataset.dc$Count/20
  # dataset.dc$exp.d.max.diff <- dataset.dc$Exp.d.max - dataset.dc$Degree
  
  # 
  # scars.high.ratio <- dataset.dc$Scar[dataset.dc$scardeg.ratio > 3]
  # 
  # cells.in.tree <- cells.in.tree[cells.in.tree$Scar %in% scars.high.ratio, ]
  
  # dc.lm <- lm(Count ~ Degree + 0, dataset.dc)
  # dataset.dc$LM.res <- dc.lm$residuals/dataset.dc$Count
  # dataset.dc$Scaled.res <- dataset.dc$LM.res/dataset.dc$Count
  
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
  
  # Investigate difference between dataset with doublets and without
  # cells.in.tree.no.d <- cells.in.tree[cells.in.tree$Cell.type != "Doublet", ]
  # scar.connections.no.d <- connections.for.graph(cells.in.tree.no.d)
  # only.once.connections.no.d <- data.frame(t(combn(unique(cells.in.tree$Scar), 2)))
  # colnames(only.once.connections.no.d) <- c("Scar.A", "Scar.B")
  # only.once.connections.no.d <- 
  #   merge(only.once.connections.no.d, 
  #         scar.connections.no.d)
  # colnames(only.once.connections.no.d)[3:5] <- c("xnd_A", "xnd_B", "xnd_AB")
  
  
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
  # dfilter.end <- Sys.time()
  # dfilter.time <- dfilter.end - dfilter.start
  # print(dfilter.time)
  
  # Iteration start conditions ####
  # iterative.sc.start <- Sys.time()
  scar.amount <- length(unique(cells.in.tree.f$Scar))
  it.tree.building <- vector("list", scar.amount)
  tree.summary <- 
    initialize.branches(cells.in.tree.f, scar.remove = "Root",
                        size.ratio = branch.size.ratio)
  # iterative.sc.end <- Sys.time()
  # iterative.sc.time <- iterative.sc.end - iterative.sc.start
  # print(iterative.sc.time)
  
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
    
    # cs <- current.cs.component
    # graph <- current.graph
    
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
      # scar.lls.unique <- 
      #   unique(scar.lls[, c("Scar", "Degree", "Scar.count", "Expected.degree", 
      #                       "Degree.p", "Mean.p_A")])
      # scar.lls.select.unique <- 
      #   scar.lls.unique[scar.lls.unique$Mean.p_A > min.detection.rate, ]
      
      if(nrow(scar.lls.select) > 0){
        scar.remove <- scar.lls.select$Scar[1]
        difficult.scars <-
          scar.lls$Scar[scar.lls$Degree.p >
                          scar.lls$Degree.p[scar.lls$Scar == scar.remove]]
        # if(length(difficult.scars) > 0){
        #   print(paste("Weak scar", difficult.scars[1], "found"))
        #   weak.scar.found <- T
        #   weak.scar <- difficult.scars[1]
        #   # print(cat("Scars", difficult.scars, "will be difficult to place"))
        # }
      }else{
        # How to relate this to weak scars?
        print("No scar above minimum detection rate. Taking best scar under minimum detection rate")
        scar.remove <- scar.lls$Scar[1]
      }
      
      it.tree.element <- list(Scar = scar.remove,
                              LLS = scar.lls,
                              LLS.select = scar.lls.select)
      # ,
      #                         LLS.unique = scar.lls.unique,
      #                         LLS.select.unique = scar.lls.select.unique)
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
        # print(paste("Testing difficult scar", difficult.scar))
        
        # remove difficult scar from (sub)graph, see if it is still connected
        difficult.subgraph <-
          graph.with[[which(unlist(lapply(graph.with,
                                          function(x) difficult.scar %in% names(V(x)))))]]
        this.scar.weak <-
          !is.connected(delete.vertices(difficult.subgraph,
                        which(names(V(difficult.subgraph)) == difficult.scar)))
        if(this.scar.weak){
          # print("This scar weak!")
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
# Placeable.main = sum(correct.cell.placement$Main),
# Placeable.off.main = sum(!correct.cell.placement$Main),
# )
# rm(actually.not.conflicting, cells.in.tree.flagged, correct.conflicting,
#    placeable.cells)

# tree.summary.out.2 <- tree.summary.old.pc[, c("Node.2", "Depth", "Main", "Node", "Size")]
# colnames(tree.summary.out.2)[1] <- "Scar"


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
    # collapsed.nodename <- tree.summary.collapse$Node[index]
    # collapsed.nnp.target <- 
    tree.summary.collapse$Node.1[tree.summary.collapse$Node.1 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse$Node.2[tree.summary.collapse$Node.2 %in% 
                                   c(single.name, downstream.name)] <- 
      collapsed.name
    tree.summary.collapse <- tree.summary.collapse[-index, ]
  }else{
    # if(sum(is.na(tree.summary.collapse$Nnp)) == nrow(tree.summary.collapse)){
    #   tree.summary.collapse$Nnp[index] <- 1
    # }else{
    #   tree.summary.collapse$Nnp[index] <- max(tree.summary.collapse$Nnp, na.rm = T) + 1
    # }
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
  # x <- as.character(node.count.working$Node[1])
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
  # cells.not.add <- 
  #   correct.cell.placement.main[!(correct.cell.placement.main$Node %in% 
  #                            unique(c(tree.plot$Child, tree.plot$Parent))),
  #                          c("Node", "Cell", "Cell.type")]
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
## No pie charts
# Without cells
# LINNAEUS.wo <- generate_tree(tree.plot)
# LINNAEUS.wo.wg <- 
#   collapsibleTree(LINNAEUS.wo, root = LINNAEUS.wo$scar, pieNode = F,
#                   pieSummary = F, fill = "fill", nodeSize = "size",
#                   width = 400, height = 400,
#                   collapsed = F, ctypes=unique(LINNAEUS.wo$Get("Cell.types")))
# htmlwidgets::saveWidget(
#   LINNAEUS.wo.wg,
#   file = "~/Documents/Projects/TOMO_scar/Images/Simulations/tree_B2_d005_LINNAEUS_tree.html")

# With cells
# LINNAEUS.with <- generate_tree(tree.plot.cells)
# collapsibleTree(LINNAEUS.with, root = LINNAEUS.wo$scar, pieNode = F,
#                 pieSummary = F, fill = "fill", nodeSize = "size",
#                 width = 800, height = 600,
#                 collapsed = F, ctypes=unique(LINNAEUS.with$Get("Cell.types")))

## With pie charts
tree.plot.cells <- 
  make.edgelist(tree.summary.collapse, node.count.cumulative.agg,
                correct.cell.placement, main.min = 5, off.main.min = 50)

tree.plot.cells.scar.blind <- tree.plot.cells
tree.plot.cells.scar.blind$Scar.acquisition <- ""
# if(sum(tree.plot.cells.scar.blind$Parent == "0") > 1){
#   tree.plot.cells.scar.blind <-
#     rbind(data.frame(Child = 0, Scar.acquisition = "", Parent = "Root", Cell.type = "NA",
#                      fill = "black", size = 1),
#           tree.plot.cells.scar.blind)
#   tree.plot.cells.scar.blind$Parent <- as.character(tree.plot.cells.scar.blind$Parent)
# }
LINNAEUS.pie <- generate_tree(tree.plot.cells.scar.blind)
# save(LINNAEUS.pie, file = "~/Dropbox/scartrace manuscript/collapsibleTrees/Trees/Z5_Ltree_pie.Robj")
# Without cells
LINNAEUS.pie.wg <-
  collapsibleTree(df = LINNAEUS.pie, root = LINNAEUS.pie$scar, pieNode = T,
                  pieSummary = T,collapsed = F,
                  width = 600, height = 500,
                  ctypes = celltype_colors$Cell.type,linkLength=50,
                  ct_colors = celltype_colors$color, angle = pi/2,
                  nodeSize_class = c(10, 20, 35), nodeSize_breaks = c(0, 50, 1000, 1e6))
LINNAEUS.pie.wg
# htmlwidgets::saveWidget(
#   LINNAEUS.pie.wg,
#   file = "~/Documents/Projects/heart_Bo/Images/Trees/tree_Hr10_LINNAEUS_pie.html")
# save(LINNAEUS.pie, file = "~/Documents/Projects/heart_Bo/Data/Trees/Hr10_tree_pie.Robj")
# Without cells but with all information
# tree.plot.cells.all <- tree.plot.cells
# tree.plot.cells.all <-
#   merge(tree.plot.cells,
#         node.count.cumulative.agg.main[, c("Node", "Freq")],
#         by.x = "Child", by.y = "Node", all.x = T)
# tree.plot.cells.all$Freq[is.na(tree.plot.cells.all$Freq)] <- ""
# tree.plot.cells.all$Scar.acquisition[tree.plot.cells.all$Freq != ""] <-
#   paste(tree.plot.cells.all$Scar.acquisition[tree.plot.cells.all$Freq != ""],
#         ", N = ", tree.plot.cells.all$Freq[tree.plot.cells.all$Freq != ""],
#         sep = "")
# LINNAEUS.pie.all.info <- generate_tree(tree.plot.cells.all)
# LINNAEUS.pie.all.info.wg <-
#   collapsibleTree(LINNAEUS.pie.all.info, root = LINNAEUS.pie.all.info$scar,
#                   pieNode = T,
#                   pieSummary = T,collapsed = F,
#                   width = 800, height = 600,
#                   ctypes = larvae.colors$Cell.type,
#                   ct_colors = larvae.colors$color,
#                   nodeSize_class = c(10, 20, 35), nodeSize_breaks = c(0, 50, 1000, 1e6))
# # htmlwidgets::saveWidget(
#   LINNAEUS.pie.all.info.wg,
#   file = "~/Documents/Projects/TOMO_scar/Images/2017_10X_10/tree_Z5_LINNAEUS_pie_scb_info.html")

# With cells
# LINNAEUS.pie.all <- generate_tree(tree.plot.cells)
# LINNAEUS.pie.all.wg <-
#   collapsibleTree(LINNAEUS.pie.all, root = LINNAEUS.pie.all$scar, pieNode = T,
#                   pieSummary = F,collapsed = F,
#                   width = 800, height = 1200,
#                   ctypes = larvae.colors$Cell.type,
#                   ct_colors = larvae.colors$color, nodeSize_sc = 1,
#                   nodeSize_class = c(10, 20, 35), nodeSize_breaks = c(0, 50, 500, 1e6))
# htmlwidgets::saveWidget(
#   LINNAEUS.pie.all.wg,
#   file = "~/Documents/Projects/TOMO_scar/Images/2017_10X_2/tree_Z2_LINNAEUS_pie_all.html")

# Extract tree ####
# print("Visualization of zoomed trees")
# parent.child.scarnodes <-
#   tree.plot.cells.scar.blind[tree.plot.cells.scar.blind$Cell.type == "NA", ]
# 
# useful.colors <- 
#   read.csv("color_table_larva.csv", 
#            stringsAsFactors = F, sep = ";")[, -1]
# colnames(useful.colors)[2] <- "Cell.type"
# 
# zoom.to <- "Lateral plate mesoderm"
# if(zoom.to == "Endoderm"){
#   colors.use <- useful.colors[, c("Cell.type", "zoom1", "color1")]
# }else if(zoom.to == "Neurectoderm"){
#   colors.use <- useful.colors[, c("Cell.type", "zoom2", "color2")]
# }else if(zoom.to == "Neural crest"){
#   colors.use <- useful.colors[, c("Cell.type", "zoom3", "color3")]
# }else if(zoom.to == "Lateral plate mesoderm"){
#   colors.use <- useful.colors[, c("Cell.type", "zoom4", "color4")]
# }
# colors.use <- colors.use[complete.cases(colors.use), ]
# colnames(colors.use)[2:3] <- c("Order", "color")
# colors.use <- colors.use[order(colors.use$Order), ]
# 
# cell.types.mini <- colors.use$Cell.type
# zoom.nodes <- 
#   unique(tree.plot.cells.scar.blind$Parent[tree.plot.cells.scar.blind$Cell.type %in%
#                                              cell.types.mini])
# zoom.parents <-
#   unique(parent.child.scarnodes$Parent[parent.child.scarnodes$Child %in%zoom.nodes])
# repeat{
#   zoom.parents.2 <-
#     unique(parent.child.scarnodes$Parent[parent.child.scarnodes$Child %in%zoom.parents])
#   if(length(zoom.parents) == length(unique(c(zoom.parents, zoom.parents.2)))){
#     rm(zoom.parents.2)
#     break}
#   zoom.parents <- unique(c(zoom.parents, zoom.parents.2))
# }
# zoom.siblings <- parent.child.scarnodes$Child[parent.child.scarnodes$Parent %in% zoom.parents]
# 
# zoom.edges <- make.edgelist(tree.summary.collapse, node.count.cumulative.agg, 
#                             correct.cell.placement,
#                             main.min = 5, off.main.min = 50)
# zoom.edges <- zoom.edges[zoom.edges$Cell.type != "NA" | zoom.edges$Child %in% zoom.siblings, ]
# zoom.edges$Scar.acquisition <- ""

# zoom.nodes <- 
#   unique(tree.plot.cells.scar.blind$Parent[tree.plot.cells.scar.blind$Cell.type %in%
#                                              cell.types.mini])
# zoom.parents <-
#   unique(parent.child.scarnodes$Parent[parent.child.scarnodes$Child %in%zoom.nodes])

# larvae.colors.zoom <- larvae.colors[larvae.colors$layer == "Neural crest", ]
# 
# cell.types.mini <- larvae.colors.zoom$Cell.type
# zoom.nodes <- 
#   unique(tree.plot.cells.scar.blind$Parent[tree.plot.cells.scar.blind$Cell.type %in%
#                                              cell.types.mini])
# zoom.parents <-
#   unique(parent.child.scarnodes$Parent[parent.child.scarnodes$Child %in%zoom.nodes])

# zoom.siblings <- parent.child.scarnodes$Child[parent.child.scarnodes$Parent %in% zoom.parents]
# zoom.edges <- rbind(parent.child.scarnodes[parent.child.scarnodes$Child %in% zoom.siblings, ],
#                     cells.add)
# LINNAEUS.zoom <- generate_tree(zoom.edges)
# LINNAEUS.pie.zoom.wg <-
#   collapsibleTree(LINNAEUS.zoom, root = LINNAEUS.pie$scar, pieNode = T,
#                   pieSummary = T,collapsed = F,
#                   width = 400, height = 400, linkLength = 50,
#                   ctypes = colors.use$Cell.type, angle = pi/2,
#                   ct_colors = colors.use$color,
#                   nodeSize_class = c(10, 20, 35), nodeSize_breaks = c(0, 50, 1000, 1e6))
# htmlwidgets::saveWidget(
#   LINNAEUS.pie.zoom.wg,
#   file = "~/Documents/Projects/TOMO_scar/Images/2017_10X_10/tree_Z5_LINNAEUS_pie_scb_lpm.html")
# sum(zoom.edges$Cell.type %in% cell.types.mini)

# colors.use$Cell.type <- 
#   factor(colors.use$Cell.type,
#          levels = colors.use$Cell.type[order(colors.use$Order)])
# zoom.colors <- colors.use$color
# names(zoom.colors) <- colors.use$Cell.type
# pdf("./Images/neural_crest_neurons_zoom_colors.pdf",
# width = 3, height = 2)
# ggplot(colors.use) +
#   geom_tile(aes(x = "", y = Cell.type, fill = Cell.type)) +
#   scale_fill_manual(values = zoom.colors) +
#   labs(x = "", y = "") +
#   guides(fill = F) +
#   theme(panel.background = element_blank(),
#         panel.border = element_blank(),
#         panel.grid.major.x = element_blank(),
#         panel.grid.major.y = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.y = element_text(size = 8))
# dev.off()

# Write parameters and statistics ####
parst.output <- data.frame(Hr10 = rbind(t(parameters), t(tree.statistics)))
# write.csv(parst.output, "./Data/Trees/Hr10_LINNAUS_par_stat.csv",
#           quote = F)

# Determine inferred scars ####
# Create scar-cell matrix with observed values (1), inferred values (0.5),
# and unobserved and uninferred values (0).
# Start with the cells, their scars and their position in the tree
# observed.tree.scars <- merge(cells.in.tree[, c("Cell", "Scar")],
#                              correct.cell.placement[,
#                                                 c("Cell", "Node")])
# observed.tree.scars <- merge(cells.in.tree[, c("Cell", "Scar")],
#                              cell.tree.position[complete.cases(cell.tree.position),
#                                                 c("Cell", "Node")])
# observed.tree.scars$Observed <- T
# Determine the full parent child structure of the tree.
# full.parent.child <- data.frame(Child = character(),
#                                 Parent = character())
# c.node <- unique(tree.summary.collapse$Node)[1]
# for(c.node in c("0", unique(tree.summary.collapse$Node))){
#   c.children <-
#     unique(tree.summary.collapse$Node[grep(paste(c.node, "_", sep = ""),
#                                            tree.summary.collapse$Node)])
#   parent.child.add <- data.frame(Child = c.children,
#                                  Parent = rep(c.node, length(c.children)))
#   full.parent.child <- rbind(full.parent.child, parent.child.add)
# }
# for(c.node in unique(scar.tree.final$Nodes$Node)){
#   c.children <-
#     unique(scar.tree.final$Nodes$Node[grep(paste(c.node, "_", sep = ""),
#                                            scar.tree.final$Nodes$Node)])
#   parent.child.add <- data.frame(Child = c.children,
#                                  Parent = rep(c.node, length(c.children)))
#   full.parent.child <- rbind(full.parent.child, parent.child.add)
# }

# Determine which scars are inferred for each cell, based on their position.
# Note that a position can include multiple scars.
# This will include upstream scars that were actually observed.
# tree.summary.nn <- data.frame(Node.2 = character(), Node = character())
# for(i in 1:nrow(tree.summary.collapse)){
#   # ts.row <- tree.summary.collapse[i, c("Node.2", "Node")]
#   if(grepl(",", tree.summary.collapse$Node.2[i])){
#     scars <- unlist(strsplit(tree.summary.collapse$Node.2[i], ","))
#     ts.add <- data.frame(Node.2 = scars,
#                          Node = tree.summary.collapse$Node[i])
#   }else{
#     ts.add <- tree.summary.collapse[i, c("Node.2", "Node")]
#   }
#   tree.summary.nn <- rbind(tree.summary.nn, ts.add)
# }
# inferred.tree.scars.1 <-
#   merge(correct.cell.placement[, c("Cell", "Node")], full.parent.child,
#         by.x = "Node", by.y = "Child")
# # Determine both scars inferred through parents and through own node.
# inferred.tree.scars.node <-
#   unique(merge(inferred.tree.scars.1[, c("Node", "Cell")], 
#                tree.summary.nn[, c("Node.2", "Node")]))
# colnames(inferred.tree.scars.node)[3] <- "Scar"
# inferred.tree.scars.parent <-
#   merge(inferred.tree.scars.1, tree.summary.nn[, c("Node.2", "Node")],
#         by.x = "Parent", by.y = "Node")
# inferred.tree.scars.parent <- inferred.tree.scars.parent[, -1]
# colnames(inferred.tree.scars.parent)[3] <- "Scar"
# inferred.tree.scars <- 
#   unique(rbind(inferred.tree.scars.node, inferred.tree.scars.parent))
# inferred.tree.scars$Inferred <- T
# inferred.tree.scars <-
#   merge(cell.tree.position[, c("Cell", "Node")], full.parent.child,
#         by.x = "Node", by.y = "Child")
# inferred.tree.scars <-
#   merge(inferred.tree.scars, scar.tree.final$Nodes,
#         by.x = "Parent", by.y = "Node")
# inferred.tree.scars$Inferred <- T
# 
# # Combine the observed and the inferred scars per cell.
# obsinf.per.cell <- data.frame(table(observed.tree.scars$Cell))
# colnames(obsinf.per.cell) <- c("Cell", "Observed")
# 
# oi.tree.scars <- merge(observed.tree.scars[, c("Cell", "Scar", "Observed")],
#                        inferred.tree.scars[, c("Cell", "Scar", "Inferred")],
#                        all = T)
# oi.tree.scars$Observed[is.na(oi.tree.scars$Observed)] <- F
# oi.tree.scars$Inferred[is.na(oi.tree.scars$Inferred)] <- F
# oi.tree.scars$Needed.inferring <- !(oi.tree.scars$Observed)
# inf.per.cell <- aggregate(oi.tree.scars$Needed.inferring,
#                           by = list(Cell = oi.tree.scars$Cell),
#                           sum)
# colnames(inf.per.cell)[2] <- "Inferred"
# obsinf.per.cell <- merge(obsinf.per.cell, inf.per.cell)
# obsinf.per.cell$Total <- obsinf.per.cell$Observed + obsinf.per.cell$Inferred
# obsinf.per.cell <- obsinf.per.cell[order(-obsinf.per.cell$Total, -obsinf.per.cell$Observed), ]
# obsinf.per.cell$Cell.order <- 1:nrow(obsinf.per.cell)


# obsinf.per.cell$Missing <- this.depth - obsinf.per.cell$Observed - obsinf.per.cell$Inferred
# obsinf.per.cell <- obsinf.per.cell[order(obsinf.per.cell$Missing, -obsinf.per.cell$Observed), ]
# obsinf.per.cell$Cell.order <- 1:nrow(obsinf.per.cell)
# 
# # Plot amount of observed, inferred and missing per cell
# obsinf.melt <- melt(obsinf.per.cell[, c(2, 3, 5)], id.vars = "Cell.order")
# colnames(obsinf.melt)[2:3] <- c("Type", "Scars")
# obsinf.melt$Type <- factor(obsinf.melt$Type,
#                            levels =c("Inferred", "Observed"))
# png("./Images/2017_10X_10/Z5_inferred_scars_LINNAEUS.png",
#     width = 5, height = 3, units = "in",
#     res = 300)
# ggplot(obsinf.melt) +
#   geom_bar(stat = "identity", width = 1,
#            aes(x = Cell.order, y = Scars, fill = Type)) +
#   scale_fill_manual(values = c("yellow", "red")) +
#   labs(x = "", fill = "") +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         panel.grid.major.x = element_blank())
# dev.off()


# Calculate cell type enrichment per node ####
# Does a cell type have more than expected cells in a node (> expected is per
# ratio of the above node). Maybe improve function output to include ratio
# observed/expected instead of all the raw numbers? Or at least change the
# column ordering.
# e.calc <- calculate.node.enrichment(node.counts = node.count.cumulative,
#                                     node.counts.agg = node.count.cumulative.agg)
