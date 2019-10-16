# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)
require(weights)

CalculateCooccurrence <- function(tree_sample){
  node_counts <- count_cumulative(tree_sample$Tree)
  
  # For each cell type and each node, if there's 0 of another cell type, calculate the chance of this given
  # the cell type's abundance in the parent node.
  # Nodes that our descendants appear in
  type_counts <- 
    node_counts[node_counts$Type_count != 0, 
                c("Cell_type", "Node", "Parent_node")]
  # Potential precursor node counts
  parent_counts <- node_counts[, c("Node", "Cell_type", "Type_count", "Parent_node")]
  colnames(parent_counts) <-  c("Node", "Precursor", "Precursor_count", "Parent_node")
  # Full node counts
  full_node_counts <- aggregate(node_counts$Type_count,
                                by = list(Node = node_counts$Node),
                                sum)
  colnames(full_node_counts)[2] <- "Total"
  # Cell type proportions
  proportions <- merge(node_counts[, c("Cell_type", "Node", "Type_count")], full_node_counts)
  proportions$Frequency <- proportions$Type_count/proportions$Total
  # Merge so that we have for each node and descendant type: potential precursors, their
  # counts in the node, the node size and the precursor frequency in the parent node.
  descendancy <- 
    merge(
      merge(
        merge(type_counts, parent_counts), 
        full_node_counts),
      proportions[, c("Node", "Cell_type", "Frequency")], by.x = c("Parent_node", "Precursor"),
      by.y = c("Node", "Cell_type"))[, c("Cell_type", "Precursor", 
                                         "Node", "Parent_node", "Precursor_count", 
                                         "Total", "Frequency")]
  descendancy$Precursor_presence_p <- 
    apply(descendancy[, c("Precursor_count", "Total", "Frequency")], 1,
          function(x){
            if(x[1] > 0){
              output <- 1
            }else{
              output <- pbinom(q = 0, size = x[2], prob = x[3], lower.tail = T)
            }
            return(output)
          }
    )
  descendancy_agg <- aggregate(descendancy$Precursor_presence_p,
            by = list(Cell_type = descendancy$Cell_type,
                      Precursor = descendancy$Precursor),
            prod)
  colnames(descendancy_agg)[3] <- "Tree_precursor_p"

  CTN <- acast(node_counts, Cell_type ~ Node, value.var = "Type_count")
  CTN <- ifelse(CTN == 0, 0, 1)
  Cooc_M <- CTN %*% t(CTN)
  Cooc_f_M <- Cooc_M/diag(Cooc_M)
  
  nodes <- node_counts[, c("Node", "Cell_type", "Type_count")]
  baseline <- nodes[nodes$Node == "nd0", c("Node", "Cell_type", "Type_count")]
  colnames(baseline)[3] <- "Total"
  nodes <- merge(nodes, baseline[, c("Cell_type", "Total")])
  nodes$Rel_freq <- nodes$Type_count/nodes$Total
  relfreqsum <-
    aggregate(nodes$Rel_freq,
              by = list(Node = nodes$Node),
              sum)
  colnames(relfreqsum)[2] <- "Relfreqsum"
  nodes <- merge(nodes, relfreqsum)
  nodes$RF_norm <- nodes$Rel_freq/nodes$Relfreqsum
  node_entropy <- aggregate(nodes$RF_norm,
                                by = list(Node = nodes$Node),
                                function(x){
                                  y <- -x*log(x, base = length(x))
                                  y[x==0] <- 0
                                  return(sum(y))
                                })
  colnames(node_entropy)[2] <- "Entropy"
  
  tree_sample$Relative_cooccurrence <- Cooc_f_M
  tree_sample$Descendancy <- descendancy
  tree_sample$Aggregated_descendancy <- descendancy_agg
  tree_sample$Node_entropy <- node_entropy
  return(tree_sample)
}

count_cumulative <- function(tree){
  # For a tree, count the cumulative cell type numbers for each node, the full
  # cumulative cell count for that node, the frequency of the cell type in that
  # node, and determine the parent node.
  
  nodes <- tree$Get('name', filterFun = function(x) {!isLeaf(x)})
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  node_cells <- data.frame(Node = character(),
                           Cell = character(),
                           Cell_type = character())
  for(n in 1:length(nodes)){
    node <- nodes[n]
    node_cells_add <- edge_list[grep(node, edge_list$from), ]
    node_cells_add <- node_cells_add[node_cells_add$Cell.type != "NA", ]
    colnames(node_cells_add) <- c("Node", "Cell", "Cell_type")
    node_cells_add$Node <- node
    node_cells <- rbind(node_cells, node_cells_add)
  }
  node_counts <- data.frame(table(node_cells$Node, node_cells$Cell_type))
  colnames(node_counts) <- c("Node", "Cell_type", "Type_count")
  node_counts$Node <- as.character(node_counts$Node)
  
  node_cocc <- node_counts[, c("Node", "Cell_type", "Type_count")]
  # node_cocc$Ccount <- NA
  # for(nl in 1:nrow(node_cocc)){
  #   n <- node_counts$Node[nl]
  #   ct <- node_counts$Cell_type[nl]
  #   node_cocc$Ccount[nl] <-
  #     sum(node_counts$Type_count[grepl(paste(n, "_|", n, "$", sep = ""), node_counts$Node) &
  #                                  node_counts$Cell_type == ct])
  # }
  node_cocc$Parent_node <-
    sapply(node_cocc$Node,
           function(x){
             if(x == "nd0"){
               parent = NA
             }else{
               y <- unlist(strsplit(x, "_"))
               parent <- paste(y[-length(y)], collapse = "_")
             }
           }
    )
  
  return(node_cocc)
}

RemapCellTypes <- function(node, reference_set){
  for(c in names(node$children)){
    type_bc <- unlist(strsplit(c, "_"))
    if(type_bc[1] == "nd0"){
      next
    }else if(c %in% reference_set$Cell_name){
      next
    }
    node$RemoveChild(c)
  }
  if(node$Cell.type != "NA"){
    cellname <- node$name
    node$Cell.type <- reference_set$Cell_type[reference_set$Cell_name == cellname]
  }
}

ZoomCellTypes <- function(node, zoom_types){
  for(c in node$children){
    type_bc <- unlist(strsplit(c$name, "_"))
    if(type_bc[1] == "nd0"){
      next
    }else if(c$Cell.type %in% zoom_types){
      next
    }
    node$RemoveChild(c$name)
  }
}

ReadTree <- function(library_name, reference_set){
  load(paste("./Data/Trees/", library_name, "_tree_pie.Robj", sep = ""))
  list_out <- list(Tree = LINNAEUS.pie)
  list_out$Tree$Do(RemapCellTypes, 
                   reference_set = reference_set)
  
  return(list_out)
}

MakePieTree <- function(tree_sample, pie_tree_name, types = NULL, ct_colors){
  tree_to_plot <- Clone(tree_sample$Tree)
  if(!is.null(types)){
    tree_to_plot$Do(ZoomCellTypes, 
                           zoom_types = types)
  }
  tree_sample$CTree <-    
    collapsibleTree(df = tree_to_plot, root = tree_to_plot$scar, pieNode = T,
                    pieSummary = T, collapsed = F,
                    width = 600, height = 500,
                    ctypes = type_colors$celltype,linkLength=60, 
                    ct_colors = ct_colors, angle = pi/2,fontSize = 0,
                    nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 

  names(tree_sample)[which(names(tree_sample) == "CTree")] <- pie_tree_name
  
  return(tree_sample)
}

# Prepare colors and cell type data ####
type_colors <- read.csv("./Data/color.table.all.sub.csv", sep = ";", stringsAsFactors = F)[, -1]

cell_types <- read.csv("./Data/all.hearts.all.cells.all.sub.sept03.csv", stringsAsFactors = F)
cell_types$Cell_type <- cell_types$immune.fibro.CM.subtypes
cell_types$Cell_name <- paste("nd", cell_types$X, sep = "")

ph_annotation <- data.frame(Celltype = type_colors$celltype, row.names = type_colors$celltype)
ph_zoom_annotation <- data.frame(Zoomtype = type_colors$celltype, row.names = type_colors$celltype)
ann_colors <-
  list(Celltype = setNames(type_colors$Color2, type_colors$celltype),
       Zoomtype = setNames(type_colors$colo1, type_colors$celltype))

zoom_types <- setdiff(type_colors$celltype[type_colors$colo1 != ""],
                      c("Ery. Duplex", "M. Duplex"))
# vertex_colors <- data.frame(celltype = zoom_types,
#                             Order = 1:(length(zoom_types)))
# vertex_colors <- merge(vertex_colors, type_colors[, c("celltype", "colo1")])
# vertex_colors <- vertex_colors[order(vertex_colors$Order), ]
# rownames(vertex_colors) <- vertex_colors$celltype

target <- "Fibroblast (nppc)"
# target <- "Fibroblast (col11a1a)"
# target <- "Fibroblast (proliferating)"
# target <- "Fibroblast (col12a1a)"
# target <- "Ventricle (ttn.2/aSMA)"
precursors_include <- c("Smooth muscle cells (Vasculature) 2", "Endocardium (Ventricle)")

# Load tree objects, append cell types, create lineage trees ####
# 3dpi
Hr10 <- ReadTree("Hr10", reference_set = cell_types[cell_types$orig.ident == "Hr10", ])
Hr11 <- ReadTree("Hr11", reference_set = cell_types[cell_types$orig.ident == "Hr11", ])
Hr12 <- ReadTree("Hr12", reference_set = cell_types[cell_types$orig.ident == "Hr12", ])
Hr24 <- ReadTree("Hr24", reference_set = cell_types[cell_types$orig.ident == "Hr24", ])
Hr26 <- ReadTree("Hr26", reference_set = cell_types[cell_types$orig.ident == "Hr26", ])
Hr27 <- ReadTree("Hr27", reference_set = cell_types[cell_types$orig.ident == "Hr27", ])
# tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr24 = Hr24, Hr26 = Hr26, Hr27 = Hr27)
# 7dpi
Hr1_ref <- cell_types[cell_types$orig.ident == "Hr1", ]
Hr1_ref$Cell_name <- paste("nd", sapply(Hr1_ref$Cell_name, function(x){unlist(strsplit(x, "_"))[2]}), sep = "")
Hr1 <- ReadTree("Hr1", reference_set = Hr1_ref)
Hr2 <- ReadTree("Hr2", reference_set = cell_types[cell_types$orig.ident %in% c("Hr2a", "Hr2b"), ])
Hr13 <- ReadTree("Hr13", reference_set = cell_types[cell_types$orig.ident == "Hr13", ])
Hr14 <- ReadTree("Hr14", reference_set = cell_types[cell_types$orig.ident == "Hr14", ])
Hr15 <- ReadTree("Hr15", reference_set = cell_types[cell_types$orig.ident == "Hr15", ])
# tree_list <- list(Hr1 = Hr1, Hr2 = Hr2, Hr13 = Hr13, Hr14 = Hr14, Hr15 = Hr15)
# 3 and 7dpi
tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr24 = Hr24, Hr26 = Hr26, Hr27 = Hr27,
                  Hr1 = Hr1, Hr2 = Hr2, Hr13 = Hr13, Hr14 = Hr14, Hr15 = Hr15)
# tree_list <- list(Hr2 = Hr2)

# Create tree visualization and zoom visualization
for(t in 1:length(tree_list)){
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Full_tree", ct_colors = type_colors$Color2)
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
                                ct_colors = type_colors$colo1)
  # htmlwidgets::saveWidget(
  #   tree_list[[t]]$Full_tree,
  #   file = paste("~/Documents/Projects/heart_Bo/Images/tree_", 
  #                names(tree_list)[t], "_LINNAEUS_pie.html", sep = ""))
  # htmlwidgets::saveWidget(
  #   tree_list[[t]]$Fibrozoom_tree,
  #   file = paste("~/Documents/Projects/heart_Bo/Images/tree_",
  #                names(tree_list)[t], "_LINNAEUS_pie_fibrozoom.html", sep = ""))
}

# Finding precursor candidates ####
# Start by calculating the correlations between a progenitor cell type and potential
# precursors. After this, investigate whether these correlations are stable when leaving
# out one of the trees. Finally, bootstrap the correlations by randomly placing the
# precursor cells over the trees to see what the null-distribution of the correlations
# is and how enriched the observed correlations are.

# ** Calculate correlations between cell types over all trees ####
# Run this over all trees, keep distances top-level, create a new tree-level in the list
comparison_list <- list(Comparison = data.frame(Tree = character(),
                                                Precursor = character(),
                        #                        Simple_distance = numeric(),
                        #                        Distance = numeric(),
                        #                        Correlation = numeric(),
                                                Type_count = integer()),
                        Node_sizes = data.frame(Size = integer()),
                        Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = length(unique(cell_types$Cell_type)))))#,
                        #PNormalized_frequencies = data.frame(matrix(nrow = 0, ncol = length(unique(cell_types$Cell_type)))))
colnames(comparison_list$Normalized_frequencies) <- unique(cell_types$Cell_type)
# colnames(comparison_list$PNormalized_frequencies) <- unique(cell_types$Cell_type)

for(t in 1:length(tree_list)){
  tree <- tree_list[[t]]$Tree
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  if(!(target %in% edge_list$Cell.type)){
    next
  }
  sample_type_count <- 
    dcast(data.frame(table(edge_list$from, edge_list$Cell.type)), Var1 ~ Var2, value.var = "Freq")
  rownames(sample_type_count) <- paste(names(tree_list)[t], sample_type_count$Var1, sep = ":")
  sample_type_count <- sample_type_count[, !(colnames(sample_type_count) %in% c("NA", "Var1"))]
  
  sample_type_cumulative <- sample_type_count
  sample_type_cumulative[,] <- NA
  for(i in 1:nrow(sample_type_cumulative)){
    node <- rownames(sample_type_cumulative)[i]
    sample_type_cumulative[i, ] <-
      colSums(sample_type_count[grep(node, rownames(sample_type_count)), ])
  }
  sample_type_nf <- sample_type_cumulative/rowSums(sample_type_cumulative)
  
  # COMPARISON WITH BELOW:
  # nodes <- tree$Get('name', filterFun = function(x) {!isLeaf(x)})
  # edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  # pure_counts <- data.frame(table(edge_list$from, edge_list$Cell.type))
  # colnames(pure_counts) <- c("Node", "Cell_type", "Type_count")
  # pure_sizes <- aggregate(pure_counts$Type_count,
  #                         by = list(Node = pure_counts$Node),
  #                         sum)
  # colnames(pure_sizes)[2] <- "Node_count"
  # pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == target, ])
  # colnames(pure_sizes)[3:4] <- c("Progeny", "Progeny_count")
  
  
  
  # Philipp-normalized frequencies (node frequencies normalized to top node frequencies)
  # sample_type_pn <- sample_type_nf
  # for(n in 1:nrow(sample_type_pn)){
  #   sample_type_pn[n, ] <- sample_type_nf[n, ]/sample_type_nf[1, ]
  # }
  
  comparison_list$Node_sizes <- rbind(comparison_list$Node_sizes,
                                      data.frame(Size = rowSums(sample_type_cumulative)))
  comparison_tree <-
    data.frame(Tree = rep(names(tree_list)[t], ncol(sample_type_count)),
               Precursor = names(sample_type_count),
               # Distance = numeric(ncol(sample_type_count)),
               # Correlation = numeric(ncol(sample_type_count)),
               Type_count = colSums(sample_type_count))
  
  # precursor_list <- list()
  # for(p in 1:ncol(sample_type_count)){#names(sample_type_count))
  #   precursor <- names(sample_type_count)[p]#"Endocardium (Ventricle)"
  #   
  #   correlations <- 
  #     data.frame(Precursor_cor = 
  #                  t(wtd.cors(sample_type_nf[[precursor]], 
  #                             sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
  #                             weight = rowSums(sample_type_cumulative))),
  #                Both_cor = t(wtd.cors(sample_type_nf[[precursor]] + sample_type_nf[[target]], 
  #                                      sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
  #                                      weight = rowSums(sample_type_cumulative))))
  #   comparison_tree$Distance[p] <- dist(t(correlations))
  #   comparison_tree$Correlation[p] <- cor(correlations$Precursor_cor, correlations$Both_cor)
  #   precursor_list[[p]] <- correlations
  #   names(precursor_list)[p] <- precursor
  # }
  
  # comparison_list[[length(comparison_list) + 1]] <- precursor_list
  # names(comparison_list)[length(comparison_list)] <- names(tree_list)[t]
  comparison_list$Comparison <- rbind(comparison_list$Comparison, comparison_tree)
  
  sample_type_nf[setdiff(names(comparison_list$Normalized_frequencies), names(sample_type_nf))] <- NA
  comparison_list$Normalized_frequencies[setdiff(names(sample_type_nf), 
                                                 names(comparison_list$Normalized_frequencies))] <- NA
  # sample_type_pn[setdiff(names(comparison_list$PNormalized_frequencies), names(sample_type_pn))] <- NA
  # comparison_list$PNormalized_frequencies[setdiff(names(sample_type_pn), 
  #                                                names(comparison_list$PNormalized_frequencies))] <- NA
  
  comparison_list$Normalized_frequencies <- 
    rbind(comparison_list$Normalized_frequencies, sample_type_nf)
  # comparison_list$PNormalized_frequencies <-
  #   rbind(comparison_list$PNormalized_frequencies, sample_type_pn)
}
comparison_list$Normalized_frequencies[is.na(comparison_list$Normalized_frequencies)] <- 0
comparison_list$Normalized_frequencies <- 
  comparison_list$Normalized_frequencies[, colSums(comparison_list$Normalized_frequencies) > 0]
# comparison_list$PNormalized_frequencies[is.na(comparison_list$PNormalized_frequencies)] <- 0
# comparison_list$PNormalized_frequencies <- 
#   comparison_list$PNormalized_frequencies[, colSums(comparison_list$PNormalized_frequencies) > 0]
comparison_list$Node_sizes$Weight <- comparison_list$Node_sizes$Size/sum(comparison_list$Node_sizes$Size)
comparison_list$All_trees_prec <- list()
comparison_list$All_trees_distances <- 
  data.frame(Precursor = names(comparison_list$Normalized_frequencies),
             #Weighted_cor = numeric(ncol(comparison_list$Normalized_frequencies)),
             #Reverse_weighted_cor = numeric(ncol(comparison_list$Normalized_frequencies)),
             #PN_weighted_cor = numeric(ncol(comparison_list$Normalized_frequencies)),
             Weighted_cor_progpos = numeric(ncol(comparison_list$Normalized_frequencies)))#,
             #Distance = numeric(ncol(comparison_list$Normalized_frequencies)))

for(p in 1:ncol(comparison_list$Normalized_frequencies)){
  precursor <- names(comparison_list$Normalized_frequencies)[p]
  
  if(precursor == target){
    comparison_list$All_trees_distances$Weighted_cor_progpos[p] <- -2
    next
  }

  # comparison_list$All_trees_distances$Weighted_cor[p] <-
  #   wtd.cors(comparison_list$Normalized_frequencies[[precursor]],
  #          comparison_list$Normalized_frequencies[[precursor]] + comparison_list$Normalized_frequencies[[target]],
  #          weight = comparison_list$Node_sizes$Size)
  
  # comparison_list$All_trees_distances$Reverse_weighted_cor[p] <-
  #   wtd.cors(comparison_list$Normalized_frequencies[[target]],
  #            comparison_list$Normalized_frequencies[[precursor]] + comparison_list$Normalized_frequencies[[target]],
  #            weight = comparison_list$Node_sizes$Size)
  
  # comparison_list$All_trees_distances$PN_weighted_cor[p] <-
  #   wtd.cors(comparison_list$PNormalized_frequencies[[precursor]],
  #            comparison_list$PNormalized_frequencies[[precursor]] + comparison_list$PNormalized_frequencies[[target]],
  #            weight = comparison_list$Node_sizes$Size)
  
  progpos <- comparison_list$Normalized_frequencies[comparison_list$Normalized_frequencies[[target]] > 0, ]
  progpos_weights <- comparison_list$Node_sizes[comparison_list$Normalized_frequencies[[target]] > 0, ]
  comparison_list$All_trees_distances$Weighted_cor_progpos[p] <-
    wtd.cors(progpos[[precursor]],
             progpos[[target]],
             weight = progpos_weights$Size)
  
  # correlations <-
  #   data.frame(
  #     Precursor_cor =
  #       t(wtd.cors(comparison_list$Normalized_frequencies[[precursor]],
  #                  comparison_list$Normalized_frequencies[, !(names(comparison_list$Normalized_frequencies) %in% c(precursor, target))],
  #                  weight = comparison_list$Node_sizes$Size)),
  #     Both_cor = 
  #       t(wtd.cors(comparison_list$Normalized_frequencies[[precursor]] + comparison_list$Normalized_frequencies[[target]],
  #                  comparison_list$Normalized_frequencies[, !(names(comparison_list$Normalized_frequencies) %in% c(precursor, target))],
  #                  weight = comparison_list$Node_sizes$Size)))
  # comparison_list$All_trees_distances$Distance[p] <- dist(t(correlations))
  
  # comparison_list$All_trees_prec[[p]] <- correlations
  # names(comparison_list$All_trees_prec)[p] <- precursor
}  


# View(comparison_list$All_trees_distances)
# Add #trees that were used, #trees that had precursor candidate cell type, average
# cells per tree used.
# View(comparison_list$Comparison)
tree_presence_precursors <- data.frame(table(comparison_list$Comparison$Precursor))
colnames(tree_presence_precursors) <- c("Precursor", "Trees")
mean_precursor_count <- aggregate(comparison_list$Comparison$Type_count,
                                  by = list(Precursor = comparison_list$Comparison$Precursor),
                                  mean)
colnames(mean_precursor_count)[2] <- "Mean_count"
precursor_ranking <-
  merge(comparison_list$All_trees_distances,
        merge(tree_presence_precursors, mean_precursor_count))
precursor_ranking <- precursor_ranking[order(-precursor_ranking$Weighted_cor_progpos), ]
precursor_ranking$Precursor <- factor(precursor_ranking$Precursor, precursor_ranking$Precursor)

# png(paste("./Images/", target, "_precursors_PN_7dpi.png", sep = ""))
print(
  ggplot(precursor_ranking[precursor_ranking$Precursor != target, ]) +
    geom_bar(aes(x = Precursor, y = Weighted_cor_progpos, fill = Precursor, 
                 alpha = ifelse(Mean_count < 10, 0.5, 1)), stat = "identity") +
    scale_fill_manual(values = ann_colors$Celltype) +
    labs(title = paste(target, ""), y = "Fitness") +
    theme(legend.position = "none",
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
)
# dev.off()
potential_prec <- 
  union(precursor_ranking$Precursor[precursor_ranking$Weighted_cor_progpos > 0.5 & 
                                precursor_ranking$Mean_count >= 10],
        precursors_include)
print(potential_prec)

# Plot correlations of two top-scoring precursors and two low-scoring precursors
for(i in c(1:2, nrow(precursor_ranking) - 1, nrow(precursor_ranking))){
  pot_prec <- as.character(precursor_ranking$Precursor[i])
  plot_freqs <- 
    data.frame(comparison_list$Normalized_frequencies)[,make.names(c(pot_prec, target))]
  colnames(plot_freqs) <- c("Precursor", "Progenitor")
  plot_freqs$Size <- comparison_list$Node_sizes$Size
  # png(paste("./Images/", pot_prec, "precursor_to_", target, "_7dpi.png", sep = ""))
  print(
    ggplot() +
      geom_point(data = plot_freqs[plot_freqs$Progenitor == 0, ], 
                 aes(x = Precursor, y = Progenitor, size = Size), color = "grey") +
      geom_point(data = plot_freqs[plot_freqs$Progenitor > 0, ], 
                 aes(x = Precursor, y = Progenitor, size = Size)) +
      labs(x = pot_prec, y = target) +
      theme(legend.position = "none")
  )
  # dev.off()
}

# ** Calculate one-(tree-)out correlations ####

# ** Bootstrap correlation ####
# If I randomly place cells in the tree, what is the correlation with the progeny?
# Select a potential precursor; how many cells in each tree? Randomly place that amount
# of cells in the tree, depending on 'pure' (i.e. not cumulative) node counts when progeny
# is removed. Then calculate correlation with progeny.
# Total node size after sampling is node size before sampling without progeny, with added
# sampled progeny. This means precursor frequencies and node weights change as well.
# Do note that for proper sampling, we need to use the 'pure' (i.e. non-cumulative) node
# sizes, but for correlations we do need to talk about cumulative nodes.

# Loop over trees to calculate node counts per cell type
sample_tree_list <- list()
for(t in 1:length(tree_list)){
  print(paste("Tree", names(tree_list)[t]))
  # Step 0: Pure node counts with and without progeny
  tree <- Clone(tree_list[[t]]$Tree) # Clone(tree_list$Hr26$Tree)
  
  nodes <- tree$Get('name', filterFun = function(x) {!isLeaf(x)})
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  # if(sum(edge_list$Cell.type == target) == 0){
  #   next
  # }
  # if(sum(edge_list$Cell.type == bootstrap_precursor) == 0){
  #   next
  # }
  # sample_tree_list[[length(sample_tree_list) + 1]] <- list()
  # sample_tree_list[[length(sample_tree_list)]]$Tree <- tree
  # names(sample_tree_list)[length(sample_tree_list)] <- names(tree_list)[length(sample_tree_list)]
  pure_counts <- data.frame(table(edge_list$from, edge_list$Cell.type))
  colnames(pure_counts) <- c("Node", "Cell_type", "Type_count")
  pure_counts <- pure_counts[pure_counts$Cell_type != "NA", ]
  
  sample_tree_list[[t]] <- list(Pure_counts = pure_counts)
  names(sample_tree_list)[t] <- names(tree_list)[t]
  # pure_sizes <- aggregate(pure_counts$Type_count,
  #                         by = list(Node = pure_counts$Node),
  #                         sum)
  # colnames(pure_sizes)[2] <- "Node_count"
  # pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == target, ])
  # colnames(pure_sizes)[3:4] <- c("Progeny", "Progeny_count")
  # 
  # # Step I: potential precursor and how many cells?
  # bootstrap_cells <- sum(pure_counts$Count[pure_counts$Cell_type == bootstrap_precursor])
  # 
  # # Step II: Pure node sizes with and without precursor
  # pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == bootstrap_precursor, ])
  # colnames(pure_sizes)[5:6] <- c("Precursor", "Precursor_count")
  # pure_sizes$Sample_chance <- 
  #   (pure_sizes$Node_count - pure_sizes$Progeny_count)/(sum(pure_sizes$Node_count) - sum(pure_sizes$Progeny_count))
  # pure_sizes$Size_no_prec <- pure_sizes$Node_count - pure_sizes$Precursor_count
  # 
  # for(n in 1:nrow(pure_sizes)){
  #   node <- as.character(pure_sizes$Node[n])
  #   
  #   pure_sizes$CNode_count[n] <- sum(pure_sizes$Node_count[grep(node, pure_sizes$Node)])
  #   pure_sizes$CSize_no_prec[n] <- sum(pure_sizes$Size_no_prec[grep(node, pure_sizes$Node)])
  #   pure_sizes$CProgeny_count[n] <- sum(pure_sizes$Progeny_count[grep(node, pure_sizes$Node)])
  # }
  # 
  # sample_tree_list[[length(sample_tree_list)]]$Pure_sizes <- pure_sizes
  # sample_tree_list[[length(sample_tree_list)]]$Pure_counts <- pure_counts
}

# Testing influence of trees
# full_s_tree_list <- sample_tree_list
# sample_tree_list <- full_s_tree_list[4]

# Step III: Perform sampling and calculate correlations
samples <- 100
sampled_wcors <- data.frame(matrix(nrow = samples, ncol = nrow(precursor_ranking) - 1))
colnames(sampled_wcors) <- precursor_ranking$Precursor[-nrow(precursor_ranking)]

set.seed(42)
bootstrap_list <- list()
for(p in 1:(nrow(precursor_ranking) - 1)){
  bootstrap_precursor <- as.character(precursor_ranking$Precursor[p])
  print(paste("Potential precursor:", bootstrap_precursor))
  
  # Step IIIa: Randomly place cells in tree; calculate resulting cumulative precursor frequencies
  all_tree_sample_frequencies <- matrix(nrow = 0, ncol = samples)
  all_tree_progeny_frequencies <- matrix(nrow = 0, ncol = samples)
  all_tree_sampled_node_sizes <- matrix(nrow = 0, ncol = samples)
  for(t in 1:length(sample_tree_list)){
    pure_counts <- sample_tree_list[[t]]$Pure_counts
    if(#!(bootstrap_precursor %in% pure_counts$Cell_type) |
       !(target %in% pure_counts$Cell_type)){
      next
    }
    
    # pure_sizes_needed <- pure_sizes[, c("Node", "Precursor_count", "Sample_chance", "CSize_no_prec", "CProgeny_count")]
    pure_sizes <- aggregate(pure_counts$Type_count,
                            by = list(Node = pure_counts$Node),
                            sum)
    colnames(pure_sizes)[2] <- "Node_count"
    pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == target, ])
    colnames(pure_sizes)[3:4] <- c("Progeny", "Progeny_count")
    
    # Step I: potential precursor and how many cells?
    bootstrap_cells <- sum(pure_counts$Type_count[pure_counts$Cell_type == bootstrap_precursor])
    
    if(bootstrap_cells == 0){
      # sample_counts_cumulative <- matrix(0, ncol = samples, nrow = nrow(pure_sizes))
      # rownames(sample_counts_cumulative) <- pure_sizes$Node
      # colnames(sample_counts_cumulative) <- 1:samples
      # pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == bootstrap_precursor, ])
      # colnames(pure_sizes)[5:6] <- c("Precursor", "Precursor_count")
      pure_sizes$Precursor <- bootstrap_precursor
      pure_sizes$Precursor_count <- 0
    }else{
      # Step II: Pure node sizes with and without precursor
      pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == bootstrap_precursor, ])
      colnames(pure_sizes)[5:6] <- c("Precursor", "Precursor_count")
    }
    pure_sizes$Sample_chance <- 
      (pure_sizes$Node_count - pure_sizes$Progeny_count)/(sum(pure_sizes$Node_count) - sum(pure_sizes$Progeny_count))
    pure_sizes$Size_no_prec <- pure_sizes$Node_count - pure_sizes$Precursor_count
    
    for(n in 1:nrow(pure_sizes)){
      node <- as.character(pure_sizes$Node[n])
      
      # pure_sizes$CNode_count[n] <- sum(pure_sizes$Node_count[grep(node, pure_sizes$Node)])
      pure_sizes$CSize_no_prec[n] <- sum(pure_sizes$Size_no_prec[grep(node, pure_sizes$Node)])
      pure_sizes$CProgeny_count[n] <- sum(pure_sizes$Progeny_count[grep(node, pure_sizes$Node)])
    }
    
    # pure_sizes <- sample_tree_list[[t]]$Pure_sizes
    
    # pure_sizes_needed <- pure_sizes[, c("Node", "Precursor_count", "Sample_chance", "CSize_no_prec", "CProgeny_count")]
    if(bootstrap_cells == 0){
      sample_counts <- matrix(0, ncol = samples, nrow = nrow(pure_sizes))
      rownames(sample_counts) <- pure_sizes$Node
      colnames(sample_counts) <- 1:samples
    }else{
      # Sample
      sample_outcomes <- data.frame(Sample = rep(1:samples, each = sum(pure_sizes$Precursor_count)),
                                    Node_sampled = sample(x = as.character(pure_sizes$Node), size = samples * sum(pure_sizes$Precursor_count), 
                                                          replace = T, prob = pure_sizes$Sample_chance))
      sample_outcomes$Node_sampled <- factor(sample_outcomes$Node_sampled, levels = as.character(pure_sizes$Node))
      # Put outcomes in table
      sample_counts <- acast(data.frame(table(sample_outcomes$Sample, sample_outcomes$Node_sampled)),
                             Var2 ~ Var1, value.var = "Freq")
    }
    # Calculate cumulative sampled values
    sample_counts_cumulative <- sample_counts
    
      for(n in 1:nrow(sample_counts_cumulative)){
        node <- as.character(rownames(sample_counts_cumulative)[n])
        
        sample_counts_cumulative[n, ] <- colSums(sample_counts[grep(node, rownames(sample_counts)), , drop = F])
      }
    
    # Calculate node sizes with sampled precursor
    rownames(pure_sizes) <- pure_sizes$Node
    new_node_sizes <- 
      sample_counts_cumulative + 
      pure_sizes[rownames(sample_counts_cumulative), ]$CSize_no_prec
    rownames(new_node_sizes) <- paste(names(sample_tree_list)[t], rownames(new_node_sizes), sep = ":")
    # Calculate resulting sampled frequencies
    sample_frequencies <- sample_counts_cumulative/new_node_sizes
    rownames(sample_frequencies) <- paste(names(sample_tree_list)[t], rownames(sample_frequencies), sep = ":")
    # Calculate resulting progeny frequencies
    progeny_frequencies <- pure_sizes$CProgeny_count/new_node_sizes
    # rownames(progeny_frequencies) <- paste(names(sample_tree_list)[t], rownames(progeny_frequencies), sep = ":")
    
    # Bind to samples from other trees
    all_tree_sampled_node_sizes <- rbind(all_tree_sampled_node_sizes, 
                                         new_node_sizes)
    all_tree_sample_frequencies <- rbind(all_tree_sample_frequencies, sample_frequencies)
    all_tree_progeny_frequencies <- rbind(all_tree_progeny_frequencies, progeny_frequencies)
  }
  
  bootstrap_list[[p]] <- list(Node_sizes = all_tree_sampled_node_sizes,
                              Precursor_sampled_freqs = all_tree_sample_frequencies,
                              Progeny_frequencies = all_tree_progeny_frequencies)
  names(bootstrap_list)[p] <- bootstrap_precursor
  
  # Step IIIb: Calculate weighted progeny-node correlation
  for(s in 1:samples){
    sampled_wcors[s, p] <- 
      wtd.cors(all_tree_sample_frequencies[(all_tree_progeny_frequencies[, s] != 0), s], 
               all_tree_progeny_frequencies[(all_tree_progeny_frequencies[, s] != 0), s], 
               weight = all_tree_sampled_node_sizes[(all_tree_progeny_frequencies[, s] != 0), s])
    # x <- data.frame(Prec = all_tree_sample_frequencies[, s],
    #                 Prog = all_tree_progeny_frequencies[, s],
    #                 Weight = all_tree_sampled_node_sizes[, s])
    # ggplot(x[x$Prog != 0, ]) +
    #   geom_point(aes(x = Prec, y = Prog, size = Weight))
    # wtd.cors(x$Prec[x$Prog != 0], x$Prog[x$Prog != 0], x$Weight[x$Prog != 0])
    # 
    # y <- progpos[, c(target, precursor)]
    # colnames(y) <- c("Prog", "Prec")
    # ggplot(y) +
    #   geom_point(aes(x = Prec, y = Prog))
    # cor(y$Prec, y$Prog)
    # wtd.cors(all_tree_sample_frequencies[, s], 
      #          all_tree_progeny_frequencies[, s], 
      #          weight = all_tree_sampled_node_sizes[, s])
    
    # all_tree_sample_frequencies[(all_tree_sample_frequencies[, s] != 0), s]
    
  }
}

precursor_ranking$Cor_bootstrap_p <- NA
precursor_ranking$Z_bootstrap <- NA
for(p in 1:(nrow(precursor_ranking) - 1)){
  # Calculate how many samples have a higher correlation than the measured one
  precursor_ranking$Cor_bootstrap_p[p] <-
    sum(sampled_wcors[, p] > precursor_ranking$Weighted_cor_progpos[p])/samples
  # Calculate z-distance to mean correlation
  precursor_ranking$Z_bootstrap[p] <-
    (precursor_ranking$Weighted_cor_progpos[p] - mean(sampled_wcors[, p]))/sd(sampled_wcors[, p])
}


data_cors <- comparison_list$Normalized_frequencies[, c(target, "M (apoeb)")]
colnames(data_cors) <- c("F_nppc", "M_apoeb")
dc <- data_cors[data_cors$F_nppc != 0,]
dc_b1 <- data.frame(F_nppc = bootstrap_list$`M (apoeb)`$Progeny_frequencies[, 1],
                    M_apoeb = bootstrap_list$`M (apoeb)`$Precursor_sampled_freqs[, 1])
dc_b1 <- dc_b1[dc_b1$F_nppc != 0, ]

ggplot(dc) +
  geom_point(aes(x = F_nppc, y = M_apoeb))
ggplot(dc_b1) +
  geom_point(aes(x = F_nppc, y = M_apoeb))

cor(dc)
cor(dc_b1)
cor(dc[grep("Hr24", rownames(dc)), ])
cor(dc_b1[grep("Hr24", rownames(dc_b1)), ])
cor(dc[grep("Hr26", rownames(dc)), ])
cor(dc_b1[grep("Hr26", rownames(dc_b1)), ])
cor(dc[grep("Hr27", rownames(dc)), ])
cor(dc_b1[grep("Hr27", rownames(dc_b1)), ])
cor(dc[grep("Hr1", rownames(dc)), ])
cor(dc_b1[grep("Hr1", rownames(dc_b1)), ])
cor(dc[grep("Hr2", rownames(dc)), ])
cor(dc_b1[grep("Hr2", rownames(dc_b1)), ])
cor(dc[grep("Hr14", rownames(dc)), ])
cor(dc_b1[grep("Hr14", rownames(dc_b1)), ])

dc_b1 <- data.frame(F_nppc = bootstrap_list$`M (apoeb)`$Progeny_frequencies[, 1],
  M_apoeb = bootstrap_list$`M (apoeb)`$Precursor_sampled_freqs[, 1])
dc_b1 <- dc_b1[dc_b1$F_nppc != 0, ]
# sampled_wcors_smcv2 <- sampled_wcors
# sampled_wcors_endoV <- sampled_wcors



# sampled_wcors_endfrzbV <- sampled_wcors
# ggplot(sampled_wcors_endoV) +
#   geom_histogram(aes(x = Wt_cor))
# mean(sampled_wcors_endoV$Wt_cor)
# sd(sampled_wcors_endoV$Wt_cor)
# 

# Analyse evidence against precursors ####
full_descendancy <- data.frame(Cell_type = character(),
                               Precursor = character(),
                               Node = character(),
                               Parent_node = character(),
                               Precursor_count = numeric(),
                               Total = numeric(),
                               Frequency = numeric(),
                               Precursor_presence_p = numeric(),
                               Entropy = numeric(),
                               Tree = character())
for(t in 1:length(tree_list)){
  tree_list[[t]] <- CalculateCooccurrence(tree_list[[t]])
  descendancy_add <- tree_list[[t]]$Descendancy
  descendancy_add$Tree <- names(tree_list)[t]
  descendancy_add <- merge(descendancy_add, tree_list[[t]]$Node_entropy)
  full_descendancy <- rbind(full_descendancy, descendancy_add)
}

full_descendancy_plot <- full_descendancy[full_descendancy$Cell_type == target &
                                            full_descendancy$Precursor %in% potential_prec, ]
full_descendancy_plot$Precursor <- factor(full_descendancy_plot$Precursor, levels = precursor_ranking$Precursor)
full_descendancy_plot$Log10p_bottom <- 
  sapply(full_descendancy_plot$Precursor_presence_p, function(x) max(log10(x), -10))

# png(paste("./Images/", target, "precursor_node_probs_7dpi.png", sep = ""),
    # width = 768, height = 256 * ceiling(length(potential_prec)/5))
print(
  # ggplot(full_descendancy[full_descendancy$Cell_type == zoom_to[i] &
  ggplot(full_descendancy_plot) +
    geom_jitter(aes(x = 0, y = Log10p_bottom, color = Precursor), height = 0) +
    # scale_y_continuous(limits = c(-10.1, 0.1)) + 
    # labs(title = paste(zoom_to[i], "precursors"),
    labs(title = paste(target, "precursors"),
         y = "Probability (log10)", x = "") +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    scale_color_manual(values = ann_colors$Celltype) +
    facet_wrap(~Precursor, ncol = 5) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 12, face = "bold"))
)
# dev.off()


