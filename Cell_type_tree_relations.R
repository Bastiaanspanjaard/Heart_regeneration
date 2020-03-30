# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)
require(weights)
require(data.table)

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

# library_name = "H5"
# reference_set
# load("./Data/Trees/A5_Ltree_pie.Robj")
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
                    width = 1000, height = 500,
                    ctypes = type_colors$celltype,linkLength=60, 
                    ct_colors = ct_colors, angle = pi/2,fontSize = 0,
                    nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 

  names(tree_sample)[which(names(tree_sample) == "CTree")] <- pie_tree_name
  
  return(tree_sample)
}

# Prepare colors and cell type data ####
# type_colors <- read.csv("./Data/color.table.all.sub.csv", sep = ";", stringsAsFactors = F)[, -1]

# cell_types <- read.csv("./Data/all.hearts.all.cells.all.sub.sept03.csv", stringsAsFactors = F)
cell_types <- fread("./Data/celltypes_zoom_allcells.csv",
                    header = T)
# cell_types$Cell_type <- cell_types$zoom.subtypes #cell_types$immune.fibro.CM.subtypes
cell_types$Cell_name <- paste("nd", cell_types$V1, sep = "")

type_rename <- data.frame(zoom.subtypes = unique(cell_types$zoom.subtypes))
type_rename$Cell_type <- as.character(type_rename$zoom.subtypes)
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("CM Atrium", "CM Ventricle (ttn.2)", "CM Ventricle", "CM Atrium (ttn.2)",
                          "CM (Proliferating)")] <- "Cardiomyocytes"
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("Bl.ves.EC (apnln)", "Bl.ves.EC (plvapb)", 
                          "Bl.ves.EC (lyve1)")] <- "Blood vessel EC"
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("Fibroblast (acta2)")] <- "Fibroblast"
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("Ery Duplex")] <- NA
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("Immune Duplex")] <- NA
type_rename$Cell_type[type_rename$zoom.subtypes %in% 
                        c("Fibroblast (proliferating)")] <- "Fibroblast (col11a1a)"
type_rename <- type_rename[complete.cases(type_rename), ]
cell_types <- merge(cell_types, type_rename, by = "zoom.subtypes")

type_colors <- data.frame(celltype = unique(type_rename$Cell_type))
type_colors$Color <-
  c("#B79F00", "#7CAE00", "#e41a1c", "#00B4F0",
    "#f781bf", "#FF64B0", "#00BA38", "#b2df8a",
    "#F8766D", "#ff7f00", "#4daf4a", "#b15928",
    "#F564E3", "#bebada", "#377eb8", "#984ea3",
    "#C77CFF", "#b2b223", "#ffff33", "#DE8C00")
# for(clr in c('Epicardium A', 'Epicardium V', 'Fibroblast', 'Fibroblast (cfd)',
#     'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)',
#     'Fibroblast (cxcl12a)', 'Perivascular cells')){
#   print(type_colors$Color[type_colors$Cell_type == clr])
# }
# type_colors <- 
#   type_colors[!(type_colors$celltype %in% 
#                   c("CM (Ery Duplex)", "CM (Immune Duplex)", "Dead cells", "Ery. Duplex",
#                     "M. duplex")), ]
# type_colors$Cell_type <-
#   c("Immune Cells", "Blood vessel EC", "Blood vessel EC", "Blood vessel EC",
#     "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", 
#     "Cardiomyocytes", "Cardiomyocytes", "Cardiomyocytes", "Cardiac neurons",
#     "Endocardium", "Endocardium", "Endocardium", "Endocardium",
#     "Endocardium", "Endocardium", )
# type_colors$celltype


# ph_annotation <- data.frame(Celltype = type_colors$celltype, row.names = type_colors$celltype)
# ph_zoom_annotation <- data.frame(Zoomtype = type_colors$celltype, row.names = type_colors$celltype)
ann_colors <-
  list(Celltype = setNames(type_colors$Color, type_colors$Cell_type))


# ann_colors <-
#   list(Celltype = setNames(type_colors$Color2, type_colors$celltype),
#        Zoomtype = setNames(type_colors$colo1, type_colors$celltype))

# zoom_types <- setdiff(type_colors$celltype[type_colors$colo1 != ""],
#                       c("Ery. Duplex", "M. Duplex"))
# vertex_colors <- data.frame(celltype = zoom_types,
#                             Order = 1:(length(zoom_types)))
# vertex_colors <- merge(vertex_colors, type_colors[, c("celltype", "colo1")])
# vertex_colors <- vertex_colors[order(vertex_colors$Order), ]
# rownames(vertex_colors) <- vertex_colors$celltype

# target <- "Fibroblast (nppc)"
# target <- "Fibroblast (col11a1a)"
target <- "Fibroblast-like cells"
# target <- "Perivascular cells"
# target <- "Fibroblast (proliferating)"
# target <- "Fibroblast (col12a1a)"
# target <- "Ventricle (ttn.2/aSMA)"
precursors_include_incl_target <- #NULL
  # c("Fibroblast (cxcl12a)", "Fibroblast", "Epicardium V",
  #   "Perivascular cells",
  #   "Fibroblast (col11a1a)", "Epicardium A", "Fibroblast (col12a1a)")
  c("Fibroblast-like cells", "Endocardium", "Fibroblast (spock3)")
precursors_include <- setdiff(precursors_include_incl_target, target)

# Load tree objects, append cell types, create lineage trees ####
# Control
H5 <- ReadTree("H5", reference_set = cell_types[cell_types$orig.ident == "H5", ])
H6 <- ReadTree("H6", reference_set = cell_types[cell_types$orig.ident == "H6", ])
H7 <- ReadTree("H7", reference_set = cell_types[cell_types$orig.ident == "H7", ])
# tree_list <- list(H5 = H5, H6 = H6, H7 = H7)
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
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Full_tree", ct_colors = type_colors$Color)
  # tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
  #                               ct_colors = type_colors$colo1)
  htmlwidgets::saveWidget(
    tree_list[[t]]$Full_tree,
    file = paste("~/Documents/Projects/heart_Bo/Images/tree_",
                 names(tree_list)[t], "_LINNAEUS_pie.html", sep = ""))
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

analysis_stats <-
  data.frame(Tree = character(),
           Included = logical(),
           Cell_type = character(),
           Count = integer())

for(t in 1:length(tree_list)){
  tree <- tree_list[[t]]$Tree
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  analysis_stats_add <-
    data.frame(table(edge_list$Cell.type[edge_list$Cell.type != "NA"]))
  colnames(analysis_stats_add) <- c("Cell_type", "Count")
  analysis_stats_add$Tree <- names(tree_list)[t]
  if(!(target %in% edge_list$Cell.type) | sum(edge_list$Cell.type == target) < 20){
    analysis_stats_add$Included <- F
    analysis_stats <- rbind(analysis_stats, analysis_stats_add)
    next
  }
  analysis_stats_add$Included <- T
  analysis_stats <- rbind(analysis_stats, analysis_stats_add)
  
  sample_type_count <- 
    reshape2::dcast(data.frame(table(edge_list$from, edge_list$Cell.type)), Var1 ~ Var2, value.var = "Freq")
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
             Weighted_cor_progpos = numeric(ncol(comparison_list$Normalized_frequencies)),
             # Weighted_boot_cor = numeric(ncol(comparison_list$Normalized_frequencies)))#,
             Weighted_cor_se = numeric(ncol(comparison_list$Normalized_frequencies)))#,
             #Distance = numeric(ncol(comparison_list$Normalized_frequencies)))
comparison_list$All_trees_distances[, unique(analysis_stats$Tree[analysis_stats$Included])] <- 0

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
  weighted_cor <- data.frame(wtd.cor(progpos[[precursor]],
                                     progpos[[target]],
                                     weight = progpos_weights$Size))#, bootse = T))
  comparison_list$All_trees_distances$Weighted_cor_progpos[p] <- weighted_cor$correlation
  # comparison_list$All_trees_distances$Weighted_boot_cor[p] <- weighted_cor$bootcor
  comparison_list$All_trees_distances$Weighted_cor_se[p] <- weighted_cor$std.err
  
  for(t in unique(analysis_stats$Tree[analysis_stats$Included])){
    progpos_oneout <- progpos[!grepl(paste(t, ":", sep = ""), rownames(progpos)), ]
    weights_oneout <- progpos_weights[!grepl(paste(t, ":", sep = ""), rownames(progpos_weights)), ]
    comparison_list$All_trees_distances[p, t] <-
      wtd.cors(progpos_oneout[[precursor]],
             progpos_oneout[[target]],
             weight = weights_oneout$Size)
    }
  
  # comparison_list$All_trees_distances$Weighted_cor_progpos[p] <-
  #   wtd.cors(progpos[[precursor]],
  #            progpos[[target]],
  #            weight = progpos_weights$Size)
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
# print(
#   ggplot(precursor_ranking[precursor_ranking$Precursor != target, ]) +
#     geom_bar(aes(x = Precursor, y = Weighted_cor_progpos, fill = Precursor, 
#                  alpha = ifelse(Mean_count < 10, 0.5, 1)), stat = "identity") +
#     scale_fill_manual(values = ann_colors$Celltype) +
#     labs(title = paste(target, ""), y = "Fitness") +
#     theme(legend.position = "none",
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank())
# )
# dev.off()

# Plot correlations of two top-scoring precursors and two low-scoring precursors
i <- 2
for(i in c(1:2, nrow(precursor_ranking) - 1, nrow(precursor_ranking))){
  pot_prec <- as.character(precursor_ranking$Precursor[i])
  plot_freqs <-
    data.frame(comparison_list$Normalized_frequencies)[,make.names(c(pot_prec, target))]
  colnames(plot_freqs) <- c("Precursor", "Progenitor")
  plot_freqs$Size <- comparison_list$Node_sizes$Size
  # png(paste("./Images/", pot_prec, "precursor_to_", target, "_37dpi_all.png", sep = ""))
  print(
    ggplot() +
      geom_point(data = plot_freqs[plot_freqs$Progenitor == 0, ],
                 aes(x = Precursor, y = Progenitor), size = 6) + #, size = Size)) + #, color = "grey") +
      geom_point(data = plot_freqs[plot_freqs$Progenitor > 0, ],
                 aes(x = Precursor, y = Progenitor), size = 6) + #, size = Size)) +
      labs(x = pot_prec, y = target) +
      theme(text = element_text(size = 36))#+
      # theme(legend.position = "none")
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
samples <- 1000
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
       !(target %in% pure_counts$Cell_type) | sum(pure_counts$Type_count[pure_counts$Cell_type == target]) < 20){
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

example_sample <- sampled_wcors[, 2, drop = F]
colnames(example_sample) <- "Wcor"

# Example plot
# png("./Images/Ventricular_endocardium_to_nppc_fibroblast_bootcor.png",
#     width = 600, height = 400)
ggplot(example_sample) +
  geom_histogram(aes(x = Wcor)) +
  geom_vline(xintercept = mean(example_sample$Wcor), col = "red", size = 3) +
  geom_vline(xintercept = 
               precursor_ranking$Weighted_cor_progpos[precursor_ranking$Precursor == "Endocardium (Ventricle)"],
             col = "blue", size = 3) +
  labs(x = "Correlation", y = "Count") +
  theme(text = element_text(size = 36))
# dev.off()

precursor_ranking$Wcor_plus095 <- precursor_ranking$Weighted_cor_progpos + 2 * precursor_ranking$Weighted_cor_se
precursor_ranking$Wcor_minus095 <- precursor_ranking$Weighted_cor_progpos - 2 * precursor_ranking$Weighted_cor_se

precursor_ranking$Cor_bootstrap_p <- NA
precursor_ranking$Z_bootstrap <- NA
precursor_ranking$Mean_cor <- NA
precursor_ranking$Plus095 <- NA
precursor_ranking$Minus095 <- NA
for(p in 1:(nrow(precursor_ranking) - 1)){
  # Average bootstrap correlation
  precursor_ranking$Mean_cor[p] <- mean(sampled_wcors[, p])
  higher_cors <- sampled_wcors[sampled_wcors[, p] > mean(sampled_wcors[, p]), p]
  precursor_ranking$Plus095[p] <- quantile(sort(higher_cors), probs = 0.95)
  lower_cors <- sampled_wcors[sampled_wcors[, p] < mean(sampled_wcors[, p]), p]
  precursor_ranking$Minus095[p] <- quantile(sort(lower_cors), probs = 0.05)
  # Calculate how many samples have a higher correlation than the measured one
  precursor_ranking$Cor_bootstrap_p[p] <-
    sum(sampled_wcors[, p] > precursor_ranking$Weighted_cor_progpos[p])/samples
  # Calculate z-distance to mean correlation
  precursor_ranking$Z_bootstrap[p] <-
    (precursor_ranking$Weighted_cor_progpos[p] - mean(sampled_wcors[, p]))/sd(sampled_wcors[, p])
}
# precursor_ranking_plot <- precursor_ranking[(precursor_ranking$Weighted_cor_progpos > 0 &
#                                               precursor_ranking$Cor_bootstrap_p < 0.001) |
#                                               precursor_ranking$Precursor %in% precursors_include, ]
potential_prec <-
  union(precursor_ranking$Precursor[precursor_ranking$Cor_bootstrap_p < 0.001 &
                                      precursor_ranking$Weighted_cor_progpos > 0],
        precursors_include)
# potential_prec <- precursor_ranking$Precursor[precursor_ranking$Precursor != target]
# potential_prec <- potential_prec[grep("Duplex", potential_prec, invert = T)]

# OLD
# precursor_ranking_plot <- precursor_ranking[(precursor_ranking$Weighted_cor_progpos > 0 &
#                                                precursor_ranking$Wcor_minus095 > precursor_ranking$Plus095) |
#                                               precursor_ranking$Precursor %in% precursors_include, ]
precursor_ranking_plot <- precursor_ranking[precursor_ranking$Precursor %in% potential_prec, ]
prp_m <- reshape2::melt(precursor_ranking_plot[, c("Precursor", "Weighted_cor_progpos", "Mean_cor")])
prp_m$Type_alpha <- ifelse(prp_m$variable == "Weighted_cor_progpos", 1, 0.9)
prp_bebar <- precursor_ranking_plot[, c("Precursor", "Minus095", "Plus095")]
prp_m <- merge(prp_m, prp_bebar)
prp_cebar <- precursor_ranking_plot[, c("Precursor", "Wcor_minus095", "Wcor_plus095")]
prp_m <- merge(prp_m, prp_cebar)

prp_m$Eplus <- ifelse(prp_m$variable == "Weighted_cor_progpos",
                      prp_m$Wcor_plus095, prp_m$Plus095)
prp_m$Eminus <- ifelse(prp_m$variable == "Weighted_cor_progpos",
                      prp_m$Wcor_minus095, prp_m$Minus095)

# prp_m$Minus099[prp_m$variable == "Weighted_cor_progpos"] <- NA
# prp_m$Plus099[prp_m$variable == "Weighted_cor_progpos"] <- NA
prp_m$Eminus[prp_m$variable == "Weighted_cor_progpos"] <- NA
prp_m$Eplus[prp_m$variable == "Weighted_cor_progpos"] <- NA


# pdf(paste("./Images/", target, "_all_precursors_with_bootstrap_37dpi.pdf", sep = ""),
    #width = 1400, height = 768)
png(paste("./Images/", target, "_all_precursors_with_bootstrap_37dpi_legend.png", sep = ""),
    type = "quartz", width = 900, height = 768)
# png("./Images/Fibroblast_col11_sources.png",
    # width = 960, height = 480)
print(
  ggplot(prp_m) +
    geom_bar(aes(x = Precursor, y = value, fill = Precursor, group = variable, alpha = Type_alpha),
             stat = "identity", position = "dodge") +
    # geom_errorbar(aes(x = Precursor, ymin = Minus099, ymax = Plus099, group = variable), 
    #               position = position_dodge(width = 1), width = 0.5) +
    geom_errorbar(aes(x = Precursor, ymin = Eminus, ymax = Eplus, group = variable), 
                  position = position_dodge(width = 1), width = 0.5) +#, size = 2) +
    scale_fill_manual(values = ann_colors$Celltype) +
    scale_alpha(range = c(0.35, 1), guide = F) +
    labs(title = paste(target, ""), y = "Correlation", fill = "Source") +
    theme(#legend.position = "none",
          text = element_text(size = 36),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
)
dev.off()

# potential_prec <- 
#   union(precursor_ranking$Precursor[precursor_ranking$Cor_bootstrap_p < 0.001 &
#                                       precursor_ranking$Weighted_cor_progpos > 0],
#         precursors_include)
# union(precursor_ranking$Precursor[precursor_ranking$Weighted_cor_progpos > 0.5 & 
#                               precursor_ranking$Mean_count >= 10],
#       precursors_include)
# print(potential_prec)

# precursor_ranking_plot <- precursor_ranking_plot[order(-precursor_ranking_plot$Z_bootstrap), ]
# precursor_ranking_plot$Precursor <- factor(precursor_ranking_plot$Precursor, 
#                                            levels = precursor_ranking_plot$Precursor)

# ggplot(precursor_ranking_plot) +
#   geom_bar(aes(x = Precursor, y = Z_bootstrap, fill = Precursor,
#                alpha = 1 - Cor_bootstrap_p), stat = "identity") +
#   scale_fill_manual(values = ann_colors$Celltype) +
#   labs(title = paste(target, ""), y = "Fitness") +
#   theme(legend.position = "none",
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())

# print(
#   ggplot(precursor_ranking[precursor_ranking$Precursor != target, ]) +
#     geom_bar(aes(x = Precursor, y = Weighted_cor_progpos, fill = Precursor,
#                  alpha = ifelse(Mean_count < 10, 0.5, 1)), stat = "identity") +
#     scale_fill_manual(values = ann_colors$Celltype) +
#     labs(title = paste(target, ""), y = "Fitness") +
#     theme(legend.position = "none",
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank())
# )


# make.names("M (apoeb)")
# precursor_candidate <- "Fibroblasts"
# # 
# data_cors <- comparison_list$Normalized_frequencies[, c(target, precursor_candidate)]
# colnames(data_cors) <- c("Target", "Precursor")
# dc <- data_cors[data_cors$Target != 0,]
# dc_b1 <- data.frame(Target = bootstrap_list[[precursor_candidate]]$Progeny_frequencies[, 2],
#                     Precursor = bootstrap_list[[precursor_candidate]]$Precursor_sampled_freqs[, 2])
# dc_b1 <- dc_b1[dc_b1$Target != 0, ]

# ggplot(dc) +
#   geom_point(aes(x = F_nppc, y = M_apoeb))
# ggplot(dc_b1) +
#   geom_point(aes(x = F_nppc, y = M_apoeb))
# 
# for(t in 1:length(sample_tree_list)){
#   counts <- sample_tree_list[[t]]$Pure_counts
#   print(paste(names(sample_tree_list)[t], "-", target, ":", sum(counts$Type_count[counts$Cell_type == target]),
#               precursor_candidate, ":", sum(counts$Type_count[counts$Cell_type== precursor_candidate])))
# }

# cor(dc)
# cor(dc_b1)
# cor(dc[grep("Hr24", rownames(dc)), ])
# cor(dc_b1[grep("Hr24", rownames(dc_b1)), ])
# cor(dc[grep("Hr26", rownames(dc)), ])
# cor(dc_b1[grep("Hr26", rownames(dc_b1)), ])
# cor(dc[grep("Hr27", rownames(dc)), ])
# cor(dc_b1[grep("Hr27", rownames(dc_b1)), ])
# cor(dc[grep("Hr1", rownames(dc)), ])
# cor(dc_b1[grep("Hr1", rownames(dc_b1)), ])
# cor(dc[grep("Hr2", rownames(dc)), ])
# cor(dc_b1[grep("Hr2", rownames(dc_b1)), ])
# cor(dc[grep("Hr14", rownames(dc)), ])
# cor(dc_b1[grep("Hr14", rownames(dc_b1)), ])

# dc_n_b1 <- cbind(dc, dc_b1)[, c(1, 2, 4)]
# dc_n_b1$Tree <- sapply(rownames(dc_n_b1), function(x){unlist(strsplit(x, ":"))[1]})
# ggplot(dc_n_b1) +
#   geom_point(aes(x = Target, y = Precursor, color = Tree)) +
#   geom_segment(aes(x = Target, xend = Target, y = Precursor, yend = Precursor.1, color = Tree),
#                arrow = arrow(length=unit(0.30,"cm"))) +
#   labs(x = target, y = precursor_candidate)
# for(t in unique(dc_n_b1$Tree)){
#   print(
#     paste("Tree", t, "observed correlation:",
#           cor(dc_n_b1$Target[dc_n_b1$Tree == t], dc_n_b1$Precursor[dc_n_b1$Tree == t]),
#           "bootstrapped correlation:",
#           cor(dc_n_b1$Target[dc_n_b1$Tree == t], dc_n_b1$Precursor.1[dc_n_b1$Tree == t]))
#   )        
# }
# 
# tree_to_look_at <- "Hr26"
# 
# precursor_candidate <- "Fibroblasts"
# bootstrapping_cors <-
#   data.frame(Precursor = rep(precursor_candidate, samples),
#              Sample = 1:samples,
#              Bootstrapped_cor = numeric(1000))
# data_cors <- data.frame(comparison_list$Normalized_frequencies[, c(target, precursor_candidate)],
#                         comparison_list$Node_sizes)
# colnames(data_cors)[1:2] <- c("Target", "Precursor")
# dc <- data_cors[data_cors$Target != 0 & grepl(paste(tree_to_look_at, ":", sep = ""), rownames(data_cors)),]
# 
# for(b in 1:samples){
#   dc_boot <- data.frame(Node_sizes = bootstrap_list[[precursor_candidate]]$Node_sizes[, b],
#                         Prec_freq = bootstrap_list[[precursor_candidate]]$Precursor_sampled_freqs[, b],
#                         Prog_freq = bootstrap_list[[precursor_candidate]]$Progeny_frequencies[, b])
#   dc_boot <- dc_boot[dc_boot$Prog_freq != 0 & grepl(paste(tree_to_look_at, ":", sep = ""), rownames(dc_boot)), ]
#   bootstrapping_cors$Bootstrapped_cor[b] <-
#     wtd.cors(dc_boot$Prec_freq, dc_boot$Prog_freq, weight = dc_boot$Node_sizes)
# }
# # 
# # png(paste("./Images/", tree_to_look_at, "_", precursor_candidate, "_to_", target, ".png", sep = ""),
# #     width = 738)
# ggplot(bootstrapping_cors) +
#   geom_histogram(aes(x = Bootstrapped_cor), binwidth = 0.05) +
#   geom_vline(xintercept = mean(bootstrapping_cors$Bootstrapped_cor), color = "red") +
#   geom_vline(xintercept = wtd.cors(dc$Target, dc$Precursor, weight = dc$Size), color = "blue") +
#   labs(title = paste(tree_to_look_at, ":", precursor_candidate, "to", target),
#        x = "Bootstrapped correlation")
# # dev.off()
# labs(title = paste("All trees:", precursor_candidate, "to", target))

# 
# dc_b1 <- data.frame(F_nppc = bootstrap_list$`M (apoeb)`$Progeny_frequencies[, 1],
#   M_apoeb = bootstrap_list$`M (apoeb)`$Precursor_sampled_freqs[, 1])
# dc_b1 <- dc_b1[dc_b1$F_nppc != 0, ]
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
  edge_list <- ToDataFrameNetwork(tree_list[[t]]$Tree, "Cell.type")
  if(!(target %in% edge_list$Cell.type) | sum(edge_list$Cell.type == target) < 20){
    next
  }
  
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

x <- length(potential_prec)
# png(paste("./Images/", target, "precursor_node_probs_37dpi.png", sep = ""),
# width = 768, height = 256 * ceiling(length(potential_prec)/5))
# png(paste("./Images/", target, "precursor_node_probs_37dpi.png", sep = ""),
#     width = max(512, 384 * ceiling(x/5)), height = 256 * min(x, min(3, x)))
print(
  # ggplot(full_descendancy[full_descendancy$Cell_type == zoom_to[i] &
  ggplot(full_descendancy_plot) +
    geom_jitter(aes(x = 0, y = Log10p_bottom, color = Precursor), height = 0, size = 3) +
    # scale_y_continuous(limits = c(-10.1, 0.1)) + 
    # labs(title = paste(zoom_to[i], "precursors"),
    labs(title = paste(target, "precursors"),
         y = "Probability (log10)", x = "") +
    theme(axis.text.x = element_text(angle = 90),
          legend.position = "none") +
    scale_color_manual(values = ann_colors$Celltype) +
    facet_wrap(~Precursor, nrow = 3) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text.x = element_text(size = 12, face = "bold"))
)
# dev.off()



# Calculate distance between cell types in a tree ####
t <- 1
tree <- tree_list[[t]]$Tree
# ** Compute background tree distance distribution ####
# Calculate node distances
edge_list <- ToDataFrameNetwork(tree, "Cell.type")
nodes <- setdiff(unique(c(edge_list$from, edge_list$to[edge_list$Cell.type == "NA"])), "nd0")
node_distances <- data.frame(t(combn(nodes, 2)))
colnames(node_distances) <- c("Node_1", "Node_2")
node_distances$Distance <- 
  apply(node_distances[, c("Node_1", "Node_2")], 1,
        function(x){
          v_1 <- unlist(strsplit(as.character(x[1]), "_"))
          v_2 <- unlist(strsplit(as.character(x[2]), "_"))
          # For comparison, add "-1"'s to the shortest vector
          v_1f <- c(v_1, rep(-1, max(length(v_2) - length(v_1), 0)))
          v_2f <- c(v_2, rep(-1, max(length(v_1) - length(v_2), 0)))
          # Distance = d(str_1 - mcra) + d(str_2 - mcra) =
          #            d(str_1) - d(mcra) + d(str_2) - d(mcra) =
          #            d(str_1) + d(str_2) - 2*d(mcra)
          return(length(v_1) + length(v_2) - 2*sum(cumprod(v_1f == v_2f)))
        })
node_distances <- rbind(node_distances,
                        data.frame(Node_1 = nodes,
                                   Node_2 = nodes,
                                   Distance = 0))
# Calculate node chances
pure_sizes <- aggregate(pure_counts$Type_count,
                        by = list(Node = pure_counts$Node),
                        sum)
colnames(pure_sizes)[2] <- "Node_count"
pure_sizes <- pure_sizes[pure_sizes$Node != "nd0", ]
pure_sizes$Chance <- pure_sizes$Node_count/sum(pure_sizes$Node_count)
# Calculate distance chances
node_distances <- merge(node_distances, pure_sizes[, c("Node", "Chance")],
                        by.x = "Node_1", by.y = "Node")
colnames(node_distances)[4] <- "Chance_1"
node_distances <- merge(node_distances, pure_sizes[, c("Node", "Chance")],
                        by.x = "Node_2", by.y = "Node")
colnames(node_distances)[5] <- "Chance_2"
node_distances$Distance_chance <- node_distances$Chance_1 * node_distances$Chance_2
node_distances$Distance_chance[node_distances$Node_1 != node_distances$Node_2] <-
  2*node_distances$Distance_chance[node_distances$Node_1 != node_distances$Node_2]
# Calculate distance distribution and average
distance_distribution <- 
  aggregate(node_distances$Distance_chance,
          by = list(Distance = node_distances$Distance),
          sum)
colnames(distance_distribution)[2] <- "Chance"
ggplot(distance_distribution) +
  geom_bar(aes(x = Distance, y = Chance), stat = "identity") +
  labs(title = paste("Distance distribution for tree ", names(tree_list)[t],
                     ", average ", 
                     round(weighted.mean(distance_distribution$Distance, distance_distribution$Chance), 2),
                     sep = ""))
       
# ** Calculate distance distribution for cell type pairs ####
cell_type_distances <-
  data.frame(t(combn(setdiff(unique(edge_list$Cell.type), "NA"), 2)))
colnames(cell_type_distances) <- c("Type_1", "Type_2")
cell_type_distances <-
  rbind(cell_type_distances,
        data.frame(Type_1 = setdiff(unique(edge_list$Cell.type), "NA"),
           Type_2 = setdiff(unique(edge_list$Cell.type), "NA")))
cell_type_distances$Distance <- -1
cell_type_distances$Measurements <- 0
# d <- 1126
for(d in 1:nrow(cell_type_distances)){
  print(paste(d, " out of ", nrow(cell_type_distances), sep = ""))
  type_1 <- cell_type_distances$Type_1[d]
  type_2 <- cell_type_distances$Type_2[d]
  # Calculate distances - node 1, node 2, distance. Note the number of comparisons is m*n for
  # two different cell types, but 0.5*m(m-1) for inter-cell type distances.
  if(type_1 == type_2){
    type_positions <- edge_list[edge_list$Cell.type == type_1 &
                                    edge_list$from != "nd0", ]
    colnames(type_positions)[1:2] <- c("Node", "Cell")
    if(nrow(type_positions) > 1){
      type_distances <- data.frame(t(combn(type_positions$Cell, 2)))
      colnames(type_distances) <- c("Cell_1", "Cell_2")
      type_distances <- merge(type_distances, type_positions[, c("Node", "Cell")],
                              by.x = "Cell_1", by.y = "Cell")
      colnames(type_distances)[3] <- "Node_1"
      type_distances <- merge(type_distances, type_positions[, c("Node", "Cell")],
                              by.x = "Cell_2", by.y = "Cell")
      colnames(type_distances)[4] <- "Node_2"
      type_distances <- merge(type_distances, node_distances[, c("Node_1", "Node_2", "Distance")])
      cell_type_distances$Distance[d] <- mean(type_distances$Distance)
      cell_type_distances$Measurements[d] <- nrow(type_distances)
    }
  }else{
    type_1_positions <- edge_list[edge_list$Cell.type == type_1 &
                                    edge_list$from != "nd0", ]
    colnames(type_1_positions) <- c("Node_1", "Cell_1", "Type_1")
    type_2_positions <- edge_list[edge_list$Cell.type == type_2 &
                                    edge_list$from != "nd0", ]
    colnames(type_2_positions) <- c("Node_2", "Cell_2", "Type_2")
    type_distances <- merge(type_1_positions, type_2_positions)
    type_distances <- merge(type_distances, node_distances[, c("Node_1", "Node_2", "Distance")])
    cell_type_distances$Distance[d] <- mean(type_distances$Distance)
    cell_type_distances$Measurements[d] <- nrow(type_distances)
  }
}  

# ** Compare distances for cell type pairs to background distribution by calculating z-values ####
mean_dist <- weighted.mean(distance_distribution$Distance, distance_distribution$Chance)
distance_distribution$x_m <- distance_distribution$Distance - mean_dist
dist_sd <- sqrt(sum(distance_distribution$Chance * distance_distribution$x_m^2)/nrow(distance_distribution))
cell_type_distances$Z_distance <- 
  sqrt(cell_type_distances$Measurements) * (cell_type_distances$Distance - mean_dist)/dist_sd
View(cell_type_distances[cell_type_distances$Type_1 == cell_type_distances$Type_2, ])

