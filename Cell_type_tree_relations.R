# Dependencies ####
source("./Scripts/HR_library.R")

# Parameters ####
tree_path <- "./Data/Trees/Sbcf/"
inclusion_limit <- 10 # How many cells are needed of a certain type to include a tree for
# that type.

# Prepare colors and cell type data ####
cell_annotations <- read.csv("./Data/final_metadata_Tmacromerged_2.csv", stringsAsFactors = F)
# Make table of existing cell types and frequencies
celltype_frequencies <- data.frame(table(cell_annotations$Cell_type), stringsAsFactors = F)
colnames(celltype_frequencies)[1] <- c("Cell_type")
celltype_frequencies$Cell_type <- as.character(celltype_frequencies$Cell_type)

# Read in cell type colors
celltype_colors_in <- read.csv("./Data/Cell_type_colors_2.csv")
celltype_colors <- setNames(celltype_colors_in$color, celltype_colors_in$Cell.type)
rm(celltype_colors_in)

# Prepare trees ####
# tree_list_in <- InitializeTreesFromDisk(tree_path, libraries, cell_annotations)
# saveRDS(tree_list_in, "./Data/Trees/Tree_list_oneEndo.rds")
tree_list_in <- readRDS("./Data/Trees/Tree_list_oneEndo.rds")

# Export all trees in a universal format
# full_edge_list <- data.frame(from = character(),
#                              to = character(),
#                              Cell.type = character(),
#                              heart = character())
# for(t in names(tree_list_in)){
#   edge_list_add <- (tree_list_in[[t]])$Edge_list
#   edge_list_add$heart <- t
#   full_edge_list <- rbind(full_edge_list, edge_list_add)
# Unfortunately, the ToNewick script does not work very well and it stops at the first large tree, H8.
# Writing a custom ToNewick function is absolutely possible but would take a while (estimated one day)
# and all the information is in the edgelists.
# cat(ToNewick((tree_list_in[[t]])$Tree),
#       file = paste("~/Documents/Projects/heart_Bo/Data/Trees/Sbcf/", t, "_Newick.txt"))
# }
# write.csv(full_edge_list, file = "~/Documents/Projects/heart_Bo/Data/Trees/All_sbcf_edgelists.csv",
#           quote = F, row.names = F)

# write.table(ToNewick(tree_list_in$H5$Tree), 
#             file = "~/Documents/Projects/heart_Bo/Data/Trees/Sbcf/H5_Newick.txt", sep = "",
#             quote = F, row.names = F)
# cat(ToNewick(tree_list_in$H5$Tree),
#     file = "~/Documents/Projects/heart_Bo/Data/Trees/Sbcf/H5_Newick.txt")
# test_newick_tree <- ToNewick(tree_list_in$H5$Tree)

# Create tree visualization and zoom visualization ####
#- only run if needed because this takes quite some time.
# t <- 1
# for(t in 1:length(tree_list_in)){
#   tree_list_in[[t]]$Pie_tree <-
#     collapsibleTree(df = tree_list_in[[t]]$Tree, root =tree_list_in[[t]]$Tree$scar, pieNode = T,
#                     pieSummary = T, collapsed = F,
#                     width = 1000, height = 500,
#                     ctypes = names(celltype_colors), linkLength=60,
#                     ct_colors = as.character(celltype_colors), angle = pi/2,fontSize = 0,
#                     nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6))
# htmlwidgets::saveWidget(
#   tree_list_in[[t]]$Pie_tree,
#   file = paste("~/Documents/Projects/heart_Bo/Images/Trees/tree_",
#                names(tree_list_in)[t], "_LINNAEUS_pie.html", sep = ""))
# }

# Determine node characteristics to annotate trees ####
t <- 7 # total size 25; 1: H5, 23: Hr25, 6: Hr2, 7: Hr3
this_tree <- tree_list_in[[t]]$metadata
print(this_tree)
this_tree_node_celltypes <- 
  data.table(tree_list_in[[t]]$Node_type_counts)
this_tree_node_celltypes_cumulative <- this_tree_node_celltypes
this_tree_node_celltypes_cumulative$Type_count <- NA
for(i in 1:nrow(this_tree_node_celltypes_cumulative)){
  node <- this_tree_node_celltypes_cumulative$Node[i]
  cell_type <- this_tree_node_celltypes_cumulative$Cell_type[i]
  
  this_tree_node_celltypes_cumulative$Type_count[i] <-
    sum(this_tree_node_celltypes$Type_count[grepl(node, this_tree_node_celltypes$Node) &
                                              this_tree_node_celltypes$Cell_type == cell_type])
}
this_tree_node_celltypes <- this_tree_node_celltypes[this_tree_node_celltypes$Type_count > 0, ]
this_tree_node_celltypes_cumulative <- 
  this_tree_node_celltypes_cumulative[this_tree_node_celltypes_cumulative$Type_count > 0, ]
# this_tree_type_sizes <- 
#   this_tree_node_celltypes[, .(Type_total = sum(Type_count)), by = list(Cell_type = Cell_type)]
# this_tree_node_celltypes_cumulative <-
#   merge(this_tree_node_celltypes_cumulative, this_tree_type_sizes)
this_tree_node_sizes_cumulative <- 
  this_tree_node_celltypes_cumulative[, .(Node_count = sum(Type_count)), by = list(Node = Node)]
this_tree_node_celltypes_cumulative <- merge(this_tree_node_celltypes_cumulative, this_tree_node_sizes_cumulative)

# Finding cell types linked by lineage ####
# Determine trees and cell types to include
timepoint <- "7"
cluster_number <- 7
trees_to_include <- 
  unlist(lapply(tree_list_in, 
                function(x){x$metadata$Name[x$metadata$dpi == timepoint]}))
types_force_include <- NULL
trees_to_include <- trees_to_include[!is.na(trees_to_include)]
tree_list <- tree_list_in[trees_to_include]
tree_type_counts <- celltype_frequencies[, "Cell_type", drop = F]
for(t in 1:length(tree_list)){
  tree_name <- names(tree_list)[t]
  tree_list[[t]]$Tree <- Clone(tree_list_in[[tree_name]]$Tree)
  # Cells of each type per tree: sum up over all nodes, output has columns cell type, count columns for each tree
  # tree_counts <- 
  this_tree_counts <- aggregate(tree_list[[t]]$Node_type_counts$Type_count,
            by = list(tree_list[[t]]$Node_type_counts$Cell_type),
            sum)
  colnames(this_tree_counts) <- c("Cell_type", tree_name)
  tree_type_counts <- merge(tree_type_counts, this_tree_counts, all = T)
}
tree_type_counts[is.na(tree_type_counts)] <- 0
included_types <- 
  tree_type_counts$Cell_type[apply(tree_type_counts[, -1], 1, function(x) sum(x >= inclusion_limit) > 0)]

# Calculate cell-type node frequencies and node sizes
comparison_list <- ExtractTreeFrequencies(tree_list)
baseline_comparison_list <- comparison_list

# Scatterplot of frequencies highlighting different trees
# cor_type_1 <- "Fibroblasts (nppc)"
# cor_type_2 <- "Endocardium (Ventricle)"
# df_ct <- merge(comparison_list$Normalized_frequencies, comparison_list$Node_sizes, by = 0)
# df_ct <- df_ct[, c("Row.names", cor_type_1, cor_type_2, "Size")]
# colnames(df_ct) <- c("Node", "CT1", "CT2", "Node_size")
# df_ct$Tree <- sapply(df_ct$Node,
#                      function(x) unlist(strsplit(x, ":"))[1])
# df_ct$Tree <- factor(df_ct$Tree, 
#                      levels = c("Hr1", "Hr2", "Hr6", "Hr7",
#                                 "Hr13", "Hr14", "Hr15"))
# # pdf("./Images/Fibnppc_endoV_7dpi_corplot_trees.pdf",
# #     width = 5, height = 4, useDingbats = F)
# ggplot(df_ct) +
#   geom_point(aes(x = CT1, y = CT2, color = Tree), size = 4) +#, size = Node_size)) +
#   labs(x = cor_type_1, y = cor_type_2) +
#   scale_color_brewer(palette = 'Set1') +
#   theme(panel.grid.minor = element_blank())#,
# # title = paste("7dpi node frequencies all trees, cor ", round(cor(df_ct$CT1, df_ct$CT2), 2), sep = ""))
# # dev.off()

# Scatterplot of frequencies including only nodes that are over inclusion limit.
# norm_freq <- comparison_list$Normalized_frequencies
# inc_exc <- comparison_list$Normalized_frequencies # Same structure as norm_freq
# inc_exc[,] <- F # Set all entries to F
# for(t in 1:length(tree_list)){
#   tree_name <- names(tree_list)[t]
#   tree_type_include <- aggregate(tree_list[[t]]$Node_type_counts$Type_count,
#                                  by = list(Cell_type = tree_list[[t]]$Node_type_counts$Cell_type),
#                                  sum) # Count cell type numbers in tree
#   tree_type_include$Include <- tree_type_include$x >= inclusion_limit # Set included celltypes to T in dataframe
#   inc_exc[grepl(tree_name, rownames(inc_exc)), tree_type_include$Cell_type] <- # Nodes in tree, cell types in tree (cells not in tree remain F)
#     matrix(rep(tree_type_include$Include, each = sum(grepl(tree_name, rownames(inc_exc)))),
#            nrow = sum(grepl(tree_name, rownames(inc_exc))),
#            ncol = nrow(tree_type_include)) # Create matrix of (nodes in tree) x (cell types in tree), with values T for cell types that
#   # have enough abundance.
# }
# x_plot_name <- "Epicardium (Ventricle)" #
# y_plot_name <- "Fibroblasts (col12a1a)" #"Fibroblast-like cells"
# x_plot_col <- which(colnames(norm_freq) == x_plot_name)
# y_plot_col <- which(colnames(norm_freq) == y_plot_name)
# nf_plot <- cbind(norm_freq[, c(x_plot_col, y_plot_col)], inc_exc[, c(x_plot_col, y_plot_col)])
# colnames(nf_plot) <- c("Xp", "Yp", "Inc_X", "Inc_Y")
# nf_plot <- nf_plot[nf_plot$Inc_X & nf_plot$Inc_Y, ]
# # pdf("./Images/Fibcol12_epiV_3dpi_corplot.pdf",
# #     width = 4, height = 4, useDingbats = F)
# ggplot(nf_plot) +
#   geom_point(aes(x = Xp, y = Yp), size = 4) +
#   labs(x = x_plot_name, y = y_plot_name) +
#   theme(text = element_text(size = 16),
#         panel.grid.minor = element_blank())
# # dev.off()

# Calculate correlations
type_tree_cor <- CalculateTypeCorrelations(tree_list, comparison_list, 
                                           force_include = types_force_include, inclusion_limit = inclusion_limit)
cor_present <- type_tree_cor[rownames(type_tree_cor) != "Dead cells", colnames(type_tree_cor) != "Dead cells"]
# Remove cell types with largest amount of NaNs consecutively; in case of ties, remove the cell
# type with the lowest number of cells.
while(T){
  type_keep_decider <- data.frame(Number_NaNs = apply(cor_present, 1, function(x) {sum(is.nan(x))}))
  type_keep_decider$Cell_type <- rownames(type_keep_decider)
  if(max(type_keep_decider$Number_NaNs) == 0){break}
  type_keep_decider <- merge(type_keep_decider, 
                             data.frame(Cell_type = tree_type_counts$Cell_type, 
                                        Type_frequency = rowSums(tree_type_counts[, -1])))
  type_keep_decider <- type_keep_decider[order(-type_keep_decider$Number_NaNs, type_keep_decider$Type_frequency), ]
  types_keep <- type_keep_decider$Cell_type[-1]
  cor_present <- cor_present[types_keep, types_keep]
}
baseline_type_correlations <- cor_present

# Cluster on correlations
norm_freq <- comparison_list$Normalized_frequencies
norm_freq <- norm_freq[, colnames(cor_present)]
cor_dist <- as.dist(1 - (cor_present + 1)/2)
cell_type_cor_clust <- hclust(cor_dist)
gaps_determinant <- data.frame(Node = rownames(norm_freq))
gaps_determinant$Tree <-
  sapply(gaps_determinant$Node, function(x) unlist(strsplit(as.character(x), ":"))[1])
gaps_row <- cumsum(table(gaps_determinant$Tree)[unique(gaps_determinant$Tree)])
# gaps_row <- cumsum(gaps_row)
# gaps_row <- gaps_row[unique(gaps_determinant$Tree)]
# pdf("./Images/Cell_type_relations_pheatmap_3dpi_min10_oneEndo_proper_treegaps.pdf")
pheatmap::pheatmap(t(norm_freq), cluster_cols = F, cluster_rows = cell_type_cor_clust,
                   cutree_rows = cluster_number, gaps_col = gaps_row,
                   border_color = NA, breaks = 0:100/100,
                   show_colnames = T)
# dev.off()
# pdf("./Images/Cell_type_cell_type_relations_pheatmap_7dpi_min10_oneEndo.pdf",
#     width = 10)
pheatmap(cor_present,
         cluster_cols = cell_type_cor_clust,
         cluster_rows = cell_type_cor_clust,
         cutree_cols = cluster_number, cutree_rows = cluster_number,
         treeheight_col = 0, show_colnames = F, fontsize = 16)
# dev.off()
baseline_cor_clustering <- cell_type_cor_clust

# Validate lineage linked cell types by downsampling ####
sampling_seed <- 1
cluster_number <- 7
sample_fraction <- 0.5
iterations <- 1000
cell_type_correspondence <- matrix(0, nrow = length(included_types), ncol = length(included_types),
                                   dimnames = list(included_types, included_types))

set.seed(sampling_seed)
for(iter in 1:iterations){
  print(iter)
  comparison_list <- ExtractTreeFrequencies(tree_list, sample_fraction, sampling_seed)
  
  # Calculate correlations
  type_tree_cor <- CalculateTypeCorrelations(tree_list, comparison_list, 
                                             force_include = types_force_include, 
                                             inclusion_limit = inclusion_limit)
  cor_present <- type_tree_cor[rownames(type_tree_cor) != "Dead cells", colnames(type_tree_cor) != "Dead cells"]
  # Remove cell types with largest amount of NaNs consecutively; in case of ties, remove the cell
  # type with the lowest number of cells.
  while(T){
    type_keep_decider <- data.frame(Number_NaNs = apply(cor_present, 1, function(x) {sum(is.nan(x))}))
    type_keep_decider$Cell_type <- rownames(type_keep_decider)
    if(max(type_keep_decider$Number_NaNs) == 0){break}
    type_keep_decider <- merge(type_keep_decider, 
                               data.frame(Cell_type = tree_type_counts$Cell_type, 
                                          Type_frequency = rowSums(tree_type_counts[, -1])))
    type_keep_decider <- type_keep_decider[order(-type_keep_decider$Number_NaNs, type_keep_decider$Type_frequency), ]
    types_keep <- type_keep_decider$Cell_type[-1]
    cor_present <- cor_present[types_keep, types_keep]
  }
  
  # Cluster on correlations
  norm_freq <- comparison_list$Normalized_frequencies
  norm_freq <- norm_freq[, colnames(cor_present)] #colnames(norm_freq) %in% colnames(cor_present)]
  cor_dist <- as.dist(1 - (cor_present + 1)/2)
  cell_type_cor_clust <- hclust(cor_dist)

  # Record which cell types are together in a cluster 
  this_ds_hclust <- data.frame(cutree(cell_type_cor_clust, k = cluster_number))
  colnames(this_ds_hclust)[1] <- "DS" 
  this_ds_hclust$Cell_type <- factor(row.names(this_ds_hclust), levels = included_types)
  this_ds_hclust$Presence <- 1
  type_clusters <- reshape2::acast(this_ds_hclust, Cell_type ~ DS, value.var = "Presence", fill = 0, drop = F)
  type_clusters_add <- type_clusters %*% t(type_clusters)
  type_clusters_add <- type_clusters_add[included_types, included_types]
  cell_type_correspondence <- cell_type_correspondence + type_clusters_add
}

# Summarize resampling results in plot
cell_type_correspondence <- cell_type_correspondence[rowSums(cell_type_correspondence) > 0, 
                                                     colSums(cell_type_correspondence) > 0]
# pdf(paste("./Images/Cell_type_relations_pheatmap_", timepoint, "dpi_", cluster_number, "cluster_", sample_fraction, "subsample", sep = ""))
# pdf("./Images/Cell_type_relations_pheatmap_3dpi_min10_oneEndo_1000resample_05rate.pdf")
pheatmap(cell_type_correspondence, cutree_rows = cluster_number, cutree_cols = cluster_number,
         treeheight_col = 0, show_colnames = F)
# dev.off()

min_occurrence <- cell_type_correspondence
min_occurrence[,] <- NA
# cell_type_correspondence_norm <- cell_type_correspondence/diag(cell_type_correspondence)
for(i in 1:ncol(cell_type_correspondence)){
  for(j in 1:ncol(cell_type_correspondence)){
    min_occurrence[i, j] <- min(diag(cell_type_correspondence)[i], diag(cell_type_correspondence)[j])
  }
}
cell_type_correspondence_norm <- cell_type_correspondence/min_occurrence

# pdf("./Images/Cell_type_relations_pheatmap_7dpi_min10_oneEndo_1000resample_05rate_orig_clustering_norm.pdf")
pheatmap(cell_type_correspondence_norm[baseline_cor_clustering$labels, baseline_cor_clustering$labels], 
         cluster_cols = baseline_cor_clustering, cluster_rows = baseline_cor_clustering,
         cutree_rows = cluster_number, cutree_cols = cluster_number,
         treeheight_col = 0, show_colnames = F)
# dev.off()

# Bootstrapping tree nodes to assess correlation confidence intervals ####
## Example plots for Hr1
# Prune nodes nd0_5, nd0_6, nd0_2_2, nd0_2_4, nd0_2_1_2.
example_node_sampling_tree <- Clone(tree_list_in[[5]]$Tree)
Prune(example_node_sampling_tree, 
      function(x){!(x$name %in% c("nd0_5", "nd0_6", "nd0_2_2", 
                                  "nd0_2_4", "nd0_2_1_2"))})
example_node_sampling_pie_tree <-
  collapsibleTree(df = example_node_sampling_tree, 
                  root = example_node_sampling_tree$scar, pieNode = T,
                  pieSummary = T, collapsed = F,
                  width = 1000, height = 500,
                  ctypes = names(celltype_colors), linkLength=60,
                  ct_colors = as.character(celltype_colors), angle = pi/2,fontSize = 0,
                  nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6))
example_node_sampling_pie_tree
# Recalculate node frequencies.
example_edge_list <- ToDataFrameNetwork(example_node_sampling_tree, "Cell.type")

example_node_type_counts <- 
  data.frame(table(example_edge_list$from, example_edge_list$Cell.type))
colnames(example_node_type_counts) <- c("Node", "Cell_type", "Type_count")
example_node_type_counts$Node <- as.character(example_node_type_counts$Node)
example_node_type_counts$Cell_type <- 
  as.character(example_node_type_counts$Cell_type)
example_node_type_counts <- 
  example_node_type_counts[example_node_type_counts$Cell_type != "NA", ]
example_node_type_counts <- 
  example_node_type_counts[example_node_type_counts$Node != "nd0_2_1_1", ]

example_tree_node_counts_df <- 
  reshape2::dcast(example_node_type_counts, Node ~ Cell_type, 
                  value.var = "Type_count")
rownames(example_tree_node_counts_df) <- example_tree_node_counts_df$Node
example_tree_node_counts_df <- 
  example_tree_node_counts_df[, -which(colnames(example_tree_node_counts_df) == "Node")]
#example_tree_node_counts_df <- 
#  example_tree_node_counts_df[, colSums(example_tree_node_counts_df) >= inclusion_limit, drop = F]

example_tree_node_cumul_counts <- example_tree_node_counts_df
example_tree_node_cumul_counts[, ] <- NA
for(i in 1:nrow(example_tree_node_cumul_counts)){
  node <- rownames(example_tree_node_cumul_counts)[i]
  example_tree_node_cumul_counts[i, ] <-
    colSums(example_tree_node_counts_df[grep(paste(node, "(_|$)", sep = ""), 
                                             rownames(example_tree_node_counts_df)), ])
}
example_tree_node_cumul_freqs <- 
  example_tree_node_cumul_counts/rowSums(example_tree_node_cumul_counts)

# Plot node frequencies.
pheat_example_tree_node_cumul_freqs <- 
  example_tree_node_cumul_freqs[, colnames(example_tree_node_cumul_freqs) %in% cell_type_cor_clust$labels]
pheat_example_tree_node_cumul_freqs$Add <- 0
colnames(pheat_example_tree_node_cumul_freqs)[ncol(pheat_example_tree_node_cumul_freqs)] <-
  setdiff(cell_type_cor_clust$labels, colnames(example_tree_node_cumul_freqs))
# pdf("./Images/Cell_type_relations_downsampling_example_pheatmap_Hr1.pdf")
pheatmap(t(pheat_example_tree_node_cumul_freqs[, baseline_cor_clustering$labels]), 
         cluster_cols = F, 
         cluster_rows = baseline_cor_clustering,
         cutree_rows = cluster_number, 
         border_color = NA, breaks = 0:100/100)
# dev.off()
# Adjust, save heatmap and put into supp. figs.

# tree_list_in[[t]]$Pie_tree <-
#   collapsibleTree(df = tree_list_in[[t]]$Tree, root =tree_list_in[[t]]$Tree$scar, pieNode = T,
#                   pieSummary = T, collapsed = F,
#                   width = 1000, height = 500,
#                   ctypes = names(celltype_colors), linkLength=60,
#                   ct_colors = as.character(celltype_colors), angle = pi/2,fontSize = 0,
#                   nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6))
# htmlwidgets::saveWidget(
#   tree_list_in[[t]]$Pie_tree,
#   file = paste("~/Documents/Projects/heart_Bo/Images/Trees/tree_", 
#                names(tree_list_in)[t], "_LINNAEUS_pie.html", sep = ""))

# Procedure: randomly remove a fraction of leaf nodes and their upstream nodes (selecting clones); then
# recalculate correlations.
# Set timepoint, fraction to remove, number of bootstrap runs
timepoint <- "7"
cluster_number <- 7
node_sampling_fraction <- 0.5
node_sampling_seed <- 42
bootstrap_runs <- 1000
cell_type_correspondence <- matrix(0, nrow = length(included_types), ncol = length(included_types),
                                   dimnames = list(included_types, included_types))

# Initialize tree list, summary dataframe and start bootstrap runs
trees_to_include <- 
  unlist(lapply(tree_list_in, 
                function(x){x$metadata$Name[x$metadata$dpi == timepoint]}))
trees_to_include <- trees_to_include[!is.na(trees_to_include)]
tree_list <- tree_list_in[trees_to_include]
clone_bootstraps <- data.frame(Type_1 = character(), Type_2 = character())
nf_ac_st_bootstraps <- data.frame(Type_1 = character(), Type_2 = character())

set.seed(node_sampling_seed)
for(strap in 1:bootstrap_runs){
  print(paste("Bootstrap", strap))
  
  # Randomly select leaf nodes for removal; calculate resulting cumulative cell numbers
  # and correlations. 
  # cumul_freqs_start_time <- Sys.time()
  cumul_freqs_novert_alltrees <- data.frame()
  all_tree_node_sizes <- data.frame()
  
  # Loop over trees to get cumulative frequencies and weights for selected nodes.
  for(t in 1:length(tree_list)){
    # Identify the nodes that we keep: make dataframe of nodes, identify leaf nodes and 
    # randomly include or remove them. Then remove nodes upstream of removed leaf nodes.
    nodes <- 
      data.frame(Node = paste(tree_list[[t]]$metadata$Name, 
                              unique(tree_list[[t]]$Node_type_counts$Node), sep = ":"), 
                 stringsAsFactors = F)
    nodes$Leaf <-
      sapply(nodes$Node,
             function(x) sum(grepl(paste(x, "(_|$)", sep = ""), nodes$Node)) == 1)
    nodes$Roll <- runif(nrow(nodes))
    # For each node: if we can find an underlying or equal leaf that remains included -> keep the node.
    # x <- nodes$Node[1]
    nodes$Include <- sapply(nodes$Node, 
                            function(x){
                              underlying_equal_nodes <- grepl(paste(x, "(_|$)", sep = ""), nodes$Node)
                              made_roll <- nodes$Roll <= node_sampling_fraction
                              return(sum(underlying_equal_nodes & nodes$Leaf & made_roll) > 0)
                            })
    # Skip tree for this iteration if by any chance no clones are left.
    if(sum(nodes$Include) == 0){
      print(paste("All clones removed from tree", tree_list[[t]]$metadata$Name, "; skipping tree this iteration."))
      next
    }
    
    # Determine node count matrix and remove cell types under inclusion limit
    tree_node_counts <- tree_list[[t]]$Node_type_counts
    tree_node_counts$Node <- paste(tree_list[[t]]$metadata$Name, tree_node_counts$Node, sep = ":")
    tree_node_counts <- merge(tree_node_counts, nodes[nodes$Include, ])
    tree_node_counts_df <- 
      reshape2::dcast(tree_node_counts, Node ~ Cell_type, value.var = "Type_count")
    rownames(tree_node_counts_df) <- tree_node_counts_df$Node
    tree_node_counts_df <- tree_node_counts_df[, -which(colnames(tree_node_counts_df) == "Node")]
    tree_node_counts_df <- tree_node_counts_df[, colSums(tree_node_counts_df) >= inclusion_limit, drop = F]
    if(ncol(tree_node_counts_df) < 2){
      print(paste("Less than two types over or equal inclusion limit for tree", tree_list[[t]]$metadata$Name,
                  ". Skipping tree this iteration."))
      next
    }
    
    # Calculate cumulative node frequencies
    tree_node_cumul_counts <- tree_node_counts_df
    tree_node_cumul_counts[, ] <- NA
    for(i in 1:nrow(tree_node_cumul_counts)){
      node <- rownames(tree_node_cumul_counts)[i]
      tree_node_cumul_counts[i, ] <-
        colSums(tree_node_counts_df[grep(paste(node, "(_|$)", sep = ""), rownames(tree_node_counts_df)), ])
    }
    tree_node_cumul_freqs <- tree_node_cumul_counts/rowSums(tree_node_cumul_counts)
    
    # Prune vertical branches
    # If a parent has only one kid, remove the kid. I.e.: keep only nodes with siblings.
    # So for each node, remove the last indicator (split on _, join back without last indicator),
    # then count nodes that are exactly one level down.
    nonvertical_nodes <-
      sapply(rownames(tree_node_cumul_freqs),
             function(x){
               # Exception case for head of tree.
               if(!grepl("_", x)){
                 return(TRUE)
               }else{
                 y <- unlist(strsplit(x, "_"))
                 sparent <- paste(y[-length(y)], collapse = "_")
                 siblings_and_self <- grepl(paste(sparent, "_[0-9]+$", sep = ""), rownames(tree_node_cumul_freqs))
                 return(sum(siblings_and_self) > 1)
               }
             })
    tree_node_cumul_freqs_novert <- tree_node_cumul_freqs[nonvertical_nodes, ]
    tree_node_cumul_counts_novert <- tree_node_cumul_counts[nonvertical_nodes, ]
    
    # Add cumulative frequencies and node sizes to dataframe for all trees
    if(ncol(cumul_freqs_novert_alltrees) == 0){
      cumul_freqs_novert_alltrees <- tree_node_cumul_freqs_novert
    }else{
      cumul_freqs_novert_alltrees[setdiff(names(tree_node_cumul_freqs_novert), names(cumul_freqs_novert_alltrees))] <- NA
      tree_node_cumul_freqs_novert[setdiff(names(cumul_freqs_novert_alltrees), names(tree_node_cumul_freqs_novert))] <- NA
      cumul_freqs_novert_alltrees <- rbind(cumul_freqs_novert_alltrees, tree_node_cumul_freqs_novert)
    }
    
    # Record node sizes for weighted correlation computation
    tree_node_sizes <- data.frame(Node = rownames(tree_node_cumul_counts_novert),
                                  Size = rowSums(tree_node_cumul_counts_novert))
    all_tree_node_sizes <- rbind(all_tree_node_sizes, tree_node_sizes)
  }
  
  # Set NA's to zero (no cells of a type in a tree should translate to frequencies of 0)
  cumul_freqs_novert_alltrees[is.na(cumul_freqs_novert_alltrees)] <- 0
  # Remove cell types that aren't present at all
  cumul_freqs_novert_alltrees <- cumul_freqs_novert_alltrees[, colSums(cumul_freqs_novert_alltrees) > 0]
  
  # Ensure weights and cumulative frequencies have the same order
  all_tree_node_sizes <- all_tree_node_sizes[rownames(cumul_freqs_novert_alltrees), ]
  # cumul_freqs_end_time <- Sys.time()
  # print(paste("Cumulative frequencies time:", cumul_freqs_end_time - cumul_freqs_start_time))

  # Calculate weighted correlations (this will automatically ignore NAs, but will return NaN if there are
  # no or only one node to correlate on). Remove Dead cells.
  # weighted_cors_start_time <- Sys.time()
  cor_present <- wtd.cors(cumul_freqs_novert_alltrees, weight = all_tree_node_sizes$Size)
  cor_present <- cor_present[rownames(cor_present) != "Dead cells", colnames(cor_present) != "Dead cells"]
  
  # Remove cell types with largest amount of NaNs consecutively; in case of ties, remove the cell
  # type with the lowest number of cells.
  while(T){
    type_keep_decider <- data.frame(Number_NaNs = apply(cor_present, 1, function(x) {sum(is.nan(x))}))
    type_keep_decider$Cell_type <- rownames(type_keep_decider)
    if(max(type_keep_decider$Number_NaNs) == 0){break}
    type_keep_decider <- merge(type_keep_decider, 
                               data.frame(Cell_type = tree_type_counts$Cell_type, 
                                          Type_frequency = rowSums(tree_type_counts[, -1])))
    type_keep_decider <- type_keep_decider[order(-type_keep_decider$Number_NaNs, type_keep_decider$Type_frequency), ]
    types_keep <- type_keep_decider$Cell_type[-1]
    cor_present <- cor_present[types_keep, types_keep]
  }
  # weighted_cors_end_time <- Sys.time()
  # print(paste("Weighted correlations time:", weighted_cors_end_time - weighted_cors_start_time))

  # Log correlations into summary dataframe.
  # logging_start_time <- Sys.time()
  this_clone_bootstrap <-
    reshape2::melt(cor_present)
  colnames(this_clone_bootstrap) <- c("Type_1", "Type_2", paste("Correlation", strap, sep = "_"))
  clone_bootstraps <- merge(clone_bootstraps, this_clone_bootstrap, all = T)
  # logging_end_time <- Sys.time()
  # print(paste("Logging time:", logging_end_time - logging_start_time))
  
  # Cluster on correlations
  cor_dist <- as.dist(1 - (cor_present + 1)/2)
  cell_type_cor_clust <- hclust(cor_dist)
  
  # Record which cell types are together in a cluster 
  this_ds_hclust <- data.frame(cutree(cell_type_cor_clust, k = cluster_number))
  colnames(this_ds_hclust)[1] <- "DS" 
  this_ds_hclust$Cell_type <- factor(row.names(this_ds_hclust), levels = included_types)
  this_ds_hclust$Presence <- 1
  type_clusters <- reshape2::acast(this_ds_hclust, Cell_type ~ DS, value.var = "Presence", fill = 0, drop = F)
  type_clusters_add <- type_clusters %*% t(type_clusters)
  type_clusters_add <- type_clusters_add[included_types, included_types]
  cell_type_correspondence <- cell_type_correspondence + type_clusters_add
  
  # HERE
  # Bootstrapping conditional probabilities
  nf_leaf <- cumul_freqs_novert_alltrees
  nf_leaf$Leaf <-
    sapply(row.names(nf_leaf),
           function(x) sum(grepl(x, row.names(nf_leaf))) == 1)
  nf_leaf <- nf_leaf[nf_leaf$Leaf, c(ncol(nf_leaf), 1:(ncol(nf_leaf) - 1))]
  
  nf_ac_st <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                     dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
  nf_ac_st_counts <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                            dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
  for(i in 1:(ncol(nf_leaf) - 1)){
    for(j in 1:(ncol(nf_leaf) - 1)){
      nf_ac_st[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)/sum(nf_leaf[, j + 1] > 0)
    }
  }
  
  # Log conditional probabilities.
  nf_ac_st_long <- reshape2::melt(nf_ac_st)
  colnames(nf_ac_st_long) <- c("Type_1", "Type_2", paste("Conditional_prob", strap, sep = "_"))
  nf_ac_st_bootstraps <- merge(nf_ac_st_bootstraps, nf_ac_st_long, all = T)
}

# 4) Analyze bootstrap results.
# While the correlations can sometimes change quite wildly, it is very hard to say if this makes
# any difference in the clustering outcome.
# example_type_1 <- "Fibroblasts (nppc)" #"Macrophages" #"Epicardium (Ventricle)" #"T-cells"
# example_type_2 <- "Endocardium (Ventricle)" #"Valve fibroblasts" #Monocytes" #"Fibroblasts (const.)" 
# bootstrap_example_histogram <-
#   data.frame(t(clone_bootstraps[clone_bootstraps$Type_1 == example_type_1 &
#                                   clone_bootstraps$Type_2 == example_type_2, -c(1:2)]))
# colnames(bootstrap_example_histogram) <- "B.cors"
# ggplot(bootstrap_example_histogram) +
#   geom_histogram(aes(x = B.cors)) +
#   geom_vline(xintercept = baseline_type_correlations[example_type_1, example_type_2], color = "red") +
#   labs(title = paste(example_type_1, "-", example_type_2),
#        x = "Clone-bootstrapped correlations")
# sd(bootstrap_example_histogram$B.cors, na.rm = T)

# Summarize resampling results in plot
cell_type_correspondence <- cell_type_correspondence[rowSums(cell_type_correspondence) > 0, 
                                                     colSums(cell_type_correspondence) > 0]
min_occurrence <- cell_type_correspondence
min_occurrence[,] <- NA
# cell_type_correspondence_norm <- cell_type_correspondence/diag(cell_type_correspondence)
for(i in 1:ncol(cell_type_correspondence)){
  for(j in 1:ncol(cell_type_correspondence)){
    min_occurrence[i, j] <- min(diag(cell_type_correspondence)[i], diag(cell_type_correspondence)[j])
  }
}
cell_type_correspondence_norm <- cell_type_correspondence/min_occurrence

# pdf("./Images/Cell_type_relations_pheatmap_7dpi_min10_oneEndo_1000_clone_downsample_baseline_ordering_norm.pdf", width = 10)
pheatmap(cell_type_correspondence_norm[baseline_cor_clustering$labels, baseline_cor_clustering$labels], 
         cluster_cols = baseline_cor_clustering, cluster_rows = baseline_cor_clustering, 
         cutree_rows = cluster_number, cutree_cols = cluster_number,
         treeheight_col = 0, show_colnames = F, fontsize = 16)
# dev.off()

cell_type_correspondence_long <-
  reshape2::melt(cell_type_correspondence)
colnames(cell_type_correspondence_long) <- c("Type_1", "Type_2", "Coclustering_count")
original_clustering <- data.frame(cutree(baseline_cor_clustering, k = cluster_number))
colnames(original_clustering)[1] <- "Original_hclust" 
original_clustering$Cell_type <- factor(row.names(original_clustering), levels = included_types)
original_clustering$Presence <- 1
orig_clust_matrix <- reshape2::acast(original_clustering, Cell_type ~ Original_hclust, value.var = "Presence", fill = 0)
orig_coclustering <- reshape2::melt(orig_clust_matrix %*% t(orig_clust_matrix))
colnames(orig_coclustering) <- c("Type_1", "Type_2", "Original_hclust")
clone_bootstrap_cluster_assessment <- merge(cell_type_correspondence_long, orig_coclustering)

# Targeted search for potential source cell types ####
# Focus on tree leaf nodes
nf_leaf <- baseline_comparison_list$Normalized_frequencies
nf_leaf$Leaf <-
  sapply(row.names(nf_leaf),
         function(x) sum(grepl(x, row.names(nf_leaf))) == 1)
nf_leaf <- nf_leaf[nf_leaf$Leaf, c(ncol(nf_leaf), 1:(ncol(nf_leaf) - 1))]

ex_type_1 <- "Epicardium (Ventricle)"
ex_type_2 <- "Fibroblasts (col12a1a)"
sum(nf_leaf[, ex_type_1] > 0)
sum(nf_leaf[, ex_type_2] > 0)
sum((nf_leaf[, ex_type_2] > 0) && nf_leaf[, ex_type_1] > 0)

# Calculate conditional probabilities: the chance of seeing (row cell type) when observing
# (column cell type)
nf_ac_st <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                   dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
nf_ac_st_counts <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                   dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
for(i in 1:(ncol(nf_leaf) - 1)){
  for(j in 1:(ncol(nf_leaf) - 1)){
    nf_ac_st[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)/sum(nf_leaf[, j + 1] > 0)
    nf_ac_st_counts[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)
  }
}
# View(nf_ac_st_counts[c("Endocardium (Ventricle)", "Fibroblasts (nppc)"), c("Endocardium (Ventricle)", "Fibroblasts (nppc)")])
# View(nf_ac_st[c("Endocardium (Ventricle)", "Fibroblasts (nppc)"), c("Endocardium (Ventricle)", "Fibroblasts (nppc)")])
# 
# View(nf_ac_st_counts[c("Epicardium (Ventricle)", "Fibroblasts (col12a1a)"), c("Epicardium (Ventricle)", "Fibroblasts (col12a1a)")])
# View(nf_ac_st[c("Epicardium (Ventricle)", "Fibroblasts (col12a1a)"), c("Epicardium (Ventricle)", "Fibroblasts (col12a1a)")])


# Find sources for indicated targets and plot
start_targets <-  "Fibroblasts (spock3)"
  # "Fibroblasts (nppc)" #c("Fibroblasts (nppc)", "Fibroblasts (spock3)", "Valve fibroblasts") #"Fibroblasts (cfd)" #"Valve fibroblasts" # c("Fibroblasts (nppc)", "Fibroblasts (spock3)") # # #c("Fibroblasts (col12a1a)")
current_targets <- start_targets
nf_ac_current <- nf_ac_st[, current_targets, drop = F]
current_targets <- names(which(rowSums(nf_ac_current >= 0.8) > 0))
ac_tograph <- nf_ac_st[current_targets, current_targets, drop = F]
# ac_tograph_2 <- ac_tograph

nodenames <- 1:nrow(ac_tograph)
longnames <- rownames(ac_tograph)
graph1<-qgraph(ac_tograph, diag = F,
               layout= "spring", edge.color = "black",
               labels = F, 
               minimum = 0, cut = 0.8, legend = F,
               color = as.character(celltype_colors[rownames(ac_tograph)]),
               borders = F,
               # filename = "~/Documents/Projects/heart_Bo/Images/Asymmetric_coinc_network_7dpi_to_spock3_all_conn",
               # filetype = "eps",
               vsize=15, esize = 20, asize = 9, width = 1.5, height = 1.5)
# Plot for arrow thickness legends
qgraph_arrow_legend <- ac_tograph
qgraph_arrow_legend[,] <- 0
diag(qgraph_arrow_legend) <- 1
qgraph_arrow_legend[2, 1] <- 0.8
qgraph_arrow_legend[3, 1] <- 0.9
qgraph_arrow_legend[4, 1] <- 1
graph1<-qgraph(qgraph_arrow_legend, diag = F,
               layout= "spring", edge.color = "black",
               labels = F,
               minimum = 0.79999, cut = 0.8, legend = F,
               borders = F,
               # filename = "~/Documents/Projects/heart_Bo/Images/coinc_network_legend",
               # filetype = "eps",
               vsize=15, esize = 20, asize = 9, width = 1.5, height = 1.5)

# Compare with bootstrapping
comparison_variable <- "Fibroblasts (nppc)"
comparison_condition <- "Endocardium (Ventricle)"

comparison_bootstrap_values <- 
  data.frame(Conditional_prob_boot =
               as.vector(t(nf_ac_st_bootstraps[nf_ac_st_bootstraps$Type_1 == comparison_variable &
                                                     nf_ac_st_bootstraps$Type_2 == comparison_condition, -c(1:2)])))
# png(paste("./Images/Conditional_probability_bootstrap_7dpi_",
#           comparison_variable, "_when_", comparison_condition, ".png"),
#     width = 750, height = 400)
ggplot(comparison_bootstrap_values) +
  geom_histogram(aes(x = Conditional_prob_boot)) +
  geom_vline(xintercept = nf_ac_st[comparison_variable, comparison_condition], color = "red", size = 3) +
  scale_x_continuous(limits = c(-0.05, 1.05)) +
  labs(x = "Conditional probabilities",
       title = paste(comparison_variable, "|", comparison_condition))
# dev.off()
