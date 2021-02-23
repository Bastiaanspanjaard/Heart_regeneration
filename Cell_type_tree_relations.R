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

# Finding cell types linked by lineage ####
# Determine trees and cell types to include
timepoint <- "7"
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

# Scatterplot of frequencies highlighting different trees
cor_type_1 <- "Fibroblasts (nppc)"
cor_type_2 <- "Endocardium (Ventricle)"
df_ct <- merge(comparison_list$Normalized_frequencies, comparison_list$Node_sizes, by = 0)
df_ct <- df_ct[, c("Row.names", cor_type_1, cor_type_2, "Size")]
colnames(df_ct) <- c("Node", "CT1", "CT2", "Node_size")
df_ct$Tree <- sapply(df_ct$Node,
                     function(x) unlist(strsplit(x, ":"))[1])
df_ct$Tree <- factor(df_ct$Tree, 
                     levels = c("Hr1", "Hr2", "Hr6", "Hr7",
                                "Hr13", "Hr14", "Hr15"))
# pdf("./Images/Fibnppc_endoV_7dpi_corplot_trees.pdf",
#     width = 5, height = 4, useDingbats = F)
ggplot(df_ct) +
  geom_point(aes(x = CT1, y = CT2, color = Tree), size = 4) +#, size = Node_size)) +
  labs(x = cor_type_1, y = cor_type_2) +
  scale_color_brewer(palette = 'Set1') +
  theme(panel.grid.minor = element_blank())#,
# title = paste("7dpi node frequencies all trees, cor ", round(cor(df_ct$CT1, df_ct$CT2), 2), sep = ""))
# dev.off()

# Scatterplot of frequencies including only nodes that are over inclusion limit.
norm_freq <- comparison_list$Normalized_frequencies
inc_exc <- comparison_list$Normalized_frequencies # Same structure as norm_freq
inc_exc[,] <- F # Set all entries to F
for(t in 1:length(tree_list)){
  tree_name <- names(tree_list)[t]
  tree_type_include <- aggregate(tree_list[[t]]$Node_type_counts$Type_count,
                                 by = list(Cell_type = tree_list[[t]]$Node_type_counts$Cell_type),
                                 sum) # Count cell type numbers in tree
  tree_type_include$Include <- tree_type_include$x >= inclusion_limit # Set included celltypes to T in dataframe
  inc_exc[grepl(tree_name, rownames(inc_exc)), tree_type_include$Cell_type] <- # Nodes in tree, cell types in tree (cells not in tree remain F)
    matrix(rep(tree_type_include$Include, each = sum(grepl(tree_name, rownames(inc_exc)))),
           nrow = sum(grepl(tree_name, rownames(inc_exc))),
           ncol = nrow(tree_type_include)) # Create matrix of (nodes in tree) x (cell types in tree), with values T for cell types that
  # have enough abundance.
}
x_plot_name <- "Epicardium (Ventricle)" #
y_plot_name <- "Fibroblasts (col12a1a)" #"Fibroblast-like cells"
x_plot_col <- which(colnames(norm_freq) == x_plot_name)
y_plot_col <- which(colnames(norm_freq) == y_plot_name)
nf_plot <- cbind(norm_freq[, c(x_plot_col, y_plot_col)], inc_exc[, c(x_plot_col, y_plot_col)])
colnames(nf_plot) <- c("Xp", "Yp", "Inc_X", "Inc_Y")
nf_plot <- nf_plot[nf_plot$Inc_X & nf_plot$Inc_Y, ]
# pdf("./Images/Fibcol12_epiV_3dpi_corplot.pdf",
#     width = 4, height = 4, useDingbats = F)
ggplot(nf_plot) +
  geom_point(aes(x = Xp, y = Yp), size = 4) +
  labs(x = x_plot_name, y = y_plot_name) +
  theme(text = element_text(size = 16),
        panel.grid.minor = element_blank())
# dev.off()

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

# Cluster on correlations
norm_freq <- comparison_list$Normalized_frequencies
norm_freq <- norm_freq[, colnames(cor_present)]
cor_dist <- as.dist(1 - (cor_present + 1)/2)
cell_type_cor_clust <- hclust(cor_dist)
gaps_determinant <- data.frame(Node = rownames(norm_freq))
gaps_determinant$Tree <-
  sapply(gaps_determinant$Node, function(x) unlist(strsplit(as.character(x), ":"))[1])
gaps_row <- cumsum(table(gaps_determinant$Tree))
# pdf("./Images/Cell_type_relations_pheatmap_3dpi_min10_oneEndo.pdf")
pheatmap::pheatmap(t(norm_freq), cluster_cols = F, cluster_rows = cell_type_cor_clust,
                   cutree_rows = cluster_number, gaps_col = gaps_row,
                   border_color = NA, breaks = 0:100/100,
                   show_colnames = F)
# dev.off()
# pdf("./Images/Cell_type_cell_type_relations_pheatmap_7dpi_min10_oneEndo.pdf",
#     width = 10)
pheatmap(cor_present,
         cluster_cols = cell_type_cor_clust,
         cluster_rows = cell_type_cor_clust,
         cutree_cols = cluster_number, cutree_rows = cluster_number,
         treeheight_col = 0, show_colnames = F, fontsize = 16)
# dev.off()

# Validate lineage linked cell types by downsampling ####
sampling_seed <- 1
cluster_number <- 7
sample_fraction <- 0.5
iterations <- 1000
cell_type_correspondence <- matrix(0, nrow = length(included_types), ncol = length(included_types),
                                   dimnames = list(included_types, included_types))

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
cell_type_correspondence <- cell_type_correspondence[rowSums(cell_type_correspondence) > 0, colSums(cell_type_correspondence) > 0]
# pdf(paste("./Images/Cell_type_relations_pheatmap_", timepoint, "dpi_", cluster_number, "cluster_", sample_fraction, "subsample", sep = ""))
# pdf("./Images/Cell_type_relations_pheatmap_3dpi_min10_oneEndo_1000resample_05rate.pdf")
pheatmap(cell_type_correspondence, cutree_rows = cluster_number, cutree_cols = cluster_number,
         treeheight_col = 0, show_colnames = F)
# dev.off()

# Targeted search for potential source cell types ####
# Focus on tree leaf nodes
nf_leaf <- comparison_list$Normalized_frequencies
nf_leaf$Leaf <-
  sapply(row.names(nf_leaf),
         function(x) sum(grepl(x, row.names(nf_leaf))) == 1)
nf_leaf <- nf_leaf[nf_leaf$Leaf, c(ncol(nf_leaf), 1:(ncol(nf_leaf) - 1))]

# Calculate asymmetric coincidences for source(row) - target(column)
nf_ac_st <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                   dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
for(i in 1:(ncol(nf_leaf) - 1)){
  for(j in 1:(ncol(nf_leaf) - 1)){
    nf_ac_st[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)/sum(nf_leaf[, j + 1] > 0)
  }
}

# Find sources for indicated targets and plot
start_targets <- c("Fibroblasts (spock3)","Fibroblasts (nppc)", "Valve fibroblasts")
current_targets <- start_targets
nf_ac_current <- nf_ac_st[, current_targets]#, drop = F]
current_targets <- names(which(rowSums(nf_ac_current >= 0.8) > 0))
ac_tograph <- nf_ac_st[current_targets, current_targets]
nodenames <- 1:nrow(ac_tograph)
longnames <- rownames(ac_tograph)
graph1<-qgraph(ac_tograph, diag = F,
               layout= "spring", edge.color = "black",
               labels = F, 
               minimum = 0.79999, cut = 0.8, legend = F,
               color = as.character(celltype_colors[rownames(ac_tograph)]),
               borders = F,
               # filename = "~/Documents/Projects/heart_Bo/Images/Asymmetric_coinc_network_7dpi_to_Fibnppcspock3fiblike",
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

