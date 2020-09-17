# Dependencies ####
source("./Scripts/HR_library.R")

# Parameters ####
tree_path <- "./Data/Trees/Sbcf/"
target <- "Fibroblast (col11a1a)"
precursors_include_incl_target <- NULL
precursors_include <- setdiff(precursors_include_incl_target, target)
inclusion_limit <- 10 # How many cells are needed of a certain type to include a tree for
# that type.
correlation_type <- "Symmetric" # "Asymmetric": only correlate over the nodes with nonzero
# target cell type. "Symmetric": correlate over all nodes.
uniform_conversion <- F
uniform_conversion_chance <- 0.5
chances <- c(0, 0.1, 0.2, 0.3)
simulation_seed <- 37 # -> 0.16 # 69 -> 0.02 # 1337 -> 0.71 top # 6 -> 0.62 top # 420 -> 0.65, z-top 
# Old list without weighted mean conversion chance: 6 -> 0.52 # 1337 -> 0.52 non-significant #420 -> 0.57 # 69 -> rho = 0.51 # 37 -> rho = 0.3
zoom_types <- c("Fibroblast", "Fibroblast (cfd)", "Fibroblast (col11a1a)", 
                "Fibroblast (col12a1a)", "Fibroblast (mpeg1.1)", "Epicardium (Atrium)", 
                "Epicardium (Ventricle)", "Fibroblast (proliferating)", "Perivascular cells",
                "Fibroblast (cxcl12a)", "Fibroblast (nppc)", "Fibroblast (spock3)",
                "Fibroblast-like cells")
# 'Fibroblast', 'Fibroblast (cfd)', 'Fibroblast (col11a1a)', 'Fibroblast (col12a1a)', 
# 'Fibroblast (cxcl12a)', 'Fibroblast (mpeg1.1)', 'Fibroblast (nppc)', 'Fibroblast (spock3)',
# 'Fibroblast-like cells', 'Epicardium (Atrium)', 'Epicardium (Ventricle)',
# 'Fibroblast (proliferating)', 'Perivascular cells'

# Prepare colors and cell type data ####
# Cell type decisions as follows: all T-cell types and all macrophage types become T-cells/macrophages. 
# Merge endocardium frzb and do not split endocardium into 1 and 2.
# Merge cardiomyocytes A and V but keep the regular, ttn.2 and proliferating subtypes.
cell_annotations <- read.csv("./Data/final_metadata.csv", stringsAsFactors = F)
colnames(cell_annotations)[1] <- "Cell"
cell_annotations$Cell_type <- cell_annotations$lineage.ident
cell_type_renaming <-
  data.frame(Original_cell_type = sort(as.character(unique(cell_annotations$Cell_type))), stringsAsFactors = F)
cell_type_renaming$New_type_name <-
  c("B-cells", "Bl.ves.EC (apnln)", "Bl.ves.EC (lyve1)", 
    "Bl.ves.EC (plvapb)", "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2)",
    "Cardiomyocytes (ttn.2)", "Cardiomyocytes (Atrium)", "Cardiomyocytes (Ventricle)",
    "Dead cells", "Endocardium (Atrium)", "Endocardium (Ventricle)", "Endocardium (frzb)",
    "Endocardium (frzb)", "Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblasts (const.)",
    "Fibroblasts (cfd)", "Fibroblasts (col11a1a)", "Fibroblasts (col12a1a)", "Fibroblasts (cxcl12a)",
    "Fibroblasts (mpeg1.1)", "Fibroblasts (nppc)", "Fibroblasts (proliferating)", "Fibroblasts (spock3)",
    "Valve fibroblasts", "Macrophages", "Macrophages", "Macrophages",
    "Macrophages", "Macrophages", "Macrophages", 
    "Macrophages", "Macrophages", "Macrophages",
    "Macrophages", "Monocytes", "Myelin cells", "Neuronal cells", "Neutrophils",
    "Perivascular cells", "Proliferating cells", "Smooth muscle cells", "T-cells",
    "T-cells", "T-cells")
cell_annotations$Cell_type <-
  cell_type_renaming$New_type_name[match(cell_annotations$Cell_type, cell_type_renaming$Original_cell_type)]
# cell_annotations$Cell_type[cell_annotations$Cell_type == "Fibroblast-like cells"] <-
#   "Valve fibroblasts"
# cell_annotations$Cell_type[grepl("T-cell", cell_annotations$Cell_type)] <- "T-cells"
# cell_annotations$Cell_type[grepl("Macrophage", cell_annotations$Cell_type)] <- "Macrophages"
# cell_annotations$Cell_type[grepl("Endocardium frzb", cell_annotations$Cell_type)] <- "Endocardium (frzb)"
# cell_annotations$Cell_type[grepl("ttn", cell_annotations$Cell_type)] <- "Cardiomyocytes (ttn.2)"
cell_annotations <- cell_annotations[cell_annotations$Cell_type != "Dead cells", ]
cell_annotations$Cell_name <- paste("nd", cell_annotations$Cell, sep = "")
# write.csv(cell_annotations, "./Data/final_metadata_Tmacromerged_2.csv", row.names = F, quote = F)

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
# Start by calculating the correlations between a progenitor cell type and potential
# precursors. After this, investigate whether these correlations are stable when leaving
# out one of the trees. Finally, bootstrap the correlations by randomly placing the
# precursor cells over the trees to see what the null-distribution of the correlations
# is and how enriched the observed correlations are.

# Input: tree list, target
# The first loop here is over trees to create data frames of counts and frequencies per tree.
# Can we do most of this in the initialization? We'd need most of this for every later step and
# it also integrates well with creating edge lists and count tables.
# The second loop does the actual correlations, using the normalized frequency and count tables.
# Output: a dataframe with correlations between a progenitor and all other progenitors. Update
# this to a list (dataframe for every progenitor) and also keep a condensed long list with only
# the correlations, but between every two cell types.

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

# Calculate cell-type lineage correlations and cluster cell types
comparison_list <- ExtractTreeFrequencies(tree_list)

# Plot example correlations
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
norm_freq <- norm_freq[, colnames(cor_present)] #colnames(norm_freq) %in% colnames(cor_present)]
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

# Validate clustering: downsample trees, calculate lineage correlations and cluster cell types
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

# Run this over all trees, keep distances top-level, create a new tree-level in the list
set.seed <- 1
cluster_number <- 7
sample_fraction <- 1
cell_type_correspondence <- matrix(0, nrow = length(included_types), ncol = length(included_types),
                                   dimnames = list(included_types, included_types))

for(iter in 1:1000){
  print(iter)
  
  # Calculate cell type frequencies per node (comparison_list) 
  # and record cell type frequencies per tree (analysis_stats)
  comparison_list <- list(Comparison = data.frame(Tree = character(),
                                                  Precursor = character(),
                                                  Type_count = integer()),
                          Node_sizes = data.frame(Size = integer()),
                          Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = nrow(celltype_frequencies))),
                          Frequencies = data.frame(matrix(nrow = 0, ncol = nrow(celltype_frequencies))))
  colnames(comparison_list$Frequencies) <- celltype_frequencies$Cell_type
  colnames(comparison_list$Normalized_frequencies) <- celltype_frequencies$Cell_type
  analysis_stats <-
    data.frame(Tree = character(),
               Included = logical(),
               Cell_type = character(),
               Count = integer())

  # Extract cell type frequencies
  for(t in 1:length(tree_list)){
    edge_list_full <- tree_list[[t]]$Edge_list
    # Cell.type == "NA" are nodes; leave those in, sample over the others.
    cell_list <- edge_list_full[edge_list_full$Cell.type != "NA", ]
    cell_list_s <- cell_list[sample.int(nrow(cell_list), round(size = sample_fraction * nrow(cell_list))), ]
    edge_list <- rbind(edge_list_full[edge_list_full$Cell.type == "NA", ], cell_list_s)
    analysis_stats_add <-
      data.frame(table(edge_list$Cell.type[edge_list$Cell.type != "NA"]))
    colnames(analysis_stats_add) <- c("Cell_type", "Count")
    analysis_stats_add$Tree <- names(tree_list)[t]
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
    
    comparison_list$Node_sizes <- rbind(comparison_list$Node_sizes,
                                        data.frame(Size = rowSums(sample_type_cumulative)))
    comparison_tree <-
      data.frame(Tree = rep(names(tree_list)[t], ncol(sample_type_count)),
                 Precursor = names(sample_type_count),
                 Type_count = colSums(sample_type_count))
    
    comparison_list$Comparison <- rbind(comparison_list$Comparison, comparison_tree)
    
    # Add cumulative counts and normalized frequencies
    sample_type_cumulative[setdiff(names(comparison_list$Frequencies), 
                                   names(sample_type_cumulative))] <- NA
    comparison_list$Frequencies[setdiff(names(sample_type_cumulative), 
                                        names(comparison_list$Frequencies))] <- NA
    comparison_list$Frequencies <- 
      rbind(comparison_list$Frequencies, sample_type_cumulative)
    sample_type_nf[setdiff(names(comparison_list$Normalized_frequencies), names(sample_type_nf))] <- NA
    comparison_list$Normalized_frequencies[setdiff(names(sample_type_nf), 
                                                   names(comparison_list$Normalized_frequencies))] <- NA
    comparison_list$Normalized_frequencies <- 
      rbind(comparison_list$Normalized_frequencies, sample_type_nf)
  }
  comparison_list$Frequencies[is.na(comparison_list$Frequencies)] <- 0
  comparison_list$Frequencies <- 
    comparison_list$Frequencies[, colSums(comparison_list$Frequencies) > 0]
  comparison_list$Normalized_frequencies[is.na(comparison_list$Normalized_frequencies)] <- 0
  comparison_list$Normalized_frequencies <- 
    comparison_list$Normalized_frequencies[, colSums(comparison_list$Normalized_frequencies) > 0]
  comparison_list$Node_sizes$Weight <- comparison_list$Node_sizes$Size/sum(comparison_list$Node_sizes$Size)

  # Collapse non-splitting branches - SEPARATE FUNCTION:
  # 1) To find a non-splitting branch, identify nodes with exactly one child. Those are the parent nodes in
  # a non-splitting branch. We can recognize those nodes by counting the nodes whose name contains their name;
  # the nodes we are looking for have exactly two (themselves and their only child).
  # 2) Remove the child (note that weights and correlations are using the cumulative counts).
  # 3) Search for non-splitting branches again.
  repeat{
    nodes_list <- data.frame(Node = rownames(comparison_list$Normalized_frequencies))
    nodes_list$One_child <- 
      sapply(nodes_list$Node, function(x) sum(grepl(x, nodes_list$Node))) == 2
    if(sum(nodes_list$One_child) == 0){break}
    nodes_list$In_one_child_parent <- 
      sapply(nodes_list$Node, function(x){
        sum(sapply(nodes_list$Node[nodes_list$One_child], function(y){grepl(y, x)}))
      })
    nodes_list$Last_only_child <-
      !nodes_list$One_child & nodes_list$In_one_child_parent
    comparison_list$Frequencies <- 
      comparison_list$Frequencies[!(rownames(comparison_list$Frequencies) %in% nodes_list$Node[nodes_list$Last_only_child]), ]
    comparison_list$Normalized_frequencies <- 
      comparison_list$Normalized_frequencies[!(rownames(comparison_list$Normalized_frequencies) %in% nodes_list$Node[nodes_list$Last_only_child]), ]
    comparison_list$Node_sizes <- 
      comparison_list$Node_sizes[!(rownames(comparison_list$Node_sizes) %in% nodes_list$Node[nodes_list$Last_only_child]), ]
  }

  # Optional plot of correlations - BREAK OUT OF FUNCTION, uses comparison_list$Normalized frequencies +
  # comparison_list$Node_sizes
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

  # Calculate correlations
  # comparison_list$All_trees_distances <- 
  #   data.frame(Precursor = names(comparison_list$Normalized_frequencies),
  #              Weighted_cor_progpos = numeric(ncol(comparison_list$Normalized_frequencies)),
  #              Weighted_cor_se = numeric(ncol(comparison_list$Normalized_frequencies)))#,
  # comparison_list$All_trees_distances[, unique(analysis_stats$Tree[analysis_stats$Included])] <- 0
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
  norm_freq <- norm_freq[, colnames(cor_present)] #colnames(norm_freq) %in% colnames(cor_present)]
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

  # For resampling: record which cell types are together in a cluster and 
  this_ds_hclust <- data.frame(cutree(cell_type_cor_clust, k = cluster_number))
  colnames(this_ds_hclust)[1] <- "DS" #paste("DS", sample_fraction, iter, sep = "_")
  this_ds_hclust$Cell_type <- factor(row.names(this_ds_hclust), levels = included_types)
  this_ds_hclust$Presence <- 1
  type_clusters <- reshape2::acast(this_ds_hclust, Cell_type ~ DS, value.var = "Presence", fill = 0, drop = F)
  type_clusters_add <- type_clusters %*% t(type_clusters)
  type_clusters_add <- type_clusters_add[included_types, included_types]
  cell_type_correspondence <- cell_type_correspondence + type_clusters_add
  
  # START Making an example correlation plot. BREAK OUT OF FUNCTION
  # inc_exc <- norm_freq # Same structure as norm_freq
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
  # y_plot_name <- "Fibroblast (col12a1a)" #"Fibroblast-like cells" 
  # x_plot_col <- which(colnames(norm_freq) == x_plot_name)
  # y_plot_col <- which(colnames(norm_freq) == y_plot_name)
  # nf_plot <- cbind(norm_freq[, c(x_plot_col, y_plot_col)], inc_exc[, c(x_plot_col, y_plot_col)])
  # colnames(nf_plot) <- c("Xp", "Yp", "Inc_X", "Inc_Y")
  # nf_plot <- nf_plot[nf_plot$Inc_X & nf_plot$Inc_Y, ]
  # 
  # # pdf("./Images/Fibcol12_epiV_3dpi_corplot.pdf",
  # #     width = 4, height = 4, useDingbats = F)
  # ggplot(nf_plot) +
  #   geom_point(aes(x = Xp, y = Yp), size = 4) +
  #   labs(x = x_plot_name, y = y_plot_name) +
  #   theme(text = element_text(size = 16),
  #         panel.grid.minor = element_blank())
  # # dev.off()
  # END making example correlation plot
}

# Summarize resampling results in plot
cell_type_correspondence <- cell_type_correspondence[rowSums(cell_type_correspondence) > 0, colSums(cell_type_correspondence) > 0]
# pdf(paste("./Images/Cell_type_relations_pheatmap_", timepoint, "dpi_", cluster_number, "cluster_", sample_fraction, "subsample", sep = ""))
# pdf("./Images/Cell_type_relations_pheatmap_3dpi_min10_oneEndo_1000resample_05rate.pdf")
pheatmap(cell_type_correspondence, cutree_rows = cluster_number, cutree_cols = cluster_number,
         treeheight_col = 0, show_colnames = F)
# dev.off()

 
nf_leaf <- comparison_list$Normalized_frequencies
nf_leaf$Leaf <-
  sapply(row.names(nf_leaf),
         function(x) sum(grepl(x, row.names(nf_leaf))) == 1)
nf_leaf <- nf_leaf[nf_leaf$Leaf, c(ncol(nf_leaf), 1:(ncol(nf_leaf) - 1))]
#
# # Calculate asymmetric coincidences for source(row) - target(column)
nf_ac_st <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
                   dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
for(i in 1:(ncol(nf_leaf) - 1)){
  for(j in 1:(ncol(nf_leaf) - 1)){
    nf_ac_st[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)/sum(nf_leaf[, j + 1] > 0)
  }
}

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

