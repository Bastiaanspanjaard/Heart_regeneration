# Dependencies ####
source("./Scripts/HR_library.R")

# Parameters ####
tree_path <- "./Data/Trees/Sbcf/"
# target <- "Harpocytes"
# target <- "Fibroblast (nppc)"
target <- "Fibroblast (col11a1a)"
# target <- "Fibroblast-like cells"
# target <- "Perivascular cells"
# target <- "Fibroblast (proliferating)"
# target <- "Fibroblast (col12a1a)"
# target <- "Ventricle (ttn.2/aSMA)"
precursors_include_incl_target <- NULL
  # c("Fibroblast (cxcl12a)", "Fibroblast", "Epicardium V",
  #   "Perivascular cells",
  #   "Fibroblast (col11a1a)", "Epicardium A", "Fibroblast (col12a1a)")
  # c("Fibroblast-like cells", "Endocardium", "Fibroblast (spock3)")
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

# Simulation of target/source process ####
# Uses as input: tree_list_real, list of trees, fictional source and target, simulation seed.
# Way of simulating conversion (+chances), celltypes
# Creates as output: tree list with modified cell types (also in edge list and abundancies!)

# Split a cell type in two to see if these show up as related
# Pick a set of trees (Hr1, Hr2, Hr6?).
# tree_list_subset <- tree_list[c("Hr1", "Hr2", "Hr6")]
tree_list <- list(Hr1 = list(Tree = Clone(tree_list_real$Hr1$Tree)), 
                  Hr2 = list(Tree = Clone(tree_list_real$Hr2$Tree)), 
                  Hr6 = list(Tree = Clone(tree_list_real$Hr6$Tree)))
# Calculate cell type numbers per tree - Endocardium (V) seems like a good split candidate.
# aggregate(sample_tree_list$Hr6$Pure_counts$Type_count,
#           by = list(Cell_type = sample_tree_list$Hr6$Pure_counts$Cell_type),
#           sum)
fictional_source <- "Endocardium 1 (V)"
fictional_target <- "Harpocytes"
if(!(fictional_target %in% celltype_frequencies$Cell_type)){
  celltype_frequencies <- rbind(celltype_frequencies, c(fictional_target, 0))
}

# Extract all endocardium cells, randomly rename half of them, put cells "back in tree".
set.seed(simulation_seed)
for(t in 1:length(tree_list)){
  edge_list <- ToDataFrameNetwork(tree_list[[t]]$Tree, "Cell.type")
  node_chances <- data.frame(table(edge_list$from[edge_list$Cell.type != "NA"]))
  colnames(node_chances) <- c("Node", "Count")

  # node_chances <- data.frame(Node = unique(edge_list$from))
  if(uniform_conversion){
    # For uniform chance for all cells
    node_chances$Conversion_chance <- uniform_conversion_chance
  }else{
    # For separate chances per node
    node_chances$Leaf <- sapply(node_chances$Node, function(x) sum(Progenyl(x, node_chances$Node)) == 0)
    # Assign chances to leaf nodes
    node_chances$Conversion_chance <- 
      sapply(node_chances$Leaf, 
             function(x){
               if(x){
                 sample(chances, 1)
               }else(return(NA))
             }
      )  
    # Cumulate chances for non-leaf nodes. For leaf nodes, below formula returns the conversion chance of the
    # leaf. The na.rm is necessary because the node itself (if not a leaf) will have NA as conversion
    # probability.
    # NOTE this is scaled by cell number per node, and actually be done over multiple rounds
    # to take into account that some nodes may have more cells than suggested by the cells in their
    # progeny nodes. Formula below is the previous formula.
    # node_chances$Conversion_chance <-
    #   sapply(node_chances$Node,
    #          function(x) mean(node_chances$Conversion_chance[grepl(x, node_chances$Node)], na.rm = T))
    while(sum(is.na(node_chances$Conversion_chance))){
    node_chances$Conversion_chance <-
      sapply(node_chances$Node,
             function(x) {
               node_progeny_bool <- grepl(x, node_chances$Node)
               # Find nodes that do not yet have a conversion chance while all progeny does.
               if(sum(is.na(node_chances$Conversion_chance[node_progeny_bool])) == 1){
                 # Calculate the weighted mean of progeny nodes. Note that this will include
                 # intermediate and leaf nodes in the weighting, and also note that this is
                 # correct.
                 return(
                   weighted.mean(x = node_chances$Conversion_chance[node_progeny_bool],
                                 w = node_chances$Count[node_progeny_bool],
                                 na.rm = T)
                 )
               }else{
                 return(node_chances$Conversion_chance[node_chances$Node == x])
               }
             }
      )
    }
  }
  
  # Transform cells
  new_cell_types <- edge_list[edge_list$Cell.type != "NA", ]
  new_cell_types <- merge(new_cell_types, node_chances, by.x = "from", by.y = "Node")
  new_cell_types$Roll <- runif(sum(nrow(new_cell_types)))
  new_cell_types$New_cell_type <- 
    apply(new_cell_types[, c("Cell.type", "Conversion_chance", "Roll")], 1,
          function(x){
            if(x[1] == fictional_source){
              if(as.numeric(x[3]) < as.numeric(x[2])){
                return(fictional_target)
              }else{
                return(fictional_source)
              }
            }else{
              return(x[1])
            }
          })
  tree_list[[t]]$Tree$Do(function(node) {node$Cell.type = new_cell_types$New_cell_type[new_cell_types$to == node$name]}, filterFun = isLeaf)
  # Also chance edge lists and cell type abundancies.
}

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

# Finding precursor candidates ####
# Start by calculating the correlations between a progenitor cell type and potential
# precursors. After this, investigate whether these correlations are stable when leaving
# out one of the trees. Finally, bootstrap the correlations by randomly placing the
# precursor cells over the trees to see what the null-distribution of the correlations
# is and how enriched the observed correlations are.

# ** Calculate correlations between cell types over all trees ####
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
  # c("B-cells", "Bl.ves.EC (apnln)", "Bl.ves.EC (lyve1)", "Bl.ves.EC (plvapb)", 
  #   "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2)", "Cardiomyocytes (ttn.2)",
  #   "Cardiomyocytes A", "Cardiomyocytes V", "Endocardium (frzb)", "Endocardium (A)", 
  #   "Endocardium (V)", "Epicardium (Atrium)", "Epicardium (Ventricle)",
  #   "Fibroblast", "Fibroblast (cfd)", "Fibroblast (col11a1a)", "Fibroblast (col12a1a)",
  #   "Fibroblast (cxcl12a)", "Fibroblast (mpeg1.1)", "Fibroblast (nppc)", "Fibroblast (proliferating)",
  #   "Fibroblast (spock3)", "Fibroblast-like cells", "Macrophages", "Myelin cells",
  #   "Neuronal cells", "Neutrophils", "Perivascular cells", "Proliferating cells", "Smooth muscle cells",           
  #   "T-cells", "Monocytes")
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
    # NOTE: we may need this for the target-based correlations
    # if(!(target %in% edge_list$Cell.type) | sum(edge_list$Cell.type == target) < 20){
    #   analysis_stats_add$Included <- F
    #   analysis_stats <- rbind(analysis_stats, analysis_stats_add)
    #   next
    # }
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

  # Collapse non-splitting branches:
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

  # Optional plot of correlations
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
  # ggplot(df_ct[df_ct$Tree %in% c("Hr1", "Hr2", "Hr6"), ]) +
  #   geom_point(aes(x = CT1, y = CT2, color = Tree, size = Node_size)) +
  #   labs(x = cor_type_1, y = cor_type_2, 
  #        title = paste("7dpi node frequencies deep injuries, cor ", 
  #                      round(cor(df_ct$CT1[df_ct$Tree %in% c("Hr1", "Hr2", "Hr6")], 
  #                                df_ct$CT2[df_ct$Tree %in% c("Hr1", "Hr2", "Hr6")]), 2), sep = ""))
  # # cor(df_ct$CT1, df_ct$CT2)
  # # Leaf nodes:
  # df_ct$Leaf <-
  #   sapply(df_ct$Node,
  #          function(x) sum(grepl(x, df_ct$Node)) == 1)
  # # trees_in <- c("Hr1", "Hr2", "Hr6")
  # df_ctl <- df_ct[df_ct$Leaf, ]
  # # Symmetrical coincidence
  # sum(df_ctl$CT1 > 0 & df_ctl$CT2 > 0)/sum(df_ctl$CT1 > 0 | df_ctl$CT2 > 0)
  # # Asymmetrical coincidence (removes shallow injury)
  # sum(df_ctl$CT1 > 0 & df_ctl$CT2 > 0)/sum(df_ctl$CT1 > 0)
  # 
  # nf_leaf <- comparison_list$Normalized_frequencies
  # nf_leaf$Leaf <-
  #   sapply(row.names(nf_leaf),
  #          function(x) sum(grepl(x, row.names(nf_leaf))) == 1)
  # nf_leaf <- nf_leaf[nf_leaf$Leaf, c(ncol(nf_leaf), 1:(ncol(nf_leaf) - 1))]
  #
  # # Calculate asymmetric coincidences for source(row) - target(column)
  # nf_ac_st <- matrix(NA, nrow = ncol(nf_leaf) - 1, ncol = ncol(nf_leaf) - 1,
  #                    dimnames = list(colnames(nf_leaf)[-1], colnames(nf_leaf)[-1]))
  # for(i in 1:(ncol(nf_leaf) - 1)){
  #   for(j in 1:(ncol(nf_leaf) - 1)){
  #     nf_ac_st[i,j] <- sum(nf_leaf[, i + 1] > 0 & nf_leaf[, j + 1] > 0)/sum(nf_leaf[, j + 1] > 0)
  #   }
  # }
  # Number of occupied leaves:
  # colSums(nf_leaf[, -1] > 0)
  
  # start_targets <- c("Fibroblast (spock3)","Fibroblast (nppc)", "Valve fibroblasts")
  # current_targets <- start_targets
  # nf_ac_current <- nf_ac_st[, current_targets, drop = F]
  # current_targets <- names(which(rowSums(nf_ac_current >= 0.8) > 0))
  # OLD METHOD
  # new_targets <- names(which(rowSums(nf_ac_current >= 0.8) > 0))
  # while(!setequal(new_targets, current_targets)){
  #   current_targets <- new_targets
  #   nf_ac_current <- nf_ac_st[, current_targets, drop = F]
  #   new_targets <- names(which(rowSums(nf_ac_current >= 0.8) > 0))
  # }
  # END OLD METHOD
  # ac_tograph <- nf_ac_st[current_targets, current_targets]
  # nodenames <- 1:nrow(ac_tograph)
  # longnames <- rownames(ac_tograph)

  # graph1<-qgraph(ac_tograph, diag = F,
  #                layout= "spring", edge.color = "black",
  #                labels = F, #nodenames, nodeNames = longnames, 
  #                minimum = 0.79999, cut = 0.8, legend = F,
  #                color = as.character(celltype_colors[rownames(ac_tograph)]),
  #                borders = F, 
  #                filename = "~/Documents/Projects/heart_Bo/Images/Asymmetric_coinc_network_7dpi_to_Fibnppcspock3fiblike",
  #                filetype = "eps",
  #                vsize=15, esize = 20, asize = 9, width = 1.5, height = 1.5)
  # Plot for arrow thickness legends
  # qgraph_arrow_legend <- ac_tograph
  # qgraph_arrow_legend[,] <- 0
  # diag(qgraph_arrow_legend) <- 1
  # qgraph_arrow_legend[2, 1] <- 0.8
  # qgraph_arrow_legend[3, 1] <- 0.9
  # qgraph_arrow_legend[4, 1] <- 1
  # graph1<-qgraph(qgraph_arrow_legend, diag = F,
  #                layout= "spring", edge.color = "black",
  #                labels = F, #nodenames, nodeNames = longnames, 
  #                minimum = 0.79999, cut = 0.8, legend = F,
  #                # color = as.character(celltype_colors[rownames(ac_tograph)]),
  #                borders = F, 
  #                filename = "~/Documents/Projects/heart_Bo/Images/coinc_network_legend",
  #                filetype = "eps",
  #                vsize=15, esize = 20, asize = 9, width = 1.5, height = 1.5)
  
  # Can the number of nodes be explained by the number of cells?
  # Test against random distribution model, distributing cells evenly
  # over nodes, weighted by node size of course.
  # size for each leaf node, number of cells per type in tree-specific leaf nodes.
  # Probability of #nodes is the probability of non-occupancy (i.e. occupancy 0)
  # f_leaf <- comparison_list$Frequencies
  # f_leaf$Leaf <-
  #   sapply(row.names(f_leaf),
  #          function(x) sum(grepl(x, row.names(f_leaf))) == 1)
  # f_leaf <- f_leaf[f_leaf$Leaf, c(ncol(f_leaf), 1:(ncol(f_leaf) - 1))]
  # f_leaf$Size <- rowSums(f_leaf[, -1])
  # f_leaf$Tree <- sapply(rownames(f_leaf), function(x) unlist(strsplit(x, ":"))[1])
  # f_leaf$Node <- row.names(f_leaf)
  # f_leaf <- f_leaf[, c(1, ncol(f_leaf) - 1, ncol(f_leaf), 2:(ncol(f_leaf) - 2))]
  
  # Cell type, number of nodes occupied, total nodes, min p-value, prod p-value
  # cell_type_occupancy_p <-
  #   data.frame(Cell_type = colnames(comparison_list$Frequencies),
  #              Occupied_nodes = colSums(comparison_list$Frequencies > 0),
  #              Min_p = NA,
  #              Prod_p = NA)
  # cell_type_occupancy_p$Min_p[1] <- 3
  
  # for(c in 1:nrow(cell_type_occupancy_p)){
  #   dist_type <- as.character(cell_type_occupancy_p$Cell_type[c]) #"Fibroblast (nppc)"
  #   dist_type_leaves <- f_leaf[, c("Leaf", "Size", "Tree", "Node", dist_type)]
  #   dist_type_leaves <- 
  #     merge(
  #       dist_type_leaves,
  #       merge(aggregate(dist_type_leaves[, 5],
  #                       by = list(Tree = dist_type_leaves$Tree),
  #                       sum),
  #             aggregate(dist_type_leaves$Size,
  #                       by = list(Tree = dist_type_leaves$Tree),
  #                       sum), by = "Tree"))
  #   colnames(dist_type_leaves)[5:7] <- c("Type_count", "Total_type_count", "Total_tree_count")
  #   dist_type_leaves$P_count_lowertail <-
  #     apply(dist_type_leaves[, c("Type_count", "Size", "Total_type_count", "Total_tree_count")], 1,
  #           function(x) pbinom(q = x[1], size = x[2], prob = x[3]/x[4], lower.tail = TRUE))
  #   cell_type_occupancy_p$Min_p[c] <-
  #     min(dist_type_leaves$P_count_lowertail[dist_type_leaves$Type_count == 0])
  #   cell_type_occupancy_p$Prod_p[c] <-
  #     prod(dist_type_leaves$P_count_lowertail[dist_type_leaves$Type_count == 0])
  # }
  
  # Calculate correlations
  comparison_list$All_trees_distances <- 
    data.frame(Precursor = names(comparison_list$Normalized_frequencies),
               Weighted_cor_progpos = numeric(ncol(comparison_list$Normalized_frequencies)),
               Weighted_cor_se = numeric(ncol(comparison_list$Normalized_frequencies)))#,
  comparison_list$All_trees_distances[, unique(analysis_stats$Tree[analysis_stats$Included])] <- 0
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
  # pheatmap::pheatmap(t(norm_freq), cluster_cols = F, cluster_rows = cell_type_cor_clust,
  #                    cutree_rows = cluster_number, gaps_col = gaps_row,
  #                    border_color = NA, breaks = 0:100/100,
  #                    show_colnames = F)
  # dev.off()
  
  # pheatmap((cor_present + 1)/2, cutree_cols = cluster_number,
  #          cluster_cols = cell_type_cor_clust)
  # pdf("./Images/Cell_type_cell_type_relations_pheatmap_7dpi_min10_oneEndo.pdf",
  #     width = 10)
  # pheatmap(cor_present,
  #          cluster_cols = cell_type_cor_clust,
  #          cluster_rows = cell_type_cor_clust,
  #          cutree_cols = cluster_number, cutree_rows = cluster_number,
  #          treeheight_col = 0, show_colnames = F, fontsize = 16)
  # dev.off()
  # plot(cell_type_cor_clust)
  
  # For resampling: record which cell types are together in a cluster and 
  this_ds_hclust <- data.frame(cutree(cell_type_cor_clust, k = cluster_number))
  colnames(this_ds_hclust)[1] <- "DS" #paste("DS", sample_fraction, iter, sep = "_")
  this_ds_hclust$Cell_type <- factor(row.names(this_ds_hclust), levels = included_types)
  this_ds_hclust$Presence <- 1
  type_clusters <- reshape2::acast(this_ds_hclust, Cell_type ~ DS, value.var = "Presence", fill = 0, drop = F)
  type_clusters_add <- type_clusters %*% t(type_clusters)
  type_clusters_add <- type_clusters_add[included_types, included_types]
  cell_type_correspondence <- cell_type_correspondence + type_clusters_add
  
  # START Making an example correlation plot
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

hclust_grouping <- data.frame(cutree(cell_type_cor_clust, k = cluster_number))
colnames(hclust_grouping)[1] <- "Hclust"
sum(rownames(hclust_grouping) != rownames(cor_present)) # Just to check if the hclust ordering is correct
hclust_grouping$Nodename <- 1:nrow(hclust_grouping)
hclust_list <- list()
for(hc in unique(hclust_grouping$Hclust)){
  hclust_list[[length(hclust_list) + 1]] <- hclust_grouping$Nodename[hclust_grouping$Hclust == hc]
  names(hclust_list)[[length(hclust_list)]] <- hc
}
nodenames <- 1:nrow(cor_present)
longnames <- rownames(cor_present)

graph1<-qgraph(cor_present, diag = F,
               layout= "spring",
               labels = nodenames, nodeNames = longnames,
               sampleSize = nrow(cor_present),
               groups=hclust_list, minimum = 0.2, 
               vsize=5, esize = 4.5,
               filename = "~/Documents/Projects/heart_Bo/Images/Cor_network_7dpi_forceinclude_min10_article_legendonly", filetype = "eps",
               border.width = 2,#cut=0, maximum=.45, border.width=1.5,
               width = 5.5, height = 5.5, legend = T,
               legend.cex = 0.2, label.cex = 2)


for(p in 1:ncol(comparison_list$Normalized_frequencies)){
  precursor <- names(comparison_list$Normalized_frequencies)[p]
  
  if(precursor == target){
    comparison_list$All_trees_distances$Weighted_cor_progpos[p] <- -2
    next
  }

  # Select nodes with target cell type (turn off for symmetric).
  progpos <- comparison_list$Normalized_frequencies[comparison_list$Normalized_frequencies[[target]] > 0, ]
  progpos_weights <- comparison_list$Node_sizes[comparison_list$Normalized_frequencies[[target]] > 0, ]
  weighted_cor <- data.frame(wtd.cor(progpos[[precursor]],
                                     progpos[[target]],
                                     weight = progpos_weights$Size))#, bootse = T))
  comparison_list$All_trees_distances$Weighted_cor_progpos[p] <- weighted_cor$correlation
  comparison_list$All_trees_distances$Weighted_cor_se[p] <- weighted_cor$std.err
  
  # One-out correlations: leave out one tree and recalculate the correlations.
  for(t in unique(analysis_stats$Tree[analysis_stats$Included])){
    progpos_oneout <- progpos[!grepl(paste(t, ":", sep = ""), rownames(progpos)), ]
    weights_oneout <- progpos_weights[!grepl(paste(t, ":", sep = ""), rownames(progpos_weights)), ]
    comparison_list$All_trees_distances[p, t] <-
      wtd.cors(progpos_oneout[[precursor]],
             progpos_oneout[[target]],
             weight = weights_oneout$Size)
    }
}  

# Add #trees that were used, #trees that had precursor candidate cell type, average
# cells per tree used.
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
#     scale_fill_manual(values = celltype_colors) +
#     labs(title = paste(target, ""), y = "Fitness") +
#     theme(legend.position = "none",
#           axis.ticks.x = element_blank(),
#           axis.text.x = element_blank())
# )
# dev.off()

# Plot correlations of two top-scoring precursors and two low-scoring precursors
# for(i in c(1:2, nrow(precursor_ranking) - 1, nrow(precursor_ranking))){
#   pot_prec <- as.character(precursor_ranking$Precursor[i])
#   plot_freqs <-
#     data.frame(comparison_list$Normalized_frequencies)[,make.names(c(pot_prec, target))]
#   colnames(plot_freqs) <- c("Precursor", "Progenitor")
#   plot_freqs$Size <- comparison_list$Node_sizes$Size
  # png(paste("./Images/", pot_prec, "precursor_to_", target, "_37dpi_all.png", sep = ""))
  # print(
    # ggplot() +
    #   geom_point(data = plot_freqs[plot_freqs$Progenitor == 0, ],
    #              aes(x = Precursor, y = Progenitor), size = 6) + #, size = Size)) + #, color = "grey") +
    #   geom_point(data = plot_freqs[plot_freqs$Progenitor > 0, ],
    #              aes(x = Precursor, y = Progenitor), size = 6) + #, size = Size)) +
    #   labs(x = pot_prec, y = target) +
    #   theme(text = element_text(size = 36))#+
      # theme(legend.position = "none")
  # )
  # dev.off()
# }

# ** Bootstrap correlation ####
# Input: tree list, but we could do this with edge lists and node counts as well.
# Again, the first loop calculates tree statistics.
# Second loop runs over all potential precursors to randomize each of them 
# Output: for the second loop, a dataframe of calculated correlations, and then
# afterwards, quantiles etc for the correlations.

# If I randomly place cells in the tree, what is the correlation with the progeny?
# Select a potential precursor; how many cells in each tree? Randomly place that amount
# of cells in the tree, depending on 'pure' (i.e. not cumulative) node counts when progeny
# is removed. Then calculate correlation with progeny.
# Total node size after sampling is node size before sampling without progeny, with added
# sampled progeny. This means precursor frequencies and node weights change as well.
# Do note that for proper sampling, we need to use the 'pure' (i.e. non-cumulative) node
# sizes, but for correlations we do need to talk about cumulative nodes.

# Loop over trees to calculate node counts per cell type
# Note this looks a little weird since it also calculates for trees that have too few target cells.
# But those trees are not used in the actual bootstrapping.
sample_tree_list <- list()
for(t in 1:length(tree_list)){
  print(paste("Tree", names(tree_list)[t]))
  # Pure node counts with and without progeny
  tree <- Clone(tree_list[[t]]$Tree) # Clone(tree_list$Hr26$Tree)
  
  nodes <- tree$Get('name', filterFun = function(x) {!isLeaf(x)})
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  pure_counts <- data.frame(table(edge_list$from, edge_list$Cell.type))
  colnames(pure_counts) <- c("Node", "Cell_type", "Type_count")
  pure_counts <- pure_counts[pure_counts$Cell_type != "NA", ]
  
  sample_tree_list[[t]] <- list(Pure_counts = pure_counts)
  names(sample_tree_list)[t] <- names(tree_list)[t]
}

# Perform sampling and calculate correlations
samples <- 1000
sampled_wcors <- data.frame(matrix(nrow = samples, ncol = nrow(precursor_ranking) - 1))
colnames(sampled_wcors) <- precursor_ranking$Precursor[-nrow(precursor_ranking)]

set.seed(42)
bootstrap_list <- list()
for(p in 1:(nrow(precursor_ranking) - 1)){
  bootstrap_precursor <- as.character(precursor_ranking$Precursor[p])
  print(paste("Potential precursor:", bootstrap_precursor))
  
  # Randomly place cells in tree; calculate resulting cumulative precursor frequencies
  all_tree_sample_frequencies <- matrix(nrow = 0, ncol = samples)
  all_tree_progeny_frequencies <- matrix(nrow = 0, ncol = samples)
  all_tree_sampled_node_sizes <- matrix(nrow = 0, ncol = samples)
  for(t in 1:length(sample_tree_list)){
    pure_counts <- sample_tree_list[[t]]$Pure_counts
    if(!(target %in% pure_counts$Cell_type) | sum(pure_counts$Type_count[pure_counts$Cell_type == target]) < 20){
      next
    }
    
    pure_sizes <- aggregate(pure_counts$Type_count,
                            by = list(Node = pure_counts$Node),
                            sum)
    colnames(pure_sizes)[2] <- "Node_count"
    pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == target, ])
    colnames(pure_sizes)[3:4] <- c("Progeny", "Progeny_count")
    
    # Potential precursor and how many cells?
    bootstrap_cells <- sum(pure_counts$Type_count[pure_counts$Cell_type == bootstrap_precursor])
    
    if(bootstrap_cells == 0){
      pure_sizes$Precursor <- bootstrap_precursor
      pure_sizes$Precursor_count <- 0
    }else{
      # Pure node sizes with and without precursor
      pure_sizes <- merge(pure_sizes, pure_counts[pure_counts$Cell_type == bootstrap_precursor, ])
      colnames(pure_sizes)[5:6] <- c("Precursor", "Precursor_count")
    }
    pure_sizes$Sample_chance <- 
      (pure_sizes$Node_count - pure_sizes$Progeny_count)/(sum(pure_sizes$Node_count) - sum(pure_sizes$Progeny_count))
    pure_sizes$Size_no_prec <- pure_sizes$Node_count - pure_sizes$Precursor_count
    
    for(n in 1:nrow(pure_sizes)){
      node <- as.character(pure_sizes$Node[n])
      pure_sizes$CSize_no_prec[n] <- sum(pure_sizes$Size_no_prec[grep(node, pure_sizes$Node)])
      pure_sizes$CProgeny_count[n] <- sum(pure_sizes$Progeny_count[grep(node, pure_sizes$Node)])
    }
    
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
  
  # Calculate weighted progeny-node correlation
  for(s in 1:samples){
    sampled_wcors[s, p] <- 
      wtd.cors(all_tree_sample_frequencies[(all_tree_progeny_frequencies[, s] != 0), s], 
               all_tree_progeny_frequencies[(all_tree_progeny_frequencies[, s] != 0), s], 
               weight = all_tree_sampled_node_sizes[(all_tree_progeny_frequencies[, s] != 0), s])
  }
}

example_sample <- sampled_wcors[, 2, drop = F]
colnames(example_sample) <- "Wcor"

# Example plot
# png("./Images/Ventricular_endocardium_to_nppc_fibroblast_bootcor.png",
#     width = 600, height = 400)
# ggplot(example_sample) +
#   geom_histogram(aes(x = Wcor)) +
#   geom_vline(xintercept = mean(example_sample$Wcor), col = "red", size = 3) +
#   geom_vline(xintercept = 
#                precursor_ranking$Weighted_cor_progpos[precursor_ranking$Precursor == "Endocardium (Ventricle)"],
#              col = "blue", size = 3) +
#   labs(x = "Correlation", y = "Count") +
#   theme(text = element_text(size = 36))
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

prp_m$Eminus[prp_m$variable == "Weighted_cor_progpos"] <- NA
prp_m$Eplus[prp_m$variable == "Weighted_cor_progpos"] <- NA


# pdf(paste("./Images/", target, "_all_precursors_with_bootstrap_37dpi.pdf", sep = ""),
    #width = 1400, height = 768)
# png(paste("./Images/", target, "_all_precursors_with_bootstrap_37dpi_legend.png", sep = ""),
#     type = "quartz", width = 900, height = 768)
# png("./Images/Fibroblast_col11_sources.png",
    # width = 960, height = 480)
# png("./Images/Harpocytes_from_endo1V_all_precursors_with_bootstrap_Hr1n2n6_legend.png",
#     type = "quartz", width = 900, height = 768)
print(
  ggplot(prp_m) +
    geom_bar(aes(x = Precursor, y = value, fill = Precursor, group = variable, alpha = Type_alpha),
             stat = "identity", position = "dodge") +
    # geom_errorbar(aes(x = Precursor, ymin = Minus099, ymax = Plus099, group = variable), 
    #               position = position_dodge(width = 1), width = 0.5) +
    geom_errorbar(aes(x = Precursor, ymin = Eminus, ymax = Eplus, group = variable), 
                  position = position_dodge(width = 1), width = 0.5) +#, size = 2) +
    scale_fill_manual(values = celltype_colors) +
    scale_alpha(range = c(0.35, 1), guide = F) +
    labs(title = paste(target, ""), y = "Correlation", fill = "Source") +
    theme(#legend.position = "none",
          text = element_text(size = 36),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
)
# dev.off()

# precursor_ranking_nonrandom <- precursor_ranking

