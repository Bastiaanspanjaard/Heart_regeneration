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
uniform_conversion <- F
uniform_conversion_chance <- 0.5
chances <- c(0, 0.25, 0.5, 0.75, 1)
simulation_seed <- 37 # -> 0.16 # 69 -> 0.02 # 1337 -> 0.71 top # 6 -> 0.62 top # 420 -> 0.65, z-top 
# Old list without weighted mean conversion chance: 6 -> 0.52 # 1337 -> 0.52 non-significant #420 -> 0.57 # 69 -> rho = 0.51 # 37 -> rho = 0.3

# Prepare colors and cell type data ####
cell_types <- read.csv("./Data/final_metadata.csv", stringsAsFactors = F)
colnames(cell_types)[1] <- "Cell"
# Cell type decisions as follows: all T-cell types and all macrophage types become T-cells/macrophages. 
# Split endocardium into 1 and 2, A and V but merge endocardium frzb.
cell_types$Cell_type <- cell_types$lineage.ident
cell_types$Celltype[grepl("T-cell", cell_types$Cell_type)] <- "T-cells"
cell_types$Cell_type[grepl("Macrophage", cell_types$Cell_type)] <- "Macrophages"
cell_types$Cell_type[grepl("Endocardium", cell_types$Cell_type)] <- 
  cell_types$final.zoom[grepl("Endocardium", cell_types$Cell_type)]
cell_types$Cell_type[grepl("Endocardium frzb", cell_types$Cell_type)] <- "Endocardium (frzb)"
cell_types$Cell_name <- paste("nd", cell_types$Cell, sep = "")

celltypes <- data.frame(table(cell_types$Cell_type), stringsAsFactors = F)
colnames(celltypes)[1] <- c("Cell_type")
celltypes$Cell_type <- as.character(celltypes$Cell_type)

celltype_colors <- readRDS("./Data/Cell_type_colors.rds")
ann_colors <- setNames(celltype_colors$color, celltype_colors$Cell.type)

# Load tree objects, append cell types, create lineage trees ####
tree_list <- list()
# lib <- unique(libraries$Sample)[1]
for(lib in unique(libraries$Sample)){
  tree.in <- ReadTree(lib, reference_set = cell_types[cell_types$orig.ident %in% libraries$Library_name[libraries$Sample == lib],
                                                      ], tree_path = tree_path)
  tree.in <- list(metadata = list(dpi = libraries$Dpi[libraries$Sample == lib],
                               Name = lib), 
               Tree = tree.in$Tree)
  tree_list[[length(tree_list) + 1]] <- tree.in
  names(tree_list)[[length(tree_list)]] <- lib
}
tree_list_real <- tree_list

# Simulation of target/source process ####
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
if(!(fictional_target %in% celltypes$Cell_type)){
  celltypes <- rbind(celltypes, c(fictional_target, 0))
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
}

# Create tree visualization and zoom visualization - only run if needed because this takes quite some time.
# for(t in 1:length(tree_list)){
#   tree_list[[t]]$Pie_tree <- 
#     collapsibleTree(df = tree_list[[t]]$Tree, root =tree_list[[t]]$Tree$scar, pieNode = T,
#                     pieSummary = T, collapsed = F,
#                     width = 1000, height = 500,
#                     ctypes = celltype_colors$Cell.type, linkLength=60, 
#                     ct_colors = celltype_colors$color, angle = pi/2,fontSize = 0,
#                     nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 
#   # tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
#   #                               ct_colors = type_colors$colo1)
#   # htmlwidgets::saveWidget(
#   #   tree_list[[t]]$Full_tree,
#   #   file = paste("~/Documents/Projects/heart_Bo/Images/tree_",
#   #                names(tree_list)[t], "_LINNAEUS_pie.html", sep = ""))
#   # htmlwidgets::saveWidget(
#   #   tree_list[[t]]$Fibrozoom_tree,
#   #   file = paste("~/Documents/Projects/heart_Bo/Images/tree_",
#   #                names(tree_list)[t], "_LINNAEUS_pie_fibrozoom.html", sep = ""))
# }

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
                                                Type_count = integer()),
                        Node_sizes = data.frame(Size = integer()),
                        Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = nrow(celltypes))))
colnames(comparison_list$Normalized_frequencies) <- celltypes$Cell_type

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
  
  comparison_list$Node_sizes <- rbind(comparison_list$Node_sizes,
                                      data.frame(Size = rowSums(sample_type_cumulative)))
  comparison_tree <-
    data.frame(Tree = rep(names(tree_list)[t], ncol(sample_type_count)),
               Precursor = names(sample_type_count),
               Type_count = colSums(sample_type_count))
  
  comparison_list$Comparison <- rbind(comparison_list$Comparison, comparison_tree)
  
  sample_type_nf[setdiff(names(comparison_list$Normalized_frequencies), names(sample_type_nf))] <- NA
  comparison_list$Normalized_frequencies[setdiff(names(sample_type_nf), 
                                                 names(comparison_list$Normalized_frequencies))] <- NA

  comparison_list$Normalized_frequencies <- 
    rbind(comparison_list$Normalized_frequencies, sample_type_nf)
}
comparison_list$Normalized_frequencies[is.na(comparison_list$Normalized_frequencies)] <- 0
comparison_list$Normalized_frequencies <- 
  comparison_list$Normalized_frequencies[, colSums(comparison_list$Normalized_frequencies) > 0]
comparison_list$Node_sizes$Weight <- comparison_list$Node_sizes$Size/sum(comparison_list$Node_sizes$Size)
comparison_list$All_trees_prec <- list()
comparison_list$All_trees_distances <- 
  data.frame(Precursor = names(comparison_list$Normalized_frequencies),
             Weighted_cor_progpos = numeric(ncol(comparison_list$Normalized_frequencies)),
             Weighted_cor_se = numeric(ncol(comparison_list$Normalized_frequencies)))#,
# Add included trees
comparison_list$All_trees_distances[, unique(analysis_stats$Tree[analysis_stats$Included])] <- 0

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
#     scale_fill_manual(values = ann_colors$Celltype) +
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
    scale_fill_manual(values = ann_colors) +
    scale_alpha(range = c(0.35, 1), guide = F) +
    labs(title = paste(target, ""), y = "Correlation", fill = "Source") +
    theme(#legend.position = "none",
          text = element_text(size = 36),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
)
# dev.off()

# precursor_ranking_nonrandom <- precursor_ranking

