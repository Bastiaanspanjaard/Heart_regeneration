# Dependencies ####
source("./Scripts/HR_library.R")

# Parameters ####
tree_path <- "./Data/Trees/Sbcf/"
target <- "Fibroblast (nppc)"
# target <- "Fibroblast (col11a1a)"
# target <- "Fibroblast-like cells"
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

celltypes <- data.frame(table(cell_types$Cell_type))
colnames(celltypes)[1] <- c("Cell_type")

celltype_colors <- readRDS("./Data/Cell_type_colors.rds")
ann_colors <-
  list(Celltype = setNames(celltype_colors$color, celltype_colors$Cell_type))

# Load tree objects, append cell types, create lineage trees ####


# Control
tree_list <- list()
# lib <- unique(libraries$Sample)[1]
for(lib in unique(libraries$Sample)){
  tree.in <- ReadTree(lib, reference_set = cell_types[cell_types$orig.ident %in% libraries$Library_name[libraries$Sample == lib],
                                                      ], tree_path = tree_path)
  tree.in <- list(metadata = list(dpi = libraries$Dpi[libraries$Sample == lib],
                               Name = lib), 
               Tree = tree.in$Tree) #not quite yet but close
  tree_list[[length(tree_list) + 1]] <- tree.in
  names(tree_list)[[length(tree_list)]] <- lib
}

# H5 <- ReadTree("H5", reference_set = cell_types[cell_types$orig.ident == "H5", ], tree_path = tree_path)
# H6 <- ReadTree("H6", reference_set = cell_types[cell_types$orig.ident == "H6", ], tree_path = tree_path)
# H7 <- ReadTree("H7", reference_set = cell_types[cell_types$orig.ident == "H7", ], tree_path = tree_path)
# H8 <- ReadTree("H8", reference_set = cell_types[cell_types$orig.ident == "H8", ], tree_path = tree_path)
# H8 <- ReadTree("H8", reference_set = cell_types[cell_types$orig.ident %in% c("H8a", "H8v"), ], tree_path = tree_path)
# tree_list <- list(H5 = H5, H6 = H6, H7 = H7)
# 3dpi
# Hr10 <- ReadTree("Hr10", reference_set = cell_types[cell_types$orig.ident == "Hr10", ], tree_path = tree_path)
# Hr11 <- ReadTree("Hr11", reference_set = cell_types[cell_types$orig.ident == "Hr11", ], tree_path = tree_path)
# Hr12 <- ReadTree("Hr12", reference_set = cell_types[cell_types$orig.ident == "Hr12", ], tree_path = tree_path)
# Hr22 <- ReadTree("Hr22", reference_set = cell_types[cell_types$orig.ident == "Hr22", ], tree_path = tree_path)
# Hr23 <- ReadTree("Hr23", reference_set = cell_types[cell_types$orig.ident == "Hr23", ], tree_path = tree_path)
# Hr24 <- ReadTree("Hr24", reference_set = cell_types[cell_types$orig.ident == "Hr24", ], tree_path = tree_path)
# Hr25 <- ReadTree("Hr25", reference_set = cell_types[cell_types$orig.ident == "Hr25", ], tree_path = tree_path)
# Hr26 <- ReadTree("Hr26", reference_set = cell_types[cell_types$orig.ident == "Hr26", ], tree_path = tree_path)
# Hr27 <- ReadTree("Hr27", reference_set = cell_types[cell_types$orig.ident == "Hr27", ], tree_path = tree_path)
# tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr24 = Hr24, Hr26 = Hr26, Hr27 = Hr27)
# 7dpi
# Hr1_ref <- cell_types[cell_types$orig.ident == "Hr1", ]
# Hr1_ref$Cell_name <- paste("nd", sapply(Hr1_ref$Cell_name, function(x){unlist(strsplit(x, "_"))[2]}), sep = "")
# Hr1 <- ReadTree("Hr1", reference_set = Hr1_ref, tree_path = tree_path)
# Hr1 <- ReadTree("Hr1", reference_set = cell_types[cell_types$orig.ident == "Hr1", ], tree_path = tree_path)
# Hr2 <- ReadTree("Hr2", reference_set = cell_types[cell_types$orig.ident %in% c("Hr2a", "Hr2b"), ], tree_path = tree_path)
# Hr6 <- ReadTree("Hr6", reference_set = cell_types[cell_types$orig.ident %in% c("Hr6a", "Hr6v"), ], tree_path = tree_path)
# Hr7 <- ReadTree("Hr7", reference_set = cell_types[cell_types$orig.ident %in% c("Hr7a", "Hr7v"), ], tree_path = tree_path)
# Hr13 <- ReadTree("Hr13", reference_set = cell_types[cell_types$orig.ident == "Hr13", ], tree_path = tree_path)
# Hr14 <- ReadTree("Hr14", reference_set = cell_types[cell_types$orig.ident == "Hr14", ], tree_path = tree_path)
# Hr15 <- ReadTree("Hr15", reference_set = cell_types[cell_types$orig.ident == "Hr15", ], tree_path = tree_path)
# tree_list <- list(Hr1 = Hr1, Hr2 = Hr2, Hr13 = Hr13, Hr14 = Hr14, Hr15 = Hr15)
# 30dpi
# Hr3 <- ReadTree("Hr3", reference_set = cell_types[cell_types$orig.ident == "Hr3", ], tree_path = tree_path)
# Hr4 <- ReadTree("Hr4", reference_set = cell_types[cell_types$orig.ident == "Hr4", ], tree_path = tree_path)
# Hr19 <- ReadTree("Hr19", reference_set = cell_types[cell_types$orig.ident == "Hr19", ], tree_path = tree_path)
# Hr20 <- ReadTree("Hr20", reference_set = cell_types[cell_types$orig.ident == "Hr20", ], tree_path = tree_path)
# Hr21 <- ReadTree("Hr21", reference_set = cell_types[cell_types$orig.ident == "Hr21", ], tree_path = tree_path)

# 3 and 7dpi
# tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr22 = Hr22, 
#                   Hr23 = Hr23, Hr24 = Hr24, Hr25 = Hr25, Hr26 = Hr26, Hr27 = Hr27,
#                   Hr1 = Hr1, Hr2 = Hr2, Hr6 = Hr6, Hr7 = Hr7, Hr13 = Hr13, #Hr14 = Hr14, 
#                   Hr15 = Hr15)
# tree_list <- list(Hr2 = Hr2)

# Create tree visualization and zoom visualization - only run if needed because this takes quite some time.
for(t in 1:length(tree_list)){
  # tree_list[[t]] <- MakePieTree(tree_list[[t]], "Full_tree", 
  #                               ctypes = celltype_colors$Cell.type, ct_colors = celltype_colors$color)
  tree_list[[t]]$Pie_tree <- 
    collapsibleTree(df = tree_list[[t]]$Tree, root =tree_list[[t]]$Tree$scar, pieNode = T,
                    pieSummary = T, collapsed = F,
                    width = 1000, height = 500,
                    ctypes = celltype_colors$Cell.type, linkLength=60, 
                    ct_colors = celltype_colors$color, angle = pi/2,fontSize = 0,
                    nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 
  # tree_sample$CTree <-    
  #   collapsibleTree(df = tree_to_plot, root = tree_to_plot$scar, pieNode = T,
  #                   pieSummary = T, collapsed = F,
  #                   width = 1000, height = 500,
  #                   ctypes = ctypes,linkLength=60, 
  #                   ct_colors = ct_colors, angle = pi/2,fontSize = 0,
  #                   nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 
  # tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
  #                               ct_colors = type_colors$colo1)
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
                                                Type_count = integer()),
                        Node_sizes = data.frame(Size = integer()),
                        Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = length(unique(cell_types$Cell_type)))))
colnames(comparison_list$Normalized_frequencies) <- unique(cell_types$Cell_type)

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
# dev.off()

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

