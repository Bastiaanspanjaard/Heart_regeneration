# Description ####
# Functions for analysing and constructing single-cell trees

# Written by B. Spanjaard, 2020

# Dependencies ####
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(tidyr))
suppressPackageStartupMessages(require(collapsibleTree))
suppressPackageStartupMessages(require(data.tree))
suppressPackageStartupMessages(require(pheatmap))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(igraph))
suppressPackageStartupMessages(require(weights))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(qgraph))

# Parameters ####
libraries <- data.frame(Library_name = c("H5", "H6", "H7", "H8a", "H8v",
                                         "Hr1", "Hr2a", "Hr2b", "Hr3", "Hr4",
                                         "Hr6a", "Hr6v", "Hr7a", "Hr7v",
                                         "Hr10", "Hr11", "Hr12", "Hr13",
                                         "Hr14", "Hr15", "Hr19", "Hr20",
                                         "Hr21", "Hr22", "Hr23", "Hr24",
                                         "Hr25", "Hr26", "Hr27"),
                        Sample =  c("H5", "H6", "H7", "H8", "H8",
                                    "Hr1", "Hr2", "Hr2", "Hr3", "Hr4",
                                    "Hr6", "Hr6", "Hr7", "Hr7",
                                    "Hr10", "Hr11", "Hr12", "Hr13",
                                    "Hr14", "Hr15", "Hr19", "Hr20",
                                    "Hr21", "Hr22", "Hr23", "Hr24",
                                    "Hr25", "Hr26", "Hr27"),
                        Dpi = c("Ctrl", "Ctrl", "Ctrl", "Ctrl", "Ctrl",
                                "7", "7", "7", "30", "30",
                                "7", "7", "7", "7",
                                "3", "3", "3", "7",
                                "7", "7", "30", "30",
                                "30", "3", "3", "3",
                                "3", "3", "3"),
                        stringsAsFactors = F)
theme_bs <- 
  theme_update(text = element_text(size = 24),
               panel.background = element_rect(fill = "white", colour = "black"),
               panel.grid.major.y = element_line(colour = "grey"),
               panel.grid.major.x = element_line(colour = "grey"),
               legend.key = element_rect(fill = "white"),
               plot.title = element_text(hjust = 0.5))

# Functions ####
loadRData <- function(fileName){
  # Loads Rdata/Robj, and returns it so it can be named
  load(fileName)
  get(ls()[ls() != "fileName"])
}

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

CalculateTypeCorrelations <- function(tree_list, comparison_list, 
                                      inclusion_limit = 20, force_include = NULL,
                                      correlation_type = "Symmetric"){
  # Calculate correlations between cell types over trees. Correlations are calculated between
  # cell type frequencies over all nodes, with nodes weighted by cell count.
  # For each cell type, trees that have at least [inclusion_limit] cells of that type are
  # included in the correlation calculation.
  # Symmetric correlations include all nodes, asymmetric correlations only include the nodes
  # in which cells of the target type are present; those are the types in the row of the output.
  # Correlations between cell types that do not co-occur in nodes or trees will be returned
  # as NaN.
  norm_freq <- comparison_list$Normalized_frequencies
  node_weights <- comparison_list$Node_sizes
  
  # Inclusion/exclusion inc_exc: use to include trees for cell types. Boolean matrix with same row/column 
  # structure as norm_freq, entries F mean do not include in correlation calculation.
  inc_exc <- norm_freq # Same structure as norm_freq
  inc_exc[,] <- F # Set all entries to F
  # Loop over all trees. In [inc_exc], set all nodes per tree to true for each cell type that has more than
  # or equal to [inclusion_limit] cells.
  for(t in 1:length(tree_list)){
    tree_name <- names(tree_list)[t]
    tree_type_include <- aggregate(tree_list[[t]]$Node_type_counts$Type_count,
                                   by = list(Cell_type = tree_list[[t]]$Node_type_counts$Cell_type),
                                   sum) # Count cell type numbers in tree
    tree_type_include$Include <- tree_type_include$x >= inclusion_limit # Set included celltypes to T in dataframe
    tree_type_include$Include[tree_type_include$Cell_type %in% force_include] <- T # Allow for forced inclusion of cell types
    inc_exc[grepl(tree_name, rownames(inc_exc)), tree_type_include$Cell_type] <- # Nodes in tree, cell types in tree (cells not in tree remain F)
      matrix(rep(tree_type_include$Include, each = sum(grepl(tree_name, rownames(inc_exc)))), 
             nrow = sum(grepl(tree_name, rownames(inc_exc))), 
             ncol = nrow(tree_type_include)) # Create matrix of (nodes in tree) x (cell types in tree), with values T for cell types that
    # have enough abundance.
  }
  # For asymmetric correlations, we also need to know which nodes (cumulatively!) contain which cell types.
  # Type inclusion type_inc: use to include nodes for cell types in case of asymmetric correlation.
  # Structure as norm_freq, entries F mean a node does not have any cells of a particular type.
  type_inc <- comparison_list$Normalized_frequencies > 0 # This information is already in comparison_list$Normalized_frequencies
  
  # Calculate the correlations. For asymmetric correlations, the rows designate the target, meaning the
  # nodes to correlate over are selected based on whether the cell type in the row is present in the node.
  cor_output <- matrix(NA, nrow = ncol(norm_freq), ncol = ncol(norm_freq))
  rownames(cor_output) <- colnames(norm_freq)
  colnames(cor_output) <- colnames(norm_freq)
  for(i in 1:ncol(norm_freq)){
    for(j in 1:ncol(norm_freq)){
      if(correlation_type == "Symmetric"){
        nodes_include <- inc_exc[, i] & inc_exc[, j]
      }else if(correlation_type == "Asymmetric"){
        nodes_include <- inc_exc[, i] & inc_exc[, j] & type_inc[, i]
      }
      
      cor_output[i, j] <- 
        wtd.cors(norm_freq[nodes_include, i], norm_freq[nodes_include, j], 
                 weight = node_weights$Size[nodes_include])#/sum(node_weights$Size[nodes_include]))
    }
  }
  
  return(cor_output)
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

ExtractTreeFrequencies <- function(tree_list, sample_fraction = 1, sampling_seed = 42){
  # Calculate cell type frequencies per node (comparison_list) 
  comparison_list <- list(Comparison = data.frame(Tree = character(),
                                                  Precursor = character(),
                                                  Type_count = integer()),
                          Node_sizes = data.frame(Size = integer()),
                          Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = nrow(celltype_frequencies))),
                          Frequencies = data.frame(matrix(nrow = 0, ncol = nrow(celltype_frequencies))))
  colnames(comparison_list$Frequencies) <- celltype_frequencies$Cell_type
  colnames(comparison_list$Normalized_frequencies) <- celltype_frequencies$Cell_type
  
  # Extract cell type frequencies
  set.seed(sampling_seed)
  for(t in 1:length(tree_list)){
    edge_list_full <- tree_list[[t]]$Edge_list
    # Cell.type == "NA" are nodes; leave those in, sample over the others.
    cell_list <- edge_list_full[edge_list_full$Cell.type != "NA", ]
    cell_list_s <- cell_list[sample.int(nrow(cell_list), round(size = sample_fraction * nrow(cell_list))), ]
    edge_list <- rbind(edge_list_full[edge_list_full$Cell.type == "NA", ], cell_list_s)
    
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
  
  return(comparison_list)
}

InitializeTreesFromDisk <- function(tree_path, libraries, cell_annotations){
  # Load all trees found in libraries from tree_path, filter cells with annotations
  # and append cell types, calculate edge lists and node counts, put all in tree list.
  tree_list_in <- list()
  
  for(this_sample in unique(libraries$Sample)){
    # Load tree, filter cells and append cell types.
    sample_libraries <- libraries$Library_name[libraries$Sample == this_sample]
    tree.in <- 
      ReadTree(this_sample, 
               reference_set = cell_annotations[cell_annotations$orig.ident %in% sample_libraries, ],
               tree_path = tree_path)
    tree.in$Tree$Do(RemapCellTypes,
                    reference_set = cell_annotations)
    
    # Calculate edge list
    edge_list <- ToDataFrameNetwork(tree.in$Tree, "Cell.type")
    
    # Calculate cell type counts per node
    node_type_counts <- data.frame(table(edge_list$from, edge_list$Cell.type))
    colnames(node_type_counts) <- c("Node", "Cell_type", "Type_count")
    node_type_counts$Node <- as.character(node_type_counts$Node)
    node_type_counts$Cell_type <- as.character(node_type_counts$Cell_type)
    node_type_counts <- node_type_counts[node_type_counts$Cell_type != "NA", ]
    
    # Put all tree data into the tree list.
    tree.in <- list(metadata = list(dpi = libraries$Dpi[libraries$Sample == this_sample],
                                    Name = this_sample),
                    Tree = tree.in$Tree,
                    Edge_list = edge_list,
                    Node_type_counts = node_type_counts)
    tree_list_in[[length(tree_list_in) + 1]] <- tree.in
    names(tree_list_in)[[length(tree_list_in)]] <- this_sample
  }
  
  return(tree_list_in)
}

Progenyl <- function(precursor, nodes){
  # Find all progeny for a precursor node.
  grepl(precursor, nodes) & precursor != nodes
} 

RemapCellTypes <- function(node, reference_set){
  if(!("Cell_type" %in% names(reference_set))){
    stop("Expecting column named Cell_type in reference_set")
  }
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

ReadTree <- function(library_name, reference_set, tree_path){
  # load(paste(tree_path, library_name, "_tree_pie.Robj", sep = ""))
  # load(paste("./Data/Trees/Old_trees/", library_name, "_tree_pie.Robj", sep = ""))
  list_out <- list(Tree = loadRData(paste(tree_path, library_name, "_tree_pie.Robj", sep = "")))
  list_out$Tree$Do(RemapCellTypes, 
                   reference_set = reference_set)
  
  return(list_out)
}

MakePieTree <- function(tree_sample, pie_tree_name, types = NULL, ctypes, ct_colors){
  tree_to_plot <- Clone(tree_sample$Tree)
  if(!is.null(types)){
    tree_to_plot$Do(ZoomCellTypes, 
                    zoom_types = types)
  }
  tree_sample$CTree <-    
    collapsibleTree(df = tree_to_plot, root = tree_to_plot$scar, pieNode = T,
                    pieSummary = T, collapsed = F,
                    width = 1000, height = 500,
                    ctypes = ctypes,linkLength=60, 
                    ct_colors = ct_colors, angle = pi/2,fontSize = 0,
                    nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6)) 
  
  names(tree_sample)[which(names(tree_sample) == "CTree")] <- pie_tree_name
  
  return(tree_sample)
}
