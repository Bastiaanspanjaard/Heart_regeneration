# Description ####
# Functions for analysing and constructing single-cell trees

# Written by B. Spanjaard, 2018

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
