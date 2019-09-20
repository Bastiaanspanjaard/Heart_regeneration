# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)
CalculateCooccurrence <- function(tree_sample){
  node_counts <- count_cumulative(tree_sample$Tree)
  
  CTN <- acast(node_counts, Cell_type ~ Node, value.var = "Type_count")
  CTN <- ifelse(CTN == 0, 0, 1)
  Cooc_M <- CTN %*% t(CTN)
  # Occ <- data.frame(Cell_type = colnames(Cooc_M),
  #                   Occurrence = diag(Cooc_M), stringsAsFactors = F)
  # Cooc <- data.frame(t(combn(Occ$Cell_type, 2)))
  # colnames(Cooc) <- c("Cell_type", "Co_occurring_celltype")
  # Cooc$Co_occurrence_count <- Cooc_M[lower.tri(Cooc_M)]
  # Cooc <- merge(Occ, Cooc)
  # Cooc <- Cooc[, c("Cell_type", "Co_occurring_celltype", "Occurrence", "Co_occurrence_count")]
  # Cooc$Co_occurrence_freq <- Cooc$Co_occurrence_count/Cooc$Occurrence
  Cooc_f_M <- Cooc_M/diag(Cooc_M)
  
  tree_sample$Relative_cooccurrence <- Cooc_f_M
  return(tree_sample)
}

CalculateProgenitors <- function(tree_sample, zoom_types){
  Cooc_f_M <- tree_sample$Relative_cooccurrence
  prog_potential <- melt(Cooc_f_M[rownames(Cooc_f_M) %in% zoom_types, 
                                  colnames(Cooc_f_M) %in% zoom_types])
  colnames(prog_potential) <- c("Child", "Progenitor", "Cooc_freq")
  prog_potential$Pot_prog <- (prog_potential$Cooc_freq == 1)
  prog_potential_graph <- 
    simplify(graph_from_edgelist(as.matrix(prog_potential[prog_potential$Pot_prog, c("Progenitor", "Child")])))
  
  plot(prog_potential_graph)
  
  tree_sample$Progenitor_potential <- prog_potential
  tree_sample$Progenitor_graph <- prog_potential_graph
  
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
  node_cocc$Ccount <- NA
  for(nl in 1:nrow(node_cocc)){
    n <- node_counts$Node[nl]
    ct <- node_counts$Cell_type[nl]
    node_cocc$Ccount[nl] <-
      sum(node_counts$Type_count[grepl(paste(n, "_|", n, "$", sep = ""), node_counts$Node) &
                                   node_counts$Cell_type == ct])
  }
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
  
  # node_cocc <- node_occ[, c("Var1", "Var2")]
  # node_cocc$CFreq <- NA
  # nl <- 1
  # for(nl in 1:nrow(node_cocc)){
  #   n <- node_occ$Var1[nl]
  #   ct <- node_occ$Var2[nl]
  #   node_cocc$CFreq[nl] <- 
  #     sum(node_occ$Freq[grepl(n, node_occ$Var1) & node_occ$Var2 == ct])
  # }
  # 
  # node_cocc_df <- dcast(node_cocc, Var1 ~ Var2, value.var = "CFreq")
  # 
  # 
  # 
  # node_totals <- aggregate(node_counts$Type_count,
  #                          by = list(Node = node_counts$Node),
  #                          sum)
  # colnames(node_totals)[2] <- "Node_total"
  # node_counts <- merge(node_counts, node_totals)
  # 
  # node_counts$Type_frequency <-
  #   node_counts$Type_count/node_counts$Node_total
  
  # return(node_counts)
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
                      c("Ery. Duplex", "M. Duplex", "Myofibroblasts", "Fibroblast (nppc)"))
vertex_colors <- data.frame(celltype = zoom_types,
                            Order = 1:(length(zoom_types)))
vertex_colors <- merge(vertex_colors, type_colors[, c("celltype", "colo1")])
vertex_colors <- vertex_colors[order(vertex_colors$Order), ]
vertex_colors$Label <-
  c("Epi A", "Epi V", "cfd", "col11", "col12", "cxcl12", "prolif", "spock3", "F")
rownames(vertex_colors) <- vertex_colors$celltype

# New code ####
# Read in tree object and append cell types
Hr10 <- ReadTree("Hr10", reference_set = cell_types[cell_types$orig.ident == "Hr10", ])
Hr11 <- ReadTree("Hr11", reference_set = cell_types[cell_types$orig.ident == "Hr11", ])
Hr12 <- ReadTree("Hr12", reference_set = cell_types[cell_types$orig.ident == "Hr12", ])
Hr24 <- ReadTree("Hr24", reference_set = cell_types[cell_types$orig.ident == "Hr24", ])
Hr26 <- ReadTree("Hr26", reference_set = cell_types[cell_types$orig.ident == "Hr26", ])
Hr27 <- ReadTree("Hr27", reference_set = cell_types[cell_types$orig.ident == "Hr27", ])
# Create tree visualization and zoom visualization
Hr10 <- MakePieTree(Hr10, "Full_tree", ct_colors = type_colors$Color2)
Hr11 <- MakePieTree(Hr11, "Full_tree", ct_colors = type_colors$Color2)
Hr12 <- MakePieTree(Hr12, "Full_tree", ct_colors = type_colors$Color2)
Hr24 <- MakePieTree(Hr24, "Full_tree", ct_colors = type_colors$Color2)
Hr26 <- MakePieTree(Hr26, "Full_tree", ct_colors = type_colors$Color2)
Hr27 <- MakePieTree(Hr27, "Full_tree", ct_colors = type_colors$Color2)
Hr10$Full_tree
Hr11$Full_tree
Hr12$Full_tree
Hr24$Full_tree
Hr26$Full_tree
Hr27$Full_tree
htmlwidgets::saveWidget(
  Hr27$Full_tree,
  file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie.html")
Hr10 <- MakePieTree(Hr10, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr11 <- MakePieTree(Hr11, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr12 <- MakePieTree(Hr12, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr24 <- MakePieTree(Hr24, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr26 <- MakePieTree(Hr26, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr27 <- MakePieTree(Hr27, "Fibrozoom_tree", types = zoom_types, ct_colors = type_colors$colo1)
Hr10$Fibrozoom_tree
Hr11$Fibrozoom_tree
Hr12$Fibrozoom_tree
Hr24$Fibrozoom_tree
Hr26$Fibrozoom_tree
Hr27$Fibrozoom_tree
htmlwidgets::saveWidget(
  Hr27$Fibrozoom_tree,
  file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie_fibrozoom.html")

# Calculate co-occurrences
Hr10 <- CalculateCooccurrence(Hr10)
Hr11 <- CalculateCooccurrence(Hr11)
Hr12 <- CalculateCooccurrence(Hr12)
Hr24 <- CalculateCooccurrence(Hr24)
Hr26 <- CalculateCooccurrence(Hr26)
Hr27 <- CalculateCooccurrence(Hr27)
# png("./Images/Hr27_celltype_cooccurrence.png")
# pheatmap(Hr27$Relative_cooccurrence, treeheight_row = 0, treeheight_col = 0, 
#          fontsize_row = 8, fontsize_col = 8, 
#          annotation_col = ph_annotation, annotation_row = ph_annotation,
#          annotation_colors = ann_colors, annotation_legend = F)
# dev.off()
# png("./Images/Hr27_celltype_cooccurrence_fibrozoom.png")
pheatmap(Hr10$Relative_cooccurrence[rownames(Hr10$Relative_cooccurrence) %in% zoom_types, 
                  colnames(Hr10$Relative_cooccurrence) %in% zoom_types], 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)
# dev.off()
View(Hr10$Relative_cooccurrence[rownames(Hr10$Relative_cooccurrence) %in% zoom_types, 
                                colnames(Hr10$Relative_cooccurrence) %in% zoom_types])
# Create potential progenitors and pp-graph
Hr10 <- CalculateProgenitors(Hr10, zoom_types)
Hr11 <- CalculateProgenitors(Hr11, zoom_types)
Hr12 <- CalculateProgenitors(Hr12, zoom_types)
Hr24 <- CalculateProgenitors(Hr24, zoom_types)
Hr26 <- CalculateProgenitors(Hr26, zoom_types)
Hr27 <- CalculateProgenitors(Hr27, zoom_types)

plot(Hr10$Progenitor_graph)

png("./Images/Hr27_progenitor_potential_fibrozoom.png", width = 640, height = 640)
vertex_colors_graph <- vertex_colors[names(V(Hr27$Progenitor_graph)), ]
plot(Hr27$Progenitor_graph, vertex.color = vertex_colors_graph$colo1,
     vertex.size = 35, vertex.label = vertex_colors_graph$Label, vertex.label.cex = 2,
     edge.color = "black", edge.width = 2)
dev.off()

