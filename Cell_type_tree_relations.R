# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)
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

# Dependencies? ####
# tree <- tree.H5
# count_scars <- function(tree){
#   scar <- data.frame(Scar = tree$Get('scar', filterFun = function(x) {!isLeaf(x)}),
#                      stringsAsFactors = F)
#   scar$Node <- rownames(scar)
#   scar$Scar_count <-
#     sapply(scar$Scar,
#            function(x) length(unlist(strsplit(x, ","))))
#   scar$Agg_scar_count <- NA
#   for(nl in 1:nrow(scar)){
#     n <- scar$Node[nl]
#     # ct <- node_counts$Cell_type[nl]
#     scar$Agg_scar_count[nl] <-
#       sum(scar$Scar_count[grepl(paste(n, "_|", n, "$", sep = ""), scar$Node)])
#   }
#   
#   return(scar)
# }
# determine_first_split <- function(enrichments){
#   nodes <- unique(enrichments[, c("Parent_node", "Node")])
#   nodes <- nodes[complete.cases(nodes), ]
#   nodes$Split <- 
#     sapply(nodes$Parent_node,
#            function(x){
#              length(nodes$Node[nodes$Parent_node == x]) > 1
#            }
#     )
#   nodes$Depth <-
#     sapply(nodes$Parent_node,
#            function(x){
#              length(unlist(strsplit(x, "_"))) - 1
#            }
#     )
#   first_split <- nodes[nodes$Split, ]
#   first_split <- first_split[first_split$Depth == min(first_split$Depth), 
#                              c("Parent_node", "Node")]
#   
#   return(first_split)
# }
# final_clone_cond_prob <- function(tree, cell_cutoff){
#   # Calculate conditional probabilities between cell types in lowest clones.
#   # First find lowest clones (lowest leaves with >= cell_cutoff cells),
#   # then calculate conditional probabilities of all combinations of cell types
#   # over these clones.
#   # Return conditional probabilities in matrix form
#   # (chance of seeing row, given column).
#   
#   # Find lowest clones
#   ## Determine cumulative counts for cell types
#   node_counts <- count_cumulative(tree)
#   colnames(node_counts)[4] <- "Type_ccount"
#   ## Remove counts associated to daughter nodes
#   node_counts <- node_counts[node_counts$Cell_type != "", ]
#   ## Calculate cumulative node counts
#   node_totals <- aggregate(node_counts$Type_ccount,
#                            by = list(Node = node_counts$Node),
#                            sum)
#   colnames(node_totals)[2] <- "Total_ccount"
#   ## Select nodes over cutoff
#   node_totals <- node_totals[node_totals$Total_ccount >= cell_cutoff, ]
#   ## Determine which nodes are lowest in a branch
#   node_totals$Lowest <-
#     sapply(node_totals$Node,
#            function(x)
#              sum(grepl(paste(x, "_", sep = ""), node_totals$Node)) == 0)
#   ## Select only lowest branch nodes
#   lowest_node_counts <- 
#     merge(node_counts, 
#           node_totals[node_totals$Lowest, c("Node", "Total_ccount")])
#   
#   # Calculate conditional probabilities P(A|B) = P(A&B)/P(B)
#   ## Make cell_type x clone matrix
#   cell_type_presence_matrix <-
#     acast(lowest_node_counts, Cell_type ~ Node, value.var = "Type_ccount")
#   ## Count only presence/absence
#   cell_type_presence_matrix[cell_type_presence_matrix > 0] <- 1
#   ## Coincidence rate matrix P(A&B) is presence matrix squared, divided over the
#   ## number of clones.
#   cell_type_coincidence <- 
#     cell_type_presence_matrix %*% t(cell_type_presence_matrix)/
#     length(unique(lowest_node_counts$Node))
#   ## Diagonal of this matrix is P(B), make diagonal-only matrix with 1/P(B)
#   cell_type_inv_chances <- diag(1/diag(cell_type_coincidence))
#   colnames(cell_type_inv_chances) <- colnames(cell_type_coincidence)
#   rownames(cell_type_inv_chances) <- rownames(cell_type_coincidence)
#   ## If we don't see a cell type, set 1/P -> 1 to
#   ## cell_type_inv_chances[is.infinite(cell_type_inv_chances)] <- 1
#   ## Calculate conditional probabilities P(A|B)=P(A&B)/P(B):
#   cell_type_cond_probs <- cell_type_coincidence %*% cell_type_inv_chances 
#   ## Note that the format here is: chance of seeing cell type in row given 
#   ## cell type in column.
#   
#   return(cell_type_cond_probs)
# }
# lik_binom_line <- function(alpha, t, kx, ky, nx, ny){
#   # alpha <- at[1]
#   # t <- at[-1]
#   calc_df <- 
#     data.frame(t = t, kx = kx, nx = nx, ky = ky, ny = ny)
#   calc_df$nx[calc_df$nx < 0] <- 0
#   calc_df$ny[calc_df$ny < 0] <- 0
#   calc_df$p1 <- calc_df$t * cos(alpha)
#   calc_df$p2 <- calc_df$t * sin(alpha)
#   calc_df$LL1 <- 
#     apply(calc_df[, c("kx", "nx", "p1")], 1,
#           function(x){
#             k <- x[1]
#             n <- x[2]
#             p <- x[3]
#             lik <- choose(n, k) * p^k * (1 - p)^(n-k)
#             if(lik <= 0){
#               LL <- LL <- 10^6
#             }else{
#               LL <- -log10(lik)
#             }
#             # if(LL == -Inf){
#             #   LL <- 10^6
#             # }
#             return(LL)
#           }
#     )
#   calc_df$LL2 <- 
#     apply(calc_df[, c("ky", "ny", "p2")], 1,
#           function(x){
#             k <- x[1]
#             n <- x[2]
#             p <- x[3]
#             lik <- choose(n, k) * p^k * (1 - p)^(n-k)
#             if(lik <= 0){
#               LL <- LL <- 10^6
#             }else{
#               LL <- -log10(lik)
#             }
#             return(LL)
#           }
#     )
#   calc_df$LL12 <- calc_df$LL1 + calc_df$LL2
#   return(sum(calc_df$LL12))
# }
# grad_lik_binom_line <- function(alpha, t, kx, ky, nx, ny){
#   # alpha <- at[1]
#   # t <- at[-1]
#   calc_df <- 
#     data.frame(t = t, kx = kx, nx = nx, ky = ky, ny = ny)
#   calc_df$nx[calc_df$nx < 0] <- 0
#   calc_df$ny[calc_df$ny < 0] <- 0
#   calc_df$df_dt <-
#     (calc_df$nx * calc_df$t * cos(alpha) - calc_df$kx)/(calc_df$t * (1 - calc_df$t * cos(alpha)) * log(10)) +
#     (calc_df$ny * calc_df$t * sin(alpha) - calc_df$ky)/(calc_df$t * (1 - calc_df$t * sin(alpha)) * log(10))
#   # calc_df$da_dt_xcomp <-
#   #   - calc_df$t * sin(alpha) * (calc_df$nx * calc_df$t * cos(alpha) - calc_df$kx)/(calc_df$t * cos(alpha) * (1 - calc_df$t * cos(alpha)) * log(10))
#   # calc_df$da_dt_ycomp <-
#   #   calc_df$t * cos(alpha) * (calc_df$ny * calc_df$t * sin(alpha) - calc_df$ky)/(calc_df$t * sin(alpha) * (1 - calc_df$t * sin(alpha)) * log(10))
#   # da_dt <- sum(sum(calc_df$da_dt_xcomp) + sum(calc_df$da_dt_ycomp))
#   
#   # return(c(da_dt, calc_df$df_dt))
#   return(calc_df$df_dt)
# }
# NodeEnrichment <- function(tree, p_cutoff){
#   # Extract node names and composition (Node_name x Cell-type, entries are counts)
#   node_counts <- count_cumulative(tree)
#   colnames(node_counts)[4] <- "Type_ccount"
#   # For each node, which cell types are enriched with respect to parent nodes?
#   # Consider: do we need to correct for average cell type depth?
#   # Node_name x Cell_type, entries are p-values for enrichment.
#   cumulative_node_counts <- aggregate(node_counts$Type_ccount,
#                                       by = list(Node = node_counts$Node),
#                                       sum)
#   colnames(cumulative_node_counts)[2] <- "Node_total"
#   node_counts <- merge(node_counts, cumulative_node_counts)
#   parent_node_probability <- node_counts[, c("Node", "Cell_type", "Type_ccount", "Node_total")]
#   colnames(parent_node_probability)[1] <- "Parent_node"
#   parent_node_probability$Parent_ratio <- 
#     parent_node_probability$Type_ccount/parent_node_probability$Node_total
#   node_counts <- merge(node_counts, 
#                        parent_node_probability[, c("Parent_node", "Cell_type", "Parent_ratio")], 
#                        all.x = T)
#   node_counts_2 <- node_counts[complete.cases(node_counts), ]
#   node_counts_2$Binom_p <-
#     apply(node_counts_2[, c("Type_ccount", "Node_total", "Parent_ratio")], 1,
#           function(x){
#             succ <- as.numeric(x[1])
#             n <- as.numeric(x[2])
#             p <- as.numeric(x[3])
#             bt <- binom.test(succ, n, p, alternative = "two.sided")
#             return(bt$p.value)
#           })
#   node_counts_2$Bp_adj <- p.adjust(node_counts_2$Binom_p, method = "fdr")
#   node_counts_2$SD <- sqrt(node_counts_2$Node_total * node_counts_2$Parent_ratio * (1 - node_counts_2$Parent_ratio))
#   node_counts_2$Z <- (node_counts_2$Type_ccount - node_counts_2$Parent_ratio * node_counts_2$Node_total)/
#     node_counts_2$SD
#   node_counts_2$Z[is.na(node_counts_2$Z)] <- 0
#   
#   enriched_types <- unique((node_counts_2[node_counts_2$Bp_adj < p_cutoff,])$Cell_type)
#   enriched_nodes <- unique((node_counts_2[node_counts_2$Bp_adj < p_cutoff,])$Node)
#   enrichment_comp_full <- dcast(node_counts_2, Cell_type ~ Node, value.var = "Z")
#   rownames(enrichment_comp_full) <- enrichment_comp_full$Cell_type
#   enrichment_comp_full <- enrichment_comp_full[, -1]
#   enrichment_comp_full <- enrichment_comp_full[which(rownames(enrichment_comp_full) %in% enriched_types), 
#                                                which(colnames(enrichment_comp_full) %in% enriched_nodes)]
#   
#   output_list <- 
#     list(Enrichment = node_counts_2,
#          Significant_enrichment = node_counts_2[node_counts_2$Bp_adj < p_cutoff,],
#          Comparison = enrichment_comp_full)
#   
#   return(output_list)
# }
# GeneratePalette <- function(x, total_colors = 100){
#   palette_size <- 5 * total_colors
#   start_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(palette_size)
#   
#   range_min <- min(range(x$Comparison))
#   range_max <- max(range(x$Comparison))
#   min_color_size <- abs(range_min) * total_colors/(abs(range_min) + abs(range_max))
#   max_color_size <- abs(range_max) * total_colors/(abs(range_min) + abs(range_max))
#   
#   neg_color_indx <- round(quantile(1:palette_size/2, probs = seq(0, 1, 1/min_color_size)))
#   pos_color_indx <- round(quantile((palette_size/2):palette_size, probs = seq(0, 1, 1/max_color_size)))
#   
#   ph_colors <- start_palette[c(neg_color_indx, pos_color_indx)]
#   
#   return(ph_colors)
# }

# Prepare colors and cell type data ####
type_colors <- read.csv("./Data/color.table.all.sub.csv", sep = ";", stringsAsFactors = F)[, -1]
# colors_interest[2] <- "#b2b220"

cell_types <- read.csv("./Data/all.hearts.all.cells.all.sub.sept03.csv", stringsAsFactors = F)
cell_types$Cell_type <- cell_types$immune.fibro.CM.subtypes
cell_types$Cell_name <- paste("nd", cell_types$X, sep = "")

ph_annotation <- data.frame(Celltype = type_colors$celltype, row.names = type_colors$celltype)
ph_zoom_annotation <- data.frame(Zoomtype = type_colors$celltype, row.names = type_colors$celltype)
ann_colors <-
  list(Celltype = setNames(type_colors$Color2, type_colors$celltype),
       Zoomtype = setNames(type_colors$colo1, type_colors$celltype))

# cell_types$Barcode <-sapply(cell_types$index,
                            # function(x) unlist(strsplit(x, "-"))[1])
# fibroblast_subtypes <- read.csv("./Scripts/Cells_types_H5_to_Hr7_fibroblast_subtypes.csv",
                                # stringsAsFactors = F)
# fibroblast_subtypes$heart <-
#   sapply(fibroblast_subtypes$Cell,
#          function(x) unlist(strsplit(x, "_"))[2])
# fibroblast_subtypes$Cell.type[fibroblast_subtypes$Cell.type == "Pericytes?"] <- "Pericytes"
# fibroblast_subtypes$Cell.type[fibroblast_subtypes$Cell.type == ""] <- NA
# cell_types_sub <- merge(cell_types[, c("Barcode", "batch", "heart", "Celltype")],
#                         fibroblast_subtypes[, c("Barcode", "heart", "Cell.type")], all = T)
# cell_types_sub <- cell_types_sub[!is.na(cell_types_sub$Celltype), ]
# cell_types_sub$Cell_type <- ifelse(is.na(cell_types_sub$Cell.type), cell_types_sub$Celltype, cell_types_sub$Cell.type)
# cell_types_sub$Cell_name <- paste("nd", cell_types_sub$heart, "_", cell_types_sub$Barcode, sep = "")

# Load trees ####
# Load individual trees, calculate enrichments and plot ####
# H5 - 0dpi
load("./Data/Trees/H5_tree_pie.Robj")
tree.H5 <- Clone(LINNAEUS.pie) # Otherwise this assignment works as a pointer
tree.H5$Do(RemapCellTypes, 
           reference_set = cell_types[cell_types$orig.ident == "H5", ])
# For every cell type and every node, calculate whether it is co-occurring with at
# least one cell of another cell type. Output: per tree, number of nodes where cell type
# occurs with any other cell type (cell type, number of nodes, co-occurring, number of nodes). 
# Later, calculate the chance it is co-occurring with 
# another cell type but we do not see it: the chance, that, given the distribution in the mother node and
# the number of cells in this node, we have simply missed the other cell type due to
# sampling.
node_counts_H5 <- count_cumulative(tree.H5)
# Calculate 0/1 matrix of node presence per cell type, multiple with transpose, diagonal
# is number of nodes for cell type, off-diag is number of co-occurrences.
CTN <- acast(node_counts_H5, Cell_type ~ Node, value.var = "Type_count")
CTN <- ifelse(CTN == 0, 0, 1)
Cooc_M <- CTN %*% t(CTN)
Occ <- data.frame(Cell_type = colnames(Cooc_H5),
                     Occurrence = diag(Cooc_H5), stringsAsFactors = F)
Cooc <- data.frame(t(combn(Occ$Cell_type, 2)))
colnames(Cooc) <- c("Cell_type", "Co_occurring_celltype")
Cooc$Co_occurrence_count <- Cooc_H5[lower.tri(Cooc_H5)]
Cooc <- merge(Occ, Cooc)
Cooc <- Cooc[, c("Cell_type", "Co_occurring_celltype", "Occurrence", "Co_occurrence_count")]
Cooc$Co_occurrence_freq <- Cooc$Co_occurrence_count/Cooc$Occurrence
Cooc_f_M <- Cooc_M/diag(Cooc_M)

pheatmap::pheatmap(Cooc_f_M)
# How to read this: in Cooc_f_M, the row is the cell type, the col is the co-occurring celltype.
# This is the same in the heatmap: the rows are the cell types, the cols are the celltypes that
# co-occur. So for instance, the lineage of the atrial cardiomyocytes could be everything, and
# that of the epicardium could as well. However, ventricular endocardium does not often co-occur
# with smooth muscle cells (rgs5a & elastic), suggesting these types split off early. How is this
# influenced by the number of times we see the cell types and co-occurring cell types?


load("./Data/Trees/Hr27_tree_pie.Robj")
tree.Hr27 <- Clone(LINNAEUS.pie) # Otherwise this assignment works as a pointer
tree.Hr27$Do(RemapCellTypes, 
           reference_set = cell_types[cell_types$orig.ident == "Hr27", ])
node_counts_Hr27 <- count_cumulative(tree.Hr27)

CTN <- acast(node_counts_Hr27, Cell_type ~ Node, value.var = "Type_count")
CTN <- ifelse(CTN == 0, 0, 1)
Cooc_M <- CTN %*% t(CTN)
Occ <- data.frame(Cell_type = colnames(Cooc_M),
                  Occurrence = diag(Cooc_M), stringsAsFactors = F)
Cooc <- data.frame(t(combn(Occ$Cell_type, 2)))
colnames(Cooc) <- c("Cell_type", "Co_occurring_celltype")
Cooc$Co_occurrence_count <- Cooc_M[lower.tri(Cooc_M)]
Cooc <- merge(Occ, Cooc)
Cooc <- Cooc[, c("Cell_type", "Co_occurring_celltype", "Occurrence", "Co_occurrence_count")]
Cooc$Co_occurrence_freq <- Cooc$Co_occurrence_count/Cooc$Occurrence
Cooc_f_M <- Cooc_M/diag(Cooc_M)

ph_annotation <- data.frame(Celltype = type_colors$celltype, row.names = type_colors$celltype)
ph_zoom_annotation <- data.frame(Zoomtype = type_colors$celltype, row.names = type_colors$celltype)
ann_colors <-
  list(Celltype = setNames(type_colors$Color2, type_colors$celltype),
       Zoomtype = setNames(type_colors$colo1, type_colors$celltype))
# zoom_colors <-
  # list(Celltype = setNames(type_colors$colo1, type_colors$celltype)

png("./Images/Hr27_celltype_cooccurrence.png")
pheatmap::pheatmap(Cooc_f_M, treeheight_row = 0, treeheight_col = 0, 
                   fontsize_row = 8, fontsize_col = 8, 
                   annotation_col = ph_annotation, annotation_row = ph_annotation,
                   annotation_colors = ann_colors, annotation_legend = F)
dev.off()
# Plot zoomed as well
png("./Images/Hr27_celltype_cooccurrence_fibrozoom.png")
pheatmap::pheatmap(Cooc_f_M[rownames(Cooc_f_M) %in% type_colors$celltype[type_colors$colo1 != ""], 
                            colnames(Cooc_f_M) %in% type_colors$celltype[type_colors$colo1 != ""]], 
                   treeheight_row = 0, treeheight_col = 0, 
                   fontsize_row = 12, fontsize_col = 12, 
                   annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
                   annotation_colors = ann_colors, annotation_legend = F)
dev.off()

# Plot tree ####
# Full tree
LINNAEUS.pie.wg <-
  collapsibleTree(df = tree.Hr27, root = tree.Hr27$scar, pieNode = T,
                  pieSummary = T, collapsed = F,
                  width = 600, height = 500,
                  ctypes = type_colors$celltype,linkLength=60, 
                  ct_colors = type_colors$Color2, angle = pi/2,
                  nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6))
LINNAEUS.pie.wg
htmlwidgets::saveWidget(
  LINNAEUS.pie.wg,
  file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie.html")

# Zoom
# zoom.edges <- make.edgelist(tree.summary.collapse, node.count.cumulative.agg, 
#                             correct.cell.placement,
#                             main.min = 5, off.main.min = 50)
# zoom.edges <- zoom.edges[zoom.edges$Cell.type != "NA" | zoom.edges$Child %in% zoom.siblings, ]
# zoom.edges$Scar.acquisition <- ""

tree.Hr27.fibrozoom <- Clone(tree.Hr27)
tree.Hr27.fibrozoom$Do(ZoomCellTypes, 
                       zoom_types = type_colors$celltype[type_colors$colo1 != ""])


# LINNAEUS.zoom <- generate_tree(zoom.edges)
LINNAEUS.pie.zoom.wg <-
  collapsibleTree(df = tree.Hr27.fibrozoom, root = tree.Hr27.fibrozoom$scar, pieNode = T,
                  pieSummary = T,collapsed = F,
                  width = 600, height = 500, linkLength = 60,
                  ctypes = type_colors$celltype, angle = pi/2,
                  ct_colors = type_colors$colo1,
                  nodeSize_class = c(20, 30, 50), nodeSize_breaks = c(0, 50, 1000, 1e6))
LINNAEUS.pie.zoom.wg
htmlwidgets::saveWidget(
  LINNAEUS.pie.zoom.wg,
  file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie_fibrozoom.html")

# Potential progenitor graph ####
zoom_types <- setdiff(type_colors$celltype[type_colors$colo1 != ""],
                      c("Ery. Duplex", "M. Duplex", "Myofibroblasts", "Fibroblast (nppc)"))
prog_potential <- melt(Cooc_f_M[rownames(Cooc_f_M) %in% zoom_types, 
                                        colnames(Cooc_f_M) %in% zoom_types])
colnames(prog_potential) <- c("Child", "Progenitor", "Cooc_freq")
prog_potential$Pot_prog <- (prog_potential$Cooc_freq == 1)
prog_potential_graph <- 
  graph_from_edgelist(as.matrix(prog_potential[prog_potential$Pot_prog, c("Progenitor", "Child")]))
vertex_colors <- data.frame(celltype = names(V(prog_potential_graph)),
                            Order = 1:(length(V(prog_potential_graph))))
vertex_colors <- merge(vertex_colors, type_colors[, c("celltype", "colo1")])
vertex_colors <- vertex_colors[order(vertex_colors$Order), ]
vertex_colors$Label <-
  c("Epi A", "Epi V", "cfd", "col11", "col12", "cxcl12", "prolif", "spock3", "F")
plot(prog_potential_graph, vertex.color = vertex_colors$colo1)
png("./Images/Hr27_progenitor_potential_fibrozoom.png", width = 640, height = 640)
plot(simplify(prog_potential_graph), vertex.color = vertex_colors$colo1,
     vertex.size = 35, vertex.label = vertex_colors$Label, vertex.label.cex = 2,
     edge.color = "black", edge.width = 2)
dev.off()

# Make list of trees for the different datasets that includes per dataset: tree, zoom,
# Cooc_f_M, prog_potential + pp-graph.

# Read in tree object and append cell types
# Create tree visualization and zoom visualization
# Calculate co-occurrences
# Optionally plot heatmap of co-occurrences
# Create potential progenitors and pp-graph