# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)

# tree_sample <- Hr27
CalculateCooccurrence <- function(tree_sample){
  node_counts <- count_cumulative(tree_sample$Tree)
  
  # For each cell type and each node, if there's 0 of another cell type, calculate the chance of this given
  # the cell type's abundance in the parent node.
  # Nodes that our descendants appear in
  type_counts <- 
    node_counts[node_counts$Type_count != 0, 
                c("Cell_type", "Node", "Parent_node")]
  # Potential precursor node counts
  parent_counts <- node_counts[, c("Node", "Cell_type", "Ccount", "Parent_node")]
  colnames(parent_counts) <-  c("Node", "Precursor", "Precursor_count", "Parent_node")
  # Full node counts
  full_node_counts <- aggregate(node_counts$Ccount,
                                by = list(Node = node_counts$Node),
                                sum)
  colnames(full_node_counts)[2] <- "Total"
  # Cell type proportions
  proportions <- merge(node_counts[, c("Cell_type", "Node", "Ccount")], full_node_counts)
  proportions$Frequency <- proportions$Ccount/proportions$Total
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
  
  tree_sample$Relative_cooccurrence <- Cooc_f_M
  tree_sample$Descendancy <- descendancy
  tree_sample$Aggregated_descendancy <- descendancy_agg
  return(tree_sample)
}

CalculateProgenitors <- function(tree_sample, zoom_to, zoom_from){
  Cooc_f_M <- tree_sample$Relative_cooccurrence
  prog_potential <- melt(Cooc_f_M[rownames(Cooc_f_M) %in% zoom_to, 
                                  colnames(Cooc_f_M) %in% zoom_from])
  colnames(prog_potential) <- c("Child", "Progenitor", "Cooc_freq")
  prog_potential$Pot_prog <- (prog_potential$Cooc_freq == 1)
  prog_potential_graph <- 
    simplify(graph_from_edgelist(as.matrix(prog_potential[prog_potential$Pot_prog, c("Progenitor", "Child")])))
  
  # plot(prog_potential_graph)
  
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

NodeEnrichment <- function(tree_sample, p_cutoff){
  tree <- tree_sample$Tree
  # Extract node names and composition (Node_name x Cell-type, entries are counts)
  node_counts <- count_cumulative(tree)
  colnames(node_counts)[4] <- "Type_ccount"
  # For each node, which cell types are enriched with respect to parent nodes?
  # Consider: do we need to correct for average cell type depth?
  # Node_name x Cell_type, entries are p-values for enrichment.
  cumulative_node_counts <- aggregate(node_counts$Type_ccount,
                                      by = list(Node = node_counts$Node),
                                      sum)
  colnames(cumulative_node_counts)[2] <- "Node_total"
  node_counts <- merge(node_counts, cumulative_node_counts)
  parent_node_probability <- node_counts[, c("Node", "Cell_type", "Type_ccount", "Node_total")]
  colnames(parent_node_probability)[1] <- "Parent_node"
  parent_node_probability$Parent_ratio <- 
    parent_node_probability$Type_ccount/parent_node_probability$Node_total
  node_counts <- merge(node_counts, 
                       parent_node_probability[, c("Parent_node", "Cell_type", "Parent_ratio")], 
                       all.x = T)
  node_counts_2 <- node_counts[complete.cases(node_counts), ]
  node_counts_2$Binom_p <-
    apply(node_counts_2[, c("Type_ccount", "Node_total", "Parent_ratio")], 1,
          function(x){
            succ <- as.numeric(x[1])
            n <- as.numeric(x[2])
            p <- as.numeric(x[3])
            bt <- binom.test(succ, n, p, alternative = "two.sided")
            return(bt$p.value)
          })
  node_counts_2$Bp_adj <- p.adjust(node_counts_2$Binom_p, method = "fdr")
  node_counts_2$SD <- sqrt(node_counts_2$Node_total * node_counts_2$Parent_ratio * (1 - node_counts_2$Parent_ratio))
  node_counts_2$Z <- (node_counts_2$Type_ccount - node_counts_2$Parent_ratio * node_counts_2$Node_total)/
    node_counts_2$SD
  node_counts_2$Z[is.na(node_counts_2$Z)] <- 0
  
  enriched_types <- unique((node_counts_2[node_counts_2$Bp_adj < p_cutoff,])$Cell_type)
  enriched_nodes <- unique((node_counts_2[node_counts_2$Bp_adj < p_cutoff,])$Node)
  enrichment_comp_full <- dcast(node_counts_2, Cell_type ~ Node, value.var = "Z")
  rownames(enrichment_comp_full) <- enrichment_comp_full$Cell_type
  enrichment_comp_full <- enrichment_comp_full[, -1]
  enrichment_comp_full <- enrichment_comp_full[which(rownames(enrichment_comp_full) %in% enriched_types), 
                                               which(colnames(enrichment_comp_full) %in% enriched_nodes)]
  
  sign_enr <- node_counts_2[node_counts_2$Bp_adj < p_cutoff,]
  enr_m <- acast(sign_enr[sign_enr$Z > 0, ], 
                      Cell_type ~ Node, value.var = "Z")
  enr_m[is.na(enr_m)] <- 0
  enr_m[enr_m > 0] <- 1
  co_enrich <- enr_m %*% t(enr_m)
  
  tree_sample$Enrichment <- node_counts_2
  tree_sample$Significant_enrichment <- node_counts_2[node_counts_2$Bp_adj < p_cutoff,]
  tree_sample$Co_enrichment <- co_enrich
  
  # output_list <- 
  #   list(Enrichment = node_counts_2,
  #        Significant_enrichment = node_counts_2[node_counts_2$Bp_adj < p_cutoff,],
  #        Comparison = enrichment_comp_full)
  
  # return(output_list)
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
zoom_to <- c("Fibroblast (col11a1a)", "Fibroblast (col12a1a)", "Fibroblast (proliferating)",
             "Fibroblasts (spock3)")
zoom_from <- setdiff(zoom_types, zoom_to)

# Load tree objects, append cell types, create lineage trees ####
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
# htmlwidgets::saveWidget(
#   Hr27$Full_tree,
#   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie.html")
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
# htmlwidgets::saveWidget(
#   Hr27$Fibrozoom_tree,
#   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie_fibrozoom.html")

# Analyse lineage potential ####
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
pheatmap(Hr27$Relative_cooccurrence[rownames(Hr27$Relative_cooccurrence) %in% zoom_to, 
                  colnames(Hr27$Relative_cooccurrence) %in% zoom_from], 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)
# dev.off()

agg_desc_Hr10 <- Hr10$Aggregated_descendancy[Hr10$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr10$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr10$Tree <- "Hr10"
agg_desc_Hr11 <- Hr11$Aggregated_descendancy[Hr11$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr11$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr11$Tree <- "Hr11"
agg_desc_Hr12 <- Hr12$Aggregated_descendancy[Hr12$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr12$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr12$Tree <- "Hr12"
agg_desc_Hr24 <- Hr24$Aggregated_descendancy[Hr24$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr24$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr24$Tree <- "Hr24"
agg_desc_Hr26 <- Hr26$Aggregated_descendancy[Hr26$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr26$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr26$Tree <- "Hr26"
agg_desc_Hr27 <- Hr27$Aggregated_descendancy[Hr27$Aggregated_descendancy$Cell_type %in% zoom_to &
                                               Hr27$Aggregated_descendancy$Precursor %in% zoom_from, ]
agg_desc_Hr27$Tree <- "Hr27"
agg_desc_trees <- rbind(agg_desc_Hr10, agg_desc_Hr11, agg_desc_Hr12,
                        agg_desc_Hr24, agg_desc_Hr26, agg_desc_Hr27)
agg_desc <- aggregate(agg_desc_trees$Tree_precursor_p,
                      by = list(Cell_type = agg_desc_trees$Cell_type,
                                Precursor = agg_desc_trees$Precursor),
                      prod)
colnames(agg_desc)[3] <- "p"
agg_d_cast <- acast(agg_desc, Cell_type ~ Precursor, value.var = "p")
# png("./Images/Potential_precursors_3dpi_p_product.png")
pheatmap(agg_d_cast, 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)
# dev.off()

# View(Hr10$Descendancy[Hr10$Descendancy$Precursor == "Fibroblasts" & Hr10$Descendancy$Cell_type == "Fibroblast (col11a1a)", ])

# Show dot plots for a descendant cell type, plotting the probabilities for all precursor cell types in all
# nodes
descendancy_Hr10 <- Hr10$Descendancy
descendancy_Hr10$Tree <- "Hr10"
descendancy_Hr11 <- Hr11$Descendancy
descendancy_Hr11$Tree <- "Hr11"
descendancy_Hr12 <- Hr12$Descendancy
descendancy_Hr12$Tree <- "Hr12"
descendancy_Hr24 <- Hr24$Descendancy
descendancy_Hr24$Tree <- "Hr24"
descendancy_Hr26 <- Hr26$Descendancy
descendancy_Hr26$Tree <- "Hr26"
descendancy_Hr27 <- Hr27$Descendancy
descendancy_Hr27$Tree <- "Hr27"
full_descendancy <- rbind(descendancy_Hr10, descendancy_Hr11, descendancy_Hr12,
                          descendancy_Hr24, descendancy_Hr26, descendancy_Hr27)
# png("./Images/Col11fib_3dpi_potential_precursors.png", width = 1366, height = 768)
ggplot(full_descendancy[full_descendancy$Cell_type == "Fibroblast (col11a1a)", ]) +
  geom_jitter(aes(x = 0, y = Precursor_presence_p, color = Precursor), size = 2) +
  labs(title = "Fibroblast (col11a1a) precursors",
       x = "", y = "Probability") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_color_manual(values = ann_colors$Celltype) +
  facet_wrap(~Precursor) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.5, 1))
# dev.off()
# This plot is very interesting because it shows something that is still
# unclear about this analysis. Many cell types look like they could be precursors to
# col11-fibroblasts, whereas such lineage relationship would not biologically make
# a lot of sense, e.g. erythrocytes. Erythrocytes look good
# because their detection rate is so low - note we are testing here what cannot be
# true, not what is maybe true. Erythrocytes are hard to call because of their low detection
# rate, so they're not showing up as a clear impossibility. This does not mean they are
# a good possibility.

# png("./Images/Col12fib_3dpi_potential_precursors_zoom.png", width = 683, height = 256)
ggplot(full_descendancy[full_descendancy$Cell_type == "Fibroblast (col12a1a)" &
                          full_descendancy$Precursor %in% zoom_from, ]) +
  geom_jitter(aes(x = 0, y = log10(Precursor_presence_p), color = Precursor)) +
  labs(title = "Fibroblast (col12a1a) precursors",
       y = "Probability (log10)", x = "") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_color_manual(values = ann_colors$Zoomtype) +
  facet_wrap(~Precursor, nrow = 1) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"))
# dev.off()  

# png("./Images/3dpi_potential_precursors_zoom.png", width = 683, height = 768)
ggplot(full_descendancy[full_descendancy$Cell_type %in% zoom_to &
                          full_descendancy$Precursor %in% zoom_from, ]) +
  geom_jitter(aes(x = 0, y = log10(Precursor_presence_p), color = Precursor), size = 2) +
  labs(title = "3dpi precursors",
       y = "Probability (log10)", x = "") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_color_manual(values = ann_colors$Zoomtype) +
  facet_grid(Cell_type~Precursor) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"))
# dev.off()

# Create potential progenitors and pp-graph (old code, not sure is useful anymore)
Hr10 <- CalculateProgenitors(Hr10, zoom_to, zoom_from)
Hr11 <- CalculateProgenitors(Hr11, zoom_to, zoom_from)
Hr12 <- CalculateProgenitors(Hr12, zoom_to, zoom_from)
Hr24 <- CalculateProgenitors(Hr24, zoom_to, zoom_from)
Hr26 <- CalculateProgenitors(Hr26, zoom_to, zoom_from)
Hr27 <- CalculateProgenitors(Hr27, zoom_to, zoom_from)

# plot(Hr10$Progenitor_graph)

# png("./Images/Hr26_progenitor_potential_fibrozoom.png", width = 640, height = 640)
vertex_colors_graph <- vertex_colors[names(V(Hr10$Progenitor_graph)), ]
plot(Hr26$Progenitor_graph, vertex.color = vertex_colors_graph$colo1,
     vertex.size = 35, vertex.label = vertex_colors_graph$Label, vertex.label.cex = 2,
     edge.color = "black", edge.width = 2)
# dev.off()
# Do plot for only nodes enriched in descendant cell type?

# Co-enrichment ####
enrichment_Hr10 <- NodeEnrichment(Hr10$Tree, p_cutoff = 0.01)
enr_Hr10_m <- acast(enrichment_Hr10$Significant_enrichment[enrichment_Hr10$Significant_enrichment$Z > 0, ], 
                   Cell_type ~ Node)
enr_Hr10_m[is.na(enr_Hr10_m)] <- 0
enr_Hr10_m[enr_Hr10_m > 0] <- 1
co_enrich <- enr_Hr10_m %*% t(enr_Hr10_m)

Hr10 <- NodeEnrichment(Hr10, p_cutoff = 0.01)
Hr11 <- NodeEnrichment(Hr11, p_cutoff = 0.01)
Hr12 <- NodeEnrichment(Hr12, p_cutoff = 0.01)
Hr24 <- NodeEnrichment(Hr24, p_cutoff = 0.01)
Hr26 <- NodeEnrichment(Hr26, p_cutoff = 0.01)
Hr27 <- NodeEnrichment(Hr27, p_cutoff = 0.01)
Hr10$Co_enrichment[rownames(Hr10$Co_enrichment) %in% zoom_to, ]
Hr10_coenrm <- melt(Hr10$Co_enrichment)
colnames(Hr10_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr10_coenrm$Tree <- "Hr10"
Hr10_coenrm <- Hr10_coenrm[Hr10_coenrm$Cell_type %in% zoom_to & 
                             Hr10_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]
Hr11_coenrm <- melt(Hr11$Co_enrichment)
colnames(Hr11_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr11_coenrm$Tree <- "Hr11"
Hr11_coenrm <- Hr11_coenrm[Hr11_coenrm$Cell_type %in% zoom_to & 
                             Hr11_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]
Hr12_coenrm <- melt(Hr12$Co_enrichment)
colnames(Hr12_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr12_coenrm$Tree <- "Hr12"
Hr12_coenrm <- Hr12_coenrm[Hr12_coenrm$Cell_type %in% zoom_to & 
                             Hr12_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]
Hr24_coenrm <- melt(Hr24$Co_enrichment)
colnames(Hr24_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr24_coenrm$Tree <- "Hr24"
Hr24_coenrm <- Hr24_coenrm[Hr24_coenrm$Cell_type %in% zoom_to & 
                             Hr24_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]
Hr26_coenrm <- melt(Hr26$Co_enrichment)
colnames(Hr26_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr26_coenrm$Tree <- "Hr26"
Hr26_coenrm <- Hr26_coenrm[Hr26_coenrm$Cell_type %in% zoom_to & 
                             Hr26_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]
Hr27_coenrm <- melt(Hr27$Co_enrichment)
colnames(Hr27_coenrm) <- c("Cell_type", "Precursor", "Co-enrichment")
Hr27_coenrm$Tree <- "Hr27"
Hr27_coenrm <- Hr27_coenrm[Hr27_coenrm$Cell_type %in% zoom_to & 
                             Hr27_coenrm$Precursor %in% setdiff(type_colors$celltype, zoom_to), ]

coenrichments_all <- rbind(Hr10_coenrm, Hr11_coenrm, Hr12_coenrm,
                           Hr24_coenrm, Hr26_coenrm, Hr27_coenrm)
coenrichments_agg <-
  aggregate(coenrichments_all$`Co-enrichment`,
            by = list(Cell_type = coenrichments_all$Cell_type,
                      Precursor = coenrichments_all$Precursor),
            sum)
