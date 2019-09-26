# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)

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
  
  nodes <- node_counts[, c("Node", "Cell_type", "Ccount")]
  baseline <- nodes[nodes$Node == "nd0", c("Node", "Cell_type", "Ccount")]
  colnames(baseline)[3] <- "Total"
  nodes <- merge(nodes, baseline[, c("Cell_type", "Total")])
  nodes$Rel_freq <- nodes$Ccount/nodes$Total
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
tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr24 = Hr24, Hr26 = Hr26, Hr27 = Hr27)
# Create tree visualization and zoom visualization
for(t in 1:length(tree_list)){
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Full_tree", ct_colors = type_colors$Color2)
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
                                ct_colors = type_colors$colo1)
}
# htmlwidgets::saveWidget(
#   tree_list$Hr27$Full_tree,
#   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie.html")
# htmlwidgets::saveWidget(
#   tree_list$Hr27$Fibrozoom_tree,
#   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie_fibrozoom.html")

# Analyse lineage potential ####
agg_desc_trees <- data.frame(Cell_type = character(),
                             Precursor = character(),
                             Tree_precursor_p = numeric(),
                             Tree = character())
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
  descendancy_add <- tree_list[[t]]$Descendancy
  descendancy_add$Tree <- names(tree_list)[t]
  descendancy_add <- merge(descendancy_add, tree_list[[t]]$Node_entropy)
  full_descendancy <- rbind(full_descendancy, descendancy_add)
  
  agg_desc_add <- 
    tree_list[[t]]$Aggregated_descendancy[tree_list[[t]]$Aggregated_descendancy$Cell_type %in% zoom_to &
                                            tree_list[[t]]$Aggregated_descendancy$Precursor %in% zoom_from, ]
  agg_desc_add$Tree <- names(tree_list)[t]
  agg_desc_trees <- rbind(agg_desc_trees, agg_desc_add)
}
agg_desc <- aggregate(agg_desc_trees$Tree_precursor_p,
                      by = list(Cell_type = agg_desc_trees$Cell_type,
                                Precursor = agg_desc_trees$Precursor),
                      prod)
colnames(agg_desc)[3] <- "p"
agg_d_cast <- acast(agg_desc, Cell_type ~ Precursor, value.var = "p")

# png("./Images/Hr27_celltype_cooccurrence.png")
# pheatmap(tree_list$Hr27$Relative_cooccurrence, treeheight_row = 0, treeheight_col = 0,
#          fontsize_row = 8, fontsize_col = 8,
#          annotation_col = ph_annotation, annotation_row = ph_annotation,
#          annotation_colors = ann_colors, annotation_legend = F)
# dev.off()
# png("./Images/Hr27_celltype_cooccurrence_fibrozoom.png")
pheatmap(tree_list$Hr27$Relative_cooccurrence[rownames(tree_list$Hr27$Relative_cooccurrence) %in% zoom_to, 
                  colnames(tree_list$Hr27$Relative_cooccurrence) %in% zoom_from], 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)
# dev.off()

# png("./Images/Potential_precursors_3dpi_p_product.png")
pheatmap(agg_d_cast, 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)
# dev.off()

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

ggplot(full_descendancy[full_descendancy$Cell_type %in% zoom_to &
                          full_descendancy$Precursor %in% zoom_from, ]) +
  geom_jitter(aes(x = Entropy, y = log10(Precursor_presence_p), color = Precursor), size = 2) +
  labs(title = "3dpi precursors",
       y = "Probability (log10)", x = "Node entropy") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_color_manual(values = ann_colors$Zoomtype) +
  facet_grid(Cell_type~Precursor) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size = 12, face = "bold"))

ggplot(full_descendancy[full_descendancy$Cell_type == "Fibroblast (col11a1a)", ]) +
  geom_point(aes(x = Entropy, y = Precursor_presence_p, color = Precursor), size = 2) +
  labs(title = "Fibroblast (col11a1a) precursors",
       x = "Node entropy", y = "Probability") +
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none") +
  scale_color_manual(values = ann_colors$Celltype) +
  facet_wrap(~Precursor) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 12, face = "bold")) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

all_node_entropies <- unique(full_descendancy[, c("Node", "Tree", "Entropy")])

ggplot(all_node_entropies) +
  geom_histogram(aes(x = Entropy))

ggplot(full_descendancy[full_descendancy$Entropy > 0.5, ]) + 
  geom_histogram(aes(x = Precursor_presence_p), binwidth = 0.01)
full_descendancy$weighed_p <- full_descendancy$Precursor_presence_p^(1 - full_descendancy$Entropy)

ggplot(full_descendancy[full_descendancy$Cell_type %in% zoom_to &
                          full_descendancy$Precursor %in% zoom_from, ]) +
  geom_jitter(aes(x = 0, y = log10(weighed_p), color = Precursor), size = 2) +
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

agg_wd <- aggregate(full_descendancy$weighed_p,
                    by = list(Cell_type = full_descendancy$Cell_type,
                              Precursor = full_descendancy$Precursor),
                    prod)
colnames(agg_wd)[3] <- "wp"

agg_wd_cast <- acast(agg_wd[agg_wd$Cell_type %in% zoom_to &
                             agg_wd$Precursor %in% zoom_from, ], Cell_type ~ Precursor, value.var = "wp")
pheatmap(agg_wd_cast, 
         treeheight_row = 0, treeheight_col = 0, 
         fontsize_row = 12, fontsize_col = 12, 
         annotation_col = ph_zoom_annotation, annotation_row = ph_zoom_annotation,
         annotation_colors = ann_colors, annotation_legend = F)

# Node occupancy correlation ####
nodes <- count_cumulative(tree_list$Hr10$Tree)[, c("Node", "Cell_type", "Ccount")]
baseline <- nodes[nodes$Node == "nd0", c("Node", "Cell_type", "Ccount")]
colnames(baseline)[3] <- "Total"
nodes <- merge(nodes, baseline[, c("Cell_type", "Total")])
nodes$Rel_freq <- nodes$Ccount/nodes$Total
nodes_c <- acast(nodes[nodes$Node != "nd0", ], Node ~ Cell_type, value.var = "Rel_freq")
nodes_c <- nodes_c[, colSums(nodes_c) != 0]
cnc <- cor(nodes_c)
View(cnc[, colnames(cnc) %in% zoom_to])
pheatmap(cor(nodes_c))

nodes <- node_counts[, c("Node", "Cell_type", "Ccount")]
baseline <- nodes[nodes$Node == "nd0", c("Node", "Cell_type", "Ccount")]
colnames(baseline)[3] <- "Total"
nodes <- merge(nodes, baseline[, c("Cell_type", "Total")])
nodes$Rel_freq <- nodes$Ccount/nodes$Total
relfreqsum <-
  aggregate(nodes$Rel_freq,
            by = list(Node = nodes$Node),
            sum)
colnames(relfreqsum)[2] <- "Relfreqsum"
nodes <- merge(nodes, relfreqsum)
nodes$RF_norm <- nodes$Rel_freq/nodes$Relfreqsum
