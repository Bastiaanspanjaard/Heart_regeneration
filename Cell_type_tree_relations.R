# Dependencies ####
source("../Devtree/Scripts/Devtree_library.R")
library(reshape2)
require(pheatmap)
require(RColorBrewer)
require(igraph)

# tree_sample <- tree_list[[1]]
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

# node <- tree
# reference_set <- cells_in_tree
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

# Hr1 <- ReadTree("Hr1", reference_set = Hr1_ref)
# library_name <- "Hr1"
# reference_set <- Hr1_ref
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
                      c("Ery. Duplex", "M. Duplex"))
vertex_colors <- data.frame(celltype = zoom_types,
                            Order = 1:(length(zoom_types)))
vertex_colors <- merge(vertex_colors, type_colors[, c("celltype", "colo1")])
vertex_colors <- vertex_colors[order(vertex_colors$Order), ]
vertex_colors$Label <-
  c("Epi A", "Epi V", "cfd", "col11", "col12", "cxcl12", "prolif", "spock3", "F")
rownames(vertex_colors) <- vertex_colors$celltype
zoom_to <- c("Fibroblast (col11a1a)", "Fibroblast (col12a1a)", "Fibroblast (proliferating)",
             "Fibroblast (nppc)") #"Fibroblasts (spock3)", 
zoom_from <- setdiff(zoom_types, zoom_to)

# Load tree objects, append cell types, create lineage trees ####
# 3dpi
# Hr10 <- ReadTree("Hr10", reference_set = cell_types[cell_types$orig.ident == "Hr10", ])
# Hr11 <- ReadTree("Hr11", reference_set = cell_types[cell_types$orig.ident == "Hr11", ])
# Hr12 <- ReadTree("Hr12", reference_set = cell_types[cell_types$orig.ident == "Hr12", ])
# Hr24 <- ReadTree("Hr24", reference_set = cell_types[cell_types$orig.ident == "Hr24", ])
# Hr26 <- ReadTree("Hr26", reference_set = cell_types[cell_types$orig.ident == "Hr26", ])
# Hr27 <- ReadTree("Hr27", reference_set = cell_types[cell_types$orig.ident == "Hr27", ])
# tree_list <- list(Hr10 = Hr10, Hr11 = Hr11, Hr12 = Hr12, Hr24 = Hr24, Hr26 = Hr26, Hr27 = Hr27)
# 7dpi
Hr1_ref <- cell_types[cell_types$orig.ident == "Hr1", ]
Hr1_ref$Cell_name <- paste("nd", sapply(Hr1_ref$Cell_name, function(x){unlist(strsplit(x, "_"))[2]}), sep = "")
Hr1 <- ReadTree("Hr1", reference_set = Hr1_ref)
Hr2 <- ReadTree("Hr2", reference_set = cell_types[cell_types$orig.ident %in% c("Hr2a", "Hr2b"), ])
Hr13 <- ReadTree("Hr13", reference_set = cell_types[cell_types$orig.ident == "Hr13", ])
Hr14 <- ReadTree("Hr14", reference_set = cell_types[cell_types$orig.ident == "Hr14", ])
Hr15 <- ReadTree("Hr15", reference_set = cell_types[cell_types$orig.ident == "Hr15", ])
tree_list <- list(Hr1 = Hr1, Hr2 = Hr2, Hr13 = Hr13, Hr14 = Hr14, Hr15 = Hr15)

# Create tree visualization and zoom visualization
for(t in 1:length(tree_list)){
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Full_tree", ct_colors = type_colors$Color2)
  tree_list[[t]] <- MakePieTree(tree_list[[t]], "Fibrozoom_tree", types = zoom_types, 
                                ct_colors = type_colors$colo1)
  htmlwidgets::saveWidget(
    tree_list[[t]]$Full_tree,
    file = paste("~/Documents/Projects/heart_Bo/Images/tree_", 
                 names(tree_list)[t], "_LINNAEUS_pie.html", sep = ""))
  htmlwidgets::saveWidget(
    tree_list[[t]]$Fibrozoom_tree,
    file = paste("~/Documents/Projects/heart_Bo/Images/tree_",
                 names(tree_list)[t], "_LINNAEUS_pie_fibrozoom.html", sep = ""))
  # htmlwidgets::saveWidget(
  #   tree_list$Hr27$Full_tree,
  #   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie.html")
  # htmlwidgets::saveWidget(
  #   tree_list$Hr27$Fibrozoom_tree,
  #   file = "~/Documents/Projects/heart_Bo/Images/tree_Hr27_LINNAEUS_pie_fibrozoom.html")
}


# Find precursor suspects ####
nodes <- data.frame(Node = character(),
                    Tree = character(),
                    Cell_type = character(),
                    Ccount = integer(),
                    Cell_type_total = integer(),
                    Node_total = integer(),
                    Node_rf = numeric(),
                    Cell_type_rf = numeric())
for(t in 1:length(tree_list)){
  nodes_add <- count_cumulative(tree_list[[t]]$Tree)[, c("Node", "Cell_type", "Type_count")]
  nodes_add$Cell_type <- as.character(nodes_add$Cell_type)
  nodes_add$Tree <- names(tree_list)[t]
  
  full_node_counts <-
    aggregate(nodes_add$Type_count, by = list(Node = nodes_add$Node), sum)
  colnames(full_node_counts)[2] <- c("Node_total")
  nodes_add <- merge(nodes_add, full_node_counts)
  nodes_add$Node_rf <- nodes_add$Type_count/nodes_add$Node_total
  
  baseline <- nodes_add[nodes_add$Node == "nd0", c("Node", "Cell_type", "Type_count")]
  colnames(baseline)[3] <- "Cell_type_total"
  
  nodes_add <- merge(nodes_add, baseline[, c("Cell_type", "Cell_type_total")])
  nodes_add$Cell_type_rf <- nodes_add$Type_count/nodes_add$Cell_type_total
  
  nodes <- rbind(nodes, nodes_add)
}
nodes$Treenode <- paste(nodes$Tree, nodes$Node, sep = "_")

nodes_nfc <- acast(nodes[nodes$Node != "nd0", ], Treenode ~ Cell_type, value.var ="Node_rf")
nodes_nfc[is.na(nodes_nfc)] <- 0

ggplot(data.frame(nodes_nfc)) +
  geom_point(aes_string(x = "Endocardium..nppc..V", y = "Fibroblast..nppc."))
ggplot(data.frame(nodes_nfc)) +
  geom_point(aes_string(x = "Endocardium..Ventricle.", y = "Fibroblast..nppc."))
# png("./Images/Col11fib_3dpi_myofibr_node_prop_example.png", width = 768, height = 768)
# ggplot(data.frame(nodes_nfc)) +
#   geom_point(aes_string(x = "Myofibroblasts", y = "Fibroblast..col11a1a."), size = 3) +
#   labs(title = "Proportion of node", x = "Myofibroblasts", y = "Fibroblast (col11a1a)")
# dev.off()

cortest_node_freqs <- data.frame(Type_1 = character(),
                                 Type_2 = character(),
                                 ct_p = numeric(),
                                 ct_padj = numeric())
for(i in 1:ncol(nodes_nfc)){
  cortest_node_add <- data.frame(Type_1 = character(ncol(nodes_nfc)),
                                 Type_2 = character(ncol(nodes_nfc)),
                                 ct_p = numeric(ncol(nodes_nfc)),
                                 ct_padj = numeric(ncol(nodes_nfc)), stringsAsFactors = F)
  for(j in 1:ncol(nodes_nfc)){
    ct <- cor.test(x = nodes_nfc[, i], y = nodes_nfc[, j])
    cortest_node_add[j, 1:3] <- c(colnames(nodes_nfc)[i], 
                                  colnames(nodes_nfc)[j], ct$p.value)
  }
  cortest_node_add$ct_padj <- p.adjust(cortest_node_add$ct_p, method = "fdr")
  cortest_node_freqs <- 
    rbind(cortest_node_add, cortest_node_freqs)
}

for(i in 1:length(zoom_to)){
  cortest_partplot <- cortest_node_freqs[cortest_node_freqs$Type_1 == zoom_to[i] &
                                           cortest_node_freqs$Type_2 != zoom_to[i], ]
  cortest_partplot <- cortest_partplot[order(cortest_partplot$ct_padj), ]
  cortest_partplot$Type_2 <- factor(cortest_partplot$Type_2, levels = cortest_partplot$Type_2)
  # png(paste("./Images/", zoom_to[i], "precursor_suspects_7dpi.png", sep = ""), width = 768, height = 768)
  print(
    ggplot(cortest_partplot[cortest_partplot$ct_padj < 0.01, ]) +
      geom_bar(aes(x = Type_2, y = -log10(ct_padj), fill = Type_2), stat = "identity") +
      scale_fill_manual(values = ann_colors$Celltype) +
      labs(fill = "", x = "", y = "p (-log10)", title = paste(zoom_to[i], "precursor suspects")) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  )
  # dev.off()
}

# Analyse evidence against precursors ####
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

# png("./Images/Col11fib_3dpi_potential_precursors.png", width = 1366, height = 768)
# ggplot(full_descendancy[full_descendancy$Cell_type == "Fibroblast (col11a1a)", ]) +
#   geom_jitter(aes(x = 0, y = Precursor_presence_p, color = Precursor), size = 2) +
#   labs(title = "Fibroblast (col11a1a) precursors",
#        x = "", y = "Probability") +
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = "none") +
#   scale_color_manual(values = ann_colors$Celltype) +
#   facet_wrap(~Precursor) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_text(size = 12, face = "bold")) +
#   scale_y_continuous(breaks = c(0, 0.5, 1))
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
# ggplot(full_descendancy[full_descendancy$Cell_type == "Fibroblast (col12a1a)" &
#                           full_descendancy$Precursor %in% zoom_from, ]) +
#   geom_jitter(aes(x = 0, y = log10(Precursor_presence_p), color = Precursor)) +
#   labs(title = "Fibroblast (col12a1a) precursors",
#        y = "Probability (log10)", x = "") +
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = "none") +
#   scale_color_manual(values = ann_colors$Zoomtype) +
#   facet_wrap(~Precursor, nrow = 1) +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_text(size = 12, face = "bold"))
# dev.off()  

for(i in 1:length(zoom_to)){
  potential_prec <- cortest_node_freqs$Type_2[cortest_node_freqs$Type_1 == zoom_to[i] &
                                       cortest_node_freqs$Type_2 != zoom_to[i] &
                                       cortest_node_freqs$ct_padj < 0.01]
  # png(paste("./Images/", zoom_to[i], "precursor_node_probs_7dpi.png", sep = ""),
  #     width = 768, height = 256 * ceiling(length(potential_prec)/5))
  print(
    ggplot(full_descendancy[full_descendancy$Cell_type == zoom_to[i] &
                              full_descendancy$Precursor %in% potential_prec, ]) +
      geom_jitter(aes(x = 0, y = log10(Precursor_presence_p), color = Precursor)) +
      labs(title = paste(zoom_to[i], "precursors"),
           y = "Probability (log10)", x = "") +
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none") +
      scale_color_manual(values = ann_colors$Celltype) +
      facet_wrap(~Precursor, ncol = 5) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text.x = element_text(size = 12, face = "bold"))
  )
  # dev.off()
}

# View(nodes_nfc[, colnames(nodes_nfc) %in% c("Fibroblast (nppc)", "Endocardium (Ventricle)")])

# Distance between A and A+B ####
target <- "Fibroblast (nppc)"
# target <- "Fibroblast (col12a1a)"
# Run this over all trees, keep distances top-level, create a new tree-level in the list
comparison_list <- list(Comparison = data.frame(Tree = character(),
                                                Precursor = character(),
                                                Simple_distance = numeric(),
                                                Distance = numeric(),
                                                Correlation = numeric(),
                                                Type_count = integer()),
                        Node_sizes = data.frame(Size = integer()),
                        Normalized_frequencies = data.frame(matrix(nrow = 0, ncol = length(unique(cell_types$Cell_type)))))
colnames(comparison_list$Normalized_frequencies) <- unique(cell_types$Cell_type)
for(t in 1:length(tree_list)){
  tree <- tree_list[[t]]$Tree
  
  edge_list <- ToDataFrameNetwork(tree, "Cell.type")
  if(!(target %in% edge_list$Cell.type)){
    next
  }
  sample_type_count <- 
    dcast(data.frame(table(edge_list$from, edge_list$Cell.type)), Var1 ~ Var2, value.var = "Freq")
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
  # new_nodesizes <- 
  comparison_list$Node_sizes <- rbind(comparison_list$Node_sizes,
                                      data.frame(Size = rowSums(sample_type_cumulative)))
  # rownames(sample_type_nf) <- paste(names(tree_list)[t], rownames(sample_type_nf), sep = ":")
  
  # comparison_list_tree <- 
    # list(Distances = 
  comparison_tree <- 
    data.frame(Tree = rep(names(tree_list)[t], ncol(sample_type_count)),
               Precursor = names(sample_type_count),
               Distance = numeric(ncol(sample_type_count)),
               Correlation = numeric(ncol(sample_type_count)),
               Type_count = colSums(sample_type_count))
  
  precursor_list <- list()
  for(p in 1:ncol(sample_type_count)){#names(sample_type_count))
    precursor <- names(sample_type_count)[p]#"Endocardium (Ventricle)"
    
    correlations <- 
      data.frame(Precursor_cor = 
                   t(wtd.cors(sample_type_nf[[precursor]], 
                              sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
                              weight = rowSums(sample_type_cumulative))),
                 Both_cor = t(wtd.cors(sample_type_nf[[precursor]] + sample_type_nf[[target]], 
                                       sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
                                       weight = rowSums(sample_type_cumulative))))
    comparison_tree$Distance[p] <- dist(t(correlations))
    comparison_tree$Correlation[p] <- cor(correlations$Precursor_cor, correlations$Both_cor)
    precursor_list[[p]] <- correlations
    names(precursor_list)[p] <- precursor
  }
  
  comparison_list[[length(comparison_list) + 1]] <- precursor_list
  names(comparison_list)[length(comparison_list)] <- names(tree_list)[t]
  comparison_list$Comparison <- rbind(comparison_list$Comparison, comparison_tree)
  
  sample_type_nf[setdiff(names(comparison_list$Normalized_frequencies), names(sample_type_nf))] <- NA
  comparison_list$Normalized_frequencies[setdiff(names(sample_type_nf), 
                                                 names(comparison_list$Normalized_frequencies))] <- NA
  
  comparison_list$Normalized_frequencies <- 
    rbind(comparison_list$Normalized_frequencies, sample_type_nf)
}
comparison_list$Normalized_frequencies[is.na(comparison_list$Normalized_frequencies)] <- 0
comparison_list$Node_sizes$Weight <- comparison_list$Node_sizes$Size/sum(comparison_list$Node_sizes$Size)
# node_weights <- comparison_list$Comparison$Type_count
comparison_list$All_trees_prec <- list()
comparison_list$All_trees_distances <- 
  data.frame(Precursor = names(comparison_list$Normalized_frequencies),
             Distance = numeric(ncol(comparison_list$Normalized_frequencies)))

View(comparison_list$Node_sizes)
for(p in 1:ncol(comparison_list$Normalized_frequencies)){
  precursor <- names(comparison_list$Normalized_frequencies)[p]
  
  if(precursor == target){
    comparison_list$All_trees_distances$Distance[p] <- 0
    next
  }
  full_node_distance <-
    data.frame(Precursor = comparison_list$Normalized_frequencies[, names(comparison_list$Normalized_frequencies) == precursor],
               Prec_prog = rowSums(comparison_list$Normalized_frequencies[, names(comparison_list$Normalized_frequencies) 
                                                                          %in% c(precursor, target)]))
  
  comparison_list$All_trees_prec[[p]] <- full_node_distance
  names(comparison_list$All_trees_prec)[p] <- precursor
  
  comparison_list$All_trees_distances$Distance[p] <-
    sqrt(sum((comparison_list$Node_sizes$Weight * (full_node_distance$Precursor - full_node_distance$Prec_prog))^2))
}  
  
# Precursor
  comparison_list$Normalized_frequencies[, names(comparison_list$Normalized_frequencies) == precursor]
  
  # Progenitor
  rowSums(comparison_list$Normalized_frequencies[, names(comparison_list$Normalized_frequencies) 
                                                 %in% c(precursor, target)])
  
  
  correlations <- 
    data.frame(Precursor_cor = 
                 t(wtd.cors(sample_type_nf[[precursor]], 
                            sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
                            weight = rowSums(sample_type_cumulative))),
               Both_cor = t(wtd.cors(sample_type_nf[[precursor]] + sample_type_nf[[target]], 
                                     sample_type_nf[, !(names(sample_type_nf) %in% c(precursor, target))], 
                                     weight = rowSums(sample_type_cumulative))))
  comparison_tree$Distance[p] <- dist(t(correlations))
  comparison_tree$Correlation[p] <- cor(correlations$Precursor_cor, correlations$Both_cor)
  precursor_list[[p]] <- correlations
  names(precursor_list)[p] <- precursor
}



ppp <- "Endocardium (Ventricle)"

ggplot(comparison_list[[ppp]]$Correlations) +
  geom_point(aes(x = Precursor_cor, y = Both_cor)) +
  labs(title = paste(ppp, "giving rise to", target))


ggplot(comparison_list$`Endocardium (nppc) V`$Correlations) +
  geom_point(aes(x = Precursor_cor, y = Both_cor))
ggplot(comparison_list$`Endocardium (nppc) V`$Correlations) +
  geom_point(aes(x = Precursor_cor, y = Both_cor))



cor_samples$Correlation[s] <- cor(sample_type_nf[[split_type_0]], sample_type_nf[[split_type_1]])
cor_samples$Wt_correlation[s] <- 
  wtd.cors(sample_type_nf[[split_type_0]], 
           sample_type_nf[[split_type_1]], weight = rowSums(sample_type_cumulative))
cor_samples$Double_correlation[s] <- 
  cor(t(cor(sample_type_nf[[split_type_0]], 
            sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))])),
      t(cor(sample_type_nf[[split_type_1]], 
            sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))])))
cor_samples$Double_wt_correlation[s] <- 
  cor(t(wtd.cors(sample_type_nf[[split_type_0]], 
                 sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))],
                 weight = rowSums(sample_type_cumulative))),
      t(wtd.cors(sample_type_nf[[split_type_1]], 
                 sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))],
                 weight = rowSums(sample_type_cumulative))))


# correlation between A and other non-A,B types
# correlation between A + B and other non-A,B types
# distance between correlations
  

# Correlation bootstrap ####
# Randomly split a cell type in two within a tree: create two artificial sub types by assigning
# all cells in a type a random 'A' or 'B'. Then calculate the tree-based correlation between those
# two.
require(weights)
cell_type_split <- "Endocardium (Ventricle)"
tree <- Clone(tree_list$Hr2$Tree)
samples <- 1000

cor_samples <- data.frame(Sample = 1:samples,
                          Correlation = numeric(samples),
                          Wt_correlation = numeric(samples),
                          Double_correlation = numeric(samples),
                          Double_wt_correlation = numeric(samples))
edge_list <- ToDataFrameNetwork(tree, "Cell.type")
# node_totals <- data.frame(table(edge_list$from[edge_list$Cell.type != "NA"]))
# node_totals$Node_total <- NA
set.seed(555)
for(s in 1:samples){
  # NEW
  edge_list$New_cell_type <- 
    sapply(edge_list$Cell.type,
           function(x){
             ifelse(x == cell_type_split,
                    paste(x, rbinom(1, 1, 0.5), sep = "_"), x) #,
                    # ifelse(x == "NA", "NA", "Other"))
  })
  split_type_0 <- paste(cell_type_split, 0, sep = "_")
  split_type_1 <- paste(cell_type_split, 1, sep = "_")
  # nodes <- edge_list$to[edge_list$Cell.type == "NA"] # Intentionally leaves out nd0
  sample_type_count <- 
    dcast(data.frame(table(edge_list$from, edge_list$New_cell_type)), Var1 ~ Var2, value.var = "Freq")
  sample_type_count$Var1 <- as.character(sample_type_count$Var1)
  # NEW
  sample_type_count <- sample_type_count[, colnames(sample_type_count) != "NA"]
  rownames(sample_type_count) <- sample_type_count$Var1
  sample_type_count <- sample_type_count[, -1]
  
  sample_type_cumulative <- sample_type_count
  sample_type_cumulative[,] <- NA
  for(i in 1:nrow(sample_type_cumulative)){
    node <- rownames(sample_type_cumulative)[i]
    sample_type_cumulative[i, ] <-
      colSums(sample_type_count[grep(node, rownames(sample_type_count)), ])
  }
  sample_type_nf <- sample_type_cumulative/rowSums(sample_type_cumulative)
  cor_samples$Correlation[s] <- cor(sample_type_nf[[split_type_0]], sample_type_nf[[split_type_1]])
  cor_samples$Wt_correlation[s] <- 
    wtd.cors(sample_type_nf[[split_type_0]], 
             sample_type_nf[[split_type_1]], weight = rowSums(sample_type_cumulative))
  cor_samples$Double_correlation[s] <- 
    cor(t(cor(sample_type_nf[[split_type_0]], 
              sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))])),
        t(cor(sample_type_nf[[split_type_1]], 
              sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))])))
  cor_samples$Double_wt_correlation[s] <- 
    cor(t(wtd.cors(sample_type_nf[[split_type_0]], 
              sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))],
              weight = rowSums(sample_type_cumulative))),
        t(wtd.cors(sample_type_nf[[split_type_1]], 
              sample_type_nf[, !(colnames(sample_type_nf) %in% c(split_type_0, split_type_1))],
              weight = rowSums(sample_type_cumulative))))
  
  # OLD
  # sample_type_count$Fibroblasts_0_cum <- NA
  # sample_type_count$Fibroblasts_1_cum <- NA
  # sample_type_count$Other_cum <- NA
  # for(i in 1:nrow(sample_type_count)){
  #   node <- sample_type_count$Var1[i]
  #   sample_type_count[i, c("Fibroblasts_0_cum", "Fibroblasts_1_cum", "Other_cum")] <-
  #     colSums(sample_type_count[grep(node, sample_type_count$Var1), c(2, 3, 5)])
  # }
  # sample_type_count$Node_total <- 
  #   sample_type_count$Fibroblasts_0_cum + sample_type_count$Fibroblasts_1_cum +
  #   sample_type_count$Other_cum
  # sample_type_count$F_0_rf <- sample_type_count$Fibroblasts_0_cum/sample_type_count$Node_total
  # sample_type_count$F_1_rf <- sample_type_count$Fibroblasts_1_cum/sample_type_count$Node_total
  # cor_samples$Correlation[s] <- cor(sample_type_count$F_0_rf, sample_type_count$F_1_rf)
  # cor_samples$Wt_correlation[s] <- 
  #   wtd.cors(sample_type_count$F_0_rf, sample_type_count$F_1_rf, weight = sample_type_count$Node_total)
  # OLDER
  # cell_division_start <- Sys.time()
  # cells_in_tree <- data.frame(Cell_name = tree$Get('name', filterFun = function(x) {isLeaf(x)}),
  #                             stringsAsFactors = F)
  # cells_in_tree <- merge(cells_in_tree, cell_types)
  # cells_in_tree$Cell_type <- 
  #   sapply(cells_in_tree$Cell_type,
  #          function(x){
  #            ifelse(x == cell_type_split,
  #                   paste(x, rbinom(1, 1, 0.5), sep = "_"),
  #                   x)})
  # tree$Do(RemapCellTypes, reference_set = cells_in_tree)
  # cell_division_time <- Sys.time() - cell_division_start
  # 
  # counting_start <- Sys.time()
  # nodes <- count_cumulative(tree)#[, c("Node", "Cell_type", "Ccount")]
  # nodes$Cell_type <- as.character(nodes$Cell_type)
  # counting_time <- Sys.time() - counting_start
  # 
  # correlation_start <- Sys.time()
  # full_node_counts <-
  #   aggregate(nodes$Type_count, by = list(Node = nodes$Node), sum)
  # colnames(full_node_counts)[2] <- c("Node_total")
  # nodes <- merge(nodes, full_node_counts)
  # nodes$Node_rf <- nodes$Type_count/nodes$Node_total
  # 
  # # Now calculate correlation.
  # nodes_nfc <- dcast(nodes[nodes$Node != "nd0", ], Node ~ Cell_type, value.var ="Node_rf")
  # nodes_nfc[is.na(nodes_nfc)] <- 0
  # cor_samples$Correlation[s] <- cor(nodes_nfc$Fibroblasts_0, nodes_nfc$Fibroblasts_1)
  # correlation_time <- Sys.time() - correlation_start
  # 
  # print(paste("Sample step", s, "with", cell_division_time, ",", counting_time, "and", correlation_time))
}
ggplot(cor_samples) +
  geom_histogram((aes(x = Correlation)))
ggplot(cor_samples) +
  geom_histogram((aes(x = Wt_correlation)))
ggplot(cor_samples) +
  geom_histogram((aes(x = Double_correlation)))
ggplot(cor_samples) +
  geom_histogram((aes(x = Double_wt_correlation)))
