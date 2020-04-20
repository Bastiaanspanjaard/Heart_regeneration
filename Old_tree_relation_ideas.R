# Parts of old cell type tree relation analysis. Probably no longer needed.

# # Analyse evidence against precursors ####
# full_descendancy <- data.frame(Cell_type = character(),
#                                Precursor = character(),
#                                Node = character(),
#                                Parent_node = character(),
#                                Precursor_count = numeric(),
#                                Total = numeric(),
#                                Frequency = numeric(),
#                                Precursor_presence_p = numeric(),
#                                Entropy = numeric(),
#                                Tree = character())
# for(t in 1:length(tree_list)){
#   tree_list[[t]] <- CalculateCooccurrence(tree_list[[t]])
#   edge_list <- ToDataFrameNetwork(tree_list[[t]]$Tree, "Cell.type")
#   if(!(target %in% edge_list$Cell.type) | sum(edge_list$Cell.type == target) < 20){
#     next
#   }
#   
#   descendancy_add <- tree_list[[t]]$Descendancy
#   descendancy_add$Tree <- names(tree_list)[t]
#   descendancy_add <- merge(descendancy_add, tree_list[[t]]$Node_entropy)
#   full_descendancy <- rbind(full_descendancy, descendancy_add)
# }
# 
# full_descendancy_plot <- full_descendancy[full_descendancy$Cell_type == target &
#                                             full_descendancy$Precursor %in% potential_prec, ]
# full_descendancy_plot$Precursor <- factor(full_descendancy_plot$Precursor, levels = precursor_ranking$Precursor)
# full_descendancy_plot$Log10p_bottom <- 
#   sapply(full_descendancy_plot$Precursor_presence_p, function(x) max(log10(x), -10))
# 
# x <- length(potential_prec)
# # png(paste("./Images/", target, "precursor_node_probs_37dpi.png", sep = ""),
# # width = 768, height = 256 * ceiling(length(potential_prec)/5))
# # png(paste("./Images/", target, "precursor_node_probs_37dpi.png", sep = ""),
# #     width = max(512, 384 * ceiling(x/5)), height = 256 * min(x, min(3, x)))
# print(
#   # ggplot(full_descendancy[full_descendancy$Cell_type == zoom_to[i] &
#   ggplot(full_descendancy_plot) +
#     geom_jitter(aes(x = 0, y = Log10p_bottom, color = Precursor), height = 0, size = 3) +
#     # scale_y_continuous(limits = c(-10.1, 0.1)) + 
#     # labs(title = paste(zoom_to[i], "precursors"),
#     labs(title = paste(target, "precursors"),
#          y = "Probability (log10)", x = "") +
#     theme(axis.text.x = element_text(angle = 90),
#           legend.position = "none") +
#     scale_color_manual(values = ann_colors$Celltype) +
#     facet_wrap(~Precursor, nrow = 3) +
#     theme(axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           strip.text.x = element_text(size = 12, face = "bold"))
# )
# # dev.off()
# 
# 
# 
# # Calculate distance between cell types in a tree ####
# t <- 1
# tree <- tree_list[[t]]$Tree
# # ** Compute background tree distance distribution ####
# # Calculate node distances
# edge_list <- ToDataFrameNetwork(tree, "Cell.type")
# nodes <- setdiff(unique(c(edge_list$from, edge_list$to[edge_list$Cell.type == "NA"])), "nd0")
# node_distances <- data.frame(t(combn(nodes, 2)))
# colnames(node_distances) <- c("Node_1", "Node_2")
# node_distances$Distance <- 
#   apply(node_distances[, c("Node_1", "Node_2")], 1,
#         function(x){
#           v_1 <- unlist(strsplit(as.character(x[1]), "_"))
#           v_2 <- unlist(strsplit(as.character(x[2]), "_"))
#           # For comparison, add "-1"'s to the shortest vector
#           v_1f <- c(v_1, rep(-1, max(length(v_2) - length(v_1), 0)))
#           v_2f <- c(v_2, rep(-1, max(length(v_1) - length(v_2), 0)))
#           # Distance = d(str_1 - mcra) + d(str_2 - mcra) =
#           #            d(str_1) - d(mcra) + d(str_2) - d(mcra) =
#           #            d(str_1) + d(str_2) - 2*d(mcra)
#           return(length(v_1) + length(v_2) - 2*sum(cumprod(v_1f == v_2f)))
#         })
# node_distances <- rbind(node_distances,
#                         data.frame(Node_1 = nodes,
#                                    Node_2 = nodes,
#                                    Distance = 0))
# # Calculate node chances
# pure_sizes <- aggregate(pure_counts$Type_count,
#                         by = list(Node = pure_counts$Node),
#                         sum)
# colnames(pure_sizes)[2] <- "Node_count"
# pure_sizes <- pure_sizes[pure_sizes$Node != "nd0", ]
# pure_sizes$Chance <- pure_sizes$Node_count/sum(pure_sizes$Node_count)
# # Calculate distance chances
# node_distances <- merge(node_distances, pure_sizes[, c("Node", "Chance")],
#                         by.x = "Node_1", by.y = "Node")
# colnames(node_distances)[4] <- "Chance_1"
# node_distances <- merge(node_distances, pure_sizes[, c("Node", "Chance")],
#                         by.x = "Node_2", by.y = "Node")
# colnames(node_distances)[5] <- "Chance_2"
# node_distances$Distance_chance <- node_distances$Chance_1 * node_distances$Chance_2
# node_distances$Distance_chance[node_distances$Node_1 != node_distances$Node_2] <-
#   2*node_distances$Distance_chance[node_distances$Node_1 != node_distances$Node_2]
# # Calculate distance distribution and average
# distance_distribution <- 
#   aggregate(node_distances$Distance_chance,
#           by = list(Distance = node_distances$Distance),
#           sum)
# colnames(distance_distribution)[2] <- "Chance"
# ggplot(distance_distribution) +
#   geom_bar(aes(x = Distance, y = Chance), stat = "identity") +
#   labs(title = paste("Distance distribution for tree ", names(tree_list)[t],
#                      ", average ", 
#                      round(weighted.mean(distance_distribution$Distance, distance_distribution$Chance), 2),
#                      sep = ""))
#        
# # ** Calculate distance distribution for cell type pairs ####
# cell_type_distances <-
#   data.frame(t(combn(setdiff(unique(edge_list$Cell.type), "NA"), 2)))
# colnames(cell_type_distances) <- c("Type_1", "Type_2")
# cell_type_distances <-
#   rbind(cell_type_distances,
#         data.frame(Type_1 = setdiff(unique(edge_list$Cell.type), "NA"),
#            Type_2 = setdiff(unique(edge_list$Cell.type), "NA")))
# cell_type_distances$Distance <- -1
# cell_type_distances$Measurements <- 0
# # d <- 1126
# for(d in 1:nrow(cell_type_distances)){
#   print(paste(d, " out of ", nrow(cell_type_distances), sep = ""))
#   type_1 <- cell_type_distances$Type_1[d]
#   type_2 <- cell_type_distances$Type_2[d]
#   # Calculate distances - node 1, node 2, distance. Note the number of comparisons is m*n for
#   # two different cell types, but 0.5*m(m-1) for inter-cell type distances.
#   if(type_1 == type_2){
#     type_positions <- edge_list[edge_list$Cell.type == type_1 &
#                                     edge_list$from != "nd0", ]
#     colnames(type_positions)[1:2] <- c("Node", "Cell")
#     if(nrow(type_positions) > 1){
#       type_distances <- data.frame(t(combn(type_positions$Cell, 2)))
#       colnames(type_distances) <- c("Cell_1", "Cell_2")
#       type_distances <- merge(type_distances, type_positions[, c("Node", "Cell")],
#                               by.x = "Cell_1", by.y = "Cell")
#       colnames(type_distances)[3] <- "Node_1"
#       type_distances <- merge(type_distances, type_positions[, c("Node", "Cell")],
#                               by.x = "Cell_2", by.y = "Cell")
#       colnames(type_distances)[4] <- "Node_2"
#       type_distances <- merge(type_distances, node_distances[, c("Node_1", "Node_2", "Distance")])
#       cell_type_distances$Distance[d] <- mean(type_distances$Distance)
#       cell_type_distances$Measurements[d] <- nrow(type_distances)
#     }
#   }else{
#     type_1_positions <- edge_list[edge_list$Cell.type == type_1 &
#                                     edge_list$from != "nd0", ]
#     colnames(type_1_positions) <- c("Node_1", "Cell_1", "Type_1")
#     type_2_positions <- edge_list[edge_list$Cell.type == type_2 &
#                                     edge_list$from != "nd0", ]
#     colnames(type_2_positions) <- c("Node_2", "Cell_2", "Type_2")
#     type_distances <- merge(type_1_positions, type_2_positions)
#     type_distances <- merge(type_distances, node_distances[, c("Node_1", "Node_2", "Distance")])
#     cell_type_distances$Distance[d] <- mean(type_distances$Distance)
#     cell_type_distances$Measurements[d] <- nrow(type_distances)
#   }
# }  
# 
# # ** Compare distances for cell type pairs to background distribution by calculating z-values ####
# mean_dist <- weighted.mean(distance_distribution$Distance, distance_distribution$Chance)
# distance_distribution$x_m <- distance_distribution$Distance - mean_dist
# dist_sd <- sqrt(sum(distance_distribution$Chance * distance_distribution$x_m^2)/nrow(distance_distribution))
# cell_type_distances$Z_distance <- 
#   sqrt(cell_type_distances$Measurements) * (cell_type_distances$Distance - mean_dist)/dist_sd
# View(cell_type_distances[cell_type_distances$Type_1 == cell_type_distances$Type_2, ])
# 
