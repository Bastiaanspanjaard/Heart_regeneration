require(data.table)
require(igraph)
require(reshape2)
require(ggplot2)

# Load data ####
ctrl <- fread("./Data/Ligand_receptor/LigRec_Ctrl_zoom_significant_means.txt")
for(col in names(ctrl)[-(1:12)]){
  set(ctrl, j = col, value = as.numeric(ctrl[[col]]))
}
inj_3dpi <- fread("./Data/Ligand_receptor/LigRec_3dpi_zoom_significant_means.txt")
for(col in names(inj_3dpi)[-(1:12)]){
  set(inj_3dpi, j = col, value = as.numeric(inj_3dpi[[col]]))
}
inj_7dpi <- fread("./Data/Ligand_receptor/LigRec_7dpi_zoom_significant_means.txt")
for(col in names(inj_7dpi)[-(1:12)]){
  set(inj_7dpi, j = col, value = as.numeric(inj_7dpi[[col]]))
}
ctrl_all <- fread("./Data/Ligand_receptor/LigRec_Ctrl_zoom_means.txt")
inj_3dpi_all <- fread("./Data/Ligand_receptor/LigRec_3dpi_zoom_means.txt")
inj_7dpi_all <- fread("./Data/Ligand_receptor/LigRec_7dpi_zoom_means.txt")
# ctrl_ps <- fread("./Data/Ligand_receptor/LigRec_Ctrl_zoom_pvalues.txt")
ctrl_dec <- fread("./Data/Ligand_receptor/LigRec_Ctrl_zoom_deconvoluted.txt")
dec_3dpi <- fread("./Data/Ligand_receptor/LigRec_3dpi_zoom_deconvoluted.txt")
dec_7dpi <- fread("./Data/Ligand_receptor/LigRec_7dpi_zoom_deconvoluted.txt")

cell_types <- fread("./Data/celltypes_zoom_allcells.csv",
                    header = T)
cell_types <- cell_types[, c("V1", "time", "zoom.subtypes")]
colnames(cell_types) <- c("Cell", "Time", "Cell_type")
cell_type_freqs <- cell_types[, .(Count = .N), by = list(Time, Cell_type)]
cell_type_freqs <- cell_type_freqs[cell_type_freqs$Time %in% c("Ctrl", "3dpi", "7dpi")]

# Which genes can we expect to find?
human_zfish_translation <- fread("./Data/GRCz11_GRCh38_orthologues.txt")
wnt_ra_gene_translation <- 
  merge(wnt_in, human_zfish_translation[human_zfish_translation$`Human gene name` != "", ], 
        by.x = "gene", by.y = "Gene name", all.x = T)

interactions_profiles <- 
fread("/Users/bastiaanspanjaard/Documents/Projects/heart_Bo/Scripts/CellPhoneDB/cpdb-venv/lib/python3.6/site-packages/cellphonedb/src/core/data/interaction_input.csv")

# Compare interactions ####
interactions <- data.frame(Pair = ctrl$interacting_pair[apply(ctrl[, -(1:12)], 1, function(x){sum(!is.na(x))>0})])
interactions$Freq_ctrl <- apply(ctrl[, -(1:12)], 1, function(x){sum(!is.na(x))})[apply(ctrl[, -(1:12)], 1, function(x){sum(!is.na(x))>0})]
interactions_3dpi <- data.frame(Pair = inj_3dpi$interacting_pair[apply(inj_3dpi[, -(1:12)], 1, function(x){sum(!is.na(x))>0})])
interactions_3dpi$Freq_3dpi <- apply(inj_3dpi[, -(1:12)], 1, function(x){sum(!is.na(x))})[apply(inj_3dpi[, -(1:12)], 1, function(x){sum(!is.na(x))>0})]
interactions_7dpi <- data.frame(Pair = inj_7dpi$interacting_pair[apply(inj_7dpi[, -(1:12)], 1, function(x){sum(!is.na(x))>0})])
interactions_7dpi$Freq_7dpi <- apply(inj_7dpi[, -(1:12)], 1, function(x){sum(!is.na(x))})[apply(inj_7dpi[, -(1:12)], 1, function(x){sum(!is.na(x))>0})]
interactions <-
  merge(interactions,
        merge(interactions_3dpi, interactions_7dpi, all = T), all = T)
interactions[is.na(interactions)] <- 0
interactions$ctrl_3dpi_ratio <- interactions$Freq_3dpi/(interactions$Freq_ctrl + interactions$Freq_3dpi)
interactions$ctrl_7dpi_ratio <- interactions$Freq_7dpi/(interactions$Freq_ctrl + interactions$Freq_7dpi)

# Wrangle into long form frame ####
# ctrl_subcol <- colnames(ctrl)[c(1:2, 13:ncol(ctrl))]
# ctrl_int <- melt(ctrl[, ..ctrl_subcol],
#                  id.vars = c("id_cp_interaction", "interacting_pair"),
#                  measure.vars = ctrl_subcol[-(1:2)],
#                  variable.name = "Type_pair", value.name = "Ctrl_mean")
# ctrl_int <- ctrl_int[!is.na(ctrl_int$Ctrl_mean)]
# ctrl_int$Type_pair <- as.character(ctrl_int$Type_pair)
# ctrl_int$Partner_1 <- sapply(ctrl_int$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
# ctrl_int$Partner_2 <- sapply(ctrl_int$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})
# ctrl_int <- merge(ctrl_int, cell_type_freqs[Time == "Ctrl", c("Cell_type", "Count")], by.x = "Partner_1", by.y = "Cell_type")
# colnames(ctrl_int)[7] <- "Ctrl_count_1"
# ctrl_int <- merge(ctrl_int, cell_type_freqs[Time == "Ctrl", c("Cell_type", "Count")], by.x = "Partner_2", by.y = "Cell_type")
# colnames(ctrl_int)[8] <- "Ctrl_count_2"

ctrl_all_subcol <- colnames(ctrl_all)[c(2, 12:ncol(ctrl_all))]
ctrl_all_int <- melt(ctrl_all[, ..ctrl_all_subcol],
                 id.vars = c("interacting_pair"),
                 measure.vars = ctrl_all_subcol[-1],
                 variable.name = "Type_pair", value.name = "Ctrl_mean")
ctrl_all_int <- ctrl_all_int[!is.na(ctrl_all_int$Ctrl_mean)]
ctrl_all_int$Type_pair <- as.character(ctrl_all_int$Type_pair)
ctrl_all_int$Partner_1 <- sapply(ctrl_all_int$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
ctrl_all_int$Partner_2 <- sapply(ctrl_all_int$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})

# subcol_3dpi <- colnames(inj_3dpi)[c(1:2, 13:ncol(inj_3dpi))]
# int_3dpi <- melt(inj_3dpi[, ..subcol_3dpi],
#                  id.vars = c("id_cp_interaction", "interacting_pair"),
#                  measure.vars = subcol_3dpi[-(1:2)],
#                  variable.name = "Type_pair", value.name = "Mean_3dpi")
# int_3dpi <- int_3dpi[!is.na(int_3dpi$Mean_3dpi)]
# int_3dpi$Type_pair <- as.character(int_3dpi$Type_pair)
# int_3dpi$Partner_1 <- sapply(int_3dpi$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
# int_3dpi$Partner_2 <- sapply(int_3dpi$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})
# int_3dpi <- merge(int_3dpi, cell_type_freqs[Time == "3dpi", c("Cell_type", "Count")], by.x = "Partner_1", by.y = "Cell_type")
# colnames(int_3dpi)[7] <- "Threedpi_count_1"
# int_3dpi <- merge(int_3dpi, cell_type_freqs[Time == "3dpi", c("Cell_type", "Count")], by.x = "Partner_2", by.y = "Cell_type")
# colnames(int_3dpi)[8] <- "Threedpi_count_2"

subcol_3dpi_all <- colnames(inj_3dpi_all)[c(2, 12:ncol(inj_3dpi_all))]
int_3dpi_all <- melt(inj_3dpi_all[, ..subcol_3dpi_all],
                 id.vars = c("interacting_pair"),
                 measure.vars = subcol_3dpi_all[-1],
                 variable.name = "Type_pair", value.name = "Mean_3dpi")
int_3dpi_all <- int_3dpi_all[!is.na(int_3dpi_all$Mean_3dpi)]
int_3dpi_all$Type_pair <- as.character(int_3dpi_all$Type_pair)
int_3dpi_all$Partner_1 <- sapply(int_3dpi_all$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
int_3dpi_all$Partner_2 <- sapply(int_3dpi_all$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})

# subcol_7dpi <- colnames(inj_7dpi)[c(1:2, 13:ncol(inj_7dpi))]
# int_7dpi <- melt(inj_7dpi[, ..subcol_7dpi],
#                  id.vars = c("id_cp_interaction", "interacting_pair"),
#                  measure.vars = subcol_7dpi[-(1:2)],
#                  variable.name = "Type_pair", value.name = "Mean_7dpi")
# int_7dpi <- int_7dpi[!is.na(int_7dpi$Mean_7dpi)]
# int_7dpi$Type_pair <- as.character(int_7dpi$Type_pair)
# int_7dpi$Partner_1 <- sapply(int_7dpi$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
# int_7dpi$Partner_2 <- sapply(int_7dpi$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})
# int_7dpi <- merge(int_7dpi, cell_type_freqs[Time == "7dpi", c("Cell_type", "Count")], by.x = "Partner_1", by.y = "Cell_type")
# colnames(int_7dpi)[7] <- "Sevendpi_count_1"
# int_7dpi <- merge(int_7dpi, cell_type_freqs[Time == "7dpi", c("Cell_type", "Count")], by.x = "Partner_2", by.y = "Cell_type")
# colnames(int_7dpi)[8] <- "Sevendpi_count_2"

subcol_7dpi_all <- colnames(inj_7dpi_all)[c(2, 12:ncol(inj_7dpi_all))]
int_7dpi_all <- melt(inj_7dpi_all[, ..subcol_7dpi_all],
                     id.vars = c("interacting_pair"),
                     measure.vars = subcol_7dpi_all[-1],
                     variable.name = "Type_pair", value.name = "Mean_7dpi")
int_7dpi_all <- int_7dpi_all[!is.na(int_7dpi_all$Mean_7dpi)]
int_7dpi_all$Type_pair <- as.character(int_7dpi_all$Type_pair)
int_7dpi_all$Partner_1 <- sapply(int_7dpi_all$Type_pair, function(x){unlist(strsplit(x, "\\|"))[1]})
int_7dpi_all$Partner_2 <- sapply(int_7dpi_all$Type_pair, function(x){unlist(strsplit(x, "\\|"))[2]})

ints <-
  merge(
    merge(ctrl_all_int[, c(2, 1, 4, 6)], int_3dpi_all[, c(2, 1, 4, 6)], 
          by = c("Partner_1", "Partner_2", "interacting_pair"), all = T),
    int_7dpi_all[, c(2, 1, 4, 6)], 
    by = c("Partner_1", "Partner_2", "interacting_pair"), all = T)
  # merge(
  #   merge(ctrl_int[, c(2, 1, 3, 4, 6)], int_3dpi[, c(2, 1, 3, 4, 6)], 
  #         by = c("Partner_1", "Partner_2", "interacting_pair", "id_cp_interaction"), all = T),
  #         int_7dpi[, c(2, 1, 3, 4, 6)], by = c("Partner_1", "Partner_2", "interacting_pair", "id_cp_interaction"), all = T)
# Extract the interactions that have NA in one of the three value columns, merge with the melted full interactions to get the
# (non-significant) value. Also mark that the values are non-significant.
# ints$Sig_ctrl <- ifelse(is.na(ints$Ctrl_mean), F, T)
# ints$Sig_3dpi <- ifelse(is.na(ints$Mean_3dpi), F, T)
# ints$Sig_7dpi <- ifelse(is.na(ints$Mean_7dpi), F, T)
# ints <- ints[, -(4:6)]

# ints <- merge(ints, ctrl_all_int[, -2], all.x = T)
# ints <- merge(ints, int_3dpi_all[, -2], all.x = T)
# ints <- merge(ints, int_7dpi_all[, -2], all.x = T)

ints[is.na(ints)] <- 0
# ggplot(ints) +
#   geom_point(aes(x = Ctrl_mean, y = Mean_3dpi))
# Draw a straight line, slope ~1, y-intercept 0.75, pick all interactions above that line. Enriched for?

# ggplot(ints) +
#   geom_point(aes(x = Ctrl_mean, y = Mean_7dpi))
# Why do we see repeating structures? If legitimate, probably caused by interactions between cell type group A and 
# cell types a, b, c.

# Overview of interacting pairs
ints_nod <- ints[!(ints$Partner_1 %in% c("Ery Duplex", "Immune Duplex")) & !(ints$Partner_2 %in% c("Ery Duplex", "Immune Duplex")), ]
collagen_like_interaction_pairs <- unique(ints_nod$interacting_pair[grepl("CAM|COL|FN1", ints_nod$interacting_pair)])
ints_nod_nocoll <- ints_nod[!(ints_nod$interacting_pair %in% collagen_like_interaction_pairs), ]

ints_r_nod_nocoll <- ints_nod_nocoll
colnames(ints_r_nod_nocoll)[1:2] <- c("Partner_2", "Partner_1")
ints_full_nod_nocoll <- rbind(ints_nod_nocoll, ints_r_nod_nocoll)
ints_full_nod_nocoll <- 
  ints_full_nod_nocoll[!(ints_full_nod_nocoll$Partner_1 %in% c("Ery Duplex", "Immune Duplex")) & 
                         !(ints_full_nod_nocoll$Partner_2 %in% c("Ery Duplex", "Immune Duplex")), ]
# collagen_like_interaction_pairs <- unique(ints_full$interacting_pair[grepl("CAM|COL|FN1", ints_full$interacting_pair)])
# ints_full_nocoll <- ints_full[!(ints_full$interacting_pair %in% collagen_like_interaction_pairs), ]

# int_sum <- aggregate(ints_full$Mean_7dpi, #Ctrl_mean, Mean_3dpi, Mean_7dpi,
#                      by = list(Partner_1 = ints_full$Partner_1, Partner_2 = ints_full$Partner_2),
#                      sum)
# int_sum_all_m <- data.frame(acast(int_sum, Partner_1 ~ Partner_2, value.var = "x"))
# colnames(int_sum_all_m) <- rownames(int_sum_all_m)
# int_sum_all_m[is.na(int_sum_all_m)] <- 0
int_sum_nod_nocoll <- aggregate(ints_full_nod_nocoll$Mean_3dpi, #Ctrl_mean, Mean_3dpi, Mean_7dpi,
                     by = list(Partner_1 = ints_full_nod_nocoll$Partner_1, Partner_2 = ints_full_nod_nocoll$Partner_2),
                     sum)
int_sum_nod_nocoll_all_m <- data.frame(acast(int_sum_nod_nocoll, Partner_1 ~ Partner_2, value.var = "x"))
colnames(int_sum_nod_nocoll_all_m) <- rownames(int_sum_nod_nocoll_all_m)
int_sum_nod_nocoll_all_m[is.na(int_sum_nod_nocoll_all_m)] <- 0
# pheatmap::pheatmap(int_sum_all_m, treeheight_col = 0, treeheight_row = 0, breaks = 0:100)
# png("./Images/Interaction_sum_nocoll_7dpi.png",
#     width = 560, height = 560)
pheatmap::pheatmap(int_sum_nod_nocoll_all_m,
                   treeheight_col = 0, treeheight_row = 0,
                   fontsize = 16,
                   breaks = seq(0, 24, length.out = 101))
# dev.off()

ggplot(ints_nod_nocoll[ints_nod_nocoll$Partner_1 == "Fibroblast (nppc)" & ints_nod_nocoll$Partner_2 == "Fibroblast (nppc)", ]) +
  geom_point(aes(x = Ctrl_mean, y = Mean_3dpi))

# Interaction #1
int_1 <- data.frame(acast(ints_nod_nocoll[ints_nod_nocoll$interacting_pair == "MDK_LRP1", ],
               Partner_1 ~ Partner_2, value.var = "Mean_3dpi", fill = 0))
unique(ints_nod_nocoll$Partner_2[ints_nod_nocoll$interacting_pair == "MDK_LRP1"])
colnames(int_1) <- unique(ints_nod_nocoll$Partner_2[ints_nod_nocoll$interacting_pair == "MDK_LRP1"])
# png("./Images/MDK_LRP1_3dpi.png",
#     width = 560, height = 560)
pheatmap::pheatmap(int_1, main = "MDK_LRP1 at 3dpi",
                   treeheight_col = 0, treeheight_row = 0,
                   fontsize = 16)
# dev.off()

nod_nocoll_wcounts <- merge(ints_nod_nocoll, interactions[, c(1:4)], by.x = "interacting_pair", by.y = "Pair")
nod_nocoll_wcounts$FC_3dpi_ctrl <- log2(nod_nocoll_wcounts$Mean_3dpi + 0.1) - log2(nod_nocoll_wcounts$Ctrl_mean + 0.1)
ggplot(nod_nocoll_wcounts[nod_nocoll_wcounts$Partner_1 == "Fibroblast (nppc)", ]) +
  geom_point(aes(x = Freq_3dpi, y = FC_3dpi_ctrl))

ggplot(interactions[interactions$Freq_3dpi < 50, ]) +
  geom_histogram(aes(x = Freq_3dpi), binwidth = 1)

low_freq_ints <- as.character(interactions$Pair[interactions$Freq_3dpi <= 10])
int_lf <- ints_full_nod_nocoll[ints_full_nod_nocoll$interacting_pair %in% low_freq_ints, ]
int_lf_sum <- aggregate(int_lf$Mean_3dpi, #Ctrl_mean, Mean_3dpi, Mean_7dpi,
                                by = list(Partner_1 = int_lf$Partner_1, Partner_2 = int_lf$Partner_2),
                                sum)#function(x) sum(x > 0))
int_lf_sum_m <- data.frame(acast(int_lf_sum, Partner_1 ~ Partner_2, value.var = "x", fill = 0))
colnames(int_lf_sum_m) <- rownames(int_lf_sum_m)
pheatmap::pheatmap(int_lf_sum_m,# main = "MDK_LRP1 at 3dpi",
                   treeheight_col = 0, treeheight_row = 0,
                   fontsize = 16)
require(igraph)
int_lf_sum_graph <- graph_from_data_frame(int_lf_sum[int_lf_sum$x > 0, ], directed = F)
int_lf_sum_graph <- graph_from_adjacency_matrix(as.matrix(int_lf_sum_m), weighted = T, mode = "undirected")
int_lf_sum_graph_main <- decompose(int_lf_sum_graph)[[1]]
plot(int_lf_sum_graph_main)
E(int_lf_sum_graph)$width <- E(int_lf_sum_graph)$weight

adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), nc=10)

interaction_counts <- aggregate(ints_nod_nocoll$Mean_3dpi,
                                by = list(interacting_pair = ints_nod_nocoll$interacting_pair),
                                function(x) sum(x > 0))

# Do any of these cell types interact with anything uniquely at all?
# Strength of ligand/receptor - 0 and 3 from deconvolved 
# Number of interaction partners @ 3 - from ints_nod_nocoll, nonzero 3dpi column
# Difference between 0 and 3dpi - calculated from strengths

# Interacting pairs from ints_nod_nocoll (nonzero 3dpi column), nppc as partner 1 and 2.
# Look up minimum expression value in nppc fibroblasts per interacting pair.
nppc_fib_p1 <- ints_nod[ints_nod$Partner_1 == "Fibroblast (nppc)", ]

# Get dataframe with all 0 and 3 interactions, with relevant gene names.
nppc_fib_interactions_0 <- unique(ctrl_dec[, c("gene_name", "id_cp_interaction", "Fibroblast (nppc)")]) #[, c("gene_name", "Fibroblast (nppc)")])
colnames(nppc_fib_interactions_0)[3] <- "Mean_ctrl" 
nppc_fib_interactions_3 <- unique(dec_3dpi[, c("gene_name", "id_cp_interaction", "Fibroblast (nppc)")])
colnames(nppc_fib_interactions_3)[3] <- "Mean_3dpi"
nppc_fib_interactions <- merge(nppc_fib_interactions_0, nppc_fib_interactions_3, all = T)
nppc_fib_interactions[is.na(nppc_fib_interactions)] <- 0

nppc_fib_p1 <- ints_nod[ints_nod$Partner_1 == "Fibroblast (nppc)", ]
nppc_fib_p1$P1_mol <- sapply(nppc_fib_p1$interacting_pair,
                             function(x){unlist(strsplit(x, "_"))[1]})
nppc_fib_p1_3dpi <- aggregate(nppc_fib_p1$Mean_3dpi,
          by = list(gene_name = nppc_fib_p1$P1_mol),
          function(x) sum(x > 0))
colnames(nppc_fib_p1_3dpi)[2] <- "P1_count"
nppc_fib_p2 <- ints_nod[ints_nod$Partner_2 == "Fibroblast (nppc)", ]
nppc_fib_p2$P2_mol <- sapply(nppc_fib_p2$interacting_pair,
                             function(x){unlist(strsplit(x, "_"))[2]})
nppc_fib_p2_3dpi <- aggregate(nppc_fib_p2$Mean_3dpi,
                              by = list(gene_name = nppc_fib_p2$P2_mol),
                              function(x) sum(x > 0))
colnames(nppc_fib_p2_3dpi)[2] <- "P2_count"
nppc_fib_3dpi <- merge(nppc_fib_p1_3dpi, nppc_fib_p2_3dpi, all = T)
nppc_fib_3dpi[is.na(nppc_fib_3dpi)] <- 0
nppc_fib_3dpi$Partners_3dpi <- nppc_fib_3dpi$P1_count + nppc_fib_3dpi$P2_count


# Interactions of interesting cell types ####
cell_type <- "Fibroblast-like cells" # "CM Ventricle (ttn.2)" #"Perivascular cells"
# ctrl_all_int[grep(cell_type, ctrl_all_int$Type_pair), ]

ctrl_all_int[, ':=' (int_signal_ctrl = sum(Ctrl_mean)), by = interacting_pair]
ctrl_all_int$Rel_ctrl_strength <- ctrl_all_int$Ctrl_mean^2/ctrl_all_int$int_signal_ctrl
int_3dpi_all[, ':=' (int_signal_3dpi = sum(Mean_3dpi)), by = interacting_pair]
int_3dpi_all$Rel_3dpi_strength <- int_3dpi_all$Mean_3dpi^2/int_3dpi_all$int_signal_3dpi
int_7dpi_all[, ':=' (int_signal_7dpi = sum(Mean_7dpi)), by = interacting_pair]
int_7dpi_all$Rel_7dpi_strength <- int_7dpi_all$Mean_7dpi^2/int_7dpi_all$int_signal_7dpi

cell_spec_ints <- ctrl_all_int[ctrl_all_int$Partner_1 == cell_type | ctrl_all_int$Partner_2 == cell_type, ] #[grep(cell_type, ctrl_all_int$Type_pair), ]
cell_spec_ints <- 
  merge(cell_spec_ints, 
        int_3dpi_all[int_3dpi_all$Partner_1 == cell_type | int_3dpi_all$Partner_2 == cell_type, ], 
        all = T)
cell_spec_ints <- 
  merge(cell_spec_ints, 
        int_7dpi_all[int_7dpi_all$Partner_1 == cell_type | int_7dpi_all$Partner_2 == cell_type, ], all = T)

ctrl_sig <- ctrl_all_int[, .(int_signal_ctrl = sum(Ctrl_mean)), by = interacting_pair]

# OLD ####
# deconvoluted <- fread("./Scripts/CellPhoneDB/out/H5_test/deconvoluted.txt")
# means <- fread("./Scripts/CellPhoneDB/out/H5_test/means.txt")
# p_values <- fread("./Scripts/CellPhoneDB/out/H5_test/pvalues.txt")
# significant_means <- fread("./Scripts/CellPhoneDB/out/H5_test/significant_means.txt")
# 
# # Find receptor-ligand networks (receptors and ligands that are found to interact together)
# interactions <- data.frame(Interaction = as.character(significant_means$interacting_pair),
#                            stringsAsFactors = F)
# interactions$A <- sapply(interactions$Interaction,
#                          function(x){
#                            unlist(strsplit(x, "_"))[1]
#                          })
# interactions$B <- sapply(interactions$Interaction,
#                          function(x){
#                            unlist(strsplit(x, "_"))[2]
#                          })
# interaction_graph <- graph_from_data_frame(interactions[, c("A", "B")], directed = T)
# interaction_components <- decompose(interaction_graph, mode = "weak")
# plot(interaction_components[[9]])
# 
# # Maybe take components 9 as example? It has 7 genes that interact, not too many, not too few.
# # Later we can also do discovery the other way around: what interacts with our cell type of interest?
# comp_ints <- interactions$Interaction[interactions$A %in% names(V(interaction_components[[9]]))]
# comp_means <- means[means$interacting_pair %in% comp_ints, ]
# comp_mm <- melt(comp_means[, c(2, 12:ncol(comp_means))])
# #data.frame(means[means$interacting_pair %in% comp_ints, ])
# type_interactions <- colnames(comp_means)[c(rep(T, 11), colSums(comp_means[, -c(1:11)]) > 0)]
# comp_means <- comp_means[, type_interactions]
# 
# # Number of interactions per cell type ####
# cell_type_int_freq <- data.frame(Cell_types = colnames(significant_means)[-(1:12)],
#                                  Frequency = apply(significant_means[, -(1:12)], 2, function(x) {sum(!is.na(x))}),
#                                  stringsAsFactors = F)
# cell_type_int_freq$Partner_1 <-
#   sapply(cell_type_int_freq$Cell_types,
#          function(x){unlist(strsplit(x, "\\|"))[1]})
# cell_type_int_freq$Partner_2 <-
#   sapply(cell_type_int_freq$Cell_types,
#          function(x){unlist(strsplit(x, "\\|"))[2]})
# ggplot(cell_type_int_freq) +
#   geom_tile(aes(x = Partner_1, y = Partner_2, fill = Frequency)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# # Sum this up over both partners (invert partners and if inversion != original, add Frequency).
# # Different scale; grey to red?
# 
# # x <- "Blood vessel endothelia (apln)|Blood vessel endothelia (apln)"
# # strsplit(x, "\\|")
# # Number of connections per interaction type ####
# interaction_connections <-
#   data.frame(Interaction = significant_means$interacting_pair,
#              Frequency = apply(significant_means[, -(1:12)], 1, function(x) {sum(!is.na(x))}))
