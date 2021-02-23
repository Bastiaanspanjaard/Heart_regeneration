# Requirements ####
require(data.table)
require(ggplot2)
require(pheatmap)
source("./Scripts/HR_library.R")

make_incidence_matrix <- function(sigmeans){
  incidences_long <- data.table(Interacting_pair = colnames(sigmeans)[-c(1:12)],
                                Interactions = apply(sigmeans[, -c(1:12)], 2, function(x) sum(!is.na(x))))
  incidences_long$From <-
    sapply(incidences_long$Interacting_pair,
           function(x) unlist(strsplit(x, "|", fixed = T))[1])
  incidences_long$To <-
    sapply(incidences_long$Interacting_pair,
           function(x) unlist(strsplit(x, "|", fixed = T))[2])
  incidence_dt <- dcast(incidences_long, From ~ To, value.var = "Interactions")
  incidence_matrix <- as.matrix(incidence_dt[, -1])
  rownames(incidence_matrix) <- incidence_dt$From
  
  return(incidence_matrix)
}

make_sigmeans_numeric <- function(sigmeans){
  for(col in names(sigmeans)[-(1:12)]){
    set(sigmeans, j = col, value = as.numeric(sigmeans[[col]]))
  }
  return(sigmeans)
}

# Note that the below is a different dataframe than the 'libraries' in HR_library.R - that contains
# only samples with lineage tracing, whereas this contains all samples with scRNA-seq.
sc_libraries <- fread("./Data/HR_setnames.csv")[, -1]

# Load data ####
ctrl <- fread("./Data/Control_alliance_zoom/significant_means.txt")
ctrl <- make_sigmeans_numeric(ctrl)
noinhib_3dpi <- fread("./Data/3dpi_noinhib_alliance_zoom/significant_means.txt")
noinhib_3dpi <- make_sigmeans_numeric(noinhib_3dpi)
inhib_3dpi <- fread("./Data/3dpi_inhib_alliance_zoom/significant_means.txt")
inhib_3dpi <- make_sigmeans_numeric(inhib_3dpi)
noinhib_7dpi <- fread("./Data/7dpi_noinhib_alliance_zoom/significant_means.txt")
noinhib_7dpi <- make_sigmeans_numeric(noinhib_7dpi)
inhib_7dpi <- fread("./Data/7dpi_inhib_alliance_zoom/significant_means.txt")
inhib_7dpi <- make_sigmeans_numeric(inhib_7dpi)

celltype_colors_in <- fread("./Data/Cell_type_colors_2.csv")
celltype_colors <- setNames(celltype_colors_in$color, celltype_colors_in$Cell.type)

cells <- fread("./Data/final_metadata_Tmacromerged_2.csv")
cells <- merge(cells[, c("Cell_name", "Cell_type", "orig.ident")], sc_libraries, by.x = "orig.ident", by.y = "heart")

count_cutoff <- 50
proportion_cutoff <- 0.01

ctrl_celltype_counts <- cells[cells$dpi == 0, .(Count = .N), by = Cell_type]
ctrl_celltype_counts$Proportion <- ctrl_celltype_counts$Count/sum(ctrl_celltype_counts$Count)
ctrl_celltype_counts <- 
  ctrl_celltype_counts[ctrl_celltype_counts$Count > count_cutoff | ctrl_celltype_counts$Proportion > proportion_cutoff, ]

noinhib_3dpi_celltype_counts <- cells[cells$dpi == 3 & cells$inhib %in% c("no", "DMSO"), .(Count = .N), by = Cell_type]
noinhib_3dpi_celltype_counts$Proportion <- noinhib_3dpi_celltype_counts$Count/sum(noinhib_3dpi_celltype_counts$Count)
noinhib_3dpi_celltype_counts <- 
  noinhib_3dpi_celltype_counts[noinhib_3dpi_celltype_counts$Count > count_cutoff | noinhib_3dpi_celltype_counts$Proportion > proportion_cutoff, ]

inhib_3dpi_celltype_counts <- cells[cells$dpi == 3 & cells$inhib %in% c("IWR1"), .(Count = .N), by = Cell_type]
inhib_3dpi_celltype_counts$Proportion <- inhib_3dpi_celltype_counts$Count/sum(inhib_3dpi_celltype_counts$Count)
inhib_3dpi_celltype_counts <- 
  inhib_3dpi_celltype_counts[inhib_3dpi_celltype_counts$Count > count_cutoff | inhib_3dpi_celltype_counts$Proportion > proportion_cutoff, ]

noinhib_7dpi_celltype_counts <- cells[cells$dpi == 7 & cells$inhib %in% c("no", "DMSO"), .(Count = .N), by = Cell_type]
noinhib_7dpi_celltype_counts$Proportion <- noinhib_7dpi_celltype_counts$Count/sum(noinhib_7dpi_celltype_counts$Count)
noinhib_7dpi_celltype_counts <- 
  noinhib_7dpi_celltype_counts[noinhib_7dpi_celltype_counts$Count > count_cutoff | noinhib_7dpi_celltype_counts$Proportion > proportion_cutoff, ]

inhib_7dpi_celltype_counts <- cells[cells$dpi == 7 & cells$inhib %in% c("IWR1"), .(Count = .N), by = Cell_type]
inhib_7dpi_celltype_counts$Proportion <- inhib_7dpi_celltype_counts$Count/sum(inhib_7dpi_celltype_counts$Count)
inhib_7dpi_celltype_counts <- 
  inhib_7dpi_celltype_counts[inhib_7dpi_celltype_counts$Count > count_cutoff | inhib_7dpi_celltype_counts$Proportion > proportion_cutoff, ]

# Make incidence matrix ####
incidence_ctrl <- make_incidence_matrix(ctrl)
incidence_ctrl <- incidence_ctrl[ctrl_celltype_counts$Cell_type, ctrl_celltype_counts$Cell_type]
incidence_inhib_3dpi <- make_incidence_matrix(inhib_3dpi)
incidence_inhib_3dpi <- incidence_inhib_3dpi[inhib_3dpi_celltype_counts$Cell_type, inhib_3dpi_celltype_counts$Cell_type]
incidence_noinhib_3dpi <- make_incidence_matrix(noinhib_3dpi)
incidence_noinhib_3dpi <- incidence_noinhib_3dpi[noinhib_3dpi_celltype_counts$Cell_type, noinhib_3dpi_celltype_counts$Cell_type]
incidence_inhib_7dpi <- make_incidence_matrix(inhib_7dpi)
incidence_inhib_7dpi <- incidence_inhib_7dpi[inhib_7dpi_celltype_counts$Cell_type, inhib_7dpi_celltype_counts$Cell_type]
incidence_noinhib_7dpi <- make_incidence_matrix(noinhib_7dpi)
incidence_noinhib_7dpi <- incidence_noinhib_7dpi[noinhib_7dpi_celltype_counts$Cell_type, noinhib_7dpi_celltype_counts$Cell_type]

# Count outgoing interactions ####
outgoing_ctrl <- data.table(From = rownames(incidence_ctrl),
                            Count = rowMeans(incidence_ctrl))
outgoing_ctrl$From <- factor(outgoing_ctrl$From, 
                             levels = outgoing_ctrl$From[order(outgoing_ctrl$Count, decreasing = T)])
outgoing_noinhib_3dpi <- data.table(From = rownames(incidence_noinhib_3dpi),
                                    Count = rowMeans(incidence_noinhib_3dpi))
outgoing_noinhib_3dpi$From <- factor(outgoing_noinhib_3dpi$From, 
                                     levels = outgoing_noinhib_3dpi$From[order(outgoing_noinhib_3dpi$Count, decreasing = T)])
outgoing_noinhib_7dpi <- data.table(From = rownames(incidence_noinhib_7dpi),
                                    Count = rowMeans(incidence_noinhib_7dpi))
outgoing_noinhib_7dpi$From <- factor(outgoing_noinhib_7dpi$From, 
                                     levels = outgoing_noinhib_7dpi$From[order(outgoing_noinhib_7dpi$Count, decreasing = T)])
outgoing_inhib_3dpi <- data.table(From = rownames(incidence_inhib_3dpi),
                                  Count = rowMeans(incidence_inhib_3dpi))
outgoing_inhib_3dpi$From <- factor(outgoing_inhib_3dpi$From, 
                                   levels = outgoing_inhib_3dpi$From[order(outgoing_inhib_3dpi$Count, decreasing = T)])
outgoing_inhib_7dpi <- data.table(From = rownames(incidence_inhib_7dpi),
                                  Count = rowMeans(incidence_inhib_7dpi))
outgoing_inhib_7dpi$From <- factor(outgoing_inhib_7dpi$From, 
                                   levels = outgoing_inhib_7dpi$From[order(outgoing_inhib_7dpi$Count, decreasing = T)])

max_out_count <- max(outgoing_ctrl$Count, outgoing_inhib_3dpi$Count, outgoing_inhib_7dpi$Count,
                 outgoing_noinhib_3dpi$Count, outgoing_noinhib_7dpi$Count)
max_out_plot <- 10 * ceiling(max_out_count/10)

# pdf("./Images/Ctrl_alliance_CPDB_outgoing_counts.pdf")
ggplot(outgoing_ctrl) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_out_plot)) +
  labs(title = paste("Outgoing interactions control")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Noinhib_3dpi_alliance_CPDB_outgoing_counts.pdf")
ggplot(outgoing_noinhib_3dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_out_plot)) +
  labs(title = paste("Outgoing interactions 3dpi")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Noinhib_7dpi_alliance_CPDB_outgoing_counts.pdf")
ggplot(outgoing_noinhib_7dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_out_plot)) +
  labs(title = paste("Outgoing interactions 7dpi")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Inhib_3dpi_alliance_CPDB_outgoing_counts.pdf")
ggplot(outgoing_inhib_3dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_out_plot)) +
  labs(title = paste("Outgoing interactions 3dpi wnt-inhibited")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Inhib_7dpi_alliance_CPDB_outgoing_counts.pdf")
ggplot(outgoing_inhib_7dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_out_plot)) +
  labs(title = paste("Outgoing interactions 7dpi wnt-inhibited")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# Count incoming interactions ####
incoming_ctrl <- data.table(From = colnames(incidence_ctrl),
                            Count = colMeans(incidence_ctrl))
incoming_ctrl$From <- factor(incoming_ctrl$From, 
                             levels = incoming_ctrl$From[order(incoming_ctrl$Count, decreasing = T)])
incoming_noinhib_3dpi <- data.table(From = colnames(incidence_noinhib_3dpi),
                                    Count = colMeans(incidence_noinhib_3dpi))
incoming_noinhib_3dpi$From <- factor(incoming_noinhib_3dpi$From, 
                                     levels = incoming_noinhib_3dpi$From[order(incoming_noinhib_3dpi$Count, decreasing = T)])
incoming_noinhib_7dpi <- data.table(From = colnames(incidence_noinhib_7dpi),
                                    Count = colMeans(incidence_noinhib_7dpi))
incoming_noinhib_7dpi$From <- factor(incoming_noinhib_7dpi$From, 
                                     levels = incoming_noinhib_7dpi$From[order(incoming_noinhib_7dpi$Count, decreasing = T)])
incoming_inhib_3dpi <- data.table(From = colnames(incidence_inhib_3dpi),
                                  Count = colMeans(incidence_inhib_3dpi))
incoming_inhib_3dpi$From <- factor(incoming_inhib_3dpi$From, 
                                   levels = incoming_inhib_3dpi$From[order(incoming_inhib_3dpi$Count, decreasing = T)])
incoming_inhib_7dpi <- data.table(From = colnames(incidence_inhib_7dpi),
                                  Count = colMeans(incidence_inhib_7dpi))
incoming_inhib_7dpi$From <- factor(incoming_inhib_7dpi$From, 
                                   levels = incoming_inhib_7dpi$From[order(incoming_inhib_7dpi$Count, decreasing = T)])

max_incoming_count <- max(incoming_ctrl$Count, incoming_inhib_3dpi$Count, incoming_inhib_7dpi$Count,
                     incoming_noinhib_3dpi$Count, incoming_noinhib_7dpi$Count)
max_incoming_plot <- 10 * ceiling(max_incoming_count/10)

# pdf("./Images/Ctrl_alliance_CPDB_incoming_counts.pdf")
ggplot(incoming_ctrl) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_incoming_plot)) +
  labs(title = paste("Incoming interactions control")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Noinhib_3dpi_alliance_CPDB_incoming_counts.pdf")
ggplot(incoming_noinhib_3dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_incoming_plot)) +
  labs(title = paste("Incoming interactions 3dpi")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Noinhib_7dpi_alliance_CPDB_incoming_counts.pdf")
ggplot(incoming_noinhib_7dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_incoming_plot)) +
  labs(title = paste("Incoming interactions 7dpi")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Inhib_3dpi_alliance_CPDB_incoming_counts.pdf")
ggplot(incoming_inhib_3dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_incoming_plot)) +
  labs(title = paste("Incoming interactions 3dpi wnt-inhibited")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# pdf("./Images/Inhib_7dpi_alliance_CPDB_incoming_counts.pdf")
ggplot(incoming_inhib_7dpi) +
  geom_bar(aes(x = From, y = Count, fill = From), stat = "identity") +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(limits = c(0, max_incoming_plot)) +
  labs(title = paste("Incoming interactions 7dpi wnt-inhibited")) +
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 10),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24)) 
# dev.off()

# Interaction heatmaps ####
pheatmap_break_max <- max(incidence_ctrl, incidence_noinhib_3dpi, incidence_noinhib_7dpi,
                          incidence_inhib_3dpi, incidence_inhib_7dpi)
# 101 breaks from 0 to max value:
incidence_breaks <- pheatmap_break_max * (0:100)/100

# pdf("./Images/Ctrl_alliance_CPDB_interaction_counts.pdf")
pheatmap(incidence_ctrl, clustering_method = "ward.D2", breaks = incidence_breaks,
         main = "Interaction counts control")
# dev.off()

# pdf("./Images/Noinhib_3dpi_alliance_CPDB_interaction_counts.pdf")
pheatmap(incidence_noinhib_3dpi, clustering_method = "ward.D2", breaks = incidence_breaks,
         main = "Interaction counts 3dpi")
# dev.off()

# pdf("./Images/Noinhib_7dpi_alliance_CPDB_interaction_counts.pdf")
pheatmap(incidence_noinhib_7dpi, clustering_method = "ward.D2", breaks = incidence_breaks,
         main = "Interaction counts 7dpi")
# dev.off()

# pdf("./Images/Inhib_3dpi_alliance_CPDB_interaction_counts.pdf")
pheatmap(incidence_inhib_3dpi, clustering_method = "ward.D2", breaks = incidence_breaks,
         main = "Interaction counts 3dpi wnt-inhibited")
# dev.off()

# pdf("./Images/Inhib_7dpi_alliance_CPDB_interaction_counts.pdf")
pheatmap(incidence_inhib_7dpi, clustering_method = "ward.D2", breaks = incidence_breaks,
         main = "Interaction counts 7dpi wnt-inhibited")
# dev.off()

colnames(noinhib_3dpi[, 1:13])
plvapb_columns <- c(1:12, grep("Bl\\.ves\\.EC [(]plvapb[)]$", colnames(noinhib_3dpi)))
plvapb_interactions <- 
  noinhib_3dpi[, ..plvapb_columns]
plvapb_interactions <-
  plvapb_interactions[apply(plvapb_interactions[, -(1:12)], 1, function(x) sum(!is.na(x)) > 0), ]

View(t(noinhib_3dpi[noinhib_3dpi$interacting_pair %in% c("VEGFA_FLT1", "VEGFA_KDR"), ]))

# OLD ####
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

# Count interactions between cell types ####
# Count number of significant interactions between two cell types (ctrl, inj_3pi, inj_7dpi)
count_interactions <- function(sig_means){
  # Returns matrix of source (row) - target (column) interaction counts.
  interaction_freq <- data.table(Cell_types = colnames(sig_means)[-(1:12)],
                         Frequency = apply(sig_means[, -(1:12)], 2, function(x) {sum(!is.na(x))}))
  interaction_freq$Partner_1 <-
    sapply(interaction_freq$Cell_types,
           function(x){unlist(strsplit(x, "\\|"))[1]})
  interaction_freq$Partner_2 <-
    sapply(interaction_freq$Cell_types,
           function(x){unlist(strsplit(x, "\\|"))[2]})
  interaction_freq_matrix <- acast(interaction_freq, Partner_1 ~ Partner_2, value.var = "Frequency")
  # interaction_freq_matrix <- interaction_freq_matrix + t(interaction_freq_matrix)
  # diag(interaction_freq_matrix) <- diag(interaction_freq_matrix)/2
  return(interaction_freq_matrix)
}

interaction_counts_ctrl <- count_interactions(ctrl)
pheatmap::pheatmap(interaction_counts_ctrl, clustering_method = "ward.D2")
interaction_counts_3dpi <- count_interactions(inj_3dpi)
pheatmap::pheatmap(interaction_counts_3dpi, clustering_method = "ward.D2")
interaction_counts_7dpi <- count_interactions(inj_7dpi)
pheatmap::pheatmap(interaction_counts_7dpi, clustering_method = "ward.D2")

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

# OLDER ####
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
