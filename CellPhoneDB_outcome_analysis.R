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
