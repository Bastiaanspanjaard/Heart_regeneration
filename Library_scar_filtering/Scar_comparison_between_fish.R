# Dependencies/parameters ####
source("~/Documents/Projects/TOMO_scar/Scripts/linnaeus-scripts/scar_helper_functions.R")
min.presence <- 2 # Cells that have to have a scar for it to be considered
# present in an organism

# Load data ####
Hr10.scars <- read.csv("./Data/scars/Hr10_used_scars.csv",
                       stringsAsFactors = F)

# Z1.scars <- read.csv("./Data/2017_10X_1/Z1_used_scars_7Larvae.csv",
#                      stringsAsFactors = F)
# Z2.scars <- read.csv("./Data/2017_10X_2/Z2_used_scars_7Larvae.csv",
#                      stringsAsFactors = F)
# Z3.scars <- read.csv("./Data/2017_10X_2/Z3_used_scars_7Larvae.csv",
#                      stringsAsFactors = F)
# Z4.scars <- read.csv("./Data/2017_10X_10_CR/Z4_used_scars_7Larvae.csv",
#                      stringsAsFactors = F)
# Z5.scars <- read.csv("./Data/2017_10X_10_CR/Z5_used_scars_7Larvae.csv",
#                      stringsAsFactors = F)
# A5.scars <- read.csv("./Data/2017_10X_7/A5_used_scars_3Adults.csv",
#                      stringsAsFactors = F)
# A6.scars <- read.csv("./Data/2017_10X_6/A6_used_scars_3Adults.csv",
#                      stringsAsFactors = F)
# A7.scars <- read.csv("./Data/2018_10X_1/A7_used_scars_3Adults.csv",
#                      stringsAsFactors = F)

# Count presence of scars ####
Hr10.us <- data.frame(table(Hr10.scars$Sequence))
colnames(Hr10.us) <- c("Sequence", "Freq.Hr10")

# Z1.us <- data.frame(table(Z1.scars$Sequence)) 
# colnames(Z1.us) <- c("Sequence", "Freq.Z1")
# Z2.us <- data.frame(table(Z2.scars$Sequence)) 
# colnames(Z2.us) <- c("Sequence", "Freq.Z2")
# Z3.us <- data.frame(table(Z3.scars$Sequence)) 
# colnames(Z3.us) <- c("Sequence", "Freq.Z3")
# Z4.us <- data.frame(table(Z4.scars$Sequence)) 
# colnames(Z4.us) <- c("Sequence", "Freq.Z4")
# Z5.us <- data.frame(table(Z5.scars$Sequence)) 
# colnames(Z5.us) <- c("Sequence", "Freq.Z5")
# A5.us <- data.frame(table(A5.scars$Sequence))
# colnames(A5.us) <- c("Sequence", "Freq.A5")
# A6.us <- data.frame(table(A6.scars$Sequence))
# colnames(A6.us) <- c("Sequence", "Freq.A6")
# A7.us <- data.frame(table(A7.scars$Sequence))
# colnames(A7.us) <- c("Sequence", "Freq.A7")

# Merge all
unique.scars <- Hr10.us
# unique.scars <-
#   merge(
#     merge(
#       merge(Z1.us, Z2.us, all = T),
#       merge(Z3.us, Z4.us, all = T), all = T),
#     merge(
#       merge(Z5.us, A5.us, all = T),
#       merge(A6.us, A7.us, all = T), all = T),
#     all = T)
unique.scars[is.na(unique.scars)] <- 0
unique.scars$Presence <- unique.scars$Freq.Hr10 >= min.presence
# unique.scars$Presence <- apply(unique.scars[, -1], 1,
#                                function(x) sum(x >= min.presence))
table(unique.scars$Presence)

all.CIGARs <- unique(Hr10.scars[, c("Sequence", "CIGAR")])
# all.CIGARs <- unique(rbind(Z1.scars[, c("Sequence", "CIGAR")], 
#                            Z2.scars[, c("Sequence", "CIGAR")],
#                            Z3.scars[, c("Sequence", "CIGAR")],
#                            Z4.scars[, c("Sequence", "CIGAR")],
#                            Z5.scars[, c("Sequence", "CIGAR")],
#                            A5.scars[, c("Sequence", "CIGAR")],
#                            A6.scars[, c("Sequence", "CIGAR")],
#                            A7.scars[, c("Sequence", "CIGAR")]))
all.CIGARs <- all.CIGARs[!duplicated(all.CIGARs$Sequence), ]
unique.scars <- merge(unique.scars, all.CIGARs)

# Compare presence with probabilities ####
unique.scars$Sequence.short <- 
  substr(unique.scars$Sequence, 3 + sc.primer.length, sc.primer.length + 53)

# Load scar probabilities. Because bulk scar sequencing is done differently than
# single-cell scar sequencing, some single cell scars cannot be assigned a
# probability because they cannot be observed in bulk sequencing. We filter 
# these out.
scar.probabilities <- read.csv("~/Documents/Projects/TOMO_scar/Data/scar_Probs.csv",
                               stringsAsFactors = F)
# scar.probabilities <- read.csv("./Data/scar_Probs.csv",
#                                stringsAsFactors = F)
unique.scars <- merge(unique.scars, scar.probabilities[, c("Sequence", "p", "Embryos")],
                      by.x = "Sequence.short", by.y = "Sequence", all.x = T)
unique.scars$Embryos[is.na(unique.scars$Embryos)] <- 0
unique.scars$p[is.na(unique.scars$p)] <- min(unique.scars$p, na.rm = T)/100
unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

# ggplot(unique.scars[unique.scars$Presence == 5, ]) +
#   geom_histogram(aes(x = p)) +
#   scale_x_log10(limits = c(-0.01, 0.2)) +
#   labs(title = "Scars shared between five embryos")
# ggplot(unique.scars[unique.scars$Presence == 4, ]) +
#   geom_histogram(aes(x = p)) +
#   scale_x_log10(limits = c(-0.01, 0.2)) +
#   labs(title = "Scars shared between four embryos")
# ggplot(unique.scars[unique.scars$Presence == 3, ]) +
#   geom_histogram(aes(x = p)) +
#   scale_x_log10(limits = c(-0.01, 0.2)) +
#   labs(title = "Scars shared between three embryos")
# ggplot(unique.scars[unique.scars$Presence == 2, ]) +
#   geom_histogram(aes(x = p)) +
#   scale_x_log10(limits = c(-0.01, 0.2)) +
#   labs(title = "Scars shared between two embryos")
# ggplot(unique.scars[unique.scars$Presence == 1, ]) +
#   geom_histogram(aes(x = p)) +
#   scale_x_log10(limits = c(-0.01, 0.2)) +
#   labs(title = "Scars unique to one embryo")

# ggplot(unique.scars[unique.scars$p < 0.01, ]) +
#   geom_histogram(aes(x = Presence))
# ggplot(unique.scars[unique.scars$p < 0.001, ]) +
#   geom_histogram(aes(x = Presence))
# sum(unique.scars$Presence > 1 & unique.scars$p < 0.01)


# Compute fraction of unscarred RFP ####
# sample.columnames <- 
#   paste("Freq", c("Z1", "Z2", "Z3", "Z4", "Z5", "A5", "A6", "A7"), sep = ".")
# RFP.counts <- data.frame(t(unique.scars[unique.scars$Scar == "0:75M", sample.columnames]))
# colnames(RFP.counts)[1] <- "RFP.count"
# RFP.counts$Total <- colSums(unique.scars[, sample.columnames])
# RFP.counts$Percentage <- 100 * RFP.counts$RFP.count/RFP.counts$Total

# Count total and unique scars per library ####
# unique.scars.nowt <- unique.scars[-1, ]
# scars.Z1 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z1 >= min.presence]
# scars.Z2 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z2 >= min.presence]
# scars.Z3 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z3 >= min.presence]
# scars.Z4 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z4 >= min.presence]
# scars.Z5 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z5 >= min.presence]
# scars.A5 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.A5 >= min.presence]
# scars.A6 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.A6 >= min.presence]
# scars.A7 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.A7 >= min.presence]
# 
# sum(scars.A7 %in% unique.scars.nowt$Scar[unique.scars.nowt$Presence == 1])

# Create Venn diagram ####
# require(VennDiagram)
# venn.fill <- RColorBrewer::brewer.pal(3, "Set1")
# venn.diagram(list(Z2 = scars.Z2, ZZ4 = scars.Z4, #Z1 = scars.Z1, 3 = scars.Z3, 
#                   Z5 = scars.Z5),
#              filename = "./Images/Z245_scar_venn.png",
#              fontfamily = "sans", main.fontfamily = "sans",
#              imagetype = "png", fill = venn.fill)

# Write output ####
Hr10.scars.compared <- merge(Hr10.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
write.csv(Hr10.scars.compared, file = "./Data/scars/Hr10_scars_compared.csv",
          quote = F, row.names = F)

# Z1.scars.compared <- 
#   merge(Z1.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# write.csv(Z1.scars.compared, file = "./Data/2017_10X_1/Z1_scars_compared.csv",
#           quote = F, row.names = F)
# sum(Z1.scars.compared$Presence > 0 & Z1.scars.compared$Presence < 3 & Z1.scars.compared$p < 0.001)
# Z2.scars.compared <- 
#   merge(Z2.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# sum(Z2.scars.compared$Presence > 0 & Z2.scars.compared$Presence < 2 & Z2.scars.compared$p < 0.001)
# write.csv(Z2.scars.compared, file = "./Data/2017_10X_2/Z2_scars_compared.csv",
#           quote = F, row.names = F)
# Z3.scars.compared <- 
#   merge(Z3.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# sum(Z3.scars.compared$Presence > 0 & Z3.scars.compared$Presence < 2 & Z3.scars.compared$p < 0.001)
# write.csv(Z3.scars.compared, file = "./Data/2017_10X_2/Z3_scars_compared.csv",
#           quote = F, row.names = F)
# Z4.scars.compared <- 
#   merge(Z4.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# sum(Z4.scars.compared$Presence > 0 & Z4.scars.compared$Presence < 2 & Z4.scars.compared$p < 0.001)
# write.csv(Z4.scars.compared, file = "./Data/2017_10X_10_CR/Z4_scars_compared.csv",
#           quote = F, row.names = F)
# Z5.scars.compared <- 
#   merge(Z5.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# sum(Z5.scars.compared$Presence > 0 & Z5.scars.compared$Presence < 2 & Z5.scars.compared$p < 0.001)
# write.csv(Z5.scars.compared, file = "./Data/2017_10X_10_CR/Z5_scars_compared.csv",
#           quote = F, row.names = F)
# A5.scars.compared <- 
#   merge(A5.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# write.csv(A5.scars.compared, file = "./Data/2017_10X_7/A5_scars_compared.csv",
#           quote = F, row.names = F)
# A6.scars.compared <- 
#   merge(A6.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# write.csv(A6.scars.compared, file = "./Data/2017_10X_6/A6_scars_compared.csv",
#           quote = F, row.names = F)
# A7.scars.compared <- 
#   merge(A7.scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
# sum(A7.scars.compared$Presence > 0 & A7.scars.compared$Presence < 2 & A7.scars.compared$p < 0.001)
# write.csv(A7.scars.compared, file = "./Data/2018_10X_1/A7_scars_compared.csv",
#           quote = F, row.names = F)
