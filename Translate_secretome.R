# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("biomaRt")
library("biomaRt")
library(data.table)
library(rjson)

# Alliance orthology ####
orthology <- fromJSON(file="~/Documents/Projects/heart_Bo/Data/ORTHOLOGY-ALLIANCE-JSON_COMBINED_37.json")
orthology_DT <- data.table(DR_name = unlist(lapply(orthology$data, function(x) {x$Gene1Symbol})),
                          HS_name = unlist(lapply(orthology$data, function(x) {x$Gene2Symbol})),
                          DR_ID = unlist(lapply(orthology$data, function(x) {x$Gene1ID})),
                          HS_ID = unlist(lapply(orthology$data, function(x) {x$Gene2ID})),
                          Species_1 = unlist(lapply(orthology$data, function(x) {x$Gene1SpeciesName})),
                          Species_2 = unlist(lapply(orthology$data, function(x) {x$Gene2SpeciesName})),
                          AlgorithmsMatch = unlist(lapply(orthology$data, function(x) {x$AlgorithmsMatch})),
                          OutOfAlgorithms = unlist(lapply(orthology$data, function(x) {x$OutOfAlgorithms})),
                          IsBestScore = unlist(lapply(orthology$data, function(x) {x$IsBestScore})),
                          IsBestRevScore = unlist(lapply(orthology$data, function(x) {x$IsBestRevScore})))
DR_to_HS <- orthology_DT[Species_1 == "Danio rerio" & Species_2 == "Homo sapiens" & IsBestScore == "Yes", ]
DR_to_HS_unique <- DR_to_HS[!duplicated(DR_to_HS$DR_name), ]
View(DR_to_HS[duplicated(DR_to_HS$DR_name), ])
# Orthology has, for zebrafish, gene IDs from zfin, and for human, gene IDs from Hugo, I think?
# Allow multiple orthologues in zebrafish, but not multiple orthologues in human - it would be unclear what to do
# with counts on one zebrafish gene that maps to multiple human genes.
# write.table(DR_to_HS_unique, "~/Documents/Projects/heart_Bo/Data/Alliance_gene_name_conversion.scsv",
#             quote = F, row.names = F, sep = ";")

HS <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
HS_attributes <- listAttributes(HS)
HS_hgnc_IDs <-
  getBM(mart = HS, attributes=c("ensembl_gene_id", "external_gene_name", "hgnc_id", "chromosome_name"))
# Note there are genes with alternative IDs (but same names) that are on patches - remove these.
HS_chromosomes <- c(1:22, "X", "Y", "MT")
HS_hgnc_IDs_chromosomes <- HS_hgnc_IDs[HS_hgnc_IDs$chromosome_name %in% HS_chromosomes &
                                         HS_hgnc_IDs$hgnc_id != "", ]
# This leaves a few genes with duplicated names and hugo IDs - pick one.
# HS_hgnc_IDs_chromosomes <- HS_hgnc_IDs_chromosomes[!duplicated(HS_hgnc_IDs_chromosomes[, c("external_gene_name", "hgnc_id")]), ]

DR <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="drerio_gene_ensembl")
DR_attributes <- listAttributes(DR)
DR_zfin_IDs <-
  getBM(mart = DR, attributes=c("ensembl_gene_id", "external_gene_name", "zfin_id_id", "chromosome_name"))
DR_zfin_IDs$zfin_id <- paste("ZFIN:", DR_zfin_IDs$zfin_id_id, sep = "")
DR_chromosomes <- c(1:25, "MT")
DR_zfin_IDs_chromosomes <- DR_zfin_IDs[DR_zfin_IDs$chromosome_name %in% DR_chromosomes &
                                         DR_zfin_IDs$zfin_id_id != "", ]
# This leaves a few genes with duplicated names and hugo IDs - pick one.
# DR_zfin_IDs_chromosomes <- DR_zfin_IDs_chromosomes[!duplicated(DR_zfin_IDs_chromosomes[, c("external_gene_name", "zfin_id")]), ]

DR_to_HS_Ensembl <-
  merge(DR_to_HS_unique[, c("DR_name", "HS_name", "DR_ID", "HS_ID")], 
        HS_hgnc_IDs_chromosomes[, c("hgnc_id", "ensembl_gene_id")], by.x = "HS_ID", by.y = "hgnc_id")
colnames(DR_to_HS_Ensembl)[5] <- "Human.gene.stable.ID"
DR_to_HS_Ensembl <-
  merge(DR_to_HS_Ensembl,
        DR_zfin_IDs_chromosomes[, c("zfin_id", "ensembl_gene_id")], by.x = "DR_ID", by.y = "zfin_id")
colnames(DR_to_HS_Ensembl)[c(3, 4, 6)] <- c("Gene.name", "Human.gene.name", "Gene.stable.ID")
# write.table(DR_to_HS_Ensembl, "~/Documents/Projects/heart_Bo/Data/Alliance_gene_name_conversion.scsv",
#             quote = F, row.names = F, sep = ";")

# Secretome ####
human_secretome <- fread("~/Documents/Projects/heart_Bo/Data/Human_secretome.tsv")
human_secretome_DR <- merge(human_secretome, DR_to_HS_unique, by.x = "Gene", by.y = "HS_name", all.x = T)
# The names in the secretome are outdated. We need the conversion between Hugo IDs and Ensembl IDs to match.

# human_secretome_hgnc <- merge(human_secretome, HS_hgnc_IDs, by.x = "Ensembl", by.y = "ensembl_gene_id")
human_secretome_hgnc <- merge(human_secretome, HS_hgnc_IDs_chromosomes, by.x = "Ensembl", by.y = "ensembl_gene_id")
human_secretome_hgnc <- human_secretome_hgnc[human_secretome_hgnc$hgnc_id!="", ]
DR_secretome <- merge(human_secretome_hgnc, DR_to_HS_unique, by.x = "hgnc_id", by.y = "HS_ID", all.x = T)
DR_secretome <- DR_secretome[, c("Ensembl", "Gene", "Gene description", "Evidence", "DR_name", "DR_ID")]
DR_secretome <- DR_secretome[!is.na(DR_secretome$DR_name), ]
write.table(DR_secretome, "~/Documents/Projects/heart_Bo/Data/Alliance_secretome_gene_names_noDRduplicates.scsv",
            quote = F, row.names = F, sep = ";")

