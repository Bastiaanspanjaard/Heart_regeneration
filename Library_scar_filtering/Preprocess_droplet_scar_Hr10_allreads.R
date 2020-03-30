# Description ####
# Pre-process scar output files: merge with succesful whole-transcriptome
# barcodes, remove missequenced scars per cell, validate outcome by hand.

# NB This version works for multiple libraries from the same organism.

# Parameters ####
log2.cutoff <- 3

# Dependencies ####
source("~/Documents/Projects/TOMO_scar/Scripts/linnaeus-scripts/scar_helper_functions.R")

# Load data ####
# Scars
scars.in <- read.csv("./Data/scars/Hr10_scar_filtered_scars.csv",
                       stringsAsFactors = F, sep = "\t")
scars.in$Library <- "Hr10"
# scars.in.2 <- read.csv("./Data/2017_10X_6/all_16/H6_scar_filtered_scars.csv",
#                      stringsAsFactors = F, sep = "\t")
# scars.in.2$Library <- "H6"
# scars.in.3 <- read.csv("./Data/2017_10X_6/all_16/P6_scar_filtered_scars.csv",
#                        stringsAsFactors = F, sep = "\t")
# scars.in.3$Library <- "P6"
# scars.in <- rbind(scars.in.1, scars.in.2, scars.in.3)
scars.in$Cell <- paste(scars.in$Library, scars.in$Barcode, sep = "_")

# Select only cells that exist in the mRNA data
wt.all.cells <- read.csv("./Data/final_metadata.csv",
                           stringsAsFactors = F, sep = ",")
colnames(wt.all.cells)[1] <- "Cell"
wt.all.cells$lineage.ident[grepl("Endocardium", wt.all.cells$lineage.ident)] <- 
  wt.all.cells$final.zoom[grepl("Endocardium", wt.all.cells$lineage.ident)]
ery.cells <- read.csv("./Data/final_erythrocytes.csv", stringsAsFactors = F, sep = ",")
colnames(ery.cells) <- c("Cell", "orig.ident", "lineage.ident")

wt.cells <- rbind(wt.all.cells[wt.all.cells$orig.ident == "Hr10", c("Cell", "orig.ident", "lineage.ident")],
                  ery.cells[ery.cells$orig.ident == "Hr10", ])
colnames(wt.cells)[2:3] <- c("Library", "Cell.type")

wt.cells$Barcode <-
  sapply(wt.cells$Cell,
         function(x){unlist(strsplit(x, "_"))[2]})

# wt.all.cells.2 <- read.csv("./Data/Adult_data/Adults567_Seurat_batch_r_out_heart_subset_types.csv",
#                            stringsAsFactors = F, sep = ",")
# wt.cells.2 <- wt.all.cells.2[wt.all.cells.2$Library == "H6", ]
# wt.all.cells.3 <- read.csv("./Data/Adult_data/Adults567_Seurat_batch_r_out_pancreas_subset_types.csv",
#                            stringsAsFactors = F, sep = ",")
# wt.cells.3 <- wt.all.cells.3[wt.all.cells.3$Library == "P6", ]
# wt.cells <- rbind(wt.cells.1, wt.cells.2, wt.cells.3)

scars.unfiltered <- merge(scars.in, wt.cells[, c("Cell", "Library", "Barcode", "Cell.type")])
scars.unfiltered$Scar.id <- 1:nrow(scars.unfiltered)
scars.unfiltered$Keep <- T
scars.unfiltered$Pair <- "With"
cells <- unique(scars.unfiltered$Cell)

# Read in all >1 read sequences (all UMIs)
all.scars.g1 <- read.table("./Data/scars/Hr10_scar_reads_over1.txt",
                             sep = "\t", stringsAsFactors = F)
colnames(all.scars.g1)[-1] <- c("Barcode", "UMI", "Location", "Sequence")
all.scars.g1$Library <- "Hr10"
# all.scars.g1.2 <- read.table("./Data/2017_10X_7/P5_scar_reads_over1.txt",
#                              sep = "\t", stringsAsFactors = F)
# colnames(all.scars.g1.2)[-1] <- c("Barcode", "UMI", "Location", "Sequence")
# all.scars.g1.2$Library <- "P5"
# all.scars.g1.3 <- read.table("./Data/2017_10X_7/B5_scar_reads_over1.txt",
#                              sep = "\t", stringsAsFactors = F)
# colnames(all.scars.g1.3)[-1] <- c("Barcode", "UMI", "Location", "Sequence")
# all.scars.g1.3$Library <- "B5"
# all.scars.g1 <- rbind(all.scars.g1.1, all.scars.g1.2, all.scars.g1.3)
all.scars.g1$Cell <- paste(all.scars.g1$Library, all.scars.g1$Barcode, sep = "_")

all.scars.g1$V11 <- trimws(all.scars.g1$V1)
all.scars.g1$Reads <- sapply(all.scars.g1$V11,
                             function(x) unlist(strsplit(x, " "))[1])
all.scars.g1$CIGAR <- sapply(all.scars.g1$V11,
                             function(x) unlist(strsplit(x, " "))[2])
all.barcode.UMIs <-
  data.frame(table(all.scars.g1$Sequence, all.scars.g1$Cell))
all.barcode.UMIs <- all.barcode.UMIs[all.barcode.UMIs$Freq > 0, ]
colnames(all.barcode.UMIs) <- c("Sequence", "Cell", "UMIs")

# Filter scars per cell ####
# Go through all cells one by one. Within each cell, calculate the Hamming
# distances between the sequences; if these are one or two and the read 
# difference is high (log2(read1/read2) > log2.cutoff), flag the scar for 
# removal.
for(c in 1:length(cells)){
  cell <- cells[c]
  cell.scars <- scars.unfiltered[scars.unfiltered$Cell == cell, ]
  if(nrow(cell.scars) < 2){next}
  cell.scars <- cell.scars[order(-cell.scars$Reads), ]
  n.scars <- nrow(cell.scars)
  cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")
  
  # Identify all sequences that have to be removed.
  for(k in 1:length(cell.seqdist)){
    if(cell.seqdist[k] >2){next}
    else{
      # Get scar indices. Calculate log2 ratio between two reads. If too high,
      # flag lowest scar (Keep = F).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
        scars.unfiltered$Keep[scars.unfiltered$Scar.id == cell.scars$Scar.id[low.scar]] <-
          F
      }
    }
  }
}

# Remove scars that have HD 1 to scars with much higher reads.
scars.filter.1 <- scars.unfiltered[scars.unfiltered$Keep, ]

# Go throught the barcodes and calculate the distances between the scars again. 
# If the HD is one and the read difference is low, mark the scars as suspect. 
# If the HD is two and the read difference is high, mark the scars as suspect.
for(c in 1:length(cells)){
  cell <- cells[c]
  cell.scars <- scars.filter.1[scars.filter.1$Cell == cell, ]
  if(nrow(cell.scars) < 2){next}
  cell.scars <- cell.scars[order(-cell.scars$Reads), ]
  n.scars <- nrow(cell.scars)
  cell.seqdist <- stringdistmatrix(cell.scars$Sequence, method = "hamming")

  for(k in 1:length(cell.seqdist)){
    if(cell.seqdist[k] >2){next}
    else if(cell.seqdist[k] == 1){
      # Get scar indices. Flag both (paste(Pair, other
      # scar)).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      # If the scars are legitimate scars that are known to be close in HD,
      # skip this pair.
      # if((cell.scars$Sequence[high.scar] %in% allowed.close.scars) &
      #    (cell.scars$Sequence[low.scar] %in% allowed.close.scars)){next}
      scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
        paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
              scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
              sep = ".")
      scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
        paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
              scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
              sep = ".")
      }
    else if(cell.seqdist[k] == 2){
      # Get scar indices. Calculate log2 ratio between two reads. If too high,
      # flag both (paste(Pair, other scar)).
      scar.indices <- get.dist.index(k, n.scars)
      high.scar <- min(scar.indices)
      low.scar <- max(scar.indices)
      if(log2(cell.scars$Reads[high.scar]/cell.scars$Reads[low.scar]) > log2.cutoff){
        scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]] <-
          paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                sep = ".ratio.")
        scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]] <-
          paste(scars.filter.1$Pair[scars.filter.1$Scar.id == cell.scars$Scar.id[high.scar]],
                scars.filter.1$Scar.id[scars.filter.1$Scar.id == cell.scars$Scar.id[low.scar]],
                sep = ".ratio.")
      }
    }
  }
}

# Assess and remove suspect scars ####
# First determine which scars may need to be filtered out - occurring more than
# once but less than their possible parent scar in the same cell.
scars.assess <- scars.filter.1[scars.filter.1$Pair != "With", ]
scars.filter.2 <- scars.filter.1[, 1:7]
scars.assess$Incidence <-
  sapply(scars.assess$Sequence,
         function(x) sum(scars.filter.1$Sequence == x))
scars.assess$Pair.incidence <-
  sapply(scars.assess$Pair,
         function(x){
           pairs <- unlist(strsplit(x, "[.]"))[-1]
           return(min(scars.assess$Incidence[scars.assess$Scar.id %in% pairs]))}
  )
scars.assess$Min.pair.incidence <-
  apply(scars.assess[, c("Incidence", "Pair.incidence")], 1, min)
scars.assess.2 <- scars.assess[scars.assess$Min.pair.incidence > 1, ]

# Determine whether scars also occur in cells without the parent sequence 
# (criterion 1)
seq.freq <- data.frame(table(scars.assess.2$Sequence))
colnames(seq.freq)[1] <- "Sequence"
scars.assess.2 <- merge(scars.assess.2, seq.freq)
scars.assess.2$Crit.1 <- 
  !(scars.assess.2$Min.pair.incidence == scars.assess.2$Freq)

# Determine whether scars occur with more than one UMI (criterion 2)
scars.assess.2 <- merge(scars.assess.2, all.barcode.UMIs)
scars.assess.2$Crit.2 <- (scars.assess.2$UMIs > 1)

# Determine whether scars have at least one UMI unrelated to the UMIs of the
# parent sequence (criterion 3).
# Select scars to test - the ones that have the lowest incidence in the full
# dataset.
scars.assess.2.low.inc <- 
  scars.assess.2[scars.assess.2$Incidence <= scars.assess.2$Pair.incidence, ]
# If criteria 1 and 2 are both false, the scars are removed automatically
scars.assess.2.out <- 
  scars.assess.2.low.inc[!scars.assess.2.low.inc$Crit.1 & 
                           !scars.assess.2.low.inc$Crit.2, ]
# Else if one of criteria 1 and 2 are true, the scars are maybe kept, if
# they have at least one UMI that is at least 2 HD away from its suspected
# parent scar UMIs.
scars.assess.2.maybe <- 
  scars.assess.2.low.inc[scars.assess.2.low.inc$Crit.1 | 
                           scars.assess.2.low.inc$Crit.2, ]

# Loop over all scars to calculate the number of UMIs they have that are more
# than 1 HD away from UMIs of their 'parent' scar
scars.assess.2.maybe$Crit.3 <- NA
for(c.s in 1:nrow(scars.assess.2.maybe)){
  print(paste(c.s, "out of", nrow(scars.assess.2.maybe)))
  c.sequence <- scars.assess.2.maybe$Sequence[c.s]
  c.cell <- scars.assess.2.maybe$Cell[c.s]
  c.scar <- scars.assess.2.maybe$Scar.id[c.s]
  c.parent.sequence <- 
    scars.assess$Sequence[grepl(paste(c.scar, "\\.", sep = ""), scars.assess$Pair) |
                              grepl(paste(c.scar, "$", sep = ""), scars.assess$Pair)]
  
  validate.seq <- 
    all.scars.g1[all.scars.g1$Cell == c.cell &
                   all.scars.g1$Sequence %in% c.sequence, ]
  validate.parent <- 
    all.scars.g1[all.scars.g1$Cell == c.cell &
                   all.scars.g1$Sequence %in% c.parent.sequence, ]
  validate.seq$Min.HD <-
    sapply(validate.seq$UMI,
           function(x) min(stringdistmatrix(x, validate.parent$UMI, 
                                            method = "hamming")))
  scars.assess.2.maybe$Crit.3[c.s] <- sum(validate.seq$Min.HD > 1)
}
scars.assess.2.maybe$Out <-ifelse(scars.assess.2.maybe$Crit.3 ==0 , T, F)

# Remove sequencing error scars from the dataset
sequencing.error.scars <-
  c(scars.assess.2.out$Scar.id, 
    scars.assess.2.maybe$Scar.id[scars.assess.2.maybe$Out])
scars.output <- scars.filter.1[!(scars.filter.1$Scar.id %in% sequencing.error.scars),
                               c(1:6, 8)]

# Write filtered results ####
# write.csv(scars.output, "./Data/2017_10X_1/Z1_preprocessed_scars_allreads_7Larvae.csv",
#           quote = F, row.names = F)

# Undersequencing filter ####
# We use two final filters for the scar data. The first is an undersequencing
# filter; this removes all scars whose read number is under a cutoff. This
# cutoff is determined based on the histogram of reads per observed scar. 

# Undersequencing filter
ggplot(scars.output[scars.output$Library == "Hr10", ]) +
  geom_histogram(aes(x = Reads), binwidth = 10) +
  scale_x_continuous(limits = c(-10, 3000))
min.scar.reads <- 500
scars.output.2 <- scars.output[scars.output$Reads >= min.scar.reads, ]
# ggplot(scars.output[scars.output$Library == "H6", ]) +
#   geom_histogram(aes(x = Reads), binwidth = 10) +
#   scale_x_continuous(limits = c(-10, 1000))
# min.scar.reads.H6 <- 250
# scars.output.2.H6 <- scars.output[scars.output$Reads >= min.scar.reads.H6 & 
#                                     scars.output$Library == "H6", ]
# ggplot(scars.output[scars.output$Library == "P6", ]) +
#   geom_histogram(aes(x = Reads), binwidth = 10) +
#   scale_x_continuous(limits = c(-10, 2000))
# min.scar.reads.P6 <- 500
# scars.output.2.P6 <- scars.output[scars.output$Reads >= min.scar.reads.P6 & 
#                                     scars.output$Library == "P6", ]
# scars.output.2 <- rbind(scars.output.2.B6, scars.output.2.H6, scars.output.2.P6)

# Doublet filter ####
# The second filter is a doublet/bleedthrough filter that removes cells that have
# more scars than a cutoff set by looking at all cells of that type. This
# cutoff is currently set the same for all cell types.

# Remove cells that have too many scars
cell.scar.count <-
  data.frame(table(scars.output.2$Cell))
colnames(cell.scar.count) <- c("Cell", "Scars")
cell.scar.count$Cell <- as.character(cell.scar.count$Cell)
cell.scar.count <- merge(cell.scar.count, wt.cells[, c("Cell", "Cell.type")])

maximum.scars <- data.frame(Cell.type = unique(cell.scar.count$Cell.type),
                            Maximum = NA)

for(c.row in 1:nrow(maximum.scars)){
  c.cell.type <- maximum.scars$Cell.type[c.row]
  print(
    ggplot(cell.scar.count[cell.scar.count$Cell.type == c.cell.type, ]) +
      geom_histogram(aes(x = Scars), binwidth = 1) +
      labs(title = c.cell.type)
  )
  maximum.scars$Maximum[c.row] <-
    readline(prompt = paste("Max number scars for cell in cell.type ", 
                            c.cell.type, "? ", sep = ""))
}
maximum.scars$Maximum <- as.integer(maximum.scars$Maximum)

# postscript("./Images/Hr10_cell_type_scar_counts_with_cutoff.eps",
#     width = 7, height = 4.5)
ggplot() +
  geom_histogram(data = cell.scar.count, aes(x = Scars), binwidth = 1) +
  geom_vline(data = maximum.scars, aes(xintercept = Maximum + 0.5), color = "red") +
  facet_wrap(~ Cell.type) +
  labs(y = "Count") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 4))
# dev.off()

# maximum.scars <- read.csv("./Data/2017_10X_6/A6_max_scars_3Adults.csv",
#                           stringsAsFactors = F)
cell.scar.count <- merge(cell.scar.count, maximum.scars)

cells.too.many.scars <-
  cell.scar.count$Cell[cell.scar.count$Scars > cell.scar.count$Maximum]
scars.output.2 <- scars.output.2[!(scars.output.2$Cell %in% cells.too.many.scars), ]

# Write final output ####
write.csv(scars.output.2, "./Data/scars/Hr10_used_scars.csv",
          row.names = F, quote = F)
write.csv(maximum.scars, "./Data/scars/Hr10_max_scars.csv")
