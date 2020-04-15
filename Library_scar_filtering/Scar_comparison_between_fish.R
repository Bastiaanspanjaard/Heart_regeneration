# Dependencies/parameters ####
source("~/Documents/Projects/TOMO_scar/Scripts/linnaeus-scripts/scar_helper_functions.R")
min.presence <- 2 # Cells that have to have a scar for it to be considered
# present in an organism

# Load data ####
libraries <- data.frame(Library_name = c("H5", "H6", "H7", "H8a", "H8v",
                                         "Hr1", "Hr2a", "Hr2b", "Hr3", "Hr4",
                                         "Hr6a", "Hr6v", "Hr7a", "Hr7v",
                                         "Hr10", "Hr11", "Hr12", "Hr13",
                                         "Hr14", "Hr15", "Hr19", "Hr20",
                                         "Hr21", "Hr22", "Hr23", "Hr24",
                                         "Hr25", "Hr26", "Hr27"),
                        Sample =  c("H5", "H6", "H7", "H8", "H8",
                                    "Hr1", "Hr2", "Hr2", "Hr3", "Hr4",
                                    "Hr6", "Hr6", "Hr7", "Hr7",
                                    "Hr10", "Hr11", "Hr12", "Hr13",
                                    "Hr14", "Hr15", "Hr19", "Hr20",
                                    "Hr21", "Hr22", "Hr23", "Hr24",
                                    "Hr25", "Hr26", "Hr27"))

# Loop over sample names. If sample name has >1 library, load all libraries
# and rbind them. Outcome: list of samples, each with its scar dataframe.
scar_list <- list()
for(sample_name in unique(libraries$Sample)){
  scars_to_list <- data.frame()
  for(library_name in libraries$Library_name[libraries$Sample == sample_name]){
    scars_to_list <- rbind(scars_to_list, 
                           read.csv(paste("./Data/scars/Seurat_bc_filter/", library_name, "_used_scars.csv", sep = ""),
                                                   stringsAsFactors = F))
  }
  scar_list[[length(scar_list) + 1]] <- list(scars = scars_to_list)
  names(scar_list)[length(scar_list)] <- sample_name
}

# Count presence of scars, collect all scars, their frequencies per sample, and CIGARs ####
unique.scars <- data.frame(Sequence = character())
all.CIGARs <- data.frame()
for(sample_number in 1:length(scar_list)){
  sample_scars <- scar_list[[sample_number]]$scars
  scar_frequencies <- data.frame(table(sample_scars$Sequence))
  colnames(scar_frequencies) <- c("Sequence", paste("Freq.", names(scar_list)[sample_number], sep = ""))
  scar_list[[sample_number]]$Frequencies <- scar_frequencies
  unique.scars <- merge(unique.scars, scar_frequencies, all = T)
  all.CIGARs <- unique(rbind(all.CIGARs, sample_scars[, c("Sequence", "CIGAR")]))
}
unique.scars[is.na(unique.scars)] <- 0
unique.scars$Presence <- apply(unique.scars[, -1], 1,
                               function(x) sum(x >= min.presence))
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
unique.scars <- merge(unique.scars, scar.probabilities[, c("Sequence", "p", "Embryos")],
                      by.x = "Sequence.short", by.y = "Sequence", all.x = T)
unique.scars$Embryos[is.na(unique.scars$Embryos)] <- 0
unique.scars$p[is.na(unique.scars$p)] <- min(unique.scars$p, na.rm = T)/100
unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

# Write output ####
for(sample_number in 1:length(scar_list)){
  scars_out <- merge(scar_list[[sample_number]]$scars, unique.scars[, c("Sequence", "Presence", "p", "Embryos", "Scar")])
  print(cat("Sample", names(scar_list)[sample_number], "has", nrow(scars_out), "scars in", length(unique(scars_out$Cell)), "cells"))
  write.csv(scars_out,
        file = paste("./Data/scars/Seurat_bc_filter/", names(scar_list)[sample_number], "_scars_compared.csv", sep = ""),
        quote = F, row.names = F)
}

