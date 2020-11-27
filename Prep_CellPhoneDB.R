require(Matrix)
require(rhdf5)

h5_to_sparse_count <- function(h5_counts){
  sparse_count_matrix <- Matrix::sparseMatrix(  
    i = h5_counts$indices,
    p = h5_counts$indptr,
    x = h5_counts$data,
    dimnames = list(
      h5_counts$features$id,
      h5_counts$barcodes
    ),
    dims = h5_counts$shape,
    index1 = FALSE
  )
  
  return(sparse_count_matrix)
}

# Vector of libraries to combine to automate loading further one (loop over array to load, convert to sparse, rename cells)
# Control
# libraries_incl <- c("H5", "H6", "H7", "H8a", "H8v")
# output_prefix <- "Control"
# 3dpi
# libraries_incl <- c("Hr10", "Hr11", "Hr12", "Hr22", "Hr23",
#                     "Hr24", "Hr25", "Hr26", "Hr27", "Hr28",
#                     "Hr29", "Hr34", "Hr35")
# output_prefix <- "3dpi"
# 7dpi
# libraries_incl <- c("Hr1", "Hr2a", "Hr2b", "Hr6a", "Hr6v",
                    # "Hr7a", "Hr7v", "Hr8", "Hr9", "Hr13",
                    # "Hr14", "Hr15", "Hr30", "Hr31", "Hr32",
                    # "Hr33")
# output_prefix <- "7dpi"

# Load and concatenate data ####
cell_types <- read.csv("/local/Bastiaan/Projects/heart_Bo/Data/final_metadata_Tmacromerged_2.csv", stringsAsFactors = F)
# cell_types <- read.csv("/local/Bo/Remap_allhearts/new_clustering_inhib3/celltypes_zoom_allcells.csv", stringsAsFactors = F)
# cell_types <- cell_types[cell_types$orig.ident %in% libraries_incl, c("X", "orig.ident", "work.ident")] # For rough cell types
# cell_types <- cell_types[cell_types$orig.ident %in% libraries_incl, c("X", "orig.ident", "zoom.subtypes")] # For zoomed cell types
# colnames(cell_types) <- c("Cell", "Library", "cell_type")
HR_setnames <- read.csv("/local/Bastiaan/Projects/heart_Bo/Data/HR_setnames.csv", stringsAsFactors = F)[, -1]
orthologue_list_old <- read.csv("/local/Bastiaan/Projects/heart_Bo/Data/GRCz11_GRCh38_orthologues.txt", sep = "\t")
orthologue_list <- 
  read.csv("/local/Bastiaan/Projects/heart_Bo/Data/Alliance_gene_name_conversion.scsv", 
           sep = ";", stringsAsFactors = F)[, c("Gene.stable.ID", "Human.gene.stable.ID", "Human.gene.name", "Gene.name")]

timepoint <- 0
inhib <- F

for(timepoint in c(0, 3, 7)){
  for(inhib in c(T, F)){
    if(timepoint == 0 && inhib){next}
    if(inhib){
      inhib_v <- "IWR1"
    }else{
      inhib_v <- c("no", "DMSO")
    }
    output_prefix <- paste("Time_", timepoint, "_inhib_", inhib, sep = "")
    print(output_prefix)
    libraries_incl <- HR_setnames$heart[HR_setnames$dpi == timepoint & HR_setnames$inhib %in% inhib_v]
    print(libraries_incl)
    # l <- 1
    for(l in 1:length(libraries_incl)){
      clib <- libraries_incl[l]
      count_data_in <- 
        h5read(paste("/local/Bastiaan/Projects/heart_Bo/velocyto/", clib, "_v3Dr11/outs/filtered_feature_bc_matrix.h5", sep = ""),
               name = "matrix")
      sparse_count <- h5_to_sparse_count(count_data_in)
      print(dim(sparse_count))
      sparse_count@Dimnames[[2]] <- 
        paste(clib, sapply(sparse_count@Dimnames[[2]], 
                           function(x){unlist(strsplit(x, "-"))[1]}), sep = "_")
      sparse_count <- sparse_count[, colnames(sparse_count) %in% cell_types$Cell[cell_types$orig.ident == clib]]
      print(dim(sparse_count))
      
      # Read count matrices have all gene names as observables. Restrict this to only those genes that have an
      # orthologue.
      if(l == 1){
        ensemble_genes <- data.frame(zGene = sparse_count@Dimnames[[1]],
                                     Order = 1:length(sparse_count@Dimnames[[1]]))
        ensemble_genes <- merge(ensemble_genes, orthologue_list[, c("Gene.stable.ID", "Human.gene.stable.ID")],
                                by.x = "zGene", by.y = "Gene.stable.ID", all.x = T)
        colnames(ensemble_genes)[3] <- "Gene"
        ensemble_genes$Gene[ensemble_genes$zGene %in% ensemble_genes$zGene[duplicated(ensemble_genes$zGene)]] <- ""
        ensemble_genes$Gene[is.na(ensemble_genes$Gene)] <- ""
        ensemble_genes <- unique(ensemble_genes)  
        ensemble_genes <- ensemble_genes[order(ensemble_genes$Order), ]
        
        sparse_count <- sparse_count[rownames(sparse_count) %in% ensemble_genes$zGene[ensemble_genes$Gene != ""], ]
        print(dim(sparse_count))
        
        sparse_count_combined <- sparse_count
      }else{
        sparse_count <- sparse_count[rownames(sparse_count) %in% ensemble_genes$zGene[ensemble_genes$Gene != ""], ]
        print(dim(sparse_count))
        
        sparse_count_combined <- cbind(sparse_count_combined, sparse_count)
      }
    }
    count_matrix_combined <- as.matrix(sparse_count_combined)
    count_matrix_combined <- cbind(data.frame(Gene = ensemble_genes$Gene[ensemble_genes$Gene != ""]),
                          count_matrix_combined)
    count_matrix_combined[, -1] <- 10000 * t(t(count_matrix_combined[, -1])/colSums(count_matrix_combined[, -1]))
    count_matrix_combined[, -1] <- log1p(count_matrix_combined[, -1])

    write.table(count_matrix_combined,
                file = paste("/local/Bastiaan/Projects/heart_Bo/Data/", output_prefix, "alliance_ln_expr.txt", sep = ""),
                sep = "\t", quote = F, row.names = F)
    cell_types_here <- cell_types[cell_types$Cell %in% colnames(count_matrix_combined), c("Cell", "Cell_type")]
    cell_types_order <- data.frame(Cell = colnames(count_matrix_combined), 
                                  Order = 1:ncol(count_matrix_combined))
    cell_types_here <- merge(cell_types_here, cell_types_order)
    # cell_types_here <- cell_types_here[as.character(colnames(count_matrix_combined)), ]
    write.table(cell_types_here[, c("Cell", "Cell_type")],
      #cell_types[cell_types$Cell %in% colnames(count_matrix_combined), c("Cell", "Cell_type")], 
                file = paste("/local/Bastiaan/Projects/heart_Bo/Data/", output_prefix, "alliance_zoom_celltypes.txt", sep = ""), 
                sep = "\t", quote = F, row.names = F) # For zoomed cell types
  }
}

# count_matrix_test <- as.matrix(sparse_count)
# For the combined 3dpi datasets, the sparse matrix is 1Gb and de-sparsing a ~200 Mb
# sparse matrix gives a 3.3Gb regular matrix, i.e. a 15-fold size increase.

# Convert single-cell expression matrix to orthologues ###
# ensemble_genes <- data.frame(zGene = sparse_count_combined@Dimnames[[1]],
#                              Order = 1:length(sparse_count_combined@Dimnames[[1]]))
# ensemble_genes <- merge(ensemble_genes, orthologue_list[, c("Gene.stable.ID", "Human.gene.stable.ID")],
#                         by.x = "zGene", by.y = "Gene.stable.ID", all.x = T)

# Remove duplicate orthologues (at least for now) by setting their human gene name to "" (which will be removed from
# the count data later on)
# colnames(ensemble_genes)[3] <- "Gene"
# ensemble_genes$Gene[ensemble_genes$zGene %in% ensemble_genes$zGene[duplicated(ensemble_genes$zGene)]] <- ""
# ensemble_genes$Gene[is.na(ensemble_genes$Gene)] <- ""
# ensemble_genes <- unique(ensemble_genes)  
# ensemble_genes <- ensemble_genes[order(ensemble_genes$Order), ]
# count_matrix_combined <- cbind(data.frame(Gene = ensemble_genes$Gene[ensemble_genes$Gene != ""]),
#                       count_matrix_combined)

# Log-normalize counts ####
# count_matrix_combined <- 
#   count_matrix_combined[count_matrix_combined$Gene != "", 
#                         colnames(count_matrix_combined) %in% c("Gene", cell_types$Cell)]
# ln_counts_combined <- count_matrix_combined
count_matrix_combined[, -1] <- 10000 * t(t(count_matrix_combined[, -1])/colSums(count_matrix_combined[, -1]))
count_matrix_combined[, -1] <- log1p(count_matrix_combined[, -1])

# Save for input in CellPhoneDB ####
# write.table(count_matrix_combined, 
#             file = paste("/local/Bastiaan/Projects/heart_Bo/Data/", output_prefix, "_ln_expr.txt", sep = ""), 
#             sep = "\t", quote = F, row.names = F)
# write.table(cell_types[cell_types$Cell %in% colnames(count_matrix_combined), c("Cell", "cell_type")], 
#             file = paste("/local/Bastiaan/Projects/heart_Bo/Data/", output_prefix, "_celltypes.txt", sep = ""), 
#             sep = "\t", quote = F, row.names = F) # For rough cell types
write.table(cell_types[cell_types$Cell %in% colnames(count_matrix_combined), c("Cell", "cell_type")], 
            file = paste("/local/Bastiaan/Projects/heart_Bo/Data/", output_prefix, "_zoom_celltypes.txt", sep = ""), 
            sep = "\t", quote = F, row.names = F) # For zoomed cell types
