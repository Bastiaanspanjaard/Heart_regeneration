metadata.in <- read.csv("./Data/final_metadata.csv", stringsAsFactors = F)
colnames(metadata.in)[1] <- "Cell"

# Cell type decisions as follows: all T-cell types and all macrophage types become T-cells/macrophages. 
# Split endocardium into 1 and 2, A and V but merge endocardium frzb.
metadata.in$Cell.type <- metadata.in$lineage.ident
metadata.in$Cell.type[grepl("T-cell", metadata.in$Cell.type)] <- "T-cells"
metadata.in$Cell.type[grepl("Macrophage", metadata.in$Cell.type)] <- "Macrophages"
metadata.in$Cell.type[grepl("Endocardium", metadata.in$Cell.type)] <- 
  metadata.in$final.zoom[grepl("Endocardium", metadata.in$Cell.type)]
metadata.in$Cell.type[grepl("Endocardium frzb", metadata.in$Cell.type)] <- "Endocardium (frzb)"

celltypes <- data.frame(table(metadata.in$Cell.type))
colnames(celltypes)[1] <- c("Cell.type")

celltype_colors <- data.frame(lapply(readRDS("./Data/col.table.rds")[, ], as.character), stringsAsFactors = F)
celltype_colors <- celltype_colors[celltype_colors$set1 != "Perivascular cells" | is.na(celltype_colors$set1), ]
celltype_colors$Cell.type <- celltype_colors$setCM
celltype_colors$Cell.type[is.na(celltype_colors$Cell.type)] <- celltype_colors$setFibro[is.na(celltype_colors$Cell.type)]
# For the remaining cell types, take the color according to set1 if the color has not been taken yet.
celltype_colors$Cell.type[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)] <- 
  celltype_colors$set1[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)]


endo_add <- celltype_colors[celltype_colors$Cell.type %in% c("Endocardium (V)", "Endocardium (A)"), ]
endo_add$Cell.type[endo_add$Cell.type == "Endocardium (V)"] <- 
  "Endocardium 2 (V)"
endo_add$Cell.type[endo_add$Cell.type == "Endocardium (A)"] <- 
  "Endocardium 2 (A)"
celltype_colors$Cell.type[celltype_colors$Cell.type == "Endocardium (V)"] <- 
  "Endocardium 1 (V)"
celltype_colors$Cell.type[celltype_colors$Cell.type == "Endocardium (A)"] <- 
  "Endocardium 1 (A)"
celltype_colors <- rbind(celltype_colors, endo_add)
saveRDS(celltype_colors, file = "./Data/Cell_type_colors.rds")