# Cell type renaming ####
# Cell type decisions as follows: all T-cell types and all macrophage types become T-cells/macrophages. 
# Merge endocardium frzb and do not split endocardium into 1 and 2.
# Merge cardiomyocytes A and V but keep the regular, ttn.2 and proliferating subtypes.
cell_annotations <- read.csv("./Data/final_metadata.csv", stringsAsFactors = F)
colnames(cell_annotations)[1] <- "Cell"
cell_annotations$Cell_type <- cell_annotations$lineage.ident
cell_type_renaming <-
  data.frame(Original_cell_type = sort(as.character(unique(cell_annotations$Cell_type))), stringsAsFactors = F)
cell_type_renaming$New_type_name <-
  c("B-cells", "Bl.ves.EC (apnln)", "Bl.ves.EC (lyve1)", 
    "Bl.ves.EC (plvapb)", "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2)",
    "Cardiomyocytes (ttn.2)", "Cardiomyocytes (Atrium)", "Cardiomyocytes (Ventricle)",
    "Dead cells", "Endocardium (Atrium)", "Endocardium (Ventricle)", "Endocardium (frzb)",
    "Endocardium (frzb)", "Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblasts (const.)",
    "Fibroblasts (cfd)", "Fibroblasts (col11a1a)", "Fibroblasts (col12a1a)", "Fibroblasts (cxcl12a)",
    "Fibroblasts (mpeg1.1)", "Fibroblasts (nppc)", "Fibroblasts (proliferating)", "Fibroblasts (spock3)",
    "Valve fibroblasts", "Macrophages", "Macrophages", "Macrophages",
    "Macrophages", "Macrophages", "Macrophages", 
    "Macrophages", "Macrophages", "Macrophages",
    "Macrophages", "Monocytes", "Myelin cells", "Neuronal cells", "Neutrophils",
    "Perivascular cells", "Proliferating cells", "Smooth muscle cells", "T-cells",
    "T-cells", "T-cells")
cell_annotations$Cell_type <-
  cell_type_renaming$New_type_name[match(cell_annotations$Cell_type, cell_type_renaming$Original_cell_type)]
cell_annotations <- cell_annotations[cell_annotations$Cell_type != "Dead cells", ]
cell_annotations$Cell_name <- paste("nd", cell_annotations$Cell, sep = "")
# write.csv(cell_annotations, "./Data/final_metadata_Tmacromerged_2.csv", row.names = F, quote = F)

# Cell type colors ####
celltype_colors <- data.frame(lapply(readRDS("./Data/col.table.rds")[, ], as.character), stringsAsFactors = F)
celltype_colors$Cell.type <- celltype_colors$set1
celltype_colors$Cell.type[celltype_colors$Cell.type == "Cardiomyocytes (ttn.2) V"] <- NA
celltype_colors$Cell.type[is.na(celltype_colors$Cell.type)] <- celltype_colors$setFibro[is.na(celltype_colors$Cell.type)]
# For the remaining cell types, take the color according to set1 if the color has not been taken yet.
celltype_colors$Cell.type[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)] <- 
  celltype_colors$set1[!(celltype_colors$set1 %in% celltype_colors$Cell.type) & !is.na(celltype_colors$set1) & is.na(celltype_colors$Cell.type)]
celltype_colors <- celltype_colors[celltype_colors$Cell.type != "Fibroblasts", ]
# This is actually the right color for neuronal cells
neuronal_add <- celltype_colors[celltype_colors$Cell.type == "Myelin cells" & !is.na(celltype_colors$set1), ]
neuronal_add$Cell.type <- "Neuronal cells"
celltype_colors <- rbind(celltype_colors, neuronal_add)
cell_type_renaming <-
  data.frame(Original_cell_type = sort(as.character(celltype_colors$Cell.type)), stringsAsFactors = F)
cell_type_renaming$New_type_name <-
  c("B-cells", "Bl.ves.EC (apnln)", "Bl.ves.EC (lyve1)", 
    "Bl.ves.EC (plvapb)", "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2)",
    "Cardiomyocytes (Atrium)", "Cardiomyocytes (Ventricle)",
    "Endocardium (Atrium)", "Endocardium (frzb)", "Endocardium (Ventricle)", 
    "Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblasts (const.)",
    "Fibroblasts (cfd)", "Fibroblasts (col11a1a)", "Fibroblasts (col12a1a)", "Fibroblasts (cxcl12a)",
    "Fibroblasts (mpeg1.1)", "Fibroblasts (nppc)", "Fibroblasts (proliferating)", "Fibroblasts (spock3)",
    "Valve fibroblasts", "Macrophages", "Monocytes", "Myelin cells", "Neuronal cells", "Neutrophils",
    "Perivascular cells", "Proliferating cells", "Smooth muscle cells", "T-cells")
celltype_colors$Cell.type <-
  cell_type_renaming$New_type_name[match(celltype_colors$Cell.type, cell_type_renaming$Original_cell_type)]
# write.csv(celltype_colors[!is.na(celltype_colors$Cell.type), c("Cell.type", "color")],
#           "./Data/Cell_type_colors_2.csv", quote = F, row.names = F)
