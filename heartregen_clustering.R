# global functions ####
combine_plots <- function(...,
                          title.text = NULL,
                          title.color = "black",
                          title.size = 16,
                          title.vjust = 0.5,
                          title.hjust = 0.5,
                          title.fontface = "bold",
                          caption.text = NULL,
                          caption.color = "black",
                          caption.size = 10,
                          caption.vjust = 0.5,
                          caption.hjust = 0.5,
                          caption.fontface = "plain",
                          sub.text = NULL,
                          sub.color = "black",
                          sub.size = 12,
                          sub.vjust = 0.5,
                          sub.hjust = 0.5,
                          sub.fontface = "plain",
                          sub.x = 0.5,
                          sub.y = 0.5,
                          sub.vpadding = ggplot2::unit(1, "lines"),
                          sub.angle = 0,
                          sub.lineheight = 0.9,
                          title.rel.heights = c(0.1, 1.2),
                          caption.rel.heights = c(1.2, 0.1),
                          title.caption.rel.heights = c(0.1, 1.2, 0.1)) {
  # preparing the basic plot
  plot <- cowplot::plot_grid(...)
  
  # preparing the title
  if (!is.null(title.text)) {
    title <-
      cowplot::ggdraw() +
      cowplot::draw_label(
        label = title.text,
        colour = title.color,
        size = title.size,
        vjust = title.vjust,
        hjust = title.hjust,
        fontface = title.fontface
      )
  }
  # preparing the caption
  if (!is.null(caption.text)) {
    caption <-
      cowplot::ggdraw() +
      cowplot::draw_label(
        label = caption.text,
        colour = caption.color,
        size = caption.size,
        vjust = caption.vjust,
        hjust = caption.hjust,
        fontface = caption.fontface
      )
  }
  
  # combining the basic plot with the either title or caption or title and
  # caption
  if (!is.null(title.text)) {
    if (!is.null(caption.text)) {
      # if both title and caption are needed
      plot <-
        cowplot::plot_grid(title,
                           plot,
                           caption,
                           ncol = 1,
                           rel_heights = title.caption.rel.heights
        )
    } else {
      # if only title is needed
      plot <-
        cowplot::plot_grid(title,
                           plot,
                           ncol = 1,
                           rel_heights = title.rel.heights
        )
    }
  } else if (!is.null(caption.text)) {
    # if only caption is needed
    plot <-
      cowplot::plot_grid(plot,
                         caption,
                         ncol = 1,
                         rel_heights = caption.rel.heights
      )
  }
  
  # finally adding sub if it's needed
  if (!is.null(sub.text)) {
    plot <-
      cowplot::ggdraw(
        cowplot::add_sub(
          plot = plot,
          label = sub.text,
          x = sub.x,
          y = sub.y,
          vpadding = sub.vpadding,
          colour = sub.color,
          size = sub.size,
          vjust = sub.vjust,
          hjust = sub.hjust,
          fontface = sub.fontface,
          angle = sub.angle,
          lineheight = sub.lineheight
        )
      )
  }
  
  # return the final, combined plot
  return(plot)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}



## load and integrate data in one seurat object, calculate mito reads, filter cells, annotate cell types. ####

library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(reshape2)

#mito.genes names
mito.genes <- read.table("/local/users/Bastiaan/Projects/heart_Bo/Data/mito.genes.vs.txt",sep = ",")
mito.genes <- mito.genes$V3
mito.genes <- as.character(mito.genes)

# load all data
batches <- c("H5","H6","H7","H8a","H8v","Hr1","Hr2","Hr2","Hr3","Hr4","Hr6a","Hr6v","Hr7a","Hr7v","Hr8",
             "Hr9","Hr10","Hr11","Hr12","Hr13","Hr14","Hr15","Hr19","Hr20","Hr21","Hr22","Hr23",
             "Hr24","Hr25","Hr26","Hr27","Hr28","Hr29","Hr30","Hr31","Hr32","Hr33","Hr34","Hr35")
gather.data <- list()
data.stats <- numeric()

for (i in batches) {
  gather.data[[i]] <- Read10X_h5(filename = paste0(i,"_filtered_matrix.h5") )
  RFP.index <- grep(pattern = "^RFP", x = rownames(gather.data[[i]]), value = FALSE) # Select row indices and not ERCC names 
  gather.data[[i]] <- gather.data[[i]][-RFP.index, ]
  gather.data[[i]] <- CreateSeuratObject(counts = gather.data[[i]],
                                         min.cells = 3, min.features = 150,
                                         project = i)
  mito.genes.use <- setdiff(mito.genes,setdiff(mito.genes,rownames(gather.data[[i]][["RNA"]])))
  gather.data[[i]][["percent.mito"]] <- PercentageFeatureSet(object = gather.data[[i]], features = mito.genes.use)
  gather.data[[i]] <- subset(x = gather.data[[i]], subset = percent.mito < 25 & nFeature_RNA < 4100)
  data.stats[i] <- length(colnames(x = gather.data[[i]]))
}

# stats
print(data.stats)
# merge all data
all.hearts <- merge(gather.data[[1]], y = gather.data[-1], add.cell.ids = batches, project = "allheart")

#remove raw data to free space
rm(gather.data)

#Annotate metadata
all.hearts@meta.data$time <- NA
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H5","H6","H7","H8a","H8v"),]$time <- "Ctrl"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr10","Hr11","Hr12","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr28","Hr29","Hr34","Hr35"),]$time <- "3dpi"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr6a","Hr6v","Hr7a","Hr7v","Hr30","Hr31","Hr32","Hr33"),]$time <- "7dpi"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr3","Hr4","Hr19","Hr20","Hr21"),]$time <- "30dpi"
all.hearts@meta.data$time <- factor(x = all.hearts@meta.data$time, levels = c("Ctrl","3dpi","7dpi","30dpi"))

all.hearts@meta.data$AV <- "Wholeheart"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H8a","Hr6a","Hr7a"),]$AV <- "Atrium"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H8v","Hr6v","Hr7v","Hr25"),]$AV <- "Ventricle"

all.hearts@meta.data$inhib <- "NULL"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr28","Hr30"),]$inhib <- "DMSO"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr29","Hr31","Hr32","Hr33","Hr34","Hr35"),]$inhib <- "IWR1"

ncol(all.hearts)

# add annotated cell type names with filtering out erythrocytes and dead cells
cell.id.table <- read.csv("cell.id.table.csv")
all.hearts <- AddMetaData(all.hearts,cell.id.table)
all.hearts <- subset(all.hearts,cells = rownames(all.hearts@meta.data)[!is.na(all.hearts@meta.data$celltypes)],)

## normalize, pca, cluster and plot umap
all.hearts <- NormalizeData(object = all.hearts)
all.hearts <- FindVariableFeatures(object = all.hearts, selection.method = "vst",nfeatures = 3500)
all.hearts <- ScaleData(object = all.hearts,vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito"))
all.hearts <- RunPCA(object = all.hearts, features = VariableFeatures(object = all.hearts),npcs = 100)
ElbowPlot(object = all.hearts,ndims = 100)

all.hearts <- FindNeighbors(object = all.hearts, dims = 1:50)
all.hearts <- FindClusters(object = all.hearts,resolution = 4)
all.hearts <- RunUMAP(all.hearts, dims = 1:15, verbose = F)

# annotate cell types, annotation is based on differental gene expression of each cluster
#annotate cell types without ery ####
ct <- data.frame(cluster = 0:102, Cell.type = 0:102)
ct[ct$cluster %in% c(87,101,18,94,83,78,18,42,9,45,82,44,11,77,46,20,59,36,93,40,58,62,95,97,43,38,16),]$Cell.type <- "Immune Cells"
ct[ct$cluster %in% c(57),]$Cell.type <- "B-cells"
ct[ct$cluster %in% c(89,79,64),]$Cell.type <- "Neutrophils"
ct[ct$cluster %in% c(44,46,11,77,20,36,93,40,59),]$Cell.type <- "Macrophages"
ct[ct$cluster %in% c(94,83,87,78,101,18,42,9,82,45),]$Cell.type <- "Macrophages"
ct[ct$cluster %in% c(62,95),]$Cell.type <- "Macrophages"
ct[ct$cluster %in% c(58,43,38,16,97),]$Cell.type <- "T-cells"
ct[ct$cluster %in% c(34,32,29,75,73,23,24,3,84,53,54,60,22,86,71),]$Cell.type <-"Fibroblasts"
ct[ct$cluster %in% c(17,74,68,1,19,69,81,31,92,98),]$Cell.type <-"Smooth muscle cells"
ct[ct$cluster %in% c(50,49,99),]$Cell.type <-"Bl.ves.EC (plvapb)"
ct[ct$cluster %in% c(66),]$Cell.type <-"Bl.ves.EC (lyve1)"
ct[ct$cluster %in% c(100,48,47),]$Cell.type <-"Bl.ves.EC (apnln)"
ct[ct$cluster %in% c(0,6,27,51,72,13,15,70,37),]$Cell.type <-"Endocardium 1 (A)"
ct[ct$cluster %in% c(2,5,28,67,33,7,21,91),]$Cell.type <-"Endocardium 1 (V)"
ct[ct$cluster %in% c(70,37),]$Cell.type <-"Endocardium 2 (A)"
ct[ct$cluster %in% c(21,91),]$Cell.type <-"Endocardium 2 (V)"
ct[ct$cluster %in% c(15),]$Cell.type <-"Endocardium frzb (A)"
ct[ct$cluster %in% c(55),]$Cell.type <-"Endocardium frzb (V)"
ct[ct$cluster %in% c(41,61,10,26),]$Cell.type <-"Cardiomyocytes A"
ct[ct$cluster %in% c(30,56,12,14,8,63,52,76),]$Cell.type <-"Cardiomyocytes V"
ct[ct$cluster %in% c(4,39,35),]$Cell.type <-"Cardiomyocytes (ttn.2) A"
ct[ct$cluster %in% c(80),]$Cell.type <-"Cardiomyocytes (ttn.2) V"
ct[ct$cluster %in% c(90),]$Cell.type <-"Cardiomyocytes (proliferating)"
ct[ct$cluster %in% c(65),]$Cell.type <-"Perivascular cells"
ct[ct$cluster %in% c(25,96),]$Cell.type <-"Fibroblast-like cells"
ct[ct$cluster %in% c(102),]$Cell.type <-"Neuronal cells"
ct[ct$cluster %in% c(88),]$Cell.type <-"Myelin cells"
ct[ct$cluster %in% c(85),]$Cell.type <-"Proliferating cells"
final.all.hearts$"first.line.annotation" <- final.all.hearts$"seurat_clusters"
final.all.hearts@meta.data$first.line.annotation <- plyr::mapvalues(final.all.hearts@meta.data$first.line.annotation, from =ct$cluster, to = ct$Cell.type)

final.all.hearts <- SetIdent(final.all.hearts,value = "first.line.annotation")

#final.all.hearts <- SetIdent(final.all.hearts,cells = rownames(niche@meta.data),value = niche@meta.data$work.ident )
#final.all.hearts <- SetIdent(final.all.hearts,cells = rownames(immune@meta.data),value = immune@meta.data$work.ident )

#plot umap
DimPlot(all.hearts,group.by = "celltypes",label = F)


