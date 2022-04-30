# Requirements ####
require(Seurat)
require(ggplot2)
library(cowplot)
library(reshape2)
library(pheatmap)
library(ggbeeswarm)

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

# Define colors for cell types ####
colors <- readRDS(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/color_scheme_seurat.rds")
col.table <- colors

set1 <- c("Epicardium (A)","Macrophages","T-cells","Neutrophils","B-cells",
  "Monocytes","Fibroblasts","Epicardium (V)",
  "Valve fibroblasts","Perivascular cells","Endocardium (A)","Smooth muscle cells",
  "Endocardium (V)","Cardiomyocytes (ttn.2) A","Cardiomyocytes V","Cardiomyocytes A", # TEST - DOES ADDING THE "A" to CM (ttn.2) not break anything?
  "Endocardium (frzb)","Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Bl.ves.EC (lyve1)","Cardiomyocytes (ttn.2) V",
  "Proliferating cells","Cardiomyocytes (proliferating)","Myelin cells")
names(set1) <- c("creamyellow2","winered","pink","syn-magenta","syn-magentalight",
                 "darkpink","yellow","darkyellow",
                 "yellowochre","darkpurple","blue","ocker",
                 "darkblue","green","darkgreen","middlegreen",
                 "mattblue","orange","mattorange","darkorange","green2",
                 "white","viridian","cream")
set1 <- data.frame(name = names(set1),set1 = set1)
col.table <- merge(x = col.table,y = set1,by.x = "name",by.y = "name",all.x = T)  

setCM <- c("Cardiomyocytes (ttn.2) V","Cardiomyocytes (ttn.2) A",
           "Cardiomyocytes V","Cardiomyocytes A",
           "Cardiomyocytes (proliferating)")
names(setCM) <- c("green2","green","darkgreen",
  "middlegreen","viridian")

setCM <- data.frame(name = names(setCM),setCM = setCM)
col.table <- merge(x = col.table,y = setCM,by.x = "name",by.y = "name",all.x = T)

setFibro <- c("Fibroblasts (const.)","Epicardium (V)","Epicardium (A)",
             "Fibroblasts (cxcl12a)","Fibroblasts (col11a1a)","Fibroblasts (cfd)",
             "Valve fibroblasts","Fibroblasts (col12a1a)","Perivascular cells",
             "Fibroblasts (nppc)","Fibroblasts (spock3)","Fibroblasts (mpeg1.1)",
             "Fibroblasts (prolif.)")
names(setFibro) <-      c("syn-lightyellow","darkyellow","creamyellow2",
                          "lightpink","lightblue","green2",
                          "yellowochre","syn-red","darkpurple",
                          "purple","lightpurple","white2",
                          "darkblue2")

setFibro <- data.frame(name = names(setFibro),setFibro = setFibro)
col.table <- merge(x = col.table,y = setFibro,by.x = "name",by.y = "name",all.x = T)


setImmune <- c("Macrophage","Macrophage (il1b)","Macrophage (apoeb)",
    "T-cells","Neutrophils","Macrophage (proliferating)",
    "Macrophage (CM duplex)","Macrophage (Endothelia duplex)","Macrophage (epdl)",
    "Macrophage (cd59)","B-cells","Monocytes",
    "Macrophage (Ery duplex)","Macrophage (Fibroblast duplex)","T-cells (il4/13)",
    "T-cells (proliferating)","Dead cells")
names(setImmune) <- c("yellow","syn-red","green",
                     "orange","darkblue","ocker",
                     "leadwhite","leadwhite","pink",
                     "winered","blue","darkgreen",
                     "leadwhite","white2","darkpurple2",
                     "darkorange","leadwhite")

setImmune <- data.frame(name = names(setImmune),setImmune = setImmune)
col.table <- merge(x = col.table,y = setImmune,by.x = "name",by.y = "name",all.x = T)

# saveRDS(col.table,file = "/local/Bo/Remap_allhearts/final_clustering/col.table.rds")


# Get all plot relevant objects ####

## allhearts
load(file = "/local/Bo/Remap_allhearts/final_clustering/final.all.hearts.noErynoDuplex.Robj")
load(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/CM2.Robj")
final.all.hearts$plot.ident <- final.all.hearts$lineage.ident
# NEW
final.all.hearts@meta.data$plot.ident <- plyr::mapvalues(final.all.hearts@meta.data$plot.ident, 
   from =c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
           "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
           "Fibroblast","Smooth muscle cells","Cardiomyocytes (proliferating)",
           "B-cells","Macrophage (CM duplex)","Macrophage (Fibroblast duplex)",
           "Macrophage (Endothelia duplex)","Proliferating cells","Fibroblast (mpeg1.1)",
           "Fibroblast (cfd)","Monocytes","Macrophage (il1b)",
           "T-cells (il4/13)","Macrophages","T-cells",
           "Macrophage (proliferating)","Macrophage (cd59)","Neutrophils",
           "Macrophage (epdl)","T-cells (proliferating)","Macrophage (Ery duplex)",
           "Macrophage (apoeb)","Dead cells","Fibroblast-like cells",
           "Endocardium (frzb)","Fibroblast (nppc)","Fibroblast (spock3)",
           "Myelin cells","Neuronal cells","Perivascular cells",
           "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblast (cxcl12a)",
           "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblast (col12a1a)",
          "Fibroblast (col11a1a)","Bl.ves.EC (lyve1)","Fibroblast (proliferating)")     , 
  to = c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
         "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
         "Fibroblasts","Smooth muscle cells","Cardiomyocytes (proliferating)",
         "B-cells","Macrophages","Macrophages",
         "Macrophages","Proliferating cells","Fibroblasts",
         "Fibroblasts","Monocytes","Macrophages",
         "T-cells","Macrophages","T-cells",
         "Macrophages","Macrophages","Neutrophils",
         "Macrophages","T-cells","Macrophages",
         "Macrophages","Macrophages","Fibroblast-like cells",
         "Endocardium (frzb)","Fibroblasts","Fibroblasts",
         "Myelin cells","Neuronal cells","Perivascular cells",
         "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblasts",
         "Epicardium (A)","Epicardium (V)","Fibroblasts",
         "Fibroblasts","Bl.ves.EC (lyve1)","Fibroblasts")  )
# OLD
# final.all.hearts@meta.data$plot.ident <- 
#   plyr::mapvalues(final.all.hearts@meta.data$plot.ident, 
#                   from =c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
#                           "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
#                           "Fibroblast","Smooth muscle cells","Cardiomyocytes (proliferating)",
#                           "B-cells","Macrophage (CM duplex)","Macrophage (Fibroblast duplex)",
#                           "Macrophage (Endothelia duplex)","Proliferating cells","Fibroblast (mpeg1.1)",
#                           "Fibroblast (cfd)","Monocytes","Macrophage (il1b)",
#                           "T-cells (il4/13)","Macrophages","T-cells",
#                           "Macrophage (proliferating)","Macrophage (cd59)","Neutrophils",
#                           "Macrophage (epdl)","T-cells (proliferating)","Macrophage (Ery duplex)",
#                           "Macrophage (apoeb)","Dead cells","Fibroblast-like cells",
#                           "Endocardium (frzb)","Fibroblast (nppc)","Fibroblast (spock3)",
#                           "Myelin cells","Neuronal cells","Perivascular cells",
#                           "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblast (cxcl12a)",
#                           "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblast (col12a1a)",
#                           "Fibroblast (col11a1a)","Bl.ves.EC (lyve1)","Fibroblast (proliferating)")     , 
#                   to = c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
#                          "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
#                          "Fibroblasts","Smooth muscle cells","Cardiomyocytes (proliferating)",
#                          "B-cells","Macrophages","Macrophages",
#                          "Macrophages","Proliferating cells","Fibroblasts",
#                          "Fibroblasts","Monocytes","Macrophages",
#                          "T-cells","Macrophages","T-cells",
#                          "Macrophages","Macrophages","Neutrophils",
#                          "Macrophages","T-cells","Macrophages",
#                          "Macrophages","Macrophages","Fibroblast-like cells",
#                          "Endocardium (frzb)","Fibroblasts","Fibroblasts",
#                          "Myelin cells","Neuronal cells","Perivascular cells",
#                          "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblasts",
#                          "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblasts",
#                          "Fibroblasts","Bl.ves.EC (lyve1)","Fibroblasts")  )
# plot_grid(
#   DimPlot(final.all.hearts,group.by = "plot.ident", label = T)+NoLegend(),
#   DimPlot(final.all.hearts, label = T,group.by = "seurat_clusters")+NoLegend()
# )
final.all.hearts <- SetIdent(final.all.hearts,value = "plot.ident")
# final.all.hearts$plot.ident2 <- plyr::mapvalues(final.all.hearts$plot.ident2, from ="Cardiomyocytes (ttn.2)", to = "Cardiomyocytes (ttn.2) A")
final.all.hearts <- SetIdent(final.all.hearts, cells = WhichCells(final.all.hearts,expression = seurat_clusters %in% c(36,73) ),value = "Epicardium (Atrium)")

final.all.hearts <- SetIdent(final.all.hearts, 
                             cells = WhichCells(CM,expression = work.ident == "Cardiomyocytes (proliferating)" ),
                             value = "Cardiomyocytes (proliferating)")
final.all.hearts <- SetIdent(final.all.hearts, 
                             cells = WhichCells(CM,expression = work.ident == "Cardiomyocytes (ttn.2) V" ),
                             value = "Cardiomyocytes (ttn.2) V")


# DimPlot(final.all.hearts,label = T)
# final.all.hearts <- StashIdent(final.all.hearts,save.name = "plot.ident2")
final.all.hearts[["plot.ident2"]] <- Idents(object = final.all.hearts)
# final.all.hearts$plot.ident2 <- plyr::mapvalues(final.all.hearts$plot.ident2, from ="Cardiomyocytes (ttn.2)", to = "Cardiomyocytes (ttn.2) A")

# BiocManager::install("SingleCellExperiment")
# final.all.hearts.small <- DietSeurat(final.all.hearts)
# final.all.hearts.sce <- as.SingleCellExperiment(final.all.hearts.small)
# final.all.hearts.sce@assays@data@listData$logcounts <- NULL
# final.all.hearts.sce@colData@listData <- final.all.hearts.sce@colData@listData[c("orig.ident", "time", "morphine", "AV", "inhib", "plot.ident2")]
# save(final.all.hearts.sce, file = "/local/users/Bastiaan/Projects/heart_Bo/Data/SCE_fullobject.Robj")

all.heart.norm.counts <- data.frame(table(final.all.hearts@meta.data$orig.ident))

# DimPlot(final.all.hearts, group.by = "plot.ident2",label = F,
#         cols = colors$color[match(c("viridian","creamyellow2","winered","pink","syn-magenta","syn-magentalight",
#                                     "darkpink","yellow","darkyellow",
#                                     "yellowochre","darkpurple","blue","ocker",
#                                     "darkblue","green","darkgreen","middlegreen",
#                                     "mattblue","orange","mattorange","darkorange","green2",
#                                     "white","cream","cream"),colors$name)] )



# cols = c("Epicardium (Atrium)","Macrophages","T-cells","Neutrophils","B-cells",
#           "Monocytes","Fibroblasts","Epicardium (Ventricle)",
#           "Fibroblast-like cells","Perivascular cells","Endocardium (A)","Smooth muscle cells",
#           "Endocardium (V)","Cardiomyocytes (ttn.2)","Cardiomyocytes V","Cardiomyocytes A",
#           "Endocardium (frzb)","Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Bl.ves.EC (lyve1)","Cardiomyocytes (ttn.2) V",
#           "Proliferating cells","Myelin cells","Cardiomyocytes (proliferating)","Neuronal cells"))+NoLegend()


## CM
# load(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/CM.Robj")
# final.all.hearts <- SetIdent(final.all.hearts,value = "lineage.ident")
# CM <- subset(final.all.hearts,idents = c("Cardiomyocytes A","Cardiomyocytes V","Cardiomyocytes (ttn.2) A","Cardiomyocytes (ttn.2) V","Cardiomyocytes (proliferating)"))
# select.cells <- CellSelector(DimPlot(CM)+NoLegend())
# CM <- subset(CM,cells = select.cells)
# CM <- StashIdent(CM,save.name = "work.ident")
# CM$work.ident <- as.character(CM$work.ident)

# select.cells <- CellSelector(DimPlot(CM)+NoLegend())
CM <- SetIdent(CM,value = "work.ident")
CM <- SetIdent(CM,cells = WhichCells(final.all.hearts,expression = seurat_clusters == 6), value = "Cardiomyocytes (ttn.2) V")
CM <- StashIdent(CM,save.name = "work.ident2")
#save(CM,file = "final_clustering/CM2.Robj")

# DimPlot(CM,group.by = "work.ident", #"new.ident", 
#         label = T,cols = colors[match(c("green2","green","darkgreen",
#                                         "middlegreen","viridian"),colors$name),]$color,pt.size = 1.2)


## endo
final.all.hearts <- SetIdent(final.all.hearts,value = "final.zoom")
endo <- subset(final.all.hearts,idents = c("Endocardium 1 (V)","Endocardium 1 (A)","Endocardium 2 (V)","Endocardium 2 (A)","Endocardium frzb (A)","Endocardium frzb (V)"))
#select.cells <- CellSelector(DimPlot(endo)+NoLegend())
endo <- subset(endo,cells = select.cells)
endo <- StashIdent(endo,save.name = "work.ident")
endo$work.ident <- as.character(endo$work.ident)
endo$lineage.ident <- as.character(endo$lineage.ident)

# DimPlot(endo,label = T,group.by = "lineage.ident",
#         cols=colors[match(c("lightblue","darkpurple","darkblue"
#                        ),colors$name),]$color,pt.size = 1)



## niche
load(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/niche.noEry.Robj")
# DimPlot(niche,group.by = "work.ident",label = T)
niche@meta.data$work.ident <- 
  plyr::mapvalues(niche@meta.data$work.ident, 
                  from = c("Fibroblast", "Fibroblast (cfd)", "Fibroblast (cxcl12a)", "Fibroblast (spock3)",
                          "Fibroblast-like cells", "Fibroblast (col12a1a)", "Fibroblast (col11a1a)",
                          "Epicardium (Ventricle)", "Epicardium (Atrium)", "Fibroblast (mpeg1.1)",
                          "Fibroblast (nppc)", "Fibroblast (proliferating)", "Perivascular cells"), 
                  to = c("Fibroblasts (const.)", "Fibroblasts (cfd)", "Fibroblasts (cxcl12a)", "Fibroblasts (spock3)",
                         "Valve fibroblasts", "Fibroblasts (col12a1a)", "Fibroblasts (col11a1a)",
                         "Epicardium (V)", "Epicardium (A)", "Fibroblasts (mpeg1.1)",
                         "Fibroblasts (nppc)", "Fibroblasts (prolif.)", "Perivascular cells")  )
 

## immune
load(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/immune.noEry.Robj")

# Plot not in manuscript
# DimPlot(immune,group.by = "work.ident",label = T)+NoLegend()
# DimPlot(immune,label = F,group.by = "work.ident",
#         cols=colors[match(
#           c("yellow","syn-red","green",
#             "orange","darkblue","ocker",
#             "leadwhite","leadwhite","pink",
#             "winered","lightblue","darkgreen",
#             "leadwhite","white2","darkpurple2",
#             "darkorange","leadwhite")
#           ,colors$name),]$color,pt.size = 1)

# Plot UMAPs ####
# Plot 1b and time-resolved subsets
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Hearts_all_UMAP.png",
#     width = 960, height = 960)
DimPlot(final.all.hearts,label = F, pt.size = 0.75,
        cols = col.table$color[match(levels(final.all.hearts$plot.ident2),col.table$set1)]) + NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
# dev.off()
hearts.ctrl <- subset(x = final.all.hearts, subset = time == "Ctrl")
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Hearts_ctrl_subset_UMAP.png",
#     width = 960, height = 960)
DimPlot(hearts.ctrl, label = F, pt.size = 0.75,
        cols = col.table$color[match(levels(hearts.ctrl$plot.ident2),col.table$set1)]) + NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
# dev.off()
hearts.3dpi <- subset(x = final.all.hearts, subset = time == "3dpi")
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Hearts_3dpi_subset_UMAP.png",
#     width = 960, height = 960)
DimPlot(hearts.3dpi, label = F, pt.size = 0.75,
        cols = col.table$color[match(levels(hearts.3dpi$plot.ident2),col.table$set1)]) + NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
# dev.off()
hearts.7dpi <- subset(x = final.all.hearts, subset = time == "7dpi")
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Hearts_7dpi_subset_UMAP.png",
#     width = 960, height = 960)
DimPlot(hearts.7dpi, label = F, pt.size = 0.75,
        cols = col.table$color[match(levels(hearts.7dpi$plot.ident2),col.table$set1)]) + NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
# dev.off()
hearts.30dpi <- subset(x = final.all.hearts, subset = time == "30dpi")
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Hearts_30dpi_subset_UMAP.png",
#     width = 960, height = 960)
DimPlot(hearts.30dpi, label = F, pt.size = 0.75,
        cols = col.table$color[match(levels(hearts.30dpi$plot.ident2),col.table$set1)]) + NoLegend() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
# dev.off()

# 2c
DimPlot(niche,label = F,group.by = "work.ident",
        cols=colors[match(
          #consistent:
          c("syn-lightyellow","darkyellow","creamyellow2",
            "lightpink","lightblue","green2",
            "yellowochre","syn-red","darkpurple",
            "purple","lightpurple","white2",
            "darkblue2")
          ,colors$name),]$color,pt.size = 1)

# Gene expression - UNUSED? ####
# FeaturePlot(final.all.hearts, features = c("vangl2", "prickle1a", "fhl2a"))
# VlnPlot(final.all.hearts, features = c("yap1"),
#         cols = col.table$color[match(levels(final.all.hearts$plot.ident2),col.table$set1)])
# 
# CalculateTimeAverages <- function(seuratobject, features){
#   time_average <-
#     AverageExpression(seuratobject, features = features,
#                       add.ident = "time")
#   time_average_long <- melt(cbind(time_average$RNA,
#                                     data.frame(Gene = rownames(time_average$RNA),
#                                                stringsAsFactors = F)))
#   colnames(time_average_long)[3] <- "Expression"
#   time_average_long$Cell_type <-
#     sapply(as.character(time_average_long$variable),
#            function(x){
#              unlist(strsplit(x, "_"))[1]
#            })
#   time_average_long$Time <-
#     sapply(as.character(time_average_long$variable),
#            function(x){
#              unlist(strsplit(x, "_"))[2]
#            })
#   return(time_average_long)
# }
# 
# wwtr1_time_av <- CalculateTimeAverages(final.all.hearts, features = c("wwtr1"))
# yap1_time_av <- CalculateTimeAverages(final.all.hearts, features = c("yap1"))
# yap1_time_av$Time <- factor(yap1_time_av$Time, levels = c("Ctrl", "3dpi", "7dpi", "30dpi"))
# ggplot(yap1_time_av) +
#   geom_point(aes(x = Time, y = Expression, color = Cell_type)) +
#   facet_wrap(~Cell_type)
# 
# 
# niche_col12_time_average <-
#   AverageExpression(niche, features = c("col12a1a", "col12a1b", "serpine1"),
#                     add.ident = "time")
# niche_time_averages <- melt(cbind(niche_col12_time_average$RNA,
#                                   data.frame(Gene = rownames(niche_col12_time_average$RNA),
#                                              stringsAsFactors = F)))
# colnames(niche_time_averages)[3] <- "Expression"
# niche_time_averages$Cell_type <-
#   sapply(as.character(niche_time_averages$variable),
#          function(x){
#            unlist(strsplit(x, "_"))[1]
#          })
# niche_time_averages$Time <-
#   sapply(as.character(niche_time_averages$variable),
#          function(x){
#            unlist(strsplit(x, "_"))[2]
#          })
# 
# FeaturePlot(final.all.hearts, features = c("nppc", "klf1"))
# DotPlot(final.all.hearts, features = c("nppc", "klf1"))
# VlnPlot(final.all.hearts, features = c("nppc", "klf1"))
# DimPlot(final.all.hearts[final.all.hearts@meta.data$time == "7dpi"])
# sevendpi.hearts <- subset(final.all.hearts, subset = time == "7dpi")
# DimPlot(sevendpi.hearts)
# VlnPlot(sevendpi.hearts, features = c("nppc", "klf1", "klf2", "klf4"))
# 
# FeaturePlot(final.all.hearts, features = c("fli1a", "rspo1", "rspo2", "rspo3", "rspo4"))
# 
# VlnPlot(niche, features = c("col12a1a", "col12a1b", "serpine1"),
#         group.by = "work.ident",
#         #idents = c("Fibroblasts (const.)", "Fibroblasts (col12a1a)", "Fibroblasts (col11a1a)"),
#         split.by = "time")
# niche <- SetIdent(niche, value = "work.ident")
# niche_col12_time_average <-
#   AverageExpression(niche, features = c("col12a1a", "col12a1b", "serpine1"),
#                   add.ident = "time")
# niche_time_averages <- melt(cbind(niche_col12_time_average$RNA,
#                             data.frame(Gene = rownames(niche_col12_time_average$RNA),
#                                        stringsAsFactors = F)))
# colnames(niche_time_averages)[3] <- "Expression"
# niche_time_averages$Cell_type <-
#   sapply(as.character(niche_time_averages$variable),
#        function(x){
#          unlist(strsplit(x, "_"))[1]
#        })
# niche_time_averages$Time <-
#   sapply(as.character(niche_time_averages$variable),
#          function(x){
#            unlist(strsplit(x, "_"))[2]
#          })

final.all.hearts <- SetIdent(final.all.hearts, value = "lineage.ident2")
final.all.hearts@meta.data$time_ident <-
  paste(final.all.hearts@meta.data$lineage.ident2,
        final.all.hearts@meta.data$time, sep = "_")
activated_endo_diffgenes <-
  FindMarkers(final.all.hearts, ident.1 = "Endocardium (V)_7dpi",
              ident.2 = "Endocardium (V)_Ctrl", group.by = "time_ident")
endo7_endoctrl_diffgenes <- activated_endo_diffgenes # 254
nppc_diffgenes <-
  FindMarkers(final.all.hearts, ident.1 = "Fibroblast (nppc)_7dpi",
              ident.2 = "Endocardium (V)_Ctrl", group.by = "time_ident")
nppc_endoctrl_diffgenes <- nppc_diffgenes #790
nppc_activated_endo_diffgenes <-
  FindMarkers(final.all.hearts, ident.1 = "Fibroblast (nppc)_7dpi",
              ident.2 = "Endocardium (V)_7dpi", group.by = "time_ident")
nppc_endo7_diffgenes <- nppc_activated_endo_diffgenes #772

endo_up <- rownames(endo7_endoctrl_diffgenes)[endo7_endoctrl_diffgenes$avg_logFC > 0 &
                                     endo7_endoctrl_diffgenes$p_val < 0.01] #83

nppc_up <- rownames(nppc_endoctrl_diffgenes)[nppc_endoctrl_diffgenes$avg_logFC > 0 &
                                    nppc_endoctrl_diffgenes$p_val < 0.01] #477

intersect(endo_up, nppc_up) #21
write.csv(nppc_endoctrl_diffgenes[nppc_endoctrl_diffgenes$avg_logFC > 0 &
                                    nppc_endoctrl_diffgenes$p_val < 0.01, ],
          "/local/users/Bastiaan/Projects/heart_Bo/Data/nppc_endoctrl_up.csv",
          quote = F, row.names = T)
write.csv(nppc_endoctrl_diffgenes[nppc_endoctrl_diffgenes$avg_logFC < 0 &
                                    nppc_endoctrl_diffgenes$p_val < 0.01, ],
          "/local/users/Bastiaan/Projects/heart_Bo/Data/nppc_endoctrl_down.csv",
          quote = F, row.names = T)
write.csv(endo7_endoctrl_diffgenes[endo7_endoctrl_diffgenes$avg_logFC > 0 &
                                     endo7_endoctrl_diffgenes$p_val < 0.01, ],
          "/local/users/Bastiaan/Projects/heart_Bo/Data/endo7_endoctrl_up.csv",
          quote = F, row.names = T)

endo_emt_sc <- 
  subset(x = final.all.hearts, 
         subset = lineage.ident2 %in% c("Endocardium (V)", "Fibroblast (nppc)"))
endmt_genes <- c("acta2", "cdh5", "pecam1", "cdh1",
                 "vwf", "cdh2", "vim", "icn", "fn1a",
                 "fn1b", "fap", "cnn1b", "tagln",
                 "vcana", "vcanb", "snai1a", "snai1b", "snai2", "snai3")
endmt_genes %in% endo_emt_sc@assays$RNA@data@Dimnames[[1]]

grep("pecam", endo_emt_sc@assays$RNA@data@Dimnames[[1]])
endo_emt_sc@assays$RNA@data@Dimnames[[1]][3127]
endo_emt <- 
  AverageExpression(endo_emt_sc, 
                    features = endmt_genes,
                    add.ident = "time")
endo_emt_time_averages <- melt(cbind(endo_emt$RNA,
                                  data.frame(Gene = rownames(endo_emt$RNA),
                                             stringsAsFactors = F)))
colnames(endo_emt_time_averages)[3] <- "Expression"
endo_emt_time_averages$Cell_type <-
  sapply(as.character(endo_emt_time_averages$variable),
         function(x){
           unlist(strsplit(x, "_"))[1]
         })
endo_emt_time_averages$Time <-
  sapply(as.character(endo_emt_time_averages$variable),
         function(x){
           unlist(strsplit(x, "_"))[2]
         })

# Pie charts 1b ####
pie.numbers <- list()
pie.numbers[["Ctrl"]] <-  summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "Ctrl",]$big.ident) / sum(summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "Ctrl",]$big.ident))
pie.numbers[["3dpi"]] <- summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "3dpi",]$big.ident) / sum(summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "3dpi",]$big.ident))
pie.numbers[["7dpi"]] <- summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "7dpi",]$big.ident) / sum(summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "7dpi",]$big.ident))
pie.numbers[["30dpi"]] <- summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "30dpi",]$big.ident) / sum(summary(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "30dpi",]$big.ident))

type_freqs_df <-
  data.frame(table(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "3dpi", ][, c("orig.ident", "big.ident")]))
dataset_freqs_df <-
  data.frame(table(final.all.hearts@meta.data[final.all.hearts@meta.data$time == "3dpi", ][, c("orig.ident")]))
type_freqs_df <-
  merge(type_freqs_df, dataset_freqs_df,
        by.x = "orig.ident", by.y = "Var1")
type_freqs_df$Ratio_freq <- type_freqs_df$Freq.x/type_freqs_df$Freq.y

pie.numbers <- data.frame(Ctrl = pie.numbers[["Ctrl"]],
                          dpi3 = pie.numbers[["3dpi"]],
                          dpi7 = pie.numbers[["7dpi"]],
                          dpi30 = pie.numbers[["30dpi"]])
pie.numbers <- t(pie.numbers)

pie.numbers <- data.frame(table(final.all.hearts@meta.data$time, final.all.hearts@meta.data$big.ident))
colnames(pie.numbers) <- c("Time","Cell.type","Freq")

for (i in unique(pie.numbers$Time)) {
  pie.numbers$ratio[pie.numbers$Time == i] <-  
    pie.numbers$Freq[pie.numbers$Time == i] / sum(pie.numbers$Freq[pie.numbers$Time == i])
}

pie.numbers <- pie.numbers[complete.cases(pie.numbers),]
pie.numbers$Time <- factor(pie.numbers$Time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

pie.numbers$Cell.type <-  plyr::mapvalues(pie.numbers$Cell.type,
            from =c("Endocardium","Smooth muscle cells","Fibroblasts",
            "Cardiomyocytes","Macrophages","T-cells",
            "Fibroblast-like cells","Bl.ves.EC (apnln)", 
            "Bl.ves.EC (plvapb)", "Bl.ves.EC (lyve1)", "B-cells",
            "Neutrophils","Perivascular cells","Proliferating cells",
            "Myelin cells","Neuronal cells"), 
            to = c("Endocardium","Smooth muscle cells","Fibroblasts",
                   "Cardiomyocytes","Macrophages","T-cells",
                   "Valve Fibroblasts","Blood ves. cells",
                   "Blood ves. cells","Blood ves. cells","B-cells",
                   "Neutrophils","Perivascular cells","Proliferating cells",
                   "Neuronal/Myelin cells","Neuronal/Myelin cells"))
pie.numbers$Cell.type <- 
  factor(pie.numbers$Cell.type,
         levels =c("Fibroblasts","Macrophages","Neutrophils","T-cells","B-cells",
                   "Blood ves. cells", "Endocardium","Cardiomyocytes",
                   "Smooth muscle cells","Perivascular cells","Valve Fibroblasts",
                   "Proliferating cells","Neuronal/Myelin cells"))

pie.plot <- list()
for (i in as.character(unique(pie.numbers$Time))) {
  pie.plot[[i]] <- 
    plot_grid(
      ggplot(pie.numbers[pie.numbers$Time == i,], aes(x=Time,y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_stack(), stat="identity") +
        coord_polar("y", start=0)+
        scale_fill_manual(values = c("#B6A030","#9d3d58","#a981b5","#b18393","#c5adcc","#bd6a47",
                                     "#3e70ab","#779538","#8a6a4b","#525566","#e18c25","grey","grey")   )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ggtitle(i)+
        theme_bw()
    )
}

combine_plots(plotlist = pie.plot)

# Significance test for pie numbers ####
# ALSO CONTAINS BARPLOTS OF CELL TYPE FREQUENCIES
pie.test <- data.frame(table(final.all.hearts@meta.data$orig.ident, 
                             final.all.hearts@meta.data$big.ident))
colnames(pie.test) <- c("Library","Cell.type","Freq")
pie.test <- merge(x=pie.test, y= unique(data.frame(Library=final.all.hearts@meta.data$orig.ident,
                                                      split=final.all.hearts@meta.data$AV,
                                                      time=final.all.hearts@meta.data$time,
                                                      inhib=final.all.hearts@meta.data$inhib)),
                    by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(pie.test$Library))) {
  pie.test$norm.f[pie.test$Library == i] <- 
    sum(pie.test$Freq[pie.test$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- pie.test[pie.test$split == "Wholeheart",]
x <- x[x$inhib == "NULL",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

#CM.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
pie.test.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()

for (i in unique(pie.test.norm$Cell.type)) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(pie.test.norm[pie.test.norm$Cell.type == i,], aes(x=time, y=ratio, fill = Cell.type ) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        #scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all cells") +
        NoLegend()+
        ggtitle(i)
    )
}

# pdf(file = "final_clustering/finalplots/Rplot01.pdf")
combine_plots(plotlist = barplots)
# dev.off()


#test
test.results <- data.frame(dpi3 = rep(NA,length(unique(x$Cell.type))),
                           dpi7 = rep(NA,length(unique(x$Cell.type))),
                           dpi30 = rep(NA,length(unique(x$Cell.type))),
                           row.names = unique(x$Cell.type)) 
for (i in unique(x$Cell.type)) {
  z <- x[x$Cell.type == i,]
  z.ctrl <- z[z$time == "Ctrl",]
  z.3dpi <- z[z$time == "3dpi",]
  z.7dpi <- z[z$time == "7dpi",]
  z.30dpi <- z[z$time == "30dpi",]
  test <- t.test(z.ctrl$Freq,z.3dpi$Freq)
  test.results[i,]$dpi3 <-  test$p.value
  test <- t.test(z.ctrl$Freq,z.7dpi$Freq)
  test.results[i,]$dpi7 <-  test$p.value
  test <- t.test(z.ctrl$Freq,z.30dpi$Freq)
  test.results[i,]$dpi30 <-  test$p.value
}
test.results < 0.05
# write.csv(test.results,file = "final_clustering/t_test_results.csv",quote = F)

# z <- x[x$Cell.type == "Neutrophils",]
# z.ctrl <- z[z$time == "Ctrl",]
# z.3dpi <- z[z$time == "3dpi",]
# z.7dpi <- z[]
# test <- t.test(z.ctrl$Freq,z.3dpi$Freq)
# test$p.value

# A-V split plot 1d ####
## CM
CM.counts <- data.frame(table(CM@meta.data$orig.ident, CM@meta.data$work.ident))
colnames(CM.counts) <- c("Library","Cell.type","Freq")
CM.counts2 <- merge(x=CM.counts, y= unique(data.frame(Library=CM@meta.data$orig.ident,
                                                      split=CM@meta.data$AV,
                                                      time=CM@meta.data$time,
                                                      inhib=CM@meta.data$inhib)),
                    by.x = "Library", by.y = "Library",all.y = F)

y <- CM.counts2[CM.counts2$inhib == "NULL",]
y <- y[y$Cell.type %in% c("Cardiomyocytes V","Cardiomyocytes A","Cardiomyocytes (ttn.2) V","Cardiomyocytes (ttn.2) A"),]
for (i in unique(y$Library)) {
  y$ratio[y$Library == i] <-  y$Freq[y$Library == i] / sum(y$Freq[y$Library == i])
}
y$Cell.type <- factor(y$Cell.type, levels =   (levels(y$Cell.type)) )

plot.av <- 
  ggplot(y[y$Library %in% c("H8a",  "H8v", "Hr6a", "Hr6v", "Hr7a","Hr7v"),], aes(x=Library,y=ratio,fill=Cell.type)) +
  geom_bar(stat = "identity",position = "stack") +
  ylab("% of all cells") + 
  scale_fill_manual(values = col.table[match(levels(as.factor(as.character(y$Cell.type))),col.table$setCM),]$color) +
  #scale_y_continuous(limits = c(0,100)) +
  #scale_x_discrete(labels = c("H5","H6","H7","Hr10","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr11","Hr12","Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr16","Hr17","Hr18","Hr3","Hr4","Hr19","Hr20","Hr21","Hr5","H8a","H8v","Hr6a","Hr6v","Hr7a","Hr7v") ) +
  theme_bw()
plot.av

# z <- y[y$Library %in% c("H8a",  "H8v", "Hr6a", "Hr6v", "Hr7a","Hr7v"),]
# z <- y[y$Cell.type %in% c("Cardiomyocytes (ttn.2) A","Cardiomyocytes (ttn.2) V"),]
# plot.av3 <- 
#   ggplot(z[z$Library %in% c("H8a",  "H8v", "Hr6a", "Hr6v", "Hr7a","Hr7v"),], aes(x=Library,y=norm.ratio,fill=Cell.type)) +
#   geom_bar(stat = "identity",position = "stack") +
#   ylab("% of all cells") + 
#   scale_fill_manual(values = rev(brewer_pal(type = "qual",palette = "Set2")(length(unique(z$Cell.type))))) +
#   #scale_y_continuous(limits = c(0,100)) +
#   #scale_x_discrete(labels = c("H5","H6","H7","Hr10","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr11","Hr12","Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr16","Hr17","Hr18","Hr3","Hr4","Hr19","Hr20","Hr21","Hr5","H8a","H8v","Hr6a","Hr6v","Hr7a","Hr7v") ) +
#   theme_bw()
# plot_grid(plot.av,plot.av2,plot.av3)
# plot.av

## endo
endo.counts <- data.frame(table(endo@meta.data$orig.ident, endo@meta.data$lineage.ident))
colnames(endo.counts) <- c("Library","Cell.type","Freq")
endo.counts2 <- merge(x=endo.counts, y= unique(data.frame(Library=endo@meta.data$orig.ident,
                                                      split=endo@meta.data$AV,
                                                      time=endo@meta.data$time,
                                                      inhib=endo@meta.data$inhib)),
                    by.x = "Library", by.y = "Library",all.y = F)

y <- endo.counts2[endo.counts2$inhib == "NULL",]
y <- y[y$Cell.type %in% c("Endocardium (V)","Endocardium (A)"),]
for (i in unique(y$Library)) {
  y$ratio[y$Library == i] <-  y$Freq[y$Library == i] / sum(y$Freq[y$Library == i])
}
y$Cell.type <- factor(y$Cell.type, levels = unique(y$Cell.type) )
plot.av <- 
  ggplot(y[y$Library %in% c("H8a",  "H8v", "Hr6a", "Hr6v", "Hr7a","Hr7v"),], aes(x=Library,y=ratio,fill=Cell.type)) +
  geom_bar(stat = "identity",position = "stack") +
  ylab("% of all cells") + 
  scale_fill_manual(values =  colors[match(c("lightblue","darkblue"
  ),colors$name),]$color) +
  #scale_y_continuous(limits = c(0,100)) +
  #scale_x_discrete(labels = c("H5","H6","H7","Hr10","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr11","Hr12","Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr16","Hr17","Hr18","Hr3","Hr4","Hr19","Hr20","Hr21","Hr5","H8a","H8v","Hr6a","Hr6v","Hr7a","Hr7v") ) +
  theme_bw()
plot.av


## epicard
epi.counts <- data.frame(table(final.all.hearts@meta.data$orig.ident, final.all.hearts@meta.data$plot.ident))
colnames(epi.counts) <- c("Library","Cell.type","Freq")
epi.counts <- merge(x=epi.counts, y= unique(data.frame(Library=final.all.hearts@meta.data$orig.ident,
                                                          split=final.all.hearts@meta.data$AV,
                                                          time=final.all.hearts@meta.data$time,
                                                          inhib=final.all.hearts@meta.data$inhib)),
                      by.x = "Library", by.y = "Library",all.y = F)

y <- epi.counts[epi.counts$inhib == "NULL",]
y <- y[y$Cell.type %in% c("Epicardium (V)","Epicardium (A)"),]
for (i in unique(y$Library)) {
  y$ratio[y$Library == i] <-  y$Freq[y$Library == i] / sum(y$Freq[y$Library == i])
}
y$Cell.type <- factor(y$Cell.type, levels = unique(y$Cell.type) )
plot.av <- 
  ggplot(y[y$Library %in% c("H8a",  "H8v", "Hr6a", "Hr6v", "Hr7a","Hr7v"),], aes(x=Library,y=ratio,fill=Cell.type)) +
  geom_bar(stat = "identity",position = "stack") +
  ylab("% of all cells") + 
  scale_fill_manual(values = col.table[match(unique(y$Cell.type),col.table$set1),]$color) +
  #scale_y_continuous(limits = c(0,100)) +
  #scale_x_discrete(labels = c("H5","H6","H7","Hr10","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr11","Hr12","Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr16","Hr17","Hr18","Hr3","Hr4","Hr19","Hr20","Hr21","Hr5","H8a","H8v","Hr6a","Hr6v","Hr7a","Hr7v") ) +
  theme_bw()
plot.av


# Time course 2a, S4c, 3a without inhib ####
## CM
CM.counts <- data.frame(table(CM@meta.data$orig.ident, CM@meta.data$lineage.ident))
colnames(CM.counts) <- c("Library","Cell.type","Freq")
CM.counts2 <- merge(x=CM.counts, y= unique(data.frame(Library=CM@meta.data$orig.ident,
                                                      split=CM@meta.data$AV,
                                                      time=CM@meta.data$time,
                                                      inhib=CM@meta.data$inhib)),
                    by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(CM.counts2$Library))) {
  CM.counts2$norm.f[CM.counts2$Library == i] <- sum(CM.counts2$Freq[CM.counts2$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- CM.counts2[CM.counts2$split == "Wholeheart",]
x <- x[x$inhib == "NULL",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

#CM.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
# Fig. 2a
CM.norm.CM <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in c("Cardiomyocytes V", "Cardiomyocytes (ttn.2) V", "Cardiomyocytes (proliferating)")){
  #            "Cardiomyocytes (ttn.2) A", "Cardiomyocytes A")){#unique(CM.norm.CM$Cell.type)) {
  barplots[[i]] <- 
    plot_grid(
      ggplot()+ 
        geom_bar(data = CM.norm.CM[CM.norm.CM$Cell.type == i,], 
                 aes(x=time, y=ratio, fill = Cell.type), 
                 position=position_dodge(), stat="identity") +
        geom_errorbar(data = CM.norm.CM[CM.norm.CM$Cell.type == i,],
                      aes(x=time, ymin=ratio-se, ymax=ratio+se), size = 2,
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        geom_beeswarm(data = x[x$Cell.type == i, ], size = 4,
                   aes(x = time, y = ratio), priority = "ascending", cex = 3) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        # ylab("Normalized ratio to all Cardiomyocytes") +
        labs(x = "", y = "") +
        NoLegend() +
        theme(line = element_line(size = 2),
              panel.grid = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(),
              axis.text.x = element_blank(),
              axis.ticks = element_line(size = 2),
              axis.ticks.length = unit(8,"pt"),
              axis.text.y = element_text(size = 36, family = "Helvetica",
                                         face = "bold"),
              plot.margin = unit(c(10, 60, 0, 0), "pt"))
      #+
        # ggtitle(i)
    )
}
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/CM_dynamics_2a.png",
#     height = 480, width = 1780)
combine_plots(plotlist = barplots, nrow = 1)
# dev.off()

# ggplot(x[x$Cell.type %in% c("Cardiomyocytes V", "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2) V",
#                             "Cardiomyocytes (ttn.2) A", "Cardiomyocytes A"), ]) +
#   geom_bar(aes(x = Library, y = ratio, fill = Cell.type), stat = "identity") +
#   facet_wrap(~time)
# 
# ggplot(x[x$Cell.type %in% c("Cardiomyocytes V", "Cardiomyocytes (proliferating)", "Cardiomyocytes (ttn.2) V",
#                             "Cardiomyocytes (ttn.2) A", "Cardiomyocytes A"), ]) +
#   geom_bar(aes(x = Library, y = Freq, fill = Cell.type), stat = "identity") +
#   facet_wrap(~time)


## niche
fibro.counts <- data.frame(table(niche@meta.data$orig.ident, niche@meta.data$work.ident))
colnames(fibro.counts) <- c("Library","Cell.type","Freq")
fibro.counts <- merge(x=fibro.counts, y= unique(data.frame(Library=niche@meta.data$orig.ident,
                                                      split=niche@meta.data$AV,
                                                      time=niche@meta.data$time,
                                                      inhib=niche@meta.data$inhib)),
                    by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(fibro.counts$Library))) {
  fibro.counts$norm.f[fibro.counts$Library == i] <- sum(fibro.counts$Freq[fibro.counts$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- fibro.counts[fibro.counts$split == "Wholeheart",]
x <- x[x$inhib == "NULL",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

niche.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
niche.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
# Fig S4c
barplots <- list()
for (i in c("Fibroblasts (const.)", "Fibroblasts (prolif.)", "Fibroblasts (col11a1a)", "Fibroblasts (col12a1a)",
            "Fibroblasts (cxcl12a)", "Fibroblasts (cfd)", "Fibroblasts (spock3)", "Fibroblasts (nppc)",
            "Epicardium (V)", "Epicardium (A)")) {
  barplots[[i]] <- 
    ggplot(niche.norm[niche.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), size = 2,
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    geom_beeswarm(data = x[x$Cell.type == i, ], size = 4,
                  aes(x = time, y = ratio), priority = "ascending", cex = 3) +
    scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
    labs(x = "", y = "") +
    ggtitle(i) +
    NoLegend() +
    theme(line = element_line(size = 2),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 32,
                                    family = "Helvetica", face = "bold"),
          axis.line = element_line(),
          axis.text.x = element_text(size = 32, family = "Helvetica",
                                     face = "bold", angle = 45, hjust = 0.5, vjust = 0.6),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(8,"pt"),
          axis.text.y = element_text(size = 32, family = "Helvetica",
                                     face = "bold"),
          plot.margin = unit(c(10, 50, 0, 0), "pt"))
}
for (i in c("Perivascular cells","Valve fibroblasts")) {
  barplots[[i]] <- 
    ggplot(niche.norm.allheart[niche.norm.allheart$Cell.type == i,], aes(x=time, y=norm.ratio, fill = Cell.type) )+ 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se), size = 2,
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    geom_beeswarm(data = x[x$Cell.type == i, ], size = 4,
                  aes(x = time, y = norm.ratio), priority = "ascending", cex = 3) +
    scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
    labs(x = "", y = "") +
    ggtitle(i) +
    NoLegend() +
    theme(line = element_line(size = 2),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 32,
                                    family = "Helvetica", face = "bold"),
          axis.line = element_line(),
          axis.text.x = element_text(size = 32, family = "Helvetica",
                                     face = "bold", angle = 45, hjust = 0.5, vjust = 0.6),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(8,"pt"),
          axis.text.y = element_text(size = 32, family = "Helvetica",
                                     face = "bold"),
          plot.margin = unit(c(10, 50, 0, 0), "pt"))
}

# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Niche_dynamics_S4c.png",
#     height = 1200, width = 1800)
combine_plots(plotlist = barplots, ncol = 4, align = "v", axis = "l")
# dev.off()


for (i in setdiff(unique(niche.norm$Cell.type),c("Perivascular cells","Valve fibroblasts")) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(niche.norm[niche.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all Fibroblasts") +
        NoLegend()+
        ggtitle(i)
    )
}
# combine_plots(plotlist = barplots)

# For 3a, combine selected fibroblast barplots with barplots for perivascular cells/valve fibroblasts
barplots <- list()
for (i in c("Fibroblasts (const.)", "Fibroblasts (col11a1a)",
            "Fibroblasts (nppc)", "Fibroblasts (col12a1a)")) {
  barplots[[i]] <- 
    ggplot(niche.norm[niche.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se), size = 2,
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    geom_beeswarm(data = x[x$Cell.type == i, ], size = 4,
                  aes(x = time, y = ratio), priority = "ascending", cex = 3) +
    scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
    labs(x = "", y = "") +
    ggtitle(i) +
    NoLegend() +
    theme(line = element_line(size = 2),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 42,
                                    family = "Helvetica", face = "bold"),
          axis.line = element_line(),
          axis.text.x = element_text(size = 36, family = "Helvetica",
                                     face = "bold"),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(8,"pt"),
          axis.text.y = element_text(size = 36, family = "Helvetica",
                                     face = "bold"),
          plot.margin = unit(c(10, 50, 0, 0), "pt"))
}
for (i in c("Perivascular cells","Valve fibroblasts")) {
  barplots[[i]] <- 
    ggplot(niche.norm.allheart[niche.norm.allheart$Cell.type == i,], aes(x=time, y=norm.ratio, fill = Cell.type) )+ 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se), size = 2,
                  width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
    geom_beeswarm(data = x[x$Cell.type == i, ], size = 4,
                  aes(x = time, y = norm.ratio), priority = "ascending", cex = 3) +
    scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
    labs(x = "", y = "") +
    ggtitle(i) +
    NoLegend() +
    theme(line = element_line(size = 2),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 42,
                                    family = "Helvetica", face = "bold"),
          axis.line = element_line(),
          axis.text.x = element_text(size = 36, family = "Helvetica",
                                     face = "bold"),
          axis.ticks = element_line(size = 2),
          axis.ticks.length = unit(8,"pt"),
          axis.text.y = element_text(size = 36, family = "Helvetica",
                                     face = "bold"),
          plot.margin = unit(c(10, 50, 0, 0), "pt"))
}

# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Niche_dynamics_3a.png",
#     height = 1900, width = 1080)
combine_plots(plotlist = barplots, ncol = 2, align = "v", axis = "l")
# dev.off()

## immune
immune.counts <- data.frame(table(immune@meta.data$orig.ident, immune@meta.data$work.ident))
colnames(immune.counts) <- c("Library","Cell.type","Freq")
immune.counts <- merge(x=immune.counts, y= unique(data.frame(Library=immune@meta.data$orig.ident,
                                                           split=immune@meta.data$AV,
                                                           time=immune@meta.data$time,
                                                           inhib=immune@meta.data$inhib)),
                      by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(immune.counts$Library))) {
  immune.counts$norm.f[immune.counts$Library == i] <- sum(immune.counts$Freq[immune.counts$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- immune.counts[immune.counts$split == "Wholeheart",]
x <- x[x$inhib == "NULL",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

#immune.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
immune.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in (unique(immune.norm$Cell.type)) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(immune.norm[immune.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all immune cells") +
        NoLegend()+
        ggtitle(i)
    )
}
combine_plots(plotlist = barplots)

#barplots <- barplots[order(names(barplots))]

# pdf(file = "final_clustering/finalplots/all_dynamics.pdf",onefile = T,width = 30,height = 30)
combine_plots(plotlist = barplots)
# dev.off()


# NOT USED? Time course with inhib all ####
## CM
CM.counts <- data.frame(table(CM@meta.data$orig.ident, CM@meta.data$lineage.ident))
colnames(CM.counts) <- c("Library","Cell.type","Freq")
CM.counts2 <- merge(x=CM.counts, y= unique(data.frame(Library=CM@meta.data$orig.ident,
                                                      split=CM@meta.data$AV,
                                                      time=CM@meta.data$is.inhib
                                                      )),
                    by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(CM.counts2$Library))) {
  CM.counts2$norm.f[CM.counts2$Library == i] <- sum(CM.counts2$Freq[CM.counts2$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- CM.counts2[CM.counts2$split == "Wholeheart",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi"))

#CM.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
CM.norm.CM <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in unique(CM.norm.CM$Cell.type)) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(CM.norm.CM[CM.norm.CM$Cell.type == i,], aes(x=time, y=ratio, fill = Cell.type ) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all Cardiomyocytes") +
        NoLegend()+
        ggtitle(i)
    )
}
combine_plots(plotlist = barplots)

## niche
fibro.counts <- data.frame(table(niche@meta.data$orig.ident, niche@meta.data$work.ident))
colnames(fibro.counts) <- c("Library","Cell.type","Freq")
fibro.counts <- merge(x=fibro.counts, y= unique(data.frame(Library=niche@meta.data$orig.ident,
                                                           split=niche@meta.data$AV,
                                                           time=niche@meta.data$is.inhib
                                                           )),
                      by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(fibro.counts$Library))) {
  fibro.counts$norm.f[fibro.counts$Library == i] <- sum(fibro.counts$Freq[fibro.counts$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- fibro.counts[fibro.counts$split == "Wholeheart",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi"))

niche.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
niche.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in setdiff(unique(niche.norm$Cell.type),c("Perivascular cells","Fibroblast-like cells")) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(niche.norm[niche.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all Fibroblasts") +
        NoLegend()+
        ggtitle(i)
    )
}
combine_plots(plotlist = barplots)

for (i in c("Perivascular cells","Fibroblast-like cells")) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(niche.norm.allheart[niche.norm.allheart$Cell.type == i,], aes(x=time, y=norm.ratio, fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all cells") +
        NoLegend()+
        ggtitle(i)
    )
}
#combine_plots(plotlist = barplots)

## immune
immune.counts <- data.frame(table(immune@meta.data$orig.ident, immune@meta.data$work.ident))
colnames(immune.counts) <- c("Library","Cell.type","Freq")
immune.counts <- merge(x=immune.counts, y= unique(data.frame(Library=immune@meta.data$orig.ident,
                                                             split=immune@meta.data$AV,
                                                             time=immune@meta.data$is.inhib
                                                             )),
                       by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(immune.counts$Library))) {
  immune.counts$norm.f[immune.counts$Library == i] <- sum(immune.counts$Freq[immune.counts$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- immune.counts[immune.counts$split == "Wholeheart",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi"))

#immune.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
immune.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
#barplots <- list()
for (i in (unique(immune.norm$Cell.type)) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(immune.norm[immune.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all immune cells") +
        NoLegend()+
        ggtitle(i)
    )
}
# combine_plots(plotlist = barplots)

## all
all.count <- data.frame(table(final.all.hearts@meta.data$orig.ident, final.all.hearts@meta.data$big.ident))
colnames(all.count) <- c("Library","Cell.type","Freq")
all.count <- merge(x=all.count, y= unique(data.frame(Library=final.all.hearts@meta.data$orig.ident,
                                                             split=final.all.hearts@meta.data$AV,
                                                             time=final.all.hearts@meta.data$is.inhib
)),
by.x = "Library", by.y = "Library",all.y = F)
for (i in as.character(unique(all.count$Library))) {
  all.count$norm.f[all.count$Library == i] <- sum(all.count$Freq[all.count$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- all.count[all.count$split == "Wholeheart",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
# for (i in unique(x$Library)) {
#   x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq  
# }

x$norm.ratio <- x$ratio * x$norm.f
x$time <- factor(x$time,levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi"))

#immune.norm.allheart <- summarySE(data = x,measurevar = "norm.ratio" ,groupvars = c("time","Cell.type"))
all.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in (unique(all.norm$Cell.type)) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(all.norm[all.norm$Cell.type == i,], aes(x=time, y=ratio,fill = Cell.type) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        #scale_fill_manual(values = col.table[which(col.table == i,arr.ind=TRUE)[which.max(which(col.table == i,arr.ind=TRUE)[,2]),1],]$color )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all cells") +
        NoLegend()+
        ggtitle(i)
    )
}

barplots <- barplots[order(names(barplots))]

# pdf(file = "all_dynamics_inhib_new.pdf",onefile = T,width = 40,height = 30)
combine_plots(plotlist = barplots)
# dev.off()




# NOT USED? Experimental timecourse####
all.counts <- data.frame(table(final.all.hearts@meta.data$orig.ident, final.all.hearts@meta.data$lineage.ident))
colnames(all.counts) <- c("Library","Cell.type","Freq")
all.counts <- merge(x=all.counts, y= unique(data.frame(Library=final.all.hearts@meta.data$orig.ident,
                                                             split=final.all.hearts@meta.data$AV,
                                                             time=final.all.hearts@meta.data$time,
                                                             inhib=final.all.hearts@meta.data$inhib)),
                       by.x = "Library", by.y = "Library",all.y = F)

x <- all.counts[all.counts$split == "Wholeheart",]
x <- x[x$inhib == "NULL",]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
}
x$time <- factor(x$time,levels = c("Ctrl","3dpi","7dpi","30dpi"))

all.norm <- summarySE(data = x,measurevar = "ratio" ,groupvars = c("time","Cell.type"))
barplots <- list()
for (i in (unique(all.norm$Cell.type)) ) {
  barplots[[i]] <- 
    plot_grid(
      ggplot(all.norm[all.norm$Cell.type == i,], aes(x=time, y=ratio) )+ 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) +
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ylab("Normalized ratio to all cells") +
        NoLegend()+
        ggtitle(i)
    )
}

barplots <- barplots[order(names(barplots))]

# pdf(file = "final_clustering/all_dynamics.pdf",onefile = T,width = 30,height = 30)
combine_plots(plotlist = barplots)
# dev.off()




# Timecourse with inhib 6c, 6d, S26 1&2 ####

## CM
# CM.names <- Idents(all.hearts)
# CM.names <- CM.names[names(CM.names) %in% rownames(CM@meta.data)]
# 
# CM.v2 <- names(CM.names[CM.names == "Cardiomyocytes (V) 2"])
# 
# CM <- SetIdent(CM,value = "lineage.ident")
# CM <- SetIdent(CM,cells = CM.v2,value = "Cardiomyocytes (ttn.2) V")
# CM <- StashIdent(CM, save.name = "new.ident")
# save(CM,file = "/local/Bo/Remap_allhearts/final_clustering/CM.noEry.Robj")

# NEW
CM.counts <- data.frame(table(CM@meta.data$orig.ident, CM@meta.data$lineage.ident))
colnames(CM.counts) <- c("Library","Cell.type","Freq")
CM.counts <- merge(x=CM.counts, y= unique(data.frame(Library=CM@meta.data$orig.ident,
                                                     split=CM@meta.data$AV,
                                                     time=CM@meta.data$time,
                                                     inhib=CM@meta.data$inhib)),
                   by.x = "Library", by.y = "Library",all.y = F)

# normalized to all CM
x <- CM.counts[CM.counts$split %in% c("Wholeheart"),]
x <- x[x$time %in% c("7dpi","3dpi"),]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
  #x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq)
}
x$time <- as.character(x$time)
x[x$Library == "Hr31",]$time <- "inhib 7dpi"
x[x$Library == "Hr32",]$time <- "inhib 7dpi"
x[x$Library == "Hr33",]$time <- "inhib 7dpi"
x[x$Library == "Hr29",]$time <- "inhib 3dpi"
x[x$Library == "Hr34",]$time <- "inhib 3dpi"
x[x$Library == "Hr35",]$time <- "inhib 3dpi"

x <- x[x$time %in% c("7dpi","3dpi","inhib 3dpi","inhib 7dpi"),]
x <- x[x$Cell.type %in% c("Cardiomyocytes (ttn.2) V"), ]#,"Cardiomyocytes (proliferating)"),]

x.numbers <- summarySE(data = x,measurevar = c("ratio"), groupvars = c("time","Cell.type") )
x.numbers$time <- factor(x.numbers$time, levels = c("3dpi","inhib 3dpi",
                                                    "7dpi","inhib 7dpi"))
# 6c
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/CM_dediff_inhib_dynamics_6c.png",
#     height = 560, width = 400)
ggplot(x.numbers, aes(x=time,y=ratio, alpha=time, fill = Cell.type)) +
  geom_bar(stat = "identity",width = 1) +#position = "dodge") +
  scale_alpha_manual(values = c(1,0.7,1,0.7))+
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),size = 2,
                width=.2,position=position_dodge(.9))  +
  geom_beeswarm(data = x, size = 4,
                aes(x = time, y = ratio), priority = "ascending", 
                cex = 3) +
  scale_fill_manual(values = rev(col.table[col.table$setCM %in% unique(as.character(x.numbers$Cell.type)),]$color) ) +
  labs(x = "", y = "") +
  NoLegend() +
  theme(line = element_line(size = 2),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(8,"pt"),
        axis.text.y = element_text(size = 36, family = "Helvetica",
                                   face = "bold"))
# dev.off()

## niche
fibro.counts <- data.frame(table(niche@meta.data$orig.ident, niche@meta.data$work.ident),
                           stringsAsFactors = F)
colnames(fibro.counts) <- c("Library","Cell.type","Freq")
fibro.counts <- merge(x=fibro.counts, y= unique(data.frame(Library=niche@meta.data$orig.ident,
                                                     split=niche@meta.data$AV,
                                                     time=niche@meta.data$time,
                                                     inhib=niche@meta.data$inhib)),
                   by.x = "Library", by.y = "Library",all.y = F)

for (i in as.character(unique(fibro.counts$Library))) {
  fibro.counts$norm.f[fibro.counts$Library == i] <- sum(fibro.counts$Freq[fibro.counts$Library == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

x <- fibro.counts[fibro.counts$split %in% c("Wholeheart"),]
x <- x[x$time %in% c("7dpi","3dpi"),]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
  #x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq)
}
x$time <- as.character(x$time)
x[x$Library == "Hr31",]$time <- "inhib 7dpi"
x[x$Library == "Hr32",]$time <- "inhib 7dpi"
x[x$Library == "Hr33",]$time <- "inhib 7dpi"
x[x$Library == "Hr29",]$time <- "inhib 3dpi"
x[x$Library == "Hr34",]$time <- "inhib 3dpi"
x[x$Library == "Hr35",]$time <- "inhib 3dpi"

x$norm.ratio <- x$ratio * x$norm.f
x1 <- x[x$Cell.type %in% c("Perivascular cells", "Valve fibroblasts",
                           "Fibroblasts (nppc)", "Fibroblasts (spock3)"),]
x1$Cell.type <- factor(x1$Cell.type,
                       levels = c("Perivascular cells", "Valve fibroblasts",
                                  "Fibroblasts (nppc)", "Fibroblasts (spock3)"))
x1$time <- factor(x1$time, levels = c("3dpi","inhib 3dpi",
                                      "7dpi","inhib 7dpi"))

niche.norm.allheart <- summarySE(data = x1,measurevar = c("norm.ratio"), 
                                 groupvars = c("time","Cell.type") )
niche.norm.allheart$time <- factor(niche.norm.allheart$time, levels = c("3dpi","inhib 3dpi",
                                                    "7dpi","inhib 7dpi"))

# 6d
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Periv_and_endo_niche_inhib_dynamics_6d.png",
#     height = 450, width = 1065)
ggplot(niche.norm.allheart, 
       aes(x=Cell.type,y=norm.ratio)) +
  geom_bar(aes(alpha=time, fill = Cell.type),
           stat = "identity",position=position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se, group = time),
                width=.2,position=position_dodge(width = 0.9)) +
  geom_beeswarm(data = x1, 
             aes(x = Cell.type, y = norm.ratio, group = time),
             size = 2, dodge.width = 0.9) +
  scale_alpha_manual(values = c(1,0.7,1,0.7))+
  scale_fill_manual(values = setNames(col.table$color[!is.na(col.table$setFibro)], col.table$setFibro[!is.na(col.table$setFibro)])) + 
  labs(x = "", y = "") +
  NoLegend() +
  theme(line = element_line(size = 2),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size = 2),
        axis.ticks.length = unit(8,"pt"),
        axis.text.y = element_text(size = 36, family = "Helvetica",
                                   face = "bold"))
# dev.off()

x2 <- x[!x$Cell.type %in% c("Valve fibroblasts","Perivascular cells"),]
x2$Cell.type <- factor(x2$Cell.type,
                       levels = col.table[col.table$setFibro %in% unique(as.character(niche.norm$Cell.type)),]$setFibro)
x2$time <- factor(x2$time, levels = c("3dpi","inhib 3dpi",
                                      "7dpi","inhib 7dpi"))
niche.norm <- summarySE(data = x2,measurevar = c("ratio"), groupvars = c("time","Cell.type") )
niche.norm$time <- factor(niche.norm$time, levels = c("3dpi","inhib 3dpi",
                                                                        "7dpi","inhib 7dpi"))
niche.norm$Cell.type <- factor(niche.norm$Cell.type ,levels = col.table[col.table$setFibro %in% unique(as.character(niche.norm$Cell.type)),]$setFibro)

# S26 1
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Niche_ratio_fib_inhib_dynamics_S26_1.png",
#     height = 448, width = 2220)
ggplot(niche.norm, aes(x=Cell.type,y=ratio)) +
  geom_bar(aes(alpha=time, fill = Cell.type),
           stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se, group = time),
                width=.2,position=position_dodge(.9)) +
  geom_beeswarm(data = x2, cex = 0.5,
                aes(x = Cell.type, y = ratio, group = time),
                size = 2, dodge.width = 0.9) +
  scale_alpha_manual(values = c(1,0.7,1,0.7))+
  scale_fill_manual(values = setNames(col.table$color[!is.na(col.table$setFibro)], col.table$setFibro[!is.na(col.table$setFibro)])) + 
  labs(x = "", y = "ratio to all fibroblasts") +
  NoLegend() +
  theme(line = element_line(size = 2),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(8,"pt"),
        axis.title.y = element_text(size = 37, family = "Helvetica",
                                    face = "bold"),
        axis.text.y = element_text(size = 36, family = "Helvetica",
                                   face = "bold"))
# dev.off()

x$Cell.type <- factor(as.character(x$Cell.type),
                       levels = col.table[col.table$setFibro %in% unique(as.character(x$Cell.type)),]$setFibro)
x$time <- factor(x$time, levels = c("3dpi","inhib 3dpi",
                                      "7dpi","inhib 7dpi"))

niche.norm2 <- summarySE(data = x,measurevar = c("norm.ratio"), groupvars = c("time","Cell.type") )

# S26 2
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/Niche_ratio_all_inhib_dynamics_S26_2.png",
#     height = 448, width = 2220)
ggplot(niche.norm2, aes(x=Cell.type,y=norm.ratio)) +
  geom_bar(aes(alpha=time, fill = Cell.type), 
           stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se, group = time),
                width=.2, position=position_dodge(.9)) +
  geom_beeswarm(data = x, cex = 0.5,
                aes(x = Cell.type, y = norm.ratio, group = time),
                size = 2, dodge.width = 0.9) +
  scale_alpha_manual(values = c(1,0.7,1,0.7))+
  scale_fill_manual(values = setNames(col.table$color[!is.na(col.table$setFibro)], col.table$setFibro[!is.na(col.table$setFibro)])) + 
  labs(x = "", y = "ratio to all cells") +
  NoLegend() +
  theme(line = element_line(size = 2),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(8,"pt"),
        axis.title.y = element_text(size = 37, family = "Helvetica",
                                    face = "bold"),
        axis.text.y = element_text(size = 36, family = "Helvetica",
                                   face = "bold"))
# dev.off()

# Endo ##

endo.counts <- data.frame(table(final.all.hearts@meta.data$orig.ident, final.all.hearts@meta.data$lineage.ident2))
colnames(endo.counts) <- c("Library","Cell.type","Freq")
endo.counts <- merge(x=endo.counts, y= unique(data.frame(Library=final.all.hearts@meta.data$orig.ident,
                                                     split=final.all.hearts@meta.data$AV,
                                                     time=final.all.hearts@meta.data$time,
                                                     inhib=final.all.hearts@meta.data$inhib)),
                   by.x = "Library", by.y = "Library",all.y = F)

# normalized to all CM
x <- endo.counts[endo.counts$split %in% c("Wholeheart"),]
for (i in unique(x$Library)) {
  x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / sum(x$Freq[x$Library == i])
  #x$ratio[x$Library == i] <-  x$Freq[x$Library == i] / (x[x$Library == i,][x[x$Library == i,]$Cell.type == "CM_A",]$Freq)
}
x$time <- as.character(x$time)
x[x$Library == "Hr31",]$time <- "inhib 7dpi"
x[x$Library == "Hr32",]$time <- "inhib 7dpi"
x[x$Library == "Hr33",]$time <- "inhib 7dpi"
x[x$Library == "Hr29",]$time <- "inhib 3dpi"
x[x$Library == "Hr34",]$time <- "inhib 3dpi"
x[x$Library == "Hr35",]$time <- "inhib 3dpi"

x <- x[x$Cell.type %in% c("Endocardium (frzb)","Endocardium (A)","Endocardium (V)"),]

x.numbers <- summarySE(data = x,measurevar = c("ratio"), groupvars = c("time","Cell.type") )
x.numbers$time <- factor(x.numbers$time, levels = c("Ctrl","3dpi","inhib 3dpi",
                                                    "7dpi","inhib 7dpi","30dpi"))
ggplot(x.numbers, aes(x=Cell.type,y=ratio, alpha=time, fill = Cell.type)) +
  geom_bar(stat = "identity",position = "dodge") +
  #scale_fill_manual(values = as.character(col.table$Color ) ) +
  scale_alpha_manual(values = c(1,1,0.7,1,0.7,1))+
  geom_errorbar(aes(ymin=ratio-se, ymax=ratio+se),width=.2,position=position_dodge(.9)) +
  scale_fill_manual(values = rev(col.table[col.table$set1 %in% unique(as.character(x.numbers$Cell.type)),]$color) ) +
  #scale_y_continuous(limits = c(0,2.8)) +
  ylab("ratio to all fibroblasts") +
  theme_bw()





# NOT USED? Ecm analysis ####
anno <-  rbind(
  data.frame (gene = c(rownames(final.all.hearts)[grep("^col[0-9]",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fn1",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fbln",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("frem",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^matn",rownames(final.all.hearts))],
                       "elnb"                                                 ),Function = "Structure",color = "#1b9e77"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^lama",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^lamb",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^lamc",rownames(final.all.hearts))]),Function = "Structure",color = "#1b9e77"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^mmp",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^cts",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^prss",rownames(final.all.hearts))]),Function = "ECM breakdown",color = "#d95f02"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^itg",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^nid",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("flrt",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fln",rownames(final.all.hearts))]),Function =  "ECM interaction",color = "#7570b3"),
  data.frame (gene = c("tnc","dcn","tnn"),Function =  "ECM interaction",color = "#7570b3")
)
rownames(anno) <- anno$gene
anno_col <- as.character(unique(anno[,3]))
names(anno_col) <- as.character(unique(anno[,2]))
#saveRDS(anno,file = "final_clustering/ecm.genelist.rds")

final.all.hearts@meta.data$big.ident <- final.all.hearts@meta.data$work.ident3
final.all.hearts@meta.data$big.ident <- plyr::mapvalues(final.all.hearts@meta.data$big.ident, 
   from =c("Cardiomyocytes A","Cardiomyocytes V","Cardiomyocytes (ttn.2) A",
           "Endocardium 1 (V)","Endocardium 1 (A)","Smooth muscle cells",
           "Cardiomyocytes (ttn.2) V","Cardiomyocytes (proliferating)","B-cells",
           "Fibroblasts","T-cells","Macrophages",
           "Proliferating cells","Neutrophils","Endocardium 2 (A)",
           "Fibroblast-like cells","Endocardium frzb (A)","Endocardium frzb (V)",
           "Endocardium 2 (V)","Myelin cells","Neuronal cells",
           "Perivascular cells","Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)",
           "Bl.ves.EC (lyve1)"), 
   to = c("Cardiomyocytes","Cardiomyocytes","Cardiomyocytes",
          "Endocardium","Endocardium","Smooth muscle cells",
          "Cardiomyocytes","Cardiomyocytes","B-cells",
          "Fibroblasts","T-cells","Macrophages",
          "Proliferating cells","Neutrophils","Endocardium",
          "Fibroblast-like cells","Endocardium","Endocardium",
          "Endocardium","Myelin cells","Neuronal cells",
          "Perivascular cells","Bl.ves.EC","Bl.ves.EC",
          "Bl.ves.EC"))

final.all.hearts <- SetIdent(final.all.hearts,value = "big.ident")

ecm <- DotPlot(final.all.hearts, features = unique(c(c("adcyap1a","wnt2ba","clu","clu","fshb","fshb","hpx","hpx","zgc:172271","zgc:172271","ngfb","ngfb","ngfb","ngfb","ngfb","ngfb","tfa","tfa","uts2b","adcyap1b","notum1a","htra1a","ins","pyyb","pyyb","npy","lamc1","bmp15","gdf3","wnt6a","apoeb","apoeb","apoa4b.1","gcgb","gcgb","wnt2","ccl38a.4","ccl38a.4","ccl38.1","ccl38.1","ccl38.6","ccl38.6","nog2","loxl2b","loxl2b","gnrh2","wnt7aa","inha","zgc:195023","nppb","nppa","pyya","epoa","epoa","spx","bglap","bmp3","plat","plat","nppc","pros1","shha","ghrh","agr2","insl5a","wnt4a","gdf10b","defbl1","gcga","vipb","scg3","fibinb","galn","urp1","igf1","apoa4b.3","wnt16"))))

ecm <- DotPlot(final.all.hearts, features = as.character(anno$gene),group.by = "big.ident")
ecm <- ecm$data
ecm <- ecm[,c(1,2,3,4)]
#ecm$value <- ecm$avg.exp*ecm$pct.exp*100
ecm <- acast(ecm, id~features.plot, value.var = "avg.exp")
ecm <- ecm[!rownames(ecm) %in% c("Myelin cells","Neuronal cells","Proliferating cells"),]
#ecm.tidy <- ecm.tidy[1:20,]
ecm <- scale(ecm)

ecm <- ecm[c("Fibroblasts","Perivascular cells","Fibroblast-like cells",
              "Smooth muscle cells","Bl.ves.EC",
              "Endocardium","Cardiomyocytes","Macrophages","Neutrophils","T-cells","B-cells"),]

pheatmap(ecm,cluster_rows = T ,cluster_cols = T)
                 #annotation_col = anno[,c(2),drop = F],
                 #annotation_row = ct[,c(2),drop = F],
                 
                 #annotation_colors =  list(Function = anno_col)[1])

structure.genes <- anno[anno$Function == "Structure",]$gene
structure.Sum <- ecm[,colnames(ecm) %in% structure.genes]
structure.Sum <- structure.Sum > 0
structure.Sum <- rowSums(structure.Sum)
structure.Sum <- as.data.frame(structure.Sum)
structure.Sum$names <- rownames(structure.Sum)
structure.Sum$names <- factor(structure.Sum$names, levels = c("Fibroblasts","Perivascular cells","Fibroblast-like cells",
             "Smooth muscle cells","Bl.ves.EC",
             "Endocardium","Cardiomyocytes","Macrophages","Neutrophils","T-cells","B-cells"))
ggplot(structure.Sum,aes(x = names,y = structure.Sum,fill = names))+
  geom_bar(stat = "identity")+
  theme_bw() +
  RotatedAxis()+
  ylab("Number of ECM structure genes") +
  NoLegend() +
  scale_fill_manual(values = c("#B6A030" ,"#525566" ,"#e18c25", "#8a6a4b" ,"#a77745" , "#3e70ab","#3f5346",     
                   "#9d3d58", "#a981b5", "#b18393", "#c5adcc"))

breakdown.genes <- anno[anno$Function == "ECM breakdown",]$gene
breakdown.Sum <- ecm[,colnames(ecm) %in% breakdown.genes]
breakdown.Sum <- breakdown.Sum > 0
breakdown.Sum <- rowSums(breakdown.Sum)
breakdown.Sum <- as.data.frame(breakdown.Sum)
breakdown.Sum$names <- rownames(breakdown.Sum)
breakdown.Sum$names <- factor(breakdown.Sum$names, levels = c("Fibroblasts","Perivascular cells","Fibroblast-like cells",
                                                              "Smooth muscle cells","Bl.ves.EC (plvapb)","Bl.ves.EC (lyve1)","Bl.ves.EC (apnln)",
                                                              "Endocardium","Cardiomyocytes","Macrophages","Neutrophils","T-cells","B-cells"))
ggplot(breakdown.Sum,aes(x = names,y = breakdown.Sum,fill = names))+
  geom_bar(stat = "identity")+
  theme_bw() +
  RotatedAxis()+
  ylab("Number of ECM breakdown genes") +
  NoLegend() +
  scale_fill_manual(values = c("#B6A030" ,"#525566" ,"#e18c25", "#8a6a4b" ,"#a77745" ,"#B68239", "#bd6a47", "#7b7f51","#3f5346",     
                               "#9d3d58", "#a981b5", "#b18393", "#c5adcc"))


# ECM genes within niche 2d ####
anno <-  rbind(
  data.frame (gene = c(rownames(final.all.hearts)[grep("^col[0-9]",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fn1",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fbln",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("frem",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^matn",rownames(final.all.hearts))],
                       "elnb"                                                 ),Function = "Structure",color = "#1b9e77"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^lama",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^lamb",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^lamc",rownames(final.all.hearts))]),Function = "Structure",color = "#1b9e77"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^mmp",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^cts",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^prss",rownames(final.all.hearts))]),Function = "ECM breakdown",color = "#d95f02"),
  data.frame (gene = c(rownames(final.all.hearts)[grep("^itg",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^nid",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("flrt",rownames(final.all.hearts))],
                       rownames(final.all.hearts)[grep("^fln",rownames(final.all.hearts))]),Function =  "ECM interaction",color = "#7570b3"),
  data.frame (gene = c("tnc","dcn","tnn"),Function =  "ECM interaction",color = "#7570b3")
)
rownames(anno) <- anno$gene
anno_col <- as.character(unique(anno[,3]))
names(anno_col) <- as.character(unique(anno[,2]))

ecm.niche <- DotPlot(niche, features = unique(c(c("adcyap1a","wnt2ba","clu","clu","fshb","fshb","hpx","hpx","zgc:172271","zgc:172271","ngfb","ngfb","ngfb","ngfb","ngfb","ngfb","tfa","tfa","uts2b","adcyap1b","notum1a","htra1a","ins","pyyb","pyyb","npy","lamc1","bmp15","gdf3","wnt6a","apoeb","apoeb","apoa4b.1","gcgb","gcgb","wnt2","ccl38a.4","ccl38a.4","ccl38.1","ccl38.1","ccl38.6","ccl38.6","nog2","loxl2b","loxl2b","gnrh2","wnt7aa","inha","zgc:195023","nppb","nppa","pyya","epoa","epoa","spx","bglap","bmp3","plat","plat","nppc","pros1","shha","ghrh","agr2","insl5a","wnt4a","gdf10b","defbl1","gcga","vipb","scg3","fibinb","galn","urp1","igf1","apoa4b.3","wnt16"))))

ecm.niche <- DotPlot(niche, features = setdiff(as.character(anno$gene),c("prss60.3","prss60.2","mmp24","mmp19","col9a1b","mmp28","mmp25b","mmp16b")),group.by = "work.ident")
ecm.niche <- ecm.niche$data
ecm.niche <- ecm.niche[,c(1,2,3,4)]
#ecm$value <- ecm$avg.exp*ecm$pct.exp*100
ecm.niche <- acast(ecm.niche, id~features.plot, value.var = "avg.exp")
#ecm.niche <- ecm.niche[!rownames(ecm.niche) %in% c("Fibroblast (spock3)","Fibroblast (mpeg1.1)","Epicardium (Atrium)","Epicardium (Ventricle)"),]
#ecm.tidy <- ecm.tidy[1:20,]
ecm.niche <- scale(ecm.niche)

# 2d
pheatmap(ecm.niche,cluster_rows = T ,cluster_cols = T, 
                 annotation_col = anno[,c(2),drop = F],
                 #annotation_row = ct[,c(2),drop = F],
                 annotation_colors =  list(Function = anno_col)[1])

# S5: How much of niche clustering is determined by ECM genes? ####
# Count number of differentially expressed ECM and non-ECM genes
marker.niche <- 
  read.csv(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/marker.niche.csv",
           stringsAsFactors = F)
marker.niche.ecm <- merge(marker.niche, anno[, c("gene", "Function")], all.x = T)
marker.niche.ecm$ECM <- ifelse(is.na(marker.niche.ecm$Function), F, T)
ECM.gene.fraction <- data.frame(acast(data.frame(table(marker.niche.ecm[, c("cluster", "ECM")])), cluster ~ ECM, value.var = "Freq"))
colnames(ECM.gene.fraction)[colnames(ECM.gene.fraction) == "TRUE."] <- "ECM"
colnames(ECM.gene.fraction)[colnames(ECM.gene.fraction) == "FALSE."] <- "Other"
ECM.gene.fraction$Fraction_ECM <- ECM.gene.fraction$ECM/(ECM.gene.fraction$ECM + ECM.gene.fraction$Other)
range(ECM.gene.fraction$Fraction_ECM)
# Between 0 and 20% of differentially expressed genes are ECM.

# Calculate portion of differentially expressed genes by ECM and non-ECM
niche <- SetIdent(niche,value = "work.ident")
niche_averages <- (AverageExpression(niche))$RNA
niche_averages_ECM <- niche_averages[rownames(niche_averages) %in% anno$gene, ]
niche_averages_Other <- niche_averages[!(rownames(niche_averages) %in% anno$gene), ]

niche_expressions <- merge(data.frame(ECM = colSums(niche_averages_ECM)),
                           data.frame(Other = colSums(niche_averages_Other)),
                           by = "row.names")
niche_expressions$Fraction_ECM_expressed <-
  niche_expressions$ECM/(niche_expressions$ECM + niche_expressions$Other)
range(niche_expressions$Fraction_ECM_expressed)
# Between 0.4 and 6% of expression is ECM-related.

# Remove ECM genes and recluster
niche_test <- niche
niche_test <- NormalizeData(object = niche_test)
niche_test <- FindVariableFeatures(object = niche_test, selection.method = "mean.var.plot",mean.cutoff = c(0.02,3.5) ,dispersion.cutoff = c(0.6,10))
niche_test <- ScaleData(object = niche_test,vars.to.regress = c("orig.ident","nFeature_RNA","nCount_RNA","percent.mito"))
niche_test <- RunPCA(object = niche_test, features = VariableFeatures(object = niche_test),npcs = 100)
niche_test <- FindNeighbors(object = niche_test,force.recalc = T, dims = 1:20)
niche_test <- FindClusters(object = niche_test,resolution = 1)
niche_test <- RunUMAP(niche_test, dims = 1:20)
DimPlot(object = niche_test, reduction = "umap",group.by = "RNA_snn_res.1",label = TRUE,pt.size = 1) + NoLegend()
DimPlot(object = niche_test, reduction = "umap", group.by = "work.ident")
# Full reproduction of niche object clustering + UMAP.

niche_no_ECM <- subset(x = niche, features = setdiff(niche@assays$RNA@counts@Dimnames[[1]], anno$gene))
niche_no_ECM <- NormalizeData(object = niche_no_ECM)
niche_no_ECM <- FindVariableFeatures(object = niche_no_ECM, selection.method = "mean.var.plot",mean.cutoff = c(0.02,3.5) ,dispersion.cutoff = c(0.6,10))
niche_no_ECM <- ScaleData(object = niche_no_ECM,vars.to.regress = c("orig.ident","nFeature_RNA","nCount_RNA","percent.mito"))
niche_no_ECM <- RunPCA(object = niche_no_ECM, features = VariableFeatures(object = niche_no_ECM),npcs = 100)
niche_no_ECM <- FindNeighbors(object = niche_no_ECM,force.recalc = T, dims = 1:20)
niche_no_ECM <- FindClusters(object = niche_no_ECM,resolution = 1)
niche_no_ECM <- RunUMAP(niche_no_ECM, dims = 1:20)
# S5a
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/No_ECM_niche_clustering.png")
DimPlot(object = niche_no_ECM, reduction = "umap", group.by = "RNA_snn_res.1",label = TRUE,pt.size = 1) + #NoLegend() +
  scale_x_reverse() + scale_y_reverse()
# dev.off()

# Clusters and UMAP look very similar
niche_clustering <- niche@meta.data[, "work.ident", drop = F]
colnames(niche_clustering) <- "Clusters_all"
niche_clustering$Cell <- rownames(niche_clustering)
niche_no_ECM_clustering <- niche_no_ECM@meta.data[, "RNA_snn_res.1", drop = F]
colnames(niche_no_ECM_clustering) <- "Clusters_no_ECM"
niche_no_ECM_clustering$Cell <- rownames(niche_no_ECM_clustering)

ECM_or_not_cluster_compare <- merge(niche_clustering, niche_no_ECM_clustering)
# Majority voting of the new clusters:
cluster_voting <- data.frame(Clusters_no_ECM = unique(ECM_or_not_cluster_compare$Clusters_no_ECM),
                             Fibro_subtype_no_ECM = "", stringsAsFactors = F)
for(i in 1:nrow(cluster_voting)){
  cluster_assignment_table <- 
    table(ECM_or_not_cluster_compare$Clusters_all[ECM_or_not_cluster_compare$Clusters_no_ECM == cluster_voting$Clusters_no_ECM[i]])
  majority_vote <- names(which(cluster_assignment_table == max(cluster_assignment_table)))
  if(length(majority_vote) > 1){
  cluster_voting$Fibro_subtype_no_ECM[i] <- paste(majority_vote, collapse = ",")
  }else{
    cluster_voting$Fibro_subtype_no_ECM[i] <- majority_vote
  }
}
test_merge <- niche_no_ECM@meta.data
test_merge$Cell_name <- row.names(test_merge)
test_merge <- merge(test_merge, cluster_voting, by.x = "RNA_snn_res.1", by.y = "Clusters_no_ECM")
rownames(test_merge) <- test_merge$Cell_name
test_merge <- test_merge[, -which(colnames(test_merge)=="Cell_name")]
test_merge$Fibro_subtype_no_ECM <- factor(test_merge$Fibro_subtype_no_ECM, levels =
                                            levels(niche@meta.data$work.ident))
# 
niche_no_ECM@meta.data <- test_merge
niche_colors <- colors[match(
  c("syn-lightyellow","darkyellow","creamyellow2",
    "lightpink","lightblue","green2",
    "yellowochre","syn-red","darkpurple",
    "purple","lightpurple","white2",
    "darkblue2")
  ,colors$name),]$color
# png("/local/users/Bastiaan/Projects/heart_Bo/Images/No_ECM_niche_clustering_majority_voting.png")
DimPlot(object = niche_no_ECM, reduction = "umap", group.by = "Fibro_subtype_no_ECM",label = FALSE,
        cols=niche_colors,pt.size = 1) + NoLegend() +
  scale_x_reverse() + scale_y_reverse()
# dev.off()

ECM_or_not_cluster_compare <- merge(ECM_or_not_cluster_compare, cluster_voting)
ECM_or_not_cluster_compare_square <- 
  acast(data.frame(table(ECM_or_not_cluster_compare[, c("Clusters_all", "Fibro_subtype_no_ECM")])),
        Clusters_all ~ Fibro_subtype_no_ECM, value.var = "Freq") # Clusters_all are rows, Clusters_no_ECM are columns
ECM_or_not_cluster_compare_freq <- t(t(ECM_or_not_cluster_compare_square)/colSums(ECM_or_not_cluster_compare_square))
colSums(ECM_or_not_cluster_compare_freq)
ECM_or_not_cluster_compare_freq <- ECM_or_not_cluster_compare_freq[colnames(ECM_or_not_cluster_compare_freq), ]

# S5b
# pdf("/local/users/Bastiaan/Projects/heart_Bo/Images/No_ECM_cluster_voting.pdf")
pheatmap(ECM_or_not_cluster_compare_freq, cluster_rows = F, cluster_cols = F)
# dev.off()

apply(ECM_or_not_cluster_compare_freq, 2, max)

# NOT USED? Secreted selected ####
secreted <- DotPlot(final.all.hearts,features = c("ghrh","bglap","pyya","nppb","nppc","wnt4a","adcyap1a","epoa","pros1","notum1a","urp1","hpx","apoa4b.1","ccl38a.4","htra1a","igf1"),group.by = "lineage.ident2")
secreted <- DotPlot(final.all.hearts, features = unique(c(c("adcyap1a","wnt2ba","clu","clu","fshb","fshb","hpx","hpx","zgc:172271","zgc:172271","ngfb","ngfb","ngfb","ngfb","ngfb","ngfb","tfa","tfa","uts2b","adcyap1b","notum1a","htra1a","ins","pyyb","pyyb","npy","lamc1","bmp15","gdf3","wnt6a","apoeb","apoeb","apoa4b.1","gcgb","gcgb","wnt2","ccl38a.4","ccl38a.4","ccl38.1","ccl38.1","ccl38.6","ccl38.6","nog2","loxl2b","loxl2b","gnrh2","wnt7aa","inha","zgc:195023","nppb","nppa","pyya","epoa","epoa","spx","bglap","bmp3","plat","plat","nppc","pros1","shha","ghrh","agr2","insl5a","wnt4a","gdf10b","defbl1","gcga","vipb","scg3","fibinb","galn","urp1","igf1","apoa4b.3","wnt16"))),group.by = "lineage.ident2")

secreted <- secreted$data
secreted <- secreted[,c(1,2,3,4)]
secreted <- acast(secreted, id~features.plot, value.var = "avg.exp")
secreted <- secreted[!rownames(secreted) %in% c("Myelin cells","Neuronal cells","Dead cells"),]
secreted <- scale(secreted)

pheatmap(secreted,cluster_rows = T ,cluster_cols = T)
#annotation_col = anno[,c(2),drop = F],
#annotation_row = ct[,c(2),drop = F],

# NOT USED? Secreted non-curated ####
secretome <- read.csv("/local/Bo/Remap_allhearts/final_clustering/secreted/resultsTable_noncur_genes.csv")
secreted <- DotPlot(final.all.hearts, features = as.character(unique(secretome$V2)) , group.by = "trashid",split.by  = "inhib")
                    
secreted <- secreted$data
secreted <- secreted[,c(1,2,3,4)]
secreted <- acast(secreted, id~features.plot, value.var = "avg.exp")
secreted <- secreted[!rownames(secreted) %in% c("Neuronal cells_NULL","Myelin cells_NULL","Dead cells_NULL",rownames(secreted)[grep("IWR1$",rownames(secreted))]) ,]
secreted.test <- secreted[,!colSums(secreted) == 0]
rownames(secreted.test) <- as.character(sapply(rownames(secreted.test),  function(x) unlist(strsplit(x,"_"))[1] ))
secreted.test<- scale(secreted.test)
#wnt.all.2 <- scale(wnt.all.2)
pheatmap(secreted.test)

# NOT USED? Secreted time ####

sec.time <- DotPlot(final.all.hearts,group.by = "time",features = c("wnt4a","clu","igf1","loxl2b","bglap","ghrh","lamc1","tfa","ngfb","gcga","plat","nog2","insl5a","galn","gnrh2","fibinb"))
sec.time <- sec.time$data
sec.time <- sec.time[,c(1,2,3,4)]
sec.time <- acast(sec.time, id~features.plot, value.var = "avg.exp")
sec.time <- scale(sec.time)

pheatmap(sec.time,cluster_rows = F ,cluster_cols = T)

sec.time <- DotPlot(final.all.hearts,group.by = "time",features = as.character(unique(secretome$V2)))
sec.time <- sec.time$data
sec.time <- sec.time[,c(1,2,3,4)]
sec.time <- acast(sec.time, id~features.plot, value.var = "avg.exp")
sec.time <- sec.time[,!colSums(sec.time) == 0]
sec.time <- scale(sec.time)

pheatmap(sec.time,cluster_rows = F ,cluster_cols = T)
#annotation_colors =  list(Function = anno_col)[1])

# NOT USED? Selected genes ####
selected <- DotPlot(final.all.hearts,group.by = "lineage.ident2",features = c("nrg1","nrg2a","nppa","nppb","tcf3","vegfc","emilin2a","cxcl8a","cxcr1","cxcr2"))
selected <- selected$data
selected <- selected[,c(1,2,3,4)]
selected <- acast(selected, id~features.plot, value.var = "avg.exp")
selected <- selected[!rownames(selected) %in% c("Myelin cells","Neuronal cells","Dead cells"),]
selected <- scale(selected)

pheatmap(selected,cluster_rows = T ,cluster_cols = T)


selected <- DotPlot(final.all.hearts,group.by = "lineage.ident2",features = c("vegfaa","vegfab","flt1","kdr","notch2","notch3","jag2a","jag2b","dll4"))
selected <- selected$data
selected <- selected[,c(1,2,3,4)]
selected <- acast(selected, id~features.plot, value.var = "avg.exp")
selected <- selected[!rownames(selected) %in% c("Myelin cells","Neuronal cells","Dead cells"),]
selected <- scale(selected)

pheatmap(selected,cluster_rows = T ,cluster_cols = F)


# Genes resposible for clustering 2a-2 ####
final.all.hearts <- SetIdent(final.all.hearts,value = "plot.ident")
# marker.final.all.hearts <- FindAllMarkers(final.all.hearts,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
# write.csv(marker.final.all.hearts, file = "final_clustering/finalplots/marker.final.all.hearts.csv",
#           quote = F, row.names = F)
marker.final.all.hearts <- 
  read.csv(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/marker.final.all.hearts.csv",
           stringsAsFactors = F)

# marker.fibro <- FindMarkers(final.all.hearts,ident.1 = "Fibroblasts",only.pos = T)

# niche <- SetIdent(niche, value = "work.ident")
# marker.niche <- FindAllMarkers(niche,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
# write.csv(marker.niche, file = "final_clustering/finalplots/marker.niche.csv",
#           quote = F, row.names = F)
marker.niche <- 
  read.csv(file = "/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/marker.niche.csv",
           stringsAsFactors = F)
View(marker.niche[marker.niche$cluster == "Fibroblast", ])

final.all.hearts$dotplot.ident <- final.all.hearts$lineage.ident
final.all.hearts@meta.data$dotplot.ident <- 
  plyr::mapvalues(final.all.hearts@meta.data$dotplot.ident, 
                  from =c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
                          "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
                          "Fibroblast","Smooth muscle cells","Cardiomyocytes (proliferating)",
                          "B-cells","Macrophage (CM duplex)","Macrophage (Fibroblast duplex)",
                          "Macrophage (Endothelia duplex)","Proliferating cells","Fibroblast (mpeg1.1)",
                          "Fibroblast (cfd)","Monocytes","Macrophage (il1b)",
                          "T-cells (il4/13)","Macrophages","T-cells",
                          "Macrophage (proliferating)","Macrophage (cd59)","Neutrophils",
                          "Macrophage (epdl)","T-cells (proliferating)","Macrophage (Ery duplex)",
                          "Macrophage (apoeb)","Dead cells","Fibroblast-like cells",
                          "Endocardium (frzb)","Fibroblast (nppc)","Fibroblast (spock3)",
                          "Myelin cells","Neuronal cells","Perivascular cells",
                          "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblast (cxcl12a)",
                          "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblast (col12a1a)",
                          "Fibroblast (col11a1a)","Bl.ves.EC (lyve1)","Fibroblast (proliferating)")     , 
                  to = c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
                         "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
                         "Fibroblasts (const.)","Smooth muscle cells","Cardiomyocytes (proliferating)",
                         "B-cells","Macrophages","Macrophages",
                         "Macrophage","Proliferating cells","Fibroblasts (mpeg1.1)",
                         "Fibroblasts (cfd)","Monocytes","Macrophages",
                         "T-cells","Macrophages","T-cells",
                         "Macrophages","Macrophages","Neutrophils",
                         "Macrophages","T-cells","Macrophages",
                         "Macrophages","Dead cells","Valve fibroblasts",
                         "Endocardium (frzb)","Fibroblasts (nppc)","Fibroblasts (spock3)",
                         "Myelin cells","Neuronal cells","Perivascular cells",
                         "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblasts (cxcl12a)",
                         "Epicardium (A)","Epicardium (V)","Fibroblasts (col12a1a)",
                         "Fibroblasts (col11a1a)","Bl.ves.EC (lyve1)","Fibroblasts (prolif.)"))
final.all.hearts <- SetIdent(final.all.hearts, value = "dotplot.ident")
DotPlot(final.all.hearts, 
        idents = c("Cardiomyocytes A","Cardiomyocytes V" ,"Cardiomyocytes (ttn.2) A",
                   "Cardiomyocytes (ttn.2) V","Endocardium (V)","Endocardium (A)",
                   "Fibroblasts (const.)","Smooth muscle cells","Cardiomyocytes (proliferating)",
                   "Proliferating cells","Fibroblasts (mpeg1.1)",
                   "Fibroblasts (cfd)", "Valve fibroblasts",
                   "Endocardium (frzb)","Fibroblasts (nppc)","Fibroblasts (spock3)",
                   "Perivascular cells",
                   "Bl.ves.EC (apnln)","Bl.ves.EC (plvapb)","Fibroblasts (cxcl12a)",
                   "Epicardium (A)","Epicardium (V)","Fibroblasts (col12a1a)",
                   "Fibroblasts (col11a1a)","Bl.ves.EC (lyve1)","Fibroblasts (prolif.)"),
        features = setdiff(marker.niche$gene[marker.niche$cluster == "Fibroblast"],
                           marker.final.all.hearts$gene[marker.final.all.hearts$cluster %in% c("Cardiomyocytes A", "Cardiomyocytes V")])) +
  # scale_y_discrete(breaks = waiver(), labels = , limits) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


immune <- SetIdent(immune, value = "work.ident")
# marker.immune <- FindAllMarkers(immune,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T)
# write.csv(marker.immune, file = "final_clustering/finalplots/marker.immune.csv",
#           quote = F, row.names = F)

CM <- SetIdent(CM, value = "work.ident2")
# marker.CM <- FindAllMarkers(CM,logfc.threshold = 0.3,min.pct = 0.1,only.pos = T)
# write.csv(marker.CM, file = "final_clustering/finalplots/CM.csv",
          # quote = F, row.names = F)
marker.CM <- read.csv("/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/CM.csv",
                      stringsAsFactors = F)

endo <- SetIdent(endo,value = "lineage.ident")
# marker.endo <- FindAllMarkers(endo,logfc.threshold = 0.5,min.pct = 0.1,only.pos = T,group)
# write.csv(marker.endo, file = "final_clustering/finalplots/endo.csv",
#           quote = F, row.names = F)

# Draw a heatmap of all cells for these marker genes
## allhearts
## niche
genes <- character()
for (i in unique(marker.final.all.hearts$cluster)) {
  genes <- append(genes,head(marker.final.all.hearts[marker.final.all.hearts$cluster %in% i,],300)$gene )
}

heatmap.data <- DotPlot(final.all.hearts, features = unique(genes)
)
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
heatmap.data <- heatmap.data[!rownames(heatmap.data) %in% c("Dead cells"),]
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- scale(heatmap.data)

pheatmap(t(heatmap.data),cluster_rows = F ,cluster_cols = F
)

## CM
CM <- SetIdent(CM, value = "work.ident2")
# marker.CM <- FindAllMarkers(CM,logfc.threshold = 0.3, only.pos = T)
#write.csv(marker.CM, file = "final_clustering/finalplots/CM.alternative.csv",
#          quote = F, row.names = F)
marker.CM <- read.csv("/data/junker/users/Bo/HeartRegen_paper/final_clustering/finalplots/CM.alternative.csv", stringsAsFactors = F)
  # FindAllMarkers(CM,logfc.threshold = 0.3, only.pos = T)

anno <-  rbind(
  data.frame (gene = c("bves","ttn.1","ttn.2","myom2a","myom1b","ndrg4","nppa","cav1","synpo2lb"),Function = "ttn2",color = "red"),
  data.frame (gene = c("aldoaa","eno3","pgam2","tpi1b","mdh2","idh2","atp5pd"),Function = "atp",color = "green"),
  data.frame (gene = c("pcna"),Function = "prolif",color = "blue")
  
)

# anno <-  rbind(
#   data.frame (gene = c("bves"), Function = "bves", color = "red"),
#            
#   data.frame (gene = c("ttn.1"),Function = "ttn.1",color = "green"),
#   data.frame (gene = c("ttn.2"),Function = "ttn.2",color = "blue"),
#   data.frame (gene = c("nppa"), Function = "nppa", color = "yellow")
#   
# )

rownames(anno) <- anno$gene
anno_col <- as.character(unique(anno[,3]))
names(anno_col) <- as.character(unique(anno[,2]))
#saveRDS(anno,file = "final_clustering/ecm.genelist.rds")
genes <- character()
for (i in c("Cardiomyocytes (ttn.2) V","Cardiomyocytes V","Cardiomyocytes (proliferating)") ) {
  genes <- append(genes,head(marker.CM[marker.CM$cluster %in% i,],200)$gene )
}

heatmap.data <- DotPlot(CM, features = unique(genes))
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
# heatmap.data <- heatmap.data[rownames(heatmap.data) %in% c("Cardiomyocytes (ttn.2) V","Cardiomyocytes V","Cardiomyocytes (proliferating)"),]
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- scale(heatmap.data)
# pdf(file = "/local/users/Bastiaan/Projects/heart_Bo/Images/CM_heatmap_2a.pdf",
#     width = 28, height = 5)
pheatmap(heatmap.data,cluster_rows = F ,cluster_cols = F,
         annotation_col = anno[,c(2),drop = F],
         #annotation_row = ct[,c(2),drop = F],
         
         annotation_colors =  list(Function = anno_col)[1])
# dev.off()
pheatmap(heatmap.data,cluster_rows = F ,cluster_cols = F,
         annotation_col = anno[,c(2),drop = F],
         #annotation_row = ct[,c(2),drop = F],
         
         annotation_colors =  list(Function = anno_col)[1])


#CM all 
genes <- character()
for (i in unique(marker.CM$cluster)) {
  genes <- append(genes,head(marker.CM[marker.CM$cluster %in% i,],200)$gene )
}

heatmap.data <- DotPlot(CM, features = unique(genes))
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- scale(heatmap.data)
# 2a-2
pdf("/local/users/Bastiaan/Projects/heart_Bo/Images/CM_expressions_2a_2.pdf")
pheatmap(heatmap.data,cluster_rows = F ,cluster_cols = F)
dev.off()

## endo
genes <- character()
for (i in unique(marker.endo$cluster)) {
  genes <- append(genes,head(marker.endo[marker.endo$cluster %in% i,],30)$gene )
}

heatmap.data <- DotPlot(endo, features = genes
                 )
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
#heatmap.data <- heatmap.data[!rownames(heatmap.data) %in% c("Cardiomyocytes (proliferating)"),]
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- scale(heatmap.data)


pheatmap(t(heatmap.data),cluster_rows = F ,cluster_cols = F
)

## niche
marker.niche <- read.csv("../../final_clustering/finalplots/marker.niche.csv")
genes <- character()
for (i in unique(marker.niche$cluster)) {
  genes <- append(genes,as.character(head(marker.niche[marker.niche$cluster %in% i,],30)$gene ))
}

heatmap.data <- DotPlot(niche, features = unique(genes)
)
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- heatmap.data[rownames(heatmap.data) %in% c("Epicardium (Ventricle)","Epicardium (Atrium)"), ]
heatmap.data <- scale(heatmap.data)
pheatmap(t(heatmap.data),cluster_rows = F ,cluster_cols = F
)

## immune
genes <- character()
for (i in unique(marker.immune$cluster)) {
  genes <- append(genes,head(marker.immune[marker.immune$cluster %in% i,],30)$gene )
}

heatmap.data <- DotPlot(immune, features = unique(genes)
)
heatmap.data <- heatmap.data$data
heatmap.data <- heatmap.data[,c(1,2,3,4)]
heatmap.data <- acast(heatmap.data, id~features.plot, value.var = "avg.exp")
heatmap.data <- heatmap.data[!rownames(heatmap.data) %in% c("Macrophage (Fibroblast duplex)","Macrophage (Ery duplex)","Macrophage (Endothelia duplex)","Macrophage (CM duplex)","Dead cells"),]
#ecm.tidy <- ecm.tidy[1:20,]
heatmap.data <- scale(heatmap.data)


pheatmap(t(heatmap.data),cluster_rows = F ,cluster_cols = F
)


# Niche gene expressions ####
png("/local/users/Bastiaan/Projects/heart_Bo/Images/ECM_genes_niche.png")
FeaturePlot(niche, features = c("twist1a", "prrx1a", "snai2", "fli1a", "id2b", "col1a1a"))
dev.off()
FeaturePlot(final.all.hearts, features = c("fli1a", "id2b", "col1a1a"))
DotPlot(final.all.hearts, features = c("twist1a", "prrx1a", "snai2", "fli1a", "id2b", "col1a1a"))
FeaturePlot(niche, features = c("pcna", "mki67"))
DotPlot(final.all.hearts, features = c("pcna", "mki67"))


#deep vs. shallow ####

deep <-  as.data.frame(table(niche@meta.data$orig.ident,niche@meta.data$work.ident))
colnames(deep) <- c("orig.ident","Cell.type","Freq")
deep <- merge(deep,as.data.frame(unique(niche@meta.data[,c(1,5)])),by.x = "orig.ident",by.y = "orig.ident")
for (i in unique(deep$orig.ident)) {
  deep$ratio <-  deep[deep$orig.ident == i,]$Freq / sum(deep[deep$orig.ident == i,]$Freq)
}
for (i in as.character(unique(deep$orig.ident))) {
  deep$norm.f[deep$orig.ident == i] <- sum(deep$Freq[deep$orig.ident == i]) / all.heart.norm.counts[all.heart.norm.counts$Var1 == i,]$Freq
}

deep$norm.ratio <- deep$ratio * deep$norm.f   


deep <- deep[deep$time == "7dpi",]
deep$orig.ident <- factor(deep$orig.ident, levels = c("Hr1","Hr2a","Hr2b","Hr6a","Hr6v","Hr7a","Hr7v","Hr8","Hr9","Hr13","Hr14","Hr15","Hr30","Hr31","Hr32","Hr33"  ) )
co <- is.na(match(col.table$setFibro, levels(deep$Cell.type) ))
co <- col.table[!co,]
co2 <- match(col.table$setFibro, levels(deep$Cell.type) )
co2 <- co2[!is.na(co2)]
co2 <- match(1:13,co2)
co <- co[co2,]

      ggplot(deep, aes(x=orig.ident,y=ratio,fill = Cell.type) )+ 
        geom_bar(stat = "identity")+
        scale_fill_manual(values = co$color )+
        
             #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ggtitle("Deep vs shallow")+
        theme_bw()
  plot_grid(
      ggplot(deep[deep$Cell.type == "Fibroblast (nppc)",], aes(x=orig.ident,y=norm.ratio,fill = Cell.type) )+ 
        geom_bar(stat = "identity")+
        scale_fill_manual(values = "#763e69" )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ggtitle("Deep vs shallow")+
        theme_bw(),
      ggplot(deep[deep$Cell.type == "Fibroblast-like cells",], aes(x=orig.ident,y=norm.ratio,fill = Cell.type) )+ 
        geom_bar(stat = "identity")+
        scale_fill_manual(values = "#763e69" )+
        #coord_cartesian(ylim = c(0, max(niche.norm.fibro[niche.norm.fibro$Cell.type == i,]$ratio)+0.15         ),) +  
        ggtitle("Deep vs shallow")+
        theme_bw(),
      
      
      
      
  )

  plot(deep[deep$Cell.type == "Fibroblast (nppc)",]$norm.ratio,deep[deep$Cell.type == "Fibroblast (spock3)",]$norm.ratio)
  plot(deep[deep$Cell.type == "Fibroblast (nppc)",]$ratio,deep[deep$Cell.type == "Fibroblast (spock3)",]$ratio)
  
combine_plots(plotlist = pie.plot)


#wnt and other stuff expression in all cells ####
load("/local/Bo/Remap_allhearts/final_clustering/final.all.hearts.noErynoDuplex.Robj")

genes.of.interest <- "axin2"

goi.all.hearts <- DotPlot(final.all.hearts, features = genes.of.interest, group.by = "big.ident",split.by = c("is.inhib"),cols = rep("red",100))
goi.all.hearts <- goi.all.hearts$data
goi.all.hearts <- goi.all.hearts[,c(1,3,4)]
goi.all.hearts$Cell.type <- 
  sapply(goi.all.hearts$id,  function(x) unlist(strsplit(as.character(x),"_"))[1]  )
goi.all.hearts$time <- 
  sapply(goi.all.hearts$id,  function(x) unlist(strsplit(as.character(x),"_"))[2]  )
goi.all.hearts <- goi.all.hearts[!goi.all.hearts$time == "60dpi",]
goi.all.hearts <- goi.all.hearts[!goi.all.hearts$Cell.type %in% c("Neuronal cells","Myelin cells"),]
goi.all.hearts$time <- factor(goi.all.hearts$time, levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi") )

time.res <- summarySE(goi.all.hearts,measurevar = "avg.exp",groupvars = c("features.plot","time"))

ggplot(time.res,aes(y=avg.exp,x = features.plot,fill = time))+
  geom_bar(stat = "identity",position = "dodge")+
  geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) 

type.res <- goi.all.hearts[goi.all.hearts$time %in% c("Ctrl","3dpi"),]
# type.res <- summarySE(goi.all.hearts,measurevar = "avg.exp",groupvars = c("features.plot","Cell.type") )
# ggplot(type.res,aes(y=avg.exp,x = features.plot,fill = Cell.type))+
#   geom_bar(stat = "identity",position = "dodge")+
#   geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) 


goi.plot <- list()
for (i in genes.of.interest) {
  goi.plot[[i]] <- 
    ggplot(type.res[type.res$features.plot == i,], aes(x=Cell.type,y=avg.exp,fill = time,alpha = time)) +
    geom_bar(stat = "identity",position = "dodge") +
    scale_fill_manual(values = c("black","red") ) +
    scale_alpha_manual( values = c(1,1,0.5,1,0.5,1,1))+
    #geom_errorbar(aes(ymin=norm.ratio-se, ymax=norm.ratio+se),width=.2,position=position_dodge(.9)) +
    #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
    #scale_y_continuous(limits = c(0,0.3)) +
    ylab(paste0(i," expression")) +
    ggtitle(i)
  
} 
#pdf(file = "new_clustering_inhib3/goiPlots.pdf",onefile = T,width = 16,height = 6)
for (p in 1:length(goi.plot)) {
  print(goi.plot[[p]])
}
# dev.off()

# time resolved genes after wnt inhibiton (hypo) ####
hypo.genes <- unique(c("crtc3","irs2a","ldhbb","hif1aa",rownames(final.all.hearts)[grep("hif",rownames(final.all.hearts))],rownames(final.all.hearts)[grep("^ldh",rownames(final.all.hearts))],"ckmb","epoa","epob","mylz3",
                "egln3", "irs2a","irs2b", "crtc3", "camk2g2", "ncam2","serp2","pim1","inhbb","opn5"))

hypo.genes <- "hif1aa"

final.all.hearts <- SetIdent(final.all.hearts,value = "non.myocytes")
final.all.hearts <- SetIdent(final.all.hearts,value = "myocytes",cells = WhichCells(final.all.hearts, expression = lineage.ident %in% c("Cardiomyocytes (ttn.2) A","Cardiomyocytes (ttn.2) V","Cardiomyocytes A","Cardiomyocytes V","Cardiomyocytes (proliferating)") ))


#final.all.hearts@meta.data$is.inhib <- as.character(final.all.hearts@meta.data$time)
#final.all.hearts@meta.data[final.all.hearts@meta.data$orig.ident %in% c("Hr29","Hr34","Hr35"),]$is.inhib <- "3dpiinhib"
#final.all.hearts@meta.data[final.all.hearts@meta.data$orig.ident %in% c("Hr31","Hr32","Hr33"),]$is.inhib <- "7dpiinhib"

hypo.genesall.hearts <- DotPlot(final.all.hearts, features = hypo.genes, group.by  = "orig.ident",cols = rep("red",100))
hypo.genesall.hearts <- hypo.genesall.hearts$data
hypo.genesall.hearts <- hypo.genesall.hearts[,c(1,3,4)]
hypo.genesall.hearts$Cell.type <- 
  sapply(hypo.genesall.hearts$id,  function(x) unlist(strsplit(as.character(x),"_"))[1]  )
hypo.genesall.hearts$sample <- 
  sapply(hypo.genesall.hearts$id,  function(x) unlist(strsplit(as.character(x),"_"))[2]  )
hypo.genesall.hearts <- hypo.genesall.hearts[!hypo.genesall.hearts$time == "60dpi",]
hypo.genesall.hearts <- hypo.genesall.hearts[!hypo.genesall.hearts$time == "15dpi",]
#hypo.genesall.hearts <- hypo.genesall.hearts[hypo.genesall.hearts$Cell.type %in% c("Cardiomyocytes","Endocardium","Fibroblasts"),]
#hypo.genesall.hearts$time <- factor(hypo.genesall.hearts$time, levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi") )

hypo.genesall.hearts$time <- NA
hypo.genesall.hearts[hypo.genesall.hearts$id %in% c("H5","H6","H7","H8a","H8v"),]$time <- "Ctrl"
hypo.genesall.hearts[hypo.genesall.hearts$id %in% c("Hr10","Hr11","Hr12","Hr22","Hr23","Hr24","Hr25","Hr26","Hr27","Hr28","Hr29","Hr34","Hr35"),]$time <- "3dpi"
hypo.genesall.hearts[hypo.genesall.hearts$id %in% c("Hr1","Hr2a","Hr2b","Hr8","Hr9","Hr13","Hr14","Hr15","Hr6a","Hr6v","Hr7a","Hr7v","Hr30","Hr31","Hr32","Hr33"),]$time <- "7dpi"
hypo.genesall.hearts[hypo.genesall.hearts$id %in% c("Hr3","Hr4","Hr19","Hr20","Hr21"),]$time <- "30dpi"
hypo.genesall.hearts[hypo.genesall.hearts$id %in%  c("Hr29","Hr34","Hr35"),]$time <- "3dpiinhib"
hypo.genesall.hearts[hypo.genesall.hearts$id %in%  c("Hr31","Hr32","Hr33"),]$time <- "7dpiinhib"

hypo.genesall.hearts$time <- factor(hypo.genesall.hearts$time, levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi") )

hypo.genesall.hearts.stats <- summarySE(hypo.genesall.hearts,measurevar = "avg.exp",groupvars = c("time") )

hypo.genes.plot <- list()
hypo.genes.plot[[1]] <- 
    combine_plots(
    ggplot(hypo.genesall.hearts.stats[hypo.genesall.hearts.stats$Cell.type == "myocytes",], aes(x=time,y=avg.exp,alpha = time)) +
    geom_bar(stat = "identity",position = "dodge") +
    #scale_fill_manual(values = c("black","black","black","black","grey","black") ) +
    geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) +
    #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
    theme(legend.position = "none")  +
    scale_alpha_manual( values = c(1,1,0.5,1,0.5,1,1))+
    #ylim(0,50)+
    #scale_y_continuous(limits = c(0,0.3)) +
    ylab(paste0("Average scaled hif1aa expression")) +
    ggtitle("In myocytes") ) 
  
hypo.genes.plot[[2]] <- 
  combine_plots(
    ggplot(hypo.genesall.hearts.stats[hypo.genesall.hearts.stats$Cell.type == "non.myocytes",], aes(x=time,y=avg.exp,alpha = time)) +
      geom_bar(stat = "identity",position = "dodge") +
      #scale_fill_manual(values = c("black","black","black","black","grey","black") ) +
      geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) +
      #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
      theme(legend.position = "none")  +
      scale_alpha_manual( values = c(1,1,0.5,1,0.5,1,1))+
      #ylim(0,50)+
      #scale_y_continuous(limits = c(0,0.3)) +
      ylab(paste0("Average scaled hif1aa expression")) +
      ggtitle("In non-myocytes") ) 


ggplot(hypo.genesall.hearts.stats, aes(x=time,y=avg.exp,alpha = time)) +
  geom_bar(stat = "identity",position = "dodge") +
  #scale_fill_manual(values = c("black","black","black","black","grey","black") ) +
  geom_errorbar(aes(ymin=avg.exp-se, ymax=avg.exp+se),width=.2,position=position_dodge(.9)) +
  #geom_text(aes(label=time),position=position_dodge(width=0.9), vjust=0.5,hjust = -0.25, angle = 90) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5,face = "bold")) +
  theme(legend.position = "none")  +
  scale_alpha_manual( values = c(1,1,0.5,1,0.5,1,1))+
  #ylim(0,50)+
  #scale_y_continuous(limits = c(0,0.3)) +
  ylab(paste0("Average scaled hif1aa expression")) +
  ggtitle("In all cells")

combine_plots(plotlist = hypo.genes.plot)
pdf(file = "/local/Bo/Remap_allhearts/final_clustering/finalplots/nppcPlots.pdf",onefile = T,width = 8,height = 4.5)
for (p in 1:length(hypo.genes.plot)) {
  print(hypo.genes.plot[[p]])
}
dev.off()



# wnt genes heatmap ####
# Collapsing cell types into larger types
final.all.hearts@meta.data$trashid <- final.all.hearts@meta.data$plot.ident2
final.all.hearts@meta.data$trashid <- 
  plyr::mapvalues(final.all.hearts@meta.data$trashid, 
                  from = c("Cardiomyocytes A", "Cardiomyocytes V", "Cardiomyocytes (ttn.2) V", "Cardiomyocytes (ttn.2) A",
                           "Cardiomyocytes (proliferating)"),
                  to = rep("Cardiomyocytes", 5))
final.all.hearts@meta.data$trashid <- 
  plyr::mapvalues(final.all.hearts@meta.data$trashid, 
                  from = c("Endocardium (A)", "Endocardium (V)", "Endocardium (frzb)"),
                  to = rep("Endocardium", 3))
final.all.hearts@meta.data$trashid <- 
  plyr::mapvalues(final.all.hearts@meta.data$trashid, 
                  from = c("Epicardium (Atrium)", "Epicardium (Ventricle)", "Fibroblast-like cells", "Monocytes"),
                  to = c("Fibroblasts", "Fibroblasts", "Valve fibroblasts", "Macrophages"))


# final.all.hearts@meta.data$trashid <- plyr::mapvalues(final.all.hearts@meta.data$trashid, 
#                                                       from = as.character(unique(final.all.hearts$trashid)),
#                                                       to = c(
#                                                         "Cardiomyocytes (Atrium)","Cardiomyocytes (Ventricle)","Cardiomyocytes (dediff.) A","Endocardium (Ventricle)","Endocardium (Atrium)","Fibroblasts (const.)",                 
#                                                         "Smooth muscle cells","Cardiomyocytes (dediff) V","Cardiomyocytes (prolif.)","B-cells","Marcrophages","Proliferating cells",         
#                                                         "Fibroblasts (mpeg1.1)","Fibroblasts (cfd)","Monocytes","T-cells","T-cells","Neutrophils",
#                                                         "T-cells (proliferating)","Dead cells","Valve Fibroblasts","Endocardium (frzb)","Fibroblasts (nppc)","Fibroblasts (spock3)",           
#                                                         "Myelin cells","Neuronal cells","Perivascular cells","Endothelial cells (apnln)","Endothelial cells (plvapb)","Fibroblasts (cxcl12a)",        
#                                                         "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblasts (col12a1a)","Fibroblasts (col11a1a)","Endothelial cells (lyve1)","Fibroblasts (proliferating)")
# )
domini <-  rbind(
  data.frame (gene = c("wnt1","wnt2","wnt2ba","wnt2bb","wnt3","wnt3a","wnt4","wnt5a","wnt5b","wnt6a","wnt6b","wnt7aa","wnt7ab","wnt7ba","wnt7bb","wnt8a","wnt8b","wnt9a","wnt9b","wnt10a","wnt10b","wnt11","wnt11f2","wnt16","wif1","waif2"),Function = "wnt-ligand",color = "#1b9d77"),
  data.frame (gene = c("fzd1","fzd2","fzd3a","fzd3b","fzd4","fzd6","fzd7a","fzd7b","fzd8a","fzd8b","fzd9a","fzd9b","fzd10","frzb","adgra2","nid1a","nid1b","szl","lrp5","lrp6","ror1","ror2","ryk"),Function = "wnt-receptor",color = "#d76126"),
  data.frame (gene = c("dkk1a","dkk1b","dkk2","dkk3a","dkk3b","sfrp1a","sfrp1b","sfrp2","sfrp2l","sfrp4","sfrp5"),Function = "wnt-modulator",color = "#7670b3")
  #data.frame (gene = c("aldh1a2"),Function = "RA-ligand",color = "#7670b3"),
  #data.frame (gene = c("cyp26b1","stra6"),Function = "RA-receptor",color = "#7670b3")
)
rownames(domini) <- domini$gene
domini_col <- as.character(unique(domini[,3]))
names(domini_col) <- as.character(unique(domini[,2]))

wnt.all <- DotPlot(final.all.hearts, features = as.character(domini$gene),group.by = "trashid",split.by = "inhib") #"big.ident"
#wnt.all <- DotPlot(final.all.hearts, features = as.character(domini$gene),group.by = "big.ident")

wnt.all <- wnt.all$data
wnt.all <- wnt.all[,c(1,2,3,4)]
#ecm$value <- ecm$avg.exp*ecm$pct.exp*100
wnt.all <- acast(wnt.all, id~features.plot, value.var = "avg.exp")
wnt.all <- wnt.all[!rownames(wnt.all) %in% c("Neuronal cells_NULL","Myelin cells_NULL",
                                             rownames(wnt.all)[grep("IWR1$",rownames(wnt.all))]) ,]


#wnt.all <- wnt.all[!rownames(wnt.all) %in% c(paste0(c("Myelin cells","Neuronal cells","Macrophage (CM duplex)","Macrophage (Fibroblast duplex)","Macrophage (Ery duplex)","Dead cells","Fibroblast (mpeg1.1)"),"_IWR1"),
#paste0(c("Myelin cells","Neuronal cells","Macrophage (CM duplex)","Macrophage (Fibroblast duplex)","Macrophage (Ery duplex)","Dead cells","Fibroblast (mpeg1.1)"),"_NULL")),]
#wnt.all.1 <- wnt.all[rownames(wnt.all) %in% paste0(unique(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1]  )),"_NULL"),]
#wnt.all.2 <- wnt.all[rownames(wnt.all) %in% paste0(unique(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1]  )),"_IWR1"),]

#wnt.all <- wnt.all[c(paste0(unique(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1]  )),"_NULL"),paste0(unique(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1]  )),"_IWR1")),]
wnt.all<- scale(wnt.all)
#wnt.all.2 <- scale(wnt.all.2)
wnt.heatmap.list <- list()
pheatmap(wnt.all,cluster_rows = T ,cluster_cols = T,
                              annotation_col = domini[,c(2),drop = F],
                              #annotation_row = ct[,c(2),drop = F],
                              annotation_colors =  list(Function = domini_col)[1])
wnt.heatmap.list[[2]] <- plot_grid(plotlist = pheatmap(wnt.all,cluster_rows = T ,cluster_cols = T,
                                                       annotation_col = domini[,c(2),drop = F],
                                                       #annotation_row = ct[,c(2),drop = F],
                                                       annotation_colors =  list(Function = domini_col)[1])[4])

plot_grid(plotlist = wnt.heatmap.list)

#wnt in subtypes
domini <-  rbind(
  data.frame (gene = c("wnt1","wnt2","wnt2ba","wnt2bb","wnt3","wnt3a","wnt4","wnt5a","wnt5b","wnt6a","wnt6b","wnt7aa","wnt7ab","wnt7ba","wnt7bb","wnt8a","wnt8b","wnt9a","wnt9b","wnt10a","wnt10b","wnt11","wnt11f2","wnt16","wif1","waif2"),Function = "wnt-ligand",color = "#1b9d77"),
  data.frame (gene = c("fzd1","fzd2","fzd3a","fzd3b","fzd4","fzd6","fzd7a","fzd7b","fzd8a","fzd8b","fzd9a","fzd9b","fzd10","frzb","adgra2","nid1a","nid1b","szl","lrp5","lrp6","ror1","ror2","ryk"),Function = "wnt-receptor",color = "#d76126"),
  data.frame (gene = c("dkk1a","dkk1b","dkk2","dkk3a","dkk3b","sfrp1a","sfrp1b","sfrp2","sfrp2l","sfrp4","sfrp5"),Function = "wnt-modulator",color = "#7670b3")
  #data.frame (gene = c("aldh1a2"),Function = "RA-ligand",color = "#7670b3"),
  #data.frame (gene = c("cyp26b1","stra6"),Function = "RA-receptor",color = "#7670b3")
)
rownames(domini) <- domini$gene
domini_col <- as.character(unique(domini[,3]))
names(domini_col) <- as.character(unique(domini[,2]))

# wnt.all <- DotPlot(final.all.hearts, features = as.character(domini$gene),group.by = "trashid",split.by = "inhib") #group.by = "trashid"
wnt.all <- DotPlot(final.all.hearts, features = as.character(domini$gene),group.by = "big.ident",split.by = "inhib")

wnt.all <- wnt.all$data
wnt.all <- wnt.all[,c(1,2,3,4)]
#ecm$value <- ecm$avg.exp*ecm$pct.exp*100
wnt.all <- acast(wnt.all, id~features.plot, value.var = "avg.exp")
wnt.all <- wnt.all[!rownames(wnt.all) %in% c("Neuronal cells_NULL","Myelin cells_NULL","Dead cells_NULL",
                                             rownames(wnt.all)[grep("IWR1$",rownames(wnt.all))]) ,]
rownames(wnt.all) <- as.character(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1] ))
wnt.all<- scale(wnt.all)
#wnt.all.2 <- scale(wnt.all.2)
#wnt.heatmap.list <- list()
pheatmap(wnt.all,cluster_rows = T ,cluster_cols = T,
         annotation_col = domini[,c(2),drop = F],
         #annotation_row = ct[,c(2),drop = F],
         annotation_colors =  list(Function = domini_col)[1])

wnt_genes_hclust <- hclust(wnt.all)

niche_annotations <- niche@meta.data
niche_annotations$Niche_plus_anno <- niche_annotations$work.ident
niche_annotations <- niche_annotations[, "Niche_plus_anno", drop = F]
niche_annotations$Niche_plus_anno <- as.character(niche_annotations$Niche_plus_anno)

test_all_annotations <- final.all.hearts@meta.data
test_all_annotations <- merge(test_all_annotations, niche_annotations, all.x = T,
                              by = "row.names")
rownames(test_all_annotations) <- test_all_annotations$Row.names
test_all_annotations <- test_all_annotations[, 2:ncol(test_all_annotations)]
test_all_annotations$Niche_plus_anno[is.na(test_all_annotations$Niche_plus_anno)] <-
  as.character(test_all_annotations$plot.ident2[is.na(test_all_annotations$Niche_plus_anno)])
final.all.hearts@meta.data <- test_all_annotations

wnt.all <- DotPlot(final.all.hearts, features = as.character(domini$gene),group.by = "Niche_plus_anno",split.by = "inhib")
wnt.all <- wnt.all$data
wnt.all <- wnt.all[,c(1,2,3,4)]
wnt.all <- acast(wnt.all, id~features.plot, value.var = "avg.exp")
wnt.all <- wnt.all[!rownames(wnt.all) %in% c("Neuronal cells_NULL","Myelin cells_NULL","Dead cells_NULL",rownames(wnt.all)[grep("IWR1$",rownames(wnt.all))]) ,]
rownames(wnt.all) <- as.character(sapply(rownames(wnt.all),  function(x) unlist(strsplit(x,"_"))[1] ))
wnt.all<- scale(wnt.all)
pheatmap(wnt.all,cluster_rows = T ,cluster_cols = T,
         annotation_col = domini[,c(2),drop = F],
         annotation_colors =  list(Function = domini_col)[1])


wnt.niche <- DotPlot(niche, features = as.character(domini$gene),group.by = "work.ident",split.by = "inhib")

wnt.niche <- wnt.niche$data
wnt.niche <- wnt.niche[,c(1,2,3,4)]
wnt.niche <- acast(wnt.niche, id~features.plot, value.var = "avg.exp")
wnt.niche <- wnt.niche[!rownames(wnt.niche) %in% c("Neuronal cells_NULL","Myelin cells_NULL","Dead cells_NULL",rownames(wnt.niche)[grep("IWR1$",rownames(wnt.niche))]) ,]
rownames(wnt.niche) <- as.character(sapply(rownames(wnt.niche),  function(x) unlist(strsplit(x,"_"))[1] ))
wnt.niche<- scale(wnt.niche)
wnt.niche <- wnt.niche[, !is.na(colSums(wnt.niche))] # Remove any genes that got NA after scaling (i.e. no variable expression)
# pdf("/local/users/Bastiaan/Projects/heart_Bo/Images/Wnt_niche.pdf")
pheatmap(wnt.niche,cluster_rows = T ,cluster_cols = T,
         annotation_col = domini[,c(2),drop = F],
         annotation_colors =  list(Function = domini_col)[1])
# dev.off()

#change cell types ####
final.all.hearts@meta.data$trashid <- final.all.hearts@meta.data$lineage.ident2
final.all.hearts@meta.data$trashid <- plyr::mapvalues(final.all.hearts@meta.data$trashid, 
                                                      from = c("Macrophage (apoeb)","Macrophage (cd59)","Macrophage (CM duplex)","Macrophage (Endothelia duplex)","Macrophage (epdl)","Macrophage (Ery duplex)","Macrophage (Fibroblast duplex)","Macrophage (il1b)","Macrophage (proliferating)","Macrophages"),
                                                      to = rep("Marcrophages",10))


final.all.hearts@meta.data$trashid <- plyr::mapvalues(final.all.hearts@meta.data$trashid, 
                                                      from = as.character(unique(final.all.hearts$trashid)),
                                                      to = c(
                                                        "Cardiomyocytes (Atrium)","Cardiomyocytes (Ventricle)","Cardiomyocytes (dediff.) A","Endocardium (Ventricle)","Endocardium (Atrium)","Fibroblasts (const.)",                 
                                                        "Smooth muscle cells","Cardiomyocytes (dediff) V","Cardiomyocytes (prolif.)","B-cells","Marcrophages","Proliferating cells",         
                                                        "Fibroblasts (mpeg1.1)","Fibroblasts (cfd)","Monocytes","T-cells","T-cells","Neutrophils",
                                                        "T-cells (proliferating)","Dead cells","Valve Fibroblasts","Endocardium (frzb)","Fibroblasts (nppc)","Fibroblasts (spock3)",           
                                                        "Myelin cells","Neuronal cells","Perivascular cells","Endothelial cells (apnln)","Endothelial cells (plvapb)","Fibroblasts (cxcl12a)",        
                                                        "Epicardium (Atrium)","Epicardium (Ventricle)","Fibroblasts (col12a1a)","Fibroblasts (col11a1a)","Endothelial cells (lyve1)","Fibroblasts (proliferating)")
                                                      )



# Fibroblast motivation 2b ####
final.all.hearts <- SetIdent(final.all.hearts,value = "big.ident")

genes <- c("aldh1a2","frzb","dkk3b","igf2b","pdgfrb","tgfb1a","shha","jak1","notch1b","tbx18","bmp4")
genes <- c("aldh1a2","frzb","dkk3b","tbx18","igf2b","pdgfrb","tgfb1a","axin2","nrg1")

moti <- AverageExpression(final.all.hearts,features = "nrg1") #nrg1, aldh1a2, fn1a to reproduce 2B
moti <- moti$RNA
moti <- moti[,c(1:14)]
moti$non.Fibroblasts <- rowMeans(moti[,!colnames(moti) == "Fibroblasts"])
moti$ratio <- moti$Fibroblasts / moti$non.Fibroblasts
#moti$ratio.log <- log2(moti$ratio)

moti2 <- data.frame(gene = c(rownames(moti),rownames(moti)),
                    value = c(moti[,3],moti[,15]),
                    mode = c(rep("Fibro",1),rep("non.Fibro",1)) )   
  
test <- DotPlot(niche,feature = "acta2",group.by = "time")  
test <- test$data
ggplot(test,aes(x = factor(test$id,levels = c("Ctrl","3dpi","7dpi","30dpi")), y= avg.exp))+
  geom_bar(stat = "identity")+
  theme_bw()

test <- DotPlot(niche,feature = "acta2",split.by = "time",cols = rep("red",12))  $data
test <- test[test$id %in% c("Fibroblast (cfd)_Ctrl","Fibroblast (cfd)_3dpi","Fibroblast (cfd)_7dpi","Fibroblast (cfd)_30dpi"), ]

ggplot(test,aes(x = factor(test$id,levels = c("Fibroblast (cfd)_Ctrl","Fibroblast (cfd)_3dpi","Fibroblast (cfd)_7dpi","Fibroblast (cfd)_30dpi")), y= avg.exp))+
  geom_bar(stat = "identity")+
  theme_bw()

ggplot(moti,aes(x = factor(rownames(moti),levels = c("aldh1a2","frzb","dkk3b","nrg1","igf2b","pdgfrb","tgfb1a")  )  , y= ratio.log))+
  geom_bar(stat = "identity")+
  theme_bw()

ggplot(moti2,aes(x = factor(gene,levels = c("aldh1a2","frzb","dkk3b","tbx18","igf2b","pdgfrb","tgfb1a","axin2","nrg1")  ), y= value,alpha = mode))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_alpha_manual(values = c(1,0.5))+
  theme_bw()
ggplot(moti2,aes(x = factor(gene ), y= value,alpha = mode))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_alpha_manual(values = c(1,0.5))+
  theme_bw()

#pd
FeaturePlot(final.all.hearts,features = c("nrg1"),pt.size = 1,cols = c("#f5f5f5","#9d3d58"))
VlnPlot(final.all.hearts, features = c("nrg1", "aldh1a2", "fn1a")) # In rebuttal
plot_grid(
  FeaturePlot(final.all.hearts, features = c("nrg1"), order = T) +
    coord_flip() + scale_y_reverse(),
  FeaturePlot(final.all.hearts, features = c("aldh1a2"), order = T) +
    coord_flip() + scale_y_reverse(),
  FeaturePlot(final.all.hearts, features = c("fn1a"), order = T) +
    coord_flip() + scale_y_reverse(),
  nrow = 1
)

FeaturePlot(final.all.hearts, features = c("nrg1", "aldh1a2", "fn1a"), order = T) +
  coord_flip()

moti3 <- (AverageExpression(final.all.hearts, features = c("nrg1", "aldh1a2", "fn1a")))$RNA


niche <- SetIdent(niche,value = "work.ident")

moti <- AverageExpression(niche,features = c("axin2","wnt8a"))
moti <- moti$RNA
moti <- as.data.frame(t(moti))

ggplot(moti,aes(x = rownames(moti) , y= axin2))+
  geom_bar(stat = "identity")+
  theme_bw() +
  RotatedAxis()


#pd
FeaturePlot(niche,features = c("axin2"),pt.size = 1,cols = c("#f5f5f5","#9d3d58"))

valve.plot <- list()

test.matrix <- DotPlot(final.all.hearts,features = c("cxcl12b","cxcr4a"),group.by = "big.ident")
test.matrix <- test.matrix$data

test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells","T-cells","B-cells"),]

test.matrix <- scale(test.matrix)
valve.plot[["one"]] <- plot_grid(plotlist = pheatmap(test.matrix,cluster_rows = F,cluster_cols = T)[4])


test.matrix <- DotPlot(final.all.hearts,features = c("snai2","snai1a","snai1b","twist1a","twist1b"),group.by = "lineage.ident2")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells"),]
#test.matrix <- scale(test.matrix)
test.matrix <- test.matrix[rownames(test.matrix) %in% c("Fibroblast (col11a1a)","Fibroblast (col12a1a)","Fibroblast (proliferating)","Fibroblast (nppc)","Fibroblast (cxcl12a)","Epicardium (Atrium)","Epicardium (Ventricle)","Perivascular cells","Bl.ves.EC (plvapb)","Fibroblast (nppc)","Fibroblast (spock3)","Fibroblast-like cells"),]
test.matrix <- scale(test.matrix)
pheatmap(test.matrix,cluster_rows = T,cluster_cols = F)

test.matrix <- DotPlot(final.all.hearts,features = unique(c("nppc","serpine1","kcne4","lepb","vmp1","serpinh1b","krt18","col1a1a","clic2","crlf1a","zgc:158343","hspb8","tfpia","ptgs2a","ctsla","sele","rhoca","aldh1a2","vwf","zfand2a","si:ch211-248e11.2","efna1b","hspa5","foxp4","sox11a","si:ch211-222l21.1","phlda2","crip2","cdh5","CR383676.1","sult2st3","txn","col1a1b","si:ch211-153b23.5","col5a1","col1a2","ramp2","bzw1b","glulb","dcn","si:dkey-7j14.6","edn2","rdh10b","tgif1","fstl1b","col5a2a","hsp90b1","ackr3b","trib3","arl4ab","si:ch211-160j14.3","postna","manf","slc38a2","inhbaa","net1","herpud1","cd276","foxc1b","calr3a","mbd2","zfand5a","adgrl4","cnn3a","errfi1a","sox7","ctsk","slit2","cyp3c1","tspan35","mvp","matn4","tmem45a","il6st","tuba8l3","bmp6","bmp16","iqgap2","il13ra1","tp53inp1","egfl7","f2rl2","chl1a","si:ch211-33e4.2","arf5","pdia6","FP102018.1","grinab","swap70b","serpinh1a","nrarpa","bambia","cdkn1bb","shroom4","efemp2b","fli1a","hmcn2.1","adcy2b","prnprs3","s1pr1","il34","twist1b","si:dkey-219e21.2","epha2a","fkbp1b","smad3a","ptgs2b","gng2","fam129ba","pde4ba","sptlc2b","si:dkey-23i12.7","plpp3","ptprk","s100u","cx43.4","fhl3a","si:ch211-18i17.2","tie1","tspan7","adarb1a","tcf12","ndel1b","angptl2b","slc38a5a","nuak2","nos1apa.1","agrn","smad7","si:ch211-14k19.8","p4ha2","aplnra","calua","eva1ba","zgc:85975","ecscr","zswim5","ms4a17a.11","fndc3ba","cxcl20","CR847566.1","sulf2b","prcp","nos1apa","colec12","selenon","CABZ01030107.1","nid2a","zgc:113363","ptprb",
                                                            "CU469362.1","ccdc80","si:dkey-19b23.8","c1s","hyal2b","bmp7a","si:ch211-284f22.3","tusc3","tram2","ergic3","sema6e","mdka","slit3","mras","zgc:77086","ptprua","cilp","fem1b","foxc1a","cpxm1b","rnf185","sat1a.2","pawr","ddx61","plod1a","yap1","cd63","hapln1b","cd82a","fkbp7","tgfb1a","cxcl8a","dusp5","ugdh","hspb1","lysmd3","mfap2","sik1","irf2","capgb","txndc5","anxa4","si:dkeyp-27c8.2","maff","f10","rasgef1bb","zgc:92242","si:ch1073-459j12.1","calr","igf2b","cdipt","krt8","tuba8l","ackr4b","tmem88a","ctgfa","gchfr","nucb1","tmem88b","tuba1c.1","rgcc","mxra8b","btg3","nck1b","hbegfa","tmem176","tspan4b","zfand5b","glsb","sdf2l1","cwc15","klf2b","tmem33","pnp5a","fkbp9","cbx7a","sparc","bach1b","DNAJA4","acsbg2","ostc","arrdc3b","ptp4a2a","zfp36l1b","fosl2","sept2","si:ch73-63e15.2","id2b","tmed9","itpkcb","oaz1a","dap1b","s100v2","dnajb11","serp1","snu13b","ptmab","lman2","si:ch73-335l21.4","zgc:195173","jpt1b","col6a1","ahsa1b","rpz5","anxa13","crebzf","si:ch211-15b10.6","atf4b","si:dkey-25e12.3","gmfb","si:dkey-177p2.6","si:ch211-71m22.1","psmb7","prelid1a","hsph1","pdia4","kdelr2a","hyou1","perp","cacybp","plk3","zgc:66160","srgn","tfg","kdelr2b","oser1","fhl1b","cxcl18b","p4ha1b","prkcsh","anxa2a","si:dkey-4p15.3","si:ch211-156j16.1","tmem50a","vcp","tmed10","arf2a","zgc:77849","adam10a","aqp1a.1","pfn2l","cd81a","slc25a33","im:7152348","h3f3d","CU929418.1","hsp90ab1","zgc:92066","canx","tspan36","tob1b","gstm.1",
                                                            "slc20a1b","eif4a1a","fkbp1aa","ddx3a","ppp1r10","elmsan1b","cyr61","uchl1","tpt1","crip1","spcs1","arpc1a","nppb","prelid3b","gapdhs","pfn2","chchd2","nupr1","bcl2l10","si:dkey-42i9.4","brdt","tmed2","midn","tnfrsf9a","si:busm1-57f23.1","CABZ01079011.1","lamp2","cct7","hspbp1","dad1","pdia3","chac1","mycb","ppib","selenof","pgam1a","sdcbp2","rcan1a","fkbp4","ip6k2b","tmem258","hnrnpa0l.1","cox6a1","pim1","rbm4.3","sec61b","rtn4a","hint1","eif5a2","tmem176l.1","FO704736.1","ssr3","dynll1","ssr4","atp1b1a","hspa8","ddx5","fosl1a","erh","brd2a","eif1b","hspe1","rpl22l1","rps6","hsp70.3","jdp2b","btg1","hsp90aa1.2","hnrnpa0b","mcl1a","ppiaa","rps12","rplp2l","rpl3","mhc1uba","zgc:153867","si:ch211-202a12.4","eef1b2","icn","msna","s100a10b","zgc:198419","fn1b","dcn","col1a1a","si:ch1073-459j12.1","postnb","mmp2","col1a2","col1a1b","col5a1","dpt","col6a1","gstm.3","c6.1","c4","mfap5","htra1b","col6a3","mdka","edil3a","pmp22a","col6a2","frzb","angptl1a","si:ch211-105c13.3","sparc","c1qtnf5","mxra8a","tnfaip6","olfml3b","nid1b","krt15","rbp4","c6","tmsb2","spock3","angptl7","clu","ccl25b","icn","thbs3a","col12a1a","mxra8b","mgp","sostdc1a","col5a2a","colec12","nid1a","f3b","fstl1b","cyp26b1","pcolcea","anxa2a","pdgfrl","krt94","glis2a","krt4","tbx18","abi3bpb","tfa","prss23","timp2a","mfap2","fibinb","tcf21","lama5","cxcl18b","zgc:195173","aebp1","serpinf1","tmem176","slco2b1","pcolce2b","c1qtnf2","scara3","spon1b","s100a10b",
                                                            "ackr4b","col18a1b","mmp14a","thbs2a","foxp4","lum","igfbp5b","dkk3b","hapln1a","rnd3b","cd81a","s100a10a","boc","cdh2","ifitm1","epas1b","appa","kirrel1a","pdgfra","f10","her6","hsd11b2","crip1","si:dkey-8k3.2","serpinh1b","yap1","antxr1c","acvrl1","angptl2b","si:ch211-106h4.12","nupr1","endouc","tmem88b","lrp1aa","si:ch211-133n4.4","s100w","podxl","egr1","zfp36l1b","ckba","fkbp7","FP102018.1","eppk1","si:ch1073-291c23.2","plpp3","s100v2","cxcl12a","ablim1b","fxyd6l","akap12b","aldh1a2","cyp2ad3","col11a1b","dap1b","adh8a","cygb1","thbs1b","tcf12","enc1","lxn","si:ch211-137i24.10","s100u","im:7152348","zcchc24","si:ch73-46j18.5","krt18","cnn3a","fn1a","efemp2a","anxa1a","lama4","cav1","si:dkey-7j14.6","qkia","serpina10a","cavin2b","si:busm1-57f23.1","mmel1","tmem176l.1","ppib","faua","rpl37.1","wls","socs3a","add3a","zfp36l2","cemip","rplp2l","crip2","cxcl14","krt8","id3","si:ch73-138n13.1","si:ch211-156j16.1","CR936442.1","jpt1b","tfpia","arl6ip5b","CR383676.1","mmp14b","ctsk","selenop","RPL41","si:ch211-222l21.1","hdlbpa","cebpd","pfn2l","bambia","si:ch211-195m9.3","si:ch211-198c19.3","pleca","rps28","tsc22d3","kcne4","adam10a","calua","steap4","rpl5b","acin1a","her9","serpine1","spaca4l","si:ch73-1a9.3","aldh9a1a.1","p4hb","txn","CR753876.1")),
                       group.by = "lineage.ident2")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells"),]
#test.matrix <- scale(test.matrix)
#test.matrix <- test.matrix[rownames(test.matrix) %in% c("Fibroblast-like cells","Fibroblast (col11a1a)","Fibroblast (col12a1a)","Fibroblast (proliferating)","Fibroblast (nppc)","Fibroblast (cxcl12a)","Epicardium (Atrium)","Epicardium (Ventricle)","Perivascular cells","Bl.ves.EC (plvapb)","Fibroblast (nppc)","Fibroblast (spock3)"),]
test.matrix <- scale(test.matrix)
pheatmap(test.matrix,cluster_rows = T,cluster_cols = T)





#timesplitt
test.matrix <- DotPlot(final.all.hearts,features = c("axin2","wif1"),group.by = "lineage.ident2",cols = rep("red",100))
test.matrix <- test.matrix$data

test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[rev(test),]


test.matrix <- scale(test.matrix)
pheatmap(test.matrix,cluster_rows = F,cluster_cols = F)

test <- order(test.matrix)
pheatmap(test,cluster_rows = F,cluster_cols = F)

FeaturePlot(niche,features = c("wif1"),
            pt.size = 1,cols = c("#e4e4e4","red"))

# jpeg(filename = "../finalplots/literature.genes.jpeg",height = 3000,width = 3000)
FeaturePlot(niche,features = c("postnb","col11a1a","col12a1a","fn1a","elnb","sfrp1a","sfrp1b","wif1","dkk3b","frzb","acta2","scxa","tagln","prrx1b","nrg1","snai2","tgfb1a","tgfb3","fn1a"),
              pt.size = 1,cols = c("#e4e4e4","red"))
# dev.off()

# jpeg(filename = "../finalplots/literature.genes2.jpeg",height = 600,width = 600)
FeaturePlot(niche,features = c("fn1a"),
            pt.size = 1,cols = c("#e4e4e4","red"))
# dev.off()

# nppc ####
test.matrix <- DotPlot(final.all.hearts,features = "nppc",group.by = "is.inhib")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells"),]
test.matrix <- as.data.frame(t(test.matrix))
test.matrix <- as.data.frame(t(test.matrix))
nppc <- ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y=V1))+
  geom_bar(stat = "identity") + RotatedAxis() + ggtitle("nppc")

test.matrix <- DotPlot(final.all.hearts,features = "fli1a",group.by = "is.inhib")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells"),]
test.matrix <- as.data.frame(t(test.matrix))
test.matrix <- as.data.frame(t(test.matrix))
fli1a <- ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y=V1))+
  geom_bar(stat = "identity") + RotatedAxis()+ ggtitle("fli1a")

test.matrix <- DotPlot(final.all.hearts,features = "pdgfrb",group.by = "is.inhib")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Dead cells","Myelin cells","Neuronal cells"),]
test.matrix <- as.data.frame(t(test.matrix))
test.matrix <- as.data.frame(t(test.matrix))
pdgfrb <- ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y=V1))+
  geom_bar(stat = "identity") + RotatedAxis()+ ggtitle("pdgfrb")

plot_grid(nppc,fli1a,pdgfrb)

# NOT USED? ####
test.matrix <- DotPlot(niche,features = c("tcf21"),group.by = "work.ident")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- as.data.frame(t(test.matrix))
test.matrix <- as.data.frame(t(test.matrix))

plot_grid(FeaturePlot(niche,"tcf21"),
DimPlot(niche,label = T,group.by = "work.ident",
        cols=colors[match(
          #consistent:
          c("syn-lightyellow","darkyellow","creamyellow2",
            "lightpink","lightblue","green2",
            "yellowochre","syn-red","darkpurple",
            "purple","lightpurple","white2",
            "darkblue2")
          ,colors$name),]$color,pt.size = 1),
pheatmap(test.matrix,cluster_rows = F,cluster_cols = F)[[4]],
ggplot(test.matrix,aes(x=rownames(test.matrix),y=tcf21))+
  geom_bar(stat = "identity") + RotatedAxis()
)

##### NOT USED OLD plot fig2e ####

test.matrix <- DotPlot(final.all.hearts,features = rev(c("nrg1","fn1a","aldh1a2","stra6","cxcl12b","cxcr4a")),group.by = "trashid",split.by = "inhib")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Neuronal cells_NULL","Myelin cells_NULL","Dead cells_NULL",rownames(test.matrix)[grep("IWR1$",rownames(test.matrix))]) ,]
rownames(test.matrix) <- as.character(sapply(rownames(test.matrix),  function(x) unlist(strsplit(x,"_"))[1] ))
test.matrix<- scale(test.matrix)

pheatmap(test.matrix,cluster_cols = F)

# genes of interest after inhib for fig 4 ####
test.matrix <- DotPlot(final.all.hearts,features = c("nrg1","fn1a","stra6","aldh1a2","cxcr4a","cxcl12b"),group.by = "lineage.ident2",split.by = "inhib")
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- test.matrix[!rownames(test.matrix) %in% c("Neuronal cells","Endocardium (frzb)"),]

test.matrix1 <- test.matrix[rownames(test.matrix) %in% paste0(c("Fibroblast (col11a1a)","Fibroblast (col12a1a)","Fibroblast (proliferating)","Fibroblast (nppc)","Fibroblast (cxcl12a)","Epicardium (Atrium)","Epicardium (Ventricle)","Perivascular cells","Bl.ves.EC (plvapb)","Fibroblast (nppc)","Fibroblast (spock3)","Fibroblast-like cells"),"_NULL"),]
test.matrix2 <- test.matrix[rownames(test.matrix) %in% paste0(c("Fibroblast (col11a1a)","Fibroblast (col12a1a)","Fibroblast (proliferating)","Fibroblast (nppc)","Fibroblast (cxcl12a)","Epicardium (Atrium)","Epicardium (Ventricle)","Perivascular cells","Bl.ves.EC (plvapb)","Fibroblast (nppc)","Fibroblast (spock3)","Fibroblast-like cells"),"_IWR1"),]


test.matrix <- cbind(test.matrix1,test.matrix2)
test.matrix <- scale(test.matrix)

plot_grid(pheatmap(test.matrix1,cluster_cols = F,cluster_rows = F)[[4]],pheatmap(test.matrix2,cluster_cols = F,cluster_rows = F)[[4]])
pheatmap(test.matrix,cluster_rows = F,cluster_cols = F)

test.matrix <- DotPlot(final.all.hearts,features = c("nrg1","fn1a","stra6","aldh1a2","cxcr4a","cxcl12b"),group.by = "is.inhib",cols = rep("red",100))
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- as.data.frame(test.matrix)

plot_grid(ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = cxcl12b))+
            geom_bar(stat = "identity" ) + RotatedAxis(),
          ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = cxcr4a))+
            geom_bar(stat = "identity" ) + RotatedAxis(),
          ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = aldh1a2))+
            geom_bar(stat = "identity" ) + RotatedAxis(),
          ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = stra6))+
            geom_bar(stat = "identity" ) + RotatedAxis(),
          ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = fn1a))+
            geom_bar(stat = "identity" ) + RotatedAxis(),
          ggplot(test.matrix,aes(x=factor(rownames(test.matrix),levels = c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")),y = nrg1))+
            geom_bar(stat = "identity" ) + RotatedAxis())



test.matrix <- DotPlot(final.all.hearts,features = c("nrg1","fn1a","stra6","aldh1a2","cxcr4a","cxcl12b"),group.by = "lineage.ident2",split.by = "is.inhib",cols = rep("red",100))
test.matrix <- test.matrix$data
test.matrix <- test.matrix[,c(1,2,3,4)]
test.matrix <- acast(test.matrix, id~features.plot, value.var = "avg.exp")
test.matrix <- as.data.frame(test.matrix)
plot.matrix <- test.matrix[rownames(test.matrix) %in% rownames(test.matrix)[grep("Fibroblast \\(nppc\\)",rownames(test.matrix))],]
plot.matrix <- test.matrix[rownames(test.matrix) %in% rownames(test.matrix)[grep("Perivascular cells",rownames(test.matrix))],]

plot.matrix <- as.data.frame(plot.matrix)
colnames(plot.matrix) <- "value"

plot.matrix2 <- test.matrix[rownames(test.matrix) %in% rownames(test.matrix)[grep("Epicardium \\(Ventricle\\)",rownames(test.matrix))],]
plot.matrix2 <- test.matrix[rownames(test.matrix) %in% rownames(test.matrix)[grep("Fibroblast \\(col11a1a\\)",rownames(test.matrix))],]

plot.matrix2 <- as.data.frame(plot.matrix2)
colnames(plot.matrix2) <- "value"

ggplot(plot.matrix2,aes(x=factor(rownames(plot.matrix),levels = paste0("Perivascular cells_",c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")) ),y = value))+
  geom_bar(stat = "identity" ) +RotatedAxis()

ggplot(plot.matrix2,aes(x=factor(rownames(plot.matrix2),levels = paste0("Fibroblast (col11a1a)_",c("Ctrl","3dpi","3dpiinhib","7dpi","7dpiinhib","30dpi")) ),y = value))+
  geom_bar(stat = "identity" ) +RotatedAxis()



## library info ####
library(dplyr)
lib.info <- final.all.hearts@meta.data[,c(1,2,3)]
lib.info <- group_by(lib.info,orig.ident)
export <- summarise(lib.info,nUMI = round(mean(nCount_RNA)),nGene = round(mean(nFeature_RNA) ))
export <- as.data.frame(export)
write.csv(export,file = "../finalplots/suppTable1.csv")

lib.info <- morphcomp@meta.data[,c(1,2,3)]
lib.info <- group_by(lib.info,orig.ident)
export <- summarise(lib.info,nUMI = round(mean(nCount_RNA)),nGene = round(mean(nFeature_RNA) ))
export <- as.data.frame(export)

# compare diff.genes in different conditions ####
save(final.all.hearts,file = "/data/junker/users/Bo/final_seurat_objs/final.all.hearts.Robj")
final.all.hearts <- SetIdent(final.all.hearts,value = "default")
final.all.hearts <- SetIdent(final.all.hearts,value = "Fibroblast (nppc)",cells = WhichCells(final.all.hearts, expression = lineage.ident == "Fibroblast (nppc)" ))
final.all.hearts <- SetIdent(final.all.hearts,value = "Endocardium (preinjury)",
                             cells = intersect(WhichCells(final.all.hearts, expression = lineage.ident %in% c("Endocardium (frzb)","Endocardium (V)")),
                                    WhichCells(final.all.hearts, expression = is.inhib == "Ctrl")))
final.all.hearts <- SetIdent(final.all.hearts,value = "Endocardium (3dpi)",
                             cells = intersect(WhichCells(final.all.hearts, expression = lineage.ident %in% c("Endocardium (frzb)","Endocardium (V)")),
                                               WhichCells(final.all.hearts, expression = is.inhib == "3dpi")))
final.all.hearts <- SetIdent(final.all.hearts,value = "Endocardium (7dpi)",
                             cells = intersect(WhichCells(final.all.hearts, expression = lineage.ident %in% c("Endocardium (frzb)","Endocardium (V)")),
                                               WhichCells(final.all.hearts, expression = is.inhib == "7dpi")))

final.all.hearts <- StashIdent(final.all.hearts, save.name = "pseudobulk")



nppc.vs.preinjuredEC <- FindMarkers(final.all.hearts,ident.1 = "Fibroblast (nppc)", ident.2 = "Endocardium (preinjury)",logfc.threshold = 0.1)
nppc.vs.dpi3EC <- FindMarkers(final.all.hearts,ident.1 = "Fibroblast (nppc)", ident.2 = "Endocardium (3dpi)",logfc.threshold = 0.1)
dpi3EC.vs.preinjuredEC <- FindMarkers(final.all.hearts,ident.1 = "Endocardium (3dpi)", ident.2 = "Endocardium (preinjury)",logfc.threshold = 0.1)

nppc.vs.dpi7EC <- FindMarkers(final.all.hearts,ident.1 = "Fibroblast (nppc)", ident.2 = "Endocardium (7dpi)",logfc.threshold = 0.1)
dpi7EC.vs.preinjuredEC <- FindMarkers(final.all.hearts,ident.1 = "Endocardium (7dpi)", ident.2 = "Endocardium (preinjury)",logfc.threshold = 0.1)



write.csv(nppc.vs.preinjuredEC,file = "nppc.vs.ctrlEC.csv")
write.csv(nppc.vs.dpi3EC,file = "nppc.vs.dpi3EC.csv")
write.csv(dpi3EC.vs.preinjuredEC,file = "dpi3EC.vs.ctrlEC.csv")


nppc.vs.preinjuredEC.high <- rownames(nppc.vs.preinjuredEC[nppc.vs.preinjuredEC$avg_logFC>0.25,])
nppc.vs.dpi3EC.high <- rownames(nppc.vs.dpi3EC[nppc.vs.dpi3EC$avg_logFC>0.25,])
dpi3EC.vs.preinjuredEC.high <- rownames(dpi3EC.vs.preinjuredEC[dpi3EC.vs.preinjuredEC$avg_logFC>0.25,])

intersect(nppc.vs.preinjuredEC.high,dpi3EC.vs.preinjuredEC.high)

nppc.vs.dpi7EC.high <- rownames(nppc.vs.dpi7EC[nppc.vs.dpi7EC$avg_logFC>0.25,])
dpi7EC.vs.preinjuredEC.high <- rownames(dpi7EC.vs.preinjuredEC[dpi7EC.vs.preinjuredEC$avg_logFC>0.25,])

myCol <- c( "#7e3b65","#7e3b65" ,"#406dac")
venn.diagram( x = list(nppc.vs.preinjuredEC.high,nppc.vs.dpi3EC.high, dpi3EC.vs.preinjuredEC.high),
    category.names = c("Upregulated genes in nppc fibroblasts \n compared to uninjured endocardium" , "Upregulated genes in nppc fibroblasts \n compared to 3dpi endocardium","Upregulated genes in 3dpi endocardium \n compared to uninjured endocardium"),
    filename = 'venn_diagramm2.png',
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 700 , 
    width = 700 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    #lwd = 7,
    lty = 'blank',
    fill = myCol,
    alpha = 0.5,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-15, 15,0),
    cat.dist = c(0.1, 0.1,0.06),
    cat.fontfamily = "sans",
    cat.col = myCol,
    #rotation = 1
)


# npp
