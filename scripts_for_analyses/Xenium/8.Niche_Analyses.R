library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(Polychrome)
library(dplyr)
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(Seurat.object.assay.version = "v5")
load("~/Projects/Human_Kidney/Atlas_V2/color_factors_v2-clusters.robj")


###Spatial image plots
xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025.Rds")
order <- c("POD","PEC","PT","PT-S1/2","PT-S3","aPT2","aPT1","frPT","cycPT","DTL2",
           "DTL1","DTL3","ATL","M-TAL","C/M-TAL","TAL","C-TAL","aTAL1","aTAL2",
           "frTAL","DCT","aDCT","CNT","aCNT","CCD-PC","OMCD-PC","IMCD","PapE",
           "CCD-IC-A","OMCD-IC-A","IC-B","EC-GC","EC-AA","EC-DVR","M-EC-PTC",
           "EC-GC/PTC","EC-PTC","infEC-PTC","EC-EA","EC-V","EC-AVR","EC-LYM",
           "cycEC","IM-FIB","OM-FIB","C-FIB","FIB","C-FIB-PATH","C-MYOF","pvFIB-RSPO3+",
           "pvFIB-PI16+","pvMYOF","MC","VSMC","VSMC/P","B","PL","T","NK",
           "MAST","resMAC-LYVE1+","MAC","MON","cDC2","pDC","cDC1","cycIMM")
Idents(xenium.obj) <- "v2.subclass.l3"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = order)
table(Idents(xenium.obj))



epi.cols <- setNames(c("#F1085C","#783FC1","#0000FF","#B1CC71","#005300","#FFFFFF",
                       "#886C00","#FFD300","#14F9FF"),
                     c("POD","PEC","PT-S1/2","PT-S2","PT-S3","DTL1",     
                       "DTL2","DTL3","ATL"))

ImageDimPlot(xenium.obj, fov = "fov",
             cells = WhichCells(xenium.obj, 
                                idents = c("POD","PEC","PT-S1/2","PT-S3","DTL1",     
                                           "DTL2","ATL")),
             cols = epi.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##niches.k = 10
xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov",
                              group.by = "v2.subclass.l3",
                              niches.k = 10, neighbors.k = 25)

Idents(xenium.obj) <- "niches"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = c(1:10))
table(Idents(xenium.obj))

saveRDS(xenium.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.Rds")
#xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.Rds")
write.table(xenium.obj@meta.data, file="Human_Kidney_Xenium_Full_Metadata_Niches-k10_04142025.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


##niches.k = 12
xenium.obj <- BuildNicheAssay(object = xenium.obj, fov = "fov",
                              group.by = "v2.subclass.l3",
                              niches.k = 12, neighbors.k = 25)

Idents(xenium.obj) <- "niches"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = c(1:12))
table(Idents(xenium.obj))

saveRDS(xenium.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k12.Rds")
#xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k12.Rds")
write.table(xenium.obj@meta.data, file="Human_Kidney_Xenium_Full_Metadata_Niches-k12_04142025.txt", sep = "\t", row.names=TRUE, col.names=TRUE)


###Visualizations with k = 10 niches
xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.Rds")

Glasbey = glasbey.colors(21)
swatch(Glasbey)

niche.cols <- setNames(c("#442288", "#71c8a5", "#8C4F7D", "#FED23F", "#FFCC99", "#d10040", "#B5D33D", "#005C31", "#1036b5","#FF5733"),
                       c(1:10))

ImageDimPlot(xenium.obj, fov = "fov",  
             size = 0.5, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = niche.cols)
table(xenium.obj$niches,xenium.obj$v2.subclass.l3)


p1 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(1)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 1") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p2 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(2)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 2") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p3 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(3)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 3") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p4 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(4)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 4") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(5)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 5") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p6 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(6)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 6") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p7 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(7)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 7") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p8 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(8)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 8") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p9 <- ImageDimPlot(xenium.obj, fov = "fov",  
                   size = 0.05, dark.background = T,
                   cells = WhichCells(xenium.obj, idents = c(9)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 9") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p10 <- ImageDimPlot(xenium.obj, fov = "fov",  
                    size = 0.05, dark.background = T,
                    cells = WhichCells(xenium.obj, idents = c(10)), border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niche 10") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file="Niche_Plots/Xenium_Individual_Niche_Plots.pdf", width=8,height=8)
(p2 | p3 ) / (p5 | p6) / (p7 | p8) / (p9 | p10)
dev.off()


pdf(file="Niche_Plots/Xenium_Combined_Niche_Plots.pdf", width=8,height=5)
ImageDimPlot(xenium.obj, fov = "fov",  
             size = 0.05, dark.background = T, border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niches") + NoLegend() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
pdf(file="Niche_Plots/Xenium_Combined_Niche_Plots_B.pdf", width=8,height=5)
ImageDimPlot(xenium.obj, fov = "fov",  
             size = 0.05, dark.background = T, border.size = NA) +
  scale_fill_manual(values = niche.cols) + ggtitle("Niches") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


# Niche by Cell Type
ct_var <- xenium.obj$v2.subclass.l3
niche_var <- xenium.obj$niches

a <- table(ct_var, niche_var) %>%
  as.data.frame()

pdf(file="Niche_Plots/Xenium_Niche_Composition_subclassl3.pdf", width=8,height=5)
a <- table(ct_var, niche_var) %>%
  as.data.frame() %>%
  ggplot(aes(x = niche_var, y = Freq, fill = ct_var)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values = v2.scl3.cols) +
  labs(x = "Cell Niche", y = "Cell Type Proportion", fill = "", title = "niche_name_var") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))
print(a)
dev.off()

# Niche by Sample
sample_var <- xenium.obj$condition_level1
niche_var <- xenium.obj$niches
pdf(file="Niche_Plots/Xenium_Niche_Distribution_Condition.pdf", width=4,height=5)
a <- table(sample_var, niche_var) %>%
  as.data.frame() %>%
  ggplot(aes(x = sample_var, y = Freq, fill = niche_var)) +
  geom_col(position = position_fill()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = niche.cols) +
  theme_classic(base_size = 14) +
  labs(x = "", y = "Niche Proportions", title = "niche_name_var") +
  theme(legend.title = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        plot.title = element_text(hjust = 0.5))
print(a)
dev.off()




# Cell type proportions across niches
path_ct_prop_df <- as.data.frame(prop.table(table(xenium.obj$v2.subclass.l3,xenium.obj$niches), 1))
names(path_ct_prop_df) <- c("CT", "Niche", "Proportion")

ct_order_4plot <- c(
  "OMCD-IC-A","OMCD-PC",
  "aTAL1","aTAL2","C-TAL",
  "M-TAL","C/M-TAL","TAL","DTL2","DTL1","M-EC-PTC","frTAL","OM-FIB","EC-AVR",
  "IM-FIB","IMCD","ATL","DTL3","PapE",
  "IC-B","CNT","aCNT","CCD-PC","CCD-IC-A",
  "pvFIB-RSPO3+","B","MAST","MAC","pDC","cDC1","PL","pvFIB-PI16+","T","cDC2","frPT","C-MYOF","MON","resMAC-LYVE1+","NK","cycIMM",
  "EC-AA","VSMC","pvMYOF","EC-DVR","VSMC/P","infEC-PTC","EC-LYM",
  "cycPT","aPT2",
  "DCT","aPT1","PT","PT-S1/2","PT-S3","EC-PTC","C-FIB-PATH","aDCT","C-FIB","FIB","EC-GC/PTC","EC-V",
  
  "EC-GC","POD","MC","cycEC","PEC","EC-EA"
                    
                    
)
color_list_4plot <- v2.scl3.cols
theme_angle <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

path_ct_prop_plot <- path_ct_prop_df %>%
  filter(Niche %in% c(1:10)) %>%
  mutate(Niche = ordered(Niche, levels = 1:10)) %>%
  mutate(CT = ordered(CT, levels = ct_order_4plot)) %>%
  filter(Proportion >= 0.01) %>%
  ggplot(aes(x = reorder(CT, dplyr::desc(CT)),
             y = Niche,
             size = Proportion, fill = CT)) +
  geom_point(shape = 21, color = "black", stroke = 0.2) +
  theme_bw(base_size = 6.5) +
  scale_fill_manual(values = color_list_4plot) +
  theme(axis.text = element_text(color = "black", size = 5),
        strip.text = element_text(color = "white", size = 5),
        axis.text.y = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
        plot.title = element_text(hjust = 0.5, size = 6),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box.background = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.background = element_blank(),
        plot.margin = margin(5.5, 5.5, 8, 8, "pt"),
        legend.margin = margin(1, 1, 1, 1, "pt"),
        legend.key.width = unit(2, "mm"),
        legend.key.height = unit(2, "mm")) +
  theme_angle +
  labs(title = "Cell Type Proportions Across Niches") +
  guides(fill = "none",
         size = guide_legend(title.theme = element_text(size = 5),
                             size = 2,
                             nrow = 1, ncol = 6, title.hjust = 0.4,
                             label.position = "bottom",
                             title.position = "top",
                             label.theme = element_text(size = 5))) +
  scale_size_continuous(limits = c(0, 1), range = c(0, 4), breaks = seq(0.2, 1, 0.2)) +
  facet_grid(Niche~., space = "free", scales = "free")


pdf(file="Niche_Plots/Xenium_Cell_Type_Proportions_across_Niches.pdf", width=9,height=5)
path_ct_prop_plot
dev.off()



####Save h5ad Object
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(future)
plan("multisession", workers = 10)
library(anndataR)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.Rds")

DefaultAssay(xenium.obj) <- "Spatial"
xenium.obj[["Spatial"]]$counts <- as(xenium.obj[["Spatial"]]$counts, "dgCMatrix")
xenium.obj[["niche"]] <- NULL
xenium.obj[["X3782"]] <- NULL
xenium.obj[["X3723"]] <- NULL
xenium.obj[["X3916"]] <- NULL
xenium.obj[["X3778"]] <- NULL

xenium.obj$segmentation_method <- NULL
coords <- GetTissueCoordinates(xenium.obj, image = "fov")

Xadata <- from_Seurat(
  xenium.obj, x_mapping = NULL,
  layers_mapping = list(counts = "counts")
)
Xadata$obsm$spatial <- cbind(array(coords$x),array(coords$y))
Xadata$var_names <- rownames(xenium.obj)
Xadata$write_h5ad("Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.h5ad")
#mv to ~/hsKidAt/blake_LTS/Atlas_V2/xenium/
