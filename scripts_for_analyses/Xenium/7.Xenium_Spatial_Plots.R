library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(Polychrome)

options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(Seurat.object.assay.version = "v5")
load("~/Projects/Human_Kidney/Atlas_V2/color_factors_v2-clusters.robj")

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


Glasbey = glasbey.colors(32)
swatch(Glasbey)

levels(Idents(xenium.obj))
epi.cols <- setNames(c("#F1085C","#783FC1","#0000FF","#005300","#FFFFFF",
                       "#886C00","#FFD300","#14F9FF"),
                     c("POD","PEC","PT-S1/2","PT-S3","DTL1",     
                       "DTL2","DTL3","ATL"))


Idents(xenium.obj) <- "v2.subclass.l3"

pdf("Image_Plots/Epi_PT-TL_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("POD","PEC","PT-S1/2","PT-S3","DTL1",     
                                                       "DTL2","DTL3","ATL")),
             cols = epi.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


epi.cols <- setNames(c("#14F9FF","gray20","gray60","#FFACFD","#FF00B6"),
                     c("aPT2","PT-S1/2","PT-S3","aPT1","frPT"))

Idents(xenium.obj) <- "v2.subclass.l3"

pdf("Image_Plots/Epi_aPT_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("aPT2","PT-S1/2","PT-S3","aPT1","frPT"
                                            )),
             cols = epi.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


epi.cols <- setNames(c("#F1085C","#783FC1","#0000FF","#B1CC71","#FFFFFF",
                       "#886C00","#FFD300","#14F9FF","#FFACFD"),
                     c("M-TAL","C/M-TAL","C-TAL","DCT","CNT",     
                       "CCD-PC","OMCD-PC","IMCD","PapE"))

pdf("Image_Plots/Epi_TAL-CD_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("M-TAL","C/M-TAL","C-TAL","DCT","CNT",     
                                                       "CCD-PC","OMCD-PC","IMCD","PapE")),
             cols = epi.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


ec.cols <- setNames(c("#F1085C","#FFFFFF","#0000FF","#005300","#720055","#FFACFD",
                      "#FFD300","#14F9FF"),
                    c("EC-GC","EC-AA","EC-DVR","EC-PTC",
                      "EC-AVR","infEC-PTC","P-EC-AVR","EC-LYM"))

pdf("Image_Plots/EC_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("EC-GC","EC-AA","EC-DVR","EC-PTC",
                                                       "EC-AVR","infEC-PTC","EC-LYM")),
             cols = ec.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


str.cols <- setNames(c("#886C00","#0000FF","#14F9FF","#005300","#B1CC71","#FFFFFF",
                       "#F1085C","#FFACFD","#FFD300","#783FC1","#FFACFD"),
                     c("IM-FIB","OM-FIB","pvFIB-RSPO3+",
                       "C-FIB","pvFIB-PI16+","pvMYOF","MC","REN","VSMC/P",
                       "VSMC","Ad"))

pdf("Image_Plots/STR_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("IM-FIB","OM-FIB","pvFIB-RSPO3+",
                                                       "C-FIB","pvFIB-PI16+","pvMYOF","MC","VSMC/P",
                                                       "VSMC")),
             cols = str.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


rc.cols <- setNames(c("#F1085C","#783FC1","#0000FF","#FFFFFF","#B1CC71","#14F9FF","#005300"),
                    c("POD","PEC","MC","REN","VSMC","EC-GC","EC-AA"))

pdf("Image_Plots/RC_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("POD","PEC","MC","VSMC","EC-GC","EC-AA")),
             cols = rc.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


imm.cols <- setNames(c("#FFD300","white","#F1085C","#720055","#005300","#FFACFD","#886C00"),
                     c("B","PL","T","NK","resMAC-LYVE1+","MAC","cDC2"))

pdf("Image_Plots/Imm_Plot_04012025.pdf", width = 16, height = 8)
ImageDimPlot(xenium.obj, fov = "fov", cells = WhichCells(xenium.obj, 
                                            idents = c("B","PL","T","NK","resMAC-LYVE1+","MAC","cDC2")),
             cols = imm.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

