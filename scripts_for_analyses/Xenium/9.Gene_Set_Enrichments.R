library(Seurat)
library(SeuratData)
library(Matrix)
library(ggplot2)
library(patchwork)
library(BPCells)
library(SeuratExtend)
library(dplyr)
options(max.print = 12, spe = "human")
setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")

xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025_Niche-k10.Rds")
xenium.obj <- UpdateSeuratObject(xenium.obj)

#GSEA analysis For Niches 6 (frPT, CKD), 8 (early aPT, AKI) and 9 (healthy PT)
xen.obj.sub <- subset(xenium.obj, niches %in% c("6","8","9"))
xen.obj.sub[["Spatial"]]$counts <- as(xen.obj.sub[["Spatial"]]$counts, Class = "dgCMatrix")
xen.obj.sub <- NormalizeData(xen.obj.sub, assay = "Spatial", normalization.method = "LogNormalize", scale.factor = 10000)
xen.obj.sub <- GeneSetAnalysis(xen.obj.sub, genesets = hall50$human, assay = "Spatial", nCores = 4)
matr <- xen.obj.sub@misc$AUCell$genesets

pdf(file='GSEA/HKAv2_Xenium_Niche-6-8-9_GSEA_Waterfall_Plot.pdf',width=7,height=5)
WaterfallPlot(matr, f = xen.obj.sub$niches, 
              ident.1 = "9",
              ident.2 = "6",
              len.threshold = 2) + DarkTheme()
WaterfallPlot(matr, f = xen.obj.sub$niches, 
              ident.1 = "8",
              ident.2 = "6",
              len.threshold = 2) + DarkTheme()

dev.off()



###CGI-less associated GSEA in PT cells only
table(xenium.obj$v2.subclass.l1)
xen.obj.sub <- subset(xenium.obj, v2.subclass.l1 %in% c("PT"))
xen.obj.sub[["Spatial"]]$counts <- as(xen.obj.sub[["Spatial"]]$counts, Class = "dgCMatrix")
xen.obj.sub <- NormalizeData(xen.obj.sub, assay = "Spatial", normalization.method = "LogNormalize", scale.factor = 10000)
DefaultAssay(xen.obj.sub) <- "Spatial"
xen.obj.sub$v2.subclass.l2[xen.obj.sub$v2.subclass.l2 %in% c("frPT")] <- "aPT"
xen.obj.sub$v2.subclass.l2[xen.obj.sub$v2.subclass.l2 %in% c("PT-S1/2")] <- "PT"
xen.obj.sub$v2.subclass.l2[xen.obj.sub$v2.subclass.l2 %in% c("PT-S3")] <- "PT"
table(xen.obj.sub$v2.subclass.l2)
#aPT      cycPT     PT 
#14330    2738      104212 

genes <- read.csv(file = "DMR_DEG_CGI_GeneSets.csv")
genes <- genes[genes$target %in% rownames(xen.obj.sub),]
#Slit into list by source
gene_list <- split(
  x = genes$target,  
  f = genes$source 
)

#GSEA analysis
xen.obj.sub <- GeneSetAnalysis(xen.obj.sub, genesets = gene_list, assay = "Spatial", nCores = 4)
matr <- xen.obj.sub@misc$AUCell$genesets

Idents(xen.obj.sub) <- "v2.subclass.l2"
pdf(file='GSEA/HKAv2_Xenium_PT_GSEA_Waterfall_Plot_DMR_DEG_CGI_genesets.pdf',width=7,height=4)
WaterfallPlot(matr, f = xen.obj.sub$v2.subclass.l2, 
              ident.1 = "aPT",
              ident.2 = "PT",
              len.threshold = 2) + DarkTheme()
dev.off()
pdf(file='GSEA/HKAv2_Xenium_PT_GSEA_Waterfall_Plot_DMR_DEG_CGI_genesets_B.pdf',width=5,height=2)
WaterfallPlot(matr, f = xen.obj.sub$v2.subclass.l2, 
              ident.1 = "aPT",
              ident.2 = "PT",
              len.threshold = 2) #+ DarkTheme()
dev.off()

pdf(file='GSEA/HKAv2_Xenium_PT_GSEA_Waterfall_Plot_Top20-DMR_DEG_CGI_genes.pdf',width=7,height=5)
WaterfallPlot(
  xen.obj.sub, 
  group.by = "v2.subclass.l2", 
  features = genes$target,
  ident.1 = "aPT", 
  ident.2 = "PT", 
  length = "logFC",
  log.base = "2",
  top.n = 20) + DarkTheme()
dev.off()
pdf(file='GSEA/HKAv2_Xenium_PT_GSEA_Waterfall_Plot_Top20-DMR_DEG_CGI_genes_B.pdf',width=6,height=3)
WaterfallPlot(
  xen.obj.sub, 
  group.by = "v2.subclass.l2", 
  features = genes$target,
  ident.1 = "aPT", 
  ident.2 = "PT", 
  length = "logFC",
  log.base = "2",
  top.n = 20) #+ DarkTheme()
dev.off()

