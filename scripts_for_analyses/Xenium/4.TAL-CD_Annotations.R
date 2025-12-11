library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(pagoda2)
library(future)
plan("multisession", workers = 10)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB


###TAL-CD Clusters
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection.Rds")
cl.order <- c(
  1, # MD
  3, # C-TAL
  6, # IC-B
  9, # CNT-PC
  11, # aCNT
  12, # DCT1
  15, # C-TAL
  18, # CNT
  20, # PC
  22, # aTAL
  24, # frDCT
  34, # frTAL
  39, # OMCD-IC-A
  40, # C-TAL
  46, # aTAL1
  52, # IMCD
  57, # M-TAL
  59, # aCNT
  60, # PapE
  61 # PapE
  
)

KB <- subset(KB, idents = cl.order)

DefaultAssay(KB) <- "RNA"
KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- SketchData(
  object = KB,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

##Integrate sketch data set
DefaultAssay(KB) <- "sketch"

KB <- NormalizeData(KB)
KB <- FindVariableFeatures(KB)
KB <- ScaleData(KB)
KB <- RunPCA(KB)

KB <- IntegrateLayers(
  object = KB, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca", k.anchor = 20,
  verbose = FALSE
)
KB <- RunUMAP(KB, reduction = "integrated.rpca", dims = 1:50)

DimPlot(
  KB,
  reduction = "umap",
  group.by = c("assay"),
  combine = TRUE, label.size = 2
)
DimPlot(
  KB,
  reduction = "umap",
  group.by = c("rpca_k100_infomap"),
  combine = TRUE, label.size = 2
)

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_TAL-CD.Rds")



###Clustering on sketch data set
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_TAL-CD.Rds")

##Join Layers and Annotate
KB[["sketch"]] <- JoinLayers(KB[["sketch"]])
KB[["sketch"]] <- as(object = KB[["sketch"]], Class = "Assay")
KB <- NormalizeData(KB)


###cluster using Pagoda2 instead
countMatrix <- KB[["sketch"]]$counts
dim(countMatrix)
#36588 100000

# Generate a new pagoda2 object
p2 <- Pagoda2$new(x = countMatrix, n.cores = 6, trim=10)

# Add in integrated PCA embeddings
p2$reductions$PCA <- Embeddings(object = KB[["integrated.rpca"]])

# Generate K-nearest neighbour graph
p2$makeKnnGraph(k = 200, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k200infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k200infomap, col.name = "rpca_k200_infomap")
DimPlot(KB, group.by = "rpca_k200_infomap", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()

##Project to full data set
meta <- KB@meta.data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_TAL-CD.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k200_infomap"], col.name = "rpca_k200_infomap_dt")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k200_infomap_dt = "rpca_k200_infomap_dt"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "v2.subclass.l3", alpha = 0.1)
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k200_infomap_dt", label = TRUE,
        alpha = 0.1) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_TAL-CD.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_TAL-CD.Rds")


#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
levels(Idents(object = KB)) <- paste("CL", levels(Idents(object = KB)), sep = "")
celltype <- Idents(object = KB)
table(Idents(KB))

getMaxIdents <- function(seurat.obj, celltype, ref.column, query.column) {
  max.idents <- do.call(rbind, lapply(levels(celltype), function(ct) {
    print(paste("Running for cell type:", ct))
    cells <- WhichCells(seurat.obj, idents = ct)
    seurat.obj.sub <- seurat.obj@meta.data[cells,]
    top2 <- function(x) {
      if (length(x) < 2) x <- c(x, rep(NA, 2 - length(x)))
      x
    }
    ref.table <- table(seurat.obj.sub[[ref.column]])
    query.table <- table(seurat.obj.sub[[query.column]])
    if (length(ref.table) == 0 || length(query.table) == 0) {
      return(data.frame(ref.cluster=NA, ref.Freq=NA, ref.Total=0, query.subclass=NA, query.Freq=NA, query.Total=0, int.cluster=ct))
    }
    ref.top2 <- top2(tail(names(sort(ref.table)), 2))
    query.top2 <- top2(tail(names(sort(query.table)), 2))
    Idents.called <- data.frame(
      ref = prop.table(ref.table)[ref.top2],
      ref.Total = sum(ref.table),
      query = prop.table(query.table)[query.top2],
      query.Total = sum(query.table)
    )
    if (any(is.na(Idents.called))) {
      Idents.called[is.na(Idents.called)] <- NA
    }
    Idents.called$cluster <- ct
    names(Idents.called) <- c("ref.cluster", "ref.Freq", "ref.Total", "query.subclass", "query.Freq", "query.Total", "int.cluster")
    Idents.called
  }))
  return(max.idents)
}
max.idents <- getMaxIdents(KB, celltype, "v2.subclass.l1", "v2.subclass.l3")
max.idents


##Check marker genes and annotate
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_TAL-CD.Rds")
Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))

KB <- JoinLayers(KB)
KB <- NormalizeData(KB)
table(KB$assay)
xenium.obj <- subset(KB, assay %in% c("Xenium 5K"))
table(Idents(xenium.obj))

##Check marker genes
epi.markers <- c("CASR","SLC12A1","UMOD",#"EGF","ESRRB",                  #TAL
                 "ANK2","CLCNKA","KCTD16","BMP6","RIPOR3","CLDN14",      #M-TAL
                 "PHACTR1","SLCO3A1",#"CXCL12","CNTN1","CABP1",           #TAL-A
                 "KCNMB2",#"RGS6",                                        #C/M-TAL-A
                 "ENOX1","CALCR",#"RBM20","PDE3A",                        #C-TAL-A
                 #"DACH1","LRMDA",                                        #TAL-B
                 
                 "TENM4","FGF14",#"PXDNL","GRM8",                         #C/M-TAL-B
                 "KCNJ10","TMEM52B","CLDN16","TMEM207","JAG1",           #TAL-B
                 "SCN7A","COL8A1","LINGO2",                              #C-TAL-B                          
                 "BBOX1","NOS1","ROBO2","TMPRSS4",                       #MD
                 
                 "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR",              #aTAL1
                 "CD44","DIAPH3","SYT16","HAVCR1",                       #aTAL1                                      
                 
                 "ITGB6","NRP1","TFPI",                                  #aTAL
                 
                 "HIF1A","ADAMTS1","DNAH5",                              #aTAL2
                 
                 "ITGB8","PROM1","ARHGAP26",#"RNF144B","TMPRSS4","RHEX",  #frTAL
                 "EZH2","CENPP","MKI67",#"TOP2A",                        #cycling
                 
                 "SLC12A3","CNNM2","KLHL3","TRPM6",                       #DCT
                 "GRAMD1B","ADAMTS17","ZNF385D",                          #DCT1
                 
                 "UNC5C","HMCN1","FSTL4","TRPV5",                         #DCT2
                 "CACNA1D", "ACSL4",                                      #frDCT
                 "FGF13","IGF2BP2","FAM155A","NRG1",                      #aDCT 
                 
                 "SLC8A1",                                                #DCT2 / CNT
                 "HSD11B2","CALB1","ANGPT1","SNTG1",                      #CNT
                 #"CTSD",                                                  #dCNT
                 "RAPGEF5","DLGAP1","BIRC3",                              #aCNT
                 
                 "GATA3","AQP2","PAPPA",                                  #PC
                 "SCNN1G","SCNN1B",
                 "SGCD","STC1",                                           #CNT-PC
                 "FGF12","PTCSC3","CNTNAP3B",                             #CCD/OMCD PC
                 "MCTP1","CPAMD8","GABRA2",                               #OMCD-PC
                 "GREB1L",                                                #OMCD-PC/IMCD
                 "SLC14A2","HS3ST5","TENM3",#"TGFBR3","AKR1C1","FAT3",    #IMCD
                 #"AOC1","TFPI2","LCN2",                                   #dIMCD
                 "ADIRF","SNCG","DHRS2","GPX2","TP63",#"FXYD3",           #PapE                              #PapE
                 
                 "ATP6V0D2", "ATP6V1C2", "CLNK",                          #IC
                 "SLC26A7", "SLC4A1",                                     #IC-A                                   
                 "HS6ST3","NRXN3", "NXPH2", "LEF1",                       #CCD-IC-A
                 
                 "FAM184B","ADTRP","AQP6","STAP1",                        #OMCD-IC-A
                 
                 "SLC4A9", "SLC26A4", "INSRR", "TLDC2",                #IC-B
                 
                 "CKB","COX8A",#"PEBP1","UQCRB",                         #Injury
                 "CST3","DEFB1","SPP1","IGFBP7","CLU"
)     


DotPlot(KB, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()
#DotPlot(KB, features = unique(int.markers.1), dot.scale = 8) + RotatedAxis()
#DotPlot(KB, features = unique(int.markers.2), dot.scale = 8) + RotatedAxis()
DotPlot(stereo.obj, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()

DimPlot(KB, reduction = "umap", group.by = "v2.subclass.l3", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "rpca_k200_infomap_dt", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "assay", label = TRUE,
        alpha = 0.1) + NoLegend()

CL47_markers <- FindMarkers(object = xenium.obj, ident.1 = 47, 
                            max.cells.per.ident = 500, only.pos = TRUE)


cl.order <- c(
  36, #M-TAL
  24, #M-TAL
  13, #C/M-TAL
  7, #TAL
  16, #C-TAL
  6, #C-TAL
  15, #aTAL1
  20, #aTAL1
  1, #aTAL2
  21, #aTAL2
  40, #aTAL2
  5, #frTAL
  30, #DCT
  22, #DCT
  14, #DCT
  25, #aDCT
  29, #aDCT
  11, #CNT
  28, #CNT
  35, #CNT
  4, #aCNT
  10, #aCNT
  23, #aCNT
  9, #CCD-PC
  17, #CCD-PC
  39, #CCD-PC
  8, #CCD-PC
  33, #OMCD-PC
  34, #OMCD-PC
  32, #IMCD
  26, #PapE
  37, #PapE
  38, #PapE
  2, #CCD-IC-A
  19, #CCD-IC-A
  27, #CCD-IC-A
  31, #OMCD-IC-A
  3, #IC-B
   
  12, #ambiguous
  18 #ambiguous
  
  )

Idents(KB) <- "rpca_k200_infomap_dt"
Idents(KB) <- factor(Idents(KB), levels = cl.order)
DotPlot(KB, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()
Idents(xenium.obj) <- "rpca_k200_infomap_dt"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = cl.order)
DotPlot(xenium.obj, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()


#Spatial localization
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
kss <- AddMetaData(kss, metadata = KB$rpca_k200_infomap_dt, col.name = "rpca_k200_infomap_dt")
Idents(kss) <- "rpca_k200_infomap_dt"

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(36, #M-TAL
                                                24 #M-TAL
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(13 #C/M-TAL
                                                
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(7 #TAL
                                                
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(16, #C-TAL
                                                6 #C-TAL
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(15, #aTAL1
                                                20 #aTAL1
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(1, #aTAL2
                                                21, #aTAL2
                                                40 #aTAL2
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(5 #frTAL
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(30, #DCT
                                                22, #DCT
                                                14 #DCT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(25, #aDCT
                                                29 #aDCT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(11, #CNT
                                                28, #CNT
                                                35 #CNT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(4, #aCNT
                                                10, #aCNT
                                                23 #aCNT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(9, #CCD-PC
                                                17, #CCD-PC
                                                39, #CCD-PC
                                                8 #CCD-PC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(33, #OMCD-PC
                                                34 #OMCD-PC
                                                )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(32 #IMCD

             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(26, #PapE
                                                37, #PapE
                                                38 #PapE
                                                )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(2, #CCD-IC-A
                                                19, #CCD-IC-A
                                                27 #CCD-IC-A
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(31 #OMCD-IC-A

             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(3 #IC-B
                                                )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(18 #ambiguous
             )))


#save metadata
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_TAL-CD_Metadata.Rds")


