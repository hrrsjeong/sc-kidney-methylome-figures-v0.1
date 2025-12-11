library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(pagoda2)
library(future)
plan("multisession", workers = 10)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB


###PT-TL Clusters
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection.Rds")
cl.order <- c(
  2, # frPT-S1/S2
  16, # aPT-S3
  19, # PT
  21, # POD
  25, # cycPT
  28, # PT-S1
  35, # PEC
  37, # PT-S1/S2
  45, # PT-S3
  54, # ATL
  55, # PT-S1
  58 # DTL2
  
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

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_PT-TL.Rds")



###Clustering on sketch data set
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_PT-TL.Rds")

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
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_PT-TL.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k200_infomap"], col.name = "rpca_k200_infomap_pt")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k200_infomap_pt = "rpca_k200_infomap_pt"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "v2.subclass.l3", alpha = 0.1)
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k200_infomap_pt", label = TRUE,
        alpha = 0.1) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_PT-TL.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_PT-TL.Rds")


#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap_pt"
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
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_PT-TL.Rds")
Idents(KB) <- "rpca_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))

KB <- JoinLayers(KB)
KB <- NormalizeData(KB)
table(KB$assay)
xenium.obj <- subset(KB, assay %in% c("Xenium 5K"))
table(Idents(xenium.obj))

##Check marker genes
epi.markers <- c("NPHS1","NPHS2",#"ST6GALNAC3","PODXL",                         #POD
                 "PTGDS","BST2","TPPP3","CHI3L1",                              #dPOD
                 "ALDH1A2","CFH",#"FAM155A","CLDN1",                            #PEC
                 
                 "LRP2","CUBN",                                                #PT
                 "SLC5A12","SLC22A6",#"HNF4A",                                  #S1/S2
                 "RALYL","PCDH15","PRODH2",#"SLC22A8",                          #S1
                 "SLC34A1","ANKS1B","SLC5A10","SLC5A11",                       #S2                                  
                 "SLC5A8","GPM6A","SLC22A24","SLC7A13",                        #PT-S3
                 
                 "IGFBP7","SPP1","ITGB8","CDH6","TMEM178B","ALPK2","HAVCR1",
                 "ITGB3",
                 "CST3","CLU","VIM","PIGR",#"APOE",                            #aPT2
                 "IL32","SOX4","VCAM1",#"MMP7","SOX9","CCL2",                   #aPT2
                 "DCC",                                                        #aPT1
                 "GDA","GLIS1",                                                #aPT-S1/S2
                 "DLGAP1","PROM1",                                             #frPTS1/S2
                 "APBB1IP","ROBO2","COL23A1","MEG3",                           #frPTS1/S2
                 "LSAMP","KCNIP1","NRXN3","WNT2B",                             #frPTS1/S2
                 "KCTD16","SPON1",                                             #aPT-S3
                 "NRG3","FAM189A1","DTNA","KITLG","GRM8",                      #frPT-S3
                 "TOP2A","MKI67",                                              #cycling
                 "EGR1","FOS",                                                 #dPT
                 
                 "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",                       #TL
                 "ABCA13",                                                     #DTL
                 "AQP1", "UNC5D","LRRC4C","DENND2A","SLCO1A2",                 #DTL2
                 "IRX3", "SERPINA1",                                           #aDTL2
                 "SIM2",                                                       #DTL1/3/ATL
                 "ADGRL3","JAG1","SMAD9","ID1",                                #DTL1
                 "SLC14A2","FAM169A","SMOC2",                                  #DTL3
                 "ABCA4","BCL6","AKR1B1","SH3GL3",                             #DTL3/ATL
                 "CLCNKA","PROX1","CACNA1C","COLEC10",                         #ATL
                 "SOD3","HSPA2","CSRP2",                                       #dATL
                 "CST3","APOE","GATM","ALDOB","CLU","DEFB1",                   #Injury
                 "SPP1","IGFBP7","CLU"
)     


DotPlot(KB, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()
#DotPlot(KB, features = unique(int.markers.1), dot.scale = 8) + RotatedAxis()
#DotPlot(KB, features = unique(int.markers.2), dot.scale = 8) + RotatedAxis()
DotPlot(stereo.obj, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()

DimPlot(KB, reduction = "umap", group.by = "v2.subclass.l3", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "rpca_k200_infomap_pt", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "assay", label = TRUE,
        alpha = 0.1) + NoLegend()

#CL47_markers <- FindMarkers(object = xenium.obj, ident.1 = 47, 
#                            max.cells.per.ident = 500, only.pos = TRUE)


cl.order <- c(
  5, #POD
  14, #POD
  22, #POD
  37, #POD
  6, #PEC
  
  1, # PT
  33, #PT
  38, #PT
  8, #PT-S1/S2
  10, #PT-S1/S2
  18, #PT-S1/S2
  12, #PT-S1/S2
  7, #PT-S3
  35, #aPT2
  25, #aPT2
  9, #aPT
  15, #aPT
  36, #aPT
  17, #aPT
  3, #aPT
  13, #aPT
  2, #frPT
  4, #frPT
  11, #frPT
    
  16, #cycPT
  
  23, #DTL2
  20, #DTL1
  30, #DTL1
  27, #DTL3
  32, #DTL3
  26, #ATL
  28, #ATL
  
  24, #cycEC
  
  29, #ambiguous (too few cells)
  31, #ambiguous
  19, #ambiguous
  21, #FIB multiplet
  34 #ambiguous
  
)

Idents(KB) <- "rpca_k200_infomap_pt"
Idents(KB) <- factor(Idents(KB), levels = cl.order)
DotPlot(KB, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()
Idents(xenium.obj) <- "rpca_k200_infomap_pt"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = cl.order)
DotPlot(xenium.obj, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()


#Spatial localization
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
kss <- AddMetaData(kss, metadata = KB$rpca_k200_infomap_pt, col.name = "rpca_k200_infomap_pt")
Idents(kss) <- "rpca_k200_infomap_pt"

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(5, #POD
                                                14, #dPOD
                                                22, #dPOD
                                                37 #dPOD
             )))
ImageDimPlot(kss, cols = c("red"), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(6 #PEC
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(1, # PT
                                                33, #PT
                                                38 #PT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(8, #PT-S1/S2
                                                10, #PT-S1/S2
                                                18, #PT-S1/S2
                                                12 #PT-S1/S2
                                                
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(7 #PT-S3
                                                
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(35, #aPT2
                                                25 #aPT2
                                                             )))

ImageDimPlot(kss, cols = rep("red",6), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(9, #aPT
                                                15, #aPT
                                                36, #aPT
                                                17, #aPT
                                                3, #aPT
                                                13 #aPT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(2, #frPT
                                                4, #frPT
                                                11 #frPT
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(16 #cycPT
            )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(23 #DTL2

             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(20, #DTL1
                                                30 #DTL1
                                                
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(29 #dDTL
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(27, #DTL3
                                                32 #DTL3
                                                
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(26, #ATL
                                                28 #ATL
                                                
                                                
             )))


ImageDimPlot(kss, cols = rep("red",5), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(31, #ambiguous
                                                19, #ambiguous
                                                21, #FIB multiplet?
                                                24, #cycEC
                                                34 #ambiguous?
                                                
                                                
                                                
             )))


#save metadata
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_PT-TL_Metadata.Rds")

