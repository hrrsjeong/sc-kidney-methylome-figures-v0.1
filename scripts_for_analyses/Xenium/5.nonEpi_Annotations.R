library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(pagoda2)
library(future)
plan("multisession", workers = 10)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB


###Non-Epi Clusters
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection.Rds")
cl.order <- c(
  4, # T
  5, # resMAC-LYVE1+
  7, # C-FIB
  8, # cDC2
  10, # B
  13, # EC-LYM
  14, # PL
  17, # ncMON
  23, # cDC1
  26, # cycT
  27, # C-FIB
  30, # C-EC-PTC
  31, # pvFIB-RSPO3+
  32, # EC-AVR
  33, # EC-GC
  36, # EC-PTC
  38, # C-FIB
  41, # VSMC
  42, # MC
  43, # EC-DVR
  44, # EC-PTC
  47, # pvMYOF
  48, # MAST
  49, # OM-FIB
  50, # VSMC/P
  51, # moFAM
  53, # IM-FIB
  56, # OM-FIB
  62 # cycMAC
  
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

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_nonEpi.Rds")



###Clustering on sketch data set
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_nonEpi.Rds")

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
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_nonEpi.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k200_infomap"], col.name = "rpca_k200_infomap_nonEpi")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")


KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k200_infomap_nonEpi = "rpca_k200_infomap_nonEpi"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k200_infomap_nonEpi"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "v2.subclass.l3", alpha = 0.1)
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k200_infomap_nonEpi", label = TRUE,
        alpha = 0.1) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_nonEpi.Rds")
#KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_nonEpi.Rds")


#Identify max predicted subclass per cluster
Idents(KB) <- "rpca_k200_infomap_nonEpi"
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
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_nonEpi.Rds")
Idents(KB) <- "rpca_k200_infomap_nonEpi"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))
table(Idents(KB))

KB <- JoinLayers(KB)
KB <- NormalizeData(KB)
table(KB$assay)
xenium.obj <- subset(KB, assay %in% c("Xenium 5K"))
table(Idents(xenium.obj))

##Check marker genes
int.markers.1 <- c(
  "PECAM1","PTPRB","FLT1",                                   #Broad EC
  "EMCN","HECW2","ITGA8",                                    #EC-GC
  "EHD3","RXFP1", "MGP","SOST",                              #EC-GC      
  #'PTCHD4',"ZMAT3",                                          #aEC-GC
  #"AQP1","FILIP1","H19","ESM1","SLC45A4",                    #EC-GC-FILIP1+
  
  "PDE3A","SULF1","NKAIN2","NOS1",                           #EC-AA
  "ADAMTS6","MCTP1","PALMD","SLC14A1",#"ITIH5",               #EC-DVR  
  #"LYPD6B","EPHA3","STXBP6","CP",
  
  "CEACAM1","PLVAP","DNASE1L3",                              #PTC/AVR
  "COL15A1","FABP5","ALPL", #"CD36",                         #M-EC-PTC
  "GPM6A","NR2F2",                                           #PTC/AVR
  "ZNF385D","RANBP3L","EDIL3",                               #EC-AVR
  "MX2","RSAD2","ISG15","IFIT1",                             #iaEC-AVR
  
  "SLCO2A1",                                                
  "VWF","RYR3","ADGRG6","CPE","TRABD2B",                     #EC-V
  "ADAMTSL1","CMKLR1",                                       #EC-V/EC-PCV
  "DOK6",                                                    #EC-PCV
  "NAV3","OSMR","SYCP2L",                                    #C-EC-PTC
  #"AFAP1L1","USP31","MYO1B","LAMA4","NETO2",                 #angEC-PTC
  "SLC6A6","FUT8","ATP13A3","AFF3",                          #EC-EA
  #"IFITM3","HLA-DRA","CAVIN2","CCL14","CA4",                 #dEC-PTC
  
  'ICAM1',"TNFAIP3",'CCL2','SELE',"CXCL2",                   #infEC-PTC
  "VCAM1",                                                   #inf/iaEC-PTC
  "CXCL10","GBP1","CXCL11","CTSS",                           #iaEC-PTC
  "MMRN1","CD36","TBX1","PROX1",                             #EC-LYM
  "TOP2A","MKI67","CENPF"                                    #cycling
  
  )
int.markers.2 <- c(
  "DCN","C7","PDGFRA",                                     #Pan FIB
  "SYT1","TNC","IGFBP2","RARRES2",                         #Pan Medullary FIB
  "CA8","HPSE2","GABRG3",                                  #IM/OM-FIB
  "MCTP2","SNAP25","BDKRB2",                               #IM-FIB
  "KIF26B","FREM1","KCNQ3","LOXHD1",                       #OM & C/M-FIB
  "SPP1","TIMP1",                                          #dM-FIB
  #"COL1A2","COL3A1","COL1A1",                              #dFIB & MYOF
  "ADAMTSL1","KCNIP1","ADARB2",                            #C/M-FIB
  "RYR2","ZNF536","SEMA3D","ACTG2",# "PAMR1",              #IM-pvMYOF
  
  "NEGR1","LAMA2","ABCA8","MEG3",                           #Pan cortical FIB
  "CCN1","CCN2","ELL2","SAMHD1","SLC2A3",                   #C-FIB (interstitial fib)
  "GRIN2A","EMID1",                                         #C-FIB-PATH
  "SELENOP","LUM","CXCL12",'GGT5',"ECRG4",                  #C-FIB-OSMRlo
  "OSMR","SOD2","UGCG","IL1R1","CCL2",                      #C-FIB-OSMRhi
  "RELB","CXCL10","CCL19",                                  #C-FIB-OSMRhi
  "SULF1","GLI2","NTM","INHBA","FAP","POSTN",               #C-MYOF
  #"SPARC","BGN","VIM",                                      #dFIB
  
  "FLRT2","COL12A1","FGF14",                                #Pan pvFIB
  "PDZRN4",'IGF1','ADAMTS3',"RSPO3","WNT5B",                #pvFIB-RSPO3+
  "C3","EBF2","SFRP2","CD34","PI16",                        #pvFIB-PI16+
  "ITGBL1","PLD5","CNTN5",                                  #pvFIB & pvMYOF
  "MGAT4C","EPHA3",                                         #pvFIB
  "ITGA8",                                                  #pvMYOF & VSMC
  "ADGRL3","MYH11","ACTA2",'KCNMA1',"PCDH7",                #pvMYOF
  "PRUNE2","MYOCD",#"SYNPO2",'MACROD2',                     #pvMYOF
  
  
  "PDGFRB","SLCO3A1",                                       #Pan VSM markers
  "SAMD12","ROBO1","PIEZO2",                                #MC & REN
  "DAAM2","GATA3","DCC","POSTN","IL1RL1",                   #MC
  "ROBO2","REN","KCNQ5","SLCO2A1","SPOCK3",                 #REN
  "MYH11","NTRK3",'RGS6','KCNMA1',"ADRA1A","MCAM",          #VSMC
  "NOTCH3",                                                 #VSMC &VSMC/P
  "RGS5",'PLCB4',"FRMD3",                                   #VSMC/P
  "ADGRB3","SLC38A11","C2CD2","DNAH11","SLC6A1",            #M-VSMC/P
  'PDE1C',"STEAP4","RXFP1",                                 #VSMC/P (cortical)
  #'FLNA', 'TAGLN',"ACTA2","VIM","ACTB",                     #dVSMC
  "CD36","PLIN1","LPL","ADIPOQ","FABP4"                     #Ad
  
)

int.markers.3 <- c(
  "PTPRC",                                             #Broad Immune
  "BANK1","MS4A1","CD37","CD79A",                      #B Cells
  "IGKC","XBP1","MZB1","JCHAIN",#"SDC1",               #PL Cells
  "CD96","CD247","BCL11B","THEMIS",                    #T
  "INPP4B","TRAC","CD3D",                              #T
  "IL7R",                                              #T/MAIT/ILC3
  "LEF1","CD4","SELL",                                 #Naïve Th
  "SLC4A10","KLRB1","CCR6",                            #MAIT
  "PCDH9","TOX2","KIT","RORC",                         #ILC3
  "IKZF2","RTKN2","IL2RA","CTLA4",#"FOXP3",            #T-REG
  "CD8A",                                              #CD8+
  "GZMK",                                              #CD8+ TEM/TRM
  "CCL5","SAMD3","GZMA","CCL4",                        #CD8+ & NK 
  "NKG7","KLRD1","GNLY","GZMB","CX3CR1",               #CD8+ TEM/TEMRA & NK
  "GZMH",                                              #CD8+ TEM/TEMRA
  "TXK","KLRF1","NCAM1",#"PDGFD",                      #NK
  
  
  "HBB","HBA2","HBA1",                                #Erythrocyte
  "CPA3","IL18R1","TPSB2","TPSAB1",                   #MAST
  "CD163","MSR1","CSF1R","CD14",                      #MAC
  "MRC1","F13A1","STAB1","CD163L1","LYVE1",           #resMAC-LYVE1+
  "HLA-DPA1","C1QA","C1QB",                           #resMAC-HLAIIhi
  "HIF1A","NAMPT","PLAUR","ITGAX","HBEGF","OSM",      #moMAC-HBEGF+
  "PSTPIP2","CXCL10","CXCL9","CCL2","CCL3","IL1B",    #moMAC-CXCL10+
  "GPNMB","SPP1","APOC1","PLA2G7","CD68","CAPG",      #moFAM
  "HMOX1","TREM2",                                    #moFAM
  "C3","KCNQ3","ADGRB3","VASH1","CX3CR1",             #moMAC-C3+
  "CLEC10A","FCER1A","CD1C",                          #cDC2
  "TCF7L2","COTL1","FCGR3A",                          #ncMON
  "FCN1",                                             #MON/ncMON
  "VCAN","LYZ","CD36",                                #MON
  "LAMP3","SLCO5A1","CCR7",#"EBI3","CCL19",           #mDC
  "WDFY4","CADM1","CLEC9A","BATF3",                   #cDC1
  "BCL11A","CLEC4C","IL3RA","PLD4",#"LILRA4",         #pDC
  "S100A9","FCGR3B","S100A8","IFITM2",                #N
  "TOP2A","MKI67",                                    #cycling
  "NRXN1","GRIK2","CDH19","NCAM2"                     #SC/NEU
  
)


#DotPlot(KB, features = unique(epi.markers), dot.scale = 8) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.1), dot.scale = 8) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.2), dot.scale = 8) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.3), dot.scale = 8) + RotatedAxis()

DimPlot(KB, reduction = "umap", group.by = "v2.subclass.l3", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "rpca_k200_infomap_nonEpi", label = TRUE,
        alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap", group.by = "assay", label = TRUE,
        alpha = 0.1) + NoLegend()

CL47_markers <- FindMarkers(object = xenium.obj, ident.1 = 47, 
                            max.cells.per.ident = 500, only.pos = TRUE)


cl.order <- c(
  9, #EC-GC
  28, #EC-AA
  27, #EC-DVR
  34, #M-EC-PTC
  35, #EC-PTC/EC-GC
  31, #EC-PTC
  43, #EC-PTC
  45, #EC-PTC
  40, #infEC-PTC
  36, #EC-EA
  11, #EC-V
  12, #EC-AVR
  18, #EC-LYM
  
  41, #IM-FIB
  44, #OM-FIB 
  42, #OM-FIB 
  14, #C-FIB
  5, #C-FIB
  24, #C-FIB
  19, #FIB
  37, #FIB
  17, #C-FIB-PATH
  26, #C-FIB-PATH
  21, #C-MYOF
  20, #pvFIB-RSPO3+
  38, #pvFIB-PI16+
  30, #pvMYOF
  16, #MC
  29, #REN
  33, #VSMC
  15, #VSMC/P
  
  4, #B
  10, #PL
  1, #T
  22, #NK/NKT
  23, #MAST
  6, #resMAC-LYVE1+
  7, #resMAC-LYVE1+
  25, #MAC
  13, #MON
  3, #cDC2
  2, #pDC
  8, #cDC1
  39, #cycIMM
  
  32, #ambiguous, possible SC/NEU
  
  46 #ambiguous, only 1 cell
  
)

Idents(KB) <- "rpca_k200_infomap_nonEpi"
Idents(KB) <- factor(Idents(KB), levels = cl.order)
DotPlot(KB, features = unique(int.markers.1), dot.scale = 8) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.2), dot.scale = 8) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.3), dot.scale = 8) + RotatedAxis()
Idents(xenium.obj) <- "rpca_k200_infomap_nonEpi"
Idents(xenium.obj) <- factor(Idents(xenium.obj), levels = cl.order)
DotPlot(xenium.obj, features = unique(int.markers.1), dot.scale = 8) + RotatedAxis()
DotPlot(xenium.obj, features = unique(int.markers.2), dot.scale = 8) + RotatedAxis()
DotPlot(xenium.obj, features = unique(int.markers.3), dot.scale = 8) + RotatedAxis()


#Spatial localization
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
kss <- AddMetaData(kss, metadata = KB$rpca_k200_infomap_nonEpi, col.name = "rpca_k200_infomap_nonEpi")
Idents(kss) <- "rpca_k200_infomap_nonEpi"

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(9 #EC-GC
                                                
                                                             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(28 #EC-AA
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(27 #EC-DVR

             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(34 #M-EC-PTC
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(43, #EC-GC
                                                45, #EC-GC
                                                31 #EC-PTC
                                                
                                                #35 #EC-PTC/EC-GC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(40 #infEC-PTC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(36 #EC-EA
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(11 #EC-V
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(12 #EC-AVR
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(18 #EC-LYM
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(41 #IM-FIB
                                                )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(44, #OM-FIB (dM-FIB)
                                                42 #OM-FIB
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(14, #C-FIB
                                                5, #C-FIB
                                                24 #C-FIB
                                                
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(37, #FIB
                                                19 #FIB
               
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(17, #C-FIB-PATH
                                                26 #C-FIB-PATH
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(21 #C-MYOF

             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(20 #pvFIB-RSPO3+
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(38 #pvFIB-PI16+
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(30 #pvMYOF
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(16 #MC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(29 #REN
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(33 #VSMC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(15 #VSMC/P
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(4 #B
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(10 #PL
             )))

ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(1 #T
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(22 #NK/NKT
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(23 #MAST
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(6, #resMAC-LYVE1+
                                                7 #resMAC-LYVE1+
                                                             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(25 #MAC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(13 #MON
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(3 #cDC2
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(2 #pDC
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(8 #cDC1
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(39 #cycIMM
             )))
ImageDimPlot(kss, cols = rep("red",4), size = 0.2, fov = "fov",
             cells = WhichCells(kss, idents = c(32 #SC/NEU
             )))


#save metadata
meta <- KB@meta.data
saveRDS(meta, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_nonEpi_Metadata.Rds")


