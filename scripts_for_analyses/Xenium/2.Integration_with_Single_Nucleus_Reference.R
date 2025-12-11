library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(pagoda2)
library(future)
plan("multisession", workers = 10)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

###Integrate with sc Reference
xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")

##Load Reference
ref <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/scratch/Kidney_AtlasV2_Seurat_11012024.rds")
Idents(ref) <- "v2.subclass.l3"
table(Idents(ref))
ref <- subset(ref, downsample = 5000)

###Merge objects
xenium.obj[["RNA"]] <- xenium.obj[["Spatial"]]
DefaultAssay(xenium.obj) <- "RNA"
xenium.obj[["fov"]] <- NULL
xenium.obj[["Spatial"]] <- NULL
xenium.obj[["X3782"]] <- NULL
xenium.obj[["X3723"]] <- NULL
xenium.obj[["X3916"]] <- NULL
xenium.obj[["X3778"]] <- NULL
xenium.obj[["RNA"]]$counts

ref[["RNA"]]$data <- NULL
ref[["pca"]] <- NULL
ref[["umap"]] <- NULL
ref[["RNA"]]$counts


KB <- merge(x = ref, y = c(xenium.obj))
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
  group.by = c("assay", "v2.subclass.l1"),
  combine = TRUE, label.size = 2
)

DimPlot(
  KB,
  reduction = "umap",
  group.by = c("v2.subclass.l2"),
  label = TRUE, label.size = 2
) + NoLegend()

saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge.Rds")


###Clustering on sketch data set
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge.Rds")

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
p2$makeKnnGraph(k = 100, type = 'PCA', center = T, weight.type = 'none', n.cores = 6, distance = 'cosine')

# Identify clusters using the infomap.community method
p2$getKnnClusters(method = igraph::infomap.community, type = 'PCA', name = 'infomap')

#Add pagoda2 clusters 
k100infomap <- p2$clusters$PCA$infomap
KB <- AddMetaData(KB, metadata = k100infomap, col.name = "rpca_k100_infomap")
DimPlot(KB, group.by = "rpca_k100_infomap", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(KB, group.by = "v2.subclass.l3", reduction = "umap", label = TRUE) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_B.Rds")



#Identify max predicted subclass per cluster
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge_B.Rds")
Idents(KB) <- "rpca_k100_infomap"
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
max.idents[1:50,]
max.idents[51:100,]
max.idents[101:150,]



###Cluster Annotations from sc ref overlap and marker genes
#1 - MD
#2 - frPT-S1/S2
#3 - C-TAL
#4 - T
#5 - resMAC-LYVE1+
#6 - IC-B
#7 - C-FIB
#8 - cDC2
#9 - CNT-PC
#10 - B
#11 - aCNT
#12 - DCT1
#13 - EC-LYM
#14 - PL
#15 - C-TAL
#16 - aPT-S3
#17 - ncMON
#18 - CNT
#19 - PT
#20 - PC
#21 - POD
#22 - aTAL
#23 - cDC1
#24 - frDCT
#25 - cycPT
#26 - cycT
#27 - C-FIB
#28 - PT-S1
#29 - ambiguous
#30 - C-EC-PTC
#31 - pvFIB-RSPO3+
#32 - EC-AVR
#33 - EC-GC
#34 - frTAL
#35 - PEC
#36 - EC-PTC
#37 - PT-S1/S2
#38 - C-FIB
#39 - OMCD-IC-A
#40 - C-TAL
#41 - VSMC
#42 - MC
#43 - EC-DVR
#44 - EC-PTC
#45 - PT-S3
#46 - aTAL1
#47 - pvMYOF
#48 - MAST
#49 - OM-FIB
#50 - VSMC/P
#51 - moFAM
#52 - IMCD
#53 - IM-FIB
#54 - ATL
#55 - PT-S1
#56 - OM-FIB
#57 - M-TAL
#58 - DTL2
#59 - aCNT
#60 - PapE
#61 - PapE
#62 - cycMAC
#63 - ambiguous


##Check marker genes
Idents(KB) <- "rpca_k100_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

epi.markers <- c(
    "NPHS1","NPHS2",#"ST6GALNAC3","PODXL",                         #POD
    #"PTGDS","BST2","TPPP3","CHI3L1",                              #dPOD
    "ALDH1A2","CFH",#"FAM155A","CLDN1",                            #PEC
    
    "LRP2","CUBN",                                                #PT
    "SLC5A12","SLC22A6",#"HNF4A",                                  #S1/S2
    "RALYL","PCDH15","PRODH2",#"SLC22A8",                          #S1
    "SLC34A1","ANKS1B","SLC5A10","SLC5A11",                       #S2                                  
    "SLC5A8","GPM6A","SLC22A24","SLC7A13",                        #PT-S3
    
    #"IGFBP7","SPP1","ITGB8","CDH6","TMEM178B","ALPK2","HAVCR1",
    "ITGB3",
    "CST3","CLU","VIM","PIGR",#"APOE",                            #aPT2
    "IL32","SOX4","VCAM1",#"MMP7","SOX9","CCL2",                   #aPT2
    "DCC",                                                        #aPT1
    "GDA","GLIS1",                                                #aPT-S1/S2
    "DLGAP1","PROM1",                                             #frPTS1/S2
    #"APBB1IP","ROBO2","COL23A1","MEG3",                           #frPTS1/S2
    #"LSAMP","KCNIP1","NRXN3","WNT2B",                             #frPTS1/S2
    "KCTD16","SPON1",                                             #aPT-S3
    "NRG3","FAM189A1","DTNA","KITLG","GRM8",                      #frPT-S3
    "TOP2A","MKI67",                                              #cycling
    #"EGR1","FOS",                                                 #dPT
    
    "PAX8-AS1","SLC44A5","CRYAB","TACSTD2",                       #TL
    #"ABCA13",                                                     #DTL
    "AQP1", "UNC5D","LRRC4C","DENND2A","SLCO1A2",                 #DTL2
    #"IRX3", "SERPINA1",                                           #aDTL2
    "SIM2",                                                       #DTL1/3/ATL
    "ADGRL3","JAG1","SMAD9","ID1",                                #DTL1
    "SLC14A2","FAM169A","SMOC2",                                  #DTL3
    #"ABCA4","BCL6","AKR1B1","SH3GL3",                             #DTL3/ATL
    "CLCNKA","PROX1","CACNA1C","COLEC10",                         #ATL
    #"SOD3","HSPA2","CSRP2",                                       #dATL
    #"CST3","APOE","GATM","ALDOB","CLU","DEFB1",                   #Injury
    #"SPP1","IGFBP7","CLU",
    
    
    "CASR","SLC12A1","UMOD",#"EGF","ESRRB",                  #TAL
    "ANK2","CLCNKA","KCTD16","BMP6","RIPOR3","CLDN14",      #M-TAL
    "PHACTR1","SLCO3A1",#"CXCL12","CNTN1","CABP1",           #TAL-A
    "KCNMB2",#"RGS6",                                        #C/M-TAL-A
    "ENOX1","CALCR",#"RBM20","PDE3A",                        #C-TAL-A
    #"DACH1","LRMDA",                                        #TAL-B
    
    "TENM4","FGF14",#"PXDNL","GRM8",                         #C/M-TAL-B
    "KCNJ10","TMEM52B","CLDN16","TMEM207","JAG1",           #TAL-B
    #"SCN7A","COL8A1","LINGO2",                              #C-TAL-B                          
    "BBOX1","NOS1","ROBO2","TMPRSS4",                       #MD
    
    "CREB5","COL4A1","ITGA3","HIVEP3","CYTOR",              #aTAL1
    #"CD44","DIAPH3","SYT16","HAVCR1",                       #aTAL1                                      
    
    "ITGB6","NRP1","TFPI",                                  #aTAL
    
    "HIF1A","ADAMTS1","DNAH5",                              #aTAL2
    
    "ITGB8","PROM1","ARHGAP26",#"RNF144B","TMPRSS4","RHEX",  #frTAL
    #"EZH2","CENPP","MKI67",#"TOP2A",                        #cycling
    
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
    #"FGF12","PTCSC3","CNTNAP3B",                             #CCD/OMCD PC
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
  #"MX2","RSAD2","ISG15","IFIT1",                             #iaEC-AVR
  
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
  "TOP2A","MKI67","CENPF",                                    #cycling
  
  
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
  #"RELB","CXCL10","CCL19",                                  #C-FIB-OSMRhi
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
int.markers.2 <- c(
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


DotPlot(KB, features = unique(epi.markers), dot.scale = 8, idents = 1:30) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.1), dot.scale = 8, idents = 1:30) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.2), dot.scale = 8, idents = 1:30) + RotatedAxis()

DotPlot(KB, features = unique(epi.markers), dot.scale = 8, idents = 31:63) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.1), dot.scale = 8, idents = 31:63) + RotatedAxis()
DotPlot(KB, features = unique(int.markers.2), dot.scale = 8, idents = 31:63) + RotatedAxis()


CL29.markers <- FindMarkers(KB, ident.1 = 29, min.pct = 0.25, only.pos = TRUE)
CL63.markers <- FindMarkers(KB, ident.1 = 63, min.pct = 0.25, only.pos = TRUE)





##Project clusters to full data set
meta <- KB@meta.data
KB <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge.Rds")
KB <- AddMetaData(KB, metadata = meta[,"rpca_k100_infomap"], col.name = "rpca_k100_infomap")
KB <- ProjectIntegration(object = KB, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")

KB <- ProjectData(object = KB, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
                  full.reduction = "integrated.rpca.full", dims = 1:50, refdata = list(rpca_k100_infomap = "rpca_k100_infomap"))

KB <- RunUMAP(KB, reduction = "integrated.rpca.full", dims = 1:50, reduction.name = "umap.full",
              reduction.key = "UMAP_full_")

Idents(KB) <- "rpca_k100_infomap"
Idents(KB) <- factor(Idents(KB), levels = 1:length(levels(Idents(KB))))

DimPlot(KB, reduction = "umap.full", group.by = "v2.subclass.l3", alpha = 0.1) + NoLegend()
DimPlot(KB, reduction = "umap.full", group.by = "rpca_k100_infomap", label = TRUE,
        alpha = 0.1) + NoLegend()


saveRDS(KB, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection.Rds")

