library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(dplyr)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB

###Spatial Object
kss <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")

#Update annotations - Proximal
meta <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_PT-TL_Metadata.Rds")
meta$rpca_k200_infomap_pt <- paste0("P_", meta$rpca_k200_infomap_pt)

cl.meta <- read.delim("Human_Kidney_Xenium_Cluster_Table_042025.txt")
emc <- c("v2.clusters","v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2",
         "v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$rpca_k200_infomap_pt, cl.meta$rpca_clusters)]
}
colnames(meta)
order <- c("library","nCount_Spatial","nFeature_Spatial",
           "control_probe_counts","genomic_control_counts","control_codeword_counts",
           "unassigned_codeword_counts","deprecated_codeword_counts","cell_area",                 
           "nucleus_area","nucleus_count","segmentation_method",
           "source","assay","experiment","patient","specimen",
           "condition_level3","condition_level2","condition_level1","condition",                                          
           "age_binned","sex","race","tissue_type",
           "v2.clusters",
           "v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2",
           "v2.state.l1","v2.class","v2.structure")

meta.pt <- meta[meta$assay %in% c("Xenium 5K"),order] 


#Update annotations - Distal
meta <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_TAL-CD_Metadata.Rds")
meta$rpca_k200_infomap_dt <- paste0("D_", meta$rpca_k200_infomap_dt)

cl.meta <- read.delim("Human_Kidney_Xenium_Cluster_Table_042025.txt")
emc <- c("v2.clusters","v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2",
         "v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$rpca_k200_infomap_dt, cl.meta$rpca_clusters)]
}
colnames(meta)
meta.dt <- meta[meta$assay %in% c("Xenium 5K"),order] 



#Update annotations - NonEpi
meta <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object_RNA-Merge-Projection_nonEpi_Metadata.Rds")
meta$rpca_k200_infomap_nonEpi <- paste0("O_", meta$rpca_k200_infomap_nonEpi)

cl.meta <- read.delim("Human_Kidney_Xenium_Cluster_Table_042025.txt")
emc <- c("v2.clusters","v2.subclass.full","v2.subclass.l3","v2.subclass.l2","v2.subclass.l1","v2.state.l2",
         "v2.state.l1","v2.class","v2.structure")
for(i in 1:length(emc)){
  meta[[emc[i]]] <- cl.meta[,emc[i]][match(meta$rpca_k200_infomap_nonEpi, cl.meta$rpca_clusters)]
}
colnames(meta)
meta.o <- meta[meta$assay %in% c("Xenium 5K"),order] 

meta <- rbind(meta.pt, meta.dt, meta.o)
dim(meta)

table(meta$v2.subclass.l3)

#update metadata 
kss <- subset(kss, cells = rownames(meta))
kss@meta.data <- meta[rownames(kss@meta.data),]

#remove ambiguous
kss <- subset(kss, v2.clusters %in% c(1:115))
table(kss$v2.subclass.l3)
kss
#An object of class Seurat 
#5001 features across 383771 samples within 1 assay 
#Active assay: Spatial (5001 features, 0 variable features)
#1 layer present: counts
#5 spatial fields of view present: X3782 X3723 X3916 X3778 fov

saveRDS(kss, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_Full_Object_04-2025.Rds")




###Final annotation/alignment Stats
meta <- xenium.obj@meta.data
stats <- cbind(
  data.frame(meta %>%
               group_by(library) %>%
               summarise_at(vars(nCount_Spatial,nFeature_Spatial), list(mean))),
  data.frame(meta %>%
               group_by(library) %>%
               tally())
)
rownames(stats) <- stats$library
stats <- stats[,-c(1,4)]
stats  

#Determine number of cell types
df <- t(table(xenium.obj$library, xenium.obj$v2.subclass.l3))
cols <- vector()
exps <- colnames(df)

for(i in exps){
  cols[i] <- length(df[which(df[,i] > 5),i]) 
}
stats$n_clusters <- cols[rownames(stats)]

write.table(stats, file="QC_Plots/Xenium_post-clustering_Stats_04012025.txt", sep = "\t", row.names=TRUE, col.names=TRUE)
