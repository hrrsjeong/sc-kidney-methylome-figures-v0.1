library(Seurat)
library(SeuratObject)
library(ggplot2)
library(BPCells)
library(future)
plan("multisession", workers = 10)
library(anndataR)

setwd("~/Projects/Human_Kidney/Atlas_V2/xenium/")
options(future.globals.maxSize = 70 * 1024^3)  # Set limit to 70 GiB


ts <- c("3782_Xen14","3723_Xen15","3916_Xen15","3778_Xen14")

ts.names <- c("3782","3723","3916","3778")
names(ts.names) <- ts

x <- "3782_Xen14"

xo.list <- lapply(ts, function(x) {
  print(paste("Running for Tissue Sample:", x))
  
  adata <- read_h5ad(paste0("~/Projects/Human_Kidney/Atlas_V2/xenium/",x,"_obj.h5ad"), to = "HDF5AnnData")
  obj <- to_Seurat(adata, layers_mapping = list(counts = "counts"),
                   assay_name = "Spatial")
  
  coords <- adata$obsm$spatial
  rownames(coords) <- adata$obs_names
  colnames(coords) <- c("x","y")
  
  expr <- obj[["Spatial"]]$counts
  expr <- convert_matrix_type(expr, type = "uint32_t")
  write_matrix_dir(
    mat = expr, 
    dir = paste0("~/hsKidAt/blake_LTS/Atlas_V2/xenium/bpcells_counts_",x),
    overwrite = TRUE)
  expr <- open_matrix_dir(dir = paste0("~/hsKidAt/blake_LTS/Atlas_V2/xenium/bpcells_counts_",x))                 
  
  obj[["Spatial"]]$counts <- expr
  coords <- coords[rownames(coords) %in% colnames(obj),]
  dim(coords)
  obj <- subset(obj, cells = rownames(coords))
  
  centroids = CreateCentroids(coords)
  fov.name <- ts.names[x]
  fov <- CreateFOV(
    centroids, 
    type = "centroids",
    assay = "Spatial",
    key = Key(fov.name, quiet = TRUE)
  )
  fov <- fov[Cells(obj)]
  obj[[fov.name]] <- fov
  return(obj) 
})

ts.names <- c("3782","3723","3916","3778")
names(xo.list) <- ts.names
xenium.obj <- merge(x = xo.list[[1]],
             y = c(xo.list[[2]],xo.list[[3]],xo.list[[4]]))  
xenium.obj
xenium.obj[["Spatial"]] <- JoinLayers(xenium.obj[["Spatial"]])


###Save object
counts <- xenium.obj[["Spatial"]]$counts
write_matrix_dir(
  mat = counts,
  dir = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Spatial_counts_003282025",
  overwrite = TRUE
)
counts <- open_matrix_dir(dir = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Spatial_counts_003282025")
xenium.obj[["Spatial"]]$counts <- counts
xenium.obj
#An object of class Seurat 
#5001 features across 439033 samples within 1 assay 
#Active assay: Spatial (5001 features, 0 variable features)
#1 layer present: counts
#4 spatial fields of view present: X3782 X3723 X3916 X3778


###make combined fov
coords1 <- GetTissueCoordinates(xenium.obj, image = "X3782")
coords2 <- GetTissueCoordinates(xenium.obj, image = "X3723")
coords3 <- GetTissueCoordinates(xenium.obj, image = "X3916")
coords4 <- GetTissueCoordinates(xenium.obj, image = "X3778")

min.y1 <- min(coords1$y)
max.y1 <- max(coords1$y)
min.x1 <- min(coords1$x)
max.x1 <- max(coords1$x)
min.y2 <- min(coords2$y)
max.y2 <- max(coords2$y)
min.x2 <- min(coords2$x)
max.x2 <- max(coords2$x)
min.y3 <- min(coords3$y)
max.y3 <- max(coords3$y)
min.x3 <- min(coords3$x)
max.x3 <- max(coords3$x)
min.y4 <- min(coords4$y)
max.y4 <- max(coords4$y)
min.x4 <- min(coords4$x)
max.x4 <- max(coords4$x)

coords2$y <- coords2$y + 8000
coords3$y <- coords3$y + 8000
coords3$x <- coords3$x + 7000
coords4$x <- coords4$x + 7000

coords <- rbind(coords1, coords2, coords3, coords4)
rownames(coords) <- coords$cell
coords <- coords[rownames(coords) %in% colnames(xenium.obj),]
dim(coords)

centroids = CreateCentroids(coords[,c("x", "y")])
fov <- CreateFOV(
  centroids, 
  type = "centroids",
  assay = "Spatial",
  key = Key("fov", quiet = TRUE)
)

fov <- fov[Cells(xenium.obj)]
xenium.obj[["fov"]] <- fov

epi.cols <- setNames(c("#F1085C","#783FC1","#0000FF","#B1CC71","#005300","#FFFFFF",
                       "#886C00","#FFD300","#14F9FF"),
                     c("POD","PEC","PT-S1","PT-S2","PT-S3","DTL1",     
                       "DTL2","DTL3","ATL"))

Idents(xenium.obj) <- "v2.subclass.l2"
ImageDimPlot(xenium.obj, fov = "fov",
             cells = WhichCells(xenium.obj, 
                                idents = c("POD","PEC","PT-S1","PT-S2","PT-S3","DTL1",     
                                           "DTL2","ATL")),
             cols = epi.cols[levels(Idents(xenium.obj))], size = 0.2, border.size = NA) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


saveRDS(xenium.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
#xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")


#update metadata
meta <- xenium.obj@meta.data
coords1 <- GetTissueCoordinates(xenium.obj, image = "X3782")
coords2 <- GetTissueCoordinates(xenium.obj, image = "X3723")
coords3 <- GetTissueCoordinates(xenium.obj, image = "X3916")
coords4 <- GetTissueCoordinates(xenium.obj, image = "X3778")

meta$library <- "3782_Xen14"
meta$library[match(coords2$cell, rownames(meta))] <- "3723_Xen15"
meta$library[match(coords3$cell, rownames(meta))] <- "3916_Xen15"
meta$library[match(coords4$cell, rownames(meta))] <- "3778_Xen14"
table(meta$library)

exp.meta <- read.delim("Human_Kidney_Xenium_Exp_Meta_03282025.txt")
emc <- c("source","assay","experiment","patient","specimen",
         "condition_level3","condition_level2","condition_level1","condition",
         "age_binned","sex","race","tissue_type")

for(i in 1:length(emc)){
  meta[[emc[i]]] <- exp.meta[,emc[i]][match(meta$library, exp.meta$library)]
}
xenium.obj@meta.data <- meta
to.use <- c("library", "nCount_Spatial", "nFeature_Spatial", "cell_id",
            "control_probe_counts", "genomic_control_counts", "control_codeword_counts",
            "unassigned_codeword_counts","deprecated_codeword_counts", "cell_area",
            "nucleus_area", "nucleus_count","segmentation_method",
            "source", "assay", "experiment", "patient", "specimen",
            "condition_level3", "condition_level2", "condition_level1", "condition",
            "age_binned",  "sex", "race","tissue_type","v2.subclass.l1", "v2.subclass.l2"
            )
xenium.obj@meta.data <- xenium.obj@meta.data[,to.use]

Idents(xenium.obj) <- "library"

saveRDS(xenium.obj, file = "~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
#xenium.obj <- readRDS("~/hsKidAt/blake_LTS/Atlas_V2/xenium/Human_Kidney_Xenium_03-2025_Object.Rds")
