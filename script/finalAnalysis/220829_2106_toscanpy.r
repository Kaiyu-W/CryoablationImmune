source("/mnt/c/Users/PC/Documents/GitHub/My_scRNA_pipeline/utilities.r")
setwd("~/Desktop/Cryo/")
TNK <- readRDS("TNK_3end.rds")

T_sub <- TNK
T_sub$IFN_Final <- unfactor(T_sub$IFN_Final)
T_sub$IFN_Final_Cluster <- unfactor(T_sub$IFN_Final_Cluster)
library(SeuratDisk)
T_sub@commands <- list()
T_sub@assays$RNA@data <- T_sub@assays$RNA@counts
T_sub <- NormalizeData(T_sub, normalization.method = "LogNormalize", scale.factor = 10000)
T_sub <- FindVariableFeatures(T_sub, selection.method = "vst", nfeatures = 2500)
T_sub@assays$RNA@scale.data <- new("Assay")@scale.data
DefaultAssay(T_sub) <-'RNA'
T_sub@assays$SCT <- NULL

# X <- new('Assay', counts = T_sub@assays$RNA@counts)
# T_sub@assays$RNA@data <- new("Assay")@data
# T_sub@assays$RNA@scale.data <- 
# T_sub@assays$RNA <- X

if (file.exists("T_3end.h5Seurat")) file.remove("T_3end.h5Seurat")
if (file.exists("T_3end.h5ad")) file.remove("T_3end.h5ad")
SaveH5Seurat(T_sub, filename = "T_3end.h5Seurat")
Convert("T_3end.h5Seurat", dest = "h5ad")
