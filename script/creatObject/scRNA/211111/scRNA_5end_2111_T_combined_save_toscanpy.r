source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/211111/grouped_analysis/")
load("Cryo_T_grouped.rda")
load("../combined_analysis/Cryo_T_real.rda")
setwd("..")

Idents(T_sub_SCT) <- 'T_sub_MainCluster'
T_sub_SCT <- my_AddMeta(T_sub_SCT, T1$MainCluster_anno, Replace = T, allow_mismatch = T)
Idents(T_sub_SCT) <- 'T1_MainCluster_anno'
T_sub_SCT <- my_AddMeta(T_sub_SCT, T2$MainCluster_anno, Replace = T, allow_mismatch = T)
T_sub_SCT$MainCluster_anno_byGrouped <- T_sub_SCT$T2_MainCluster_anno

my_plotDim(T_sub_SCT, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, 
           group.by = 'MainCluster_anno_byGrouped', 
           title = '1/2wk T-cells Clusters by Grouped', 
           split.by = 'MainCluster_anno_byGrouped')

my_DotPlot_split(T_sub_SCT, features = T_markers_tmp2, group.by = 'MainCluster_anno_byGrouped') + RotatedAxis()


T_combined_meta <- T_sub_SCT@meta.data
save(T_combined_meta, file = "T_combined_meta.rda")

T_sub <- T_sub_SCT
library(SeuratDisk)
T_sub@commands <- list()
T_sub@assays$RNA@data <- T_sub@assays$RNA@counts
T_sub <- NormalizeData(T_sub, normalization.method = "LogNormalize", scale.factor = 10000)
T_sub@assays$RNA@scale.data <- new("Assay")@scale.data
DefaultAssay(T_sub) <-'RNA'
T_sub@assays$SCT <- NULL
T_sub@reductions <- list()
T_sub@graphs <- list()

# X <- new('Assay', counts = T_sub@assays$RNA@counts)
# T_sub@assays$RNA@data <- new("Assay")@data
# T_sub@assays$RNA@scale.data <- 
# T_sub@assays$RNA <- X

file.remove("T.h5Seurat")
file.remove("T.h5ad")
SaveH5Seurat(T_sub, filename = "T.h5Seurat")
Convert("T.h5Seurat", dest = "h5ad")
