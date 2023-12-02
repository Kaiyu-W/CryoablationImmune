
library(Seurat)
library(SeuratDisk)

setwd("E://Cryo-TCR/server/210924/")
source("../auto/utilities.r")
load("Combined_analysis.rda")
load("Combined_analysis_Myeloid_real.rda")
load("Combined_analysis_TNK_real.rda")
load("Combined_analysis_T_CD4_8.rda")


# cell type defined
Cryo_merge_temp <- my_AddMeta(Cryo_merge, T_sub_real$T_sub_MainCluster2, Replace = T)
Idents(Cryo_merge_temp) <- 'T_sub_real_T_sub_MainCluster2'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, Myeloid_sub_real$MainCluster2_new, Replace = T)
Cryo_merge_temp$Myeloid_sub_real_MainCluster2_new[Cryo_merge_temp$Myeloid_sub_real_MainCluster2_new == 'T'] <- 'epithelial'
Idents(Cryo_merge_temp) <- 'Myeloid_sub_real_MainCluster2_new'
Cryo_merge_temp$EndMainCluster <- Cryo_merge_temp$Myeloid_sub_real_MainCluster2_new

my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 1.5, label.size = 6, group.by = 'EndMainCluster', title = 'Cryo_merge Cluster')
#
# Cryo_merge_temp@assays$SCT <- NULL
# Cryo_merge_temp@assays$RNA@data <- matrix()

file.remove("Cryo_merge.h5Seurat")
file.remove("Cryo_merge.h5ad")
SaveH5Seurat(Cryo_merge_temp, filename = "Cryo_merge.h5Seurat")
Convert("Cryo_merge.h5Seurat", dest = "h5ad")

