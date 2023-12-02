source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load("Combined_analysis_TNK_real.rda")
load("T_sub_real_meta_220630.rda")
T_sub_real@meta.data <- T_meta_220630

## split Naive cells into CD4/CD8
FeaturePlot(T_sub_real, c("Cd4","Cd8a","Cd8b1"), cells = grep("T_Naive", T_sub_real$CD8_sub_MainCluster_paper))

my_Heatmap(T_sub_real, 'CD8_sub_MainCluster_paper', c("Cd4","Cd8a","Cd8b1"), 
           cells = grep("T_Naive", T_sub_real$CD8_sub_MainCluster_paper), use_pheatmap = F,
           cluster_rows = T, cluster_cols = T, border_color = NA)

# some cells has no expression of Cd4/Cd8
# unsupervised clustering:
Idents(T_sub_real) <- 'CD8_sub_MainCluster_paper'
T_sub_real <- FindSubCluster(T_sub_real, 'T_Naive', 'RNA_snn', 'Naive', 0.2)

my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Naive', 
           cells = grep("^T_Naive", T_sub_real$Naive))

my_Heatmap(T_sub_real, 'Naive', c("Cd4","Cd8a","Cd8b1"), 
           cells = grep("^T_Naive", T_sub_real$Naive), use_pheatmap = F,
           cluster_rows = T, cluster_cols = F, border_color = NA)

T_sub_real$Naive[T_sub_real$Naive == 'T_Naive_0'] = 'CD4_Naive'
T_sub_real$Naive[T_sub_real$Naive == 'T_Naive_1'] = 'CD8_Naive'
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'Naive')


## Find I-IFN in CD4
T_CD4 <- readRDS("D://Ji_Wangxue/T_CD4_sub.rds")
T_CD4_9_barcodes <- sapply(
    names(Idents(T_CD4))[Idents(T_CD4) == 9],
    FUN = function(x) {
        if (grepl("_1$",x)) {
            paste0("Old_Cryo_", sub("_1$","",x))
        } else {
            paste0("Old_NonCryo_", sub("_2$","",x))
        }
    }, USE.NAMES = F)

my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'Naive',
           cells.highlight = colnames(T_sub_real)[colnames(T_sub_real) %in% T_CD4_9_barcodes]
)
FeaturePlot(T_sub_real, feature = c("Cd4","Cd8a","Cd8b1"), 
            cells = colnames(T_sub_real)[colnames(T_sub_real) %in% T_CD4_9_barcodes])
FeaturePlot(T_sub_real, feature = c("Cd4","Cd8a","Cd8b1"), 
            cells = colnames(T_sub_real)[T_sub_real$CD8_sub_MainCluster_paper == 'I-IFN'])

# Actually I-IFN group involves that of both CD4 and CD8
my_Heatmap(T_sub_real, 'Naive', c("Cd4","Cd8a","Cd8b1"), 
           cells = colnames(T_sub_real)[T_sub_real$CD8_sub_MainCluster_paper == 'I-IFN'], use_pheatmap = F,
           cluster_rows = T, cluster_cols = T, border_color = NA, show_rownames = T)


# split I-IFN into CD4/CD8 and rename
Idents(T_sub_real) <- 'Naive'
IFN <- subset(T_sub_real, idents = 'I-IFN')

IFN <- my_process_seurat(IFN)
my_plotDim(IFN)
FeaturePlot(IFN, feature = c("Cd4","Cd8a","Cd8b1"))
IFN <- FindClusters(IFN, resolution = 2)
my_plotDim(IFN)
my_Heatmap(IFN, 'RNA_snn_res.2', c("Cd4","Cd8a","Cd8b1"), use_pheatmap = F,
           cluster_rows = T, cluster_cols = F, border_color = NA, show_rownames = T)
my_DotPlot_split(IFN, feature = c("Cd4","Cd8a","Cd8b1"))

anno_list <- list(
    "Tfh_I-IFN" = c(0,5,8),
    "CD8_Tem_I-IFN" = c(1,2,3,4,6,7)
)

IFN <- my_ClusterAnnote(
    Object = IFN, 
    anno_list = anno_list,
    meta_slot = 'RNA_snn_res.2',
    anno_slot = 'IFN_CD48',
    DefaultIdent_changeIf = TRUE,
    Duplicated_removeIf = TRUE
)
my_plotDim(IFN, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'IFN_CD48')

Idents(T_sub_real) <- "Naive"
T_sub_real <- my_AddMeta(T_sub_real, IFN$IFN_CD48, Replace = T)
T_sub_real$IFN_Final <- as.factor(T_sub_real$IFN_IFN_CD48)
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final')
Idents(T_sub_real) <- "IFN_Final"

my_DotPlot_split(T_sub_real, feature = T_markers_tmp)


## Generate Cluster-C
Idents(T_sub_real) <- "RNA_0.3"
tmp_obj <- subset(T_sub_real, ident = c("C3","C7","C4","C11"))
Idents(T_sub_real) <- "IFN_Final"
T_sub_real <- my_AddMeta(T_sub_real, tmp_obj$RNA_0.3, Replace = T)
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'tmp_obj_RNA_0.3')

# anno_list2 <- list(
#     'C2' = 'CD4_Naive',
#     'C3' = 'CD8_Naive',
#     'C1' = 'CD8_Tem',
#     'C15' = 'CD8_Tem_I-IFN',
#     'C5' = c('C4', 'CD8_Tex'),
#     'C10' = 'CD8_Tex_Proliferative',
#     'C4' = 'C3',
#     'C11' = 'NK_Proliferative',
#     'C12' = 'T_others',
#     'C7' = 'Tfh',
#     'C16' = 'Tfh_I-IFN',
#     'C6' = 'Th1',
#     'C14' = 'Th17',
#     'C9' = 'Treg',
#     'C13' = 'C11',
#     'C8' = 'C7'
# )
anno_list2 <- list(
    'C1' = 'CD4_Naive',
    'C2' = 'CD8_Naive',
    'C3' = 'CD8_Tem',
    'C4' = 'CD8_Tem_I-IFN',
    'C5' = c('C4', 'CD8_Tex'),
    'C6' = 'CD8_Tex_Proliferative',
    'C7' = 'C3',
    'C8' = 'NK_Proliferative',
    'C9' = 'T_others',
    'C10' = 'Tfh',
    'C11' = 'Tfh_I-IFN',
    'C12' = 'Th1',
    'C13' = 'Th17',
    'C14' = 'Treg',
    'C15' = 'C11',
    'C16' = 'C7'
)
T_sub_real <- my_ClusterAnnote(
    Object = T_sub_real, 
    anno_list = anno_list2,
    meta_slot = 'tmp_obj_RNA_0.3',
    anno_slot = 'IFN_Final_Cluster',
    DefaultIdent_changeIf = TRUE,
    Duplicated_removeIf = TRUE
)
T_sub_real$IFN_Final_Cluster <- factor(unfactor(T_sub_real$IFN_Final_Cluster), levels = paste0("C", 1:16))
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final_Cluster')


# new meta saved in :
#     T_sub_real$IFN_Final
#     T_sub_real$IFN_Final_Cluster
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final') + theme(legend.position = 'none')
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final_Cluster')  + theme(legend.position = 'none')
my_DotPlot_split(T_sub_real, feature = T_markers_tmp, group.by = 'IFN_Final') + RotatedAxis()
my_DotPlot_split(T_sub_real, feature = T_markers_tmp, group.by = 'IFN_Final_Cluster') + RotatedAxis()

## Stats
my_CountCluster(T_sub_real, group1 = 'IFN_Final', group2 = 'orig.ident2')
my_CountCluster(T_sub_real, group1 = 'IFN_Final_Cluster', group2 = 'orig.ident2')


## Save
T_meta_220713 <- T_sub_real@meta.data
save(T_meta_220713, file = "T_sub_real_meta_220713.rda")

## To scanpy
T_sub <- T_sub_real
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
T_sub@reductions <- list()
T_sub@graphs <- list()

# X <- new('Assay', counts = T_sub@assays$RNA@counts)
# T_sub@assays$RNA@data <- new("Assay")@data
# T_sub@assays$RNA@scale.data <- 
# T_sub@assays$RNA <- X

if (file.exists("T_sub_real.h5Seurat")) file.remove("T_sub_real.h5Seurat")
if (file.exists("T_sub_real.h5ad")) file.remove("T_sub_real.h5ad")
SaveH5Seurat(T_sub, filename = "T_sub_real.h5Seurat")
Convert("T_sub_real.h5Seurat", dest = "h5ad")
