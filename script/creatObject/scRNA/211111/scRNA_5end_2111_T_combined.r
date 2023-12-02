source("E://Cryo-TCR/server/auto/utilities.r")
library(reshape)
library(RColorBrewer)
setwd("E://Cryo-TCR/server/211111/")

Cryo_data_1 <- Read10X_h5("Cryo_1wk/filtered_feature_bc_matrix.h5")
NonCryo_data_1 <- Read10X_h5("NonCryo_1wk/filtered_feature_bc_matrix.h5")
Cryo_data_2 <- Read10X_h5("Cryo_2wk/filtered_feature_bc_matrix.h5")
NonCryo_data_2 <- Read10X_h5("NonCryo_2wk/filtered_feature_bc_matrix.h5")

colnames(Cryo_data_1) <- paste0("Cryo1wk_", colnames(Cryo_data_1))
colnames(NonCryo_data_1) <- paste0("NonCryo1wk_", colnames(NonCryo_data_1))
colnames(Cryo_data_2) <- paste0("Cryo2wk_", colnames(Cryo_data_2))
colnames(NonCryo_data_2) <- paste0("NonCryo2wk_", colnames(NonCryo_data_2))

Cryo_1 <- CreateSeuratObject(counts = Cryo_data_1,
                             project = "Cryo1wk",
                             min.cells = 3,
                             min.features = 200)
NonCryo_1 <- CreateSeuratObject(counts = NonCryo_data_1,
                                project = "NonCryo1wk",
                                min.cells = 3,
                                min.features = 200)
Cryo_2 <- CreateSeuratObject(counts = Cryo_data_2,
                             project = "Cryo2wk",
                             min.cells = 3,
                             min.features = 200)
NonCryo_2 <- CreateSeuratObject(counts = NonCryo_data_2,
                                project = "NonCryo2wk",
                                min.cells = 3,
                                min.features = 200)

Cryo_1@meta.data$orig.ident <- Cryo_1@project.name
NonCryo_1@meta.data$orig.ident <- NonCryo_1@project.name
Cryo_2@meta.data$orig.ident <- Cryo_2@project.name
NonCryo_2@meta.data$orig.ident <- NonCryo_2@project.name

Cryo_combined <- merge(Cryo_1, c(NonCryo_1, Cryo_2, NonCryo_2))
Cryo_combined

Cryo_combined[["percent.mt"]] <- PercentageFeatureSet(Cryo_combined, pattern = "^mt-")
Idents(Cryo_combined) <- 'orig.ident'
Cryo_combined$orig.ident_group <- Cryo_combined$orig.ident
Cryo_combined$orig.ident_group <- sub("^Cryo.*$", "Cryo", Cryo_combined$orig.ident_group)
Cryo_combined$orig.ident_group <- sub("^NonCryo.*$", "NonCryo", Cryo_combined$orig.ident_group)
Cryo_combined$orig.ident_time <- Cryo_combined$orig.ident
Cryo_combined$orig.ident_time <- sub("^.*1wk$", "1wk", Cryo_combined$orig.ident_time)
Cryo_combined$orig.ident_time <- sub("^.*2wk$", "2wk", Cryo_combined$orig.ident_time)

# quality filter
Cryo_merge <- subset(Cryo_combined, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA > 1000 & nCount_RNA < 20000 & percent.mt < 5)
Cryo_merge

# process
DefaultAssay(Cryo_merge) <- 'RNA'
Cryo_merge <- NormalizeData(Cryo_merge, normalization.method = "LogNormalize", scale.factor = 10000)
Cryo_merge <- ScaleData(Cryo_merge, assay = "RNA", features = rownames(Cryo_merge))

Cryo_merge <- SCTransform(Cryo_merge)
DefaultAssay(Cryo_merge) <- 'SCT'
Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = 2500)
Cryo_merge <- RunPCA(Cryo_merge)

Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10)
Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10)

Cryo_merge  <- FindClusters(Cryo_merge, resolution = 1:9/100)
Cryo_merge  <- FindClusters(Cryo_merge, resolution = 1:10/10)
clustree(Cryo_merge)

Cryo_merge  <- FindClusters(Cryo_merge, resolution = 0.6)
Idents(Cryo_merge) <- 'SCT_snn_res.0.3'
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge R=0.3', group.by = 'SCT_snn_res.0.3')
my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge R=0.3', group.by = 'SCT_snn_res.0.3')

my_DotPlot_split(Cryo_merge, features = my_MarkersList[1:7]) + RotatedAxis()
my_violin(Cryo_merge, features = unlist(my_MarkersList[1:7]), pt.size = 0, mode = 'mtx')

my_DotPlot_split(Cryo_merge, features = celltype_marker) + RotatedAxis()

# cluster 8/14/15 includes some non-T cells; 9/12 chaos
DoHeatmap(Cryo_merge, features = unlist(my_MarkersList[1:7]), cells = which(Cryo_merge$SCT_snn_res.0.3 %in% c(8,9,12,14,15)), slot = 'scale.data', assay = 'RNA')
DoHeatmap(Cryo_merge, features = c('Ptprc','Pecam1','Cd3e','Cd4','Cd8a','Cd8b1','Itgam','Cd14'), cells = which(Cryo_merge$SCT_snn_res.0.3 %in% c(8,9,12,14,15)), slot = 'scale.data', assay = 'RNA')
DoHeatmap(Cryo_merge, features = c('Ptprc','Pecam1','Cd3e','Cd4','Cd8a','Cd8b1','Itgam','Cd14'), slot = 'scale.data', assay = 'RNA')

# cell filter
T_sub <- subset(Cryo_merge, ident = setdiff(0:15, c(8,14,15)))

# ReCluster
DefaultAssay(T_sub) <- 'RNA'
T_sub <- NormalizeData(T_sub, normalization.method = "LogNormalize", scale.factor = 10000)
T_sub <- ScaleData(T_sub, assay = "RNA", features = rownames(T_sub))

T_sub <- SCTransform(T_sub)
DefaultAssay(T_sub) <- 'SCT'
T_sub <- FindVariableFeatures(T_sub, selection.method = "vst", nfeatures = 2500)
T_sub <- RunPCA(T_sub)

T_sub <- FindNeighbors(T_sub, dims = 1:12)
T_sub <- RunUMAP(T_sub, dims = 1:10)
T_sub <- RunTSNE(T_sub, dims = 1:10)

T_sub  <- FindClusters(T_sub, resolution = 1:9/100)
T_sub  <- FindClusters(T_sub, resolution = 1:10/10)
clustree(T_sub)

my_plotDim(T_sub, reduction = "umap", label = F, pt.size = 0.1, label.size = 5, group.by = 'orig.ident')

Idents(T_sub) <- 'SCT_snn_res.0.5'
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T_sub R=0.5', group.by = 'SCT_snn_res.0.5')
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T_sub R=0.5', group.by = 'SCT_snn_res.0.5')

my_DotPlot_split(T_sub, features = my_MarkersList[1:7]) + RotatedAxis()

my_DotPlot_split(T_sub, features = T_markers_tmp) + RotatedAxis()
my_DotPlot_split(T_sub, features = T_markers_tmp[1:2]) + RotatedAxis()
my_violin(T_sub, features = unlist(T_markers_tmp), pt.size = 0, mode = 'mtx')
my_violin(T_sub, features = unlist(T_markers_tmp[1:2]), pt.size = 0, mode = 'mtx')
my_violin(T_sub, features = unlist(NK_marker), pt.size = 0, mode = 'mtx')

my_violin(T_sub, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = c(3,8,12,13))
my_violin(T_sub, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = c(0,1,2,4,5,6,7,9,10,11))
my_DotPlot_split(T_sub, features = T_marker$CD4, idents = c(3,8,12,13)) + RotatedAxis()
my_DotPlot_split(T_sub, features = T_marker$CD8, idents = c(0,1,2,4,5,6,7,9,10,11)) + RotatedAxis()

# Cluster 3/12: CD4/CD8
mtx_tmp <- T_sub@assays$RNA@data[unlist(T_markers_tmp[1:2]), which(T_sub$SCT_snn_res.0.5 %in% c(3,12))]
pheatmap::pheatmap(mtx_tmp, show_rownames = T, show_colnames = F)

gene_T <- unlist(T_markers_tmp[1:2])
data_T <- mtx_tmp
data_stat_T <- unlist(apply(data_T, 2, function(x) {
    if(x[1]>0) 
        0
    else if (x[2]>0 || x[3]>0)
        1
    else
        0
}))
names_bc <- names(data_stat_T)[data_stat_T == 1]
T_sub$sub3_12 <- unfactor(T_sub$SCT_snn_res.0.5)
T_sub@meta.data$sub3_12[which(rownames(T_sub@meta.data) %in% names_bc)] <- 14
rm(gene_T, data_T, data_stat_T, names_bc)
Idents(T_sub) <- 'sub3_12'
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T_sub R=0.5', group.by = 'sub3_12', cells = grep(14,T_sub$sub3_12))

# Cluster 1: CD4/CD8
mtx_tmp <- T_sub@assays$RNA@data[unlist(T_markers_tmp[1:2]), which(T_sub$SCT_snn_res.0.5 %in% c(1))]
pheatmap::pheatmap(mtx_tmp, show_rownames = T, show_colnames = F)

gene_T <- unlist(T_markers_tmp[1:2])
data_T <- mtx_tmp
data_stat_T <- unlist(apply(data_T, 2, function(x) {
    if(x[2]>0 || x[3]>0) 
        0
    else if (x[1]>0)
        1
    else
        0
}))
names_bc <- names(data_stat_T)[data_stat_T == 1]
T_sub$sub1 <- T_sub$sub3_12
T_sub@meta.data$sub1[which(rownames(T_sub@meta.data) %in% names_bc)] <- 15
rm(gene_T, data_T, data_stat_T, names_bc)
Idents(T_sub) <- 'sub1'
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T_sub R=0.5', group.by = 'sub1', cells = grep(15,T_sub$sub1))

# Cluster 3 : Th1/fh
Idents(T_sub) <- 'sub1'
T_sub <- FindSubCluster(T_sub, cluster = '3', graph.name = 'SCT_snn', resolution = 0.1, subcluster.name = 'sub3')
Idents(T_sub) <- 'sub3'
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = 'sub3')
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'sub3')
my_DotPlot_split(T_sub, features = T_marker$CD4, idents = c(unique(T_sub$sub3)[grep("^3_", unique(T_sub$sub3))],8,12,13,15)) + RotatedAxis()
my_violin(T_sub, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = c(unique(T_sub$sub3)[grep("^3_", unique(T_sub$sub3))],8,12,13,15))

# Cluster 12 : Th1/2/reg
DoHeatmap(T_sub, features = unlist(T_marker$CD4), cells = grep("^12_", T_sub$sub12)) # Tregs also express the markers of Th1/2

# markers
my_violin(T_sub, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = c("3_0","3_1",8,12,13,15))
my_violin(T_sub, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = c(0,1,2,4,5,6,7,9,10,11,14))
my_DotPlot_split(T_sub, features = T_marker$CD4, idents = c("3_0","3_1",8,12,13,15)) + RotatedAxis()
my_DotPlot_split(T_sub, features = T_marker$CD8, idents = c(0,1,2,4,5,6,7,9,10,11,14)) + RotatedAxis()

# cell type identify
T_sub$T_sub_MainCluster <- as.factor(T_sub$sub3)
T_sub_MainCluster_df <- rbind(
    data.frame(x='T_Naive',y=c(13)),
    data.frame(x='Treg',y=c(12)),
    data.frame(x='Th17',y=c(8)),
    data.frame(x='Th1',y=c("3_0",15)),
    data.frame(x='Tfh',y=c("3_1")),
    data.frame(x='CD8_Tem',y=c(4,5,9)),
    data.frame(x='CD8_Tea',y=c(7)),
    data.frame(x='CD8_Tex',y=c(6,10,11)),
    data.frame(x='T_others',y=c(0,1,2,14))
)
T_sub_MainCluster_df <- T_sub_MainCluster_df[order(T_sub_MainCluster_df$y, decreasing = F),]
if(all(T_sub_MainCluster_df$y == levels(T_sub$T_sub_MainCluster))) {
    levels(T_sub$T_sub_MainCluster) <- T_sub_MainCluster_df$x
}
Idents(T_sub) <- 'T_sub_MainCluster'
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster")
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster")

# Cluster T_others(0/1/2/14) : Tem/T_Naive/Tea
T_others <- subset(T_sub, ident = 'T_others')
# Re-ReCluster
DefaultAssay(T_others) <- 'RNA'
# T_others <- NormalizeData(T_others, normalization.method = "LogNormalize", scale.factor = 10000)
# T_others <- ScaleData(T_others, assay = "RNA", features = rownames(T_others))
T_others <- SCTransform(T_others)
DefaultAssay(T_others) <- 'SCT'
T_others <- FindVariableFeatures(T_others, selection.method = "vst", nfeatures = 1000)
T_others <- RunPCA(T_others)
T_others <- FindNeighbors(T_others, dims = 1:12)
T_others <- RunUMAP(T_others, dims = 1:10)
T_others <- RunTSNE(T_others, dims = 1:10)
T_others  <- FindClusters(T_others, resolution = 1:9/100)
T_others  <- FindClusters(T_others, resolution = 1:10/10)
clustree(T_others)
my_plotDim(T_others, reduction = "umap", label = F, pt.size = 0.1, label.size = 5, group.by = 'orig.ident')
# 0.5
T_others  <- FindClusters(T_others, resolution = 0.5)
my_plotDim(T_others, reduction = "umap", label = F, pt.size = 0.1, label.size = 5)
my_CountCluster(T_others, group2 = 'sub3')
my_DotPlot_split(T_others, features = T_marker$CD8) + RotatedAxis()
my_violin(T_others, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx')
# cell type identify
T_others$T_others_MainCluster <- as.factor(T_others$SCT_snn_res.0.5)
T_others$T_others <- T_others$SCT_snn_res.0.5
T_others_MainCluster <- rbind(
    data.frame(x='CD8_Tem',y=c(9)),
    data.frame(x='CD8_Tea',y=c(2,3,4,5,8)),
    data.frame(x='CD8_Naive',y=c(0,1,7)),
    data.frame(x='CD8_Tex',y=c(6))
)
T_others_MainCluster <- T_others_MainCluster[order(T_others_MainCluster$y, decreasing = F),]
if(all(T_others_MainCluster$y == levels(T_others$T_others_MainCluster))) {
    levels(T_others$T_others_MainCluster) <- T_others_MainCluster$x
}
Idents(T_others) <- 'T_others_MainCluster'
my_plotDim(T_others, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_others_MainCluster")
my_plotDim(T_others, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_others_MainCluster")
T_sub$T_others_T_others_MainCluster <- NULL
T_sub$T_others_T_others <- NULL

T_sub <- my_AddMeta(T_sub, T_others$T_others_MainCluster, Replace = T)
T_sub <- my_AddMeta(T_sub, T_others$T_others, Replace = T)
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_others_T_others_MainCluster")
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_others_T_others_MainCluster")
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_others_T_others")
my_CountCluster(T_sub, group1 = 'sub3', group2 = 'T_others_T_others')

# Cluster Th1("3_0",15) : Th1/Th2
Th1 <- subset(T_sub, ident = 'Th1')
# Re-ReCluster
DefaultAssay(Th1) <- 'RNA'
# Th1 <- NormalizeData(Th1, normalization.method = "LogNormalize", scale.factor = 10000)
# Th1 <- ScaleData(Th1, assay = "RNA", features = rownames(Th1))
Th1 <- SCTransform(Th1)
DefaultAssay(Th1) <- 'SCT'
Th1 <- FindVariableFeatures(Th1, selection.method = "vst", nfeatures = 1000)
Th1 <- RunPCA(Th1)
Th1 <- FindNeighbors(Th1, dims = 1:12)
Th1 <- RunUMAP(Th1, dims = 1:10)
Th1 <- RunTSNE(Th1, dims = 1:10)
#
Th1  <- FindClusters(Th1, resolution = 1:9/100)
Th1  <- FindClusters(Th1, resolution = 1:10/10)
clustree(Th1)
my_plotDim(Th1, reduction = "umap", label = F, pt.size = 0.1, label.size = 5, group.by = 'orig.ident')
Th1  <- FindClusters(Th1, resolution = 0.5)
my_plotDim(Th1, reduction = "umap", label = F, pt.size = 0.1, label.size = 5)
my_CountCluster(Th1, group2 = 'sub3')
my_DotPlot_split(Th1, features = T_marker$CD4) + RotatedAxis()
my_DotPlot_split(Th1, features = c(T_marker$Th,T_marker$Cyto_T)) + RotatedAxis()
my_DotPlot_split(Th1, features = T_marker$CD8) + RotatedAxis()
my_violin(Th1, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx')
# cell type identify
Th1$Th1_MainCluster <- as.factor(Th1$SCT_snn_res.0.5)
Th1$Th1 <- Th1$SCT_snn_res.0.5
Th1_MainCluster <- rbind(
    data.frame(x='Tfh',y=c(0,3)),
    data.frame(x='Th1',y=c(1,2,5)),
    data.frame(x='Th2',y=c(4))
)
Th1_MainCluster <- Th1_MainCluster[order(Th1_MainCluster$y, decreasing = F),]
if(all(Th1_MainCluster$y == levels(Th1$Th1_MainCluster))) {
    levels(Th1$Th1_MainCluster) <- Th1_MainCluster$x
}
Idents(Th1) <- 'Th1_MainCluster'
my_plotDim(Th1, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "Th1_MainCluster")
my_plotDim(Th1, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "Th1_MainCluster")
T_sub$Th1_Th1_MainCluster <- NULL
T_sub$Th1_Th1 <- NULL

Idents(T_sub) <- 'T_others_T_others_MainCluster'
T_sub <- my_AddMeta(T_sub, Th1$Th1_MainCluster, Replace = T)
T_sub <- my_AddMeta(T_sub, Th1$Th1, Replace = T)
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "Th1_Th1_MainCluster")
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "Th1_Th1_MainCluster")
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "Th1_Th1")
my_CountCluster(T_sub, group1 = 'sub3', group2 = 'Th1_Th1')

# integrate clusters
my_plotDim(T_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "sub3")

T_sub_tmp <- T_sub
T_sub_tmp$tmp <- as.factor(T_sub_tmp$sub3)
Idents(T_sub_tmp) <- 'tmp'
T_sub_tmp$T_others_T_others_MainCluster <- NULL
T_sub_tmp <- my_AddMeta(T_sub_tmp, T_others$T_others_MainCluster, Replace = T)
Idents(T_sub_tmp) <- 'T_others_T_others_MainCluster'
T_sub_tmp$Th1_Th1_MainCluster <- NULL
T_sub_tmp <- my_AddMeta(T_sub_tmp, Th1$Th1_MainCluster, Replace = T)
my_plotDim(T_sub_tmp, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "Th1_Th1_MainCluster")

T_sub_tmp$T_sub_MainCluster2 <- factor(T_sub_tmp$Th1_Th1_MainCluster)
T_sub_MainCluster2_df <- data.frame(x = 0:17,
                                    y = c('Tfh','Th1','Th2','3_1',4:13,'CD8_Naive','CD8_Tea','CD8_Tem','CD8_Tex')
                                    )
T_sub_MainCluster2_df <- T_sub_MainCluster2_df[order(T_sub_MainCluster2_df$y, decreasing = F),]
if(all(T_sub_MainCluster2_df$y == levels(T_sub_tmp$T_sub_MainCluster2))) {
    levels(T_sub_tmp$T_sub_MainCluster2) <- T_sub_MainCluster2_df$x
    T_sub_tmp$T_sub_MainCluster2 <- factor(T_sub_tmp$T_sub_MainCluster2, levels = 0:17)
}
Idents(T_sub_tmp) <- 'T_sub_MainCluster2'
my_plotDim(T_sub_tmp, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster2")
my_plotDim(T_sub_tmp, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster2")

T_sub_tmp$T_sub_MainCluster <- factor(T_sub_tmp$T_sub_MainCluster2)
T_sub_MainCluster_df <- rbind(
    data.frame(x='CD4_Naive',y=c(13)),
    data.frame(x='CD8_Naive',y=c(14)),
    data.frame(x='Treg',y=c(12)),
    data.frame(x='Th17',y=c(8)),
    data.frame(x='Th1',y=c(1)),
    data.frame(x='Th2',y=c(2)),
    data.frame(x='Tfh',y=c(0,3)),
    data.frame(x='CD8_Tem',y=c(4,5,9,16)),
    data.frame(x='CD8_Tea',y=c(7,15)),
    data.frame(x='CD8_Tex',y=c(6,10,11,17))
)
T_sub_MainCluster_df <- T_sub_MainCluster_df[order(T_sub_MainCluster_df$y, decreasing = F),]
if(all(T_sub_MainCluster_df$y == levels(T_sub_tmp$T_sub_MainCluster))) {
    levels(T_sub_tmp$T_sub_MainCluster) <- T_sub_MainCluster_df$x
}
Idents(T_sub_tmp) <- 'T_sub_MainCluster'
my_plotDim(T_sub_tmp, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster")
my_plotDim(T_sub_tmp, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster")

my_DotPlot_split(T_sub_tmp, features = T_markers_tmp) + RotatedAxis()
my_violin(T_sub_tmp, features = unlist(T_markers_tmp), pt.size = 0, mode = 'mtx')

T_sub$T_sub_MainCluster <- T_sub_tmp$T_sub_MainCluster
T_sub$T_sub_MainCluster2 <- T_sub_tmp$T_sub_MainCluster2

# ProjecTIL project
library(ProjecTILs)
ref_TILAtlas <- readRDS("E://Cryo-TCR/TILAtlas/ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref_TILAtlas, label = T, cols = refCols)

Idents(T_sub) <- 'T_sub_MainCluster'
DefaultAssay(T_sub) <- "SCT"
query_T <- make.projection(T_sub, ref = ref_TILAtlas, skip.normalize = T, fast.mode = T, filter.cells = F)
DimPlot(query_T, group.by = 'T_sub_MainCluster', cols = c(refCols,"black"))
DimPlot(query_T, group.by = 'T_sub_MainCluster2', reduction = 'umap', label = T)

plot.projection(ref_TILAtlas, query_T)

# query.projected <- cellstate.predict(ref = ref_TILAtlas, query = query_T)
# table(query.projected$functional.cluster)
# plot.statepred.composition(ref_TILAtlas, query.projected, metric = "Percent")
# plot.states.radar(ref_TILAtlas, query = query.projected, min.cells = 30)

markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(T_sub, features = markers,stack = T, flip = T, assay = "RNA")

##############################
if (!dir.exists("combined_analysis"))
    dir.create('combined_analysis')
setwd('combined_analysis')

# T_sub_SCT <- T_sub
# T_sub_SCT@assays$RNA@data <- T_sub_SCT@assays$RNA@data[1:10,1:10]
# T_sub_SCT@assays$RNA@scale.data <- T_sub_SCT@assays$RNA@scale.data[1:10,1:10]
# save(T_sub_SCT, file = "Cryo_T_real.rda")

# load("Cryo_T_real.rda")
# T_sub <- T_sub_SCT
# rm(T_sub_SCT)
# DefaultAssay(T_sub) <- 'RNA'
# T_sub <- NormalizeData(T_sub, normalization.method = "LogNormalize", scale.factor = 10000)
# T_sub <- ScaleData(T_sub, assay = "RNA", features = rownames(T_sub))
# DefaultAssay(T_sub) <- 'SCT'

##############################
# output result

# T DimPlot
T_dim_umap <- my_plotDim(T_sub, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster', title = 'T Cluster')
ggsave(filename = "T_dim_umap.pdf", plot = T_dim_umap)

T_dim_tsne <- my_plotDim(T_sub, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster', title = 'T Cluster')
ggsave(filename = "T_dim_tsne.pdf", plot = T_dim_tsne)

# T VlnPlot
T_Vln <- my_violin(T_sub, features = unlist(T_markers_tmp2), pt.size = 0, mode = 'mtx', group.by = 'T_sub_MainCluster')
ggsave(filename = "T_vln.pdf", plot = T_Vln)

# T DotPlot
T_Dot <- my_DotPlot_split(T_sub, features = T_markers_tmp2, group.by = 'T_sub_MainCluster') + RotatedAxis()
ggsave(filename = "T_Dot.pdf", plot = T_Dot, width = 16, height = 7)

# T annotation DimPlot
T_dim_umap2 <- my_plotDim(T_sub, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster2', title = 'T Cluster')
ggsave(filename = "T_anno_dim_umap.pdf", plot = T_dim_umap2)

T_dim_tsne2 <- my_plotDim(T_sub, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster2', title = 'T Cluster')
ggsave(filename = "T_anno_dim_tsne.pdf", plot = T_dim_tsne2)

# T annotation VlnPlot
T_Vln2 <- my_violin(T_sub, features = unlist(T_markers_tmp2), pt.size = 0, mode = 'mtx', group.by = 'T_sub_MainCluster2')
ggsave(filename = "T_anno_vln.pdf", plot = T_Vln2)

# T annotation DotPlot
T_Dot2 <- my_DotPlot_split(T_sub, features = T_markers_tmp2, group.by = 'T_sub_MainCluster2') + RotatedAxis()
ggsave(filename = "T_anno_Dot.pdf", plot = T_Dot2, width = 16, height = 7)

# T count
T_count <- as.data.frame(my_CountCluster(T_sub, group1 = 'T_sub_MainCluster', group2 = 'orig.ident')[-11, 5:8])
T_count$cell <- rownames(T_count)
T_count <- reshape::melt(T_count, id.vars = 'cell')
T_count$variable <- sub("_percentage$", "", T_count$variable)
T_count$variable <- sapply(T_count$variable, function(x) sub("^NonCryo", "Non-CA-", sub("^Cryo", "CA-", x)))
T_count$variable <- factor(T_count$variable, levels = c("Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk"))
T_Plot <- ggplot(T_count, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(T_count)/4 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "T_count.pdf", plot = T_Plot)

# T annotation count
T_count2 <- as.data.frame(my_CountCluster(T_sub, group1 = 'T_sub_MainCluster2', group2 = 'orig.ident')[-19, 5:8])
T_count2$cell <- rownames(T_count2)
T_count2 <- reshape::melt(T_count2, id.vars = 'cell')
T_count2$variable <- sub("_percentage$", "", T_count2$variable)
T_count2$variable <- sapply(T_count2$variable, function(x) sub("^NonCryo", "Non-CA-", sub("^Cryo", "CA-", x)))
T_count2$variable <- factor(T_count2$variable, levels = c("Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk"))
T_count2$cell <- factor(T_count2$cell, levels = 0:17)
T_Plot2 <- ggplot(T_count2, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cluster") +
    scale_fill_manual(values=c("grey75", colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nrow(T_count2)/4 - 1))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "T_anno_count.pdf", plot = T_Plot2)




