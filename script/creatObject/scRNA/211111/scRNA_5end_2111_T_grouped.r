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

merge_object <- function(Cryo, NonCryo) {
    Cryo@meta.data$orig.ident <- Cryo@project.name
    NonCryo@meta.data$orig.ident <- NonCryo@project.name
    Cryo_combined <- merge(Cryo, NonCryo)
    Cryo_combined[["percent.mt"]] <- PercentageFeatureSet(Cryo_combined, pattern = "^mt-")
    Idents(Cryo_combined) <- 'orig.ident'
    Cryo_combined$orig.ident_group <- Cryo_combined$orig.ident
    Cryo_combined$orig.ident_group <- sub("^Cryo.*$", "Cryo", Cryo_combined$orig.ident_group)
    Cryo_combined$orig.ident_group <- sub("^NonCryo.*$", "NonCryo", Cryo_combined$orig.ident_group)
    Cryo_combined$orig.ident_time <- Cryo_combined$orig.ident
    Cryo_combined$orig.ident_time <- sub("^.*1wk$", "1wk", Cryo_combined$orig.ident_time)
    Cryo_combined$orig.ident_time <- sub("^.*2wk$", "2wk", Cryo_combined$orig.ident_time)
    return(Cryo_combined)
}
Cryo_combined1 <- merge_object(Cryo_1, NonCryo_1)
Cryo_combined2 <- merge_object(Cryo_2, NonCryo_2)
Cryo_combined1
Cryo_combined2

QC_vis <- function(Cryo_combined) {
    # visualization
    print(my_plotQC(Cryo_combined, pt.size = 0))
    # print(my_plotFC(Cryo_combined, combinedIf = F, pt.size = 0.5))
    hist(Cryo_combined$nFeature_RNA, breaks = 500, xaxt = 'n')
    axis(1, seq(0,10000,200))
    hist(Cryo_combined$nCount_RNA, breaks = 500, xaxt = 'n')
    axis(1, seq(0,200000,2000))
    hist(Cryo_combined$percent.mt, breaks = 500, xaxt = 'n')
    axis(1, seq(0,100,1))
}
QC_vis(Cryo_combined1)
QC_vis(Cryo_combined2)

# quality filter
Cryo_merge1 <- subset(Cryo_combined1, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA > 800 & nCount_RNA < 15000 & percent.mt < 8)
Cryo_merge1
Cryo_merge2 <- subset(Cryo_combined2, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & nCount_RNA > 800 & nCount_RNA < 15000 & percent.mt < 8)
Cryo_merge2

################################################################################
process_obj <- function(Cryo_merge, nVariableFeatures = 2500, SCT = FALSE) {
    DefaultAssay(Cryo_merge) <- 'RNA'
    if (SCT) {
        Cryo_merge <- SCTransform(Cryo_merge)
        DefaultAssay(Cryo_merge) <- 'SCT'
    } else {
        Cryo_merge <- NormalizeData(Cryo_merge, normalization.method = "LogNormalize", scale.factor = 10000)
        Cryo_merge <- ScaleData(Cryo_merge, assay = "RNA", features = rownames(Cryo_merge))
    }
    Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = nVariableFeatures)
    Cryo_merge <- RunPCA(Cryo_merge)
    Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
    Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10)
    Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10)
}

# process for data clean (remove non-T)
Cryo_merge1 <- process_obj(Cryo_merge1, 2500, FALSE)
Cryo_merge2 <- process_obj(Cryo_merge2, 2500, FALSE)

#
Cryo_merge1  <- FindClusters(Cryo_merge1, resolution = c(1:9/100,1:10/10))
clustree(Cryo_merge1)
Cryo_merge2  <- FindClusters(Cryo_merge2, resolution = c(1:9/100,1:10/10))
clustree(Cryo_merge2)

# Cryo_merge1 RNA_snn_res.0.05 C1: non-T
Cryo_merge1  <- FindClusters(Cryo_merge1, resolution = 0.05)
Idents(Cryo_merge1) <- 'RNA_snn_res.0.05'
my_plotDim(Cryo_merge1, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge1 R=0.05', group.by = 'RNA_snn_res.0.05')
my_plotDim(Cryo_merge1, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge1 R=0.05', group.by = 'RNA_snn_res.0.05')
my_DotPlot_split(Cryo_merge1, features = my_MarkersList[1:7]) + RotatedAxis()
my_violin(Cryo_merge1, features = unlist(my_MarkersList[1:7]), pt.size = 0, mode = 'mtx')
DoHeatmap(Cryo_merge1, features = c('Ptprc','Pecam1','Cd3e','Cd4','Cd8a','Cd8b1','Itgam','Cd14'), assay = 'RNA')

# Cryo_merge2 RNA_snn_res.0.05 C3/5: non-T
Cryo_merge2  <- FindClusters(Cryo_merge2, resolution = 0.06)
Idents(Cryo_merge2) <- 'RNA_snn_res.0.06'
my_plotDim(Cryo_merge2, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge2 R=0.05', group.by = 'RNA_snn_res.0.06')
my_plotDim(Cryo_merge2, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'Cryo_merge2 R=0.05', group.by = 'RNA_snn_res.0.06')
my_DotPlot_split(Cryo_merge2, features = my_MarkersList[1:7]) + RotatedAxis()
my_violin(Cryo_merge2, features = unlist(my_MarkersList[1:7]), pt.size = 0, mode = 'mtx')
DoHeatmap(Cryo_merge2, features = c('Ptprc','Pecam1','Cd3e','Cd4','Cd8a','Cd8b1','Itgam','Cd14'), assay = 'RNA')

# cell filter
T1 <- subset(Cryo_merge1, ident = c(0,2:4))
T2 <- subset(Cryo_merge2, ident = c(0:2,4))

# ReCluster
T1 <- process_obj(T1, 3000, SCT = FALSE)
T2 <- process_obj(T2, 3000, SCT = FALSE)


################################################################################
# process for T1:
T1 <- FindClusters(T1, resolution = c(1:9/100, 1:10/10))
clustree(T1, prefix = "RNA_snn_res.")

my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 5, group.by = 'orig.ident')

Idents(T1) <- 'RNA_snn_res.0.6'
my_plotDim(T1, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T1 R=0.6', group.by = 'RNA_snn_res.0.6')
my_plotDim(T1, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T1 R=0.6', group.by = 'RNA_snn_res.0.6')

my_DotPlot_split(T1, features = my_MarkersList[1:7]) + RotatedAxis()

# cluster 13: B
# CD4: 3/9/11
# CD8: 0/1/2/4/5/6/7/8/10/12
T1$MainCluster <- T1$RNA_snn_res.0.6
levels(T1$MainCluster) <- c(rep("CD8",3),"CD4",rep("CD8",5),"CD4","CD8","CD4","CD8","B")
Idents(T1) <- 'MainCluster'

T1_real <- subset(T1, ident = c("CD4","CD8"))
my_plotDim(T1_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T1 R=0.6', group.by = 'RNA_snn_res.0.6')
my_plotDim(T1_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T1 R=0.6', group.by = 'MainCluster')

Idents(T1_real) <- 'RNA_snn_res.0.6'
my_DotPlot_split(T1_real, features = T_markers_tmp) + RotatedAxis()
my_violin(T1_real, features = unlist(T_markers_tmp), pt.size = 0, mode = 'mtx')

my_violin(T1_real, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = c(3,9,11))
my_violin(T1_real, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = c(0,1,2,4,5,6,7,8,10,12))
my_DotPlot_split(T1_real, features = T_marker$CD4, idents = c(3,9,11)) + RotatedAxis()
my_DotPlot_split(T1_real, features = T_marker$CD8, idents = c(0,1,2,4,5,6,7,8,10,12)) + RotatedAxis()

# 11: Treg
# 9: Naive/T17
# 3: Tfh/Th1/Th2

# 1/6/10: Tex
# 0/4/5/7: Tem
# 2/8/12: Tea

DoHeatmap(T1_real, features = unlist(T_marker$CD4), cells = which(T1_real$RNA_snn_res.0.6 %in% c(3,9) )) # Tregs also express the markers of Th1/2

my_plotDim(T1_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T1 R=0.6', group.by = 'RNA_snn_res.0.6')

# Recluster for 3/9/11
T1_real$tmp2 <- factor(T1_real$RNA_snn_res.0.6)
levels(T1_real$tmp2) <- c(0:8,3,10,3,12)
Idents(T1_real) <- "tmp2"
T_sub3 <- subset(T1_real, ident = 3)
T_sub3 <- process_obj(T_sub3, 2000, SCT = FALSE)
for (i in colnames(T_sub3@meta.data)) if(length(grep("^RNA", i)) == 1) T_sub3[[i]] <- NULL
T_sub3 <- FindClusters(T_sub3, resolution = 1:5/10)
clustree::clustree(T_sub3)
Idents(T_sub3) <- "RNA_snn_res.0.3"
my_plotDim(T_sub3, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub3 R=0.3', group.by = 'RNA_snn_res.0.3')

my_DotPlot_split(T_sub3, features = my_MarkersList[1:7]) + RotatedAxis() # 5: CD8
my_DotPlot_split(T_sub3, features = T_marker$CD4) + RotatedAxis() # 4:Naive; 3:Treg; 1:Tfh; 6:Th2; 7: Th17; 0:Th1; 2:T_others
my_DotPlot_split(T_sub3, features = T_marker$CD8) + RotatedAxis() # 5:CD8_Tex

T_sub3 <- FindSubCluster(T_sub3, cluster = "2", subcluster.name = "sub2", resolution = 0.3, graph.name = "RNA_snn")
my_plotDim(T_sub3, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub3 sub2 R=0.3', group.by = 'sub2')
my_plotDim(T_sub3, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'T_sub3 sub2 R=0.3', group.by = 'sub2')
Idents(T_sub3) <- "sub2"

T_sub3$MainCluster_tmp <- factor(T_sub3$sub2)
levels(T_sub3$MainCluster_tmp) <- c("Th1", "Tfh", "Th1", "Th17", "Treg", "T_Naive", "CD8_Tex", "Th2", "Th17")
my_plotDim(T_sub3, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_tmp', title = 'T_sub3 Cluster')

T1_real <- my_AddMeta(T1_real, T_sub3$MainCluster_tmp, Replace = T)
my_plotDim(T1_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub3_MainCluster_tmp')
my_plotDim(T1_real, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub3_MainCluster_tmp')


# Recluster for 0,1,2,4,5,6,7,8,10,12

Idents(T1_real) <- "RNA_snn_res.0.6"
T_sub8 <- subset(T1_real, ident = c(0,1,2,4,5,6,7,8,10,12))
T_sub8 <- process_obj(T_sub8, 2000, SCT = FALSE)
for (i in colnames(T_sub8@meta.data)) if(length(grep("^RNA", i)) == 1) T_sub8[[i]] <- NULL
T_sub8 <- FindClusters(T_sub8, resolution = 1:5/10)
clustree::clustree(T_sub8)

Idents(T_sub8) <- "RNA_snn_res.0.3"
my_plotDim(T_sub8, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub8 R=0.3', group.by = 'RNA_snn_res.0.3')

my_DotPlot_split(T_sub8, features = my_MarkersList[1:7]) + RotatedAxis() # 8: NKT(Cd3e/Klrb1c/Klrk1)
my_DotPlot_split(T_sub8, features = T_marker$CD8) + RotatedAxis() # 1/4/6:CD8_Tex; 0/3:CD8_Tem; 2/5/7:CD8_Tea

T_sub8$MainCluster_tmp <- factor(T_sub8$RNA_snn_res.0.3)
levels(T_sub8$MainCluster_tmp) <- c("CD8_Tem", "CD8_Tex", "CD8_Tea", "CD8_Tem", "CD8_Tex", "CD8_Tea", "CD8_Tex", "CD8_Tea", "NKT")
my_plotDim(T_sub8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_tmp', title = 'T_sub8 Cluster')
my_plotDim(T_sub8, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_tmp', title = 'T_sub8 Cluster')

Idents(T1_real) <- "T_sub3_MainCluster_tmp"
T1_real <- my_AddMeta(T1_real, T_sub8$MainCluster_tmp, Replace = T)
my_plotDim(T1_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub8_MainCluster_tmp')
my_plotDim(T1_real, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub8_MainCluster_tmp')

# change Clusters
T1_real$MainCluster_anno <- T1_real$T_sub8_MainCluster_tmp
T1_real$MainCluster_C <- factor(T1_real$MainCluster_anno)
levels(T1_real$MainCluster_C) <- paste0("C", seq(levels(T1_real$MainCluster_C))) # need to be changed
Idents(T1) <- "MainCluster"
T1 <- my_AddMeta(T1, T1_real$MainCluster_anno, Replace = T)
T1$MainCluster_anno <- T1$T1_real_MainCluster_anno
T1$MainCluster_C <- factor(T1$MainCluster_anno)
levels(T1$MainCluster_C) <- paste0("C", seq(levels(T1$MainCluster_C)))

my_plotDim(T1_real, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_C')
my_plotDim(T1_real, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_anno')
my_plotDim(T1, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_anno')
my_plotDim(T1, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_C')

my_CountCluster(T1_real, group1 = "MainCluster_anno", group2 = 'orig.ident')
my_CountCluster(T1, group1 = "MainCluster_anno", group2 = 'orig.ident')


################################################################################
# process for T2:
T2 <- FindClusters(T2, resolution = c(1:9/100, 1:10/10))
clustree(T2, prefix = "RNA_snn_res.")

my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 5, group.by = 'orig.ident')

cluster_X <- 'RNA_snn_res.0.4'
Idents(T2) <- cluster_X
my_plotDim(T2, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T2 R=0.3', group.by = cluster_X)
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T2 R=0.3', group.by = cluster_X)

my_DotPlot_split(T2, features = my_MarkersList[1:7]) + RotatedAxis() 
my_DotPlot_split(T2, features = T_markers_tmp[1:2]) + RotatedAxis() # CD4:3/6; CD8:0/1/2/5/7/9/12; Chaos:4/8/10/11
my_violin(T2, features = unlist(my_MarkersList[1:7]), mode = 'mtx', pt.size = 0) + RotatedAxis()

# Recluster for 4/8/10/11
Idents(T2) <- "RNA_snn_res.0.4"
T_sub7 <- subset(T2, ident = c(4,8,10,11))
T_sub7 <- process_obj(T_sub7, 1000, SCT = FALSE)
for (i in colnames(T_sub7@meta.data)) if(length(grep("^RNA", i)) == 1) T_sub7[[i]] <- NULL
T_sub7 <- FindClusters(T_sub7, resolution = 1:5/10)
clustree::clustree(T_sub7)

Idents(T_sub7) <- "RNA_snn_res.0.5"
my_plotDim(T_sub7, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub7 R=0.5', group.by = 'RNA_snn_res.0.5')
FeaturePlot(T_sub7, features = c("Cd4", "Cd8a", "Cd8b1"))

my_DotPlot_split(T_sub7, features = my_MarkersList[1:7]) + RotatedAxis() # 7:NKT; 1/3/5/6:CD8; 2/4:CD4; 0:Chaos

T_sub7 <- FindSubCluster(T_sub7, cluster = '0', graph.name = 'RNA_snn', resolution = 0.2, subcluster.name = 'sub0')
Idents(T_sub7) <- "sub0"
my_plotDim(T_sub7, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub7 R=0.5', group.by = 'sub0')

# 2/4/0_1:CD4; 1/3/5/6/0_0:CD8; 7:NKT
T_sub7$MainCluster_tmp <- factor(T_sub7$sub0)
levels(T_sub7$MainCluster_tmp) <- c("CD8", "CD4", "CD8", "CD4", "CD8", "CD4", "CD8", "CD8", "NKT")
my_plotDim(T_sub7, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_tmp', title = 'T_sub7 Cluster')

T2 <- my_AddMeta(T2, T_sub7$MainCluster_tmp, Replace = T)
my_plotDim(T2, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub7_MainCluster_tmp')
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub7_MainCluster_tmp')
Idents(T2) <- "T_sub7_MainCluster_tmp"

my_DotPlot_split(T2, features = T_markers_tmp[1:2]) + RotatedAxis() # CD4:3/6/CD4; CD8:0/1/2/5/7/9/12/CD8

T2$CD48 <- factor(T2$T_sub7_MainCluster_tmp)
levels(T2$CD48) <- c("CD8","CD8","CD8","CD8","CD4","CD8","CD4","CD8","CD8","CD4","CD8","NKT")
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'CD48')


Idents(T2) <- "CD48"
# CD4:
T2_CD4 <- subset(T2, ident = "CD4")
T2_CD4 <- process_obj(T2_CD4, 2000, SCT = FALSE)

T2_CD4 <- FindClusters(T2_CD4, resolution = c(1:9/100, 1:10/10))
clustree(T2_CD4, prefix = "RNA_snn_res.")

my_plotDim(T2_CD4, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = 'orig.ident')

cluster_X <- 'RNA_snn_res.0.6'
Idents(T2_CD4) <- cluster_X
my_plotDim(T2_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T2_CD4 R=0.3', group.by = cluster_X)
my_plotDim(T2_CD4, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'T2_CD4 R=0.3', group.by = cluster_X)

my_DotPlot_split(T2_CD4, features = my_MarkersList[1:7]) + RotatedAxis() # 4:CD8
my_violin(T2_CD4, features = my_MarkersList[1:7], pt.size = 0, mode = 'mtx')

my_violin(T2_CD4, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = setdiff(0:8,4))
my_violin(T2_CD4, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = c(4))
my_DotPlot_split(T2_CD4, features = T_marker$CD4, idents = setdiff(0:8,4)) + RotatedAxis() # 5:Naive; 1/8:Tfh; 7:Treg; 3/6:Th17; 0/2:Th1
my_DotPlot_split(T2_CD4, features = T_marker$CD8, idents = c(4)) + RotatedAxis()
DoHeatmap(T2_CD4, features = c('Cd4','Cd8a','Cd8b1'), cells = which(T2_CD4$RNA_snn_res.0.6 == 4))
mtx_tmp <- T2_CD8@assays$RNA@data[c('Cd4','Cd8a','Cd8b1'), T2_CD4$RNA_snn_res.0.6 == 4]
pheatmap::pheatmap(mtx_tmp, show_colnames = F) # 4 all CD8

# change
T2_CD4$MainCluster_tmp <- T2_CD4$RNA_snn_res.0.6
levels(T2_CD4$MainCluster_tmp) <- c("CD4","CD4","CD4","CD4","CD8","CD4","CD4","CD4","CD4")
Idents(T2) <- "CD48"
T2 <- my_AddMeta(T2, T2_CD4$MainCluster_tmp, Replace = T)

# CD8:
T2_CD8 <- subset(T2, ident = "CD8")
T2_CD8 <- process_obj(T2_CD8, 2000, SCT = FALSE)

T2_CD8 <- FindClusters(T2_CD8, resolution = c(1:9/100, 1:10/10))
clustree(T2_CD8, prefix = "RNA_snn_res.")

my_plotDim(T2_CD8, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = 'orig.ident')

cluster_X <- 'RNA_snn_res.0.5'
Idents(T2_CD8) <- cluster_X
my_plotDim(T2_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T2_CD8 R=0.3', group.by = cluster_X)
my_plotDim(T2_CD8, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = 'T2_CD8 R=0.3', group.by = cluster_X)

my_DotPlot_split(T2_CD8, features = my_MarkersList[1:7]) + RotatedAxis()
my_violin(T2_CD8, features = unlist(my_MarkersList[1:7]), pt.size = 0, mode = 'mtx')

DoHeatmap(T2_CD8, features = c('Cd4','Cd8a','Cd8b1'), cells = which(T2_CD8$RNA_snn_res.0.5 %in% c(9,10,13,14)))
mtx_tmp <- T2_CD8@assays$RNA@data[c('Cd4','Cd8a','Cd8b1'), T2_CD8$RNA_snn_res.0.5 %in% c(9,10,13,14)]
names_bc_CD4 <- colnames(mtx_tmp)[apply(mtx_tmp, 2, function(x) if (x[1]>0 & x[2]==0 & x[3]==0) 1 else 0) == 1]
T2_CD8$tmp <- unfactor(T2_CD8$RNA_snn_res.0.5)
T2_CD8$tmp[names_bc_CD4] <- "CD4"
Idents(T2_CD8) <- "tmp"
T2_CD8$tmp[T2_CD8$tmp!="CD4"] <- "CD8"

# change
Idents(T2) <- "T2_CD4_MainCluster_tmp"
T2 <- my_AddMeta(T2, T2_CD8$tmp, Replace = T)
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = "T2_CD8_tmp")

T2$CD48_2 <- factor(T2$T2_CD8_tmp)

Idents(T2) <- "CD48_2"
# ReCluster for CD4
T2_CD4_2 <- subset(T2, ident = "CD4")
T2_CD4_2 <- process_obj(T2_CD4_2, 2000, SCT = FALSE)

T2_CD4_2 <- FindClusters(T2_CD4_2, resolution = c(1:9/100, 1:10/10))
clustree(T2_CD4_2, prefix = "RNA_snn_res.")

my_plotDim(T2_CD4_2, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = 'orig.ident')

cluster_X <- 'RNA_snn_res.0.5'
Idents(T2_CD4_2) <- cluster_X
my_plotDim(T2_CD4_2, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T2_CD4_2 R=0.3', group.by = cluster_X)
my_plotDim(T2_CD4_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T2_CD4_2 R=0.3', group.by = cluster_X)

# my_DotPlot_split(T2_CD4_2, features = my_MarkersList[1:7]) + RotatedAxis()
# my_violin(T2_CD4_2, features = my_MarkersList[1:7], pt.size = 0, mode = 'mtx')

my_violin(T2_CD4_2, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = setdiff(0:8,4))
my_DotPlot_split(T2_CD4_2, features = T_marker$CD4) + RotatedAxis() # 4:Naive; 1/7/8:Tfh; 5:Treg; 3/6:Th17; 0/2:Th1; chaos:0/1/2/3:Th1/Th2

# Recluster for 0/1/2/3
T2_CD4_sub0 <- subset(T2_CD4_2, ident = c(0:3))
T2_CD4_sub0 <- process_obj(T2_CD4_sub0, 1000, SCT = FALSE)

T2_CD4_sub0 <- FindClusters(T2_CD4_sub0, resolution = c(1:9/100, 1:10/10))
clustree(T2_CD4_sub0, prefix = "RNA_snn_res.")

cluster_X <- 'RNA_snn_res.0.5'
Idents(T2_CD4_sub0) <- cluster_X
my_plotDim(T2_CD4_sub0, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = cluster_X)
my_DotPlot_split(T2_CD4_sub0, features = T_marker$CD4) + RotatedAxis() # 2:Tfh; 5:Th2; 3:Th17; 0/1/4:Th1

T2_CD4_sub0$MainCluster_tmp <- T2_CD4_sub0$RNA_snn_res.0.5
levels(T2_CD4_sub0$MainCluster_tmp) <- c("Th1","Th1","Tfh","Th17","Th1","Th2")
my_plotDim(T2_CD4_sub0, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_tmp")

# change
T2_CD4_2 <- my_AddMeta(T2_CD4_2, T2_CD4_sub0$MainCluster_tmp, Replace = T)
my_plotDim(T2_CD4_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "T2_CD4_sub0_MainCluster_tmp")
T2_CD4_2$MainCluster_anno <- factor(T2_CD4_2$T2_CD4_sub0_MainCluster_tmp)
T2_CD4_2$MainCluster_C <- factor(T2_CD4_2$T2_CD4_sub0_MainCluster_tmp)
levels(T2_CD4_2$MainCluster_C) <- paste0("C", seq(levels(T2_CD4_2$MainCluster_C)))
levels(T2_CD4_2$MainCluster_anno) <- c("CD4_Naive", "Treg", "Th17", "Tfh", "Tfh", "Tfh", "Th1", "Th17", "Th2")

my_plotDim(T2_CD4_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")
my_plotDim(T2_CD4_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")


Idents(T2) <- "CD48_2"
# ReCluster for CD8
T2_CD8_2 <- subset(T2, ident = "CD8")
T2_CD8_2 <- process_obj(T2_CD8_2, 2000, SCT = FALSE)

T2_CD8_2 <- FindClusters(T2_CD8_2, resolution = c(1:9/100, 1:10/10))
clustree(T2_CD8_2, prefix = "RNA_snn_res.")

my_plotDim(T2_CD8_2, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = 'orig.ident')

cluster_X <- 'RNA_snn_res.0.5'
Idents(T2_CD8_2) <- cluster_X
my_plotDim(T2_CD8_2, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T2_CD8_2 R=0.5', group.by = cluster_X)
my_plotDim(T2_CD8_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, title = 'T2_CD8_2 R=0.5', group.by = cluster_X)

my_DotPlot_split(T2_CD8_2, features = my_MarkersList[1:7]) + RotatedAxis() # 10:NKT

my_violin(T2_CD8_2, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = setdiff(0:13,10))
my_DotPlot_split(T2_CD8_2, features = T_marker$CD8, idents = setdiff(0:13,10)) + RotatedAxis() # 0/8/11:CD8_Tea; 1/2/3/4/6:CD8_Tem; 5/7/9/12/13:CD8_Tex

T2_CD8_2$MainCluster_anno <- factor(T2_CD8_2$RNA_snn_res.0.5)
T2_CD8_2$MainCluster_C <- factor(T2_CD8_2$RNA_snn_res.0.5)
levels(T2_CD8_2$MainCluster_C) <- c(paste0("C", 1:10), "NKT", paste0("C", 11:13))
levels(T2_CD8_2$MainCluster_anno) <- c("CD8_Tea", "CD8_Tem", "CD8_Tem", "CD8_Tem", "CD8_Tem", "CD8_Tex", "CD8_Tem", "CD8_Tex", "CD8_Tea", "CD8_Tex", "NKT", "CD8_Tea", "CD8_Tex","CD8_Tex")

my_plotDim(T2_CD8_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")
my_plotDim(T2_CD8_2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")


# change Cluster
# MainCluster_anno
Idents(T2) <- "CD48_2"
T2 <- my_AddMeta(T2, T2_CD4_2$MainCluster_anno, Replace = T)
Idents(T2) <- "T2_CD4_2_MainCluster_anno"
T2 <- my_AddMeta(T2, T2_CD8_2$MainCluster_anno, Replace = T)
T2$MainCluster_anno <- T2$T2_CD8_2_MainCluster_anno
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")
my_plotDim(T2, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")
# MainCluster
Idents(T2) <- "CD48_2"
T2 <- my_AddMeta(T2, T2_CD4_2$MainCluster_C, Replace = T)
T2$T2_CD4_2_MainCluster_C[grep("^C[0-9]$", T2$T2_CD4_2_MainCluster_C)] <- paste0("CD4-", T2$T2_CD4_2_MainCluster_C[grep("^C[0-9]$", T2$T2_CD4_2_MainCluster_C)])
Idents(T2) <- "T2_CD4_2_MainCluster_C"
T2 <- my_AddMeta(T2, T2_CD8_2$MainCluster_C, Replace = T)
T2$T2_CD8_2_MainCluster_C[grep("^C[0-9]+$", T2$T2_CD8_2_MainCluster_C)] <- paste0("CD8-", T2$T2_CD8_2_MainCluster_C[grep("^C[0-9]+$", T2$T2_CD8_2_MainCluster_C)])
T2$MainCluster_C <- T2$T2_CD8_2_MainCluster_C
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")
my_plotDim(T2, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")


################################################################################
# all result
T1
T2

my_plotDim(T1, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")
my_plotDim(T1, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")

my_CountCluster(T1, group1 = "MainCluster_C")
my_CountCluster(T1, group1 = "MainCluster_anno")

my_plotDim(T2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_C")
my_plotDim(T2, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = "MainCluster_anno")

my_CountCluster(T2, group1 = "MainCluster_C")
my_CountCluster(T2, group1 = "MainCluster_anno")

# if (!dir.exists("grouped_analysis"))
#     dir.create("grouped_analysis")
# save(T1, T2, file = "grouped_analysis/Cryo_T_grouped.rda")
load(file = "grouped_analysis/Cryo_T_grouped.rda")








# ProjecTIL project
library(ProjecTILs)
ref_TILAtlas <- readRDS("E://Cryo-TCR/TILAtlas/ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref_TILAtlas, label = T, cols = refCols)

Idents(T1) <- 'MainCluster_anno'
DefaultAssay(T1) <- "RNA"
query_1 <- make.projection(T1, ref = ref_TILAtlas, skip.normalize = T, fast.mode = T, filter.cells = F)
# DimPlot(query_1, group.by = 'MainCluster_anno', cols = c(refCols,"black"))
DimPlot(query_1, group.by = 'MainCluster_anno', reduction = 'umap', label = T)
plot.projection(ref_TILAtlas, query_1)

Idents(T2) <- 'MainCluster_anno'
DefaultAssay(T2) <- "RNA"
query_2 <- make.projection(T2, ref = ref_TILAtlas, skip.normalize = T, fast.mode = T, filter.cells = F)
# DimPlot(query_2, group.by = 'MainCluster_anno', cols = c(refCols,"black"))
DimPlot(query_2, group.by = 'MainCluster_anno', reduction = 'umap', label = T)
plot.projection(ref_TILAtlas, query_2)

# # query.projected <- cellstate.predict(ref = ref_TILAtlas, query = query_T)
# # table(query.projected$functional.cluster)
# # plot.statepred.composition(ref_TILAtlas, query.projected, metric = "Percent")
# # plot.states.radar(ref_TILAtlas, query = query.projected, min.cells = 30)
# 
# markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
# VlnPlot(T_sub, features = markers,stack = T, flip = T, assay = "RNA")


##############################
# output result
setwd(dir = "grouped_analysis/")

# T1
# T DimPlot
T1_dim_umap <- my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '1wk T-cells Clusters')
ggsave(filename = "T1_dim_umap.pdf", plot = T1_dim_umap)

T1_dim_tsne <- my_plotDim(T1, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '1wk T-cells Clusters')
ggsave(filename = "T1_dim_tsne.pdf", plot = T1_dim_tsne)

# T VlnPlot
T1_Vln <- my_violin(T1, features = unlist(T_markers_tmp3), pt.size = 0, mode = 'mtx', group.by = 'MainCluster_C')
ggsave(filename = "T1_vln.pdf", plot = T1_Vln)

# T DotPlot
T1_Dot <- my_DotPlot_split(T1, features = T_markers_tmp3, group.by = 'MainCluster_C') + RotatedAxis()
ggsave(filename = "T1_Dot.pdf", plot = T1_Dot, width = 16, height = 7)

# T annotation DimPlot
T1_dim_umap2 <- my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '1wk T-cells Clusters')
ggsave(filename = "T1_anno_dim_umap.pdf", plot = T1_dim_umap2)

T1_dim_tsne2 <- my_plotDim(T1, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '1wk T-cells Clusters')
ggsave(filename = "T1_anno_dim_tsne.pdf", plot = T1_dim_tsne2)

# T annotation VlnPlot
T1_Vln2 <- my_violin(T1, features = unlist(T_markers_tmp3), pt.size = 0, mode = 'mtx', group.by = 'MainCluster_anno')
ggsave(filename = "T1_anno_vln.pdf", plot = T1_Vln2)

# T annotation DotPlot
T1_Dot2 <- my_DotPlot_split(T1, features = T_markers_tmp3, group.by = 'MainCluster_anno') + RotatedAxis()
ggsave(filename = "T1_anno_Dot.pdf", plot = T1_Dot2, width = 16, height = 7)

# # T count
# T_count <- as.data.frame(my_CountCluster(T_sub, group1 = 'T_sub_MainCluster', group2 = 'orig.ident')[-11, 5:8])
# T_count$cell <- rownames(T_count)
# T_count <- reshape::melt(T_count, id.vars = 'cell')
# T_count$variable <- sub("_percentage$", "", T_count$variable)
# T_count$variable <- sapply(T_count$variable, function(x) sub("^NonCryo", "Non-CA-", sub("^Cryo", "CA-", x)))
# T_count$variable <- factor(T_count$variable, levels = c("Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk"))
# T_Plot <- ggplot(T_count, aes(x = variable, y = value, fill = cell)) +
#     geom_bar(width = 0.9, stat = "identity") +
#     xlab("Group") +
#     ylab("Proportion, %") +
#     labs(fill="Cell Type") +
#     scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(T_count)/4 - 1, "Paired"))) +
#     theme_bw() +
#     theme(legend.text=element_text(size=8), 
#           axis.title.x=element_text(size=12), 
#           axis.title.y=element_text(size=12), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     guides(fill = guide_legend(reverse = FALSE))
# ggsave(filename = "T_count.pdf", plot = T_Plot)
# 
# # T annotation count
# T_count2 <- as.data.frame(my_CountCluster(T_sub, group1 = 'T_sub_MainCluster2', group2 = 'orig.ident')[-19, 5:8])
# T_count2$cell <- rownames(T_count2)
# T_count2 <- reshape::melt(T_count2, id.vars = 'cell')
# T_count2$variable <- sub("_percentage$", "", T_count2$variable)
# T_count2$variable <- sapply(T_count2$variable, function(x) sub("^NonCryo", "Non-CA-", sub("^Cryo", "CA-", x)))
# T_count2$variable <- factor(T_count2$variable, levels = c("Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk"))
# T_count2$cell <- factor(T_count2$cell, levels = 0:17)
# T_Plot2 <- ggplot(T_count2, aes(x = variable, y = value, fill = cell)) +
#     geom_bar(width = 0.9, stat = "identity") +
#     xlab("Group") +
#     ylab("Proportion, %") +
#     labs(fill="Cluster") +
#     scale_fill_manual(values=c("grey75", colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nrow(T_count2)/4 - 1))) +
#     theme_bw() +
#     theme(legend.text=element_text(size=8), 
#           axis.title.x=element_text(size=12), 
#           axis.title.y=element_text(size=12), 
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(), 
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black")) +
#     guides(fill = guide_legend(reverse = FALSE))
# ggsave(filename = "T_anno_count.pdf", plot = T_Plot2)

# T2
# T DimPlot
T2_dim_umap <- my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '2wk T-cells Clusters')
ggsave(filename = "T2_dim_umap.pdf", plot = T2_dim_umap)

T2_dim_tsne <- my_plotDim(T2, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '2wk T-cells Clusters')
ggsave(filename = "T2_dim_tsne.pdf", plot = T2_dim_tsne)

# T VlnPlot
T2_Vln <- my_violin(T2, features = unlist(T_markers_tmp3), pt.size = 0, mode = 'mtx', group.by = 'MainCluster_C')
ggsave(filename = "T2_vln.pdf", plot = T2_Vln)

# T DotPlot
T2_Dot <- my_DotPlot_split(T2, features = T_markers_tmp3, group.by = 'MainCluster_C') + RotatedAxis()
ggsave(filename = "T2_Dot.pdf", plot = T2_Dot, width = 16, height = 7)

# T annotation DimPlot
T2_dim_umap2 <- my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '2wk T-cells Clusters')
ggsave(filename = "T2_anno_dim_umap.pdf", plot = T2_dim_umap2)

T2_dim_tsne2 <- my_plotDim(T2, reduction = "tsne", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '2wk T-cells Clusters')
ggsave(filename = "T2_anno_dim_tsne.pdf", plot = T2_dim_tsne2)

# T annotation VlnPlot
T2_Vln2 <- my_violin(T2, features = unlist(T_markers_tmp3), pt.size = 0, mode = 'mtx', group.by = 'MainCluster_anno')
ggsave(filename = "T2_anno_vln.pdf", plot = T2_Vln2)

# T annotation DotPlot
T2_Dot2 <- my_DotPlot_split(T2, features = T_markers_tmp3, group.by = 'MainCluster_anno') + RotatedAxis()
ggsave(filename = "T2_anno_Dot.pdf", plot = T2_Dot2, width = 16, height = 7)



