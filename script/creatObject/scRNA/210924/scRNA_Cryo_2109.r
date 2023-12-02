library(Seurat)
source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")

Cryo_data <- Read10X_h5("Cryo/filtered_feature_bc_matrix.h5")
NonCryo_data <- Read10X_h5("NonCryo/filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("Cryo_", colnames(Cryo_data))
colnames(NonCryo_data) <- paste0("NonCryo_", colnames(NonCryo_data))
Cryo <- CreateSeuratObject(counts = Cryo_data,
                           project = "Cryo",
                           min.cells = 3,
                           min.features = 200)
NonCryo <- CreateSeuratObject(counts = NonCryo_data,
                              project = "NonCryo",
                              min.cells = 3,
                              min.features = 200)
Cryo_merge_raw <- merge(Cryo, NonCryo)
Cryo_merge_raw[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge_raw, pattern = "^mt-")

# visualization
#######################
vln_cryo_m <- my_plotQC(Cryo_merge_raw)
FC_cryo_m <- my_plotFC(Cryo_merge_raw)

hist(Cryo_merge_raw$nFeature_RNA, breaks = 500, xaxt = 'n')
axis(1, seq(0,10000,50))
hist(Cryo_merge_raw$nCount_RNA, breaks = 500, xaxt = 'n')
axis(1, seq(0,200000,1000))
hist(Cryo_merge_raw$percent.mt, breaks = 500, xaxt = 'n')
axis(1, seq(0,100,1))
#######################

Cryo_merge <- subset(Cryo_merge_raw, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 17)
Cryo_merge

Cryo_merge <- SCTransform(Cryo_merge)
# Cryo_merge <- NormalizeData(Cryo_merge)
# Cryo_merge <- ScaleData(Cryo_merge)

Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = 2000)
Cryo_merge <- RunPCA(Cryo_merge)

Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.01)

Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10)
Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10)

# visualization
#######################
my_plotDim(Cryo_merge, reduction = "pca", label = T, pt.size = 1, label.size = 5, title = "Cluster")
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "Cluster")
my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = "Cluster")

my_plotDim(Cryo_merge, reduction = "pca", label = T, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
#######################

Cryo_merge  <- FindClusters(Cryo_merge, resolution = seq(1,10,1)/100)
clustree(Cryo_merge)

Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.01)
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "R=0.01")

Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.4)
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "R=0.4")

DotPlot(Cryo_merge, features = my_MarkersList[[8]]) + RotatedAxis()
DotPlot(Cryo_merge, features = my_MarkersList[-8]) + RotatedAxis()

# R=0.4 cluster 8 has chaos cell-type markers
Cluster8Rdot4_markers <- FindMarkers(Cryo_merge, ident.1 = 8, group.by = 'SCT_snn_res.0.4')
Cluster8Rdot4_markers_list <- my_Markers2df_1Cluster(Cluster8Rdot4_markers, ntop = 200, positive = T, logFC_threshold = 0.25)
my_GO(Cluster8Rdot4_markers_list, ont = 'BP', show = 30, Simplify = T)
Cluster8Rdot4_markers_list

FeaturePlot(Cryo_merge, features = my_MarkersList$T, cells = which(Cryo_merge$SCT_snn_res.0.4=='8'))
FeaturePlot(Cryo_merge, features = my_MarkersList$B, cells = which(Cryo_merge$SCT_snn_res.0.4=='8'))
FeaturePlot(Cryo_merge, features = my_MarkersList$Myeloid, cells = which(Cryo_merge$SCT_snn_res.0.4=='8'))
FeaturePlot(Cryo_merge, features = Cluster8Rdot4_markers_list[1:12])
# cluster 8 should be Myeloid

# choose result of R=0.01
Cryo_merge$MainCluster <- Cryo_merge$SCT_snn_res.0.01
levels(Cryo_merge$MainCluster) <- c('T/NK', 'Myeloid', 'B', 'Myeloid')
Idents(Cryo_merge) <- 'MainCluster'
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "MainCluster")

DotPlot(Cryo_merge, features = my_MarkersList[[8]]) + RotatedAxis()
custom_markers_main <- list(T = 'Cd3e', B = c('Cd19', 'Cd79a'), NK = c('Klrk1', 'Klrb1c'), Myeloid = c('Itgam', 'Cd14'))
my_DotPlot_split(Cryo_merge, features = custom_markers_main, group.by = 'SCT_snn_res.0.01', split.by = 'orig.ident') + RotatedAxis()
my_violin(Cryo_merge, features = unlist(custom_markers_main), pt.size = 0, group.by = 'SCT_snn_res.0.01', split.by = 'orig.ident', mode = 'mtx')

########################
# subcluster for T/NK:

my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "MainCluster")
TNK <- subset(Cryo_merge, ident = 'T/NK')
DefaultAssay(TNK) <- 'RNA'
TNK <- SCTransform(TNK)
TNK <- FindVariableFeatures(TNK, selection.method = "vst", nfeatures = 3000)
TNK <- RunPCA(TNK)
TNK <- RunUMAP(TNK, dims = 1:10)
TNK <- RunTSNE(TNK, dims = 1:10)
TNK <- FindNeighbors(TNK, dims = 1:12)
TNK <- FindClusters(TNK, resolution = 0.01)

# cluster
TNK  <- FindClusters(TNK, resolution = seq(1,10,1)/100)
TNK  <- FindClusters(TNK, resolution = seq(1,5,1)/10)
clustree(TNK)

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'TNK Cryo_NonCryo')
my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(paste0("0", 1:9), 1:3)))
Idents(TNK) <- 'SCT_snn_res.0.3'

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'TNK R=0.3')
my_DotPlot_split(TNK, features = my_MarkersList[3:4]) + RotatedAxis()
my_violin(TNK, features = unlist(my_MarkersList[3:4]), pt.size = 0, split.by = 'orig.ident', mode = 'mtx')
FeaturePlot(TNK, features = unlist(my_MarkersList[3:4]))

# NKT? No
FeaturePlot(TNK, features = NKT_marker)
my_DotPlot_split(TNK, features = NKT_marker) + RotatedAxis()
my_violin(TNK, features = NKT_marker, pt.size = 0, split.by = 'orig.ident', mode = 'mtx')

# Cluster 0 chaos?
TNK <- FindSubCluster(TNK, cluster = '0', graph.name = 'SCT_snn', subcluster.name = 'New0_R0.2', resolution = 0.2)
my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New0_R0.2")
my_DotPlot_split(TNK, features = my_MarkersList[3], group.by = "New0_R0.2") + RotatedAxis()
my_violin(TNK, features = unlist(my_MarkersList[3]), pt.size = 0, group.by = "New0_R0.2", mode = 'mtx')

# Cluster 7 chaos?
Idents(TNK) <- 'New0_R0.2'
mtx_c7 <- as.matrix(TNK@assays$SCT@data)[unlist(my_MarkersList[3:4])[-c(2,8)], grep("^7", TNK$New0_R0.2)]
res_heat <- pheatmap::pheatmap(
    mtx_c7, 
    show_colnames = F, 
    cluster_rows = F)
col_cluster <- cutree(res_heat$tree_col, k = 3)
pheatmap::pheatmap(
    rbind(mtx_c7, col_cluster), 
    show_colnames = F, 
    cluster_rows = F)
# 1 is CD8, 2 is CD4, 3 is NK
col_cluster_res <- sapply(col_cluster, function(x) switch(
    x,
    '7_8',
    '7_4',
    '7_nk'
))
TNK$New7 <- TNK$New0_R0.2
TNK$New7[names(col_cluster_res)] <- col_cluster_res

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New7")
Seurat::DoHeatmap(TNK, group.by = 'New7', features = unlist(my_MarkersList[3:4]), cells = grep("^7", TNK$New7))
my_violin(TNK, features = unlist(my_MarkersList[3:4]), pt.size = 0, group.by = "New7", mode = 'mtx')

# Cluster 3 chaos?
Idents(TNK) <- 'New0_R0.2'
mtx_c3 <- as.matrix(TNK@assays$SCT@data)[my_MarkersList[[3]][-c(1:2)], grep("^3", TNK$New0_R0.2)]
celltype_3 <- apply(mtx_c3, 2, function(x) {
    if (x[1] >= max(x[2:3]) | all(x[2:3] == 0)) 
        '3_4'
    else
        '3_8'
})
temp <- sapply(celltype_3, function(x) grep(x, c('3_4', '3_8')))
pheatmap::pheatmap(
    rbind(mtx_c3, temp), 
    show_colnames = F, 
    cluster_rows = F)

TNK$New3 <- TNK$New7
TNK$New3[names(celltype_3)] <- celltype_3

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New3")
Seurat::DoHeatmap(TNK, group.by = 'New3', features = unlist(my_MarkersList[3]), cells = grep("^3", TNK$New3))
my_violin(TNK, features = unlist(my_MarkersList[3]), pt.size = 0, group.by = "New3", mode = 'mtx')

# Cluster 8 chaos?
Cluster8Rdot3_markers <- FindMarkers(TNK, ident.1 = 8, group.by = 'SCT_snn_res.0.3')
Cluster8Rdot3_markers_list <- my_Markers2df_1Cluster(Cluster8Rdot3_markers, ntop = 50, positive = T, logFC_threshold = 0.25)
my_GO(Cluster8Rdot3_markers_list, ont = 'BP', show = 30, Simplify = T)
Cluster8Rdot3_markers_list

# cell type identify
TNK$TNK_MainCluster <- as.factor(TNK$New3)
TNK_MainCluster_df <- rbind(
    data.frame(x='NK',y=c(5,6,'7_nk')),
    data.frame(x='T_CD4',y=c('0_0',2,'3_4',8,'7_4')),
    data.frame(x='T_CD8',y=c('0_1',1,4,9,'3_8','7_8'))
)
TNK_MainCluster_df <- TNK_MainCluster_df[order(TNK_MainCluster_df$y, decreasing = F),]
if(all(TNK_MainCluster_df$y == levels(TNK$TNK_MainCluster)))
    levels(TNK$TNK_MainCluster) <- TNK_MainCluster_df$x
Idents(TNK) <- 'TNK_MainCluster'
my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "TNK_MainCluster")

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New3")
custom_markers_TNKmain <- list(T = 'Cd3e', NK = c('Klrk1', 'Klrb1c'), T_CD4 = 'Cd4', T_CD8 = c('Cd8a', 'Cd8b1'))
my_DotPlot_split(TNK, features = custom_markers_TNKmain, group.by = 'New3') + RotatedAxis()
my_violin(TNK, features = unlist(custom_markers_TNKmain), pt.size = 0, group.by = 'New3', mode = 'mtx')

# count
my_CountCluster(Cryo_merge, 'MainCluster')
my_CountCluster(TNK, trend_factor = 4607/4751)

# # Enrichment
# Markers_merge_list9 <- FindAllMarkers(TNK, only.pos = T)
# Markers_merge_df9 <- my_Markers2df_multiple(Markers_merge_list9, logFC_threshold = 0, positive = T, n_top = 100)
# temp_BP_Cryo <- my_GO(Markers_merge_df9$Cluster_T_CD4, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
# my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
# temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
# temp_BP_Cryo_df$Description

# subcluster for T_CD8
T_CD8 <- subset(TNK, ident = 'T_CD8')
DefaultAssay(T_CD8) <- 'RNA'
T_CD8 <- SCTransform(T_CD8)
T_CD8 <- FindVariableFeatures(T_CD8, selection.method = "vst", nfeatures = 3000)
T_CD8 <- RunPCA(T_CD8)
T_CD8 <- RunUMAP(T_CD8, dims = 1:10)
T_CD8 <- RunTSNE(T_CD8, dims = 1:10)
T_CD8 <- FindNeighbors(T_CD8, dims = 1:12)

T_CD8  <- FindClusters(T_CD8, resolution = seq(1,20)/10)
clustree(T_CD8)

my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'T_CD8 Cryo_NonCryo')
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(1:9)))

Idents(T_CD8) <- 'SCT_snn_res.0.7'

my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_CD8 R=0.7')
my_DotPlot_split(T_CD8, features = T_marker$CD8) + RotatedAxis()
my_violin(T_CD8, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx')
FeaturePlot(T_CD8, features = 'Mki67')

T_CD8$T_CD8_MainCluster <- as.factor(T_CD8$SCT_snn_res.0.7)
T_CD8_MainCluster_df <- rbind(
    data.frame(x='CD8_Central_Memory',y=c(6)),
    data.frame(x='CD8_Effector_Memory',y=c(1,2)),
    data.frame(x='CD8_Effector',y=c(3,4,5,7,8)),
    data.frame(x='CD8_Naive',y=c(0))
)
T_CD8_MainCluster_df <- T_CD8_MainCluster_df[order(T_CD8_MainCluster_df$y, decreasing = F),]
if (all(T_CD8_MainCluster_df$y == levels(T_CD8$T_CD8_MainCluster)))
    levels(T_CD8$T_CD8_MainCluster) <- T_CD8_MainCluster_df$x

T_CD8$T_CD8_ExCluster <- as.factor(T_CD8$SCT_snn_res.0.7)
T_CD8_ExCluster_df <- rbind(
    data.frame(x='CD8_Central_Memory',y=c(6)),
    data.frame(x='CD8_Effector_Memory',y=c(1,2)),
    data.frame(x='CD8_Effector',y=c(4)),
    data.frame(x='CD8_Effector_Exhausted',y=c(3,5,7)),
    data.frame(x='CD8_Effector_Proliferating',y=c(8)),
    data.frame(x='CD8_Naive',y=c(0))
)
T_CD8_ExCluster_df <- T_CD8_ExCluster_df[order(T_CD8_ExCluster_df$y, decreasing = F),]
if (all(T_CD8_ExCluster_df$y == levels(T_CD8$T_CD8_ExCluster)))
    levels(T_CD8$T_CD8_ExCluster) <- T_CD8_ExCluster_df$x

my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD8_MainCluster")
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD8_ExCluster")

my_CountCluster(T_CD8, 'T_CD8_MainCluster', trend_factor = 4607/4751)
my_CountCluster(T_CD8, 'T_CD8_ExCluster', trend_factor = 4607/4751)

Idents(TNK) <- 'TNK_MainCluster'
TNK <- my_AddMeta(TNK, new_ident = T_CD8$T_CD8_MainCluster, Replace = T)
TNK <- my_AddMeta(TNK, new_ident = T_CD8$T_CD8_ExCluster, Replace = T)
# my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD8_T_CD8_ExCluster")


# subcluster for T_CD4
T_CD4 <- subset(TNK, ident = 'T_CD4')
DefaultAssay(T_CD4) <- 'RNA'
T_CD4 <- SCTransform(T_CD4)
T_CD4 <- FindVariableFeatures(T_CD4, selection.method = "vst", nfeatures = 3000)
T_CD4 <- RunPCA(T_CD4)
T_CD4 <- RunUMAP(T_CD4, dims = 1:10)
T_CD4 <- RunTSNE(T_CD4, dims = 1:10)
T_CD4 <- FindNeighbors(T_CD4, dims = 1:12)

T_CD4  <- FindClusters(T_CD4, resolution = seq(1,20)/10)
clustree(T_CD4)

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'T_CD4 Cryo_NonCryo')
my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(1:9)))

Idents(T_CD4) <- 'SCT_snn_res.0.5'

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_CD4 R=0.5')
my_DotPlot_split(T_CD4, features = T_marker$CD4) + RotatedAxis()
my_violin(T_CD4, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx')
FeaturePlot(T_CD4, features = 'Mki67')

T_CD4$T_CD4_MainCluster <- as.factor(T_CD4$SCT_snn_res.0.5)
T_CD4_MainCluster_df <- rbind(
    data.frame(x='CD4_Naive',y=c(0)),
    data.frame(x='CD4_Treg',y=c(3)),
    data.frame(x='CD4_Th1',y=c(1,4)),
    data.frame(x='CD4_Th2',y=c(6)),
    data.frame(x='CD4_Th17',y=c(5)),
    data.frame(x='CD4_Tfh',y=c(2))
)
T_CD4_MainCluster_df <- T_CD4_MainCluster_df[order(T_CD4_MainCluster_df$y, decreasing = F),]
if (all(T_CD4_MainCluster_df$y == levels(T_CD4$T_CD4_MainCluster)))
    levels(T_CD4$T_CD4_MainCluster) <- T_CD4_MainCluster_df$x

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD4_MainCluster")

my_CountCluster(T_CD4, 'T_CD4_MainCluster', trend_factor = 4607/4751)

Idents(TNK) <- 'T_CD8_T_CD8_ExCluster'
TNK <- my_AddMeta(TNK, new_ident = T_CD4$T_CD4_MainCluster, Replace = T)

my_plotDim(TNK, reduction = "umap", label = T, pt.size = 2, label.size = 10, group.by = "T_CD4_T_CD4_MainCluster", title = "T MainCluster")

my_CountCluster(TNK, 'T_CD4_T_CD4_MainCluster', trend_factor = 4607/4751, Heatmap = T)

TNK$FineCluster <- TNK$T_CD4_T_CD4_MainCluster
Idents(TNK) <- 'FineCluster'

# Enrichment
Markers_TNK_list <- FindAllMarkers(TNK, only.pos = T)
Markers_TNK_df <- my_Markers2df_multiple(Markers_TNK_list, logFC_threshold = 0.25, positive = T, n_top = 50)

temp_GO <- function(genelist) {
    temp_BP <- my_GO(genelist, return_plot = F, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
    plot <- my_GO(temp_BP, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
    temp_BP_df <- temp_BP@result[grep("interferon", temp_BP@result$Description), ]
    if (length(temp_BP_df$Description) == 0)
        return(plot)
    ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
    ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
    ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
    index <- which(genelist%in%ifn_gene)
    
    res <- list(pathway = temp_BP_df$Description, gene_rank = index)
    plot
    return(res)
}
# temp_GO(Markers_TNK_df$Cluster_CD4_Naive)
# temp_GO(Markers_TNK_df$Cluster_CD8_Naive)
temp_GO(Markers_TNK_df$Cluster_CD8_Effector_Memory)
temp_GO(Markers_TNK_df$Cluster_CD8_Central_Memory)


########################
# library(MultiK)
# multik <- MultiK(T_CD8, reps=100, resolution = seq(0.05, 2, 0.05))
# DiagMultiKPlot(multik$k, multik$consensus) #7 13
# clusters <- getClusters(seu = T_CD8, optK = 13)
# pval <- CalcSigClust(T_CD8, clusters$clusters)
# PlotSigClust(T_CD8, clusters$clusters, pval)
########################
########################after changing strategy
#CD8:
# library(MultiK)
# multik <- MultiK(T_CD8, reps=100, resolution = seq(0.05, 2, 0.05))
# DiagMultiKPlot(multik$k, multik$consensus) #7 13
clustree(T_CD8)
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "orig.ident")
# 7 : 0.6
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "SCT_snn_res.0.6")
# enrich clusters: 0/1/4/5
my_CountCluster(T_CD8, 'SCT_snn_res.0.6', trend_factor = 4607/4751)

# 13: 1.6
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "SCT_snn_res.1.6")
# enrich clusters: 1/3/4/6/7/8/9/11
my_CountCluster(T_CD8, 'SCT_snn_res.1.6', trend_factor = 4607/4751)

## Enrichment
# 0.6
Idents(T_CD8) <- 'SCT_snn_res.0.6'
Markers_T_CD8_0.6_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_0.6_df <- my_Markers2df_multiple(Markers_T_CD8_0.6_list, logFC_threshold = 0.25, positive = T, n_top = 100)
temp_GO(Markers_T_CD8_0.6_df$Cluster_0) # alpha
temp_GO(Markers_T_CD8_0.6_df$Cluster_1)
temp_GO(Markers_T_CD8_0.6_df$Cluster_4) # significant ! alpha/beta
temp_GO(Markers_T_CD8_0.6_df$Cluster_5)
# 1.6
Idents(T_CD8) <- 'SCT_snn_res.1.6'
Markers_T_CD8_1.6_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_1.6_df <- my_Markers2df_multiple(Markers_T_CD8_1.6_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD8_1.6_df$Cluster_1)
temp_GO(Markers_T_CD8_1.6_df$Cluster_3)
temp_GO(Markers_T_CD8_1.6_df$Cluster_4)
temp_GO(Markers_T_CD8_1.6_df$Cluster_6) # significant !
temp_GO(Markers_T_CD8_1.6_df$Cluster_7) # alpha
temp_GO(Markers_T_CD8_1.6_df$Cluster_8)
temp_GO(Markers_T_CD8_1.6_df$Cluster_9)
temp_GO(Markers_T_CD8_1.6_df$Cluster_11)
# orig.ident
Idents(T_CD8) <- 'orig.ident'
Markers_T_CD8_0_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_0_df <- my_Markers2df_multiple(Markers_T_CD8_0_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD8_0_df$Cluster_Cryo)

#CD4:
# library(MultiK)
# multik <- MultiK(T_CD4, reps=100, resolution = seq(0.05, 2, 0.05))
# DiagMultiKPlot(multik$k, multik$consensus) #9 4 3
clustree(T_CD4)
my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "orig.ident")
# 9 : 1.2-1.4
my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "SCT_snn_res.1.2")
# enrich clusters: 0/1/10
my_CountCluster(T_CD4, 'SCT_snn_res.1.2', trend_factor = 4607/4751)

## Enrichment
# 0.6
Idents(T_CD4) <- 'SCT_snn_res.1.2'
Markers_T_CD4_1.2_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_1.2_df <- my_Markers2df_multiple(Markers_T_CD4_1.2_list, logFC_threshold = 0.25, positive = T, n_top = 100)
temp_GO(Markers_T_CD4_1.2_df$Cluster_0)
temp_GO(Markers_T_CD4_1.2_df$Cluster_1)
temp_GO(Markers_T_CD4_1.2_df$Cluster_10) # significant ! alpha/beta
# orig.ident
Idents(T_CD4) <- 'orig.ident'
Markers_T_CD4_0_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_0_df <- my_Markers2df_multiple(Markers_T_CD4_0_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD4_0_df$Cluster_Cryo)
########################

save(Cryo_merge, TNK, T_CD4, T_CD8, file = 'Single_analysis.rda')
