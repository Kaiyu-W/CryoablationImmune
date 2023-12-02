source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/211111/grouped_analysis//")
load("Cryo_T_grouped.rda")

# 2wk
T2_dim_umap <- my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '2wk T-cells Clusters')
T2_dim_umap2 <- my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '2wk T-cells Clusters')

# Klrb1c Fcer1g
my_DotPlot_split(T2, features = c("Klrb1c","Fcer1g")) + RotatedAxis()
my_DotPlot_split(T2, features = c("Klrb1c","Fcer1g"), group.by='MainCluster_anno') + RotatedAxis()

FeaturePlot(T2, features = c("Klrb1c","Fcer1g"), blend = T)
FeaturePlot(T2, features = c("Klrb1c","Fcer1g","Mki67","Pdcd1",unlist(my_MarkersList[3:4])), 
            cells = grep("^NKT$", T2$MainCluster_anno))

my_plotDim(T2, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', 
           title = '2wk T-cells Clusters', cells = grep("^NKT$", T2$MainCluster_anno))

anno_df <- T2@meta.data[,'MainCluster_anno',drop=F]
anno_df<- anno_df[order(anno_df$MainCluster_anno),,drop = F]

pheatmap::pheatmap(mtx, cluster_cols = F, cluster_rows = F, annotation_col = anno_df, show_colnames = F)

my_pHeatmap <- function(
    seurat_obj, 
    group.by, 
    genes,
    slot = c('data', 'scale.data', 'counts'), 
    default.assay = c("active.assay", names(seurat_obj@assays)), 
    cell = NULL, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    ...
) {
    match.arg(default.assay)
    match.arg(slot)
    default.assay <- default.assay[1]
    slot <- slot[1]
    if (default.assay == "active.assay")
        default.assay <- DefaultAssay(seurat_obj)
    
    data <- methods::slot(seurat_obj@assays[[default.assay]], slot)
    
    anno_df <- seurat_obj@meta.data[, group.by, drop = F]
    anno_df<- anno_df[order(anno_df[, 1, drop = T]), , drop = F]
    genes_in <- genes[genes %in% rownames(data)]
    mtx <- data[genes_in, rownames(anno_df)]
    
    pheatmap::pheatmap(
        mtx, 
        color = color,
        cluster_cols = cluster_cols, 
        cluster_rows = cluster_rows, 
        annotation_col = anno_df, 
        show_colnames = show_colnames,
        ...
    )
}
my_pHeatmap(T2, 
            c('MainCluster_anno', 'MainCluster_C'), 
            c("Klrb1c","Fcer1g","Mki67","Pdcd1",unlist(my_MarkersList[3:4])), 
            color = colorRampPalette(colors = c('Blue',"grey",'Red'))(100)
)

# 1wk
T1_dim_umap <- my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_C', title = '2wk T-cells Clusters')
T1_dim_umap2 <- my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', title = '2wk T-cells Clusters')

# Klrb1c Fcer1g
my_DotPlot_split(T1, features = c("Klrb1c","Fcer1g")) + RotatedAxis()
my_DotPlot_split(T1, features = c("Klrb1c","Fcer1g"), group.by='MainCluster_anno') + RotatedAxis()

FeaturePlot(T1, features = c("Klrb1c","Fcer1g"), blend = T)
FeaturePlot(T1, features = c("Klrb1c","Fcer1g","Mki67","Pdcd1",unlist(my_MarkersList[3:4])), 
            cells = grep("^NKT$", T1$MainCluster_anno))

my_plotDim(T1, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'MainCluster_anno', 
           title = '2wk T-cells Clusters', cells = grep("^NKT$", T1$MainCluster_anno))

anno_df <- T1@meta.data[,'MainCluster_anno',drop=F]
anno_df<- anno_df[order(anno_df$MainCluster_anno),,drop = F]
mtx <- T1@assays$RNA@data[c("Klrb1c","Fcer1g","Mki67", "Pdcd1"), rownames(anno_df)]
pheatmap::pheatmap(mtx, cluster_cols = F, cluster_rows = F, annotation_col = anno_df, show_colnames = F)

# count

Count_T2 <- rbind(my_CountCluster(T2, 'MainCluster_anno')[c("NKT","sum"),1:2],'CD8'= my_CountCluster(T2, 'CD48')['CD8',1:2])
Count_T1 <- rbind(my_CountCluster(T1, 'MainCluster_anno')[c("NKT","sum"),1:2],'CD8'= my_CountCluster(T1, 'MainCluster')['CD8',1:2])

Countperc_T2_all <- apply(Count_T2, 2, function(x) x/x[2] * 100)
Countperc_T1_all <- apply(Count_T1, 2, function(x) x/x[2] * 100)

Countperc_T2_CD8 <- apply(Count_T2, 2, function(x) x/x[3] * 100)
Countperc_T1_CD8 <- apply(Count_T1, 2, function(x) x/x[3] * 100)

# # combined
# T_dim_umap <- my_plotDim(T_sub_SCT, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster', title = 'T Cluster')
# T_dim_umap2 <- my_plotDim(T_sub_SCT, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster2', title = 'T Cluster')
# 
# # Klrb1c Fcer1g
# my_DotPlot_split(T_sub_SCT, features = c("Klrb1c","Fcer1g")) + RotatedAxis()
# my_DotPlot_split(T_sub_SCT, features = c("Klrb1c","Fcer1g"), group.by='T_sub_MainCluster') + RotatedAxis()
# 
# FeaturePlot(T_sub_SCT, features = c("Klrb1c","Fcer1g"), blend = T)
# FeaturePlot(T_sub_SCT, features = c("Klrb1c","Fcer1g","Mki67","Pdcd1"), cells = grep("^CD8_Tea", T_sub_SCT$T_sub_MainCluster))
# FeaturePlot(T_sub_SCT, features = unlist(my_MarkersList[3:4]), cells = grep("^CD8_Tea", T_sub_SCT$T_sub_MainCluster))
# 
# my_plotDim(T_sub_SCT, reduction = "umap", label = F, pt.size = 0.1, label.size = 10, group.by = 'T_sub_MainCluster', 
#            title = '2wk T-cells Clusters', cells = grep("^CD8_Tea", T_sub_SCT$T_sub_MainCluster))
