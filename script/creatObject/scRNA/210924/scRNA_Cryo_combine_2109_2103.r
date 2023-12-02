
source("E://Cryo-TCR/server/auto/utilities.r")

temp_GO <- function(genelist, return_ntop = NULL, pathway_name = 'interferon') {
    temp_BP <- my_GO(genelist, return_plot = F, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
    plot <- my_GO(temp_BP, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
    temp_BP_df <- temp_BP@result[grep(pathway_name, temp_BP@result$Description), ]
    index_1 <- grep(pathway_name, temp_BP@result$Description)
    if (length(temp_BP_df$Description) == 0) {
        if (is.null(return_ntop)) {
            return(plot) 
        } else {
            plot
            res_top <- temp_BP@result$Description[1:return_ntop]
            res <- list(top_pathway = res_top)
            return(res)
        }
    }
    ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
    ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
    ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
    index_2 <- which(genelist%in%ifn_gene)
    
    res <- list(pathway = temp_BP_df$Description, pathway_rank = index_1, gene = ifn_gene, gene_rank = index_2)
    if (!is.null(return_ntop)) {
        res_top <- temp_BP@result$Description[1:return_ntop]
        res <- append(res, list(top_pathway = res_top))
    }
    plot
    return(res)
}
temp_GO_detail <- function(genelist) {
    temp_BP <- my_GO(genelist, return_plot = F, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
    my_GO(temp_BP, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
    temp_BP_df <- temp_BP@result[grep("interferon", temp_BP@result$Description), ]
    temp_BP_df <- temp_BP_df[-grep("gamma", temp_BP_df$Description), ]
    
    temp_BP_df$geneID <- sapply(temp_BP_df$geneID, function(geneID) {
        ifn_gene_id <- paste0(geneID, collapse = "/")
        ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
        ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
        paste(ifn_gene, collapse = '/')
    }, USE.NAMES = F)
    
    temp_BP_df
}

load(file = "C://celldex_reference/Database_form_celldex.rda") # reference data
###########################################1. batch between first and second one without CD45- cells
###
# Integret by merge function
# or FindIntegrationAnchors with IntegrateData
###
###########################################import raw_data

setwd("D://Ji_Wangxue/")
Cryo_data_1 <- Read10X_h5("Cryo_filtered_feature_bc_matrix.h5")
NonCryo_data_1 <- Read10X_h5("NonCryo_filtered_feature_bc_matrix.h5")
colnames(Cryo_data_1) <- paste0("Old_Cryo_", colnames(Cryo_data_1))
colnames(NonCryo_data_1) <- paste0("Old_NonCryo_", colnames(NonCryo_data_1))
Cryo_1 <- CreateSeuratObject(counts = Cryo_data_1,
                             project = "Old_Cryo", 
                             min.cells = 3, 
                             min.features = 200)
NonCryo_1 <- CreateSeuratObject(counts = NonCryo_data_1, 
                                project = "Old_NonCryo", 
                                min.cells = 3, 
                                min.features = 200)
Cryo_1@meta.data$orig.ident <- Cryo_1@project.name
NonCryo_1@meta.data$orig.ident <- NonCryo_1@project.name
Cryo_merge_1 <- merge(Cryo_1, NonCryo_1)

setwd("E://Cryo-TCR/server/210924/")
Cryo_data_2 <- Read10X_h5("Cryo/filtered_feature_bc_matrix.h5")
NonCryo_data_2 <- Read10X_h5("NonCryo/filtered_feature_bc_matrix.h5")
colnames(Cryo_data_2) <- paste0("New_Cryo_", colnames(Cryo_data_2))
colnames(NonCryo_data_2) <- paste0("New_NonCryo_", colnames(NonCryo_data_2))
Cryo_2 <- CreateSeuratObject(counts = Cryo_data_2,
                             project = "New_Cryo",
                             min.cells = 3,
                             min.features = 200)
NonCryo_2 <- CreateSeuratObject(counts = NonCryo_data_2,
                                project = "New_NonCryo",
                                min.cells = 3,
                                min.features = 200)
Cryo_2@meta.data$orig.ident <- Cryo_2@project.name
NonCryo_2@meta.data$orig.ident <- NonCryo_2@project.name
Cryo_merge_2 <- merge(Cryo_2, NonCryo_2)

Cryo_combined <- merge(Cryo_merge_1, Cryo_merge_2)
Cryo_combined[["percent.mt"]] <- PercentageFeatureSet(Cryo_combined, pattern = "^mt-")
Idents(Cryo_combined) <- 'orig.ident'

Cryo_combined$orig.ident2 <- Cryo_combined$orig.ident
Cryo_combined$orig.ident2 <- sub("^.*_Cryo$", "Cryo", Cryo_combined$orig.ident2)
Cryo_combined$orig.ident2 <- sub("^.*_NonCryo$", "NonCryo", Cryo_combined$orig.ident2)
Cryo_combined$orig.ident3 <- Cryo_combined$orig.ident
Cryo_combined$orig.ident3 <- sub("^New_.*$", "New", Cryo_combined$orig.ident3)
Cryo_combined$orig.ident3 <- sub("^Old_.*$", "Old", Cryo_combined$orig.ident3)

# visualization
#######################
vln_cryo_m <- my_plotQC(Cryo_combined, pt.size = 0)
FC_cryo_m <- my_plotFC(Cryo_combined)

hist(Cryo_combined$nFeature_RNA, breaks = 500, xaxt = 'n')
axis(1, seq(0,10000,50))
hist(Cryo_combined$nCount_RNA, breaks = 500, xaxt = 'n')
axis(1, seq(0,200000,1000))
hist(Cryo_combined$percent.mt, breaks = 500, xaxt = 'n')
axis(1, seq(0,100,1))
#######################
Cryo_merge <- subset(Cryo_combined, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA < 50000 & percent.mt < 17)
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

my_plotDim(Cryo_merge, reduction = "pca", label = F, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
my_plotDim(Cryo_merge, reduction = "umap", label = F, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
#######################
Cryo_merge  <- FindClusters(Cryo_merge, resolution = seq(1,10,1)/100)
clustree(Cryo_merge)

Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.06)
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "R=0.06")

DotPlot(Cryo_merge, features = my_MarkersList[[8]]) + RotatedAxis()
DotPlot(Cryo_merge, features = my_MarkersList[-8]) + RotatedAxis()
custom_markers_Cryomain <- list(T = c('Cd3e', 'Cd247'), NK = c('Klrk1', 'Klrb1c'), B = c('Cd19', 'Cd79a'), Myeloid = c('Itgam', 'Cd14'), Mast = c('Mcpt1', 'Mcpt2'))
my_DotPlot_split(Cryo_merge, features = custom_markers_Cryomain) + RotatedAxis()
my_violin(Cryo_merge, features = unlist(custom_markers_Cryomain), pt.size = 0, mode = 'mtx')

# cell type identify
Cryo_merge$MainCluster <- Cryo_merge$SCT_snn_res.0.06
MainCluster <- rbind(
    data.frame(x='T',y=c(0,4,6,9)),
    data.frame(x='NK',y=c(7)),
    data.frame(x='B',y=c(2)),
    data.frame(x='Myeloid',y=c(1,3,8)),
    data.frame(x='Mast',y=c(5))
)
MainCluster <- MainCluster[order(MainCluster$y, decreasing = F),]
if(all(MainCluster$y == levels(Cryo_merge$MainCluster)))
    levels(Cryo_merge$MainCluster) <- MainCluster$x
Idents(Cryo_merge) <- 'MainCluster'
my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "MainCluster")

# count
my_CountCluster(Cryo_merge, 'MainCluster', 'orig.ident2')
my_CountCluster(Cryo_merge, 'MainCluster', 'orig.ident3')

# enrichment
Idents(Cryo_merge) <- 'orig.ident2'
Markers_MainCluster_list <- FindAllMarkers(Cryo_merge, only.pos = T)
Markers_MainCluster_df <- my_Markers2df_multiple(Markers_MainCluster_list, logFC_threshold = 0.25, positive = T, n_top = 120)
GO_cryo <- temp_GO(Markers_MainCluster_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_MainCluster_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]


########################
# subcluster for T:
Idents(Cryo_merge) <- 'MainCluster'
T_sub <- subset(Cryo_merge, ident = 'T')
DefaultAssay(T_sub) <- 'RNA'
T_sub <- SCTransform(T_sub)
T_sub <- FindVariableFeatures(T_sub, selection.method = "vst", nfeatures = 3000)
T_sub <- RunPCA(T_sub)
T_sub <- RunUMAP(T_sub, dims = 1:10)
T_sub <- RunTSNE(T_sub, dims = 1:10)
T_sub <- FindNeighbors(T_sub, dims = 1:12)
T_sub <- FindClusters(T_sub, resolution = 0.01)

# cluster
T_sub  <- FindClusters(T_sub, resolution = seq(1,10,1)/100)
T_sub  <- FindClusters(T_sub, resolution = seq(1,5,1)/10)
clustree(T_sub)

my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'T_sub Cryo_NonCryo')
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(paste0("0", 1:9), 1:3)))
Idents(T_sub) <- 'SCT_snn_res.0.4'

my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub R=0.4')
my_DotPlot_split(T_sub, features = my_MarkersList[3:4]) + RotatedAxis()
my_violin(T_sub, features = unlist(my_MarkersList[3:4]), pt.size = 0, mode = 'mtx')

# Cluster 6 chaos?
T_sub <- FindSubCluster(T_sub, cluster = '6', graph.name = 'SCT_snn', subcluster.name = 'New6_R0.2', resolution = 0.2)
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New6_R0.2")
my_DotPlot_split(T_sub, features = my_MarkersList[3], group.by = "New6_R0.2") + RotatedAxis()
my_violin(T_sub, features = c("Ptprc", unlist(my_MarkersList[3:4]), "Sell"), pt.size = 0, group.by = "New6_R0.2", mode = 'mtx')

# Cluster 1 chaos?
Idents(T_sub) <- 'New6_R0.2'
mtx_c1 <- as.matrix(T_sub@assays$SCT@data)[my_MarkersList[[3]][-c(1:2)], grep("^1$", T_sub$New6_R0.2)]
celltype_1 <- apply(mtx_c1, 2, function(x) {
    if (x[1] >= max(x[2:3]) | all(x[2:3] == 0)) 
        '1_4'
    else
        '1_8'
})
temp <- sapply(celltype_1, function(x) grep(x, c('1_4', '1_8')))
pheatmap::pheatmap(
    rbind(mtx_c1, temp), 
    show_colnames = F, 
    cluster_rows = F)

T_sub$New1 <- T_sub$New6_R0.2
T_sub$New1[names(celltype_1)] <- celltype_1
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "New1")
my_violin(T_sub, features = c("Ptprc", unlist(my_MarkersList[3:4]), "Sell"), pt.size = 0, group.by = "New1", mode = 'mtx')

# cell type identify
T_sub$T_sub_MainCluster <- as.factor(T_sub$New1)
T_sub_MainCluster_df <- rbind(
    data.frame(x='Non_Immune',y=c(12)),
    data.frame(x='NKT',y=c(5,11)),
    data.frame(x='T_CD4',y=c(3,4,7,'6_0','6_2','1_4')),
    data.frame(x='T_CD8',y=c(0,2,8,9,'6_1','1_8')),
    data.frame(x='T',y=c(10, '6_3'))
)
T_sub_MainCluster_df <- T_sub_MainCluster_df[order(T_sub_MainCluster_df$y, decreasing = F),]
if(all(T_sub_MainCluster_df$y == levels(T_sub$T_sub_MainCluster)))
    levels(T_sub$T_sub_MainCluster) <- T_sub_MainCluster_df$x
Idents(T_sub) <- 'T_sub_MainCluster'
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "T_sub_MainCluster")

# count
my_CountCluster(Cryo_merge, 'MainCluster', 'orig.ident2')
my_CountCluster(T_sub, 'T_sub_MainCluster', 'orig.ident2', trend_factor = 8649/9021)

# enrichment
# T_CD4/8 cannot enrich for interferon-a/b
Idents(T_sub) <- 'T_sub_MainCluster'
Markers_T_sub_list <- FindAllMarkers(T_sub, only.pos = T)
Markers_T_sub_df <- my_Markers2df_multiple(Markers_T_sub_list, logFC_threshold = 0.25, positive = T, n_top = 100)
temp_GO(Markers_T_sub_df$Cluster_T_CD8)
temp_GO(Markers_T_sub_df$Cluster_T_CD4)

Idents(T_sub) <- 'orig.ident2'
Markers_T_sub_list2 <- FindAllMarkers(T_sub, only.pos = T)
Markers_T_sub_df1 <- my_Markers2df_multiple(Markers_T_sub_list2, logFC_threshold = 0, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_sub_df1$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_sub_df1$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]

########################
# save(Cryo_merge, Cryo_combined, T_sub, file = 'Combined_analysis.rda')
load(file = 'Combined_analysis.rda')

# subcluster for T_CD8:
Idents(T_sub) <- 'T_sub_MainCluster'
T_CD8 <- subset(T_sub, ident = 'T_CD8')
DefaultAssay(T_CD8) <- 'RNA'
T_CD8 <- SCTransform(T_CD8)
T_CD8 <- FindVariableFeatures(T_CD8, selection.method = "vst", nfeatures = 3000)
T_CD8 <- RunPCA(T_CD8)
T_CD8 <- RunUMAP(T_CD8, dims = 1:10)
T_CD8 <- RunTSNE(T_CD8, dims = 1:10)
T_CD8 <- FindNeighbors(T_CD8, dims = 1:12)

T_CD8  <- FindClusters(T_CD8, resolution = seq(1,20)/10)
clustree(T_CD8)

# library(MultiK)
# multik_CD8 <- MultiK(T_CD8, reps=100, resolution = seq(0.05, 2, 0.05))
# DiagMultiKPlot(multik_CD8$k, multik_CD8$consensus) #15

my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'T_CD8 Cryo_NonCryo')
my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(1:9)))

Idents(T_CD8) <- 'SCT_snn_res.0.9'

my_plotDim(T_CD8, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_CD8 R=0.9')
my_DotPlot_split(T_CD8, features = T_marker$CD8) + RotatedAxis()
my_violin(T_CD8, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx')
FeaturePlot(T_CD8, features = 'Mki67')
my_CountCluster(T_CD8, group2 = 'orig.ident2', trend_factor = 8649/9021)

T_CD8$T_CD8_MainCluster <- as.factor(T_CD8$SCT_snn_res.0.9)
T_CD8_MainCluster_df <- rbind(
    data.frame(x='CD8_Central_Memory',y=c(10)),
    data.frame(x='CD8_Effector_Memory',y=c(0,3,6,7,8)),
    data.frame(x='CD8_Effector',y=c(1,4,5,9)),
    data.frame(x='CD8_Naive',y=c(2))
)
T_CD8_MainCluster_df <- T_CD8_MainCluster_df[order(T_CD8_MainCluster_df$y, decreasing = F),]
if (all(T_CD8_MainCluster_df$y == levels(T_CD8$T_CD8_MainCluster)))
    levels(T_CD8$T_CD8_MainCluster) <- T_CD8_MainCluster_df$x

T_CD8$T_CD8_ExCluster <- as.factor(T_CD8$SCT_snn_res.0.9)
T_CD8_ExCluster_df <- rbind(
    data.frame(x='CD8_Central_Memory',y=c(10)),
    data.frame(x='CD8_Effector_Memory',y=c(0,3,6,7)),
    data.frame(x='CD8_Effector_Exhausted',y=c(1,4,9)),
    data.frame(x='CD8_Effector_Memory_Exhausted',y=c(8)),
    data.frame(x='CD8_Effector_Proliferating',y=c(5)),
    data.frame(x='CD8_Naive',y=c(2))
)
T_CD8_ExCluster_df <- T_CD8_ExCluster_df[order(T_CD8_ExCluster_df$y, decreasing = F),]
if (all(T_CD8_ExCluster_df$y == levels(T_CD8$T_CD8_ExCluster)))
    levels(T_CD8$T_CD8_ExCluster) <- T_CD8_ExCluster_df$x

my_plotDim(T_CD8, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "T_CD8_MainCluster")
my_plotDim(T_CD8, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "T_CD8_ExCluster")

my_CountCluster(T_CD8, 'T_CD8_MainCluster', 'orig.ident2', trend_factor = 8649/9021)
my_CountCluster(T_CD8, 'T_CD8_ExCluster', 'orig.ident2', trend_factor = 8649/9021)

# Idents(TNK) <- 'TNK_MainCluster'
# TNK <- my_AddMeta(TNK, new_ident = T_CD8$T_CD8_MainCluster, Replace = T)
# TNK <- my_AddMeta(TNK, new_ident = T_CD8$T_CD8_ExCluster, Replace = T)
# # my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD8_T_CD8_ExCluster")

# enrichment
Idents(T_CD8) <- 'T_CD8_MainCluster'
Markers_T_CD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_df <- my_Markers2df_multiple(Markers_T_CD8_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector_Memory)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Naive)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Central_Memory)

Idents(T_CD8) <- 'T_CD8_ExCluster'
Markers_T_CD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_df <- my_Markers2df_multiple(Markers_T_CD8_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector_Memory)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector_Memory_Exhausted)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector_Exhausted)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Effector_Proliferating)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Naive)
temp_GO(Markers_T_CD8_df$Cluster_CD8_Central_Memory)

Idents(T_CD8) <- 'SCT_snn_res.0.9'
Markers_T_CD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_df <- my_Markers2df_multiple(Markers_T_CD8_list, logFC_threshold = 0.25, positive = T, n_top = 200)
for (i in 0:10) {
    cat('Cluster_', i, "\n", sep = "")
    print(temp_GO(Markers_T_CD8_df[, paste0('Cluster_', i)]))
}
my_plotDim(T_CD8, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "orig.ident2")


Cluster_7_GO_df <- temp_GO_detail(Markers_T_CD8_df$Cluster_7)
Cluster_8_GO_df <- temp_GO_detail(Markers_T_CD8_df$Cluster_8)
write.csv(Cluster_7_GO_df, file = 'T_CD8_Cluster_7_GO_df.csv')
write.csv(Cluster_8_GO_df, file = 'T_CD8_Cluster_8_GO_df.csv')

# my_DotPlot_split(T_CD8, features = T_marker$CD8, idents = c('7', '8')) + RotatedAxis()
# my_violin(T_CD8, features = unlist(T_marker$CD8), pt.size = 0, idents = c('7', '8'), mode = 'mtx')


# subcluster for T_CD4:
Idents(T_sub) <- 'T_sub_MainCluster'
T_CD4 <- subset(T_sub, ident = 'T_CD4')
DefaultAssay(T_CD4) <- 'RNA'
T_CD4 <- SCTransform(T_CD4)
T_CD4 <- FindVariableFeatures(T_CD4, selection.method = "vst", nfeatures = 3000)
T_CD4 <- RunPCA(T_CD4)
T_CD4 <- RunUMAP(T_CD4, dims = 1:10)
T_CD4 <- RunTSNE(T_CD4, dims = 1:10)
T_CD4 <- FindNeighbors(T_CD4, dims = 1:12)

T_CD4  <- FindClusters(T_CD4, resolution = seq(1,20)/10)
clustree(T_CD4)

# library(MultiK)
# multik_CD8 <- MultiK(T_CD4, reps=100, resolution = seq(0.05, 2, 0.05))
# DiagMultiKPlot(multik_CD8$k, multik_CD8$consensus) #15

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident', title = 'T_CD4 Cryo_NonCryo')
my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(1:9)))

Idents(T_CD4) <- 'SCT_snn_res.0.5'

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_CD4 R=0.5')
my_DotPlot_split(T_CD4, features = T_marker$CD4) + RotatedAxis()
my_violin(T_CD4, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx')
FeaturePlot(T_CD4, features = 'Mki67')
my_CountCluster(T_CD4, group2 = 'orig.ident2', trend_factor = 8649/9021)

T_CD4_temp <- my_ReferenceCluster(T_CD4, reference = DB_Mi, label_fine = T, celltype_list = c("CD4+ T cells"))
my_plotDim(T_CD4_temp, pt.size = 1, reduction = 'umap', group.by = 'DB_Mi')

my_DotPlot_split(T_CD4, features = T_marker$CD4, idents = c('0', '4')) + RotatedAxis()
my_violin(T_CD4, features = unlist(T_marker$CD4), pt.size = 0, idents = c('0', '4'), mode = 'mtx')

T_CD4$T_CD4_MainCluster <- as.factor(T_CD4$SCT_snn_res.0.5)
T_CD4_MainCluster_df <- rbind(
    data.frame(x='CD4_Th1',y=c(0,5)),
    data.frame(x='CD4_Th2',y=c(6)),
    data.frame(x='CD4_Treg',y=c(3)),
    data.frame(x='CD4_Tfh',y=c(1,4)),
    data.frame(x='CD4_Naive',y=c(2))
)
T_CD4_MainCluster_df <- T_CD4_MainCluster_df[order(T_CD4_MainCluster_df$y, decreasing = F),]
if (all(T_CD4_MainCluster_df$y == levels(T_CD4$T_CD4_MainCluster)))
    levels(T_CD4$T_CD4_MainCluster) <- T_CD4_MainCluster_df$x

my_plotDim(T_CD4, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD4_MainCluster")
my_CountCluster(T_CD4, 'T_CD4_MainCluster', 'orig.ident2', trend_factor = 8649/9021)

# Idents(TNK) <- 'TNK_MainCluster'
# TNK <- my_AddMeta(TNK, new_ident = T_CD4$T_CD4_MainCluster, Replace = T)
# # my_plotDim(TNK, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_CD4_T_CD4_MainCluster")

# enrichment
Idents(T_CD4) <- 'T_CD4_MainCluster'
Markers_T_CD4_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_df <- my_Markers2df_multiple(Markers_T_CD4_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_T_CD4_df$Cluster_CD4_Th1)
temp_GO(Markers_T_CD4_df$Cluster_CD4_Th2)
temp_GO(Markers_T_CD4_df$Cluster_CD4_Treg)
temp_GO(Markers_T_CD4_df$Cluster_CD4_Tfh)
temp_GO(Markers_T_CD4_df$Cluster_CD4_Naive)

Idents(T_CD4) <- 'SCT_snn_res.0.5'
Markers_T_CD4_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_df <- my_Markers2df_multiple(Markers_T_CD4_list, logFC_threshold = 0.25, positive = T, n_top = 100)
for (i in 0:6) {
    cat('Cluster_', i, "\n", sep = "")
    print(temp_GO(Markers_T_CD4_df[, paste0('Cluster_', i)]))
}
my_plotDim(T_CD4, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "orig.ident2")

# enrich between case-control only in either T_CD8 or T_CD4
# CD8
Idents(T_CD8) <- 'orig.ident2'
Markers_T_CD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_df <- my_Markers2df_multiple(Markers_T_CD8_list, logFC_threshold = 0.25, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_CD8_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_CD8_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]
# CD4
Idents(T_CD4) <- 'orig.ident2'
Markers_T_CD4_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_df <- my_Markers2df_multiple(Markers_T_CD4_list, logFC_threshold = 0.25, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_CD4_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_CD4_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]

# save(T_CD4, T_CD8, file = 'Combined_analysis_T_CD4_8.rda')
load(file = 'Combined_analysis_T_CD4_8.rda')

########################
# subcluster for Myeloid:
Idents(Cryo_merge) <- 'MainCluster'
Myeloid_sub <- subset(Cryo_merge, ident = c('Myeloid', 'Mast'))
DefaultAssay(Myeloid_sub) <- 'RNA'
Myeloid_sub <- SCTransform(Myeloid_sub)
Myeloid_sub <- FindVariableFeatures(Myeloid_sub, selection.method = "vst", nfeatures = 3000)
Myeloid_sub <- RunPCA(Myeloid_sub)
Myeloid_sub <- RunUMAP(Myeloid_sub, dims = 1:10)
Myeloid_sub <- RunTSNE(Myeloid_sub, dims = 1:10)
Myeloid_sub <- FindNeighbors(Myeloid_sub, dims = 1:12)
Myeloid_sub <- FindClusters(Myeloid_sub, resolution = 0.01)

# cluster
Myeloid_sub  <- FindClusters(Myeloid_sub, resolution = seq(1,10,1)/100)
Myeloid_sub  <- FindClusters(Myeloid_sub, resolution = seq(1,10,1)/10)
clustree(Myeloid_sub)

my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident2', title = 'Myeloid_sub Cryo_NonCryo')
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(paste0("0", 1:9), 1:3)))

Idents(Myeloid_sub) <- 'SCT_snn_res.0.1'
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Myeloid_sub R=0.1')
my_DotPlot_split(Myeloid_sub, features = Myeloid_marker) + RotatedAxis()
my_violin(Myeloid_sub, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx')

Myeloid_sub$MainCluster <- as.factor(Myeloid_sub$SCT_snn_res.0.1)
levels(Myeloid_sub$MainCluster) <- c('M', 'M', 'Mast', 'DC', 'M', 'DC')
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "MainCluster")

# Idents(Myeloid_sub) <- 'SCT_snn_res.1'
# my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Myeloid_sub R=1')
# my_DotPlot_split(Myeloid_sub, features = Myeloid_marker) + RotatedAxis()
# my_violin(Myeloid_sub, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx')
# # Cluster 9 chaos?
# Markers_9_list <- FindMarkers(Myeloid_sub, ident.1 = '9', only.pos = T)
# Markers_9_df <- my_Markers2df_1Cluster(Markers_9_list, logFC_threshold = 0.25, positive = T, ntop = 100)
# my_GO(Markers_9_df, ont = 'BP', Simplify = T, return_res = F, show = 30, font.size = 8)

Idents(Myeloid_sub) <- 'MainCluster'
# subcluster for DC:
DC_sub <- subset(Myeloid_sub, ident = 'DC')
DefaultAssay(DC_sub) <- 'RNA'
DC_sub <- SCTransform(DC_sub)
DC_sub <- FindVariableFeatures(DC_sub, selection.method = "vst", nfeatures = 2000)
DC_sub <- RunPCA(DC_sub)
DC_sub <- RunUMAP(DC_sub, dims = 1:10)
DC_sub <- RunTSNE(DC_sub, dims = 1:10)
DC_sub <- FindNeighbors(DC_sub, dims = 1:12)
DC_sub <- FindClusters(DC_sub, resolution = 0.2)

my_plotDim(DC_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident2', title = 'DC_sub Cryo_NonCryo')

Idents(DC_sub) <- 'SCT_snn_res.0.2'
my_plotDim(DC_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'DC_sub R=0.2')
my_DotPlot_split(DC_sub, features = Myeloid_marker[6:8]) + RotatedAxis()
my_violin(DC_sub, features = unlist(Myeloid_marker[6:8]), pt.size = 0, mode = 'mtx')

DC_sub$MainCluster <- as.factor(DC_sub$SCT_snn_res.0.2)
levels(DC_sub$MainCluster) <- c('cDC', 'cDC', 'cDC', 'pDC')
my_plotDim(DC_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "MainCluster")

Myeloid_sub <- my_AddMeta(Myeloid_sub, DC_sub$MainCluster, Replace = T)
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'DC_sub_MainCluster', title = 'DC_sub Cryo_NonCryo')

Idents(Myeloid_sub) <- 'MainCluster'
# subcluster for M(Macro/Mono):
M_sub <- subset(Myeloid_sub, ident = 'M')
DefaultAssay(M_sub) <- 'RNA'
M_sub <- SCTransform(M_sub)
M_sub <- FindVariableFeatures(M_sub, selection.method = "vst", nfeatures = 2000)
M_sub <- RunPCA(M_sub)
M_sub <- RunUMAP(M_sub, dims = 1:10)
M_sub <- RunTSNE(M_sub, dims = 1:10)
M_sub <- FindNeighbors(M_sub, dims = 1:12, k.param = 50)
M_sub <- FindClusters(M_sub, resolution = 1)

my_plotDim(M_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident2', title = 'M_sub Cryo_NonCryo')

# mtx_temp <- as.matrix(M_sub@assays$SCT@data)[c("Adgre1", "Cd14", 'S100a8', 'S100a9'), ]
# pheatmap::pheatmap(
#     mtx_temp,
#     show_colnames = F, 
#     cluster_rows = T, cluster_cols = T)

Idents(M_sub) <- 'SCT_snn_res.1'
my_plotDim(M_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'M_sub R=1')
my_DotPlot_split(M_sub, features = Myeloid_marker[c(2,3,4,5,9)]) + RotatedAxis()
my_violin(M_sub, features = unlist(Myeloid_marker[c(2,3,4,5,9)]), pt.size = 0, mode = 'mtx')
FeaturePlot(M_sub, features = unlist(Myeloid_marker[c(2,3,4,5,9)]))

my_DotPlot_split(M_sub, features = my_MarkersList$T) + RotatedAxis()
my_violin(M_sub, features = my_MarkersList$T, ident = '7', pt.size = 0, mode = 'mtx')
FeaturePlot(M_sub, features = my_MarkersList$T)
my_DotPlot_split(M_sub, features = T_marker$CD8) + RotatedAxis()
my_violin(M_sub, features = unlist(T_marker$CD8), ident = '7', pt.size = 0, mode = 'mtx')

M_sub$MainCluster <- as.factor(M_sub$SCT_snn_res.1)
levels(M_sub$MainCluster) <- c('Macrophage_M2', 'Macrophage_M1', 'Monocyte', 'Monocyte', 'Monocyte', 'Macrophage_M1', 'Macrophage_M1', 'T_CD8(EM_ex)', 'Monocyte')
my_plotDim(M_sub, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "MainCluster")

Idents(Myeloid_sub) <- 'DC_sub_MainCluster'
Myeloid_sub <- my_AddMeta(Myeloid_sub, M_sub$MainCluster, Replace = T)
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'M_sub_MainCluster', title = 'M_sub Cryo_NonCryo')

Idents(Myeloid_sub) <- 'M_sub_MainCluster'
my_CountCluster(Myeloid_sub, group2 = 'orig.ident2', trend_factor = 8649/9021)

# enrichment
Idents(Myeloid_sub) <- 'orig.ident2'
Markers_Myeloid_sub_list <- FindAllMarkers(Myeloid_sub, only.pos = T)
Markers_Myeloid_sub_df <- my_Markers2df_multiple(Markers_Myeloid_sub_list, logFC_threshold = 0.25, positive = T, n_top = 120)
GO_cryo <- temp_GO(Markers_Myeloid_sub_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_Myeloid_sub_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]

Idents(Myeloid_sub) <- 'M_sub_MainCluster'
Markers_Myeloid_sub_list <- FindAllMarkers(Myeloid_sub, only.pos = T)
Markers_Myeloid_sub_df <- my_Markers2df_multiple(Markers_Myeloid_sub_list, logFC_threshold = 0.25, positive = T, n_top = 200)
temp_GO(Markers_Myeloid_sub_df$Cluster_Monocyte)
temp_GO(Markers_Myeloid_sub_df$Cluster_cDC)
temp_GO(Markers_Myeloid_sub_df$Cluster_Mast)
temp_GO(Markers_Myeloid_sub_df$Cluster_Macrophage_M1)
temp_GO(Markers_Myeloid_sub_df$Cluster_Macrophage_M2)
temp_GO(Markers_Myeloid_sub_df$Cluster_pDC)

#####################################
# PLOT FOR pDC
temp_pDC <- my_GO(Markers_Myeloid_sub_df$Cluster_pDC, return_plot = F, ont = 'BP', Simplify = T, return_res = T, type = 'bar', 
                  pAdjustMethod = 'BH', pvalue = 0.101, qvalue = 0.1)
barplot(temp_pDC, showCategory = 33)

temp_pDC1 <- my_GO(Markers_Myeloid_sub_df$Cluster_pDC, return_plot = F, ont = 'BP', Simplify = T, return_res = T, type = 'bar', 
                  pAdjustMethod = 'BH', pvalue = 0.05, qvalue = 0.2)
temp_pDC1@result <- temp_pDC1@result[grep("interferon",temp_pDC1@result$Description),]
temp_pDC@result <- temp_pDC1@result
temp_pDC@pvalueCutoff <- 0.2
temp_pDC@qvalueCutoff <- 0.2
dotplot(temp_pDC, showCategory = 20)
#####################################

pDC_interferon_gene <- c("Tlr7","Flt3","Ptprs","Il12a","Tnfrsf13c","Kynu")
FeaturePlot(Myeloid_sub, features = pDC_interferon_gene, pt.size = 2)
my_DotPlot_split(Myeloid_sub, features = pDC_interferon_gene) + RotatedAxis()
my_violin(Myeloid_sub, features = pDC_interferon_gene, pt.size = 0, mode = 'mtx')

Idents(Myeloid_sub) <- 'SCT_snn_res.1'
my_CountCluster(Myeloid_sub, group2 = 'orig.ident2')
my_plotDim(Myeloid_sub, reduction = "umap", label = F, pt.size = 1, label.size = 5, group.by = "orig.ident2")
my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, title = 'Myeloid_sub R=1')
Markers_Myeloid_sub_list <- FindAllMarkers(Myeloid_sub, only.pos = T)
Markers_Myeloid_sub_df <- my_Markers2df_multiple(Markers_Myeloid_sub_list, logFC_threshold = 0.25, positive = T, n_top = 200)
for (i in 0:16) {
    cat('Cluster_', i, "\n", sep = "")
    print(temp_GO(Markers_Myeloid_sub_df[, paste0('Cluster_', i)]))
}

cDC_GO <- temp_GO(Markers_Myeloid_sub_df$Cluster_12, pathway_name = 'interferon-alpha')
cDC_interferon_gene <- c("Havcr2", "Mmp12",  "Flt3")
FeaturePlot(Myeloid_sub, features = cDC_interferon_gene, pt.size = 2)
my_DotPlot_split(Myeloid_sub, features = cDC_interferon_gene) + RotatedAxis()
my_violin(Myeloid_sub, features = cDC_interferon_gene, pt.size = 0, mode = 'mtx')

# save(Myeloid_sub, DC_sub, M_sub, file = 'Combined_analysis_Myeloid.rda')
load(file = 'Combined_analysis_Myeloid.rda')

########################
# cell type defined
Cryo_merge_temp <- my_AddMeta(Cryo_merge, T_CD4$T_CD4_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD4_T_CD4_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, T_CD8$T_CD8_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD8_T_CD8_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, Myeloid_sub$M_sub_MainCluster, Replace = T)
Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster[Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster == 'T_CD8(EM_ex)'] <- 'CD8_Effector_Memory'
Idents(Cryo_merge_temp) <- 'Myeloid_sub_M_sub_MainCluster'

my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 1.5, label.size = 6, group.by = 'Myeloid_sub_M_sub_MainCluster', title = 'Cryo_merge Cluster')
my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 1.5, label.size = 6, group.by = 'orig.ident2', title = 'Cryo_merge Cryo_NonCryo')
my_CountCluster(Cryo_merge_temp, group2 = 'orig.ident2')
########################


########################PCA analysis
T_CD8 <- RunPCA(T_CD8)
tmp_PCA <- sort(T_CD8@reductions$pca@feature.loadings[,1], decreasing = T)
tmp_pos <- tmp_PCA[1:30]
tmp_neg <- rev(tmp_PCA)[1:30]
PCA_top30 <- names(c(tmp_pos, tmp_neg))
PCA_GO <- my_GO(PCA_top30, ont = 'BP', type = 'dot', Simplify = T, return_plot = T, return_res = T)
PCA_GO_ifn <- tmp_GO@result[grep("interferon", tmp_GO@result$Description), ]
PCA_GO_ifn_gene <- bitr(strsplit(PCA_GO_ifn$geneID, split = "/")[[1]], fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
PCA_top30
PCA_GO_ifn_gene
# "Hmgb2" "Hmgb1"
# positive regulation of interferon-beta production

mtx_tmp0 <- T_CD8@assays$SCT@data
mtx_tmp0 <- t(apply(mtx_tmp0, 1, scale))
colnames(mtx_tmp0) <- colnames(T_CD8@assays$SCT@data)
mtx_tmp <- mtx_tmp0[PCA_top30, ]
anno_col <- T_CD8$orig.ident2
anno_col <- sort(anno_col, decreasing = T)
anno_col <- as.data.frame(anno_col)
pheatmap::pheatmap(mtx_tmp, 
                   color = grDevices::colorRampPalette(colors = c("blue","white","red"))(100),
                   annotation_col = anno_col, 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   show_colnames = F, 
                   show_rownames = T)

T_sub <- RunPCA(T_sub)
tmp_PCA <- sort(T_sub@reductions$pca@feature.loadings[,1], decreasing = T)
tmp_pos <- tmp_PCA[1:30]
tmp_neg <- rev(tmp_PCA)[1:30]
PCA_top30 <- names(c(tmp_pos, tmp_neg))
PCA_GO <- my_GO(PCA_top30, ont = 'BP', type = 'dot', Simplify = T, return_plot = T, return_res = T)
PCA_GO_ifn <- tmp_GO@result[grep("interferon", tmp_GO@result$Description), ]
PCA_GO_ifn_gene <- bitr(strsplit(PCA_GO_ifn$geneID, split = "/")[[1]], fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
PCA_top30
PCA_GO_ifn_gene
# "Hmgb2" "Hmgb1"
# positive regulation of interferon-beta production

mtx_tmp0 <- T_sub@assays$SCT@data
mtx_tmp0 <- t(apply(mtx_tmp0, 1, scale))
colnames(mtx_tmp0) <- colnames(T_sub@assays$SCT@data)
mtx_tmp <- mtx_tmp0[PCA_top30, ]
anno_col <- T_sub$orig.ident2
anno_col <- sort(anno_col, decreasing = T)
anno_col <- as.data.frame(anno_col)
pheatmap::pheatmap(mtx_tmp, 
                   color = grDevices::colorRampPalette(colors = c("blue","white","red"))(100),
                   annotation_col = anno_col, 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   show_colnames = F, 
                   show_rownames = T)
########################


########################Communication
cellchatCreate <- function(object, object.identity, DBtype = c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor")[1]) {
    library(CellChat)
    library(Seurat)
    if (class(object)[1] == "Seurat") {
        data.input <- object@assays$RNA@data
    } else {
        data.input <- object
    }
    objectx <- createCellChat(data.input)
    objectx <- addMeta(objectx, meta = object.identity, meta.name = "labels")
    objectx <- setIdent(objectx, ident.use = "labels")
    CellChatDB <- CellChatDB.mouse
    CellChatDB.use <- subsetDB(CellChatDB, search = DBtype)
    objectx@DB <- CellChatDB.use
    objectx <- subsetData(objectx)  
    future::plan("multiprocess", workers = 4)
    objectx <- identifyOverExpressedGenes(objectx)  #很耗时
    objectx <- identifyOverExpressedInteractions(objectx)
    objectx <- projectData(objectx, PPI.mouse)  #小耗时
    objectx <- computeCommunProb(objectx)  #很耗时
    objectx <- computeCommunProbPathway(objectx)
    objectx <- aggregateNet(objectx)
    
    return(objectx)
}
cellchat_obj <- cellchatCreate(Cryo_merge_temp, 
                               Cryo_merge_temp@meta.data[, c('Myeloid_sub_M_sub_MainCluster', 'orig.ident', 'orig.ident2', 'orig.ident3')]
                               )
save(cellchat_obj, file = "Combined_CellChat_obj.rda")
cellchat_obj@netP$pathways
groupSize <- table(cellchat_obj@idents)
# lapply(
#     cellchat_obj@netP$pathways, 
#     function(x) 
#         netVisual_aggregate(cellchat_obj, 
#                             signaling = x, 
#                             layout = "circle", 
#                             pt.title=20, 
#                             vertex.label.cex = 1.7, 
#                             vertex.size = groupSize)
# )
netVisual_aggregate(cellchat_obj, 
                    signaling = cellchat_obj@netP$pathways[8], 
                    layout = "circle", 
                    pt.title=20, 
                    vertex.label.cex = 1.7, 
                    vertex.size = groupSize)
netVisual_aggregate(cellchat_obj, 
                    signaling = cellchat_obj@netP$pathways[11], 
                    layout = "circle", 
                    pt.title=20, 
                    vertex.label.cex = 1.7, 
                    vertex.size = groupSize)
netVisual_aggregate(cellchat_obj, 
                    signaling = cellchat_obj@netP$pathways[18], 
                    layout = "circle", 
                    pt.title=20, 
                    vertex.label.cex = 1.7, 
                    vertex.size = groupSize)
netVisual_aggregate(cellchat_obj, 
                    signaling = cellchat_obj@netP$pathways[20], 
                    layout = "circle", 
                    pt.title=20, 
                    vertex.label.cex = 1.7, 
                    vertex.size = groupSize)
########################