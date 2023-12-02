source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/211111/")

# Cryo_data_1 <- Read10X_h5("Cryo_1wk/filtered_feature_bc_matrix.h5")
# NonCryo_data_1 <- Read10X_h5("NonCryo_1wk/filtered_feature_bc_matrix.h5")
# Cryo_data_2 <- Read10X_h5("Cryo_2wk/filtered_feature_bc_matrix.h5")
# NonCryo_data_2 <- Read10X_h5("NonCryo_2wk/filtered_feature_bc_matrix.h5")
# 
# colnames(Cryo_data_1) <- paste0("Cryo1wk_", colnames(Cryo_data_1))
# colnames(NonCryo_data_1) <- paste0("NonCryo1wk_", colnames(NonCryo_data_1))
# colnames(Cryo_data_2) <- paste0("Cryo2wk_", colnames(Cryo_data_2))
# colnames(NonCryo_data_2) <- paste0("NonCryo2wk_", colnames(NonCryo_data_2))
# 
# Cryo_1 <- CreateSeuratObject(counts = Cryo_data_1,
#                              project = "Cryo1wk",
#                              min.cells = 3,
#                              min.features = 200)
# NonCryo_1 <- CreateSeuratObject(counts = NonCryo_data_1,
#                                 project = "NonCryo1wk",
#                                 min.cells = 3,
#                                 min.features = 200)
# Cryo_2 <- CreateSeuratObject(counts = Cryo_data_2,
#                              project = "Cryo2wk",
#                              min.cells = 3,
#                              min.features = 200)
# NonCryo_2 <- CreateSeuratObject(counts = NonCryo_data_2,
#                                 project = "NonCryo2wk",
#                                 min.cells = 3,
#                                 min.features = 200)
# 
# Cryo_1@meta.data$orig.ident <- Cryo_1@project.name
# NonCryo_1@meta.data$orig.ident <- NonCryo_1@project.name
# Cryo_2@meta.data$orig.ident <- Cryo_2@project.name
# NonCryo_2@meta.data$orig.ident <- NonCryo_2@project.name
# 
# Cryo_combined <- merge(Cryo_1, c(NonCryo_1, Cryo_2, NonCryo_2))
# 
# Cryo_combined[["percent.mt"]] <- PercentageFeatureSet(Cryo_combined, pattern = "^mt-")
# Idents(Cryo_combined) <- 'orig.ident'
# Cryo_combined$orig.ident_group <- Cryo_combined$orig.ident
# Cryo_combined$orig.ident_group <- sub("^Cryo.*$", "Cryo", Cryo_combined$orig.ident_group)
# Cryo_combined$orig.ident_group <- sub("^NonCryo.*$", "NonCryo", Cryo_combined$orig.ident_group)
# Cryo_combined$orig.ident_time <- Cryo_combined$orig.ident
# Cryo_combined$orig.ident_time <- sub("^.*1wk$", "1wk", Cryo_combined$orig.ident_time)
# Cryo_combined$orig.ident_time <- sub("^.*2wk$", "2wk", Cryo_combined$orig.ident_time)
# 
# # visualization
# #######################
# vln_cryo_m <- my_plotQC(Cryo_combined, pt.size = 0)
# FC_cryo_m <- my_plotFC(Cryo_combined)
# 
# hist(Cryo_combined$nFeature_RNA, breaks = 500, xaxt = 'n')
# axis(1, seq(0,10000,50))
# hist(Cryo_combined$nCount_RNA, breaks = 500, xaxt = 'n')
# axis(1, seq(0,200000,1000))
# hist(Cryo_combined$percent.mt, breaks = 500, xaxt = 'n')
# axis(1, seq(0,100,1))
# #######################
# Cryo_merge <- subset(Cryo_combined, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA < 50000 & percent.mt < 17)
# Cryo_merge
# 
# Cryo_merge <- SCTransform(Cryo_merge)
# # Cryo_merge <- NormalizeData(Cryo_merge)
# # Cryo_merge <- ScaleData(Cryo_merge)
# 
# Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = 2000)
# Cryo_merge <- RunPCA(Cryo_merge)
# 
# Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
# Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.01)
# 
# Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10)
# Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10)
# 
# # visualization
# #######################
# my_plotDim(Cryo_merge, reduction = "pca", label = T, pt.size = 1, label.size = 5, title = "Cluster")
# my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "Cluster")
# my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, label.size = 5, title = "Cluster")
# 
# my_plotDim(Cryo_merge, reduction = "pca", label = F, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
# my_plotDim(Cryo_merge, reduction = "umap", label = F, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
# my_plotDim(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, group.by = 'orig.ident', title = 'Cryo_NonCryo')
# #######################
# Cryo_merge  <- FindClusters(Cryo_merge, resolution = seq(1,10,1)/100)
# clustree(Cryo_merge)
# 
# Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.06)
# my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "R=0.06")
# 
# DotPlot(Cryo_merge, features = my_MarkersList[[8]]) + RotatedAxis()
# DotPlot(Cryo_merge, features = my_MarkersList[-8]) + RotatedAxis()
# custom_markers_Cryomain <- list(T = c('Cd3e', 'Cd247'), NK = c('Klrk1', 'Klrb1c'), B = c('Cd19', 'Cd79a'), Myeloid = c('Itgam', 'Cd14'), Mast = c('Mcpt1', 'Mcpt2'))
# my_DotPlot_split(Cryo_merge, features = custom_markers_Cryomain) + RotatedAxis()
# my_violin(Cryo_merge, features = unlist(custom_markers_Cryomain), pt.size = 0, mode = 'mtx')
# 
# # cell type identify
# Cryo_merge$MainCluster <- Cryo_merge$SCT_snn_res.0.06
# MainCluster <- rbind(
#     data.frame(x='T',y=c(0,4,6,9)),
#     data.frame(x='NK',y=c(7)),
#     data.frame(x='B',y=c(2)),
#     data.frame(x='Myeloid',y=c(1,3,8)),
#     data.frame(x='Mast',y=c(5))
# )
# MainCluster <- MainCluster[order(MainCluster$y, decreasing = F),]
# if(all(MainCluster$y == levels(Cryo_merge$MainCluster)))
#     levels(Cryo_merge$MainCluster) <- MainCluster$x
# Idents(Cryo_merge) <- 'MainCluster'
# my_plotDim(Cryo_merge, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = "MainCluster")
# 
# # count
# my_CountCluster(Cryo_merge, 'MainCluster', 'orig.ident2')
# my_CountCluster(Cryo_merge, 'MainCluster', 'orig.ident3')


library(ProjecTILs)
ref_TILAtlas <- readRDS("E://Cryo-TCR/TILAtlas/ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref_TILAtlas,label = T, cols = refCols)

load("Cryo_T_analysis.rda")
DimPlot(Cryo_merge_T)
query_T <- make.projection(Cryo_merge_T, ref = ref_TILAtlas, skip.normalize = T, fast.mode = T)
refCols2 <- c(refCols, "black", "yellow", "green")
DimPlot(query_T, group.by = 'T_CD8_CD8_Cluster', cols = refCols2)
DimPlot(query_T, group.by = 'T_CD8_CD8_Cluster', reduction = 'pca')
query_T <- FindNeighbors(query_T, dims = 1:12)
query_T <- FindClusters(query_T, resolution = 1)


ClusterSCT4Query <- function(Object, ref_obj, resolution, obj_assay = 'SCT') {
    name_tmp <- paste0(obj_assay, "_snn_res.", resolution)
    cat(name_tmp, "\n")
    DefaultAssay(Object) <- obj_assay
    Object <- FindClusters(Object, resolution = resolution)
    # ref_obj@meta.data[, name_tmp] <- Object@meta.data[, name_tmp]
    eval(parse(text = paste0("ref_obj$", name_tmp, " <- Object$", name_tmp)))
    print(DimPlot(ref_obj, group.by = name_tmp, reduction = 'umap'))
    return(ref_obj)
}
ClusterSCT4Query(Object = Cryo_merge_T, ref_obj = query_T, resolution = 1)
ClusterSCT4Query(Object = Cryo_merge_T, ref_obj = query_T, resolution = 0.5)
ClusterSCT4Query(Object = Cryo_merge_T, ref_obj = query_T, resolution = 1)

Cryo_merge_T <- FindClusters(Cryo_merge_T, resolution = 1)
query_T$SCT_snn_res.1 <- Cryo_merge_T$SCT_snn_res.1
DimPlot(query_T, group.by = 'SCT_snn_res.1', reduction = 'umap')

DimPlot(query_T, group.by = 'integrated_snn_res.1', reduction = 'umap')

markers <- c("Cd4","Cd8a","Sell","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(query_T, features=markers,stack = T,flip = T,assay = "SCT", group.by = 'integrated_snn_res.1')
VlnPlot(query_T, features=markers,stack = T,flip = T,assay = "SCT", group.by = 'T_CD8_CD8_Cluster')
VlnPlot(ref_TILAtlas, features=markers,stack = T,flip = T,assay = "RNA")


plot.projection(ref_TILAtlas, query_T)

query.projected <- cellstate.predict(ref=ref_TILAtlas, query=query_T)
table(query.projected$functional.cluster)
plot.statepred.composition(ref_TILAtlas, query.projected,metric = "Percent")

plot.states.radar(ref_TILAtlas, query=query.projected, min.cells=30)

query.control <- subset(query.projected, subset=`orig.ident_group` == 'Cryo')
query.perturb <- subset(query.projected, subset=`orig.ident_group` == 'NonCryo')

plot.states.radar(ref_TILAtlas, query=list("Control" = query.control, "Query" = query.perturb))

discriminantGenes <- find.discriminant.genes(ref=ref_TILAtlas, query=query.perturb, query.control=query.control, state="CD8_Tex", query.assay = 'RNA')
discriminantGenes
discriminantGenes[discriminantGenes$p_val_adj < 0.05 & discriminantGenes$avg_log2FC > 0, ]
