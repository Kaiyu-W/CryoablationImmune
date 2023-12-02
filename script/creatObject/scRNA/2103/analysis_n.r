setwd("D://Ji_Wangxue/")
library(Seurat)
library(ggplot2)

cluster_number <- function(cluster_object, Percentage = FALSE, Round = FALSE) {
  df <- data.frame(
    cluster=cluster_object@active.ident, 
    project=cluster_object$orig.ident
  )
  fren <- table(df)
  fren_perc <- as.matrix(fren)
  fren_perc[,1] <- fren_perc[,1] / sum(fren_perc[,1]) * (100 ** Percentage)
  fren_perc[,2] <- fren_perc[,2] / sum(fren_perc[,2]) * (100 ** Percentage)
  fren_prop <- round(fren_perc[,1] / fren_perc[,2], 3)
  if (Round)
    fren_perc <- round(fren_perc, ifelse(Percentage, 1, 3))
  res <- cbind(fren, fren_perc, fren_prop)
  res <- rbind(res, apply(res, 2, sum))
  res[nrow(res), ncol(res)] <- NA
  names(dimnames(res)) <- c("cluster", "project")
  dimnames(res)[[1]][nrow(res)] <- "sum" 
  dimnames(res)[[2]] <- c("Cryo_count", "NonCryo_count", 
                          "Cryo_precentage", "NonCryo_percentage", "trend")
  res
}

#load data
Cryo_merge <- readRDS("Cryo_merge_sct.rds")

#Clustering
Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.05)
DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
DimPlot(Cryo_merge, reduction = "tsne", label = T, pt.size = 1)
table(Cryo_merge@active.ident)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.3)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.6)
DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
table(Cryo_merge@active.ident)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.9)
DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
table(Cryo_merge@active.ident)
#Dimension Reduction : UMAP/TSNE
# reticulate::py_install(packages = 'umap-learn')
Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10, seed.use = 0)
Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10, seed.use = 0)

#Visualization

DimPlot(Cryo_merge, reduction = "pca")
DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
DimPlot(Cryo_merge, reduction = "tsne", label = T, pt.size = 1)

#3D-Visualization
source('C://Users/PC/Desktop/3DPlot_seurat.r')
Cryo_merge_3d <- RunUMAP(Cryo_merge, dims = 1:10, seed.use = 0, n.components = 3L)
#Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.05)
DimPlot_3d(Cryo_merge_3d, meta = 'seurat_clusters', size = 4)
FeaturePlot_3d(Cryo_merge_3d, gene = "Cd3e", Colors = c("blue", "red"), cutoff = c(2.5, 'p99')[1], size = 4)


#Marker genes
Markers_cryo_m <- FindAllMarkers(Cryo_merge, test.use = 't')
# If there isn't batch effect between samples, output genes from function FindMarkers can be clusters markers
# However, if batch effect exists, should use function FindConservedMarkers to find overlap marker genes between sample batches
# ConservedMarkers <- lapply(
#   c("T","NK","B","Myeloid","Proliferating"), 
#   function(x) FindConservedMarkers(Cryo_merge, ident.1 = x, grouping.var = 'orig.ident')
#   )
# saveRDS(ConservedMarkers, file = 'ConservedMarkers.rds')

trans_df_markers <- function(df_input, n_top) {
  cluster <- levels(df_input$cluster)
  res <- rep(NA, n_top*length(cluster))
  dim(res) <- c(n_top, length(cluster))
  dimnames(res) <- list(1:n_top, paste0("Cluster", cluster))
  for (i in 1:length(cluster)) {
    cluster_name <- cluster[i]
    cluster_index <- which(df_input$cluster == cluster_name)[1:n_top]
    cluster_gene <- df_input$gene[cluster_index] 
    res[, i] <- cluster_gene
  }
  return(res)
}

trans_df_conservedmarkers <- function(df_input, n_top) {
  cluster <- levels(df_input$cluster)
  res <- rep(NA, n_top*length(cluster))
  dim(res) <- c(n_top, length(cluster))
  dimnames(res) <- list(1:n_top, paste0("Cluster", cluster))
  for (i in 1:length(cluster)) {
    cluster_name <- cluster[i]
    cluster_index <- which(df_input$cluster == cluster_name)[1:n_top]
    cluster_gene <- df_input$gene[cluster_index] 
    res[, i] <- cluster_gene
  }
  return(res)
}

Markers_cryo_merge_df <- trans_df_markers(Markers_cryo_m, 1000)

Markers <- list()
Markers$Leukocytes <- "Ptprc"
Markers$Lymphoid_Cells <- ""
Markers$B_Cells <- "Cd19"
Markers$Plasma_Cells <- c("Tnfrsf17", "Sdc1")
Markers$NK_Cells <- c("Klrb1c", "Ncr1", "Klrk1", "Cd3e", "Cd247")
Markers$NKT_Cells <- c("Cd3e", "Cd247", "Klrb1c")
Markers$T_Cells <- c("Cd3e", "Cd247")
Markers$Myeloid_Cells <- "Itgam"
Markers$Monocytes <- c("Cd14", "Adgre1")
Markers$Macrophages <- "Adgre1"
Markers$cDCs <- c("Itgax", "H2-Ab1")
Markers$pDCs <- c("Siglech", "Bst2", "Cd83")
Markers$MDSCs <- c("Arg1","Arg2")
Markers$Granulocytes <- ""

#rank
#merge
markers_vector <- character()
for (i in names(Markers)) markers_vector <- append(markers_vector, Markers[[i]])
markers_vector <- unique(markers_vector)
markers_vector <- markers_vector[markers_vector != ""]
markers_rank <- matrix(nrow = length(markers_vector), 
                       ncol = ncol(Markers_cryo_merge_df), 
                       dimnames = list(
                         markers_vector, 
                         colnames(Markers_cryo_merge_df)))
for (i in 1:ncol(Markers_cryo_merge_df)) {
  index <- sapply(markers_vector, function(x) which(Markers_cryo_merge_df[, i] %in% x))
  index <- unlist(index)
  markers_rank[rownames(markers_rank) %in% names(index), i] <- unname(index)
  markers_rank[-which(rownames(markers_rank) %in% names(index)), i] <- 1000
}

heatmap(t(1000 - markers_rank), Rowv = NA, Colv = NA)


###Plot markers in umap:
###Merge
#B_Cells:CD19
FeaturePlot(Cryo_merge, Markers$B_Cells, reduction = 'umap', pt.size = 1)
#Plasma_Cells
for (i in Markers$Plasma_Cells) print(FeaturePlot(Cryo_merge, i, reduction = 'umap', pt.size = 1))
#NK_Cells
for (i in Markers$NK_Cells[1:2]) print(FeaturePlot(Cryo_merge, i, reduction = 'umap', pt.size = 1))
FeaturePlot(Cryo_merge, "Klrk1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd3e", cols=c("blue", "white"), reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd247", cols=c("blue", "white"), reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Klrk1", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd3e", cols=c("blue", "white"), reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd247", cols=c("blue", "white"), reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo_merge, c("Klrk1","Cd247"), pt.size = 1, reduction = 'umap', blend = T)
FeaturePlot(Cryo_merge, c("Klrk1","Cd247"), pt.size = 1, reduction = 'tsne', blend = T)
#NKT_Cells
FeaturePlot(Cryo_merge, "Klrb1c", pt.size = 1)
FeaturePlot(Cryo_merge, "Cd3e", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd247", reduction = 'umap', pt.size = 1)
#T_Cells
FeaturePlot(Cryo_merge, "Cd3e", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Cd247", reduction = 'umap', pt.size = 1)
#CD4
FeaturePlot(Cryo_merge, "Cd4", reduction = 'umap', pt.size = 1)
#Cd8
FeaturePlot(Cryo_merge, "Cd8a", reduction = 'umap', pt.size = 1)
#Myeloid_Cells
FeaturePlot(Cryo_merge, Markers$Myeloid_Cells, reduction = 'umap', pt.size = 1)
#Monocytes
FeaturePlot(Cryo_merge, "Cd14", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Adgre1", cols=c("blue", "white"), reduction = 'umap', pt.size = 0.1)
  #M1:Cd80,Cd86,Nos2
  FeaturePlot(Cryo_merge, "Cd80", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo_merge, "Cd86", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo_merge, "Nos2", reduction = 'umap', pt.size = 1)
  #M2:Cd163,Mrc1,Arg1/Arg2
  FeaturePlot(Cryo_merge, "Cd163", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo_merge, "Mrc1", reduction = 'umap', pt.size = 1)
#Macrophages
FeaturePlot(Cryo_merge, "Adgre1", reduction = 'umap', pt.size = 1)
#cDCs
FeaturePlot(Cryo_merge, "Itgax", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "H2-Ab1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, c("H2-Ab1",'Itgax'), reduction = 'umap', pt.size = 1, blend = T, combine = F)
#pDCs
FeaturePlot(Cryo_merge, 'Siglech', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Bst2', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, c("Siglech",'Bst2'), reduction = 'umap', pt.size = 1, blend = T, cols = c("lightgrey",'red','blue'), combine = F)
FeaturePlot(Cryo_merge, Markers$pDCs[3], reduction = 'umap', pt.size = 1)
#MDSCs
for (i in Markers$MDSCs) print(FeaturePlot(Cryo_merge, i, reduction = 'umap', pt.size = 1))
#Granulocytes
FeaturePlot(Cryo_merge, 'Ifng', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Il6', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Ccl2', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Mpo', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Tnf', reduction = 'umap', pt.size = 1)

#
FeaturePlot(Cryo_merge, "Ptprc", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Epcam", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Pecam1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Mme", reduction = 'umap', pt.size = 1)

#CD4 T_cells : Il7r
FeaturePlot(Cryo_merge, "Il7r", reduction = 'umap', pt.size = 1)
#CD14+ Monocytes
FeaturePlot(Cryo_merge, "Cd14", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, "Il7r", reduction = 'umap', pt.size = 1)
#B cells : Ms4a1
FeaturePlot(Cryo_merge, 'Ms4a1', reduction = 'umap', pt.size = 1)
#CD8 T cells:	CD8a
FeaturePlot(Cryo_merge, 'CD8a', reduction = 'umap', pt.size = 1)
#NK cells: GNLY, NKG7
FeaturePlot(Cryo_merge, '', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Nkg7', reduction = 'umap', pt.size = 1)
#FCGR3A+ Monocytes: FCGR3A, MS4A7
FeaturePlot(Cryo_merge, 'Fcgr3', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Ms4a7', reduction = 'umap', pt.size = 1)
#Dendritic Cells: FCER1A, CST3
FeaturePlot(Cryo_merge, 'Fcer1a', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo_merge, 'Cst3', reduction = 'umap', pt.size = 1)
#Megakaryocytes:	PPBP
FeaturePlot(Cryo_merge, 'Ppbp', reduction = 'umap', pt.size = 1)

####
# FeaturePlot(Cryo_merge, Markers_cryo_merge_df[1:4,5], reduction = 'umap', pt.size = 1)


## dotplot for markers to clusters
# marker_gene <- c(
#   "Fcer1a","Tpsb2","Cpa3",'Cd63','Alox5ap',
#   "Cd79a","Ebf1",
#   "Ctsh","Ly86","Lyz2",
#   "Cebpb","Cxcl2","Il1rn","Il1b","Il1r2","Ets2","Ccr1","Cd14",
#   "Gzma","Irf8","Nrarp","Klre1","Txk","Car2",
#   "Cd3g","Trbc2","Cd3d","Trac","Cd3e",
#   "Stmn1","Tuba1b","Mki67","Hmgn2","Birc5","Pclaf"
# )
marker_gene <- list(
  Proliferating = c("Stmn1","Tuba1b","Mki67","Hmgn2","Birc5","Pclaf"),
  B = c("Cd79a","Ebf1"),
  NK = c("Gzma","Irf8","Nrarp","Klre1","Txk","Car2"),
  Myeloid = c("Fcer1a","Tpsb2","Cpa3",'Cd63','Alox5ap',"Ctsh","Ly86","Lyz2","Cebpb","Cxcl2","Il1rn","Il1b","Il1r2","Ets2","Ccr1","Cd14"),
  T = c("Cd3g","Trbc2","Cd3d","Trac","Cd3e")
)

#define clusters to cell types
levels(Cryo_merge$seurat_clusters) <- c(
  "T",
  "Myeloid",
  "NK",
  "Myeloid",
  "Myeloid",
  "B",
  "Proliferating"
  )
Cryo_merge <- SetIdent(Cryo_merge, value = "seurat_clusters")

##plot
DimPlot(Cryo_merge, reduction = "umap", group.by = 'seurat_clusters', label = T, pt.size = 1)
DotPlot(Cryo_merge, features = marker_gene) + RotatedAxis()

##count clusters
table(Cryo_merge$seurat_clusters)
# df <- data.frame(
#   cluster=Cryo_merge$seurat_clusters, 
#   project=Cryo_merge$orig.ident
# )
# fren <- table(df)
# fren_perc <- as.matrix(fren)
# fren_perc[,1] <- fren_perc[,1] / sum(fren_perc[,1])
# fren_perc[,2] <- fren_perc[,2] / sum(fren_perc[,2])
# fren
# fren_perc
cluster_number <- function(cluster_object, Percentage = FALSE, Round = FALSE) {
  df <- data.frame(
    cluster=cluster_object@active.ident, 
    project=cluster_object$orig.ident
  )
  fren <- table(df)
  fren_perc <- as.matrix(fren)
  fren_perc[,1] <- fren_perc[,1] / sum(fren_perc[,1]) * (100 ** Percentage)
  fren_perc[,2] <- fren_perc[,2] / sum(fren_perc[,2]) * (100 ** Percentage)
  fren_prop <- round(fren_perc[,1] / fren_perc[,2], 3)
  if (Round)
    fren_perc <- round(fren_perc, ifelse(Percentage, 1, 3))
  res <- cbind(fren, fren_perc, fren_prop)
  res <- rbind(res, apply(res, 2, sum))
  res[nrow(res), ncol(res)] <- NA
  names(dimnames(res)) <- c("cluster", "project")
  dimnames(res)[[1]][nrow(res)] <- "sum" 
  dimnames(res)[[2]] <- c("Cryo_count", "NonCryo_count", 
                          "Cryo_precentage", "NonCryo_percentage", "trend")
  res
}

Cryo_m_number <- cluster_number(Cryo_merge, Percentage = T, Round = T)
Cryo_m_number

# ##sub-clustering data preparation
# B_merge <- subset(Cryo_merge, idents = "B")
# T_merge <- subset(Cryo_merge, idents = "T")
# Myeloid_merge <- subset(Cryo_merge, idents = "Myeloid")
# NK_merge <- subset(Cryo_merge, idents = "NK")
# Proliferating_merge <- subset(Cryo_merge, idents = "Proliferating")
# 
# #pca shows much heterogeneity of macro/mono/-phils
# DimPlot(Cryo_merge, reduction = "pca", group.by = "SCT_snn_res.0.05", label = T)
# 
# #T
# T_merge <- SCTransform(T_merge, assay = "RNA", new.assay.name = "SCT")
# T_merge <- FindVariableFeatures(T_merge, selection.method = "vst", nfeatures = 2000)
# T_merge <- RunPCA(T_merge)
# T_merge <- RunUMAP(T_merge, dims = 1:10, seed.use = 1)
# T_merge <- RunTSNE(T_merge, dims = 1:10, seed.use = 0)
# T_merge <- FindNeighbors(T_merge, dims = 1:12)
# T_merge <- FindClusters(T_merge, resolution = 0.2)
# DimPlot(T_merge, reduction = "umap", label = T, pt.size = 1)
# DimPlot(T_merge, reduction = "umap", group.by = "orig.ident", label = T, pt.size = 1)
# DimPlot(T_merge, reduction = "tsne", label = T, pt.size = 1)
# 
# T_count <- as.matrix(T_merge@assays$RNA@counts)
# T_count <- unname(unlist(apply(T_count, 2, sum)))
# hist(T_count, breaks = 100)
# 
# Markers_T_m_wilcox <- FindAllMarkers(T_merge, test.use = 'wilcox')
# Markers_T_m_ttest <- FindAllMarkers(T_merge, test.use = 't')
# #write
# #
# 
# #NK
# NK_merge <- SCTransform(NK_merge, assay = "RNA", new.assay.name = "SCT")
# NK_merge <- FindVariableFeatures(NK_merge, selection.method = "vst", nfeatures = 2000)
# NK_merge <- RunPCA(NK_merge)
# NK_merge <- RunUMAP(NK_merge, dims = 1:10, seed.use = 0)
# NK_merge <- RunTSNE(NK_merge, dims = 1:10, seed.use = 0)
# NK_merge <- FindNeighbors(NK_merge, dims = 1:12)
# NK_merge <- FindClusters(NK_merge, resolution = 0.3)
# DimPlot(NK_merge, reduction = "umap", label = T, pt.size = 1)
# DimPlot(NK_merge, reduction = "umap", group.by = "orig.ident", label = T, pt.size = 1)
# DimPlot(NK_merge, reduction = "tsne", label = T, pt.size = 1)
# 
# NK_count <- as.matrix(NK_merge@assays$RNA@counts)
# NK_count <- unname(unlist(apply(NK_count, 2, sum)))
# hist(NK_count, breaks = 100)
# 
# Markers_NK_m_wilcox <- FindAllMarkers(NK_merge, test.use = 'wilcox')
# Markers_NK_m_ttest <- FindAllMarkers(NK_merge, test.use = 't')
# 
# #Myeloid
# Myeloid_merge <- SCTransform(Myeloid_merge, assay = "RNA", new.assay.name = "SCT")
# Myeloid_merge <- FindVariableFeatures(Myeloid_merge, selection.method = "vst", nfeatures = 2000)
# Myeloid_merge <- RunPCA(Myeloid_merge)
# Myeloid_merge <- RunUMAP(Myeloid_merge, dims = 1:10, seed.use = 0)
# Myeloid_merge <- RunTSNE(Myeloid_merge, dims = 1:10, seed.use = 0)
# Myeloid_merge <- FindNeighbors(Myeloid_merge, dims = 1:12)
# Myeloid_merge <- FindClusters(Myeloid_merge, resolution = 0.1)
# DimPlot(Myeloid_merge, reduction = "umap", label = T, pt.size = 1)
# DimPlot(Myeloid_merge, reduction = "umap", group.by = "orig.ident", label = T, pt.size = 1)
# DimPlot(Myeloid_merge, reduction = "tsne", label = T, pt.size = 1)
# 
# Myeloid_count <- as.matrix(Myeloid_merge@assays$RNA@counts)
# Myeloid_count <- unname(unlist(apply(Myeloid_count, 2, sum)))
# hist(Myeloid_count, breaks = 1000)
# #Not norm in raw counts
# 
# Markers_Myeloid_m_wilcox <- FindAllMarkers(Myeloid_merge, test.use = 'wilcox')
# Markers_Myeloid_m_ttest <- FindAllMarkers(Myeloid_merge, test.use = 't')
# 
# #
# trans_df_markers <- function(df_input, n_top) {
#   cluster <- levels(df_input$cluster)
#   res <- rep(NA, n_top*length(cluster))
#   dim(res) <- c(n_top, length(cluster))
#   dimnames(res) <- list(1:n_top, paste0("Cluster", cluster))
#   for (i in 1:length(cluster)) {
#     cluster_name <- cluster[i]
#     cluster_index <- which(df_input$cluster == cluster_name)[1:n_top]
#     cluster_gene <- df_input$gene[cluster_index] 
#     res[, i] <- cluster_gene
#   }
#   return(res)
# }
# 
# Markers_T_m_wilcox <- trans_df_markers(Markers_T_m_wilcox, 100)
# Markers_T_m_ttest <- trans_df_markers(Markers_T_m_ttest, 100)
# Markers_NK_m_wilcox <- trans_df_markers(Markers_NK_m_wilcox, 100)
# Markers_NK_m_ttest <- trans_df_markers(Markers_NK_m_ttest, 100)
# Markers_Myeloid_m_wilcox <- trans_df_markers(Markers_Myeloid_m_wilcox, 100)
# Markers_Myeloid_m_ttest <- trans_df_markers(Markers_Myeloid_m_ttest, 100)
# 
# write.table(Markers_T_m_wilcox, file = "Markers_T_wilcox.csv", sep = ",", row.names = T, col.names = T, quote = F)
# write.table(Markers_T_m_ttest, file = "Markers_T_ttest.csv", sep = ",", row.names = T, col.names = T, quote = F)
# write.table(Markers_NK_m_wilcox, file = "Markers_NK_wilcox.csv", sep = ",", row.names = T, col.names = T, quote = F)
# write.table(Markers_NK_m_ttest, file = "Markers_NK_ttest.csv", sep = ",", row.names = T, col.names = T, quote = F)
# write.table(Markers_Myeloid_m_wilcox, file = "Markers_Myeloid_wilcox.csv", sep = ",", row.names = T, col.names = T, quote = F)
# write.table(Markers_Myeloid_m_ttest, file = "Markers_Myeloid_ttest.csv", sep = ",", row.names = T, col.names = T, quote = F)
# ##
# 
# ##
# DimPlot(Myeloid_merge, reduction = "umap", label = T, pt.size = 1) + labs(title = "Myeloid_merge") + DimPlot(Myeloid_merge, reduction = "umap", group.by = "orig.ident", pt.size = 1)
# DimPlot(T_merge, reduction = "umap", label = T, pt.size = 1) + labs(title = "T_merge") + DimPlot(T_merge, reduction = "umap", group.by = "orig.ident", pt.size = 1)
# DimPlot(NK_merge, reduction = "umap", label = T, pt.size = 1) + labs(title = "NK_merge") + DimPlot(NK_merge, reduction = "umap", group.by = "orig.ident", pt.size = 1)
# table(T_merge$seurat_clusters)
# table(NK_merge$seurat_clusters)
# table(Myeloid_merge$seurat_clusters)
# 
# T_number <- cluster_number(T_merge, T, T)
# NK_number <- cluster_number(NK_merge, T, T)
# Myeloid_number <- cluster_number(Myeloid_merge, T, T)
# 
# T_number
# NK_number
# Myeloid_number
# 
# FeaturePlot(Myeloid_merge, Markers_Myeloid_m_ttest[1:4,1], reduction = 'umap', pt.size = 1)
# Markers_T_m_ttest[1:20, ]
# Markers_NK_m_ttest[1:20, ]
# Markers_Myeloid_m_wilcox[1:20, ]

###T subclustering:
Cryo_merge <- readRDS("Cryo_merge_sct.rds")
levels(Cryo_merge$seurat_clusters) <- c(
  "T",
  "Myeloid",
  "NK",
  "Myeloid",
  "Myeloid",
  "B",
  "Proliferating"
)
Cryo_merge <- SetIdent(Cryo_merge, value = "seurat_clusters")
T_sub <- FindSubCluster(Cryo_merge, 
                        cluster = "T", 
                        subcluster.name = "T", 
                        graph.name = "SCT_snn",
                        resolution = 0.1)
T_sub <- SetIdent(T_sub, value = "T")
# DimPlot(T_sub, reduction = "umap", label = T)
# DotPlot(T_sub, features = c("Cd3e", "Cd4", "Il7r", "Cd8a")) + RotatedAxis()
# FeaturePlot(T_sub, features = c("Cd4", "Cd8a"), pt.size = 1)
T_sub <- FindSubCluster(T_sub,
                        cluster = "T_2", 
                        subcluster.name = "T_2", 
                        graph.name = "SCT_snn",
                        resolution = 0.2)
T_sub <- SetIdent(T_sub, value = "T_2")
# DimPlot(T_sub, reduction = "umap", label = T)
# DotPlot(T_sub, features = c("Cd3e", "Cd4", "Il7r", "Cd8a")) + RotatedAxis()
T_sub <- FindSubCluster(T_sub,
                        cluster = "Proliferating", 
                        subcluster.name = "Proliferating", 
                        graph.name = "SCT_snn",
                        resolution = 0.5)
T_sub <- SetIdent(T_sub, value = "Proliferating")
# DimPlot(T_sub, reduction = "umap", label = T)
# DotPlot(T_sub, features = c("Cd3e", "Cd4", "Il7r", "Cd8a")) + RotatedAxis()
# FeaturePlot(T_sub, features = "Cd4", reduction = "umap", pt.size = 1, label = T, label.size = 2)
# FeaturePlot(T_sub, features = "Il7r", reduction = "umap", pt.size = 1, label = T, label.size = 2)
# FeaturePlot(T_sub, features = "Cd8a", reduction = "umap", pt.size = 1, label = T, label.size = 2)
##
Proliferating_sub <- subset(SetIdent(T_sub, value = "seurat_clusters"), idents = "Proliferating")
Proliferating_sub <- SetIdent(Proliferating_sub, value = "Proliferating")

DimPlot(Proliferating_sub, reduction = "umap", group.by = "Proliferating", label = T)
DotPlot(Proliferating_sub, features = c("Cd3e", "Cd4", "Cd8a")) + RotatedAxis()

Proliferating_sub <- SCTransform(Proliferating_sub, assay = "RNA", new.assay.name = "SCT")
Proliferating_sub <- FindVariableFeatures(Proliferating_sub, selection.method = "vst", nfeatures = 2000)
Proliferating_sub <- RunPCA(Proliferating_sub)
Proliferating_sub <- RunUMAP(Proliferating_sub, dims = 1:10, seed.use = 1)
DimPlot(Proliferating_sub, reduction = "umap", label = T, pt.size = 2)
DotPlot(Proliferating_sub, features = c("Cd3e", "Cd4", "Cd8a")) + RotatedAxis()
FeaturePlot(Proliferating_sub, features = c("Cd3e", "Cd4", "Cd8a"), reduction = "umap", pt.size = 1, label = T, label.size = 2)

Pro_markers_3 <- FindMarkers(object = Proliferating_sub, ident.1 = "Proliferating_3")
Pro_markers_3_vector <- rownames(Pro_markers_3[Pro_markers_3$p_val_adj<0.05, ])[order(Pro_markers_3$avg_log2FC[Pro_markers_3$p_val_adj<0.05], decreasing = T)]
FeaturePlot(Proliferating_sub, features = Pro_markers_3_vector[1:6])
#NK: Klra7 Klra4 Gzma
##
levels(T_sub@active.ident) <- c(
  "Myeloid",
  "NK",
  "T_CD4",
  "T_CD8",
  "Proliferating_T_CD4",
  "B",
  "T_CD8",
  "Proliferating_T_CD8",
  "T_CD4",
  "Proliferating_NK",
  "Proliferating_T_CD8"
)
DimPlot(T_sub, reduction = "umap", label = T)
DotPlot(T_sub, features = c("Cd3e", "Cd4", "Cd8a")) + RotatedAxis()
DotPlot(T_sub, features = marker_gene) + RotatedAxis()
#saveRDS(T_sub, file = "T_sub.rds")

T_sub_number <- cluster_number(T_sub, Percentage = T, Round = T)
T_sub_number

##
T_sub <- readRDS("T_sub.rds")
##

#subclustering CD4+ and CD+8 T_cells:
#
T_CD4 <- subset(T_sub, idents = c("T_CD4", "Proliferating_T_CD4"))
T_CD8 <- subset(T_sub, idents = c("T_CD8", "Proliferating_T_CD8"))
T_CD4 <- SCTransform(T_CD4, assay = "RNA", new.assay.name = "SCT")
T_CD4 <- FindVariableFeatures(T_CD4, selection.method = "vst", nfeatures = 2000)
T_CD4 <- RunPCA(T_CD4)
T_CD4 <- RunUMAP(T_CD4, dims = 1:10, seed.use = 0)
T_CD4 <- RunTSNE(T_CD4, dims = 1:10, seed.use = 0)
T_CD4 <- FindNeighbors(T_CD4, dims = 1:12)
T_CD8 <- SCTransform(T_CD8, assay = "RNA", new.assay.name = "SCT")
T_CD8 <- FindVariableFeatures(T_CD8, selection.method = "vst", nfeatures = 2000)
T_CD8 <- RunPCA(T_CD8)
T_CD8 <- RunUMAP(T_CD8, dims = 1:10, seed.use = 0)
T_CD8 <- RunTSNE(T_CD8, dims = 1:10, seed.use = 0)
T_CD8 <- FindNeighbors(T_CD8, dims = 1:12)

T_CD4 <- FindClusters(T_CD4, resolution = 0.8)
T_CD8 <- FindClusters(T_CD8, resolution = 0.8)

# library(clustree)
# T_CD4 <- FindClusters(T_CD4, resolution = seq(2,15,2)/10)
# T_CD8 <- FindClusters(T_CD8, resolution = seq(2,15,2)/10)
# clustree(T_CD4)
# clustree(T_CD8)
# #r=0.8 may be good

T_CD4 <- SetIdent(T_CD4, value = "SCT_snn_res.0.8")
T_CD8 <- SetIdent(T_CD8, value = "SCT_snn_res.0.8")
# DimPlot(T_CD4, reduction = "umap", label = T, label.size = 5, pt.size = 2) + labs(title = "T_CD4")
# DimPlot(T_CD4, reduction = "umap", label = T, label.size = 5, pt.size = 2, group.by = "orig.ident") + labs(title = "T_CD4")
# DimPlot(T_CD8, reduction = "umap", label = T, label.size = 5, pt.size = 2) + labs(title = "T_CD8")
# DimPlot(T_CD8, reduction = "umap", label = T, label.size = 5, pt.size = 2, group.by = "orig.ident") + labs(title = "T_CD8")
# FeaturePlot(T_CD4, features = "Pdcd1")
# FeaturePlot(T_CD8, features = "Mki67")
# DotPlot(T_CD4, features = "Pdcd1") + RotatedAxis()

#exhausted cells£ºPD1:Pdcd1/CTLA4:Ctla4/TIM-3:Havcr2
#T_CD4:cluster1/6/9
Exhausted_markers <- c("Pdcd1", "Ctla4", "Havcr2")
FeaturePlot(T_CD4, Exhausted_markers)
DotPlot(T_CD4, features = Exhausted_markers, group.by = "orig.ident") + RotatedAxis()
DotPlot(T_CD4, features = Exhausted_markers) + RotatedAxis()
FeaturePlot(T_CD8, Exhausted_markers)
DotPlot(T_CD8, features = Exhausted_markers, group.by = "orig.ident") + RotatedAxis()

#Treg:Foxp3
#T_CD4:cluster2
Treg_marker <- "Foxp3"
FeaturePlot(T_CD4, c("Foxp3","Mki67"))
DotPlot(T_CD4, features = Treg_marker) + RotatedAxis()

#naive: CD62L: Sell

# Markers_T_CD4 <- FindAllMarkers(T_CD4)
# Markers_T_CD8 <- FindAllMarkers(T_CD8)
# # save(Markers_T_CD4, file = "Markers_T_CD4_raw.rda")
# # save(Markers_T_CD8, file = "Markers_T_CD8_raw.rda")
# # Markers_T_CD4 <- load("Markers_T_CD4_raw.rda")
# # Markers_T_CD8 <- load("Markers_T_CD8_raw.rda")
# Markers_T_CD4_df <- trans_df_markers(Markers_T_CD4, 100)
# Markers_T_CD8_df <- trans_df_markers(Markers_T_CD8, 100)

T_CD4_number <- cluster_number(T_CD4, Percentage = T, Round = T)
T_CD8_number <- cluster_number(T_CD8, Percentage = T, Round = T)
T_CD4_number
T_CD8_number

# write.csv(Markers_T_CD4_df, file = "Markers_T_CD4.csv")
# write.csv(Markers_T_CD8_df, file = "Markers_T_CD8.csv")
# pdf("T_CD4_CD8.pdf")
# DimPlot(T_CD4, reduction = "umap", label = T, label.size = 5, pt.size = 1) + labs(title = "T_CD4")
# DimPlot(T_CD4, reduction = "umap", label = T, label.size = 5, pt.size = 1, group.by = "orig.ident") + labs(title = "T_CD4")
# FeaturePlot(T_CD4, features = "Mki67")
# DimPlot(T_CD8, reduction = "umap", label = T, label.size = 5, pt.size = 1) + labs(title = "T_CD8")
# DimPlot(T_CD8, reduction = "umap", label = T, label.size = 5, pt.size = 1, group.by = "orig.ident") + labs(title = "T_CD8")
# FeaturePlot(T_CD8, features = "Mki67")
# dev.off()

#####NK:
NK_sub <- subset(T_sub, idents = c("NK","Proliferating_NK"))
# DimPlot(NK_sub)
NK_sub <- SCTransform(NK_sub, assay = "RNA", new.assay.name = "SCT")
NK_sub <- FindVariableFeatures(NK_sub, selection.method = "vst", nfeatures = 2000)
NK_sub <- RunPCA(NK_sub)
NK_sub <- RunUMAP(NK_sub, dims = 1:10, seed.use = 0)
NK_sub <- FindNeighbors(NK_sub, dims = 1:12)

# NK_sub <- FindClusters(NK_sub, resolution = seq(2,15,2)/10)
# library(clustree)
# clustree(NK_sub)

NK_sub <- FindClusters(NK_sub, resolution = 0.6)
# DimPlot(NK_sub, label = T, label.size = 5)
# DimPlot(NK_sub, group.by = "orig.ident")
# FeaturePlot(NK_sub, "Mki67")
NK_sub_number <- cluster_number(NK_sub, Percentage = T, Round = T)
NK_sub_number

# NK_markers <- FindAllMarkers(NK_sub)
# NK_markers_df <- trans_df_markers(NK_markers, 50)
# # save(NK_markers, file = "NK_markers.rda")
# write.csv(NK_markers_df, file = "NK_markers_df.csv")

#####Myeloid:
Myeloid_sub <- subset(T_sub, idents = c("Myeloid"))
# DimPlot(Myeloid_sub)
Myeloid_sub <- SCTransform(Myeloid_sub, assay = "RNA", new.assay.name = "SCT")
Myeloid_sub <- FindVariableFeatures(Myeloid_sub, selection.method = "vst", nfeatures = 2000)
Myeloid_sub <- RunPCA(Myeloid_sub)
Myeloid_sub <- RunUMAP(Myeloid_sub, dims = 1:10, seed.use = 0)
Myeloid_sub <- FindNeighbors(Myeloid_sub, dims = 1:12)

# Myeloid_sub <- FindClusters(Myeloid_sub, resolution = seq(2,15,2)/10)
# library(clustree)
# clustree(Myeloid_sub)

Myeloid_sub <- FindClusters(Myeloid_sub, resolution = 0.6)
# DimPlot(Myeloid_sub, label = T, label.size = 5)
# DimPlot(Myeloid_sub, group.by = "orig.ident")
# FeaturePlot(Myeloid_sub, "Mki67")
Myeloid_sub_number <- cluster_number(Myeloid_sub, Percentage = T, Round = T)
Myeloid_sub_number

# Myeloid_markers <- FindAllMarkers(Myeloid_sub)
# Myeloid_markers_df <- trans_df_markers(Myeloid_markers, 50)
# # save(Myeloid_markers, file = "Myeloid_markers.rda")
# write.csv(Myeloid_markers_df, file = "Myeloid_markers_df.csv")

##Neutrophil:Ly-6G
FeaturePlot(Myeloid_sub, "Ly6g", pt.size = 2)
DotPlot(Myeloid_sub, features = "Ly6g")
#M2: CD206:Mrc1
FeaturePlot(Myeloid_sub, "Mrc1", pt.size = 2)
DotPlot(Myeloid_sub, features = "Mrc1")
#M1: 	Ccr7
FeaturePlot(Myeloid_sub, "Ccr7", pt.size = 2)
DotPlot(Myeloid_sub, features = "Ccr7")

#MDSC-M: Ly-6C : non-significant
FeaturePlot(Myeloid_sub, "Ly6c1", pt.size = 2)
DotPlot(Myeloid_sub, features = "Ly6c1")

# saveRDS(T_CD4, file = "T_CD4_sub.rds")
# saveRDS(T_CD8, file = "T_CD8_sub.rds")
# saveRDS(NK_sub, file = "NK_sub.rds")
# saveRDS(Myeloid_sub, file = "Myeloid_sub.rds")

##
T_CD4 <- readRDS("T_CD4_sub.rds")
T_CD8 <- readRDS("T_CD8_sub.rds")
NK_sub <- readRDS("NK_sub.rds")
Myeloid_sub <- readRDS("Myeloid_sub.rds")
Markers_T_CD4_df <- read.csv(file = "Markers_T_CD4.csv", row.names = 1)
Markers_T_CD8_df <- read.csv(file = "Markers_T_CD8.csv", row.names = 1)
NK_markers_df <- read.csv(file = "NK_markers_df.csv", row.names = 1)
Myeloid_markers_df <- read.csv(file = "Myeloid_markers_df.csv", row.names = 1)
T_CD4_number <- cluster_number(T_CD4, Percentage = T, Round = T)
T_CD8_number <- cluster_number(T_CD8, Percentage = T, Round = T)
NK_sub_number <- cluster_number(NK_sub, Percentage = T, Round = T)
Myeloid_sub_number <- cluster_number(Myeloid_sub, Percentage = T, Round = T)
T_CD4_number
T_CD8_number
NK_sub_number
Myeloid_sub_number

## GO / KEGG annotation
library(clusterProfiler)
library(enrichplot)
options(connectionObserver = NULL)
OrgDb = get("org.Mm.eg.db") # mouse database

GO_cluster <- function(marker_df, cluster, ntop = 50, ont = "ALL", type = 'dot', pvalue = 0.05, qvalue = 0.2, show = 20, title = NULL, font.size = 7) {
  message('ont can be one of "BP", "MF", "CC" and "ALL"(default)')
  message('type can be some of "bar", "dot", "cnet" and "emap"')
  GeneSet <- marker_df[,cluster][seq(ntop)]
  title <- ifelse(is.null(title), cluster, title)
  genes <- bitr(GeneSet, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  ego <- enrichGO(gene = genes$ENTREZID, ont = ont, OrgDb = OrgDb, pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (nrow(ego) != 0) {
    if ('bar' %in% type) print(barplot(ego, showCategory = show, title = title, font.size = font.size))
    if ('dot' %in% type) print(dotplot(ego, showCategory = show, title = title, font.size = font.size))
    if ('cnet' %in% type) print(cnetplot(ego, showCategory = show))
    if ('emap' %in% type) {
      ego1 <- pairwise_termsim(ego, method = "JC", semData = NULL, showCategory = show)
      emapplot(ego1, showCategory = show)
      }
    }
}

GO_sim <- function(Geneset, ont = "ALL", type = 'dot', pvalue = 0.05, qvalue = 0.2, show = 20, font.size = 7) {
  message('ont can be one of "BP", "MF", "CC" and "ALL"(default)')
  message('type can be some of "bar", "dot", "cnet" and "emap"')
  genes <- bitr(Geneset, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  ego <- enrichGO(gene = genes$ENTREZID, ont = ont, OrgDb = OrgDb, pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (nrow(ego) != 0) {
    if ('bar' %in% type) print(barplot(ego, showCategory = show, font.size = font.size))
    if ('dot' %in% type) print(dotplot(ego, showCategory = show, font.size = font.size))
    if ('cnet' %in% type) print(cnetplot(ego, showCategory = show))
    if ('emap' %in% type) {
      ego1 <- pairwise_termsim(ego, method = "JC", semData = NULL, showCategory = show)
      emapplot(ego1, showCategory = show)
      }
    }
}

KEGG_cluster <- function(marker_df, cluster, ntop = 50, type = c('bar', 'dot'), pvalue = 0.05, qvalue = 0.2, show = 20, title = NULL, font.size = 7) {
  GeneSet <- marker_df[,cluster][seq(ntop)]
  title <- ifelse(is.null(title), cluster, title)
  genes <- bitr(GeneSet, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  kk <- enrichKEGG(gene = genes$ENTREZID, organism  = 'mmu', pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (nrow(kk) != 0) {
    if ('bar' %in% type) print(barplot(kk, showCategory = show, title = title, font.size = font.size))
    if ('dot' %in% type) print(dotplot(kk, showCategory = show, title = title, font.size = font.size))
  }
}

KEGG_sim <- function(Geneset, type = c('bar', 'dot'), pvalue = 0.05, qvalue = 0.2, show = 20, font.size = 7) {
  genes <- bitr(Geneset, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  kk <- enrichKEGG(gene = genes$ENTREZID, organism  = 'mmu', pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (nrow(kk) != 0) {
    if ('bar' %in% type) print(barplot(kk, showCategory = show, font.size = font.size))
    if ('dot' %in% type) print(dotplot(kk, showCategory = show, font.size = font.size))
  }
}

#T_CD4:2/
T_CD4_diff_cluster <- rownames(T_CD4_number)[-nrow(T_CD4_number)][T_CD4_number[-nrow(T_CD4_number),"trend"]>1]
T_CD8_diff_cluster <- rownames(T_CD8_number)[-nrow(T_CD8_number)][T_CD8_number[-nrow(T_CD8_number),"trend"]>1]
NK_diff_cluster <- rownames(NK_sub_number)[-nrow(NK_sub_number)][NK_sub_number[-nrow(NK_sub_number),"trend"]>1]
Myeloid_diff_cluster <- rownames(Myeloid_sub_number)[-nrow(Myeloid_sub_number)][Myeloid_sub_number[-nrow(Myeloid_sub_number),"trend"]>1]

##
ont = "BP"
show = 20
type = "dot"
ntop = 100
pdf("T_CD4_annotation_BP.pdf")
for (i in T_CD4_diff_cluster) {
  i <- paste0("Cluster", i)
  GO_cluster(Markers_T_CD4_df, i, ntop = ntop, ont = ont, type = type, show = show, title = paste(i, "GO", ont))
  KEGG_cluster(Markers_T_CD4_df, i, ntop = ntop, type = type, show = show, title = paste(i, "KEGG"))
}
dev.off()
pdf("T_CD8_annotation_BP.pdf")
for (i in T_CD8_diff_cluster) {
  i <- paste0("Cluster", i)
  GO_cluster(Markers_T_CD8_df, i, ntop = ntop, ont = ont, type = type, show = show, title = paste(i, "GO", ont))
  KEGG_cluster(Markers_T_CD8_df, i, ntop = ntop, type = type, show = show, title = paste(i, "KEGG"))
}
dev.off()
pdf("NK_annotation_BP.pdf")
for (i in NK_diff_cluster) {
  i <- paste0("Cluster", i)
  GO_cluster(NK_markers_df, i, ntop = ntop, ont = ont, type = type, show = show, title = paste(i, "GO", ont))
  KEGG_cluster(NK_markers_df, i, ntop = ntop, type = type, show = show, title = paste(i, "KEGG"))
}
dev.off()
pdf("Myeloid_annotation_BP.pdf")
for (i in Myeloid_diff_cluster) {
  i <- paste0("Cluster", i)
  GO_cluster(Myeloid_markers_df, i, ntop = ntop, ont = ont, type = type, show = show, title = paste(i, "GO", ont))
  KEGG_cluster(Myeloid_markers_df, i, ntop = ntop, type = type, show = show, title = paste(i, "KEGG"))
}
dev.off()

# 
# GO_cluster(Markers_T_CD4_df, 9, ntop = 100, ont = "BP", type = "dot", show = 50, title = paste(9, "GO"))
# KEGG_cluster(Markers_T_CD4_df, 9, ntop = 100, type = "dot", show = 50, title = paste(9, "KEGG"))
# 
# GO_cluster(Markers_T_CD8_df, 4, ntop = 100, type = "dot", show = 50, title = paste(4, "GO"))
# KEGG_cluster(Markers_T_CD8_df, 4, ntop = 100, type = "dot", show = 50, title = paste(4, "KEGG"))
# GO_cluster(Markers_T_CD8_df, 2, ntop = 100, type = "dot", show = 50, title = paste(2, "GO"))
# KEGG_cluster(Markers_T_CD8_df, 2, ntop = 100, type = "dot", show = 50, title = paste(2, "KEGG"))
# 
# GO_cluster(NK_markers_df, 1, ntop = 100, type = "dot", show = 50, title = paste(1, "GO"))
# KEGG_cluster(NK_markers_df, 1, ntop = 100, type = "dot", show = 50, title = paste(1, "KEGG"))
# 
# GO_cluster(Myeloid_markers_df, 0, ntop = 100, type = "dot", show = 50, title = paste(0, "GO"))
# KEGG_cluster(Myeloid_markers_df, 0, ntop = 100, type = "dot", show = 50, title = paste(0, "KEGG"))
# GO_cluster(Myeloid_markers_df, 2, ntop = 100, type = "dot", show = 50, title = paste(2, "GO"))
# KEGG_cluster(Myeloid_markers_df, 2, ntop = 100, type = "dot", show = 50, title = paste(2, "KEGG"))


## GSEA annotation analysis:
my_fc <- function(data, 
                  meta, 
                  LOG = T,
                  average_method = c("arithmetic", "geometric")[1]
) {
  if (class(data)=="dgCMatrix")
    data <- as.matrix(data)
  N_cluster <- length(unique(meta))
  meta_levels <- sort(unique(meta))
  fc <- apply(data, 1,
              function(x, Log = LOG) {
                fc <- c()
                for (i in 1:N_cluster) {
                  index <- meta == meta_levels[i]
                  if (average_method == "arithmetic") {
                    target <- mean(x[index])
                    non_target <- mean(x[!index])
                  }
                  if (average_method == "geometric") {
                    target <- expm1(mean(log1p(x[index])))
                    non_target <- expm1(mean(log1p(x[!index])))
                  }
                  fc_temp <- target / non_target
                  fc[i] <- ifelse(test = Log, yes = log(fc_temp), no = fc_temp)
                }
                return(fc)
              }
              )
  rownames(fc) <- meta_levels
  res <- t(fc)
  return(res)
}

GSEA_anno <- function(FC, symbols, anno = "GO", ont = "All", pAdjustMethod = "BH", show = 20, title = NULL, font.size = 12) {
  message('The parameters first and second "FC" and "symbols" are respectively the input FoldChange vector and corresponding gene symbols')
  message('pAdjustMethod can be one of "holm", "hochberg", "hommel", "bonferroni", "BH"(default), "BY", "fdr" and "none"')
  message('Annotation method can be some of "GO"(default), "KEGG" and "WP"')
  message('If GO_annotation, ont can be one of "BP", "MF", "CC" and "ALL"(default)')
  names(FC) <- symbols
  genes <- symbols
  genes <- bitr(genes, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  FC <- FC[names(FC) %in% genes$SYMBOL]
  names(FC) <- as.character(genes$ENTREZID)
  FC <- FC[!is.infinite(FC)]
  FC <- sort(FC, decreasing = T)
  if ("GO" %in% anno) {
    Gse_GO <- gseGO(FC, ont = ont, OrgDb = OrgDb, pAdjustMethod = pAdjustMethod)
    if (nrow(Gse_GO@result) != 0) {
      Dot_Gse_GO <- dotplot(Gse_GO, split=".sign", showCategory = show, title = paste("GSE_GO", title), font.size = font.size) + facet_wrap(~.sign, scales = "free")
      print(Dot_Gse_GO)
    }
  }
  if ("KEGG" %in% anno) {
    Gse_KEGG <- gseKEGG(FC, organism = "mmu", pAdjustMethod = pAdjustMethod)
    if (nrow(Gse_KEGG@result) != 0) {
      Dot_Gse_KEGG <- dotplot(Gse_KEGG, split=".sign", showCategory = show, title = paste("GSE_KEGG", title), font.size = font.size) + facet_wrap(~.sign, scales = "free")
      print(Dot_Gse_KEGG)
    }
  }
  if ("WP" %in% anno) {
    Gse_WikiPathway <- gseWP(FC, organism = "Mus musculus")
    if (nrow(Gse_WikiPathway@result) != 0) {
      Dot_Gse_WikiPathway <- dotplot(Gse_WikiPathway, showCategory = show, title = paste("GSE_WikiPathway", title), font.size= font.size)
      print(Dot_Gse_WikiPathway)
    }
  }
}

T_CD4_fc <- my_fc(data = T_CD4@assays$SCT@data, 
                  meta = T_CD4@active.ident, 
                  LOG = F,
                  average_method = "geometric")
T_CD8_fc <- my_fc(data = T_CD8@assays$SCT@data, 
                  meta = T_CD8@active.ident, 
                  LOG = F,
                  average_method = "geometric")
NK_fc <- my_fc(data = NK_sub@assays$SCT@data, 
               meta = NK_sub@active.ident, 
               LOG = F,
               average_method = "geometric")
Myeloid_fc <- my_fc(data = Myeloid_sub@assays$SCT@data, 
                    meta = Myeloid_sub@active.ident, 
                    LOG = F,
                    average_method = "geometric")

pdf("GSEA_annotation.pdf")
GSEA_anno(T_CD4_fc[,'9'], rownames(T_CD4_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "T_CD4_cluster9")
GSEA_anno(T_CD8_fc[,'2'], rownames(T_CD8_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "T_CD8_cluster2")
GSEA_anno(T_CD8_fc[,'4'], rownames(T_CD8_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "T_CD8_cluster4")
GSEA_anno(NK_fc[,'1'], rownames(NK_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "NK_cluster1")
GSEA_anno(NK_fc[,'7'], rownames(NK_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "NK_cluster7")
GSEA_anno(Myeloid_fc[,'0'], rownames(Myeloid_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Myeloid_cluster0")
GSEA_anno(Myeloid_fc[,'2'], rownames(Myeloid_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Myeloid_cluster2", font.size = 3)
dev.off()


############
#### some code from that of the second seq
VlnPlot(T_CD8, features = unlist(list(
  Activated = c("Cd69", "Il2ra"),
  Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
  Naive = c("Sell"),
  Effector_Memory = c("Il7r"),
  Effector = c("Cd44"),
  Exhausted = c("Pdcd1", "Ctla4", "Lag3", "Havcr2")
)), group.by = 'orig.ident')

DotPlot(T_CD8, features = list(
  Activated = c("Cd69", "Il2ra"),
  Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
  Naive = c("Sell"),
  Effector_Memory = c("Il7r"),
  Effector = c("Cd44"),
  Exhausted = c("Pdcd1", "Ctla4", "Lag3", "Havcr2")
), group.by = 'orig.ident')

VlnPlot(T_CD4, features = unlist(list(
  Th1 = 'Tbx21',
  Th2 = 'Gata3',
  Th17 = 'Il17a',
  Treg = c('Foxp3', 'Il2ra')
)), group.by = 'orig.ident')

DotPlot(T_CD4, features = list(
  Th1 = 'Tbx21',
  Th2 = 'Gata3',
  Th17 = 'Il17a',
  Treg = c('Foxp3', 'Il2ra')
), group.by = 'orig.ident')


#CD8
T_CD8_temp <- SetIdent(T_CD8, value = 'orig.ident')
Markers_T_CD8_list <- FindAllMarkers(T_CD8_temp, only.pos = T)
Markers_T_CD8_df <- trans_df_markers(Markers_T_CD8_list, 120)

VlnPlot(T_CD8_temp, features = Markers_T_CD8_df[1:30, 1])
DotPlot(T_CD8_temp, features = Markers_T_CD8_df[1:30, 1]) + RotatedAxis()

T_CD8_GOBP <- my_GO(Markers_T_CD8_df[1:120, 1], ont = "BP", Simplify = T, return_res = T, font.size = 15)
T_CD8_GOMF <- my_GO(Markers_T_CD8_df[1:120, 1], ont = "MF", Simplify = T, return_res = T, font.size = 15)
T_CD8_GOCC <- my_GO(Markers_T_CD8_df[1:120, 1], ont = "CC", Simplify = T, return_res = T, font.size = 15)
T_CD8_GOBP@result[which(T_CD8_GOBP@result$Description %in% grep("interferon", T_CD8_GOBP@result$Description, value = T)), ]
IFN_genes <- unique(unlist(strsplit(T_CD8_GOBP@result[which(T_CD8_GOBP@result$Description %in% grep("interferon", T_CD8_GOBP@result$Description, value = T)), ]$geneID, split="/")))
IFN_genes <- bitr(IFN_genes, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
FeaturePlot(T_CD8_temp, features = IFN_genes, pt.size = 1)
VlnPlot(T_CD8_temp, features = IFN_genes)

T_CD8_fc <- my_fc(T_CD8_temp, Feature = 'Variable', average_method = 'geometric')

my_GSEA(T_CD8_fc[,1], rownames(T_CD8_fc), c("GO", "KEGG"), "BP", "BH", 200, "Cryo", Simplify = T)

#CD4
T_CD4_temp <- SetIdent(T_CD4, value = 'orig.ident')
Markers_T_CD4_list <- FindAllMarkers(T_CD4_temp, only.pos = T)
Markers_T_CD4_df <- trans_df_markers(Markers_T_CD4_list, 120)

VlnPlot(T_CD4_temp, features = Markers_T_CD4_df[1:30, 1])
DotPlot(T_CD4_temp, features = Markers_T_CD4_df[1:30, 1]) + RotatedAxis()

T_CD4_GOBP <- my_GO(Markers_T_CD4_df[1:120, 1], ont = "BP", Simplify = T, return_res = T, font.size = 15)
T_CD4_GOMF <- my_GO(Markers_T_CD4_df[1:120, 1], ont = "MF", Simplify = T, return_res = T, font.size = 15)
T_CD4_GOCC <- my_GO(Markers_T_CD4_df[1:120, 1], ont = "CC", Simplify = T, return_res = T, font.size = 15)
T_CD4_GOBP@result[which(T_CD4_GOBP@result$Description %in% grep("interferon", T_CD4_GOBP@result$Description, value = T)), ]
IFN_genes <- unique(unlist(strsplit(T_CD4_GOBP@result[which(T_CD4_GOBP@result$Description %in% grep("interferon", T_CD4_GOBP@result$Description, value = T)), ]$geneID, split="/")))
IFN_genes <- bitr(IFN_genes, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
FeaturePlot(T_CD4_temp, features = IFN_genes, pt.size = 1)
VlnPlot(T_CD4_temp, features = IFN_genes)

T_CD4_fc <- my_fc(T_CD4_temp, Feature = 'Variable', average_method = 'geometric')

my_GSEA(T_CD4_fc[,1], rownames(T_CD4_fc), c("GO", "KEGG"), "BP", "BH", 200, "Cryo", Simplify = T)
