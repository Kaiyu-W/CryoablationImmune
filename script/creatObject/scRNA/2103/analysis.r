setwd("D://Ji_Wangxue/")
library(Seurat)
library(ggplot2)
Cryo_data <- Read10X_h5("Cryo_filtered_feature_bc_matrix.h5")
NonCryo_data <- Read10X_h5("NonCryo_filtered_feature_bc_matrix.h5")
Cryo <- CreateSeuratObject(counts = Cryo_data, 
                           project = "Cryo", 
                           min.cells = 3, 
                           min.features = 200)
NonCryo <- CreateSeuratObject(counts = NonCryo_data, 
                              project = "NonCryo", 
                              min.cells = 3, 
                              min.features = 200)
Cryo_merge <- merge(Cryo, NonCryo)

Cryo
NonCryo
Cryo[["percent.mt"]] <- PercentageFeatureSet(Cryo, pattern = "^mt-")
NonCryo[["percent.mt"]] <- PercentageFeatureSet(NonCryo, pattern = "^mt-")

Cryo_merge[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge, pattern = "^mt-")

#visualization
vln_cryo <- VlnPlot(Cryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Cryo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cryo, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot3 <- FeatureScatter(Cryo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vln_cryo
plot1 + plot2 + plot3
vln_noncryo <- VlnPlot(NonCryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot4 <- FeatureScatter(NonCryo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(NonCryo, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(NonCryo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vln_noncryo
plot4 + plot5 + plot6

vln_cryo_m <- VlnPlot(Cryo_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1_m <- FeatureScatter(Cryo_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_m <- FeatureScatter(Cryo_merge, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot3_m <- FeatureScatter(Cryo_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vln_cryo_m
plot1_m + plot2_m + plot3_m

#subset
Cryo <- subset(Cryo, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 10)
NonCryo <- subset(NonCryo, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 15)
Cryo
NonCryo

Cryo_merge <- subset(Cryo_merge, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 10)
Cryo_merge

##
vln_cryo_1 <- VlnPlot(Cryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1_1 <- FeatureScatter(Cryo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(Cryo, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot3_1 <- FeatureScatter(Cryo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vln_cryo_1
plot1_1 + plot2_1 + plot3_1
vln_noncryo_1 <- VlnPlot(NonCryo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot4_1 <- FeatureScatter(NonCryo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5_1 <- FeatureScatter(NonCryo, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot6_1 <- FeatureScatter(NonCryo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
vln_noncryo_1
plot4_1 + plot5_1 + plot6_1
##

#Normolize
Cryo <- SCTransform(Cryo)
NonCryo <- SCTransform(NonCryo)
Cryo_merge <- SCTransform(Cryo_merge)
#or
#Cryo <- NormalizeData(Cryo, normalization.method = "LogNormalize", scale.factor = 10000)
#Cryo <- ScaleData(Cryo, assay = "RNA", features = rownames(Cryo))
#NonCryo <- NormalizeData(NonCryo, normalization.method = "LogNormalize", scale.factor = 10000)
#NonCryo <- ScaleData(NonCryo, assay = "RNA", features = rownames(NonCryo))

#Variable
Cryo <- FindVariableFeatures(Cryo, selection.method = "vst", nfeatures = 2000)
NonCryo <- FindVariableFeatures(NonCryo, selection.method = "vst", nfeatures = 2000)
Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = 2000)
top10_cryo <- head(VariableFeatures(Cryo), 10)
plot7 <- VariableFeaturePlot(Cryo)
plot8 <- LabelPoints(plot = plot7, points = top10_cryo, repel = TRUE)
plot7 + plot8
plot8
top10_cryo
top10_noncryo <- head(VariableFeatures(NonCryo), 10)
plot9 <- VariableFeaturePlot(NonCryo)
plot10 <- LabelPoints(plot = plot9, points = top10_noncryo, repel = TRUE)
plot9 + plot10
plot10
top10_noncryo

top10_cryo_m <- head(VariableFeatures(Cryo_merge), 10)
plot7_m <- VariableFeaturePlot(Cryo_merge)
plot8_m <- LabelPoints(plot = plot7_m, points = top10_cryo_m, repel = TRUE)
plot7_m + plot8_m
plot8_m
top10_cryo_m

#Dimension Reduction : PCA
Cryo <- RunPCA(Cryo, features = VariableFeatures(object = Cryo))
#Cryo <- RunPCA(Cryo, npcs = 100)
NonCryo <- RunPCA(NonCryo, features = VariableFeatures(object = NonCryo))
Cryo_merge <- RunPCA(Cryo_merge)

pca_cryo <- DimPlot(Cryo, reduction = "pca") + labs(title = "Cryo")
pca_noncryo <- DimPlot(NonCryo, reduction = "pca") + labs(title = "NonCryo")
pca_cryo_m <- DimPlot(Cryo_merge, reduction = "pca") + labs(title = "Cryo_merge")
#DimHeatmap(Cryo, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(NonCryo, dims = 1, cells = 500, balanced = TRUE)
#DimHeatmap(Cryo_merge, dims = 1, cells = 500, balanced = TRUE)
VizDimLoadings(Cryo, dims = 1:2, reduction = "pca")
VizDimLoadings(NonCryo, dims = 1:2, reduction = "pca")
VizDimLoadings(Cryo_merge, dims = 1:2, reduction = "pca")
###optional(not sct)
#Cryo <- JackStraw(Cryo, num.replicate = 100)
#Cryo <- ScoreJackStraw(Cryo, dims = 1:20)
#JackStrawPlot(Cryo, dims = 1:20)
#ElbowPlot(Cryo)
#NonCryo <- JackStraw(NonCryo, num.replicate = 100)
#NonCryo <- ScoreJackStraw(NonCryo, dims = 1:20)
#JackStrawPlot(NonCryo, dims = 1:20)
#ElbowPlot(NonCryo)
#Cryo_merge <- JackStraw(Cryo_merge, num.replicate = 100)
#Cryo_merge <- ScoreJackStraw(Cryo_merge, dims = 1:20)
#JackStrawPlot(Cryo_merge, dims = 1:20)
#ElbowPlot(Cryo_merge)

#Clustering
Cryo <- FindNeighbors(Cryo, dims = 1:12)
Cryo <- FindClusters(Cryo, resolution = 0.3)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
Cryo <- FindClusters(Cryo, resolution = 0.4)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
Cryo <- FindClusters(Cryo, resolution = 0.5)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
Cryo <- FindClusters(Cryo, resolution = 0.6)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
DimPlot(Cryo, reduction = "tsne", label = T, pt.size = 1)
table(Cryo@active.ident)
Cryo <- FindClusters(Cryo, resolution = 0.7)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
#r=0.6/0.7 return same results
Cryo <- FindClusters(Cryo, resolution = 0.8)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
Cryo <- FindClusters(Cryo, resolution = 0.9)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
DimPlot(Cryo, reduction = "tsne", label = T, pt.size = 1)
table(Cryo@active.ident)
Cryo <- FindClusters(Cryo, resolution = 1)
DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
#r=0.9/0.1 return same results


NonCryo <- FindNeighbors(NonCryo, dims = 1:12)
NonCryo <- FindClusters(NonCryo, resolution = 0.3)


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
Cryo <- RunUMAP(Cryo, dims = 1:10, seed.use = 0)
Cryo <- RunTSNE(Cryo, dims = 1:10, seed.use = 0)
NonCryo <- RunUMAP(NonCryo, dims = 1:10, seed.use = 0)
NonCryo <- RunTSNE(NonCryo, dims = 1:10, seed.use = 0)
Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10, seed.use = 0)
Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10, seed.use = 0)
#Save RData
#save(Cryo, NonCryo, file="Cryo_nonCryo_sct.rda")
#save(Cryo, NonCryo, file="Cryo_nonCryo_vst.rda")
#saveRDS(Cryo_merge, file="Cryo_merge_sct.rds")
#load('Cryo_nonCryo_xxx.rda')

#Visualization
pca_cryo_1 <- DimPlot(Cryo, reduction = "pca")
pca_noncryo_1 <- DimPlot(NonCryo, reduction = "pca")
umap_cryo <- DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1)
tsne_cryo <- DimPlot(Cryo, reduction = "tsne", label = T, pt.size = 1)
umap_noncryo <- DimPlot(NonCryo, reduction = "umap", label = T, pt.size = 1)
tsne_noncryo <- DimPlot(NonCryo, reduction = "tsne", label = T, pt.size = 1)

umap_cryo_1 <- DimPlot(Cryo, reduction = "umap", label = T, pt.size = 1, repel = F)
tsne_cryo_1 <- DimPlot(Cryo, reduction = "tsne", label = T, pt.size = 1)
umap_noncryo_1 <- DimPlot(NonCryo, reduction = "umap", label = T, pt.size = 1)
tsne_noncryo_1 <- DimPlot(NonCryo, reduction = "tsne", label = T, pt.size = 1)
markers <- FindAllMarkers(Cryo, test.use = "t")

DimPlot(Cryo_merge, reduction = "pca")
DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
DimPlot(Cryo_merge, reduction = "tsne", label = T, pt.size = 1)

pca_cryo_1
pca_noncryo_1
umap_cryo
tsne_cryo
umap_noncryo
tsne_noncryo


#Marker genes
Markers_cryo <- FindAllMarkers(Cryo,test.use = 't')
Markers_noncryo <- FindAllMarkers(NonCryo,test.use = 't')
Markers_cryo_m <- FindAllMarkers(Cryo_merge,test.use = 't')

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

Markers_cryo_df <- trans_df_markers(Markers_cryo, 1000)
Markers_noncryo_df <- trans_df_markers(Markers_noncryo, 100)
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
markers_vector <- character()
for (i in names(Markers)) markers_vector <- append(markers_vector, Markers[[i]])
markers_vector <- unique(markers_vector)
markers_vector <- markers_vector[markers_vector != ""]
markers_rank <- matrix(nrow = length(markers_vector), 
                       ncol = ncol(Markers_cryo_df), 
                       dimnames = list(
                         markers_vector, 
                         colnames(Markers_cryo_df)))
for (i in 1:ncol(Markers_cryo_df)) {
  index <- sapply(markers_vector, function(x) which(Markers_cryo_df[, i] %in% x))
  index <- unlist(index)
  markers_rank[rownames(markers_rank) %in% names(index), i] <- unname(index)
  markers_rank[-which(rownames(markers_rank) %in% names(index)), i] <- 1000
}

heatmap(t(1000 - markers_rank), Rowv = NA, Colv = NA)

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
#Leukocytes:CD45
FeaturePlot(Cryo, Markers$Leukocytes, reduction = 'umap', cols=c("lightgrey", "#ff0000", "#00ff00"), pt.size = 1)
FeaturePlot(Cryo, Markers$Leukocytes, reduction = 'tsne', cols=c("lightgrey", "#ff0000", "#00ff00"), pt.size = 1)
#B_Cells:CD19
FeaturePlot(Cryo, Markers$B_Cells, reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, Markers$B_Cells, reduction = 'tsne', pt.size = 1)
#Plasma_Cells
for (i in Markers$Plasma_Cells) print(FeaturePlot(Cryo, i, reduction = 'umap', pt.size = 1))
for (i in Markers$Plasma_Cells) print(FeaturePlot(Cryo, i, reduction = 'tsne', pt.size = 1))
#NK_Cells
for (i in Markers$NK_Cells[1:2]) print(FeaturePlot(Cryo, i, reduction = 'umap', pt.size = 1))
for (i in Markers$NK_Cells[1:2]) print(FeaturePlot(Cryo, i, reduction = 'tsne', pt.size = 1))
FeaturePlot(Cryo, "Klrk1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd3e", cols=c("blue", "white"), reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd247", cols=c("blue", "white"), reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Klrk1", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "Cd3e", cols=c("blue", "white"), reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "Cd247", cols=c("blue", "white"), reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, c("Klrk1","Cd247"), pt.size = 1, reduction = 'umap', blend = T)
FeaturePlot(Cryo, c("Klrk1","Cd247"), pt.size = 1, reduction = 'tsne', blend = T)
#NKT_Cells
FeaturePlot(Cryo, "Klrb1c", pt.size = 1)
FeaturePlot(Cryo, "Cd3e", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd247", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd3e", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "Cd247", reduction = 'tsne', pt.size = 1)
#T_Cells
FeaturePlot(Cryo, "Cd3e", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd247", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd3e", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "Cd247", reduction = 'tsne', pt.size = 1)
  #CD4
  FeaturePlot(Cryo, "Cd4", reduction = 'umap', pt.size = 1)
  #Cd8
  FeaturePlot(Cryo, "Cd8a", reduction = 'umap', pt.size = 1)
#Myeloid_Cells
FeaturePlot(Cryo, Markers$Myeloid_Cells, reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, Markers$Myeloid_Cells, reduction = 'tsne', pt.size = 1)
#Monocytes
FeaturePlot(Cryo, "Cd14", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Cd14", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "Adgre1", cols=c("blue", "white"), reduction = 'umap', pt.size = 0.1)
FeaturePlot(Cryo, "Adgre1", cols=c("blue", "white"), reduction = 'tsne', pt.size = 0.1)
  #M1:Cd80,Cd86,Nos2
  FeaturePlot(Cryo, "Cd80", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo, "Cd86", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo, "Nos2", reduction = 'umap', pt.size = 1)
  #M2:Cd163,Mrc1,Arg1/Arg2
  FeaturePlot(Cryo, "Cd163", reduction = 'umap', pt.size = 1)
  FeaturePlot(Cryo, "Mrc1", reduction = 'umap', pt.size = 1)
#Macrophages
FeaturePlot(Cryo, "Adgre1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Adgre1", reduction = 'tsne', pt.size = 1)
#cDCs
FeaturePlot(Cryo, "Itgax", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "H2-Ab1", reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, "Itgax", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "H2-Ab1", reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, "H2-Ea", pt.size = 1)
#pDCs
FeaturePlot(Cryo, 'Siglech', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Bst2', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Siglech', reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, 'Bst2', reduction = 'tsne', pt.size = 1)
FeaturePlot(Cryo, Markers$pDCs[3], reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, Markers$pDCs[3], reduction = 'tsne', pt.size = 1)
#MDSCs
for (i in Markers$MDSCs) print(FeaturePlot(Cryo, i, reduction = 'umap', pt.size = 1))
for (i in Markers$MDSCs) print(FeaturePlot(Cryo, i, reduction = 'tsne', pt.size = 1))
#Granulocytes
FeaturePlot(Cryo, 'Ifng', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Il6', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Ccl2', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Mpo', reduction = 'umap', pt.size = 1)
FeaturePlot(Cryo, 'Tnf', reduction = 'umap', pt.size = 1)

###
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
FeaturePlot(Cryo_merge, Markers_cryo_merge_df[1:4,5], reduction = 'umap', pt.size = 1)


Markers_cryo_merge_df[1:10,]

