setwd("E://Cryo-TCR/server/auto/")
source("utilities.r")
###
# Integret by merge function
# or FindIntegrationAnchors with IntegrateData
###
###########################################import raw_data
setwd("D://Ji_Wangxue/")
Cryo_merge1 <- readRDS("T_sub.rds")
Cryo_merge1@meta.data <- Cryo_merge1@meta.data[,c(1:6,8,12)]
colnames(Cryo_merge1@meta.data)[7] <- 'Cryo_merge_1'
levels(Cryo_merge1@meta.data$Proliferating) <- c('B', 'Myeloid', 'NK', 'Proliferating_T_CD8', 'Proliferating_T_CD8', 'Proliferating_T_CD4', 'Proliferating_NK', 'T_CD8', 'T_CD4', 'T_CD4', 'T_CD8')
Cryo_merge1@meta.data$Proliferating <- as.factor(Cryo_merge1@meta.data$Proliferating)
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"] <- paste0("Cryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"])
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"] <- paste0("NonCryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"])
Cryo_merge1_meta.data <- Cryo_merge1@meta.data

Cryo_data <- Read10X_h5("Cryo_filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo_filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("NonCryo_", colnames(NonCryo_data))
Cryo <- CreateSeuratObject(counts = Cryo_data, 
                           project = "Cryo", 
                           min.cells = 0, 
                           min.features = 200)
NonCryo <- CreateSeuratObject(counts = NonCryo_data, 
                              project = "NonCryo", 
                              min.cells = 0, 
                              min.features = 200)
Cryo_merge1 <- merge(Cryo, NonCryo)
Cryo_merge1[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge1, pattern = "^mt-")
Cryo_merge1 <- subset(Cryo_merge1, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 10)
Cryo_merge1@meta.data <- Cryo_merge1_meta.data
Cryo_merge1@meta.data$orig.ident <- paste0("Old_", Cryo_merge1@meta.data$orig.ident)

setwd("E://Cryo-TCR/server/auto/")
Cryo_merge2_300 <- readRDS("scRNA_Cryo_2_300.rds")
Cryo_merge2_300@meta.data <- Cryo_merge2_300@meta.data[,c(1:7,18)]
colnames(Cryo_merge2_300@meta.data)[8]<-'Cluster_Cryo_merge2_300'
Cryo_merge2_300_meta.data <- Cryo_merge2_300@meta.data

Cryo_merge2_500 <- readRDS("scRNA_Cryo_2_500.rds")
Cryo_merge2_500@meta.data <- Cryo_merge2_500@meta.data[,c(1:7,18,20:24)]
colnames(Cryo_merge2_500@meta.data)[8]<-'Cluster_Cryo_merge2_500'
Cryo_merge2_500_meta.data <- Cryo_merge2_500@meta.data

Cryo_data <- Read10X_h5("Cryo/filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo/filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("NonCryo_", colnames(NonCryo_data))
Cryo <- CreateSeuratObject(counts = Cryo_data,
                           project = "Cryo",
                           min.cells = 0,
                           min.features = 200)
NonCryo <- CreateSeuratObject(counts = NonCryo_data,
                              project = "NonCryo",
                              min.cells = 0,
                              min.features = 200)
Cryo_merge2 <- merge(Cryo, NonCryo)
Cryo_merge2[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge2, pattern = "^mt-")

Cryo_merge2_300 <- subset(Cryo_merge2, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
Cryo_merge2_300@meta.data <- Cryo_merge2_300_meta.data
Cryo_merge2_300@meta.data$orig.ident <- paste0("New_", Cryo_merge2_300@meta.data$orig.ident)

Cryo_merge2_500 <- subset(Cryo_merge2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
Cryo_merge2_500 <- Cryo_merge2_500[, colnames(Cryo_merge2_500) %in% rownames(Cryo_merge2_500_meta.data)]
Cryo_merge2_500@meta.data <- Cryo_merge2_500_meta.data
Cryo_merge2_500@meta.data$orig.ident <- paste0("New_", Cryo_merge2_500@meta.data$orig.ident)

####################################
################# merge
Cryo_merge3_300 <- merge(Cryo_merge1, Cryo_merge2_300, add.cell.ids = c("Old","New"))
Cryo_merge3_500 <- merge(Cryo_merge1, Cryo_merge2_500, add.cell.ids = c("Old","New"))

Cryo_merge3_300 <- SCTransform(Cryo_merge3_300)
Cryo_merge3_300 <- FindVariableFeatures(Cryo_merge3_300, selection.method = "vst", nfeatures = 2000)
Cryo_merge3_300 <- RunPCA(Cryo_merge3_300)
Cryo_merge3_300 <- RunUMAP(Cryo_merge3_300, dims = 1:10)
Cryo_merge3_300 <- RunTSNE(Cryo_merge3_300, dims = 1:10)
Cryo_merge3_300 <- FindNeighbors(Cryo_merge3_300, dims = 1:12)
Cryo_merge3_300 <- FindClusters(Cryo_merge3_300, resolution = 0.04)

Cryo_merge3_500 <- SCTransform(Cryo_merge3_500)
Cryo_merge3_500 <- FindVariableFeatures(Cryo_merge3_500, selection.method = "vst", nfeatures = 2000)
Cryo_merge3_500 <- RunPCA(Cryo_merge3_500)
Cryo_merge3_500 <- RunUMAP(Cryo_merge3_500, dims = 1:10)
Cryo_merge3_500 <- RunTSNE(Cryo_merge3_500, dims = 1:10)
Cryo_merge3_500 <- FindNeighbors(Cryo_merge3_500, dims = 1:12)
Cryo_merge3_500 <- FindClusters(Cryo_merge3_500, resolution = 0.04)

###################

DimPlot(Cryo_merge3_300, reduction = "pca", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_300, reduction = "umap", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_300, reduction = "tsne", group.by = 'orig.ident', pt.size = 1)

DimPlot(Cryo_merge3_300, reduction = "pca", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_300, reduction = "umap", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_300, reduction = "tsne", label = T, label.size = 10, pt.size = 1)

DimPlot(Cryo_merge3_300, label = T, label.size = 5, pt.size = 1, group.by = 'Proliferating')
DimPlot(Cryo_merge3_300, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_1')
DimPlot(Cryo_merge3_300, label = T, label.size = 5, pt.size = 1, group.by = 'TCR')
DimPlot(Cryo_merge3_300, label = T, label.size = 5, pt.size = 1, group.by = 'Cluster_Cryo_merge2_300')

VlnPlot(Cryo_merge3_300, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Cryo_merge3_300, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')

Cryo2_cluster2 <- read.table(file = "temp_barcode.txt")$x
Cryo2_cluster2 <- colnames(Cryo_merge3_300) %in% Cryo2_cluster2
Cryo_merge3_300$cluster2 <- Cryo2_cluster2
table(Cryo_merge3_300@active.ident, Cryo2_cluster2)

################# 500
DimPlot(Cryo_merge3_500, reduction = "pca", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "umap", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "tsne", group.by = 'orig.ident', pt.size = 1)

DimPlot(Cryo_merge3_500, reduction = "pca", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "umap", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "tsne", label = T, label.size = 10, pt.size = 1)

DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Proliferating')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_1')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'TCR')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cluster_Cryo_merge2_500')

VlnPlot(Cryo_merge3_500, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Cryo_merge3_500, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')

Cryo2_cluster2 <- read.table(file = "temp_barcode.txt")$x
Cryo2_cluster2 <- colnames(Cryo_merge3_500) %in% Cryo2_cluster2
Cryo_merge3_500$cluster2 <- Cryo2_cluster2
table(Cryo_merge3_500@active.ident, Cryo2_cluster2)