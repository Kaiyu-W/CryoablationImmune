setwd("E://Cryo-TCR/server/auto/")
source("utilities.r")

load(file = "C://celldex_reference/Database_form_celldex.rda") # reference data
###########################################1. batch between first and second one without CD45- cells
###
# Integret by merge function
# or FindIntegrationAnchors with IntegrateData
###
###########################################import raw_data
setwd("D://Ji_Wangxue/")
Cryo_merge1 <- readRDS("T_sub.rds")
Cryo_merge1@meta.data <- Cryo_merge1@meta.data[,c(1:6,8,12)]
colnames(Cryo_merge1@meta.data)[7] <- 'Cryo_merge_1'
Cryo_merge1@meta.data$Proliferating <- as.factor(Cryo_merge1@meta.data$Proliferating)
levels(Cryo_merge1@meta.data$Proliferating) <- c('B', 'Myeloid', 'NK', 'Proliferating_T_CD8', 'Proliferating_T_CD8', 'Proliferating_T_CD4', 'Proliferating_NK', 'T_CD8', 'T_CD4', 'T_CD4', 'T_CD8')
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"] <- paste0("Old_Cryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"])
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"] <- paste0("Old_NonCryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"])
Cryo_merge1_meta.data <- Cryo_merge1@meta.data
rownames(Cryo_merge1_meta.data) <- sub("_[12]$", "", rownames(Cryo_merge1_meta.data))

Cryo_data <- Read10X_h5("Cryo_filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("Old_Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo_filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("Old_NonCryo_", colnames(NonCryo_data))
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

#

setwd("E://Cryo-TCR/server/auto/")
Cryo_merge2_500 <- readRDS("scRNA_Cryo_2_500.rds")
Cryo_merge2_500@meta.data <- Cryo_merge2_500@meta.data[,c(1:7,18,20:24)]
colnames(Cryo_merge2_500@meta.data)[8]<-'Cluster_Cryo_merge2_500'
Cryo_merge2_500_meta.data <- Cryo_merge2_500@meta.data
rownames(Cryo_merge2_500_meta.data) <- paste0("New_", rownames(Cryo_merge2_500_meta.data))

Cryo_data <- Read10X_h5("Cryo/filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("New_Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo/filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("New_NonCryo_", colnames(NonCryo_data))
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

Cryo_merge2_500 <- subset(Cryo_merge2, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 15)
Cryo_merge2_500 <- Cryo_merge2_500[, colnames(Cryo_merge2_500) %in% rownames(Cryo_merge2_500_meta.data)]
Cryo_merge2_500@meta.data <- Cryo_merge2_500_meta.data
Cryo_merge2_500@meta.data$orig.ident <- paste0("New_", Cryo_merge2_500@meta.data$orig.ident)

PTPRC <- Cryo_merge2_500@assays$RNA@counts['Ptprc', ]
PTPRC <- PTPRC > 0
Cryo_merge2_500$PTPRC <- PTPRC
Cryo_merge_CD45 <- SetIdent(Cryo_merge2_500, value = 'PTPRC')
Cryo_merge_CD45 <- subset(Cryo_merge_CD45, ident = "TRUE")
Cryo_merge2_500$PTPRC <- NULL
Cryo_merge_CD45$PTPRC <- NULL

Cryo_merge_CD45_500_Cluster_T_sub <- read.csv(file = 'Cryo_merge_CD45_500_Cluster_T_sub.csv')
Cryo_merge_CD45@meta.data$Cryo_merge_CD45_500_Cluster_T_sub <- Cryo_merge_CD45_500_Cluster_T_sub$x

##############################################################
################## Integret by merge and SCTransform functions
# Cryo_merge3_500 <- merge(Cryo_merge1, Cryo_merge_CD45, add.cell.ids = c("Old","New"))
Cryo_merge3_500 <- merge(Cryo_merge1, Cryo_merge_CD45)
Cryo_merge3_500[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge3_500, pattern = "^mt-")

Cryo_merge3_500 <- SCTransform(Cryo_merge3_500, vars.to.regress = "percent.mt")
Cryo_merge3_500 <- FindVariableFeatures(Cryo_merge3_500, selection.method = "vst", nfeatures = 2000)
Cryo_merge3_500 <- RunPCA(Cryo_merge3_500)
Cryo_merge3_500 <- RunUMAP(Cryo_merge3_500, dims = 1:10)
Cryo_merge3_500 <- RunTSNE(Cryo_merge3_500, dims = 1:10)
Cryo_merge3_500 <- FindNeighbors(Cryo_merge3_500, dims = 1:12)
Cryo_merge3_500 <- FindClusters(Cryo_merge3_500, resolution = 0.09)

#################
DimPlot(Cryo_merge3_500, reduction = "pca", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "umap", group.by = 'orig.ident', pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "tsne", group.by = 'orig.ident', pt.size = 1)

DimPlot(Cryo_merge3_500, reduction = "pca", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "umap", label = T, label.size = 10, pt.size = 1)
DimPlot(Cryo_merge3_500, reduction = "tsne", label = T, label.size = 10, pt.size = 1)

DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Proliferating', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Proliferating', reduction = 'tsne')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_1', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_1', reduction = 'tsne')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'TCR')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_CD45_500_Cluster_T_sub', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'Cryo_merge_CD45_500_Cluster_T_sub', reduction = 'tsne')

l1 <- Cryo_merge3_500$Proliferating
l1[which(l1 == 'T_CD4')] <- 'CD4+'
l1[which(l1 == 'T_CD8')] <- 'CD8+'
l2 <- Cryo_merge3_500$Cryo_merge_CD45_500_Cluster_T_sub
cluster_combined <- c(l1, l2)[!is.na(c(l1, l2))]
Cryo_merge3_500 <- my_AddMeta(Cryo_merge3_500, cluster_combined)

DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'cluster_combined', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'SCT_snn_res.0.09', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, reduction = 'umap')

cluster_combined_main <- cluster_combined
cluster_combined_main[grep("^CD8+", cluster_combined_main)] <- 'CD8+'
cluster_combined_main[grep("^Proliferating", cluster_combined_main)] <- 'Proliferating'
Cryo_merge3_500 <- my_AddMeta(Cryo_merge3_500, cluster_combined_main)
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'cluster_combined_main', reduction = 'umap')
DimPlot(Cryo_merge3_500, label = F, pt.size = 1, group.by = 'cluster_combined_main', reduction = 'umap')

# saveRDS(Cryo_merge3_500, file = "Cryo_merge3_500_merge_sctransform.rds")
Cryo_merge3_500 <- readRDS('Cryo_merge3_500_merge_sctransform.rds')

DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
VlnPlot(Cryo_merge3_500, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
VlnPlot(Cryo_merge3_500, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'cluster_combined_main')
VlnPlot(Cryo_merge3_500, features = my_MarkersList$universal_markers, group.by = 'cluster_combined_main')
FeaturePlot(Cryo_merge3_500, features = unlist(my_MarkersList[-8]), slot = 'data')
DotPlot(Cryo_merge3_500, features = my_MarkersList[-8], group.by = 'cluster_combined_main', scale = F)
VlnPlot(Cryo_merge3_500, features = unlist(my_MarkersList[-8]), group.by = 'cluster_combined_main')

# Cryo_merge3_500 <- FindClusters(Cryo_merge3_500, resolution = c(1:9/100,1:10/10))
# clustree(Cryo_merge3_500)
DimPlot(Cryo_merge3_500, label = T, label.size = 5, pt.size = 1, group.by = 'SCT_snn_res.0.09', reduction = 'umap')

temp <- my_CountCluster(Cryo_merge3_500, 'SCT_snn_res.0.09')
str(temp)
temp <- as.data.frame(temp)
temp$new <- temp$New_Cryo_count + temp$New_NonCryo_count
temp$old <- temp$Old_Cryo_count + temp$Old_NonCryo_count
temp$Cryo <- temp$New_Cryo_count + temp$Old_Cryo_count 
temp$NonCryo <- temp$New_NonCryo_count + temp$Old_NonCryo_count 
temp$new_percent <- 100 * temp$new / temp$new[9]
temp$old_percent <- 100 * temp$old / temp$old[9]
temp$Cryo_percent <- 100 * temp$Cryo / temp$Cryo[9]
temp$NonCryo_percent <- 100 * temp$NonCryo / temp$NonCryo[9]
temp

temp_obj <- SetIdent(Cryo_merge3_500, value = 'SCT_snn_res.0.09')
Markers_merge_list <- FindAllMarkers(temp_obj, slot = 'data', only.pos = T)
Markers_merge_df <- my_Markers2df_multiple(Markers_merge_list,logFC_threshold = 0, positive = T, n_top = 100)
temp_BP_1 <- my_GO(Markers_merge_df$Cluster_1, return_plot = T, ont = 'BP', Simplify = T, return_res = T)
temp_BP_0 <- my_GO(Markers_merge_df$Cluster_0, return_plot = T, ont = 'BP', Simplify = T, return_res = T)


temp_BP_1@result$Description[grep("inter", temp_BP_1@result$Description)]
temp_BP_0@result$Description[grep("inter", temp_BP_0@result$Description)]
save(Markers_merge_list, Markers_merge_df, file = 'Markers_merge_sctransform.rda')

temp_obj$orig.ident2 <- temp_obj$orig.ident
temp_obj$orig.ident2 <- sub("^.*_Cryo$", "Cryo", temp_obj$orig.ident2)
temp_obj$orig.ident2 <- sub("^.*_NonCryo$", "NonCryo", temp_obj$orig.ident2)
DimPlot(temp_obj, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident2', reduction = 'umap')
temp_obj <- SetIdent(temp_obj, value = 'orig.ident2')
Markers_merge_list2 <- FindAllMarkers(temp_obj, slot = 'data', only.pos = T)
Markers_merge_df2 <- my_Markers2df_multiple(Markers_merge_list2,logFC_threshold = 0, positive = T, n_top = 50)
temp_BP_Cryo <- my_GO(Markers_merge_df2$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T)
temp_BP_NonCryo <- my_GO(Markers_merge_df2$Cluster_NonCryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T)
temp_BP_Cryo@result$Description[grep("inter", temp_BP_Cryo@result$Description)]
temp_BP_NonCryo@result$Description[grep("inter", temp_BP_NonCryo@result$Description)]
temp_BP_Cryo@result[grep("inter", temp_BP_Cryo@result$Description),]
ifn_gene_id <- '66141/54123/100038882/80876/56045/66141/69550/80876'
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL

temp_obj2 <- SetIdent(temp_obj, value = 'cluster_combined_main')
temp_obj2 <- subset(temp_obj2, idents = c('CD4+', 'CD8+', 'NK', 'Proliferating', 'Treg'), )
temp_obj2 <- SetIdent(temp_obj2, value = 'orig.ident2')
Markers_merge_list3 <- FindAllMarkers(temp_obj2, slot = 'data', only.pos = T)
Markers_merge_df3 <- my_Markers2df_multiple(Markers_merge_list3,logFC_threshold = 0, positive = T, n_top = 30)
temp_BP_Cryo <- my_GO(Markers_merge_df3$Cluster_Cryo, return_plot = T, return_res = T)
temp_BP_NonCryo <- my_GO(Markers_merge_df3$Cluster_NonCryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T)
temp_BP_Cryo@result$Description[grep("inter", temp_BP_Cryo@result$Description)]
temp_BP_NonCryo@result$Description[grep("inter", temp_BP_NonCryo@result$Description)]
temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description),]
ifn_gene_id <- '69550/15900/20304/100038882/15900'
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL

temp_obj2 <- SetIdent(temp_obj2, value = 'cluster_combined')
Markers_merge_list4 <- FindAllMarkers(temp_obj2, slot = 'data', only.pos = T)
Markers_merge_df4 <- my_Markers2df_multiple(Markers_merge_list4,logFC_threshold = 0, positive = T, n_top = 200)
temp_BP_CD8_cyto <- my_GO(Markers_merge_df4$'Cluster_CD8+_cyto', return_plot = T, return_res = T, ont = 'BP', Simplify = T)
temp_BP_CD8_cyto_prolifer <- my_GO(Markers_merge_df4$'Cluster_CD8+_cyto_prolifer', return_plot = T, ont = 'BP', Simplify = T, return_res = T)

temp_BP_CD8_cyto@result$Description[grep("inter", temp_BP_CD8_cyto@result$Description)]
temp_BP_CD8_cyto_prolifer@result$Description[grep("inter", temp_BP_CD8_cyto_prolifer@result$Description)]
temp_BP_CD8_cyto@result[grep("interferon", temp_BP_CD8_cyto@result$Description),]
ifn_gene_id <- '20304/11465/74117/22793/16423/22352/68713'
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df4$'Cluster_CD8+_cyto'%in%ifn_gene)

VlnPlot(temp_obj2, c("Ifitm3","Irf7","Isg15","Ifitm2","Samhd1","Bst2"), group.by = 'orig.ident2')
DotPlot(temp_obj2, features = c("Ifitm3","Irf7","Isg15","Ifitm2","Samhd1","Bst2"), group.by = 'orig.ident2')
DimPlot(temp_obj2, label = T, label.size = 5, pt.size = 1, reduction = 'umap')

temp_obj2 <- SCTransform(temp_obj2, vars.to.regress = "percent.mt")
temp_obj2 <- FindVariableFeatures(temp_obj2, selection.method = "vst", nfeatures = 2000)
temp_obj2 <- RunPCA(temp_obj2)
temp_obj2 <- RunUMAP(temp_obj2, dims = 1:10)
temp_obj2 <- RunTSNE(temp_obj2, dims = 1:10)
temp_obj2 <- FindNeighbors(temp_obj2, dims = 1:12)
temp_obj2 <- FindClusters(temp_obj2, resolution = 0.2)

clustree(temp_obj2)
VlnPlot(temp_obj2, features = c("nFeature_SCT", "nCount_SCT", "percent.mt"), ncol = 3, group.by = 'cluster_combined')
VlnPlot(temp_obj2, features = unlist(my_MarkersList[-8])[19:28], ncol = 3, group.by = 'cluster_combined', slot = 'data')

DimPlot(temp_obj2, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
DimPlot(temp_obj2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'orig.ident2')
DimPlot(temp_obj2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'cluster_combined_main')
my_CountCluster(temp_obj2, 'cluster_combined', 'orig.ident')
my_CountCluster(temp_obj2, group2 = 'orig.ident2')
xx <- my_rank_avgExp(temp_obj2, 'cluster_combined_main', assay = 'SCT', slot = 'data')
FeaturePlot(temp_obj2, xx[1:5,5])
##############################################################
##############################################################
################## Integret by RunHarmony
library(harmony)

seuratObject <- RunHarmony(temp_obj, 'orig.ident', lambda = .1, verbose = T, reduction = 'pca')
## Harmony cell embeddings
harmony_embedding <- Seurat::Embeddings(seuratObject, 'harmony')
harmony_embedding[seq_len(5), seq_len(5)]
## Harmony gene loadings
harmony_loadings <- Seurat::Loadings(seuratObject, 'harmony')
harmony_loadings[seq_len(5), seq_len(5)]

p1 <- Seurat::DimPlot(seuratObject, reduction = 'harmony',
                      group.by = 'orig.ident')
p2 <- Seurat::VlnPlot(seuratObject, features = 'harmony_1',
                      group.by = 'orig.ident')
cowplot::plot_grid(p1, p2)


seuratObject <- RunUMAP(seuratObject, dims = 1:10, reduction = 'harmony')
seuratObject <- RunTSNE(seuratObject, dims = 1:10, reduction = 'harmony')
seuratObject <- FindNeighbors(seuratObject, dims = 1:12, reduction = 'harmony')
seuratObject <- FindClusters(seuratObject, resolution = 0.09)
DimPlot(seuratObject, reduction = 'umap', group.by = 'cluster_combined_main', label = T)
DimPlot(seuratObject, reduction = 'umap', label = T)

##############################################################
##############################################################
################## Integret by FindIntegrationAnchors with IntegrateData functions

Cryo.list <- list(Cryo_merge1, Cryo_merge_CD45)
Cryo.list <- lapply(X = Cryo.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Cryo.list)
# Cryo.anchors <- FindIntegrationAnchors(object.list = Cryo.list, anchor.features = features, reduction = 'cca')

Cryo.list <- lapply(X = Cryo.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
})
Cryo.anchors2 <- FindIntegrationAnchors(
    object.list = Cryo.list, 
    reduction = 'rpca'
    )

# this command creates an 'integrated' data assay
# Cryo.combined <- IntegrateData(anchorset = Cryo.anchors)
# DefaultAssay(Cryo.combined) <- "integrated"

Cryo.combined <- IntegrateData(anchorset = Cryo.anchors2, new.assay.name = 'integrated_rpca')
DefaultAssay(Cryo.combined) <- "integrated_rpca"

##############
# mnn
Cryo.list <- list(Cryo_merge1, Cryo_merge_CD45)
Cryo.list <- lapply(X = Cryo.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
Cryo_fastMNN <- RunFastMNN(object.list = Cryo.list)
DefaultAssay(Cryo_fastMNN) <- 'mnn.reconstructed'
Cryo_fastMNN <- RunUMAP(Cryo_fastMNN, reduction = "mnn", dims = 1:10)
Cryo_fastMNN <- RunTSNE(Cryo_fastMNN, reduction = "mnn", dims = 1:10)
Cryo_fastMNN <- FindNeighbors(Cryo_fastMNN, dims = 1:12, reduction = 'mnn')
Cryo_fastMNN <- FindClusters(Cryo_fastMNN, resolution = 0.09)

l1 <- Cryo_fastMNN$Proliferating
l1[which(l1 == 'T_CD4')] <- 'CD4+'
l1[which(l1 == 'T_CD8')] <- 'CD8+'
l2 <- Cryo_fastMNN$Cryo_merge_CD45_500_Cluster_T_sub
cluster_combined <- c(l1, l2)[!is.na(c(l1, l2))]
cluster_combined_main <- cluster_combined
cluster_combined_main[grep("^CD8+", cluster_combined_main)] <- 'CD8+'
cluster_combined_main[grep("^Proliferating", cluster_combined_main)] <- 'Proliferating'
Cryo_fastMNN <- my_AddMeta(Cryo_fastMNN, cluster_combined)
Cryo_fastMNN <- my_AddMeta(Cryo_fastMNN, cluster_combined_main)
DimPlot(Cryo_fastMNN, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined")
DimPlot(Cryo_fastMNN, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined_main")
DimPlot(Cryo_fastMNN, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident', reduction = 'umap')
DimPlot(Cryo_fastMNN, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
DimPlot(Cryo_fastMNN, label = T, label.size = 5, pt.size = 1, group.by = 'cluster_combined', reduction = 'umap', split.by = 'orig.ident')

Cryo_fastMNN$orig.ident2 <- Cryo_fastMNN$orig.ident
Cryo_fastMNN$orig.ident2 <- sub("^.*_Cryo$", "Cryo", Cryo_fastMNN$orig.ident2)
Cryo_fastMNN$orig.ident2 <- sub("^.*_NonCryo$", "NonCryo", Cryo_fastMNN$orig.ident2)
DimPlot(Cryo_fastMNN, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident2', reduction = 'umap')
DimPlot(Cryo_fastMNN_temp, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
###########################

# Run the standard workflow for visualization and clustering
Cryo.combined <- ScaleData(Cryo.combined)
Cryo.combined <- RunPCA(Cryo.combined)
Cryo.combined <- RunUMAP(Cryo.combined, reduction = "pca", dims = 1:10)
Cryo.combined <- RunTSNE(Cryo.combined, reduction = "pca", dims = 1:10)
Cryo.combined <- FindNeighbors(Cryo.combined, dims = 1:12)
Cryo.combined <- FindClusters(Cryo.combined, resolution = 0.09)

# Visualization
p1 <- DimPlot(Cryo.combined, reduction = "umap", split.by = "orig.ident")
p2 <- DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

l1 <- Cryo.combined$Proliferating
l1[which(l1 == 'T_CD4')] <- 'CD4+'
l1[which(l1 == 'T_CD8')] <- 'CD8+'
l2 <- Cryo.combined$Cryo_merge_CD45_500_Cluster_T_sub
cluster_combined <- c(l1, l2)[!is.na(c(l1, l2))]
cluster_combined_main <- cluster_combined
cluster_combined_main[grep("^CD8+", cluster_combined_main)] <- 'CD8+'
cluster_combined_main[grep("^Treg", cluster_combined_main)] <- 'CD4+'
cluster_combined_main[grep("^Proliferating", cluster_combined_main)] <- 'Proliferating'
Cryo.combined <- my_AddMeta(Cryo.combined, cluster_combined)
Cryo.combined <- my_AddMeta(Cryo.combined, cluster_combined_main)
Cryo.combined$orig.ident2 <- Cryo.combined$orig.ident
Cryo.combined$orig.ident2 <- sub("^.*_Cryo$", "Cryo", Cryo.combined$orig.ident2)
Cryo.combined$orig.ident2 <- sub("^.*_NonCryo$", "NonCryo", Cryo.combined$orig.ident2)
Cryo.combined$orig.ident3 <- Cryo.combined$orig.ident
Cryo.combined$orig.ident3 <- sub("^New_.*$", "New", Cryo.combined$orig.ident3)
Cryo.combined$orig.ident3 <- sub("^Old_.*$", "Old", Cryo.combined$orig.ident3)

DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined")
DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined_main")
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident2', reduction = 'umap')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident3', reduction = 'umap')
VlnPlot(Cryo.combined, features = c('nCount_RNA','nFeature_RNA'), group.by = 'cluster_combined_main', pt.size = 0)

Cryo.combined_temp <- SetIdent(Cryo.combined, value = 'orig.ident2')
Markers_merge_list5 <- FindAllMarkers(Cryo.combined_temp, only.pos = T)
Markers_merge_df5 <- my_Markers2df_multiple(Markers_merge_list5,logFC_threshold = 0, positive = T, n_top = 140)
temp_BP_Cryo <- my_GO(Markers_merge_df5$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df5$Cluster_Cryo%in%ifn_gene)

FeaturePlot(Cryo.combined, ifn_gene)
VlnPlot(Cryo.combined, features = ifn_gene, group.by = 'cluster_combined_main', pt.size = 0)
VlnPlot(Cryo.combined, features = ifn_gene, group.by = 'orig.ident2', pt.size = 0)
DotPlot(Cryo.combined, features = ifn_gene, group.by = 'cluster_combined_main') + RotatedAxis()

DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'cluster_combined', reduction = 'umap', split.by = 'orig.ident')

annotation_x <- data.frame(type = Cryo.combined$cluster_combined, group = Cryo.combined$orig.ident2, row.names = rownames(Cryo.combined@meta.data))
annotation_x <- annotation_x[order(annotation_x$group), , drop = F]
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
gap_x <- sapply(unique(annotation_x$type), function(x) which(annotation_x$type == x)[1])[-1] -1
pheatmap::pheatmap(
    Cryo.combined@assays$integrated@data[ifn_gene, rownames(annotation_x)], 
    show_colnames = F, show_rownames = T, 
    annotation_col = annotation_x, #color = c("#EFF3FF","#6BAED6","#6BAED6","#6BAED6","#3182BD"),
    cluster_cols = F, cluster_rows = T, gaps_col = gap_x)

avg_ifn <- AverageExpression(Cryo.combined, features = ifn_gene, assays = 'integrated_rpca', group.by = 'cluster_combined', slot = 'scale.data')[[1]]
pheatmap::pheatmap(avg_ifn, show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, angle_col = 45)

my_DotPlot_split(Cryo.combined, features = ifn_gene, split.by = 'orig.ident2', cols = c('red', 'blue'), dot.scale = 10, group.by = 'cluster_combined') + RotatedAxis()
VlnPlot(Cryo.combined, features = ifn_gene, slot = 'data', group.by = 'cluster_combined', pt.size = 0) 

# clustering
Cryo.combined <- FindClusters(Cryo.combined, resolution = 0.4)
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', split.by = 'orig.ident')
# That means our previous clustering result correct.
# saveRDS(Cryo.combined, file = 'scRNA_Cryo_combined_rpca.rds')
Cryo.combined <- readRDS(file = 'scRNA_Cryo_combined_rpca.rds')

# Sub-clustering of T/NK/Proliferating
T_sub2 <- SetIdent(Cryo.combined, value = 'cluster_combined_main')
T_sub2 <- subset(T_sub2, idents = c('CD4+', 'CD8+', 'NK', 'Proliferating'))
T_sub2 <- SetIdent(T_sub2, value = 'orig.ident2')

T_sub2 <- ScaleData(T_sub2)
T_sub2 <- RunPCA(T_sub2)
T_sub2 <- RunUMAP(T_sub2, reduction = "pca", dims = 1:10)
T_sub2 <- RunTSNE(T_sub2, reduction = "pca", dims = 1:10)
T_sub2 <- FindNeighbors(T_sub2, dims = 1:12)
T_sub2 <- FindClusters(T_sub2, resolution = 0.1)

DimPlot(T_sub2, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined")
DimPlot(T_sub2, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined_main")
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, group.by = 'orig.ident2', reduction = 'umap')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, group.by = 'orig.ident', reduction = 'umap')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
DotPlot(T_sub2, features = my_MarkersList[-8][c('T','NK')], dot.scale = 10)
DotPlot(T_sub2, features = T_marker[[5]], dot.scale = 10, idents = c('0','1','3'))
VlnPlot(T_sub2, features = unlist(T_marker[[5]]), idents = c('0','1','3'), pt.size = 0)
# T_sub2 <- RenameIdents(T_sub2, '0'='CD8+', '1'='CD8+', '2'='CD4+', '3'='CD8+', '4'='NK')
T_sub2$sub_T <- T_sub2$integrated_rpca_snn_res.0.1
levels(T_sub2$sub_T) <- c('CD8_cyto', 'CD8_naive', 'CD4+', 'CD8_cyto', 'NK')
T_sub2 <- SetIdent(T_sub2, value = 'sub_T')

DimPlot(T_sub2, label = T, label.size = 10, pt.size = 1, reduction = 'umap')
my_CountCluster(Cryo.combined, group1 = 'cluster_combined_main', group2 = 'orig.ident2')
my_CountCluster(T_sub2, group2 = 'orig.ident2', trend_factor = 10073 / 14277)
VlnPlot(T_sub2, features = c('nCount_RNA','nFeature_RNA'), pt.size = 0)

# CD8+ 
T_sub2$sub_T_main <- T_sub2$integrated_rpca_snn_res.0.1
levels(T_sub2$sub_T_main) <- c('CD8+', 'CD8+', 'CD4+', 'CD8+', 'NK')
T_sub2 <- SetIdent(T_sub2, value = 'sub_T_main')
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD8+', subcluster.name = 'CD8.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD8+', subcluster.name = 'CD8.1', graph.name = 'integrated_rpca_snn', resolution = 1)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD8+', subcluster.name = 'CD8.1.5', graph.name = 'integrated_rpca_snn', resolution = 1.5)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD8+', subcluster.name = 'CD8.2', graph.name = 'integrated_rpca_snn', resolution = 2)
clustree(T_sub2, prefix = 'CD8.')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.0.5')
my_DotPlot_split(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.0.5', idents = 'CD8+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.0.5', idents = 'CD8+')
VlnPlot(T_sub2, features = unlist(T_marker[[5]]), group.by = 'CD8.0.5', idents = 'CD8+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.1')
my_DotPlot_split(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.1', idents = 'CD8+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.1', idents = 'CD8+')
VlnPlot(T_sub2, features = unlist(T_marker[[5]]), group.by = 'CD8.1', idents = 'CD8+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.1.5')
my_DotPlot_split(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.1.5', idents = 'CD8+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.1.5', idents = 'CD8+')
VlnPlot(T_sub2, features = unlist(T_marker[[5]]), group.by = 'CD8.1.5', idents = 'CD8+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.2')
my_DotPlot_split(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.2', idents = 'CD8+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.2', idents = 'CD8+')
VlnPlot(T_sub2, features = unlist(T_marker[[5]]), group.by = 'CD8.2', idents = 'CD8+', pt.size = 0, slot = 'data')
DotPlot(T_sub2, features = my_MarkersList$marker_gene_cryo1$Proliferating, dot.scale = 8, group.by = 'CD8.2', idents = 'CD8+')
VlnPlot(T_sub2, features = unlist(my_MarkersList$marker_gene_cryo1$Proliferating), group.by = 'CD8.2', idents = 'CD8+', pt.size = 0, slot = 'data')

############### count
Cryo.combined_temp_Old <- subset(SetIdent(Cryo.combined, value = 'orig.ident3'), ident = 'Old')
T_sub2_temp_Old <- subset(SetIdent(T_sub2, value = 'orig.ident3'), ident = 'Old')
my_CountCluster(T_sub2, group1 = 'CD8.2', group2 = 'orig.ident2', trend_factor = 10073 / 14277)
DimPlot(T_sub2_temp_Old, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.2')
my_CountCluster(Cryo.combined_temp_Old, group1 = 'cluster_combined_main', group2 = 'orig.ident2')
my_CountCluster(T_sub2_temp_Old, group1 = 'CD8.2', group2 = 'orig.ident2', trend_factor = 3678 / 3974)
############### 

T_sub2$CD8_detail <- as.factor(T_sub2$CD8.2)
levels(T_sub2$CD8_detail) <- c('CD4+', 
                                 paste0("CD8_", c('Naive', 
                                                  'Effector', 'Naive', 'Cytotoxic', 
                                                  'Effector', 'Central_Memory', 'Central_Memory', 
                                                  'Effector_Memory', 'Cytotoxic', 'Central_Memory', 
                                                  'Naive', 'Naive', 'Effector_Memory', 
                                                  'Naive', 'Effector_Memory', 'Effector', 
                                                  'Effector_Memory', 'Naive', 'Naive', 
                                                  'Effector', 'Central_Memory', 'Cytotoxic'
                                                  )
                                 ), 'NK')
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail')
T_sub2 <- SetIdent(T_sub2, value = 'CD8_detail')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident2')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident3')

# CD4+
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD4+', subcluster.name = 'CD4.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD4+', subcluster.name = 'CD4.1', graph.name = 'integrated_rpca_snn', resolution = 1)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD4+', subcluster.name = 'CD4.1.5', graph.name = 'integrated_rpca_snn', resolution = 1.5)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'CD4+', subcluster.name = 'CD4.2', graph.name = 'integrated_rpca_snn', resolution = 2)
clustree(T_sub2, prefix = 'CD4.')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4.0.5', cells = grep("CD4+", T_sub2$CD4.0.5))
my_DotPlot_split(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.0.5', idents = 'CD4+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.0.5', idents = 'CD4+') + RotatedAxis()
VlnPlot(T_sub2, features = unlist(T_marker[[3]]), group.by = 'CD4.0.5', idents = 'CD4+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4.1', cells = grep("CD4+", T_sub2$CD4.1))
my_DotPlot_split(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.1', idents = 'CD4+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.1', idents = 'CD4+') + RotatedAxis()
VlnPlot(T_sub2, features = unlist(T_marker[[3]]), group.by = 'CD4.1', idents = 'CD4+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4.1.5', cells = grep("CD4+", T_sub2$CD4.1.5))
my_DotPlot_split(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.1.5', idents = 'CD4+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.1.5', idents = 'CD4+') + RotatedAxis()
VlnPlot(T_sub2, features = unlist(T_marker[[3]]), group.by = 'CD4.1.5', idents = 'CD4+', pt.size = 0, slot = 'data')

DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4.2', cells = grep("CD4+", T_sub2$CD4.2))
my_DotPlot_split(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.2', idents = 'CD4+', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = T_marker[[3]], dot.scale = 8, group.by = 'CD4.2', idents = 'CD4+') + RotatedAxis()
VlnPlot(T_sub2, features = unlist(T_marker[[3]]), group.by = 'CD4.2', idents = 'CD4+', pt.size = 0, slot = 'data')

T_sub2$CD4_detail <- as.factor(T_sub2$CD4.1)
levels(T_sub2$CD4_detail) <- c(c('Treg', 'Th1', 'T_FH', 'CD4_Effector_Memory',
                                 'Th1','Th17', 'T_FH', 'Treg', 
                                 'Th17', 'Treg', 'Th17'),
                               "CD8_Central_Memory", "CD8_Cytotoxic", "CD8_Effector", "CD8_Effector_Memory", "CD8_Naive", 'NK')
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail')
T_sub2 <- SetIdent(T_sub2, value = 'CD4_detail')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident2')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident3')

# NK
# all NKs are cytotoxic , except NKT (not cytotoxic but activated)
T_sub2 <- FindSubCluster(T_sub2, cluster = 'NK', subcluster.name = 'NK.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK.0.5', cells = grep("NK", T_sub2$CD4.0.5))
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'TCR', cells = grep("NK", T_sub2$CD4.0.5))
my_DotPlot_split(T_sub2, features = NK_marker, dot.scale = 8, group.by = 'NK.0.5', idents = 'NK', split.by = 'orig.ident2', cols = c('red', 'blue'))
DotPlot(T_sub2, features = NK_marker, dot.scale = 8, group.by = 'NK.0.5', idents = 'NK')
VlnPlot(T_sub2, features = c(unlist(NK_marker), unlist(T_marker$T), "Tcra"), group.by = 'NK.0.5', idents = 'NK', pt.size = 0, slot = 'data')
T_sub2$TCR2 <- T_sub2$TCR
T_sub2$TCR2[is.na(T_sub2$TCR2)] <- 0
T_sub2$TCR2[T_sub2$TCR2=='T_Cryo'] <- 1
T_sub2$TCR2[T_sub2$TCR2=='T_NonCryo'] <- 1
T_sub2$TCR2[T_sub2$TCR2=='Non_T'] <- 0
T_sub2$TCR2 <- as.numeric(T_sub2$TCR2)
VlnPlot(T_sub2, features = 'TCR2', group.by = 'NK.0.5', idents = 'NK', pt.size = 0)

T_sub2$NK_detail <- as.factor(T_sub2$NK.0.5)
levels(T_sub2$NK_detail) <- c(levels(T_sub2$NK_detail)[1:6], "NK", "NK", "NKT", "NK", "NK", levels(T_sub2$NK_detail)[12:15])
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK_detail', cells = grep("NK", T_sub2$CD4.0.5))
DimPlot(T_sub2, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK_detail')

my_CountCluster(T_sub2, group1 = 'NK_detail', group2 = 'orig.ident2', trend_factor = 10073 / 14277)

# sub-cluster over
T_sub2$sub_T_detail <- T_sub2$NK_detail
T_sub2 <- SetIdent(T_sub2, value = 'sub_T_detail')
# saveRDS(T_sub2, file = 'scRNA_Cryo_combined_rpca_T_sub2.rds')
T_sub2 <- readRDS(file = 'scRNA_Cryo_combined_rpca_T_sub2.rds')

my_CountCluster(Cryo.combined, group1 = 'cluster_combined_main', group2 = 'orig.ident2')
my_CountCluster(T_sub2, group2 = 'orig.ident2', trend_factor = 10073 / 14277)

##
Cryo.combined_temp_Old <- subset(SetIdent(Cryo.combined, value = 'orig.ident3'), ident = 'Old')
T_sub2_temp_Old <- subset(SetIdent(T_sub2, value = 'orig.ident3'), ident = 'Old')
my_CountCluster(T_sub2, group2 = 'orig.ident2', trend_factor = 10073 / 14277)
DimPlot(T_sub2_temp_Old, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'sub_T_detail')
my_CountCluster(Cryo.combined_temp_Old, group1 = 'cluster_combined_main', group2 = 'orig.ident2')
my_CountCluster(T_sub2_temp_Old, group1 = 'sub_T_detail', group2 = 'orig.ident2', trend_factor = 3678 / 3974)


Cryo_merge1_temp_temp <- readRDS("D://Ji_Wangxue/T_sub.rds")
Cryo_merge1_temp <- SetIdent(Cryo_merge1, value = 'Cryo_merge_1') 
Cryo_merge1_temp <- my_AddMeta(Cryo_merge1_temp, T_sub2_temp_Old$sub_T_detail, Replace = T)
Cryo_merge1_temp@reductions$umap <- Cryo_merge1_temp_temp@reductions$umap
rownames(Cryo_merge1_temp@reductions$umap@cell.embeddings) <- sapply(
    rownames(Cryo_merge1_temp@reductions$umap@cell.embeddings),
    function(x) {
        if(length(grep("_1$", x))==1)
            paste0("Old_Cryo_", sub("_1$", "", x))
        else
            paste0("Old_NonCryo_", sub("_2$", "", x))
    }
)
DimPlot(Cryo_merge1_temp, label = F, label.size = 5, pt.size = 1, group.by = 'T_sub2_temp_Old_sub_T_detail')
my_CountCluster(Cryo_merge1_temp, group1 = 'T_sub2_temp_Old_sub_T_detail', group2 = 'orig.ident')
##

DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', split.by = 'orig.ident2')
DimPlot(T_sub2, label = F, label.size = 5, pt.size = 1, reduction = 'umap', split.by = 'orig.ident3')
DotPlot(T_sub2, features = my_MarkersList[-8][c('T','NK')], dot.scale = 10)
VlnPlot(T_sub2, features = unlist(my_MarkersList[-8][c('T','NK')]), pt.size = 0)

# T_sub_temp <- T_sub2
# T_sub_temp <- my_ReferenceCluster(T_sub_temp, reference = DB_Mi, label_fine = F)
# DimPlot(T_sub_temp, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'DB_Mi.y')
# #
# Cryo.combined_temp <- SetIdent(Cryo.combined, value = 'cluster_combined_main')
# Cryo.combined_temp <- my_AddMeta(Cryo.combined_temp, new_ident = T_sub_temp$DB_Mi.y, Replace = T)
# DimPlot(Cryo.combined_temp, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'T_sub_temp_DB_Mi.y')

##########
###########################################
# pseudotime
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
cds <- as.cell_data_set(T_sub2)
cds <- cluster_cells(cds, reduction_method = 'UMAP')
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

max.avp <- which(integrated.sub$sub_T_detail=='CD8_Naive')[200]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE, label_roots = F)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
FeaturePlot(integrated.sub, "monocle3_pseudotime")
############################################
# Cryo-NonCryo
T_sub2 <- SetIdent(T_sub2, value = 'orig.ident2')
Markers_merge_list6 <- FindAllMarkers(T_sub2, only.pos = T)
Markers_merge_df6 <- my_Markers2df_multiple(Markers_merge_list6,logFC_threshold = 0, positive = T, n_top = 44)
temp_BP_Cryo <- my_GO(Markers_merge_df6$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description
write.csv(Markers_merge_df6, file = 'temp.csv')

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df6$Cluster_Cryo%in%ifn_gene)

# annotation_x <- data.frame(type = T_sub2$sub_T_detail, group = T_sub2$orig.ident2, row.names = rownames(T_sub2@meta.data))
# annotation_x <- annotation_x[order(annotation_x$group), , drop = F]
# pheatmap::pheatmap(T_sub2@assays$integrated_rpca@scale.data[ifn_gene, rownames(annotation_x)], 
#                    show_colnames = F, show_rownames = T, 
#                    annotation_col = annotation_x, 
#                    cluster_cols = T, cluster_rows = T)

# Among Clusters
T_sub2 <- SetIdent(T_sub2, value = 'sub_T_detail')
Markers_merge_list7 <- FindAllMarkers(T_sub2, only.pos = T)
Markers_merge_df7 <- my_Markers2df_multiple(Markers_merge_list7,logFC_threshold = 0, positive = T, n_top = 200)
write.csv(Markers_merge_df7, file = 'temp.csv')

res <- list()
for (i in colnames(Markers_merge_df7)) {
    temp_BP <- my_GO(
        Markers_merge_df7[[i]], return_plot = T, ont = 'BP', Simplify = T, return_res = T, 
        type = 'bar', font.size = 18, show = 30, title = i)
    temp_BP_df <- temp_BP@result[grep("interferon", temp_BP@result$Description), ]
    
    ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
    ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
    ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
    
    if (length(ifn_gene) != 0)
        res[[i]] <- list(GO = temp_BP, df = temp_BP_df, gene = ifn_gene, rank = which(Markers_merge_df7[[i]]%in%ifn_gene))
}

temp_function <- function(i) {
    cat(i, "\n")
    my_GO(res[[i]]$GO, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, title = i)
    print(res[[i]]$df)
    cat("\n")
}

cat("\n")
for (i in names(res)) {
    xx <- sapply(res[[i]]$df$Description, function(x) {
        if(length(grep("alpha", x))==1)
            a<- 'alpha'
        if(length(grep("beta", x))==1)
            a<- 'beta'
        if(length(grep("gamma", x))==1)
            a<- 'gamma'
        if(length(grep("type I", x))==1)
            a<- 'type I'
        if(length(grep("response", x))==1)
            b <- 'response'
        if(length(grep("production", x))==1)
            b <- 'production'
        ifelse(a == 'type I', 'type I signaling', paste(a, b, sep = "-"))
    })
    cat(
        sub("Cluster_","", i),
        "\t",
        paste(xx, collapse = " | ")
        ,"\n")
}

for (i in seq(res))
    temp_function(names(res)[i])

ifn_high <- unique(unlist(lapply(res, function(x) x$gene)))
annotation_x <- data.frame(type = T_sub2$sub_T_detail, group = T_sub2$orig.ident2, row.names = rownames(T_sub2@meta.data))
annotation_x <- annotation_x[order(annotation_x$group), , drop = F]
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
gap_x <- sapply(unique(annotation_x$type), function(x) which(annotation_x$type == x)[1])[-1] -1
pheatmap::pheatmap(T_sub2@assays$integrated_rpca@data[ifn_high, rownames(annotation_x)], 
                   show_colnames = F, show_rownames = T, 
                   annotation_col = annotation_x, color = c("#EFF3FF","#6BAED6","#6BAED6","#6BAED6","#3182BD"),
                   cluster_cols = F, cluster_rows = T, gaps_col = gap_x)

avg_ifn <- AverageExpression(T_sub2, features = ifn_high, assays = 'integrated_rpca', group.by = 'sub_T_detail', slot = 'data')[[1]]
pheatmap::pheatmap(log1p(avg_ifn), show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, angle_col = 45, cutree_cols = 2, fontsize_row = 8)

DotPlot(T_sub2, features = ifn_high, split.by = 'orig.ident2', cols = c('red', 'blue'), dot.scale = 10) + RotatedAxis()
VlnPlot(T_sub2, features = ifn_high[1:16], slot = 'data', pt.size = 0, group.by = 'orig.ident2') 
VlnPlot(T_sub2, features = ifn_high[17:32], slot = 'data', pt.size = 0, group.by = 'orig.ident') 
VlnPlot(T_sub2, features = ifn_high[33:43], slot = 'data', pt.size = 0, group.by = 'orig.ident') 

############
plot_avg <- function(object, ident_name, genes = ifn_gene, data_slot = 'integrated_rpca') {
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())
    temp_cells <- subset(object, idents = ident_name)
    Idents(temp_cells) <- "orig.ident2"
    
    avg.temp_cells <- as.data.frame(log1p(AverageExpression(temp_cells, verbose = FALSE)[[data_slot]]))
    avg.temp_cells$gene <- rownames(avg.temp_cells)
    
    genes.to.label <- genes
    p1 <- ggplot(avg.temp_cells, aes(Cryo, NonCryo)) + geom_point() + ggtitle("CD8_Effector")
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
    print(p1)
}
for (i in seq(unique(T_sub2@active.ident))) {
    plot_avg(T_sub2, unique(T_sub2@active.ident)[i])
}

FeaturePlot(T_sub2, features = ifn_gene[1:3], split.by = "orig.ident2", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(T_sub2, features = ifn_gene[1:3], split.by = "orig.ident2", group.by = "sub_T_detail", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)
############
############
#TCR
CryoTCR_filtered_contig_annotations <- read.csv('../../data/TCR_Cryo_Cellranger/CryoTCR/outs/filtered_contig_annotations.csv')
CryoTCR_filtered_contig_annotations$barcode <- paste0("New_Cryo_", CryoTCR_filtered_contig_annotations$barcode)
NonCryoTCR_filtered_contig_annotations <- read.csv('../../data/TCR_Cryo_Cellranger/NonCryoTCR/outs/filtered_contig_annotations.csv')
NonCryoTCR_filtered_contig_annotations$barcode <- paste0("New_Cryo_", NonCryoTCR_filtered_contig_annotations$barcode)
CryoTCR_filtered_contig_annotations <- CryoTCR_filtered_contig_annotations[CryoTCR_filtered_contig_annotations$barcode%in%colnames(T_sub2), ]
CryoTCR_filtered_contig_annotations$cluster <- T_sub2$sub_T_detail[CryoTCR_filtered_contig_annotations$barcode]
NonCryoTCR_filtered_contig_annotations <- NonCryoTCR_filtered_contig_annotations[NonCryoTCR_filtered_contig_annotations$barcode%in%colnames(T_sub2), ]
NonCryoTCR_filtered_contig_annotations$cluster <- T_sub2$sub_T_detail[NonCryoTCR_filtered_contig_annotations$barcode]

############

###########################################
###########################################2. batch between first and second one with CD45- cells
###
# since 3' scRNA-seq for visualization, not strictly filtered, so 300 and no-CD45 filter.
###########################################import raw_data
setwd("D://Ji_Wangxue/")
Cryo_merge1 <- readRDS("T_sub.rds")
Cryo_merge1@meta.data <- Cryo_merge1@meta.data[,c(1:6,8,12)]
colnames(Cryo_merge1@meta.data)[7] <- 'Cryo_merge_1'
Cryo_merge1@meta.data$Proliferating <- as.factor(Cryo_merge1@meta.data$Proliferating)
levels(Cryo_merge1@meta.data$Proliferating) <- c('B', 'Myeloid', 'NK', 'Proliferating_T_CD8', 'Proliferating_T_CD8', 'Proliferating_T_CD4', 'Proliferating_NK', 'T_CD8', 'T_CD4', 'T_CD4', 'T_CD8')
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"] <- paste0("Old_Cryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="Cryo"])
rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"] <- paste0("Old_NonCryo_", rownames(Cryo_merge1@meta.data)[Cryo_merge1@meta.data$orig.ident=="NonCryo"])
Cryo_merge1_meta.data <- Cryo_merge1@meta.data
rownames(Cryo_merge1_meta.data) <- sub("_[12]$", "", rownames(Cryo_merge1_meta.data))

Cryo_data <- Read10X_h5("Cryo_filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("Old_Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo_filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("Old_NonCryo_", colnames(NonCryo_data))
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
#
setwd("E://Cryo-TCR/server/auto/")
Cryo_merge2_300 <- readRDS("scRNA_Cryo_2_300.rds")
Cryo_merge2_300@meta.data <- Cryo_merge2_300@meta.data[,1:7]
Cryo_merge2_300_meta.data <- Cryo_merge2_300@meta.data
rownames(Cryo_merge2_300_meta.data) <- paste0("New_", rownames(Cryo_merge2_300_meta.data))

Cryo_data <- Read10X_h5("Cryo/filtered_feature_bc_matrix.h5")
colnames(Cryo_data) <- paste0("New_Cryo_", colnames(Cryo_data))
NonCryo_data <- Read10X_h5("NonCryo/filtered_feature_bc_matrix.h5")
colnames(NonCryo_data) <- paste0("New_NonCryo_", colnames(NonCryo_data))
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
# Cryo_merge2_300 <- Cryo_merge2_300[, colnames(Cryo_merge2_300) %in% rownames(Cryo_merge2_300_meta.data)]
Cryo_merge2_300@meta.data <- Cryo_merge2_300_meta.data
Cryo_merge2_300@meta.data$orig.ident <- paste0("New_", Cryo_merge2_300@meta.data$orig.ident)

Cryo_merge2 <- Cryo_merge2_300
Cryo_merge_CD45_500_Cluster_T_sub <- read.csv(file = 'Cryo_merge_CD45_500_Cluster_T_sub.csv', row.names = 1)
rownames(Cryo_merge_CD45_500_Cluster_T_sub) <- paste0("New_", rownames(Cryo_merge_CD45_500_Cluster_T_sub))
colnames(Cryo_merge_CD45_500_Cluster_T_sub) <- 'Cryo_merge_CD45_500_Cluster_T_sub'
Cryo_merge2 <- my_AddMeta(Cryo_merge2, Cryo_merge_CD45_500_Cluster_T_sub)
#
# Integret
Cryo.list <- list(Cryo_merge1, Cryo_merge2)
Cryo.list <- lapply(X = Cryo.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = Cryo.list)

Cryo.list <- lapply(X = Cryo.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
})
Cryo.anchors2 <- FindIntegrationAnchors(
    object.list = Cryo.list, 
    reduction = 'rpca'
)

Cryo.combined <- IntegrateData(anchorset = Cryo.anchors2, new.assay.name = 'integrated_rpca')
DefaultAssay(Cryo.combined) <- "integrated_rpca"
#
# Run the standard workflow for visualization and clustering
Cryo.combined <- ScaleData(Cryo.combined)
Cryo.combined <- RunPCA(Cryo.combined)
Cryo.combined <- RunUMAP(Cryo.combined, reduction = "pca", dims = 1:10)
Cryo.combined <- RunTSNE(Cryo.combined, reduction = "pca", dims = 1:10)
Cryo.combined <- FindNeighbors(Cryo.combined, dims = 1:12)
Cryo.combined <- FindClusters(Cryo.combined, resolution = 0.1)

# Visualization
p1 <- DimPlot(Cryo.combined, reduction = "umap", split.by = "orig.ident")
p2 <- DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

l1 <- Cryo.combined$Proliferating
l1[which(l1 == 'T_CD4')] <- 'CD4+'
l1[which(l1 == 'T_CD8')] <- 'CD8+'
l2 <- Cryo.combined$Cryo_merge_CD45_500_Cluster_T_sub
cluster_combined <- c(l1, l2)[c(l1, l2)!='Others']
cluster_combined <- cluster_combined[!is.na(cluster_combined)]
cluster_combined_main <- cluster_combined
cluster_combined_main[grep("^CD8+", cluster_combined_main)] <- 'CD8+'
cluster_combined_main[grep("^Treg", cluster_combined_main)] <- 'CD4+'
cluster_combined_main[grep("^Proliferating", cluster_combined_main)] <- 'Proliferating'
Cryo.combined <- my_AddMeta(Cryo.combined, cluster_combined)
Cryo.combined <- my_AddMeta(Cryo.combined, cluster_combined_main)
Cryo.combined$orig.ident2 <- Cryo.combined$orig.ident
Cryo.combined$orig.ident2 <- sub("^.*_Cryo$", "Cryo", Cryo.combined$orig.ident2)
Cryo.combined$orig.ident2 <- sub("^.*_NonCryo$", "NonCryo", Cryo.combined$orig.ident2)
Cryo.combined$orig.ident3 <- Cryo.combined$orig.ident
Cryo.combined$orig.ident3 <- sub("^New_.*$", "New", Cryo.combined$orig.ident3)
Cryo.combined$orig.ident3 <- sub("^Old_.*$", "Old", Cryo.combined$orig.ident3)

DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined")
DimPlot(Cryo.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = "cluster_combined_main")
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident2', reduction = 'umap')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'orig.ident3', reduction = 'umap')
VlnPlot(Cryo.combined, features = c('nCount_RNA','nFeature_RNA'), group.by = 'cluster_combined_main', pt.size = 0)

Cryo.combined_temp <- SetIdent(Cryo.combined, value = 'orig.ident2')
Markers_merge_list8 <- FindAllMarkers(Cryo.combined_temp, only.pos = T)
Markers_merge_df8 <- my_Markers2df_multiple(Markers_merge_list8,logFC_threshold = 0, positive = T, n_top = 110)
temp_BP_Cryo <- my_GO(Markers_merge_df8$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df8$Cluster_Cryo%in%ifn_gene)

FeaturePlot(Cryo.combined, ifn_gene)
DotPlot(Cryo.combined, features = ifn_gene, group.by = 'cluster_combined_main') + RotatedAxis()
plots1 <- VlnPlot(Cryo.combined, features = ifn_gene[1:4], split.by = "orig.ident2", group.by = "cluster_combined_main", 
                 pt.size = 0, combine = FALSE)
plots2 <- VlnPlot(Cryo.combined, features = ifn_gene[5:8], split.by = "orig.ident2", group.by = "cluster_combined_main", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots1, ncol = 1)
wrap_plots(plots = plots2, ncol = 1)

DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, group.by = 'cluster_combined', reduction = 'umap', split.by = 'orig.ident')

annotation_x <- data.frame(type = Cryo.combined$cluster_combined, group = Cryo.combined$orig.ident2, row.names = rownames(Cryo.combined@meta.data))
annotation_x <- annotation_x[order(annotation_x$group), , drop = F]
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
gap_x <- sapply(unique(annotation_x$type), function(x) which(annotation_x$type == x)[1])[-1] -1
pheatmap::pheatmap(
    Cryo.combined@assays$integrated@data[ifn_gene, rownames(annotation_x)], 
    show_colnames = F, show_rownames = T, 
    annotation_col = annotation_x, #color = c("#EFF3FF","#6BAED6","#6BAED6","#6BAED6","#3182BD"),
    cluster_cols = F, cluster_rows = T, gaps_col = gap_x)

avg_ifn <- AverageExpression(Cryo.combined, features = ifn_gene, assays = 'integrated_rpca', group.by = 'cluster_combined', slot = 'scale.data')[[1]]
pheatmap::pheatmap(avg_ifn, show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, angle_col = 45)

my_DotPlot_split(Cryo.combined, features = ifn_gene, split.by = 'orig.ident2', cols = c('red', 'blue'), dot.scale = 10, group.by = 'cluster_combined') + RotatedAxis()

# clustering
Cryo.combined <- FindClusters(Cryo.combined, resolution = 1)

DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', split.by = 'orig.ident')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'TCR')
# That means our previous clustering result correct.
# But with group 'Others' existing, re-clustering required.
DotPlot(Cryo.combined, features = my_MarkersList[-8][c(2:7,10,11,12,15,16)], dot.scale = 10) + RotatedAxis()

Cryo.combined$Re_cluster_all <- Cryo.combined$integrated_rpca_snn_res.1
levels(Cryo.combined$Re_cluster_all) <- c('Unknown', 'T_CD8', 'Unknown', 'B', 'T_CD8',
                                          'B', 'T_CD4', 'Myeloid', 'T_CD8',
                                          'Myeloid', 'T_CD8', 'T_CD8', 'Myeloid', 
                                          'T_CD8', 'NK', 'T_CD4', 'Myeloid',
                                          'Myeloid', 'B', 'Mast', 'B', 
                                          'Plasma', 'T_CD8', 'T_CD4', 'pDC',
                                          'B', 'Neutrophil')
Cryo.combined <- SetIdent(Cryo.combined, value = 'Re_cluster_all')
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap')

Re_cluster_all_main <- as.character(Cryo.combined$Re_cluster_all)
Re_cluster_all_main[Re_cluster_all_main == 'pDC'] <- 'Myeloid'
Re_cluster_all_main[Re_cluster_all_main == 'Mast'] <- 'Myeloid'
Re_cluster_all_main[Re_cluster_all_main == 'Neutrophil'] <- 'Myeloid'
Re_cluster_all_main[Re_cluster_all_main == 'T_CD8'] <- 'T'
Re_cluster_all_main[Re_cluster_all_main == 'T_CD4'] <- 'T'
Re_cluster_all_main[Re_cluster_all_main == 'Plasma'] <- 'B'
Cryo.combined$Re_cluster_all_main <- Re_cluster_all_main
DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Re_cluster_all_main')

my_DotPlot_split(Cryo.combined, features = list(
                                       T = c("Cd4", "Cd8a", "Cd8b1"),
                                       NK = c("Klrb1c", "Ncr1", "Klrk1"),
                                       Myeloid = c("Itgam", "Cd14", "Aif1"),
                                       B = c("Cd19", "Cd79a", "Jchain")
                                       ), 
        dot.scale = 10, group.by = 'Re_cluster_all_main', split.by = 'orig.ident3', cols = c('red', 'blue')) + RotatedAxis()
DotPlot(Cryo.combined, features = list(
                            Neutrophils = c("S100a8", "S100a9"),
                            pDC = c("Cox6a2", "Siglech"),
                            Plasma = "Jchain", 
                            Mast = c('Mcpt1', 'Mcpt2'),
                            NK = c("Klrb1c", "Ncr1", "Klrk1"),
                            Myeloid = c("Itgam", "Cd14", "Aif1"),
                            T_CD4 = "Cd4",
                            B = c("Cd19", "Cd79a"),
                            T_CD8 = c("Cd8a", "Cd8b1")
                            ),
        dot.scale = 10, group.by = 'Re_cluster_all') + RotatedAxis()


# saveRDS(Cryo.combined, file = 'scRNA_Cryo_combined_rpca_nofilter.rds')
Cryo.combined <- readRDS(file = 'scRNA_Cryo_combined_rpca_nofilter.rds')

# subcluster for T_CD8:
######subcluster
# Cryo.combined <- SetIdent(Cryo.combined, value = 'Re_cluster_all')
# Cryo.combined <- FindSubCluster(Cryo.combined, cluster = 'T_CD8', subcluster.name = 'CD8.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)
# Cryo.combined <- FindSubCluster(Cryo.combined, cluster = 'T_CD8', subcluster.name = 'CD8.1', graph.name = 'integrated_rpca_snn', resolution = 1)
# clustree(Cryo.combined, prefix = 'CD8.')
# 
# DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.0.5', cells = grep("T_", Cryo.combined$Re_cluster_all))
# my_DotPlot_split(Cryo.combined, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.0.5', idents = 'T_CD8', split.by = 'orig.ident2', cols = c('red', 'blue'))
# DotPlot(Cryo.combined, features = T_marker[[5]], dot.scale = 10, group.by = 'CD8.0.5', idents = 'T_CD8') + RotatedAxis()
# VlnPlot(Cryo.combined, features = unlist(T_marker[[5]]), group.by = 'CD8.0.5', idents = 'T_CD8', pt.size = 0, slot = 'data')
# DotPlot(Cryo.combined, features = my_MarkersList$marker_gene_cryo1$Proliferating, dot.scale = 10, group.by = 'CD8.0.5', idents = 'T_CD8')
# VlnPlot(Cryo.combined, features = unlist(my_MarkersList$marker_gene_cryo1$Proliferating), group.by = 'CD8.0.5', idents = 'T_CD8', pt.size = 0, slot = 'data')
# 
# DimPlot(Cryo.combined, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8.1', cells = grep("T_", Cryo.combined$Re_cluster_all))
# my_DotPlot_split(Cryo.combined, features = T_marker[[5]], dot.scale = 8, group.by = 'CD8.1', idents = 'T_CD8', split.by = 'orig.ident2', cols = c('red', 'blue'))
# DotPlot(Cryo.combined, features = T_marker[[5]], dot.scale = 10, group.by = 'CD8.1', idents = 'T_CD8') + RotatedAxis()
# VlnPlot(Cryo.combined, features = unlist(T_marker[[5]]), group.by = 'CD8.1', idents = 'T_CD8', pt.size = 0, slot = 'data')
# DotPlot(Cryo.combined, features = my_MarkersList$marker_gene_cryo1$Proliferating, dot.scale = 10, group.by = 'CD8.1', idents = 'T_CD8')
# VlnPlot(Cryo.combined, features = unlist(my_MarkersList$marker_gene_cryo1$Proliferating), group.by = 'CD8.1', idents = 'T_CD8', pt.size = 0, slot = 'data')
######
######sub-object cluster
T_CD8 <- subset(Cryo.combined, ident = 'T_CD8')
T_CD8.list <- SplitObject(T_CD8, split.by = "orig.ident3")
T_CD8.list <- lapply(X = T_CD8.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA' 
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x
})
features <- SelectIntegrationFeatures(object.list = T_CD8.list, nfeatures = 3000)
T_CD8.list <- lapply(X = T_CD8.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
    x
})
T_CD8.anchors <- FindIntegrationAnchors(
    object.list = T_CD8.list, 
    reduction = 'rpca'
)

T_CD8 <- IntegrateData(anchorset = T_CD8.anchors, new.assay.name = 'integrated_rpca')
DefaultAssay(T_CD8) <- "integrated_rpca"
T_CD8 <- ScaleData(T_CD8)
T_CD8 <- RunPCA(T_CD8)
T_CD8 <- RunUMAP(T_CD8, reduction = "pca", dims = 1:10)
T_CD8 <- RunTSNE(T_CD8, reduction = "pca", dims = 1:10)
T_CD8 <- FindNeighbors(T_CD8, dims = 1:12)
T_CD8 <- FindClusters(T_CD8, resolution = 0.1)
T_CD8 <- FindClusters(T_CD8, resolution = 1)

DimPlot(T_CD8, reduction = "umap", label = T, label.size = 10, pt.size = 1, group.by = 'orig.ident')
DimPlot(T_CD8, reduction = "umap", label = T, label.size = 10, pt.size = 1)

DotPlot(T_CD8, features = T_marker[[5]], dot.scale = 10) + RotatedAxis()
DotPlot(T_CD8, features = my_MarkersList$marker_gene_cryo1$Proliferating, dot.scale = 10)
DotPlot(T_CD8, features = list(NK = my_MarkersList$marker_gene_cryo1$NK, T = c("Cd4", "Cd8a", "Cd8b1")), dot.scale = 10)
VlnPlot(T_CD8, features = unlist(T_marker[[5]]), slot = 'data', pt.size = 0)

# # cluster 0/16
# DotPlot(T_CD8, features = list(NK = my_MarkersList$universal_markers, T = c("Cd4", "Cd8a", "Cd8b1")), dot.scale = 10)
# FeaturePlot(T_CD8, 'Cd44')
# Temp_diff <- FindMarkers(T_CD8, ident.1 = '16', only.pos = T)
# Temp_diff_df <- my_Markers2df_1Cluster(Temp_diff, 100)
# Temp_diff_df
# my_GO(Temp_diff_df, ont = 'BP', return_plot = T, return_res = F, Simplify = T)

############### 
T_CD8$CD8_detail <- as.factor(T_CD8$integrated_rpca_snn_res.1)
CD8_df <- rbind(
    data.frame(x='CD8_Naive',y=c(2,6,8)),
    data.frame(x='CD8_Central_Memory',y=c(9,12)),
    data.frame(x='CD8_Effector_Memory',y=c(1,5)),
    data.frame(x='CD8_Effector',y=c(4,10,15)),
    data.frame(x='CD8_CTL',y=c(7,14)),
    data.frame(x='CD8_CTL_Il7ra',y=c(0,16)),
    data.frame(x='NKT',y=11),
    data.frame(x='T_CD4',y=c(3,13))
)
CD8_df <- CD8_df[order(CD8_df$y, decreasing = F),]
levels(T_CD8$CD8_detail) <- CD8_df$x
DimPlot(T_CD8, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail')

# proliferating
T_CD8$CD8_detail_proliferating <- as.factor(T_CD8$integrated_rpca_snn_res.1)
CD8_df_proliferating <- rbind(
    data.frame(x='NKT',y=11),
    data.frame(x='T_CD4',y=c(3,13)),
    data.frame(x='CD8_proliferating',y=c(4,7,10,14)),
    data.frame(x='CD8_others',y=c(0:16)[-(c(4,7,10,11,14,3,13)+1)])
)
CD8_df_proliferating <- CD8_df_proliferating[order(CD8_df_proliferating$y, decreasing = F),]
levels(T_CD8$CD8_detail_proliferating) <- CD8_df_proliferating$x
DimPlot(T_CD8, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail_proliferating')
# Exhausted
T_CD8$CD8_detail_Exhausted <- as.factor(T_CD8$integrated_rpca_snn_res.1)
CD8_detail_Exhausted <- rbind(
    data.frame(x='NKT',y=11),
    data.frame(x='T_CD4',y=c(3,13)),
    data.frame(x='CD8_Exhausted',y=c(1,4,5,12,15)),
    data.frame(x='CD8_Naive',y=c(2,6,8)),
    data.frame(x='CD8_others',y=c(0:16)[-(c(c(1,4,5,12,13,15),c(2,3,6,8),11)+1)])
)
CD8_detail_Exhausted <- CD8_detail_Exhausted[order(CD8_detail_Exhausted$y, decreasing = F),]
levels(T_CD8$CD8_detail_Exhausted) <- CD8_detail_Exhausted$x
DimPlot(T_CD8, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail_Exhausted')
############### count
my_CountCluster(T_CD8, group1 = 'integrated_rpca_snn_res.1', group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_CountCluster(T_CD8, group1 = 'integrated_rpca_snn_res.1', group2 = 'orig.ident', trend_factor = 15446 / 26228)

my_CountCluster(T_CD8, group1 = 'CD8_detail', group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_CountCluster(T_CD8, group1 = 'CD8_detail_proliferating', group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_CountCluster(T_CD8, group1 = 'CD8_detail_Exhausted', group2 = 'orig.ident2', trend_factor = 15446 / 26228)
############### 
############### cluster 10/12
Temp_diff_10 <- FindMarkers(T_CD8, ident.1 = '10', only.pos = T, group.by = 'integrated_rpca_snn_res.1')
Temp_diff_df_10 <- my_Markers2df_1Cluster(Temp_diff_10, 200)
Temp_diff_df_10
my_GO(Temp_diff_df_10, ont = 'All', return_plot = T, return_res = F, Simplify = F, show = 30)

Temp_diff_12 <- FindMarkers(T_CD8, ident.1 = '12', only.pos = T, group.by = 'integrated_rpca_snn_res.1')
Temp_diff_df_12 <- my_Markers2df_1Cluster(Temp_diff_12, 200)
Temp_diff_df_12
my_GO(Temp_diff_df_12, ont = 'All', return_plot = T, return_res = F, Simplify = F, show = 30)
############### 
T_CD8 <- SetIdent(T_CD8, value = 'CD8_detail')
# T_CD8 <- SetIdent(T_CD8, value = 'CD8_detail_Exhausted')
# T_CD8 <- SetIdent(T_CD8, value = 'CD8_detail_proliferating')
DimPlot(T_CD8, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident2')
DimPlot(T_CD8, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident')
DimPlot(T_CD8, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail', split.by = 'orig.ident3')

Cryo.combined <- SetIdent(Cryo.combined, value = 'Re_cluster_all')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD8@meta.data[, 'CD8_detail', drop=F], Replace = T)
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD8@meta.data[, 'CD8_detail_Exhausted', drop=F], Replace = T)
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD8@meta.data[, 'CD8_detail_proliferating', drop=F], Replace = T)
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail',  cells = grep("CD", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail_Exhausted',  cells = grep("CD", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD8_detail_proliferating',  cells = grep("CD", Cryo.combined$CD8_detail))

marker_temp <- T_marker[[5]]; marker_temp$Proliferating <- 'Mki67'; marker_temp$T <- c('Cd4', 'Cd8a', 'Cd8b1')
DotPlot(T_CD8, features = marker_temp, dot.scale = 10, group.by = 'CD8_detail') + RotatedAxis()
DotPlot(T_CD8, features = marker_temp, dot.scale = 10, group.by = 'CD8_detail_Exhausted') + RotatedAxis()
DotPlot(T_CD8, features = marker_temp, dot.scale = 10, group.by = 'CD8_detail_proliferating') + RotatedAxis()

# subcluster for T_CD4:
######sub-object cluster
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD8_detail')
T_CD4 <- subset(Cryo.combined, ident = 'T_CD4')
T_CD4.list <- SplitObject(T_CD4, split.by = "orig.ident3")
T_CD4.list <- lapply(X = T_CD4.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA' 
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x
})
features <- SelectIntegrationFeatures(object.list = T_CD4.list, nfeatures = 3000)
T_CD4.list <- lapply(X = T_CD4.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
    x
})
T_CD4.anchors <- FindIntegrationAnchors(
    object.list = T_CD4.list, 
    reduction = 'rpca'
)

T_CD4 <- IntegrateData(anchorset = T_CD4.anchors, new.assay.name = 'integrated_rpca')
DefaultAssay(T_CD4) <- "integrated_rpca"
T_CD4 <- ScaleData(T_CD4)
T_CD4 <- RunPCA(T_CD4)
T_CD4 <- RunUMAP(T_CD4, reduction = "pca", dims = 1:10)
T_CD4 <- RunTSNE(T_CD4, reduction = "pca", dims = 1:10)
T_CD4 <- FindNeighbors(T_CD4, dims = 1:12)
for (i in grep("integrated_rpca_snn_res", colnames(T_CD4@meta.data), value = T)) {
    T_CD4[[i]] <- NULL
}
T_CD4 <- FindClusters(T_CD4, resolution = 0.3)
my_DotPlot_split(T_CD4, features = T_marker[[3]], dot.scale = 10, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.0.3') + RotatedAxis()
my_DotPlot_split(T_CD4, features = list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK), dot.scale = 10, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.0.3') + RotatedAxis()
VlnPlot(T_CD4, features = unlist(list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK)), pt.size = 0, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.0.3') + RotatedAxis()
my_CountCluster(T_CD4, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
DimPlot(T_CD4, reduction = "umap", label = T, label.size = 10, pt.size = 1)

T_CD4$T_CD48 <- as.factor(T_CD4$integrated_rpca_snn_res.0.3)
levels(T_CD4$T_CD48) <- c(rep('CD4', 4), 'T_unknown', 'CD4', 'CD8')
T_CD4 <- SetIdent(T_CD4, value = 'T_CD48')
DimPlot(T_CD4, reduction = "umap", label = T, label.size = 10, pt.size = 1)
VlnPlot(T_CD4, features = c(unlist(T_marker[[5]]),'Mki67'), pt.size = 0, split.by = 'orig.ident3', group.by = 'T_CD48', idents = 'CD8') + RotatedAxis()
levels(T_CD4$T_CD48) <- c('CD4', 'T_unknown', 'CD8_Effector_Memory')
T_CD4$T_CD48_2 <- T_CD4$T_CD48
levels(T_CD4$T_CD48_2) <- c('CD4', 'T_unknown', 'CD8_Exhausted')

T_CD4 <- SetIdent(T_CD4, value = 'T_CD48')
DimPlot(T_CD4, reduction = "umap", label = F, pt.size = 1, group.by = 'T_CD48')
# T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.0.1', graph.name = 'integrated_rpca_snn', resolution = 0.1)
# T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.0.2', graph.name = 'integrated_rpca_snn', resolution = 0.2)
# T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.0.3', graph.name = 'integrated_rpca_snn', resolution = 0.3)
T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.0.4', graph.name = 'integrated_rpca_snn', resolution = 0.4)
T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)
# T_CD4 <- FindSubCluster(T_CD4, cluster = 'CD4', subcluster.name = 'CD4.1', graph.name = 'integrated_rpca_snn', resolution = 1)
# clustree(T_CD4, prefix = 'CD4.')
T_CD4 <- SetIdent(T_CD4, value = 'CD4.0.4')
my_CountCluster(T_CD4, group1 = 'CD4.0.4', group2 = 'orig.ident3')
DimPlot(T_CD4, reduction = "umap", label = T, label.size = 10, pt.size = 1)
my_CountCluster(T_CD4, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_DotPlot_split(T_CD4, features = T_marker[[3]], dot.scale = 10, split.by = 'orig.ident3', group.by = 'CD4.0.4') + RotatedAxis()
DotPlot(T_CD4, features = T_marker[[3]], dot.scale = 10, group.by = 'CD4.0.4', idents = paste0('CD4_', 0:6)) + RotatedAxis()
VlnPlot(T_CD4, features = unlist(T_marker[[3]]), group.by = 'CD4.0.4', slot = 'scale.data', idents = paste0('CD4_', 0:6), pt.size = 0) + RotatedAxis()
VlnPlot(T_CD4, features = unlist(T_marker[[3]]), group.by = 'CD4.0.4', slot = 'scale.data', idents = paste0('CD4_', 0:6), pt.size = 0, split.by = 'orig.ident2') + RotatedAxis()
DotPlot(T_CD4, features = T_marker$CD4$Th17, dot.scale = 10, group.by = 'CD4.0.4', idents = paste0('CD4_', 0:6)) + RotatedAxis()

T_CD4 <- SetIdent(T_CD4, value = 'CD4.0.5')
my_CountCluster(T_CD4, group1 = 'CD4.0.5', group2 = 'orig.ident3')
DimPlot(T_CD4, reduction = "umap", label = T, label.size = 5, pt.size = 1)
my_CountCluster(T_CD4, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_DotPlot_split(T_CD4, features = T_marker[[3]], dot.scale = 10, split.by = 'orig.ident3', group.by = 'CD4.0.5') + RotatedAxis()
DotPlot(T_CD4, features = T_marker[[3]], dot.scale = 10, group.by = 'CD4.0.5', idents = paste0('CD4_', 0:7)) + RotatedAxis()
VlnPlot(T_CD4, features = unlist(T_marker[[3]]), group.by = 'CD4.0.5', slot = 'scale.data', idents = paste0('CD4_', 0:7), pt.size = 0) + RotatedAxis()
VlnPlot(T_CD4, features = unlist(T_marker[[3]]), group.by = 'CD4.0.5', slot = 'scale.data', idents = paste0('CD4_', 0:7), pt.size = 0, split.by = 'orig.ident2') + RotatedAxis()
DotPlot(T_CD4, features = list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK), dot.scale = 10, group.by = 'CD4.0.5') + RotatedAxis()

DotPlot(T_CD4, features = my_MarkersList$marker_gene_cryo1$Proliferating, dot.scale = 10, group.by = 'CD4.0.5') + RotatedAxis()

# # cluster 2/7
# Temp_diff <- FindMarkers(T_CD4, ident.1 = 'CD4_2', only.pos = T)
# Temp_diff_df <- my_Markers2df_1Cluster(Temp_diff, 200)
# Temp_diff_df
# my_GO(Temp_diff_df, ont = 'BP', return_plot = T, return_res = F, Simplify = T)
# 
# Temp_diff <- FindMarkers(T_CD4, ident.1 = 'CD4_7', only.pos = T)
# Temp_diff_df <- my_Markers2df_1Cluster(Temp_diff, 200)
# Temp_diff_df
# my_GO(Temp_diff_df, ont = 'BP', return_plot = T, return_res = F, Simplify = T)

############### 
T_CD4$CD4_detail <- as.factor(T_CD4$'CD4.0.5')
levels(T_CD4$CD4_detail) <- c('CD4_Naive', 'CD4_Follicular_helper', 'CD4_Sh2d1a', 'CD4_Helper_1', 
                              'CD4_Helper_17', 'CD4_Regulatory', 'CD4_Follicular_helper', 'CD4_Naive',
                              'CD8_Effector_Memory', "T_unknown")
T_CD4$CD4_detail_Exhausted <- T_CD4$CD4_detail
levels(T_CD4$CD4_detail_Exhausted) <- sub("CD8_Effector_Memory", "CD8_Exhausted", levels(T_CD4$CD4_detail_Exhausted))
T_CD4$CD4_detail_proliferating <- T_CD4$CD4_detail
levels(T_CD4$CD4_detail_proliferating) <- sub("CD8_Effector_Memory", "CD8_others", levels(T_CD4$CD4_detail_proliferating))

DimPlot(T_CD4, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail')

T_CD4 <- SetIdent(T_CD4, value = 'CD4_detail')
my_CountCluster(T_CD4, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
DimPlot(T_CD4, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident2')
DimPlot(T_CD4, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident')
DimPlot(T_CD4, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail', split.by = 'orig.ident3')

Cryo.combined <- SetIdent(Cryo.combined, value = 'CD8_detail')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD4@meta.data[, 'CD4_detail', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD8_detail_Exhausted')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD4@meta.data[, 'CD4_detail_Exhausted', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD8_detail_proliferating')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = T_CD4@meta.data[, 'CD4_detail_proliferating', drop=F], Replace = T)
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail',  cells = grep("CD", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail_Exhausted',  cells = grep("CD4", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'CD4_detail_proliferating',  cells = grep("CD4", Cryo.combined$CD8_detail))


# subcluster for NK:
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD4_detail')
NK <- subset(Cryo.combined, ident = c('NK', 'NKT'))
NK.list <- SplitObject(NK, split.by = "orig.ident3")
NK.list <- lapply(X = NK.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA' 
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x
})
features <- SelectIntegrationFeatures(object.list = NK.list, nfeatures = 3000)
NK.list <- lapply(X = NK.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
    x
})
NK.anchors <- FindIntegrationAnchors(
    object.list = NK.list, 
    reduction = 'rpca'
)

NK <- IntegrateData(anchorset = NK.anchors, new.assay.name = 'integrated_rpca')
DefaultAssay(NK) <- "integrated_rpca"
NK <- ScaleData(NK)
NK <- RunPCA(NK)
NK <- RunUMAP(NK, reduction = "pca", dims = 1:10)
NK <- RunTSNE(NK, reduction = "pca", dims = 1:10)
NK <- FindNeighbors(NK, dims = 1:12)
for (i in grep("integrated_rpca_snn_res", colnames(NK@meta.data), value = T)) {
    NK[[i]] <- NULL
}
NK <- FindClusters(NK, resolution = 1)

VlnPlot(NK, features = unlist(list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK)), pt.size = 0, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.1') + RotatedAxis()
my_CountCluster(NK, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
DimPlot(NK, reduction = "umap", label = T, label.size = 10, pt.size = 1)
DotPlot(NK, features = list(NK = NK_marker$NK, T = c("Cd3e", "Cd4", "Cd8a", "Cd8b1")), dot.scale = 10)

NK$NK_detail <- as.factor(NK$integrated_rpca_snn_res.1)
levels(NK$NK_detail) <- c(rep('NK', 4), rep('NKT', 2), rep('NK', 2), 'NKT')
NK <- SetIdent(NK, value = 'NK_detail')
DimPlot(NK, reduction = "umap", label = T, label.size = 10, pt.size = 1)
my_CountCluster(NK, group2 = 'orig.ident2', trend_factor = 15446 / 26228)
DotPlot(NK, features = list(NK = NK_marker$NK, T = c("Cd3e", "Cd4", "Cd8a", "Cd8b1")), dot.scale = 10)
VlnPlot(NK, features = unlist(list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK)), pt.size = 0, split.by = 'orig.ident3', group.by = 'NK_detail') + RotatedAxis()

NK$NK_detail_Exhausted <- NK$NK_detail
NK$NK_detail_proliferating <- NK$NK_detail

Cryo.combined <- SetIdent(Cryo.combined, value = 'CD4_detail')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = NK@meta.data[, 'NK_detail', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD4_detail_Exhausted')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = NK@meta.data[, 'NK_detail_Exhausted', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'CD4_detail_proliferating')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = NK@meta.data[, 'NK_detail_proliferating', drop=F], Replace = T)
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK_detail',  cells = grep("NK", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK_detail_Exhausted',  cells = grep("NK", Cryo.combined$CD8_detail))
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'NK_detail_proliferating',  cells = grep("NK", Cryo.combined$CD8_detail))


# subcluster for Myeloid:
Cryo.combined <- SetIdent(Cryo.combined, value = 'NK_detail')
Myeloid <- subset(Cryo.combined, ident = c('Myeloid', 'Mast', 'pDC', 'Neutrophil'))
Myeloid.list <- SplitObject(Myeloid, split.by = "orig.ident3")
Myeloid.list <- lapply(X = Myeloid.list, FUN = function(x) {
    DefaultAssay(x) <- 'RNA' 
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    x
})
features <- SelectIntegrationFeatures(object.list = Myeloid.list, nfeatures = 3000)
Myeloid.list <- lapply(X = Myeloid.list, FUN = function(x) {
    x <- ScaleData(x, features = features)
    x <- RunPCA(x, features = features)
    x
})
Myeloid.anchors <- FindIntegrationAnchors(
    object.list = Myeloid.list, 
    reduction = 'rpca'
)

Myeloid <- IntegrateData(anchorset = Myeloid.anchors, new.assay.name = 'integrated_rpca')
DefaultAssay(Myeloid) <- "integrated_rpca"
Myeloid <- ScaleData(Myeloid)
Myeloid <- RunPCA(Myeloid)
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:10)
Myeloid <- RunTSNE(Myeloid, reduction = "pca", dims = 1:10)
Myeloid <- FindNeighbors(Myeloid, dims = 1:12)
for (i in grep("integrated_rpca_snn_res", colnames(Myeloid@meta.data), value = T)) {
    Myeloid[[i]] <- NULL
}
Myeloid <- FindClusters(Myeloid, resolution = 0.1)
# clustree(Myeloid)
DimPlot(Myeloid, reduction = "umap", label = T, label.size = 5, pt.size = 1)
DimPlot(Myeloid, reduction = "umap", label = F, label.size = 10, pt.size = 1, group.by = 'orig.ident')
DimPlot(Myeloid, reduction = "umap", label = F, label.size = 10, pt.size = 1, group.by = 'NK_detail')

VlnPlot(Myeloid, features = unlist(list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK)), pt.size = 0, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.1') + RotatedAxis()
my_CountCluster(Myeloid, group2 = 'orig.ident2', trend_factor = 15446 / 26228)

# DotPlot(Myeloid, features = list(T = c('Cd3e', 'Cd247'), CD4 = 'Cd4', CD8 = c('Cd8a', 'Cd8b1'), NK = NK_marker$NK), dot.scale = 10) + RotatedAxis()
DotPlot(Myeloid, features = Myeloid_marker, dot.scale = 10) + RotatedAxis()
VlnPlot(Myeloid, features = c('nFeature_RNA', 'nCount_RNA'), pt.size = 0)
VlnPlot(Myeloid, features = unlist(Myeloid_marker), pt.size = 0, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.0.1')
VlnPlot(Myeloid, features = c("S100a8", "S100a9"), pt.size = 0, split.by = 'orig.ident3', group.by = 'integrated_rpca_snn_res.0.1')

Myeloid$Myeloid_detail_0 <- as.factor(Myeloid$integrated_rpca_snn_res.0.1)
levels(Myeloid$Myeloid_detail_0) <- c(rep('others', 4), 'Mast', 'pDC', 'Neutrophils', 'pDC_Ly6c1')
Myeloid <- SetIdent(Myeloid, value = 'Myeloid_detail_0')
DimPlot(Myeloid, reduction = "umap", label = T, label.size = 5, pt.size = 1)

# Myeloid_temp <- my_ReferenceCluster(Myeloid, reference = DB_Mrna, label_fine = F, celltype_list = c('Monocytes', 'Granulocytes', 'Dendritic cells', 'Macrophages'))
# DimPlot(Myeloid_temp, reduction = "umap", label = F, label.size = 10, pt.size = 1, group.by = 'DB_Mrna')
# Myeloid_temp <- my_ReferenceCluster(Myeloid, reference = DB_Ig, label_fine = F, celltype_list = c('Monocytes', 'DC', 'Neutrophils', 'Macrophages', 'Eosinophils', 'Mast cells'))
# DimPlot(Myeloid_temp, reduction = "umap", label = F, label.size = 10, pt.size = 1, group.by = 'DB_Ig')

Myeloid <- FindSubCluster(Myeloid, cluster = 'others', subcluster.name = 'others.0.5', graph.name = 'integrated_rpca_snn', resolution = 0.5)

DimPlot(Myeloid, reduction = "umap", label = T, label.size = 5, pt.size = 1, group.by = 'others.0.5')
DimPlot(Myeloid, reduction = "tsne", label = T, label.size = 5, pt.size = 1, group.by = 'others.0.5')

Myeloid <- SetIdent(Myeloid, value = 'others.0.5')
DotPlot(Myeloid, features = Myeloid_marker[-c(1, 7, 10, 13:15)], dot.scale = 10, idents = paste0('others_', 0:11)) + RotatedAxis()
VlnPlot(Myeloid, features = c('nFeature_RNA', 'nCount_RNA'), pt.size = 0)
VlnPlot(Myeloid, features = unlist(Myeloid_marker[-c(1, 7, 10, 13:15)]), pt.size = 0, split.by = 'orig.ident3', group.by = 'others.0.5', idents = paste0('others_', 0:11))

Myeloid_marker_temp <- unlist(Myeloid_marker[-c(1, 7, 10, 13:15)])
annotation_col_x <- data.frame(sub_celltype = sub("[0-9]$", "", names(Myeloid_marker_temp)),
                               row.names = unname(Myeloid_marker_temp))
annotation_col_x <- annotation_col_x[order(annotation_col_x$sub_celltype), , drop = F]
annotation_col_x$celltype <- c(rep('cDC', 3), rep('Macrophages', 8), 'Monocyte', rep('Neutrophils', 3))
avg_Myeloid <- AverageExpression(Myeloid, features = rownames(annotation_col_x), slot = 'scale.data', assays = 'integrated_rpca', group.by = 'others.0.5')[[1]]
avg_Myeloid <- t(avg_Myeloid[, grep('others', colnames(avg_Myeloid))])
pheatmap::pheatmap(avg_Myeloid, show_colnames = T, show_rownames = T, cluster_cols = F, cluster_rows = T, angle_col = 45, annotation_col = annotation_col_x)


Myeloid.DB_Ig.ref <- my_ReferenceCluster(Myeloid, reference = DB_Ig, label_fine = F, celltype_list = c('Monocytes', 'DC', 'Neutrophils', 'Macrophages', 'Eosinophils', 'Mast cells'), forPlot_likelihood = T)

my_ReferenceScore(
    Myeloid.DB_Ig.ref, Myeloid$others.0.5, avg = T, avg_method = 'arithmetic',
    ref_order = c('Mast.cells', 'Neutrophils', 'Monocytes', 'DC', 'Macrophages'),  
    meta_order = 'others_', cluster_rows = T)
my_ReferenceScore(
    Myeloid.DB_Ig.ref, Myeloid$others.0.5, avg = T, avg_method = 'geometric',
    ref_order = c('Mast.cells', 'Neutrophils', 'Monocytes', 'DC', 'Macrophages'),  
    meta_order = 'others_', cluster_rows = T)
my_ReferenceScore(
    Myeloid.DB_Ig.ref, Myeloid$others.0.5, avg = T, avg_method = 'harmonic',
    ref_order = c('Mast.cells', 'Neutrophils', 'Monocytes', 'DC', 'Macrophages'),  
    meta_order = 'others_', cluster_rows = T)

my_ReferenceLabelStat(Myeloid.DB_Ig.ref, Myeloid$others.0.5, Heatmap = T)


Myeloid$Myeloid_detail <- as.factor(Myeloid$others.0.5)
levels(Myeloid$Myeloid_detail) <- c('Mast', 'Neutrophils', 
                                    'Monocyte', 'Monocyte', 'Neutrophils', 'cDC_activated', 'Macrophage', 'Monocyte', 'cDC_activated', 'Neutrophils', 'Macrophage', 'Neutrophils', 'cDC', 'Neutrophils',
                                    'pDC', 'pDC_Ly6c1')
Myeloid$Myeloid_detail_Exhausted <- Myeloid$Myeloid_detail
Myeloid$Myeloid_detail_proliferating <- Myeloid$Myeloid_detail
Myeloid$Myeloid_detail_main <- Myeloid$Myeloid_detail
levels(Myeloid$Myeloid_detail_main) <- c('Mast', 'Neutrophils', 'Monocyte', 'cDC', 'Macrophage', 'cDC', 'pDC', 'pDC')

DimPlot(Myeloid, label = T, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail')
DimPlot(Myeloid, label = T, label.size = 5, pt.size = 1, reduction = 'tsne', group.by = 'Myeloid_detail')
DimPlot(Myeloid, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail_main')
Myeloid <- SetIdent(Myeloid, value = 'Myeloid_detail')
my_CountCluster(Myeloid, group1 = 'Myeloid_detail', group2 = 'orig.ident2', trend_factor = 15446 / 26228)
my_CountCluster(Myeloid, group1 = 'Myeloid_detail_main', group2 = 'orig.ident2', trend_factor = 15446 / 26228)


Cryo.combined <- SetIdent(Cryo.combined, value = 'NK_detail')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = Myeloid@meta.data[, 'Myeloid_detail', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'NK_detail_Exhausted')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = Myeloid@meta.data[, 'Myeloid_detail_Exhausted', drop=F], Replace = T)
Cryo.combined <- SetIdent(Cryo.combined, value = 'NK_detail_proliferating')
Cryo.combined <- my_AddMeta(Cryo.combined, new_ident = Myeloid@meta.data[, 'Myeloid_detail_proliferating', drop=F], Replace = T)
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail')
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail_Exhausted')
DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail_proliferating')

Cryo.combined <- SetIdent(Cryo.combined, value = 'Myeloid_detail')

CountCluster_heatmap(object = Cryo.combined, group1 = 'Myeloid_detail', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_Exhausted', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_Exhausted', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_proliferating', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_proliferating', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

DimPlot(Cryo.combined, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail', cells = grep('Old', Cryo.combined$orig.ident3))

# subset Old
Cryo.Old <- subset(SetIdent(Cryo.combined, value = 'orig.ident3'), ident = "Old")
DimPlot(Cryo.Old, label = F, label.size = 5, pt.size = 1, reduction = 'umap', group.by = 'Myeloid_detail')

CountCluster_heatmap(object = Cryo.Old, group1 = 'Myeloid_detail', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

CountCluster_heatmap(Cryo.Old, group1 = 'Myeloid_detail_Exhausted', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_Exhausted', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

CountCluster_heatmap(Cryo.Old, group1 = 'Myeloid_detail_proliferating', group2 = 'orig.ident2')
# CountCluster_heatmap(Cryo.combined, group1 = 'Myeloid_detail_proliferating', group2 = 'orig.ident2', trend_factor = (15446-1886)/(26228-7483))

##############################################################
# saveRDS(Cryo.combined, file = 'scRNA_Cryo_combined_rpca_nofilter.rds')
Cryo.combined <- readRDS(file = 'scRNA_Cryo_combined_rpca_nofilter.rds')
# Cryo.Old <- subset(SetIdent(Cryo.combined, value = 'orig.ident3'), ident = "Old")

VlnPlot(Cryo.combined, features = c('nFeature_RNA', 'nCount_RNA'), pt.size = 0)
VlnPlot(Cryo.Old, features = c('nFeature_RNA', 'nCount_RNA'), pt.size = 0, group.by = 'Myeloid_detail')

CountCluster_plot(Cryo.combined, 'Myeloid_detail', xlable_y = -150)
CountCluster_plot(Cryo.Old, 'Myeloid_detail', xlable_y = -15)
#
############################################
############################################ GO analysis
# Cryo-NonCryo with unknown
Cryo.combined_forGO <- SetIdent(Cryo.combined, value = 'Myeloid_detail')

Cryo.combined_forGO <- SetIdent(Cryo.combined_forGO, value = 'orig.ident2')
Markers_merge_list9 <- FindAllMarkers(Cryo.combined_forGO, only.pos = T)
Markers_merge_df9 <- my_Markers2df_multiple(Markers_merge_list9,logFC_threshold = 0, positive = T, n_top = 110)
temp_BP_Cryo <- my_GO(Markers_merge_df9$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description
write.csv(Markers_merge_df9, file = 'temp.csv')

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df9$Cluster_Cryo%in%ifn_gene)

annotation_x <- data.frame(type = Cryo.combined_forGO$Myeloid_detail, group = Cryo.combined_forGO$orig.ident2, row.names = rownames(Cryo.combined_forGO@meta.data))
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
pheatmap::pheatmap(Cryo.combined_forGO@assays$integrated_rpca@scale.data[ifn_gene, rownames(annotation_x)],
                   show_colnames = F, show_rownames = T,
                   annotation_col = annotation_x,
                   cluster_cols = F, cluster_rows = T)

# Cryo-NonCryo without unknown
Cryo.combined_forGO2 <- SetIdent(Cryo.combined, value = 'Myeloid_detail')
Cryo.combined_forGO2 <- subset(Cryo.combined_forGO2, idents = levels(Cryo.combined_forGO2)[levels(Cryo.combined_forGO2)!='Unknown'])
Cryo.combined_forGO2 <- SetIdent(Cryo.combined_forGO2, value = 'orig.ident2')

Markers_merge_list10 <- FindAllMarkers(Cryo.combined_forGO2, only.pos = T)
Markers_merge_df10 <- my_Markers2df_multiple(Markers_merge_list10,logFC_threshold = 0, positive = T, n_top = 105)
temp_BP_Cryo <- my_GO(Markers_merge_df10$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description
write.csv(Markers_merge_df10, file = 'temp.csv')

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df9$Cluster_Cryo%in%ifn_gene)

annotation_x <- data.frame(type = Cryo.combined_forGO2$Myeloid_detail, group = Cryo.combined_forGO2$orig.ident2, row.names = rownames(Cryo.combined_forGO2@meta.data))
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
pheatmap::pheatmap(Cryo.combined_forGO2@assays$integrated_rpca@scale.data[ifn_gene, rownames(annotation_x)],
                   show_colnames = F, show_rownames = T,
                   annotation_col = annotation_x,
                   cluster_cols = F, cluster_rows = T)

# Cryo-NonCryo without unknown, Old
Cryo.Old_forGO <- SetIdent(Cryo.Old, value = 'Myeloid_detail')
Cryo.Old_forGO <- SetIdent(Cryo.Old_forGO, value = 'orig.ident2')

Markers_merge_list11 <- FindAllMarkers(Cryo.Old_forGO, only.pos = T)
Markers_merge_df11 <- my_Markers2df_multiple(Markers_merge_list11,logFC_threshold = 0, positive = T, n_top = 115)
temp_BP_Cryo <- my_GO(Markers_merge_df10$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = T, type = 'bar')
my_GO(temp_BP_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30)
temp_BP_Cryo_df <- temp_BP_Cryo@result[grep("interferon", temp_BP_Cryo@result$Description), ]
temp_BP_Cryo_df$Description
write.csv(Markers_merge_df11, file = 'temp.csv')

ifn_gene_id <- paste0(temp_BP_Cryo_df$geneID, collapse = "/")
ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
which(Markers_merge_df10$Cluster_Cryo%in%ifn_gene)

annotation_x <- data.frame(type = Cryo.Old_forGO$Myeloid_detail, group = Cryo.Old_forGO$orig.ident2, row.names = rownames(Cryo.Old_forGO@meta.data))
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
Cryo.Old_forGO <- ScaleData(Cryo.Old_forGO, assay = 'RNA')
pheatmap::pheatmap(Cryo.Old_forGO@assays$RNA@scale.data[ifn_gene, rownames(annotation_x)],
                   show_colnames = F, show_rownames = T,
                   annotation_col = annotation_x,
                   cluster_cols = F, cluster_rows = T)

# Among Clusters
# Cryo.combined_forGO for Combined;  Cryo.Old_forGO for Old
Cryo.Old_forGO <- SetIdent(Cryo.Old, value = 'Myeloid_detail')

Cryo.Old_forGO <- SetIdent(Cryo.Old_forGO, value = 'Myeloid_detail')
Markers_merge_list11 <- FindAllMarkers(Cryo.Old_forGO, only.pos = T)
Markers_merge_df11 <- my_Markers2df_multiple(Markers_merge_list11,logFC_threshold = 0.25, positive = T, n_top = 200)
write.csv(Markers_merge_df11, file = 'temp.csv')

res <- list()
for (i in colnames(Markers_merge_df11)) {
    temp_BP <- my_GO(
        Markers_merge_df11[[i]], return_plot = F, ont = 'BP', Simplify = T, return_res = T, 
        type = 'bar', font.size = 18, show = 30, title = i)
    temp_BP_df <- temp_BP@result[grep("interferon", temp_BP@result$Description), ]
    
    ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
    ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
    ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
    
    if (length(ifn_gene) != 0)
        res[[i]] <- list(GO = temp_BP, df = temp_BP_df, gene = ifn_gene, rank = which(Markers_merge_df11[[i]]%in%ifn_gene))
}

temp_function <- function(i) {
    cat(i, "\n")
    my_GO(res[[i]]$GO, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, title = i)
    print(res[[i]]$df)
    cat("\n")
}

cat("\n")
for (i in names(res)) {
    xx <- sapply(res[[i]]$df$Description, function(x) {
        if(length(grep("alpha", x))==1)
            a<- 'alpha'
        if(length(grep("beta", x))==1)
            a<- 'beta'
        if(length(grep("gamma", x))==1)
            a<- 'gamma'
        if(length(grep("type I", x))==1)
            a<- 'type I'
        if(length(grep("response", x))==1)
            b <- 'response'
        if(length(grep("production", x))==1)
            b <- 'production'
        ifelse(a == 'type I', 'type I signaling', paste(a, b, sep = "-"))
    })
    cat(
        sub("Cluster_","", i),
        "\t",
        paste(xx, collapse = " | ")
        ,"\n")
}

for (i in seq(res))
    temp_function(names(res)[i])

ifn_high <- unique(unlist(lapply(res, function(x) x$gene)))
annotation_x <- data.frame(type = Cryo.Old_forGO$Myeloid_detail, group = Cryo.Old_forGO$orig.ident2, row.names = rownames(Cryo.Old_forGO@meta.data))
annotation_x <- annotation_x[order(annotation_x$group), , drop = F]
annotation_x <- annotation_x[order(annotation_x$type), , drop = F]
gap_x <- sapply(unique(annotation_x$type), function(x) which(annotation_x$type == x)[1])[-1] -1
pheatmap::pheatmap(Cryo.Old_forGO@assays$integrated_rpca@data[ifn_high, rownames(annotation_x)], 
                   show_colnames = F, show_rownames = T, 
                   annotation_col = annotation_x, color = c("#EFF3FF","#BDD7E7","#6BAED6","#3182BD"),
                   cluster_cols = F, cluster_rows = T, gaps_col = gap_x)

# avg_ifn <- AverageExpression(Cryo.Old_forGO, features = ifn_high, assays = 'integrated_rpca', group.by = 'sub_T_detail', slot = 'data')[[1]]
# pheatmap::pheatmap(log1p(avg_ifn), show_colnames = T, show_rownames = T, cluster_cols = T, cluster_rows = T, angle_col = 45, cutree_cols = 2, fontsize_row = 8)

DotPlot(Cryo.Old_forGO, features = ifn_high, split.by = 'orig.ident2', cols = c('red', 'blue'), dot.scale = 8) + RotatedAxis()
VlnPlot(Cryo.Old_forGO, features = c("Ifitm1", "Ifitm2", "Plscr1"), slot = 'data', pt.size = 0, group.by = 'orig.ident') 

############
plot_avg <- function(object, ident_name, genes = ifn_high, data_slot = 'integrated_rpca') {
    library(ggplot2)
    library(cowplot)
    theme_set(theme_cowplot())
    temp_cells <- subset(object, idents = ident_name)
    Idents(temp_cells) <- "orig.ident2"
    
    avg.temp_cells <- as.data.frame(log1p(AverageExpression(temp_cells, verbose = FALSE)[[data_slot]]))
    avg.temp_cells$gene <- rownames(avg.temp_cells)
    
    genes.to.label <- genes
    p1 <- ggplot(avg.temp_cells, aes(Cryo, NonCryo)) + geom_point() + ggtitle(ident_name)
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
    print(p1)
}
for (i in seq(unique(Cryo.Old_forGO@active.ident))) {
    plot_avg(Cryo.Old_forGO, unique(Cryo.Old_forGO@active.ident)[i])
}

FeaturePlot(T_sub2, features = ifn_gene[1:3], split.by = "orig.ident2", max.cutoff = 3, 
            cols = c("grey", "red"))
plots <- VlnPlot(T_sub2, features = ifn_gene[1:3], split.by = "orig.ident2", group.by = "sub_T_detail", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

# Among Clusters, with main cell type of Old 
Cryo.Old_forGO <- SetIdent(Cryo.Old, value = 'Myeloid_detail')
levels(Cryo.Old_forGO$Myeloid_detail) <- sub("^CD4.*", "CD4", levels(Cryo.Old_forGO$Myeloid_detail))
levels(Cryo.Old_forGO$Myeloid_detail) <- sub("^CD8.*", "CD8", levels(Cryo.Old_forGO$Myeloid_detail))
levels(Cryo.Old_forGO$Myeloid_detail) <- sub("^pDC.*", "pDC", levels(Cryo.Old_forGO$Myeloid_detail))
levels(Cryo.Old_forGO$Myeloid_detail) <- sub("^cDC.*", "cDC", levels(Cryo.Old_forGO$Myeloid_detail))
Cryo.Old_forGO <- SetIdent(Cryo.Old_forGO, value = 'Myeloid_detail')

Markers_merge_list12 <- FindAllMarkers(Cryo.Old_forGO, only.pos = T)
Markers_merge_df12 <- my_Markers2df_multiple(Markers_merge_list12,logFC_threshold = 0.25, positive = T, n_top = 200)

res <- list()
for (i in colnames(Markers_merge_df12)) {
    temp_BP <- my_GO(
        Markers_merge_df12[[i]], return_plot = T, ont = 'BP', Simplify = T, return_res = T, 
        type = 'bar', font.size = 18, show = 30, title = i)
    temp_BP_df <- temp_BP@result[grep("interferon", temp_BP@result$Description), ]
    
    ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
    ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
    ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
    
    if (length(ifn_gene) != 0)
        res[[i]] <- list(GO = temp_BP, df = temp_BP_df, gene = ifn_gene, rank = which(Markers_merge_df12[[i]]%in%ifn_gene))
}
cat("\n")
for (i in names(res)) {
    xx <- sapply(res[[i]]$df$Description, function(x) {
        if(length(grep("alpha", x))==1)
            a<- 'alpha'
        if(length(grep("beta", x))==1)
            a<- 'beta'
        if(length(grep("gamma", x))==1)
            a<- 'gamma'
        if(length(grep("type I", x))==1)
            a<- 'type I'
        if(length(grep("response", x))==1)
            b <- 'response'
        if(length(grep("production", x))==1)
            b <- 'production'
        ifelse(a == 'type I', 'type I signaling', paste(a, b, sep = "-"))
    })
    cat(
        sub("Cluster_","", i),
        "\t",
        paste(xx, collapse = " | ")
        ,"\n")
}
############
############
#TCR
CryoTCR_filtered_contig_annotations <- read.csv('../../data/TCR_Cryo_Cellranger/CryoTCR/outs/filtered_contig_annotations.csv')
CryoTCR_filtered_contig_annotations$barcode <- paste0("New_Cryo_", CryoTCR_filtered_contig_annotations$barcode)
NonCryoTCR_filtered_contig_annotations <- read.csv('../../data/TCR_Cryo_Cellranger/NonCryoTCR/outs/filtered_contig_annotations.csv')
NonCryoTCR_filtered_contig_annotations$barcode <- paste0("New_Cryo_", NonCryoTCR_filtered_contig_annotations$barcode)
CryoTCR_filtered_contig_annotations <- CryoTCR_filtered_contig_annotations[CryoTCR_filtered_contig_annotations$barcode%in%colnames(T_sub2), ]
CryoTCR_filtered_contig_annotations$cluster <- T_sub2$sub_T_detail[CryoTCR_filtered_contig_annotations$barcode]
NonCryoTCR_filtered_contig_annotations <- NonCryoTCR_filtered_contig_annotations[NonCryoTCR_filtered_contig_annotations$barcode%in%colnames(T_sub2), ]
NonCryoTCR_filtered_contig_annotations$cluster <- T_sub2$sub_T_detail[NonCryoTCR_filtered_contig_annotations$barcode]

############
