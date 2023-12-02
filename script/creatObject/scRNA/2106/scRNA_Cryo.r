setwd("E://Cryo-TCR/server/auto/")
library(Seurat)
library(ggplot2)
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
Cryo_merge <- merge(Cryo, NonCryo)
Cryo_merge[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge, pattern = "^mt-")

#visualization
# vln_cryo_m <- VlnPlot(Cryo_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# plot1_m <- FeatureScatter(Cryo_merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2_m <- FeatureScatter(Cryo_merge, feature1 = "nFeature_RNA", feature2 = "percent.mt")
# plot3_m <- FeatureScatter(Cryo_merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# vln_cryo_m
# plot1_m + plot2_m + plot3_m

Cryo_merge <- subset(Cryo_merge, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 15)
Cryo_merge

Cryo_merge <- SCTransform(Cryo_merge)
Cryo_merge <- FindVariableFeatures(Cryo_merge, selection.method = "vst", nfeatures = 2000)
Cryo_merge <- RunPCA(Cryo_merge)

Cryo_merge <- FindNeighbors(Cryo_merge, dims = 1:12)
Cryo_merge <- FindClusters(Cryo_merge, resolution = 0.05)

Cryo_merge <- RunUMAP(Cryo_merge, dims = 1:10)
Cryo_merge <- RunTSNE(Cryo_merge, dims = 1:10)

# DimPlot(Cryo_merge, reduction = "pca") + labs(title = "Cryo_merge")
# DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1)
# DimPlot(Cryo_merge, reduction = "tsne", label = T, pt.size = 1)
# 
# DimPlot(Cryo_merge, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident')
# DimPlot(Cryo_merge, reduction = "tsne", label = T, pt.size = 1, group.by = 'orig.ident')
# 
saveRDS(Cryo_merge, file="scRNA_Cryo_2.rds")

# ####### merge process
# 
# setwd("D://Ji_Wangxue/")
# library(Seurat)
# library(ggplot2)
# Cryo_merge <- readRDS("Cryo_merge_sct.rds")
# setwd("E://Cryo-TCR/server/auto/")
# Cryo_merge2 <- readRDS("scRNA_Cryo_auto.rds")
# Cryo_merge@meta.data$orig.ident<-paste0(Cryo_merge@meta.data$orig.ident,"_old")
# Cryo_merge2@meta.data$orig.ident<-paste0(Cryo_merge2@meta.data$orig.ident,"_new")
# Cryo_merge3 <- merge(Cryo_merge, Cryo_merge2)
# 
# Cryo_merge3[["percent.mt"]] <- PercentageFeatureSet(Cryo_merge3, pattern = "^mt-")
# Cryo_merge3 <- subset(Cryo_merge3, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 10)
# Cryo_merge3 <- SCTransform(Cryo_merge3)
# Cryo_merge3 <- FindVariableFeatures(Cryo_merge3, selection.method = "vst", nfeatures = 2000)
# Cryo_merge3 <- RunPCA(Cryo_merge3)
# Cryo_merge3 <- RunUMAP(Cryo_merge3, dims = 1:10)
# Cryo_merge3 <- RunTSNE(Cryo_merge3, dims = 1:10)
# 
# DimPlot(Cryo_merge3, reduction = "pca", group.by = 'orig.ident')
# DimPlot(Cryo_merge3, reduction = "umap", group.by = 'orig.ident')
# DimPlot(Cryo_merge3, reduction = "tsne", group.by = 'orig.ident')
# 
# str(Cryo_merge3,3)
# 
# Cryo_merge3 <- FindNeighbors(Cryo_merge3, dims = 1:12)
# Cryo_merge3 <- FindClusters(Cryo_merge3, resolution = 0.3)
# 
# DimPlot(Cryo_merge3, reduction = "pca")
# DimPlot(Cryo_merge3, reduction = "umap")
# DimPlot(Cryo_merge3, reduction = "tsne")

##
library(clustree)

setwd("E://Cryo-TCR/server/auto/")
Cryo_merge2 <- readRDS("scRNA_Cryo_2.rds")

Cryo_merge2  <- FindClusters(Cryo_merge2, resolution = seq(1,10,1)/10)
clustree(Cryo_merge2)
# 0.1
Cryo_merge2 <- SetIdent(Cryo_merge2, value = 'SCT_snn_res.0.1')

DimPlot(Cryo_merge2, reduction = "pca")
DimPlot(Cryo_merge2, reduction = "umap", label = T, label.size = 5, group.by = 'SCT_snn_res.0.1')
DimPlot(Cryo_merge2, reduction = "tsne", label = T, label.size = 5, group.by = 'SCT_snn_res.0.1')

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

trans_df_markers <- function(df_input, n_top, logFC_threshold = 0.25, positive = TRUE) {
  index <- if (positive) df_input$avg_log2FC >= logFC_threshold else df_input$avg_log2FC <= -logFC_threshold
  df_input <- df_input[index, ]
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


Cryo_m_number <- cluster_number(Cryo_merge2, Percentage = T, Round = T)
Cryo_m_number

# Markers_cryo_m <- FindAllMarkers(Cryo_merge2, test.use = 't')
# Markers_cryo_m_w <- FindAllMarkers(Cryo_merge2, test.use = 'wilcox')
# save(Markers_cryo_m, Markers_cryo_m_w, file = "Markers_cryo_m.rda")
load(file = "Markers_cryo_m.rda")

Markers_cryo_merge_df <- trans_df_markers(Markers_cryo_m, 1000, 0)
Markers_cryo_merge_df_w <- trans_df_markers(Markers_cryo_m_w, 1000, 0)

for (i in c(1:3,5:9)) print(FeaturePlot(Cryo_merge2, Markers_cryo_merge_df_w[1:12,i], reduction = 'umap'))
FeaturePlot(Cryo_merge2, Markers_cryo_merge_df_w[1,4], reduction = 'umap')

universal_markers <- c('Ptprc','Epcam', 'Pecam1')
B_markers <- "Cd19"
T_markers <- c("Cd3e", "Cd247", "Cd4", "Cd8a")
NK_markers <- c("Klrb1c", "Ncr1", "Klrk1")
Myeloid_markers <- c("Itgam", "Adgre1", "Cd14", "Itgax","Cd83")
marker_gene_1 <- list(
  Proliferating = c("Stmn1","Tuba1b","Mki67","Hmgn2","Birc5","Pclaf"),
  B = c("Cd79a","Ebf1"),
  NK = c("Gzma","Irf8","Nrarp","Klre1","Txk","Car2"),
  Myeloid = c("Fcer1a","Tpsb2","Cpa3",'Cd63','Alox5ap',"Ctsh","Ly86","Lyz2","Cebpb","Cxcl2","Il1rn","Il1b","Il1r2","Ets2","Ccr1","Cd14"),
  T = c("Cd3g","Trbc2","Cd3d","Trac","Cd3e")
)
Texhausted_markers <- c("Pdcd1", "Ctla4", "Havcr2")
Treg_marker <- "Foxp3"

FeaturePlot(Cryo_merge2, universal_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, B_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, T_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, NK_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, Myeloid_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, marker_gene_1, reduction = 'umap')
FeaturePlot(Cryo_merge2, Texhausted_markers, reduction = 'umap')
FeaturePlot(Cryo_merge2, Treg_marker, reduction = 'umap')

DotPlot(Cryo_merge2, features = marker_gene_1, group.by = 'SCT_snn_res.0.1', dot.scale = 10)
DotPlot(Cryo_merge2, features = list(universal_markers=universal_markers,
                                     B_markers=B_markers,
                                     T_markers=T_markers,
                                     NK_markers=NK_markers,
                                     Myeloid_markers=Myeloid_markers,
                                     Texhausted_markers=Texhausted_markers,
                                     Treg_marker=Treg_marker
                                     ),
        group.by = 'SCT_snn_res.0.1', dot.scale = 10)

Panglao_Marker <- readRDS("E://Cryo_markers_from_PanglaoDB.rds")
Panglao_Marker_filter <- lapply(Panglao_Marker, function(x) x[x%in%rownames(Cryo_merge2)])

Panglao_Marker_temp1 <- Panglao_Marker_filter[c(1,4,5,6,10,11,13,15,16,17,19,20,24)]


DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[1])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[2])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[3])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[4])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[5])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[6])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[7])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[8])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[9])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[10])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[11])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[12])
DotPlot(Cryo_merge2, features = Panglao_Marker_temp1[13])


overlap_markers <- function(markers_DB_list, markers_cluster_diff, duplicate.rm = T) {
  temp <- matrix(nrow=length(markers_DB_list), 
                 ncol=ncol(markers_cluster_diff), 
                 dimnames = list(names(markers_DB_list), colnames(markers_cluster_diff)))
  temp2 <- as.list(colnames(markers_cluster_diff))
  names(temp2) <- colnames(markers_cluster_diff)
  for (j in colnames(temp)) {
    x <- temp[,j]
    y <- names(x)
    temp3 <- as.list(y)
    names(temp3) <- y
    for (i in y) {
      p <- unname(markers_cluster_diff[, j])
      q <- unname(markers_DB_list[[i]])
      pq <- intersect(p, q)
      temp[i,j] <- length(pq)
      temp3[[i]] <- pq
    }
    temp4 <- unique(unlist(temp3)[duplicated(unlist(temp3))])
    temp3 <- lapply(temp3, function(z) z[!z%in%temp4] )
    temp3$'duplicated' <- temp4
    temp2[[j]] <- temp3
  }
  return(list(heat=temp,detail=temp2))
}

overlap_res <- overlap_markers(Panglao_Marker_temp1, Markers_cryo_merge_df_w[1:50,])
overlap_res_t <- overlap_markers(Panglao_Marker_temp1, Markers_cryo_merge_df[1:50,])

library(pheatmap)
# wilcox
pheatmap(overlap_res[[1]], cluster_rows = F, cluster_cols = F, angle_col = 315, fontsize = 15)
for (i in 1:ncol(overlap_res[[1]])) print(DotPlot(Cryo_merge2, features = overlap_res[[2]][[i]], ))
#t
pheatmap(overlap_res_t[[1]], cluster_rows = F, cluster_cols = F, angle_col = 315, fontsize = 15)
for (i in 1:ncol(overlap_res[[1]])) print(DotPlot(Cryo_merge2, features = overlap_res[[2]][[i]], ))


DimPlot(Cryo_merge2, reduction = "umap", label = T)
# according to PanglaoDB, annotation is following:
# split 0 into 0 and 9
levels(Cryo_merge2$SCT_snn_res.0.2) <- c(2,3,1,4,0,5,9,6,7,6,8,1,2,3)
Cryo_merge2 <- SetIdent(Cryo_merge2, value = 'SCT_snn_res.0.2')
# save(Cryo_merge2, file = "scRNA_Cryo_2_modified.rds")

# Markers_cryo_m_2 <- FindAllMarkers(Cryo_merge2, test.use = 't')
# Markers_cryo_m_w2 <- FindAllMarkers(Cryo_merge2, test.use = 'wilcox')
# # save(Markers_cryo_m_2, Markers_cryo_m_w2, file = "Markers_cryo_m2.rda")
load("Markers_cryo_m2.rda")
Markers_cryo_merge_df_w2 <- trans_df_markers(Markers_cryo_m_w2, 1000)
Markers_cryo_merge_df_2 <- trans_df_markers(Markers_cryo_m_2, 1000)

overlap_res_t2 <- overlap_markers(Panglao_Marker_temp1, Markers_cryo_merge_df_2[1:48,])
overlap_res2 <- overlap_markers(Panglao_Marker_temp1, Markers_cryo_merge_df_w2[1:48,])


pheatmap(overlap_res_t2[[1]], cluster_rows = F, cluster_cols = F, angle_col = 45, fontsize = 15)

# GSVA
library(GSVA)
library(GSEABase)
GSV_path <- "C://DeepTL_dev-main/pathway/GSVA/"
gmt_file <- paste0(GSV_path,"c5.go.v7.4.symbols1.gmt")
geneset <- getGmt(gmt_file)
# exp_X_poisson <- Cryo_merge2@assays$RNA@data
# exp_X_poisson <- as.matrix(exp_X_poisson)
# exp_X_poisson <- exp_X_poisson[rownames(Cryo_merge2@assays$SCT@scale.data),]
# rownames(exp_X_poisson) <- toupper(rownames(exp_X_poisson))
exp_X_gaussian <- Cryo_merge2@assays$SCT@data
exp_X_gaussian <- as.matrix(exp_X_gaussian)
exp_X_gaussian <-exp_X_gaussian[rownames(Cryo_merge2@assays$SCT@scale.data),]
rownames(exp_X_gaussian) <- toupper(rownames(exp_X_gaussian))
## GSVA analysis
# # GSVA_res_poisson <- gsva(exp_X_poisson, geneset, min.sz=2, max.sz=1000, verbose=TRUE, method="gsva", kcdf="Poisson")
# GSVA_res_gaussian <- gsva(exp_X_gaussian, geneset, min.sz=2, max.sz=1000, verbose=TRUE, method="gsva", kcdf="Gaussian")
# save(GSVA_res_gaussian, file = "GSVA_res_gaussian_0.1.rda")
load("GSVA_res_gaussian_0.1.rda")

cluster_r0.1 <- Idents(Cryo_merge2)
cluster_r0.1 <- sort(cluster_r0.1)
cluster_r0.1_anno <- as.data.frame(cluster_r0.1)

ifn_ppathway <- grep(toupper("interferon"), rownames(GSVA_res_gaussian), value = T)
ifn_mtx <- GSVA_res_gaussian[ifn_ppathway,names(cluster_r0.1)]
pheatmap::pheatmap(ifn_mtx, show_colnames = F, show_rownames = T, cluster_rows = F, cluster_cols = F, annotation_col = cluster_r0.1_anno)



## differential analysis between Cryo and NonCryo
# Cryo_merge2_temp <- SetIdent(Cryo_merge2, value = 'orig.ident')
# Markers_cryo_group <- FindAllMarkers(Cryo_merge2_temp, test.use = 't')
# Markers_cryo_group_w <- FindAllMarkers(Cryo_merge2_temp, test.use = 'wilcox')
# # save(Markers_cryo_group, Markers_cryo_group_w, file = "Markers_cryo_group.rda")
Markers_cryo_merge_df_group <- trans_df_markers(Markers_cryo_group, 35)
Markers_cryo_merge_df_group_w <- trans_df_markers(Markers_cryo_group_w, 35)


# GO R=0.1
library(clusterProfiler)
library(enrichplot)
options(connectionObserver = NULL)
library("org.Mm.eg.db")
OrgDb = get("org.Mm.eg.db")
library("goProfiles")
GO <- function(Geneset, 
               ont = "ALL", 
               type = 'dot', 
               pvalue = 0.05, 
               qvalue = 0.2, 
               show = 20, 
               font.size = 10, 
               title = NULL, 
               Simplify = FALSE, 
               return_plot = FALSE, 
               return_res = FALSE
               ) {
  message('ont can be one of "BP", "MF", "CC" and "ALL"(default)')
  message('type can be some of "bar", "dot", "cnet" and "emap"')
  if (Simplify & ont == "ALL") stop('simplify function only applies to a single ontology!')
  genes <- bitr(Geneset, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
  ego <- enrichGO(gene = genes$ENTREZID, ont = ont, OrgDb = OrgDb, pvalueCutoff = pvalue, qvalueCutoff = qvalue)
  if (nrow(ego) != 0) {
    if (Simplify) ego <- simplify(ego)
    if ('bar' %in% type) res <- barplot(ego, showCategory = show, font.size = font.size, title = paste("GO", title))
    if ('dot' %in% type) res <- dotplot(ego, showCategory = show, font.size = font.size, title = paste("GO", title))
    if ('cnet' %in% type) res <- cnetplot(ego, showCategory = show, circular = T, colorEdge = T, categorySize = 'p.adjust', title = paste("GO", title))
    if ('emap' %in% type) {
      ego1 <- pairwise_termsim(ego, method = "JC", semData = NULL, showCategory = show)
      res <- emapplot(ego1, showCategory = show)
    }
    if ('go' %in% type) {
      if (ont == "ALL") {
        message("Ontology being ALL cannot goplot! Transfer into one of the three!")
      } else {
        res <- goplot(ego, showCategory = show, color = "p.adjust", layout = "sugiyama", geom = "text")
      }
    }
    if ('top' %in% type) {
      if (ont == "ALL") {
        message("Ontology being ALL cannot plotGOgraph! Transfer into one of the three!")
      } else {
        res <- plotGOgraph(ego, firstSigNodes = show, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
      }
    }
    print(res)
    if (return_plot) return(res)
    if (return_res) return(ego)
  } else {
    message("No output information")
    return(ego)
  }
}
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

## GO
GO(Markers_cryo_merge_df_w[1:100,1], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,1], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,1], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,2], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,2], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,2], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,3], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,3], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,3], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,5], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,5], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,5], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,6], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,6], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,6], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,7], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,7], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,7], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,8], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,8], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,8], ont = "MF", Simplify = T, title = "MF")

GO(Markers_cryo_merge_df_w[1:100,9], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_w[1:100,9], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_w[1:100,9], ont = "MF", Simplify = T, title = "MF")

# Cryo vs NonCryo
GO(Markers_cryo_merge_df_group_w[1:31,1], ont = "BP", Simplify = T, title = "BP")
GO(Markers_cryo_merge_df_group_w[1:31,1], ont = "CC", Simplify = T, title = "CC")
GO(Markers_cryo_merge_df_group_w[1:31,1], ont = "MF", Simplify = T, title = "MF")


## GSEA
# Cryo_merge2_fc <- my_fc(data = Cryo_merge2@assays$SCT@data, 
#                   meta = Cryo_merge2@meta.data$SCT_snn_res.0.1, 
#                   LOG = F,
#                   average_method = "geometric")
# # save(Cryo_merge2_fc, file = "Cryo_merge2_fc_0.1.rda")
load("Cryo_merge2_fc_0.1.rda")

pdf("GSEA_annotation.pdf")
GSEA_anno(Cryo_merge2_fc[,'0'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster0", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'1'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster1", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'2'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster2", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'3'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster3", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'4'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster4", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'5'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster5", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'6'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster6", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'7'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster7", font.size = 4)
GSEA_anno(Cryo_merge2_fc[,'8'], rownames(Cryo_merge2_fc), c("GO", "KEGG"), "ALL", "bonferroni", 200, "Cryo_merge2_cluster8", font.size = 4)
dev.off()

## T cell filter from TCR barcodes:
T_barcodes_Cryo <- read.table("Cryo_T_barcode.txt", header = F)[[1]]
T_barcodes_NonCryo <- read.table("NonCryo_T_barcode.txt", header = F)[[1]]
T_barcodes_Cryo <- paste0("Cryo_", T_barcodes_Cryo)
T_barcodes_NonCryo <- paste0("NonCryo_", T_barcodes_NonCryo)

barcodes_temp <- data.frame(Cryo_merge2=colnames(Cryo_merge2))
barcodes_temp$TCR <- "Non_T"
barcodes_temp$TCR[barcodes_temp$Cryo_merge2%in%T_barcodes_Cryo] <- "T_Cryo"
barcodes_temp$TCR[barcodes_temp$Cryo_merge2%in%T_barcodes_NonCryo] <- "T_NonCryo"

Cryo_merge2@meta.data$TCR <- barcodes_temp$TCR

DimPlot(Cryo_merge2, reduction = "umap", group.by = 'TCR')
DimPlot(Cryo_merge2, reduction = "tsne", group.by = 'TCR')
