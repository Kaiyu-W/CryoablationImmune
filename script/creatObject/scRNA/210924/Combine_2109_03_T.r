source("E://Cryo-TCR/server/auto/utilities.r")
load(file = "C://celldex_reference/Database_form_celldex.rda") # reference data
setwd("E://Cryo-TCR/server/210924/")
load(file = 'Combined_analysis.rda')
load(file = 'Combined_analysis_T_CD4_8.rda')
load(file = 'Combined_analysis_Myeloid.rda')

Cryo_merge_temp <- my_AddMeta(Cryo_merge, T_CD4$T_CD4_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD4_T_CD4_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, T_CD8$T_CD8_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD8_T_CD8_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, Myeloid_sub$M_sub_MainCluster, Replace = T)
Cryo_merge_temp$MainCluster[Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster == 'T_CD8(EM_ex)'] <- 'T'
Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster[Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster == 'T_CD8(EM_ex)'] <- 'CD8_Effector_Memory'
Idents(Cryo_merge_temp) <- 'Myeloid_sub_M_sub_MainCluster'

my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 1.5, label.size = 6, group.by = 'MainCluster', title = 'Cryo_merge MainCluster')
my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 0.5, label.size = 6, group.by = 'Myeloid_sub_M_sub_MainCluster', title = 'Cryo_merge Cluster')
my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 0.5, label.size = 6, group.by = 'orig.ident2', title = 'Cryo_merge Cryo_NonCryo')

my_CountCluster(Cryo_merge_temp, group1 = 'MainCluster', group2 = 'orig.ident2')

# T visualization
Idents(Cryo_merge_temp) <- 'MainCluster'
T_sub <- subset(Cryo_merge_temp, ident = c('T', 'NK'))
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'Myeloid_sub_M_sub_MainCluster', title = 'T_sub Cryo_NonCryo')

# ReCluster
DefaultAssay(T_sub) <- 'RNA'
T_sub <- SCTransform(T_sub)
T_sub <- FindVariableFeatures(T_sub, selection.method = "vst", nfeatures = 2500)
T_sub <- RunPCA(T_sub)
T_sub <- RunUMAP(T_sub, dims = 1:10)
T_sub <- RunTSNE(T_sub, dims = 1:10)
T_sub <- FindNeighbors(T_sub, dims = 1:12)
T_sub <- FindClusters(T_sub, resolution = 0.01)

# cluster
T_sub  <- FindClusters(T_sub, resolution = 1:9/100)
T_sub  <- FindClusters(T_sub, resolution = 1:10/10)
clustree(T_sub)

# my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident2', title = 'T_sub Cryo_NonCryo')
# my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(paste0("0", 1:9), 1:3)))
# my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = paste0('SCT_snn_res.0.', 1:9))

# select R=0.2
Idents(T_sub) <- 'SCT_snn_res.0.2'
my_plotDim(T_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub R=0.2')
my_DotPlot_split(T_sub, features = my_MarkersList[1:4]) + RotatedAxis()

# cluster 9 is epithelial
T_sub_real <- subset(T_sub, ident = 0:8)
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1)

# Re-ReCluster
DefaultAssay(T_sub_real) <- 'RNA'
# T_sub_real <- SCTransform(T_sub_real)
T_sub_real <- NormalizeData(T_sub_real, normalization.method = "LogNormalize", scale.factor = 10000)
T_sub_real <- ScaleData(T_sub_real, assay = "RNA", features = rownames(T_sub_real))
T_sub_real <- FindVariableFeatures(T_sub_real, selection.method = "vst", nfeatures = 2500)
T_sub_real <- RunPCA(T_sub_real)
T_sub_real <- RunUMAP(T_sub_real, dims = 1:10)
T_sub_real <- RunTSNE(T_sub_real, dims = 1:10)
T_sub_real <- FindNeighbors(T_sub_real, dims = 1:12)
T_sub_real <- FindClusters(T_sub_real, resolution = 0.2)

# cluster
T_sub_real  <- FindClusters(T_sub_real, resolution = 1:9/100)
T_sub_real  <- FindClusters(T_sub_real, resolution = 1:10/10)
clustree(T_sub_real)

# select R=0.3
Idents(T_sub_real) <- 'RNA_snn_res.0.3'
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'T_sub_real R=0.3', group.by = 'RNA_snn_res.0.3')

my_DotPlot_split(T_sub_real, features = my_MarkersList[1:7], group.by = 'RNA_snn_res.0.3') + RotatedAxis()
my_violin(T_sub_real, features = unlist(my_MarkersList[1:7]), pt.size = 0, mode = 'mtx')
my_DotPlot_split(T_sub_real, features = T_markers_tmp) + RotatedAxis()
my_violin(T_sub_real, features = unlist(T_markers_tmp), pt.size = 0, mode = 'mtx')
my_violin(T_sub_real, features = unlist(NK_marker), pt.size = 0, mode = 'mtx')
my_CountCluster(T_sub_real, group1 = 'RNA_snn_res.0.3', group2 = 'Myeloid_sub_M_sub_MainCluster')

my_violin(T_sub_real, features = unlist(T_marker$CD4), pt.size = 0, mode = 'mtx', idents = c(4,5,7,9,11))
my_violin(T_sub_real, features = unlist(T_marker$CD8), pt.size = 0, mode = 'mtx', idents = c(0,3,8,10))
my_DotPlot_split(T_sub_real, features = T_marker$CD4, idents = c(4,5,7,11)) + RotatedAxis()
my_DotPlot_split(T_sub_real, features = T_marker$CD8, idents = c(0,3,8,10)) + RotatedAxis()

# ProjecTIL_markers <- list(
#     Th1 = c('Ifngr1', 'Fasl'),
#     Tfh = c("Cxcr5", "Tox", "Slamf6"),
#     Treg = c("Foxp3")
# )
# my_DotPlot_split(T_sub_real, features = ProjecTIL_markers, idents = c(4,5,7,9,11)) + RotatedAxis()
# my_violin(T_sub_real, features = unlist(ProjecTIL_markers), pt.size = 0, mode = 'mtx', idents = c(4,5,7,9,11))

# cell type identify
T_sub_real$T_sub_MainCluster <- as.factor(T_sub_real$RNA_snn_res.0.3)
T_sub_real$T_sub_MainCluster2 <- as.factor(T_sub_real$RNA_snn_res.0.3)
T_sub_MainCluster_df <- rbind(
    data.frame(x='NK',y=c(2,6)),
    data.frame(x='T_Naive',y=c(1)),
    data.frame(x='T_CD4',y=c(4,5,7,11)),
    data.frame(x='T_CD8',y=c(0,3,8,10)),
    data.frame(x='T_others',y=c(9))
)
T_sub_MainCluster_df2 <- rbind(
    data.frame(x='NK',y=c(2,6)),
    data.frame(x='T_Naive',y=c(1)),
    data.frame(x='Treg',y=c(7)),
    data.frame(x='Th17',y=c(11)),
    data.frame(x='Th1',y=c(4)),
    data.frame(x='Tfh',y=c(5)),
    data.frame(x='CD8_Tem',y=c(0)),
    data.frame(x='CD8_Tex',y=c(3,8,10)),
    data.frame(x='T_others',y=c(9))
)
T_sub_MainCluster_df <- T_sub_MainCluster_df[order(T_sub_MainCluster_df$y, decreasing = F),]
T_sub_MainCluster_df2 <- T_sub_MainCluster_df2[order(T_sub_MainCluster_df2$y, decreasing = F),]
if(all(T_sub_MainCluster_df$y == levels(T_sub_real$T_sub_MainCluster))) {
    levels(T_sub_real$T_sub_MainCluster) <- T_sub_MainCluster_df$x
    levels(T_sub_real$T_sub_MainCluster2) <- T_sub_MainCluster_df2$x
}
Idents(T_sub_real) <- 'T_sub_MainCluster'
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster")
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = "T_sub_MainCluster2")

# count
my_CountCluster(T_sub_real, group1 = 'RNA_snn_res.0.3', group2 = 'orig.ident2')
my_CountCluster(T_sub_real, group1 = 'T_sub_MainCluster', group2 = 'orig.ident2')
my_CountCluster(T_sub_real, group1 = 'T_sub_MainCluster2', group2 = 'orig.ident2')

# ProjecTIL project
library(ProjecTILs)
ref_TILAtlas <- readRDS("E://Cryo-TCR/TILAtlas/ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref_TILAtlas, label = T, cols = refCols)

T_sub_real2 <- subset(T_sub_real, ident = c("T_CD8","T_Naive","T_CD4","T_others"))
query_T <- make.projection(T_sub_real2, ref = ref_TILAtlas, skip.normalize = T, fast.mode = T, filter.cells = F)
DimPlot(query_T, group.by = 'T_sub_MainCluster', cols = refCols)
DimPlot(query_T, group.by = 'RNA_snn_res.0.3', reduction = 'umap', label = T)
DimPlot(query_T, group.by = 'T_sub_MainCluster2', reduction = 'umap', label = T)

DimPlot(query_T, group.by = 'RNA_snn_res.0.3', reduction = 'umap', label = T, cells = grep("^10$", query_T$RNA_snn_res.0.3))
plot.projection(ref_TILAtlas, query_T)

query.projected <- cellstate.predict(ref = ref_TILAtlas, query = query_T)
table(query.projected$functional.cluster)
plot.statepred.composition(ref_TILAtlas, query.projected, metric = "Percent")
plot.states.radar(ref_TILAtlas, query = query.projected, min.cells = 30)

markers <- c("Cd4","Cd8a","Ccr7","Tcf7","Pdcd1","Havcr2","Tox","Izumo1r","Cxcr6","Xcl1","Gzmb","Gzmk","Ifng","Foxp3")
VlnPlot(T_sub_real, features = markers,stack = T, flip = T, assay = "RNA")

# save result
save(T_sub_real, query_T, file = 'Combined_analysis_TNK_real.rda')

# save plot
T_sub_real$RNA_0.3 <- T_sub_real$RNA_snn_res.0.3
levels(T_sub_real$RNA_0.3) <- paste0("C", as.numeric(levels(T_sub_real$RNA_0.3))+1)
Idents(T_sub_real) <- "RNA_0.3"
TNK_dim <- my_plotDim(T_sub_real, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, title = 'T/NK R=0.3', group.by = 'RNA_0.3')
TNK_TNK_dot <- my_DotPlot_split(T_sub_real, features = my_MarkersList[3:4]) + RotatedAxis()
TNK_CD4_dot <- my_DotPlot_split(T_sub_real, features = T_marker$CD4, idents = paste0("C", c(4,5,7,11)+1)) + RotatedAxis()
TNK_CD8_dot <- my_DotPlot_split(T_sub_real, features = T_marker$CD8, idents = paste0("C", c(0,3,8,10)+1)) + RotatedAxis()
TNK_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, title = "T/NK Clusters", group.by = "T_sub_MainCluster2")

ggsave(filename = "TNK_dim.png", plot = TNK_dim)
ggsave(filename = "TNK_TNK_dot.png", plot = TNK_TNK_dot)
ggsave(filename = "TNK_CD4_dot.png", plot = TNK_CD4_dot, width = 10)
ggsave(filename = "TNK_CD8_dot.png", plot = TNK_CD8_dot, width = 10)
ggsave(filename = "TNK_anno_dim.png", plot = TNK_anno_dim)

ggsave(filename = "TNK_CD4_dot.pdf", plot = TNK_CD4_dot)
ggsave(filename = "TNK_CD8_dot.pdf", plot = TNK_CD8_dot)

# T/NK annotation count
T_sub_count2 <- as.data.frame(my_CountCluster(T_sub_real, group1 = 'T_sub_MainCluster2', group2 = 'orig.ident2')[-10, 1:2])
T_sub_sum2 <- my_CountCluster(Cryo_merge, group1 = 'MainCluster', group2 = 'orig.ident2')["sum", 1:2]
T_sub_count2[,1] <- T_sub_count2[,1] / T_sub_sum2[1]
T_sub_count2[,2] <- T_sub_count2[,2] / T_sub_sum2[2]
T_sub_count2[,1] <- T_sub_count2[,1] / sum(T_sub_count2[,1])
T_sub_count2[,2] <- T_sub_count2[,2] / sum(T_sub_count2[,2])
T_sub_count2
T_sub_count2$cell <- rownames(T_sub_count2)
T_sub_count2 <- reshape::melt(T_sub_count2, id.vars = 'cell')
T_sub_count2$variable <- sub("_count$", "", T_sub_count2$variable)
T_sub_count2$variable <- ifelse(T_sub_count2$variable == "Cryo", "CA", "Non-CA")
T_sub_count2$variable <- factor(T_sub_count2$variable, levels = c("Non-CA", "CA"))
T_sub_Plot2 <- ggplot(T_sub_count2, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(T_sub_count2)/2 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "TNK_anno_count.png", plot = T_sub_Plot2)

# T/NK count
T_sub_count <- as.data.frame(my_CountCluster(T_sub_real, group1 = 'RNA_0.3', group2 = 'orig.ident2')[-13, 1:2])
T_sub_sum <- my_CountCluster(Cryo_merge, group1 = 'MainCluster', group2 = 'orig.ident2')["sum", 1:2]
T_sub_count[,1] <- T_sub_count[,1] / T_sub_sum[1]
T_sub_count[,2] <- T_sub_count[,2] / T_sub_sum[2]
T_sub_count[,1] <- T_sub_count[,1] / sum(T_sub_count[,1])
T_sub_count[,2] <- T_sub_count[,2] / sum(T_sub_count[,2])
T_sub_count
T_sub_count$cell <- rownames(T_sub_count)
T_sub_count <- reshape::melt(T_sub_count, id.vars = 'cell')
T_sub_count$variable <- sub("_count$", "", T_sub_count$variable)
T_sub_count$variable <- ifelse(T_sub_count$variable == "Cryo", "CA", "Non-CA")
T_sub_count$variable <- factor(T_sub_count$variable, levels = c("Non-CA", "CA"))
T_sub_count$cell <- factor(T_sub_count$cell, levels = paste0("C", 1:12))
T_sub_Plot <- ggplot(T_sub_count, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(T_sub_count)/2 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "TNK_count.png", plot = T_sub_Plot)



# enrich between case-control only in either T_CD8 or T_CD4
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
levels(T_sub_real$T_sub_MainCluster2)
T_sub_real$MainCluster_T <- T_sub_real$T_sub_MainCluster2
levels(T_sub_real$MainCluster_T) <- c("T_CD8","T_Naive","NK","T_CD8","T_CD4","T_CD4","T_CD4","T_others","T_CD4")
DimPlot(T_sub_real, group.by = "MainCluster_T")
# T
Idents(T_sub_real) <- 'MainCluster_T'
T_T <- subset(T_sub_real, ident = c("T_CD8","T_CD4","T_Naive","T_others"))
Idents(T_T) <- 'orig.ident2'
Markers_T_T_list <- FindAllMarkers(T_T, only.pos = T)
Markers_T_T_df <- my_Markers2df_multiple(Markers_T_T_list, logFC_threshold = 0.25, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_T_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_T_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]
# CD8
Idents(T_sub_real) <- 'MainCluster_T'
T_CD8 <- subset(T_sub_real, ident = "T_CD8")
Idents(T_CD8) <- 'orig.ident2'
Markers_T_CD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_T_CD8_df <- my_Markers2df_multiple(Markers_T_CD8_list, logFC_threshold = 0.25, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_CD8_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_CD8_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]
# CD4
Idents(T_sub_real) <- 'MainCluster_T'
T_CD4 <- subset(T_sub_real, ident = "T_CD4")
Idents(T_CD4) <- 'orig.ident2'
Markers_T_CD4_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_T_CD4_df <- my_Markers2df_multiple(Markers_T_CD4_list, logFC_threshold = 0.25, positive = T, n_top = 100)
GO_cryo <- temp_GO(Markers_T_CD4_df$Cluster_Cryo, return_ntop = 50)
GO_noncryo <- temp_GO(Markers_T_CD4_df$Cluster_NonCryo, return_ntop = 50)
diff_pathway_between_Cryo_NonCryo <- setdiff(GO_cryo$top_pathway, GO_noncryo$top_pathway)
diff_pathway_between_Cryo_NonCryo[grep("interferon", diff_pathway_between_Cryo_NonCryo)]

GO_T_CA <- my_GO(Markers_T_T_df$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "CA-T cells")
GO_8_CA <- my_GO(Markers_T_CD8_df$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "CA-T_CD8 cells")
GO_4_CA <- my_GO(Markers_T_CD4_df$Cluster_Cryo, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "CA-T_CD4 cells")
GO_T_CA + theme(plot.title = element_text(size = 25))
GO_8_CA + theme(plot.title = element_text(size = 25))
GO_4_CA + theme(plot.title = element_text(size = 25))

# among Clusters
Idents(T_T) <- 'T_sub_MainCluster2'
Markers_T_T_list2 <- FindAllMarkers(T_T, only.pos = T)
Markers_T_T_df2 <- my_Markers2df_multiple(Markers_T_T_list2, logFC_threshold = 0.25, positive = T, n_top = 200)

GO_T_CD8_Tem <- my_GO(Markers_T_T_df2$Cluster_CD8_Tem, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 48, title = "T_CD8_Tem cells") + theme(plot.title = element_text(size = 25))
GO_T_Naive <- my_GO(Markers_T_T_df2$Cluster_T_Naive, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_Naive cells") + theme(plot.title = element_text(size = 25))
GO_T_CD8_Tex <- my_GO(Markers_T_T_df2$Cluster_CD8_Tex, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_CD8_Tex cells") + theme(plot.title = element_text(size = 25))
GO_T_Th1 <- my_GO(Markers_T_T_df2$Cluster_Th1, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_Th1 cells") + theme(plot.title = element_text(size = 25))
GO_T_Tfh <- my_GO(Markers_T_T_df2$Cluster_Tfh, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_Tfh cells") + theme(plot.title = element_text(size = 25))
GO_T_Treg <- my_GO(Markers_T_T_df2$Cluster_Treg, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_Treg cells") + theme(plot.title = element_text(size = 25))
GO_T_others <- my_GO(Markers_T_T_df2$Cluster_T_others, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 48, title = "T_others cells") + theme(plot.title = element_text(size = 25))
GO_T_Th17 <- my_GO(Markers_T_T_df2$Cluster_Th17, return_plot = T, ont = 'BP', Simplify = T, return_res = F, type = 'bar', font.size = 18, show = 30, nChar = 40, title = "T_Th17 cells") + theme(plot.title = element_text(size = 25))

GO_T_CD8_Tem
GO_T_others

# get the ifn-I genes
ifn_i_gene <- c()
for (i in colnames(Markers_T_T_df2)) {
    temp_BP <- my_GO(Markers_T_T_df2[, i], return_plot = F, ont = 'BP', Simplify = T, return_res = T)
    
    temp_BP_index1 <- grep('interferon', temp_BP@result$Description)
    temp_BP_index2 <- grep('gamma', temp_BP@result$Description)
    temp_BP_index <- setdiff(temp_BP_index1, temp_BP_index2)
    if (length(temp_BP_index) > 0) {
        temp_BP_df <- temp_BP@result[temp_BP_index, , drop = F]
        ifn_gene_id <- paste0(temp_BP_df$geneID, collapse = "/")
        ifn_gene <- unique(strsplit(ifn_gene_id, split = "/")[[1]])
        ifn_gene <- bitr(ifn_gene, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
        ifn_i_gene <- append(ifn_i_gene, ifn_gene)
        ifn_i_gene <- unique(ifn_i_gene)
    }
}
ifn_i_gene

Markers_T_T_ifni <- Markers_T_T_df2[1:length(ifn_i_gene), ]
for (i in colnames(Markers_T_T_df2)) {
    genes_tmp <- intersect(Markers_T_T_df2[, i], ifn_i_gene)
    if (length(genes_tmp) < length(ifn_i_gene))
        genes_tmp <- append(genes_tmp, rep(NA, length(ifn_i_gene)-length(genes_tmp)))
    Markers_T_T_ifni[, i] <- genes_tmp
    
    if (!all(is.na(genes_tmp))) {
        genes_tmp <- genes_tmp[!is.na(genes_tmp)]
        rank_tmp <- sapply(genes_tmp, function(x) which(Markers_T_T_df2[, i]==x))
        index_tmp <- order(rank_tmp, decreasing = F)
        genes_tmp <- genes_tmp[index_tmp]
        rank_tmp <- rank_tmp[index_tmp]
        cat(i, ": ", paste(genes_tmp, rank_tmp, sep = "/", collapse = ", "), "\n")
    }
}
Markers_T_T_ifni

# save(list = ls()[grep("^Markers.*list", ls())], file = "Combined_analysis_T_T_markers.rda")









#########################################cell communication
source("E://Cryo-TCR/server/auto/utilities.r")
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
    objectx <- computeCommunProbPathway(objectx, )
    objectx <- aggregateNet(objectx)
    
    return(objectx)
}

load(file = 'Combined_analysis.rda')
load(file = 'Combined_analysis_Myeloid_real.rda') # Myeloid_sub_real
load(file = 'Combined_analysis_TNK_real.rda') # T_sub_real

DimPlot(Myeloid_sub_real, group.by = "MainCluster2_new")
DimPlot(Cryo_merge, group.by = "MainCluster")
Idents(Cryo_merge) <- "MainCluster"
Cryo_merge_tmp <- my_AddMeta(Cryo_merge, Myeloid_sub_real$MainCluster2_new, Replace = T)
Idents(Cryo_merge_tmp) <- "Myeloid_sub_real_MainCluster2_new"
Cryo_merge_tmp <- my_AddMeta(Cryo_merge_tmp, T_sub_real$T_sub_MainCluster2, Replace = T)
Cryo_merge_tmp$MainCluster_detail <- as.factor(Cryo_merge_tmp$T_sub_real_T_sub_MainCluster2)
levels(Cryo_merge_tmp$MainCluster_detail)[13] <- 'epithelial'

DimPlot(Cryo_merge_tmp, group.by = "MainCluster_detail")

cellchat_obj <- cellchatCreate(Cryo_merge_tmp, 
                               Cryo_merge_tmp@meta.data[, c('MainCluster_detail', 'orig.ident2')]
)
# save(cellchat_obj, file = "Combined_CellChat_obj2.rda")

setwd("E://Cryo-TCR/server/210924/")
load(file = "Combined_CellChat_obj2.rda")

cellchat_obj@netP$pathways
groupSize <- table(cellchat_obj@idents)
C_pathways <- cellchat_obj@netP$pathways[grep("^C", cellchat_obj@netP$pathways)]
lapply(
    C_pathways,
    function(x)
        netVisual_aggregate(cellchat_obj,
                            signaling = x,
                            layout = "circle",
                            pt.title=20,
                            vertex.label.cex = 1.7,
                            vertex.weight = groupSize)
)
netVisual_aggregate(cellchat_obj, 
                    signaling = cellchat_obj@netP$pathways[1], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)

pp1 <- netVisual_aggregate(cellchat_obj, 
                    signaling = C_pathways[1], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)
pp2 <- netVisual_aggregate(cellchat_obj, 
                    signaling = C_pathways[2], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)
pp3 <- netVisual_aggregate(cellchat_obj, 
                    signaling = C_pathways[4], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)
pp4 <- netVisual_aggregate(cellchat_obj, 
                    signaling = C_pathways[6], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)
pp5 <- netVisual_aggregate(cellchat_obj, 
                    signaling = C_pathways[7], 
                    layout = "circle", 
                    pt.title = 20, 
                    vertex.label.cex = 1.7, 
                    vertex.weight = groupSize)
save_png <- function(file, plot, width = 1000, height = 1000) {
    png(filename = file, width = width, height = height)
    print(plot)
    dev.off()
}
save_png("CellChat_C1.png", pp1)
save_png("CellChat_C2.png", pp2)
save_png("CellChat_C3.png", pp3)
save_png("CellChat_C4.png", pp4)
save_png("CellChat_C5.png", pp5)
#########################################