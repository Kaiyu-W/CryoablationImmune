source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load("Combined_analysis_TNK_real.rda")
T_sub_real$RNA_0.3 <- T_sub_real$RNA_snn_res.0.3
levels(T_sub_real$RNA_0.3) <- paste0("C", as.numeric(levels(T_sub_real$RNA_0.3))+1)
Idents(T_sub_real) <- "RNA_0.3"
TNK_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'T/NK R=0.3', group.by = 'RNA_0.3')
TNK_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = "T/NK Clusters", group.by = "T_sub_MainCluster2")

# Klrb1c Fcer1g
my_DotPlot_split(T_sub_real, features = c("Klrb1c","Fcer1g")) + RotatedAxis()
my_DotPlot_split(T_sub_real, features = c("Klrb1c","Fcer1g"), group.by='T_sub_MainCluster2') + RotatedAxis()

FeaturePlot(T_sub_real, features = c("Klrb1c","Fcer1g"), blend = T)
FeaturePlot(T_sub_real, features = c("Klrb1c","Fcer1g","Mki67","Pdcd1"))


# C9 subcluster
T_sub_real <- FindSubCluster(T_sub_real, cluster = 'C9', subcluster.name = 'C9', resolution = 0.1, graph.name = 'RNA_snn')
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'C9', group.by = 'C9')
my_DotPlot_split(T_sub_real, features = c("Klrb1c","Fcer1g","Mki67","Pdcd1"), group.by='C9') + RotatedAxis()
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, title = 'C9', group.by = 'C9', cells = grep("^C9_",T_sub_real$C9))

FeaturePlot(T_sub_real, features = c("Klrb1c","Fcer1g"), blend = T, cells = grep("^C9_",T_sub_real$C9))
FeaturePlot(T_sub_real, features = c("Klrb1c","Fcer1g","Mki67", "Pdcd1"), cells = grep("^C9_",T_sub_real$C9))


anno_df <- T_sub_real@meta.data[grep("^C9_",T_sub_real$C9),'C9',drop=F]
anno_df <- anno_df[order(anno_df$C9),,drop=F]
mtx <- T_sub_real@assays$RNA@data[c("Klrb1c","Fcer1g","Mki67", "Pdcd1"), rownames(anno_df)]
pheatmap::pheatmap(mtx, cluster_cols = F, cluster_rows = F, annotation_col = anno_df, show_colnames = F)


load(file = 'Combined_analysis.rda')
Count_C9_all <- rbind(
    my_CountCluster(T_sub_real, "C9")[paste0("C9_",0:2),1:4],
    'CD8'=my_CountCluster(T_sub_real, 'T_sub_MainCluster')['T_CD8',1:4], 
    'CD45'=my_CountCluster(Cryo_merge)['sum',1:4]
    )
Countperc_C9_all <- apply(Count_C9_all, 2, function(x) x/x[5] * 100)
Countperc_C9_CD8 <- apply(Count_C9_CD8[1:4,], 2, function(x) x/x[4] * 100)

Count_C9_all_combined <- rbind(
    my_CountCluster(T_sub_real, "C9", 'orig.ident2')[paste0("C9_",0:2),1:2], 
    'CD8'=my_CountCluster(T_sub_real, 'T_sub_MainCluster', 'orig.ident2')['T_CD8',1:2],
    'CD45'=my_CountCluster(Cryo_merge, group2='orig.ident2')['sum',1:2]
    )
Countperc_C9_all_combined <- apply(Count_C9_all_combined, 2, function(x) x/x[5] * 100)
Countperc_C9_CD8_combined <- apply(Count_C9_CD8_combined[1:4,], 2, function(x) x/x[4] * 100)


# CD8 subcluster
Idents(T_sub_real) <- "T_sub_MainCluster"
CD8_sub <- subset(T_sub_real, ident=c('T_CD8','T_Naive'))
CD8_sub <- my_process_seurat(CD8_sub, nVariableFeatures = 2000, normalize = T, norm.method = 'Log', default.assay = 'RNA', mt.to.regress = T, tsne = T)

my_plotDim(CD8_sub, group.by = 'orig.ident2')
my_plotDim(CD8_sub, group.by = 'T_sub_MainCluster2', label = T)
CD8_sub <- FindClusters(CD8_sub, resolution = 0.5)
my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5)

Markers_paper <- list(
    "naive" = c("Il7r","Tcf7"), 
    "exausted" = c("Pdcd1", "Tox"), 
    "abIL TCKs" = c("Gzmb","Klrb1c"), 
    "I-IFN act" = c("Isg15", "Ifit3"),
    "Proliferative" = c("Mki67", "Top2a")
    )
my_DotPlot_split(CD8_sub, feature = Markers_paper)
my_violin(CD8_sub, feature = unlist(Markers_paper), mode = 'mtx')
my_violin(CD8_sub, feature = unlist(Markers_paper)[c(seq(1,10,2),seq(2,10,2))], mode = 'raw', ncol = 5)
FeaturePlot(CD8_sub, feature = c(unlist(Markers_paper),'Cd44','Sell'))
FeaturePlot(CD8_sub, unlist(my_MarkersList[3:4]))

anno_list <- list(
    'T_Naive' = c(0,8),
    'CD8_Tem' = c(1,4),
    'CD8_Tex_Proliferative' = 6,
    'I-IFN' = 7,
    'NK_Proliferative' = 10,
    'CD8_Tex' = c(2,3,5,9)
)

CD8_sub <- my_ClusterAnnote(
    Object = CD8_sub, 
    anno_list = anno_list,
    meta_slot = 'RNA_snn_res.0.5',
    anno_slot = 'MainCluster_paper',
    DefaultIdent_changeIf = TRUE,
    Duplicated_removeIf = TRUE
)
my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'MainCluster_paper')
Idents(T_sub_real) <- "T_sub_MainCluster2"
T_sub_real <- my_AddMeta(T_sub_real, CD8_sub$MainCluster_paper, Replace = T)
my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')

# Count:
my_CountCluster(T_sub_real, 'CD8_sub_MainCluster_paper', 'orig.ident2')

T_meta_220630 <- T_sub_real@meta.data
save(T_meta_220630, file = "T_sub_real_meta_220630.rda")
