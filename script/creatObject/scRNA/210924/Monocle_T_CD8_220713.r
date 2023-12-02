source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load("T_sub_real_meta_220713.rda")
load("Combined_analysis_TNK_real.rda")

if(!dir.exists("T_Pseudotime"))
    dir.create("T_Pseudotime")
setwd("T_Pseudotime")

T_sub_real@meta.data <- T_meta_220713
CD8_idents <- c("CD8_Naive","CD8_Tem","CD8_Tex","CD8_Tex_Proliferative",'CD8_Tem_I-IFN')

T_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'IFN_Final')
T_CD8_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, 
                             group.by = 'IFN_Final',
                             cells = T_sub_real$IFN_Final %in% CD8_idents)

ggsave(filename = "T_anno_dim.pdf", plot = T_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "T_CD8_anno_dim.pdf", plot = T_CD8_anno_dim, width = 9.81, height = 6.19)



## pseudotime: Monocle3

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

Idents(T_sub_real) <- 'IFN_Final'
T_CD8 <- subset(T_sub_real, ident = CD8_idents)
cds <- as.cell_data_set(T_CD8, group.by = 'IFN_Final')
cds <- cluster_cells(cds, cluster_method = 'louvain')
tmp <- cds@colData@listData$IFN_Final_Cluster
names(tmp) <- names(cds@clusters@listData$UMAP$clusters)
cds@clusters@listData$UMAP$clusters <- tmp
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
# cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$IFN_Final =='CD8_Naive')[1]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
integrated.sub$monocle3_pseudotime <- ifelse(
    is.infinite(integrated.sub$monocle3_pseudotime),
    0,
    integrated.sub$monocle3_pseudotime
)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

ggsave(
    'T_CD8_state_Monocle3.pdf',
    plot_cells(cds, 
               label_groups_by_cluster = F, label_leaves = F, label_branch_points = F),
    width = 9.81, height = 6.19
)
ggsave(
    'T_CD8_pseudotime_Monocle3.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", 
               label_groups_by_cluster = F, label_leaves = F, label_branch_points = F),
    width = 9.81, height = 6.19
)
ggsave(
    'T_CD8_cluster_Monocle3.pdf',
    plot_cells(cds, color_cells_by = "IFN_Final", 
               label_groups_by_cluster = F, label_leaves = F, label_branch_points = F),
    width = 9.81, height = 6.19
)

detach("package:monocle3", unload = TRUE)



## Monocle2-diffgenes CD8
cds_all_CD8 <- my_create_monocle(T_CD8)
cds_all_CD8_diffgenes <- my_featureSelect_cds(
    cds_all_CD8,
    method = "diffgenes",
    seurat_obj = T_CD8,
    FilterCondition = "p_val_adj<1e-6"
    # FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)"
)
sum(cds_all_CD8_diffgenes@featureData@data[["use_for_ordering"]])

cds_all_CD8_diffgenes <- my_process_cds(cds_all_CD8_diffgenes)
my_plotPseudo(
    cds_all_CD8_diffgenes,
    seurat_obj = T_CD8,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_CD8_diffgenes,
    seurat_obj = T_CD8,
    color_legend.position = c(.8, .2),
    color_by = 'IFN_Final',
    orig.ident = "orig.ident2"
)

plot_cell_trajectory(cds_all_CD8_diffgenes, color_by = 'IFN_Final')
plot_cell_trajectory(cds_all_CD8_diffgenes, color_by = 'State')

T_CD8 <- my_AddSeuratPseudo(T_CD8, cds_all_CD8_diffgenes, reduction_key = 'Monocle_Diffgenes')
my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final')
my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final', 
           cells = grep("I-IFN", T_CD8$IFN_Final))

ggsave(
    "T_CD8_Pseudotime_Monocle2Diffgenes.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, 
               group.by = 'IFN_Final'),
    width = 9.81, height = 6.19
)
ggsave(
    "T_CD8_Pseudotime_Monocle2Diffgenes_1IFN.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, 
               group.by = 'IFN_Final', cells = grep("I-IFN", T_CD8$IFN_Final)),
    width = 9.81, height = 6.19
)


## Sub-Cluster Monocle3
Idents(T_sub_real) <- "IFN_Final"
CD8_sub <- subset(T_sub_real, ident=c('T_CD8','T_Naive'))
CD8_sub <- my_process_seurat(CD8_sub, nVariableFeatures = 2000, normalize = T, norm.method = 'Log', default.assay = 'RNA', mt.to.regress = T, tsne = F)
CD8_sub <- FindClusters(CD8_sub, resolution = 0.5)
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
key_tmp <- levels(CD8_sub$MainCluster_paper)
Idents(CD8_sub) <- 'MainCluster_paper'
CD8_sub <- subset(CD8_sub, ident = key_tmp[!grepl('^NK', key_tmp)])
if (!file.exists('CD8_sub_220704.rda'))
    save(CD8_sub, file = 'CD8_sub_220704.rda')

CD8_sub_anno_dim <- my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_paper')
CD8_sub_dim <- my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'RNA_snn_res.0.5')
ggsave(filename = "CD8_sub_anno_dim.pdf", plot = CD8_sub_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "CD8_sub_dim.pdf", plot = CD8_sub_dim, width = 9.81, height = 6.19)

