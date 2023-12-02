source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load("T_sub_real_meta_220630.rda")
load("Combined_analysis_TNK_real.rda")

if(!dir.exists("T_Pseudotime"))
    dir.create("T_Pseudotime")
setwd("T_Pseudotime")

T_sub_real@meta.data <- T_meta_220630

T_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
# my_plotDim(T_sub_real, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
T_CD8_anno_dim <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, 
           group.by = 'CD8_sub_MainCluster_paper',
           cells = T_sub_real$CD8_sub_MainCluster_paper %in% c("T_Naive","CD8_Tem","CD8_Tex","CD8_Tex_Proliferative",'I-IFN'))

ggsave(filename = "T_anno_dim.pdf", plot = T_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "T_CD8_anno_dim.pdf", plot = T_CD8_anno_dim, width = 9.81, height = 6.19)

# pseudotime: Monocle2

# seurat
cds_all <- my_create_monocle(T_sub_real)
cds_all_seurat <- my_featureSelect_cds(
    cds_all,
    method = "seurat",
    seurat_obj = T_sub_real
)
# cds_seurat
cds_all_seurat <- my_process_cds(cds_all_seurat)
my_plotPseudo(
    cds_all_seurat,
    seurat_obj = T_sub_real,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_seurat,
    seurat_obj = T_sub_real,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

T_sub_real <- my_AddSeuratPseudo(T_sub_real, cds_all_seurat, reduction_key = 'Monocle_Seurat')
my_plotDim(T_sub_real, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
my_plotDim(T_sub_real, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(T_sub_real, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')
ggsave(
    "T_Pseudotime_Monocle2.pdf",
    my_plotDim(T_sub_real, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, 
               group.by = 'CD8_sub_MainCluster_paper'),
    width = 9.81, height = 6.19
)

# pseudotime: Monocle3

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

Idents(T_sub_real) <- 'CD8_sub_MainCluster_paper'
T_CD8 <- subset(T_sub_real, ident = c("T_Naive","CD8_Tem","CD8_Tex","CD8_Tex_Proliferative",'I-IFN'))
cds <- as.cell_data_set(T_CD8, group.by = 'CD8_sub_MainCluster_paper')
cds <- cluster_cells(cds, cluster_method = 'louvain')
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

# integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
# cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$CD8_sub_MainCluster_paper =='T_Naive')[200]
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

detach("package:monocle3", unload = TRUE)

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
    plot_cells(cds, color_cells_by = "CD8_sub_MainCluster_paper", 
               label_groups_by_cluster = F, label_leaves = F, label_branch_points = F),
    width = 9.81, height = 6.19
)


# Monocle2-Seurat CD8
cds_all_CD8 <- my_create_monocle(T_CD8)
cds_all_CD8_seurat <- my_featureSelect_cds(
    cds_all_CD8,
    method = "seurat",
    seurat_obj = T_CD8
)
cds_all_CD8_seurat <- my_process_cds(cds_all_CD8_seurat)
my_plotPseudo(
    cds_all_CD8_seurat,
    seurat_obj = T_CD8,
    orig.ident = "orig.ident2"
)
plot_cell_trajectory(cds_all_CD8_seurat, color_by = 'State')
cds_all_CD8_seurat <- my_process_cds(cds_all_CD8_seurat, root_state = 5)
my_plotPseudo(
    cds_all_CD8_seurat,
    seurat_obj = T_CD8,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

T_CD8 <- my_AddSeuratPseudo(T_CD8, cds_all_CD8_seurat, reduction_key = 'Monocle_Seurat')
my_plotDim(T_CD8, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
my_plotDim(T_CD8, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(T_CD8, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')
ggsave(
    "T_CD8_Pseudotime_Monocle2Seurat.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, 
               group.by = 'CD8_sub_MainCluster_paper'),
    width = 9.81, height = 6.19
)


# Monocle2-monocle CD8
cds_all_CD8 <- my_create_monocle(T_CD8)
cds_all_CD8_monocle <- my_featureSelect_cds(
    cds_all_CD8,
    method = "monocle",
    seurat_obj = T_CD8,
    FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit"
)
cds_all_CD8_monocle <- my_process_cds(cds_all_CD8_monocle)
my_plotPseudo(
    cds_all_CD8_monocle,
    seurat_obj = T_CD8,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_CD8_monocle,
    seurat_obj = T_CD8,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

T_CD8 <- my_AddSeuratPseudo(T_CD8, cds_all_CD8_monocle, reduction_key = 'Monocle_Monocle')
my_plotDim(T_CD8, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
my_plotDim(T_CD8, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(T_CD8, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')
ggsave(
    "T_CD8_Pseudotime_Monocle2Monocle.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, 
               group.by = 'CD8_sub_MainCluster_paper'),
    width = 9.81, height = 6.19
)


# Monocle2-diffgenes CD8
cds_all_CD8 <- my_create_monocle(T_CD8)
cds_all_CD8_diffgenes <- my_featureSelect_cds(
    cds_all_CD8,
    method = "diffgenes",
    seurat_obj = T_CD8,
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)"
)
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
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

T_CD8 <- my_AddSeuratPseudo(T_CD8, cds_all_CD8_diffgenes, reduction_key = 'Monocle_Diffgenes')
my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')
ggsave(
    "T_CD8_Pseudotime_Monocle2Diffgenes.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, 
               group.by = 'CD8_sub_MainCluster_paper'),
    width = 9.81, height = 6.19
)


# Monocle2-dpFeature CD8
cds_all_CD8 <- my_create_monocle(T_CD8)
cds_all_CD8_dpFeature <- my_featureSelect_cds(
    cds_all_CD8,
    method = "dpFeature",
    seurat_obj = T_CD8,
    DefaultMeta = "CD8_sub_MainCluster_paper"
)
cds_all_CD8_dpFeature <- my_process_cds(cds_all_CD8_dpFeature)
my_plotPseudo(
    cds_all_CD8_dpFeature,
    seurat_obj = T_CD8,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_CD8_dpFeature,
    seurat_obj = T_CD8,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

T_CD8 <- my_AddSeuratPseudo(T_CD8, cds_all_CD8_dpFeature, reduction_key = 'Monocle_dpFeature')
my_plotDim(T_CD8, reduction = "Monocle_dpFeature", label = T, pt.size = 1, label.size = 5, group.by = 'CD8_sub_MainCluster_paper')
my_plotDim(T_CD8, reduction = "Monocle_dpFeature", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(T_CD8, reduction = "Monocle_dpFeature", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')
ggsave(
    "T_CD8_Pseudotime_Monocle2dpFeature.pdf",
    my_plotDim(T_CD8, reduction = "Monocle_dpFeature", label = T, pt.size = 1, label.size = 5, 
               group.by = 'CD8_sub_MainCluster_paper'),
    width = 9.81, height = 6.19
)


T_CD8_reductions <- T_CD8@reductions
if (!file.exists('T_CD8_reductions_220704.rda'))
    save(T_CD8_reductions, file = 'T_CD8_reductions_220704.rda')

###########################
Idents(T_sub_real) <- "T_sub_MainCluster"
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


# Monocle2-Seurat CD8_sub
cds_all_CD8sub <- my_create_monocle(CD8_sub)
cds_all_CD8sub_seurat <- my_featureSelect_cds(
    cds_all_CD8sub,
    method = "seurat",
    seurat_obj = CD8_sub
    # FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit"
)
cds_all_CD8sub_seurat <- my_process_cds(cds_all_CD8sub_seurat)
my_plotPseudo(
    cds_all_CD8sub_seurat,
    seurat_obj = CD8_sub,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_CD8sub_seurat,
    seurat_obj = CD8_sub,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

CD8_sub <- my_AddSeuratPseudo(CD8_sub, cds_all_CD8sub_seurat, reduction_key = 'Monocle_Seurat')
FeaturePlot(CD8_sub, 'Monocle_Seurat_Pseudotime')
# ggsave(filename = "CD8_sub_pseudotime_Monocle2Seurat.pdf", 
#        plot = FeaturePlot(CD8_sub, 'Monocle_Seurat_Pseudotime'), 
#        width = 9.81, height = 6.19)

my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Seurat_State')
my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Seurat_State')
my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'RNA_snn_res.0.5')
my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_paper')
my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')

my_plotPseudo(
    cds_all_CD8sub_seurat,
    seurat_obj = CD8_sub,
    color_by = 'State',
    orig.ident = "orig.ident2"
)

ggsave(
    "CD8_sub_State_Monocle2Seurat.pdf",
    my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Seurat_State'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_StatePseudo_Monocle2Seurat.pdf",
    plot_cell_trajectory(cds_all_CD8sub_seurat, color_by = 'State', size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_Pseudotime_Monocle2Seurat.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_paper'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeCluster_Monocle2Seurat.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'RNA_snn_res.0.5'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeMain_Monocle2Seurat.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeNaive_Monocle2Seurat.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Seurat", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster'),
    width = 9.81, height = 6.19
)



# Monocle2-Monocle CD8_sub
cds_all_CD8sub <- my_create_monocle(CD8_sub)
cds_all_CD8sub_monocle <- my_featureSelect_cds(
    cds_all_CD8sub,
    method = "monocle",
    seurat_obj = CD8_sub,
    FilterCondition = "mean_expression >= 0.4 & dispersion_empirical >= 0.9 * dispersion_fit"
)
cds_all_CD8sub_monocle <- my_process_cds(cds_all_CD8sub_monocle)
my_plotPseudo(
    cds_all_CD8sub_monocle,
    seurat_obj = CD8_sub,
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds_all_CD8sub_monocle,
    seurat_obj = CD8_sub,
    color_legend.position = c(.8, .2),
    color_by = 'CD8_sub_MainCluster_paper',
    orig.ident = "orig.ident2"
)

plot_cell_trajectory(cds_all_CD8sub_monocle, color_by = 'State')
cds_all_CD8sub_monocle <- my_process_cds(cds_all_CD8sub_monocle, root_state = 4)
CD8_sub <- my_AddSeuratPseudo(CD8_sub, cds_all_CD8sub_monocle, reduction_key = 'Monocle_Monocle')
FeaturePlot(CD8_sub, 'Monocle_Monocle_Pseudotime')
ggsave(filename = "CD8_sub_pseudotime_Monocle2Monocle.pdf",
       plot = FeaturePlot(CD8_sub, 'Monocle_Monocle_Pseudotime'),
       width = 9.81, height = 6.19)

my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Monocle_State')
my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Monocle_State')
my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'RNA_snn_res.0.5')
my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_paper')
my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2')
my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster')

my_plotPseudo(
    cds_all_CD8sub_monocle,
    seurat_obj = CD8_sub,
    color_by = 'State',
    orig.ident = "orig.ident2"
)

ggsave(
    "CD8_sub_State_Monocle2Monocle.pdf",
    my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Monocle_Monocle_State'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_StatePseudo_Monocle2Monocle.pdf",
    plot_cell_trajectory(cds_all_CD8sub_monocle, color_by = 'State', size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_Pseudotime_Monocle2Monocle.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_paper'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeCluster_Monocle2Monocle.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'RNA_snn_res.0.5'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeMain_Monocle2Monocle.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster2'),
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_sub_PseudotimeNaive_Monocle2Monocle.pdf",
    my_plotDim(CD8_sub, reduction = "Monocle_Monocle", label = T, pt.size = 1, label.size = 5, group.by = 'T_sub_MainCluster'),
    width = 9.81, height = 6.19
)


# monocle3
library(monocle3)
# k15 leiden no-resolution
cds <- as.cell_data_set(CD8_sub, group.by = 'RNA_snn_res.0.5')
cds <- cluster_cells(cds, k = 15, cluster_method = 'leiden')
# p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
# p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
# wrap_plots(p1, p2)

# integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
# cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$MainCluster_paper =='CD8_Tem')[100]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
# plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T,
#            label_branch_points = T, label_roots = T)
# plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = T, label_leaves = T,
#            label_branch_points = T, label_roots = T)

# Set the assay back as 'integrated'
integrated.sub <- as.Seurat(cds, assay = "integrated")
integrated.sub$monocle3_pseudotime <- ifelse(
    is.infinite(integrated.sub$monocle3_pseudotime),
    0,
    integrated.sub$monocle3_pseudotime
)
# FeaturePlot(integrated.sub, "monocle3_pseudotime")

# detach("package:monocle3", unload = TRUE)

plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
FeaturePlot(integrated.sub, "monocle3_pseudotime")


ggsave(
    'CD8_sub_pseudotime_Monocle3K15.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_sub_cluster_Monocle3K15.pdf',
    plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)

# k20 louvain
cds <- cluster_cells(cds, k = 20, cluster_method = 'louvain')
cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$MainCluster_paper =='CD8_Tem')[100]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
ggsave(
    'CD8_sub_pseudotime_Monocle3K20louvain.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_sub_cluster_Monocle3K20louvain.pdf',
    plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)

# k30 leiden
cds <- cluster_cells(cds, k = 30, cluster_method = 'leiden', resolution = 0.01)
cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$MainCluster_paper =='CD8_Tem')[100]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
ggsave(
    'CD8_sub_pseudotime_Monocle3K30.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_sub_cluster_Monocle3K30.pdf',
    plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)

# k35 leiden
cds <- cluster_cells(cds, k = 35, cluster_method = 'leiden', resolution = 0.01)
cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$MainCluster_paper =='CD8_Tem')[100]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
ggsave(
    'CD8_sub_pseudotime_Monocle3K35.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_sub_cluster_Monocle3K35.pdf',
    plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)

# k40 leiden
cds <- cluster_cells(cds, k = 40, cluster_method = 'leiden', resolution = 0.01)
cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
integrated.sub <- as.Seurat(cds, assay = "integrated")
max.avp <- which(integrated.sub$MainCluster_paper =='CD8_Tem')[100]
max.avp <- rownames(integrated.sub@meta.data)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
           label_branch_points = F, label_roots = F)
ggsave(
    'CD8_sub_pseudotime_Monocle3K40.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_sub_cluster_Monocle3K40.pdf',
    plot_cells(cds, color_cells_by = "MainCluster_paper", label_cell_groups = F, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1),
    width = 9.81, height = 6.19
)
