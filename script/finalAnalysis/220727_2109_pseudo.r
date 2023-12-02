source("/mnt/c/Users/PC/Documents/GitHub/My_scRNA_pipeline/utilities.r")
setwd("~/Desktop/Cryo/")

TNK <- readRDS("TNK_3end.rds")
TNK

my_plotDim(TNK, group.by='IFN_Final', reduction = "umap", label = T, pt.size = 0.5, label.size = 5)

if(!dir.exists("T_Pseudotime")) dir.create("T_Pseudotime")
setwd("T_Pseudotime")

T_anno_dim <- my_plotDim(TNK, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, 
                         group.by = 'IFN_Final')
T_anno2_dim <- my_plotDim(TNK, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, 
                             group.by = 'IFN_Final_Cluster')

ggsave(filename = "T_anno_dim.pdf", plot = T_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "T_anno2_dim.pdf", plot = T_anno2_dim, width = 9.81, height = 6.19)

## Sub-Cluster CD8
Idents(TNK) <- "IFN_Final"
CD8_sub <- subset(TNK, ident = levels(TNK)[grep("^CD8", levels(TNK))])
CD8_sub <- my_process_seurat(
    CD8_sub, nVariableFeatures = 2000, normalize = T, norm.method = 'Log', 
    default.assay = 'RNA', mt.to.regress = F, tsne = F, dims_umap = 1:15)
CD8_sub <- RunUMAP(CD8_sub, dims = 1:15, n.neighbors = 20)
CD8_sub_anno_dim <- my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, 
                               label.size = 5, group.by = 'IFN_Final')
CD8_sub_dim <- my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, 
                          label.size = 5, group.by = 'IFN_Final_Cluster')

ggsave(filename = "CD8_sub_anno_dim.pdf", plot = CD8_sub_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "CD8_sub_dim.pdf", plot = CD8_sub_dim, width = 9.81, height = 6.19)

# 3d plot
CD8_sub <- RunUMAP(CD8_sub, dims = 1:15, n.neighbors = 20, n.components = 3)
plotly::plot_ly(
    data = as.data.frame(CD8_sub@reductions$umap@cell.embeddings), 
    x = ~UMAP_1, y= ~UMAP_2, z = ~UMAP_3,
    color = CD8_sub$IFN_Final, size = 1
    )

# Monocle3 for CD8
cds <- my_process_monocle3(
    seurat_obj = CD8_sub, group.by = 'IFN_Final', root_ident = 'CD8_Naive',
    cluster_method = 'leiden', partition_qval = 1, resolution = 0.01)

require(monocle3)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
           trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F')
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5,
           trajectory_graph_segment_size = 1, show_trajectory_graph = F)#, trajectory_graph_color = 'turquoise')

# sub-cluster for Cryo/NonCryo
CD8_sub_anno_dim2 <- my_plotDim(CD8_sub, reduction = "umap", label = T, pt.size = 1, 
           label.size = 5, group.by = 'IFN_Final', split.by = 'orig.ident2')
ggsave(filename = "CD8_sub_anno_dim_split.pdf", plot = CD8_sub_anno_dim2, width = 9.81*2, height = 6.19)

Idents(CD8_sub) <- 'orig.ident2'
CD8_sub_Cryo <- subset(CD8_sub, ident = 'Cryo')
CD8_sub_NonCryo <- subset(CD8_sub, ident = 'NonCryo')

cds_Cyro <- my_process_monocle3(
    seurat_obj = CD8_sub_Cryo, group.by = 'IFN_Final', root_ident = 'CD8_Naive',
    k = 25)
cds_NonCyro <- my_process_monocle3(
    seurat_obj = CD8_sub_NonCryo, group.by = 'IFN_Final', root_ident = 'CD8_Naive')

plot_cells(cds_Cyro, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
           trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F')
plot_cells(cds_NonCyro, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
           trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F')

# save
# all
ggsave(
    'CD8_all_cluster_Monocle3.pdf',
    plot_cells(cds, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
               trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F'),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_all_pseudotime_Monocle3.pdf',
    plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5,
               trajectory_graph_segment_size = 1, show_trajectory_graph = F),
    width = 9.81, height = 6.19
)
# Cryo
ggsave(
    'CD8_Cryo_cluster_Monocle3.pdf',
    plot_cells(cds_Cyro, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
               trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F'),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_Cryo_pseudotime_Monocle3.pdf',
    plot_cells(cds_Cyro, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5,
               trajectory_graph_segment_size = 1, show_trajectory_graph = F),
    width = 9.81, height = 6.19
)
# NonCryo
ggsave(
    'CD8_NonCryo_cluster_Monocle3.pdf',
    plot_cells(cds_NonCyro, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
               trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F'),
    width = 9.81, height = 6.19
)
ggsave(
    'CD8_NonCryo_pseudotime_Monocle3.pdf',
    plot_cells(cds_NonCyro, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5,
               trajectory_graph_segment_size = 1, show_trajectory_graph = F),
    width = 9.81, height = 6.19
)


# Monocle2 for CD8
Idents(CD8_sub) <- 'IFN_Final'
cds2_CD8 <- my_create_monocle(CD8_sub)

cds2_CD8_seurat <- my_featureSelect_cds(
    cds2_CD8,
    method = 'seurat',
    ntop = 1000,
    seurat_obj = CD8_sub
    # method = "diffgenes",
    # # FilterCondition = "p_val_adj<1e-6"
    # FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1.5)"
)
sum(cds2_CD8_seurat@featureData@data[["use_for_ordering"]])

cds2_CD8_seurat <- my_process_cds(cds2_CD8_seurat)
my_plotPseudo(
    cds2_CD8_seurat,
    seurat_obj = CD8_sub,
    color_legend.position = c(.6, .2),
    color_by = 'IFN_Final',
    orig.ident = "orig.ident2"
)
my_plotPseudo(
    cds2_CD8_seurat,
    seurat_obj = CD8_sub,
    color_legend.position = c(.6, .2),
    color_by = 'orig.ident2',
    orig.ident = "orig.ident2"
)

CD8_sub <- my_AddSeuratPseudo(CD8_sub, cds2_CD8_seurat, reduction_key = 'Monocle_Diffgenes')
# my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final')
# my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'orig.ident2', 
#            cells = grep("I-IFN", CD8_sub$IFN_Final))
# my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, label.size = 5, group.by = 'IFN_Final_Cluster', 
           # cells = grep("CD8_Tex", CD8_sub$IFN_Final))
CD8_sub_splitMonocle2 <- my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, 
                                    label.size = 5, split.by = 'IFN_Final_Cluster', group.by = 'IFN_Final')
CD8_sub_split2Monocle2 <- my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, 
                                    label.size = 5, split.by = 'IFN_Final', group.by = 'IFN_Final_Cluster')
CD8_sub_clusterMonocle2 <- my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, 
                                      label.size = 5, group.by = 'IFN_Final')
CD8_sub_cluster2Monocle2 <- my_plotDim(CD8_sub, reduction = "Monocle_Diffgenes", label = T, pt.size = 1, 
                                       label.size = 5, group.by = 'IFN_Final_Cluster')

ggsave(
    "CD8_all_clusterSplit_Monocle2.pdf",
    CD8_sub_splitMonocle2,
    width = 9.81*2, height = 6.19
)
ggsave(
    "CD8_all_clusterSplit2_Monocle2.pdf",
    CD8_sub_split2Monocle2,
    width = 9.81*2, height = 6.19
)
ggsave(
    "CD8_all_cluster_Monocle2.pdf",
    CD8_sub_clusterMonocle2,
    width = 9.81, height = 6.19
)
ggsave(
    "CD8_all_cluster2_Monocle2.pdf",
    CD8_sub_cluster2Monocle2,
    width = 9.81, height = 6.19
)


## Sub-Cluster CD4
Idents(TNK) <- "IFN_Final"
CD4_sub <- subset(TNK, ident = levels(TNK)[
    !grepl("^CD8", levels(TNK)) & 
        !grepl("^NK", levels(TNK)) & 
        !grepl("T_others", levels(TNK))
    ]
)
CD4_sub <- my_process_seurat(
    CD4_sub, nVariableFeatures = 2000, normalize = T, norm.method = 'Log', 
    default.assay = 'RNA', mt.to.regress = F, tsne = F)
CD4_sub_anno_dim <- my_plotDim(CD4_sub, reduction = "umap", label = T, pt.size = 1, 
                               label.size = 5, group.by = 'IFN_Final')
CD4_sub_dim <- my_plotDim(CD4_sub, reduction = "umap", label = T, pt.size = 1, 
                          label.size = 5, group.by = 'IFN_Final_Cluster')

ggsave(filename = "CD4_sub_anno_dim.pdf", plot = CD4_sub_anno_dim, width = 9.81, height = 6.19)
ggsave(filename = "CD4_sub_dim.pdf", plot = CD4_sub_dim, width = 9.81, height = 6.19)

# Monocle3 for Tfh
Idents(CD4_sub) <- 'IFN_Final'
CD4_Tfh <- subset(CD4_sub, ident = c('CD4_Naive', 'Tfh', 'Tfh_I-IFN'))
CD4_Tfh <- my_process_seurat(
    CD4_Tfh, nVariableFeatures = 2000, normalize = T, norm.method = 'Log', 
    default.assay = 'RNA', mt.to.regress = F, tsne = F)
CD4_Tfh_anno_dim <- my_plotDim(CD4_Tfh, reduction = "umap", label = T, pt.size = 1, 
                               label.size = 5, group.by = 'IFN_Final')
cds_Tfh <- my_process_monocle3(
    seurat_obj = CD4_Tfh, group.by = 'IFN_Final', root_ident = 'CD4_Naive'
    )

plot_cells(cds_Tfh, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
           label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
           trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F')
ggsave(
    'CD4_Tfh_cluster_Monocle3.pdf',
    plot_cells(cds_Tfh, color_cells_by = "cluster", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5, 
               trajectory_graph_segment_size = 1, trajectory_graph_color = '#4F4F4F'),
    width = 9.81, height = 6.19
)
ggsave(
    'CD4_Tfh_pseudotime_Monocle3.pdf',
    plot_cells(cds_Tfh, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = F,
               label_branch_points = F, label_roots = F, cell_size = 1, group_label_size = 5,
               trajectory_graph_segment_size = 1, show_trajectory_graph = F),
    width = 9.81, height = 6.19
)


# Monocle2 for Tfh
Idents(CD4_Tfh) <- 'IFN_Final'
cds2_Tfh <- my_create_monocle(CD4_Tfh)

cds2_Tfh_dg <- my_featureSelect_cds(
    cds2_Tfh,
    seurat_obj = CD4_Tfh,
    method = "diffgenes",
    # FilterCondition = "p_val_adj<1e-6"
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)"
)
sum(cds2_Tfh_dg@featureData@data[["use_for_ordering"]])

cds2_Tfh_dg <- my_process_cds(cds2_Tfh_dg)
my_plotPseudo(
    cds2_Tfh_dg,
    seurat_obj = CD4_Tfh,
    color_legend.position = c(.6, .2),
    color_by = 'IFN_Final',
    orig.ident = NULL
)
plot_cell_trajectory(
    cds2_Tfh_dg,
    color_by = 'Pseudotime'
)

CD4_Tfh <- my_AddSeuratPseudo(CD4_Tfh, cds2_Tfh_dg, reduction_key = 'Monocle_Diffgenes')
CD4_Tfh_splitMonocle2 <- my_plotDim(CD4_Tfh, reduction = "Monocle_Diffgenes", label = T, pt.size = 1.5, 
                                      label.size = 5, group.by = 'IFN_Final', split.by = 'orig.ident2')
CD4_Tfh_clusterMonocle2 <- my_plotDim(CD4_Tfh, reduction = "Monocle_Diffgenes", label = T, pt.size = 1.5, 
                                      label.size = 5, group.by = 'IFN_Final')

ggsave(
    "CD4_Tfh_clusterSplit_Monocle2.pdf",
    CD4_Tfh_splitMonocle2,
    width = 9.81*1.5, height = 6.19
)
ggsave(
    "CD4_Tfh_Psuedotime_Monocle2.pdf",
    plot_cell_trajectory(
        cds2_Tfh_dg,
        color_by = 'Pseudotime'
    ),
    width = 9.81, height = 6.19
)
ggsave(
    "CD4_Tfh_cluster_Monocle2.pdf",
    CD4_Tfh_clusterMonocle2,
    width = 9.81, height = 6.19
)
