# source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
source("/mnt/c/Users/PC/Documents/GitHub/My_scRNA_pipeline/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis_TNK_real.rda")
load("T_sub_real_meta_220713.rda")
T_sub_real@meta.data <- T_meta_220713
Idents(T_sub_real) <- "IFN_Final"
T_sub <- subset(T_sub_real, idents = setdiff(unique(T_sub_real$IFN_Final), c("NK", "NK_Proliferative")))
T_sub$IFN_Final <- unfactor(T_sub$IFN_Final)
T_sub$IFN_Final[T_sub$IFN_Final == "Tfh_I-IFN"] <- "CD4_T_I-IFN"
Idents(T_sub) <- "IFN_Final"

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

# color <- 'viridis'

# CD8-T cells
CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")
Idents(T_sub) <- "IFN_Final"
CD8_sub <- subset(T_sub, ident = CD8_idents)
CD8_sub <- my_process_seurat(
    CD8_sub,
    nVariableFeatures = 2000,
    normalize = T,
    norm.method = "Log",
    default.assay = "RNA",
    mt.to.regress = F,
    tsne = F,
    dims_umap = 1:15
)
CD8_sub <- RunUMAP(CD8_sub, dims = 1:15, n.neighbors = 20)

# 1.umap
CD8_sub_anno_dim <- my_plotDim(
    CD8_sub,
    reduction = "umap",
    label = T,
    pt.size = 1,
    label.size = 5,
    group.by = "IFN_Final",
    title = "CD8-T cells MainCluster"
) +
    theme(
        legend.position = c(0.7, 0.9),
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
ggsave(
    filename = "plot/scRNA3_TCD8_dim.pdf",
    plot = CD8_sub_anno_dim,
    height = 6, width = 7.69
)

# 2.Monocle3 for CD8
cds <- my_process_monocle3(
    seurat_obj = CD8_sub,
    group.by = "IFN_Final",
    root_ident = "CD8_Naive",
    cluster_method = "leiden",
    partition_qval = 1,
    resolution = 0.01
)

require(monocle3)
# cluster
CD8_monocle3_cluster <- plot_cells(
    cds,
    color_cells_by = "cluster",
    label_cell_groups = T,
    label_leaves = F,
    label_branch_points = F,
    label_roots = F,
    cell_size = 1,
    group_label_size = 5,
    trajectory_graph_segment_size = 1,
    trajectory_graph_color = "#4F4F4F"
) +
    labs(title = "CD8-T cells Trajectories") +
    theme(
        # legend.position = c(0.75,0.9),
        legend.position = "none",
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e"
        )
        # ) +
        # guides(
        #     color = guide_legend(
        #         override.aes = list(size=4)
        #     )
    )
ggsave(
    filename = "plot/scRNA3_TCD8_monocle3_cluster.pdf",
    plot = CD8_monocle3_cluster,
    height = 6, width = 7.69
)

# pseudotime
CD8_monocle3_pseudotime <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = T,
    label_leaves = F,
    label_branch_points = F,
    label_roots = F,
    cell_size = 1,
    group_label_size = 5,
    trajectory_graph_segment_size = 1,
    show_trajectory_graph = F,
    # , trajectory_graph_color = 'turquoise')
) +
    labs(title = "Pseudotime") +
    theme(
        legend.position = c(0.9, 0.9),
        # legend.position = 'none',
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        text = element_text(face = "bold", size = 12)
    ) + scale_color_viridis_c(
        option = "plasma" # 'viridis', 'plasma', 'cividis', 'inferno'
    )
CD8_monocle3_pseudotime$labels$colour <- "Value"
ggsave(
    filename = "plot/scRNA3_TCD8_monocle3_pseudotime.pdf",
    plot = CD8_monocle3_pseudotime,
    height = 6, width = 7.69
)

# sub-cluster for Cryo/NonCryo
Idents(CD8_sub) <- "orig.ident2"
CD8_sub_Cryo <- subset(CD8_sub, ident = "Cryo")
CD8_sub_NonCryo <- subset(CD8_sub, ident = "NonCryo")

cds_Cyro <- my_process_monocle3(
    seurat_obj = CD8_sub_Cryo,
    group.by = "IFN_Final",
    root_ident = "CD8_Naive",
    k = 25
)
cds_NonCyro <- my_process_monocle3(
    seurat_obj = CD8_sub_NonCryo,
    group.by = "IFN_Final",
    root_ident = "CD8_Naive"
)

CD8_monocle3_cluster_CA <- plot_cells(
    cds_Cyro,
    color_cells_by = "cluster",
    label_cell_groups = T,
    label_leaves = F,
    label_branch_points = F,
    label_roots = F,
    cell_size = 1,
    group_label_size = 5,
    trajectory_graph_segment_size = 1,
    trajectory_graph_color = "#4F4F4F"
) +
    labs(title = "CA's CD8-T cells Trajectories") +
    theme(
        # legend.position = c(0.75,0.9),
        legend.position = "none",
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e"
        )
        # ) +
        # guides(
        #     color = guide_legend(
        #         override.aes = list(size=4)
        #     )
    )
ggsave(
    filename = "plot/scRNA3_TCD8_CA_monocle3_cluster.pdf",
    plot = CD8_monocle3_cluster_CA,
    height = 6, width = 7.69
)

CD8_monocle3_cluster_NonCA <- plot_cells(
    cds_NonCyro,
    color_cells_by = "cluster",
    label_cell_groups = T,
    label_leaves = F,
    label_branch_points = F,
    label_roots = F,
    cell_size = 1,
    group_label_size = 5,
    trajectory_graph_segment_size = 1,
    trajectory_graph_color = "#4F4F4F"
) +
    labs(title = "Non-CA's CD8-T cells Trajectories") +
    theme(
        # legend.position = c(0.75,0.9),
        legend.position = "none",
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e"
        )
        # ) +
        # guides(
        #     color = guide_legend(
        #         override.aes = list(size=4)
        #     )
    )
ggsave(
    filename = "plot/scRNA3_TCD8_NonCA_monocle3_cluster.pdf",
    plot = CD8_monocle3_cluster_NonCA,
    height = 6, width = 7.69
)
