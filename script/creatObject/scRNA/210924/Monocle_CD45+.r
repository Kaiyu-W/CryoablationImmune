source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load(file = "Combined_analysis_Myeloid_real.rda")

# pseudotime: Monocole2
# all Myeloid

cds_all <- my_create_monocle(Myeloid_sub_real)
cds_all_seurat <- my_featureSelect_cds(
    cds_all,
    method = "seurat",
    seurat_obj = Myeloid_sub_real
)
cds_all_monocle <- my_featureSelect_cds(
    cds_all,
    method = "monocle",
    seurat_obj = Myeloid_sub_real,
    FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit"
)
cds_all_diffgenes <- my_featureSelect_cds(
    cds_all,
    method = "diffgenes",
    seurat_obj = Myeloid_sub_real,
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)"
)
cds_all_dpFeature <- my_featureSelect_cds(
    cds_all,
    method = "dpFeature",
    seurat_obj = Myeloid_sub_real,
    DefaultMeta = "MainCluster_new"
)

# cds_dpFeature
cds_all_dpFeature <- my_process_cds(cds_all_dpFeature)
my_plotPseudo(
    cds_all_dpFeature,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_seurat
cds_all_seurat <- my_process_cds(cds_all_seurat)
my_plotPseudo(
    cds_all_seurat,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_monocle
cds_all_monocle <- my_process_cds(cds_all_monocle)
my_plotPseudo(
    cds_all_monocle,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_diffgenes
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes)
plot_cell_trajectory(cds_all_diffgenes, color_by = "State")
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes, root_state = 2)
my_plotPseudo(
    cds_all_diffgenes,
    color_by = "Pseudotime",
    theta = -20,
    show_branch_points = F,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )
plot_cell_trajectory(
    cds_all_diffgenes,
    color_by = "MainCluster2_new",
    theta = -20,
    show_branch_points = F,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )


# pseudotime: Monocole2
# only Macro/Mono/DC
Idents(Myeloid_sub_real) <- "MainCluster_new"
cds_sub <- my_create_monocle(Myeloid_sub_real, idents = paste0("C", 2:7))

cds_sub_seurat <- my_featureSelect_cds(
    cds_sub,
    method = "seurat",
    seurat_obj = Myeloid_sub_real
    )
cds_sub_monocle <- my_featureSelect_cds(
    cds_sub,
    method = "monocle",
    seurat_obj = Myeloid_sub_real,
    FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit"
    )
cds_sub_diffgenes <- my_featureSelect_cds(
    cds_sub,
    method = "diffgenes",
    seurat_obj = Myeloid_sub_real,
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 0.5 | avg_log2FC < -0.5)"
    )
cds_sub_dpFeature <- my_featureSelect_cds(
    cds_sub,
    method = "dpFeature",
    seurat_obj = Myeloid_sub_real,
    DefaultMeta = "MainCluster_new",
    qval_threshold = 0.01
    )

# cds_dpFeature
cds_sub_dpFeature <- my_process_cds(cds_sub_dpFeature)
my_plotPseudo(
    cds_sub_dpFeature,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_seurat
cds_sub_seurat <- my_process_cds(cds_sub_seurat)
my_plotPseudo(
    cds_sub_seurat,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_monocle
cds_sub_monocle <- my_process_cds(cds_sub_monocle)
my_plotPseudo(
    cds_sub_monocle,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )

# cds_diffgenes
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes)
plot_cell_trajectory(
    cds_sub_diffgenes,
    color_by = "State",
    theta = -20,
    show_branch_points = F
    )
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes, root_state = 2)
my_plotPseudo(
    cds_sub_diffgenes,
    color_by = "Pseudotime",
    theta = -20,
    show_branch_points = F,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )
plot_cell_trajectory(
    cds_sub_diffgenes,
    color_by = "MainCluster2_new",
    theta = -20,
    show_branch_points = F,
    seurat_obj = Myeloid_sub_real,
    orig.ident = "orig.ident2"
    )


Seurat_tmp <- as.Seurat(cds_sub_diffgenes)
Myeloid_sub_real$Pseudotime <- Seurat_tmp$Pseudotime
Myeloid_sub_real$State <- unfactor(Seurat_tmp$State)
Myeloid_sub_real$State[is.na(Myeloid_sub_real$State)] <- 4
DimPlot(
    Myeloid_sub_real,
    reduction = "umap",
    group.by = "State",
    pt.size = 1
    )
DimPlot(
    Myeloid_sub_real,
    reduction = "umap",
    group.by = "Pseudotime",
    pt.size = 1
    ) +
    theme(legend.position = "none")

# getwd()
# ls()[grep("^cds", ls())]
save(list = ls()[grep("^cds", ls())], file = "monocle_cds_Myeloid.rda")

# pseudotime: Monocole2
# all T

load(file = "Combined_analysis_TNK_real.rda")
Idents(T_sub_real) <- "T_sub_MainCluster2"
T_T <- subset(
    T_sub_real, 
    ident = setdiff(
        levels(T_sub_real$T_sub_MainCluster2), 
        "NK"
        )
    )

cds_all <- my_create_monocle(T_T)
cds_all_seurat <- my_featureSelect_cds(
    cds_all,
    method = "seurat",
    seurat_obj = T_T
)
cds_all_monocle <- my_featureSelect_cds(
    cds_all,
    method = "monocle",
    seurat_obj = T_T,
    FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit"
)
deg.cluster <- Seurat::FindAllMarkers(T_T)
cds_all_diffgenes <- my_featureSelect_cds(
    cds_all,
    method = "diffgenes",
    seurat_obj = T_T,
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)"
)
# cds_all_dpFeature <- my_featureSelect_cds(
#     cds_all,
#     method = "dpFeature",
#     seurat_obj = T_T,
#     DefaultMeta = "MainCluster_new"
# )

# # cds_dpFeature
# cds_all_dpFeature <- my_process_cds(cds_all_dpFeature)
# my_plotPseudo(
#     cds_all_dpFeature,
#     seurat_obj = T_T,
#     orig.ident = "orig.ident2"
#     )

# cds_seurat
cds_all_seurat <- my_process_cds(cds_all_seurat)
my_plotPseudo(
    cds_all_seurat,
    seurat_obj = T_T,
    color_by = "T_sub_MainCluster2",
    show_branch_points = F,
    orig.ident = "orig.ident2",
    color_legend.position = c(.8, .2),
    color_guide_legend_direction = "horizontal",
    color_guide_legend_nrow = 4,
    color_guide_legend_ncol = NULL,
    )
plot_cell_trajectory(
    cds_all_seurat,
    color_by = "T_sub_MainCluster2",
    # theta = -20,
    show_branch_points = F,
    seurat_obj = T_T,
    orig.ident = "orig.ident2"
    )

# cds_monocle
cds_all_monocle <- my_process_cds(cds_all_monocle)
my_plotPseudo(
    cds_all_monocle,
    seurat_obj = T_T,
    color_by = "T_sub_MainCluster2",
    orig.ident = "orig.ident2"
    )
plot_cell_trajectory(
    cds_all_monocle,
    color_by = "T_sub_MainCluster2",
    theta = -20,
    show_branch_points = F,
    seurat_obj = T_T,
    orig.ident = "orig.ident2"
    )

# cds_diffgenes
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes)
plot_cell_trajectory(cds_all_diffgenes, color_by = "State")
plot_cell_trajectory(cds_all_diffgenes, color_by = "T_sub_MainCluster2")
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes, root_state = 6)
my_plotPseudo(
    cds_all_diffgenes,
    color_by = "Pseudotime",
    theta = -20,
    show_branch_points = F,
    seurat_obj = T_T,
    orig.ident = "orig.ident2"
    )
plot_cell_trajectory(
    cds_all_diffgenes,
    color_by = "T_sub_MainCluster2",
    theta = -20,
    show_branch_points = F,
    seurat_obj = T_T,
    orig.ident = "orig.ident2"
    )
