source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
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

# ref T cells
ref_TILAtlas <- readRDS("/mnt/e/Cryo-TCR/TILAtlas/ref_TILAtlas_mouse_v1.rds")
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref_TILAtlas, label = T, cols = refCols)

# my T cells
T_sub
query_T

# prepare data
query_T_real <- query_T
query_T_real$real <- sapply(
    rownames(query_T@meta.data),
    function(x) {
        x %in% rownames(T_sub@meta.data)
    }
)
Idents(query_T_real) <- "real"
query_T_real <- subset(query_T_real, ident = "TRUE")
if (all(rownames(T_sub@meta.data) == rownames(query_T_real@meta.data))) {
    query_T_real@meta.data <- T_sub@meta.data
}

# 1. query data plot
query_plot <- my_plotDim(
    query_T_real,
    reduction = "umap", label = T,
    pt.size = 1, label.size = 6,
    group.by = "IFN_Final",
    title = "T cells mapped by ProjecTILs"
)
query_plot <- query_plot +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "right",
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
            "CD4_Naive" = "#709ae1",
            "CD4_T_I-IFN" = "#197ec0",
            "Th1" = "#d2af81",
            "Th17" = "#bd559f",
            "Treg" = "#f7c530",
            "Tfh" = "#46732e",
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e",
            "T_others" = "#748aa6"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
query_plot
ggsave(
    filename = "plot/scRNA3_queryT_dim.pdf",
    plot = query_plot,
    height = 6, width = 7.69 + 2
)

# 2. reference data plot
ref_TILAtlas$functional.cluster <- factor(
    ref_TILAtlas$functional.cluster,
    levels = c(
        "CD4_NaiveLike",
        "Th1",
        "Treg",
        "Tfh",
        "CD8_NaiveLike",
        "CD8_EarlyActiv",
        "CD8_EffectorMemory",
        "CD8_Tpex",
        "CD8_Tex"
    )
)
Idents(ref_TILAtlas) <- "functional.cluster"
ref_cols <- c(
    "CD4_NaiveLike" = "#709ae1",
    "Th1" = "#d2af81",
    "Treg" = "#f7c530",
    "Tfh" = "#46732e",
    "CD8_NaiveLike" = "#d5e4a2",
    "CD8_EarlyActiv" = "#F8766D", #
    "CD8_EffectorMemory" = "#95cc5e",
    "CD8_Tpex" = "#A58AFF", #
    "CD8_Tex" = "#d0dfe6"
)

ref_plot <- my_plotDim(
    ref_TILAtlas,
    reduction = "umap", label = T,
    pt.size = 1, label.size = 6,
    # group.by = 'IFN_Final',
    title = "Reference T cells from ProjecTILs"
)
ref_plot <- ref_plot +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "right",
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = ref_cols
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
ref_plot
ggsave(
    filename = "plot/scRNA3_refT_dim.pdf",
    plot = ref_plot,
    height = 6, width = 7.69 + 2
)


# 3. reference data markers
aref <- list()
markers48_list <- c(
    "Cd4", "Cd8a", "Sell", "Mki67", "Gzma", "Il7r", "Cd44", "Ctla4", "Pdcd1",
    "Il21r", "Tbx21", "Cxcr3", "Rorc", "Il17a", "Foxp3", "Il2ra"
)
for (i in markers48_list) {
    aref[[i]] <- FeaturePlot(
        ref_TILAtlas,
        features = i, reduction = "umap",
        pt.size = 0.5, cols = c("#e0e7c8", "#cc340a")
    ) +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            text = element_blank(),
            plot.title = element_blank()
        ) +
        annotate(
            "text",
            x = 4.2,
            y = -7,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
refT_express <- patchwork::wrap_plots(plots = aref, ncol = 4)
ggsave(
    filename = "plot/scRNA3_refT_markers_umap.pdf",
    plot = refT_express,
    height = 6 * 2, width = 7.69 * 1.5 * 4 / 5
)

# 4. my T data markers
T_sub@reductions$umap <- query_T_real@reductions$umap
aque <- list()
markers48_list <- c(
    "Cd4", "Cd8a", "Sell", "Mki67", "Gzma", "Il7r", "Cd44", "Ctla4", "Pdcd1",
    "Il21r", "Tbx21", "Cxcr3", "Rorc", "Il17a", "Foxp3", "Il2ra"
)
for (i in markers48_list) {
    aque[[i]] <- FeaturePlot(
        T_sub,
        features = i, reduction = "umap",
        pt.size = 0.5, cols = c("#e0e7c8", "#cc340a")
    ) +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            text = element_blank(),
            plot.title = element_blank()
        ) +
        annotate(
            "text",
            x = 4,
            y = -7,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
queT_express <- patchwork::wrap_plots(plots = aque, ncol = 4)
ggsave(
    filename = "plot/scRNA3_queryT_markers_umap.pdf",
    plot = queT_express,
    height = 6 * 2, width = 7.69 * 1.5 * 4 / 5
)
