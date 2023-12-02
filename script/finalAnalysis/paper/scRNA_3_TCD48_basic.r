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

# 2. all T cells
T_sub
if (!file.exists("scRNA3_T_cellbarcode.csv")) {
    write.csv(colnames(T_sub), file = "scRNA3_T_cellbarcode.csv")
}
if (!file.exists("scRNA3_T_umap.csv")) {
    write.csv(T_sub@reductions$umap@cell.embeddings, file = "scRNA3_T_umap.csv")
}
if (!file.exists("scRNA3_T_meta.csv")) {
    write.csv(T_sub@meta.data, file = "scRNA3_T_meta.csv")
}

# 2.1 T umap
T_plot <- my_plotDim(
    T_sub,
    reduction = "umap", label = T,
    pt.size = 2, label.size = 6, repel = T,
    group.by = "IFN_Final",
    title = "T cells MainCluster"
)
T_plot <- T_plot +
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
T_plot
ggsave(
    filename = "plot/scRNA3_T_dim.pdf",
    plot = T_plot,
    height = 6, width = 7.69 + 2
)

my_CountCluster(T_sub, group2 = "orig.ident2", group1 = "active.ident")

# 2.2 T markers violin
CD4_makers <- list(
    Helper = "Cd4",
    Naive = "Sell",
    Th1 = c("Tbx21", "Cxcr3"),
    Th17 = c("Rorc", "Il17a"),
    T_FH = c("Il21r", "Pdcd1"),
    Treg = c("Foxp3", "Il2ra")
)
CD8_makers <- list(
    Cytotoxic = c("Cd8a", "Gzma"),
    Naive = c("Sell"),
    Memory = c("Il7r"),
    Effector = c("Cd44"),
    Exhausted = c("Pdcd1", "Ctla4"),
    Proliferative = c("Mki67")
)

CD4_idents <- c("CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh")
CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")
T_CD4_dot <- my_DotPlot_split(
    T_sub,
    features = CD4_makers,
    idents = CD4_idents
) + RotatedAxis()
T_CD8_dot <- my_DotPlot_split(
    T_sub,
    features = CD8_makers,
    idents = CD8_idents
) + RotatedAxis()

T_CD4_violin <- VlnPlot(
    T_sub,
    idents = CD4_idents,
    feature = unlist(CD4_makers),
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    cols = c(
        "CD4_Naive" = "#709ae1",
        "CD4_T_I-IFN" = "#197ec0",
        "Th1" = "#d2af81",
        "Th17" = "#bd559f",
        "Treg" = "#f7c530",
        "Tfh" = "#46732e"
    )
) +
    # ggsci::scale_fill_jco() +
    theme(
        aspect.ratio = 5,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.line.x = element_blank(),
        strip.text.x = element_text(angle = 45, size = 15, face = "bold"),
        # panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )
ggsave(
    filename = "plot/scRNA3_TCD4_markers_violin.pdf",
    plot = T_CD4_violin,
    height = 6, width = 7.69 + 1
)


T_CD8_violin <- VlnPlot(
    T_sub,
    idents = CD8_idents,
    feature = unlist(CD8_makers),
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    cols = c(
        "CD8_Naive" = "#d5e4a2",
        "CD8_Tem_I-IFN" = "#ff410d",
        "CD8_Tem" = "#95cc5e",
        "CD8_Tex" = "#d0dfe6",
        "CD8_Tex_Proliferative" = "#f79d1e"
    )
) +
    # ggsci::scale_fill_jco() +
    theme(
        aspect.ratio = 5,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.line.x = element_blank(),
        strip.text.x = element_text(angle = 45, size = 15, face = "bold"),
        # panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )
ggsave(
    filename = "plot/scRNA3_TCD8_markers_violin.pdf",
    plot = T_CD8_violin,
    height = 4.5, width = 7.69 + 1
)

# 2.3 T markers express
T_Marker_rough <- append(CD4_makers, CD8_makers)

a8 <- list()
markers8_list <- c("Cd8a", "Sell", "Il7r", "Pdcd1", "Gzma", "Mki67", "Cd44", "Ctla4")
for (i in markers8_list) {
    a8[[i]] <- FeaturePlot(
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
            x = 4.2,
            y = -5,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
TCD8_express <- patchwork::wrap_plots(plots = a8, ncol = 4)
ggsave(
    filename = "plot/scRNA3_TCD8_markers_umap.pdf",
    plot = TCD8_express,
    height = 6, width = 7.69 * 1.5 * 4 / 5
)

a4 <- list()
markers4_list <- c("Cd4", "Tbx21", "Rorc", "Il21r", "Foxp3", "Sell", "Cxcr3", "Il17a", "Pdcd1", "Il2ra")
for (i in markers4_list) {
    a4[[i]] <- FeaturePlot(
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
            x = 4.2,
            y = -5,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
TCD4_express <- patchwork::wrap_plots(plots = a4, ncol = 5)
ggsave(
    filename = "plot/scRNA3_TCD4_markers_umap.pdf",
    plot = TCD4_express,
    height = 6, width = 7.69 * 1.5
)

a48 <- list()
markers48_list <- c(
    "Cd4", "Cd8a", "Sell", "Mki67", "Gzma", "Il7r", "Cd44", "Ctla4", "Pdcd1",
    "Il21r", "Tbx21", "Cxcr3", "Rorc", "Il17a", "Foxp3", "Il2ra"
)
for (i in markers48_list) {
    a48[[i]] <- FeaturePlot(
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
            x = 4.2,
            y = -5,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
TCD48_express <- patchwork::wrap_plots(plots = a48, ncol = 4)
ggsave(
    filename = "plot/scRNA3_TCD48_markers_umap.pdf",
    plot = TCD48_express,
    height = 6 * 2, width = 7.69 * 1.5 * 4 / 5
)

# 2.4 T celltype count
colors_bar <- c(
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
T_sub_count_sum <- as.data.frame(my_CountCluster(T_sub, group1 = "IFN_Final", group2 = "orig.ident2")[-13, 1:2])
T_sub_count_sum2 <- my_CountCluster(T_sub, group1 = "IFN_Final", group2 = "orig.ident2")["sum", 1:2]
T_sub_count_sum[, 1] <- T_sub_count_sum[, 1] / T_sub_count_sum2[1]
T_sub_count_sum[, 2] <- T_sub_count_sum[, 2] / T_sub_count_sum2[2]
T_sub_count_sum[, 1] <- T_sub_count_sum[, 1] / sum(T_sub_count_sum[, 1])
T_sub_count_sum[, 2] <- T_sub_count_sum[, 2] / sum(T_sub_count_sum[, 2])
T_sub_count_sum
T_sub_count_sum$cell <- factor(rownames(T_sub_count_sum), levels = names(colors_bar))
T_sub_count_sum <- reshape2::melt(T_sub_count_sum, id.vars = "cell")
T_sub_count_sum$variable <- sub("_count$", "", T_sub_count_sum$variable)
T_sub_count_sum$variable <- ifelse(T_sub_count_sum$variable == "Cryo", "CA", "Non-CA")
T_sub_count_sum$variable <- factor(T_sub_count_sum$variable, levels = c("Non-CA", "CA"))
T_sub_count_Plot <- ggplot(T_sub_count_sum, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.8, stat = "identity") +
    ylab("Percentage") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "grey", size = 2)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = colors_bar
    )
ggsave(
    filename = "plot/scRNA3_T_count.pdf",
    plot = T_sub_count_Plot,
    height = 6, width = 4.5
)
