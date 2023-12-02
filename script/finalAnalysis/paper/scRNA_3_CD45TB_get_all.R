source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis.rda")
load("Combined_analysis_TNK_real.rda")
load("T_sub_real_meta_220713.rda")
T_sub_real@meta.data <- T_meta_220713
Idents(T_sub_real) <- "IFN_Final"
T_sub <- subset(T_sub_real, idents = setdiff(unique(T_sub_real$IFN_Final), c("NK", "NK_Proliferative")))
T_sub$IFN_Final <- unfactor(T_sub$IFN_Final)
T_sub$IFN_Final[T_sub$IFN_Final == "Tfh_I-IFN"] <- "CD4_T_I-IFN"
Idents(T_sub) <- "IFN_Final"
rm("T_sub_real")

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

#####
# 1. all CD45+ cells:
Cryo_merge

# 1.1 CD45+ umap
Cryo_merge$MainCluster <- factor(unfactor(Cryo_merge$MainCluster), levels = c("T", "B", "NK", "Myeloid", "Mast"))
CD45plus_plot <- my_plotDim(
    Cryo_merge,
    reduction = "umap", label = T,
    pt.size = 1.5, label.size = 6, repel = T,
    group.by = "MainCluster",
    title = "CD45+ MainCluster"
) +
    theme(
        legend.position = c(0.8, 0.2),
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
            "T" = "#c0abcb",
            "B" = "#4f8ab2",
            "NK" = "#de9393",
            "Myeloid" = "#a1c4d6",
            "Mast" = "#d86762"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
ggsave(
    filename = "plot/scRNA3_CD45plus_dim.pdf",
    plot = CD45plus_plot,
    height = 6, width = 7.69
)

my_CountCluster(Cryo_merge, group2 = "orig.ident2", group1 = "active.ident")

# 1.2 CD45+ markers express
CD45_Marker_rough <- list(
    T = c("Cd3e", "Cd247"),
    NK = c("Klrk1", "Klrb1c"),
    B = c("Cd19", "Cd79a"),
    Myeloid = c("Itgam", "Cd14"),
    Mast = c("Mcpt1", "Mcpt2")
)

a <- list()
markers_list <- c("Cd3e", "Klrk1", "Cd19", "Itgam", "Mcpt1", "Cd247", "Klrb1c", "Cd79a", "Cd14", "Mcpt2")
for (i in markers_list) {
    a[[i]] <- FeaturePlot(
        Cryo_merge,
        features = i, reduction = "umap",
        pt.size = 0.5, cols = c("#e0e7c8", "#cc340a")
    ) +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            text = element_blank()
        ) +
        annotate(
            "text",
            x = 9.5,
            y = -10,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
CD45plus_express <- patchwork::wrap_plots(plots = a, ncol = 5)
ggsave(
    filename = "plot/scRNA3_CD45plus_markers_umap.pdf",
    plot = CD45plus_express,
    height = 6, width = 7.69 * 1.5
)


# 1.3 CD45+ markers violin
Cryo_merge$MainCluster <- factor(unfactor(Cryo_merge$MainCluster), levels = rev(c("T", "B", "NK", "Myeloid", "Mast")))
CD45plus_violin <- VlnPlot(
    Cryo_merge,
    feature = unlist(CD45_Marker_rough),
    group.by = "MainCluster",
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    cols = c(
        "T" = "#c0abcb",
        "B" = "#4f8ab2",
        "NK" = "#de9393",
        "Myeloid" = "#a1c4d6",
        "Mast" = "#d86762"
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
    filename = "plot/scRNA3_CD45plus_markers_violin.pdf",
    plot = CD45plus_violin,
    height = 6, width = 7.69
)

# 1.4 CD45+ celltype count
colors_bar <- c(
    "T" = "#c0abcb",
    "B" = "#4f8ab2",
    "NK" = "#de9393",
    "Myeloid" = "#a1c4d6",
    "Mast" = "#d86762"
)
CD45plus_count_sum <- as.data.frame(my_CountCluster(Cryo_merge, group1 = "MainCluster", group2 = "orig.ident2")[-6, 1:2])
CD45plus_count_sum2 <- my_CountCluster(Cryo_merge, group1 = "MainCluster", group2 = "orig.ident2")["sum", 1:2]
CD45plus_count_sum[, 1] <- CD45plus_count_sum[, 1] / CD45plus_count_sum2[1]
CD45plus_count_sum[, 2] <- CD45plus_count_sum[, 2] / CD45plus_count_sum2[2]
CD45plus_count_sum[, 1] <- CD45plus_count_sum[, 1] / sum(CD45plus_count_sum[, 1])
CD45plus_count_sum[, 2] <- CD45plus_count_sum[, 2] / sum(CD45plus_count_sum[, 2])
CD45plus_count_sum
CD45plus_count_sum$cell <- factor(rownames(CD45plus_count_sum), levels = names(colors_bar))
CD45plus_count_sum <- reshape2::melt(CD45plus_count_sum, id.vars = "cell")
CD45plus_count_sum$variable <- sub("_count$", "", CD45plus_count_sum$variable)
CD45plus_count_sum$variable <- ifelse(CD45plus_count_sum$variable == "Cryo", "CA", "Non-CA")
CD45plus_count_sum$variable <- factor(CD45plus_count_sum$variable, levels = c("Non-CA", "CA"))
CD45plus_count_Plot <- ggplot(CD45plus_count_sum, aes(x = variable, y = value, fill = cell)) +
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
        axis.ticks.y = element_line(color = "grey", size = 2, lineend = "square"),
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
    filename = "plot/scRNA3_CD45plus_count.pdf",
    plot = CD45plus_count_Plot,
    height = 6, width = 3.4
)
#####

#####
# 2. all T cells
T_sub

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
            "Th17" = "#f05c3b",
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


T_plot2 <- my_plotDim(
    T_sub,
    reduction = "umap", label = F,
    pt.size = 2, label.size = 6, repel = T,
    group.by = "orig.ident2",
    title = "T cells Case-Control"
)
T_plot2 <- T_plot2 +
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
            "NonCryo" = "#00bfc4",
            "Cryo" = "#f8766d"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
T_plot2
ggsave(
    filename = "plot/scRNA3_T_dim_casecontrol.pdf",
    plot = T_plot2,
    height = 6, width = 7.69 + 2
)

my_CountCluster(T_sub, group2 = "orig.ident2", group1 = "active.ident")

T_sub$IFN_highlight <- ifelse(T_sub$IFN_Final %in% c('CD8_Tem_I-IFN', 'CD4_T_I-IFN'), 'I-IFN-associated', 'others')
idx_tmp <- T_sub$IFN_highlight != 'others'
T_sub$IFN_highlight[idx_tmp] <- paste(T_sub$IFN_highlight[idx_tmp], T_sub$orig.ident2[idx_tmp], sep = "-")
T_sub$IFN_highlight <- factor(
    T_sub$IFN_highlight, 
    levels = c("I-IFN-associated-NonCryo","I-IFN-associated-Cryo","others")
)
levels(T_sub$IFN_highlight) <- c("I-IFN-associated Non-CA","I-IFN-associated CA","others")

df_plot3 <- as.data.frame(T_sub@reductions$umap@cell.embeddings)
df_plot3$IFN_highlight <- T_sub@meta.data[rownames(df_plot3), 'IFN_highlight']
df_plot3 <- df_plot3[order(df_plot3$IFN_highlight, decreasing = F), ]

df_plot3_IFN_CA <- df_plot3[df_plot3$IFN_highlight == 'I-IFN-associated CA',]
df_plot3_IFN_NonCA <- df_plot3[df_plot3$IFN_highlight == 'I-IFN-associated Non-CA',]

# T_plot3 <- my_plotDim(
#     T_sub,
#     # cells.highlight = colnames(T_sub)[T_sub$IFN_Final %in% c('CD8_Tem_I-IFN', 'CD4_T_I-IFN')],
#     # sizes.highlight = 2,
#     reduction = "umap", label = F,
#     pt.size = 2, label.size = 6, repel = T,
#     group.by = "IFN_highlight",
#     title = "I-IFN assiociated T cells"
# )
T_plot3 <- ggplot(aes(x = UMAP_1, y = UMAP_2, color = IFN_highlight), data = df_plot3) +
    geom_point(size = 2) + 
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = IFN_highlight), data = df_plot3_IFN_CA, size = 2) + 
    geom_point(aes(x = UMAP_1, y = UMAP_2, color = IFN_highlight), data = df_plot3_IFN_NonCA, size = 2) + 
    labs(title = 'I-IFN assiociated T cells') + 
    cowplot::theme_cowplot() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "right",
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(
        values = rev(c(
            "I-IFN-associated Non-CA" = "#00bfc4",
            "I-IFN-associated CA" = "#f8766d",
            "others" = 'grey'
        ))
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
# T_plot3
ggsave(
    filename = "plot/scRNA3_T_dim_highlight_I-IFN.pdf",
    plot = T_plot3,
    height = 6, width = 7.69 + 2
)


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
        "Th17" = "#f05c3b",
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
    "Th17" = "#f05c3b",
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
        axis.ticks.y = element_line(color = "grey", size = 2, lineend = "square"),
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
#####

#####
# 3. CD8 T cells I-IFN signiture
source("~/Desktop/Github/My_scRNA_pipeline/Pathway_score.r")
CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")
gmt_file <- "/mnt/e/Cryo-TCR/server/210924/custom_signatures_for_GSEA.gmt"
geneset <- getGmt(gmt_file)
DefaultAssay(T_sub) <- "SCT"
T_sub_CD8 <- subset(T_sub, ident = c(setdiff(CD8_idents, "CD8_Naive")))
my_CountCluster(T_sub_CD8, group2 = "orig.ident2", group1 = "active.ident")

T_sub_CD8 <- geneset_score(T_sub_CD8, geneset[1:3], method = "average", slot = "data", highly_variable = T)
T_sub_CD8[["IFN.gamma"]] <- T_sub_CD8$HALLMARK_INTERFERON_GAMMA_RESPONSE
T_sub_CD8[["Effector"]] <- T_sub_CD8$GSE9650_EFFECTOR_GENESET
T_sub_CD8[["Cytolytic"]] <- T_sub_CD8$CYTOLYTIC_SCORE

df_signiture <- T_sub_CD8[[c("IFN.gamma", "Effector", "Cytolytic", "orig.ident2", "IFN_Final")]]
df_signiture <- df_signiture[df_signiture$IFN_Final %in% c(setdiff(CD8_idents, "CD8_Naive")), ]

df <- data.frame(
    signature = unlist(lapply(
        c("IFN.gamma", "Effector", "Cytolytic"),
        function(x) rep(x, nrow(df_signiture))
    )),
    scale_value = unlist(lapply(
        c("IFN.gamma", "Effector", "Cytolytic"),
        function(x) {
            y <- df_signiture[, x, drop = T] / max(df_signiture[, x, drop = T])
            y
        }
    )),
    bc = rep(rownames(df_signiture), 3)
)
df$group <- as.character(T_sub_CD8$orig.ident2[df$bc])
df$group[df$group == "Cryo"] <- "CA"
df$group[df$group == "NonCryo"] <- "Non-CA"
df$group <- factor(df$group, levels = c("Non-CA", "CA"))

CD8_signiture_violin <- ggplot(df, aes(x = signature, y = scale_value)) +
    geom_violin(
        aes(fill = group),
        adjust = 1,
        linetype = 0
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif"
    ) +
    ggtitle("Signitures in CD8 T cells") +
    ylab("Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = c(0.85, 0.08),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin
ggsave(
    filename = "plot/scRNA3_TCD8_signiture_violin.pdf",
    plot = CD8_signiture_violin,
    height = 6, width = 4
)

df$group <- factor(df$group, levels = c("CA", "Non-CA"))
CD8_signiture_boxplot <- ggplot(df, aes(x = signature, y = scale_value)) +
    geom_boxplot(
        # aes(fill = group, color = group),
        # color = rep(rev(c(
        #     "Non-CA" = "#55a0fb",
        #     "CA" = "#ff8080"
        # )),3),
        aes(fill = group),
        size = 1,
        outlier.color = "#a1aab7",
        outlier.size = 0.5
    ) +
    geom_jitter(color = "grey", alpha = 0.5) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 8,
        angle = 90, vjust = 1.3
    ) +
    ggtitle("Signitures in CD8 T cells") +
    xlab("Signitures") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_text(size = 8, color = "#363a38"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        # strip.text.x = element_text(angle = 45, size = 12, face = 'bold'),
        legend.position = c(0.85, 0.08),
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    coord_flip()
CD8_signiture_boxplot
ggsave(
    filename = "plot/scRNA3_TCD8_signiture_boxplot.pdf",
    plot = CD8_signiture_boxplot,
    height = 6, width = 7
)

### single
df_s <- data.frame(
    signature = unlist(lapply(
        c("IFN.gamma", "Effector", "Cytolytic"),
        function(x) rep(x, nrow(df_signiture))
    )),
    scale_value = unlist(lapply(
        c("IFN.gamma", "Effector", "Cytolytic"),
        function(x) {
            df_signiture[, x, drop = T]
        }
    )),
    bc = rep(rownames(df_signiture), 3)
)
df_s$group <- as.character(T_sub_CD8$orig.ident2[df_s$bc])
df_s$group[df_s$group == "Cryo"] <- "CA"
df_s$group[df_s$group == "NonCryo"] <- "Non-CA"
df_s$group <- factor(df_s$group, levels = c("Non-CA", "CA"))
#
df_s1 <- df_s[df_s$signature == unique(df_s$signature)[1], ]
CD8_signiture_violin_s1 <- ggplot(df_s1, aes(x = signature, y = scale_value)) +
    geom_violin(
        aes(fill = group)
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s1$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s1
ggsave(
    filename = "plot/scRNA3_TCD8_signiture1_violin.pdf",
    plot = CD8_signiture_violin_s1,
    height = 3, width = 3
)
#
df_s2 <- df_s[df_s$signature == unique(df_s$signature)[2], ]
CD8_signiture_violin_s2 <- ggplot(df_s2, aes(x = signature, y = scale_value)) +
    geom_violin(
        aes(fill = group)
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s2$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s2
ggsave(
    filename = "plot/scRNA3_TCD8_signiture2_violin.pdf",
    plot = CD8_signiture_violin_s2,
    height = 3, width = 3
)
#
df_s3 <- df_s[df_s$signature == unique(df_s$signature)[3], ]
CD8_signiture_violin_s3 <- ggplot(df_s3, aes(x = signature, y = scale_value)) +
    geom_violin(
        aes(fill = group)
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s3$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s3
ggsave(
    filename = "plot/scRNA3_TCD8_signiture3_violin.pdf",
    plot = CD8_signiture_violin_s3,
    height = 3, width = 3
)
# boxplot
#
df_s1b <- df_s[df_s$signature == unique(df_s$signature)[1], ]
CD8_signiture_violin_s1b <- ggplot(df_s1b, aes(x = signature, y = scale_value)) +
    geom_boxplot(
        aes(fill = group),
        outlier.size = 0.3,
        outlier.alpha = 0.2
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s1b$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s1b
ggsave(
    filename = "plot/scRNA3_TCD8_signiture1_boxplot.pdf",
    plot = CD8_signiture_violin_s1b,
    height = 3, width = 3
)
#
df_s2b <- df_s[df_s$signature == unique(df_s$signature)[2], ]
CD8_signiture_violin_s2b <- ggplot(df_s2b, aes(x = signature, y = scale_value)) +
    geom_boxplot(
        aes(fill = group),
        outlier.size = 0.3,
        outlier.alpha = 0.2
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s2b$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s2b
ggsave(
    filename = "plot/scRNA3_TCD8_signiture2_boxplot.pdf",
    plot = CD8_signiture_violin_s2b,
    height = 3, width = 3
)
#
df_s3b <- df_s[df_s$signature == unique(df_s$signature)[3], ]
CD8_signiture_violin_s3b <- ggplot(df_s3b, aes(x = signature, y = scale_value)) +
    geom_boxplot(
        aes(fill = group),
        outlier.size = 0.3,
        outlier.alpha = 0.2
    ) +
    stat_compare_means(
        aes(group = group),
        method = "t.test",
        method.args = list(alternative = "less"),
        label = "p.signif",
        size = 4
    ) +
    ggtitle(paste(unique(df_s3b$signature), "signiture")) +
    ylab("Signiture Score") +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "black"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 8, color = "#363a38"),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
        legend.position = "right",
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent")
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    )
CD8_signiture_violin_s3b
ggsave(
    filename = "plot/scRNA3_TCD8_signiture3_boxplot.pdf",
    plot = CD8_signiture_violin_s3b,
    height = 3, width = 3
)
#####

#####
# 4. T cells GO
T_sub
CD4_idents <- c("CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh")
CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")

makeTable <- function(GO_result, width = 60) {
    merge <- function(vec, len, sep = "\n") {
        if (length(vec) < len) {
            paste(vec, collapse = "")
        } else {
            x1 <- paste(vec[1:len], collapse = "")
            x2 <- merge(vec[(len + 1):length(vec)], len)
            paste(x1, x2, sep = sep)
        }
    }

    Description <- sapply(GO_result$Description, function(x) {
        y <- strsplit(x, split = "")[[1]]
        y[1] <- toupper(y[1])
        paste(y, collapse = "")
    })
    if (sum(nchar(Description) > width) == 0) {
        Description <- sapply(Description, function(x) {
            ifelse(
                nchar(x) < width,
                paste0(paste(rep(" ", width - nchar(x)), collapse = ""), x),
                x
            )
        })
    }
    Description <- sapply(Description, function(x) {
        if (nchar(x) <= width) {
            x
        } else {
            y <- strsplit(x, "")[[1]]
            merge(y, width, "\n")
        }
    })

    Description <- factor(
        Description,
        levels = rev(Description)
    )
    pvalue <- GO_result$p.adjust
    df <- data.frame(
        name = Description,
        pvalue = -log10(pvalue)
    )
    return(df)
}

# 4.1 CA_vs_Non-CA with all T cells
Idents(T_sub) <- "orig.ident2"
Markers_T_list <- FindAllMarkers(T_sub, only.pos = T)
Markers_T_df <- my_Markers2df_multiple(
    Markers_T_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 120
)

GO_T <- my_GO(
    Markers_T_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_T_res <- GO_T@result
df_for_plot <- makeTable(GO_T_res[1:25, ], 60)
T_GO_bar <- ggplot(
    df_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with all T cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(rep(0.5,4),'cm')
    )
T_GO_bar
ggsave(
    filename = "plot/scRNA3_T_GO_bar.pdf",
    plot = T_GO_bar,
    height = 6, width = 9
)

# 4.2 CA_vs_Non-CA with CD4-T cells
Idents(T_sub) <- "IFN_Final"
T_CD4 <- subset(T_sub, ident = CD4_idents)
Idents(T_CD4) <- "orig.ident2"
Markers_TCD4_list <- FindAllMarkers(T_CD4, only.pos = T)
Markers_TCD4_df <- my_Markers2df_multiple(
    Markers_TCD4_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 120
)

GO_TCD4 <- my_GO(
    Markers_TCD4_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_TCD4_res <- GO_TCD4@result
df4_for_plot <- makeTable(GO_TCD4_res[1:25, ], 60) # filter
TCD4_GO_bar <- ggplot(
    df4_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with CD4-T cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
TCD4_GO_bar
ggsave(
    filename = "plot/scRNA3_TCD4_GO_bar.pdf",
    plot = TCD4_GO_bar,
    height = 6, width = 9
)

# 4.3 CA_vs_Non-CA with CD8-T cells
Idents(T_sub) <- "IFN_Final"
T_CD8 <- subset(T_sub, ident = CD8_idents)
Idents(T_CD8) <- "orig.ident2"
Markers_TCD8_list <- FindAllMarkers(T_CD8, only.pos = T)
Markers_TCD8_df <- my_Markers2df_multiple(
    Markers_TCD8_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 120
)

GO_TCD8 <- my_GO(
    Markers_TCD8_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_TCD8_res <- GO_TCD8@result
df8_for_plot <- makeTable(GO_TCD8_res[1:25, ], 60) # filter

TCD8_GO_bar <- ggplot(
    df8_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with CD8-T cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,4)),'cm')
    )
TCD8_GO_bar
ggsave(
    filename = "plot/scRNA3_TCD8_GO_bar.pdf",
    plot = TCD8_GO_bar,
    height = 6, width = 9
)

# 4.4 CA_vs_Non-CA with CD8-Tem cells
Idents(T_sub) <- "IFN_Final"
Tem_CD8 <- subset(T_sub, ident = c("CD8_Tem_I-IFN", "CD8_Tem"))
Idents(Tem_CD8) <- "orig.ident2"
Markers_Tem_CD8_list <- FindAllMarkers(Tem_CD8, only.pos = T)
Markers_Tem_CD8_df <- my_Markers2df_multiple(
    Markers_Tem_CD8_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 85
)

GO_Tem_CD8 <- my_GO(
    Markers_Tem_CD8_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_Tem_CD8_res <- GO_Tem_CD8@result
dfem_for_plot <- makeTable(GO_Tem_CD8_res[1:25, ], 60) # filter

Tem_CD8_GO_bar <- ggplot(
    dfem_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with CD8-Tem cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,4)),'cm')
    )
Tem_CD8_GO_bar
ggsave(
    filename = "plot/scRNA3_TemCD8_GO_bar.pdf",
    plot = Tem_CD8_GO_bar,
    height = 6, width = 9
)

# 4.5 CA_vs_Non-CA with T-others cells
Idents(T_sub) <- "IFN_Final"
T_others <- subset(T_sub, ident = "T_others")
Idents(T_others) <- "orig.ident2"
Markers_T_others_list <- FindAllMarkers(T_others, only.pos = T)
Markers_T_others_df <- my_Markers2df_multiple(
    Markers_T_others_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 100
)

GO_Tothers <- my_GO(
    Markers_T_others_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_Tothers_res <- GO_Tothers@result
dfo_for_plot <- makeTable(GO_Tothers_res[1:25, ], 60) # filter

Tothers_GO_bar <- ggplot(
    dfo_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with T-others cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,4)),'cm')
    )
Tothers_GO_bar
ggsave(
    filename = "plot/scRNA3_Tothers_GO_bar.pdf",
    plot = Tothers_GO_bar,
    height = 6, width = 9
)

# 4.6 CA_vs_Non-CA with B cells
Idents(Cryo_merge) <- "MainCluster"
B <- subset(Cryo_merge, ident = "B")
Idents(B) <- "orig.ident2"
Markers_B_list <- FindAllMarkers(B, only.pos = T)
Markers_B_df <- my_Markers2df_multiple(
    Markers_B_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 60
)

GO_B <- my_GO(
    Markers_B_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_B_res <- GO_B@result
dfB_for_plot <- makeTable(GO_B_res[1:25, ], 55) # filter

B_GO_bar <- ggplot(
    dfB_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with B cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
B_GO_bar
ggsave(
    filename = "plot/scRNA3_B_GO_bar.pdf",
    plot = B_GO_bar,
    height = 6, width = 9
)
#####
