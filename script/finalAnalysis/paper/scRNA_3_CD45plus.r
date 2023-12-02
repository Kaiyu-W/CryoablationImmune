source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis.rda")

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

# 1. all CD45+ cells:
Cryo_merge

if (!file.exists("scRNA3_cellbarcode.csv")) {
    write.csv(colnames(Cryo_merge), file = "scRNA3_cellbarcode.csv")
}
if (!file.exists("scRNA3_umap.csv")) {
    write.csv(Cryo_merge@reductions$umap@cell.embeddings, file = "scRNA3_umap.csv")
}
if (!file.exists("scRNA3_meta.csv")) {
    write.csv(Cryo_merge@meta.data, file = "scRNA3_meta.csv")
}
# save gene expression
if (!dir.exists("plot/scRNA3_CD45plus_matrix") & !file.exists("plot/scRNA3_CD45plus_matrix.zip")) {
    dir.create("plot/scRNA3_CD45plus_matrix")
    # mtx
    Matrix::writeMM(Cryo_merge@assays$RNA@counts, "plot/scRNA3_CD45plus_matrix/matrix.mtx")
    system("gzip plot/scRNA3_CD45plus_matrix/matrix.mtx")
    # scales::number_bytes((file.size("plot/scRNA3_CD45plus_matrix/matrix.mtx.gz")))
    # barcode
    write.table(
        colnames(Cryo_merge@assays$RNA@counts), 
        file = "plot/scRNA3_CD45plus_matrix/barcodes.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F
    )
    system("gzip -f plot/scRNA3_CD45plus_matrix/barcodes.tsv")
    # feature
    feature_list = rownames(Cryo_merge@assays$RNA@counts)
    feature_df = read.table("/mnt/e/Cryo-TCR/server/210924/Cryo/filtered_feature_bc_matrix/features.tsv.gz", sep = '\t')
    feature_df$V2 <- make.unique(feature_df$V2)
    feature_df <- feature_df[match(feature_list, feature_df$V2), ]
    write.table(
        feature_df, 
        file = "plot/scRNA3_CD45plus_matrix/features.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F,
        sep = '\t'
    )
    system("gzip -f plot/scRNA3_CD45plus_matrix/features.tsv")
}

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
