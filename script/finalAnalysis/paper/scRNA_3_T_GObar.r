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
