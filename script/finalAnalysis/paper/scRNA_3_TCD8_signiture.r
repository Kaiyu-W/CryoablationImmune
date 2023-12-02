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

CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")

# 3. CD8 T cells I-IFN signiture
library(ggpubr)
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

# violin + boxplot
plot_violinBox <- function(
    df, 
    jitter_alpha = 0.1, jitter_size = 0.1, jitter_width = 0.2, boxplot_width = 0.2,
    outlier.size = 0, outlier.shape = NA, outlier.alpha = 0.1
) {
    p <- ggviolin(
        data = df,
        x = "group",
        y = "scale_value",
        fill = 'group',
        palette = c("#55a0fb","#ff8080"),
        width = 0.9,
        trim = TRUE,
        draw_quantiles = NULL,
        error.plot = 'crossbar',
        add = c("boxplot", "jitter"), 
        add.params  = list(
            color = 'black',
            fill = 'white',
            alpha = jitter_alpha,
            size = jitter_size,
            jitter = jitter_width,
            width = boxplot_width
        ),
        scale = 'width'
    ) +
        # scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
        ggtitle(paste(unique(df$signature), "signiture")) +
        ylab("Signiture Score") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12, color = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black"),
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(
                color = "black", size = 12
            ),
            axis.line.x.bottom = element_line(),
            axis.line.y.left = element_line(),
            strip.text.x = element_text(angle = 45, size = 12, face = "bold"),
            legend.position = "none",
            # legend.text = element_text(color = "black", size = 12),
            # legend.title = element_blank(),
            # legend.background = element_rect(fill = "transparent"),
            panel.background = element_rect(fill = "transparent")
        ) + stat_compare_means(
            aes(group = group),
            method = "t.test",
            method.args = list(alternative = "less"),
            label = "p.signif",
            size = 4,
            comparisons = list(1:2)
        )
        # theme(
        #     legend.position = "none",
        #     axis.text.y = element_text(color = "black", size = 12),
        #     axis.text.x = element_text(
        #         color = "black", size = 12,
        #         angle = angle, vjust = x_vjust, hjust = x_hjust
        #     ),
        #     axis.title.y = element_text(size = 12),
        #     axis.title.x = element_blank(),
        #     axis.ticks.x = element_blank(),
        #     axis.ticks.y = element_line(color = "black")
        # )
    for (i in seq(length(p$layers))) {
        if (class(p$layers[[i]]$geom)[1] == 'GeomBoxplot') {
            p$layers[[i]]$geom_params$outlier.size <- outlier.size
            p$layers[[i]]$geom_params$outlier.shape <- outlier.shape
            p$layers[[i]]$geom_params$outlier.alpha <- outlier.alpha
            p$layers[[2]]$aes_params$alpha <- 1
        }
    }
    p
}

CD8_signiture_viobox_s1 <- plot_violinBox(df_s1)
CD8_signiture_viobox_s1
ggsave(
    filename = "plot/scRNA3_TCD8_signiture1_violinbox.pdf",
    plot = CD8_signiture_viobox_s1,
    height = 4, width = 3.5
)

CD8_signiture_viobox_s2 <- plot_violinBox(df_s2)
CD8_signiture_viobox_s2
ggsave(
    filename = "plot/scRNA3_TCD8_signiture2_violinbox.pdf",
    plot = CD8_signiture_viobox_s2,
    height = 4, width = 3.5
)

CD8_signiture_viobox_s3 <- plot_violinBox(df_s3)
CD8_signiture_viobox_s3
ggsave(
    filename = "plot/scRNA3_TCD8_signiture3_violinbox.pdf",
    plot = CD8_signiture_viobox_s3,
    height = 4, width = 3.5
)
# 


# # train the linear model for classfication
my_linearModel_SignitureCluster <- function(obj, gene_list, labels, group.by, slot = "data", assay = "RNA") {
    GEX <- methods::slot(obj@assays[[assay]], slot)
    gene_list <- intersect(gene_list, rownames(GEX))
    data <- as.data.frame(t(GEX[gene_list, ]))
    data$label <- as.factor(ifelse(
        obj@meta.data[colnames(GEX), group.by] %in% labels,
        1,
        0
    ))
    my_formula <- paste("label ~", paste(gene_list, collapse = " + "))
    fit <- eval(parse(text = paste0("glm(", my_formula, ", data = data, family = binomial())")))
    fit
}
my_linearModel_SignitureScore.default <- function(obj, fit, use_predict = FALSE, slot = "data", assay = "RNA") {
    GEX <- methods::slot(obj@assays[[assay]], slot)
    data <- as.data.frame(t(GEX))
    res <- if (use_predict) {
        predict(fit, newdata = data, type = "response")
    } else {
        coeff <- fit$coefficients
        coeff <- coeff[names(coeff) != "(Intercept)"]
        apply(
            data[, names(coeff), drop = F],
            1,
            function(x) sum(x * coeff)
        )
    }
    return(res)
}
my_linearModel_SignitureScore <- function(...) {
    res <- my_linearModel_SignitureScore.default(...)
    gc()
    return(res)
}
# my_plotScore <- function(
#     obj, score_col, group.by,
#     idents = NULL, label = NULL,
#     fill_col = "red", jitter = TRUE,
#     outlier.shape = NA, outlier.size = 1, outlier.alpha = 0.1,
#     angle = 45, x_vjust = 1, x_hjust = 1,
#     use_ggpubr = F, comparisons_list = NULL,
#     use_kw = F, kw_x = 1, kw_y = 0.9
# ) {
#     df <- obj@meta.data[, c(score_col, group.by)]
#     rownames(df) <- 1:nrow(df)

#     if (!is.null(idents)) {
#         df <- df[df[[group.by]] %in% idents, , drop = F]
#         df[[group.by]] <- factor(df[[group.by]], levels = idents)
#     }
#     df$tmp1 <- df[[group.by]]
#     df$tmp2 <- df[[score_col]]

#     if (use_ggpubr) {
#         if (!require(ggpubr)) stop("No ggpubr!")
#         p <- ggboxplot(
#             data = df,
#             x = "tmp1",
#             y = "tmp2",
#             fill = fill_col,
#             outlier.shape = outlier.shape,
#             outlier.size = outlier.size,
#             outlier.alpha = outlier.alpha
#         ) +
#             scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
#             ylab(sub("^ ", "", paste(label, "Score"))) +
#             theme_classic() +
#             theme(
#                 legend.position = "none",
#                 axis.text.y = element_text(color = "black", size = 12),
#                 axis.text.x = element_text(
#                     color = "black", size = 12,
#                     angle = angle, vjust = x_vjust, hjust = x_hjust
#                 ),
#                 axis.title.y = element_text(size = 12),
#                 axis.title.x = element_blank(),
#                 axis.ticks.x = element_blank(),
#                 axis.ticks.y = element_line(color = "black")
#             )
#         if (use_kw) {
#             p <- p + stat_compare_means(
#                 label.y = kw_y, label.x = kw_x,
#                 method = "kruskal.test"
#             )
#         }
#         if (!is.null(comparisons_list)) {
#             p <- p + stat_compare_means(
#                 label = "p.signif",
#                 comparisons = comparisons_list
#             )
#         }
#     } else {
#         p <- ggplot(df, aes(x = tmp1, y = tmp2)) +
#             geom_boxplot(
#                 outlier.shape = outlier.shape,
#                 fill = fill_col,
#                 outlier.size = outlier.size,
#                 outlier.alpha = outlier.alpha
#             ) +
#             scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
#             ylab(sub("^ ", "", paste(label, "Score"))) +
#             theme_classic() +
#             theme(
#                 legend.position = "none",
#                 axis.text.y = element_text(color = "black", size = 12),
#                 axis.text.x = element_text(
#                     color = "black", size = 12,
#                     angle = angle, vjust = x_vjust, hjust = x_hjust
#                 ),
#                 axis.title.y = element_text(size = 12),
#                 axis.title.x = element_blank(),
#                 axis.ticks.x = element_blank(),
#                 axis.ticks.y = element_line(color = "black")
#             )
#     }
#     if (jitter) p + geom_jitter(size = 0.01, alpha = 0.1) else p
# }
my_plotScore_violin <- function(
    obj, score_col, group.by,
    idents = NULL, label = NULL,
    fill_col = "red",
    palette = NULL,
    draw_quantiles = FALSE, 
    draw_boxplot = FALSE,
    boxplot_width = 0.2,
    jitter_alpha = 0.1,
    jitter_width = 0.2,
    jitter_size = 0.1,
    outlier.shape = NA, outlier.size = 1, outlier.alpha = 0.1,
    angle = 45, x_vjust = 1, x_hjust = 1,
    comparisons_list = NULL,
    use_kw = F, kw_x = 1, kw_y = 0.9
) {
    df <- obj@meta.data[, c(score_col, group.by)]
    rownames(df) <- 1:nrow(df)

    if (!is.null(idents)) {
        df <- df[df[[group.by]] %in% idents, , drop = F]
        df[[group.by]] <- factor(df[[group.by]], levels = idents)
    }
    df$tmp1 <- df[[group.by]]
    df$tmp2 <- df[[score_col]]

    if (!require(ggpubr)) stop("No ggpubr!")
    p <- ggviolin(
        data = df,
        x = "tmp1",
        y = "tmp2",
        fill = fill_col,
        palette = palette,
        width = 0.9,
        trim = TRUE,
        draw_quantiles = if (draw_quantiles) c(0.25,0.5,0.75) else NULL,
        error.plot = 'crossbar',
        add = if (draw_boxplot) c("boxplot", "jitter") else "none", 
        add.params  = list(
            color = 'black',
            fill = 'white',
            alpha = jitter_alpha,
            size = jitter_size,
            jitter = jitter_width,
            width = boxplot_width
        ),
        scale = 'width'
    ) +
        scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
        ylab(sub("^ ", "", paste(label, "Score"))) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(
                color = "black", size = 12,
                angle = angle, vjust = x_vjust, hjust = x_hjust
            ),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black")
        )
    if (use_kw) {
        p <- p + stat_compare_means(
            label.y = kw_y, label.x = kw_x,
            method = "kruskal.test"
        )
    }
    if (!is.null(comparisons_list)) {
        p <- p + stat_compare_means(
            label = "p.signif",
            comparisons = comparisons_list
        )
    }

    for (i in seq(length(p$layers))) {
        if (class(p$layers[[i]]$geom)[1] == 'GeomBoxplot') {
            p$layers[[i]]$geom_params$outlier.size <- outlier.size
            p$layers[[i]]$geom_params$outlier.shape <- outlier.shape
            p$layers[[i]]$geom_params$outlier.alpha <- outlier.alpha
            p$layers[[2]]$aes_params$alpha <- 1
        }
    }
    p
}
paper_gene_list <- list(
    cytotoxic = c('Nkg7','Ccl5','Gzma','Gzmk','Ccl4','Cst7','Itm2c','Ifng'),
    exhaustion = c('Lag3','Gzmb','Havcr2','Ptms','Cxcl13','Vcam1','Prf1','Tnfrsf9','Tigit','Pdcd1'),
    naive = c('Ccr7','Tcf7','Lef1','Actn1','S1pr1','Mal','Il7r','Plac8','Spint2','Bach2')
)
T_sub_CD8plus <- subset(T_sub, ident = CD8_idents)
T_sub_CD8$orig.ident <- factor(
    T_sub_CD8$orig.ident,
    levels = c("New_NonCryo", "New_Cryo", "Old_NonCryo", "Old_Cryo")
)
levels(T_sub_CD8$orig.ident) <- c("New\nNonCryo", "New\nCryo", "Old\nNonCryo", "Old\nCryo")
T_sub_CD8$orig.ident2 <- factor(
    T_sub_CD8$orig.ident2,
    levels = c("NonCryo", "Cryo")
)
levels(T_sub_CD8$orig.ident2) <- c('Non-CA', 'CA')
fit1p <- my_linearModel_SignitureCluster(
    T_sub_CD8plus, paper_gene_list$naive,
    label = c("CD8_Naive"),
    group.by = "IFN_Final", assay = "SCT"
)
T_sub_CD8plus$Na_score_p <- my_linearModel_SignitureScore(T_sub_CD8plus, fit1p, use_predict = T, assay = "SCT")
fit2p <- my_linearModel_SignitureCluster(
    T_sub_CD8plus, paper_gene_list$exhaustion,
    label = c("CD8_Tex", "CD8_Tex_Proliferative"),
    group.by = "IFN_Final", assay = "SCT"
)
T_sub_CD8$Ex_score_p <- my_linearModel_SignitureScore(T_sub_CD8, fit2p, use_predict = T, assay = "SCT")
fit3p <- my_linearModel_SignitureCluster(
    T_sub_CD8plus, paper_gene_list$cytotoxic,
    label = c("CD8_Tem_I-IFN", "CD8_Tem"),
    group.by = "IFN_Final", assay = "SCT"
)
T_sub_CD8$Cy_score_p <- my_linearModel_SignitureScore(T_sub_CD8, fit3p, use_predict = T, assay = "SCT")
p_Exp_g1 <- my_plotScore_violin(
    T_sub_CD8[, T_sub_CD8$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score_p", "orig.ident",
    label = "Exhaustion",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080','#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)
p_Exp_g2 <- my_plotScore_violin(
    T_sub_CD8[, T_sub_CD8$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score_p", "orig.ident2",
    label = "Exhaustion",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)
p_Cyp_g1 <- my_plotScore_violin(
    T_sub_CD8,
    "Cy_score_p", "orig.ident",
    label = "Cytotoxic",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080','#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)
p_Cyp_g2 <- my_plotScore_violin(
    T_sub_CD8,
    "Cy_score_p", "orig.ident2",
    label = "Cytotoxic",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)
p_Nap_g1 <- my_plotScore_violin(
    T_sub_CD8plus[, T_sub_CD8plus$IFN_Final %in% c(
        "CD8_Naive"
    )],
    "Na_score_p", "orig.ident",
    label = "Naive",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080','#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)

p_Nap_g2 <- my_plotScore_violin(
    T_sub_CD8plus[, T_sub_CD8plus$IFN_Final %in% c(
        "CD8_Naive"
    )],
    "Na_score_p", "orig.ident2",
    label = "Naive",
    fill_col = "tmp1",
    palette = c('#55a0fb','#ff8080'),
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    draw_boxplot = TRUE, boxplot_width = 0.15, jitter_width = 0.15,
    comparisons_list = list(1:2,3:4),
    use_kw = F
)
ggsave(
    filename = "plot/scRNA3_Score_Cytotoxic_CD8+_byPaper_violinBoxplot1.pdf",
    plot = p_Cyp_g1,
    height = 4, width = 5
)
ggsave(
    filename = "plot/scRNA3_Score_Cytotoxic_CD8+_byPaper_violinBoxplot2.pdf",
    plot = p_Cyp_g2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA3_Score_Naive_CD8+_byPaper_violinBoxplot1.pdf",
    plot = p_Nap_g1,
    height = 4, width = 5
)
ggsave(
    filename = "plot/scRNA3_Score_Naive_CD8+_byPaper_violinBoxplot2.pdf",
    plot = p_Nap_g2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA3_Score_Exhaustion_CD8+_byPaper_violinBoxplot1.pdf",
    plot = p_Exp_g1,
    height = 4, width = 5
)
ggsave(
    filename = "plot/scRNA3_Score_Exhaustion_CD8+_byPaper_violinBoxplot2.pdf",
    plot = p_Exp_g2,
    height = 4, width = 3.5
)
