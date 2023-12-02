setwd("/mnt/e/Cryo-TCR/bulk/res_counts")
# source("../../server/auto/utilities.r")
library(ComplexHeatmap)

if (file.exists("~/Desktop/Cryo/Paper_Result_Plot/bulk_deseq2_data.rds")) {
    data <- readRDS("~/Desktop/Cryo/Paper_Result_Plot/bulk_deseq2_data.rds")
} else {
    library(DESeq2)
    count <- read.table("genome_counts.tsv", sep = "\t", header = T, row.names = 1)
    idx <- read.table("transcriptome_counts.tsv", sep = "\t", header = T, row.names = 1)

    # remove 3
    count$`M.Cryo.3` <- NULL
    colnames(count)[3] <- "M.Cryo.3"

    data <- DESeqDataSetFromMatrix(
        countData = count,
        colData = data.frame(
            row.names = colnames(count),
            condition = factor(c(rep("CA", 3), rep("Non-CA", 3)), levels = c("CA", "Non-CA"))
        ),
        design = ~condition
    )
    #
    # data <- DESeqDataSetFromMatrix(
    #     countData = idx[, -1],
    #     colData = data.frame(
    #         row.names = colnames(idx[, -1]),
    #         condition = factor(c(rep('Cryo', 4), rep('NonCryo', 3)), levels = c('Cryo', 'NonCryo'))
    #     ),
    #     design = ~ condition)

    data <- data[rowSums(counts(data)) > 1, ]
    rl_data <- rlog(object = data, blind = F)

    data <- estimateSizeFactors(data)
    # normalized_counts <- counts(data, normalized = T)

    data <- DESeq(data)

    saveRDS(data, "~/Desktop/Cryo/Paper_Result_Plot/bulk_deseq2_data.rds")
}

setFont <- function(p, font_family = "Arial", legend_font_family = "Bold") {
    library(showtext)
    font_add("Arial", "/usr/share/fonts/windows10/arial.ttf")
    font_add("Bold", "/usr/share/fonts/windows10/arialbd.ttf")
    showtext_auto()

    p@column_names_param$gp$fontfamily <- font_family # 列名称, no use("M.Nonc.1" "M.Nonc.2" "M.Nonc.3" "M.Cryo.1" "M.Cryo.2" "M.Cryo.3")
    p@column_names_param$gp$fontsize <- 12

    p@row_names_param$gp$fontfamily <- font_family # 行名称(gene names)
    p@row_names_param$gp$fontsize <- 12

    p@row_title_param$gp$fontfamily <- font_family # 行标题(Up/Down)
    p@row_title_param$gp$fontsize <- 13

    p@column_title_param$gp$fontfamily <- font_family # 列标题(CA/Non-CA)
    p@column_title_param$gp$fontsize <- 13

    p@matrix_legend_param <- gpar( # 矩阵legend(Z-score)
        title_gp = gpar(fontfamily = legend_font_family, fontsize = 10),
        labels_gp = gpar(fontfamily = font_family, fontsize = 10)
    )

    p@left_annotation@anno_list[["Significant"]]@legend_param <- gpar( # 行分类legend (Significant)
        title_gp = gpar(fontfamily = legend_font_family, fontsize = 10),
        labels_gp = gpar(fontfamily = font_family, fontsize = 10)
    )

    p@top_annotation@anno_list[["Group"]]@name_param$gp$fontfamily <- font_family # 列分类标题 (Group)
    p@top_annotation@anno_list[["Group"]]@name_param$gp$fontsize <- 12

    return(p)
}

#############################
############################# GSVA for pathway
library(GSVA)
library(GSEABase)
red <- "#ff5656" # "red"
blue <- "#5858ff" # "blue"
RWB_color <- grDevices::colorRampPalette(colors = c(blue, "white", red))(100)
GSV_path <- "/mnt/c/DeepTL_dev-main/pathway/GSVA/"
gmt_file <- paste0(GSV_path, "c5.go.v7.4.symbols1.gmt")
geneset <- getGmt(gmt_file)

exp_X_Poisson <- data@assays@data$counts
rownames(exp_X_Poisson) <- toupper(rownames(exp_X_Poisson))
## GSVA analysis by samples
GSVA_res_Poisson <- gsva(
    exp_X_Poisson,
    geneset,
    min.sz = 10,
    max.sz = 1000,
    verbose = TRUE,
    method = "gsva",
    kcdf = "Poisson"
)

GSVA_diff <- apply(GSVA_res_Poisson, 1, function(x) wilcox.test(x[1:3], x[4:6], alternative = "greater")$p.value)
# GSVA_diff_padj <- p.adjust(GSVA_diff, method = 'BH')
GSVA_diff <- sort(GSVA_diff, decreasing = F)
# GSVA_diff_padj <- sort(GSVA_diff_padj, decreasing = F)
pathway_up <- names(GSVA_diff)[GSVA_diff <= 0.05]
GSVA_scale <- t(apply(GSVA_res_Poisson, 1, scale))
colnames(GSVA_scale) <- colnames(GSVA_res_Poisson)

# custom
pathway_custom0 <- c(
    "chemokine_production",
    "response_to_chemokine",
    "cytokine_production",
    "response_to_interferon_gamma",
    "response_to_interferon_alpha",
    "response_to_interferon_beta",
    "interferon_alpha_production",
    "interferon_beta_production",
    "interferon_gamma_production",
    "inflammatory_response"
)
pathway_custom <- sapply(pathway_custom0, function(x) names(GSVA_diff)[grep(paste0("^GO.._", x, "$"), names(GSVA_diff), ignore.case = T)], USE.NAMES = F)
pathway_pvalue_t <- apply(GSVA_res_Poisson[pathway_custom, ], 1, function(x) t.test(x[1:3], x[4:6], alternative = "greater")$p.value)
pathway_up_custom <- unlist(sapply(pathway_custom, function(x) pathway_up[grep(paste0("^", x, "$"), pathway_up, ignore.case = T)], USE.NAMES = F))

custom_scale_mtx <- GSVA_scale[pathway_custom, c(4:6, 1:3)]
rownames(custom_scale_mtx) <- gsub("_", " ", tolower(sub("^GO.._", "", rownames(custom_scale_mtx))))

annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(custom_scale_mtx)
)

column_split <- data.frame(
    group = factor(
        c(rep("Non-CA", 3), rep("CA", 3)),
        levels = c("Non-CA", "CA")
    ),
    row.names = colnames(custom_scale_mtx)
)
pvalue_vec <- unlist(lapply(
    pathway_pvalue_t,
    function(x) {
        if (x > 0.051) {
            "Not significant"
        } else if (x > 0.01 & x <= 0.051) {
            "P-value < 0.05"
        } else if (x > 0.001 & x <= 0.01) {
            "P-value < 0.01"
        } else if (x <= 0.001) "P-value < 0.001"
    }
))
pp_pathway_custom <- ComplexHeatmap::Heatmap(
    custom_scale_mtx, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = T, 
    show_column_dend = F,
    cluster_rows = T,
    show_row_dend = T,
    show_column_names = F,
    show_row_names = T,
    column_split = column_split,
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    row_names_gp = gpar(fontsize = 12),
    top_annotation = HeatmapAnnotation(
        Group = factor(
            annotation_col$Group,
            levels = c("Non-CA", "CA")
        ),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange")),
        show_legend = F
    ),
    left_annotation = rowAnnotation(
        Significant = factor(
            pvalue_vec,
            levels = c("Not significant", "P-value < 0.05", "P-value < 0.01", "P-value < 0.001", "P-value < 0.0001")
        ),
        col = list(
            Significant = c(
                "Not significant" = "#D3D3D3",
                "P-value < 0.05" = "#B4EEB4",
                "P-value < 0.01" = "#3CB371",
                "P-value < 0.001" = "#006400"
            )
        ),
        gp = gpar(col = "black", lwd = 1),
        show_annotation_name = F
    )
)
pp_pathway_custom

pdf("~/Desktop/Cryo/Paper_Result_Plot/plot/pathway_custom_heatmap.pdf", width = 6.65, height = 3)
draw(setFont(pp_pathway_custom))
dev.off()
