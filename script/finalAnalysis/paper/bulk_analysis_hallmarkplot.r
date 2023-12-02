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
# HallMark pathway
# GSVA for pathway
library(GSVA)
library(GSEABase)
red <- "#ff5656" # "red"
blue <- "#5858ff" # "blue"
RWB_color <- grDevices::colorRampPalette(colors = c(blue, "white", red))(100)
GSV_path <- "/mnt/c/DeepTL_dev-main/pathway/MSigDB/"
gmt_file2 <- paste0(GSV_path, "h.all.v7.4.symbols.gmt")
geneset2 <- getGmt(gmt_file2)
exp_X_Poisson <- data@assays@data$counts
rownames(exp_X_Poisson) <- toupper(rownames(exp_X_Poisson))
## GSVA analysis by samples
GSVA_res_Poisson2 <- gsva(
    exp_X_Poisson,
    geneset2,
    min.sz = 10,
    max.sz = 1000,
    verbose = TRUE,
    method = "gsva",
    kcdf = "Poisson"
)

rownames(GSVA_res_Poisson2) <- gsub("_", " ", sub("^HALLMARK_", "", rownames(GSVA_res_Poisson2)))

GSVA_res_Poisson2

pathway2_pvalue_t_up <- apply(GSVA_res_Poisson2, 1, function(x) t.test(x[1:3], x[4:6], alternative = "greater")$p.value)
pathway2_pvalue_t_down <- apply(GSVA_res_Poisson2, 1, function(x) t.test(x[1:3], x[4:6], alternative = "less")$p.value)
sum(pathway2_pvalue_t_down < 0.05)
pathway2_pvalue_t_up <- sort(pathway2_pvalue_t_up, decreasing = F)
pathway2_up_names <- names(pathway2_pvalue_t_up)[pathway2_pvalue_t_up < 0.05]
pathway2_no_names <- names(pathway2_pvalue_t_up)[!pathway2_pvalue_t_up < 0.05]
pathway2_up_names
pathway2_no_names

pathway2_pvalue_t <- apply(GSVA_res_Poisson2, 1, function(x) t.test(x[1:3], x[4:6], alternative = "greater")$p.value)
pathway2_pvalue_t <- sort(pathway2_pvalue_t, decreasing = F)

pathway2_pvalue_t <- pathway2_pvalue_t[c(pathway2_up_names, setdiff(names(pathway2_pvalue_t), pathway2_up_names))]

GSVA_scale2 <- t(apply(GSVA_res_Poisson2, 1, scale))
colnames(GSVA_scale2) <- colnames(GSVA_res_Poisson2)
GSVA_scale2 <- GSVA_scale2[names(pathway2_pvalue_t), c(4:6, 1:3)]

GSVA_scale2_no <- GSVA_scale2[pathway2_no_names, ]
GSVA_scale2_up <- GSVA_scale2[pathway2_up_names, ]

xx <- ComplexHeatmap::Heatmap(
    GSVA_scale2_no, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 2),
    cluster_columns = F,
    show_row_dend = F,
    show_column_names = F
)
GSVA_scale2_all <- rbind(
    GSVA_scale2_no[rownames(xx@matrix), ],
    GSVA_scale2_up
)
GSVA_scale2_all_number <- unlist(lapply(
    rownames(GSVA_scale2_all),
    function(X) {
        if (X %in% rownames(GSVA_scale2_no)) {
            0
        } else if (X %in% rownames(GSVA_scale2_up)) {
            1
        }
    }
))
pathway2_pvalue_t <- pathway2_pvalue_t[rownames(GSVA_scale2_all)]
pvalue_vec2 <- unlist(lapply(
    pathway2_pvalue_t,
    function(x) {
        if (x > 0.05) {
            "Not significant"
        } else if (x > 0.01 & x <= 0.05) {
            "P-value < 0.05"
        } else if (x > 0.001 & x <= 0.01) {
            "P-value < 0.01"
        } else if (x <= 0.001) "P-value < 0.001"
    }
))
row_split <- data.frame(
    group = factor(GSVA_scale2_all_number, labels = c("No", "Up")),
    row.names = rownames(GSVA_scale2_all)
)
column_split <- data.frame(
    group = factor(c(rep("Non-CA", 3), rep("CA", 3)), levels = c("Non-CA", "CA")),
    row.names = colnames(GSVA_scale2_all)
)

rownames(GSVA_scale2_all) <- tolower(rownames(GSVA_scale2_all))

annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(GSVA_scale2_all)
)
pp_pathway_custom2 <- ComplexHeatmap::Heatmap(
    GSVA_scale2_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = T, 
    show_column_dend = F,
    cluster_rows = T,
    # show_heatmap_legend = F,
    show_row_dend = T,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    column_split = column_split,
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    row_names_gp = gpar(fontsize = 12),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange")),
        show_legend = F
    ),
    left_annotation = rowAnnotation(
        Significant = factor(
            pvalue_vec2,
            levels = c(
                "Not significant",
                "P-value < 0.05",
                "P-value < 0.01",
                "P-value < 0.001",
                "P-value < 0.0001"
            )
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
pp_pathway_custom2


library(showtext)
font_add("Arial", "/usr/share/fonts/windows10/arial.ttf")
font_add("Bold", "/usr/share/fonts/windows10/arialbd.ttf")
showtext_auto()
setFont <- function(p, font_family = "Arial", legend_font_family = "Bold") {
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

    p
}

pdf("~/Desktop/Cryo/Paper_Result_Plot/plot/pathway_hallmark_heatmap.pdf", width = 7, height = 10)
draw(setFont(pp_pathway_custom2))
dev.off()
