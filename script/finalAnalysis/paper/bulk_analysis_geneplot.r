setwd("/mnt/e/Cryo-TCR/bulk/res_counts")
# source("../../server/auto/utilities.r")
library(ComplexHeatmap)

if (file.exists("~/Desktop/Cryo/Paper_Result_Plot/bulk_deseq2_data.rds")) {
    data <- readRDS("~/Desktop/Cryo/Paper_Result_Plot/bulk_deseq2_data.rds")
    count <- read.table("genome_counts.tsv", sep = "\t", header = T, row.names = 1)

    # resultsNames(data)
    res <- DESeq2::results(data, contrast = c("condition", "CA", "Non-CA"), alpha = 0.05)
    summary(res)

    res <- res[order(res$padj), ]
    # diff_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange > 1)
    diff_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
    diff_gene_deseq2_name <- rownames(diff_gene_deseq2)
    up_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange > 0))
    down_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange < 0))
    nosignificant_gene <- setdiff(rownames(res), c(up_genes, down_genes))
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
# search the chemokine and cytokine pathway genes to plot heatmap
setwd("/mnt/c/DeepTL_dev-main/pathway/GO_mus/direct-experimentevidence/")
go_bp <- read.table("GO_BP.tsv", sep = "\t", header = T)
go_cc <- read.table("GO_CC.tsv", sep = "\t", header = T)
go_mf <- read.table("GO_MF.tsv", sep = "\t", header = T)

go_bp <- read.table("GO_BP.tsv", sep = "\t", header = T)
go_cc <- read.table("GO_CC.tsv", sep = "\t", header = T)
go_mf <- read.table("GO_MF.tsv", sep = "\t", header = T)

searchGenefromPathway <- function(pathaway) {
    res_bp <- go_bp$element[grep(pathaway, go_bp$pathway, ignore.case = T)]
    res_cc <- go_cc$element[grep(pathaway, go_cc$pathway, ignore.case = T)]
    res_mf <- go_mf$element[grep(pathaway, go_mf$pathway, ignore.case = T)]
    res <- unique(c(res_bp, res_cc, res_mf))
    return(res)
}
searchGenefromPathway2 <- function(pathaway1, pathaway2) {
    res_bp <- go_bp$element[intersect(grep(pathaway1, go_bp$pathway, ignore.case = T), grep(pathaway2, go_bp$pathway, ignore.case = T))]
    res_cc <- go_cc$element[intersect(grep(pathaway1, go_cc$pathway, ignore.case = T), grep(pathaway2, go_cc$pathway, ignore.case = T))]
    res_mf <- go_mf$element[intersect(grep(pathaway1, go_mf$pathway, ignore.case = T), grep(pathaway2, go_mf$pathway, ignore.case = T))]
    res <- unique(c(res_bp, res_cc, res_mf))
    return(res)
}
searchPathwayfromGene <- function(Gene) {
    res_bp <- go_bp$pathway[grep(Gene, go_bp$element, ignore.case = T)]
    res_cc <- go_cc$pathway[grep(Gene, go_cc$element, ignore.case = T)]
    res_mf <- go_mf$pathway[grep(Gene, go_mf$element, ignore.case = T)]
    res <- unique(c(res_bp, res_cc, res_mf))
    return(res)
}
subsetMtx <- function(geneset) {
    geneset <- unique(geneset)
    mtx0 <- log1p(count_mtx[geneset[geneset %in% rownames(count_mtx)], ])
    mtx <- t(apply(mtx0, 1, scale))
    colnames(mtx) <- colnames(mtx0)
    mtx
}
chemokine <- searchGenefromPathway("chemokine")
cytokine <- searchGenefromPathway("cytokine")
interferon_beta <- searchGenefromPathway2("interferon", "beta")
interferon_alpha <- searchGenefromPathway2("interferon", "alpha")
interferon_I <- searchGenefromPathway2("interferon", "type I")
interferon <- unique(c(interferon_beta, interferon_alpha, interferon_I))
custom_genes <- c("Cxcr4", "Ccl19", "Cxcl9", "Cxcl3", "Cxcl10", "Cxcl11", "Cxcl12", "Stat1", "Stat2", "Ifnl2", "Ifnb1", "Ifna1", "Ifnar1", "Ifnar2", "Ifna2", "Ifna4", "Ifng", "Xcl1", "Ccl5", "Ccl9", "Ccr5", "Ifngr1")
count_mtx <- data@assays@data$counts

# pathway_Ifnar1 <- searchPathwayfromGene('Ifnar1')

#######
geneall <- rownames(count)

CC <- geneall[grep("^Cc[lr]", geneall)]
CXC <- geneall[grep("^Cxc[lr]", geneall)]
CX3C <- geneall[grep("^Cx3c[lr]", geneall)]
XC <- geneall[grep("^Xc[lr]", geneall)]
CX <- geneall[grep("^Cx[lr]", geneall)]
ACKR <- geneall[grep("^Ackr", geneall)]
others_chemokines <- geneall[grep("^Gpr35$", geneall)]
chemokines <- unique(c(CC, CXC, CX3C, XC, CX, ACKR, others_chemokines))
#######
#######
if (file.exists("~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines.csv") &
    file.exists("~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines_chemokines.csv")
) {
    cytokines <- read.csv("~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines.csv")$x
    cytokines_chemokines <- read.csv("~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines_chemokines.csv")$x
} else {
    cytokines_human <- scan("/mnt/e/Cryo-TCR/bulk/cytokines_human.txt", what = "c")
    cytokines_receptor_human <- scan("/mnt/e/Cryo-TCR/bulk/cytokines_receptor_human.txt", what = "c")
    cytokines_receptor_mouse <- scan("/mnt/e/Cryo-TCR/bulk/cytokines_receptor_mouse.txt", what = "c")

    library(biomaRt)
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://may2021.archive.ensembl.org/")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://may2021.archive.ensembl.org/")

    cytokines_mouse <- getLDS(
        attributes = "hgnc_symbol",
        filters = "hgnc_symbol",
        values = cytokines_human,
        mart = human,
        attributesL = "mgi_symbol",
        martL = mouse,
        uniqueRows = TRUE
    )
    cytokines_receptor_mouse2 <- getLDS(
        attributes = "hgnc_symbol",
        filters = "hgnc_symbol",
        values = cytokines_receptor_human,
        mart = human,
        attributesL = "mgi_symbol",
        martL = mouse,
        uniqueRows = TRUE
    )

    cytokines <- unique(c(cytokines_mouse$MGI.symbol, cytokines_receptor_mouse2$MGI.symbol, cytokines_receptor_mouse))
    cytokines_chemokines <- unique(c(chemokines, cytokines))
}

#######

chemokine_mtx <- subsetMtx(chemokine)
cytokine_mtx <- subsetMtx(cytokine)
interferon_mtx <- subsetMtx(interferon)
custom_mtx <- subsetMtx(custom_genes)
chemokines_mtx <- subsetMtx(chemokines)
cytokines_mtx <- subsetMtx(cytokines)
cytokines_chemokines_mtx <- subsetMtx(cytokines_chemokines)

diff_mtx <- subsetMtx(diff_gene_deseq2_name)

diff_targets <- intersect(up_genes, c(chemokine, cytokine))
difftarget_mtx <- subsetMtx(diff_targets)

diff_targets_ifn <- intersect(up_genes, interferon)
difftarget_ifn_mtx <- subsetMtx(diff_targets_ifn)

diff_targets_custom <- intersect(up_genes, custom_genes)
difftargets_custom_mtx <- subsetMtx(diff_targets_custom)

red <- "#ff5656" # "red"
blue <- "#5858ff" # "blue"
RWB_color <- grDevices::colorRampPalette(colors = c(blue, "white", red))(100)

# cytokines_chemokines
annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(cytokines_chemokines_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(cytokines_chemokines_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
cytokines_chemokines_mtx_up <- cytokines_chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
cytokines_chemokines_mtx_down <- cytokines_chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
cytokines_chemokines_mtx_no <- cytokines_chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]


xx <- ComplexHeatmap::Heatmap(
    cytokines_chemokines_mtx_no[, c(4:6, 1:3)], # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 2),
    cluster_columns = F,
    show_row_dend = F,
    show_column_names = F,
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
cytokines_chemokines_mtx_all <- rbind(
    #     cytokines_chemokines_mtx_no[rownames(xx@matrix), c(4:6, 1:3)],
    cytokines_chemokines_mtx_up[, c(4:6, 1:3)],
    cytokines_chemokines_mtx_down[, c(4:6, 1:3)]
)
cytokines_chemokines_mtx_all_number <- unlist(lapply(
    rownames(cytokines_chemokines_mtx_all),
    function(X) {
        #         if (X %in% rownames(cytokines_chemokines_mtx_no))
        #             1
        if (X %in% rownames(cytokines_chemokines_mtx_down)) {
            2
        } else if (X %in% rownames(cytokines_chemokines_mtx_up)) {
            1
        }
    }
))


cytokines_chemokines_mtx_all

gene_pvalue_t <- res[rownames(cytokines_chemokines_mtx_all), "padj", drop = F]
row_split <- data.frame(
    group = factor(cytokines_chemokines_mtx_all_number, labels = c("Up", "Down")),
    row.names = rownames(cytokines_chemokines_mtx_all)
)
column_split <- data.frame(
    group = factor(c(rep("Non-CA", 3), rep("CA", 3)), levels = c("Non-CA", "CA")),
    row.names = colnames(cytokines_chemokines_mtx_all)
)
pvalue_vec <- unlist(lapply(
    gene_pvalue_t[, 1],
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
pp_cytokines_chemokines <- ComplexHeatmap::Heatmap(
    cytokines_chemokines_mtx_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = T, 
    show_column_dend = F,
    cluster_rows = T,
    show_row_dend = T,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    column_split = column_split,
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    row_names_gp = gpar(fontsize = 15),
    top_annotation = HeatmapAnnotation(
        Group = factor(
            annotation_col$Group,
            levels = c("Non-CA", "CA")
        ),
        col = list(
            Group = c("Non-CA" = "grey", "CA" = "orange")
        ),
        show_legend = F
    ),
    left_annotation = rowAnnotation(
        Significant = factor(
            pvalue_vec,
            levels = c(
                "Not significant",
                "P-value < 0.05",
                "P-value < 0.01",
                "P-value < 0.001"
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
pp_cytokines_chemokines

pdf("~/Desktop/Cryo/Paper_Result_Plot/plot/cytokines_chemokines_heatmap.pdf", width = 5.5, height = 15)
draw(setFont(pp_cytokines_chemokines))
dev.off()
