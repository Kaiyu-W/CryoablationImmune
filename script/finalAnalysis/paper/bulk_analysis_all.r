setwd("/mnt/e/Cryo-TCR/bulk/res_counts")
# source("../../server/auto/utilities.r")
library(DESeq2)
library(ComplexHeatmap)

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

# resultsNames(data)
res <- results(data, contrast = c("condition", "CA", "Non-CA"), alpha = 0.05)
summary(res)

res <- res[order(res$padj), ]
# diff_gene_deseq2 <-subset(res, padj < 0.05 & log2FoldChange > 1)
diff_gene_deseq2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
diff_gene_deseq2_name <- rownames(diff_gene_deseq2)
up_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange > 0))
down_genes <- rownames(subset(res, padj < 0.05 & log2FoldChange < 0))
nosignificant_gene <- setdiff(rownames(res), c(up_genes, down_genes))
resdata <- merge(as.data.frame(res), as.data.frame(counts(data, normalize = TRUE)), by = "row.names", sort = FALSE)

# save_res <- resdata[resdata$Row.names %in% diff_gene_deseq2_name, -c(4:5)][order(diff_gene_deseq2$log2FoldChange, decreasing = T), ]
#
# write.csv(save_res, file = 'Up_genelist.csv', row.names = F)


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
# write.csv(cytokines, file = '~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines.csv', quote = F)
# write.csv(cytokines_chemokines, file = '~/Desktop/Cryo/Paper_Result_Plot/bulk_cytokines_chemokines.csv', quote = F)
#######

chemokine_mtx <- subsetMtx(chemokine)
cytokine_mtx <- subsetMtx(cytokine)
interferon_mtx <- subsetMtx(interferon)
custom_mtx <- subsetMtx(custom_genes)
chemokines_mtx <- subsetMtx(chemokines)
cytokines_mtx <- subsetMtx(cytokines)
cytokines_chemokines_mtx <- subsetMtx(cytokines_chemokines)

# pheatmap::pheatmap(chemokine_mtx, cluster_cols = F)
# pheatmap::pheatmap(cytokine_mtx, cluster_cols = F)
# pheatmap::pheatmap(interferon_mtx, cluster_cols = F)
# pheatmap::pheatmap(custom_mtx, cluster_cols = F)
# pheatmap::pheatmap(chemokines_mtx, cluster_cols = F)
# pheatmap::pheatmap(cytokines_mtx, cluster_cols = F)
# pheatmap::pheatmap(cytokines_chemokines_mtx, cluster_cols = F)

diff_mtx <- subsetMtx(diff_gene_deseq2_name)

diff_targets <- intersect(up_genes, c(chemokine, cytokine))
difftarget_mtx <- subsetMtx(diff_targets)

diff_targets_ifn <- intersect(up_genes, interferon)
difftarget_ifn_mtx <- subsetMtx(diff_targets_ifn)

diff_targets_custom <- intersect(up_genes, custom_genes)
difftargets_custom_mtx <- subsetMtx(diff_targets_custom)

RWB_color <- colorRampPalette(colors = c("blue", "white", "red"))(100)

# chemokine cytokine
chemokine_cytokine_mtx <- subsetMtx(c(chemokine, cytokine))

annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(chemokine_cytokine_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(chemokine_cytokine_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
chemokine_cytokine_mtx_up <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
chemokine_cytokine_mtx_down <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
chemokine_cytokine_mtx_no <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

ComplexHeatmap::Heatmap(
    chemokine_cytokine_mtx_up[, c(4:6, 1:3)], # put control in front
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

# interferon
annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(interferon_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(interferon_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
interferon_mtx_up <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
interferon_mtx_down <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
interferon_mtx_no <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

xx <- ComplexHeatmap::Heatmap(
    interferon_mtx_no[, c(4:6, 1:3)], # put control in front
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
interferon_mtx_all <- rbind(
    interferon_mtx_no[rownames(xx@matrix), c(4:6, 1:3)],
    interferon_mtx_up[, c(4:6, 1:3)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(interferon_mtx_all) %in% rownames(interferon_mtx_up)), labels = c("No", "Up")), row.names = rownames(interferon_mtx_all))
ComplexHeatmap::Heatmap(
    interferon_mtx_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split, gap = unit(0.2, "inches"),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)

# custom
annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(custom_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(custom_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
custom_mtx_up <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
custom_mtx_down <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
custom_mtx_no <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

xx <- ComplexHeatmap::Heatmap(
    custom_mtx_no[, c(4:6, 1:3)], # put control in front
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
custom_mtx_all <- rbind(
    custom_mtx_no[rownames(xx@matrix), c(4:6, 1:3)],
    custom_mtx_up[, c(4:6, 1:3)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(custom_mtx_all) %in% rownames(custom_mtx_up)), labels = c("No", "Up")), row.names = rownames(custom_mtx_all))
ComplexHeatmap::Heatmap(
    custom_mtx_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split, gap = unit(0.2, "inches"),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)

# chemokines
annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(chemokines_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(chemokines_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
chemokines_mtx_up <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
chemokines_mtx_down <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
chemokines_mtx_no <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

xx <- ComplexHeatmap::Heatmap(
    chemokines_mtx_no[, c(4:6, 1:3)], # put control in front
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
chemokines_mtx_all <- rbind(
    chemokines_mtx_no[rownames(xx@matrix), c(4:6, 1:3)],
    chemokines_mtx_up[, c(4:6, 1:3)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(chemokines_mtx_all) %in% rownames(chemokines_mtx_up)), labels = c("No", "Up")), row.names = rownames(chemokines_mtx_all))
ComplexHeatmap::Heatmap(
    chemokines_mtx_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split, gap = unit(0.2, "inches"),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)

# cytokines
annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(difftarget_mtx)
)
annotation_row <- data.frame(
    Significant = sapply(rownames(cytokines_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"),
    row.names = rownames(cytokines_mtx)
)
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
cytokines_mtx_up <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == "up"], rownames(annotation_col), drop = F]
cytokines_mtx_down <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == "down"], rownames(annotation_col), drop = F]
cytokines_mtx_no <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]


xx <- ComplexHeatmap::Heatmap(
    cytokines_mtx_no[, c(4:6, 1:3)], # put control in front
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
cytokines_mtx_all <- rbind(
    cytokines_mtx_no[rownames(xx@matrix), c(4:6, 1:3)],
    cytokines_mtx_down[, c(4:6, 1:3)],
    cytokines_mtx_up[, c(4:6, 1:3)]
)
cytokines_mtx_all_number <- unlist(lapply(
    rownames(cytokines_mtx_all),
    function(X) {
        if (X %in% rownames(cytokines_mtx_no)) {
            1
        } else if (X %in% rownames(cytokines_mtx_down)) {
            2
        } else if (X %in% rownames(cytokines_mtx_up)) {
            3
        }
    }
))
row_split <- data.frame(group = factor(cytokines_mtx_all_number, labels = c("No", "Down", "Up")), row.names = rownames(cytokines_mtx_all))
ComplexHeatmap::Heatmap(
    cytokines_mtx_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = F,
    row_split = row_split, gap = unit(0.2, "inches"),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)

# pp_cytokines <- ComplexHeatmap::Heatmap(cytokines_mtx_all, # put control in front
#                                         name = "Z-score",
#                                         col = RWB_color,
#                                         rect_gp = gpar(col = "white", lwd = 1),
#                                         cluster_columns = F,
#                                         cluster_rows = F,
#                                         show_row_dend = F,
#                                         show_column_names = F,
#                                         show_row_names = T,
#                                         row_split = row_split,
#                                         gap = unit(0.2, "inches"),
#                                         row_names_gp = gpar(fontsize = 2.5),
#                                         top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
#                                                                            col = list(Group = c("Non-CA" = "grey", "CA" = "orange")))
# )
# write.csv(row_split[order(row_split$group),,drop=F], file = 'E://Cryo-TCR/bulk/cytokines_list_rm3.csv')
# pdf("E://Cryo-TCR/bulk/cytokines_heatmap.pdf", width = 3.5, height = 15)
# draw(pp_cytokines)
# dev.off()

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
row_split <- data.frame(
    group = factor(cytokines_chemokines_mtx_all_number, labels = c("Up", "Down")),
    row.names = rownames(cytokines_chemokines_mtx_all)
)
column_split <- data.frame(
    group = factor(c(rep("Non-CA", 3), rep("CA", 3)), levels = c("Non-CA", "CA")),
    row.names = colnames(cytokines_chemokines_mtx_all)
)
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
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
pp_cytokines_chemokines

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
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange")),
        show_legend = F
    ),
    left_annotation = rowAnnotation(
        Significant = factor(pvalue_vec, levels = c("Not significant", "P-value < 0.05", "P-value < 0.01", "P-value < 0.001")),
        col = list(Significant = c(
            "Not significant" = "#D3D3D3",
            "P-value < 0.05" = "#B4EEB4",
            "P-value < 0.01" = "#3CB371",
            "P-value < 0.001" = "#006400"
        )),
        gp = gpar(col = "black", lwd = 1),
        show_annotation_name = F
    )
)
pp_cytokines_chemokines

# pdf("E://Cryo-TCR/bulk/cytokines_chemokines_heatmap.pdf", width = 6, height = 15)
# draw(pp_cytokines_chemokines)
# dev.off()


#############################
############################# GSVA for pathway
library(GSVA)
library(GSEABase)
RWB_color <- grDevices::colorRampPalette(colors = c("blue", "white", "red"))(100)
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

annotation_col <- data.frame(
    Group = c(rep("Non-CA", 3), rep("CA", 3)),
    row.names = colnames(count_mtx)
)

# interferon
pathway_interferon <- names(GSVA_diff)[grep("interferon", names(GSVA_diff), ignore.case = T)]
pathway_up_interferon <- pathway_up[grep("interferon", pathway_up, ignore.case = T)]

interferon_scale_mtx <- GSVA_scale[pathway_interferon, c(4:6, 1:3)]
rownames(interferon_scale_mtx) <- tolower(sub("^GO.._", "", rownames(interferon_scale_mtx)))
up_index_tmp <- which(tolower(sub("^GO.._", "", pathway_up_interferon)) %in% rownames(interferon_scale_mtx))
interferon_scale_mtx <- interferon_scale_mtx[c(up_index_tmp, setdiff(1:nrow(interferon_scale_mtx), up_index_tmp)), ]
row_split <- data.frame(group = factor(1:nrow(interferon_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(interferon_scale_mtx))
rownames(interferon_scale_mtx) <- gsub("_", " ", rownames(interferon_scale_mtx))
pp_pathway_interferon <- ComplexHeatmap::Heatmap(
    interferon_scale_mtx, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 3),
    #                         row_names_max_width = unit(1, "inches"),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    gap = unit(0.2, "inches"),
    row_names_gp = gpar(fontsize = 10),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
pp_pathway_interferon

# chemokine
pathway_chemokine <- names(GSVA_diff)[grep("chemokine", names(GSVA_diff), ignore.case = T)]
pathway_up_chemokine <- pathway_up[grep("chemokine", pathway_up, ignore.case = T)]

chemokine_scale_mtx <- GSVA_scale[pathway_chemokine, c(4:6, 1:3)]
rownames(chemokine_scale_mtx) <- tolower(sub("^GO.._", "", rownames(chemokine_scale_mtx)))
up_index_tmp <- which(tolower(sub("^GO.._", "", pathway_up_chemokine)) %in% rownames(chemokine_scale_mtx))
chemokine_scale_mtx <- chemokine_scale_mtx[c(up_index_tmp, setdiff(1:nrow(chemokine_scale_mtx), up_index_tmp)), ]
row_split <- data.frame(group = factor(1:nrow(chemokine_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(chemokine_scale_mtx))
rownames(chemokine_scale_mtx) <- gsub("_", " ", rownames(chemokine_scale_mtx))
pp_pathway_chemokine <- ComplexHeatmap::Heatmap(
    chemokine_scale_mtx, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 3),
    #                         row_names_max_width = unit(1, "inches"),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    gap = unit(0.2, "inches"),
    row_names_gp = gpar(fontsize = 10),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
pp_pathway_chemokine

# cytokine
pathway_cytokine <- names(GSVA_diff)[grep("cytokine", names(GSVA_diff), ignore.case = T)]
pathway_up_cytokine <- pathway_up[grep("cytokine", pathway_up, ignore.case = T)]
cytokine_scale_mtx <- GSVA_scale[pathway_cytokine, c(4:6, 1:3)]
rownames(cytokine_scale_mtx) <- tolower(sub("^GO.._", "", rownames(cytokine_scale_mtx)))
up_index_tmp <- which(tolower(sub("^GO.._", "", pathway_up_cytokine)) %in% rownames(cytokine_scale_mtx))
cytokine_scale_mtx <- cytokine_scale_mtx[c(up_index_tmp, setdiff(1:nrow(cytokine_scale_mtx), up_index_tmp)), ]
row_split <- data.frame(group = factor(1:nrow(cytokine_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(cytokine_scale_mtx))
rownames(cytokine_scale_mtx) <- gsub("_", " ", rownames(cytokine_scale_mtx))
pp_pathway_cytokine <- ComplexHeatmap::Heatmap(
    cytokine_scale_mtx, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 3),
    #                         row_names_max_width = unit(1, "inches"),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    gap = unit(0.2, "inches"),
    row_names_gp = gpar(fontsize = 10),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
pp_pathway_cytokine

# inflammatory
pathway_inflammatory <- names(GSVA_diff)[grep("inflammatory", names(GSVA_diff), ignore.case = T)]
pathway_up_inflammatory <- pathway_up[grep("inflammatory", pathway_up, ignore.case = T)]

inflammatory_scale_mtx <- GSVA_scale[pathway_inflammatory, c(4:6, 1:3)]
rownames(inflammatory_scale_mtx) <- tolower(sub("^GO.._", "", rownames(inflammatory_scale_mtx)))
up_index_tmp <- which(tolower(sub("^GO.._", "", pathway_up_inflammatory)) %in% rownames(inflammatory_scale_mtx))
inflammatory_scale_mtx <- inflammatory_scale_mtx[c(up_index_tmp, setdiff(1:nrow(inflammatory_scale_mtx), up_index_tmp)), ]
row_split <- data.frame(group = factor(1:nrow(inflammatory_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(inflammatory_scale_mtx))
rownames(inflammatory_scale_mtx) <- gsub("_", " ", rownames(inflammatory_scale_mtx))
pp_pathway_inflammatory <- ComplexHeatmap::Heatmap(
    inflammatory_scale_mtx, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 3),
    #                         row_names_max_width = unit(1, "inches"),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    show_row_names = T,
    row_split = row_split,
    gap = unit(0.2, "inches"),
    row_names_gp = gpar(fontsize = 10),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
    )
)
pp_pathway_inflammatory

# pdf("E://Cryo-TCR/bulk/pathway_temp.pdf", width = 10, height = 10)
# draw(pp_pathway_interferon)
# draw(pp_pathway_chemokine)
# draw(pp_pathway_cytokine)
# draw(pp_pathway_inflammatory)
# dev.off()

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
if (length(pathway_up_custom) == length(pathway_custom)) {
    custom_scale_mtx <- GSVA_scale[pathway_custom, c(4:6, 1:3)]
    rownames(custom_scale_mtx) <- gsub("_", " ", tolower(sub("^GO.._", "", rownames(custom_scale_mtx))))
    ComplexHeatmap::Heatmap(
        custom_scale_mtx, # put control in front
        name = "Z-score",
        col = RWB_color,
        rect_gp = gpar(col = "white", lwd = 3),
        row_names_max_width = unit(20, "inches"),
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_names = F,
        show_row_names = T,
        #                             row_split = row_split,
        gap = unit(0.2, "inches"),
        left_annotation = rowAnnotation(pvalue = pathway_pvalue_t),
        top_annotation = HeatmapAnnotation(
            Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
            col = list(Group = c("Non-CA" = "grey", "CA" = "orange"))
        )
    )
}

column_split <- data.frame(
    group = factor(c(rep("Non-CA", 3), rep("CA", 3)), levels = c("Non-CA", "CA")),
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
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        col = list(Group = c("Non-CA" = "grey", "CA" = "orange")),
        show_legend = F
    ),
    # left_annotation = rowAnnotation(pvalue = pvalue_vec)
    left_annotation = rowAnnotation(
        Significant = factor(
            pvalue_vec,
            levels = c("Not significant", "P-value < 0.05", "P-value < 0.01", "P-value < 0.001", "P-value < 0.0001")
        ),
        col = list(Significant = c(
            "Not significant" = "#D3D3D3",
            "P-value < 0.05" = "#B4EEB4",
            "P-value < 0.01" = "#3CB371",
            "P-value < 0.001" = "#006400"
        )),
        gp = gpar(col = "black", lwd = 1),
        show_annotation_name = F
    )
)
pp_pathway_custom

# pdf("E://Cryo-TCR/bulk/pathway_custom_heatmap.pdf", width = 6, height = 4.5)
# draw(pp_pathway_custom)
# dev.off()

# HallMark pathway
# GSVA for pathway
library(GSVA)
library(GSEABase)
RWB_color <- grDevices::colorRampPalette(colors = c("blue", "white", "red"))(100)
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
row_split <- data.frame(
    group = factor(GSVA_scale2_all_number, labels = c("No", "Up")),
    row.names = rownames(GSVA_scale2_all)
)
column_split <- data.frame(
    group = factor(c(rep("Non-CA", 3), rep("CA", 3)), levels = c("Non-CA", "CA")),
    row.names = colnames(GSVA_scale2_all)
)

rownames(GSVA_scale2_all) <- tolower(rownames(GSVA_scale2_all))
pp_pathway_custom2 <- ComplexHeatmap::Heatmap(
    GSVA_scale2_all, # put control in front
    name = "Z-score",
    col = RWB_color,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = T, 
    show_column_dend = F,
    cluster_rows = T,
    show_heatmap_legend = F,
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
        Significant = factor(pvalue_vec2, levels = c(
            "Not significant",
            "P-value < 0.05",
            "P-value < 0.01",
            "P-value < 0.001",
            "P-value < 0.0001"
        )),
        col = list(Significant = c(
            "Not significant" = "#D3D3D3",
            "P-value < 0.05" = "#B4EEB4",
            "P-value < 0.01" = "#3CB371",
            "P-value < 0.001" = "#006400"
        )),
        gp = gpar(col = "black", lwd = 1),
        show_annotation_name = F,
        show_legend = F
    )
)
pp_pathway_custom2

# pdf("E://Cryo-TCR/bulk/pathway_hallmark_heatmap.pdf", width = 6, height = 10)
# draw(pp_pathway_custom2)
# dev.off()
