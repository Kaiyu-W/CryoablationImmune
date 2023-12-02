setwd("E://Cryo-TCR/bulk/res_counts")
source("../../server/auto/utilities.r")

count<-read.table("genome_counts.tsv", sep = "\t", header = T, row.names = 1)
idx<-read.table("transcriptome_counts.tsv", sep = "\t", header = T, row.names = 1)

library(DESeq2)
data <- DESeqDataSetFromMatrix(
    countData = count,
    colData = data.frame(
        row.names = colnames(count),
        condition = factor(c(rep('Cryo', 4), rep('NonCryo', 3)), levels = c('Cryo', 'NonCryo'))
        ),
    design = ~ condition)
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
res <- results(data, contrast = c("condition", "Cryo", "NonCryo"), alpha = 0.05)
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
# 
# bulk_GO <- my_GO(diff_gene_deseq2_name, ont = 'BP', Simplify = T, return_plot = F, return_res = T)
# which(bulk_GO@result$Description %in% grep("interferon", bulk_GO@result$Description, value = T))
# barplot(bulk_GO, showCategory = 35, font.size = 10, title = "GO BP")
# dotplot(bulk_GO, showCategory = 35, font.size = 10, title = "GO BP")
# bulk_GO@result[which(bulk_GO@result$Description %in% grep("interferon", bulk_GO@result$Description, value = T)), ]
# # ID                  Description GeneRatio   BgRatio       pvalue     p.adjust       qvalue
# # GO:0032609 GO:0032609  interferon-gamma production    16/447 128/23328 3.316056e-09 1.345825e-07 9.613962e-08
# # GO:0034341 GO:0034341 response to interferon-gamma    12/447 137/23328 1.355453e-05 2.228902e-04 1.592226e-04
# # geneID Count
# # GO:0032609 12481/20737/21940/81897/66824/21947/27218/19260/15900/12501/50931/74131/18173/60505/16162/18578    16
# # GO:0034341                         15018/24088/18126/20296/15900/20304/16365/66141/24047/22376/20308/18173    12
# bulk_GO_IFN1 <- strsplit('12481/20737/21940/81897/66824/21947/27218/19260/15900/12501/50931/74131/18173/60505/16162/18578', split = '/')[[1]]
# bulk_GO_IFN2 <- strsplit('15018/24088/18126/20296/15900/20304/16365/66141/24047/22376/20308/18173', split = '/')[[1]]
# bulk_GO_IFN <- unique(c(bulk_GO_IFN1, bulk_GO_IFN2))
# IFN_gene <- bitr(bulk_GO_IFN, fromType = "ENTREZID", toType = 'SYMBOL', OrgDb = OrgDb)$SYMBOL
# 
# res_GOBP <- bulk_GO@result[,-c(1,8,9)]
# res_IFN <- resdata[resdata$Row.names%in%IFN_gene, ]
# write.csv(res_GOBP, file = 'GOBP_result.csv', row.names = F)
# write.csv(res_IFN, file = 'IFN_result.csv', row.names = F)
# ggsave(filename = 'GOBP_result_bar.png', plot = barplot(bulk_GO, showCategory = 35, font.size = 5, title = "GO BP"))
# ggsave(filename = 'GOBP_result_dot.png', plot = dotplot(bulk_GO, showCategory = 35, font.size = 5, title = "GO BP"))
# 
# library(pheatmap)
# library(RColorBrewer)
# 
# pheatmap::pheatmap(log(resdata[resdata[,1]%in%diff_gene_deseq2,(ncol(resdata)-6):ncol(resdata)]+1))
# 
# rl_data
# sampleDist <- dist(t(assay(rl_data)))
# 
# #poisd <- PoiClaClu::PoissonDistance(t(counts(dds))) 
# sampleDistMatrix <- as.matrix(sampleDist)  #样品间距离的矩阵
# rownames(sampleDistMatrix) <- paste0(rl_data$cell, "-", rl_data$dex)
# colnames(sampleDistMatrix) <- NULL
# head(sampleDistMatrix)  #样品间距离的数据???
# colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# pheatmap(sampleDistMatrix,
#          clustering_distance_rows=sampleDist,
#          clustering_distance_cols=sampleDist,
#          color = colors)
# 
# library(ggplot2)
# library(ggrepel)
# pcadata <- plotPCA(rl_data, intgroup = "condition", returnData=TRUE)
# percentVar <- round(100*attr(pcadata, "percentVar"), 1)
# ggplot(pcadata, aes(PC1, PC2, color=condition)) + 
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     geom_text_repel(aes(PC1, PC2,color=condition,label=colnames(rl_data)),size=3) +
#     theme_bw()
# 
# topGene <- rownames(res)[which.min(res$padj)]  #padj最小的一个基???
# plotCounts(data, gene=topGene, intgroup=c("condition"))  #画出这个基因的标准化后的表达???
# #以散点图的形式画出这个基因在各样本中的表达量
# Plot <- plotCounts(data, gene=topGene, intgroup=c("condition"), returnData=TRUE)
# ggplot(Plot, aes(x=condition, y=count, color=condition)) +  
#     scale_y_log10() +
#     geom_point(position=position_jitter(width=.1,height=0), size=3)  
# ggplot(Plot, aes(x=condition, y=count, fill=condition)) +
#     scale_y_log10() +
#     geom_dotplot(binaxis="y", stackdir="center")
# ggplot(Plot, aes(x=condition, y=count, color=condition, group=condition)) +
#     scale_y_log10() + geom_point(size=3) + geom_line()
# 
# 
# exprSet <- idx[,-1]
# exprSet <- log2(exprSet+1)
# cg <-names(tail(sort(apply(exprSet,1,sd)),1000))
# n <-t(scale(t(exprSet[cg,])))
# n[n>2] <-2
# n[n<-2] <-2
# ac <- data.frame(sample=factor(c(rep('Cryo', 4), rep('NonCryo', 3)), levels = c('Cryo', 'NonCryo')))
# rownames(ac) <- colnames(n)
# 
# ##rlog ??׼??
# rld <- rlog(data)
# ##??ȡ
# exprMatrix_rlog <- assay(rld)
# 
# #??׼????raw count ?ıȽ?
# 
# par(cex= 0.7)
# n.samples <- ncol(exprSet)
# if(n.samples > 40) 
#     par(cex = 0.5)
# cols <- rainbow(n.samples*1.2)
# par(mfrow = c(2,2))
# boxplot(exprSet, col = cols, main = "expression value",las=2)
# boxplot(exprMatrix_rlog, col = cols, main = "expression value", las = 2)
# hist(as.matrix(exprSet))
# hist(exprMatrix_rlog)
# 
# 
# ##DEseq2 ????dds
# dds <- DESeq(data)
# res <- results(dds)
# res <- res[order(res$padj), ]
# DEG <- as.data.frame(res)
# 
# ##ȥ??NA
# DEG <- na.omit(DEG)
# 
# ##??ͼ
# library(pheatmap)
# choose_gene <- head(rownames(DEG), 100)  ##50 maybe better
# choose_matrix <- exprSet[choose_gene, ]
# choose_matrix <- t(scale(t(choose_matrix)))
# pheatmap(choose_matrix, show_rownames = F, show_colnames = T, annotation_col = ac)
# 
# 
# 
# # MA plot
# library(apeglm)  
# resultsNames(data)  #??һ??Ҫshrink??ά??;shrink???ݸ??ӽ???,????һ??stat??????δ?ı?padj?????ı???foldchange
# res_shrink <- lfcShrink(data, coef="condition_NonCryo_vs_Cryo", type="apeglm") #???Ƽ?apeglm?㷨;????resultsNames(dds)?ĵ?5??ά?ȣ?coef=5??Ҳ??ֱ??""ָ??;apeglm??allow contrast??????Ҫָ??coef
# plotMA(res_shrink, ylim=c(-10,10), alpha=0.1, main="MA plot: ")
# 
# library(ggplot2)
# voldata <- resdata[!is.na(resdata$padj), ]
# 
# Up_log2FoldChange <- voldata$log2FoldChange >= 1 & voldata$padj < 0.05
# Down_log2FoldChange <- voldata$log2FoldChange <= -1 & voldata$padj < 0.05
# Not_Significant <- voldata$padj >= 0.05 | (voldata$log2FoldChange > -1 & voldata$log2FoldChange < 1)
# Significant_resdata <- unlist(apply(
#     data.frame(Up = Up_log2FoldChange, Down = Down_log2FoldChange, "Not\ Significant" = Not_Significant), 
#     1,
#     function(x)
#         names(x)[x]
#     ))
#     
# voldata$significant <- Significant_resdata
# 
# ggplot(data = voldata, aes(x = log2FoldChange, y= -1 * log10(padj))) +
#     geom_point(aes(color = significant)) +
#     scale_color_manual(values = c("#546de5", "#d2dae2", "#ff4757")) + 
#     labs(title = "Volcano Plot: ", x = expression(log[2](FC), y = expression(-log[10](padj)))) +
#     geom_hline(yintercept=1.3, linetype=4) +  #??????,????0.05????
#     geom_vline(xintercept=c(-1, 1), linetype = 4) +
#     theme_bw() + theme(panel.grid = element_blank())  #?????????߾?Ϊ?հ?


#############################
# search the chemokine and cytokine pathway genes to plot heatmap
setwd("C://DeepTL_dev-main/pathway/GO_mus/direct-experimentevidence/")
go_bp <- read.table("GO_BP.tsv", sep= '\t', header = T)
go_cc <- read.table("GO_CC.tsv", sep= '\t', header = T)
go_mf <- read.table("GO_MF.tsv", sep= '\t', header = T)

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
chemokine <- searchGenefromPathway("chemokine")
cytokine <- searchGenefromPathway("cytokine")
interferon_beta <- searchGenefromPathway2("interferon", 'beta')
interferon_alpha <- searchGenefromPathway2("interferon", 'alpha')
interferon_I <- searchGenefromPathway2("interferon", 'type I')
interferon <- unique(c(interferon_beta, interferon_alpha, interferon_I))
custom_genes <- c('Cxcr4', 'Ccl19', 'Cxcl9', 'Cxcl3', 'Cxcl10', 'Cxcl11', 'Cxcl12', 'Stat1', 'Stat2', 'Ifnl2', 'Ifnb1', 'Ifna1', 'Ifnar1', 'Ifnar2', 'Ifna2', 'Ifna4', 'Ifng', 'Xcl1', 'Ccl5', 'Ccl9', 'Ccr5', 'Ifngr1')
count_mtx <- data@assays@data$counts

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
cytokines_human <- scan("E://Cryo-TCR/bulk/cytokines_human.txt", what = 'c')
cytokines_receptor_human <- scan("E://Cryo-TCR/bulk/cytokines_receptor_human.txt", what = 'c')
cytokines_receptor_mouse <- scan("E://Cryo-TCR/bulk/cytokines_receptor_mouse.txt", what = 'c')

library(biomaRt)
human <- useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "https://may2021.archive.ensembl.org/")
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl", host = "https://may2021.archive.ensembl.org/")

cytokines_mouse <- getLDS(attributes = "hgnc_symbol", 
                          filters = "hgnc_symbol",
                          values = cytokines_human, 
                          mart = human,
                          attributesL = "mgi_symbol",
                          martL = mouse,
                          uniqueRows = TRUE)
cytokines_receptor_mouse2 <- getLDS(attributes = "hgnc_symbol",
                                    filters = "hgnc_symbol",
                                    values = cytokines_receptor_human, 
                                    mart = human,
                                    attributesL = "mgi_symbol",
                                    martL = mouse,
                                    uniqueRows = TRUE)

cytokines <- unique(c(cytokines_mouse$MGI.symbol, cytokines_receptor_mouse2$MGI.symbol, cytokines_receptor_mouse))
#######


pathway_Ifnar1 <- searchPathwayfromGene('Ifnar1')

subsetMtx <- function(geneset) {
    geneset <- unique(geneset)
    mtx0 <- log1p(count_mtx[geneset[geneset %in% rownames(count_mtx)], ])
    mtx <- t(apply(mtx0, 1, scale))
    colnames(mtx) <- colnames(mtx0)
    mtx
}
chemokine_mtx <- subsetMtx(chemokine)
cytokine_mtx <- subsetMtx(cytokine)
interferon_mtx <- subsetMtx(interferon)
custom_mtx <- subsetMtx(custom_genes)
chemokines_mtx <- subsetMtx(chemokines)
cytokines_mtx <- subsetMtx(cytokines)

pheatmap::pheatmap(chemokine_mtx, cluster_cols = F)
pheatmap::pheatmap(cytokine_mtx, cluster_cols = F)
pheatmap::pheatmap(interferon_mtx, cluster_cols = F)
pheatmap::pheatmap(custom_mtx, cluster_cols = F)

diff_mtx <- subsetMtx(diff_gene_deseq2_name)
pheatmap::pheatmap(diff_mtx, cluster_cols = F)

diff_targets <- intersect(up_genes, c(chemokine, cytokine))
difftarget_mtx <- subsetMtx(diff_targets)
pheatmap::pheatmap(difftarget_mtx, 
                   cluster_cols = F, 
                   angle_col = 45, 
                   annotation_col = data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                                               row.names = colnames(difftarget_mtx))
                   )

diff_targets_ifn <- intersect(up_genes, interferon)
difftarget_ifn_mtx <- subsetMtx(diff_targets_ifn)
pheatmap::pheatmap(difftarget_ifn_mtx, 
                   cluster_cols = F, 
                   angle_col = 45, 
                   annotation_col = data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                                               row.names = colnames(difftarget_ifn_mtx))
)

diff_targets_custom <- intersect(up_genes, custom_genes)
difftargets_custom_mtx <- subsetMtx(diff_targets_custom)
pheatmap::pheatmap(difftargets_custom_mtx, 
                   cluster_cols = F, 
                   angle_col = 45, 
                   annotation_col = data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                                               row.names = colnames(difftargets_custom_mtx))
)

chemokine_cytokine_mtx <- subsetMtx(c(chemokine, cytokine))
annotation_col <- data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                             row.names = colnames(difftarget_mtx))
annotation_row <- data.frame(Significant = sapply(rownames(chemokine_cytokine_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"), 
                             row.names = rownames(chemokine_cytokine_mtx))
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
chemokine_cytokine_mtx_up <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == 'up'], rownames(annotation_col), drop = F]
chemokine_cytokine_mtx_down <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == 'down'], rownames(annotation_col), drop = F]
chemokine_cytokine_mtx_no <- chemokine_cytokine_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]


annotation_col <- data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                             row.names = colnames(difftargets_custom_mtx))
annotation_row <- data.frame(Significant = sapply(rownames(interferon_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"), 
                             row.names = rownames(interferon_mtx))
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
interferon_mtx_up <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == 'up'], rownames(annotation_col), drop = F]
interferon_mtx_down <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == 'down'], rownames(annotation_col), drop = F]
interferon_mtx_no <- interferon_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

annotation_col <- data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                             row.names = colnames(difftarget_ifn_mtx))
annotation_row <- data.frame(Significant = sapply(rownames(custom_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"), 
                             row.names = rownames(custom_mtx))
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
custom_mtx_up <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == 'up'], rownames(annotation_col), drop = F]
custom_mtx_down <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == 'down'], rownames(annotation_col), drop = F]
custom_mtx_no <- custom_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

annotation_col <- data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                             row.names = colnames(difftarget_mtx))
annotation_row <- data.frame(Significant = sapply(rownames(chemokines_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"), 
                             row.names = rownames(chemokines_mtx))
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
chemokines_mtx_up <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == 'up'], rownames(annotation_col), drop = F]
chemokines_mtx_down <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == 'down'], rownames(annotation_col), drop = F]
chemokines_mtx_no <- chemokines_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]

annotation_col <- data.frame(Group = c(rep("Cryo", 4), rep("NonCryo", 3)), 
                             row.names = colnames(difftarget_mtx))
annotation_row <- data.frame(Significant = sapply(rownames(cytokines_mtx), function(x) if (x %in% up_genes) "up" else if (x %in% down_genes) "down" else "Not significant"), 
                             row.names = rownames(cytokines_mtx))
annotation_row <- annotation_row[order(annotation_row$Significant), , drop = F]
cytokines_mtx_up <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == 'up'], rownames(annotation_col), drop = F]
cytokines_mtx_down <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == 'down'], rownames(annotation_col), drop = F]
cytokines_mtx_no <- cytokines_mtx[rownames(annotation_row)[annotation_row$Significant == "Not significant"], rownames(annotation_col), drop = F]


RWB_color <- colorRampPalette(colors = c("blue","white","red"))(100)
pheatmap::pheatmap(chemokine_cytokine_mtx[rownames(annotation_row), rownames(annotation_col)],
                   cluster_cols = F,
                   cluster_rows = F, 
                   color = RWB_color,
                   show_colnames = F,
                   annotation_col = annotation_col,
                   annotation_row = annotation_row
                   )
library(ComplexHeatmap)
p_up <- ComplexHeatmap::pheatmap(chemokine_cytokine_mtx_up, 
                   cluster_cols = F, 
                   cluster_rows = T,
                   color = RWB_color,
                   rect_gp = gpar(col = "white", lwd = 2),
                   annotation_col = annotation_col,
)
p_down <- ComplexHeatmap::pheatmap(chemokine_cytokine_mtx_down, 
                   cluster_cols = F, 
                   cluster_rows = F,
                   color = RWB_color,
                   annotation_col = annotation_col,
)
p_no <- ComplexHeatmap::pheatmap(chemokine_cytokine_mtx_no, 
                   cluster_cols = F, 
                   cluster_rows = T,
                   color = RWB_color,
                   annotation_col = annotation_col,
)
p_up
p_down
p_no

# chemokine cytokine
ComplexHeatmap::Heatmap(chemokine_cytokine_mtx_up[, c(5:7, 1:4)], # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 2), 
                        cluster_columns = F, 
                        show_row_dend = F, 
                        show_column_names = F,
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# interferon
ComplexHeatmap::Heatmap(interferon_mtx_up[, c(5:7, 1:4)], # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 2), 
                        cluster_columns = F, 
                        show_row_dend = F, 
                        show_column_names = F,
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
xx <- ComplexHeatmap::Heatmap(interferon_mtx_no[, c(5:7, 1:4)], # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 2), 
                        cluster_columns = F, 
                        show_row_dend = F, 
                        show_column_names = F,
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
interferon_mtx_all <- rbind(
    interferon_mtx_no[rownames(xx@matrix), c(5:7, 1:4)],
    interferon_mtx_up[, c(5:7, 1:4)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(interferon_mtx_all) %in% rownames(interferon_mtx_up)), labels = c("No", "Up")), row.names = rownames(interferon_mtx_all))
ComplexHeatmap::Heatmap(interferon_mtx_all, # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 1), 
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split, gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# custom
xx <- ComplexHeatmap::Heatmap(custom_mtx_no[, c(5:7, 1:4)], # put control in front
                              name = "Z-score",
                              col = RWB_color, 
                              rect_gp = gpar(col = "white", lwd = 2), 
                              cluster_columns = F, 
                              show_row_dend = F, 
                              show_column_names = F,
                              top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                                 col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
custom_mtx_all <- rbind(
    custom_mtx_no[rownames(xx@matrix), c(5:7, 1:4)],
    custom_mtx_up[, c(5:7, 1:4)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(custom_mtx_all) %in% rownames(custom_mtx_up)), labels = c("No", "Up")), row.names = rownames(custom_mtx_all))
ComplexHeatmap::Heatmap(custom_mtx_all, # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 1), 
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split, gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# chemokines
xx <- ComplexHeatmap::Heatmap(chemokines_mtx_no[, c(5:7, 1:4)], # put control in front
                              name = "Z-score",
                              col = RWB_color, 
                              rect_gp = gpar(col = "white", lwd = 2), 
                              cluster_columns = F, 
                              show_row_dend = F, 
                              show_column_names = F,
                              top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                                 col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
chemokines_mtx_all <- rbind(
    chemokines_mtx_no[rownames(xx@matrix), c(5:7, 1:4)],
    chemokines_mtx_up[, c(5:7, 1:4)]
)
row_split <- data.frame(group = factor(as.numeric(rownames(chemokines_mtx_all) %in% rownames(chemokines_mtx_up)), labels = c("No", "Up")), row.names = rownames(chemokines_mtx_all))
ComplexHeatmap::Heatmap(chemokines_mtx_all, # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 1), 
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split, gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# cytokines
xx <- ComplexHeatmap::Heatmap(cytokines_mtx_no[, c(5:7, 1:4)], # put control in front
                              name = "Z-score",
                              col = RWB_color, 
                              rect_gp = gpar(col = "white", lwd = 2), 
                              cluster_columns = F, 
                              show_row_dend = F, 
                              show_column_names = F,
                              top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                                 col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
cytokines_mtx_all <- rbind(
    cytokines_mtx_no[rownames(xx@matrix), c(5:7, 1:4)],
    cytokines_mtx_down[, c(5:7, 1:4)],
    cytokines_mtx_up[, c(5:7, 1:4)]
)
cytokines_mtx_all_number <- unlist(lapply(
    rownames(cytokines_mtx_all),
    function(X) {
        if (X %in% rownames(cytokines_mtx_no))
            1
        else if (X %in% rownames(cytokines_mtx_down))
            2
        else if (X %in% rownames(cytokines_mtx_up))
            3
    }
))
row_split <- data.frame(group = factor(cytokines_mtx_all_number, labels = c("No", "Down", 'Up')), row.names = rownames(cytokines_mtx_all))
ComplexHeatmap::Heatmap(cytokines_mtx_all, # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 1), 
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = F, 
                        row_split = row_split, gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
pp_cytokines <- ComplexHeatmap::Heatmap(cytokines_mtx_all, # put control in front
                                        name = "Z-score",
                                        col = RWB_color, 
                                        rect_gp = gpar(col = "white", lwd = 1), 
                                        cluster_columns = F, 
                                        cluster_rows = F,
                                        show_row_dend = F, 
                                        show_column_names = F,
                                        show_row_names = T, 
                                        row_split = row_split, 
                                        gap = unit(0.2, "inches"), 
                                        row_names_gp = gpar(fontsize = 2.5),
                                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")), 
                                                                           col = list(Group = c("Non-CA" = "grey", "CA" = "orange")))
)
write.csv(row_split[order(row_split$group),,drop=F], file = 'E://Cryo-TCR/bulk/cytokines_list_rm3.csv')
pdf("E://Cryo-TCR/bulk/cytokines_heatmap.pdf", width = 3.5, height = 15)
draw(pp_cytokines)
dev.off()
#############################
#############################GSVA for pathway
library(GSVA)
library(GSEABase)
RWB_color <- grDevices::colorRampPalette(colors = c("blue","white","red"))(100)
GSV_path <- "C://DeepTL_dev-main/pathway/GSVA/"
gmt_file <- paste0(GSV_path,"c5.go.v7.4.symbols1.gmt")
geneset <- getGmt(gmt_file)

exp_X_Poisson <- data@assays@data$counts
rownames(exp_X_Poisson) <- toupper(rownames(exp_X_Poisson))
## GSVA analysis by samples
GSVA_res_Poisson <- gsva(exp_X_Poisson, 
                          geneset, 
                          min.sz = 10, 
                          max.sz = 1000, 
                          verbose = TRUE, 
                          method = "gsva", 
                          kcdf = "Poisson")

GSVA_diff <- apply(GSVA_res_Poisson, 1, function(x) wilcox.test(x[1:4], x[5:7], alternative = 'greater')$p.value)
# GSVA_diff_padj <- p.adjust(GSVA_diff, method = 'BH')
GSVA_diff <- sort(GSVA_diff, decreasing = F)
# GSVA_diff_padj <- sort(GSVA_diff_padj, decreasing = F)

pathway_up <- names(GSVA_diff)[GSVA_diff < 0.05]
pathway_interferon <- names(GSVA_diff)[grep("interferon", names(GSVA_diff), ignore.case = T)]
pathway_up_interferon <- pathway_up[grep("interferon", pathway_up, ignore.case = T)]
pathway_chemokine <- names(GSVA_diff)[grep("chemokine", names(GSVA_diff), ignore.case = T)]
pathway_up_chemokine <- pathway_up[grep("chemokine", pathway_up, ignore.case = T)]
pathway_cytokine <- names(GSVA_diff)[grep("cytokine", names(GSVA_diff), ignore.case = T)]
pathway_up_cytokine <- pathway_up[grep("cytokine", pathway_up, ignore.case = T)]
pathway_inflammatory <- names(GSVA_diff)[grep("inflammatory", names(GSVA_diff), ignore.case = T)]
pathway_up_inflammatory <- pathway_up[grep("inflammatory", pathway_up, ignore.case = T)]

GSVA_scale <- t(apply(GSVA_res_Poisson, 1, scale))
colnames(GSVA_scale) <- colnames(GSVA_res_Poisson)


ComplexHeatmap::Heatmap(GSVA_up_scale, # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 1),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = F, 
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
# interferon
interferon_scale_mtx <- GSVA_scale[pathway_interferon, c(5:7,1:4)]
rownames(interferon_scale_mtx) <- tolower(sub("^GO.._", "", rownames(interferon_scale_mtx)))
up_index_tmp <- which(rownames(interferon_scale_mtx) %in% tolower(sub("^GO.._", "", pathway_up_interferon)))
row_split <- data.frame(group = factor(1:nrow(interferon_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(interferon_scale_mtx))
row_split <- row_split[c(up_index_tmp, setdiff(1:nrow(interferon_scale_mtx), up_index_tmp)), , drop = F]
interferon_scale_mtx <- interferon_scale_mtx[c(up_index_tmp, setdiff(1:nrow(interferon_scale_mtx), up_index_tmp)), ]
rownames(interferon_scale_mtx) <- gsub("_", " ", rownames(interferon_scale_mtx))
ComplexHeatmap::Heatmap(interferon_scale_mtx,  # put control in front
                        name = "Z-score", 
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 3),
                        row_names_max_width = unit(20, "inches"),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split,
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# chemokine
chemokine_scale_mtx <- GSVA_scale[pathway_chemokine, c(5:7,1:4)]
rownames(chemokine_scale_mtx) <- tolower(sub("^GO.._", "", rownames(chemokine_scale_mtx)))
up_index_tmp <- which(rownames(chemokine_scale_mtx) %in% tolower(sub("^GO.._", "", pathway_up_chemokine)))
row_split <- data.frame(group = factor(1:nrow(chemokine_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(chemokine_scale_mtx))
row_split <- row_split[c(up_index_tmp, setdiff(1:nrow(chemokine_scale_mtx), up_index_tmp)), , drop = F]
chemokine_scale_mtx <- chemokine_scale_mtx[c(up_index_tmp, setdiff(1:nrow(chemokine_scale_mtx), up_index_tmp)), ]
rownames(chemokine_scale_mtx) <- gsub("_", " ", rownames(chemokine_scale_mtx))
ComplexHeatmap::Heatmap(chemokine_scale_mtx,  # put control in front
                        name = "Z-score", 
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 3),
                        row_names_max_width = unit(20, "inches"),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split,
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# cytokine
cytokine_scale_mtx <- GSVA_scale[pathway_cytokine, c(5:7,1:4)]
rownames(cytokine_scale_mtx) <- tolower(sub("^GO.._", "", rownames(cytokine_scale_mtx)))
up_index_tmp <- which(rownames(cytokine_scale_mtx) %in% tolower(sub("^GO.._", "", pathway_up_cytokine)))
row_split <- data.frame(group = factor(1:nrow(cytokine_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(cytokine_scale_mtx))
row_split <- row_split[c(up_index_tmp, setdiff(1:nrow(cytokine_scale_mtx), up_index_tmp)), , drop = F]
cytokine_scale_mtx <- cytokine_scale_mtx[c(up_index_tmp, setdiff(1:nrow(cytokine_scale_mtx), up_index_tmp)), ]
rownames(cytokine_scale_mtx) <- gsub("_", " ", rownames(cytokine_scale_mtx))
ComplexHeatmap::Heatmap(cytokine_scale_mtx,  # put control in front
                        name = "Z-score", 
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 3),
                        row_names_max_width = unit(20, "inches"),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split,
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
write.csv(annotation_row, file = "E://Cryo-TCR/bulk/cytokines_list.csv")

# inflammatory
inflammatory_scale_mtx <- GSVA_scale[pathway_inflammatory, c(5:7,1:4)]
rownames(inflammatory_scale_mtx) <- tolower(sub("^GO.._", "", rownames(inflammatory_scale_mtx)))
up_index_tmp <- which(rownames(inflammatory_scale_mtx) %in% tolower(sub("^GO.._", "", pathway_up_inflammatory)))
row_split <- data.frame(group = factor(1:nrow(inflammatory_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(inflammatory_scale_mtx))
row_split <- row_split[c(up_index_tmp, setdiff(1:nrow(inflammatory_scale_mtx), up_index_tmp)), , drop = F]
inflammatory_scale_mtx <- inflammatory_scale_mtx[c(up_index_tmp, setdiff(1:nrow(inflammatory_scale_mtx), up_index_tmp)), ]
rownames(inflammatory_scale_mtx) <- gsub("_", " ", rownames(inflammatory_scale_mtx))
ComplexHeatmap::Heatmap(inflammatory_scale_mtx,  # put control in front
                        name = "Z-score", 
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 3),
                        row_names_max_width = unit(20, "inches"),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split,
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)

# custom
pathway_custom0 <- c("chemokine_production", 
                    "T_helper_2_cell_cytokine_production",
                    "cytokine_mediated_signaling_pathway",
                    "cytokine_production",
                    "response_to_type_i_interferon",
                    "response_to_interferon_gamma",
                    "response_to_interferon_alpha",
                    "response_to_interferon_beta",
                    "interferon_alpha_production",
                    "interferon_gamma_production",
                    "cellular_response_to_interferon_alpha",
                    "inflammatory_response")
pathway_custom <- sapply(pathway_custom0, function(x) names(GSVA_diff)[grep(paste0("^GO.._", x, "$"), names(GSVA_diff), ignore.case = T)], USE.NAMES = F)
pathway_up_custom <- unlist(sapply(pathway_custom, function(x) pathway_up[grep(paste0("^", x, "$"), pathway_up, ignore.case = T)], USE.NAMES = F))

custom_scale_mtx <- GSVA_scale[pathway_custom, c(5:7,1:4)]
rownames(custom_scale_mtx) <- tolower(sub("^GO.._", "", rownames(custom_scale_mtx)))
up_index_tmp <- which(rownames(custom_scale_mtx) %in% tolower(sub("^GO.._", "", pathway_up_custom)))
row_split <- data.frame(group = factor(1:nrow(custom_scale_mtx) %in% up_index_tmp, labels = c("No", "Up")), row.names = rownames(custom_scale_mtx))
row_split <- row_split[c(up_index_tmp, setdiff(1:nrow(custom_scale_mtx), up_index_tmp)), , drop = F]
custom_scale_mtx <- custom_scale_mtx[c(up_index_tmp, setdiff(1:nrow(custom_scale_mtx), up_index_tmp)), ]
rownames(custom_scale_mtx) <- gsub("_", " ", rownames(custom_scale_mtx))
ComplexHeatmap::Heatmap(custom_scale_mtx,  # put control in front
                        name = "Z-score",
                        col = RWB_color, 
                        rect_gp = gpar(col = "white", lwd = 3),
                        row_names_max_width = unit(20, "inches"),
                        cluster_columns = F, 
                        cluster_rows = F,
                        show_row_dend = F, 
                        show_column_names = F,
                        show_row_names = T, 
                        row_split = row_split,
                        gap = unit(0.2, "inches"), 
                        top_annotation = HeatmapAnnotation(Group = factor(annotation_col$Group[c(5:7, 1:4)], levels = c("NonCryo", "Cryo")), 
                                                           col = list(Group = c("Cryo" = "orange", "NonCryo" = "grey")))
)
