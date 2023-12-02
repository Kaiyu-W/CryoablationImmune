# GSVA for signature

source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
# load(file = 'Combined_analysis.rda')
load(file = 'Combined_analysis_TNK_real.rda')
Idents(T_sub_real) <- 'T_sub_MainCluster2'
my_plotDim(T_sub_real, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, title = 'T/NK R=0.3', group.by = 'T_sub_MainCluster2')

# setwd("E://Cryo-TCR/server/211111/")
# load(file = "Cryo_T_analysis.rda")

# library(GSVA)
# library(GSEABase)
# RWB_color <- grDevices::colorRampPalette(colors = c("blue","white","red"))(100)
# gmt_file <- "E://Cryo-TCR/server/210924/custom_signatures_for_GSEA.gmt"
# geneset <- getGmt(gmt_file)
# exp_X_gaussian <- T_sub_real@assays$SCT@data
# exp_X_gaussian <- as.matrix(exp_X_gaussian)
# rownames(exp_X_gaussian) <- toupper(rownames(exp_X_gaussian))
# # set_of_genes <- unique(unlist(lapply(geneset@.Data, function(x) x@geneIds)))
# # set_of_genes_overlap <- intersect(set_of_genes, rownames(exp_X_gaussian))
# # exp_X_gaussian <- exp_X_gaussian[set_of_genes_overlap, ]
# 
# ## GSVA analysis by cells
# GSVA_res_gaussian <- gsva(exp_X_gaussian, 
#                           geneset, 
#                           min.sz = 0, 
#                           max.sz = 1000, 
#                           verbose = TRUE, 
#                           method = "gsva", 
#                           kcdf = "Gaussian")
# 
# my_heatmap <- function(SeuratObj, meta, mtx, scale = FALSE) {
#     RWB_color <- grDevices::colorRampPalette(colors = c("blue","white","red"))(100)
#     
#     anno_col <- SeuratObj[[meta]]
#     anno_col <- anno_col[order(anno_col[,1]), , drop = F]
#     anno_col[,1] <- unfactor(anno_col[,1])
#     mtx <- mtx[, rownames(anno_col)]
#     
#     gap_col <- c()
#     for (i in 1:(nrow(anno_col)-1)) {
#         if (anno_col[i, 1] != anno_col[i+1, 1])
#             gap_col <- append(gap_col, i)
#     }
#     
#     if (scale) {
#         mtx1 <- t(apply(mtx, 1, scale))
#         colnames(mtx1) <- colnames(mtx)
#         mtx <- mtx1
#     }
#     
#     pheatmap::pheatmap(mtx, 
#                        fontsize_row = 5,
#                        col = RWB_color,
#                        show_colnames = F,
#                        show_rownames = T,
#                        cluster_rows = F,
#                        cluster_cols = F,
#                        annotation_col = anno_col, 
#                        gaps_col = gap_col
#     )
#     
#     
# }
# 
# my_heatmap(T_sub_real, 'T_sub_MainCluster2', GSVA_res_gaussian, scale = T)
# my_heatmap(T_sub_real, 'orig.ident2', GSVA_res_gaussian, scale = T)
# 
# 
# mtx_tmp <- GSVA_res_gaussian
# 
# df <- data.frame(
#     # signature = unlist(lapply(rownames(mtx_tmp), function(x) rep(x, ncol(mtx_tmp)))),
#     signature = unlist(lapply(c("IFN-gamma", "Effector", "Cytolytic"), function(x) rep(x, ncol(mtx_tmp)))),
#     value = unlist(lapply(rownames(mtx_tmp), function(x) mtx_tmp[x, , drop = T])),
#     bc = rep(colnames(mtx_tmp), nrow(mtx_tmp))
# )
# df$group <- unfactor(T_sub_real$orig.ident2[df$bc])
# df$group[df$group == "Cryo"] <- "CA"
# df$group[df$group == "NonCryo"] <- "Non-CA"
# df$group <- factor(df$group, levels = c("Non-CA", "CA"))
# df$cluster <- T_sub_real$T_sub_MainCluster2[df$bc]
# df$TNK <- unfactor(df$cluster)
# df$TNK[grep("CD8", df$TNK)] <- "T_CD8"
# df$TNK[grep("Th", df$TNK)] <- "T_CD4"
# df$TNK[grep("Treg", df$TNK)] <- "T_CD4"
# df$TNK[grep("Tfh", df$TNK)] <- "T_CD4"
# 
# library(ggpubr)
# 
# # plot case-control of all T/NK cells
# p <- ggboxplot(df, x = 'signature', y = 'value', color = 'group', palette = 'jco', add = 'boxplot')
# p + stat_compare_means(aes(group = group), method = "t.test", method.args = list(alternative = "less")) + RotatedAxis()
# p + stat_compare_means(aes(group = group), method = "wilcox.test", method.args = list(alternative = "less")) + RotatedAxis()
# 
# # plot case-control among each T clusters
# clusters <- levels(df$cluster)
# plot_box <- function(df, cluster, test_method = "t.test", significant_format = c('p.format', 'p.signif'), angle = 90) {
#     match.arg(significant_format)
#     significant_format <- significant_format[1]
#     
#     df <- df[df$cluster == cluster, ]
#     p <- ggboxplot(df, x = 'signature', y = 'value', color = 'group', palette = 'jco', add = 'boxplot')
#     p + stat_compare_means(aes(group = group), method = test_method, method.args = list(alternative = "less"), label = significant_format) + 
#         theme(axis.text.x = element_text(angle = angle, hjust = 1)) +
#         ggtitle(cluster) + 
#         theme(plot.title = element_text(hjust = 0.5))
# }
# 
# # #9
# plot_box(df, clusters[1])
# plot_box(df, clusters[2])
# plot_box(df, clusters[3])
# plot_box(df, clusters[4])
# plot_box(df, clusters[5])
# plot_box(df, clusters[6])
# plot_box(df, clusters[7])
# plot_box(df, clusters[8])
# plot_box(df, clusters[9])
# 
# lapply(clusters, 
#        function(x) 
#            ggsave(plot_box(df, x), 
#                   file = paste0("custom_signature_", x, ".png")
#                   )
#        )
# 
# save(GSVA_res_gaussian, file = "custom_signature_GSVA_res_gaussian.rda")
# 
# # plot case-control of T_CD4/CD8
# clusters <- unique(df$TNK)
# plot_box <- function(df, cluster, test_method = "t.test", significant_format = c('p.format', 'p.signif'), angle = 90) {
#     match.arg(significant_format)
#     significant_format <- significant_format[1]
#     
#     df <- df[df$TNK == cluster, ]
#     p <- ggboxplot(df, x = 'signature', y = 'value', color = 'group', palette = 'jco', add = 'boxplot')
#     p + stat_compare_means(aes(group = group), method = test_method, method.args = list(alternative = "less"), label = significant_format) +
#         theme(axis.text.x = element_text(angle = angle, hjust = 1)) +
#         ggtitle(cluster) +
#         theme(plot.title = element_text(hjust = 0.5))
# }
# 
# plot_box(df, clusters[1])


############### by average
RWB_color <- grDevices::colorRampPalette(colors = c("blue","white","red"))(100)
gmt_file <- "E://Cryo-TCR/server/210924/custom_signatures_for_GSEA.gmt"
geneset <- getGmt(gmt_file)
exp_X_gaussian <- T_sub_real@assays$SCT@data
# exp_X_gaussian <- T_sub_real@assays$SCT@scale.data
# exp_X_gaussian <- as.matrix(exp_X_gaussian)
rownames(exp_X_gaussian) <- toupper(rownames(exp_X_gaussian))
set_of_genes <- unique(unlist(lapply(geneset@.Data, function(x) x@geneIds)))
set_of_genes_overlap <- intersect(set_of_genes, rownames(exp_X_gaussian))
exp_X_gaussian <- exp_X_gaussian[set_of_genes_overlap, ]

list_of_genes <- lapply(geneset@.Data, function(x) x@geneIds)

res_list <- list()
for (i in seq(list_of_genes)) {
    tmp <- apply(exp_X_gaussian, 2, function(x) {
        tapply(x, rownames(exp_X_gaussian) %in% list_of_genes[[i]], mean)
    })
    
    res_list[[i]] <- tmp['TRUE',]-tmp['FALSE',]
}
names(res_list) <- names(geneset)
res_mtx <- as.matrix(data.frame(res_list))
head(res_mtx)

# colnames(res_mtx)[grep("Exhausted", colnames(res_mtx), ignore.case = T)]


res_mtx <- apply(res_mtx, 2, function(x) (x-min(x))/(max(x)-min(x)))

mtx_tmp <- t(res_mtx)[1:3,]
mtx_tmp

df <- data.frame(
    # signature = unlist(lapply(rownames(mtx_tmp), function(x) rep(x, ncol(mtx_tmp)))),
    signature = unlist(lapply(c("IFN-gamma", "Effector", "Cytolytic"), function(x) rep(x, ncol(mtx_tmp)))),
    scale_value = unlist(lapply(rownames(mtx_tmp), function(x) mtx_tmp[x, , drop = T])),
    bc = rep(colnames(mtx_tmp), nrow(mtx_tmp))
)
df$group <- unfactor(T_sub_real$orig.ident2[df$bc])
df$group[df$group == "Cryo"] <- "CA"
df$group[df$group == "NonCryo"] <- "Non-CA"
df$group <- factor(df$group, levels = c("Non-CA", "CA"))
df$cluster <- T_sub_real$T_sub_MainCluster2[df$bc]
df$TNK <- unfactor(df$cluster)
df$TNK[grep("CD8", df$TNK)] <- "T_CD8"
df$TNK[grep("Th", df$TNK)] <- "T_CD4"
df$TNK[grep("Treg", df$TNK)] <- "T_CD4"
df$TNK[grep("Tfh", df$TNK)] <- "T_CD4"

library(ggpubr)

# plot case-control of all T/NK cells
p <- ggboxplot(df, x = 'signature', y = 'scale_value', color = 'group', palette = 'jco', add = 'violinplot')
p <- ggviolin(df, x = 'signature', y = 'scale_value', color = 'group', palette = 'jco', add = 'boxplot', add.params = list(size = 1))
p + stat_compare_means(aes(group = group), method = "t.test", method.args = list(alternative = "less")) + RotatedAxis()
# p + stat_compare_means(aes(group = group), method = "wilcox.test", method.args = list(alternative = "less")) + RotatedAxis()

# plot case-control among each T clusters
clusters <- levels(df$cluster)
plot_box <- function(df, cluster, test_method = "t.test", significant_format = c('p.format', 'p.signif'), angle = 90) {
    match.arg(significant_format)
    significant_format <- significant_format[1]
    
    df <- df[df$cluster == cluster, ]
    p <- ggboxplot(df, x = 'signature', y = 'scale_value', color = 'group', palette = 'jco', add = 'boxplot')
    p + stat_compare_means(aes(group = group), method = test_method, method.args = list(alternative = "less"), label = significant_format) +
        theme(axis.text.x = element_text(angle = angle, hjust = 1)) +
        ggtitle(cluster) + 
        theme(plot.title = element_text(hjust = 0.5))
}

# #9
plot_box(df, clusters[1])
plot_box(df, clusters[2])
plot_box(df, clusters[3])
plot_box(df, clusters[4])
plot_box(df, clusters[5])
plot_box(df, clusters[6])
plot_box(df, clusters[7])
plot_box(df, clusters[8])
plot_box(df, clusters[9])

lapply(clusters,
       function(x)
           ggsave(plot_box(df, x, significant_format = 'p.signif', angle = 90),
                  file = paste0("custom_signature_", x, ".png"),
                  height = 15, width = 10
                  )
       )

# plot case-control of T_CD4/CD8
clusters <- unique(df$TNK)
plot_box <- function(df, cluster, test_method = "t.test", significant_format = c('p.format', 'p.signif'), angle = 90) {
    match.arg(significant_format)
    significant_format <- significant_format[1]
    
    df <- df[df$TNK == cluster, ]
    p <- ggboxplot(df, x = 'signature', y = 'scale_value', color = 'group', palette = 'jco', add = 'boxplot')
    p + stat_compare_means(aes(group = group), method = test_method, method.args = list(alternative = "less"), label = significant_format) +
        theme(axis.text.x = element_text(angle = angle, hjust = 1)) +
        ggtitle(cluster) + 
        theme(plot.title = element_text(hjust = 0.5))
}

p1 <- plot_box(df, clusters[1], significant_format = 'p.signif', angle = 45)

plot_violin <- function(df, cluster, test_method = "t.test", significant_format = c('p.format', 'p.signif'), angle = 90) {
    match.arg(significant_format)
    significant_format <- significant_format[1]

    df <- df[df$TNK == cluster, ]
    p <- ggviolin(df, x = 'signature', y = 'scale_value', color = 'group', palette = 'jco', add = 'boxplot')
    p + stat_compare_means(aes(group = group), method = test_method, method.args = list(alternative = "less"), label = significant_format) +
        theme(axis.text.x = element_text(angle = angle, hjust = 1)) +
        ggtitle(cluster) +
        theme(plot.title = element_text(hjust = 0.5))
}

p2 <- plot_violin(df, clusters[1], significant_format = 'p.signif', angle = 45)

ggsave(p1,
       file = paste0("custom_signature_T_CD8_box.png")
)
ggsave(p2,
       file = paste0("custom_signature_T_CD8_violin.png")
)
