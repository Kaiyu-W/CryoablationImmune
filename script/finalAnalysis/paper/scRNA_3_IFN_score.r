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

source("~/Desktop/Github/My_scRNA_pipeline/Pathway_score.r")
gmt_file <- "/mnt/c/DeepTL_dev-main/pathway/GSVA/c5.go.v7.4.symbols1_ifn.gmt"
geneset <- getGmt(gmt_file)
pathway_names <- names(geneset)[grepl('RESPONSE',names(geneset))&!grepl('GAMMA',names(geneset))] 
geneset_IIFN <- geneset[pathway_names]

set_of_genes <- unique(unlist(lapply(geneset_IIFN@.Data, function(x) x@geneIds)))
set_of_genes_overlap <- Hmisc::capitalize(tolower(
    intersect(set_of_genes, toupper(rownames(T_sub@assays$SCT@data)))
))
set_of_genes_overlap_filter <- c(
    "Adar","Irf9","Ifitm3","Ifitm2","Usp18",
    "Gbp2","Ifi35","Ifit2","Ifit1","Ifit3",
    "Irf7","Isg20","Myd88","Oas2","Oas3","Ythdf2","Xaf1","Shfl",
    "Sp100","Stat1","Stat2","Bst2",
    "Zbp1","Trim56","Nlrc5","Ifitm1","Rsad2","Isg15","Eif2ak2"
)

T_sub_tmp1 <- geneset_score(T_sub, geneset_IIFN, method = "average", slot = "data", highly_variable = T)
FeaturePlot(T_sub_tmp1, features = pathway_names, pt.size = 1.5, cols = c('white', 'red'))

# dotplot
p1 <- my_DotPlot_split(T_sub_tmp1, features = pathway_names, group.by = 'CellType') + 
    theme(
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
    )
ggsave(
    filename = "plot/scRNA3_I_IFN_dotplot.pdf", 
    plot = p1,
    width = 7, height = 8
    )
# my_violin(T_sub_tmp1, features = set_of_genes_overlap, mode = 'mtx', group.by = 'CellType')


# # heatmap
# color_list <- list(
#     'CellType' = c(
#         "CD4_Naive" = "#709ae1",
#         "CD4_T_I-IFN" = "#197ec0",
#         "Th1" = "#d2af81",
#         "Th17" = "#bd559f",
#         "Treg" = "#f7c530",
#         "Tfh" = "#46732e",
#         "CD8_Naive" = "#d5e4a2",
#         "CD8_Tem_I-IFN" = "#ff410d",
#         "CD8_Tem" = "#95cc5e",
#         "CD8_Tex" = "#d0dfe6",
#         "CD8_Tex_Proliferative" = "#f79d1e",
#         "T_others" = "#748aa6"
#     )
# )
# T_sub_tmp1$CellType <- factor(
#     T_sub_tmp1$IFN_Final,
#     levels = names(color_list$CellType)
# )
# p2 <- ComplexHeatmap::pheatmap(
#     apply(
#         T_sub_tmp1@meta.data[names(sort(T_sub_tmp1$CellType)), pathway_names], 
#         2,
#         function(x) x / max(x)
#     ),
#     annotation_row = T_sub_tmp1@meta.data[names(sort(T_sub_tmp1$CellType)), 'CellType', drop = F],
#     cluster_rows = F,
#     cluster_cols = F,
#     show_rownames = F, 
#     show_colnames = T,
#     angle_col = '90',
#     col = grDevices::colorRampPalette(colors = c("blue", "white", "red"))(100),
#     annotation_colors = color_list,
#     fontsize_col = 5,
#     heatmap_legend_param = list(title = 'activity')
# )
# pdf("plot/scRNA3_I_IFN_heatmap.pdf", width = 5, height = 10)
# draw(p2)
# dev.off()

# GO:0034340 GO:0035455 GO:0035456 GO:0035457 GO:0035458
# GOBP_RESPONSE_TO_TYPE_I_INTERFERON GOBP_RESPONSE_TO_INTERFERON_ALPHA GOBP_RESPONSE_TO_INTERFERON_BETA GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA

# matrixplot
mtx <- apply(
    sapply(
        pathway_names,
        function(x) {
            tapply(T_sub_tmp1@meta.data[[x]], T_sub_tmp1$CellType, mean)
        }
    )[names(color_list$CellType), ],
    2,
    scale
)
rownames(mtx) <- names(color_list$CellType)
colnames(mtx) <- c("GO:0034340","GO:0035455","GO:0035456","GO:0035457","GO:0035458")
p3 <- ComplexHeatmap::Heatmap(
    t(mtx),
    cluster_rows = F,
    cluster_columns = F,
    show_row_names = T, 
    show_column_names = F,
    col = grDevices::colorRampPalette(colors = c("#5badc0", "white", "#c53c3d"))(100),
    top_annotation = HeatmapAnnotation(
        CellType = factor(names(color_list$CellType), levels = names(color_list$CellType)),
        col = color_list,
        show_legend = T
    ),
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    name = 'z-score'
)
pdf("plot/scRNA3_I_IFN_matrixplot.pdf", width = 8, height = 4)
draw(p3, merge_legends = T)
dev.off()


# genes
# matrixplot
p11 <- my_DotPlot_split(T_sub_tmp1, features = set_of_genes_overlap_filter, group.by = 'CellType') + 
    theme(
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
    )
mat_df <- reshape::cast(p11$data, features.plot ~ id, value = 'avg.exp')
rownames(mat_df) <- mat_df$features.plot
dn <- dimnames(mat_df)
mat_df <- as.matrix(mat_df[,2:ncol(mat_df)])
rownames(mat_df) <- dn[[1]]
colnames(mat_df) <- dn[[2]][2:length(dn[[2]])]
mat_df

# matrixplot
mtx1 <- t(apply(
    mat_df,
    1,
    scale
))
colnames(mtx1) <- colnames(mat_df)
p5 <- ComplexHeatmap::Heatmap(
    mtx1,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = T, 
    show_column_names = F,
    col = grDevices::colorRampPalette(colors = c("#5badc0", "white", "#c53c3d"))(100),
    top_annotation = HeatmapAnnotation(
        CellType = factor(names(color_list$CellType), levels = names(color_list$CellType)),
        col = color_list,
        show_legend = T
    ),
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    name = 'z-score'
)
pdf("plot/scRNA3_I_IFN_gene_matrixplot.pdf", width = 8, height = 8)
draw(p5, merge_legends = T)
dev.off()



# paper gene
genes_paper <- unique(c(
    'IFI44L','HCK','NMI','IFI16','IFIT1','EIF2AK2',
    'C3AR1','OAS1','IFI30','TLR2','MX1','PLSCR1',
    'IFITM1','PIK3R2','TRIM22','S100A8','CYBB','LYN',
    'ISG20','IFI6','XAF1','UBE2L6','STAT1','MS4A4A',
    'IFI44','ISG15','MX2','PSMB9','IFIT3','IFI27',
    'IFITM2','EGR1','IRF7','JUN','MAX','NFIL3',
    'PLSCR1','SP100','SP110','STAT1','STAT2','TCF7L2',
    'TFEC','ZNFX1'    
))
# library("biomaRt")
# # library(httr)
# # httr::set_config(config(ssl_verifypeer = 0L))
# human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://usaeast.ensembl.org")
# mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://usaeast.ensembl.org")
# genes_paper2 <- getLDS(
#     attributes = c("hgnc_symbol"), 
#     filters = "hgnc_symbol", 
#     values = genes_paper, 
#     mart = human, 
#     attributesL = c("mgi_symbol"), 
#     martL = mouse, 
#     uniqueRows = T
# )
# genes_paper_list <- genes_paper2$MGI.symbol
genes_paper_list <- c(
    "Mx2","Mx1","Lyn","Nfil3","Hck","Xaf1","Ifit3","Ifit3b",
    "Isg20","Pik3r2","Ifi30","Psmb9","Cybb","Egr1","Ube2l6","Ms4a4a",
    "Ms4a4b","Ms4a4c","Ms4a4d","Ifitm1","Ifitm7","Ifitm3","Ifitm2","Plscr2",
    "Plscr1","1700057G04Rik","Tcf7l2","Max","Oas1g","Oas1h","Oas1e","Oas1c",
    "Oas1f","Oas1b","Oas1a","Oas1d","Eif2ak2","Nmi","Ifitm1","Ifitm7",
    "Ifitm3","Ifitm2","Ifi44l","Isg15","Tlr2","Ifi44","C3ar1","Stat2",
    "Irf7","Trim75","Sp110","Stat1","Ifi207","Ifi209","Ifi203","Ifi213",
    "Ifi204","Ifi206","Ifi202b","Ifi214","Ifi208","S100a8","Jun"
)
set_of_genes_overlap2 <- intersect(genes_paper_list, rownames(T_sub@assays$SCT@data))

p12 <- my_DotPlot_split(T_sub_tmp1, features = set_of_genes_overlap2, group.by = 'CellType') + 
    theme(
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
    )
mat_df <- reshape::cast(p12$data, features.plot ~ id, value = 'avg.exp')
rownames(mat_df) <- mat_df$features.plot
dn <- dimnames(mat_df)
mat_df <- as.matrix(mat_df[,2:ncol(mat_df)])
rownames(mat_df) <- dn[[1]]
colnames(mat_df) <- dn[[2]][2:length(dn[[2]])]
# mat_df

# matrixplot
mtx2 <- t(apply(
    mat_df,
    1,
    scale
))
colnames(mtx2) <- colnames(mat_df)
p6 <- ComplexHeatmap::Heatmap(
    mtx2,
    cluster_rows = T,
    cluster_columns = F,
    show_row_names = T, 
    show_column_names = F,
    col = grDevices::colorRampPalette(colors = c("#5badc0", "white", "#c53c3d"))(100),
    top_annotation = HeatmapAnnotation(
        CellType = factor(names(color_list$CellType), levels = names(color_list$CellType)),
        col = color_list,
        show_legend = T
    ),
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    name = 'z-score'
)
pdf("plot/scRNA3_I_IFN_gene_paper_matrixplot.pdf", width = 8, height = 8)
draw(p6, merge_legends = T)
dev.off()
