source("E://Cryo-TCR/server/auto/utilities.r")
setwd("E://Cryo-TCR/server/210924/")
load("T_sub_real_meta_220630.rda")
load(file = 'Combined_analysis.rda')
load(file = 'Combined_analysis_Myeloid_real.rda')

# my_plotDim(Cryo_merge, group.by = 'MainCluster')
# my_plotDim(Myeloid_sub_real, group.by = 'MainCluster2_new')
# my_plotDim(T_sub_real, group.by = 'CD8_sub_MainCluster_paper')

Idents(Cryo_merge) <- 'MainCluster'
Cryo_merge <- my_AddMeta(Cryo_merge, new_ident = Myeloid_sub_real$MainCluster2_new, Replace = T)
Idents(Cryo_merge) <- 'Myeloid_sub_real_MainCluster2_new'
Cryo_merge <- my_AddMeta(Cryo_merge, new_ident = T_meta_220630[,'CD8_sub_MainCluster_paper',drop=F], Replace = T)
Cryo_merge$MainCluster_sub <- Cryo_merge$CD8_sub_MainCluster_paper
# my_plotDim(Cryo_merge, group.by = 'MainCluster_sub')
Idents(Cryo_merge) <- 'MainCluster_sub'


uniprot_tb <- readxl::read_excel("uniprot_cell_membrane.xlsx")
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
genes <- bitr(geneID = uniprot_tb$Entry, fromType = 'UNIPROT', toType = 'SYMBOL', OrgDb = org.Hs.eg.db)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
genes_mouse <- getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = genes$SYMBOL, 
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows=T
)
genes_membrane <- genes_mouse[,2]


length(genes_membrane)
nrow(Cryo_merge@assays$RNA@data)
gene_over <- intersect(genes_membrane, rownames(Cryo_merge@assays$RNA@data))
length(gene_over)


DefaultAssay(Cryo_merge) <- 'RNA'
Cryo_merge <- NormalizeData(Cryo_merge, normalization.method = "LogNormalize")
Cryo_merge <- ScaleData(Cryo_merge, features = rownames(Cryo_merge))

avg_counts <- my_avg(Cryo_merge, c('RNA','SCT'), gene_over, group.by = 'MainCluster_sub', slot = 'counts')
avg_data <- my_avg(Cryo_merge, c('RNA','SCT'), gene_over, group.by = 'MainCluster_sub', slot = 'data')
colnames_custom <- c(
    "B",
    "Monocyte_S100a8/9+","Monocyte","Macrophage_M1","Macrophage_M2","cDC_Itgax-","cDC_Itgax+","pDC","Mast",
    "T","T_others","T_Naive","Tfh","Th1","Th17","Treg",
    "I-IFN","CD8_Tem","CD8_Tex","CD8_Tex_Proliferative",
    "NK","NK_Proliferative"
)
avg_counts <- lapply(avg_counts, function(x)x[,colnames_custom])
avg_data <- lapply(avg_data, function(x)x[,colnames_custom])

write.csv(avg_counts$RNA, file = 'membrane_genes_expression_rawcount1.csv')
write.csv(avg_counts$SCT, file = 'membrane_genes_expression_rawcount2.csv')
write.csv(avg_data$RNA, file = 'membrane_genes_expression_norm1.csv')
write.csv(avg_data$SCT, file = 'membrane_genes_expression_norm2.csv')


# my_Heatmap <- function(
#     seurat_obj, 
#     group.by, 
#     genes,
#     slot = c('data', 'scale.data', 'counts'), 
#     default.assay = c("active.assay", names(seurat_obj@assays)), 
#     cell = NULL, 
#     cluster_rows = FALSE, 
#     cluster_cols = FALSE, 
#     show_colnames = FALSE,
#     show_rownames = FALSE,
#     color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
#     use_pheatmap = FALSE,
#     ...
# ) {
#     match.arg(default.assay)
#     match.arg(slot)
#     default.assay <- default.assay[1]
#     slot <- slot[1]
#     if (default.assay == "active.assay")
#         default.assay <- DefaultAssay(seurat_obj)
#     
#     data <- methods::slot(seurat_obj@assays[[default.assay]], slot)
#     
#     anno_df <- seurat_obj@meta.data[, group.by, drop = F]
#     anno_df<- anno_df[order(anno_df[, 1, drop = T]), , drop = F]
#     genes_in <- genes[genes %in% rownames(data)]
#     mtx <- data[genes_in, rownames(anno_df)]
#     
#     fun <- if (use_pheatmap)
#         pheatmap::pheatmap
#     else
#         ComplexHeatmap::pheatmap
#     fun(
#         mat = as.matrix(mtx), 
#         color = color,
#         cluster_cols = cluster_cols, 
#         cluster_rows = cluster_rows, 
#         annotation_col = anno_df, 
#         show_colnames = show_colnames,
#         show_rownames = show_rownames,
#         ...
#     )
# }

# pdf('tmp2.pdf', width = 25, height = 25)
# my_Heatmap(
#     seurat_obj = Cryo_merge, 
#     group.by = 'MainCluster_sub', 
#     slot = 'data',
#     genes = gene_over,
#     default.assay = 'SCT',
#     cluster_rows = T,
#     cluster_cols = T,
#     show_colnames = F,
#     show_rownames = F,
#     use_pheatmap = F
# )
# dev.off()

saveHeat <- function(df, file_pdf = 'tmp.pdf', width = 10, height = 250, ...) {
    pdf(file_pdf, width = width, height = height)
    df <- df[apply(df, 1, function(x)!all(x==0)),]
    # df <- scale(df, center = F, scale = T)
    draw(
        ComplexHeatmap::Heatmap(df, cluster_columns = F, ...)
    )
    dev.off()
}

saveHeat(avg_counts$RNA, 'membrane_genes_heatmap_rawcount1.pdf')
saveHeat(avg_counts$SCT, 'membrane_genes_heatmap_rawcount2.pdf')
saveHeat(avg_data$RNA, 'membrane_genes_heatmap_norm1.pdf')
saveHeat(avg_data$SCT, 'membrane_genes_heatmap_norm2.pdf')
