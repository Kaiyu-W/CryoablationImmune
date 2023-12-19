source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis.rda")

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

# 1. all CD45+ cells:
Cryo_merge

if (!file.exists("scRNA3_cellbarcode.csv")) {
    write.csv(colnames(Cryo_merge), file = "scRNA3_cellbarcode.csv")
}
if (!file.exists("scRNA3_umap.csv")) {
    write.csv(Cryo_merge@reductions$umap@cell.embeddings, file = "scRNA3_umap.csv")
}
if (!file.exists("scRNA3_meta.csv")) {
    write.csv(Cryo_merge@meta.data, file = "scRNA3_meta.csv")
}
# save gene expression
if (!dir.exists("plot/scRNA3_CD45plus_matrix") & !file.exists("plot/scRNA3_CD45plus_matrix.zip")) {
    dir.create("plot/scRNA3_CD45plus_matrix")
    # mtx
    Matrix::writeMM(Cryo_merge@assays$RNA@counts, "plot/scRNA3_CD45plus_matrix/matrix.mtx")
    system("gzip plot/scRNA3_CD45plus_matrix/matrix.mtx")
    # scales::number_bytes((file.size("plot/scRNA3_CD45plus_matrix/matrix.mtx.gz")))
    # barcode
    write.table(
        colnames(Cryo_merge@assays$RNA@counts), 
        file = "plot/scRNA3_CD45plus_matrix/barcodes.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F
    )
    system("gzip -f plot/scRNA3_CD45plus_matrix/barcodes.tsv")
    # feature
    feature_list = rownames(Cryo_merge@assays$RNA@counts)
    feature_df = read.table("/mnt/e/Cryo-TCR/server/210924/Cryo/filtered_feature_bc_matrix/features.tsv.gz", sep = '\t')
    feature_df$V2 <- make.unique(feature_df$V2)
    feature_df <- feature_df[match(feature_list, feature_df$V2), ]
    write.table(
        feature_df, 
        file = "plot/scRNA3_CD45plus_matrix/features.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F,
        sep = '\t'
    )
    system("gzip -f plot/scRNA3_CD45plus_matrix/features.tsv")
}

# 1.1 CD45+ umap
Cryo_merge$MainCluster <- factor(unfactor(Cryo_merge$MainCluster), levels = c("T", "B", "NK", "Myeloid", "Mast"))
CD45plus_plot <- my_plotDim(
    Cryo_merge,
    reduction = "umap", label = T,
    pt.size = 1.5, label.size = 6, repel = T,
    group.by = "MainCluster",
    title = "CD45+ MainCluster"
) +
    theme(
        legend.position = c(0.8, 0.2),
        # legend.key.size = unit(0.5, "inch"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "T" = "#c0abcb",
            "B" = "#4f8ab2",
            "NK" = "#de9393",
            "Myeloid" = "#a1c4d6",
            "Mast" = "#d86762"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
ggsave(
    filename = "plot/scRNA3_CD45plus_dim.pdf",
    plot = CD45plus_plot,
    height = 6, width = 7.69
)

my_CountCluster(Cryo_merge, group2 = "orig.ident2", group1 = "active.ident")

# 1.2 CD45+ markers express
CD45_Marker_rough <- list(
    T = c("Cd3e", "Cd247"),
    NK = c("Klrk1", "Klrb1c"),
    B = c("Cd19", "Cd79a"),
    Myeloid = c("Itgam", "Cd14"),
    Mast = c("Mcpt1", "Mcpt2")
)

a <- list()
markers_list <- c("Cd3e", "Klrk1", "Cd19", "Itgam", "Mcpt1", "Cd247", "Klrb1c", "Cd79a", "Cd14", "Mcpt2")
for (i in markers_list) {
    a[[i]] <- FeaturePlot(
        Cryo_merge,
        features = i, reduction = "umap",
        pt.size = 0.5, cols = c("#e0e7c8", "#cc340a")
    ) +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            text = element_blank()
        ) +
        annotate(
            "text",
            x = 9.5,
            y = -10,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
CD45plus_express <- patchwork::wrap_plots(plots = a, ncol = 5)
ggsave(
    filename = "plot/scRNA3_CD45plus_markers_umap.pdf",
    plot = CD45plus_express,
    height = 6, width = 7.69 * 1.5
)


# 1.3 CD45+ markers violin
Cryo_merge$MainCluster <- factor(unfactor(Cryo_merge$MainCluster), levels = rev(c("T", "B", "NK", "Myeloid", "Mast")))
CD45plus_violin <- VlnPlot(
    Cryo_merge,
    feature = unlist(CD45_Marker_rough),
    group.by = "MainCluster",
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    cols = c(
        "T" = "#c0abcb",
        "B" = "#4f8ab2",
        "NK" = "#de9393",
        "Myeloid" = "#a1c4d6",
        "Mast" = "#d86762"
    )
) +
    # ggsci::scale_fill_jco() +
    theme(
        aspect.ratio = 5,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold"),
        axis.line.x = element_blank(),
        strip.text.x = element_text(angle = 45, size = 15, face = "bold"),
        # panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none"
    )
ggsave(
    filename = "plot/scRNA3_CD45plus_markers_violin.pdf",
    plot = CD45plus_violin,
    height = 6, width = 7.69
)

# 1.4 CD45+ celltype count
colors_bar <- c(
    "T" = "#c0abcb",
    "B" = "#4f8ab2",
    "NK" = "#de9393",
    "Myeloid" = "#a1c4d6",
    "Mast" = "#d86762"
)
CD45plus_count_sum <- as.data.frame(my_CountCluster(Cryo_merge, group1 = "MainCluster", group2 = "orig.ident2")[-6, 1:2])
CD45plus_count_sum2 <- my_CountCluster(Cryo_merge, group1 = "MainCluster", group2 = "orig.ident2")["sum", 1:2]
CD45plus_count_sum[, 1] <- CD45plus_count_sum[, 1] / CD45plus_count_sum2[1]
CD45plus_count_sum[, 2] <- CD45plus_count_sum[, 2] / CD45plus_count_sum2[2]
CD45plus_count_sum[, 1] <- CD45plus_count_sum[, 1] / sum(CD45plus_count_sum[, 1])
CD45plus_count_sum[, 2] <- CD45plus_count_sum[, 2] / sum(CD45plus_count_sum[, 2])
CD45plus_count_sum
CD45plus_count_sum$cell <- factor(rownames(CD45plus_count_sum), levels = names(colors_bar))
CD45plus_count_sum <- reshape2::melt(CD45plus_count_sum, id.vars = "cell")
CD45plus_count_sum$variable <- sub("_count$", "", CD45plus_count_sum$variable)
CD45plus_count_sum$variable <- ifelse(CD45plus_count_sum$variable == "Cryo", "CA", "Non-CA")
CD45plus_count_sum$variable <- factor(CD45plus_count_sum$variable, levels = c("Non-CA", "CA"))
CD45plus_count_Plot <- ggplot(CD45plus_count_sum, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.8, stat = "identity") +
    ylab("Percentage") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "grey", size = 2, lineend = "square"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "grey", size = 2)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = colors_bar
    )
ggsave(
    filename = "plot/scRNA3_CD45plus_count.pdf",
    plot = CD45plus_count_Plot,
    height = 6, width = 3.4
)

# 1.4 to find where I-IFN signals come from
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
########### 
# by GSVA, with all cell types
library(GSVA)
library(GSEABase)
source("~/Desktop/Github/My_scRNA_pipeline/Pathway_score.r")
red <- "#ff5656" # "red"
blue <- "#5858ff" # "blue"
RWB_color <- grDevices::colorRampPalette(colors = c(blue, "white", red))(100)
GSV_path <- "/mnt/c/DeepTL_dev-main/pathway/GSVA/"
gmt_file <- paste0(GSV_path, "c5.go.v7.4.symbols1.gmt")
geneset <- getGmt(gmt_file)
pathway_names <- names(geneset)[grepl("interferon.*production", names(geneset), ignore.case = TRUE)]
geneset_IIFN <- geneset[pathway_names]

Cryo_merge_tmp1 <- Cryo_merge
Cryo_merge_tmp1 <- geneset_score(
    Cryo_merge_tmp1, geneset_IIFN, 
    slot = "data", highly_variable = T,
    method = 'GSVA', 
    gsva_method = 'ssgsea', # "gsva", "ssgsea", "zscore", "plage"
    gsva_kcdf = 'Gaussian'
)
pathway_names <- pathway_names[!grepl("gamma", pathway_names, ignore.case = TRUE)]
pathway_names <- pathway_names[!grepl("regulation", pathway_names, ignore.case = TRUE)]
p_dot <- my_DotPlot_split(
    Cryo_merge_tmp1, 
    features = pathway_names, 
    group.by = 'MainCluster', 
    cols = 'RdBu',
    scale.min = 0
) + coord_flip() + RotatedAxis() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(), 
    # legend.position = 'bottom'
) + guides(size = FALSE)
ggsave(
    filename = "plot/scRNA3_CD45plus_I-IFN_prod_dot.pdf",
    plot = p_dot,
    height = 2.5, width = 10
)

# CA_vs_Non-CA within each cell type
Idents(Cryo_merge) <- "MainCluster"
obj_list <- list()
for (g in unique(Cryo_merge$MainCluster)) {
    sub_obj <- subset(Cryo_merge, ident = g)
    Idents(sub_obj) <- "orig.ident2"
    Markers_sub_list <- FindAllMarkers(sub_obj, only.pos = T)
    Markers_sub_df <- my_Markers2df_multiple(
        Markers_sub_list,
        logFC_threshold = 0.25,
        positive = T,
        n_top = 120
    )
    GO_ <- my_GO(
        Markers_sub_df$Cluster_Cryo,
        return_plot = F, return_res = T,
        ont = "BP", Simplify = F, type = "bar", font.size = 18, show = 30, 
        pvalue = 1, qvalue = 1
    )
    obj_list[[g]] <- GO_@result[grepl('interferon', GO_@result$Description),]
}

ifn_prod_pathways <- sort(unique(unlist(lapply(
    obj_list,
    function(x) {
        y <- x$Description
        y[grepl('production',y)&!grepl('gamma',y)&!grepl('regulation',y)]
    }
))))
ifn_prod_matrix_p <- matrix(
    1, 
    nrow = length(ifn_prod_pathways), 
    ncol = length(obj_list), 
    dimnames = list(
        ifn_prod_pathways, 
        names(obj_list)
    )
)
ifn_prod_matrix_q <- matrix(
    1, 
    nrow = length(ifn_prod_pathways), 
    ncol = length(obj_list), 
    dimnames = list(
        ifn_prod_pathways, 
        names(obj_list)
    )
)
for (i in names(obj_list)) {
    res_ <- obj_list[[i]]
    for (j in ifn_prod_pathways) {
        if (j %in% res_$Description) {
            padj <- res_[res_$Description == j, 'p.adjust']
            q_v <- res_[res_$Description == j, 'qvalue']
            ifn_prod_matrix_p[j,i] <- -log10(padj)
            ifn_prod_matrix_q[j,i] <- q_v
        }
    }
}

p_mat <- ComplexHeatmap::Heatmap(
    mat = ifn_prod_matrix_p,
    name = '-Log10(p.adj)\n*>1.3',
    show_heatmap_legend = TRUE,
    cell_fun = function(j, i, x, y, w, h, col) {
        grid.text(
            ifelse(
                ifn_prod_matrix_q <= 0.2,
                ifelse(
                    ifn_prod_matrix_p >= -log10(0.001),
                    '***',
                    ifelse(
                        ifn_prod_matrix_p >= -log10(0.01),
                        '**',
                        ifelse(
                            ifn_prod_matrix_p >= -log10(0.05),
                            '*',
                            ''
                        )
                    )
                ),
                ''
            )[i, j], 
            x, y
        )
    },
    col = colorRampPalette(rev(brewer.pal(n = 7, name ="PRGn")))(100),
    rect_gp = gpar(col = "white", lwd = 1),
    gap = unit(0.05, "inches"),
    column_gap = unit(0.05, "inches"),
    row_names_gp = gpar(fontsize = 12),
    width = unit(6, "cm"),
    height = unit(4, "cm"), 
    heatmap_legend_param = list(direction = "horizontal"),
    column_names_rot = 0,
    column_names_centered = TRUE
)
p_mat

# save
library(showtext)
font_add("Arial", "/usr/share/fonts/windows10/arial.ttf")
font_add("Bold", "/usr/share/fonts/windows10/arialbd.ttf")
showtext_auto()
setFont <- function(p, font_family = "Arial", legend_font_family = "Bold") {
    p@column_names_param$gp$fontfamily <- font_family # 列名称
    p@column_names_param$gp$fontsize <- 12

    p@row_names_param$gp$fontfamily <- font_family # 行名称
    p@row_names_param$gp$fontsize <- 12

    p
}
pdf("plot/scRNA3_CD45plus_I-IFN_heatmap.pdf", width = 6, height = 5)
draw(
    setFont(p_mat),
    heatmap_legend_side = 'top'
)
dev.off()


# IIFN_prod_genes_ref <- unique(unlist(lapply(geneset[pathway_names], function(x) x@geneIds)))
# allgenes <- rownames(Cryo_merge)
# # allgenes <- Myeloid_sub_real@assays$SCT@var.features
# IIFN_prod_genes <- allgenes[which(toupper(allgenes) %in% IIFN_prod_genes_ref)]

# p_dotx <- my_DotPlot_split(
#     Cryo_merge_tmp1, 
#     features = IIFN_prod_genes, 
#     group.by = 'MainCluster',
#     # cols = 'RdBu',
#     scale.min = 0
# ) + theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(), 
#     # legend.position = 'bottom'
# ) + guides(size = FALSE)
# p1<-my_Heatmap(
#     Cryo_merge[,Cryo_merge$orig.ident2 == 'Cryo'],
#     gene=IIFN_prod_genes,group.by = 'MainCluster',slot='scale.data',
#     show_rownames = TRUE, cluster_rows = TRUE
# )
# p2<-my_Heatmap(
#     Cryo_merge[,Cryo_merge$orig.ident2 == 'NonCryo'],
#     gene=IIFN_prod_genes,group.by = 'MainCluster',slot='scale.data',
#     show_rownames = TRUE, cluster_rows = TRUE
# )
# p1+p2
###########