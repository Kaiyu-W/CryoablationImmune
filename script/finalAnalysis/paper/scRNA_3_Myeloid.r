source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
# load("Combined_analysis.rda")
load("Combined_analysis_Myeloid_real.rda")

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

# 1. all Myeloid cells
Myeloid_sub_real
if (!file.exists("scRNA3_Myeloid_cellbarcode.csv")) {
    write.csv(colnames(Myeloid_sub_real), file = "scRNA3_Myeloid_cellbarcode.csv")
}
if (!file.exists("scRNA3_Myeloid_tsne.csv")) {
    write.csv(Myeloid_sub_real@reductions$tsne@cell.embeddings, file = "scRNA3_Myeloid_tsne.csv")
}
if (!file.exists("scRNA3_Myeloid_meta.csv")) {
    write.csv(Myeloid_sub_real@meta.data, file = "scRNA3_Myeloid_meta.csv")
}

Myeloid_sub_real$MainCluster2_new <- factor(
    Myeloid_sub_real$MainCluster2_new,
    levels = rev(c(
        "Monocyte",
        "Monocyte_S100a8/9+",
        "Macrophage_M1",
        "Macrophage_M2",
        "cDC_Itgax+",
        "cDC_Itgax-",
        "pDC",
        "Mast"
    ))
)
Idents(Myeloid_sub_real) <- "MainCluster2_new"

#####
# 1.1 Myeloid tsne
Myeloid_plot <- my_plotDim(
    Myeloid_sub_real,
    reduction = "tsne", label = T,
    pt.size = 2, label.size = 6,
    group.by = "MainCluster2_new",
    title = "Myeloid MainCluster"
) +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "right",
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
            "Monocyte" = "#61baa2",
            "Monocyte_S100a8/9+" = "#709ae1",
            "Macrophage_M1" = "#e98c71",
            "Macrophage_M2" = "#d7b592",
            "cDC_Itgax+" = "#d787b0",
            "cDC_Itgax-" = "#9dcb72",
            "pDC" = "#938cb7",
            "Mast" = "#91a9cc"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
Myeloid_plot

ggsave(
    filename = "plot/scRNA3_Myeloid_dim.pdf",
    plot = Myeloid_plot,
    height = 6, width = 7.69 + 2
)

my_CountCluster(Myeloid_sub_real, group2 = "orig.ident2", group1 = "MainCluster2_new")
#####

#####
# 1.2 Myeloid markers violin
Myeloid_markers <- list(
    Monocytes = c("Cd14", "Fcgr3"),
    Macrophages = c("Itgam", "Adgre1"),
    "M1/cDC" = c("Cd86", "Cd80"),
    M2 = "Mrc1",
    cDC = c("Itgax", "H2-K1", "H2-D1"),
    DCs_activated = "Cd83",
    pDC = c("Bst2", "Siglech"),
    Mast = c("Fcer1a", "Fcer1g", "Kit"),
    Neutrophils = c("S100a8", "S100a9", "Cdk5")
)

Myeloid_dot <- my_DotPlot_split(
    Myeloid_sub_real,
    features = Myeloid_markers
) + RotatedAxis()

Myeloid_violin <- VlnPlot(
    Myeloid_sub_real,
    feature = unlist(Myeloid_markers),
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    cols = c(
        "Monocyte" = "#61baa2",
        "Monocyte_S100a8/9+" = "#709ae1",
        "Macrophage_M1" = "#e98c71",
        "Macrophage_M2" = "#d7b592",
        "cDC_Itgax+" = "#d787b0",
        "cDC_Itgax-" = "#9dcb72",
        "pDC" = "#938cb7",
        "Mast" = "#91a9cc"
    )
) +
    # ggsci::scale_fill_jco() +
    theme(
        aspect.ratio = 8,
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
    filename = "plot/scRNA3_Myeloid_markers_violin.pdf",
    plot = Myeloid_violin,
    height = 6, width = 7.69 + 6
)
#####

#####
# 1.3 Myeloid markers express

am <- list()
# Myeloid_markers
for (i in unlist(Myeloid_markers)) {
    am[[i]] <- FeaturePlot(
        Myeloid_sub_real,
        features = i, reduction = "tsne",
        pt.size = 0.5, cols = c("#e0e7c8", "#cc340a")
    ) +
        theme(
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.title = element_blank(),
            text = element_blank(),
            plot.title = element_blank()
        ) +
        annotate(
            "text",
            x = 31,
            y = -42,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
Myeloid_express <- patchwork::wrap_plots(plots = am, ncol = 4)
ggsave(
    filename = "plot/scRNA3_Myeloid_markers_tsne.pdf",
    plot = Myeloid_express,
    height = 6 * 2, width = 7.69 * 1.5 * 4 / 5
)
#####

#####
# 1.4 Myeloid celltype count
colors_bar <- c(
    "Monocyte" = "#61baa2",
    "Monocyte_S100a8/9+" = "#709ae1",
    "Macrophage_M1" = "#e98c71",
    "Macrophage_M2" = "#d7b592",
    "cDC_Itgax+" = "#d787b0",
    "cDC_Itgax-" = "#9dcb72",
    "pDC" = "#938cb7",
    "Mast" = "#91a9cc"
)
Myeloid_count_sum <- as.data.frame(my_CountCluster(Myeloid_sub_real, group1 = "MainCluster2_new", group2 = "orig.ident2")[-9, 1:2])
Myeloid_count_sum2 <- my_CountCluster(Myeloid_sub_real, group1 = "MainCluster2_new", group2 = "orig.ident2")["sum", 1:2]
Myeloid_count_sum[, 1] <- Myeloid_count_sum[, 1] / Myeloid_count_sum2[1]
Myeloid_count_sum[, 2] <- Myeloid_count_sum[, 2] / Myeloid_count_sum2[2]
Myeloid_count_sum[, 1] <- Myeloid_count_sum[, 1] / sum(Myeloid_count_sum[, 1])
Myeloid_count_sum[, 2] <- Myeloid_count_sum[, 2] / sum(Myeloid_count_sum[, 2])
Myeloid_count_sum
Myeloid_count_sum$cell <- factor(rownames(Myeloid_count_sum), levels = names(colors_bar))
Myeloid_count_sum <- reshape2::melt(Myeloid_count_sum, id.vars = "cell")
Myeloid_count_sum$variable <- sub("_count$", "", Myeloid_count_sum$variable)
Myeloid_count_sum$variable <- ifelse(Myeloid_count_sum$variable == "Cryo", "CA", "Non-CA")
Myeloid_count_sum$variable <- factor(Myeloid_count_sum$variable, levels = c("Non-CA", "CA"))
Myeloid_count_Plot <- ggplot(Myeloid_count_sum, aes(x = variable, y = value, fill = cell)) +
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
    filename = "plot/scRNA3_Myeloid_count.pdf",
    plot = Myeloid_count_Plot,
    height = 6, width = 4.5
)
#####

#####
# 2. Myeloid GO/GSEA

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

# CA_vs_Non-CA with Myeloid cells
Idents(Myeloid_sub_real) <- "orig.ident2"
Markers_Myeloid_list <- FindAllMarkers(Myeloid_sub_real, only.pos = T)
Markers_Myeloid_df <- my_Markers2df_multiple(
    Markers_Myeloid_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 120
)

GO_Myeloid <- my_GO(
    Markers_Myeloid_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_Myeloid_res <- GO_Myeloid@result
dfM_for_plot <- makeTable(GO_Myeloid_res[1:25, ], 55) # filter

Myeloid_GO_bar <- ggplot(
    dfM_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("CA_vs_Non-CA with Myeloid cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
Myeloid_GO_bar
ggsave(
    filename = "plot/scRNA3_Myeloid_GO_bar.pdf",
    plot = Myeloid_GO_bar,
    height = 6, width = 9
)

########### to find where I-IFN signals come from
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

Myeloid_sub_real_tmp1 <- Myeloid_sub_real
Myeloid_sub_real_tmp1$MainCluster2_new <- factor(
    Myeloid_sub_real_tmp1$MainCluster2_new,
    levels = rev(levels(Myeloid_sub_real_tmp1$MainCluster2_new))
)
Myeloid_sub_real_tmp1 <- geneset_score(
    Myeloid_sub_real_tmp1, geneset_IIFN, 
    slot = "data", highly_variable = F,
    method = 'GSVA', 
    gsva_method = 'ssgsea', # "gsva", "ssgsea", "zscore", "plage"
    gsva_kcdf = 'Gaussian'
)
pathway_names <- pathway_names[grepl("(alpha|beta|gamma)", pathway_names, ignore.case = TRUE)]
pathway_names <- pathway_names[!grepl("regulation.*gamma", pathway_names, ignore.case = TRUE)]
p_dot <- my_DotPlot_split(
    Myeloid_sub_real_tmp1, 
    features = pathway_names, 
    group.by = 'MainCluster2_new', 
    cols = 'RdBu'
) + coord_flip() + RotatedAxis() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
)
ggsave(
    filename = "plot/scRNA3_Myeloid_I-IFN_prod_dot.pdf",
    plot = p_dot,
    height = 5, width = 15
)


# different cell types within Myeloid cells
Idents(Myeloid_sub_real) <- "MainCluster2_new"
Markers_Myeloid_list <- FindAllMarkers(Myeloid_sub_real, only.pos = T)
Markers_Myeloid_df <- my_Markers2df_multiple(
    Markers_Myeloid_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 120
)
good_list <- list()
for (col in colnames(Markers_Myeloid_df)) {
    GO_ <- my_GO(
        Markers_Myeloid_df[[col]],
        return_plot = T, return_res = T,
        ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
    )
    if (any(grepl('interferon', GO_@result$Description))) {
        good_list[[col]] <- GO_
        cat(col, " ", GO_@result$Description[grepl('interferon', GO_@result$Description)],"\n")
    }
}
# find where I-IFN came from
for (col in names(good_list)) {
    GO_ <- good_list[[col]]
    if (any(grepl('interferon.*production', GO_@result$Description))) {
        cat(col, ":\n")
        print(GO_@result$Description[grepl('interferon.*production', GO_@result$Description)])
    }
}
# maybe pDC
GO_pDC <- good_list[['Cluster_pDC']]
dfp_for_plot <- makeTable(GO_pDC@result[1:20, ], 53) # filter

pDC_GO_bar <- ggplot(
    dfp_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("pDC cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
pDC_GO_bar
ggsave(
    filename = "plot/scRNA3_pDC_GO_bar.pdf",
    plot = pDC_GO_bar,
    height = 5, width = 9
)

# genes
GO_pDC@result$geneSymbol <- sapply(
    strsplit(GO_pDC@result$geneID, split = "/"),
    function(x)
        paste(
            bitr(
                x, 
                fromType = 'ENTREZID', 
                toType = "SYMBOL", 
                OrgDb = OrgDb
            )$SYMBOL,
            collapse = ","
        )
)
pDC_genes <- unique(unlist(strsplit(
    GO_pDC@result$geneSymbol[grepl('interferon.*production', GO_pDC@result$Description)],
    split = ","
))) # "Tlr7"  "Ptprs" "Flt3"
write.csv(GO_pDC@result, file = "plot/scRNA3_pDC_GO_result.csv")
for (gene in pDC_genes) {
    gene_plot <- FeaturePlot(
        Myeloid_sub_real,
        features = gene,
        reduction = "tsne", # label = T,
        pt.size = 1, label.size = 6,
        order = TRUE, ncol = 1
    )
    ggsave(
        filename = paste0("plot/scRNA3_pDC_",gene,".pdf"),
        plot = gene_plot,
        height = 6, width = 6
    )
}

# gsea
Idents(Myeloid_sub_real) <- "MainCluster2_new"
Markers_Myeloid_list <- FindAllMarkers(
    Myeloid_sub_real, only.pos = F, logfc.threshold = 0,
)
good_gsea_list <- list()
for (col in unique(Markers_Myeloid_list$cluster)) {
    Markers_x_list <- Markers_Myeloid_list[Markers_Myeloid_list$cluster == col, ]
    x_gsea <- my_GSEA(
        Markers_x_list$avg_log2FC, Markers_x_list$gene, 
        ont = 'BP',
        pAdjustMethod = 'BH', # "holm", "hochberg", "hommel", "bonferroni", "BH"(default), "BY", "fdr" and "none"
        return_res = TRUE, pvalueCutoff = 0.1,
        eps = 0,
        nPermSimple = nrow(Markers_x_list)
    )
    x_gsea_IDs <- x_gsea$Gse_GO@result$ID[
        grepl('interferon.*production', x_gsea$Gse_GO@result$Description)
    ]
    if (length(x_gsea_IDs) > 0) {
        good_gsea_list[[col]] <- x_gsea
        cat(col, " ", x_gsea_IDs, "\n")
        plot(gseaplot2(
            x_gsea$Gse_GO, x_gsea_IDs, 
            color = colorspace::rainbow_hcl(4), subplots=1:2,
            pvalue_table = TRUE
        ))
    }
}

# find where I-IFN came from
for (col in names(good_gsea_list)) {
    Gsea_ <- good_gsea_list[[col]]
    if (any(grepl('interferon.*production', Gsea_$Gse_GO@result$Description))) {
        cat(col, ":\n")
        print(Gsea_$Gse_GO@result$Description[grepl('interferon.*production', Gsea_$Gse_GO@result$Description)])
    }
}
Gsea_Mast <- good_gsea_list[['Mast']]
res_ <- Gsea_Mast$Gse_GO@result[Gsea_Mast$Gse_GO@result$p.adjust<0.05,]
gsea_plot <- gseaplot2(
    Gsea_Mast$Gse_GO, res_$ID[grepl('interferon.*production', res_$Description)], 
    color = colorspace::rainbow_hcl(4), subplots=1:2,
    pvalue_table = TRUE
)
ggsave(
    filename = paste0("plot/scRNA3_Mast_gsea.pdf"),
    plot = gsea_plot,
    height = 6, width = 12
)
# # genes
# res_interferon <- res_[grepl('interferon.*production', res_$Description),]
# res_interferon$geneSymbol <- sapply(
#     strsplit(res_interferon$core_enrichment, split = "/"),
#     function(x)
#         paste(
#             bitr(
#                 x, 
#                 fromType = 'ENTREZID', 
#                 toType = "SYMBOL", 
#                 OrgDb = OrgDb
#             )$SYMBOL,
#             collapse = ","
#         )
# )
# Mast_genes <- unique(unlist(strsplit(
#     res_interferon$geneSymbol,
#     split = ","
# )))
# my_DotPlot_split(Myeloid_sub_real, features = Mast_genes, split.by = 'orig.ident2')
# FeaturePlot(Myeloid_sub_real, Mast_genes)

# CA_vs_Non-CA within each cell type of Myeloid cells
Idents(Myeloid_sub_real) <- "MainCluster2_new"
good_obj_list <- list()
for (g in unique(Myeloid_sub_real$MainCluster2_new)) {
    sub_obj <- subset(Myeloid_sub_real, ident = g)
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
        return_plot = T, return_res = T,
        ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
    )
    if (any(grepl('interferon', GO_@result$Description))) {
        good_obj_list[[g]] <- GO_
        cat(g, " ", GO_@result$Description[grepl('interferon', GO_@result$Description)],"\n")
    }
}
# maybe Mast
GO_Mast <- good_obj_list[['Mast']]
dfm_for_plot <- makeTable(GO_Mast@result[1:20, ], 55) # filter

Mast_GO_bar <- ggplot(
    dfm_for_plot,
    aes(
        y = name,
        x = pvalue
    )
) +
    geom_col(orientation = "y", fill = "#4682b4") +
    geom_vline(
        xintercept = c(-log10(0.05)),
        linetype = 8,
        color = "white",
        size = 0.7
    ) +
    ggtitle("pDC cells") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 0.2),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
Mast_GO_bar
ggsave(
    filename = "plot/scRNA3_Mast_GO_bar.pdf",
    plot = Mast_GO_bar,
    height = 5, width = 9
)
# genes
GO_Mast@result$geneSymbol <- sapply(
    strsplit(GO_Mast@result$geneID, split = "/"),
    function(x)
        paste(
            bitr(
                x, 
                fromType = 'ENTREZID', 
                toType = "SYMBOL", 
                OrgDb = OrgDb
            )$SYMBOL,
            collapse = ","
        )
)
Mast_genes <- unique(unlist(strsplit(
    GO_Mast@result$geneSymbol[grepl(
        '(interferon.*production|type I interferon signaling pathway)', 
        GO_Mast@result$Description
    )],
    split = ","
))) # "Irf7"  "Isg15" 
# "Ifi27"  "Ifitm3"
write.csv(GO_Mast@result, file = "plot/scRNA3_Mast_GO_result.csv")
for (gene in Mast_genes) {
    gene_plot <- FeaturePlot(
        Myeloid_sub_real,
        features = gene,
        reduction = "tsne", # label = T,
        pt.size = 1, label.size = 6,
        order = TRUE, ncol = 1
    )
    ggsave(
        filename = paste0("plot/scRNA3_Mast_",gene,".pdf"),
        plot = gene_plot,
        height = 6, width = 6
    )
}
ggsave(
    filename = paste0("plot/scRNA3_Mast_GOgene_violin.pdf"),
    plot = my_violin(Myeloid_sub_real, features = Mast_genes, split.by = 'orig.ident2', mode='mtx'),
)
ggsave(
    filename = paste0("plot/scRNA3_Mast_GOgene_dot.pdf"),
    plot = my_DotPlot_split(Myeloid_sub_real, features = Mast_genes, split.by = 'orig.ident2')
)
###########


#####

#####
# 3. Myeloid Pseudotime (Mono-Macro-DC)
my_AddSeuratPseudo <- function(
    seurat_obj,
    cds_obj,
    reduction_key = "Monocle",
    theta = 0) {
    # @reductions
    if (reduction_key %in% names(seurat_obj@reductions)) {
        warning(reduction_key, " has existed in seurat_obj@reductions! It will be covered.")
    }
    seurat_obj@reductions[[reduction_key]] <- new("DimReduc")
    if (theta == 0) {
        seurat_obj@reductions[[reduction_key]]@cell.embeddings <- t(cds_obj@reducedDimS)
    } else {
        dimPseudotime <- cds_obj@reducedDimS
        theta0 <- theta / 180 * pi # theta: degree; theta0:radian
        wMtx <- sin(theta0) * matrix(c(0, 1, -1, 0), 2, 2) + cos(theta0) * matrix(c(1, 0, 0, 1), 2, 2) # Rotation matrix
        dimPseudotime1 <- wMtx %*% dimPseudotime
        seurat_obj@reductions[[reduction_key]]@cell.embeddings <- t(dimPseudotime1)
    }
    colnames(seurat_obj@reductions[[reduction_key]]@cell.embeddings) <- paste(reduction_key, 1:2, sep = "_")
    seurat_obj@reductions[[reduction_key]]@key <- paste0(reduction_key, "_")
    message("Add $", reduction_key, " into @reduction!")

    # @meta.data
    for (i in c("State", "Pseudotime")) {
        if (i %in% names(cds_obj@phenoData@data)) {
            meta_key <- paste(reduction_key, i, sep = "_")
            seurat_obj@meta.data[[meta_key]] <- cds_obj@phenoData@data[[i]]
            message("Add $", meta_key, " into @meta.data!")
        }
    }

    return(seurat_obj)
}


Idents(Myeloid_sub_real) <- "MainCluster_new"
Myeloid_sub_real_sub <- subset(Myeloid_sub_real, idents = paste0("C", 2:7))
cds_sub <- my_create_monocle(Myeloid_sub_real, idents = paste0("C", 2:7))

cds_sub_diffgenes <- my_featureSelect_cds(
    cds_sub,
    method = "diffgenes",
    seurat_obj = Myeloid_sub_real_sub,
    FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 0.5 | avg_log2FC < -0.5)"
)
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes)
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes, root_state = 2)

Myeloid_sub_real_sub <- my_AddSeuratPseudo(
    Myeloid_sub_real_sub,
    cds_sub_diffgenes,
    reduction_key = "Monocle_Diffgenes",
    theta = -20
)
Idents(Myeloid_sub_real_sub) <- "MainCluster2_new"

# 3.1 Pseudotime reduction dimension
Myeloid_Pseudo <- my_plotDim(
    Myeloid_sub_real_sub,
    reduction = "Monocle_Diffgenes",
    label = T,
    pt.size = 2,
    label.size = 6,
    group.by = "MainCluster2_new",
    title = "Mono-Macro-DC Trajectories"
) +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        legend.position = c(0.1, 0.2),
        legend.text = element_text(face = "bold", size = 12),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
        values = c(
            "Monocyte" = "#61baa2",
            "Macrophage_M1" = "#e98c71",
            "Macrophage_M2" = "#d7b592",
            "cDC_Itgax+" = "#d787b0",
            "cDC_Itgax-" = "#9dcb72",
            "pDC" = "#938cb7"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )

ggsave(
    filename = "plot/scRNA3_Myeloid_pseudo_dim.pdf",
    plot = Myeloid_Pseudo,
    height = 6, width = 7.69
)

# 3.2 Pseudotime value
Myeloid_Pseudo_value1 <- FeaturePlot(
    Myeloid_sub_real_sub,
    features = "Monocle_Diffgenes_Pseudotime",
    reduction = "Monocle_Diffgenes",
    label = T,
    pt.size = 2,
    label.size = 6
)
Myeloid_Pseudo_value1$labels$colour <- "Pseudotime"
Myeloid_Pseudo_value1 <- Myeloid_Pseudo_value1 +
    xlab("Component 1") +
    ylab("Component 2") +
    scale_x_continuous(position = "bottom") +
    scale_y_continuous(position = "right") +
    theme(
        legend.position = c(0.1, 0.2),
        # legend.key.size = unit(0.5, "inch"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # axis.line = element_blank(),
        plot.title = element_blank(),
        text = element_text(size = 12)
    )

dimPseudotime2 <- Myeloid_sub_real_sub@reductions$Monocle_Diffgenes@cell.embeddings
colnames(dimPseudotime2) <- c("x", "y")
dimPseudotime2 <- as.data.frame(dimPseudotime2)
dimPseudotime2$orig.ident <- factor(
    sapply(
        rownames(dimPseudotime2),
        function(x) {
            Myeloid_sub_real_sub@meta.data[x, "orig.ident2"]
        }
    ),
    levels = c("NonCryo", "Cryo")
)
levels(dimPseudotime2$orig.ident) <- c("Non-CA", "CA")

Myeloid_Pseudo_value3 <- ggplot(
    dimPseudotime2,
    aes(x = x, color = orig.ident, fill = orig.ident)
) +
    geom_density(
        mapping = aes(
            x = x, color = orig.ident, linetype = orig.ident
        ),
        alpha = 0.1,
        show.legend = F
    ) +
    scale_color_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "Non-CA" = "longdash",
            "CA" = "solid"
        )
    ) +
    theme_classic() +
    ylab("Density") +
    scale_y_reverse(position = "right") +
    # scale_x_continuous(position = "bottom") +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

Myeloid_Pseudo_value2 <- ggplot(
    dimPseudotime2,
    aes(y = y, color = orig.ident, fill = orig.ident)
) +
    geom_density(
        mapping = aes(
            y = y, color = orig.ident, linetype = orig.ident
        ),
        alpha = 0.1,
        show.legend = F
    ) +
    scale_color_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "Non-CA" = "longdash",
            "CA" = "solid"
        )
    ) +
    theme_classic() +
    xlab("Density") +
    # scale_x_reverse() +
    scale_y_continuous(position = "left") +
    theme(
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )

p_for_legend <- ggplot(
    dimPseudotime2,
    aes(x = x, color = orig.ident, fill = orig.ident)
) +
    geom_density(
        mapping = aes(
            x = x, color = orig.ident, linetype = orig.ident
        ),
        alpha = 0.1,
        show.legend = T
    ) +
    scale_color_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_fill_manual(
        values = c(
            "Non-CA" = "#55a0fb",
            "CA" = "#ff8080"
        )
    ) +
    scale_linetype_manual(
        values = c(
            "Non-CA" = "longdash",
            "CA" = "solid"
        )
    ) +
    guides(
        color = guide_legend(title = "Group"),
        fill = guide_legend(title = "Group"),
        linetype = guide_legend(title = "Group")
    ) +
    theme_classic() +
    theme(
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 12)
    )
require(cowplot)
Myeloid_Pseudo_value4 <- plot_grid(
    get_legend(p_for_legend), NULL,
    ncol = 1
)

Myeloid_Pseudo_value <- plot_grid(
    Myeloid_Pseudo_value1,
    Myeloid_Pseudo_value2,
    Myeloid_Pseudo_value3,
    get_legend(p_for_legend),
    ncol = 2,
    nrow = 2,
    align = "hv",
    axis = "rltb",
    rel_widths = c(4, 1),
    rel_heights = c(3.5, 1)
)

ggsave(
    filename = "plot/scRNA3_Myeloid_pseudo_stat.pdf",
    plot = Myeloid_Pseudo_value,
    height = 6, width = 7.69
)
#####
