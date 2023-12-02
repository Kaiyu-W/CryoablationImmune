library(immunarch)
source("/mnt/e/Cryo-TCR/server/auto/utilities.r")

# 2111 data
setwd("/mnt/e/Cryo-TCR/server/211111/TCR_for_immunarch/")
for (i in dir()) {
    if (grepl("wk.txt$", i)) {
        file_copy <- paste0("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/", i)
        if (!file.exists(file_copy)) {
            file.copy(i, file_copy)
            print(i)
        }
    }
}
# 2106 data
setwd("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/")
for (i in c("CryoTCR", "NonCryoTCR")) {
    file_raw <- paste0("/mnt/e/Cryo-TCR/data/202106_scRNATCR/scTCR/", i, "/outs/filtered_contig_annotations.csv")
    file_copy <- paste0(i, "_1wk_2.txt")
    if (!file.exists(file_copy)) {
        file.copy(file_raw, file_copy)
        print(file_raw)
    }
}

if (!file.exists("T_5end.rds")) {
    # load T5_2111
    # get matrix from server
    datas_2111 <- paste0("/home/kaiyu/Desktop/Cryo/2111_data/", c("Cryo1wk", "Cryo2wk", "NonCryo1wk", "NonCryo2wk"), "_filtered_feature_bc_matrix.h5")
    CA1_2111_mtx <- Seurat::Read10X_h5(datas_2111[1])
    CA2_2111_mtx <- Seurat::Read10X_h5(datas_2111[2])
    NonCA1_2111_mtx <- Seurat::Read10X_h5(datas_2111[3])
    NonCA2_2111_mtx <- Seurat::Read10X_h5(datas_2111[4])
    colnames(CA1_2111_mtx) <- paste0("Cryo1wk2111_", colnames(CA1_2111_mtx))
    colnames(CA2_2111_mtx) <- paste0("Cryo2wk2111_", colnames(CA2_2111_mtx))
    colnames(NonCA1_2111_mtx) <- paste0("NonCryo1wk2111_", colnames(NonCA1_2111_mtx))
    colnames(NonCA2_2111_mtx) <- paste0("NonCryo2wk2111_", colnames(NonCA2_2111_mtx))

    CA1_2111 <- Seurat::CreateSeuratObject(CA1_2111_mtx, min.cells = 0, min.features = 0)
    CA2_2111 <- Seurat::CreateSeuratObject(CA2_2111_mtx, min.cells = 0, min.features = 0)
    NonCA1_2111 <- Seurat::CreateSeuratObject(NonCA1_2111_mtx, min.cells = 0, min.features = 0)
    NonCA2_2111 <- Seurat::CreateSeuratObject(NonCA2_2111_mtx, min.cells = 0, min.features = 0)
    Cryo_2111 <- merge(CA1_2111, c(CA2_2111, NonCA1_2111, NonCA2_2111))

    # load T5_2106
    # get matrix from server
    datas_2106 <- paste0("/home/kaiyu/Desktop/Cryo/2106_data/", c("Cryo", "NonCryo"), "_filtered_feature_bc_matrix.h5")
    CA_2106_mtx <- Seurat::Read10X_h5(datas_2106[1])
    NonCA_2106_mtx <- Seurat::Read10X_h5(datas_2106[2])
    colnames(CA_2106_mtx) <- paste0("Cryo1wk2106_", colnames(CA_2106_mtx))
    colnames(NonCA_2106_mtx) <- paste0("NonCryo1wk2106_", colnames(NonCA_2106_mtx))
    CA_2106 <- Seurat::CreateSeuratObject(CA_2106_mtx, min.cells = 0, min.features = 0)
    NonCA_2106 <- Seurat::CreateSeuratObject(NonCA_2106_mtx, min.cells = 0, min.features = 0)
    Cryo_2106 <- merge(CA_2106, NonCA_2106)

    # load meta and umap
    meta_data <- read.csv("/home/kaiyu/Desktop/Cryo/T5_2106_2111_obs.csv", row.names = 1)
    X_umap <- read.csv("/home/kaiyu/Desktop/Cryo/T5_2106_2111_X_umap_3end.csv", row.names = 1)
    rownames(X_umap) <- rownames(meta_data)
    colnames(X_umap) <- c("UMAP_1", "UMAP_2")

    # filter cells
    Cryo_tmp <- merge(Cryo_2106, Cryo_2111)
    T5_sub <- Cryo_tmp[, rownames(Cryo_tmp@meta.data) %in% rownames(meta_data)]
    T5_sub@meta.data <- meta_data[rownames(T5_sub@meta.data), ]
    T5_sub$IFN_Final[T5_sub$IFN_Final == "Tfh_I-IFN"] <- "CD4_T_I-IFN"
    T5_sub$orig.ident_time <- sub("^.*Cryo", "", T5_sub$orig.ident)

    T5_sub <- my_process_seurat(T5_sub, nVariableFeatures = 2500, norm.method = "SCT")
    T5_sub@reductions$umap@cell.embeddings <- as.matrix(
        X_umap[rownames(T5_sub@reductions$umap@cell.embeddings), ]
    )
    Idents(T5_sub) <- "IFN_Final"
    if (!file.exists("T_5end.rds")) saveRDS(T5_sub, file = "T_5end.rds")
} else {
    T5_sub <- readRDS("T_5end.rds")
}

if (!file.exists("T_5end.h5ad")) {
    library(SeuratDisk)
    # if (file.exists("T_5end.h5Seurat")) file.remove("T_5end.h5Seurat")
    # if (file.exists("T_5end.h5ad")) file.remove("T_5end.h5ad")
    T5_sub_tmp <- T5_sub
    T5_sub_tmp@assays$SCT@scale.data <- matrix()
    SaveH5Seurat(T5_sub_tmp, filename = "T_5end.h5Seurat", overwrite = T)
    Convert("T_5end.h5Seurat", dest = "h5ad", assay = "SCT", overwrite = T)
}

if (!file.exists("T5_cellbarcode.csv")) {
    write.csv(colnames(T5_sub), file = "T5_cellbarcode.csv")
}
if (!file.exists("T5_umap.csv")) {
    write.csv(T5_sub@reductions$umap@cell.embeddings, file = "T5_umap.csv")
}
if (!file.exists("T5_meta.csv")) {
    write.csv(T5_sub@meta.data, file = "T5_meta.csv")
}
# save gene expression
if (!dir.exists("plot/scRNA5_T_matrix") & !file.exists("plot/scRNA5_T_matrix.zip")) {
    dir.create("plot/scRNA5_T_matrix")
    # mtx
    Matrix::writeMM(T5_sub@assays$RNA@counts, "plot/scRNA5_T_matrix/matrix.mtx")
    system("gzip plot/scRNA5_T_matrix/matrix.mtx")
    # scales::number_bytes((file.size("plot/scRNA3_CD45plus_matrix/matrix.mtx.gz")))
    # barcode
    write.table(
        colnames(T5_sub@assays$RNA@counts), 
        file = "plot/scRNA5_T_matrix/barcodes.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F
    )
    system("gzip -f plot/scRNA5_T_matrix/barcodes.tsv")
    # feature
    feature_list = rownames(T5_sub@assays$RNA@counts)
    feature_df = read.table("/mnt/e/Cryo-TCR/server/211111/Cryo_1wk/filtered_feature_bc_matrix/features.tsv.gz", sep = '\t')
    feature_df$V2 <- make.unique(feature_df$V2)
    feature_df <- feature_df[match(feature_list, feature_df$V2), ]
    write.table(
        feature_df, 
        file = "plot/scRNA5_T_matrix/features.tsv", 
        quote = F, 
        row.names = F, 
        col.names = F,
        sep = '\t'
    )
    system("gzip -f plot/scRNA5_T_matrix/features.tsv")
}

#####
# plot dim
T_plot <- my_plotDim(
    T5_sub,
    reduction = "umap", label = T,
    pt.size = 0.5, label.size = 6, repel = T,
    group.by = "IFN_Final",
    title = "T cells MainCluster"
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
            "CD4_Naive" = "#709ae1",
            "CD4_T_I-IFN" = "#197ec0",
            "Th1" = "#d2af81",
            "Th17" = "#bd559f",
            "Treg" = "#f7c530",
            "Tfh" = "#46732e",
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e",
            "T_others" = "#748aa6"
        )
        # values = sapply(
        #     c(
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
        #     ), function(x) {
        #         alpha(alpha = 0.7, colour = x)
        #     }
        # )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
T_plot
ggsave(
    filename = "plot/scRNA5_T_dim.pdf",
    plot = T_plot,
    height = 6, width = 7.69 + 2
)

my_CountCluster(T5_sub, group2 = "orig.ident", group1 = "active.ident")

# plot dim by each times
T_plot1 <- my_plotDim(
    T5_sub,
    reduction = "umap", label = T,
    pt.size = 1, label.size = 6, repel = T,
    group.by = "IFN_Final", cells = T5_sub$orig.ident_time == "1wk",
    title = "T cells MainCluster at 1wk"
) +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "none",
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
            "CD4_Naive" = "#709ae1",
            "CD4_T_I-IFN" = "#197ec0",
            "Th1" = "#d2af81",
            "Th17" = "#bd559f",
            "Treg" = "#f7c530",
            "Tfh" = "#46732e",
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e",
            "T_others" = "#748aa6"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
T_plot2 <- my_plotDim(
    T5_sub,
    reduction = "umap", label = T,
    pt.size = 1, label.size = 6, repel = T,
    group.by = "IFN_Final", cells = T5_sub$orig.ident_time == "2wk",
    title = "T cells MainCluster at 2wk"
) +
    # ggsci::scale_color_tron() +
    # ggsci::scale_color_simpsons() +
    theme(
        # legend.position = c(0.75,0.25),
        legend.position = "none",
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
            "CD4_Naive" = "#709ae1",
            "CD4_T_I-IFN" = "#197ec0",
            "Th1" = "#d2af81",
            "Th17" = "#bd559f",
            "Treg" = "#f7c530",
            "Tfh" = "#46732e",
            "CD8_Naive" = "#d5e4a2",
            "CD8_Tem_I-IFN" = "#ff410d",
            "CD8_Tem" = "#95cc5e",
            "CD8_Tex" = "#d0dfe6",
            "CD8_Tex_Proliferative" = "#f79d1e",
            "T_others" = "#748aa6"
        )
    ) +
    guides(
        color = guide_legend(
            override.aes = list(size = 4)
        )
    )
T_plot1 + T_plot2
ggsave(
    filename = "plot/scRNA5_T_dim_byTimes.pdf",
    plot = T_plot1 + T_plot2,
    height = 6, width = 14
)

# T celltype count
colors_bar <- c(
    "CD4_Naive" = "#709ae1",
    "CD4_T_I-IFN" = "#197ec0",
    "Th1" = "#d2af81",
    "Th17" = "#bd559f",
    "Treg" = "#f7c530",
    "Tfh" = "#46732e",
    "CD8_Naive" = "#d5e4a2",
    "CD8_Tem_I-IFN" = "#ff410d",
    "CD8_Tem" = "#95cc5e",
    "CD8_Tex" = "#d0dfe6",
    "CD8_Tex_Proliferative" = "#f79d1e",
    "T_others" = "#748aa6"
)
T_sub_count_sum <- as.data.frame(my_CountCluster(T5_sub, group1 = "IFN_Final", group2 = "orig.ident")[-13, 1:4])
T_sub_count_sum2 <- my_CountCluster(T5_sub, group1 = "IFN_Final", group2 = "orig.ident")["sum", 1:4]
T_sub_count_sum_tmp <- mapply(function(x, y) x / y, T_sub_count_sum, T_sub_count_sum2)
rownames(T_sub_count_sum_tmp) <- rownames(T_sub_count_sum)
T_sub_count_sum <- as.data.frame(T_sub_count_sum_tmp)
T_sub_count_sum$cell <- factor(rownames(T_sub_count_sum), levels = names(colors_bar))
T_sub_count_sum <- reshape2::melt(T_sub_count_sum, id.vars = "cell")
T_sub_count_sum$variable <- sub("_count$", "", T_sub_count_sum$variable)

T_sub_count_sum$variable <- sub("NonCryo", "Non-CA-", T_sub_count_sum$variable)
T_sub_count_sum$variable <- sub("Cryo", "CA-", T_sub_count_sum$variable)
T_sub_count_sum$variable <- sub("CA-", "CA\n", T_sub_count_sum$variable)
T_sub_count_sum$variable <- factor(
    T_sub_count_sum$variable,
    levels = c("Non-CA\n1wk", "Non-CA\n2wk", "CA\n1wk", "CA\n2wk")
)
T_sub_count_Plot <- ggplot(T_sub_count_sum, aes(x = variable, y = value, fill = cell)) +
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
        axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
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
    filename = "plot/scRNA5_T_count.pdf",
    plot = T_sub_count_Plot,
    height = 6.2, width = 6
)

# T markers violin
CD4_makers <- list(
    Helper = "Cd4",
    Naive = "Sell",
    Th1 = c("Tbx21", "Cxcr3"),
    Th17 = c("Rorc", "Il17a"),
    T_FH = c("Il21r", "Pdcd1"),
    Treg = c("Foxp3", "Il2ra")
)
CD8_makers <- list(
    Cytotoxic = c("Cd8a", "Gzma"),
    Naive = c("Sell"),
    Memory = c("Il7r"),
    Effector = c("Cd44"),
    Exhausted = c("Pdcd1", "Ctla4"),
    Proliferative = c("Mki67")
)

CD4_idents <- c("CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh")
CD8_idents <- c("CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative")
T_CD4_dot <- my_DotPlot_split(
    T5_sub,
    features = CD4_makers,
    idents = CD4_idents
) + RotatedAxis()
T_CD8_dot <- my_DotPlot_split(
    T5_sub,
    features = CD8_makers,
    idents = CD8_idents
) + RotatedAxis()

T5_sub$IFN_Final <- factor(
    as.character(T5_sub$IFN_Final),
    levels = rev(c(CD4_idents, CD8_idents, "T_others"))
)
T_CD4_violin <- VlnPlot(
    T5_sub,
    idents = CD4_idents,
    feature = unlist(CD4_makers),
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    group.by = "IFN_Final",
    cols = c(
        "CD4_Naive" = "#709ae1",
        "CD4_T_I-IFN" = "#197ec0",
        "Th1" = "#d2af81",
        "Th17" = "#bd559f",
        "Treg" = "#f7c530",
        "Tfh" = "#46732e"
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
    filename = "plot/scRNA5_TCD4_markers_violin.pdf",
    plot = T_CD4_violin,
    height = 6, width = 7.69 + 1
)

T_CD8_violin <- VlnPlot(
    T5_sub,
    idents = CD8_idents,
    feature = unlist(CD8_makers),
    pt.size = 0,
    stack = TRUE,
    flip = F,
    fill.by = "ident",
    group.by = "IFN_Final",
    cols = c(
        "CD8_Naive" = "#d5e4a2",
        "CD8_Tem_I-IFN" = "#ff410d",
        "CD8_Tem" = "#95cc5e",
        "CD8_Tex" = "#d0dfe6",
        "CD8_Tex_Proliferative" = "#f79d1e"
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
    filename = "plot/scRNA5_TCD8_markers_violin.pdf",
    plot = T_CD8_violin,
    height = 4.5, width = 7.69 + 1
)

# T markers express
a48 <- list()
markers48_list <- c(
    "Cd4", "Cd8a", "Sell", "Mki67", "Gzma", "Il7r", "Cd44", "Ctla4", "Pdcd1",
    "Il21r", "Tbx21", "Cxcr3", "Rorc", "Il17a", "Foxp3", "Il2ra"
)
for (i in markers48_list) {
    a48[[i]] <- FeaturePlot(
        T5_sub,
        features = i, reduction = "umap",
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
            x = 3.7,
            y = -5,
            label = paste0("italic(", i, ")"),
            parse = T,
            color = "black",
            size = 8
        )
}
TCD48_express <- patchwork::wrap_plots(plots = a48, ncol = 4)
ggsave(
    filename = "plot/scRNA5_TCD48_markers_umap.pdf",
    plot = TCD48_express,
    height = 6 * 2, width = 7.69 * 1.5 * 4 / 5, dpi = 72
)
#####

#####
# load TCR data
meta_info <- data.frame(
    Sample = c(
        "CryoTCR_1wk_2", "CryoTCR_1wk", "CryoTCR_2wk",
        "NonCryoTCR_1wk_2", "NonCryoTCR_1wk", "NonCryoTCR_2wk"
    ),
    group = c("CA", "CA", "CA", "Non-CA", "Non-CA", "Non-CA"),
    time = paste0((c(1, 1, 2, 1, 1, 2)), "wk"),
    prefix = c(
        "Cryo1wk2106", "Cryo1wk2111", "Cryo2wk2111",
        "NonCryo1wk2106", "NonCryo1wk2111", "NonCryo2wk2111"
    ),
    prefix2 = c(
        "Cryo1wk", "Cryo1wk", "Cryo2wk",
        "NonCryo1wk", "NonCryo1wk", "NonCryo2wk"
    )
)
if (!dir.exists("TCR_files")) {
    dir.create("_TCR_files")
    for (idx in seq(rownames(meta_info))) {
        prefix <- meta_info$Sample[idx]
        file <- paste0(prefix, ".txt")
        prefix_replace <- meta_info$prefix[idx]

        data <- read.csv(file)
        data$barcode <- paste0(prefix_replace, "_", data$barcode)
        data$contig_id <- paste0(prefix_replace, "_", data$contig_id)
        data$raw_clonotype_id <- paste0(data$raw_clonotype_id, "_", prefix_replace)
        data$raw_consensus_id <- sub("consensus", paste0(prefix_replace, "_consensus"), data$raw_consensus_id)
        message(file, "\t", nrow(data))
        file_out <- paste0("_TCR_files/", file)
        write.table(data, file = file_out, quote = F, row.names = F, sep = ",")
    }
    file.copy("_TCR_files", "TCR_files", recursive = T)
}
immdata_10x <- immunarch::repLoad("./_TCR_files", .mode = "paired")
immdata_10x$meta <- meta_info
TCR <- immdata_10x

Barcode_RNA <- as.data.frame(T5_sub$IFN_Final)
colnames(Barcode_RNA) <- "Cluster"
Barcode_list_RNA <- list()
for (i in seq(TCR$data)) {
    prefix_tmp <- paste0("^", TCR$meta$prefix2[i])
    name <- names(TCR$data)[i]
    Barcode_list_RNA[[name]] <- Barcode_RNA[grep(prefix_tmp, rownames(Barcode_RNA)), , drop = F]
}
for (i in names(TCR$data)) {
    barcodes_tmp <- rownames(Barcode_list_RNA[[i]])
    TCR$data[[i]] <- select_barcodes(TCR$data[[i]], barcodes_tmp)
}

Barcode_RNA$Clone <- "No"
Barcode_RNA$count <- 0
Barcode_RNA$prop <- 0
Barcode_RNA$aa <- ""
for (i in names(TCR$data)) {
    df <- TCR$data[[i]][, c("Clones", "Proportion", "Barcode", "CDR3.aa")]
    group <- paste(TCR$meta$group, TCR$meta$time, sep = "-")[TCR$meta$Sample == i]
    for (j in 1:nrow(df)) {
        count <- df[j, 1]
        prop <- df[j, 2]
        clonotype <- paste0("C", j, "_", count, "_", group)
        aa <- df[j, 4]
        btmp <- strsplit(df[j, 3], split = ";")[[1]]
        Barcode_RNA[btmp, 2] <- clonotype
        Barcode_RNA[btmp, 3] <- count
        Barcode_RNA[btmp, 4] <- prop
        Barcode_RNA[btmp, 5] <- aa
    }
}

TCR$meta <- meta_info[c(5, 6, 2, 3), ]
rownames(TCR$meta) <- seq(rownames(TCR$meta))

TCR$data[["CryoTCR_1wk"]] <- rbind(TCR$data[["CryoTCR_1wk"]], TCR$data[["CryoTCR_1wk_2"]])
TCR$data[["NonCryoTCR_1wk"]] <- rbind(TCR$data[["NonCryoTCR_1wk"]], TCR$data[["NonCryoTCR_1wk_2"]])
TCR$data[["CryoTCR_1wk_2"]] <- NULL
TCR$data[["NonCryoTCR_1wk_2"]] <- NULL
TCR$data <- TCR$data[c("NonCryoTCR_1wk", "NonCryoTCR_2wk", "CryoTCR_1wk", "CryoTCR_2wk")]
names(TCR$data) <- paste(TCR$meta$group, TCR$meta$time, sep = "-")[match(TCR$meta$Sample, names(TCR$data))]

Barcode_list_RNA[["CryoTCR_1wk"]] <- rbind(Barcode_list_RNA[["CryoTCR_1wk"]], Barcode_list_RNA[["CryoTCR_1wk_2"]])
Barcode_list_RNA[["NonCryoTCR_1wk"]] <- rbind(Barcode_list_RNA[["NonCryoTCR_1wk"]], Barcode_list_RNA[["NonCryoTCR_1wk_2"]])
Barcode_list_RNA[["CryoTCR_1wk_2"]] <- NULL
Barcode_list_RNA[["NonCryoTCR_1wk_2"]] <- NULL
Barcode_list_RNA <- Barcode_list_RNA[c("NonCryoTCR_1wk", "NonCryoTCR_2wk", "CryoTCR_1wk", "CryoTCR_2wk")]
names(Barcode_list_RNA) <- names(TCR$data)

for (i in 1:4) {
    btm <- unlist(lapply(TCR$data[[i]]$Barcode, function(x) strsplit(x, split = ";")[[1]]))
    cat(i, sum(duplicated(btm)), "\n", sep = " ")
    if (sum(duplicated(btm)) != 0) {
        print(TCR$data[[i]][grep(btm[duplicated(btm)], TCR$data[[i]]$Barcode), ])
    }
}
# TCR$data[[4]] <- TCR$data[[4]][-670, ]
TCR$data[[4]] <- TCR$data[[4]][-5944, ]

for (i in names(TCR$data)) {
    barcodes_tmp <- Barcode_list_RNA[[i]]
    TCR$data[[i]]$Celltype <- sapply(TCR$data[[i]]$Barcode, function(x) {
        btmp <- strsplit(x, split = ";")[[1]]
        btmp2 <- barcodes_tmp[btmp, ]
        btmp3 <- table(btmp2)
        btmp4 <- paste0(paste(names(btmp3), btmp3, sep = ":"), collapse = ",")
        btmp4
    })
}


# plot V(D)J usage preference
plot_usage <- function(
    obj,
    only_group = NULL,
    gene_stats = c("musmus.trav", "musmus.traj", "musmus.trbv", "musmus.trbj", "musmus.trbd")[1],
    rank_col = "all",
    grid = TRUE,
    title = NULL,
    change_colname = NULL,
    return_plot = TRUE,
    return_res = FALSE,
    ...) {
    # change_colname is a vector
    # that refers to post-alternative names with custom orders,
    # with name that refers to pre-alternative names
    if (!require(immunarch)) stop("Not install R package immunarch!")

    obj_data <- if (is.null(only_group)) obj$data else obj$data[[only_group]]
    imm_gu <- geneUsage(obj_data, gene_stats, .norm = T, .ambig = "exc")
    imm_gu[is.na(imm_gu)] <- 0

    if (!is.null(rank_col)) {
        rank_method <-
            if (is.null(rank_col)) {
                NULL
            } else if (rank_col == "all") {
                function(x) mean(x[!is.na(x)])
            } else if (rank_col %in% colnames(imm_gu)) {
                function(x) x[colnames(imm_gu)[-1] == rank_col]
            } else if (rank_col %in% change_colname) {
                function(x) x[sapply(colnames(imm_gu)[-1], function(x) change_colname[x]) == rank_col]
            } else {
                stop("ERROR rank_col! if only_group is set as non-NULL, avoid this parameters (NULL) or set it as 'all' (default).")
            }

        mean_tmp <- apply(imm_gu[, -1], 1, rank_method)
        index_tmp <- order(mean_tmp, decreasing = T)
        imm_gu <- imm_gu[index_tmp, ]
        imm_gu[, 1] <- factor(imm_gu[, 1, drop = T], levels = imm_gu[, 1, drop = T])
    }

    if (is.null(only_group)) {
        if (!is.null(change_colname)) {
            colnames(imm_gu) <- c(
                colnames(imm_gu)[1],
                sapply(
                    colnames(imm_gu)[-1],
                    function(x) {
                        change_colname[x]
                    }
                )
            )
            imm_gu <- imm_gu[, c(colnames(imm_gu)[1], change_colname)]
        }
    }

    pp <- vis(imm_gu, .grid = grid, .title = title, ...)

    if (return_plot & return_res) print(pp)
    if (return_plot & !return_res) {
        return(pp)
    }
    if (return_res) {
        return(imm_gu)
    }
}

plot_usage_multigene <- function(
    obj,
    only_group = NULL,
    gene_stats = c("musmus.trav", "musmus.traj", "musmus.trbv", "musmus.trbj", "musmus.trbd")[1:2],
    rank_col = "all",
    grid = TRUE,
    title = NULL,
    change_colname = NULL,
    return_plot = TRUE,
    return_res = FALSE,
    ...) {
    if (length(gene_stats) <= 1) {
        stop("Need more than 2 gene_stats! Use function plot_usage when only 1 gene_stat.")
    }

    imm_gu_list <- lapply(
        seq(gene_stats),
        function(x) {
            plot_usage(
                obj,
                only_group = only_group,
                gene_stats = gene_stats[x],
                rank_col = rank_col,
                grid = FALSE,
                title = NULL,
                change_colname = change_colname,
                return_plot = FALSE,
                return_res = TRUE
            )
        }
    )

    imm_gu <- imm_gu_list[[1]]
    for (i in 2:length(imm_gu_list)) {
        imm_gu <- rbind(imm_gu, imm_gu_list[[i]])
    }

    pp <- vis(imm_gu, .grid = grid, .title = title, ...)

    if (return_plot & return_res) print(pp)
    if (return_plot & !return_res) {
        return(pp)
    }
    if (return_res) {
        return(imm_gu)
    }
}

TcrAVJplot <- plot_usage_multigene(
    TCR,
    gene_stats = c("musmus.trav", "musmus.traj"),
    rank_col = "Non-CA-1wk",
    title = "Mus.TcrA-VJ",
    change_colname = NULL,
    .add.layer = list(
        theme(
            axis.text.x = element_text(vjust = 0.5, size = 5),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15)
        ),
        # geom_vline(xintercept = c(20.5,22.5),linetype=2,col = 'red')
        geom_rect(aes(xmin = 0, xmax = 75.5, ymin = 0, ymax = Inf), fill = "#e1a3b0", alpha = .2),
        geom_rect(aes(xmin = 75.5, xmax = 114.5, ymin = 0, ymax = Inf), fill = "#72c8b5", alpha = .2)
    ),
    return_res = F,
    return_plot = T
)

TcrBVDJplot <- plot_usage_multigene(
    TCR,
    gene_stats = c("musmus.trbv", "musmus.trbd", "musmus.trbj"),
    rank_col = "Non-CA-1wk",
    title = "Mus.TcrB-VDJ",
    change_colname = NULL,
    .add.layer = list(
        theme(
            axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15)
        ),
        # geom_vline(xintercept = c(20.5,22.5),linetype=2,col = 'red')
        geom_rect(aes(xmin = 0, xmax = 20.5, ymin = 0, ymax = Inf), fill = "#e1a3b0", alpha = .2),
        geom_rect(aes(xmin = 20.5, xmax = 22.5, ymin = 0, ymax = Inf), fill = "#ffffff", alpha = .2),
        geom_rect(aes(xmin = 22.5, xmax = 33.5, ymin = 0, ymax = Inf), fill = "#72c8b5", alpha = .2)
    ),
    return_res = F,
    return_plot = T
)

TcrAVJplot_value <- plot_usage_multigene(
    TCR,
    gene_stats = c("musmus.trav", "musmus.traj"),
    rank_col = "Non-CA-1wk",
    title = "Mus.TcrA-VJ",
    change_colname = NULL,
    .add.layer = list(
        theme(
            axis.text.x = element_text(vjust = 0.5, size = 5),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15)
        ),
        # geom_vline(xintercept = c(20.5,22.5),linetype=2,col = 'red')
        geom_rect(aes(xmin = 0, xmax = 75.5, ymin = 0, ymax = Inf), fill = "#e1a3b0", alpha = .2),
        geom_rect(aes(xmin = 75.5, xmax = 114.5, ymin = 0, ymax = Inf), fill = "#72c8b5", alpha = .2)
    ),
    return_res = T,
    return_plot = F
)
TcrBVDJplot_value <- plot_usage_multigene(
    TCR,
    gene_stats = c("musmus.trbv", "musmus.trbd", "musmus.trbj"),
    rank_col = "Non-CA-1wk",
    title = "Mus.TcrB-VDJ",
    change_colname = NULL,
    .add.layer = list(
        theme(
            axis.text.x = element_text(vjust = 0.5, size = 10),
            axis.title.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 15)
        ),
        # geom_vline(xintercept = c(20.5,22.5),linetype=2,col = 'red')
        geom_rect(aes(xmin = 0, xmax = 20.5, ymin = 0, ymax = Inf), fill = "#e1a3b0", alpha = .2),
        geom_rect(aes(xmin = 20.5, xmax = 22.5, ymin = 0, ymax = Inf), fill = "#ffffff", alpha = .2),
        geom_rect(aes(xmin = 22.5, xmax = 33.5, ymin = 0, ymax = Inf), fill = "#72c8b5", alpha = .2)
    ),
    return_res = T,
    return_plot = F
)

# Hmisc::rcorr(TcrBVDJplot_value$`Non-CA-1wk`, TcrBVDJplot_value$`CA-1wk`, type = 'p')$r
# Hmisc::rcorr(TcrBVDJplot_value$`Non-CA-2wk`, TcrBVDJplot_value$`CA-2wk`, type = 'p')$r
# Hmisc::rcorr(TcrAVJplot_value$`Non-CA-1wk`, TcrAVJplot_value$`CA-1wk`, type = 'p')$r
# Hmisc::rcorr(TcrAVJplot_value$`Non-CA-2wk`, TcrAVJplot_value$`CA-2wk`, type = 'p')$r
for (i in c("TcrAVJplot_value", "TcrBVDJplot_value")) {
    tmp <- get(i)
    res1 <- eval(parse(text = "Hmisc::rcorr(tmp$`Non-CA-1wk`, tmp$`CA-1wk`, type = 'p')"))
    res2 <- eval(parse(text = "Hmisc::rcorr(tmp$`Non-CA-2wk`, tmp$`CA-2wk`, type = 'p')"))
    cat(i, " 1wk: PCC:", res1$r[2, 1], " pvalue:", res1$P[2, 1], "\n", sep = "")
    cat(i, " 2wk: PCC:", res2$r[2, 1], " pvalue:", res2$P[2, 1], "\n", sep = "")
}

ggsave(TcrAVJplot, filename = "plot/scRNA5_TcrAVJplot.pdf", width = 15, height = 7)
ggsave(TcrBVDJplot, filename = "plot/scRNA5_TcrBVDJplot.pdf", width = 10, height = 7)
#####

#####
# plot TCR components between clusters
T5_sub$clone_TCR <- Barcode_RNA$Clone
T5_sub$TCR_if <- ifelse(T5_sub$clone_TCR == "No", "No_TCR", "TCR")


df_TCR_if1 <- table(T5_sub@meta.data[, c("TCR_if", "IFN_Final")])
df_TCR_if2 <- reshape2::melt(df_TCR_if1)
df_TCR_if3 <- reshape2::melt(apply(df_TCR_if1, 2, function(x) 100 * x / sum(x)))
df_TCR_if2$IFN_Final <- factor(df_TCR_if2$IFN_Final, levels = rev(names(colors_bar)))
df_TCR_if3$IFN_Final <- factor(df_TCR_if3$IFN_Final, levels = rev(names(colors_bar)))

pal <- RColorBrewer::brewer.pal(ncol(df_TCR_if1) - 1, "Paired")
p_TCR <- ggplot(df_TCR_if2, aes(x = IFN_Final, y = value, fill = TCR_if)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Cell type") +
    ylab("Cell count") +
    labs(title = "TCR detected by scTCR-seq in T cells of scRNA-seq") +
    scale_fill_manual(values = c("grey75", pal)) +
    theme_bw() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    coord_flip()


p_TCRp <- ggplot(df_TCR_if3, aes(x = IFN_Final, y = value, fill = TCR_if)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Cell type") +
    ylab("Prop, %") +
    labs(title = "TCR detected by scTCR-seq in T cells of scRNA-seq") +
    scale_fill_manual(values = c("grey75", pal)) +
    theme_bw() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    coord_flip()


ggsave(p_TCR, filename = "plot/scRNA5_TCR_detected_byCluster.pdf", width = 8, height = 4)
ggsave(p_TCRp, filename = "plot/scRNA5_TCR_detected_prob_byCluster.pdf", width = 8, height = 4)




# plot TCR Clonal expansion between clusters
Clone_n_ge_2 <- names(table(Barcode_RNA$Clone))[table(Barcode_RNA$Clone) > 1]
T5_sub$clone_TCR_ge_2 <- sapply(
    Barcode_RNA$Clone,
    function(x) {
        ifelse(x %in% Clone_n_ge_2, x, "No")
    }
)
T5_sub$TCR_ge_2 <- ifelse(T5_sub$clone_TCR_ge_2 == "No", "No_TCR", "TCR")
df_TCR_if4 <- sapply(
    unique(T5_sub$orig.ident),
    function(x) {
        mtx_tmp <- table(T5_sub@meta.data[T5_sub$orig.ident == x, c("TCR_ge_2", "IFN_Final")])
        mtx_tmp[2, ]
    }
)
# df_TCR_if4 <- df_TCR_if4 / nrow(T5_sub@meta.data)
# df_TCR_if4 <- apply(df_TCR_if4, 2, function(x) x / sum(x))
df_TCR_if4 <- sapply(
    colnames(df_TCR_if4),
    function(x) {
        df_TCR_if4[, x] / sum(T5_sub$orig.ident == x)
    }
)
# re-order & rename
colnames(df_TCR_if4) <- c("CA-1wk", "Non-CA-1wk", "CA-2wk", "Non-CA-2wk")
df_TCR_if4 <- df_TCR_if4[
    c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    ),
    c(2:1, 4, 3) # put control in front
]
annotation_col <- data.frame(
    Group = c("Non-CA", "CA", "Non-CA", "CA"),
    Time = c("1wk", "1wk", "2wk", "2wk")
)

RWB_color <- grDevices::colorRampPalette(
    colors = c("#f0f2f4", "#f4bda5", "#f03a3a")
)(50)
p_TCR_mtx <- ComplexHeatmap::Heatmap(
    df_TCR_if4,
    name = "Clonal expansion",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 2),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    # column_names_rot = 45,
    heatmap_legend_param = list(
        at = seq(0, 0.4, 0.1),
        direction = "horizontal",
        labels = paste0(seq(0, 40, 10), "%"),
        title_position = "topcenter"
    ),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        Time = factor(annotation_col$Time, levels = c("1wk", "2wk")),
        col = list(
            Group = c("Non-CA" = "grey", "CA" = "orange"),
            Time = c("1wk" = "#80c498", "2wk" = "#c076a6")
        )
    )
)

pdf("plot/scRNA5_TCR_clonalExpansion_byCluster_matrixplot.pdf", width = 4, height = 8)
draw(p_TCR_mtx, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# plot TCR Clonal-Size between clusters
Clone_size <- sapply(
    table(Barcode_RNA$Clone),
    function(x) {
        c(x == 1, x == 2, x >= 3)
    }
)
T5_sub$clone_TCR_size <- sapply(
    Barcode_RNA$Clone,
    function(x) {
        ifelse(x == "No", "=0", c("=1", "=2", ">=3")[Clone_size[, x]])
    }
)
df_TCR_if5 <- table(T5_sub@meta.data[, c("orig.ident", "clone_TCR_size")])
df_TCR_if5 <- df_TCR_if5[, 2:4]
# re-order & rename
rownames(df_TCR_if5) <- c("CA\n1wk", "CA\n2wk", "Non-CA\n1wk", "Non-CA\n2wk")
df_TCR_if5 <- df_TCR_if5[
    c(3, 1, 4, 2) # put control in front
    ,
]
df_TCR_if5_all <- apply(df_TCR_if5, 1, function(x) x / sum(x) * 100)
df_TCR_if5_Group <- apply(
    apply(
        df_TCR_if5, 2, function(x) {
            c(
                "Non-CA" = unname(x[1] + x[3]),
                "CA" = unname(x[2] + x[4])
            )
        }
    ), 1, function(x) x / sum(x) * 100
)
df_TCR_if5_Time <- apply(
    apply(
        df_TCR_if5, 2, function(x) {
            c(
                "1wk" = unname(x[1] + x[2]),
                "2wk" = unname(x[3] + x[4])
            )
        }
    ), 1, function(x) x / sum(x) * 100
)

# proportion barplot
colors_bar2 <- c(
    "=0" = "#f8eef5",
    "=1" = "#e5bbd9",
    "=2" = "#ca77b2",
    ">=3" = "#97447f"
)
df2_TCR_if5_all <- reshape2::melt(df_TCR_if5_all)
p_TCR_bar <- ggplot(df2_TCR_if5_all, aes(x = orig.ident, y = value, fill = clone_TCR_size)) +
    geom_bar(width = 0.8, stat = "identity") +
    ylab("Fraction of clone size (%)") +
    labs(fill = "Clone size") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "grey", size = 2)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = colors_bar2
    )
ggsave(
    filename = "plot/scRNA5_TCR_clonalSize_bar.pdf",
    plot = p_TCR_bar,
    height = 6, width = 4.5
)

df2_TCR_if5_Group <- reshape2::melt(df_TCR_if5_Group)
p_TCR_bar_Group <- ggplot(df2_TCR_if5_Group, aes(x = Var2, y = value, fill = clone_TCR_size)) +
    geom_bar(width = 0.8, stat = "identity") +
    ylab("Fraction of clone size (%)") +
    labs(fill = "Clone size") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "grey", size = 2)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = colors_bar2
    )
ggsave(
    filename = "plot/scRNA5_TCR_clonalSize_byGroup_bar.pdf",
    plot = p_TCR_bar_Group,
    height = 6, width = 3
)

df2_TCR_if5_Time <- reshape2::melt(df_TCR_if5_Time)
p_TCR_bar_Time <- ggplot(df2_TCR_if5_Time, aes(x = Var2, y = value, fill = clone_TCR_size)) +
    geom_bar(width = 0.8, stat = "identity") +
    ylab("Fraction of clone size (%)") +
    labs(fill = "Clone size") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "grey", size = 2)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = colors_bar2
    )
ggsave(
    filename = "plot/scRNA5_TCR_clonalSize_byTime_bar.pdf",
    plot = p_TCR_bar_Time,
    height = 6, width = 3
)


# plot TCR Clonal expansion between clusters by Clone-Size
df_TCR_if6 <- lapply(c("=2", ">=3"), function(y) {
    sapply(
        unique(T5_sub$orig.ident),
        function(x) {
            mtx_tmp <- table(
                T5_sub@meta.data[
                    T5_sub$orig.ident == x & T5_sub$clone_TCR_size == y,
                    c("TCR_ge_2", "IFN_Final")
                ]
            )
            mtx_tmp["TCR", ]
        }
    )
})

df_TCR_if6 <- lapply(1:2, function(y) {
    sapply(
        colnames(df_TCR_if6[[y]]),
        function(x) {
            df_TCR_if6[[y]][, x] / sum(T5_sub$orig.ident == x)
        }
    )
})
# re-order & rename
colnames(df_TCR_if6[[1]]) <- c("CA-1wk", "Non-CA-1wk", "CA-2wk", "Non-CA-2wk")
colnames(df_TCR_if6[[2]]) <- c("CA-1wk", "Non-CA-1wk", "CA-2wk", "Non-CA-2wk")
df_TCR_if6[[1]] <- df_TCR_if6[[1]][
    c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    ),
    c(2:1, 4, 3) # put control in front
]
df_TCR_if6[[2]] <- df_TCR_if6[[2]][
    c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    ),
    c(2:1, 4, 3) # put control in front
]
annotation_col <- data.frame(
    Group = c("Non-CA", "CA", "Non-CA", "CA"),
    Time = c("1wk", "1wk", "2wk", "2wk")
)

RWB_color <- grDevices::colorRampPalette(
    colors = c("#f0f2f4", "#f4bda5", "#f03a3a")
)(50)
p_TCR_mtx_2 <- ComplexHeatmap::Heatmap(
    df_TCR_if6[[1]],
    name = "Clonal expansion",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 2),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    # column_names_rot = 45,
    heatmap_legend_param = list(
        at = seq(0, 0.3, 0.1),
        direction = "horizontal",
        labels = paste0(seq(0, 30, 10), "%"),
        title_position = "topcenter"
    ),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        Time = factor(annotation_col$Time, levels = c("1wk", "2wk")),
        col = list(
            Group = c("Non-CA" = "grey", "CA" = "orange"),
            Time = c("1wk" = "#80c498", "2wk" = "#c076a6")
        )
    )
)
p_TCR_mtx_3 <- ComplexHeatmap::Heatmap(
    df_TCR_if6[[2]],
    name = "Clonal expansion",
    col = RWB_color,
    rect_gp = gpar(col = "white", lwd = 2),
    cluster_columns = F,
    cluster_rows = F,
    show_row_dend = F,
    show_column_names = F,
    # column_names_rot = 45,
    heatmap_legend_param = list(
        at = seq(0, 0.15, 0.05),
        direction = "horizontal",
        labels = paste0(seq(0, 15, 5), "%"),
        title_position = "topcenter"
    ),
    top_annotation = HeatmapAnnotation(
        Group = factor(annotation_col$Group, levels = c("Non-CA", "CA")),
        Time = factor(annotation_col$Time, levels = c("1wk", "2wk")),
        col = list(
            Group = c("Non-CA" = "grey", "CA" = "orange"),
            Time = c("1wk" = "#80c498", "2wk" = "#c076a6")
        )
    )
)

pdf("plot/scRNA5_TCR_clonalExpansion_byCluster_matrixplot_size2.pdf", width = 4, height = 8)
draw(p_TCR_mtx_2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
pdf("plot/scRNA5_TCR_clonalExpansion_byCluster_matrixplot_size3+.pdf", width = 4, height = 8)
draw(p_TCR_mtx_3, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


# exhaustion / cytotoxic / Treg / Naïve score by Clone Size
my_cor <- function(GEX, gene) {
    gene_exp <- GEX[gene, ]
    # gene_exp_others <- GEX[rownames(GEX) != gene, ]
    gene_exp_others <- GEX
    pb <- txtProgressBar(style = 3)
    iii <- 0
    n <- nrow(gene_exp_others)
    res_mtx <- apply(
        gene_exp_others,
        1,
        function(x) {
            res <- Hmisc::rcorr(gene_exp, x, type = "p")
            iii <<- iii + 1
            setTxtProgressBar(pb, iii / n)
            c("r" = res$r[2, 1], "p" = res$P[2, 1])
        }
    )
    close(pb)
    t(res_mtx)
}
my_select_topCor <- function(df, n, p_thre = 0.05) {
    idx0 <- unlist(apply(df, 1, function(x) is.na(x[1]) | is.nan(x[1] | is.na(x[2]) | is.nan(x[2]))))
    df <- as.data.frame(df[!idx0, , drop = F])

    idx <- df$p <= p_thre
    if (sum(idx) == 0) stop("No gene passed p-value threshold!")
    df <- df[idx, , drop = F]
    df <- df[order(df$r, decreasing = T), , drop = F]
    ntop <- min(n, nrow(df))
    df[1:ntop, ]
}
my_plotCor <- function(df, fill_col = "red") {
    df$gene <- factor(
        rownames(df),
        levels = rev(rownames(df))
    )
    ggplot(df, aes(x = r, y = gene)) +
        geom_col(color = "black", fill = fill_col, width = 0.8) +
        scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
        xlab("Correlation") +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(face = "italic", color = "black", size = 12),
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_line(color = "black")
        )
}

# get correlation and select genes
Idents(T5_sub) <- "IFN_Final"
# 1. Exhaustion: LAG3
Exhaustion_cells <- subset(T5_sub, idents = c("CD8_Tex", "CD8_Tex_Proliferative"))
LAG3_cor <- my_cor(Exhaustion_cells@assays$SCT@data, "Lag3")
LAG3_cor_df <- my_select_topCor(LAG3_cor, 20)
# 2. Cytotoxic: NKG7
Cytotoxic_cells <- subset(T5_sub, idents = c("CD8_Tem", "CD8_Tem_I-IFN"))
NKG7_cor <- my_cor(Cytotoxic_cells@assays$SCT@data, "Nkg7")
NKG7_cor_df <- my_select_topCor(NKG7_cor, 20)
# 3. Naive: CCR7
Naive_cells <- subset(T5_sub, idents = c("CD4_Naive", "CD8_Naive"))
CCR7_cor <- my_cor(Naive_cells@assays$SCT@data, "Ccr7")
CCR7_cor_df <- my_select_topCor(CCR7_cor, 20)
# 4. Tregs: FOXP3
Tregs_cells <- subset(T5_sub, idents = c("Treg"))
Tregs_cor <- my_cor(Tregs_cells@assays$SCT@data, "Foxp3")
Tregs_cor_df <- my_select_topCor(Tregs_cor, 20)
# plot
score_list <- list(
    "Exhaustion" = c("LAG3", "#d0dfe6"),
    "Cytotoxic" = c("NKG7", "#95cc5e"),
    "Naive" = c("CCR7", "#709ae1"),
    "Tregs" = c("Tregs", "#f7c530")
)
for (i in names(score_list)) {
    fill_color <- score_list[[i]][2]
    gene <- score_list[[i]][1]
    eval(parse(text = paste0(
        "ggsave(
            filename = 'plot/scRNA5_TCR_", i, "Score_cor.pdf',
            plot = my_plotCor(", gene, "_cor_df, fill_color),
            height = 4, width = 3
        )"
    )))
}

# train the linear model for classfication
my_linearModel_SignitureCluster <- function(obj, gene_list, labels, group.by, slot = "data", assay = "RNA") {
    GEX <- methods::slot(obj@assays[[assay]], slot)
    gene_list <- intersect(gene_list, rownames(GEX))
    data <- as.data.frame(t(GEX[gene_list, ]))
    data$label <- as.factor(ifelse(
        obj@meta.data[colnames(GEX), group.by] %in% labels,
        1,
        0
    ))
    my_formula <- paste("label ~", paste(gene_list, collapse = " + "))
    fit <- eval(parse(text = paste0("glm(", my_formula, ", data = data, family = binomial())")))
    fit
}
my_linearModel_SignitureScore.default <- function(obj, fit, use_predict = FALSE, slot = "data", assay = "RNA") {
    GEX <- methods::slot(obj@assays[[assay]], slot)
    data <- as.data.frame(t(GEX))
    res <- if (use_predict) {
        predict(fit, newdata = data, type = "response")
    } else {
        coeff <- fit$coefficients
        coeff <- coeff[names(coeff) != "(Intercept)"]
        apply(
            data[, names(coeff), drop = F],
            1,
            function(x) sum(x * coeff)
        )
    }
    return(res)
}
my_linearModel_SignitureScore <- function(...) {
    res <- my_linearModel_SignitureScore.default(...)
    gc()
    return(res)
}
my_plotScore <- function(
    obj, score_col, group.by,
    idents = NULL, label = NULL,
    fill_col = "red", jitter = TRUE,
    outlier.shape = NA, outlier.size = 1, outlier.alpha = 0.1,
    angle = 45, x_vjust = 1, x_hjust = 1,
    use_ggpubr = F, comparisons_list = NULL,
    use_kw = F, kw_x = 1, kw_y = 0.9
) {
    df <- obj@meta.data[, c(score_col, group.by)]
    rownames(df) <- 1:nrow(df)

    if (!is.null(idents)) {
        df <- df[df[[group.by]] %in% idents, , drop = F]
        df[[group.by]] <- factor(df[[group.by]], levels = idents)
    }
    df$tmp1 <- df[[group.by]]
    df$tmp2 <- df[[score_col]]

    if (use_ggpubr) {
        if (!require(ggpubr)) stop("No ggpubr!")
        p <- ggboxplot(
            data = df,
            x = "tmp1",
            y = "tmp2",
            fill = fill_col,
            outlier.shape = outlier.shape,
            outlier.size = outlier.size,
            outlier.alpha = outlier.alpha
        ) +
            scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
            ylab(sub("^ ", "", paste(label, "Score"))) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.text.y = element_text(color = "black", size = 12),
                axis.text.x = element_text(
                    color = "black", size = 12,
                    angle = angle, vjust = x_vjust, hjust = x_hjust
                ),
                axis.title.y = element_text(size = 12),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(color = "black")
            )
        if (use_kw) {
            p <- p + stat_compare_means(
                label.y = kw_y, label.x = kw_x,
                method = "kruskal.test"
            )
        }
        if (!is.null(comparisons_list)) {
            p <- p + stat_compare_means(
                label = "p.signif",
                comparisons = comparisons_list
            )
        }
    } else {
        p <- ggplot(df, aes(x = tmp1, y = tmp2)) +
            geom_boxplot(
                outlier.shape = outlier.shape,
                fill = fill_col,
                outlier.size = outlier.size,
                outlier.alpha = outlier.alpha
            ) +
            scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
            ylab(sub("^ ", "", paste(label, "Score"))) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.text.y = element_text(color = "black", size = 12),
                axis.text.x = element_text(
                    color = "black", size = 12,
                    angle = angle, vjust = x_vjust, hjust = x_hjust
                ),
                axis.title.y = element_text(size = 12),
                axis.title.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(color = "black")
            )
    }
    if (jitter) p + geom_jitter(size = 0.01, alpha = 0.1) else p
}
my_plotScore_violin <- function(
    obj, score_col, group.by,
    idents = NULL, label = NULL,
    fill_col = "red",
    palette = NULL,
    draw_quantiles = FALSE, 
    draw_boxplot = FALSE,
    boxplot_width = 0.1,
    jitter_alpha = 0.01,
    jitter_width = 0.1,
    jitter_size = 0.1,
    outlier.shape = NA, outlier.size = 1, outlier.alpha = 0.1,
    angle = 45, x_vjust = 1, x_hjust = 1,
    comparisons_list = NULL,
    use_kw = F, kw_x = 1, kw_y = 0.9
) {
    df <- obj@meta.data[, c(score_col, group.by)]
    rownames(df) <- 1:nrow(df)

    if (!is.null(idents)) {
        df <- df[df[[group.by]] %in% idents, , drop = F]
        df[[group.by]] <- factor(df[[group.by]], levels = idents)
    }
    df$tmp1 <- df[[group.by]]
    df$tmp2 <- df[[score_col]]

    if (!require(ggpubr)) stop("No ggpubr!")
    p <- ggviolin(
        data = df,
        x = "tmp1",
        y = "tmp2",
        fill = fill_col,
        palette = palette,
        width = 0.9,
        trim = TRUE,
        draw_quantiles = if (draw_quantiles) c(0.25,0.5,0.75) else NULL,
        error.plot = 'crossbar',
        add = if (draw_boxplot) c("boxplot", "jitter") else "none", 
        add.params  = list(
            color = 'black',
            fill = 'white',
            alpha = jitter_alpha,
            size = jitter_size,
            jitter = jitter_width,
            width = boxplot_width
        ),
        scale = 'width'
    ) +
        scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0.05, 0.05)) +
        ylab(sub("^ ", "", paste(label, "Score"))) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(
                color = "black", size = 12,
                angle = angle, vjust = x_vjust, hjust = x_hjust
            ),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(color = "black")
        )
    if (use_kw) {
        p <- p + stat_compare_means(
            label.y = kw_y, label.x = kw_x,
            method = "kruskal.test"
        )
    }
    if (!is.null(comparisons_list)) {
        p <- p + stat_compare_means(
            label = "p.signif",
            comparisons = comparisons_list
        )
    }

    for (i in seq(length(p$layers))) {
        if (class(p$layers[[i]]$geom)[1] == 'GeomBoxplot') {
            p$layers[[i]]$geom_params$outlier.size <- outlier.size
            p$layers[[i]]$geom_params$outlier.shape <- outlier.shape
            p$layers[[i]]$geom_params$outlier.alpha <- outlier.alpha
            p$layers[[2]]$aes_params$alpha <- 1
        }
    }
    p
}

Idents(T5_sub) <- "IFN_Final"
Train_cells <- subset(T5_sub, idents = c(
    "CD4_Naive", "CD8_Naive", "Treg",
    "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex", "CD8_Tex_Proliferative"
))
# 1. Exhaustion: LAG3
fit1 <- my_linearModel_SignitureCluster(
    Train_cells, rownames(LAG3_cor_df),
    label = c("CD8_Tex", "CD8_Tex_Proliferative"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Ex_score <- my_linearModel_SignitureScore(T5_sub, fit1, use_predict = T, assay = "SCT")
p_Ex <- my_plotScore(
    T5_sub, "Ex_score", "IFN_Final",
    idents = c("CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"),
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 4), c(1, 5), c(2, 4), c(2, 5), c(3, 4), c(3, 5), 4:5),
    use_kw = T, kw_y = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion.pdf",
    plot = p_Ex,
    height = 4, width = 5
)
# 2. Cytotoxic: NKG7
fit2 <- my_linearModel_SignitureCluster(
    Train_cells, rownames(NKG7_cor_df),
    label = c("CD8_Tem_I-IFN", "CD8_Tem"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Cy_score <- my_linearModel_SignitureScore(T5_sub, fit2, use_predict = T, assay = "SCT")
p_Cy <- my_plotScore(
    T5_sub, "Cy_score", "IFN_Final",
    idents = c("CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"),
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(2, 5), c(3, 4), c(3, 5)),
    use_kw = T, kw_y = 1.7
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic.pdf",
    plot = p_Cy,
    height = 4, width = 5
)
# 3. Naive: CCR7
fit3 <- my_linearModel_SignitureCluster(
    Train_cells, rownames(CCR7_cor_df),
    label = c("CD4_Naive", "CD8_Naive"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Na_score <- my_linearModel_SignitureScore(T5_sub, fit3, use_predict = T, assay = "SCT")
p_Na <- my_plotScore(
    T5_sub, "Na_score", "IFN_Final",
    idents = c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ),
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(
        c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(1, 6),
        c(7, 8), c(7, 9), c(7, 10), c(7, 11)
    ),
    use_kw = T, kw_y = 1.9, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive.pdf",
    plot = p_Na,
    height = 4, width = 7
)
# 4. Tregs: FOXP3
fit4 <- my_linearModel_SignitureCluster(
    Train_cells, rownames(Tregs_cor_df),
    label = c("Treg"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Tr_score <- my_linearModel_SignitureScore(T5_sub, fit4, use_predict = T, assay = "SCT")
p_Tr <- my_plotScore(
    T5_sub, "Tr_score", "IFN_Final",
    idents = c("CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg"),
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 5), c(2, 5), c(3, 5), c(4, 5)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg.pdf",
    plot = p_Tr,
    height = 4, width = 5.5
)
# 5. by Clone-Size
# 5.1 Exhaustion
p_Ex_size <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score", "clone_TCR_size",
    idents = c("=1", "=2", ">=3"),
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(2, 3), c(1, 3)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_byCloneSize.pdf",
    plot = p_Ex_size,
    height = 4, width = 3.5
)
# 5.2 Cytotoxic
p_Cy_size <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Cy_score", "clone_TCR_size",
    idents = c("=1", "=2", ">=3"),
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(2, 3), c(1, 3)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_byCloneSize.pdf",
    plot = p_Cy_size,
    height = 4, width = 3.5
)
# 5.3 Naive
p_Na_size <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Na_score", "clone_TCR_size",
    label = "Naive",
    idents = c("=1", "=2", ">=3"),
    fill_col = "#709ae1",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(2, 3), c(1, 3)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_byCloneSize.pdf",
    plot = p_Na_size,
    height = 4, width = 3.5
)
# 5.4 Treg
p_Tr_size <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    )],
    "Tr_score", "clone_TCR_size",
    idents = c("=1", "=2", ">=3"),
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(2, 3), c(1, 3)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_byCloneSize.pdf",
    plot = p_Tr_size,
    height = 4, width = 3.5
)

# ggpubr::compare_means(Tr_score~clone_TCR_size, T5_sub@meta.data[
#         T5_sub$IFN_Final %in% c(
#             "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
#         ) & T5_sub$clone_TCR_size %in% c("=1", "=2", ">=3"),
#     ])

# by Clone-Size & Case-Control
T5_sub$orig_ident_tmp <- factor(
    T5_sub$orig.ident,
    levels = c("NonCryo1wk", "Cryo1wk", "NonCryo2wk", "Cryo2wk")
)
levels(T5_sub$orig_ident_tmp) <- c("Non-CA\n1wk", "CA\n1wk", "Non-CA\n2wk", "CA\n2wk")
# 5.1 Exhaustion
p_Ex_size_1 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1"],
    "Ex_score", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
p_Ex_size_2 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2"],
    "Ex_score", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
p_Ex_size_3 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3"],
    "Ex_score", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size1.pdf",
    plot = p_Ex_size_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size2.pdf",
    plot = p_Ex_size_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size3.pdf",
    plot = p_Ex_size_3,
    height = 4, width = 3.5
)
p_Ex_size123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Ex_score", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
p_Ex_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size123.pdf",
    plot = p_Ex_size123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_all.pdf",
    plot = p_Ex_size_all,
    height = 4, width = 3.5
)
# 5.2 Cytotoxic
p_Cy_size_1 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1"],
    "Cy_score", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Cy_size_2 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2"],
    "Cy_score", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Cy_size_3 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3"],
    "Cy_score", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size1.pdf",
    plot = p_Cy_size_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size2.pdf",
    plot = p_Cy_size_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size3.pdf",
    plot = p_Cy_size_3,
    height = 4, width = 3.5
)
p_Cy_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Cy_score", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Cy_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Cy_score", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_all.pdf",
    plot = p_Cy_size_all,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size123.pdf",
    plot = p_Cy_size_123,
    height = 4, width = 3.5
)
# 5.3 Naive
p_Na_size_1 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1"],
    "Na_score", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Na_size_2 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2"],
    "Na_score", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Na_size_3 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3"],
    "Na_score", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size1.pdf",
    plot = p_Na_size_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size2.pdf",
    plot = p_Na_size_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size3.pdf",
    plot = p_Na_size_3,
    height = 4, width = 3.5
)
p_Na_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Na_score", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Na_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Na_score", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
# ggpubr::compare_means(Na_score~orig_ident_tmp, T5_sub@meta.data[
#         T5_sub$IFN_Final %in% c(
#         "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
#         "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
#     )& T5_sub$clone_TCR_size != "=0", ], alter = 'greater')
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size123.pdf",
    plot = p_Na_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_all.pdf",
    plot = p_Na_size_all,
    height = 4, width = 3.5
)
# 5.3 Treg
p_Tr_size <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    )],
    "Tr_score", "clone_TCR_size",
    idents = c("=1", "=2", ">=3"),
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(2, 3), c(1, 3)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_byCloneSize.pdf",
    plot = p_Tr_size,
    height = 4, width = 3.5
)

p_Tr_size_1 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=1"],
    "Tr_score", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Tr_size_2 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=2"],
    "Tr_score", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Tr_size_3 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == ">=3"],
    "Tr_score", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size1.pdf",
    plot = p_Tr_size_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size2.pdf",
    plot = p_Tr_size_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size3.pdf",
    plot = p_Tr_size_3,
    height = 4, width = 3.5
)
p_Tr_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Tr_score", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Tr_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    )],
    "Tr_score", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size123.pdf",
    plot = p_Tr_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_all.pdf",
    plot = p_Tr_size_all,
    height = 4, width = 3.5
)


# exhaustion / cytotoxic / Treg / Naïve score by paper genes by Case-Control
paper_gene_list <- list(
    cytotoxic = c('Nkg7','Ccl5','Gzma','Gzmk','Ccl4','Cst7','Itm2c','Ifng'),
    exhaustion = c('Lag3','Gzmb','Havcr2','Ptms','Cxcl13','Vcam1','Prf1','Tnfrsf9','Tigit','Pdcd1'),
    naive = c('Ccr7','Tcf7','Lef1','Actn1','S1pr1','Mal','Il7r','Plac8','Spint2','Bach2'),
    Treg = c('Foxp3','Il2ra','Tnfrsf18','Il1r2','Tnfrsf4','Ccr8','Layn','Ctla4','Tigit')
)
# 1. train
# 1. Exhaustion: LAG3
fit1p <- my_linearModel_SignitureCluster(
    Train_cells, paper_gene_list$exhaustion,
    label = c("CD8_Tex", "CD8_Tex_Proliferative"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Ex_score_p <- my_linearModel_SignitureScore(T5_sub, fit1p, use_predict = T, assay = "SCT")
p_Exp <- my_plotScore(
    T5_sub, "Ex_score_p", "IFN_Final",
    idents = c("CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"),
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 4), c(1, 5), c(2, 4), c(2, 5), c(3, 4), c(3, 5), 4:5),
    use_kw = T, kw_y = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_byPaper.pdf",
    plot = p_Exp,
    height = 4, width = 5
)
# 2. Cytotoxic: NKG7
fit2p <- my_linearModel_SignitureCluster(
    Train_cells, paper_gene_list$cytotoxic,
    label = c("CD8_Tem_I-IFN", "CD8_Tem"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Cy_score_p <- my_linearModel_SignitureScore(T5_sub, fit2p, use_predict = T, assay = "SCT")
p_Cyp <- my_plotScore(
    T5_sub, "Cy_score_p", "IFN_Final",
    idents = c("CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"),
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 2), c(1, 3), c(2, 3), c(2, 4), c(2, 5), c(3, 4), c(3, 5)),
    use_kw = T, kw_y = 1.7
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_byPaper.pdf",
    plot = p_Cyp,
    height = 4, width = 5
)
# 3. Naive: CCR7
fit3p <- my_linearModel_SignitureCluster(
    Train_cells, paper_gene_list$naive,
    label = c("CD4_Naive", "CD8_Naive"),
    group.by = "IFN_Final", assay = "RNA"
)
T5_sub$Na_score_p <- my_linearModel_SignitureScore(T5_sub, fit3p, use_predict = T, assay = "RNA")
p_Nap <- my_plotScore(
    T5_sub, "Na_score_p", "IFN_Final",
    idents = c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ),
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(
        c(1, 2), c(1, 3), c(1, 4), c(1, 5), c(1, 6),
        c(7, 8), c(7, 9), c(7, 10), c(7, 11)
    ),
    use_kw = T, kw_y = 1.9, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_byPaper.pdf",
    plot = p_Nap,
    height = 4, width = 7
)
# 4. Tregs: FOXP3
fit4p <- my_linearModel_SignitureCluster(
    Train_cells, paper_gene_list$Treg,
    label = c("Treg"),
    group.by = "IFN_Final", assay = "SCT"
)
T5_sub$Tr_score_p <- my_linearModel_SignitureScore(T5_sub, fit4p, use_predict = T, assay = "SCT")
p_Trp <- my_plotScore(
    T5_sub, "Tr_score_p", "IFN_Final",
    idents = c("CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg"),
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    use_ggpubr = T, comparisons_list = list(c(1, 5), c(2, 5), c(3, 5), c(4, 5)),
    use_kw = T, kw_y = 1.4
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_byPaper.pdf",
    plot = p_Trp,
    height = 4, width = 5.5
)
# 2.score
###########################################################################
# 5.1 Exhaustion
p_Exp_size123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
p_Exp_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size123_byPaper.pdf",
    plot = p_Exp_size123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_all_byPaper.pdf",
    plot = p_Exp_size_all,
    height = 4, width = 3.5
)
p_Exp_size123 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
p_Exp_size_all <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = "#d0dfe6",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.2, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size123_byPaper_violin.pdf",
    plot = p_Exp_size123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_all_byPaper_violin.pdf",
    plot = p_Exp_size_all,
    height = 4, width = 3.5
)
# 5.2 Cytotoxic
p_Cyp_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Cyp_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_all_byPaper.pdf",
    plot = p_Cyp_size_all,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size123_byPaper.pdf",
    plot = p_Cyp_size_123,
    height = 4, width = 3.5
)
p_Cyp_size_123 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Cyp_size_all <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = "#95cc5e",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_all_byPaper_violin.pdf",
    plot = p_Cyp_size_all,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size123_byPaper_violin.pdf",
    plot = p_Cyp_size_123,
    height = 4, width = 3.5
)
# 5.3 Naive
p_Nap_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Nap_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
# ggpubr::compare_means(Na_score_p~orig_ident_tmp, T5_sub@meta.data[
#         T5_sub$IFN_Final %in% c(
#         "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
#         "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
#     ),], alter = 'less')
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size123_byPaper.pdf",
    plot = p_Nap_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_all_byPaper.pdf",
    plot = p_Nap_size_all,
    height = 4, width = 3.5
)
p_Nap_size_123 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Nap_size_all <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    )],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = "#709ae1",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size123_byPaper_violin.pdf",
    plot = p_Nap_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_all_byPaper_violin.pdf",
    plot = p_Nap_size_all,
    height = 4, width = 3.5
)
# # only CD4/8 Naive
# p_Nap_size_123 <- my_plotScore(
#     T5_sub[, T5_sub$IFN_Final %in% c(
#         "CD4_Naive", 
#         "CD8_Naive"
#     ) & T5_sub$clone_TCR_size != "=0"],
#     "Na_score_p", "orig_ident_tmp",
#     label = "Naive",
#     fill_col = "#709ae1",
#     jitter = F,
#     angle = 0,
#     x_vjust = 0.5, x_hjust = 0.5,
#     use_ggpubr = T, comparisons_list = list(1:2, 3:4),
#     use_kw = T, kw_y = 1.25, kw_x = 1.5
# )
# p_Nap_size_all <- my_plotScore(
#     T5_sub[, T5_sub$IFN_Final %in% c(
#         "CD4_Naive", 
#         "CD8_Naive"
#     )],
#     "Na_score_p", "orig_ident_tmp",
#     label = "Naive",
#     fill_col = "#709ae1",
#     jitter = F,
#     angle = 0,
#     x_vjust = 0.5, x_hjust = 0.5,
#     use_ggpubr = T, comparisons_list = list(1:2, 3:4),
#     use_kw = T, kw_y = 1.25, kw_x = 1.5
# )
# # ggpubr::compare_means(Na_score_p~orig_ident_tmp, T5_sub@meta.data[
# #         T5_sub$IFN_Final %in% c(
# #         "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
# #         "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
# #     ),], alter = 'less')
# ggsave(
#     filename = "plot/scRNA5_TCR_Score_Naive_Size123_byPaper_onlyNaive.pdf",
#     plot = p_Nap_size_123,
#     height = 4, width = 3.5
# )
# ggsave(
#     filename = "plot/scRNA5_TCR_Score_Naive_all_byPaper_onlyNaive.pdf",
#     plot = p_Nap_size_all,
#     height = 4, width = 3.5
# )
# 5.4 Treg
p_Trp_size_123 <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Trp_size_all <- my_plotScore(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    )],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size123_byPaper.pdf",
    plot = p_Trp_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_all_byPaper.pdf",
    plot = p_Trp_size_all,
    height = 4, width = 3.5
)
p_Trp_size_123 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size != "=0"],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
p_Trp_size_all <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    )],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = "#f7c530",
    jitter = F,
    draw_boxplot = T, draw_quantiles = F,
    boxplot_width = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    use_ggpubr = T, comparisons_list = list(1:2, 3:4),
    use_kw = T, kw_y = 1.25, kw_x = 1.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size123_byPaper_violin.pdf",
    plot = p_Trp_size_123,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_all_byPaper_violin.pdf",
    plot = p_Trp_size_all,
    height = 4, width = 3.5
)
###########################################################################
# by Clone-Size & Case-Control
T5_sub$orig_ident_tmp <- factor(
    T5_sub$orig.ident,
    levels = c("NonCryo1wk", "Cryo1wk", "NonCryo2wk", "Cryo2wk")
)
levels(T5_sub$orig_ident_tmp) <- c("Non-CA\n1wk", "CA\n1wk", "Non-CA\n2wk", "CA\n2wk")
# 5.1 Exhaustion
p_Exp_size_1_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '1wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.3,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Exp_size_2_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '1wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Exp_size_3_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '1wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Exp_size_1_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '2wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Exp_size_2_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '2wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Exp_size_3_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '2wk'],
    "Ex_score_p", "orig_ident_tmp",
    label = "Exhaustion",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size1_byPaper_1wk.pdf",
    plot = p_Exp_size_1_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size2_byPaper_1wk.pdf",
    plot = p_Exp_size_2_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size3_byPaper_1wk.pdf",
    plot = p_Exp_size_3_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size1_byPaper_2wk.pdf",
    plot = p_Exp_size_1_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size2_byPaper_2wk.pdf",
    plot = p_Exp_size_2_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Exhaustion_Size3_byPaper_2wk.pdf",
    plot = p_Exp_size_3_2,
    height = 4, width = 3.5
)
# 5.2 Cytotoxic
p_Cyp_size_1_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '1wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.3,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Cyp_size_2_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '1wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Cyp_size_3_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '1wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Cyp_size_1_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '2wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Cyp_size_2_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '2wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Cyp_size_3_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '2wk'],
    "Cy_score_p", "orig_ident_tmp",
    label = "Cytotoxic",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size1_byPaper_1wk.pdf",
    plot = p_Cyp_size_1_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size2_byPaper_1wk.pdf",
    plot = p_Cyp_size_2_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size3_byPaper_1wk.pdf",
    plot = p_Cyp_size_3_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size1_byPaper_2wk.pdf",
    plot = p_Cyp_size_1_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size2_byPaper_2wk.pdf",
    plot = p_Cyp_size_2_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Cytotoxic_Size3_byPaper_2wk.pdf",
    plot = p_Cyp_size_3_2,
    height = 4, width = 3.5
)
# 5.3 Naive
p_Nap_size_1_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '1wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.3,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Nap_size_2_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '1wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Nap_size_3_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '1wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Nap_size_1_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '2wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Nap_size_2_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '2wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Nap_size_3_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem", "CD8_Tem_I-IFN", "CD8_Tex_Proliferative", "CD8_Tex"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '2wk'],
    "Na_score_p", "orig_ident_tmp",
    label = "Naive",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size1_byPaper_1wk.pdf",
    plot = p_Nap_size_1_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size2_byPaper_1wk.pdf",
    plot = p_Nap_size_2_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size3_byPaper_1wk.pdf",
    plot = p_Nap_size_3_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size1_byPaper_2wk.pdf",
    plot = p_Nap_size_1_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size2_byPaper_2wk.pdf",
    plot = p_Nap_size_2_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Naive_Size3_byPaper_2wk.pdf",
    plot = p_Nap_size_3_2,
    height = 4, width = 3.5
)
# 5.3 Tr4g
p_Trp_size_1_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '1wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.3,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Trp_size_2_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '1wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Trp_size_3_1 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '1wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Trp_size_1_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=1" & T5_sub$orig.ident_time == '2wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.05,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Trp_size_2_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == "=2" & T5_sub$orig.ident_time == '2wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
p_Trp_size_3_2 <- my_plotScore_violin(
    T5_sub[, T5_sub$IFN_Final %in% c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh"
    ) & T5_sub$clone_TCR_size == ">=3" & T5_sub$orig.ident_time == '2wk'],
    "Tr_score_p", "orig_ident_tmp",
    label = "Treg",
    fill_col = 'tmp1',
    palette = c('#55a0fb','#ff8080'),
    draw_quantiles = F, draw_boxplot = T, 
    boxplot_width = 0.1, jitter_alpha = 0.1,
    angle = 0,
    x_vjust = 0.5, x_hjust = 0.5,
    comparisons_list = list(1:2), use_kw = F
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size1_byPaper_1wk.pdf",
    plot = p_Trp_size_1_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size2_byPaper_1wk.pdf",
    plot = p_Trp_size_2_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size3_byPaper_1wk.pdf",
    plot = p_Trp_size_3_1,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size1_byPaper_2wk.pdf",
    plot = p_Trp_size_1_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size2_byPaper_2wk.pdf",
    plot = p_Trp_size_2_2,
    height = 4, width = 3.5
)
ggsave(
    filename = "plot/scRNA5_TCR_Score_Treg_Size3_byPaper_2wk.pdf",
    plot = p_Trp_size_3_2,
    height = 4, width = 3.5
)


# TCR clone top distribution
data_list_top <- lapply(
    TCR$data,
    function(x) {
        x$top <- 1:nrow(x)
        tmp1 <- lapply(x$Barcode, function(y) strsplit(y, split = ";")[[1]])
        x$n <- sapply(tmp1, length)
        df <- data.frame(
            "cb" = unlist(tmp1),
            "ntop" = unlist(apply(x[, c("n", "top")], 1, function(y) rep(y[2], y[1]))),
            "CDR3.aa" = unlist(apply(x[, c("n", "CDR3.aa")], 1, function(y) rep(y[2], y[1]))),
            "freq" = unlist(apply(x[, c("n", "Proportion")], 1, function(y) rep(y[2], y[1]))),
            "count" = unlist(apply(x[, c("n", "Clones")], 1, function(y) rep(y[2], y[1])))
        )
        df$cluster <- T5_sub@meta.data[df$cb, "IFN_Final"]
        df
    }
)
intervals <- c(10, 30, 100)
top_clone_byCluster_list <- lapply(
    data_list_top,
    function(x) {
        y <- tapply(x$ntop, x$cluster, table)
        res <- list()
        for (i in seq(length(intervals) + 1)) {
            if (i == 1) {
                i_min <- 0
                i_max <- intervals[1]
            } else if (i == length(intervals) + 1) {
                i_min <- intervals[length(intervals)]
                i_max <- Inf
            } else {
                i_min <- intervals[i - 1]
                i_max <- intervals[i]
            }
            res[[i]] <- c(i_min, i_max)
        }
        mtx <- sapply(
            y,
            function(z) sapply(res, function(k) sum(z[names(z) > k[1] & names(z) <= k[2]]))
        )
        rownames(mtx) <- c(intervals, Inf)
        as.data.frame(mtx)
    }
)
top_clone_byCluster_list <- lapply(
    top_clone_byCluster_list,
    function(x) {
        df <- x[, colnames(x) != "T_others"]
        df2 <- reshape2::melt(
            apply(df, 2, function(y) y / sum(y) * 100)
        )
        df2$Var2 <- factor(df2$Var2, levels = c(
            "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
            "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
        ))
        df2$Var1 <- factor(df2$Var1, levels = c(
            intervals, Inf
        ))
        levels(df2$Var1) <- c(paste("<=", intervals), "Others")
        df2
    }
)
colors_bar3 <- c(
    "<= 10" = "#e84518",
    "<= 30" = "#f2d06f",
    "<= 100" = "#b2dbf1",
    "Others" = "#26429e"
)
for (i in names(top_clone_byCluster_list)) {
    p_TCR_bar_top <- ggplot(top_clone_byCluster_list[[i]], aes(x = Var2, y = value, fill = Var1)) +
        geom_bar(width = 0.95, stat = "identity", color = "black") +
        ylab("TCR repertoire (%)") +
        labs(fill = "Top Clones") +
        theme_bw() +
        theme(
            legend.text = element_text(color = "black", size = 12),
            legend.position = "right",
            legend.title = element_text(color = "black", size = 12),
            axis.title.x = element_blank(),
            axis.title.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.line = element_blank(),
            axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.background = element_blank(),
            panel.border = element_rect(color = "grey", size = 2)
        ) +
        guides(fill = guide_legend(reverse = FALSE)) +
        scale_fill_manual(
            values = rev(colors_bar3)
        ) +
        scale_y_continuous(expand = c(0.01, 0, 0.025, 0))
    ggsave(
        filename = paste0("plot/scRNA5_TCR_topBar_", i, ".pdf"),
        plot = p_TCR_bar_top,
        height = 6, width = 6
    )
}

# TCR clone `novel, expanded and contracted` distribution
# 1. 1wk Case vs Control
# names(data_list_top) # 1 / 3
all_aa_1 <- unique(c(data_list_top[[1]][[3]], data_list_top[[3]][[3]]))
all_aa_df1 <- data.frame(
    aa = all_aa_1,
    case = sapply(all_aa_1, function(x) {
        idx <- data_list_top[[3]]$CDR3.aa == x
        if (any(idx)) data_list_top[[3]][which(idx)[1], "freq"] else 0
    }),
    control = sapply(all_aa_1, function(x) {
        idx <- data_list_top[[1]]$CDR3.aa == x
        if (any(idx)) data_list_top[[1]][which(idx)[1], "freq"] else 0
    }),
    case_count = sapply(all_aa_1, function(x) {
        idx <- data_list_top[[3]]$CDR3.aa == x
        if (any(idx)) data_list_top[[3]][which(idx)[1], "count"] else 0
    })
)
all_aa_df1$fc <- apply(all_aa_df1[, c("case", "control")], 1, function(x) {
    if (x[1] == 0) -Inf else if (x[2] == 0) Inf else if (x[1] > x[2]) x[1] / x[2] else x[2] / x[1]
})
all_aa_df1$type <- factor(
    apply(all_aa_df1[, c("case", "control", "fc", "case_count")], 1, function(x) {
        if (x[1] > 0 & x[2] == 0) {
            if (x[4] > 1) "Novel" else "Persistent"
        } else if ((x[3] < exp(0.5) & x[3] >= 1) | x[1] == 0) {
            "Persistent"
        } else if (x[1] > x[2]) {
            "Expanded"
        } else {
            "Contracted"
        }
    }),
    levels = c("Novel", "Persistent", "Expanded", "Contracted")
)
colors_bar4 <- c(
    "Novel" = "#279a2a",
    "Persistent" = "#bfbcbb",
    "Expanded" = "#bb291e",
    "Contracted" = "#4b50aa"
)
p_aaDot1 <- ggplot(all_aa_df1, aes(x = control, y = case, fill = type, col = type, shape = type)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = 5) +
    ylab("Clone frequency of CA") +
    xlab("Clone frequency of Non-CA") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 12),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect("transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5, lineend = "square"),
        axis.ticks.x = element_line(color = "black", size = 0.5, lineend = "square"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1)
    ) +
    scale_fill_manual(
        values = colors_bar4
    ) +
    scale_color_manual(
        values = colors_bar4
    ) +
    scale_shape_manual(values = c(17, 20, 15, 18))
ggsave(
    filename = paste0("plot/scRNA5_TCR_aaType_dot_1wk.pdf"),
    plot = p_aaDot1,
    height = 4, width = 6
)

# 2. 2wk Case vs Control
# names(data_list_top) # 2 / 4
all_aa_2 <- unique(c(data_list_top[[2]][[3]], data_list_top[[4]][[3]]))
all_aa_df2 <- data.frame(
    aa = all_aa_2,
    case = sapply(all_aa_2, function(x) {
        idx <- data_list_top[[4]]$CDR3.aa == x
        if (any(idx)) data_list_top[[4]][which(idx)[1], "freq"] else 0
    }),
    control = sapply(all_aa_2, function(x) {
        idx <- data_list_top[[2]]$CDR3.aa == x
        if (any(idx)) data_list_top[[2]][which(idx)[1], "freq"] else 0
    }),
    case_count = sapply(all_aa_2, function(x) {
        idx <- data_list_top[[4]]$CDR3.aa == x
        if (any(idx)) data_list_top[[4]][which(idx)[1], "count"] else 0
    })
)
all_aa_df2$fc <- apply(all_aa_df2[, c("case", "control")], 1, function(x) {
    if (x[1] == 0) -Inf else if (x[2] == 0) Inf else if (x[1] > x[2]) x[1] / x[2] else x[2] / x[1]
})
all_aa_df2$type <- factor(
    apply(all_aa_df2[, c("case", "control", "fc", "case_count")], 1, function(x) {
        if (x[1] > 0 & x[2] == 0) {
            if (x[4] > 1) "Novel" else "Persistent"
        } else if ((x[3] < exp(0.5) & x[3] >= 1) | x[1] == 0) {
            "Persistent"
        } else if (x[1] > x[2]) {
            "Expanded"
        } else {
            "Contracted"
        }
    }),
    levels = c("Novel", "Persistent", "Expanded", "Contracted")
)
p_aaDot2 <- ggplot(all_aa_df2, aes(x = control, y = case, fill = type, col = type, shape = type)) +
    geom_point(size = 3) +
    geom_abline(intercept = 0, slope = 1, linetype = 5) +
    ylab("Clone frequency of CA") +
    xlab("Clone frequency of Non-CA") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 12),
        legend.position = c(0.8, 0.6),
        legend.background = element_rect("transparent"),
        legend.title = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.line = element_blank(),
        axis.ticks.y = element_line(color = "black", size = 0.5, lineend = "square"),
        axis.ticks.x = element_line(color = "black", size = 0.5, lineend = "square"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1)
    ) +
    scale_fill_manual(
        values = colors_bar4
    ) +
    scale_color_manual(
        values = colors_bar4
    ) +
    scale_shape_manual(values = c(17, 20, 15, 18))
ggsave(
    filename = paste0("plot/scRNA5_TCR_aaType_dot_2wk.pdf"),
    plot = p_aaDot2,
    height = 4, width = 5
)

# 3. merge into T5_sub
data_list_1wk <- rbind(data_list_top[[1]], data_list_top[[3]])
data_list_2wk <- rbind(data_list_top[[2]], data_list_top[[4]])
T5_sub$aa <- "No"
T5_sub@meta.data[data_list_1wk$cb, "aa"] <- data_list_1wk$CDR3.aa
T5_sub@meta.data[data_list_2wk$cb, "aa"] <- data_list_2wk$CDR3.aa
T5_sub$aa_type <- "No"

all_aa_df1$type <- as.character(all_aa_df1$type)
all_aa_df1 <- rbind(all_aa_df1, "No" = c("", "", "", "", "", "No"))
all_aa_df2$type <- as.character(all_aa_df2$type)
all_aa_df2 <- rbind(all_aa_df2, "No" = c("", "", "", "", "", "No"))
idx_from1 <- T5_sub$orig.ident_time == "1wk"
aa_type_from1 <- all_aa_df1[T5_sub$aa[idx_from1], "type"]
aa_type_from2 <- all_aa_df2[T5_sub$aa[!idx_from1], "type"]

T5_sub$aa_type[idx_from1] <- aa_type_from1
T5_sub$aa_type[!idx_from1] <- aa_type_from2

# 4. barplot
type_clone_byCluster_list <- lapply(
    lapply(
        unique(T5_sub$orig.ident),
        function(x) {
            y <- T5_sub@meta.data[T5_sub$orig.ident == x, ]
            table(y$aa_type, y$IFN_Final)
        }
    ),
    function(x) {
        df <- x[, colnames(x) != "T_others"]
        df2 <- reshape2::melt(
            apply(df, 1, function(y) y / sum(y) * 100)
        )
        df2$Var1 <- factor(df2$Var1, levels = c(
            "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
            "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
        ))
        df2$Var2 <- factor(df2$Var2, levels = c(
            c("No", "Persistent", "Contracted", "Expanded", "Novel")
        ))
        df2
    }
)
names(type_clone_byCluster_list) <- unique(T5_sub$orig.ident)

for (i in names(type_clone_byCluster_list)) {
    p_TCR_bar_aatype <- ggplot(type_clone_byCluster_list[[i]], aes(x = Var2, y = value, fill = Var1)) +
        geom_bar(width = 0.95, stat = "identity", color = "white") +
        ylab("Fraction of Cells (%)") +
        labs(fill = "Cell Type") +
        theme_bw() +
        theme(
            legend.text = element_text(color = "black", size = 12),
            legend.position = "right",
            legend.title = element_text(color = "black", size = 12),
            axis.title.x = element_blank(),
            axis.title.y = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            axis.line = element_blank(),
            axis.ticks.y = element_line(color = "grey", size = 1, lineend = "square"),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.background = element_blank(),
            panel.border = element_rect(color = "grey", size = 2)
        ) +
        guides(fill = guide_legend(reverse = FALSE)) +
        scale_fill_manual(
            values = rev(colors_bar)
        ) +
        scale_y_continuous(expand = c(0.01, 0, 0.025, 0))
    ggsave(
        filename = paste0("plot/scRNA5_TCR_aatypeBar_", i, ".pdf"),
        plot = p_TCR_bar_aatype,
        height = 6, width = if (i == "NonCryo2wk") 4.2 else 6
    )
}


######
Barcode_RNA$group <- sapply(
    strsplit(rownames(Barcode_RNA), split = "_"),
    function(x) x[1]
)
top3_prop_list <- tapply(
    Barcode_RNA$prop,
    Barcode_RNA$group,
    function(x) sort(unique(x), decreasing = T)[1:3]
)
top3_prop_list <- top3_prop_list[c(2:3, 5:6)]
top3_clone_df <- as.data.frame(sapply(
    seq(top3_prop_list),
    function(x) {
        name <- names(top3_prop_list)[x]
        df <- Barcode_RNA[Barcode_RNA$group == name, ]
        sort(unique(df$Clone[df$prop %in% top3_prop_list[[x]]]), decreasing = F)
    }
))
colnames(top3_clone_df) <- c("CA\n1wk", "CA\n2wk", "Non-CA\n1wk", "Non-CA\n2wk")
top3_aa_df <- as.data.frame(apply(
    top3_clone_df,
    2,
    function(x) {
        sapply(x, function(y) {
            unique(Barcode_RNA$aa[Barcode_RNA$Clone == y])
        })
    }
))

top3_aaType_df <- apply(top3_clone_df, 2, function(x) {
    sapply(x, function(y) {
        unique(T5_sub$aa_type[rownames(Barcode_RNA)[Barcode_RNA$Clone == y]])
    })
})
top3_aaClone_list <- apply(top3_aa_df, 2, function(x) {
    sapply(x, function(y) {
        unique(Barcode_RNA$Clone[Barcode_RNA$aa == y])
    })
})
# filter
for (i in seq(top3_aaClone_list)) {
    time_tmp <- sub("^.*\n", "", names(top3_aaClone_list)[i])
    for (j in seq(top3_aaClone_list[[i]])) {
        top3_aaClone_list[[i]][[j]] <- top3_aaClone_list[[i]][[j]][grep(time_tmp, top3_aaClone_list[[i]][[j]])]
    }
}
top3_aaClone_list[[1]] <- top3_aaClone_list[[1]][c(2, 1, 3)]
top3_aaType_df[, 1] <- top3_aaType_df[c(2, 1, 3), 1]
top3_aa_df[, 1] <- top3_aa_df[c(2, 1, 3), 1]
top3_clone_df[, 1] <- top3_clone_df[c(2, 1, 3), 1]

T5_sub$orig.ident_group <- ifelse(grepl("NonCryo", T5_sub$orig.ident), "Non-CA", "CA")
top3_CellCount_list <- lapply(top3_aaClone_list, function(x) {
    sapply(x, function(y) {
        Clones <- y
        res <- table(T5_sub@meta.data[T5_sub$clone_TCR %in% Clones, c("IFN_Final", "orig.ident_group")])
        if (ncol(res) == 1) {
            new_col <- setdiff(c("CA", "Non-CA"), colnames(res))
            res <- cbind(res, tmp = rep(0, nrow(res)))
            colnames(res)[2] <- new_col
        }
        res
    }, simplify = F)
})
top3_CellCount_list <- top3_CellCount_list[c(3:4, 1:2)]

top3_CellCount_list2 <- lapply(1:4, function(x) {
    cell_sum <- sum(TCR$data[[x]]$Clones)
    lapply(
        top3_CellCount_list[[x]],
        function(y) {
            apply(y, 2, function(z) z / cell_sum * 10000)
        }
    )
})
names(top3_CellCount_list2) <- names(top3_CellCount_list)
top3_CellCount_list <- top3_CellCount_list2
top3_CellCount_list <- top3_CellCount_list[c(3:4, 1:2)]

top3_CellCount_df <- list()
for (i in 1:ncol(top3_aa_df)) {
    tmp <- do.call(rbind, lapply(
        1:3,
        function(x) {
            y <- reshape2::melt(top3_CellCount_list[[i]][[x]])
            y$Clone <- names(top3_CellCount_list[[i]])[x]
            colnames(y)[1:2] <- c("Var1", "Var2")
            y
        }
    ))
    tmp$Sample <- names(top3_CellCount_list)[i]
    tmp$Clone <- factor(tmp$Clone, levels = names(top3_CellCount_list[[i]]))
    top3_CellCount_df[[i]] <- tmp
}
top3_CellCount_df <- do.call(rbind, top3_CellCount_df)

top3_CellCount_df$group <- factor(
    paste(top3_CellCount_df$Var2, top3_CellCount_df$Clone, sep = "_"),
    levels = c(
        paste("CA", levels(top3_CellCount_df$Clone), sep = "_"),
        paste("Non-CA", levels(top3_CellCount_df$Clone), sep = "_")
    )
)
top3_CellCount_df$Var1 <- factor(
    top3_CellCount_df$Var1,
    levels = rev(levels(top3_CellCount_df$Var1))
)

top3_CellCount_df$Type <- top3_CellCount_df$Clone
levels(top3_CellCount_df$Type) <- sapply(
    levels(top3_CellCount_df$Type),
    function(x) top3_aaType_df[top3_aa_df == x]
)


p_bar_list <- list()
for (i in unique(top3_CellCount_df$Sample)[c(3, 1, 4, 2)]) {
    df_plot <- top3_CellCount_df[top3_CellCount_df$Sample == i, ]
    df_plot$group <- factor(
        df_plot$group,
        levels = levels(df_plot$group)[c(rbind(13:24, 1:12))]
    )
    df_plot_text <- df_plot[, c("Type", "group")]
    df_plot_text <- df_plot_text[!duplicated(df_plot_text), ]
    df_plot_text$value <- sapply(df_plot_text$group, function(x) max(df_plot$value[df_plot$group == x]))
    df_plot_text <- df_plot_text[df_plot_text$value > 0, ]

    max_value_idx <- df_plot$value %in% sort(df_plot_text$value, decreasing = T)[1:3]
    df_plot$Type <- ""
    df_plot$Type[max_value_idx] <- sapply(strsplit(
        as.character(df_plot_text$Type)[order(df_plot_text$value, decreasing = T)][1:3],
        split = ""
    ), function(x) x[1])

    geom_text(mapping = aes(label = Type, x = group, y = value), data = df_plot_text)

    p_bar_list[[i]] <- ggplot(df_plot, aes(x = group, y = value, fill = Var1)) +
        geom_col() +
        geom_text(aes(label = Type), vjust = 1) +
        theme_classic() +
        theme(
            legend.text = element_text(color = "black", size = 12),
            legend.position = "right",
            legend.title = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", size = 12),
            axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
            axis.ticks.x = element_line(color = "black", size = 0.5),
            axis.title.y = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.ticks.y = element_line(color = "black", size = 0.5)
        ) +
        guides(fill = guide_legend(title = "Cell Type", reverse = FALSE)) +
        scale_fill_manual(values = rev(colors_bar)) +
        labs(x = i, y = "Relative clone proportion (x 10^4)") +
        scale_y_continuous(expand = c(0.01, 0), limits = c(0, 100)) +
        scale_x_discrete(labels = rep(c("Non-CA", "CA"), 3))
}
p_bar_list[2:4] <- lapply(
    p_bar_list[2:4],
    function(x) {
        x + theme(
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    }
)
ggsave(
    filename = paste0("plot/scRNA5_TCR_top3Bar_byCluster_bySample.pdf"),
    plot = patchwork::wrap_plots(p_bar_list, ncol = 4, guides = "collect"),
    height = 6, width = 9
)


# shared TCR clones' matrix-plot
# by Sample
plot_shareTCR_mtx_ident <- function(
    ident, 
    ats = seq(0, 0.06, 0.02), 
    labels = seq(0, 6, 2), 
    cell_or_clone = "clone"
) {
    CA_1_meta <- T5_sub@meta.data[T5_sub$orig.ident == ident, ]
    shared_clone_mtx <- matrix(
        0,
        nrow = length(unique(CA_1_meta$IFN_Final)),
        ncol = length(unique(CA_1_meta$aa)),
        dimnames = list(
            unique(CA_1_meta$IFN_Final),
            unique(CA_1_meta$aa)
        )
    )
    tmp <- tapply(CA_1_meta$aa, CA_1_meta$IFN_Final, table)
    for (i in names(tmp)) {
        shared_clone_mtx[i, names(tmp[[i]])] <- tmp[[i]]
    }
    shared_clone_mtx <- shared_clone_mtx[, colnames(shared_clone_mtx) != "No"]
    Nshared_clone_mtx <- sapply(
        rownames(shared_clone_mtx),
        function(x) {
            sapply(
                rownames(shared_clone_mtx),
                function(y) {
                    if (cell_or_clone == "clone") {
                        z1 <- shared_clone_mtx[x, ] > 0
                        z2 <- shared_clone_mtx[y, ] > 0
                        sum(z1 & z2) / sum(z1 | z2)
                    } else {
                        z1 <- shared_clone_mtx[x, ]
                        z2 <- shared_clone_mtx[y, ]
                        idx <- z1 > 0 & z2 > 0
                        (sum(z1[idx]) + sum(z2[idx])) / (sum(z1) + sum(z2))
                    }
                }
            )
        }
    )
    cell_order <- c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    )
    Nshared_clone_mtx <- Nshared_clone_mtx[cell_order, cell_order]
    Nshared_clone_mtx[Nshared_clone_mtx == 1] <- min(
        1,
        max(Nshared_clone_mtx[Nshared_clone_mtx != max(Nshared_clone_mtx)]) * 2
    )

    p_shared <- ComplexHeatmap::Heatmap(
        Nshared_clone_mtx,
        name = "Fraction shared clones",
        col = grDevices::colorRampPalette(
            colors = c("#2f53a1", "white")
        )(50),
        rect_gp = gpar(col = "white", lwd = 0.5),
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_names = T,
        show_row_names = T,
        heatmap_legend_param = list(
            at = ats,
            direction = "horizontal",
            labels = paste0(labels, "%"),
            title_position = "topcenter",
            legend_width = unit(3, "in"),
            title_gp = gpar(fontsize = 12)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
                ifelse(
                    i == j,
                    1,
                    ifelse(
                        Nshared_clone_mtx[i, j] == 0,
                        0,
                        sprintf("%.3f", Nshared_clone_mtx[i, j])
                    )
                ),
                x,
                y,
                gp = gpar(fontsize = 10)
            )
        }
    )

    pdf(paste0("plot/scRNA5_TCR_sharedClone_byCluster_matrixplot_", ident, "_by", cell_or_clone, ".pdf"), width = 8, height = 9)
    draw(p_shared, heatmap_legend_side = "bottom")
    dev.off()
}
plot_shareTCR_mtx_ident("Cryo1wk", ats = c(0,0.1,1), labels = c(0,10,100))
plot_shareTCR_mtx_ident("Cryo2wk", ats = c(0,0.1,1), labels = c(0,10,100))
plot_shareTCR_mtx_ident("NonCryo1wk", ats = c(0,0.1,1), labels = c(0,10,100))
plot_shareTCR_mtx_ident("NonCryo2wk", ats = c(0,0.1,1), labels = c(0,10,100))

plot_shareTCR_mtx_ident("Cryo1wk", cell_or_clone = "cell", ats = c(0,0.1,1), labels = c(0,10,100))
plot_shareTCR_mtx_ident("Cryo2wk", cell_or_clone = "cell", ats = c(0,0.1,0.2,1), labels = c(0,10,20,100))
plot_shareTCR_mtx_ident("NonCryo1wk", cell_or_clone = "cell", ats = c(0,0.1,0.2,1), labels = c(0,10,20,100))
plot_shareTCR_mtx_ident("NonCryo2wk", cell_or_clone = "cell", ats = c(0,0.1,0.2,0.3,0.4,1), labels = c(0,10,20,30,40,100))


# by Time for Case-Control
plot_shareTCR_mtx_time <- function(
    time, 
    ats = seq(0, 0.06, 0.02), 
    labels = seq(0, 6, 2),
    cell_or_clone = "clone"
) {
    tmp_meta <- T5_sub@meta.data[T5_sub$orig.ident_time == time, ]
    idx <- grepl("NonCryo", tmp_meta$orig.ident)
    tmp_meta_case <- tmp_meta[!idx, ]
    tmp_meta_control <- tmp_meta[idx, ]

    shared_clone_mtx_1 <- matrix( # control
        0,
        nrow = length(unique(tmp_meta$IFN_Final)),
        ncol = length(unique(tmp_meta$aa)),
        dimnames = list(
            unique(tmp_meta$IFN_Final),
            unique(tmp_meta$aa)
        )
    )
    shared_clone_mtx_2 <- shared_clone_mtx_1 # case

    tmp1 <- tapply(tmp_meta_control$aa, tmp_meta_control$IFN_Final, table)
    for (i in names(tmp1)) {
        shared_clone_mtx_1[i, names(tmp1[[i]])] <- tmp1[[i]]
    }
    shared_clone_mtx_1 <- shared_clone_mtx_1[, colnames(shared_clone_mtx_1) != "No"]

    tmp2 <- tapply(tmp_meta_case$aa, tmp_meta_case$IFN_Final, table)
    for (i in names(tmp2)) {
        shared_clone_mtx_2[i, names(tmp2[[i]])] <- tmp2[[i]]
    }
    shared_clone_mtx_2 <- shared_clone_mtx_2[, colnames(shared_clone_mtx_2) != "No"]


    Nshared_clone_mtx <- sapply(
        rownames(shared_clone_mtx_1),
        function(x) {
            sapply(
                rownames(shared_clone_mtx_2),
                function(y) {
                    if (cell_or_clone == "clone") {
                        z1 <- shared_clone_mtx_1[x, ] > 0
                        z2 <- shared_clone_mtx_2[y, ] > 0
                        sum(z1 & z2) / sum(z1 | z2)
                    } else {
                        z1 <- shared_clone_mtx_1[x, ]
                        z2 <- shared_clone_mtx_2[y, ]
                        idx <- z1 > 0 & z2 > 0
                        (sum(z1[idx]) + sum(z2[idx])) / (sum(z1) + sum(z2))
                    }
                }
            )
        }
    )
    cell_order <- c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    )
    Nshared_clone_mtx <- Nshared_clone_mtx[cell_order, cell_order]

    p_shared <- ComplexHeatmap::Heatmap(
        Nshared_clone_mtx,
        name = "Fraction shared clones",
        col = grDevices::colorRampPalette(
            colors = c("#eff7d7", "#bbe913")
        )(50),
        rect_gp = gpar(col = "white", lwd = 0.5),
        cluster_columns = F,
        cluster_rows = F,
        show_row_dend = F,
        show_column_names = T,
        show_row_names = T,
        column_title = "Non-CA",
        row_title = "CA",
        heatmap_legend_param = list(
            at = ats,
            direction = "horizontal",
            labels = paste0(labels, "%"),
            title_position = "topcenter",
            legend_width = unit(3, "in"),
            title_gp = gpar(fontsize = 12)
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(
                ifelse(
                    Nshared_clone_mtx[i, j] == 0,
                    0,
                    sprintf("%.3f", Nshared_clone_mtx[i, j])
                ),
                x,
                y,
                gp = gpar(fontsize = 10)
            )
        }
    )

    pdf(paste0("plot/scRNA5_TCR_sharedClone_byCluster_matrixplot_", time, "_by", cell_or_clone, ".pdf"), width = 8, height = 9)
    draw(p_shared, heatmap_legend_side = "bottom")
    dev.off()
}
plot_shareTCR_mtx_time("1wk", ats = seq(0, 0.02, 0.005), labels = seq(0, 2, 0.5))
plot_shareTCR_mtx_time("2wk", ats = seq(0, 0.02, 0.005), labels = seq(0, 2, 0.5))

plot_shareTCR_mtx_time("1wk", cell_or_clone = "cell", ats = seq(0, 0.09, 0.03), labels = seq(0, 9, 3))
plot_shareTCR_mtx_time("2wk", cell_or_clone = "cell", ats = seq(0, 0.09, 0.03), labels = seq(0, 9, 3))

# Transition among Cluster by different samples
# too much edge
# only CD8+
compute_shareTCR_mtx_ident <- function(
    ident,
    cell_or_clone = "clone",
    cell_order = c(
        "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
        "T_others"
    )
) {
    CA_1_meta <- T5_sub@meta.data[T5_sub$orig.ident == ident, ]
    shared_clone_mtx <- matrix(
        0,
        nrow = length(unique(CA_1_meta$IFN_Final)),
        ncol = length(unique(CA_1_meta$aa)),
        dimnames = list(
            unique(CA_1_meta$IFN_Final),
            unique(CA_1_meta$aa)
        )
    )
    tmp <- tapply(CA_1_meta$aa, CA_1_meta$IFN_Final, table)
    for (i in names(tmp)) {
        shared_clone_mtx[i, names(tmp[[i]])] <- tmp[[i]]
    }
    shared_clone_mtx <- shared_clone_mtx[, colnames(shared_clone_mtx) != "No"]
    Nshared_clone_mtx <- sapply(
        rownames(shared_clone_mtx),
        function(x) {
            sapply(
                rownames(shared_clone_mtx),
                function(y) {
                    if (cell_or_clone == "clone") {
                        z1 <- shared_clone_mtx[x, ] > 0
                        z2 <- shared_clone_mtx[y, ] > 0
                        sum(z1 & z2) / sum(z1 | z2)
                    } else {
                        z1 <- shared_clone_mtx[x, ]
                        z2 <- shared_clone_mtx[y, ]
                        idx <- z1 > 0 & z2 > 0
                        (sum(z1[idx]) + sum(z2[idx])) / (sum(z1) + sum(z2))
                    }
                }
            )
        }
    )
    # cell_order <- c(
    #     "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
    #     "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
    #     "T_others"
    # )
    Nshared_clone_mtx <- Nshared_clone_mtx[cell_order, cell_order]
    Nshared_clone_mtx[Nshared_clone_mtx == 1] <- 0

    Nshared_clone_mtx
}
# write into .csv for later cytoscape's process 
write.csv(
    reshape2::melt(compute_shareTCR_mtx_ident('Cryo1wk', cell_order = c(
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
    ))),
    file = 'plot/scRNA5_TCR_Transition_byCluster_CA1wk.csv',
    quote = F
)
write.csv(
    reshape2::melt(compute_shareTCR_mtx_ident('Cryo2wk', cell_order = c(
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
    ))),
    file = 'plot/scRNA5_TCR_Transition_byCluster_CA2wk.csv',
    quote = F
)
write.csv(
    reshape2::melt(compute_shareTCR_mtx_ident('NonCryo1wk', cell_order = c(
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
    ))),
    file = 'plot/scRNA5_TCR_Transition_byCluster_NonCA1wk.csv',
    quote = F
)
write.csv(
    reshape2::melt(compute_shareTCR_mtx_ident('NonCryo2wk', cell_order = c(
        "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative"
    ))),
    file = 'plot/scRNA5_TCR_Transition_byCluster_NonCA2wk.csv',
    quote = F
)


####
Hig_list <- unique(Barcode_RNA$Clone[Barcode_RNA$prop > 0.005])
Hig_cell <- names(T5_sub$IFN_Final[T5_sub$clone_TCR %in% Hig_list])

df_high1 <- table(T5_sub@meta.data[Hig_cell, c("clone_TCR", "IFN_Final")])
df_high2 <- reshape2::melt(df_high1)
df_high3 <- reshape2::melt(apply(df_high1, 1, function(x) 100 * x / sum(x)))
df_high2$IFN_Final <- factor(
    df_high2$IFN_Final,
    levels = rev(names(colors_bar))[rev(names(colors_bar)) %in% df_high2$IFN_Final]
)
df_high2$clone_TCR <- factor(
    as.character(df_high2$clone_TCR),
    levels = rev(c(
        "C1_56_Non-CA-1wk", "C1_34_Non-CA-1wk", "C2_31_Non-CA-1wk", "C3_25_Non-CA-1wk", "C4_24_Non-CA-1wk",
        "C1_100_Non-CA-2wk", "C2_95_Non-CA-2wk", "C3_66_Non-CA-2wk", "C4_63_Non-CA-2wk",
        "C1_28_CA-1wk", "C2_27_CA-1wk", "C3_26_CA-1wk",
        "C1_70_CA-2wk", "C2_56_CA-2wk", "C3_42_CA-2wk"
    ))
)
levels(df_high2$clone_TCR)[11:15] <- c(
    "C5_24_Non-CA-1wk", "C4_25_Non-CA-1wk", "C3_31_Non-CA-1wk", "C2_34_Non-CA-1wk", "C1_56_Non-CA-1wk"
)
df_high3$IFN_Final <- factor(df_high3$IFN_Final, levels = rev(names(colors_bar))[rev(names(colors_bar)) %in% df_high3$IFN_Final])
df_high3$clone_TCR <- factor(
    as.character(df_high3$clone_TCR),
    levels = rev(c(
        "C1_56_Non-CA-1wk", "C1_34_Non-CA-1wk", "C2_31_Non-CA-1wk", "C3_25_Non-CA-1wk", "C4_24_Non-CA-1wk",
        "C1_100_Non-CA-2wk", "C2_95_Non-CA-2wk", "C3_66_Non-CA-2wk", "C4_63_Non-CA-2wk",
        "C1_28_CA-1wk", "C2_27_CA-1wk", "C3_26_CA-1wk",
        "C1_70_CA-2wk", "C2_56_CA-2wk", "C3_42_CA-2wk"
    ))
)
levels(df_high3$clone_TCR)[11:15] <- c(
    "C5_24_Non-CA-1wk", "C4_25_Non-CA-1wk", "C3_31_Non-CA-1wk", "C2_34_Non-CA-1wk", "C1_56_Non-CA-1wk"
)


p_top <- ggplot(df_high2, aes(x = clone_TCR, y = value, fill = IFN_Final)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("TCR clonotype") +
    ylab("Cell count") +
    labs(title = "Cell composition of High-Freq Clonotype", fill = "Cell type") +
    scale_fill_manual(values = colors_bar[levels(df_high2$IFN_Final)]) +
    theme_bw() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    guides(fill = guide_legend(reverse = T)) +
    coord_flip()

p_topp <- ggplot(df_high3, aes(x = clone_TCR, y = value, fill = IFN_Final)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("TCR clonotype") +
    ylab("Prop, %") +
    labs(title = "Cell composition of High-Freq Clonotype", fill = "Cell type") +
    scale_fill_manual(values = colors_bar[levels(df_high3$IFN_Final)]) +
    # theme_bw() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15)
    ) +
    guides(fill = guide_legend(reverse = T)) +
    coord_flip()

ggsave(p_top, filename = "plot/scRNA5_TCR_topClonotype_byCluster.pdf", width = 8, height = 4)
ggsave(p_topp, filename = "plot/scRNA5_TCR_topClonotype_prob_byCluster.pdf", width = 8, height = 4)
#####

#####
# plot common clonotype among clusters
pr.nt <- pubRep(TCR$data, "nt", .verbose = F, .coding = T)
pr.aa <- pubRep(TCR$data, "aa", .coding = T, .verbose = F)
pr.nt_aa <- pubRep(TCR$data, "nt+aa", .verbose = F, .coding = T)

pr.nt

data_x <- pr.nt
idxs <- 3:6
mtx <- matrix(
    0,
    nrow = 4, ncol = 4,
    dimnames = list(
        colnames(data_x)[idxs],
        colnames(data_x)[idxs]
    )
)
for (idx in 1:nrow(data_x)) {
    data_tmp <- data_x[idx, ..idxs]
    # data_tmp[is.na(data_tmp)] <- 0
    nSample <- data_x$Samples[idx]
    if (nSample == 1) {
        idx2 <- which(!is.na(data_tmp))
        mtx[idx2, idx2] <- mtx[idx2, idx2] + 1
    } else if (nSample == 2) {
        idx2 <- which(!is.na(data_tmp))
        mtx[idx2[1], idx2[2]] <- mtx[idx2[1], idx2[2]] + 1
        mtx[idx2[2], idx2[1]] <- mtx[idx2[2], idx2[1]] + 1
    } else {
        idx2 <- which(!is.na(data_tmp))
        for (idx21 in idx2) {
            for (idx22 in setdiff(idx2, idx21)) {
                mtx[idx21, idx22] <- mtx[idx21, idx22] + 1
            }
        }
    }
}
mtx
nrow(data_x)
mtx2 <- mtx
mtx2[c(2:4, 7:8, 12)] <- 0

# library(circlize)
# set.seed(3)
# chordDiagram(mtx, symmetric = F, keep.diagonal = T)
# circos.clear()

tmp <- reshape2::melt(mtx2)
tmp <- tmp[tmp$value != 0, ]
tmp$label <- seq(rownames(tmp))
colnames(tmp)[3] <- "value2"
tmp2 <- reshape2::melt(
    tmp,
    measure.vars = c("Var1", "Var2")
)
tmp2$value <- factor(
    tmp2$value,
    levels = c("Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk")
)
tmp2$variable <- factor(
    tmp2$variable,
    levels = rev(c("Var1", "Var2"))
)
# tmp2 <- tmp2[order(tmp2$value), ]
library(ggalluvial)
tmp2
p_TCRshared2 <- ggplot(
    tmp2,
    aes(
        x = variable,
        stratum = value,
        alluvium = label,
        y = value2,
        fill = value,
        label = value
    )
) +
    geom_flow(
        stat = "alluvium",
        lode.guidance = "frontback",
        color = c(
            "#66c2a5", "#66c2a5",
            "#fc8d62",
            "#8da0cb", "#8da0cb", "#8da0cb",
            "#e78ac3", "#e78ac3", "#e78ac3", "#e78ac3",
            rep("white", 10)
        ),
        linetype = 1
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    geom_stratum() +
    # geom_text(stat = "stratum", size = 3) +
    labs(title = "TCR clonotypes shared among Case-Control") +
    ylab("Clonotype number") +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.3, size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text = element_text(size = 12)
    )
p_TCRshared2$layers[[1]]$geom$default_aes$size <- 1
p_TCRshared2$layers[[1]]$geom$default_aesalpha <- 0.5
# factor(
#     c(
#         'Non-CA-1wk' = '#66c2a5',
#         'Non-CA-2wk' = '#fc8d62',
#         'CA-1wk' = '#e78ac3',
#         'CA-2wk' = '#8da0cb'
#     )[tmp2$value],
#     levels = c('#66c2a5','#fc8d62','#e78ac3','#8da0cb')
# )

ggsave(
    filename = "plot/scRNA5_TCR_sharedClonotype_byGroup.pdf",
    plot = p_TCRshared2,
    height = 6, width = 4
)

library(ggvenn)
data_venn <- as.list(as.data.frame(
    !t(apply(pr.nt[, c(4, 3, 5, 6)], 1, is.na))
))
data_venn2 <- lapply(data_venn, which)
p_venn <- ggvenn(
    data_venn2,
    show_elements = FALSE,
    show_percentage = TRUE,
    digits = 1,
    fill_color = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb"),
    fill_alpha = 0.7,
    stroke_color = "white",
    stroke_alpha = 1,
    stroke_size = 1,
    stroke_linetype = "solid",
    set_name_color = c("#66c2a5", "#fc8d62", "#e78ac3", "#8da0cb"),
    set_name_size = 5,
    text_color = "black",
    text_size = 3.5
) +
    scale_x_continuous(limits = c(-3, 3)) +
    labs(title = "TCR clonotypes shared among Case-Control") +
    theme(
        plot.title = element_text(size = 15, hjust = 0.5, vjust = -20)
    )

ggsave(
    filename = "plot/scRNA5_TCR_sharedClonotype_byGroup_venn.pdf",
    plot = p_venn,
    height = 6, width = 8
)

#####

#####
# # TCR expansion between clusters splited by Case-Control
# color_4 <- c('#e5e5e5', '#87cee7', '#ec7b8b', '#595757')
#
# head(Barcode_RNA)
# Barcode_RNA$Group <- sub("_[ACTGN]*-\\d$", "", rownames(Barcode_RNA))
# Barcode_RNA_list <- lapply(unique(Barcode_RNA$Group), function(x) Barcode_RNA[Barcode_RNA$Group == x, ])
# names(Barcode_RNA_list) <- unique(Barcode_RNA$Group)
# Barcode_RNA_list2 <- a
# Barcode_RNA_list[[1]]
#

#####
