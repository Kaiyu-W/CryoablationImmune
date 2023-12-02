source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/TCR_files/")

# load data
#####
if (file.exists("T5_TCR.rds")) {
    T5_sub <- readRDS(file = "T5_TCR.rds")
} else {
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
    meta_conga <- read.csv("bicluster_obs_df.csv", row.names = 1)
    X_umap <- read.csv("/home/kaiyu/Desktop/Cryo/T5_2106_2111_X_umap_3end.csv", row.names = 1)
    rownames(X_umap) <- rownames(meta_data)
    colnames(X_umap) <- c("UMAP_1", "UMAP_2")

    # filter cells
    Cryo_tmp <- merge(Cryo_2106, Cryo_2111)
    T5_sub <- Cryo_tmp[, rownames(Cryo_tmp@meta.data) %in% rownames(meta_conga)]
    T5_sub@meta.data <- cbind(
        meta_data[rownames(T5_sub@meta.data), ],
        meta_conga[rownames(T5_sub@meta.data), ]
    )
    T5_sub$IFN_Final[T5_sub$IFN_Final == "Tfh_I-IFN"] <- "CD4_T_I-IFN"
    T5_sub$orig.ident_time <- sub("^.*Cryo", "", T5_sub$orig.ident)
    T5_sub <- my_process_seurat(T5_sub, nVariableFeatures = 2500, norm.method = "SCT")
    # T5_sub@reductions$umap <- new('DimReduc')
    T5_sub@reductions$umap@cell.embeddings <- as.matrix(
        X_umap[colnames(T5_sub), ]
    )
    Idents(T5_sub) <- "IFN_Final"
    if (!file.exists("T5_TCR.rds")) saveRDS(T5_sub, file = "T5_TCR.rds")
}
T5_sub$bicluster <- paste(T5_sub$IFN_Final, T5_sub$TCR_clusters, sep = "-")
as.data.frame(sort(table(T5_sub$bicluster)))
#####

# load meta from python
#####
# R to plot
setwd("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/TCR_files/")
bicluster_count <- read.csv("bicluster_count_df.csv", row.names = 1)
bicluster_count_scale <- read.csv("bicluster_count_scale_df.csv", row.names = 1)
bicluster_count_raw <- read.csv("bicluster_count_raw_df.csv", row.names = 1)
bicluster_count_scale_raw <- read.csv("bicluster_count_scale_raw_df.csv", row.names = 1)
#####

# define function
#####
cluster_level <- c(
    "CD4_Naive", "CD4_T_I-IFN", "Th1", "Th17", "Treg", "Tfh",
    "CD8_Naive", "CD8_Tem_I-IFN", "CD8_Tem", "CD8_Tex", "CD8_Tex_Proliferative",
    "T_others"
)
my_CountCluster_clonesize <- function(
    cluster_object,
    group1 = "active.ident",
    group2 = "orig.ident",
    Heatmap = FALSE,
    trend_factor = 0,
    Percentage = TRUE,
    Round = TRUE) {
    # if trend_factor be not used, do not assign this variable
    # trend_factor means the other sum value and the present sum value
    # should be considered together.
    meta <- cluster_object@meta.data
    if (!all(c("clone_sizes", "conga_scores") %in% colnames(meta))) {
        stop("clone_sizes and conga_scores didn't exist!")
    }
    meta$active.ident <- cluster_object@active.ident
    df <- meta[, c(group1, group2)]
    fren <- table(df)
    for (x in colnames(fren)) {
        for (y in rownames(fren)) {
            idx <- meta[[group2]] == x & meta[[group1]] == y & meta[["conga_scores"]] <= 1
            fren[y, x] <- sum(meta[idx, "clone_sizes"])
        }
    }
    fren_perc <- as.matrix(fren)
    fren_perc <- apply(
        X = fren_perc,
        MARGIN = 2,
        FUN = function(x) x / sum(x) * (100**Percentage)
    )

    # if (Round)
    #     fren_perc <- round(fren_perc, ifelse(Percentage, 1, 3))
    res <- cbind(fren, fren_perc)
    if (ncol(fren_perc) == 2) {
        trend_factor_all <- ifelse(
            trend_factor == 0,
            1,
            trend_factor * sum(res[, 1]) / sum(res[, 2])
        )
        fren_prop <- round(fren_perc[, 1] / fren_perc[, 2] * trend_factor_all, 3)
        res <- cbind(res, fren_prop)
    }
    res <- rbind(res, apply(res, 2, sum))
    if (ncol(res) %% 2 == 1) {
        res[nrow(res), ncol(res)] <- NA
    }
    names(dimnames(res)) <- c(group1, group2)
    dimnames(res)[[1]][nrow(res)] <- "sum"
    dimnames(res)[[2]] <- c(
        paste0(colnames(fren), "_count"),
        paste0(colnames(fren), "_percentage"),
        "trend"
    )[1:ncol(res)]

    if (Round) {
        res[, 3:4] <- round(res[, 3:4], ifelse(Percentage, 1, 3))
    }
    if (Heatmap) {
        if (is.null(dim(res))) {
            message("result with ONE dimension cannot process heatmap!")
            return(res)
        }
        heatmap_temp <- res[-nrow(res), grep("_percentage", colnames(res), value = TRUE)]
        print(
            pheatmap::pheatmap(
                heatmap_temp,
                angle_col = 45,
                cluster_rows = FALSE,
                cluster_cols = FALSE
            )
        )
    }
    return(res)
}
myMelt <- function(df) {
    df_bicluster <- reshape2::melt(
        df,
        measure.vars = c("NonCryo1wk", "NonCryo2wk", "Cryo1wk", "Cryo2wk")
    )
    df_bicluster$x <- paste(
        df_bicluster$level_0, df_bicluster$level_1,
        sep = "_"
    )
    levels(df_bicluster$variable) <- c(
        "Non-CA-1wk", "Non-CA-2wk", "CA-1wk", "CA-2wk"
    )
    df_bicluster$level_0 <- factor(
        df_bicluster$level_0,
        levels = cluster_level
    )
    df_bicluster$value <- as.numeric(df_bicluster$value)
    df_bicluster
}
myMelt_combine <- function(df, by = c("time", "group")[1]) {
    df_bicluster <- reshape2::melt(
        df,
        measure.vars = c("NonCryo1wk", "NonCryo2wk", "Cryo1wk", "Cryo2wk")
    )
    df_bicluster$x <- paste(
        df_bicluster$level_0, df_bicluster$level_1,
        sep = "_"
    )
    df_bicluster$variable <- sapply(
        df_bicluster$variable,
        function(x) {
            ifelse(
                by == "time",
                ifelse(grepl("1wk", x), "1wk", "2wk"),
                ifelse(grepl("Non", x), "Non-CA", "CA")
            )
        }
    )
    df_bicluster$variable <- factor(
        df_bicluster$variable,
        levels = if (by == "time") c("1wk", "2wk") else c("Non-CA", "CA")
    )
    df_bicluster$level_0 <- factor(
        df_bicluster$level_0,
        levels = cluster_level
    )
    df_bicluster$value <- as.numeric(df_bicluster$value)
    df_bicluster
}
pbar <- function(df, y_title = "Cell count", x_title = "TCR cluster", by = NA, percent = F) {
    color_bar <- if (is.na(by)) {
        c(
            "Non-CA-1wk" = "#66c2a5",
            "Non-CA-2wk" = "#fc8d62",
            "CA-1wk" = "#e78ac3",
            "CA-2wk" = "#8da0cb"
        )
    } else if (by == "time") {
        c(
            "1wk" = "#66c2a5",
            "2wk" = "#fc8d62"
        )
    } else if (by == "group") {
        c(
            "Non-CA" = "#e78ac3",
            "CA" = "#8da0cb"
        )
    }
    df$value <- as.numeric(df$value)
    if (percent) df$value <- df$value * 100
    ggplot(df, aes(x = level_1, y = value, fill = variable)) +
        geom_bar(width = 1, stat = "identity") +
        facet_grid(. ~ level_0, scales = "free") +
        ylab(y_title) +
        xlab(x_title) +
        labs(title = "Group composition of each bi-cluster of significant CoNGA score", fill = "Group") +
        scale_fill_manual(
            values = rev(color_bar)
        ) +
        theme_bw() +
        theme(
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            axis.text.x = element_text(
                size = 6, angle = 45, hjust = 1, vjust = 1
            ),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            plot.title = element_text(hjust = 0.5, size = 15)
        ) +
        guides(fill = guide_legend(reverse = T)) +
        scale_y_continuous(expand = c(0.01, 0.01))
}
#####

# plot
#####
setwd("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/")
# T celltype count
#####
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
T_sub_count_sum <- as.data.frame(my_CountCluster_clonesize(T5_sub, group1 = "IFN_Final", group2 = "orig.ident")[-13, 1:4])
T_sub_count_sum2 <- my_CountCluster_clonesize(T5_sub, group1 = "IFN_Final", group2 = "orig.ident")["sum", 1:4]
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
    filename = "plot/scRNA5_T_count_cellbyconga.pdf",
    plot = T_sub_count_Plot,
    height = 6.2, width = 6
)

T_sub_count_sum <- as.data.frame(my_CountCluster(T5_sub[, T5_sub$conga_scores <= 1], group1 = "IFN_Final", group2 = "orig.ident")[-13, 1:4])
T_sub_count_sum2 <- my_CountCluster(T5_sub[, T5_sub$conga_scores <= 1], group1 = "IFN_Final", group2 = "orig.ident")["sum", 1:4]
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
    filename = "plot/scRNA5_T_count_clonebyconga.pdf",
    plot = T_sub_count_Plot,
    height = 6.2, width = 6
)


# # all T:
# NonCryo2wk    17888
# Cryo1wk       12752
# Cryo2wk       11432
# NonCryo1wk    11088
# # TCR T:
# Cryo1wk       6552
# NonCryo1wk    6315
# NonCryo2wk    5671
# Cryo2wk       5201
# # conga T:
# my_CountCluster_clonesize(T5_sub, group1 = 'IFN_Final', group2 = 'orig.ident')['sum',1:4]
# Cryo1wk_count    Cryo2wk_count NonCryo1wk_count NonCryo2wk_count
# 2843             2257             2863             2306
df_for_bar <- matrix(
    c(
        11088, 17888, 12752, 11432,
        6315, 5671, 6552, 5201,
        2863, 2306, 2843, 2257
    ),
    nrow = 4,
    dimnames = list(
        c("Non-CA\n1wk", "Non-CA\n2wk", "CA\n1wk", "CA\n2wk"),
        c("All", "TCR exist", "CoNGA score <= 1")
    )
)
df_for_bar_minus <- t(apply(
    df_for_bar,
    1,
    function(x) {
        c(x[1] - x[2], x[2] - x[3], x[3])
    }
))
df_for_bar_2 <- reshape2::melt(
    df_for_bar_minus,
    measure.vars = rownames(df_for_bar)
)
df_for_bar_2$Var1 <- factor(
    df_for_bar_2$Var1,
    levels = rev(c("Non-CA\n1wk", "Non-CA\n2wk", "CA\n1wk", "CA\n2wk"))
)
p_countbar <- ggplot(df_for_bar_2, aes(x = Var1, y = value, fill = Var2)) +
    geom_bar(width = 0.9, stat = "identity") +
    ylab("Number of Cells") +
    theme_bw() +
    coord_flip() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = "#363a38", size = 8),
        axis.text.y = element_text(color = "#363a38", size = 8),
        axis.line = element_blank(),
        axis.ticks.x = element_line(color = "#363a38", size = 0.5, lineend = "square"),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
        panel.border = element_rect(color = "#363a38", size = 1, linetype = 4)
    ) +
    guides(fill = guide_legend(reverse = FALSE)) +
    scale_fill_manual(
        values = rev(c(
            "All" = "#595757",
            "TCR exist" = "#87cee7",
            "CoNGA score <= 1" = "#ec7b8b"
        ))
    ) +
    scale_y_continuous(
        position = "right",
        expand = expansion(c(0, 0.05)),
        limits = c(0, 18000),
        breaks = seq(0, 20000, 3000)
    )

ggsave(
    filename = "plot/scRNA5_T_count_cellfilterbygroup.pdf",
    plot = p_countbar,
    height = 2.5, width = 6
)
#####
# TCR cluster count
#####
df_bicluster <- myMelt(bicluster_count)
for (i in as.character(unique(df_bicluster$variable))) {
    for (j in unique(df_bicluster$level_1)) {
        for (k in as.character(unique(df_bicluster$level_0))) {
            if (!any(df_bicluster$variable == i & df_bicluster$level_1 == j & df_bicluster$level_0 == k)) {
                df_bicluster <- rbind(df_bicluster, c(k, j, i, 0, paste(k, j, sep = "_")))
            }
        }
    }
}
p_bicluster <- pbar(df_bicluster)
ggsave(
    filename = "plot/scRNA5_T_cellcount_bycluster_byconga.pdf",
    plot = p_bicluster,
    height = 8, width = 16
)

df_bicluster_scale <- myMelt(bicluster_count_scale)
for (i in as.character(unique(df_bicluster_scale$variable))) {
    for (j in unique(df_bicluster_scale$level_1)) {
        for (k in as.character(unique(df_bicluster_scale$level_0))) {
            if (!any(df_bicluster_scale$variable == i & df_bicluster_scale$level_1 == j & df_bicluster_scale$level_0 == k)) {
                df_bicluster_scale <- rbind(df_bicluster_scale, c(k, j, i, 0, paste(k, j, sep = "_")))
            }
        }
    }
}
p_bicluster_scale <- pbar(df_bicluster_scale, y_title = "Cell Percent (%)", percent = T)
ggsave(
    filename = "plot/scRNA5_T_cellratio_bycluster_byconga.pdf",
    plot = p_bicluster_scale,
    height = 8, width = 16
)
#####
# pie plot
#####
T5_sub$conga_sig <- T5_sub$conga_scores <= 1 * T5_sub$clone_sizes
df_conga_sig <- lapply(
    unique(T5_sub$orig.ident),
    function(x) {
        xx <- tapply(
            T5_sub@meta.data[T5_sub$orig.ident == x, 'conga_sig'],
            T5_sub@meta.data[T5_sub$orig.ident == x, 'IFN_Final'],
            function(y) c(sum(y) / length(y), 1 - sum(y) / length(y))
        )
        xx2 <- as.data.frame(lapply(xx, function(x) x), check.names = F)
        rownames(xx2) <- c('sig', 'no-sig')
        xx3 <- reshape2::melt(t(xx2))
        xx3$sample <- x
        xx3
    }
)
names(df_conga_sig) <- unique(T5_sub$orig.ident)
df_conga_sig <- do.call(rbind, df_conga_sig)
rownames(df_conga_sig) <- 1:nrow(df_conga_sig)
df_conga_sig$sample <- factor(df_conga_sig$sample, levels = rev(c('NonCryo1wk','NonCryo2wk','Cryo1wk','Cryo2wk')))
levels(df_conga_sig$sample) <- rev(c('Non-CA\n1wk','Non-CA\n2wk','CA\n1wk','CA\n2wk'))
df_conga_sig$x <- paste(df_conga_sig$sample, df_conga_sig$Var1, sep = "-")

p_pie0 <- ggplot(
    df_conga_sig,
    aes(y = sample, x = Var1, group = x)
) +
    geom_jjPointPie(
        aes(pievar = value, fill = Var2),
        width = 0.8,
        color = "white",
        line.size = 0.5
    ) +
    # coord_fixed(clip = 'off') +
    xlab("GEX clsuter") +
    ylab("Group") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
    ) +
    scale_edge_fill_manual(
        values = c(
            "sig" = "#f8766d",
            "no-sig" = "#bebebe"
        )
    )
p_pie0
ggsave(
    filename = "plot/scRNA5_T_cellcount_bySig_byconga.pdf",
    plot = p_pie0,
    height = 5, width = 6
)

#####
tcr_clusters <- c("0", "1_av13", "2_AV6", "3_AV12", "4_bv13", "5_AV9", "6_AV16", "7_BV1", "8_AV8", "9_BV20", "10_AV7")
tcr_colorbars_tab10 <- c(
    "#1f77b4",
    "#ff7f0e",
    "#279e68",
    "#d62728",
    "#aa40fc",
    "#8c564b",
    "#e377c2",
    "#b5bd61",
    "#17becf",
    "#aec7e8",
    "#ffbb78"
)
tcr_colorbars_tab20 <- c(
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5"
)
tcr_colorbars <- tcr_colorbars_tab20[1:11]
names(tcr_colorbars) <- tcr_clusters
library(jjPlot)
df_bicluster_scale <- myMelt(bicluster_count_scale)
df_bicluster_scale$level_1 <- factor(
    df_bicluster_scale$level_1,
    levels = tcr_clusters
)
p_pie <- ggplot(
    df_bicluster_scale,
    aes(x = level_0, y = level_1, group = x)
) +
    geom_jjPointPie(
        aes(pievar = value, fill = variable),
        width = 0.8,
        color = "white",
        line.size = 0.5
    ) +
    # coord_fixed(clip = 'off') +
    xlab("GEX clsuter") +
    ylab("TCR cluster") +
    theme_bw() +
    theme(
        legend.text = element_text(color = "black", size = 8),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.title.x = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
    ) +
    scale_edge_fill_manual(
        values = c(
            # 'Non-CA-1wk' = '#66c2a5',
            # 'Non-CA-2wk' = '#fc8d62',
            # 'CA-1wk' = '#e78ac3',
            # 'CA-2wk' = '#8da0cb'
            "Non-CA-1wk" = "#f8766d",
            "Non-CA-2wk" = "#bed780",
            "CA-1wk" = "#80dfe2",
            "CA-2wk" = "#e3beff"
        )
    )
# + guides(
#     color = guide_legend(
#         override.aes = list(
#             shape = 21
#         )
#     )
# )
p_pie
gex_sum <- reshape2::melt(
    tapply(
        apply(bicluster_count[3:6], 1, sum),
        bicluster_count$level_0,
        sum
    )
)
gex_sum$Var1 <- factor(gex_sum$Var1, levels = cluster_level)
tcr_sum <- reshape2::melt(
    tapply(
        apply(bicluster_count[3:6], 1, sum),
        bicluster_count$level_1,
        sum
    )
)
tcr_sum$Var1 <- factor(tcr_sum$Var1, levels = tcr_clusters)

p_pietop <- ggplot(gex_sum, aes(x = Var1, y = value, fill = Var1)) +
    geom_bar(stat = "identity") +
    # coord_fixed(clip = 'off') +
    ylab("Number of Cells") +
    theme_classic(base_size = 11 / 2) +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        # axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
    ) +
    scale_fill_manual(values = colors_bar) +
    scale_y_continuous(
        expand = expansion(c(0, 0.05))
    )

p_pieright <- ggplot(tcr_sum, aes(x = Var1, y = value, fill = Var1)) +
    geom_bar(stat = "identity") +
    # coord_fixed(clip = 'off') +
    ylab("Number of Cells") +
    theme_classic(base_size = 11 / 2) +
    theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(color = "black", size = 8),
        # axis.text.y = element_text(color = 'black', size=8),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.background = element_blank(),
    ) +
    coord_flip() +
    scale_fill_manual(values = tcr_colorbars) +
    scale_y_continuous(
        expand = expansion(c(0, 0.05))
    )


p_pie_all <- patchwork::wrap_plots(
    list(
        p_pietop, ggplot() +
            cowplot::theme_nothing(),
        p_pie, p_pieright
    ),
    ncol = 2, nrow = 2,
    widths = c(0.8, 0.2), heights = c(0.2, 0.8)
)

ggsave(
    filename = "plot/scRNA5_T_cellcount_bybicluster_byconga.pdf",
    plot = p_pie_all,
    height = 8, width = 8
)
#####

# GO for Naive cells
#####
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
bicluster_count_naive <- bicluster_count[grepl("Naive", bicluster_count$level_0), ]

T5_sub$"orig.ident_group" <- sub("\\dwk", "", T5_sub$orig.ident)
T5_sub_sig <- subset(
    T5_sub,
    cell = colnames(T5_sub)[T5_sub$conga_scores <= 1 & T5_sub$orig.ident_time == "1wk"]
)
T5_sub_sig_CD4Naive <- subset(
    T5_sub_sig,
    cell = colnames(T5_sub_sig)[T5_sub_sig$IFN_Final == "CD4_Naive"]
)
T5_sub_sig_CD8Naive <- subset(
    T5_sub_sig,
    cell = colnames(T5_sub_sig)[T5_sub_sig$IFN_Final == "CD8_Naive"]
)
T5_sub_sig_CD4Naive
T5_sub_sig_CD8Naive
table(T5_sub_sig_CD4Naive$orig.ident_group)
table(T5_sub_sig_CD8Naive$orig.ident_group)

# VlnPlot(T5_sub_sig, features = c('Cd3e','Cd3d','Cd4','Cd8a','Cd8b1'))
# FeaturePlot(T5_sub_sig, features = c('Cd69','Cd27'),pt.size = 2)
# FeaturePlot(T5_sub_sig, features = c('Cd69','Cd44'),pt.size = 2)

# CD4 Naive
Idents(T5_sub_sig_CD4Naive) <- "orig.ident_group"
Markers_CD4Naive_list <- FindAllMarkers(T5_sub_sig_CD4Naive, only.pos = F, logfc.threshold = 0)
Markers_CD4Naive_df <- my_Markers2df_multiple(
    Markers_CD4Naive_list[Markers_CD4Naive_list$p_val_adj < 0.05, ],
    logFC_threshold = 0.15,
    positive = T,
    # n_top = 40
    n_top = 20
)
# my_plotVolcano(
#     Markers_CD4Naive_list,
#     case_group = 'Cryo',
#     control_group = 'NonCryo',
#     lfc_thre = 0.15,
#     p_thre = 0.05
# )
# my_DotPlot_split(
#     T5_sub_sig_CD4Naive,
#     feature = lapply(as.list(Markers_CD4Naive_df), na.omit)
# )
# Seurat::DoHeatmap(
#     T5_sub_sig_CD4Naive,
#     feature = unlist(lapply(as.list(Markers_CD4Naive_df), na.omit)),
#     slot = 'scale.data'
# )

GO_CD4Naive <- my_GO(
    Markers_CD4Naive_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_CD4Naive_res <- GO_CD4Naive@result
df_CD4Naive_for_plot <- makeTable(GO_CD4Naive_res[1:30, ], 70) # filter

CD4Naive_GO_bar <- ggplot(
    df_CD4Naive_for_plot,
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
    ggtitle("CA_vs_Non-CA with CD4 Naive cells with significant CoNGA score") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 1.18),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
CD4Naive_GO_bar
ggsave(
    filename = "plot/scRNA5_CD4Naive_byconga_GO_bar.pdf",
    plot = CD4Naive_GO_bar,
    height = 6 * 1.2, width = 9
)

# CD8 Naive
Idents(T5_sub_sig_CD8Naive) <- "orig.ident_group"
Markers_CD8Naive_list <- FindAllMarkers(T5_sub_sig_CD8Naive, only.pos = T)
Markers_CD8Naive_df <- my_Markers2df_multiple(
    Markers_CD8Naive_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 40
)

GO_CD8Naive <- my_GO(
    Markers_CD8Naive_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_CD8Naive_res <- GO_CD8Naive@result
df_CD8Naive_for_plot <- makeTable(GO_CD8Naive_res[1:30, ], 70) # filter

CD8Naive_GO_bar <- ggplot(
    df_CD8Naive_for_plot,
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
    ggtitle("CA_vs_Non-CA with CD8 Naive cells with significant CoNGA score") +
    theme_classic() +
    xlab("-log10 FDR") +
    theme(
        plot.title = element_text(size = 18, color = "black", hjust = 1.18),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", lineheight = 0.6),
        axis.ticks.length.y = unit(5, "pt"),
        # plot.margin=unit(c(rep(0.5,3),1),'cm')
    )
CD8Naive_GO_bar
ggsave(
    filename = "plot/scRNA5_CD8Naive_byconga_GO_bar.pdf",
    plot = CD8Naive_GO_bar,
    height = 6 * 1.2, width = 9
)
#####

# # activation signature
# source("~/Desktop/Github/My_scRNA_pipeline/Pathway_score.r")
# act_sig_genes <- list(
#     activation_signature = c("Ifng", "Ccl4", "Ccl3", "Il3", "Xcl1", "Csf2", "Gzmb", "Fabp5", "Xcl2", "Lta", "Lag3", "Tnfrsf4", "Tnfrsf9", "Pim3", "Mir155hg")
# )
# T_tmp <- geneset_score(
#     T5_sub,
#     act_sig_genes,
#     gene_upper = F, # method = 'GSVA'
# )
# T_tmp$"orig.ident_group" <- sub("\\dwk", "", T_tmp$orig.ident)

# FeaturePlot(T_tmp, "activation_signature")
# VlnPlot(T_tmp, "activation_signature")
# VlnPlot(T_tmp, "activation_signature", split.by = "orig.ident_group", pt.size = 0)

# wilcox.test(
#     T_tmp$activation_signature[T_tmp$IFN_Final == "CD4_Naive" & T_tmp$orig.ident_group == "Cryo"],
#     T_tmp$activation_signature[T_tmp$IFN_Final == "CD4_Naive" & T_tmp$orig.ident_group == "NonCryo"],
#     alternative = "greater"
# )
# wilcox.test(
#     T_tmp$activation_signature[T_tmp$IFN_Final == "CD8_Naive" & T_tmp$orig.ident_group == "Cryo"],
#     T_tmp$activation_signature[T_tmp$IFN_Final == "CD8_Naive" & T_tmp$orig.ident_group == "NonCryo"],
#     alternative = "greater"
# )


# rds_file <- "~/Desktop/Cryo/Paper_Result_Plot/T_sub_all.rds"
# if (!file.exists(rds_file)) {
#     #####
#     # load T5_2111
#     # get matrix from server
#     datas_2111 <- paste0("/home/kaiyu/Desktop/Cryo/2111_data/", c("Cryo1wk", "Cryo2wk", "NonCryo1wk", "NonCryo2wk"), "_filtered_feature_bc_matrix.h5")
#     CA1_2111_mtx <- Seurat::Read10X_h5(datas_2111[1])
#     CA2_2111_mtx <- Seurat::Read10X_h5(datas_2111[2])
#     NonCA1_2111_mtx <- Seurat::Read10X_h5(datas_2111[3])
#     NonCA2_2111_mtx <- Seurat::Read10X_h5(datas_2111[4])
#     colnames(CA1_2111_mtx) <- paste0("Cryo1wk2111_", colnames(CA1_2111_mtx))
#     colnames(CA2_2111_mtx) <- paste0("Cryo2wk2111_", colnames(CA2_2111_mtx))
#     colnames(NonCA1_2111_mtx) <- paste0("NonCryo1wk2111_", colnames(NonCA1_2111_mtx))
#     colnames(NonCA2_2111_mtx) <- paste0("NonCryo2wk2111_", colnames(NonCA2_2111_mtx))

#     CA1_2111 <- Seurat::CreateSeuratObject(CA1_2111_mtx, min.cells = 0, min.features = 0)
#     CA2_2111 <- Seurat::CreateSeuratObject(CA2_2111_mtx, min.cells = 0, min.features = 0)
#     NonCA1_2111 <- Seurat::CreateSeuratObject(NonCA1_2111_mtx, min.cells = 0, min.features = 0)
#     NonCA2_2111 <- Seurat::CreateSeuratObject(NonCA2_2111_mtx, min.cells = 0, min.features = 0)
#     Cryo_2111 <- merge(CA1_2111, c(CA2_2111, NonCA1_2111, NonCA2_2111))

#     # load T5_2106
#     # get matrix from server
#     datas_2106 <- paste0("/home/kaiyu/Desktop/Cryo/2106_data/", c("Cryo", "NonCryo"), "_filtered_feature_bc_matrix.h5")
#     CA_2106_mtx <- Seurat::Read10X_h5(datas_2106[1])
#     NonCA_2106_mtx <- Seurat::Read10X_h5(datas_2106[2])
#     colnames(CA_2106_mtx) <- paste0("Cryo1wk2106_", colnames(CA_2106_mtx))
#     colnames(NonCA_2106_mtx) <- paste0("NonCryo1wk2106_", colnames(NonCA_2106_mtx))
#     CA_2106 <- Seurat::CreateSeuratObject(CA_2106_mtx, min.cells = 0, min.features = 0)
#     NonCA_2106 <- Seurat::CreateSeuratObject(NonCA_2106_mtx, min.cells = 0, min.features = 0)
#     Cryo_2106 <- merge(CA_2106, NonCA_2106)

#     # load meta and umap
#     meta_data <- read.csv("/home/kaiyu/Desktop/Cryo/T5_2106_2111_obs.csv", row.names = 1)
#     X_umap <- read.csv("/home/kaiyu/Desktop/Cryo/T5_2106_2111_X_umap_3end.csv", row.names = 1)
#     rownames(X_umap) <- rownames(meta_data)
#     colnames(X_umap) <- c("UMAP_1", "UMAP_2")

#     # filter cells
#     Cryo_tmp <- merge(Cryo_2106, Cryo_2111)
#     T5_sub_all <- Cryo_tmp[, rownames(Cryo_tmp@meta.data) %in% rownames(meta_data)]
#     T5_sub_all@meta.data <- meta_data[rownames(T5_sub_all@meta.data), ]
#     T5_sub_all$IFN_Final[T5_sub_all$IFN_Final == "Tfh_I-IFN"] <- "CD4_T_I-IFN"
#     T5_sub_all$orig.ident_time <- sub("^.*Cryo", "", T5_sub_all$orig.ident)

#     T5_sub_all <- my_process_seurat(T5_sub_all, nVariableFeatures = 2500, norm.method = "SCT")
#     T5_sub_all@reductions$umap@cell.embeddings <- as.matrix(
#         X_umap[rownames(T5_sub_all@reductions$umap@cell.embeddings), ]
#     )
#     Idents(T5_sub_all) <- "IFN_Final"

#     #####
#     saveRDS(T5_sub_all, rds_file)
# } else {
#     T5_sub_all <- readRDS(rds_file)
# }


# T5_sub_all$tmp <- paste0(T5_sub_all$IFN_Final, "_noconga")
# Idents(T5_sub_all) <- "tmp"
# T5_sub_all <- my_AddMeta(T5_sub_all, T_tmp$IFN_Final, T)
# DimPlot(T5_sub_all, group.by = "T_tmp_IFN_Final")
# Idents(T5_sub_all) <- "T_tmp_IFN_Final"
# T_tmp2 <- subset(T5_sub_all, ident = c("CD4_Naive_noconga", "CD4_Naive", "CD8_Naive_noconga", "CD8_Naive"))

# T_tmp2 <- geneset_score(
#     T_tmp2,
#     act_sig_genes,
#     gene_upper = F, # method = 'GSVA'
# )

# VlnPlot(T_tmp2, "activation_signature", split.by = "T_tmp_IFN_Final")

# wilcox.test(
#     T_tmp2$activation_signature[T_tmp2$T_tmp_IFN_Final == "CD4_Naive"],
#     T_tmp2$activation_signature[T_tmp2$T_tmp_IFN_Final == "CD4_Naive_noconga"],
#     alternative = "greater"
# )
# wilcox.test(
#     T_tmp2$activation_signature[T_tmp2$T_tmp_IFN_Final == "CD8_Naive"],
#     T_tmp2$activation_signature[T_tmp2$T_tmp_IFN_Final == "CD8_Naive_noconga"],
#     alternative = "greater"
# )


# T5_sub_all <- geneset_score(
#     T5_sub_all,
#     act_sig_genes,
#     gene_upper = F, # method = 'GSVA'
# )

# VlnPlot(T5_sub_all, "activation_signature", split.by = "T_tmp_IFN_Final") + theme(legend.position = "none")
