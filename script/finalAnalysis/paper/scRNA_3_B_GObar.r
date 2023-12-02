source("/mnt/e/Cryo-TCR/server/auto/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis.rda")

setwd("~/Desktop/Cryo/Paper_Result_Plot/")
if (!dir.exists("plot")) dir.create("plot")

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

# 4.6 CA_vs_Non-CA with B cells
Idents(Cryo_merge) <- "MainCluster"
B <- subset(Cryo_merge, ident = "B")
Idents(B) <- "orig.ident2"
Markers_B_list <- FindAllMarkers(B, only.pos = T)
Markers_B_df <- my_Markers2df_multiple(
    Markers_B_list,
    logFC_threshold = 0.25,
    positive = T,
    n_top = 60
)

GO_B <- my_GO(
    Markers_B_df$Cluster_Cryo,
    return_plot = T, return_res = T,
    ont = "BP", Simplify = T, type = "bar", font.size = 18, show = 30
)

GO_B_res <- GO_B@result
dfB_for_plot <- makeTable(GO_B_res[1:25, ], 55) # filter

B_GO_bar <- ggplot(
    dfB_for_plot,
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
    ggtitle("CA_vs_Non-CA with B cells") +
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
B_GO_bar
ggsave(
    filename = "plot/scRNA3_B_GO_bar.pdf",
    plot = B_GO_bar,
    height = 6, width = 9
)
