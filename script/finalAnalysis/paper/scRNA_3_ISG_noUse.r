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

if (!require(AUCell)) {
    stop("No AUCell!")
}

getDEGs <- function(df, logFC_threshold = 0.25) {
    df <- subset(df, subset = avg_log2FC >= logFC_threshold)
    genes <- sub(".\\d$", "", rownames(df))
    unique(genes)
}
color_bar <- c("#050415", "#e34e64", "#fcf5b7")

# CD45+ differential genes
Idents(Cryo_merge) <- "orig.ident2"
Markers_CD45plus_list <- FindAllMarkers(Cryo_merge, only.pos = F)
Markers_CD45plus <- getDEGs(
    Markers_CD45plus_list,
    logFC_threshold = 0.25
)

# T differential genes
Idents(T_sub) <- "orig.ident2"
Markers_T_list <- FindAllMarkers(T_sub, only.pos = F)
Markers_T <- getDEGs(
    Markers_T_list,
    logFC_threshold = 0.25
)

cat(Markers_CD45plus, sep = ",") # input into interferome database
cat(Markers_T, sep = ",") # input into interferome database
# http://interferome.its.monash.edu.au/interferome/home.jspx
# get genelist from interferome
interferome_list <- c("Akap13", "Alox5ap", "Atf3", "Basp1", "Bcl2a1b", "Bcl2a1d", "Bhlhe40", "C1qb", "Ccl2", "Ccl3", "Ccl4", "Ccl6", "Ccrl2", "Cd14", "Cd3g", "Cd44", "Cd63", "Cd79b", "Cd9", "Cd9-ps", "Cdkn1a", "Cebpb", "Clec4e", "Cma1", "Cst3", "Ctla2a", "Cxcl2", "Cxcl3", "Cxcr6", "Dnaja1", "Dnajb1", "Dusp1", "Dusp5", "Eef1g", "Egr1", "Fabp5", "Fcer1g", "Fcrla", "Fosb", "Gadd45b", "Gm15470", "Gm7816", "Gm9089", "Gsn", "Gsr", "Gzma", "H2-Aa", "H2-DMa", "H2-Ob", "H3f3a", "Hdc", "Hist1h1c", "Hopx", "Hspa1a", "Hspa1b", "Hspa5", "Hsph1", "Icam1", "Icos", "Ier3", "Ifi27l2a", "Ifrd1", "Ighm", "Ikzf2", "Il1a", "Il1b", "Il1r2", "Il1rl1", "Il1rn", "Itm2b", "Jun", "Junb", "Kdm6b", "Klf6", "Lmnb1", "Ly6a", "Maf", "Marcks", "Marcksl1", "Mef2c", "Mmp9", "Mpp7", "Ms4a4c", "Mt1", "Napsa", "Ncf1", "Nfkbia", "Nfkbiz", "Nlrp3", "Nr3c1", "Nr4a3", "Odc1", "Phlda1", "Pim1", "Plek", "Ppp1r15a", "Ptgs2", "Rab8b", "Ramp3", "Rbpj", "Rgs1", "Rgs16", "Rgs2", "Rora", "Rpl10a", "Rpl13a", "Rpl35a", "Rps3", "Rpsa", "Samsn1", "Sell", "Slpi", "Srgn", "Tbc1d4", "Thbs1", "Tmem176a", "Tmem176b", "Tmsb4x", "Tnf", "Tnfaip2", "Tnfaip3", "Tnfaip8", "Trim30a", "Vegfa", "Xcl1", "Zfp36")
interferome_list2 <- c("Aebp2", "Ahnak", "Akap13", "Arhgap31", "Atf3", "AW112010", "Bcl2a1b", "Bcl2a1d", "Bhlhe40", "Btg2", "Camk2n1", "Capg", "Ccr5", "Cd3g", "Cd44", "Cd5", "Cd52", "Crem", "Cst3", "Cxcl2", "Cxcr4", "Cxcr6", "Dusp1", "Eef1g", "Fnbp1", "Fos", "Fosb", "Gadd45b", "Gm13160", "Gm13237", "Gm15470", "Got1", "Gzma", "Gzmb", "Herpud1", "Hist1h1c", "Hmgb2", "Hopx", "Hspa1a", "Hspa1b", "Hspa5", "Hsph1", "Icos", "Ifi27l2a", "Ifngr1", "Ikzf2", "Il18r1", "Il2ra", "Itga4", "Itgae", "Itm2b", "Jun", "Kdm6b", "Klf6", "Lgals3bp", "Lmnb1", "Lyz2", "Maf", "Malat1", "Malt1", "Mpp7", "Nabp1", "Nfkbia", "Nr3c1", "Odc1", "Pdcd1", "Phlda1", "Pim1", "Ptprcap", "Rab8b", "Rabgap1l", "Rbpj", "Rel", "Rgs1", "Rgs16", "Rgs2", "Rora", "Rpl10a", "Rpl13a", "Rpl35a", "Rps27a", "Rps3", "Rpsa", "Runx1", "S100a4", "Samsn1", "Sell", "Sla", "Tbc1d4", "Tmem176a", "Tmem176b", "Tmsb4x", "Tnfaip3", "Traf1", "Trps1", "Ugcg", "Xist", "Zfp36l1")
interferome_list <- interferome_list[interferome_list %in% Markers_CD45plus]
interferome_list2 <- interferome_list2[interferome_list2 %in% Markers_T]

rankings <- AUCell_buildRankings(Cryo_merge@assays$RNA@counts, plotStats = T)
cellsAUC <- AUCell_calcAUC(interferome_list, rankings)
threshold <- AUCell_exploreThresholds(cellsAUC, plotHist = T, thrP = 0.05)

cellsAUC2 <- AUCell_calcAUC(interferome_list2, rankings)
threshold2 <- AUCell_exploreThresholds(cellsAUC, plotHist = F)

if (all(
    rownames(Cryo_merge@meta.data) ==
        rownames(cellsAUC@colData@rownames)
)) {
    value <- cellsAUC@assays@data@listData$AUC[1, ]
    # thre <- threshold$geneSet$aucThr$selected
    thre <- threshold$geneSet$aucThr$thresholds["Global_k1", "threshold"]
    Cryo_merge$AUC <- value
    Cryo_merge$AUC_thre <- value * (value >= thre)

    FeaturePlot(Cryo_merge, features = "AUC_thre", cols = color_bar, pt.size = 1)
}

rankings2 <- AUCell_buildRankings(T_sub@assays$RNA@counts, plotStats = F)
cellsAUC22 <- AUCell_calcAUC(interferome_list, rankings2)
threshold22 <- AUCell_exploreThresholds(cellsAUC22, plotHist = T, thrP = 0.9)
if (all(
    rownames(T_sub@meta.data) ==
        rownames(cellsAUC22@colData@rownames)
)) {
    value <- cellsAUC22@assays@data@listData$AUC[1, ]
    # thre <- threshold22$geneSet$aucThr$selected
    thre <- threshold22$geneSet$aucThr$thresholds["Global_k1", "threshold"]
    T_sub$AUC <- value
    T_sub$AUC_thre <- value * (value >= thre)
    FeaturePlot(T_sub, features = "AUC", cols = color_bar, pt.size = 1)
}
