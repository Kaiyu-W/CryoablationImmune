library(Seurat)
library(SingleR)
library(scater)
library(SummarizedExperiment)
library(celldex) # reference RNA-seq data

library(harmony) # batch
library(SeuratWrappers) # batch, pseudotime and more

library(ggplot2) # plot
library(RColorBrewer)
library(pheatmap) # heatmap
library(ComplexHeatmap) # heatmap
library(clustree) # cluster resolution

library(clusterProfiler) # annotation analysis
library(enrichplot)
options(connectionObserver = NULL)
library("org.Mm.eg.db")
OrgDb = get("org.Mm.eg.db")
library("goProfiles")

library(GSVA) # annotation analysis
library(GSEABase)

# cluster by differential genes:
# library(clustermole)
# my_overlaps <- clustermole_overlaps(genes = my_genes_vec, species = "mm")

# sample.txt example:
# A data.frame: 16 × 7
# sample_ID   case_control    time    seq_type    sample_resource Prefix  batch
# <chr>   <chr>   <chr>   <chr>   <chr>   <chr>   <chr>
# Cryo1wk_RNA5end_2111    CA  1wk scRNA_5end  CD45CD3_pos CA_5r1w2b   batch_2111
# Cryo1wk_TCR_2111    CA  1wk scTCR   CD45CD3_pos CA_t1w2b    batch_2111
# Cryo2wk_RNA5end_2111    CA  2wk scRNA_5end  CD45CD3_pos CA_5r2w2b   batch_2111
# Cryo2wk_TCR_2111    CA  2wk scTCR   CD45CD3_pos CA_t2w2b    batch_2111
# Cryo_TCR_2106   CA  1wk scTCR   CD45_pos    CA_t1w1b    batch_2106
# Cryo_RNA3end_2103   CA  1wk scRNA_3end  CD45_pos    CA_3r1w1b   batch_2103
# Cryo_RNA3end_2109   CA  1wk scRNA_3end  CD45_pos    CA_3r1w2b   batch_2109
# Cryo_RNA5end_2106   CA  1wk scRNA_5end  CD45_pos    CA_5r1w1b   batch_2106
# NonCryo1wk_RNA5end_2111 Non_CA  1wk scRNA_5end  CD45CD3_pos NonCA_5r1w2b    batch_2111
# NonCryo1wk_TCR_2111 Non_CA  1wk scTCR   CD45CD3_pos NonCA_t1w2b batch_2111
# NonCryo2wk_RNA5end_2111 Non_CA  2wk scRNA_5end  CD45CD3_pos NonCA_5r2w2b    batch_2111
# NonCryo2wk_TCR_2111 Non_CA  2wk scTCR   CD45CD3_pos NonCA_t2w2b batch_2111
# NonCryo_RNA3end_2103    Non_CA  1wk scRNA_3end  CD45_pos    NonCA_3r1w1b    batch_2103
# NonCryo_RNA3end_2109    Non_CA  1wk scRNA_3end  CD45_pos    NonCA_3r1w2b    batch_2109
# NonCryo_RNA5end_2106    Non_CA  1wk scRNA_5end  CD45_pos    NonCA_5r1w1b    batch_2106
# NonCryo_TCR_2106    Non_CA  1wk scTCR   CD45_pos    NonCA_t1w1b batch_2106

check_datafile_exist <- function(
    IDs, 
    data_dir = "./data/output/", 
    suffix = "_filtered_feature_bc_matrix.h5"
) {
    for (ID in IDs) {
        file <- paste0(data_dir, ID, suffix)
        state <- file.exists(file)
        printout <- if (state) "checks OK!" else "did not exist!"
        print(paste(ID, printout))
    }
    if (length(IDs) == 1)
        return(state)
}

my_CreateSeuratObject_fromh5 <- function(
    ID, 
    prefix = ID, 
    case_control = NULL, 
    time = NULL, 
    seq_type = NULL, 
    sample_resource = NULL, 
    batch = NULL,
    data_dir = "./data/output/", 
    suffix = "_filtered_feature_bc_matrix.h5"
) {
    if (!check_datafile_exist(IDs = ID, data_dir = data_dir, suffix = suffix))
        stop("File did not exist! Check carefully please!")
    data_path <- paste0(data_dir, ID, suffix)
    
    data_mtx <- Read10X_h5(data_path)
    colnames(data_mtx) <- paste(prefix, colnames(data_mtx), sep = "_")

    Seurat_Obj <- CreateSeuratObject(
        counts = data_mtx,
        project = ID, 
        min.cells = 5, 
        min.features = 300
    )

    Seurat_Obj@meta.data$orig.ident <- ID
    if (length(case_control) != 0) 
        Seurat_Obj@meta.data$orig.group <- case_control
    if (length(time) != 0) 
        Seurat_Obj@meta.data$orig.time <- time
    if (length(seq_type) != 0) 
        Seurat_Obj@meta.data$orig.method <- seq_type
    if (length(sample_resource) != 0) 
        Seurat_Obj@meta.data$orig.resource <- sample_resource
    if (length(prefix) != 0) 
        Seurat_Obj@meta.data$orig.prefix <- prefix
    if (length(batch) != 0) 
        Seurat_Obj@meta.data$orig.batch <- batch

    return(Seurat_Obj)
}

my_process_seurat <- function(
    seurat_obj, 
    nVariableFeatures = 2500, 
    normalize = TRUE, 
    norm.method = c("LogNormalize", "SCT"), 
    default.assay = c("active.assay", names(seurat_obj@assays)),
    mt.to.regress = TRUE,
    tsne = TRUE # tsne run very slow when large data
) {
    match.arg(norm.method)
    norm.method <- norm.method[1] # default LogNorm
    match.arg(default.assay)
    default.assay <- default.assay[1]
    if (default.assay == "active.assay")
        default.assay <- DefaultAssay(seurat_obj)
    DefaultAssay(seurat_obj) <- default.assay

    obj_tmp <- seurat_obj

    if (normalize) {
        if (norm.method == "SCT") {
            if (mt.to.regress) {
                if (! 'percent.mt' %in% colnames(obj_tmp@meta.data))
                    stop("'percent.mt' slot did not exist in seurat_obj! Please assign that name if exists indeed!")
                obj_tmp <- SCTransform(obj_tmp, vars.to.regress = 'percent.mt')
            } else {
                obj_tmp <- SCTransform(obj_tmp)
            }    
            DefaultAssay(obj_tmp) <- 'SCT'
        } else {
            obj_tmp <- NormalizeData(obj_tmp, normalization.method = "LogNormalize")
            obj_tmp <- ScaleData(obj_tmp, features = rownames(obj_tmp))
        } 
    }
    
    if (nVariableFeatures != 0)
        obj_tmp <- FindVariableFeatures(obj_tmp, selection.method = "vst", nfeatures = nVariableFeatures)
    
    obj_tmp <- RunPCA(obj_tmp)
    obj_tmp <- FindNeighbors(obj_tmp, dims = 1:12)
    obj_tmp <- RunUMAP(obj_tmp, dims = 1:10)
    if (tsne) obj_tmp <- RunTSNE(obj_tmp, dims = 1:10)

    seurat_obj <- obj_tmp
    
    return(seurat_obj)
}

my_Integrate_Seurat <- function(
    Merge_Obj, 
    split_ident, 
    norm_method = c("LogNormalize", "SCT"),
    redu_method = c("rpca", "cca"), 
    fvf.nfeatures = 2000, 
    sif.nfeature = 2000
) {
    match.arg(norm_method)
    match.arg(redu_method)
    norm_method <- norm_method[1] # default LogNormalize
    redu_method <- redu_method[1] # default rpca
    if (norm_method == 'SCT')
        message("Process with redu_method of 'SCT' would be very slow...")

    DefaultAssay(Merge_Obj) <- "RNA"
    obj.list <- SplitObject(Merge_Obj, split.by = split_ident)
    
    my_NormFUN <- if (norm_method == "SCT") SCTransform else NormalizeData
    new.assay.name <- paste(
        "integrated", 
        redu_method, 
        ifelse(
            norm_method == "SCT", 
            'sct', 
            'log'
        ), 
        sep = "_"
    )

    # normalization
    obj.list <- lapply(
        X = obj.list, 
        FUN = function(x)
            FindVariableFeatures(
                my_NormFUN(x), 
                selection.method = "vst", 
                nfeatures = fvf.nfeatures
            )
    )

    # select features
    features <- SelectIntegrationFeatures(
        object.list = obj.list, 
        nfeatures = sif.nfeature, 
        fvf.nfeatures = fvf.nfeatures
    )

    # preSCT
    if (norm_method == "SCT")
        obj.list <- PrepSCTIntegration(
            object.list = obj.list, 
            anchor.features = features
        )
    
    # PCA
    if (redu_method == 'rpca')
        obj.list <- lapply(
            X = obj.list, 
            FUN = function(x) 
                if (norm_method == "SCT")
                    RunPCA(x, features = features)
                else
                    RunPCA(ScaleData(x, features = features), features = features)
        )
    
    # Find Anchors
    obj.anchors <- FindIntegrationAnchors(
        object.list = obj.list, 
        normalization.method = norm_method, 
        reduction = redu_method,
        anchor.features = features
    ) 

    # Integrate
    obj.integrated <- IntegrateData(
        anchorset = obj.anchors, 
        normalization.method = norm_method,
        new.assay.name = new.assay.name
    ) 

    # change assay
    DefaultAssay(obj.integrated) <- new.assay.name
    
    scale_if <- length(obj.integrated@assays[[new.assay.name]]@scale.data) != 1
    if (!scale_if)
        obj.integrated <- ScaleData(obj.integrated)
    
    return(obj.integrated)
}

my_CountCluster <- function(
    cluster_object, 
    group1 = 'active.ident', 
    group2 = 'orig.ident', 
    Heatmap = FALSE,
    trend_factor = 0,
    Percentage = TRUE, 
    Round = TRUE
) {
    # if trend_factor be not used, do not assign this variable
    # trend_factor means the other sum value and the present sum value 
    # should be considered together.
    meta <- cluster_object@meta.data
    meta$active.ident <- cluster_object@active.ident
    df <- meta[, c(group1, group2)]
    fren <- table(df)
    fren_perc <- as.matrix(fren)
    fren_perc <- apply(
        X = fren_perc, 
        MARGIN = 2,
        FUN = function(x) x / sum(x) * (100 ** Percentage)
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
        fren_prop <- round(fren_perc[,1] / fren_perc[,2] * trend_factor_all, 3)
        res <- cbind(res, fren_prop)
    }
    res <- rbind(res, apply(res, 2, sum))
    if (ncol(res) %% 2 == 1)
        res[nrow(res), ncol(res)] <- NA
    names(dimnames(res)) <- c(group1, group2)
    dimnames(res)[[1]][nrow(res)] <- "sum" 
    dimnames(res)[[2]] <- c(
        paste0(colnames(fren), "_count"), 
        paste0(colnames(fren), "_percentage"),
        "trend"
        )[1:ncol(res)]

    if (Round)
        res[, 3:4] <- round(res[, 3:4], ifelse(Percentage, 1, 3))
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

CountCluster_heatmap <- function(object, group1, group2 = 'orig.ident2', trend_factor = 1) {
    all_count_0 <- my_CountCluster(object, group1 = group1, group2 = group2, trend_factor = trend_factor)

    trend_factor_all <- ifelse(
        trend_factor == 1, 
        trend_factor,
        trend_factor * all_count_0[nrow(all_count_0), 1] / all_count_0[nrow(all_count_0), 2]
    )

    all_count <- apply(all_count_0[, 1:2], 2, function(x) { y <- x[-length(x)]; y / sum(y) })
    all_count[, 1] <- all_count[, 1] * trend_factor_all
    all_count <- apply(all_count, 1, function(x) if(max(x) == 0) x*0 else if(min(x) == 0) x / max(x) else log(x / min(x)))
    all_count <- t(all_count[, match(sort(levels(object@meta.data[, group1, drop = T])), colnames(all_count))])

    pheatmap::pheatmap(all_count, cluster_rows = F, cluster_cols = F, angle_col = 0, fontsize_row = 15, fontsize_col = 15)
    return(all_count_0)
}

CountCluster_plot <- function(object, group.by, xlable_y = -150, srt = 45) {
    count <- table(object[[group.by]])
    xlabels <- names(count)

    plot(count, xaxt = 'n')
    axis(1, at = seq(xlabels), labels = rep("", length(xlabels)), pos = 0)
    text(x = seq(xlabels), y = xlable_y, srt = srt, adj = 1, labels = xlabels, xpd = TRUE)
}

my_Markers2df_1Cluster <- function(
    df_1, 
    ntop, 
    logFC_threshold = 0.25, 
    positive = TRUE
) {
    index <- if (positive) 
                df_1$avg_log2FC >= logFC_threshold 
            else 
                df_1$avg_log2FC <= -logFC_threshold
    df_1 <- df_1[index, ]
    ntop <- min(ntop, nrow(df_1))
    if ('gene' %in% colnames(df_1))
        df_1$gene[1:ntop]
    else
        rownames(df_1)[1:ntop]
}

my_Markers2df_multiple <- function(
    df_input, 
    n_top,
    logFC_threshold = 0.25, 
    positive = TRUE
) {
    index <- if (positive) 
                df_input$avg_log2FC >= logFC_threshold 
            else 
                df_input$avg_log2FC <= -logFC_threshold
    df_input <- df_input[index, ]
    cluster <- levels(df_input$cluster)
    res <- rep(NA, n_top*length(cluster))
    dim(res) <- c(n_top, length(cluster))
    dimnames(res) <- list(1:n_top, paste("Cluster", cluster, sep = "_"))
    for (i in seq(cluster)) {
        cluster_name <- cluster[i]
        cluster_index <- which(df_input$cluster == cluster_name)[1:n_top]
        cluster_gene <- df_input$gene[cluster_index] 
        res[, i] <- cluster_gene
    }
    res <- as.data.frame(res)
    return(res)
}

my_MarkersList <- list(
    universal_markers = c('Ptprc','Epcam', 'Pecam1'),
    B = c("Cd19", "Cd79a"),
    T = c("Cd3e", "Cd247", "Cd4", "Cd8a", "Cd8b1"),
    NK = c("Klrb1c", "Ncr1", "Klrk1"),
    Myeloid = "Itgam",
    Macrophages = "Adgre1",
    Monocytes = "Cd14",
    marker_gene_cryo1 = list(
      Proliferating = c("Stmn1","Tuba1b","Mki67","Hmgn2","Birc5","Pclaf"),
      B = c("Cd79a","Ebf1"),
      NK = c("Irf8","Nrarp","Klre1","Txk","Car2"),
      Myeloid = c("Fcer1a","Tpsb2","Cpa3",'Cd63','Alox5ap',"Ctsh","Ly86","Lyz2","Cebpb","Cxcl2","Il1rn","Il1b","Il1r2","Ets2","Ccr1","Cd14"),
      T = c("Cd3g","Trbc2","Cd3d","Trac","Cd3e")
      ),
    T_exhausted = c("Pdcd1", "Ctla4", "Havcr2"),
    T_reg = "Foxp3",
    Plasma = c("Sdc1", "Jchain"), #CD138
    DCs_activated = "Cd83",
    # pro_metastatic_factors = c("S100a8", "S100a9"),
    Neutrophils = c("S100a8", "S100a9"),
    T_NK_activate = 'Cd69',
    T_NK_cytotoxic = c('Prf1', 'Gzma', 'Gzmb'),
    Mast = c('Mcpt1', 'Mcpt2', 'Mcpt4', 'Tpsb2', 'Cma1'),
    pDC = c('Cox6a2', 'Siglech', 'Klk1', 'Smim5')
    )

celltype_marker <- list(
    epithelial = "Epcam",
    endothelial= "Pecam1",
    fibroblasts= "Col3a1", 
    Myeloid = c("Cd163", "Aif1"),
    B = "Cd79a",
    Plasma = "Jchain",
    T = c("Cd3e","Cd4","Cd8a","Cd8b1"),
    Nk = "Nkg7",
    immune = "Ptprc"
    )

CD45_Markers <- list(
    immune = "Ptprc",
    epithelial = "Epcam",
    endothelial= "Pecam1",
    fibroblasts= "Col3a1", 
    
    B = c("Cd19", "Cd79a"),
    Plasma = c("Sdc1", "Jchain"),
    
    T = c("Cd3e","Cd4", "Cd8a", "Cd8b1"),
    Nk = c("Klrb1c", "Ncr1", "Klrk1", "Nkg7"),
    
    Myeloid = c("Itgam", "Cd14", "Cd163", "Aif1"),
    Pan_Granulocyte = c("Fcgr3", "Fcgr2b"),
    Neutrophils = "Cdk5",
    Eosinophils = c("Il5ra", "Ccr3"),
    Basophil = "Il3ra",
    Mast = c("Kit", "Mcpt1", "Mcpt2", "Mcpt4"),
    'Basophils/Mast' = "Fcer1a"
)

T_marker <- list(
    T = c("Cd3e", "Cd3d", "Cd3g"),
    Th = "Cd4",
    CD4 = list(
        Naive = 'Sell',
        Memory = c('Cd44', 'Cd40lg'),
        Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
        Th1 = c('Tbx21', 'Cxcr3'), # 'Il12rb1', 'Ifngr1'
        Th2 = c('Gata3', 'Il1rl1'), # 'Il4ra', 'Ccr4', 'Il17rb', 'Ptgdr2'
        # Th9 = c('Spi1', 'Il9'), 
        Th17 = c('Rorc', 'Il17a', 'Il22', 'Ccr6', 'Il23r'), # the first 2 are main
        T_FH = c('Il21r', 'Cxcr5', 'Pdcd1', 'Sh2d1a'), # 'Icos', 'Il21'
        Treg = c('Foxp3', 'Il2ra', 'Cd69') # 'Ctla4', 'Tnfrsf18'
        ),
    Cyto_T = c('Cd8a', 'Cd8b1'),
    CD8 = list(
        Activated = c("Cd69", "Il2ra"),
        Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
        Naive = c("Sell"),
        Effector_Memory = c("Il7r"),
        Effector = c("Cd44"),
        Exhausted = c("Pdcd1", "Ctla4", "Lag3", "Havcr2")
        )
    )

T_markers_tmp <- list(
    CD4 = "Cd4",
    CD8 = c('Cd8a', 'Cd8b1'),
    Naive = 'Sell',
    # CD8
    Activated = c("Cd69", "Il2ra"),
    Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
    Effector = "Cd44",
    Memory = c('Il7r', 'Cd40lg'),
    Exhausted = c("Pdcd1", "Ctla4", "Lag3", "Havcr2"),
    # CD4
    Th1 = c('Tbx21', 'Cxcr3'),
    Th2 = c('Gata3', 'Il1rl1'),
    Th17 = c('Rorc', 'Il17a'),
    T_FH = c('Il21r', 'Cxcr5', 'Sh2d1a'),
    Treg = 'Foxp3'
    )

T_markers_tmp2 <- list(
    CD4 = "Cd4",
    CD8 = c('Cd8a', 'Cd8b1'),
    Naive = 'Sell',
    # CD8
    Activated = c("Cd69", "Il2ra"),
    Cytotoxic = c('Gzma', 'Gzmb'),
    Effector = "Cd44",
    Memory = c('Il7r'),
    Exhausted = c("Pdcd1", "Ctla4", "Lag3"),
    # CD4
    Th1 = c('Cxcr3'),
    Th2 = c('Gata3'),
    Th17 = c('Il17a'),
    T_FH = c('Il21r', 'Sh2d1a'),
    Treg = 'Foxp3'
    )

T_markers_tmp3 <- list(
    CD4 = "Cd4",
    CD8 = c('Cd8a', 'Cd8b1'),
    NKT = c("Cd3e", "Klrb1c", "Klrk1"), 
    Naive = 'Sell',
    # CD8
    Activated = c("Cd69", "Il2ra"),
    Cytotoxic = c('Gzma', 'Gzmb'),
    Effector = "Cd44",
    Memory = c('Il7r'),
    Exhausted = c("Pdcd1", "Ctla4", "Lag3"),
    # CD4
    Th1 = c('Cxcr3'),
    Th2 = c('Gata3'),
    Th17 = c('Il17a'),
    T_FH = c('Il21r', 'Sh2d1a'),
    Treg = 'Foxp3'
    )

TNK_marker <- list(
    T = c("Cd3e", "Cd4", 'Cd8a'),
    NK = c("Klrb1c", "Klrk1"),
    Naive = 'Sell',
    Memory = c('Cd44', "Il7r"),
    Cytotoxic = c('Gzma', 'Gzmb', 'Prf1'),
    'Activated/Treg' = c("Cd69", "Il2ra"),
    Proliferating = 'Mki67',
    Exhausted = c("Lag3", "Havcr2"),
    'Exhausted/Tfh' = "Pdcd1", 
    'Exhausted/Treg' = "Ctla4", 
    Th1 = c('Tbx21', 'Cxcr3'), # 'Il12rb1', 'Ifngr1'
    Th2 = c('Gata3', 'Il1rl1'), # 'Il4ra', 'Ccr4', 'Il17rb', 'Ptgdr2'
    # Th9 = c('Spi1', 'Il9'), 
    Th17 = c('Rorc', 'Il17a', 'Ccr6', 'Il23r'), # the first 2 are main
    Tfh = c('Il21r', 'Cxcr5', 'Sh2d1a'), # 'Icos', 'Il21'
    Treg = 'Foxp3'
    )

NKT_marker <- c('Klrb1c', 'Slamf1', 'Slamf6', 'Tgfbr1', 'Tcra-V14', 'Tcra-J')

NK_marker <- list(
    NK = c("Klrb1c", "Ncr1", "Klrk1", "Nkg7"),
    NK_Activated = "Cd69",
    NK_Cytotoxic = c('Gzma', 'Gzmb', 'Prf1')
    )

Myeloid_marker <- list(
    Myeloid = "Itgam",
    'Macrophages/Neutrophils_med/Eosinophils_med' = "Adgre1",
    Macrophages_M1 = c('Cd86', 'Cd80', 'Nos2'),
    Macrophages_M2 = c('Cd163', 'Mrc1'),
    Monocytes = c("Cd14", "Fcgr3"),
    cDC = c('Itgax', 'H2'), 
    pDC = c('Bst2', 'Siglech'),
    # pDC = c('Bst2', 'Siglech', 'Cox6a2', 'Siglech', 'Klk1', 'Smim5'),
    DCs_activated = "Cd83",
    'M2/MDSCs' = c('Arg1', 'Arg2'),
    M_MDSCs = 'Ly6c1',
    'PMN_MDSCs/Neutrophils' = 'Ly6g',
    Neutrophils = c("S100a8", "S100a9"), # Cdk5
    Eosinophils = c('Ccr3', 'Siglecf'), # Il5ra
    'Basophils/Mast' = 'Fcer1a', # Basophils -> Il3ra
    # Mast = c('Mcpt1', 'Mcpt2', 'Mcpt4', 'Tpsb2', 'Cma1'),
    Mast = c('Fcer1g', 'Kit')
)

invitrogen_Granulocyte <- list(
    Pan = c("Fcgr3", "Fcgr2b"),
    Neutrophils = "Cdk5",
    Eosinophils = c("Il5ra", "Ccr3"),
    Basophil = "Il3ra",
    Mast = "Kit",
    'Basophils/Mast' = "Fcer1a"
)

cDC_marker <- list(
    cDC = c('Itgax', 'H2', 'Cd86', 'Cd80'),
    XCR1_CLEC9A_plus = c("Xcr1", "Clec9a"),
    Resident = 'Cd8a',
    Migratory = 'Itgae',
    Activated = 'Cd83',
    CD11b_plus = 'Itgam'
    # cDC1 = 'Thbd',
    # cDC2 = 'Cd1'
)

# cell type identify
my_ClusterAnnote <- function(
    Object, 
    anno_list, 
    meta_slot, 
    anno_slot, 
    DefaultIdent_changeIf = TRUE,
    Duplicated_removeIf = TRUE
) {
    Object$ClusterAnnote_tmp <- as.factor(Object@meta.data[[meta_slot]])
    
    ClusterAnnote_df <- NULL
    for (anno_name in names(anno_list)) {
        df_tmp <- data.frame(
            x = anno_name,
            y = anno_list[[anno_name]]
        )
        if (is.null(ClusterAnnote_df))
            ClusterAnnote_df <- df_tmp
        else
            ClusterAnnote_df <- rbind(
                ClusterAnnote_df, 
                df_tmp
            )
    }
    ClusterAnnote_df <- ClusterAnnote_df[order(ClusterAnnote_df$y, decreasing = F), ]
    
    if (all(ClusterAnnote_df$y == levels(Object$ClusterAnnote_tmp))) {
        levels(Object$ClusterAnnote_tmp) <- ClusterAnnote_df$x
    } else {
        print(levels(Object$ClusterAnnote_tmp))
        print(ClusterAnnote_df)
        stop("ERROR! Different levels!")
    } 
    
    if (anno_slot %in% colnames(Object@meta.data)) {
        colname_index1 <- colnames(Object@meta.data) == anno_slot
        if (Duplicated_removeIf)
            Object@meta.data[, colname_index1] <- NULL
        else
            stop("Duplicated anno_slot! Please set Duplicated_removeIf as TRUE to avoid unknown error.")
    }

    colname_index <- colnames(Object@meta.data) == 'ClusterAnnote_tmp'
    colnames(Object@meta.data)[colname_index] <- anno_slot
    
    if (DefaultIdent_changeIf)
        Idents(Object) <- anno_slot
    
    return(Object)
}

rename_colname <- function(input_mtx, name_pattern) {
    # name_pattern should be a vector of unique name, with regexp as well.

    colnames_input_mtx <- colnames(input_mtx)
    res_temp <- lapply(
        name_pattern, 
        function(x) {
            index <- grepl(pattern = x, colnames_input_mtx)
            name_temp <- colnames_input_mtx[index]
            name_temp <- sort(name_temp)
            name_temp
        }
    )
    names(res_temp) <- name_pattern

    res <- unname(unlist(res_temp))

    if (length(res) > length(colnames_input_mtx)) {
        warning("Colnames of input mtx/df are not involved in name_pattern !")
    } else {
        # index <- sapply(res, function(x) which(x == colnames(input_mtx)))
        # input_mtx <- input_mtx[, index, drop = F]
        input_mtx <- input_mtx[, res, drop = F]
    }

    input_mtx
}

reorder_df <- function(df, colname, name_pattern) {
    # name_pattern should be a vector of unique name, with regexp as well.

    vv <- df[, colname]
    res_temp <- lapply(
        name_pattern, 
        function(x) {
            which(grepl(pattern = x, vv))
        }
    )
    names(res_temp) <- name_pattern

    res <- unname(unlist(res_temp))

    if (length(res) > length(vv)) {
        warning("Values of the col in input mtx/df are not involved in name_pattern !")
    } else {
        index <- as.numeric(res)
        df <- df[index, , drop = F]
    }

    df
}

my_plotQC <- function(object, ...) {
    Seurat::VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, ...)
}
my_plotFC <- function(object, combinedIf = TRUE, ...) {
    plot1_m <- Seurat::FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt", ...)
    plot2_m <- Seurat::FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.mt", ...)
    plot3_m <- Seurat::FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", ...)
    if (combinedIf) {
        plot1_m + plot2_m + plot3_m
    } else {
        print(plot1_m)
        print(plot2_m)
        print(plot3_m)
    }
}
my_QC_vis <- function(object, FC_combinedIf = FALSE) {
    # visualization
    print(my_plotQC(object, pt.size = 0))
    print(my_plotFC(object, combinedIf = FC_combinedIf, pt.size = 0.5))
    hist(object$nFeature_RNA, breaks = 500, xaxt = 'n')
    axis(1, seq(0, min(10000, max(object$nFeature_RNA)), 500))
    hist(object$nCount_RNA, breaks = 500, xaxt = 'n')
    axis(1, seq(0, min(200000, max(object$nCount_RNA)), 2000))
    hist(object$percent.mt, breaks = 500, xaxt = 'n')
    axis(1, seq(0, 100, 1))
}

my_plotDim <- function(..., title = NULL) {
    # if (is.null(title) & exists('group.by'))
    #     title <- group.by
    #   if group.by exists, then seurat plot function 
    #   will add title as group.by automatically
    if (is.null(title))
        DimPlot(...)
    else
        DimPlot(...) + labs(title = title) + theme(plot.title = element_text(hjust = 0.5))
}
my_violin <- function(..., mode = 'raw') {
    mode_mode <- c('mtx', 'line', 'raw')
    index <- grep(mode, mode_mode)
    switch(
        index,
        # VlnPlot(..., stack = TRUE, flip = TRUE) + ggsci::scale_fill_jco(),
        VlnPlot(..., stack = TRUE, flip = TRUE),
        patchwork::wrap_plots(plots = VlnPlot(..., combine = FALSE), ncol = 1),
        VlnPlot(...)
        )
}
my_DotPlot_split <- function(object, split.by = NULL, cols = c('red', 'blue'), ...) {
    if (is.null(split.by)) {
        Seurat::DotPlot(object, ...)
    } else {
        p <- Seurat::DotPlot(object, split.by = split.by, cols = cols, ...)

        groupnames <- unique(object[[split.by, drop = TRUE]])
        cols <- cols
        names(cols) <- groupnames

        data.plot <- p$data
        data.plot$split.by <- unlist(lapply(data.plot$id, function(x) {
            res <- x
            while (! res %in% groupnames) {
                res <- gsub(".+?_", "", x)
            }
            return(res)
        }))

        data.plot$split.color <- cols[data.plot$split.by]

        data.plot$colors <- unlist(apply(data.plot, 1, function(x) {
            color <- x[["split.color"]]
            value <- as.numeric(x[["avg.exp.scaled"]])
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }))

        p$data$colors <- data.plot$colors

        # print(p)
        return(p)
    }
}

my_Heatmap <- function(
    seurat_obj, 
    group.by, 
    genes,
    slot = c('data', 'scale.data', 'counts'), 
    default.assay = c("active.assay", names(seurat_obj@assays)), 
    cells = NULL, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE,
    show_rownames = FALSE,
    color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    use_pheatmap = FALSE,
    border_color = NA, # NA if no border to draw
    ...
) {
    match.arg(default.assay)
    match.arg(slot)
    default.assay <- default.assay[1]
    slot <- slot[1]
    if (default.assay == "active.assay")
        default.assay <- DefaultAssay(seurat_obj)
    
    data <- methods::slot(seurat_obj@assays[[default.assay]], slot)
    anno_df <- seurat_obj@meta.data[, group.by, drop = F]

    if (!is.null(cells)) {
        data <- data[, cells, drop = F]
        anno_df <- anno_df[cells, , drop = F]
    }

    anno_df <- anno_df[order(anno_df[, 1, drop = T]), , drop = F]
    genes_in <- genes[genes %in% rownames(data)]
    mtx <- data[genes_in, rownames(anno_df)]
    
    fun <- if (use_pheatmap)
        pheatmap::pheatmap
    else
        ComplexHeatmap::pheatmap
    fun(
        mat = as.matrix(mtx), 
        color = color,
        cluster_cols = cluster_cols, 
        cluster_rows = cluster_rows, 
        annotation_col = anno_df, 
        show_colnames = show_colnames,
        show_rownames = show_rownames,
        border_color = border_color,
        ...
    )
}

my_GO <- function(
    object, 
    ont = "ALL", 
    type = 'dot', 
    pvalue = 0.05,
    qvalue = 0.2, 
    show = 20, 
    font.size = 10, 
    title = NULL, 
    Simplify = FALSE, 
    return_plot = TRUE, 
    return_res = FALSE,
    pAdjustMethod = 'BH',
    nChar = 55
) {
    message("object can be either Geneset or GO_output")
    message('ont can be one of "BP", "MF", "CC" and "ALL"(default)')
    message('type can be some of "bar", "dot", "cnet" and "emap"')
    if (!is.numeric(nChar))
        stop('nChar should be one number!')
    if (Simplify & ont == "ALL") 
        stop('simplify function only applies to a single ontology!')
    if (is.null(title))
        title <- ont

    if (class(object) == "enrichResult") {
        ego <- object
    } else {
        Geneset <- object
        Geneset <- Geneset[!is.na(Geneset)]
        genes <- bitr(Geneset, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
        ego <- enrichGO(gene = genes$ENTREZID, ont = ont, OrgDb = OrgDb, pvalueCutoff = pvalue, qvalueCutoff = qvalue, pAdjustMethod = pAdjustMethod)
    }
    
    if (nrow(ego) != 0) {
        if (Simplify) 
            ego <- clusterProfiler::simplify(ego)
        ego1 <- ego
        ego@result$Description <- sapply(
            ego@result$Description,
            function(x) {
                ifelse(
                    nchar(x) < nChar, 
                    x, 
                    paste0(paste(strsplit(x, split = "")[[1]][1:nChar], collapse = ""), "...")
                    )
            }
        )
        if ('bar' %in% type) 
            res <- barplot(ego, showCategory = show, font.size = font.size, title = paste("GO", title))
        if ('dot' %in% type) 
            res <- dotplot(ego, showCategory = show, font.size = font.size, title = paste("GO", title))
        if ('cnet' %in% type) 
            res <- cnetplot(ego, showCategory = show, circular = T, colorEdge = T, categorySize = 'p.adjust', title = paste("GO", title))
        if ('emap' %in% type) {
            ego1 <- pairwise_termsim(ego, method = "JC", semData = NULL, showCategory = show)
            res <- emapplot(ego1, showCategory = show)
        }
        if ('go' %in% type) {
            if (ont == "ALL") {
                message("Ontology being ALL cannot goplot! Transfer into one of the three!")
            } else {
                res <- goplot(ego, showCategory = show, color = "p.adjust", layout = "sugiyama", geom = "text")
            }
        }
        if ('top' %in% type) {
            if (ont == "ALL") {
                message("Ontology being ALL cannot plotGOgraph! Transfer into one of the three!")
            } else {
                res <- plotGOgraph(ego, firstSigNodes = show, useInfo = "all", sigForAll = TRUE, useFullNames = TRUE)
            }
        }
        if (return_plot & return_res) print(res)
        if (return_plot & !return_res) return(res)
        if (return_res) return(ego1)
    } else {
        message("No output information")
        return(ego)
    }
}

my_GSEA <- function(
    FC, 
    symbols, 
    anno = "GO", 
    ont = "All", 
    pAdjustMethod = "BH", 
    pvalueCutoff = 0.05,
    show = 20, 
    title = NULL, 
    font.size = 12,
    Simplify = FALSE,
    return_plot = TRUE, 
    return_res = FALSE,
    ...
) {
    message('The parameters first and second "FC" and "symbols" are respectively the input FoldChange vector and corresponding gene symbols')
    message('pAdjustMethod can be one of "holm", "hochberg", "hommel", "bonferroni", "BH"(default), "BY", "fdr" and "none"')
    message('Annotation method can be some of "GO"(default), "KEGG" and "WP"')
    message('If GO_annotation, ont can be one of "BP", "MF", "CC" and "ALL"(default)')
    if (Simplify & ont == "ALL") 
        stop('simplify function only applies to a single ontology!')
    names(FC) <- symbols
    genes <- symbols
    genes <- bitr(genes, fromType = "SYMBOL", toType = 'ENTREZID', OrgDb = OrgDb)
    FC <- FC[names(FC) %in% genes$SYMBOL]
    genes <- genes[which(genes$SYMBOL %in% names(FC)), , drop = FALSE]
    names(FC) <- as.character(genes$ENTREZID)
    FC <- FC[!is.infinite(FC)]
    FC <- FC[!is.nan(FC)]
    FC <- sort(FC, decreasing = T)
    res <- list()
    if ("GO" %in% anno) {
        Gse_GO <- gseGO(FC, ont = ont, OrgDb = OrgDb, pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, ...)
        if (nrow(Gse_GO@result) != 0) {
            if (Simplify)
                Gse_GO <- simplify(Gse_GO)
            Dot_Gse_GO <- dotplot(Gse_GO, split=".sign", showCategory = show, title = paste("GSE_GO", title), font.size = font.size) + facet_wrap(~.sign, scales = "free")
            if (return_plot)
                print(Dot_Gse_GO)
        }
        res <- append(res, list(Gse_GO = Gse_GO))
    }
    if ("KEGG" %in% anno) {
        Gse_KEGG <- gseKEGG(FC, organism = "mmu", pAdjustMethod = pAdjustMethod)
        if (nrow(Gse_KEGG@result) != 0) {
            Dot_Gse_KEGG <- dotplot(Gse_KEGG, split=".sign", showCategory = show, title = paste("GSE_KEGG", title), font.size = font.size) + facet_wrap(~.sign, scales = "free")
            if (return_plot)
                print(Dot_Gse_KEGG)
        }
        res <- append(res, list(Gse_KEGG = Gse_KEGG))
    }
    if ("WP" %in% anno) {
        Gse_WikiPathway <- gseWP(FC, organism = "Mus musculus")
        if (nrow(Gse_WikiPathway@result) != 0) {
            Dot_Gse_WikiPathway <- dotplot(Gse_WikiPathway, showCategory = show, title = paste("GSE_WikiPathway", title), font.size= font.size)
            if (return_plot)
                print(Dot_Gse_WikiPathway)
        }
        res <- append(res, list(Gse_WikiPathway = Gse_WikiPathway))
    }
    if (return_res)
        return(res)
}

my_fc <- function(
    object, 
    meta_slot = NULL, 
    Clusters = NULL,
    Feature = c("ALL", "Variable")[1],
    LOG = TRUE,
    log_func = log1p,
    average_method = c("arithmetic", "geometric")[1],
    showProcess = FALSE
) {
    # if (class(data)=="dgCMatrix")
    #     data <- as.matrix(data)
    # data <- as.matrix(object@assays[[object@active.assay]]@data)
    data <- object@assays[[object@active.assay]]@data
    if (Feature == 'Variable') {
        var_names <- object@assays[[object@active.assay]]@var.features
        data <- data[var_names, ]
    }
    if (is.null(meta_slot))
        meta <- object@active.ident
    else
        meta <- object[[meta_slot]][[1]]
    if (is.null(Clusters)) {
        meta_levels <- meta
    } else {
        meta_index <- which(meta %in% Clusters)
        meta_levels <- meta[meta_index]
        data <- data[, meta_index]
    }
    if (average_method == "arithmetic") {
        my_fc_fun <- function(In, All)
            mean(In) / (mean(All) - mean(In))
    } else if (average_method == "geometric") {
        my_mean <- function(...) expm1(mean(log1p(...)))
        my_fc_fun <- function(In, All)
            (my_mean(In) / my_mean(All)) ** (length(All) / (length(All) - length(In)))
    } else {
        stop("ERROR average method! Should be one of 'arithmetic' and 'geometric'!")
    }
            
    if (!showProcess) {
        FC <- apply(
            X = data, 
            MARGIN = 1,
            FUN = function(x) 
                tapply(x, meta_levels, function(y) my_fc_fun(In = y, All = x))
        )
    } else {
        pb <- txtProgressBar(style = 3)
        iii <- 0
        iii_sum <- nrow(data)
        FC <- apply(
            X = data, 
            MARGIN = 1,
            FUN = function(x) {
                fc <- tapply(x, meta_levels, function(y) my_fc_fun(y, x))
                setTxtProgressBar(pb, iii / iii_sum)
                iii <<- iii + 1
                fc
            }
        )
        close(pb)
    }
    FC <- if (LOG) log_func(FC) else FC
    res <- if (!is.null(dim(FC))) t(FC) else FC
    return(res)
}

my_avg <- function(...) {
    return(Seurat::AverageExpression(...))
}

my_rank_avgExp <- function(
    object, 
    meta_slot = 'orig.ident', 
    assay = 'RNA', 
    slot = 'data'
) {
    avgExp <- Seurat::AverageExpression(
        object, 
        group.by = meta_slot, 
        assays = assay, 
        slot = slot
        )[[1]]
    avgExp <- as.data.frame(avgExp)
    output <- apply(
        X = avgExp, 
        MARGIN = 2,
        FUN = function(x) {
            index <- order(x, decreasing = T)
            rownames(avgExp)[index]
            }
        )

    as.data.frame(output)
}

my_ReferenceCluster <- function(
    object, 
    reference, 
    celltype_list = NULL,
    label_fine = FALSE,
    forPlot_likelihood = FALSE
) {
    message("Input Seurat object and celldex reference object simultaneously as well as cell types list if required.\nThis process would cost plenty time.\nYou'll get the Seurat object with a new slot named by input reference.")
    message("DB_Ig_celltype_list = c('B cells', 'T cells', 'NK cells', 'DC', 'ILC', 'Macrophages', 'Mast cells', 'Monocytes', 'Neutrophils', 'Eosinophils', 'B cells, pro')")
    message("DB_Mrna_celltype_list = c('T cells', 'B cells', 'Monocytes', 'Granulocytes', 'Dendritic cells', 'NK cells', 'Macrophages')")
    message("DB_Mi_celltype_list = c('CD4+ T cells', 'B cells', 'Monocytes', 'T cells', 'CD8+ T cells', 'Neutrophils', 'NK cells', 'Dendritic cells', 'Basophils')")
    message("DB_Ice_celltype_list = c('T cells, CD4+', 'B cells', 'Monocytes', 'T cells, CD8+', 'NK cells')")
    message("DB_Be_celltype_list = c('CD4+ T-cells', 'B-cells', 'Monocytes', 'NK cells', 'CD8+ T-cells', 'Neutrophils', 'DC', 'Macrophages')")

    new_meta_name <- deparse(substitute(reference))

    rownames(reference) <- lapply(
        rownames(reference), 
        function(x) {
            y <- strsplit(x, split ="")[[1]]
            y[-1] <- tolower(y[-1])
            paste0(y, collapse = "")
        }
        )
    if (!is.null(celltype_list)) {
        if (!all(celltype_list %in% reference$label.main))
            stop("Cell-type list input didn't match the reference!")
        reference <- reference[, reference$label.main %in% celltype_list]
    }
    test_temp <- as.data.frame(object[["RNA"]]@counts)
    common_test <- intersect(rownames(test_temp), rownames(reference))
    reference.test <- reference[common_test, ]
    test_temp_forref <- test_temp[common_test, ]
    test_temp_forref.test <- SummarizedExperiment(
        assays = list(counts = test_temp_forref)
        )
    test_temp_forref.test <- logNormCounts(test_temp_forref.test)

    if (label_fine) {
        labels <- reference.test$label.fine
    } else {
        labels <- reference.test$label.main
    }

    test.main.ref <- SingleR(
        test = test_temp_forref.test, 
        ref = reference.test, 
        labels = labels
        )
    if (forPlot_likelihood) 
        return(test.main.ref)

    result_main_ref <- as.data.frame(test.main.ref$labels)
    result_main_ref$CB <- rownames(test.main.ref)
    colnames(result_main_ref) <- c(new_meta_name, 'CB')

    Barcodes <- rownames(object@meta.data)
    object@meta.data$CB <- Barcodes
    object@meta.data <- merge(
        object@meta.data, 
        result_main_ref, 
        by = "CB"
        )
    rownames(object@meta.data) <- object@meta.data$CB
    object$CB <- NULL
    object@meta.data <- object@meta.data[Barcodes, ]

    return(object)
}

my_ReferenceCelldexCreate <- function(download = FALSE) {
    message("
        For example, celldex::ImmGenData()
        DB_Ig <- celldex::ImmGenData()
        DB_Be <- celldex::BlueprintEncodeData()
        DB_Ice <- celldex::DatabaseImmuneCellExpressionData()
        DB_Hpca <- celldex::HumanPrimaryCellAtlasData()
        DB_Mi <- celldex::MonacoImmuneData()
        DB_Mrna <- celldex::MouseRNAseqData()
        DB_Nh <- celldex::NovershternHematopoieticData()
        ")
    if (!download) return(NULL)
    if (!file.exists("C://celldex_reference/Database_form_celldex.rda")) {
        DB_Ig <- celldex::ImmGenData()
        DB_Be <- celldex::BlueprintEncodeData()
        DB_Ice <- celldex::DatabaseImmuneCellExpressionData()
        DB_Hpca <- celldex::HumanPrimaryCellAtlasData()
        DB_Mi <- celldex::MonacoImmuneData()
        DB_Mrna <- celldex::MouseRNAseqData()
        DB_Nh <- celldex::NovershternHematopoieticData()
        save(DB_Ig, DB_Be, DB_Ice, DB_Hpca, DB_Mi, DB_Mrna, DB_Nh, file = "C://celldex_reference/Database_form_celldex.rda")
    }
    
}

my_ReferenceScore <- function(
    SingleR_object, 
    meta_vector_withName, 
    meta_order = NULL,
    ref_order = NULL,
    avg = FALSE, 
    avg_method = c('arithmetic', 'geometric', 'harmonic')[1], 
    cluster_cols = F, cluster_rows = F, ...
) {
    if (is.null(names(meta_vector_withName)))
        stop("Wrong meta input! Vector with no names(barcodes)!")
    meta_vector <- meta_vector_withName

    test <- SingleR_object[, 'scores', drop = F]
    test <- as.data.frame(test)
    colnames(test) <- sub("^scores.", "", colnames(test))
    meta_vector <- meta_vector[names(meta_vector) %in% rownames(test)]
    ref_order <- if (is.null(ref_order)) colnames(test) else ref_order
    ref_order <- if (all(ref_order %in% colnames(test))) ref_order else colnames(test)
    
    if (!avg) {
        annotation_col_x <- data.frame(celltype = unname(meta_vector),
                                       row.names = names(meta_vector))
        annotation_col_x <- annotation_col_x[order(annotation_col_x$celltype), , drop = F]
        
        if (!is.null(meta_order)) 
            annotation_col_x <- reorder_df(df = annotation_col_x, colname = 'celltype', name_pattern = meta_order)

        gap_x <- sapply(
            unique(annotation_col_x$celltype), 
            function(x) which(annotation_col_x$celltype == x)[1]
            )[-1] -1
         
        test.forPlot <- t(test[rownames(annotation_col_x), ref_order])
        
        Plot <- pheatmap::pheatmap(
            test.forPlot, 
            show_colnames = F, show_rownames = T, 
            annotation_col = annotation_col_x, 
            gaps_col = gap_x, 
            cluster_cols = cluster_cols, cluster_rows = cluster_rows, 
            ...
        )
    } else {
        if (!all(names(meta_vector) == rownames(test)))
            meta_vector <- meta_vector[rownames(test)]

        avg_warning <- 'Wrong input of avg_method! avg_method can be character/numeric: 1/arithmetic; 2/geometric; 3/harmonic. Default by mean(arithmetic)'
        if (is.character(avg_method)) {
            avg_method_rank <- match(avg_method, c('arithmetic', 'geometric', 'harmonic'))
            if (is.na(avg_method_rank)) {
                avg_method_rank <- 1
                warning(avg_warning)
            }
        } else if (is.numeric(avg_method) & avg_method > 0 & avg_method < 4) {
            avg_method_rank <- avg_method 
        } else { 
            avg_method_rank <- 1
            warning(avg_warning)
        }
        avg_fun <- switch(
                avg_method_rank,
                mean,
                function(x) exp(sum(log(x + 1e-5 - min(x))) / length(x)) + min(x),
                function(x) length(x) / sum(1 / x)
                )

        test.temp <- apply(
            test, 
            2, 
            function(x) {
                tapply(x, unname(meta_vector), avg_fun)
            }
        )

        # meta_order <- if (is.null(meta_order)) rownames(test.temp) else meta_order
        # meta_order <- if (all(meta_order %in% rownames(test.temp))) meta_order else rownames(test.temp)
        test.forPlot <- test.temp[, ref_order, drop = F]
        test.forPlot <- if (is.null(meta_order)) test.forPlot else t(rename_colname(t(test.forPlot), meta_order))

        Plot <- pheatmap::pheatmap(
            test.forPlot, angle_col = 45, 
            show_colnames = T, show_rownames = T, 
            cluster_cols = cluster_cols, cluster_rows = cluster_rows,
            ...
        )
    }

    return(Plot)
}

my_ReferenceLabelStat <- function(
    SingleR_object, 
    meta_vector_withName, 
    label_col = c('first.labels', 'labels', 'pruned.labels')[2],
    Heatmap = T, angle_col = 45, ...
) {
    if (is.null(names(meta_vector_withName)))
        stop("Wrong meta input! Vector with no names(barcodes)!")
    meta_vector <- meta_vector_withName

    test <- SingleR_object[, label_col, drop = F]
    test <- as.data.frame(test)
    colnames(test) <- 'label'
    test <- test[names(meta_vector), , drop = F]

    Stat_table <- tapply(test$label, meta_vector, table)
    
    row <- unique(test$label)
    col <- unique(unlist(meta_vector_withName))

    res_mtx <- matrix(0, nrow = length(row), ncol = length(col), dimnames = list(row, col))
    for (i in names(Stat_table)) {
        Stat_temp <- Stat_table[[i]]
        index <- names(Stat_temp)
        res_mtx[index, i] <- Stat_temp
    }

    if (Heatmap) {
        res_mtx_temp <- res_mtx
        res_mtx_temp <- apply(res_mtx_temp, 2, function(x) x / sum(x) * 100)
        pheatmap::pheatmap(res_mtx_temp, cutree_rows = nrow(res_mtx_temp), cutree_cols = ncol(res_mtx_temp), angle_col = angle_col, ...)
    }

    res_mtx
}

my_AddMeta <- function(
    object, 
    new_ident,
    Replace = FALSE,
    others_name = "Others",
    allow_mismatch = FALSE
) {
    message("Input new_ident should be a data.frame with one column as well as cell-barcode rownames, or a list(factor) with cell-barcode names")
    message("New meta name will be the colname of input data.frame , or object_name of inputting list(factor)")
    message("If Replace, cells with not-matched barcodes will add the mata label of object@active.ident")
    newSlotName <- deparse(substitute(new_ident))
    newSlotName <- gsub("\\$", "_", newSlotName)

    # test input barcodes
    Barcodes <- rownames(object@meta.data)
    if (class(new_ident) == 'factor') {
        newBC <- names(new_ident)
    } else if (class(new_ident) == 'data.frame') {
        newBC <- rownames(new_ident)
    } else if (class(new_ident) == 'character') {
        new_ident <- as.factor(new_ident)
        newBC <- names(new_ident)
    } else {
        stop("Not allowed format input!\nSuggest 'Factor' format!")
    }
    if (!all(newBC %in% Barcodes)) {
        if (allow_mismatch)
            newBC <- newBC[newBC %in% Barcodes]
        else
            stop("Input idents with mismatched barcodes! Set allow_mismatch as TRUE!")
    }

    # standardize input data
    object@meta.data$CB <- Barcodes
    if (class(new_ident) == 'data.frame') {
        df_temp <- new_ident[newBC, 1, drop = FALSE]
        df_temp$CB <- newBC
        message("Only add your first column of input data.frame")
    }
    if (class(new_ident) == 'factor') {
        df_temp <- data.frame(
            temp1 = as.character(unname(new_ident[newBC])),
            CB = newBC,
            row.names = newBC
            )
        colnames(df_temp) <- c(newSlotName, "CB")
    }

    # data completion
    if (nrow(df_temp) < length(Barcodes)) {
        Barcode_residual <- Barcodes[!Barcodes %in% newBC]
        if (!Replace) {
            df_temp_2 <- data.frame(
                temp1 = others_name,
                CB = Barcode_residual,
                row.names = Barcode_residual
                )
        } else {
            index_temp <- match(Barcode_residual, names(object@active.ident))
            df_temp_2 <- data.frame(
                temp1 = as.character(unname(object@active.ident[index_temp])),
                CB = Barcode_residual,
                row.names = Barcode_residual
                )
        }
        colnames(df_temp_2) <- colnames(df_temp)
        df_temp <- rbind(df_temp, df_temp_2)
    }

    # add meta slot
    object@meta.data <- merge(
        object@meta.data, 
        df_temp,
        by = "CB"
    )
    rownames(object@meta.data) <- object@meta.data$CB
    object@meta.data$CB <- NULL
    
    object@meta.data <- object@meta.data[Barcodes, ]
    return(object)    
}

#########################################
# Pseudotime
library(monocle)

my_create_monocle <- function(seurat_obj, idents = NULL, assay = "RNA") {
    if (!is.null(idents))
        seurat_obj <- subset(seurat_obj, ident = idents)
    
    expr_matrix <- seurat_obj@assays[[assay]]@counts
    p_data <- seurat_obj@meta.data 
    f_data <- data.frame(gene_short_name = row.names(expr_matrix), row.names = row.names(expr_matrix))
    
    # create cds object
    pd <- new('AnnotatedDataFrame', data = p_data) 
    fd <- new('AnnotatedDataFrame', data = f_data)
    cds <- newCellDataSet(expr_matrix,
                          phenoData = pd,
                          featureData = fd,
                          expressionFamily = negbinomial.size()
    )
    
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    return(cds)
}

my_featureSelect_cds <-  function(cds_obj, method = c("seurat", 'monocle', 'diffgenes', 'dpFeature'), seurat_obj = NULL, ...) {
    method <- match.arg(method)
    method <- method[1]
    
    switch(
        which(c("seurat", 'monocle', 'diffgenes', 'dpFeature') == method),
        message("seurat_obj should be assigned"),
        message("... should be FilterCondition, like FilterCondition = 'mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit'"),
        message("seurat_obj should be assigned and ... should be FilterCondition, like FilterCondition = 'p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)'"),
        message("... should involve 1)DefaultMeta, 2)qval_threshold (0.1, default)")
    )
    
    selected_FUN <- switch(
        which(c("seurat", 'monocle', 'diffgenes', 'dpFeature') == method),
        function() VariableFeatures(seurat_obj),
        function(FilterCondition) {
            if (!exists("FilterCondition"))
                stop("Assign FilterCondition please if use method 'monocle'!")
            disp_table <- dispersionTable(cds_obj)
            eval(parse(text = paste0("disp.genes <- unique(subset(disp_table, ", FilterCondition, ")$gene_id)")))
            disp.genes
        },
        function(FilterCondition) {
            if (is.null(seurat_obj))
                stop("Assign seurat_obj please if use method 'diffgenes'!")
            if (!exists("FilterCondition"))
                stop("Assign FilterCondition please if use method 'diffgenes'!")
            if (!exists("deg.cluster"))
                deg.cluster <<- FindAllMarkers(seurat_obj)
            eval(parse(text = paste0("express_genes <- unique(subset(deg.cluster, ", FilterCondition, ")$gene)")))
            express_genes
        },
        function(DefaultMeta, qval_threshold = 0.1, ntop = NULL) {
            if (!exists("DefaultMeta"))
                stop("Assign DefaultMeta please if use method 'dpFeature'!")
            cds_obj <- detectGenes(cds_obj, min_expr = 0.1)
            expressed_genes <- row.names(subset(fData(cds_obj), num_cells_expressed >= 10))
            diff <- differentialGeneTest(cds_obj[expressed_genes, ], fullModelFormulaStr = paste0("~", DefaultMeta), cores = 1) 
            deg <- subset(diff, qval < qval_threshold)
            deg <- deg[order(deg$qval, decreasing = F), ]
            ordergene <- rownames(deg)
            if (is.integer(ntop)) {
                if(ntop < length(ordergene))
                    ordergene[1:ntop]
            } else {
                ordergene
            }
        }
    )
    
    selected_Features <- selected_FUN(...)
    cds_res <- setOrderingFilter(cds_obj, selected_Features)
    print(plot_ordering_genes(cds_res))
    return(cds_res)
}

my_process_cds <- function(cds_obj, root_state = NULL) {
    if (is.null(root_state)) {
        cds <- cds_obj
        cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
        cds <- orderCells(cds)
    } else {
        cds <- orderCells(cds_obj, root_state = root_state)
    }
    return(cds)
}

my_plotPseudo <- function(
    cds, 
    seurat_obj, 
    theta = 0, 
    color_by = "State",
    orig.ident = "orig.ident", 
    type = c("density", "hist"), 
    size = 1, 
    show_backbone = TRUE, 
    color_legend.position = NULL, 
    color_guide_legend_direction = "horizontal",
    color_guide_legend_nrow = NULL,
    color_guide_legend_ncol = NULL,
    ...
) {
    type <- match.arg(type)
    type <- type[1]
    
    dimPseudotime <- cds@reducedDimS
    theta0 <- theta / 180 * pi # theta: degree; theta0:radian
    wMtx <- sin(theta0) * matrix(c(0,1,-1,0),2,2) + cos(theta0) * matrix(c(1,0,0,1),2,2) # Rotation matrix
    dimPseudotime1 <- wMtx %*% dimPseudotime
    dimPseudotime2 <- t(dimPseudotime1)
    colnames(dimPseudotime2) <- c("x", "y")
    dimPseudotime2 <- as.data.frame(dimPseudotime2)
    dimPseudotime2$orig.ident <- sapply(rownames(dimPseudotime2), 
                                        function(x)
                                            seurat_obj@meta.data[x, orig.ident]
    )
    packageExist <- require(ggplot2) && require(cowplot)
    if (!packageExist)
        stop("Packages ggplot2 and cowplot should be installed!")
    p1 <- plot_cell_trajectory(cds, color_by = color_by, size = size, show_backbone = show_backbone, theta = theta, ...) + 
        theme_classic()
    if (is.null(color_legend.position))
        p1 <- p1 + theme(legend.position = 'none')
    else
        p1 <- p1 + 
        theme(
            legend.position = color_legend.position,
            legend.title = element_blank(),
            legend.key.height = grid::unit(0.35, "in"), 
            legend.key = element_blank(),
            legend.background = element_rect(fill = "transparent")
        ) + 
        guides(
            color = guide_legend(
                direction = color_guide_legend_direction, 
                nrow = color_guide_legend_nrow, 
                ncol = color_guide_legend_ncol
            )
        )
    p4 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = size/4, show_backbone = FALSE, theta = theta, ...) + 
        labs(x = NULL, y = NULL) + theme_void() + theme(legend.position = 'none') 
    p2_tmp <- ggplot(dimPseudotime2, aes(x = x, color = orig.ident)) + 
        geom_density(mapping = aes(x = x, color = orig.ident, linetype = orig.ident), alpha = 0.1, show.legend = T) +
        theme_classic()
    p2 <- ggplot(dimPseudotime2, aes(x = x, color = orig.ident)) + 
        geom_density(mapping = aes(x = x, color = orig.ident, linetype = orig.ident), alpha = 0.1, show.legend = F) + 
        theme_classic() + xlab("Component 1") + ylab("Density") + 
        theme(
            axis.title.x = element_text(),
            axis.title.y = element_text(),
            axis.text.x = element_text(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    p3 <- ggplot(dimPseudotime2, aes(y = y, color = orig.ident)) + 
        geom_density(mapping = aes(y = y, color = orig.ident, linetype = orig.ident), alpha = 0.1, show.legend = F) + 
        theme_classic() + xlab("Density") + ylab("Component 2") + 
        scale_y_continuous(position = "right") +
        theme(
            axis.title.x = element_text(),
            axis.title.y = element_text(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(),
            axis.ticks = element_blank()
        )
    p22_tmp <- ggplot(dimPseudotime2, aes(x = x, color = orig.ident)) + 
        geom_freqpoly(mapping = aes(x = x, color = orig.ident, linetype = orig.ident), show.legend = T) +
        theme_classic()
    p22 <- ggplot(dimPseudotime2, aes(x = x, color = orig.ident)) + 
        geom_freqpoly(mapping = aes(x = x, color = orig.ident, linetype = orig.ident), show.legend = F) + 
        theme_classic() + xlab("Component 1") + ylab("Density") + 
        theme(
            axis.title.x = element_text(),
            axis.title.y = element_text(),
            axis.text.x = element_text(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()
        )
    p33 <- ggplot(dimPseudotime2, aes(y = y, color = orig.ident)) + 
        geom_freqpoly(mapping = aes(y = y, color = orig.ident, linetype = orig.ident), show.legend = F) + 
        theme_classic() + xlab("Density") + ylab("Component 2") + 
        scale_y_continuous(position = "right") +
        theme(
            axis.title.x = element_text(),
            axis.title.y = element_text(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(),
            axis.ticks = element_blank()
        )
    p23 <- plot_grid(NULL, get_legend(p2_tmp), ncol = 1)
    p2233 <- plot_grid(NULL, get_legend(p22_tmp), ncol = 1)

    pp <- if (type == 'density')
        plot_grid(p1, p3, p2, p23, ncol = 2, nrow = 2, align = "vh", rel_widths = c(4,1), rel_heights = c(4,1))
    else
        plot_grid(p1, p33, p22, p2233, ncol = 2, nrow = 2, align = "vh", rel_widths = c(4,1), rel_heights = c(4,1))
    
    return(pp)
}

my_AddSeuratPseudo <- function(
    seurat_obj,
    cds_obj, 
    reduction_key = 'Monocle'
) {
    # @reductions
    if (reduction_key %in% names(seurat_obj@reductions))
        warning(reduction_key, " has existed in seurat_obj@reductions! It will be covered.")
    seurat_obj@reductions[[reduction_key]] <- new('DimReduc')
    seurat_obj@reductions[[reduction_key]]@cell.embeddings <- t(cds_obj@reducedDimS)
    colnames(seurat_obj@reductions[[reduction_key]]@cell.embeddings) <- paste(reduction_key, 1:2, sep = "_")
    seurat_obj@reductions[[reduction_key]]@key <- paste0(reduction_key, "_")
    message("Add $", reduction_key, " into @reduction!")

    # @meta.data
    for (i in c('State', 'Pseudotime')) {
        if (i %in% names(cds_obj@phenoData@data)) {
            meta_key <- paste(reduction_key, i, sep = "_")
            seurat_obj@meta.data[[meta_key]] <- cds_obj@phenoData@data[[i]]
            message("Add $", meta_key, " into @meta.data!")
        }
    }

    return(seurat_obj)
}

subset_cds <- function(cds, ...) {
    cds2 <- cds[, ...]
    cds2@reducedDimS <- cds2@reducedDimS[, colnames(cds2@assayData$exprs)]
    return(cds2)
}


#########################################
# TCR/BCR
# require(immunarch)

plot_usage <- function(
    obj, 
    only_group = NULL, 
    gene_stats = c("musmus.trav", "musmus.traj", "musmus.trbv","musmus.trbj","musmus.trbd")[1], 
    rank_col = "all", 
    grid = TRUE, 
    title = NULL, 
    norm = TRUE,
    change_colname = NULL,
    return_plot = TRUE,
    return_res = FALSE,
    ...
) {
    # change_colname is a vector 
    # that refers to post-alternative names with custom orders, 
    # with name that refers to pre-alternative names
    if (!require(immunarch)) stop("Not install R package immunarch!")
    
    obj_data <- if(is.null(only_group)) obj$data else obj$data[[only_group]]
    imm_gu <- geneUsage(obj_data, gene_stats, .norm = norm, .ambig = "exc")
    imm_gu[is.na(imm_gu)] <- 0

    if (!is.null(rank_col)) {
        rank_method <- 
            if (is.null(rank_col)) 
                NULL
            else if (rank_col == "all") 
                function(x) mean(x[!is.na(x)])
            else if (rank_col %in% colnames(imm_gu)) 
                function(x) x[colnames(imm_gu)[-1] == rank_col]
            else if (rank_col %in% change_colname) 
                function(x) x[sapply(colnames(imm_gu)[-1], function(x) change_colname[x]) == rank_col]
            else
                stop("ERROR rank_col! if only_group is set as non-NULL, avoid this parameters (NULL) or set it as 'all' (default).")

        mean_tmp <- apply(imm_gu[, -1], 1, rank_method)    
        index_tmp <- order(mean_tmp, decreasing = T)
        imm_gu <- imm_gu[index_tmp, ]
        imm_gu[,1] <- factor(imm_gu[, 1, drop = T], levels = imm_gu[, 1, drop = T])
    }

    if (is.null(only_group))
        if (!is.null(change_colname)) {
            colnames(imm_gu) <- c(
                colnames(imm_gu)[1], 
                sapply(
                    colnames(imm_gu)[-1], 
                    function(x) 
                        change_colname[x]
                )
            )
            imm_gu <- imm_gu[, c(colnames(imm_gu)[1], change_colname)]
        }

    pp <- vis(imm_gu, .grid = grid, .title = title, ...)

    if (return_plot & return_res) print(pp)
    if (return_plot & !return_res) return(pp)
    if (return_res) return(imm_gu)
}

# TcrBVplot <- plot_usage(
#     TCR_meta, 
#     gene_stats = "musmus.trbv",
#     rank_col = 'NonCryoTCR_1wk', 
#     title = "Mus.TcrBV",
#     change_colname = colnames_alter,
#     return_plot = T,
#     return_res = F,
#     .add.layer = theme(axis.text.x = element_text(vjust = 0.5, size = 10), 
#                        axis.title.x = element_blank(),
#                        plot.title = element_text(hjust = 0.5)
#                        )
# )

plot_usage_multigene <- function(
    obj, 
    only_group = NULL, 
    gene_stats = c("musmus.trav", "musmus.traj", "musmus.trbv","musmus.trbj","musmus.trbd")[1:2], 
    rank_col = "all", 
    grid = TRUE, 
    title = NULL, 
    norm = TRUE,
    change_colname = NULL,
    return_plot = TRUE,
    return_res = FALSE,
    ...
) {
    if (length(gene_stats) <= 1)
        stop("Need more than 2 gene_stats! Use function plot_usage when only 1 gene_stat.")

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
                norm = norm,
                change_colname = change_colname,
                return_plot = FALSE,
                return_res = TRUE
            )
        }
    )

    imm_gu <- imm_gu_list[[1]]
    for (i in 2:length(imm_gu_list))
        imm_gu <- rbind(imm_gu, imm_gu_list[[i]])

    pp <- vis(imm_gu, .grid = grid, .title = title, ...)

    if (return_plot & return_res) print(pp)
    if (return_plot & !return_res) return(pp)
    if (return_res) return(imm_gu)
}

# TcrAVJplot_multigroup <- plot_usage_multigene(
#     TCR_meta, 
#     gene_stats = c("musmus.trav", "musmus.traj"), 
#     rank_col = 'NonCryoTCR_1wk', 
#     title = "Mus.TcrA-VJ",
#     change_colname = colnames_alter,
#     .add.layer = theme(axis.text.x = element_text(vjust = 0.5, size = 5), 
#                        axis.title.x = element_blank(),
#                        plot.title = element_text(hjust = 0.5, size = 25)
#                        ),
#     return_res = F,
#     return_plot = T
# )
# TcrAVJplot_singlegroup <- plot_usage_multigene(
#     TCR_meta, 
#     only_group = "NonCryoTCR_1wk",
#     gene_stats = c("musmus.trav", "musmus.traj"), 
#     rank_col = 'all', 
#     title = "Mus.TcrA-VJ",
#     change_colname = colnames_alter,
#     .add.layer = theme(axis.text.x = element_text(vjust = 0.5, size = 5), 
#                        axis.title.x = element_blank(),
#                        plot.title = element_text(hjust = 0.5, size = 25)
#                        ),
#     return_res = F,
#     return_plot = T
# )

plot_topClone <- function(object, clone, ...) {
    obj <- object[[clone]]
    cells <- obj$barcode
    clonotypes <- paste(unique(obj$clonotype), collapse =",")
    # p <- DimPlot(Cryo_merge_T, reduction = 'tsne', group.by = 'T_CD8_CD8_Cluster', pt.size = 0.1, 
    #         split.by = 'orig.ident', cells.highlight = cells) + labs(title = clone, subtitle = clonotypes)
    DimPlot(..., cells.highlight = cells) + labs(title = clone, subtitle = clonotypes)
}

plotNetwork <- function(mtx, count_trim = NULL, grep_select = NULL) {
    # suppressPackageStartupMessages(require(igraph))
    if (!require(igraph)) stop("Not install R package igraph!")

    ct_nmtx <- mtx
    if (!is.null(count_trim)) {
        ct_nmtx[ct_nmtx <= count_trim] <- 0
        tmp <- apply(ct_nmtx,1,sum) 
        tmp_name <- which(tmp == 0)
        ct_nmtx <- ct_nmtx[-tmp_name, -tmp_name]
    }
    if (!is.null(grep_select)) {
        tmp_index <- grep(grep_select, rownames(ct_nmtx))
        ct_nmtx <- ct_nmtx[tmp_index, tmp_index]
    }
    
    index_0 <- which(apply(ct_nmtx, 1, sum) == 0)
    if (length(index_0) > 0)
        ct_nmtx <- ct_nmtx[-index_0, -index_0]
    
    graphCt <- graph_from_adjacency_matrix(ct_nmtx, mode = "upper", weighted = T, diag = F)
    E(graphCt)$width <- E(graphCt)$weight / 2
    deg <- degree(graphCt, mode = "all")
    V(graphCt)$location <- sub("_[^_]*_.*$", "", V(graphCt)$name)
    vcolor <- c("orange","red","lightblue","tomato")
    V(graphCt)$color <- vcolor[factor(V(graphCt)$location)]

    plot(
        graphCt,
        layout = layout_in_circle, 
        vertex.size = deg,
        vertex.label.cex = .5,
        vertex.label.dist = 1,
        edge.color = "gray50",
        edge.arrow.size = .4, 
        edge.curved = .1
    )
    legend(x = -1.5,
           y = 1.5,
           levels(factor(V(graphCt)$location)),
           pch = 21,
           col = "#777777", 
           pt.bg = vcolor
    )

#     return(sort(deg, decreasing = T))
}