source("E://Cryo-TCR/server/auto/utilities.r")
load(file = "C://celldex_reference/Database_form_celldex.rda") # reference data
setwd("E://Cryo-TCR/server/210924/")
load(file = 'Combined_analysis.rda')
load(file = 'Combined_analysis_T_CD4_8.rda')
load(file = 'Combined_analysis_Myeloid.rda')

Cryo_merge_temp <- my_AddMeta(Cryo_merge, T_CD4$T_CD4_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD4_T_CD4_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, T_CD8$T_CD8_MainCluster, Replace = T)
Idents(Cryo_merge_temp) <- 'T_CD8_T_CD8_MainCluster'
Cryo_merge_temp <- my_AddMeta(Cryo_merge_temp, Myeloid_sub$M_sub_MainCluster, Replace = T)
Cryo_merge_temp$MainCluster[Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster == 'T_CD8(EM_ex)'] <- 'T'
Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster[Cryo_merge_temp$Myeloid_sub_M_sub_MainCluster == 'T_CD8(EM_ex)'] <- 'CD8_Effector_Memory'
Idents(Cryo_merge_temp) <- 'Myeloid_sub_M_sub_MainCluster'

my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 1.5, label.size = 6, group.by = 'MainCluster', title = 'Cryo_merge MainCluster')
my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 0.5, label.size = 6, group.by = 'Myeloid_sub_M_sub_MainCluster', title = 'Cryo_merge Cluster')
my_plotDim(Cryo_merge_temp, reduction = "umap", label = T, pt.size = 0.5, label.size = 6, group.by = 'orig.ident2', title = 'Cryo_merge Cryo_NonCryo')

my_CountCluster(Cryo_merge_temp, group1 = 'MainCluster', group2 = 'orig.ident2')

my_plotDim(Myeloid_sub, reduction = "umap", label = T, pt.size = 1, group.by = 'M_sub_MainCluster', title = 'M_sub_raw Cryo_NonCryo')
# Myeloid visualization
Idents(Cryo_merge_temp) <- 'MainCluster'
Myeloid_sub_real <- subset(Cryo_merge_temp, ident = c('Myeloid', 'Mast'))
Myeloid_sub_real$M_sub_MainCluster <- 0
Myeloid_sub_real$M_sub_MainCluster <- sapply(rownames(Myeloid_sub_real@meta.data), function(x) Myeloid_sub$M_sub_MainCluster[x])
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, group.by = 'M_sub_MainCluster', title = 'M_sub_real Cryo_NonCryo')

# ReCluster
DefaultAssay(Myeloid_sub_real) <- 'RNA'
Myeloid_sub_real <- SCTransform(Myeloid_sub_real)
Myeloid_sub_real <- FindVariableFeatures(Myeloid_sub_real, selection.method = "vst", nfeatures = 2000)
Myeloid_sub_real <- RunPCA(Myeloid_sub_real)
Myeloid_sub_real <- RunUMAP(Myeloid_sub_real, dims = 1:10)
Myeloid_sub_real <- RunTSNE(Myeloid_sub_real, dims = 1:10)
Myeloid_sub_real <- FindNeighbors(Myeloid_sub_real, dims = 1:12)

# cluster
Myeloid_sub_real  <- FindClusters(Myeloid_sub_real, resolution = seq(1,10,1)/100)
Myeloid_sub_real  <- FindClusters(Myeloid_sub_real, resolution = seq(1,10,1)/10)
clustree(Myeloid_sub_real)

# my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident2', title = 'Myeloid_sub_real Cryo_NonCryo')
# my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = paste0('SCT_snn_res.0.',c(paste0("0", 1:9), 1:3)))
# my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = paste0('SCT_snn_res.0.', 1:9))

# select R=0.2
Idents(Myeloid_sub_real) <- 'SCT_snn_res.0.2'
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'Myeloid_sub_real R=0.2')
my_DotPlot_split(Myeloid_sub_real, features = Myeloid_marker) + RotatedAxis()
my_violin(Myeloid_sub_real, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx')

# subcluster for DC into pDC/cDC
Myeloid_sub_real <- FindSubCluster(Myeloid_sub_real, cluster = '3', graph.name = 'SCT_snn', resolution = 0.05, subcluster.name = 'sub3')
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'sub3', title = 'Myeloid_sub R=0.2 subR=0.05')
Idents(Myeloid_sub_real) <- 'sub3'
my_violin(Myeloid_sub_real, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx')
my_violin(Myeloid_sub_real, features = unlist(cDC_marker), pt.size = 0, mode = 'mtx')
DoHeatmap(Myeloid_sub_real, features = unlist(cDC_marker), cells = which(Myeloid_sub_real$sub3 %in% c("6", "3_0")))

# try to split cDC into more, failed
# Myeloid_sub_real <- FindSubCluster(Myeloid_sub_real, cluster = '4', graph.name = 'SCT_snn', resolution = 0.1, subcluster.name = 'sub4')
# my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'sub4')
# Idents(Myeloid_sub_real) <- 'sub4'
# my_violin(Myeloid_sub_real, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx')
# FeaturePlot(Myeloid_sub_real, c('Cd163', 'Mrc1', 'Arg1', 'Arg2'))
# FeaturePlot(Myeloid_sub_real, c("Adgre1", 'Cd86', 'Cd80', 'Nos2'))

# try to find Granulocyte, no
my_DotPlot_split(Myeloid_sub_real, features = invitrogen_Granulocyte) + RotatedAxis()
my_violin(Myeloid_sub_real, features = unlist(invitrogen_Granulocyte), pt.size = 0, mode = 'mtx')
DoHeatmap(Myeloid_sub_real, features = unlist(invitrogen_Granulocyte))

# so result:
Myeloid_sub_real$MainCluster <- as.factor(Myeloid_sub_real$sub3)
Myeloid_sub_real$MainCluster2 <- as.factor(Myeloid_sub_real$sub3)
MainCluster <- rbind(
    data.frame(x='C1',y=0),
    data.frame(x='C2',y=1),
    data.frame(x='C3',y=4),
    data.frame(x='C4',y='3_0'),
    data.frame(x='C5',y=6),
    data.frame(x='C6',y='3_1'),
    data.frame(x='C7',y=c(2,5))
    # data.frame(x='C7',y=2),
    # data.frame(x='C8',y=5)
)
MainCluster2 <- rbind(
    data.frame(x='Monocyte_S100a8/9+',y=0),
    data.frame(x='Monocyte',y=1),
    data.frame(x='Macrophage',y=4),
    data.frame(x='cDC_Itgax+',y='3_0'),
    data.frame(x='cDC_Itgax-',y=6),
    data.frame(x='pDC',y='3_1'),
    data.frame(x='Mast',y=c(2,5))
    # data.frame(x='Mast',y=2),
    # data.frame(x='Mast',y=5)
)
MainCluster <- MainCluster[order(MainCluster$y, decreasing = F),]
MainCluster2 <- MainCluster2[order(MainCluster2$y, decreasing = F),]
if (all(MainCluster$y == levels(Myeloid_sub_real$MainCluster))) {
    levels(Myeloid_sub_real$MainCluster) <- MainCluster$x
    levels(Myeloid_sub_real$MainCluster2) <- MainCluster2$x
    Myeloid_sub_real$MainCluster <- factor(Myeloid_sub_real$MainCluster, levels = paste0("C", 1:7))
    Myeloid_sub_real$MainCluster2 <- factor(Myeloid_sub_real$MainCluster2, levels = c("Monocyte_S100a8/9+","Monocyte","Macrophage","cDC_Itgax+","cDC_Itgax-","pDC","Mast"))
}
Idents(Myeloid_sub_real) <- 'MainCluster'
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster', title = 'Myeloid_sub Cluster')
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster2', title = 'Myeloid_sub Cluster')
my_plotDim(Myeloid_sub_real, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster2', title = 'Myeloid_sub Cluster')

MyeloidMarkers <- list(
    Monocytes = c("Cd14", "Fcgr3"),
    Macrophages = c("Itgam", "Adgre1"),
    "M1/cDC" = c('Cd86', 'Cd80'),
    M2 = 'Mrc1',
    cDC = 'Itgax',
    # DCs_activated = "Cd83",
    pDC = c('Bst2', 'Siglech'),
    Mast = c('Fcer1a', 'Fcer1g', 'Kit'),
    Neutrophils = c("S100a8", "S100a9", "Cdk5")
)
my_violin(Myeloid_sub_real, features = unlist(MyeloidMarkers), pt.size = 0, mode = 'mtx', group.by = 'MainCluster')
my_violin(Myeloid_sub_real, features = unlist(MyeloidMarkers), pt.size = 0, mode = 'mtx', group.by = 'MainCluster2')
my_DotPlot_split(Myeloid_sub_real, features = Myeloid_marker) + RotatedAxis()
my_DotPlot_split(Myeloid_sub_real, features = Myeloid_marker, group.by = "MainCluster2") + RotatedAxis()
my_violin(Myeloid_sub_real, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx', group.by = 'MainCluster2')
my_DotPlot_split(Myeloid_sub_real, features = c("Trem2", "Retnla","C1qa", "Egr2"), group.by = "MainCluster2", split.by = 'orig.ident2') + RotatedAxis()
M12_markers <- c('Cd163', 'Mrc1', "Trem2", "Retnla","C1qa", "Egr2")
DoHeatmap(Myeloid_sub_real,features = c("Cd14", "Fcgr3","Adgre1",'Cd86', 'Cd80', 'Nos2','Cd163', 'Mrc1'), group.by = 'MainCluster2', cells = grep("^Monocyte$", Myeloid_sub_real$MainCluster2))

Idents(Myeloid_sub_real) <- 'MainCluster2'
# 0.06
Myeloid_sub_real <- FindSubCluster(Myeloid_sub_real, cluster = 'Monocyte', graph.name = 'SCT_snn', resolution = 0.06, subcluster.name = 'Monocyte_sub')
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'Monocyte_sub', title = 'Monocyte_sub Cluster')
my_DotPlot_split(Myeloid_sub_real, features = Myeloid_marker, group.by = "Monocyte_sub") + RotatedAxis()
my_violin(Myeloid_sub_real, features = unlist(Myeloid_marker), pt.size = 0, mode = 'mtx', group.by = 'Monocyte_sub')
FeaturePlot(Myeloid_sub_real, "Ly6g")
FeaturePlot(Myeloid_sub_real, "Adgre1")
Myeloid_sub_real$Monocyte_sub <- as.factor(Myeloid_sub_real$Monocyte_sub)
levels(Myeloid_sub_real$Monocyte_sub) <- c("cDC_Itgax-", "cDC_Itgax+", "Macrophage_M1", "Mast", "Macrophage_M2", "Monocyte", "Monocyte_S100a8/9+", "pDC")

##################################
Idents(Myeloid_sub_real) <- 'MainCluster2'
# ReCluster
M_sub <- subset(Myeloid_sub_real, ident = c("Monocyte", "Macrophage"))
DefaultAssay(M_sub) <- 'RNA'
# M_sub <- SCTransform(M_sub)
M_sub <- NormalizeData(M_sub, normalization.method = "LogNormalize", scale.factor = 10000)
M_sub <- ScaleData(M_sub, assay = "RNA", features = rownames(M_sub))
M_sub <- FindVariableFeatures(M_sub, selection.method = "vst", nfeatures = 1000)
M_sub <- RunPCA(M_sub)
M_sub <- RunUMAP(M_sub, dims = 1:10)
M_sub <- RunTSNE(M_sub, dims = 1:10)
M_sub <- FindNeighbors(M_sub, dims = 1:12)
M_sub  <- FindClusters(M_sub, resolution = 0.6)

# select R=0.6
Idents(M_sub) <- 'RNA_snn_res.0.6'
my_plotDim(M_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, title = 'M_sub R=0.1')
##
Idents(Myeloid_sub_real) <- 'MainCluster2'
Myeloid_sub_real_tmp <- my_AddMeta(Myeloid_sub_real, M_sub$RNA_snn_res.0.6, Replace = T)
my_plotDim(Myeloid_sub_real_tmp, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'M_sub_RNA_snn_res.0.6')
my_DotPlot_split(Myeloid_sub_real_tmp, features = Myeloid_marker, group.by = 'M_sub_RNA_snn_res.0.6') + RotatedAxis()
##
my_plotDim(M_sub, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster2')
my_DotPlot_split(M_sub, features = Myeloid_marker) + RotatedAxis()
my_violin(M_sub, features = setdiff(unlist(Myeloid_marker),"Ly6g"), pt.size = 0, mode = 'mtx')
my_CountCluster(M_sub, group2 = 'MainCluster2')

# Cluster 2/6:M2; 0/8:M1; 1/3/4/5/7:Monocyte
M_sub$M_sub <- M_sub$RNA_snn_res.0.6
levels(M_sub$M_sub) <- c("M1","Mo","M2","Mo","Mo","Mo","M2","Mo","M1")
my_plotDim(M_sub, reduction = "tsne", label = T, pt.size = 1, label.size = 5, group.by = 'M_sub')
# Myeloid_sub_real$M_sub_M_sub <- NULL
Myeloid_sub_real <- my_AddMeta(Myeloid_sub_real, new_ident = M_sub$M_sub, Replace = T)
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'M_sub_M_sub')
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster')

# add meta
Myeloid_sub_real$MainCluster_new <- as.factor(Myeloid_sub_real$MainCluster)
levels(Myeloid_sub_real$MainCluster_new) <- c("C1", paste0("C",3:8))
Myeloid_sub_real$MainCluster_new <- unfactor(Myeloid_sub_real$MainCluster_new)
Myeloid_sub_real$MainCluster_new[Myeloid_sub_real$M_sub_M_sub == "M1"] <- "C2"
Myeloid_sub_real$MainCluster_new[Myeloid_sub_real$M_sub_M_sub == "M2"] <- "C3"
Myeloid_sub_real$MainCluster_new[Myeloid_sub_real$M_sub_M_sub == "Mo"] <- "C4"
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster_new')
Myeloid_sub_real$MainCluster2_new <- as.factor(Myeloid_sub_real$MainCluster_new)
levels(Myeloid_sub_real$MainCluster2_new) <- c("Monocyte_S100a8/9+", "Macrophage_M1", "Macrophage_M2", "Monocyte", "cDC_Itgax+", "cDC_Itgax-", "pDC", "Mast")
my_plotDim(Myeloid_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 5, group.by = 'MainCluster2_new')

##################################
#### save plot
# CD45 DimPlot
CD45plus_plot <- my_plotDim(Cryo_merge_temp, reduction = "umap", label = F, pt.size = 0.2, label.size = 6, group.by = 'MainCluster', title = 'CD45+ Cluster')
ggsave(filename = "CD45plus_dim.pdf", plot = CD45plus_plot)

# CD45 count
CD45plus_count <- as.data.frame(my_CountCluster(Cryo_merge_temp, group1 = 'MainCluster', group2 = 'orig.ident2')[-6, 3:4])
CD45plus_count$cell <- rownames(CD45plus_count)
CD45plus_count <- reshape::melt(CD45plus_count, id.vars = 'cell')
CD45plus_count$variable <- sub("_percentage$", "", CD45plus_count$variable)
CD45plus_count$variable <- ifelse(CD45plus_count$variable == "Cryo", "CA", "Non-CA")
CD45plus_count$variable <- factor(CD45plus_count$variable, levels = c("Non-CA", "CA"))
CD45plus_Plot <- ggplot(CD45plus_count, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(CD45plus_count)/2 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "CD45plus_count.pdf", plot = CD45plus_Plot)

# Myeloid DimPlot
Myeloid_plot <- my_plotDim(Myeloid_sub_real, reduction = "umap", label = F, pt.size = 0.5, label.size = 10, group.by = 'MainCluster_new', title = 'Myeloid Cluster')
ggsave(filename = "Myeloid_dim.pdf", plot = Myeloid_plot)

# Myeloid VlnPlot
Myeloid_Vln <- my_violin(Myeloid_sub_real, features = unlist(MyeloidMarkers), pt.size = 0, mode = 'mtx', group.by = 'MainCluster_new')
ggsave(filename = "Myeloid_vln.pdf", plot = Myeloid_Vln)

# Myeloid count
Myeloid_count <- as.data.frame(my_CountCluster(Myeloid_sub_real, group1 = 'MainCluster_new', group2 = 'orig.ident2')[-9, 1:2])
Myeloid_sum <- my_CountCluster(Cryo_merge_temp, group1 = 'MainCluster', group2 = 'orig.ident2')["sum", 1:2]
Myeloid_count[,1] <- Myeloid_count[,1] / Myeloid_sum[1]
Myeloid_count[,2] <- Myeloid_count[,2] / Myeloid_sum[2]
Myeloid_count[,1] <- Myeloid_count[,1] / sum(Myeloid_count[,1])
Myeloid_count[,2] <- Myeloid_count[,2] / sum(Myeloid_count[,2])
Myeloid_count
Myeloid_count$cell <- rownames(Myeloid_count)
Myeloid_count <- reshape::melt(Myeloid_count, id.vars = 'cell')
Myeloid_count$variable <- sub("_count$", "", Myeloid_count$variable)
Myeloid_count$variable <- ifelse(Myeloid_count$variable == "Cryo", "CA", "Non-CA")
Myeloid_count$variable <- factor(Myeloid_count$variable, levels = c("Non-CA", "CA"))
Myeloid_Plot <- ggplot(Myeloid_count, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(Myeloid_count)/2 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "Myeloid_count.pdf", plot = Myeloid_Plot)

# Myeloid annotation DimPlot
Myeloid_plot2 <- my_plotDim(Myeloid_sub_real, reduction = "umap", label = F, pt.size = 0.5, label.size = 10, group.by = 'MainCluster2_new', title = 'Myeloid Cluster')
ggsave(filename = "Myeloid_anno_dim.pdf", plot = Myeloid_plot2)

# Myeloid annotation VlnPlot
Myeloid_Vln2 <- my_violin(Myeloid_sub_real, features = unlist(MyeloidMarkers), pt.size = 0, mode = 'mtx', group.by = 'MainCluster2_new') + 
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.8),'cm'))
ggsave(filename = "Myeloid_anno_vln.pdf", plot = Myeloid_Vln2)

# Myeloid annotation count
Myeloid_count2 <- as.data.frame(my_CountCluster(Myeloid_sub_real, group1 = 'MainCluster2_new', group2 = 'orig.ident2')[-9, 1:2])
Myeloid_sum2 <- my_CountCluster(Cryo_merge_temp, group1 = 'MainCluster', group2 = 'orig.ident2')["sum", 1:2]
Myeloid_count2[,1] <- Myeloid_count2[,1] / Myeloid_sum2[1]
Myeloid_count2[,2] <- Myeloid_count2[,2] / Myeloid_sum2[2]
Myeloid_count2[,1] <- Myeloid_count2[,1] / sum(Myeloid_count2[,1])
Myeloid_count2[,2] <- Myeloid_count2[,2] / sum(Myeloid_count2[,2])
Myeloid_count2
Myeloid_count2$cell <- rownames(Myeloid_count2)
Myeloid_count2 <- reshape::melt(Myeloid_count2, id.vars = 'cell')
Myeloid_count2$variable <- sub("_count$", "", Myeloid_count2$variable)
Myeloid_count2$variable <- ifelse(Myeloid_count2$variable == "Cryo", "CA", "Non-CA")
Myeloid_count2$variable <- factor(Myeloid_count2$variable, levels = c("Non-CA", "CA"))
Myeloid_Plot2 <- ggplot(Myeloid_count2, aes(x = variable, y = value, fill = cell)) +
    geom_bar(width = 0.9, stat = "identity") +
    xlab("Group") +
    ylab("Proportion, %") +
    labs(fill="Cell Type") +
    scale_fill_manual(values=c("grey75", RColorBrewer::brewer.pal(nrow(Myeloid_count2)/2 - 1, "Paired"))) +
    theme_bw() +
    theme(legend.text=element_text(size=8), 
          axis.title.x=element_text(size=12), 
          axis.title.y=element_text(size=12), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(reverse = FALSE))
ggsave(filename = "Myeloid_anno_count.pdf", plot = Myeloid_Plot2)
#### 

ggsave(filename = "Myeloid_dim.png", plot = Myeloid_plot)
ggsave(filename = "Myeloid_anno_dim.png", plot = Myeloid_plot2)
ggsave(filename = "Myeloid_Vln.png", plot = Myeloid_Vln)
ggsave(filename = "Myeloid_anno_count.png", plot = Myeloid_Plot2)

# save(Myeloid_sub_real, file = 'Combined_analysis_Myeloid_real.rda')

###########################################
# pseudotime: Monocle3
# 
# library(monocle3)
# library(Seurat)
# library(SeuratWrappers)
# library(ggplot2)
# library(patchwork)
# library(magrittr)
# cds <- as.cell_data_set(Myeloid_sub_real, group.by = 'Monocyte_sub')
# cds <- cluster_cells(cds, reduction_method = 'UMAP')
# p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
# p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
# wrap_plots(p1, p2)
# 
# # integrated.sub <- subset(as.Seurat(cds), monocle3_partitions == 1)
# # cds <- as.cell_data_set(integrated.sub)
# cds <- learn_graph(cds)
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
# 
# max.avp <- which(integrated.sub$MainCluster=='C2')[200]
# max.avp <- rownames(integrated.sub@meta.data)[max.avp]
# cds <- order_cells(cds, root_cells = max.avp)
# plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
#            label_branch_points = FALSE, label_roots = F)
# 
# # Set the assay back as 'integrated'
# integrated.sub <- as.Seurat(cds, assay = "integrated")
# FeaturePlot(integrated.sub, "monocle3_pseudotime")

# detach("package:monocle3", unload = TRUE)


# pseudotime: Monocole2
# all Myeloid
library(monocle)

cds_all <- my_create_monocle(Myeloid_sub_real)
cds_all_seurat <- my_featureSelect_cds(cds_all, method = "seurat", seurat_obj = Myeloid_sub_real)
cds_all_monocle <- my_featureSelect_cds(cds_all, method = "monocle", seurat_obj = Myeloid_sub_real, FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit")
cds_all_diffgenes <- my_featureSelect_cds(cds_all, method = "diffgenes", seurat_obj = Myeloid_sub_real, FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 1 | avg_log2FC < -1)")
cds_all_dpFeature <- my_featureSelect_cds(cds_all, method = "dpFeature", seurat_obj = Myeloid_sub_real, DefaultMeta = "MainCluster_new")

# cds_dpFeature
cds_all_dpFeature <- my_process_cds(cds_all_dpFeature)
my_plotPseudo(cds_all_dpFeature, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_seurat
cds_all_seurat <- my_process_cds(cds_all_seurat)
my_plotPseudo(cds_all_seurat, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_monocle
cds_all_monocle <- my_process_cds(cds_all_monocle)
my_plotPseudo(cds_all_monocle, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_diffgenes
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes)
plot_cell_trajectory(cds_all_diffgenes, color_by = "State")
cds_all_diffgenes <- my_process_cds(cds_all_diffgenes, root_state = 2)
my_plotPseudo(cds_all_diffgenes, color_by = "Pseudotime", theta = -20, show_branch_points = F, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")
plot_cell_trajectory(cds_all_diffgenes, color_by = "MainCluster2_new", theta = -20, show_branch_points = F, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")


# pseudotime: Monocole2
# only Macro/Mono/DC
Idents(Myeloid_sub_real) <- 'MainCluster_new'
cds_sub <- my_create_monocle(Myeloid_sub_real, idents = paste0("C", 2:7))

cds_sub_seurat <- my_featureSelect_cds(cds_sub, method = "seurat", seurat_obj = Myeloid_sub_real)
cds_sub_monocle <- my_featureSelect_cds(cds_sub, method = "monocle", seurat_obj = Myeloid_sub_real, FilterCondition = "mean_expression >= 0.5 & dispersion_empirical >= 1 * dispersion_fit")
cds_sub_diffgenes <- my_featureSelect_cds(cds_sub, method = "diffgenes", seurat_obj = Myeloid_sub_real, FilterCondition = "p_val_adj<1e-6 & (avg_log2FC > 0.5 | avg_log2FC < -0.5)")
cds_sub_dpFeature <- my_featureSelect_cds(cds_sub, method = "dpFeature", seurat_obj = Myeloid_sub_real, DefaultMeta = "MainCluster_new", qval_threshold = 0.01)

# cds_dpFeature
cds_sub_dpFeature <- my_process_cds(cds_sub_dpFeature)
my_plotPseudo(cds_sub_dpFeature, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_seurat
cds_sub_seurat <- my_process_cds(cds_sub_seurat)
my_plotPseudo(cds_sub_seurat, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_monocle
cds_sub_monocle <- my_process_cds(cds_sub_monocle)
my_plotPseudo(cds_sub_monocle, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# cds_diffgenes
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes)
plot_cell_trajectory(cds_sub_diffgenes, color_by = "State", theta = -20, show_branch_points = F)
cds_sub_diffgenes <- my_process_cds(cds_sub_diffgenes, root_state = 2)
my_plotPseudo(cds_sub_diffgenes, color_by = "Pseudotime", theta = -20, show_branch_points = F, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")
plot_cell_trajectory(cds_sub_diffgenes, color_by = "MainCluster2_new", theta = -20, show_branch_points = F, seurat_obj = Myeloid_sub_real, orig.ident = "orig.ident2")

# #############
# pdata <- Biobase::pData(cds)
# s.cells <- subset(pdata, State=="7") %>% rownames()
# 
# keygenes <- head(ordergene,4)
# cds_subset <- cds[keygenes,]
# ##可视化：以state/celltype/pseudotime进行
# p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
# p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "cell_type")
# p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
# plotc <- p1|p2|p3
# 
# s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","S100A4")
# p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
# p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
# p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
# plotc <- p1|p2|p3
# 
# #这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
# #如果不设置，就会用所有基因来做它们与拟时间的相关性
# Time_diff <- differentialGeneTest(cds[ordergene, ], cores = 1, fullModelFormulaStr = "~sm.ns(Pseudotime)")
# Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] #把gene放前面，也可以不改
# 
# Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
# p <- plot_pseudotime_heatmap(cds[Time_genes, ], num_clusters = 4, show_rownames = T, return_heatmap = T)
# p
# ggsave("Time_heatmap.pdf", p, width = 5, height = 10)


Seurat_tmp <- as.Seurat(cds)
Myeloid_sub_real$Pseudotime <- Seurat_tmp$Pseudotime
Myeloid_sub_real$State <- unfactor(Seurat_tmp$State)
Myeloid_sub_real$State[is.na(Myeloid_sub_real$State)] <- 4
DimPlot(Myeloid_sub_real, reduction = 'umap', group.by = 'State', pt.size = 1)
DimPlot(Myeloid_sub_real, reduction = 'umap', group.by = 'Pseudotime', pt.size = 1) + theme(legend.position = 'none')

