library(Seurat)
library(SeuratDisk)

setwd("E://Cryo-TCR/server/210924/")
load(file = 'Combined_analysis_TNK_real.rda')

T_sub_real@commands <- list()
T_sub_real@assays$SCT <- NULL
# T_sub_real@assays$RNA@data <- matrix()
T_sub_real@assays$RNA@scale.data <- matrix()

file.remove("T_sub_real.h5Seurat")
file.remove("T_sub_real.h5ad")
SaveH5Seurat(T_sub_real, filename = "T_sub_real.h5Seurat")
Convert("T_sub_real.h5Seurat", dest = "h5ad")


################################################################################
#after python
source("E://Cryo-TCR/server/auto/utilities.r")

meta <- read.csv("TNK_python_out_meta.csv", row.names = 1)

all(rownames(meta) == rownames(T_sub_real@meta.data))
T_sub_real@meta.data$GAE_snn_leiden_0.2 <- as.factor(meta$GAE_snn_leiden_0.2)
T_sub_real@meta.data$GAE_snn_leiden_0.5 <- as.factor(meta$GAE_snn_leiden_0.5)
T_sub_real$RNA_louvain_0.3 <- T_sub_real$RNA_snn_res.0.3
levels(T_sub_real$RNA_louvain_0.3) <- paste0("C", as.numeric(levels(T_sub_real$RNA_louvain_0.3))+1)
levels(T_sub_real$GAE_snn_leiden_0.2) <- paste0("C", as.numeric(levels(T_sub_real$GAE_snn_leiden_0.2))+1)

Idents(T_sub_real) <- "RNA_0.3"
TNK_dim_se <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 8, title = 'Louvain', group.by = 'RNA_louvain_0.3')
TNK_TNK_dot_se <- my_DotPlot_split(T_sub_real, features = my_MarkersList[3:4]) + RotatedAxis()
TNK_anno_dim_se <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 8, title = "T/NK Clusters", group.by = "T_sub_MainCluster2")

Idents(T_sub_real) <- "GAE_snn_leiden_0.2"
TNK_dim_gt <- my_plotDim(T_sub_real, reduction = "umap", label = T, pt.size = 1, label.size = 8, title = 'GAE')
TNK_TNK_dot_gt <- my_DotPlot_split(T_sub_real, features = my_MarkersList[3:4]) + RotatedAxis()

ggsave(filename = 'louvain_res.png', TNK_dim_se + theme(legend.position = 'none'))
ggsave(filename = 'gae_res.png', TNK_dim_gt + theme(legend.position = 'none'))
ggsave(filename = 'anno_res.png', TNK_anno_dim_se + theme(legend.position = 'none'))


count_mtx <- table(T_sub_real@meta.data[,c('RNA_louvain_0.3','GAE_snn_leiden_0.2')])
heat_mtx <- matrix(count_mtx, nrow=nrow(count_mtx), ncol = ncol(count_mtx), dimnames = dimnames(count_mtx))
heat_mtx <- heat_mtx[c(1:2,4,5,3,6:9,11,10,12),]
png("louvai_gae_count.png")
pheatmap(heat_mtx, cluster_rows = F, cluster_cols = F, color = c('white','blue','black'), legend = F)
dev.off()

count_tmp <- my_CountCluster(T_sub_real, group1 = 'GAE_snn_leiden_0.2', group2 = 'orig.ident2')
count_ratio <- apply(count_tmp[1:13,2:1],2,function(x)x/sum(x)*100)
write.csv(count_ratio, file='anno_count_ratio.csv')


################################################################################
#pathway activities visualization
RWB_color <- grDevices::colorRampPalette(colors = c("blue","white","red"))(100)
gmt_file <- "C://DeepTL_dev-main/pathway/GSVA/c5.go.v7.4.symbols1_ifn.gmt"
geneset <- getGmt(gmt_file)
load(file = 'Combined_analysis_TNK_real.rda')
exp_X <- T_sub_real@assays$SCT@data
# exp_X <- T_sub_real@assays$SCT@scale.data
# exp_X <- as.matrix(exp_X)
rownames(exp_X) <- toupper(rownames(exp_X))
set_of_genes <- unique(unlist(lapply(geneset@.Data, function(x) x@geneIds)))
set_of_genes_overlap <- intersect(set_of_genes, rownames(exp_X))
exp_X <- exp_X[set_of_genes_overlap, ]

list_of_genes <- lapply(geneset@.Data, function(x) x@geneIds)
res_list <- list()
for (i in seq(list_of_genes)) {
    if (all(!rownames(exp_X) %in% list_of_genes[[i]])) {
        res_list[[i]] <- rep(0, ncol(exp_X))
    } else {
        tmp <- apply(exp_X, 2, function(x) {
            tapply(x, rownames(exp_X) %in% list_of_genes[[i]], mean)
        })
        res_list[[i]] <- tmp['TRUE',]
    }
}
names(res_list) <- names(geneset)
res_mtx <- as.matrix(data.frame(res_list))
head(res_mtx)

T_sub_real@meta.data <- cbind(T_sub_real@meta.data, res_mtx)

pathway_names <- colnames(res_mtx)[grepl('RESPONSE',colnames(res_mtx))&!grepl('GAMMA',colnames(res_mtx))] 
FeaturePlot(T_sub_real, features = pathway_names)

p1<-FeaturePlot(T_sub_real, features = pathway_names[1], pt.size = 1.5, cols = c('white', 'red'))
p2<-FeaturePlot(T_sub_real, features = pathway_names[2], pt.size = 1.5, cols = c('white', 'red'))
p3<-FeaturePlot(T_sub_real, features = pathway_names[3], pt.size = 1.5, cols = c('white', 'red'))
p4<-FeaturePlot(T_sub_real, features = pathway_names[4], pt.size = 1.5, cols = c('white', 'red'))
# FeaturePlot(T_sub_real, features = pathway_names[5], pt.size = 1.5, cols = c('white', 'red'))

pathway_names[1:4]

ggsave(filename = 'pathway1.pdf', plot = p1)
ggsave(filename = 'pathway2.pdf', plot = p2)
ggsave(filename = 'pathway3.pdf', plot = p3)
ggsave(filename = 'pathway4.pdf', plot = p4)