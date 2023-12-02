source("/mnt/c/Users/PC/Documents/GitHub/My_scRNA_pipeline/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/211111/")
load("combined_analysis/Cryo_T_real.rda") #T_sub_SCT
setwd("~/Desktop/Cryo/")

X_umap_3end <- read.csv("T_2111_X_umap_3end.csv", row.names = 1)
T_2111_obs <- read.csv("T_2111_obs.csv", row.names = 1)

T_sub_SCT@meta.data <- T_2111_obs
rownames(X_umap_3end) <- rownames(T_sub_SCT@reductions$umap@cell.embeddings)
colnames(X_umap_3end) <- colnames(T_sub_SCT@reductions$umap@cell.embeddings)
X_umap_3end <- as.matrix(X_umap_3end)
T_sub_SCT@reductions$umap@cell.embeddings <- X_umap_3end

my_plotDim(T_sub_SCT, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'IFN_Final')
my_plotDim(T_sub_SCT, reduction = "umap", label = T, pt.size = 0.5, label.size = 5, group.by = 'MainCluster_unique_detail_fromSub')

my_DotPlot_split(T_sub_SCT, features = T_markers_tmp2, group.by = 'IFN_Final') + RotatedAxis()



#pathway activities visualization
source("~/Desktop/Github/My_scRNA_pipeline/Pathway_score.r")
gmt_file <- "/mnt/c/DeepTL_dev-main/pathway/GSVA/c5.go.v7.4.symbols1_ifn.gmt"
geneset <- getGmt(gmt_file)

exp_X <- T_sub_SCT@assays$SCT@data
rownames(exp_X) <- toupper(rownames(exp_X))
set_of_genes <- unique(unlist(lapply(geneset@.Data, function(x) x@geneIds)))
set_of_genes_overlap <- intersect(set_of_genes, rownames(exp_X))
exp_X <- exp_X[set_of_genes_overlap, ]

list_of_genes <- lapply(geneset@.Data, function(x) x@geneIds)
names(list_of_genes) <- names(geneset)
pathway_names <-names(list_of_genes)[grepl('RESPONSE',names(list_of_genes))&!grepl('GAMMA',names(list_of_genes))] 
list_of_genes <- list_of_genes[pathway_names]

T_sub_SCT <- geneset_score(T_sub_SCT, list_of_genes, 'average', 'data', TRUE, TRUE)

p1<-FeaturePlot(T_sub_SCT, features = pathway_names[1], pt.size = 1.5, cols = c('white', 'red'))
p2<-FeaturePlot(T_sub_SCT, features = pathway_names[2], pt.size = 1.5, cols = c('white', 'red'))
p3<-FeaturePlot(T_sub_SCT, features = pathway_names[3], pt.size = 1.5, cols = c('white', 'red'))
p4<-FeaturePlot(T_sub_SCT, features = pathway_names[4], pt.size = 1.5, cols = c('white', 'red'))
p5<-FeaturePlot(T_sub_SCT, features = pathway_names[5], pt.size = 1.5, cols = c('white', 'red'))

pathway_names

p_IFN_violin <- my_violin(T_sub_SCT, pathway_names, group.by = 'IFN_Final', mode = 'mtx') + theme(legend.position = 'none')
write.csv(T_sub_SCT@meta.data, file = "T_2111_meta_data_pathwayscore.csv")
