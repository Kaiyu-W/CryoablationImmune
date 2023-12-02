library(Seurat)
setwd("E://Cryo-TCR/server/211111/")
load(file = "Cryo_T_analysis.rda")

setwd("E://Cryo-TCR/server/211111/TCR_for_immunarch/")
library(immunarch)
immdata_10x <- immunarch::repLoad(".", .mode = "paired")

repExplore(immdata_10x$data, .method = "lens", .col = "aa") %>% vis()
# repExplore(immdata_10x$data, "lens", .col = "aa") %>% vis() %>% fixVis()
repExplore(immdata_10x$data, .method = "count") %>% vis()
repExplore(immdata_10x$data, .method = "volume") %>% vis()
repClonality(immdata_10x$data, "homeo") %>% vis()
repOverlap(immdata_10x$data) %>% vis()
geneUsage(immdata_10x$data[[1]]) %>% vis()
repDiversity(immdata_10x$data) %>% vis(.by = "group", .meta = immdata_10x$meta)
repDiversity(immdata_10x$data) %>% vis(.by = "time", .meta = immdata_10x$meta)


imm_pr <- repClonality(immdata_10x$data, .method = "clonal.prop")
imm_pr
imm_top <- repClonality(immdata_10x$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
imm_top
imm_rare <- repClonality(immdata_10x$data, .method = "rare")
imm_rare
imm_hom <- repClonality(immdata_10x$data,
                        .method = "homeo",
                        .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1)
)
imm_hom

vis(imm_top) + vis(imm_top, .by = "group", .meta = immdata_10x$meta)
vis(imm_rare) + vis(imm_rare, .by = "group", .meta = immdata_10x$meta)
vis(imm_hom) + vis(imm_hom, .by = c("group", "time"), .meta = immdata_10x$meta)







repDiversity(immdata_10x$data, "raref", .verbose = F) %>% vis(.log = TRUE)
trackClonotypes(immdata_10x$data, list(1, 100), .col = "nt") %>% vis()
trackClonotypes(immdata_10x$data, list("CryoTCR_1wk", 10), .col = "aa+v") %>% vis()
kmers <- getKmers(immdata_10x$data, 10)
vis(kmers)


gene_stats()
imm_gu <- geneUsage(immdata_10x$data, c("musmus.ighj"), .norm = T)

vis(imm_gu)

vis(imm_gu, .grid = T)

vis(imm_gu, .by = "group", .meta = immdata_10x$meta, .plot = "box")

imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = F)
vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 5)

imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = F)
vis(imm_cl_pca, .plot = "clust")

imm_cl_pca2 <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .k = 3, .verbose = F)
vis(imm_cl_pca2)

p1 <- vis(spectratype(immdata_10x$data[[1]], .quant = "id", .col = "nt"))
p2 <- vis(spectratype(immdata_10x$data[[1]], .quant = "count", .col = "aa+v"))
p1 + p2

p1 <- repOverlap(immdata_10x$data) %>% vis()
p2 <- repDiversity(immdata_10x$data) %>% vis()
target <- c("CARAGYLRGFDYW;CQQYGSSPLTF", "CARATSFYYFHHW;CTSYTTRTTLIF", "CARDLSRGDYFPYFSYHMNVW;CQSDDTANHVIF", "CARGFDTNAFDIW;CTAWDDSLSGVVF", "CTREDYW;CMQTIQLRTF")
p3 <- trackClonotypes(immdata_10x$data, target, .col = "aa") %>% vis()
p1 + p2 + p3


barcodes <- c("AGTAGTCAGTGTACTC-1", "GGCGACTGTACCGAGA-1", "TTGAACGGTCACCTAA-1")
new_df <- select_barcodes(immdata_10x$data[[1]], barcodes)
new_df


exp_vol <- repExplore(immdata_10x$data, .method = "volume")
p1 <- vis(exp_vol, .by = c("group"), .meta = immdata_10x$meta)
p2 <- vis(exp_vol, .by = c("group", "time"), .meta = immdata_10x$meta)
p1 + p2
