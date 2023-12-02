### from pyScenic
pyScenicDir <- "E://Cryo-TCR/server/210924/"

library(SCENIC)
library(SCopeLoomR)

pyScenicLoomFile <- file.path(pyScenicDir, "pyscenic_output.loom")
loom <- open_loom(pyScenicLoomFile, mode="r")

# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name='Regulons')
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)

exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
# clusterings <- get_clusterings_with_name(loom)

close_loom(loom)
###


pheatmap::pheatmap(regulonsAUC@assays@data$AUC,
                   color = colorRampPalette(c("blue","white","red"))(100),
                   show_rownames = T,
                   show_colnames = F)
