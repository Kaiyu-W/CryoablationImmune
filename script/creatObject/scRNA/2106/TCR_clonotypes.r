
# 1. filtered sample
setwd('E://Cryo-TCR/server/auto/')

cell_type_df <- read.csv('E://Cryo-TCR/server/auto/Cryo_merge_CD45_500_Cluster_T_sub.csv', row.names = 1)

TCR_clonotypes_celltypes <- function(Group = c('Cryo', 'NonCryo')[1]) {
    dat.file1 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/filtered_contig_annotations.csv')
    dat.file2 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/consensus_annotations.csv')
    dat.file3 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/clonotypes.csv')
    data1 <- read.csv(dat.file1, stringsAsFactors = F, head=T)
    data2 <- read.csv(dat.file2, stringsAsFactors = F, head=T)
    data3 <- read.csv(dat.file3, stringsAsFactors = F, head=T)
    
    data1$raw_consensus_id <- sub("consensus_", "consensus", data1$raw_consensus_id)
    data4 <- data1[, c('barcode', 'raw_clonotype_id', 'raw_consensus_id')]
    head(data4)
    data4$barcode <- paste(Group, data4$barcode, sep = "_")
    data4$cell_type <- sapply(data4$barcode, function(x) cell_type_df[x,1])
    data4 <- data4[!is.na(data4$cell_type), ]
    
    xxx <- tapply(data4$cell_type, data4$raw_clonotype_id, table)
    names_xxx <- names(xxx)[order(as.numeric(sub("clonotype", "", names(xxx))), decreasing = F)]
    result <- matrix(0, nrow = length(xxx), ncol = length(unique(data4$cell_type)), dimnames = list(names_xxx, unique(data4$cell_type)))
    
    for (i in 1:nrow(result)) {
        x_temp <- rownames(result)[i]
        xx_temp <- xxx[[x_temp]]
        for (j in names(xx_temp)) {
            result[i, j] <- unname(xx_temp[j])
        }
    }
    
    data5 <- data3[data3$clonotype_id %in% rownames(result), ]
    if (all(data5$clonotype_id == rownames(result))) {
        data6 <- cbind(data5, result)
    }
    
    return(data6)
}

Cryo_TCR_clonotypes <- TCR_clonotypes_celltypes(Group = 'Cryo')
NonCryo_TCR_clonotypes <- TCR_clonotypes_celltypes(Group = 'NonCryo')

head(Cryo_TCR_clonotypes)
head(NonCryo_TCR_clonotypes)

# write.csv(Cryo_TCR_clonotypes, file = 'Cryo_TCR_clonotypes_filter.csv', row.names = F, quote = F)
# write.csv(NonCryo_TCR_clonotypes, file = 'NonCryo_TCR_clonotypes_filter.csv', row.names = F, quote = F)

# pheatmap::pheatmap(as.matrix(Cryo_TCR_clonotypes[, 8:17]), cluster_rows = T, show_rownames = F, cluster_cols = T)

Cryo_TCR_clonotypes$aa_special <- sapply(
    Cryo_TCR_clonotypes$cdr3s_aa, 
    function(x) {
        if (x %in% NonCryo_TCR_clonotypes$cdr3s_aa) {
            index <- which(NonCryo_TCR_clonotypes$cdr3s_aa %in% x)
            other <- NonCryo_TCR_clonotypes$clonotype_id[index]
            other <- paste(paste("NonCryo", other, sep = "_"), collapse = ", ")
            other
        } else {
            "Cryo_specific"
        }
    }
)
NonCryo_TCR_clonotypes$aa_special <- sapply(
    NonCryo_TCR_clonotypes$cdr3s_aa, 
    function(x) {
        if (x %in% Cryo_TCR_clonotypes$cdr3s_aa) {
            index <- which(Cryo_TCR_clonotypes$cdr3s_aa %in% x)
            other <- Cryo_TCR_clonotypes$clonotype_id[index]
            other <- paste(paste("Cryo", other, sep = "_"), collapse = ", ")
            other
        } else {
            "NonCryo_specific"
        }
    }
)
Cryo_TCR_clonotypes$nt_special <- sapply(
    Cryo_TCR_clonotypes$cdr3s_nt, 
    function(x) {
        if (x %in% NonCryo_TCR_clonotypes$cdr3s_nt) {
            index <- which(NonCryo_TCR_clonotypes$cdr3s_nt %in% x)
            other <- NonCryo_TCR_clonotypes$clonotype_id[index]
            other <- paste(paste("NonCryo", other, sep = "_"), collapse = ", ")
            other
        } else {
            "Cryo_specific"
        }
    }
)
NonCryo_TCR_clonotypes$nt_special <- sapply(
    NonCryo_TCR_clonotypes$cdr3s_nt, 
    function(x) {
        if (x %in% Cryo_TCR_clonotypes$cdr3s_nt) {
            index <- which(Cryo_TCR_clonotypes$cdr3s_nt %in% x)
            other <- Cryo_TCR_clonotypes$clonotype_id[index]
            other <- paste(paste("Cryo", other, sep = "_"), collapse = ", ")
            other
        } else {
            "NonCryo_specific"
        }
    }
)

write.csv(Cryo_TCR_clonotypes, file = 'Cryo_TCR_clonotypes_filter.csv', row.names = F, quote = F)
write.csv(NonCryo_TCR_clonotypes, file = 'NonCryo_TCR_clonotypes_filter.csv', row.names = F, quote = F)

#
TCR_clonotypes_vdjtools <- function(Group = c('Cryo', 'NonCryo')[1]) {
    dat.file1 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/filtered_contig_annotations.csv')
    dat.file2 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/consensus_annotations.csv')
    dat.file3 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/clonotypes.csv')
    dat.file4 <- paste0('E://Cryo-TCR/data/TCR_Cryo_Cellranger/', Group, 'TCR/outs/airr_rearrangement.tsv')
    data1 <- read.csv(dat.file1, stringsAsFactors = F, head=T)
    data2 <- read.csv(dat.file2, stringsAsFactors = F, head=T)
    data3 <- read.csv(dat.file3, stringsAsFactors = F, head=T)
    data4 <- read.table(dat.file4, stringsAsFactors = F, head=T, sep = "\t")
    
    temp1 <- data1[, c('cdr3', 'cdr3_nt', 'v_gene', 'd_gene', 'j_gene', 'raw_consensus_id', 'raw_clonotype_id', 'contig_id')]
    temp2 <- data4[, c('sequence_id', 'v_sequence_end', 'd_sequence_start', 'd_sequence_end', 'j_sequence_start')]
    head(res)
    head(data1)
    head(data2)
    head(data3)
    str(data4)
    duplicated(data1[,c('exact_subclonotype_id','raw_consensus_id')])
}
################################################################################

# 1.1) Load the package into R:
library(immunarch)

# 1.2) Replace with the path to your processed 10x data or to the clonotypes file
file_path_Cryo <- "E://Cryo-TCR/data/TCR_Cryo_Cellranger_output/CryoTCR/"
file_path_NonCryo <- "E://Cryo-TCR/data/TCR_Cryo_Cellranger_output/NonCryoTCR/"

# 1.3) Load 10x data with repLoad
immdata_Cryo <- repLoad(file_path_Cryo)
immdata_NonCryo <- repLoad(file_path_NonCryo)

repOverlap(immdata_Cryo$data) %>% vis()
geneUsage(immdata_Cryo$data[[1]]) %>% vis()
repDiversity(immdata_Cryo$data) %>% vis(.meta = immdata_Cryo$meta) 

