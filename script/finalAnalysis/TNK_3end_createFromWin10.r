source("/mnt/c/Users/PC/Documents/GitHub/My_scRNA_pipeline/utilities.r")
setwd("/mnt/e/Cryo-TCR/server/210924/")
load("Combined_analysis_TNK_real.rda")
load("T_sub_real_meta_220713.rda")
T_sub_real@meta.data <- T_meta_220713

T_sub_real$IFN_Final_Cluster[T_sub_real$IFN_Final == 'CD8_Tem_I-IFN'] <- 'C4'
T_sub_real$IFN_Final_Cluster[T_sub_real$IFN_Final == 'CD8_Tex_Proliferative'] <- 'C6'
T_sub_real$IFN_Final_Cluster[T_sub_real$IFN_Final == 'CD8_Tem'] <- 'C3'
T_sub_real$IFN_Final_Cluster[T_sub_real$IFN_Final == 'CD8_Naive'] <- 'C2'

setwd("~/Desktop/Cryo/")
saveRDS(T_sub_real, file = 'TNK_3end.rds')
