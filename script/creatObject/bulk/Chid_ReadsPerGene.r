args=commandArgs(T)
input_dir <- args[1]
output_dir <- input_dir

id_ref <- args[2]

setwd(input_dir)

ref <- read.table(id_ref, sep = "\t", header = F, row.names = 1)
input <- read.table("ReadsPerGene_merge.tsv", sep = "\t", header = T, row.names = 1)

a <- sort(rownames(ref))
b <- sort(rownames(input))
ref1 <- ref[a, ]
input1 <- input[b, ]

output_temp <- as.data.frame(apply(input1, 2, function(x) {tapply(x, ref1, sum)}))

output <- output_temp
output$'M.Nonc.1' <- output$'M.Nonc.1' + output$'X2_M.Nonc.1'
output$'X2_M.Nonc.1' <- NULL

setwd(output_dir)
write.table(output, file = "genome_counts.tsv", sep = "\t", col.names = T, row.names = T, quote = F)

cat("genome_counts.tsv saved in", output_dir, "\n")
