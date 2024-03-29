args=commandArgs(T)
input_dir <- args[1]
output_dir <- input_dir

output_name <- 'merge_idxstats.tsv'

setwd(input_dir)
files_name <- dir()
files_name <- files_name[grep("_idxstats.tsv$", files_name)]
sample_id <- sub("_idxstats.tsv", "", files_name)

temp <- read.table(files_name[1], sep = "\t", header = T, row.names = 1)
res <- temp[, 1:2]
colnames(res) <- c("GeneLength", sample_id[1])

for (i in 2:length(files_name)) {
    temp <- read.table(files_name[i], sep = "\t", header = T, row.names = 1)
    temp <- temp[, 2, drop = F]
    colnames(temp) <- sample_id[i]

    res <- cbind(res, temp)
}

setwd(output_dir)
write.table(res, file = output_name, sep = '\t', row.names = T, col.names = T, quote = F)
message("Idxstats Merge OVER!")