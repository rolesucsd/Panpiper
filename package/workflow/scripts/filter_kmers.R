#################
# Title: filter_kmers.R
# Author- Renee Oles
# Purpose: R Script to filter and correct p values
# Date- 06/07/2022
################
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

file <- read.delim(input, header = TRUE, check.names = FALSE)
file$holm_pval <- p.adjust(file$`lrt-pvalue`, method = "holm", n = nrow(file))
file <- file[file$holm_pval < 0.05,]
write.table(file,output,quote=FALSE,col.names = TRUE,row.names = FALSE, sep="\t")
