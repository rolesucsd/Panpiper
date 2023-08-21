# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
pattern <- args[2]
output <- args[3]

# The unitig file to be filtered
file <- read.delim(input, header = TRUE, check.names = FALSE)
# The filtering threshold determined by the pattern program in pyseer
pattern <- read.delim(pattern, header=FALSE)
# The adjusted pvalue is made by subtracted the pvalue by the threshold
file$pv_adj <- file$`lrt-pvalue` - pattern[2,2]
# The unitigs are filtered by the threshold
file <-file[file$`lrt-pvalue` <= pattern[2,2],]
#file <-file[file$`lrt-pvalue` <= 0.05,]
file <- file[order(file$`lrt-pvalue`),]
write.table(file,output,quote=FALSE,col.names = TRUE,row.names = FALSE, sep="\t")
