# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
# Title: filter_kmers.R
# Author: Renee Oles
# Purpose: R Script to filter and correct p values by holm's correction
# Input: 
#	input: string
#		A string of the full path of the file to read in
#	output: string
#		A string of thefull path of the output file

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]

file <- read.delim(input, header = TRUE, check.names = FALSE)
file$bh_pv_adj <- p.adjust(file$`lrt-pvalue`, method = "BH", n = nrow(file))
file <- file[order(file$bh_pv_adj),]
write.table(file,output,quote=FALSE,col.names = TRUE,row.names = FALSE, sep="\t")
