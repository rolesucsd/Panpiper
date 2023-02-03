# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Title: AMR.R
# Author: Renee Oles
# Purpose: R Script to summarize the paths from the graph.py program
# Input: 
# 	amr: string, required
#		A string of the filename of the AMR matrix to use to make the AMR matrix
#	outfile: string, required
# 		A string of the file prefix for the heatmap and edited wide format matrix
# Output: 
# 		A matrix, wide-format file, of each sample and whether (1) or not (0) it has the given amr gene
#		A heatmap, png file file, of the presence/absence of the amr genes per sample


# Libraries
library(tidyverse)
library(RColorBrewer)
library(pheatmap)

# Input
args <- commandArgs(trailingOnly = TRUE)
amr <- read.delim(args[1])

# Output
outfile <- paste(args[2],"amr_wide.txt",sep="/")
outpng <- paste(args[2],"amr.png",sep="/")

# Convert from long to wide based off Sample
amr_wide <- unique(amr[,c(2,19)])
amr_wide$pr <- 1
amr_wide <- reshape(amr_wide, idvar = "Sample", timevar = "Gene.symbol", direction = "wide")
names(amr_wide) <- sub('^pr.', '', names(amr_wide))
amr_wide[is.na(amr_wide)] <- 0
amr_wide$Sample <- sub(".txt","",amr_wide$Sample)
rownames(amr_wide) <- amr_wide[,1]
amr_wide <- amr_wide[,-1]
write.table(amr_wide,outfile,quote=F,row.names = T,sep="\t")

amr_wide_edit <-  amr_wide[, colSums(amr_wide != 0) > 10]
amr_wide_edit    <- amr_wide_edit[,order(colnames(amr_wide_edit))]

# Create heatmap
png(filename=outpng, units="in", width=7, height=15, res=300)
pheatmap(amr_wide_edit, show_rownames = F, color=colorRampPalette(c("white", "red"))(50), cluster_cols = FALSE)
dev.off()
