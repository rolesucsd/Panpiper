#################
# Title: AMR.R
# Author- Renee Oles
# Purpose: R Script to summarize the paths from the graph.py program
# Output: A heatmap of presence/absence in an operon
# Date- 9/23/2022
################

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
write.table(amr_wide,output,quote=F,row.names = T,sep="\t")

amr_wide_edit <-  amr_wide[, colSums(amr_wide != 0) > 10]
amr_wide_edit    <- amr_wide_edit[,order(colnames(amr_wide_edit))]

# Create heatmpa
png(filename=output, units="in", width=7, height=15, res=300)
pheatmap(amr_wide_edit, fontsize_row=1, color=colorRampPalette(c("white", "red"))(50), cluster_cols = FALSE)
dev.off()
