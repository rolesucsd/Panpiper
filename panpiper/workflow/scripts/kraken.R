# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
# Title: 
# Author: Renee Oles
# Purpose: Read and summarize kraken results
# Input: 
# Output: 

library(tidyverse)

# Input
args <- commandArgs(trailingOnly = TRUE)
kraken_files <- args[1:(length(args)-1)]
outpath <- args[length(args)]

kraken_contig <- data.frame()
kraken_ag <- data.frame()
for (i in 1:length(kraken_files)){
  f <- kraken_files[i]
  print(f)
  df <- read.delim(f,header=FALSE)
  df <- df[,c(1:4)]
  f <- sub(".out","",f)
  df$sample <- f
  df_sum <- aggregate(V4 ~ V3, data=df, FUN=sum)
  df_sum$V3 <- sapply(df_sum$V3, function(x){sub("(\\w+\\s+\\w+).*", "\\1", x)})
  df_sum$V3 <- sapply(df_sum$V3, function(x){gsub("\\(.*","",x)})
  df_sum$V3 <- sapply(df_sum$V3, function(x){gsub("Phocaeicola","Bacteroides",x)})
  df_sum <- aggregate(V4 ~ V3, data=df_sum, FUN=sum)
  df_sum$perc <- sapply(df_sum$V4, function(x){100*x/sum(df_sum$V4)})
  df_sum$sample <- f
  kraken_contig <- rbind(df, kraken_contig)
  kraken_ag <- rbind(df_sum, kraken_ag)
}

colnames(kraken_ag) <- c("Species_ID","length","perc","Bin_ID")
kraken_ag$Bin_ID <- gsub(".out","",kraken_ag$Bin_ID)
kraken_ag$old_sample_name <- substr(kraken_ag$Bin_ID,1,9)
write.table(kraken_ag, paste(outpath,"kraken_ag.txt",sep="/"),quote=FALSE,sep="\t",row.names=FALSE)
kraken_cont <- unique(kraken_ag[kraken_ag$perc >= 5 & kraken_ag$perc < 95,4])
write.table(kraken_cont,paste(outpath,"kraken_cut.txt",sep="/"),quote=FALSE,sep="\t",row.names=FALSE)

