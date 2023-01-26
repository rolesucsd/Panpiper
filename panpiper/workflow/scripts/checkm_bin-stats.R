# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
# Title: checkm_bin-stats.R- program to load checkM results from one file and graph
# Author- Renee Oles
# Purpose: R Script to read the CheckM results: bin_stats.analyze.tsv and concatenated.tre

require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
input <- args
out_pref <- gsub("Quality.*", "", args[1])
output <- paste(out_pref,"Quality/CheckM/checkm_stats",sep="")

checkm <- data.frame()
for (i in input){
  line <- read.delim(i, header=FALSE)
  i <- gsub("../Assembly/Shovill/","",i)
  i <- gsub("/storage/bin_stats.analyze.tsv","",i)
  line$Sample <- i  
  checkm <- rbind(checkm,line)
}
checkm2 <- str_split_fixed(checkm$V2, ",", 15)
checkm2<-apply(checkm2,2,function(x) gsub("'","",as.character(x)))
checkm2<-apply(checkm2,2,function(x) gsub("\\{","",as.character(x)))
checkm2<-apply(checkm2,2,function(x) gsub("\\}","",as.character(x)))
checkm2<-apply(checkm2,2,function(x) sub(".*\\:","",as.character(x)))
checkm2 <- as.data.frame(checkm2)
checkm <- cbind(checkm[,3],checkm2)
colnames(checkm) <- c("Sample", "GC", "GC std", "Genome size",  "ambiguous bases",  "scaffolds",  "contigs",  "Longest scaffold",
                      "Longest contig", "N50 scaffolds",  "N50 contigs",  "Mean scaffold length", "Mean contig length",
                      "Coding density", "Translation table",  "predicted genes")
checkm$Sample <- basename(checkm$Sample)
write.table(checkm,paste(output,".txt",sep=""),quote=FALSE,col.names = TRUE,row.names = FALSE, sep="\t")