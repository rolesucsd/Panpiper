#################
# Title: checkm-log.R
# Author: Renee Oles
# Purpose: R Script to summarize the paths from the graph.py program
# Input: 
#   out_pref: string, required
#     A string of the file prefix for the output files
# Output: 
#   A data-frame, wide-format file, of each sample and summary of assembly details include completeness and contamination from log CheckM output
# TODO: remove hard-coding
################

require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
#meta <- "metadata/groups.txt"
# 
out_pref <- gsub("Quality.*", "", args[1])
output <- paste(out_pref,"Quality/CheckM/checkm_log",sep="")

file.create(paste(output,".txt",sep=""))
processFile = function(input) {
  con = file(input, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    else if (grepl("---", line)){
      while ( TRUE ) {
        line = readLines(con, n = 1)
        if ( length(line) == 0 ) {
          break
        }
        if(!grepl("Bin",line) && !grepl("---",line) & !grepl("\\[",line)){
          input <- gsub("*Quality/Assembly_filter/","",input)
          input <- gsub("/lineage.log","",input)
          line <- gsub("filter", input, line)
          write(line, file=paste(output,".txt",sep=""), append=TRUE)
        }
      }
    }
  }
  close(con)
}
for(i in args){processFile(i)}

# LOAD DATA
checkm_r <- data.frame()
df_c <- read.table(paste(output,".txt",sep=""), header=FALSE, sep="\t")
#df_c$old_sample_name <- as.numeric(substr(df_c$Sample,1,9))
colnames(df_c)[1] <- "first"
df_c <- separate(df_c, first,sep=("\\s+"), into=c("Blank","Bin_ID","Marker_lineage","Marker_lineage_ID","genomes","markers","marker_sets","0","1","2","3","4","5+","Completeness","Contamination","Strain_heterogeneity"))
df_c <- df_c[,c(2,3,14,15,16)]
df_c$Contamination <- as.numeric(df_c$Contamination)
df_c$Completeness <- as.numeric(df_c$Completeness)
df_c$Strain_heterogeneity <- as.numeric(df_c$Strain_heterogeneity)
df_c[df_c$Strain_heterogeneity == 100 & df_c$Contamination == 100, 4] <- 0
write.table(df_c, paste(output,".txt",sep=""), quote=FALSE,row.names = FALSE,col.names = TRUE,sep="\t")
write.table(df_c[df_c$Contamination < 5 & df_c$Completeness > 95 ,1], paste(output,"_filter.txt",sep=""), quote=FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
