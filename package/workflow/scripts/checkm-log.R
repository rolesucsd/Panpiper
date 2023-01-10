#################
# Check.R- program to load checkM results from one file and graph
# Author- Renee Oles
# Date- 7/2/2021
################

#install.packages("broom",repos = "https://cloud.r-project.org")
#install.packages("plyr",repos = "https://cloud.r-project.org")
#install.packages("tidyverse",repos = "https://cloud.r-project.org")
#install.packages("reshape2",repos = "https://cloud.r-project.org")
#install.packages("ggplot2",repos = "https://cloud.r-project.org")
require(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
inarg <- args
#meta <- "metadata/groups.txt"
output <- "results/Quality/checkm_log"

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
          input <- gsub("../Assembly/Shovill/","",input)
          input <- gsub("/lineage.log","",input)
          line <- gsub("filter", input, line)
          write(line, file=paste(output,".txt",sep=""), append=TRUE)
        }
      }
    }
  }
  close(con)
}
for(i in inarg){processFile(i)}

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
