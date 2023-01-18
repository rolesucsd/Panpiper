#################
# Title: checkm_bin-stats.R- program to load checkM results from one file and graph
# Author- Renee Oles
# Purpose: R Script to read the CheckM results: bin_stats.analyze.tsv and concatenated.tre
# Date- 11/12/2021
################
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
require(stringr)
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

# plot the log file with all three variables in the same same graph
require(reshape2)
data_long <- melt(checkm, id.vars=c("Sample"))  # convert to long format
data_long$value <- as.numeric(as.character(data_long$value))

require(ggplot2)
ggplot(data_long, aes(fill="fill", x=value))+
  geom_histogram(bins=20)+
  theme_bw()+
  ggtitle("CheckM")+
  xlab("Isolates")+
  ylab("")+
  theme(axis.text.x = element_text(angle=90, hjust=1), legend.position="none")+
  scale_fill_manual(values=c("#62b587"))+
  facet_wrap(~variable, scale = "free")
ggsave(paste(output,".png",sep=""), height=7,width=15)

checkm$contigs <- as.numeric(as.character(checkm$contigs))
ggplot(checkm, aes(x=contigs))+
  #geom_bar(aes(fill = cut(`# contigs`, c(0, 1000, Inf))), stat="identity")+
  geom_histogram(bins=20)+
  theme_bw()+
  ggtitle("CheckM: Number of Contigs")+
  xlab("Isolates")+
  ylab("")+
  guides(fill=guide_legend(title="Cutoff"))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=c("grey","#664697"))
ggsave(paste(output,"#contigs.png",sep=""), height=7,width=15)

checkm$`N50 contigs` <- as.numeric(as.character(checkm$`N50 contigs`))
ggplot(checkm, aes(x=`N50 contigs`))+
  #geom_bar(aes(fill = cut(`N50 (contigs)`, c(0, 100000, Inf))), stat="identity")+
  geom_histogram(bins=20)+
  theme_bw()+
  ggtitle("CheckM: N50 (contigs)")+
  xlab("Isolates")+
  ylab("")+
  guides(fill=guide_legend(title="Cutoff"))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=c("#664697","grey"))
ggsave(paste(output,"N50.png",sep=""), height=7,width=15)

checkm$`Genome size` <- as.numeric(as.character(checkm$`Genome size`))
ggplot(checkm, aes(x=`Genome size`))+
  #geom_bar(aes(fill = cut(`Genome size`, c(0, (5205140-5205140*.5), (5205140+5205140*.5), Inf))), stat="identity")+
  geom_histogram(bins=20)+
  theme_bw()+
  ggtitle("CheckM: Genome size")+
  xlab("Isolates")+
  ylab("")+
  guides(fill=guide_legend(title="Cutoff"))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values=c("grey","#664697"))
  scale_fill_manual(values=c("#664697","grey","#664697"))
ggsave(paste(output,"genome_size.png",sep=""), height=7,width=15)

checkm$GC <- as.numeric(as.character(checkm$GC))
ggplot(checkm, aes(x=GC))+
  geom_histogram(bins=20)+
  #geom_bar(aes(fill = cut(GC, c(0, (0.4319-0.4319*.3), (0.4319+0.4319*.3), Inf))), stat="identity")+
  theme_bw()+
  ggtitle("CheckM: GC")+
  xlab("Isolates")+
  ylab("")+
  guides(fill=guide_legend(title="Cutoff"))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values=c("grey","#664697"))
  scale_fill_manual(values=c("#664697","grey","#664697"))
ggsave(paste(output,"gc.png",sep=""), height=7,width=15)

checkm$`Mean contig length` <- as.numeric(as.character(checkm$`Mean contig length`))
ggplot(checkm, aes(x=`Mean contig length`))+
  geom_histogram(bins=20)+
  #geom_bar(aes(fill = cut(`Mean contig length`, c(0, 10000, Inf))), stat="identity")+
  theme_bw()+
  ggtitle("CheckM: Mean contig length")+
  xlab("Isolates")+
  ylab("")+
  guides(fill=guide_legend(title="Cutoff"))+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=c("#664697","grey"))
ggsave(paste(output,"avg_contig_len.png",sep=""), height=7,width=15)
