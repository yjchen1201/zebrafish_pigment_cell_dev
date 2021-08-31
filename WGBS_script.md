### Tool version
#cutadapt v1.10: Used to trim adapter reads 
#samtools version 1.3.1: Used to sort and downsample bam files for downstream processing
#bismark v0.16.1: Used to align WGBS reads to the genome and call CpG methylation
#bwa v0.7.15: bwa-mem used to align ATAC-seq reads to the genome
#macs2 v2.1.1: Used to call peaks for ATAC-seq
#STAR v2.5.1b: Used to align RNA-sequencing data
#stringtie v1.3.3: Used to annotate aligned reads and quantify transcripts. 
#picard v2.8.1: Used to mark and remove duplicate reads
#bedtools v2.27.1: closest, intersect, shuffle command used as specified in manuscript 
#meme v5.0.3: Used AME command to find motif enrichment and FIMO command to scan for motif presence 
#R v3.3.0: Used to analyze sequencing library results
#DSS v2.14.0 :Used to call differentially methylated regions
#DESeq2 v1.12.4: Used to call differentially expressed genes
#DiffBind v 2.2.12: Used to call differentially accessible regions/peaks
#Metascape v3.0: Used for GO enrichment analysis

### Analysis
```{bash}
###################################################### {bash script} ######################################################
## Example of BS-seq pipeline#

# Adapter trimming
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1.fastq.gz" "SAMPLE_R2.fastq.gz"

# Reads alignment
bismark --bowtie2 --skip 9 -N 1 -L 28 --score_min L,0,-0.6 "reference genome path" -1 "SAMPLE_R1_trimmed.fastq.gz" -2 "SAMPLE_R2_trimmed.fastq.gz"

# SAM->BAM, and only include reads mapped to the input bed region 
samtools view -b -L "danRer10_genome_lite.size.bed" "SAMPLE_R1_trimmed_bismark_bt2_pe.bam" > "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam"

# sort
samtools sort -o "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam" "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam"  

# Remove duplicate
java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=True REMOVE_DUPLICATES=True I="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam" O="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam" M="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.dup_removed.txt" ASSUME_SORTED=true

# sort bam file by read name
samtools sort -n -o "SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam" "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam"

# Extract the methylation information from the Bismark alignment output
bismark_methylation_extractor -p --no_overlap --comprehensive --gzip --no_header --bedgraph --counts --ignore 10 --ignore_r2 10 --report -o ./ --genome_folder <path to genome folder> "SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam"

##Bismark ON LAMBDA###
bismark --quiet <path to lambda genome> -1 "SAMPLE_R1.fastq.gz" -2 "SAMPLE_R2.fastq.gz"


###convert bismark coverage output into DSS input file format
python3 convert_DSS.py Merged_ALLRUNS_Iri_Rep1_trimmed_bismark_resorted.rmdup.bismark.cov
```
```{R}
###################################################### {R analysis} ######################################################
##QC individual WGBS Replicates##
setwd(<WGBS_folder>)
library(ggplot2)
library(RColorBrewer)
mypalette <- brewer.pal(12,"Paired")
mypalette2 <- c("#a8e6cf","#ffd3b6","#ff8b94")
s15_R1 <- read.table("15somite_NCC_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s15_R2 <- read.table("15somite_NCC_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s24_R1 <- read.table("24hpf_NCC_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s24_R2 <- read.table("24hpf_NCC_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Mel_R1 <- read.table("Mel_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Mel_R2 <- read.table("Mel_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Iri_R1 <- read.table("Iri_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Iri_R2 <- read.table("Iri_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)

s15_R1$M <- round(s15_R1$X/s15_R1$N, digits = 2)
s15_R2$M <- round(s15_R2$X/s15_R2$N, digits = 2)
s24_R1$M <- round(s24_R1$X/s24_R1$N, digits = 2)
s24_R2$M <- round(s24_R2$X/s24_R2$N, digits = 2)
Mel_R1$M <- round(Mel_R1$X/Mel_R1$N, digits = 2)
Mel_R2$M <- round(Mel_R2$X/Mel_R2$N, digits = 2)
Iri_R1$M <- round(Iri_R1$X/Iri_R1$N, digits = 2)
Iri_R2$M <- round(Iri_R2$X/Iri_R2$N, digits = 2)

# Averge methylation
filter = 5 # Using 5 as cutoff
me <- data.frame(c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"),c(mean(s15_R1[s15_R1$N >= filter,]$M),mean(s15_R2[s15_R2$N >= filter,]$M),mean(s24_R1[s24_R1$N >= filter,]$M),mean(s24_R2[s24_R2$N >= filter,]$M),mean(Mel_R1[Mel_R1$N >= filter,]$M),mean(Mel_R2[Mel_R2$N >= filter,]$M),mean(Iri_R1[Iri_R1$N >= filter,]$M),mean(Iri_R2[Iri_R2$N >= filter,]$M)))
colnames(me) <- c("Type", "Methylation")
me$Type <- factor(me$Type, levels = c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"))
p <- ggplot(data=me, aes(x=Type, y=Methylation, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Average Methylation across samples (coverage filter >=5)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "Average Methylation")+scale_y_continuous(lim = c(0,1))+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=me,aes(x=Type,y=Methylation+0.03,label=paste(round(Methylation,3)*100,"%",sep ="")),size=8)+ theme(legend.position="none")
p
mypalette <- brewer.pal(12,"Paired")

#Number of CpG covered #
co <- data.frame(c(rep("Mel-R1", 3),rep("Mel-R2", 3),rep("Iri-R1", 3),rep("Iri-R2", 3),rep("24hpf-R1",3),rep("24hpf-R2", 3),rep("s15-R1", 3),rep("s15-R2", 3)),c(nrow(Mel_R1[Mel_R1$N>0,]), nrow(Mel_R1[Mel_R1$N>=5,]), nrow(Mel_R1[Mel_R1$N>=10,]),nrow(Mel_R2[Mel_R2$N>0,]), nrow(Mel_R2[Mel_R2$N>=5,]), nrow(Mel_R2[Mel_R2$N>=10,]),nrow(Iri_R1[Iri_R1$N>0,]), nrow(Iri_R1[Iri_R1$N>=5,]), nrow(Iri_R1[Iri_R1$N>=10,]),nrow(Iri_R2[Iri_R2$N>0,]), nrow(Iri_R2[Iri_R2$N>=5,]), nrow(Iri_R2[Iri_R2$N>=10,]),nrow(s24_R1[s24_R1$N>0,]), nrow(s24_R1[s24_R1$N>=5,]), nrow(s24_R1[s24_R1$N>=10,]),nrow(s24_R2[s24_R2$N>0,]), nrow(s24_R2[s24_R2$N>=5,]), nrow(s24_R2[s24_R2$N>=10,]),nrow(s15_R1[s15_R1$N>0,]), nrow(s15_R1[s15_R1$N>=5,]), nrow(s15_R1[s15_R1$N>=10,]),nrow(s15_R2[s15_R2$N>0,]), nrow(s15_R2[s15_R2$N>=5,]), nrow(s15_R2[s15_R2$N>=10,])),c(rep(c(">0",">=5",">=10"),8)))
colnames(co) <- c("Type", "count","Filter")
co$Type <- factor(co$Type, levels = c("s15-R1","s15-R2","24hpf-R1","24hpf-R2","Mel-R1","Mel-R2","Iri-R1","Iri-R2"))
co$Filter <- factor(co$Filter, levels = c(">0",">=5",">=10"))
p <- ggplot(data=co, aes(x=Type, y=count, fill=Filter)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle(paste("CpG coverage cutoff",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "CpG count")+scale_fill_manual(values = mypalette2)+scale_color_manual(values=mypalette2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p

#Methylation distribution#
filter = 5
meth <- data.frame(c(rep("Mel-R1", nrow(Mel_R1[Mel_R1$N >= filter,])),rep("Mel-R2", nrow(Mel_R2[Mel_R2$N >= filter,])),rep("Iri-R1", nrow(Iri_R1[Iri_R1$N >= filter,])),rep("Iri-R2", nrow(Iri_R2[Iri_R2$N >= filter,])),rep("24hpf-R1", nrow(s24_R1[s24_R1$N >= filter,])),rep("24hpf-R2", nrow(s24_R2[s24_R2$N >= filter,])),rep("s15-R1", nrow(s15_R1[s15_R1$N >= filter,])),rep("s15-R2", nrow(s15_R2[s15_R2$N >= filter,]))), c(Mel_R1[Mel_R1$N >= filter,]$M,Mel_R2[Mel_R2$N >= filter,]$M,Iri_R1[Iri_R1$N >= filter,]$M,Iri_R2[Iri_R2$N >= filter,]$M,s24_R1[s24_R1$N >= filter,]$M,s24_R2[s24_R2$N >= filter,]$M,s15_R1[s15_R1$N >= filter,]$M,s15_R2[s15_R2$N >= filter,]$M))
colnames(meth) <- c("Type", "Meth")
a <-ggplot(meth, aes(Meth, fill = Type, colour = Type)) + geom_density(alpha = 0.05, adjust = 3)+ggtitle(paste("Pigment CpG methylation distribution (coverage filter >=",filter,")",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Coverage", y = "Density")+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a

#PCA plot#
filter = 5
m<-Reduce(function(x, y) merge(x, y, by=c("chr","pos")), list(Mel_R1[Mel_R1$N >= filter,c(1,2,5)],Mel_R2[Mel_R2$N >= filter,c(1,2,5)],Iri_R1[Iri_R1$N >= filter,c(1,2,5)],Iri_R2[Iri_R2$N >= filter,c(1,2,5)],s24_R1[s24_R1$N >= filter,c(1,2,5)],s24_R2[s24_R2$N >= filter,c(1,2,5)],s15_R1[s15_R1$N >= filter,c(1,2,5)],s15_R2[s15_R2$N >= filter,c(1,2,5)]))
names(m) <- c("chr","start","Mel-R1","Mel-R2","Iri-R1","Iri-R2","24hpf-R1","24hpf-R2","s15-R1","s15-R2")
m.pca <-prcomp(m[,3:10])
summary(m.pca)

#Importance of components:
#                          PC1     PC2    PC3     PC4     PC5     PC6     PC7     PC8
#Standard deviation     0.7544 0.16948 0.1023 0.07458 0.07282 0.06331 0.06319 0.05850
#Proportion of Variance 0.9025 0.04555 0.0166 0.00882 0.00841 0.00636 0.00633 0.00543
#Cumulative Proportion  0.9025 0.94805 0.9647 0.97347 0.98188 0.98824 0.99457 1.00000

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

loadings <- data.frame(m.pca$rotation, .names = row.names(m.pca$rotation))
a<-ggplot(loadings, aes(x = PC1, y = PC2))+geom_text(data=loadings, 
             mapping=aes(x = PC1, y = PC2, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC1 (Variance = 90.3%)", y = "PC2 (Variance = 4.6%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC1 vs PC2")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
b<-ggplot(loadings, aes(x = PC1, y = PC3))+geom_text(data=loadings, 
              mapping=aes(x = PC1, y = PC3, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC1 (Variance = 90.3%)", y = "PC3 (Variance = 1.7%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC1 vs PC3")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
c<-ggplot(loadings, aes(x = PC2, y = PC3))+geom_text(data=loadings, 
              mapping=aes(x = PC2, y = PC3, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC2 (Variance = 4.6%)", y = "PC3 (Variance = 1.7%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC2 vs PC3")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
multiplot(a,b,c, cols = 3)

library(corrplot)
gsc_cor <- cor(m[,3:10])
corrplot(gsc_cor, order="hclust", addrect=4, method ="number")
corrplot.mixed(gsc_cor, order="hclust", addrect=4)
corrplot(gsc_cor, order="hclust", addrect=4)

#Combine 15somite NCC Rep1 with Rep2#
s15pos_combined <- merge(s15_R1,s15_R2,by=c("chr","pos"), all = TRUE)
s15pos_combined[is.na(s15pos_combined)]<-0
s15pos_combined$N <- s15pos_combined$N.x+s15pos_combined$N.y
s15pos_combined$X <- s15pos_combined$X.x+s15pos_combined$X.y
s15pos_combined <-s15pos_combined[,c(1,2,9,10)]
write.table(s15pos_combined, "15somite_NCC_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

#Combine 24hpf NCC Rep1 with Rep2#
s24pos_combined <- merge(s24_R1,s24_R2,by=c("chr","pos"), all = TRUE)
s24pos_combined[is.na(s24pos_combined)]<-0
s24pos_combined$N <- s24pos_combined$N.x+s24pos_combined$N.y
s24pos_combined$X <- s24pos_combined$X.x+s24pos_combined$X.y
s24pos_combined <-s24pos_combined[,c(1,2,9,10)]
write.table(s24pos_combined, "24hpf_NCC_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

#Combine Mel Rep1 with Rep2#
Mel_combined <- merge(Mel_R1,Mel_R2,by=c("chr","pos"), all = TRUE)
Mel_combined[is.na(Mel_combined)]<-0
Mel_combined$N <- Mel_combined$N.x+Mel_combined$N.y
Mel_combined$X <- Mel_combined$X.x+Mel_combined$X.y
Mel_combined <-Mel_combined[,c(1,2,9,10)]
write.table(Mel_combined, "Mel_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

#Combine Iri Rep1 with Rep2#
Iri_combined <- merge(Iri_R1,Iri_R2,by=c("chr","pos"), all = TRUE)
Iri_combined[is.na(Iri_combined)]<-0
Iri_combined$N <- Iri_combined$N.x+Iri_combined$N.y
Iri_combined$X <- Iri_combined$X.x+Iri_combined$X.y
Iri_combined <-Iri_combined[,c(1,2,9,10)]
write.table(Iri_combined, "Iri_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

s15pos_combined$M <- round(s15pos_combined$X/s15pos_combined$N, digits = 2)
s24pos_combined$M <- round(s24pos_combined$X/s24pos_combined$N, digits = 2)
Mel_combined$M <- round(Mel_combined$X/Mel_combined$N, digits = 2)
Iri_combined$M <- round(Iri_combined$X/Iri_combined$N, digits = 2)


###Look at combined global profiles##
#average methylation

filter = 5 # coverage filter
me <- data.frame(c("15somite","24hpf","Mel","Iri"),c(mean(s15pos_combined[s15pos_combined$N >= filter,]$M),mean(s24pos_combined[s24pos_combined$N >= filter,]$M),mean(Mel_combined[Mel_combined$N >= filter,]$M),mean(Iri_combined[Iri_combined$N >= filter,]$M)))
colnames(me) <- c("Type", "Methylation")
me$Type <- factor(me$Type, levels = c("15somite","24hpf","Mel","Iri"))
# Plot average methylation
p1 <- ggplot(data=me, aes(x=Type, y=Methylation, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Average Methylation across samples (coverage filter >=5)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "Average Methylation")+scale_y_continuous(lim = c(0,1))+scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=me,aes(x=Type,y=Methylation+0.03,label=paste(round(Methylation,3)*100,"%",sep ="")),size=8)+ theme(legend.position="none")
p1

#Number of CpG covered #
co <- data.frame(c(rep("Mel", 3),rep("Iri", 3),rep("24hpf", 3),rep("s15", 3)),c(nrow(Mel_combined[Mel_combined$N>0,]), nrow(Mel_combined[Mel_combined$N>=5,]), nrow(Mel_combined[Mel_combined$N>=10,]),nrow(Iri_combined[Iri_combined$N>0,]), nrow(Iri_combined[Iri_combined$N>=5,]), nrow(Iri_combined[Iri_combined$N>=10,]),nrow(s24pos_combined[s24pos_combined$N>0,]), nrow(s24pos_combined[s24pos_combined$N>=5,]), nrow(s24pos_combined[s24pos_combined$N>=10,]),nrow(s15pos_combined[s15pos_combined$N>0,]), nrow(s15pos_combined[s15pos_combined$N>=5,]), nrow(s15pos_combined[s15pos_combined$N>=10,])),c(rep(c(">0",">=5",">=10"),4)))
colnames(co) <- c("Type", "count","Filter")
co$Type <- factor(co$Type, levels = c("s15","24hpf","Mel","Iri"))
co$Filter <- factor(co$Filter, levels = c(">0",">=5",">=10"))
# Plot CpG coverage cutoff
p2 <- ggplot(data=co, aes(x=Type, y=count, fill=Filter)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle(paste("CpG coverage cutoff",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "CpG count")+scale_fill_manual(values = mypalette2)+scale_color_manual(values=mypalette2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p2

#Methylation distribution
filter = 5
meth <- data.frame(c(rep("Mel", nrow(Mel_combined[Mel_combined$N >= filter,])),rep("Iri", nrow(Iri_combined[Iri_combined$N >= filter,])),rep("24hpf", nrow(s24pos_combined[s24pos_combined$N >= filter,])),rep("s15", nrow(s15pos_combined[s15pos_combined$N >= filter,]))), c(Mel_combined[Mel_combined$N >= filter,]$M,Iri_combined[Iri_combined$N >= filter,]$M,s24pos_combined[s24pos_combined$N >= filter,]$M,s15pos_combined[s15pos_combined$N >= filter,]$M))
colnames(meth) <- c("Type", "Meth")
meth$Type <- factor(meth$Type, levels = c("s15","24hpf","Mel","Iri"))
# Plot CpG methylation distribution
p3 <-ggplot(meth, aes(Meth, fill = Type, colour = Type)) + geom_density(alpha = 0.1, adjust = 4)+ggtitle(paste("Pigment CpG methylation distribution (coverage filter >=",filter,")",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Coverage", y = "Density")+scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p3

#PCA#
filter = 5
m<-Reduce(function(x, y) merge(x, y, by=c("chr","pos")), list(Mel_combined[Mel_combined$N >= filter,c(1,2,5)],Iri_combined[Iri_combined$N >= filter,c(1,2,5)],s24pos_combined[s24pos_combined$N >= filter,c(1,2,5)],s15pos_combined[s15pos_combined$N >= filter,c(1,2,5)]))
names(m) <- c("chr","start","Mel","Iri","24hpf","s15")
m.pca <-prcomp(m[,3:6])
summary(m.pca)

Importance of components:
                          PC1     PC2     PC3     PC4
Standard deviation     0.5325 0.12170 0.07764 0.06329
Proportion of Variance 0.9194 0.04802 0.01954 0.01299
Cumulative Proportion  0.9194 0.96747 0.98701 1.00000

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# PCA plot
loadings <- data.frame(m.pca$rotation, .names = row.names(m.pca$rotation))
a<-ggplot(loadings, aes(x = PC1, y = PC2))+geom_text(data=loadings, 
             mapping=aes(x = PC1, y = PC2, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC1 (Variance = 91.9%)", y = "PC2 (Variance = 4.8%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC1 vs PC2")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
b<-ggplot(loadings, aes(x = PC1, y = PC3))+geom_text(data=loadings, 
              mapping=aes(x = PC1, y = PC3, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC1 (Variance = 91.9%)", y = "PC3 (Variance = 2%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC1 vs PC3")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
c<-ggplot(loadings, aes(x = PC2, y = PC3))+geom_text(data=loadings, 
              mapping=aes(x = PC2, y = PC3, label = .names, colour = .names), fontface = "bold", size = 8) +
    labs(x = "PC2 (Variance = 4.8%)", y = "PC3 (Variance = 2%)")+ggtitle("CpG Methylation PCA plot of replicates (Cov >=5): PC2 vs PC3")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
multiplot(a,b,c, cols = 3)

library(corrplot)
gsc_cor <- cor(m[,3:6])
corrplot(gsc_cor, order="hclust", addrect=4, method ="number")
corrplot.mixed(gsc_cor, order="hclust", addrect=4)
corrplot(gsc_cor, order="hclust", addrect=4)


##CALL DMLs USING DSS##
BSobj <-makeBSseqData(list(s15pos_combined,s24pos_combined,Mel_combined,Iri_combined), c("s15","24hpf","Mel","Iri"))

# get DMRs between each two samples
dml15v24 <- DMLtest(BSobj, group1 = c("s15"), group2 = c("24hpf"), smoothing = TRUE)
dml24vMel <- DMLtest(BSobj, group1 = c("24hpf"), group2 = c("Mel"), smoothing = TRUE)
dml24vIri <- DMLtest(BSobj, group1 = c("24hpf"), group2 = c("Iri"), smoothing = TRUE)
dmlMelvIri <- DMLtest(BSobj, group1 = c("Mel"), group2 = c("Iri"), smoothing = TRUE)
dml15vMel <- DMLtest(BSobj, group1 = c("s15"), group2 = c("Mel"), smoothing = TRUE)
dml15vIri <- DMLtest(BSobj, group1 = c("s15"), group2 = c("Iri"), smoothing = TRUE)

# get DMRs using delta 0.3 and p value 0.01
DMR15v24_d30_p0.01<- callDMR(dml15v24, delta = 0.30, p.threshold= 0.01)
DMR24vMel_d30_p0.01<- callDMR(dml24vMel, delta = 0.30, p.threshold= 0.01)
DMR24vIri_d30_p0.01<- callDMR(dml24vIri, delta = 0.30, p.threshold= 0.01)
DMRMelvIri_d30_p0.01<- callDMR(dmlMelvIri, delta = 0.30, p.threshold= 0.01)
DMR15vMel_d30_p0.01<- callDMR(dml15vMel, delta = 0.30, p.threshold= 0.01)
DMR15vIri_d30_p0.01<- callDMR(dml15vIri, delta = 0.30, p.threshold= 0.01)

# save called DMRs to files
write.table(DMR15v24_d30_p0.01, "DMR15v24_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(DMR24vMel_d30_p0.01, "DMR24vMel_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(DMR24vIri_d30_p0.01, "DMR24vIri_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(DMRMelvIri_d30_p0.01, "DMRMelvIri_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(DMR15vMel_d30_p0.01, "DMR15vMel_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(DMR15vIri_d30_p0.01, "DMR15vIri_d30_p0.01_wSMOOTHING.txt", row.names = F, col.names = F, sep = "\t",quote =F)

# Merge all DMRs into master list
cat DMR15v24_d30_p0.01_wSMOOTHING.txt DMR24vIri_d30_p0.01_wSMOOTHING.txt DMR24vMel_d30_p0.01_wSMOOTHING.txt DMRMelvIri_d30_p0.01_wSMOOTHING.txt> Combined_DMRs_d30_p0.01.txt
sort -k1,1 -k2,2n Combined_DMRs_d30_p0.01.txt > Combined_DMRs_d30_p0.01.sorted.txt #delete headers manually
awk -F',' '{gsub(/"/, "", $1); print $1}' Combined_DMRs_d30_p0.01.sorted.txt >  Combined_DMRs_d30_p0.01.sorted.fixed.txt #remove quotes from column 1

# merge combined DMRs
bedtools merge -i Combined_DMRs_d30_p0.01.sorted.fixed.txt > Combined_DMRs_d30_p0.01.txt 
rm Combined_DMRs_d30_p0.01.sorted.txt
wc -l  Combined_DMRs_d30_p0.01.txt 
rm Combined_DMRs_d30_p0.01.sorted.fixed.txt
rm Combined_DMRs_d30_p0.01.sorted.txt

#pull out DMRs associated with combined_DMR
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR15v24_d30_p0.01_wSMOOTHING.txt > Combined_DMR15v24_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR24vMel_d30_p0.01_wSMOOTHING.txt > Combined_DMR24vMel_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR24vIri_d30_p0.01_wSMOOTHING.txt > Combined_DMR24vIri_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMRMelvIri_d30_p0.01_wSMOOTHING.txt > Combined_DMRMelvIri_d30_p0.01_wSMOOTHING.txt


#Graph Methylation differences and also correlation in shared DMRs
s24 <-combined_dmr_wInfo[combined_dmr_wInfo$s15vs24 !=0,] 
M<-combined_dmr_wInfo[combined_dmr_wInfo$s24vMel !=0,]
I<- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri !=0,]
MI <-combined_dmr_wInfo[combined_dmr_wInfo$MelvsIri !=0,]

ps24 <-ggplot(s24, aes(s15vs24)) + geom_density(alpha = 0.5,adjust =0.15,color = mypalette[3],fill = mypalette[3])+ggtitle(paste("DMR methylation change from 15somite to 24hpf"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ps24
pM <-ggplot(M, aes(s24vMel)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[5],fill = mypalette[5])+ggtitle(paste("DMR methylation change from 24hpf to Melanophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM
pI <-ggplot(I, aes(s24vIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[7],fill = mypalette[7])+ggtitle(paste("DMR methylation change from 24hpf to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI
pMI <-ggplot(shar, aes(MelvsIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[9],fill = mypalette[9])+ggtitle(paste("DMR methylation change from Melanophore to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI

multiplot(ps24,pI,pM,pMI, cols = 2)


Mel_Iri_shared_hyperDMR$Type <- "Shared HyperDMR"
Mel_Iri_shared_hypoDMR$Type <- "Shared HypoDMR"
Mel_specific_hyperDMR$Type <- "Mel-specific HyperDMR"
Mel_specific_hypoDMR$Type <- "Mel-specific HypoDMR"
Iri_specific_hyperDMR$Type <- "Iri-specific HyperDMR"
Iri_specific_hypoDMR$Type <- "Iri-specific HypoDMR"

mypalette3 <- brewer.pal(8,"Paired")

shared <- rbind(Mel_Iri_shared_hyperDMR,Mel_Iri_shared_hypoDMR,Mel_specific_hyperDMR,Mel_specific_hypoDMR[Mel_specific_hypoDMR$s24vMel > 0.2 | Mel_specific_hypoDMR$s24vMel < -0.2,],Iri_specific_hyperDMR,Iri_specific_hypoDMR)
ggplot(shared[shared$type == "Shared hyperDMR" | shared$type == "Shared hypoDMR",], aes(x=s24vMel, y=s24vIri,fill=Type,color=Type)) +
  geom_point(size = 0.2) + geom_jitter(data = shared[shared$type != "Shared hyperDMR" | shared$type != "Shared hypoDMR",],aes(x=s24vMel, y=s24vIri,fill=Type,color=Type),width = 0.01, height = 0.01,size = 0.2)+
  ggtitle(paste("DMR methylation change comparison from 24hpf to Pigment cells"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 20), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+labs(x = "Methylaton change from 24hpf to Melanophore", y = "Methylaton change from 24hpf to Iridophore")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title = element_text(size = 12, face = "bold"),legend.text = element_text(size = 10, face = "bold"))+ theme(legend.position="top")+ guides(color=guide_legend(override.aes=list(fill=NA,size = 6)))+ theme(legend.key=element_blank())


#Merge all samples d30
combined_dmr <- read.table("Combined_DMRs_d30_p0.01.txt",sep = "\t", header = F, stringsAsFactors = F)

#Plot DMR size distribution
combined_dmr$size <- combined_dmr$V3-combined_dmr$V2
hist(combined_dmr$size,breaks = 1000, xlim = c(0,2000),xlab = "DMR Size", main = "Combined DMR size distribution (d30)", col = "grey")

#Read DMR files
cDMR15v24_d30_p0.01 <- read.table("Combined_DMR15v24_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vMel_d30_p0.01 <- read.table("Combined_DMR24vMel_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vIri_d30_p0.01 <- read.table("Combined_DMR24vIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMRMelvIri_d30_p0.01 <- read.table("Combined_DMRMelvIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)

# Add a chromosome position column, like chr1:12345-67890
cDMR15v24_d30_p0.01$chrompos<- paste(cDMR15v24_d30_p0.01$V1,":",cDMR15v24_d30_p0.01$V2,"-",cDMR15v24_d30_p0.01$V3,sep = "")
cDMR24vMel_d30_p0.01$chrompos<- paste(cDMR24vMel_d30_p0.01$V1,":",cDMR24vMel_d30_p0.01$V2,"-",cDMR24vMel_d30_p0.01$V3,sep = "")
cDMR24vIri_d30_p0.01$chrompos<- paste(cDMR24vIri_d30_p0.01$V1,":",cDMR24vIri_d30_p0.01$V2,"-",cDMR24vIri_d30_p0.01$V3,sep = "")
cDMRMelvIri_d30_p0.01$chrompos<- paste(cDMRMelvIri_d30_p0.01$V1,":",cDMRMelvIri_d30_p0.01$V2,"-",cDMRMelvIri_d30_p0.01$V3,sep = "")

# Merge DMR files
combined_dmr_wInfo <-Reduce(function(x, y) merge(x, y, by=c("chrompos"), all = T), list(cDMR15v24_d30_p0.01[,c("V1","V2","V3","chrompos","V11")],cDMR24vMel_d30_p0.01[,c("chrompos","V11")],cDMR24vIri_d30_p0.01[,c("chrompos","V11")],cDMRMelvIri_d30_p0.01[,c("chrompos","V11")]))
colnames(combined_dmr_wInfo) <- c("chrompos","chr","start","end","s15vs24","s24vMel","s24vIri","MelvsIri")
combined_dmr_wInfo$size <- combined_dmr_wInfo$end-combined_dmr_wInfo$start
combined_dmr_wInfo[combined_dmr_wInfo == "."] <- 0
combined_dmr_wInfo$s15vs24 <- as.numeric(combined_dmr_wInfo$s15vs24)
combined_dmr_wInfo$s24vMel <- as.numeric(combined_dmr_wInfo$s24vMel)
combined_dmr_wInfo$s24vIri <- as.numeric(combined_dmr_wInfo$s24vIri)
combined_dmr_wInfo$MelvsIri <- as.numeric(combined_dmr_wInfo$MelvsIri)

combined_dmr_wInfo <- combined_dmr_wInfo[,c(2,3,4,1,9,5,6,7,8)]

combined_dmr_wInfo <- combined_dmr_wInfo[!duplicated(combined_dmr_wInfo$chrompos),]
write.table(combined_dmr_wInfo, "Combined_DMRs_d30_p0.01_wINFO.bed", row.names = F, col.names = F, sep = "\t",quote =F)

# pull out cell type-specific
Mel_specific_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <0 & combined_dmr_wInfo$s24vIri >=0,] 
Mel_specific_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri <=0,] 
# using more strigent way to call Mel hypoDMR
Mel_specific_hypoDMR_stringent <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri <=0 & combined_dmr_wInfo$MelvsIri <0,] 

Iri_specific_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >=0 & combined_dmr_wInfo$s24vIri <0,] 
Iri_specific_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <=0 & combined_dmr_wInfo$s24vIri >0,] 
# using more strigent way to call Iri hypoDMR
Iri_specific_hypoDMR_stringent <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <=0 & combined_dmr_wInfo$s24vIri >0 & combined_dmr_wInfo$MelvsIri >0,] 

Mel_Iri_shared_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <0 & combined_dmr_wInfo$s24vIri <0,] 
Mel_Iri_shared_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri >0,] 

# write hypo DMRs to files
write.table(Mel_specific_hyperDMR, "Mel_specific_hyperDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_hypoDMR, "Mel_specific_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_hypoDMR_stringent, "Mel_specific_hypoDMR_d30_p0.01_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Iri_specific_hyperDMR, "Iri_specific_hyperDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_hypoDMR, "Iri_specific_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_hypoDMR_stringent, "Iri_specific_hypoDMR_d30_p0.01_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Mel_Iri_shared_hyperDMR , "Mel_Iri_shared_hyperDMR _d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_Iri_shared_hypoDMR, "Mel_Iri_shared_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)

#Graph DMR size distribution
NCC_DMR <- combined_dmr_wInfo[combined_dmr_wInfo$s15vs24 != 0,]
M_DMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel != 0,]
I_DMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri != 0,]
MI_DMR <- combined_dmr_wInfo[combined_dmr_wInfo$MelvsIri != 0,]

ps24 <-ggplot(NCC_DMR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[4],fill = mypalette[4])+ggtitle(paste("DMR size distribution\n(15somite vs 24hpf)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM <-ggplot(M_DMR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[6],fill = mypalette[6])+ggtitle(paste("DMR size distribution\n(24hpf vs Melanophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI <-ggplot(I_DMR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[8],fill = mypalette[8])+ggtitle(paste("DMR size distribution\n(24hpf vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI <-ggplot(MI_DMR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[10],fill = mypalette[10])+ggtitle(paste("DMR size distribution\n(Melanophore vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
multiplot(ps24,pI,pM,pMI, cols = 2)


#Graph Methylation differences and also correlation in shared DMRs
s24 <-combined_dmr_wInfo[combined_dmr_wInfo$s15vs24 !=0,] 
M<-combined_dmr_wInfo[combined_dmr_wInfo$s24vMel !=0,]
I<- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri !=0,]
MI <-combined_dmr_wInfo[combined_dmr_wInfo$MelvsIri !=0,]

ps24 <-ggplot(s24, aes(s15vs24)) + geom_density(alpha = 0.5,adjust =0.15,color = mypalette[3],fill = mypalette[3])+ggtitle(paste("DMR methylation change from 15somite to 24hpf"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ps24
pM <-ggplot(M, aes(s24vMel)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[5],fill = mypalette[5])+ggtitle(paste("DMR methylation change from 24hpf to Melanophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM
pI <-ggplot(I, aes(s24vIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[7],fill = mypalette[7])+ggtitle(paste("DMR methylation change from 24hpf to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI
pMI <-ggplot(shar, aes(MelvsIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[9],fill = mypalette[9])+ggtitle(paste("DMR methylation change from Melanophore to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI

multiplot(ps24,pI,pM,pMI, cols = 2)
```
