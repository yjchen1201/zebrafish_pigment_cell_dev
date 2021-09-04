### Tool version
- cutadapt v1.10: Used to trim adapter reads 
- samtools version 1.3.1: Used to sort and downsample bam files for downstream processing
- bismark v0.16.1: Used to align WGBS reads to the genome and call CpG methylation
- bwa v0.7.15: bwa-mem used to align ATAC-seq reads to the genome
- macs2 v2.1.1: Used to call peaks for ATAC-seq
- STAR v2.5.1b: Used to align RNA-sequencing data
- stringtie v1.3.3: Used to annotate aligned reads and quantify transcripts. 
- picard v2.8.1: Used to mark and remove duplicate reads
- bedtools v2.27.1: closest, intersect, shuffle command used as specified in manuscript 
- meme v5.0.3: Used AME command to find motif enrichment and FIMO command to scan for motif presence 
- R v3.3.0: Used to analyze sequencing library results
- DSS v2.14.0 :Used to call differentially methylated regions
- DESeq2 v1.12.4: Used to call differentially expressed genes
- DiffBind v 2.2.12: Used to call differentially accessible regions/peaks
- Metascape v3.0: Used for GO enrichment analysis

### Analysis example
```{bash}
## Example ##
# Adapter trimming
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1.fastq.gz" "SAMPLE_R2.fastq.gz"
# Alignment using STAR, output is BAM file sorted by coordinate
STAR --readFilesCommand zcat --outTmpDir <out_temp_dirctory> --clip5pNbases 11 --outSAMtype BAM SortedByCoordinate --genomeDir /bar/genomes/danRer10/STAR --sjdbGTFfile "Danio_rerio.GRCz10.85.gtf" --outFileNamePrefix "SAMPLE" --readFilesIn "SAMPLE_R1_trimmed.fastq.gz" "SAMPLE_R2_trimmed.fastq.gz"

# Normalziation to RPM
STAR --runMode inputAlignmentsFromBAM --inputBAMfile "SAMPLE_sortedByCoord.bam" --outWigType bedGraph --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix "SAMPLE_sortedByCoord.out.bam"
htseq-count -f bam -r pos -s no -t transcript $i "Danio_rerio.GRCz10.85.gtf" > "SAMPLE_sortedByCoord.counts"

stringtie "SAMPLE_sortedByCoord.out.bam" -o "RefSeq.gtf" -A "RefSeq.gene.abundance.txt" -G "danRer10_RefSeq_forstringtie.gtf"

#Filter NM genes only#
python3 Filter_only_NM.py "SAMPLE_abundance.txt"
#Make list of all genes expressed in NCC and pigment cells (RefSeq)
python3 All_expressed_genelist_RefSeq.py RNA_15somiteNCC_Rep2.gene.abundance.filteredNM.txt RNA_15somiteNCC_Rep4.gene.abundance.filteredNM.txt RNA_24hpfNCC_Rep1.gene.abundance.filteredNM.txt RNA_24hpfNCC_Rep2.gene.abundance.filteredNM.txt RNA_Mel_Rep1.gene.abundance.filteredNM.txt RNA_Mel_Rep2.gene.abundance.filteredNM.txt RNA_Iri_Rep2.gene.abundance.filteredNM.txt RNA_Iri_Rep3.gene.abundance.filteredNM.txt 
```

```{R}
############################################# {R script} #############################################
#Merge all gene expression
library(ggplot2)
library(RColorBrewer)
mypalette <- brewer.pal(12,"Paired")

genelist <- read.table("<path to RNA folder>/stringtie/Combined_gene_list_RefSeq.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

s15_R1 <- read.table("<path to RNA folder>/stringtie/RNA_15somiteNCC_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s15_R2 <- read.table("<path to RNA folder>/stringtie/RNA_15somiteNCC_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_R1 <- read.table("<path to RNA folder>/stringtie/RNA_24hpfNCC_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_R2 <- read.table("<path to RNA folder>/stringtie/RNA_24hpfNCC_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Mel_R1 <- read.table("<path to RNA folder>/stringtie/RNA_Mel_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Mel_R2 <- read.table("<path to RNA folder>/stringtie/RNA_Mel_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Iri_R1 <- read.table("<path to RNA folder>/stringtie/RNA_Iri_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Iri_R2 <- read.table("<path to RNA folder>/stringtie/RNA_Iri_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

geneTPM<-Reduce(function(x, y) merge(x, y,by = "V1",all.x = T), list(s15_R1[,c(1,9)],s15_R2[,c(1,9)], s24_R1[,c(1,9)],s24_R2[,c(1,9)],Mel_R1[,c(1,9)],Mel_R2[,c(1,9)], Iri_R1[,c(1,9)],Iri_R2[,c(1,9)]))
colnames(geneTPM) <- c("gene","s15_Rep1","s15_Rep2","s24_Rep1","s24_Rep2","Mel_Rep1","Mel_Rep2","Iri_Rep1","Iri_Rep2")
geneTPM <- geneTPM[!duplicated(geneTPM$gene),]
geneTPM$s15 <- (geneTPM$s15_Rep1+geneTPM$s15_Rep2)/2
geneTPM$s24 <- (geneTPM$s24_Rep1+geneTPM$s24_Rep2)/2
geneTPM$Mel <- (geneTPM$Mel_Rep1+geneTPM$Mel_Rep2)/2
geneTPM$Iri <- (geneTPM$Iri_Rep1+geneTPM$Iri_Rep2)/2


#Use DESeq2 to find differentially expressed genes in each comparison
library(DESeq2)
setwd("<path to RNA folder>")
directory <- "<path to RNA folder>/htseq"
sampleFiles <-grep("out.counts",list.files(directory),value =TRUE)
sampleCondition <- c("15somite","15somite","24hpf","24hpf","Iri","Iri","Mel","Mel")
sampleNames <- c("s15pos_R1","s15pos_R2","s24pos_R1","s24pos_R2","Iri_R1","Iri_R2","Mel_R1","Mel_R2")
sampleTable <- data.frame(sampleName = sampleNames,fileName = sampleFiles,condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory, design= ~ condition)
dds <- ddsHTSeq[rowSums(counts(ddsHTSeq))>1,]
dds1 <- DESeq(dds)

res_15pos_24hpf <- results(dds1, contrast=c("condition","15somite","24hpf"))
res_24hpf_Mel<- results(dds1, contrast=c("condition","24hpf","Mel"))
res_24hpf_Iri<- results(dds1, contrast=c("condition","24hpf","Iri"))
res_Mel_Iri<- results(dds1, contrast=c("condition","Mel","Iri"))

res_15pos_24hpf_p0.01 <- subset(res_15pos_24hpf, padj < 0.01)
res_24hpf_Mel_p0.01 <- subset(res_24hpf_Mel, padj < 0.01)
res_24hpf_Iri_p0.01 <- subset(res_24hpf_Iri, padj < 0.01)
res_Mel_Iri_p0.01 <- subset(res_Mel_Iri, padj < 0.01)

#export the DEGs, Columns: gene  baseMean  log2FoldChange  lfcSE stat  pvalue  padj
write.table(as.data.frame(res_15pos_24hpf_p0.01),file="<path to RNA folder>/DESeq2/DEGs_15somite_24hpf_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(res_24hpf_Mel_p0.01),file="<path to RNA folder>/DESeq2/DEGs_24hpf_Mel_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(res_24hpf_Iri_p0.01),file="<path to RNA folder>/DESeq2/DEGs_24hpf_Iri_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(res_Mel_Iri_p0.01),file="<path to RNA folder>y/DESeq2/DEGs_Mel_Iri_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")


#PCA plot
rld <- rlog(dds1, blind=FALSE)
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition)) +
geom_point(size=6) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +labs(title = "PCA analysis of RNA expression")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))

#Hierchal clustering heatmap
library("RColorBrewer")#heatmap based on euclidian distance
library(pheatmap)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
clustering_distance_rows=sampleDists,
clustering_distance_cols=sampleDists,
col=colors)

#MAplots
par(mfrow=c(2,2))
plotMA(res_15pos_24hpf, main="15somite vs 24hpf", ylim=c(-6,6), alpha = 0.01)
plotMA(res_24hpf_Mel, main="24hpf vs Melanophore", ylim=c(-6,6), alpha = 0.01)
plotMA(res_24hpf_Iri, main="24hpf vs Iridophore", ylim=c(-6,6), alpha = 0.01)
plotMA(res_Mel_Iri, main="Melanophore vs Iridophore", ylim=c(-6,6), alpha = 0.01)

par(mfrow=c(2,2))
plotMA(res_15pos_24hpf, main="15somite vs 24hpf", ylim=c(-6,6), alpha = 0.001)
plotMA(res_24hpf_Mel, main="24hpf vs Melanophore", ylim=c(-6,6), alpha = 0.001)
plotMA(res_24hpf_Iri, main="24hpf vs Iridophore", ylim=c(-6,6), alpha = 0.001)
plotMA(res_Mel_Iri, main="Melanophore vs Iridophore", ylim=c(-6,6), alpha = 0.001)
```

```{bash}
############################################# {Bash script} #############################################
#MAKE master list of DEGs for comparision
python3 All_DEGs_list.py DEGs_15somite_24hpf_p0.01.txt DEGs_24hpf_Mel_p0.01.txt DEGs_24hpf_Iri_p0.01.txt DEGs_Mel_Iri_p0.01.txt
```

```{R}
############################################# {R script} #############################################
setwd("<path to RNA folder>/DESeq2")
genelist <- read.table("Combined_DEGs_list.txt",,sep = "\t", header = F,quote = "", stringsAsFactors = F)
s15_s24<- read.table("DEGs_15somite_24hpf_p0.01.txt",,sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_M <- read.table("DEGs_24hpf_Mel_p0.01.txt",,sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_I <- read.table("DEGs_24hpf_Iri_p0.01.txt",,sep = "\t", header = F,quote = "", stringsAsFactors = F)
M_I <- read.table("DEGs_Mel_Iri_p0.01.txt",,sep = "\t", header = F,quote = "", stringsAsFactors = F)

DEG <-Reduce(function(x, y) merge(x, y,by = "V1", all = T), list(genelist,s15_s24[,c(1:3)],s24_M[,c(1:3)],s24_I[,c(1:3)],M_I[,c(1:3)]))
colnames(DEG) <- c("gene","s15v24hpf_mean_expression","s15v24hpf_log2_change","s24vMel_mean_expression","s24vMel_log2_change","s24vIri_mean_expression","s24vIri_log2_change","MelvIri_mean_expression","MelvIri_log2_change")
DEG[is.na(DEG)]<-0

write.table(DEG,"DEGs_combined_samples_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")

Mel_specific_DEG_on <- DEG[DEG$s24vMel_log2_change <0 & DEG$s24vIri_log2_change >=0 & DEG$MelvIri_log2_change >0,]
Mel_specific_DEG_off <- DEG[DEG$s24vMel_log2_change >0 & DEG$s24vIri_log2_change <=0 & DEG$MelvIri_log2_change <0,]


setwd("<path to RNA folder>/DESeq2")
genelist2 <- read.table("Combined_DEGs_list_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s15_s242<- read.table("DEGs_15somite_24hpf_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_M2 <- read.table("DEGs_24hpf_Mel_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_I2 <- read.table("DEGs_24hpf_Iri_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
M_I2 <- read.table("DEGs_Mel_Iri_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

DEG2 <-Reduce(function(x, y) merge(x, y,by = "V1", all = T), list(genelist2,s15_s242[,c(1:3)],s24_M2[,c(1:3)],s24_I2[,c(1:3)],M_I2[,c(1:3)]))
colnames(DEG2) <- c("gene","s15v24hpf_mean_expression","s15v24hpf_log2_change","s24vMel_mean_expression","s24vMel_log2_change","s24vIri_mean_expression","s24vIri_log2_change","MelvIri_mean_expression","MelvIri_log2_change")
DEG2[is.na(DEG2)]<-0

DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG2_genenames) <- c("gene","genename")
DEG2<- merge(DEG2,DEG2_genenames,by = "gene",all = T)
DEG2 <- DEG2[!duplicated(DEG2$gene),]

write.table(DEG2,"DEGs_combined_samples_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")


### GOT TF list from AnimalTFDB ####
# Lots of TF missing in zebrafish. Combine with human TFs that have been converted using orthoretriever
TFlist <- read.table("Danio_rerio_transcription_factors_gene_list.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
homo_TFlist <- read.table("homo_danrer_TF.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
cofactorlist <- read.table("Danio_rerio_transcription_co-factors_gene_list.txt",header =F, ,quote = "", sep = "\t",stringsAsFactors = F)
homo_cofactorlist <- read.table("homo_danrer_cofactor.txt",header =F, ,quote = "", sep = "\t",stringsAsFactors = F)
CRMlist <- read.table("Danio_rerio_chromatin_remodeling_factors_gene_list.txt",header = F, ,quote = "", sep = "\t",stringsAsFactors = F)
homo_CRMlist <- read.table("homo_danrer_CHRF.txt",header = F, ,quote = "", sep = "\t",stringsAsFactors = F)


TFs <- c(TFlist[,1],homo_TFlist[,3])
TF <- unique(TFs)

CFs <- c(cofactorlist[,1],homo_cofactorlist[,3])
CF <- unique(CFs)

CRFs <- c(CRMlist[,1],homo_CRMlist[,3])
CRF <- unique(CRFs )

TFlist <- data.frame(TF, c(rep("TF",length(TF))),stringsAsFactors = F)
CFlist <- data.frame(CF, c(rep("Cofactor",length(CF))),stringsAsFactors = F)
CRFlist <- data.frame(CRF, c(rep("ChomatinRemodeler",length(CRF))),stringsAsFactors = F)

colnames(TFlist) <- c( "gene" ,  "type" )
colnames(CFlist) <- c( "gene" , "type" )
colnames(CRFlist) <- c( "gene"    , "type" )

all_list <- rbind(TFlist,CFlist,CRFlist, stringsAsFactors = F)

DEG2 <- read.table("DEGs_combined_samples_p0.01.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
DEG2<- merge(DEG2,all_list,by = "gene",all.x=T)
DEG2[is.na(DEG2)] <- "other"

#sometimes the name conversion is wrong so also make sure common name is included
homo_TFlist$Input.Common.Name <- tolower(homo_TFlist$Input.Common.Name)
homo_cofactorlist$Input.Common.Name <- tolower(homo_cofactorlist$Input.Common.Name)
homo_CRMlist$Input.Common.Name <- tolower(homo_CRMlist$Input.Common.Name)

DEG2[DEG2$genename %in% homo_TFlist$Input.Common.Name,]$type <- "TF"
DEG2[DEG2$genename %in% homo_cofactorlist$Input.Common.Name,]$type <- "Cofactor"
DEG2[DEG2$genename %in% homo_CRMlist$Input.Common.Name,]$type <- "ChomatinRemodeler"

write.table(DEG2,"DEGs_combined_samples_p0.01_wTFinfo.txt",col.names = T, row.names = F,quote = F, sep = "\t")

#### Identify DEGs for each cell type
DEG2 <- read.table("DEGs_combined_samples_p0.01_wTFinfo.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
DEG2<-merge(DEG2,geneTPM,by = "gene", all.x =T)##ALL GENES
DEG2<-DEG2[DEG2$s15 > 5 | DEG2$s24 > 5 |DEG2$Mel > 5 |DEG2$Iri > 5,]

str(DEG2[DEG2$s15v24hpf_log2_change >0,])

# get cell type-specific up(on) down(off) regulated DEGs
Mel_specific_DEG_on <- DEG2[DEG2$s24vMel_log2_change <0 & DEG2$s24vIri_log2_change >=0 & DEG2$MelvIri_log2_change >0,] 
Mel_specific_DEG_off <- DEG2[DEG2$s24vMel_log2_change >0 & DEG2$s24vIri_log2_change <=0 & DEG2$MelvIri_log2_change <0,]
Iri_specific_DEG_on <- DEG2[DEG2$s24vMel_log2_change >=0 & DEG2$s24vIri_log2_change <0 & DEG2$MelvIri_log2_change <0,] #351
Iri_specific_DEG_off <- DEG2[DEG2$s24vMel_log2_change <=0 & DEG2$s24vIri_log2_change >0 & DEG2$MelvIri_log2_change >0,] #241
Shared_Mel_Iri_DEG_on <- DEG2[DEG2$s24vMel_log2_change <0 & DEG2$s24vIri_log2_change <0,] #1105
Shared_Mel_Iri_DEG_off <- DEG2[DEG2$s24vMel_log2_change >0 & DEG2$s24vIri_log2_change >0,] #1327

Mel_specific_DEG_on_loose <- DEG2[DEG2$s24vMel_log2_change <0 & DEG2$s24vIri_log2_change >=0,] #523
Mel_specific_DEG_off_loose <- DEG2[DEG2$s24vMel_log2_change >0 & DEG2$s24vIri_log2_change <=0,] #445
Iri_specific_DEG_on_loose <- DEG2[DEG2$s24vMel_log2_change >=0 & DEG2$s24vIri_log2_change <0,] #730
Iri_specific_DEG_off_loose <- DEG2[DEG2$s24vMel_log2_change <=0 & DEG2$s24vIri_log2_change >0,] #594

s24_specific_DEG_on_on <- DEG2[DEG2$s15v24hpf_log2_change < 0 & DEG2$s24vMel_log2_change <=0 & DEG2$s24vIri_log2_change <=0,] #373
s24_specific_DEG_off_off <- DEG2[DEG2$s15v24hpf_log2_change > 0 & DEG2$s24vMel_log2_change >=0 & DEG2$s24vIri_log2_change >=0,] #269

s24_Iri_only_DEG_on_on <- DEG2[DEG2$s15v24hpf_log2_change < 0 & DEG2$s24vMel_log2_change >0 & DEG2$s24vIri_log2_change <=0,] #43
s24_Mel_only_DEG_on_on <- DEG2[DEG2$s15v24hpf_log2_change < 0 & DEG2$s24vMel_log2_change <=0 & DEG2$s24vIri_log2_change >0,] #74

MelvsIri_DEG_off <- DEG2[DEG2$MelvIri_log2_change < 0,]
MelvsIri_DEG_on <- DEG2[DEG2$MelvIri_log2_change > 0,]

# Save DEGs into a file
write.table(s24_Iri_only_DEG_on_on,"s24_Iri_only_DEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_Mel_only_DEG_on_on,"s24_Mel_only_DEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_on,"Mel_specific_DEG_on_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_off,"Mel_specific_DEG_off_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_on,"Iri_specific_DEG_on_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_off,"Iri_specific_DEG_off_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Shared_Mel_Iri_DEG_on,"Shared_Mel_Iri_DEG_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Shared_Mel_Iri_DEG_off,"Shared_Mel_Iri_DEG_off_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_on_loose,"Mel_specific_DEG_on_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_off_loose,"Mel_specific_DEG_off_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_on_loose,"Iri_specific_DEG_on_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_off_loose,"Iri_specific_DEG_off_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_specific_DEG_on_on,"s24_specific_DEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_specific_DEG_off_off,"s24_specific_DEG_off_off_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")


################# MAKE GENE EXPRESSION PLOTS (check to make sure samples reflect cell type biology)###########
#15somite vs 24hpf
res_15pos_24hpf $log2FoldChange <- res_15pos_24hpf $log2FoldChange*-1

plotMA(res_15pos_24hpf , main="15somite NCC vs 24hpf NCC DEGs", ylim=c(-10,10), alpha = 0.01,colNonSig = "grey", colSig = "#de7b5c")
s15gene<-"ENSDARG00000003293"
s15genename <- "sox9a"
NCCGene <- c("ENSDARG00000100398","ENSDARG00000076010","ENSDARG00000030402")
NCCGeneName <- c("pax7a","twist1b","twist1a") 
with(res_15pos_24hpf[NCCGene, ], {
  points(baseMean, log2FoldChange, col="#a200ff", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, NCCGeneName,cex =2, pos=4, col="#a200ff")
})
with(res_15pos_24hpf[s15gene, ], {
  points(baseMean, log2FoldChange, col="black", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, s15genename,cex =2, pos=4, col="black")
})

DEG_s24_Mel <- as.data.frame(res_24hpf_Mel_p0.01)
DEG_s24_Mel$gene <- rownames(res_24hpf_Mel_p0.01)
DEG_s24_Mel$log2FoldChange <- DEG_s24_Mel$log2FoldChange*-1

DEG_s24_Mel$color <-0
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_on_loose$gene,]$color <- 1
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_on$gene,]$color <- 2
DEG_s24_Mel[DEG_s24_Mel$gene %in% Shared_Mel_Iri_DEG_on$gene,]$color <- 3
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_off_loose$gene,]$color <-4
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_off$gene,]$color <- 5
DEG_s24_Mel[DEG_s24_Mel$gene %in% Shared_Mel_Iri_DEG_off$gene,]$color <- 6

#24hpf vs Mel
res_24hpf_Mel$log2FoldChange <- res_24hpf_Mel$log2FoldChange*-1

plotMA(res_24hpf_Mel, main="24hpf NCC vs Melanophore DEGs", ylim=c(-10,10), alpha = 0.01,colNonSig = "grey", colSig = "#de7b5c")
MelGene <- c("ENSDARG00000039077","ENSDARG00000006008","ENSDARG00000091298","ENSDARG00000043317")
MelGeneName <- c("tyr","dct","pmela","kita") 
NCCGene <- c("ENSDARG00000100398","ENSDARG00000076010","ENSDARG00000030402")
NCCGeneName <- c("pax7a","twist1b","twist1a") 
with(res_24hpf_Mel[MelGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, MelGeneName,cex =2, pos=4, col="dodgerblue")
})
with(res_24hpf_Mel[NCCGene, ], {
  points(baseMean, log2FoldChange, col="#a200ff", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, NCCGeneName,cex =2, pos=4, col="#a200ff")
})

#24hpf vs Iri
res_24hpf_Iri$log2FoldChange <- res_24hpf_Iri$log2FoldChange*-1

plotMA(res_24hpf_Iri, main="24hpf NCC vs Iridophore DEGs", ylim=c(-10,10), alpha = 0.01,colNonSig = "grey", colSig = "#de7b5c")

IriGene <- c("ENSDARG00000057575","ENSDARG00000042861","ENSDARG00000089334")
IriGeneName <- c("pnp4a","ltk","ednrba") 
NCCGene <- c("ENSDARG00000100398","ENSDARG00000076010","ENSDARG00000030402")
NCCGeneName <- c("pax7a","twist1b","twist1a") 
with(res_24hpf_Iri[IriGene, ], {
  points(baseMean, log2FoldChange, col="#1f9400", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, IriGeneName,cex =2, pos=4, col="#1f9400")
})
with(res_24hpf_Iri[NCCGene, ], {
  points(baseMean, log2FoldChange, col="#a200ff", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, NCCGeneName,cex =2, pos=4, col="#a200ff")
})


#Mel vs Iri
res_Mel_Iri$log2FoldChange <- res_Mel_Iri$log2FoldChange*-1

plotMA(res_Mel_Iri, main="Melanophore vs Iridophore DEGs", ylim=c(-10,10), alpha = 0.01,colNonSig = "grey", colSig = "#de7b5c")
IriGene <- c("ENSDARG00000057575","ENSDARG00000042861","ENSDARG00000089334")
IriGeneName <- c("pnp4a","ltk","ednrba") 
MelGene <- c("ENSDARG00000039077","ENSDARG00000006008","ENSDARG00000091298","ENSDARG00000043317")
MelGeneName <- c("tyr","dct","pmela","kita") 
with(res_Mel_Iri[MelGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, MelGeneName,cex =2, pos=4, col="dodgerblue")
})
with(res_Mel_Iri[IriGene, ], {
  points(baseMean, log2FoldChange, col="#1f9400", cex=3, lwd=4)
  text(baseMean+(baseMean/10), log2FoldChange, IriGeneName,cex =2, pos=4, col="#1f9400")
})


# Get transcription factors that differentially expressed

## All DEGs
DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG2_genenames) <- c("gene","genename")
DEG3 <- merge(DEG2,DEG2_genenames, by = "gene", all.x = T)
DEG3<-merge(DEG3,geneTPM, by = "gene", all.x =T) ##ALL DEGs
DEG3[is.na(DEG3)] <-0

## All genes
DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG2_genenames) <- c("gene","genename")
DEG3 <- merge(DEG2,DEG2_genenames, by = "gene", all.x = T)
DEG4<-merge(geneTPM, DEG3,by = "gene", all.x =T)##ALL GENES
DEG4[is.na(DEG4)] <-0


TF_DEG2 <- DEG3[DEG3$type == "TF",] 
TF_DEG <- merge(TF_DEG2,geneTPM,by = "gene",all.x=TRUE)
TF_DEG$Ave_RPKM <- (TF_DEG$s15+TF_DEG$s24+TF_DEG$Mel+TF_DEG$Iri)/4
TF_DEG2<- TF_DEG[,c(1:10,16)]
library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- TF_DEG2[,1]
df.OG <- TF_DEG2[,c(5,7,9)]
row.names(df.OG) <- name

#kclus <- kmeans(df.OG,20)
#split <-  kclus$cluster
#ht = Heatmap(df.OG, column_title = "Mel-specific only DAR (closing)",name= "Methylation",col = colorRamp2(c(-5, 0, 5), c("#0392cf", "#fdf498","#ee4035")), 
#    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE,split = split)
#ht

# Get Cell type-specific on/off TF DEGs
Mel_specific_DEG_on <- TF_DEG2[TF_DEG2$s24vMel_log2_change <0 & TF_DEG2$s24vIri_log2_change >=0 & TF_DEG2$MelvIri_log2_change >0,] #15
Mel_specific_DEG_off <- TF_DEG2[TF_DEG2$s24vMel_log2_change >0 & TF_DEG2$s24vIri_log2_change <=0 & TF_DEG2$MelvIri_log2_change <0,] #12
Iri_specific_DEG_on <- TF_DEG2[TF_DEG2$s24vMel_log2_change >=0 & TF_DEG2$s24vIri_log2_change <0 & TF_DEG2$MelvIri_log2_change <0,] #13
Iri_specific_DEG_off <- TF_DEG2[TF_DEG2$s24vMel_log2_change <=0 & TF_DEG2$s24vIri_log2_change >0 & TF_DEG2$MelvIri_log2_change >0,] #42
Shared_Mel_Iri_DEG_on <- TF_DEG2[TF_DEG2$s24vMel_log2_change <0 & TF_DEG2$s24vIri_log2_change <0,] #63 
Shared_Mel_Iri_DEG_off <- TF_DEG2[TF_DEG2$s24vMel_log2_change >0 & TF_DEG2$s24vIri_log2_change >0,] #276

# Get Cell type-specific on/off TF DEGs using less stringent cutoff
Mel_specific_DEG_on_loose <- TF_DEG2[TF_DEG2$s24vMel_log2_change <0 & TF_DEG2$s24vIri_log2_change >=0,] #33
Mel_specific_DEG_off_loose <- TF_DEG2[TF_DEG2$s24vMel_log2_change >0 & TF_DEG2$s24vIri_log2_change <=0,] #60
Iri_specific_DEG_on_loose <- TF_DEG2[TF_DEG2$s24vMel_log2_change >=0 & TF_DEG2$s24vIri_log2_change <0,] #42
Iri_specific_DEG_off_loose <- TF_DEG2[TF_DEG2$s24vMel_log2_change <=0 & TF_DEG2$s24vIri_log2_change >0,] #102

s24_specific_DEG_on_on <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change < 0 & TF_DEG2$s24vMel_log2_change <=0 & TF_DEG2$s24vIri_log2_change <=0,] #14
s24_specific_DEG_off_off <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change > 0 & TF_DEG2$s24vMel_log2_change >=0 & TF_DEG2$s24vIri_log2_change >=0,] #47

s24_Mel_DEG_on_on <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change < 0 & TF_DEG2$s24vMel_log2_change ==0 & TF_DEG2$s24vIri_log2_change >0 & TF_DEG2$MelvIri_log2_change >=0,] #6
s24_Iri_DEG_on_on <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change < 0 & TF_DEG2$s24vIri_log2_change ==0 & TF_DEG2$s24vMel_log2_change >0 & TF_DEG2$MelvIri_log2_change <=0,] #6

s24_specific_DEG_on_off <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change < 0 & TF_DEG2$s24vMel_log2_change >0 & TF_DEG2$s24vIri_log2_change >0,] #15
s15_specific_DEG_on_off <- TF_DEG2[TF_DEG2$s15v24hpf_log2_change > 0 & TF_DEG2$s24vMel_log2_change >=0 & TF_DEG2$s24vIri_log2_change >=0,] #15

# Record to files
write.table(Mel_specific_DEG_on,"Mel_specific_TranscriptionFactor_DEG_on_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_off,"Mel_specific_TranscriptionFactor_DEG_off_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_on,"Iri_specific_TranscriptionFactor_DEG_on_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_off,"Iri_specific_TranscriptionFactor_DEG_off_stringentParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Shared_Mel_Iri_DEG_on,"Shared_Mel_Iri_TranscriptionFactor_DEG_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Shared_Mel_Iri_DEG_off,"Shared_Mel_Iri_TranscriptionFactor_DEG_off_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_on_loose,"Mel_specific_TranscriptionFactor_DEG_on_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Mel_specific_DEG_off_loose,"Mel_specific_TranscriptionFactor_DEG_off_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_on_loose,"Iri_specific_TranscriptionFactor_DEG_on_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(Iri_specific_DEG_off_loose,"Iri_specific_TranscriptionFactor_DEG_off_looseParameters_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_specific_DEG_on_on,"s24_specific_TranscriptionFactorDEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_specific_DEG_off_off,"s24_specific_TranscriptionFactor_DEG_off_off_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_Mel_DEG_on_on,"Shared_s24_Mel_TranscriptionFactorDEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_Iri_DEG_on_on,"Shared_s24_Iri_TranscriptionFactorDEG_on_on_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s24_specific_DEG_on_off,"s24_specific_TranscriptionFactorDEG_on_OFF_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")
write.table(s15_specific_DEG_on_off,"s15_specific_TranscriptionFactorDEG_on_OFF_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")

#heatmap for TFs
TFs <- read.table("TFlist_for_heatmap.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG2_genenames) <- c("gene","genename")
TFs<- merge(TFs,DEG2_genenames, by = "gene", all.x = T)
TFs<-merge(TFs,geneTPM, by = "gene", all.x =T)
TFs[is.na(TFs)] <-0

TFs[TFs$gene == "ENSDARG00000007641",]$genename.x <- "msx1a"
TFs[TFs$gene == "ENSDARG00000007641",]$genename.y <- "msx1a"
TFs[TFs$gene == "ENSDARG00000010124",]$genename.x <- "sp5l"
TFs[TFs$gene == "ENSDARG00000010124",]$genename.y <- "sp5l"
TFs[TFs$gene == "ENSDARG00000067727",]$genename.x <- "im:7142702"
TFs[TFs$gene == "ENSDARG00000067727",]$genename.y <- "im:7142702"
TFs[TFs$gene == "ENSDARG00000073856",]$genename.x <- "znf319"
TFs[TFs$gene == "ENSDARG00000073856",]$genename.y <- "znf319"
TFs[TFs$gene == "ENSDARG00000087185",]$genename.x <- "CABZ01069006.1"
TFs[TFs$gene == "ENSDARG00000087185",]$genename.y <- "CABZ01069006.1"
TFs[TFs$gene == "ENSDARG00000091086",]$genename.x <- "alx3"
TFs[TFs$gene == "ENSDARG00000091086",]$genename.y <- "alx3"


# filter low expressed TFs and write to a file
TFs<-TFs[TFs$s15 >=3 | TFs$s24 >=3 | TFs$Mel >=3| TFs$Iri >=3,] #filter low expressed TFs
TFs$X15v24 <-TFs$X15v24*-1
TFs$X24vM <-TFs$X24vM*-1
TFs$X24vI <-TFs$X24vI*-1
TFs$MvI <-TFs$MvI*-1

TFs[TFs$type == "Pigment_OFF" & TFs$MvI > 0,]$type <- "Pigment_OFF_Iri"
TFs[TFs$type == "Pigment_OFF" & TFs$MvI < 0,]$type <- "Pigment_OFF_Mel"
TFs[TFs$type == "Pigment_ON" & TFs$MvI > 0,]$type <- "Pigment_ON_Iri"
TFs[TFs$type == "Pigment_ON" & TFs$MvI < 0,]$type <- "Pigment_OFF_Iri"

write.table(TFs,"TFlist_for_heatmap_updated.txt",col.names = T, ,quote =F, sep = "\t",row.names = F)


s15spec<-TFs[TFs$type == "s15_specific",]
s15spec<-s15spec[order(s15spec$genename.x),]
P_off<-TFs[TFs$type == "Pigment_OFF",]
P_off<-P_off[order(P_off$genename.x),]
sP_on<-TFs[TFs$type == "24_Mel_Iri_ON",]
sP_on<-sP_on[order(sP_on$genename.x),]
P_on<-TFs[TFs$type == "Pigment_ON",]
P_on<-P_on[order(P_on$genename.x),]
I_off<-TFs[TFs$type == "Iri_OFF",]
I_off<-I_off[order(I_off$genename.x),]
I_on<-TFs[TFs$type == "Iri_ON",]
I_on<-I_on[order(I_on$genename.x),]
sI_on<-TFs[TFs$type == "s24_Iri_ON" | TFs$type == "Pigment_OFF_Iri" |  TFs$type == "Pigment_ON_Iri",]
sI_on<-sI_on[order(sI_on$genename.x),]
M_off<-TFs[TFs$type == "Mel_OFF",]
M_off<-M_off[order(M_off$genename.x),]
M_on<-TFs[TFs$type == "Mel_ON",]
M_on<-M_on[order(M_on$genename.x),]
sM_on<-TFs[TFs$type == "s24_Mel_ON"| TFs$type == "Pigment_OFF_Mel" |  TFs$type == "Pigment_ON_Mel",]
sM_on<-sM_on[order(sM_on$genename.x),]

All<-rbind(s15spec,P_off,sP_on,P_on,sM_on,M_on,M_off,sI_on,I_on,I_off)
NCC <- rbind(s15spec,P_off)
M <- rbind(sP_on,P_on,sM_on,M_on,M_off)
I<-rbind(sP_on,P_on,sI_on,I_on,I_off)
```
