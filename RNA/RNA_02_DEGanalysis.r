# Plot set up
source("PlotSetUp.r")

### Differential expression analysis ###
# Use DESeq2 to find differentially expressed genes in each comparison
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

s15_s24 <- subset(res_15pos_24hpf, padj < 0.01)
s24_M <- subset(res_24hpf_Mel, padj < 0.01)
s24_I <- subset(res_24hpf_Iri, padj < 0.01)
M_I <- subset(res_Mel_Iri, padj < 0.01)

#export the DEGs, Columns: gene  baseMean  log2FoldChange  lfcSE stat  pvalue  padj
write.table(as.data.frame(s15_s24),file="DEGs_15somite_24hpf_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(s24_M),file="DEGs_24hpf_Mel_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(s24_I),file="DEGs_24hpf_Iri_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")
write.table(as.data.frame(M_I),file="DEGs_Mel_Iri_p0.01.txt", col.names = F, row.names = T,quote = F, sep = "\t")

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
#heatmap based on euclidian distance
library(pheatmap)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
pheatmap(sampleDistMatrix,
	clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

#MAplots 
par(mfrow=c(2,2))
plotMA(res_15pos_24hpf, main="15somite vs 24hpf", ylim=c(-6,6), alpha = 0.001)
plotMA(res_24hpf_Mel, main="24hpf vs Melanophore", ylim=c(-6,6), alpha = 0.001)
plotMA(res_24hpf_Iri, main="24hpf vs Iridophore", ylim=c(-6,6), alpha = 0.001)
plotMA(res_Mel_Iri, main="Melanophore vs Iridophore", ylim=c(-6,6), alpha = 0.001)
## Cell type-specific MAplots
### 15somite vs 24hpf
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
DEG_s24_Mel <- as.data.frame(s24_M)
DEG_s24_Mel$gene <- rownames(s24_M)

DEG_s24_Mel$color <-0
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_on_loose$gene,]$color <- 1
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_on$gene,]$color <- 2
DEG_s24_Mel[DEG_s24_Mel$gene %in% Shared_Mel_Iri_DEG_on$gene,]$color <- 3
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_off_loose$gene,]$color <-4
DEG_s24_Mel[DEG_s24_Mel$gene %in% Mel_specific_DEG_off$gene,]$color <- 5
DEG_s24_Mel[DEG_s24_Mel$gene %in% Shared_Mel_Iri_DEG_off$gene,]$color <- 6

### 24hpf vs Mel
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

### 24hpf vs Iri
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

### Mel vs Iri
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

# MAKE master list of DEGs for comparision
## Output saved as Combined_DEGs_list.txt. Mv to Combined_DEGs_list_p0.01.txt
system('python3 All_DEGs_list.py DEGs_15somite_24hpf_p0.01.txt DEGs_24hpf_Mel_p0.01.txt DEGs_24hpf_Iri_p0.01.txt DEGs_Mel_Iri_p0.01.txt')

setwd("<path to RNA folder>/DESeq2")
genelist <- read.table("Combined_DEGs_list_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s15_s24<- read.table("DEGs_15somite_24hpf_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_M <- read.table("DEGs_24hpf_Mel_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_I <- read.table("DEGs_24hpf_Iri_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
M_I <- read.table("DEGs_Mel_Iri_p0.01.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

DEG <-Reduce(function(x, y) merge(x, y, by = "V1", all = T), list(genelist,s15_s24[,c(1:3)],s24_M[,c(1:3)],s24_I[,c(1:3)],M_I[,c(1:3)]))
colnames(DEG) <- c("gene","s15v24hpf_mean_expression","s15v24hpf_log2_change","s24vMel_mean_expression","s24vMel_log2_change","s24vIri_mean_expression","s24vIri_log2_change","MelvIri_mean_expression","MelvIri_log2_change")
DEG[is.na(DEG)]<-0

# Add genenames
## All DEGs
DEG_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG_genenames) <- c("gene","genename")
DEG <- merge(DEG,DEG_genenames, by = "gene", all.x = T)
DEG[is.na(DEG)] <-0

write.table(DEG,"DEGs_combined_samples_p0.01.txt",col.names = T, row.names = F,quote = F, sep = "\t")

### TF analysis ####
# Obtain TF list from AnimalTFDB
## Lots of TF missing in zebrafish. Combine with human TFs that have been converted using orthoretriever
TFlist <- read.table("Danio_rerio_transcription_factors_gene_list.txt",header = T, quote = "", sep = "\t",stringsAsFactors = F)
homo_TFlist <- read.table("homo_danrer_TF.txt",header = T, quote = "", sep = "\t",stringsAsFactors = F)
cofactorlist <- read.table("Danio_rerio_transcription_co-factors_gene_list.txt",header =F, quote = "", sep = "\t",stringsAsFactors = F)
homo_cofactorlist <- read.table("homo_danrer_cofactor.txt",header =F, quote = "", sep = "\t",stringsAsFactors = F)
CRMlist <- read.table("Danio_rerio_chromatin_remodeling_factors_gene_list.txt",header = F, quote = "", sep = "\t",stringsAsFactors = F)
homo_CRMlist <- read.table("homo_danrer_CHRF.txt",header = F, quote = "", sep = "\t",stringsAsFactors = F)
## Label and combine TF, CF and CRFs.
TFs <- c(TFlist[,1],homo_TFlist[,3])
TF <- unique(TFs)

CFs <- c(cofactorlist[,1],homo_cofactorlist[,3])
CF <- unique(CFs)

CRFs <- c(CRMlist[,1],homo_CRMlist[,3])
CRF <- unique(CRFs)

TFlist <- data.frame(TF, c(rep("TF",length(TF))),stringsAsFactors = F)
CFlist <- data.frame(CF, c(rep("Cofactor",length(CF))),stringsAsFactors = F)
CRFlist <- data.frame(CRF, c(rep("ChromatinRemodeler",length(CRF))),stringsAsFactors = F)

colnames(TFlist) <- c("gene", "type")
colnames(CFlist) <- c("gene", "type")
colnames(CRFlist) <- c("gene", "type")

all_list <- rbind(TFlist,CFlist,CRFlist, stringsAsFactors = F)

## Merge DEG list with TF/CF/CRF list 
DEG<- merge(DEG,all_list,by = "gene",all.x=T)
## Label DEGs not in all_list as other
DEG[is.na(DEG)] <- "other"

# Sometimes the name conversion is wrong so also make sure common name is included
homo_TFlist$Input.Common.Name <- tolower(homo_TFlist$Input.Common.Name)
homo_cofactorlist$Input.Common.Name <- tolower(homo_cofactorlist$Input.Common.Name)
homo_CRMlist$Input.Common.Name <- tolower(homo_CRMlist$Input.Common.Name)

DEG[DEG$genename %in% homo_TFlist$Input.Common.Name,]$type <- "TF"
DEG[DEG$genename %in% homo_cofactorlist$Input.Common.Name,]$type <- "Cofactor"
DEG[DEG$genename %in% homo_CRMlist$Input.Common.Name,]$type <- "ChromatinRemodeler"

write.table(DEG,"DEGs_combined_samples_p0.01_wTFinfo.txt",col.names = T, row.names = F,quote = F, sep = "\t")


### Identify DEGs for each cell type ###
# Merge all gene expression at transcript level
## Load all datasets 
genelist <- read.table("<path to RNA folder>/stringtie/Combined_gene_list_RefSeq.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

s15_R1 <- read.table("<path to RNA folder>/stringtie/RNA_15somiteNCC_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s15_R2 <- read.table("<path to RNA folder>/stringtie/RNA_15somiteNCC_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_R1 <- read.table("<path to RNA folder>/stringtie/RNA_24hpfNCC_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_R2 <- read.table("<path to RNA folder>/stringtie/RNA_24hpfNCC_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Mel_R1 <- read.table("<path to RNA folder>/stringtie/RNA_Mel_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Mel_R2 <- read.table("<path to RNA folder>/stringtie/RNA_Mel_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Iri_R1 <- read.table("<path to RNA folder>/stringtie/RNA_Iri_Rep1.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
Iri_R2 <- read.table("<path to RNA folder>/stringtie/RNA_Iri_Rep2.gene.abundance.filteredNM.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
## Merge dataframes
geneTPM<-Reduce(function(x, y) merge(x, y,by = "V1",all.x = T), list(s15_R1[,c(1,9)],s15_R2[,c(1,9)], s24_R1[,c(1,9)],s24_R2[,c(1,9)],Mel_R1[,c(1,9)],Mel_R2[,c(1,9)], Iri_R1[,c(1,9)],Iri_R2[,c(1,9)]))
## Rename columns to match with input
colnames(geneTPM) <- c("gene","s15_Rep1","s15_Rep2","s24_Rep1","s24_Rep2","Mel_Rep1","Mel_Rep2","Iri_Rep1","Iri_Rep2")
## Remove duplicated genes
geneTPM <- geneTPM[!duplicated(geneTPM$gene),]
## Calculate average value of replicates for each sample type
geneTPM$s15 <- (geneTPM$s15_Rep1+geneTPM$s15_Rep2)/2
geneTPM$s24 <- (geneTPM$s24_Rep1+geneTPM$s24_Rep2)/2
geneTPM$Mel <- (geneTPM$Mel_Rep1+geneTPM$Mel_Rep2)/2
geneTPM$Iri <- (geneTPM$Iri_Rep1+geneTPM$Iri_Rep2)/2
## Merge transcript level counts with DEG list 
# Convert Refseq mRNA ID to Ensembl gene ID then merge DEG and geneTPM
library(biomaRt)
ensembl <- useMart("ensembl", dataset="drerio_gene_ensembl")
value <- DEG$gene
DEGbm <- getBM(attributes=c("ensembl_gene_id", "refseq_mrna"), filters = "ensembl_gene_id", values = value, mart= ensembl)
# Remove rows when no Refseq ID fetched. 
neDEGbm <- DEGbm[nchar(DEGbm$refseq_mrna) != 0, ]
neDEGbm <- neDEGbm[!duplicated(neDEGbm$ensembl_gene_id),]
colnames(neDEGbm) <- c("gene","refseq_id")
DEG <- merge(DEG, neDEGbm, by = "gene")
DEG$Ensembl <- DEG$gene
DEG$gene <- DEG$refseq_id
DEG<-merge(DEG,geneTPM,by = "gene", all.x =T)
## Keep only transcripts with TPM over 5 in any of the four samples
DEG<-DEG[DEG$s15 > 5 | DEG$s24 > 5 |DEG$Mel > 5 |DEG$Iri > 5,]
DEG$gene <- DEG$Ensembl

# get cell type-specific up(on) down(off) regulated DEGs
Mel_specific_DEG_on <- DEG[DEG$s24vMel_log2_change <0 & DEG$s24vIri_log2_change >=0 & DEG$MelvIri_log2_change >0,] 
Mel_specific_DEG_off <- DEG[DEG$s24vMel_log2_change >0 & DEG$s24vIri_log2_change <=0 & DEG$MelvIri_log2_change <0,] 
Iri_specific_DEG_off <- DEG[DEG$s24vMel_log2_change <=0 & DEG$s24vIri_log2_change >0 & DEG$MelvIri_log2_change >0,] 
Shared_Mel_Iri_DEG_on <- DEG[DEG$s24vMel_log2_change <0 & DEG$s24vIri_log2_change <0,] 
Shared_Mel_Iri_DEG_off <- DEG[DEG$s24vMel_log2_change >0 & DEG$s24vIri_log2_change >0,] 

Mel_specific_DEG_on_loose <- DEG[DEG$s24vMel_log2_change <0 & DEG$s24vIri_log2_change >=0,]
Mel_specific_DEG_off_loose <- DEG[DEG$s24vMel_log2_change >0 & DEG$s24vIri_log2_change <=0,]
Iri_specific_DEG_on_loose <- DEG[DEG$s24vMel_log2_change >=0 & DEG$s24vIri_log2_change <0,]
Iri_specific_DEG_off_loose <- DEG[DEG$s24vMel_log2_change <=0 & DEG$s24vIri_log2_change >0,] 

s24_specific_DEG_on_on <- DEG[DEG$s15v24hpf_log2_change < 0 & DEG$s24vMel_log2_change <=0 & DEG$s24vIri_log2_change <=0,] 
s24_specific_DEG_off_off <- DEG[DEG$s15v24hpf_log2_change > 0 & DEG$s24vMel_log2_change >=0 & DEG$s24vIri_log2_change >=0,] 

s24_Iri_only_DEG_on_on <- DEG[DEG$s15v24hpf_log2_change < 0 & DEG$s24vMel_log2_change >0 & DEG$s24vIri_log2_change <=0,] 
s24_Mel_only_DEG_on_on <- DEG[DEG$s15v24hpf_log2_change < 0 & DEG$s24vMel_log2_change <=0 & DEG$s24vIri_log2_change >0,] 

MelvsIri_DEG_off <- DEG[DEG$MelvIri_log2_change < 0,]
MelvsIri_DEG_on <- DEG[DEG$MelvIri_log2_change > 0,]

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

# Get transcription factors that differentially expressed
## Use the geneTPM merged DEG object
TF_DEG <- DEG[DEG$type == "TF",] 
TF_DEG$Ave_RPKM <- (TF_DEG$s15+TF_DEG$s24+TF_DEG$Mel+TF_DEG$Iri)/4
## Keep only the columns with DEseq results, genename and the newly calculated Ave_RPKM
TF_DEG<- TF_DEG[,c(1:10,26)]
name <- TF_DEG[,1]

# Get Cell type-specific on/off TF DEGs
Mel_specific_DEG_on <- TF_DEG[TF_DEG$s24vMel_log2_change <0 & TF_DEG$s24vIri_log2_change >=0 & TF_DEG$MelvIri_log2_change >0,] #15
Mel_specific_DEG_off <- TF_DEG[TF_DEG$s24vMel_log2_change >0 & TF_DEG$s24vIri_log2_change <=0 & TF_DEG$MelvIri_log2_change <0,] #12
Iri_specific_DEG_on <- TF_DEG[TF_DEG$s24vMel_log2_change >=0 & TF_DEG$s24vIri_log2_change <0 & TF_DEG$MelvIri_log2_change <0,] #13
Iri_specific_DEG_off <- TF_DEG[TF_DEG$s24vMel_log2_change <=0 & TF_DEG$s24vIri_log2_change >0 & TF_DEG$MelvIri_log2_change >0,] #42
Shared_Mel_Iri_DEG_on <- TF_DEG[TF_DEG$s24vMel_log2_change <0 & TF_DEG$s24vIri_log2_change <0,] #63 
Shared_Mel_Iri_DEG_off <- TF_DEG[TF_DEG$s24vMel_log2_change >0 & TF_DEG$s24vIri_log2_change >0,] #276

# Get Cell type-specific on/off TF DEGs using less stringent cutoff
Mel_specific_DEG_on_loose <- TF_DEG[TF_DEG$s24vMel_log2_change <0 & TF_DEG$s24vIri_log2_change >=0,] #33
Mel_specific_DEG_off_loose <- TF_DEG[TF_DEG$s24vMel_log2_change >0 & TF_DEG$s24vIri_log2_change <=0,] #60
Iri_specific_DEG_on_loose <- TF_DEG[TF_DEG$s24vMel_log2_change >=0 & TF_DEG$s24vIri_log2_change <0,] #42
Iri_specific_DEG_off_loose <- TF_DEG[TF_DEG$s24vMel_log2_change <=0 & TF_DEG$s24vIri_log2_change >0,] #102

s24_specific_DEG_on_on <- TF_DEG[TF_DEG$s15v24hpf_log2_change < 0 & TF_DEG$s24vMel_log2_change <=0 & TF_DEG$s24vIri_log2_change <=0,] #14
s24_specific_DEG_off_off <- TF_DEG[TF_DEG$s15v24hpf_log2_change > 0 & TF_DEG$s24vMel_log2_change >=0 & TF_DEG$s24vIri_log2_change >=0,] #47

s24_Mel_DEG_on_on <- TF_DEG[TF_DEG$s15v24hpf_log2_change < 0 & TF_DEG$s24vMel_log2_change ==0 & TF_DEG$s24vIri_log2_change >0 & TF_DEG$MelvIri_log2_change >=0,] #6
s24_Iri_DEG_on_on <- TF_DEG[TF_DEG$s15v24hpf_log2_change < 0 & TF_DEG$s24vIri_log2_change ==0 & TF_DEG$s24vMel_log2_change >0 & TF_DEG$MelvIri_log2_change <=0,] #6

s24_specific_DEG_on_off <- TF_DEG[TF_DEG$s15v24hpf_log2_change < 0 & TF_DEG$s24vMel_log2_change >0 & TF_DEG$s24vIri_log2_change >0,] #15
s15_specific_DEG_on_off <- TF_DEG[TF_DEG$s15v24hpf_log2_change > 0 & TF_DEG$s24vMel_log2_change >=0 & TF_DEG$s24vIri_log2_change >=0,] #15

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

# Update TF list
TFs <- read.table("TFlist_for_heatmap.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F)
DEG_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG_genenames) <- c("gene","genename")
TFs<- merge(TFs,DEG_genenames, by = "gene", all.x = T)
value <- TFs$gene
TFsbm <- getBM(attributes=c("ensembl_gene_id", "refseq_mrna"), filters = "ensembl_gene_id", values = value, mart= ensembl)
# Remove rows when no Refseq ID fetched
neTFsbm <- TFsbm[nchar(TFsbm$refseq_mrna) != 0, ]
colnames(neTFsbm) <- c("gene","refseq_id")
TFs <- merge(TFs, neTFsbm, by = "gene")
TFs$Ensembl <- TFs$gene
TFs$gene <- TFs$refseq_id
TFs<- merge(TFs,geneTPM, by = "gene", all.x =T)
TFs[is.na(TFs)] <-0
TFs$gene <- TFs$Ensembl

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
