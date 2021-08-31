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

### Analysis script
```{R}
###################################### {R script} ######################################

######## COMBINING WGBS METHYLATION DATA (d30) AND ATAC SEQ DATA (IDR DARs) ###########
#Compare DMR size vs DAR size
combined_DMR <- read.table("Combined_DMRs_d30_p0.01_wINFO.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
combined_DAR <- read.table("All_DiffBind_peaks_merged_wINFO.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
#32891 combined DMRs vs 55380  combined DARs

d <- data.frame(c(rep("DMR",nrow(combined_DMR)),rep("DAR",nrow(combined_DAR))),c(combined_DMR$V5,combined_DAR$V5))
colnames(d) <- c("Type","Size")
a <-ggplot(d, aes(Size, fill = Type, colour = Type)) + geom_density(alpha = 0.5, adjust = 1)+ggtitle(paste("Combined DMR (d30) & Combined DAR size distribution"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Size (bp)", y = "Density")+scale_x_continuous(limits = c(0,2000))+scale_fill_manual(values = mypalette2[c(1,3)])+scale_color_manual(values=mypalette2[c(1,3)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a
```
```{bash}
###################################### {Bash script} ######################################
#combine DMRs and DARs bed files to make a master list
# combined DMRs and DARs
cat Combined_DMRs_d30_p0.01_wINFO.bed All_DiffBind_peaks_merged_wINFO.bed >Combined_DMR_DAR.bed
# sort based on chrom
sort -k1,1 -k2,2n Combined_DMR_DAR.bed > Combined_DMR_DAR.sorted.bed
# merge regions
bedtools merge -i Combined_DMR_DAR.sorted.bed > Combined_DMR_DAR.bed

rm Combined_DMR_DAR.sorted.bed

# Add DMR and DAR information to the combined DMR_DAR bed file
bedtools intersect -wao -a Combined_DMR_DAR.bed -b Combined_DMRs_d30_p0.01_wINFO.bed > Combined_DMR_DAR_wDMRinfo.bed
bedtools intersect -wao -a Combined_DMR_DAR.bed -b All_DiffBind_peaks_merged_wINFO.bed > Combined_DMR_DAR_wDARinfo.bed

# Find in region peaks and count overlapps
bedtools intersect -c -a Combined_DMR_DAR.bed -b 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_15s.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_24hpf.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_Mel.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_Iri.txt

#calculate average methylation of DMARs
for i in "SAMPLE_Combined_DSS.txt"; do sed '1d' $i | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,1}' > ${i/.txt/.bed};done
bedtools intersect -wo -a Combined_DMR_DAR.bed -b 15somite_NCC_Combined_DSS.bed > Combined_DMR_DAR_wMETH_15s.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b 24hpf_NCC_Combined_DSS.bed > Combined_DMR_DAR_wMETH_24hpf.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b Mel_Combined_DSS.bed > Combined_DMR_DAR_wMETH_Mel.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b Iri_Combined_DSS.bed > Combined_DMR_DAR_wMETH_Iri.txt
```
```{R}
###################################### {R script} ######################################
# Pull out regions with methylation information from each samples
DMAR_meth15<- read.table("Combined_DMR_DAR_wMETH_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_meth24<- read.table("Combined_DMR_DAR_wMETH_24hpf.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_methM<- read.table("Combined_DMR_DAR_wMETH_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_methI<- read.table("Combined_DMR_DAR_wMETH_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

# Get regions with CpGs with >=5 coverage
DMAR_meth15 <- DMAR_meth15[DMAR_meth15$V7 >=5,]#filter CpGs with >=5 cov
DMAR_meth24 <- DMAR_meth24[DMAR_meth24$V7 >=5,]#filter CpGs with >=5 cov
DMAR_methM <- DMAR_methM[DMAR_methM$V7 >=5,]#filter CpGs with >=5 cov
DMAR_methI <- DMAR_methI[DMAR_methI$V7 >=5,]#filter CpGs with >=5 cov

# Convert to methylation percentage
DMAR_meth15$meth <- DMAR_meth15$V8/DMAR_meth15$V7*100
DMAR_meth24$meth <- DMAR_meth24$V8/DMAR_meth24$V7*100
DMAR_methM$meth <- DMAR_methM$V8/DMAR_methM$V7*100
DMAR_methI$meth <- DMAR_methI$V8/DMAR_methI$V7*100

#Aggregate CpG counts and calculate average methylation
library(plyr)
DMAR_meth15_aggregate <- ddply(DMAR_meth15,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth))
DMAR_meth24_aggregate <- ddply(DMAR_meth24,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth))
DMAR_methM_aggregate <- ddply(DMAR_methM,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth))
DMAR_methI_aggregate <- ddply(DMAR_methI,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth))

#parse data together
#  Read combined region information
DMAR_DAR<- read.table("Combined_DMR_DAR_wDARinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_DMR<- read.table("Combined_DMR_DAR_wDMRinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR15<- read.table("Combined_DMR_DAR_wIDRinfo_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR24<- read.table("Combined_DMR_DAR_wIDRinfo_24hpf.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRM<- read.table("Combined_DMR_DAR_wIDRinfo_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRI<- read.table("Combined_DMR_DAR_wIDRinfo_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

#Add chromosome position information
DMAR_DAR$chrompos<- paste(DMAR_DAR$V1,":",DMAR_DAR$V2,"-",DMAR_DAR$V3,sep = "")
DMAR_DMR$chrompos<- paste(DMAR_DMR$V1,":",DMAR_DMR$V2,"-",DMAR_DMR$V3,sep = "")
DMAR_IDR15$chrompos<- paste(DMAR_IDR15$V1,":",DMAR_IDR15$V2,"-",DMAR_IDR15$V3,sep = "")
DMAR_IDR24$chrompos<- paste(DMAR_IDR24$V1,":",DMAR_IDR24$V2,"-",DMAR_IDR24$V3,sep = "")
DMAR_IDRM$chrompos<- paste(DMAR_IDRM$V1,":",DMAR_IDRM$V2,"-",DMAR_IDRM$V3,sep = "")
DMAR_IDRI$chrompos<- paste(DMAR_IDRI$V1,":",DMAR_IDRI$V2,"-",DMAR_IDRI$V3,sep = "")

DMAR_meth15_aggregate$chrompos<- paste(DMAR_meth15_aggregate$V1,":",DMAR_meth15_aggregate$V2,"-",DMAR_meth15_aggregate$V3,sep = "")
DMAR_meth24_aggregate$chrompos<- paste(DMAR_meth24_aggregate$V1,":",DMAR_meth24_aggregate$V2,"-",DMAR_meth24_aggregate$V3,sep = "")
DMAR_methM_aggregate$chrompos<- paste(DMAR_methM_aggregate$V1,":",DMAR_methM_aggregate$V2,"-",DMAR_methM_aggregate$V3,sep = "")
DMAR_methI_aggregate$chrompos<- paste(DMAR_methI_aggregate$V1,":",DMAR_methI_aggregate$V2,"-",DMAR_methI_aggregate$V3,sep = "")

# Combine all information into one table
DMAR <-Reduce(function(x, y) merge(x, y,by = "chrompos", all.x = T), list(DMAR_DMR[,c("chrompos","V1","V2","V3","V8","V9","V10","V11","V12")],DMAR_DAR[,c("chrompos","V8","V9","V10","V11","V12")],DMAR_IDR15[,c("chrompos","V4")],DMAR_IDR24[,c("chrompos","V4")],DMAR_IDRM[,c("chrompos","V4")],DMAR_IDRI[,c("chrompos","V4")],DMAR_meth15_aggregate[,c("chrompos","AveMeth")],DMAR_meth24_aggregate[,c("chrompos","AveMeth")],DMAR_methM_aggregate[,c("chrompos","AveMeth")],DMAR_methI_aggregate[,c("chrompos","AveMeth")],DMAR_meth15_aggregate[,c("chrompos","CpGcount")],DMAR_meth24_aggregate[,c("chrompos","CpGcount")],DMAR_methM_aggregate[,c("chrompos","CpGcount")],DMAR_methI_aggregate[,c("chrompos","CpGcount")]))
colnames(DMAR) <- c("chrompos","chr","start","end","DMRsize","DMRs15vs24","DMRs24vMel","DMRs24vIri","DMRMelvIri","DARsize","DARs15vs24","DARs24vMel","DARs24vIri","DARMelvIri","IDR_s15","IDR_s24","IDR_M","IDR_I","Meth_s15","Meth_s24","Meth_Mel","Meth_Iri","CpG_s15","CpG_s24","CpG_Mel","CpG_Iri")
DMAR$DMARsize <- DMAR$end-DMAR$start
DMAR[is.na(DMAR)] <- 0
DMAR[DMAR == "."] <- 0
DMAR[DMAR == "-1"] <- 0
DMAR$DMRsize <- as.numeric(DMAR$DMRsize)
DMAR$DARsize <- as.numeric(DMAR$DARsize)
DMAR$DMRs15vs24 <- as.numeric(DMAR$DMRs15vs24 )
DMAR$DMRs24vMel <- as.numeric(DMAR$DMRs24vMel )
DMAR$DMRs24vIri <- as.numeric(DMAR$DMRs24vIri )
DMAR$DMRMelvIri <- as.numeric(DMAR$DMRMelvIri)
DMAR$DARs15vs24 <- as.numeric(DMAR$DARs15vs24 )
DMAR$DARs24vMel <- as.numeric(DMAR$DARs24vMel )
DMAR$DARs24vIri <- as.numeric(DMAR$DARs24vIri )
DMAR$DARMelvIri <- as.numeric(DMAR$DARMelvIri)

DMAR <- DMAR[!duplicated(DMAR$chrompos),]

DMAR <- DMAR[,c(2,3,4,1,27,5:26)]
DMAR$aveCpG <- apply(DMAR,1,function(x) {(as.numeric(x[24])+as.numeric(x[25])+as.numeric(x[26])+as.numeric(x[27]))/4})
DMAR$CpGdensity100bp <- apply(DMAR,1,function(x) {(as.numeric(x[28])/as.numeric(x[5]))*100})
write.table(DMAR, "All_DMAR_Combined_wINFO.bed", row.names = F, col.names = T, sep = "\t",quote =F)
```
```{bash}
########### Annotate Peaks ##########
###################################### {Bash script} ####################################
# annotatePeaks.pl is from HOMER （http://homer.ucsd.edu/homer/ngs/annotation.html）
annotatePeaks.pl All_DMAR_Combined_wINFO.bed danRer10 -annStats All_DMAR_Combined_wINFO.annstats > All_DMAR_Combined_wINFO.annotated.txt
annotatePeaks.pl All_DiffBind_peaks_merged.txt danRer10 -annStats All_DiffBind_peaks_merged.annstats > All_DiffBind_peaks_merged.annotated.txt
annotatePeaks.pl Combined_DMRs_d30_p0.01.txt danRer10 -annStats Combined_DMRs_d30_p0.01.annstats > Combined_DMRs_d30_p0.01.annotated.txt
```
```{R}
###################################### {R script} ######################################
DMAR <- read.table("All_DMAR_Combined_wINFO.bed", sep = "\t", header = T, quote = "",stringsAsFactors = F)
anno_DMAR <- read.delim("All_DMAR_Combined_wINFO.annotated.txt",skip =1, header = F, quote = "",stringsAsFactors = F)
colnames(anno_DMAR)[1] <- "chrompos"
DMAR <- merge(DMAR,anno_DMAR[,c(1,8,9)],by = "chrompos",all.x = T)
colnames(DMAR)[30] <- "Annotation"
DMAR$Annotation <- apply(DMAR, 1, function(x) {unlist(strsplit(x[30], " \\("))[1]})
DMAR$V9 <- apply(DMAR, 1, function(x) {unlist(strsplit(x[31], " \\("))[1]})
colnames(DMAR)[31] <- "Annotation2"

DMAR <- DMAR[,c(2,3,4,1,5:31)]
write.table(DMAR, "All_DMAR_Combined_wINFO.annotated.bed", row.names = F, col.names = T, sep = "\t",quote =F)
# without header
write.table(DMAR, "All_DMAR_Combined_wINFO.annotated2.bed", row.names = F, col.names = F, sep = "\t",quote =F)


#Plot size distribution of DAR, DMR, DMAR
d <- data.frame(c(rep("DMR",nrow(DMAR[DMAR$DMRsize > 0,])),rep("DAR",nrow(DMAR[DMAR$DARsize > 0,])),rep("DMAR",nrow(DMAR))),c(DMAR[DMAR$DMRsize > 0,]$DMRsize,DMAR[DMAR$DARsize > 0,]$DARsize,DMAR$DMARsize))
colnames(d) <- c("Type","Size")
a <-ggplot(d, aes(Size, fill = Type, colour = Type)) + geom_density(alpha = 0.5, adjust = 1)+ggtitle(paste("DMAR size distribution"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Size (bp)", y = "Density")+scale_x_continuous(limits = c(0,2000))+scale_fill_manual(values = mypalette2[c(1,2,3)])+scale_color_manual(values=mypalette2[c(1,2,3)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a


#Identify different classes of DMARs
Only_DARs <- DMAR[DMAR$DMRsize == 0 & DMAR$DARsize >0,] #39669
Only_DMRs <- DMAR[DMAR$DARsize == 0 & DMAR$DMRsize >0,] #17463

str(Only_DMRs[Only_DMRs$DMRs24vMel <0 | Only_DMRs$DMRs24vIri <0,])#650
str(DMAR[DMAR$DMRs24vMel <0 | DMAR$DMRs24vIri <0,])#695    ~93.5% hyperDMRs are solitary hyper DMRS and occur in regions with no peaks (650-7-9-15)/650*100 = >95%

#DMR specific
Only_DMRs_15v24_hypo <-  Only_DMRs[Only_DMRs$DMRs15vs24>0,] #287
Only_DMRs_15v24_hyper <-  Only_DMRs[Only_DMRs$DMRs15vs24<0,] #155

Only_DMRs_Mel_specific_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel >0 & Only_DMRs$DMRs24vIri <= 0,]  #6701 
Only_DMRs_Mel_specific_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel <0 & Only_DMRs$DMRs24vIri >= 0,] #210

Only_DMRs_Iri_specific_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel <=0 & Only_DMRs$DMRs24vIri > 0,] #5672
Only_DMRs_Iri_specific_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel >= 0 & Only_DMRs$DMRs24vIri < 0,] #289

Only_DMRs_M_I_shared_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel >0 & Only_DMRs$DMRs24vIri > 0,]#3886
Only_DMRs_M_I_shared_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel <0 & Only_DMRs$DMRs24vIri < 0,]#151

write.table(Only_DMRs_Mel_specific_hypo, "Solo_DMRs_Mel_specific_hypo.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DMRs_Mel_specific_hyper, "Solo_DMRs_Mel_specific_hyper.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DMRs_Iri_specific_hypo, "Solo_DMRs_Iri_specific_hypo.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DMRs_Iri_specific_hyper, "Solo_DMRs_Iri_specific_hyper.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DMRs_M_I_shared_hypo, "Solo_DMRs_Mel_Iri_shared_hypo.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DMRs_M_I_shared_hyper, "Solo_DMRs_Mel_Iri_shared_hyper.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Only_DMRs_Mel_specific_hypo_00 <- Only_DMRs_Mel_specific_hypo[Only_DMRs_Mel_specific_hypo$IDR_s24 ==0 & Only_DMRs_Mel_specific_hypo$IDR_M ==0,] #1993
Only_DMRs_Mel_specific_hypo_10 <- Only_DMRs_Mel_specific_hypo[Only_DMRs_Mel_specific_hypo$IDR_s24 >0 & Only_DMRs_Mel_specific_hypo$IDR_M ==0,] #731
Only_DMRs_Mel_specific_hypo_01 <- Only_DMRs_Mel_specific_hypo[Only_DMRs_Mel_specific_hypo$IDR_s24 ==0 & Only_DMRs_Mel_specific_hypo$IDR_M >0,] #1607
Only_DMRs_Mel_specific_hypo_11 <- Only_DMRs_Mel_specific_hypo[Only_DMRs_Mel_specific_hypo$IDR_s24 >0 & Only_DMRs_Mel_specific_hypo$IDR_M >0,] #2370 

write.table(Only_DMRs_Mel_specific_hypo_11, "Solo_DMRs_Mel_specific_hypo_alwaysOPEN.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Only_DMRs_Iri_specific_hypo_00 <- Only_DMRs_Iri_specific_hypo[Only_DMRs_Iri_specific_hypo$IDR_s24 ==0 & Only_DMRs_Iri_specific_hypo$IDR_I ==0,] #1983
Only_DMRs_Iri_specific_hypo_10 <- Only_DMRs_Iri_specific_hypo[Only_DMRs_Iri_specific_hypo$IDR_s24 >0 & Only_DMRs_Iri_specific_hypo$IDR_I ==0,] #335
Only_DMRs_Iri_specific_hypo_01 <- Only_DMRs_Iri_specific_hypo[Only_DMRs_Iri_specific_hypo$IDR_s24 ==0 & Only_DMRs_Iri_specific_hypo$IDR_I >0,] #1711
Only_DMRs_Iri_specific_hypo_11 <- Only_DMRs_Iri_specific_hypo[Only_DMRs_Iri_specific_hypo$IDR_s24 >0 & Only_DMRs_Iri_specific_hypo$IDR_I >0,] #1643

write.table(Only_DMRs_Iri_specific_hypo_11, "Solo_DMRs_Iri_specific_hypo_alwaysOPEN.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Only_DMRs_Mel_specific_hyper_00 <- Only_DMRs_Mel_specific_hyper[Only_DMRs_Mel_specific_hyper$IDR_s24 ==0 & Only_DMRs_Mel_specific_hyper$IDR_M ==0,] #193
Only_DMRs_Mel_specific_hyper_10 <- Only_DMRs_Mel_specific_hyper[Only_DMRs_Mel_specific_hyper$IDR_s24 >0 & Only_DMRs_Mel_specific_hyper$IDR_M ==0,] #10
Only_DMRs_Mel_specific_hyper_01 <- Only_DMRs_Mel_specific_hyper[Only_DMRs_Mel_specific_hyper$IDR_s24 ==0 & Only_DMRs_Mel_specific_hyper$IDR_M >0,] #0
Only_DMRs_Mel_specific_hyper_11 <- Only_DMRs_Mel_specific_hyper[Only_DMRs_Mel_specific_hyper$IDR_s24 >0 & Only_DMRs_Mel_specific_hyper$IDR_M >0,] #7 

Only_DMRs_Iri_specific_hyper_00 <- Only_DMRs_Iri_specific_hyper[Only_DMRs_Iri_specific_hyper$IDR_s24 ==0 & Only_DMRs_Iri_specific_hyper$IDR_I ==0,] #269
Only_DMRs_Iri_specific_hyper_10 <- Only_DMRs_Iri_specific_hyper[Only_DMRs_Iri_specific_hyper$IDR_s24 >0 & Only_DMRs_Iri_specific_hyper$IDR_I ==0,] #10
Only_DMRs_Iri_specific_hyper_01 <- Only_DMRs_Iri_specific_hyper[Only_DMRs_Iri_specific_hyper$IDR_s24 ==0 & Only_DMRs_Iri_specific_hyper$IDR_I >0,] #1
Only_DMRs_Iri_specific_hyper_11 <- Only_DMRs_Iri_specific_hyper[Only_DMRs_Iri_specific_hyper$IDR_s24 >0 & Only_DMRs_Iri_specific_hyper$IDR_I >0,] #9

Only_DMRs_M_I_shared_hypo_000 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 ==0 & Only_DMRs_M_I_shared_hypo$IDR_I ==0 & Only_DMRs_M_I_shared_hypo$IDR_M ==0,] #338
Only_DMRs_M_I_shared_hypo_100 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 >0 & Only_DMRs_M_I_shared_hypo$IDR_I ==0 & Only_DMRs_M_I_shared_hypo$IDR_M ==0,] #137
Only_DMRs_M_I_shared_hypo_001 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 ==0 & Only_DMRs_M_I_shared_hypo$IDR_I ==0 & Only_DMRs_M_I_shared_hypo$IDR_M > 0,] #108
Only_DMRs_M_I_shared_hypo_010 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 ==0 & Only_DMRs_M_I_shared_hypo$IDR_I >0 & Only_DMRs_M_I_shared_hypo$IDR_M == 0,] #163
Only_DMRs_M_I_shared_hypo_011 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 ==0 & Only_DMRs_M_I_shared_hypo$IDR_I >0 & Only_DMRs_M_I_shared_hypo$IDR_M > 0,] #354
Only_DMRs_M_I_shared_hypo_111 <- Only_DMRs_M_I_shared_hypo[Only_DMRs_M_I_shared_hypo$IDR_s24 >0 & Only_DMRs_M_I_shared_hypo$IDR_I >0& Only_DMRs_M_I_shared_hypo$IDR_M >0,] #2416

write.table(Only_DMRs_Iri_specific_hypo_11, "Solo_DMRs_Mel_Iri_shared_hypo_alwaysOPEN.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Only_DMRs_M_I_shared_hyper_000 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M ==0,] #110
Only_DMRs_M_I_shared_hyper_100 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 >0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M ==0,] #17
Only_DMRs_M_I_shared_hyper_001 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M > 0,] #0
Only_DMRs_M_I_shared_hyper_010 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I >0 & Only_DMRs_M_I_shared_hyper$IDR_M == 0,] #1
Only_DMRs_M_I_shared_hyper_011 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I >0 & Only_DMRs_M_I_shared_hyper$IDR_M > 0,] #2
Only_DMRs_M_I_shared_hyper_111 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 >0 & Only_DMRs_M_I_shared_hyper$IDR_I >0& Only_DMRs_M_I_shared_hyper$IDR_M >0,] #15

#DAR specific
Only_DARs_15v24_opening <-  Only_DARs[Only_DARs$DARs15vs24>0,] #287
Only_DARs_15v24_closing <-  Only_DARs[Only_DARs$DARs15vs24<0,] #155

Only_DARs_Mel_specific_open <- Only_DARs[Only_DARs$DARs24vMel >0 & Only_DARs$DARs24vIri <= 0,]  #3605
Only_DARs_Mel_specific_closed <- Only_DARs[Only_DARs$DARs24vMel <0 & Only_DARs$DARs24vIri >= 0,] #8807

str(Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$Meth_Mel <30,]) #791
str(Only_DARs_Mel_specific_closed[Only_DARs_Mel_specific_closed$Meth_Mel <30,]) #3752

Only_DARs_Iri_specific_open <- Only_DARs[Only_DARs$DARs24vMel <=0 & Only_DARs$DARs24vIri > 0,] #3645
Only_DARs_Iri_specific_closed <- Only_DARs[Only_DARs$DARs24vMel >= 0 & Only_DARs$DARs24vIri < 0,] #6043

str(Only_DARs_Iri_specific_open[Only_DARs_Iri_specific_open$Meth_Iri <30,]) #1453
str(Only_DARs_Iri_specific_closed[Only_DARs_Iri_specific_closed$Meth_Iri <30,])#911


Only_DARs_M_I_shared_open <- Only_DARs[Only_DARs$DARs24vMel >0 & Only_DARs$DARs24vIri > 0,]#589
Only_DARs_M_I_shared_closed <- Only_DARs[Only_DARs$DARs24vMel <0 & Only_DARs$DARs24vIri < 0,]#15058

write.table(Only_DARs_Mel_specific_open, "Solo_DARs_Mel_specific_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Mel_specific_closed, "Solo_DARs_Mel_specific_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Iri_specific_open, "Solo_DARs_Iri_specific_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Iri_specific_closed, "Solo_DARs_Iri_specific_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_M_I_shared_open, "Solo_DARs_Mel_Iri_shared_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_M_I_shared_closed, "Solo_DARs_Mel_Iri_shared_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Functional_Only_DARs_Mel_specific_open <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$Meth_Mel >0 & Only_DARs_Mel_specific_open$Meth_Mel <= 30 & Only_DARs_Mel_specific_open$Meth_s24 >0 & Only_DARs_Mel_specific_open$Meth_s24 <=30,]
as.data.frame(table(Functional_Only_DARs_Mel_specific_open$Annotation))

#make heatmap Mel specific DAR vs meth
library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- Only_DARs_Mel_specific_open[,4]
df.OG <- Only_DARs_Mel_specific_open[,c(21:23)]
row.names(df.OG) <- name

kclus <- kmeans(df.OG,10)
split <-  kclus$cluster
ht = Heatmap(df.OG, column_title = "Mel-specific only DAR (closing)",name= "Methylation",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE,split = split)
ht

order <- row_order(ht)

c1 <- t(t(row.names(df.OG[order[[1]],])))
c2 <- t(t(row.names(df.OG[order[[2]],])))
c3 <- t(t(row.names(df.OG[order[[3]],])))
c4 <- t(t(row.names(df.OG[order[[4]],])))
c5 <- t(t(row.names(df.OG[order[[5]],])))
c6 <- t(t(row.names(df.OG[order[[6]],])))
c7 <- t(t(row.names(df.OG[order[[7]],])))
c8 <- t(t(row.names(df.OG[order[[8]],])))
c9 <- t(t(row.names(df.OG[order[[9]],])))
c10 <- t(t(row.names(df.OG[order[[10]],])))


Only_DARs_Mel_specific_open_c1 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c1,]
Only_DARs_Mel_specific_open_c2 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c2,]
Only_DARs_Mel_specific_open_c3 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c3,]
Only_DARs_Mel_specific_open_c4 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c4,]
Only_DARs_Mel_specific_open_c5 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c5,]
Only_DARs_Mel_specific_open_c6 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c6,]
Only_DARs_Mel_specific_open_c7 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c7,]
Only_DARs_Mel_specific_open_c8 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c8,]
Only_DARs_Mel_specific_open_c9 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c9,]
Only_DARs_Mel_specific_open_c10 <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$chrompos %in% c10,]

#check if there's any annotation enrichment?
t <- data.frame(c("non-coding","Intergenic","promoter-TSS","5' UTR","exon","intron","3' UTR","TTS"))
colnames(t) <- c("Var1")
tt <-Reduce(function(x, y) merge(x, y,by = "Var1",all.x=T), list(t,
as.data.frame(table(Only_DARs_Mel_specific_open_c1$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c2$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c3$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c4$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c5$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c6$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c7$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c8$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c9$Annotation)),
as.data.frame(table(Only_DARs_Mel_specific_open_c10$Annotation))))

colnames(tt) <- c("Annotation","k1","k2","k3","k4","k5","k6","k7","k8","k9","k10")
tt[is.na(tt)]<-0
tt$k1 <- tt$k1/sum(tt$k1)*100
tt$k2 <- tt$k2/sum(tt$k2)*100
tt$k3 <- tt$k3/sum(tt$k3)*100
tt$k4 <- tt$k4/sum(tt$k4)*100
tt$k5 <- tt$k5/sum(tt$k5)*100
tt$k6 <- tt$k6/sum(tt$k6)*100
tt$k7 <- tt$k7/sum(tt$k7)*100
tt$k8 <- tt$k8/sum(tt$k8)*100
tt$k9 <- tt$k9/sum(tt$k9)*100
tt$k10 <- tt$k10/sum(tt$k10)*100

kmelt <- melt(tt,id.vars = "Annotation")
colnames(kmelt)[2]<- "Type"
kmelt$Type <- factor(kmelt$Type, levels = c("k1","k2","k3","k4","k5","k6","k7","k8","k9","k10"))
kmelt$Annotation <- factor(kmelt$Annotation, levels = c("non-coding","Intergenic","promoter-TSS","5' UTR","exon","intron","3' UTR","TTS"))


mypalette <- c(brewer.pal(9,"Pastel1"),brewer.pal(8,"Pastel2"))
p <- ggplot(kmelt, aes(x=Annotation, y=value, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Annotation of kmean cluster")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Annotation", y = "Percent of peaks")+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p



#CpG density distribution in epigenetically dynamic regions
mypalette <- brewer.pal(8,"Pastel1")
d<- data.frame(c(rep("k1",nrow(Only_DARs_Mel_specific_closed_c1)),rep("k2",nrow(Only_DARs_Mel_specific_closed_c2)),rep("k3",nrow(Only_DARs_Mel_specific_closed_c3)),rep("k4",nrow(Only_DARs_Mel_specific_closed_c4)),rep("k5",nrow(Only_DARs_Mel_specific_closed_c5)),rep("k6",nrow(Only_DARs_Mel_specific_closed_c6)),rep("k7",nrow(Only_DARs_Mel_specific_closed_c7)),rep("k8",nrow(Only_DARs_Mel_specific_closed_c8)),rep("k9",nrow(Only_DARs_Mel_specific_closed_c9)),rep("k10",nrow(Only_DARs_Mel_specific_closed_c10))),c(Only_DARs_Mel_specific_closed_c1$CpGdensity100bp, Only_DARs_Mel_specific_closed_c2$CpGdensity100bp, Only_DARs_Mel_specific_closed_c3$CpGdensity100bp, Only_DARs_Mel_specific_closed_c4$CpGdensity100bp, Only_DARs_Mel_specific_closed_c5$CpGdensity100bp, Only_DARs_Mel_specific_closed_c6$CpGdensity100bp, Only_DARs_Mel_specific_closed_c7$CpGdensity100bp, Only_DARs_Mel_specific_closed_c8$CpGdensity100bp, Only_DARs_Mel_specific_closed_c9$CpGdensity100bp, Only_DARs_Mel_specific_closed_c10$CpGdensity100bp))
colnames(d) <- c("kmean","CpGdensity100bp")
a <-ggplot(d, aes(CpGdensity100bp, fill = kmean, colour = kmean)) + geom_density(alpha = 0.2, adjust = 1)+ggtitle(paste("CpGs density (per 100bp) in solitary DARs (Mel-specific closeding)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Ave CpG Count", y = "Density")+scale_x_continuous(limits = c(0,10))+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a


#heatmap of Average dna methylation in pigment cell-specific DARs
ht_order <- c()

OG <- Only_DARs_Mel_specific_open 
name <- OG[,4]
df.OG <- OG[,c(2:16)]
row.names(df.OG) <- name
setdiff(OG$chrompos,ht_order)
OG$chrompos <- factor(OG$chrompos, levels = c(ht_order))
OG2 <- OG[order(OG$chrompos),]
name <- OG2[,1]
df.OG2 <- OG2[,c(21:23)]
row.names(df.OG2) <- name
ht = Heatmap(df.OG2, col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE)
ht


#dynamic DMARs#
dynamic_DMARs <- DMAR[DMAR$DMRsize != 0 & DMAR$DARsize !=0,] #15252

# Make DMAR file and double check data set
DMAR_DAR<- read.table("Combined_DMR_DAR_wDARinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_DMR<- read.table("Combined_DMR_DAR_wDMRinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR15<- read.table("Combined_DMR_DAR_wIDRinfo_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR24<- read.table("Combined_DMR_DAR_wIDRinfo_24hpf.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRM<- read.table("Combined_DMR_DAR_wIDRinfo_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRI<- read.table("Combined_DMR_DAR_wIDRinfo_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

DMAR_DAR$chrompos<- paste(DMAR_DAR$V1,":",DMAR_DAR$V2,"-",DMAR_DAR$V3,sep = "")
DMAR_DMR$chrompos<- paste(DMAR_DMR$V1,":",DMAR_DMR$V2,"-",DMAR_DMR$V3,sep = "")
DMAR_IDR15$chrompos<- paste(DMAR_IDR15$V1,":",DMAR_IDR15$V2,"-",DMAR_IDR15$V3,sep = "")
DMAR_IDR24$chrompos<- paste(DMAR_IDR24$V1,":",DMAR_IDR24$V2,"-",DMAR_IDR24$V3,sep = "")
DMAR_IDRM$chrompos<- paste(DMAR_IDRM$V1,":",DMAR_IDRM$V2,"-",DMAR_IDRM$V3,sep = "")
DMAR_IDRI$chrompos<- paste(DMAR_IDRI$V1,":",DMAR_IDRI$V2,"-",DMAR_IDRI$V3,sep = "")


#methylation value, filter out DMRs with less than ave 5 cov
DMAR_meth15<- read.table("Combined_DMR_DAR_wMETH_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_meth24<- read.table("Combined_DMR_DAR_wMETH_24hpf.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_methM<- read.table("Combined_DMR_DAR_wMETH_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_methI<- read.table("Combined_DMR_DAR_wMETH_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

DMAR_meth15$meth <- DMAR_meth15$V8/DMAR_meth15$V7*100
DMAR_meth24$meth <- DMAR_meth24$V8/DMAR_meth24$V7*100
DMAR_methM$meth <- DMAR_methM$V8/DMAR_methM$V7*100
DMAR_methI$meth <- DMAR_methI$V8/DMAR_methI$V7*100

library(plyr)
DMAR_meth15_aggregate <- ddply(DMAR_meth15,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth),AveCov=mean(V7))
DMAR_meth24_aggregate <- ddply(DMAR_meth24,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth),AveCov=mean(V7))
DMAR_methM_aggregate <- ddply(DMAR_methM,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth),AveCov=mean(V7))
DMAR_methI_aggregate <- ddply(DMAR_methI,.(V1,V2,V3),summarise,CpGcount=sum(V9),AveMeth=mean(meth),AveCov=mean(V7))

DMAR_meth15_aggregate$chrompos<- paste(DMAR_meth15_aggregate$V1,":",DMAR_meth15_aggregate$V2,"-",DMAR_meth15_aggregate$V3,sep = "")
DMAR_meth24_aggregate$chrompos<- paste(DMAR_meth24_aggregate$V1,":",DMAR_meth24_aggregate$V2,"-",DMAR_meth24_aggregate$V3,sep = "")
DMAR_methM_aggregate$chrompos<- paste(DMAR_methM_aggregate$V1,":",DMAR_methM_aggregate$V2,"-",DMAR_methM_aggregate$V3,sep = "")
DMAR_methI_aggregate$chrompos<- paste(DMAR_methI_aggregate$V1,":",DMAR_methI_aggregate$V2,"-",DMAR_methI_aggregate$V3,sep = "")


DMAR <-Reduce(function(x, y) merge(x, y,by = "chrompos", all.x = T), list(DMAR_DMR[,c("chrompos","V1","V2","V3","V8","V9","V10","V11","V12")],DMAR_DAR[,c("chrompos","V8","V9","V10","V11","V12")],DMAR_IDR15[,c("chrompos","V4")],DMAR_IDR24[,c("chrompos","V4")],DMAR_IDRM[,c("chrompos","V4")],DMAR_IDRI[,c("chrompos","V4")],DMAR_meth15_aggregate[,c("chrompos","AveMeth")],DMAR_meth24_aggregate[,c("chrompos","AveMeth")],DMAR_methM_aggregate[,c("chrompos","AveMeth")],DMAR_methI_aggregate[,c("chrompos","AveMeth")],DMAR_meth15_aggregate[,c("chrompos","CpGcount")],DMAR_meth24_aggregate[,c("chrompos","CpGcount")],DMAR_methM_aggregate[,c("chrompos","CpGcount")],DMAR_methI_aggregate[,c("chrompos","CpGcount")],DMAR_meth15_aggregate[,c("chrompos","AveCov")],DMAR_meth24_aggregate[,c("chrompos","AveCov")],DMAR_methM_aggregate[,c("chrompos","AveCov")],DMAR_methI_aggregate[,c("chrompos","AveCov")]))
colnames(DMAR) <- c("chrompos","chr","start","end","DMRsize","DMRs15vs24","DMRs24vMel","DMRs24vIri","DMRMelvIri","DARsize","DARs15vs24","DARs24vMel","DARs24vIri","DARMelvIri","IDR_s15","IDR_s24","IDR_M","IDR_I","Meth_s15","Meth_s24","Meth_Mel","Meth_Iri","CpG_s15","CpG_s24","CpG_Mel","CpG_Iri","AveCovCpG_s15","AveCovCpG_s24","AveCovCpG_Mel","AveCovCpG_Iri")
DMAR$DMARsize <- DMAR$end-DMAR$start
DMAR[is.na(DMAR)] <- 0
DMAR[DMAR == "."] <- 0
DMAR[DMAR == "-1"] <- 0
DMAR$DMRsize <- as.numeric(DMAR$DMRsize)
DMAR$DARsize <- as.numeric(DMAR$DARsize)
DMAR$DMRs15vs24 <- as.numeric(DMAR$DMRs15vs24 )
DMAR$DMRs24vMel <- as.numeric(DMAR$DMRs24vMel )
DMAR$DMRs24vIri <- as.numeric(DMAR$DMRs24vIri )
DMAR$DMRMelvIri <- as.numeric(DMAR$DMRMelvIri)
DMAR$DARs15vs24 <- as.numeric(DMAR$DARs15vs24 )
DMAR$DARs24vMel <- as.numeric(DMAR$DARs24vMel )
DMAR$DARs24vIri <- as.numeric(DMAR$DARs24vIri )
DMAR$DARMelvIri <- as.numeric(DMAR$DARMelvIri)

DMAR <- DMAR[!duplicated(DMAR$chrompos),]

DMAR <- DMAR[,c(2,3,4,1,31,5:30)]
DMAR$aveCpGcount <- apply(DMAR,1,function(x) {(as.numeric(x[24])+as.numeric(x[25])+as.numeric(x[26])+as.numeric(x[27]))/4})
DMAR$CpGdensity100bp <- apply(DMAR,1,function(x) {(as.numeric(x[32])/as.numeric(x[5]))*100})
DMAR$DARMelvIri <- DMAR$DARMelvIri*-1   ######### the directionality of DAR is switch for Mel v Iri for DARs!!!!###################


# export to a final bed file including all the DMAR information
write.table(DMAR, "All_DMAR_Combined_wINFO_NEW_042620.bed", row.names = F, col.names = T, sep = "\t",quote =F)


# Calculate the counts of differentially methylated and accessible regions
str(DMAR[DMAR$DMRMelvIri >0 & DMAR$DARMelvIri >0,]) #1697 hypo opening
str(DMAR[DMAR$DMRMelvIri >0 & DMAR$DARMelvIri <0,]) #2 hypo closing
str(DMAR[DMAR$DMRMelvIri <0 & DMAR$DARMelvIri >0,]) #3 hyper opening
str(DMAR[DMAR$DMRMelvIri <0 & DMAR$DARMelvIri <0,]) #838 hyper closing
str(DMAR[DMAR$DMRMelvIri <0 & DMAR$DARMelvIri ==0,]) #1422 solo hyper DMR
str(DMAR[DMAR$DMRMelvIri >0 & DMAR$DARMelvIri ==0,]) # 1918 solo hypoDMR
str(DMAR[DMAR$DMRMelvIri ==0 & DMAR$DARMelvIri >0,]) #6239 solo opening
str(DMAR[DMAR$DMRMelvIri ==0 & DMAR$DARMelvIri <0,]) #5197 solo closing


str(DMAR[DMAR$DMRs15vs24 > 0 & DMAR$DARs15vs24 == 0,]) #352 hypoDMR
str(DMAR[DMAR$DMRs15vs24 < 0 & DMAR$DARs15vs24 == 0,]) #158 hyperDMR
str(DMAR[DMAR$DMRs15vs24 == 0 & DMAR$DARs15vs24 > 0,]) #877 solo openDAR
str(DMAR[DMAR$DMRs15vs24 == 0 & DMAR$DARs15vs24 < 0,]) #8 solo closeDAR
str(DMAR[DMAR$DMRs15vs24 > 0 & DMAR$DARs15vs24 < 0,]) # 0 hypoclosing
str(DMAR[DMAR$DMRs15vs24 > 0 & DMAR$DARs15vs24 > 0,]) # 8 hypo opening

str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri == 0,]) #4378 solo shared hypoDMR
str(DMAR[DMAR$DMRs24vMel < 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri < 0 & DMAR$DARs24vIri == 0,]) #151 solo shared hyperDMR
str(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel > 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri > 0,]) #589 solo shared openDAR
str(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel < 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri < 0,]) #15116 solo shared closeDAR
str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel < 0 & DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri < 0,]) # 423 shared hypoclosing
str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel > 0 & DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri > 0,]) #1286 shared hypo opening

str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri <=0 & DMAR$DARs24vIri == 0,]) #7319 solo Mel hypoDMR 
str(DMAR[DMAR$DMRs24vMel < 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri >= 0 & DMAR$DARs24vIri == 0,]) #210 solo Mel hyperDMR
str(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel > 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri <= 0,]) #3617 solo Mel openDAR
str(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel < 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri >= 0,]) #8819 solo Mel closeDAR
str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel < 0 & DMAR$DMRs24vIri <= 0 & DMAR$DARs24vIri >= 0,]) # 140 Melhypoclosing
str(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel > 0 & DMAR$DMRs24vIri <= 0 & DMAR$DARs24vIri <= 0,]) #2446 Mel hypo opening

str(DMAR[DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri == 0 & DMAR$DMRs24vMel <=0 & DMAR$DARs24vMel == 0,]) #6319 solo Iri hypoDMR 
str(DMAR[DMAR$DMRs24vIri < 0 & DMAR$DARs24vIri == 0 & DMAR$DMRs24vMel >= 0 & DMAR$DARs24vMel == 0,]) #289 solo Iri hyperDMR
str(DMAR[DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri > 0 & DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel <= 0,]) #3657 solo Iri openDAR
str(DMAR[DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri < 0 & DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel >= 0,]) #6052 solo Iri closeDAR
str(DMAR[DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri < 0 & DMAR$DMRs24vMel <= 0 & DMAR$DARs24vMel >= 0,]) # 62 Iri hypoclosing
str(DMAR[DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri > 0 & DMAR$DMRs24vMel <= 0 & DMAR$DARs24vMel <= 0,]) #3095 Iri hypo opening


#manual order of rows in heatmap
row1 <- DMAR[DMAR$DARs15vs24>0 & DMAR$Meth_s15 >50 & DMAR$Meth_s24 >50 & DMAR$aveCpGcount > 0 & DMAR$Meth_Mel >50 &  DMAR$Meth_Iri >50,]
row1$combineMI <- (row1$Meth_Mel+ row1$Meth_I)/2
row1 <- row1[order(-row1$combineMI),]

row2 <- DMAR[DMAR$DARs15vs24>0 & DMAR$Meth_s15 >50 & DMAR$Meth_s24 >50 & DMAR$aveCpGcount > 0 & DMAR$Meth_Mel >50 &  DMAR$Meth_Iri <50,]
row2 <- row2[order(-row2$Meth_Iri),]

row3 <- DMAR[DMAR$DARs15vs24>0 & DMAR$Meth_s15 >50 & DMAR$Meth_s24 >50 & DMAR$aveCpGcount > 0 & DMAR$Meth_Mel <50 &  DMAR$Meth_Iri >50,]
row3 <- row3[order(-row3$Meth_Mel),]

row4 <- DMAR[DMAR$DARs15vs24>0 & DMAR$Meth_s15 >50 & DMAR$Meth_s24 >50 & DMAR$aveCpGcount > 0 & DMAR$Meth_Mel <50 &  DMAR$Meth_Iri <50,]
row4$combineMI <- (row4$Meth_Mel+ row4$Meth_I)/2
row4 <- row4[order(-row4$combineMI),]

row5 <- DMAR[DMAR$DARs15vs24>0 & DMAR$Meth_s15 <50 & DMAR$Meth_s24 <50 & DMAR$aveCpGcount > 0 & DMAR$Meth_Mel <50 &  DMAR$Meth_Iri <50,]
row5$combineMI <- (row5$Meth_Mel+ row5$Meth_I)/2
row5 <- row5[order(-row5$combineMI),]


roworder<- c(row1$chrompos,row2$chrompos,row3$chrompos,row4$chrompos,row5$chrompos)
rowadd<-setdiff(DMAR[DMAR$DARs15vs24>0 & DMAR$aveCpGcount > 0,]$chrompos,roworder)
roworder2 <- c(roworder,rowadd)

name <- All[,4]
df.OG2 <- All[,c("Meth_s15","Meth_s24","Meth_Mel","Meth_Iri")]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("15somite", "24hpf","Mel", "Iri")

ht4 = Heatmap(df.OG2[roworder2,], column_title = "Methylation of 15somite to 24hpf DARs",name= "Methylaton",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), cluster_rows=FALSE, cluster_columns = FALSE,show_row_names = FALSE)
ht4




### MAKE A BAR CHART WITH NUMBER OF DMRs DARs and DMARs####
library(RColorBrewer)
mypalette <- brewer.pal(12,"Paired")

DMAR_counts <- read.table("DMR_DAR_DMAR_countlist_for_figure.txt", sep = "\t", header = T, quote = "",stringsAsFactors = F)

dmrplot <- DMAR_counts[DMAR_counts$Type == "DMR",]
dmrplot$Cell <- factor(dmrplot$Cell, levels = c("s15vs24","shared","s24vMel","s24vIri","MelvIri"))
p <- ggplot(data=dmrplot, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DMR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DMRs")+scale_fill_manual(values = mypalette[c(1,2)])+scale_color_manual(values=mypalette[c(1,2)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=dmrplot,position = position_dodge(width =0.9),aes(x=Cell,y=Count+300,label=Count),fontface='bold',hjust=0.5,size=5)
p


darplot <- DMAR_counts[DMAR_counts$Type == "DAR",]
darplot$Cell <- factor(darplot$Cell, levels = c("s15vs24","shared","s24vMel","s24vIri","MelvIri"))
p <- ggplot(data=darplot, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DAR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DARs")+scale_fill_manual(values = mypalette[c(3,4)])+scale_color_manual(values=mypalette[c(3,4)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=darplot,position = position_dodge(width =0.9),aes(x=Cell,y=Count+500,label=Count),fontface='bold',hjust=0.5,size=5)
p


degplot <- DMAR_counts[DMAR_counts$Type == "DEG",]
degplot$Cell <- factor(degplot$Cell, levels = c("s15vs24","shared","s24vMel","s24vIri","MelvIri"))
p <- ggplot(data=degplot, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DEG count (filtering genes > 5 RPKM)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DEGs")+scale_fill_manual(values = mypalette[c(5,6)])+scale_color_manual(values=mypalette[c(5,6)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=degplot,position = position_dodge(width =0.9),aes(x=Cell,y=Count+50,label=Count),fontface='bold',hjust=0.5,size=5)
p


dmarplot <- DMAR_counts[DMAR_counts$Type == "DMAR" | DMAR_counts$Type == "DMAR2",]
dmarplot <- dmarplot[dmarplot$Name != "hypercloseDMAR" & dmarplot$Name != "hyperopenDMAR",]
dmarplot$Cell <- factor(dmarplot$Cell, levels = c("s15vs24","shared","s24vMel","s24vIri","MelvIri"))
dmarplot$Name <- factor(dmarplot$Name, levels = c("solohyperDMR","solohypoDMR","solocloseDAR","soloopenDAR","hypocloseDMAR","hypoopenDMAR"))

p <- ggplot(data=dmarplot, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DMAR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DMARs")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_y_continuous(lim = c(0,16000))+scale_color_manual(values=mypalette[c(1,2,3,4,9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=dmarplot,position = position_dodge(width =0.9),aes(x=Cell,y=Count+400,label=Count),fontface='bold',hjust=0.5,size=5)
p


dmarplot2 <- DMAR_counts[DMAR_counts$Type == "DMAR" | DMAR_counts$Type == "DMAR2",]
dmarplot2 <- dmarplot2[dmarplot2$Cell == "MelvIri",]
dmarplot2$Name <- factor(dmarplot2$Name, levels = c("solohyperDMR","solohypoDMR","solocloseDAR","soloopenDAR","hypocloseDMAR","hypoopenDMAR","hypercloseDMAR","hyperopenDMAR"))

p <- ggplot(data=dmarplot2, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DMAR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DMARs")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10,7,8)])+scale_y_continuous(lim = c(0,16000))+scale_color_manual(values=mypalette[c(1,2,3,4,9,10,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=dmarplot2,position = position_dodge(width =0.9),aes(x=Cell,y=Count+400,label=Count),fontface='bold',hjust=0.5,size=5)
p
```
```{bash}
###################################### {Bash script} ######################################
####### FOR methylC track on new browser ######
chr1    10542   10543   CG      0.923   -       26 <--- this is the format. combine all CpG strands so all will be +

head Reformatted_s15pos_CombinedRep_cpg_combined_DSS.txt
chr1  81    82    39    39    1
chr1  92    93    37    30    0.810811
chr1  149   150   20    20    1
chr1  234   235   23    18    0.782609
chr1  317   318   17    14    0.823529
chr1  357   358   10    1     0.1
chr1  374   375   5     4     0.8
chr1  609   610   32    29    0.90625
chr1  716   717   15    15    1
chr1  735   736   18    0     0

for i in Reformatted*s15pos*DSS.txt; do awk -v OFS="\t" '{print $1,$2,$3,"CG",$6,"+",$4}' $i > ${i/.txt/_forNewBrowser_methylCtrack.bed};done
for i in Reformatted*Mel*DSS.txt; do awk -v OFS="\t" '{print $1,$2,$3,"CG",$6,"+",$4}' $i > ${i/.txt/_forNewBrowser_methylCtrack.bed};done
for i in Reformatted*Iri*DSS.txt; do awk -v OFS="\t" '{print $1,$2,$3,"CG",$6,"+",$4}' $i > ${i/.txt/_forNewBrowser_methylCtrack.bed};done
for i in Reformatted*s24pos*DSS.txt; do awk -v OFS="\t" '{print $1,$2,$3,"CG",$6,"+",$4}' $i > ${i/.txt/_forNewBrowser_methylCtrack.bed};done

for i in *_forNewBrowser_methylCtrack.bed; do echo $i && sort -k1,1 -k2,2n $i > ${i/bed/sorted.bed};done
for i in *_forNewBrowser_methylCtrack.sorted.bed; do echo $i && bgzip $i;done
for i in *_forNewBrowser_methylCtrack.sorted.bed.gz; do echo $i && tabix -p bed $i;done

cd /bar/jjang/public_html/pigment/pigment_BS
ln -s /scratch/jjang/PIGMENT_PROJECT/Pigment_BS_ALL/DSS/*_forNewBrowser_methylCtrack.sorted.bed.gz .
ln -s /scratch/jjang/PIGMENT_PROJECT/Pigment_BS_ALL/DSS/*_forNewBrowser_methylCtrack.sorted.bed.gz.tbi .


## Use https://biit.cs.ut.ee/gprofiler/orth for gene conversion for AME ##
Some gene names still missing. Manually curate using http://www.zebrafishmine.org/template.do?name=Human_Gene_Curated_Zeb_Ortho. Use https://resources.altius.org/~jvierstra/projects/motif-clustering/ to identify cluster and DBD
MASTER-LIST-TF-orthologs.csv
MASTER-LIST-TF-orthologs_converted_ensemble_withMotifClusterInfo.csv

######### Use AME to call motif enrichment #########
#make background file of all peaks#
awk -v OFS="\t" '{print $1, $2, $3, $4}' All_IDR_peaks_Merged.txt > All_IDR_peaks_Merged_forfasta.txt
twoBitToFa -bed=All_IDR_peaks_Merged_forfasta.txt danRer10.2bit All_IDR_peaks_Merged.fa

#make fasta files of DARs
#hypoDMRs#
awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_hypoDMR_d30_p0.01.bed > Iri_specific_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Iri_specific_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Iri_specific_hypoDMR_d30_p0.01.fasta

ame --verbose 1 --oc Iri_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_hypoDMR_d30_p0.01.bed > Mel_specific_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_specific_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Mel_specific_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_hypoDMR -scoring totalhits --method fisher --control /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/All_IDR_peaks_Merged.fa Mel_specific_hypoDMR_d30_p0.01.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_hypoDMR_d30_p0.01.bed > Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Mel_Iri_shared_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_Iri_shared_hypoDMR -scoring totalhits --method fisher --control /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/All_IDR_peaks_Merged.fa Mel_Iri_shared_hypoDMR_d30_p0.01.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_hypoDMR_d30_p0.01.bed > Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Mel_Iri_shared_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_Iri_shared_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme

#DARs#
awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_opening_DAR.bed > Iri_specific_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_opening_DAR.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Iri_specific_opening_DAR.fasta
ame --verbose 1 --oc Iri_specific_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_opening_DAR.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_closing_DAR.bed > Iri_specific_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_closing_DAR.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Iri_specific_closing_DAR.fasta
ame --verbose 1 --oc Iri_specific_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_closing_DAR.fasta JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_opening_DAR.bed > Mel_specific_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_opening_DAR.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Mel_specific_opening_DAR.fasta
ame --verbose 1 --oc Mel_specific_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_opening_DAR.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_closing_DAR.bed > Mel_specific_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_closing_DAR.onlyCoord.bed /bar/genomes/danRer10/danRer10.2bit Mel_specific_closing_DAR.fasta
ame --verbose 1 --oc Mel_specific_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_closing_DAR.fasta JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_opening_DAR.bed > Mel_Iri_shared_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_opening_DAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_opening_DAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_opening_DAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_closing_DAR.bed > Mel_Iri_shared_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_closing_DAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_closing_DAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_closing_DAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme


#DMARs#
awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Iri_specific_hypo_opening.DMAR.bed > Iri_specific_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Iri_specific_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Iri_specific_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Mel_specific_hypo_opening.DMAR.bed > Mel_specific_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Mel_specific_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Mel_specific_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Mel_Iri_shared_hypo_opening.DMAR.bed > Mel_Iri_shared_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme
```
```{R}
###################################### {R script} ######################################
#Load AME findings on R
masterTF_orth_list <- read.csv("MASTER-LIST-TF-orthologs_converted_ensemble_withMotifClusterInfo.csv", header = T, sep = ",", stringsAsFactors = F)
colnames(masterTF_orth_list) <- c("motif_alt_ID","gene","genename","Cluster","DBD")

Iri_hypoDMR_AME <- read.table("Iri_hypoDMR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_hypoDMR_AME <- read.table("Mel_hypoDMR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Shared_hypoDMR_AME <- read.table("Mel_Iri_shared_hypoDMR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Iri_openDAR_AME <- read.table("Iri_specific_openDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_openDAR_AME <- read.table("Mel_specific_openDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Iri_closeDAR_AME <- read.table("Iri_specific_closeDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_closeDAR_AME <- read.table("Mel_specific_closeDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_Iri_shared_openDAR_AME <- read.table("Mel_Iri_shared_openDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_Iri_shared_closeDAR_AME <- read.table("Mel_Iri_shared_closeDAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Iri_openDMAR_AME <- read.table("Iri_specific_openDMAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_openDMAR_AME <- read.table("Mel_specific_openDMAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")
Mel_Iri_shared_openDMAR_AME <- read.table("Mel_Iri_shared_openDMAR_ame.tsv", header = T, stringsAsFactors = F, sep = "\t")

# Change ID to upper case
Iri_hypoDMR_AME$motif_alt_ID <- toupper(Iri_hypoDMR_AME$motif_alt_ID)
Mel_hypoDMR_AME$motif_alt_ID <- toupper(Mel_hypoDMR_AME$motif_alt_ID)
Shared_hypoDMR_AME$motif_alt_ID <- toupper(Shared_hypoDMR_AME$motif_alt_ID)
Iri_openDAR_AME$motif_alt_ID <- toupper(Iri_openDAR_AME$motif_alt_ID)
Mel_openDAR_AME$motif_alt_ID <- toupper(Mel_openDAR_AME$motif_alt_ID)
Mel_Iri_shared_openDAR_AME$motif_alt_ID <- toupper(Mel_Iri_shared_openDAR_AME$motif_alt_ID)
Iri_closeDAR_AME$motif_alt_ID <- toupper(Iri_closeDAR_AME$motif_alt_ID)
Mel_closeDAR_AME$motif_alt_ID <- toupper(Mel_closeDAR_AME$motif_alt_ID)
Mel_Iri_shared_closeDAR_AME$motif_alt_ID <- toupper(Mel_Iri_shared_closeDAR_AME$motif_alt_ID)
Iri_openDMAR_AME$motif_alt_ID <- toupper(Iri_openDMAR_AME$motif_alt_ID)
Mel_openDMAR_AME$motif_alt_ID <- toupper(Mel_openDMAR_AME$motif_alt_ID)
Mel_Iri_shared_openDMAR_AME$motif_alt_ID <- toupper(Mel_Iri_shared_openDMAR_AME$motif_alt_ID)

# Merge masterTF_orth_list with pigment specific/shared motifs by motif alt ID
Iri_hypoDMR_AME_genename <- merge(masterTF_orth_list, Iri_hypoDMR_AME, by = "motif_alt_ID")
Mel_hypoDMR_AME_genename <- merge(masterTF_orth_list, Mel_hypoDMR_AME, by = "motif_alt_ID")
Shared_hypoDMR_AME_genename <- merge(masterTF_orth_list, Shared_hypoDMR_AME, by = "motif_alt_ID")

Iri_openDAR_AME_genename <- merge(masterTF_orth_list, Iri_openDAR_AME, by = "motif_alt_ID")
Mel_openDAR_AME_genename <- merge(masterTF_orth_list, Mel_openDAR_AME, by = "motif_alt_ID")
Shared_openDAR_AME_genename <- merge(masterTF_orth_list, Mel_Iri_shared_openDAR_AME, by = "motif_alt_ID")
Iri_closeDAR_AME_genename <- merge(masterTF_orth_list, Iri_closeDAR_AME, by = "motif_alt_ID")
Mel_closeDAR_AME_genename <- merge(masterTF_orth_list, Mel_closeDAR_AME, by = "motif_alt_ID")
Shared_closeDAR_AME_genename <- merge(masterTF_orth_list, Mel_Iri_shared_closeDAR_AME, by = "motif_alt_ID")

Iri_openDMAR_AME_genename <- merge(masterTF_orth_list, Iri_openDMAR_AME, by = "motif_alt_ID")
Mel_openDMAR_AME_genename <- merge(masterTF_orth_list, Mel_openDMAR_AME, by = "motif_alt_ID")
Shared_openDMAR_AME_genename <- merge(masterTF_orth_list, Mel_Iri_shared_openDMAR_AME, by = "motif_alt_ID")


# Merge DEG information with pigment specific/shared DM/ARs enriched motifs 
Iri_hypoDMR_AME_genename_DEG <- merge(Iri_hypoDMR_AME_genename,DEG2_Z, by = "gene")
Mel_hypoDMR_AME_genename_DEG <- merge(Mel_hypoDMR_AME_genename,DEG2_Z, by = "gene")
Shared_hypoDMR_AME_genename_DEG <- merge(Shared_hypoDMR_AME_genename,DEG2_Z, by = "gene")
Iri_closeDAR_AME_genename_DEG <- merge(Iri_closeDAR_AME_genename,DEG2_Z, by = "gene")
Mel_closeDAR_AME_genename_DEG <- merge(Mel_closeDAR_AME_genename,DEG2_Z, by = "gene")
Shared_closeDAR_AME_genename_DEG <- merge(Shared_closeDAR_AME_genename,DEG2_Z, by = "gene")
Iri_openDAR_AME_genename_DEG <- merge(Iri_openDAR_AME_genename,DEG2_Z, by = "gene")
Mel_openDAR_AME_genename_DEG <- merge(Mel_openDAR_AME_genename,DEG2_Z, by = "gene")
Shared_openDAR_AME_genename_DEG <- merge(Shared_openDAR_AME_genename,DEG2_Z, by = "gene")
Iri_openDMAR_AME_genename_DEG <- merge(Iri_openDMAR_AME_genename,DEG2_Z, by = "gene")
Mel_openDMAR_AME_genename_DEG <- merge(Mel_openDMAR_AME_genename,DEG2_Z, by = "gene")
Shared_openDMAR_AME_genename_DEG <- merge(Shared_openDMAR_AME_genename,DEG2_Z, by = "gene")


Iri_hypoDMR_AME_genename_DEG$pvalue_log <- -log10(Iri_hypoDMR_AME_genename_DEG$p.value)
Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500 
Mel_hypoDMR_AME_genename_DEG$pvalue_log <- -log10(Mel_hypoDMR_AME_genename_DEG$p.value)
Mel_hypoDMR_AME_genename_DEG[Mel_hypoDMR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Shared_hypoDMR_AME_genename_DEG$pvalue_log <- -log10(Shared_hypoDMR_AME_genename_DEG$p.value)
Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500

Iri_openDAR_AME_genename_DEG$pvalue_log <- -log10(Iri_openDAR_AME_genename_DEG$p.value)
Iri_openDAR_AME_genename_DEG[Iri_openDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Mel_openDAR_AME_genename_DEG$pvalue_log <- -log10(Mel_openDAR_AME_genename_DEG$p.value)
Mel_openDAR_AME_genename_DEG[Mel_openDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Shared_openDAR_AME_genename_DEG$pvalue_log <- -log10(Shared_openDAR_AME_genename_DEG$p.value)
Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500

Iri_closeDAR_AME_genename_DEG$pvalue_log <- -log10(Iri_closeDAR_AME_genename_DEG$p.value)
Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Mel_closeDAR_AME_genename_DEG$pvalue_log <- -log10(Mel_closeDAR_AME_genename_DEG$p.value)
Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Shared_closeDAR_AME_genename_DEG$pvalue_log <- -log10(Shared_closeDAR_AME_genename_DEG$p.value)
Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500

Iri_openDMAR_AME_genename_DEG$pvalue_log <- -log10(Iri_openDMAR_AME_genename_DEG$p.value)
Iri_openDMAR_AME_genename_DEG[Iri_openDMAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Mel_openDMAR_AME_genename_DEG$pvalue_log <- -log10(Mel_openDMAR_AME_genename_DEG$p.value)
Mel_openDMAR_AME_genename_DEG[Mel_openDMAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500
Shared_openDMAR_AME_genename_DEG$pvalue_log <- -log10(Shared_openDMAR_AME_genename_DEG$p.value)
Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$pvalue_log == Inf,]$pvalue_log <- 500




#Iri open AME motif figure
Iri_hypoDMR_AME_filtered <-Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$s24vIri_log2_change <0 | Iri_hypoDMR_AME_genename_DEG$MelvIri_log2_change <0 & Iri_hypoDMR_AME_genename_DEG$Iri > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in iri 

Iri_AME_masterlist <- unique(rbind(Iri_hypoDMR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Iri_openDAR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Iri_openDMAR_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

Iri_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(Iri_AME_masterlist,Iri_hypoDMR_AME_filtered,Iri_openDAR_AME_filtered,Iri_openDMAR_AME_filtered))
Iri_AME_masterlist_withLogP[is.na(Iri_AME_masterlist_withLogP)] <- 0

Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),]

Iri_AME_masterlist_withLogP_all3 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x >0 & Iri_AME_masterlist_withLogP$pvalue_log.y >0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]
Iri_AME_masterlist_withLogP_all2 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x ==0 & Iri_AME_masterlist_withLogP$pvalue_log.y >0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]
Iri_AME_masterlist_withLogP_all1 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x ==0 & Iri_AME_masterlist_withLogP$pvalue_log.y ==0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]

Iri_AME_masterlist_withLogP_reordered <- rbind(Iri_AME_masterlist_withLogP_all3,Iri_AME_masterlist_withLogP_all2,Iri_AME_masterlist_withLogP_all1)
Iri_AME_masterlist_withLogP_reordered<-unique(Iri_AME_masterlist_withLogP_reordered)


library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- Iri_AME_masterlist_withLogP_reordered[,1]
df.OG2 <- -Iri_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- Iri_AME_masterlist_withLogP_reordered[,c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- Iri_AME_masterlist_withLogP_reordered[,c(8,9,10)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #12 x 3

# Order based on DBD
Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),]
name <- Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),1]
df.OG2 <- -Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),c(5,6,7)]
row.names(df.OG) <- name
#kclus <- kmeans(df.OG,10)
#split <-  kclus$cluster
#ht =Heatmap(df.OG, column_title = "TFs",name= "RPKM",col = colorRamp2(c(0,1,5,10,50), c("#e7eff4", "#c3f7dd","#c1fa9e","#fcda79","#ff526a")), 
#    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE,split = split)
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),c(9,10,11)]
row.names(df.OG) <- name
#kclus <- kmeans(df.OG,10)
#split <-  kclus$cluster
#ht =Heatmap(df.OG, column_title = "TFs",name= "RPKM",col = colorRamp2(c(0,1,5,10,50), c("#e7eff4", "#c3f7dd","#c1fa9e","#fcda79","#ff526a")), 
#    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE,split = split)
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #12 x 3



# Mel open AME motif figure
Mel_hypoDMR_AME_filtered <- Mel_hypoDMR_AME_genename_DEG[Mel_hypoDMR_AME_genename_DEG$s24vMel_log2_change <0 | Mel_hypoDMR_AME_genename_DEG$MelvIri_log2_change >0 & Mel_hypoDMR_AME_genename_DEG$Mel > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in mel
Mel_openDAR_AME_filtered <- Mel_openDAR_AME_genename_DEG[Mel_openDAR_AME_genename_DEG$s24vMel_log2_change <0 | Mel_openDAR_AME_genename_DEG$MelvIri_log2_change >0 & Mel_openDAR_AME_genename_DEG$Mel > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"),] #upreg in mel
Mel_openDMAR_AME_filtered <- Mel_openDMAR_AME_genename_DEG[Mel_openDMAR_AME_genename_DEG$s24vMel_log2_change <0 | Mel_openDMAR_AME_genename_DEG$MelvIri_log2_change >0 & Mel_openDMAR_AME_genename_DEG$Mel > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"),] #upreg in mel
Mel_AME_masterlist <- unique(rbind(Mel_hypoDMR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_openDAR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_openDMAR_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

Mel_AME_masterlist <- unique(rbind(Mel_hypoDMR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_openDAR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_openDMAR_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

Mel_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(Mel_AME_masterlist,Mel_hypoDMR_AME_filtered,Mel_openDAR_AME_filtered,Mel_openDMAR_AME_filtered))
Mel_AME_masterlist_withLogP[is.na(Mel_AME_masterlist_withLogP)] <- 0
Mel_AME_masterlist_withLogP <- unique(Mel_AME_masterlist_withLogP)

Mel_AME_masterlist_withLogP_all3 <- Mel_AME_masterlist_withLogP[Mel_AME_masterlist_withLogP$pvalue_log.x >0 & Mel_AME_masterlist_withLogP$pvalue_log.y >0 & Mel_AME_masterlist_withLogP$pvalue_log >0,]
Mel_AME_masterlist_withLogP_all2 <- Mel_AME_masterlist_withLogP[Mel_AME_masterlist_withLogP$pvalue_log.x ==0 & Mel_AME_masterlist_withLogP$pvalue_log.y >0 & Mel_AME_masterlist_withLogP$pvalue_log >0,]
Mel_AME_masterlist_withLogP_all1 <- Mel_AME_masterlist_withLogP[Mel_AME_masterlist_withLogP$pvalue_log.x ==0 & Mel_AME_masterlist_withLogP$pvalue_log.y ==0 & Mel_AME_masterlist_withLogP$pvalue_log >0,]

Mel_AME_masterlist_withLogP_reordered <- rbind(Mel_AME_masterlist_withLogP_all3,Mel_AME_masterlist_withLogP_all2,Mel_AME_masterlist_withLogP_all1)
Mel_AME_masterlist_withLogP_reordered <-unique(Mel_AME_masterlist_withLogP_reordered)


library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- Mel_AME_masterlist_withLogP_reordered[,1]
df.OG2 <- -Mel_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- Mel_AME_masterlist_withLogP_reordered[,c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- Mel_AME_masterlist_withLogP_reordered[,c(8,9,10)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #12 x 3.5

# Order based on DBD
name <- Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),1]
df.OG2 <- -Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),c(9,10,11)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #12 x 3.5



#Shared AME opening
shared_hypoDMR_up_AME_filtered <-Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$s24vIri_log2_change <0 & Shared_hypoDMR_AME_genename_DEG$s24vMel_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in both
shared_openDAR_up_AME_filtered <-Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$s24vIri_log2_change <0 & Shared_openDAR_AME_genename_DEG$s24vMel_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] # upreg in both
shared_openDMAR_up_AME_filtered <-Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$s24vIri_log2_change <0 & Shared_openDMAR_AME_genename_DEG$s24vMel_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] # upreg in both

shared_hypoDMR_up_AME_masterlist <- unique(rbind(shared_hypoDMR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)],shared_openDAR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)],shared_openDMAR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

shared_hypoDMR_up_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(shared_hypoDMR_up_AME_masterlist,shared_hypoDMR_up_AME_filtered,shared_openDAR_up_AME_filtered,shared_openDMAR_up_AME_filtered))
shared_hypoDMR_up_AME_masterlist_withLogP[is.na(shared_hypoDMR_up_AME_masterlist_withLogP)] <- 0
shared_hypoDMR_up_AME_masterlist_withLogP<- unique(shared_hypoDMR_up_AME_masterlist_withLogP)

shared_hypoDMR_up_AME_masterlist_withLogP_all3 <- shared_hypoDMR_up_AME_masterlist_withLogP[shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.x >0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.y >0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log >0,]
shared_hypoDMR_up_AME_masterlist_withLogP_all2 <- shared_hypoDMR_up_AME_masterlist_withLogP[shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.x ==0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.y >0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log >0,]
shared_hypoDMR_up_AME_masterlist_withLogP_all1 <- shared_hypoDMR_up_AME_masterlist_withLogP[shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.x ==0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log.y ==0 & shared_hypoDMR_up_AME_masterlist_withLogP$pvalue_log >0,]

shared_hypoDMR_up_AME_masterlist_withLogP_reordered <- rbind(shared_hypoDMR_up_AME_masterlist_withLogP_all3,shared_hypoDMR_up_AME_masterlist_withLogP_all2,shared_hypoDMR_up_AME_masterlist_withLogP_all1)
shared_hypoDMR_up_AME_masterlist_withLogP_reordered<-unique(shared_hypoDMR_up_AME_masterlist_withLogP_reordered)


library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- shared_hypoDMR_up_AME_masterlist_withLogP_reordered[,1]
df.OG2 <- -shared_hypoDMR_up_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- shared_hypoDMR_up_AME_masterlist_withLogP_reordered[,c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- shared_hypoDMR_up_AME_masterlist_withLogP_reordered[,c(8,9,10)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #6 x 3

# Order based on DBD
shared_hypoDMR_up_AME_masterlist_withLogP[order(shared_hypoDMR_up_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_up_AME_masterlist_withLogP[,"genename.x"]),]
name <- shared_hypoDMR_up_AME_masterlist_withLogP[order(shared_hypoDMR_up_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_up_AME_masterlist_withLogP[,"genename.x"]),1]
df.OG2 <- -shared_hypoDMR_up_AME_masterlist_withLogP[order(shared_hypoDMR_up_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_up_AME_masterlist_withLogP[,"genename.x"]),c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- shared_hypoDMR_up_AME_masterlist_withLogP[order(shared_hypoDMR_up_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_up_AME_masterlist_withLogP[,"genename.x"]),c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- shared_hypoDMR_up_AME_masterlist_withLogP[order(shared_hypoDMR_up_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_up_AME_masterlist_withLogP[,"genename.x"]),c(9,10,11)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #6 x 3


#shared open downregulated
shared_hypoDMR_down_AME_filtered <-Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$s24vIri_log2_change >0 & Shared_hypoDMR_AME_genename_DEG$s24vMel_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #downreg in both
shared_openDAR_down_AME_filtered <-Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$s24vIri_log2_change >0 & Shared_openDAR_AME_genename_DEG$s24vMel_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] # downreg in both
shared_openDMAR_down_AME_filtered <-Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$s24vIri_log2_change >0 & Shared_openDMAR_AME_genename_DEG$s24vMel_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] # downreg in both

shared_hypoDMR_down_AME_masterlist <- unique(rbind(shared_hypoDMR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)],shared_openDAR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)],shared_openDMAR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

shared_hypoDMR_down_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(shared_hypoDMR_down_AME_masterlist,shared_hypoDMR_down_AME_filtered,shared_openDAR_down_AME_filtered,shared_openDMAR_down_AME_filtered))
shared_hypoDMR_down_AME_masterlist_withLogP[is.na(shared_hypoDMR_down_AME_masterlist_withLogP)] <- 0
shared_hypoDMR_down_AME_masterlist_withLogP<-unique(shared_hypoDMR_down_AME_masterlist_withLogP)

shared_hypoDMR_down_AME_masterlist_withLogP_all3 <- shared_hypoDMR_down_AME_masterlist_withLogP[shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.x >0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.y >0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log >0,]
shared_hypoDMR_down_AME_masterlist_withLogP_all2a <- shared_hypoDMR_down_AME_masterlist_withLogP[shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.x >0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.y ==0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log >0,]
shared_hypoDMR_down_AME_masterlist_withLogP_all2b <- shared_hypoDMR_down_AME_masterlist_withLogP[shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.x ==0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.y >0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log >0,]
shared_hypoDMR_down_AME_masterlist_withLogP_all1 <- shared_hypoDMR_down_AME_masterlist_withLogP[shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.x >0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log.y ==0 & shared_hypoDMR_down_AME_masterlist_withLogP$pvalue_log ==0,]

shared_hypoDMR_down_AME_masterlist_withLogP_reordered <- rbind(shared_hypoDMR_down_AME_masterlist_withLogP_all3,shared_hypoDMR_down_AME_masterlist_withLogP_all2a,shared_hypoDMR_down_AME_masterlist_withLogP_all2b,shared_hypoDMR_down_AME_masterlist_withLogP_all1)
shared_hypoDMR_down_AME_masterlist_withLogP_reordered <- unique(shared_hypoDMR_down_AME_masterlist_withLogP_reordered)

library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- make.unique(shared_hypoDMR_down_AME_masterlist_withLogP_reordered[,1],sep=".")
df.OG2 <- -shared_hypoDMR_down_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- shared_hypoDMR_down_AME_masterlist_withLogP_reordered[,c(5,6,7)]
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- shared_hypoDMR_down_AME_masterlist_withLogP_reordered[,c(8,9,10)]
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #24x3

# Order based on DBD
shared_hypoDMR_down_AME_masterlist_withLogP[order(shared_hypoDMR_down_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_down_AME_masterlist_withLogP[,"genename.x"]),]
name <- make.unique(shared_hypoDMR_down_AME_masterlist_withLogP[order(shared_hypoDMR_down_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_down_AME_masterlist_withLogP[,"genename.x"]),1],sep=".")
df.OG2 <- -shared_hypoDMR_down_AME_masterlist_withLogP[order(shared_hypoDMR_down_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_down_AME_masterlist_withLogP[,"genename.x"]),c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- shared_hypoDMR_down_AME_masterlist_withLogP[order(shared_hypoDMR_down_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_down_AME_masterlist_withLogP[,"genename.x"]),c(5,6,7)]
row.names(df.OG) <- name
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- shared_hypoDMR_down_AME_masterlist_withLogP[order(shared_hypoDMR_down_AME_masterlist_withLogP[,"DBD"],shared_hypoDMR_down_AME_masterlist_withLogP[,"genename.x"]),c(9,10,11)]
row.names(df.OG) <- name
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #6 x 3



#check closing DARs
Iri_closeDAR_up_AME_filtered <-Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$s24vIri_log2_change <0 | Iri_closeDAR_AME_genename_DEG$MelvIri_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in iri 
Iri_closeDAR_down_AME_filtered <-Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$s24vIri_log2_change >0 | Iri_closeDAR_AME_genename_DEG$MelvIri_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #downreg in iri

Mel_closeDAR_up_AME_filtered <-Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$s24vMel_log2_change <0 | Mel_closeDAR_AME_genename_DEG$MelvIri_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in mel
Mel_closeDAR_down_AME_filtered <-Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$s24vMel_log2_change >0 | Mel_closeDAR_AME_genename_DEG$MelvIri_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #downreg in iri

Shared_closeDAR_up_AME_filtered <-Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$s24vIri_log2_change <0 & Shared_closeDAR_AME_genename_DEG$s24vMel_log2_change <0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in both
Shared_closeDAR_down_AME_filtered <-Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$s24vIri_log2_change >0 & Shared_closeDAR_AME_genename_DEG$s24vMel_log2_change >0,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #downreg in both 


closeDAR_up_AME_masterlist <- unique(rbind(Iri_closeDAR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_closeDAR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)],Shared_closeDAR_up_AME_filtered[,c(1,3,4,5,6,7,8,9)]))
closeDAR_down_AME_masterlist <- unique(rbind(Iri_closeDAR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)],Mel_closeDAR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)],Shared_closeDAR_down_AME_filtered[,c(1,3,4,5,6,7,8,9)]))


closeDAR_up_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(closeDAR_up_AME_masterlist,Iri_closeDAR_up_AME_filtered,Mel_closeDAR_up_AME_filtered ,Shared_closeDAR_up_AME_filtered))
closeDAR_up_AME_masterlist_withLogP[is.na(closeDAR_up_AME_masterlist_withLogP)] <- 0
closeDAR_up_AME_masterlist_withLogP<-unique(closeDAR_up_AME_masterlist_withLogP)


closeDAR_up_AME_masterlist_withLogP_all3 <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_up_AME_masterlist_withLogP_all2a <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_up_AME_masterlist_withLogP_all2b <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_up_AME_masterlist_withLogP_all2c <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_up_AME_masterlist_withLogP_all1a <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_up_AME_masterlist_withLogP_all1b <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_up_AME_masterlist_withLogP_all1c <- closeDAR_up_AME_masterlist_withLogP[closeDAR_up_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_up_AME_masterlist_withLogP$pvalue_log >0,]

closeDAR_up_AME_masterlist_withLogP_all3 <-closeDAR_up_AME_masterlist_withLogP_all3[order(closeDAR_up_AME_masterlist_withLogP_all3[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all3[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all2a <-closeDAR_up_AME_masterlist_withLogP_all2a[order(closeDAR_up_AME_masterlist_withLogP_all2a[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all2a[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all2b <-closeDAR_up_AME_masterlist_withLogP_all2b[order(closeDAR_up_AME_masterlist_withLogP_all2b[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all2b[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all2c <-closeDAR_up_AME_masterlist_withLogP_all2c[order(closeDAR_up_AME_masterlist_withLogP_all2c[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all2c[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all1a <-closeDAR_up_AME_masterlist_withLogP_all1a[order(closeDAR_up_AME_masterlist_withLogP_all1a[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all1a[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all1b <-closeDAR_up_AME_masterlist_withLogP_all1b[order(closeDAR_up_AME_masterlist_withLogP_all1b[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all1b[,"genename.x"]),]
closeDAR_up_AME_masterlist_withLogP_all1c <-closeDAR_up_AME_masterlist_withLogP_all1c[order(closeDAR_up_AME_masterlist_withLogP_all1c[,"DBD"],closeDAR_up_AME_masterlist_withLogP_all1c[,"genename.x"]),]

closeDAR_up_AME_masterlist_withLogP_reordered <- rbind(closeDAR_up_AME_masterlist_withLogP_all3,closeDAR_up_AME_masterlist_withLogP_all2a,closeDAR_up_AME_masterlist_withLogP_all2b,closeDAR_up_AME_masterlist_withLogP_all2c,closeDAR_up_AME_masterlist_withLogP_all1a,closeDAR_up_AME_masterlist_withLogP_all1b,closeDAR_up_AME_masterlist_withLogP_all1c)
closeDAR_up_AME_masterlist_withLogP_reordered <- unique(closeDAR_up_AME_masterlist_withLogP_reordered)

library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- make.unique(closeDAR_up_AME_masterlist_withLogP_reordered[,1],sep=".")
df.OG2 <- -closeDAR_up_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- closeDAR_up_AME_masterlist_withLogP_reordered[,c(5,6,7)]
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- closeDAR_up_AME_masterlist_withLogP_reordered[,c(9,10,11)]
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #20 x 



#closeDAR downregulated genes
closeDAR_down_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(closeDAR_down_AME_masterlist,Iri_closeDAR_down_AME_filtered,Mel_closeDAR_down_AME_filtered ,Shared_closeDAR_down_AME_filtered))
closeDAR_down_AME_masterlist_withLogP[is.na(closeDAR_down_AME_masterlist_withLogP)] <- 0
closeDAR_down_AME_masterlist_withLogP<-unique(closeDAR_down_AME_masterlist_withLogP)


closeDAR_down_AME_masterlist_withLogP_all3 <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_down_AME_masterlist_withLogP_all2a <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_down_AME_masterlist_withLogP_all2b <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_down_AME_masterlist_withLogP_all2c <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log >0,]
closeDAR_down_AME_masterlist_withLogP_all1a <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_down_AME_masterlist_withLogP_all1b <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y >0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log ==0,]
closeDAR_down_AME_masterlist_withLogP_all1c <- closeDAR_down_AME_masterlist_withLogP[closeDAR_down_AME_masterlist_withLogP$pvalue_log.x ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log.y ==0 & closeDAR_down_AME_masterlist_withLogP$pvalue_log >0,]

closeDAR_down_AME_masterlist_withLogP_all3 <-closeDAR_down_AME_masterlist_withLogP_all3[order(closeDAR_down_AME_masterlist_withLogP_all3[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all3[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all2a <-closeDAR_down_AME_masterlist_withLogP_all2a[order(closeDAR_down_AME_masterlist_withLogP_all2a[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all2a[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all2b <-closeDAR_down_AME_masterlist_withLogP_all2b[order(closeDAR_down_AME_masterlist_withLogP_all2b[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all2b[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all2c <-closeDAR_down_AME_masterlist_withLogP_all2c[order(closeDAR_down_AME_masterlist_withLogP_all2c[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all2c[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all1a <-closeDAR_down_AME_masterlist_withLogP_all1a[order(closeDAR_down_AME_masterlist_withLogP_all1a[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all1a[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all1b <-closeDAR_down_AME_masterlist_withLogP_all1b[order(closeDAR_down_AME_masterlist_withLogP_all1b[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all1b[,"genename.x"]),]
closeDAR_down_AME_masterlist_withLogP_all1c <-closeDAR_down_AME_masterlist_withLogP_all1c[order(closeDAR_down_AME_masterlist_withLogP_all1c[,"DBD"],closeDAR_down_AME_masterlist_withLogP_all1c[,"genename.x"]),]

closeDAR_down_AME_masterlist_withLogP_reordered <- rbind(closeDAR_down_AME_masterlist_withLogP_all3,closeDAR_down_AME_masterlist_withLogP_all2a,closeDAR_down_AME_masterlist_withLogP_all2b,closeDAR_down_AME_masterlist_withLogP_all2c,closeDAR_down_AME_masterlist_withLogP_all1a,closeDAR_down_AME_masterlist_withLogP_all1b,closeDAR_down_AME_masterlist_withLogP_all1c)
closeDAR_down_AME_masterlist_withLogP_reordered <- unique(closeDAR_down_AME_masterlist_withLogP_reordered)

library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
name <- make.unique(closeDAR_down_AME_masterlist_withLogP_reordered[,1],sep=".")
df.OG2 <- -closeDAR_down_AME_masterlist_withLogP_reordered[,c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

ht2 =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht2

df.OG <- closeDAR_down_AME_masterlist_withLogP_reordered[,c(5,6,7)]
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- closeDAR_down_AME_masterlist_withLogP_all1clist_withLogP_reordered[,c(9,10,11)]
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #28x3 




### check known early genes for expression in shared sox10, mitfa, tfec, mitfa
sox10_AME<-Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"), all = T), list(Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDAR_AME_genename_DEG[Iri_openDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDMAR_AME_genename_DEG[Iri_openDMAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_hypoDMR_AME_genename_DEG[Mel_hypoDMR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDAR_AME_genename_DEG[Mel_openDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDMAR_AME_genename_DEG[Mel_openDMAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$genename.x == "sox10",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]))

sox10_AME[is.na(sox10_AME)] <-0
colnames(sox10_AME)<-c("gene","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","Iri_hypoDMR","Iri_openDAR","Iri_openDMAR","Mel_hypoDMR","Mel_openDAR","Mel_openDMAR","shared_hypoDMR","shared_openDAR","shared_openDMAR","Iri_closeDAR","Mel_closeDAR","shared_closeDAR")


mitfa_AME<-Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"), all = T), list(Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDAR_AME_genename_DEG[Iri_openDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDMAR_AME_genename_DEG[Iri_openDMAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_hypoDMR_AME_genename_DEG[Mel_hypoDMR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDAR_AME_genename_DEG[Mel_openDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDMAR_AME_genename_DEG[Mel_openDMAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$genename.x == "mitfa",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]))

mitfa_AME[is.na(mitfa_AME)] <-0
colnames(mitfa_AME)<-c("gene","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","Iri_hypoDMR","Iri_openDAR","Iri_openDMAR","Mel_hypoDMR","Mel_openDAR","Mel_openDMAR","shared_hypoDMR","shared_openDAR","shared_openDMAR","Iri_closeDAR","Mel_closeDAR","shared_closeDAR")


tfec_AME<-Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"), all = T), list(Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDAR_AME_genename_DEG[Iri_openDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_openDMAR_AME_genename_DEG[Iri_openDMAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_hypoDMR_AME_genename_DEG[Mel_hypoDMR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDAR_AME_genename_DEG[Mel_openDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_openDMAR_AME_genename_DEG[Mel_openDMAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_hypoDMR_AME_genename_DEG[Shared_hypoDMR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDAR_AME_genename_DEG[Shared_openDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_openDMAR_AME_genename_DEG[Shared_openDMAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Iri_closeDAR_AME_genename_DEG[Iri_closeDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Mel_closeDAR_AME_genename_DEG[Mel_closeDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]
,Shared_closeDAR_AME_genename_DEG[Shared_closeDAR_AME_genename_DEG$genename.x == "tfec",c("genename.x","pvalue_log","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")]))

tfec_AME[is.na(tfec_AME)] <-0
colnames(tfec_AME)<-c("gene","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","Iri_hypoDMR","Iri_openDAR","Iri_openDMAR","Mel_hypoDMR","Mel_openDAR","Mel_openDMAR","shared_hypoDMR","shared_openDAR","shared_openDMAR","Iri_closeDAR","Mel_closeDAR","shared_closeDAR")


######### IDENTIFY TOP 3000 REGIONS for DMR, DAR AND DMAR #############
combined_dmr <- read.table("Combined_DMRs_d30_p0.01.txt",sep = "\t", header = F, stringsAsFactors = F)
combined_dmr$size <- combined_dmr$V3-combined_dmr$V2
hist(combined_dmr$size,breaks = 1000, xlim = c(0,2000),xlab = "DMR Size", main = "Combined DMR size distribution (d30)", col = "grey")

cDMR15v24_d30_p0.01 <- read.table("Combined_DMR15v24_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vMel_d30_p0.01 <- read.table("Combined_DMR24vMel_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vIri_d30_p0.01 <- read.table("Combined_DMR24vIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMRMelvIri_d30_p0.01 <- read.table("Combined_DMRMelvIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)

cDMR15v24_d30_p0.01$chrompos<- paste(cDMR15v24_d30_p0.01$V1,":",cDMR15v24_d30_p0.01$V2,"-",cDMR15v24_d30_p0.01$V3,sep = "")
cDMR24vMel_d30_p0.01$chrompos<- paste(cDMR24vMel_d30_p0.01$V1,":",cDMR24vMel_d30_p0.01$V2,"-",cDMR24vMel_d30_p0.01$V3,sep = "")
cDMR24vIri_d30_p0.01$chrompos<- paste(cDMR24vIri_d30_p0.01$V1,":",cDMR24vIri_d30_p0.01$V2,"-",cDMR24vIri_d30_p0.01$V3,sep = "")
cDMRMelvIri_d30_p0.01$chrompos<- paste(cDMRMelvIri_d30_p0.01$V1,":",cDMRMelvIri_d30_p0.01$V2,"-",cDMRMelvIri_d30_p0.01$V3,sep = "")

combined_dmr_wInfo <-Reduce(function(x, y) merge(x, y, by=c("chrompos"), all = T), list(cDMR15v24_d30_p0.01[,c("V1","V2","V3","chrompos","V11","V12")],cDMR24vMel_d30_p0.01[,c("chrompos","V11","V12")],cDMR24vIri_d30_p0.01[,c("chrompos","V11","V12")],cDMRMelvIri_d30_p0.01[,c("chrompos","V11","V12")]))
colnames(combined_dmr_wInfo) <- c("chrompos","chr","start","end","s15vs24","s15vs24_areastat","s24vMel","s24vMel_areastat","s24vIri","s24vIri_areastat","MelvsIri","MelvsIri_areastat")
combined_dmr_wInfo$size <- combined_dmr_wInfo$end-combined_dmr_wInfo$start
combined_dmr_wInfo[combined_dmr_wInfo == "."] <- 0
combined_dmr_wInfo$s15vs24 <- as.numeric(combined_dmr_wInfo$s15vs24)
combined_dmr_wInfo$s24vMel <- as.numeric(combined_dmr_wInfo$s24vMel)
combined_dmr_wInfo$s24vIri <- as.numeric(combined_dmr_wInfo$s24vIri)
combined_dmr_wInfo$MelvsIri <- as.numeric(combined_dmr_wInfo$MelvsIri)
combined_dmr_wInfo$s15vs24_areastat <- as.numeric(combined_dmr_wInfo$s15vs24_areastat)
combined_dmr_wInfo$s24vMel_areastat <- as.numeric(combined_dmr_wInfo$s24vMel_areastat)
combined_dmr_wInfo$s24vIri_areastat <- as.numeric(combined_dmr_wInfo$s24vIri_areastat)
combined_dmr_wInfo$MelvsIri_areastat <- as.numeric(combined_dmr_wInfo$MelvsIri_areastat)

# Find 5000 DMRs with highest areastat
Mel_hypo <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel_areastat>0 & combined_dmr_wInfo$s24vIri_areastat<=0,]
Mel_hypo <-Mel_hypo[order(-Mel_hypo$s24vMel_areastat),]
Mel_hypo <- Mel_hypo[c(1:3000),]
write.table(Mel_hypo[,c(2,3,4,1)],"Mel-specific-hypoDMR-top3000_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)

Mel_hypo <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel_areastat>0 & combined_dmr_wInfo$s24vIri_areastat<=0,]
Mel_hypo <-Mel_hypo[order(-Mel_hypo$s24vMel_areastat),]
Mel_hypo <- Mel_hypo[c(1:5000),]
write.table(Mel_hypo[,c(2,3,4,1)],"Mel-specific-hypoDMR-top5000_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)

Iri_hypo <- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri_areastat>0 & combined_dmr_wInfo$s24vMel_areastat<=0,]
Iri_hypo <-Iri_hypo[order(-Iri_hypo$s24vIri_areastat),]
Iri_hypo <- Iri_hypo[c(1:3000),]
write.table(Iri_hypo[,c(2,3,4,1)],"Iri-specific-hypoDMR-top3000_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)


Iri_hypo <- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri_areastat>0 & combined_dmr_wInfo$s24vMel_areastat<=0,]
Iri_hypo <-Iri_hypo[order(-Iri_hypo$s24vIri_areastat),]
Iri_hypo <- Iri_hypo[c(1:5000),]
write.table(Iri_hypo[,c(2,3,4,1)],"Iri-specific-hypoDMR-top5000_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)
```

```{bash}
###################################### {Bash script} ######################################
sort -k1,1 -k2,2n Mel-specific-hypoDMR-top3000_forMetascape.bed > Mel-specific-hypoDMR-top3000_forMetascape_sorted.bed
sort -k1,1 -k2,2n Mel-specific-hypoDMR-top5000_forMetascape.bed > Mel-specific-hypoDMR-top5000_forMetascape_sorted.bed
sort -k1,1 -k2,2n Iri-specific-hypoDMR-top5000_forMetascape.bed > Iri-specific-hypoDMR-top5000_forMetascape_sorted.bed
sort -k1,1 -k2,2n Iri-specific-hypoDMR-top3000_forMetascape.bed > Iri-specific-hypoDMR-top3000_forMetascape_sorted.bed
sort -k1,1 -k2,2n Danio_rerio.GRCz10.85.GENE.PROMOTER.bed > promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed


bedtools closest -a Mel-specific-hypoDMR-top3000_forMetascape_sorted.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Mel-specific-hypoDMR-top3000_forMetascape_closestpromoter.bed
bedtools closest -a Iri-specific-hypoDMR-top3000_forMetascape_sorted.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Iri-specific-hypoDMR-top3000_forMetascape_closestpromoter.bed
bedtools closest -a Mel-specific-hypoDMR-top5000_forMetascape_sorted.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Mel-specific-hypoDMR-top5000_forMetascape_closestpromoter.bed
bedtools closest -a Iri-specific-hypoDMR-top5000_forMetascape_sorted.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Iri-specific-hypoDMR-top5000_forMetascape_closestpromoter.bed


Iri_specifc_DMAR <- DMAR[DMAR$DMRs24vIri >0 & DMAR$DARs24vIri >0 & DMAR$DMRs24vMel <=0 & DMAR$DARs24vMel <=0,]
write.table(Iri_specifc_DMAR[,c(1,2,3,4)],"Iri-specific-DMAR_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)
Mel_specifc_DMAR <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DARs24vMel >0 & DMAR$DMRs24vIri <=0 & DMAR$DARs24vIri <=0,]
write.table(Mel_specifc_DMAR[,c(1,2,3,4)],"Mel-specific-DMAR_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)

bedtools closest -a Iri-specific-DMAR_forMetascape.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Iri-specific-DMAR_forMetascape_closestpromoter.bed
bedtools closest -a Mel-specific-DMAR_forMetascape.bed -b Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > Mel-specific-DMAR_forMetascape_closestpromoter.bed
```

```{R}
###################################### {R script} ######################################
##FOR DAR, pick highest fold change top 5000
Mel_specific_openDAR <-DMAR[DMAR$DARs24vMel > 0 & DMAR$DARs24vIri <=0,]
Mel_specific_openDAR <- Mel_specific_openDAR[order(-Mel_specific_openDAR$DARs24vMel),]
Mel_specific_openDAR <- Mel_specific_openDAR[c(1:5000),]
write.table(Mel_specific_openDAR[,c(1,2,3,4)],"Mel-specific-openDAR_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)

Iri_specific_openDAR <-DMAR[DMAR$DARs24vMel <= 0 & DMAR$DARs24vIri >0,]
Iri_specific_openDAR <- Iri_specific_openDAR[order(-Iri_specific_openDAR$DARs24vIri),]
Iri_specific_openDAR <- Iri_specific_openDAR[c(1:5000),]
write.table(Iri_specific_openDAR[,c(1,2,3,4)],"Iri-specific-openDAR_forMetascape.bed",sep = "\t", col.names = F, row.names = F, quote = F)
```
```{bash}
###################################### {Bash script} ######################################
bedtools closest -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Iri-specific-openDAR_forMetascape.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Iri-specific-openDAR_forMetascape_closestpromoter.bed
bedtools closest -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Mel-specific-openDAR_forMetascape.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER_sorted.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Mel-specific-openDAR_forMetascape_closestpromoter.bed
```
```{R}
###################################### {R script} ######################################
#make heatmap of methylation status in DARs
library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
All <- DMAR[DMAR$DARs15vs24>0 & DMAR$aveCpGcount > 0,]
name <- All[,4]
df.OG2 <- All[,c("Meth_s15","Meth_s24","Meth_Mel","Meth_Iri")]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("15somite vs 24hpf", "24hpf vs Mel", "24hpf vs Iri")
kclus2 <- kmeans(df.OG2,10)
split2 <-  kclus2$cluster
ht2 = Heatmap(df.OG2, column_title = "Methylation of DARs",name= "heyo",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), cluster_rows = TRUE, cluster_columns = FALSE,show_row_names = FALSE,split = split2)
ht2

library(dendextend)
row_dend = hclust(dist(df.OG2),method="ward.D")
ht3 = Heatmap(df.OG2, column_title = "Methylation of 15somite to 24hpf DARs",name= "Methylaton %",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), cluster_rows = color_branches(row_dend, k = 4), cluster_columns = FALSE,show_row_names = FALSE)
ht3


############## Promoter centric Analysis (got danRer10 promoter (1kb upstream TSS) from UCSC Table browser) #################
############ {Bash script} ################
#pull out gene info from GTF file
awk '($3=="gene"){OFS="\t"; if ($7~/+/){print $1,$4-1000,$4+500,$10}; if ($7~/-/){print $1,$5-500,$5+1000,$10}}' /bar/genomes/danRer10/Ensembl/Danio_rerio.GRCz10.85.gtf| sed 's/[";]//g;' > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed
awk '($3=="gene"){OFS="\t"; if ($7~/+/){print $1,$4,$5,$10}; if ($7~/-/){print $1,$4,$5,$10}}' /bar/genomes/danRer10/Ensembl/Danio_rerio.GRCz10.85.gtf| sed 's/[";]//g;' > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.LOCATION.bed
############ {R script} ###################
prom <- read.table("Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t", stringsAsFactors =F, header = F, quote ="")
prom[prom$V2 <0,]$V2 <- 0
prom$chrompos <- paste(prom$V1,":",prom$V2,"-",prom$V3,sep = "")

write.table(prom,"Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t",quote = F, col.names = F, row.names = F)

#overlap CpG methylation
############ {Bash script} ################
bedtools intersect -wo -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/BS_only/15somite_NCC_Combined_DSS.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_15s.txt
bedtools intersect -wo -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/BS_only/24hpf_NCC_Combined_DSS.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_24s.txt
bedtools intersect -wo -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/BS_only/Mel_Combined_DSS.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Mel.txt
bedtools intersect -wo -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/BS_only/Iri_Combined_DSS.bed > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Iri.txt

############ {R script} ################
prom_meth15<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_meth24<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_24s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_methM<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_methI<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

prom_meth15 <- prom_meth15[prom_meth15$V9 >=5,]#filter CpGs with >=5 cov
prom_meth24 <- prom_meth24[prom_meth24$V9 >=5,]#filter CpGs with >=5 cov
prom_methM <- prom_methM[prom_methM$V9 >=5,]#filter CpGs with >=5 cov
prom_methI <- prom_methI[prom_methI$V9 >=5,]#filter CpGs with >=5 cov

prom_meth15$meth <- prom_meth15$V10/prom_meth15$V9*100
prom_meth24$meth <- prom_meth24$V10/prom_meth24$V9*100
prom_methM$meth <- prom_methM$V10/prom_methM$V9*100
prom_methI$meth <- prom_methI$V10/prom_methI$V9*100

library(plyr)
prom_meth15_aggregate <- ddply(prom_meth15,.(V1,V2,V3,V4,V5),summarise,CpGcount=sum(V11),AveMeth=mean(meth))
prom_meth24_aggregate <- ddply(prom_meth24,.(V1,V2,V3,V4,V5),summarise,CpGcount=sum(V11),AveMeth=mean(meth))
prom_methM_aggregate <- ddply(prom_methM,.(V1,V2,V3,V4,V5),summarise,CpGcount=sum(V11),AveMeth=mean(meth))
prom_methI_aggregate <- ddply(prom_methI,.(V1,V2,V3,V4,V5),summarise,CpGcount=sum(V11),AveMeth=mean(meth))

m<-Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3"), all.x = T), list(prom,prom_meth15_aggregate[,c(1,2,3,6,7)],prom_meth24_aggregate[,c(1,2,3,6,7)],prom_methM_aggregate[,c(1,2,3,6,7)],prom_methI_aggregate[,c(1,2,3,6,7)]))
```
```{bash}
###################################### {Bash script} ######################################
#overlap IDR peaks
bedtools intersect -c -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks15s.bed
bedtools intersect -c -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks24s.bed
bedtools intersect -c -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksMel.bed
bedtools intersect -c -a /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/promoter_centric/ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksIri.bed
```
```{R}
###################################### {R script} ######################################
prom_IDR15<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks15s.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDR24<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks24s.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDRM<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksMel.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDRI<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksIri.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)

m_IDR<-Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3"), all.x = T), list(m,prom_IDR15[,c(1,2,3,6)],prom_IDR24[,c(1,2,3,6)],prom_IDRM[,c(1,2,3,6)],prom_IDRI[,c(1,2,3,6)]))
names(m_IDR) <- c("chr","start","end","gene","chrompos","s15_CpGcount","s15_AveMeth","s24_CpGcount","s24_AveMeth","Mel_CpGcount","Mel_AveMeth","Iri_CpGcount","Iri_AveMeth","s15_IDRpeak","s24_IDRpeak","Mel_IDRpeak","Mel_IDRpeak")

m_IDR[is.na(m_IDR)]<-0

write.table(m_IDR, "Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo.bed", sep = "\t", quote =F, col.names =T, row.names = F)
write.table(m_IDR, "Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo2.bed", sep = "\t", quote =F, col.names =F, row.names = F)

############ {Bash script} ################
bedtools intersect -wao -a Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo2.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/All_DMAR_Combined_wINFO.annotated2.bed > Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation.bed

############ {R script} ################
prom_DMAR <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
names(prom_DMAR) <- c("chr","start","end","gene","chrompos","s15_CpGcount","s15_AveMeth","s24_CpGcount","s24_AveMeth","Mel_CpGcount","Mel_AveMeth","Iri_CpGcount","Iri_AveMeth","s15_IDRpeak","s24_IDRpeak","Mel_IDRpeak","Iri_IDRpeak", "chr_DMAR","start_DMAR" , "end_DMAR"  ,    "chrompos_DMAR",  "DMARsize", "DMRsize" , "DMRs15vs24" , "DMRs24vMel", "DMRs24vIri",  "DMRMelvIri" ,"DARsize" ,   "DARs15vs24",  "DARs24vMel",
"DARs24vIri","DARMelvIri", "IDR_s15" ,"IDR_s24"  , "IDR_M"  , "IDR_I" ,    "Meth_s15"   ,"Meth_s24" , "Meth_Mel" , "Meth_Iri"   , "CpG_s15",  "CpG_s24",  "CpG_Mel" ,"CpG_Iri", "aveCpG"   , "CpGdensity100bp","Annotation", "Annotation2","prom_DMAR_overlap")
prom_DMAR[is.na(prom_DMAR)] <- 0
prom_DMAR[prom_DMAR =="."] <- 0


#add gene expression info
DEG2_genenames<- read.table("<path DESeq2_051718/DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) #from BioMart
colnames(DEG2_genenames) <- c("gene","genename")
DEG3 <- merge(DEG2,DEG2_genenames, by = "gene", all.x = T)
DEG4<-merge(geneTPM, DEG3,by = "gene", all.x =T)##ALL GENES
DEG4[is.na(DEG4)] <-0
DEG4<-DEG4[DEG4$s15 > 5 | DEG4$s24 > 5 |DEG4$Mel > 5 |DEG4$Iri > 5,] #filter genes with at least 5 rpkm

prom_DMAR2 <- merge(prom_DMAR,DEG4[,c(1:13)],by="gene",all.x = T)
prom_DMAR2 <- prom_DMAR2[,c(2,3,4,1,5:61)]
prom_DMAR2[is.na(prom_DMAR2)] <-0
write.table(prom_DMAR2,"Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation_wDEGs.bed", sep = "\t", col.names = T, row.names =F, quote =F)


prom_DMAR2 <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation_wDEGs.bed", sep = "\t", header = T, stringsAsFactors =F, quote = "")
prom_DMAR2[,c(22:46)] <- sapply(prom_DMAR2[,c(22:46)],as.numeric)

###DEG promoter centric analysis
prom_DMAR2_Mel_On <- prom_DMAR2[prom_DMAR2$s24vMel_log2_change < 0,]
prom_DMAR2_Iri_On <- prom_DMAR2[prom_DMAR2$s24vIri_log2_change < 0,]

prom_DMAR2_Mel_Off <- prom_DMAR2[prom_DMAR2$s24vMel_log2_change > 0 ,]
prom_DMAR2_Iri_Off <- prom_DMAR2[prom_DMAR2$s24vIri_log2_change > 0,]

#Mel On
prom_DMAR2_Mel_On_nochange <- prom_DMAR2_Mel_On[prom_DMAR2_Mel_On$DMRs24vMel == 0 & prom_DMAR2_Mel_On$DARs24vMel ==0,]
prom_DMAR2_Mel_On_only_hypo<- prom_DMAR2_Mel_On[prom_DMAR2_Mel_On$DMRs24vMel > 0 & prom_DMAR2_Mel_On$DARs24vMel ==0,]
prom_DMAR2_Mel_On_only_opening<- prom_DMAR2_Mel_On[prom_DMAR2_Mel_On$DMRs24vMel == 0 & prom_DMAR2_Mel_On$DARs24vMel >0,]
prom_DMAR2_Mel_On_both_hypo_open<- prom_DMAR2_Mel_On[prom_DMAR2_Mel_On$DMRs24vMel > 0 & prom_DMAR2_Mel_On$DARs24vMel >0,]

gl <- c(prom_DMAR2_Mel_On_nochange$chrompos,prom_DMAR2_Mel_On_only_hypo$chrompos,prom_DMAR2_Mel_On_only_opening$chrompos,prom_DMAR2_Mel_On_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Mel_On_others <- prom_DMAR2_Mel_On[prom_DMAR2_Mel_On$chrompos %out% gl,]

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

mypalette1 <- c("#ee4035","#ffcb85","#fdf498","#7bc043","#0392cf")

library(waffle)
vals <- c(nrow(prom_DMAR2_Mel_On_nochange),nrow(prom_DMAR2_Mel_On_only_hypo),nrow(prom_DMAR2_Mel_On_only_opening),nrow(prom_DMAR2_Mel_On_both_hypo_open),nrow(prom_DMAR2_Mel_On_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hypoDMR", "Only openDAR", "Opening DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

aa<-waffle::waffle(vals,colors = mypalette1,rows = 18,size = 0.5,title = "Mel DEG's promoter dynamics (On)")


#iri on
prom_DMAR2_Iri_On_nochange <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri == 0 & prom_DMAR2_Iri_On$DARs24vIri ==0,]
prom_DMAR2_Iri_On_only_hypo<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri > 0 & prom_DMAR2_Iri_On$DARs24vIri ==0,]
prom_DMAR2_Iri_On_only_opening<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri == 0 & prom_DMAR2_Iri_On$DARs24vIri >0,]
prom_DMAR2_Iri_On_both_hypo_open<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri > 0 & prom_DMAR2_Iri_On$DARs24vIri >0,]

gl <- c(prom_DMAR2_Iri_On_nochange$chrompos,prom_DMAR2_Iri_On_only_hypo$chrompos,prom_DMAR2_Iri_On_only_opening$chrompos,prom_DMAR2_Iri_On_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Iri_On_others <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$chrompos %out% gl,]

mypalette1 <- c("#ee4035","#ffcb85","#fdf498","#7bc043","#0392cf")

library(waffle)
vals <- c(nrow(prom_DMAR2_Iri_On_nochange),nrow(prom_DMAR2_Iri_On_only_hypo),nrow(prom_DMAR2_Iri_On_only_opening),nrow(prom_DMAR2_Iri_On_both_hypo_open),nrow(prom_DMAR2_Iri_On_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hypoDMR", "Only openDAR", "Opening DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

bb<-waffle::waffle(vals,colors = mypalette1,rows = 18,size = 0.5,title = "Iri DEG's promoter dynamics (On)")


#Mel off
prom_DMAR2_Mel_Off_nochange <- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel == 0 & prom_DMAR2_Mel_Off$DARs24vMel ==0,]
prom_DMAR2_Mel_Off_only_hypo<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel > 0 & prom_DMAR2_Mel_Off$DARs24vMel ==0,]
prom_DMAR2_Mel_Off_only_opening<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel == 0 & prom_DMAR2_Mel_Off$DARs24vMel <0,]
prom_DMAR2_Mel_Off_both_hypo_open<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel > 0 & prom_DMAR2_Mel_Off$DARs24vMel <0,]

gl <- c(prom_DMAR2_Mel_Off_nochange$chrompos,prom_DMAR2_Mel_Off_only_hypo$chrompos,prom_DMAR2_Mel_Off_only_opening$chrompos,prom_DMAR2_Mel_Off_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Mel_Off_others <- prom_DMAR2_Mel_On[prom_DMAR2_Mel_Off$chrompos %out% gl,]

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

mypalette1 <- c("#ee4035","#ffcb85","#fdf498","#7bc043","#0392cf")

library(waffle)
vals <- c(nrow(prom_DMAR2_Mel_Off_nochange),nrow(prom_DMAR2_Mel_Off_only_hypo),nrow(prom_DMAR2_Mel_Off_only_opening),nrow(prom_DMAR2_Mel_Off_both_hypo_open),nrow(prom_DMAR2_Mel_Off_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hypoDMR", "Only closingDAR", "Closing DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

cc<-waffle::waffle(vals,colors = mypalette1,rows = 18,size = 0.5,title = "Mel DEG's promoter dynamics (Off)")


#iri off
prom_DMAR2_Iri_Off_nochange <- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri == 0 & prom_DMAR2_Iri_Off$DARs24vIri ==0,]
prom_DMAR2_Iri_Off_only_hypo<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri < 0 & prom_DMAR2_Iri_Off$DARs24vIri ==0,]
prom_DMAR2_Iri_Off_only_opening<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri == 0 & prom_DMAR2_Iri_Off$DARs24vIri <0,]
prom_DMAR2_Iri_Off_both_hypo_open<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri < 0 & prom_DMAR2_Iri_Off$DARs24vIri <0,]

gl <- c(prom_DMAR2_Iri_Off_nochange$chrompos,prom_DMAR2_Iri_Off_only_hypo$chrompos,prom_DMAR2_Iri_Off_only_opening$chrompos,prom_DMAR2_Iri_Off_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Iri_Off_others <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_Off$chrompos %out% gl,]

mypalette1 <- c("#ee4035","#ffcb85","#fdf498","#7bc043","#0392cf")

library(waffle)
vals <- c(nrow(prom_DMAR2_Iri_Off_nochange),nrow(prom_DMAR2_Iri_Off_only_hypo),nrow(prom_DMAR2_Iri_Off_only_opening),nrow(prom_DMAR2_Iri_Off_both_hypo_open),nrow(prom_DMAR2_Iri_Off_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hyperDMR", "Only closingDAR", "Closing DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

dd<-waffle::waffle(vals,colors = mypalette1,rows =18,size = 0.5,title = "Iri DEG's promoter dynamics (Off)")

multiplot(aa,bb,cc,dd,cols=1)


############# Check conservation across DMRS, DMARs and DARs using deeptools #####################
Iri_DMR_all <- DMAR[DMAR$DMRs24vIri !=0,]
Mel_DMR_all <- DMAR[DMAR$DMRs24vMel !=0,]
Iri_DAR_all <- DMAR[DMAR$DARs24vIri !=0,]
Mel_DAR_all <- DMAR[DMAR$DARs24vMel !=0,]
Iri_DMAR_all <- DMAR[DMAR$DARs24vIri !=0 & DMAR$DMRs24vIri !=0,]
Mel_DMAR_all <- DMAR[DMAR$DARs24vMel !=0 & DMAR$DMRs24vMel !=0,]

write.table(Iri_DMR_all[,c(1,2,3)],"Iri_DMR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)
write.table(Mel_DMR_all[,c(1,2,3)],"Mel_DMR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)
write.table(Iri_DAR_all[,c(1,2,3)],"Iri_DAR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)
write.table(Mel_DAR_all[,c(1,2,3)],"Mel_DAR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)
write.table(Iri_DMAR_all[,c(1,2,3)],"Iri_DMAR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)
write.table(Mel_DMAR_all[,c(1,2,3)],"Mel_DMAR_all.bed",sep = "\t", col.names = F, row.names =F,quote =F)

# Shuffle the regions and use them as control
bedtools shuffle -i Iri_DAR_all.bed -g /bar/genomes/danRer10/conservation/danRer10.chrom.sizes > Iri_shuffled_genome.bed ##for control in plotting
bedtools shuffle -i Mel_DAR_all.bed -g /bar/genomes/danRer10/conservation/danRer10.chrom.sizes > Mel_shuffled_genome.bed ##for control in plotting

############ {Bash script} #############
# Download these files:
##danRer10.vertebrate.phastCons8way.bw 
##danRer10.vertebrate.phyloP8way.bw 

module load python3

computeMatrix reference-point -p 20 -S danRer10.vertebrate.phastCons8way.bw -R Iri_DMR_all.bed Iri_DAR_all.bed Iri_DMAR_all.bed Iri_shuffled_genome.bed --missingDataAsZero --referencePoint center -a 5000 -b 5000 -o Iri_phastCons_matrix.gz --outFileNameMatrix Iri_phastCons_matrix.tab --outFileSortedRegions Iri_phastCons_matrix.bed
plotHeatmap -m Iri_phastCons_matrix.gz -out Iri_phastCons_heatmap_only.png -T "Iridophore-specific region phastCons" --missingDataColor 0.8 --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 20

computeMatrix reference-point -p 20 -S danRer10.vertebrate.phastCons8way.bw -R Mel_DMR_all.bed Mel_DAR_all.bed Mel_DMAR_all.bed Mel_shuffled_genome.bed --missingDataAsZero --referencePoint center -a 5000 -b 5000 -o Mel_phastCons_matrix.gz --outFileNameMatrix Mel_phastCons_matrix.tab --outFileSortedRegions Mel_phastCons_matrix.bed
plotHeatmap -m Mel_phastCons_matrix.gz -out Mel_phastCons_heatmap_only.png -T "Melanophore-specific region phastCons" --missingDataColor 0.8 --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 20


computeMatrix reference-point -p 20 -S danRer10.vertebrate.phyloP8way.bw -R Iri_DMR_all.bed Iri_DAR_all.bed Iri_DMAR_all.bed Iri_shuffled_genome.bed --missingDataAsZero --referencePoint center -a 5000 -b 5000 -o Iri_phyloP_matrix.gz --outFileNameMatrix Iri_phyloP_matrix.tab --outFileSortedRegions Iri_phyloP_matrix.bed
plotHeatmap -m Iri_phyloP_matrix.gz -out Iri_phyloP_heatmap_only.png -T "Iridophore-specific region phyloP" --missingDataColor 0.8 --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 20

computeMatrix reference-point -p 20 -S danRer10.vertebrate.phyloP8way.bw -R Mel_DMR_all.bed Mel_DAR_all.bed Mel_DMAR_all.bed Mel_shuffled_genome.bed --missingDataAsZero --referencePoint center -a 5000 -b 5000 -o Mel_phyloP_matrix.gz --outFileNameMatrix Mel_phyloP_matrix.tab --outFileSortedRegions Mel_phyloP_matrix.bed
plotHeatmap -m Mel_phyloP_matrix.gz -out Mel_phyloP_heatmap_only.png -T "Melanophore-specific region phyloP" --missingDataColor 0.8 --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 20
```
```{R}
############ {R script} #############
####### Make WashU browser bed tracks for DMR,DAR,DMAR ###########
DMAR_new<-read.table("All_DMAR_Combined_wINFO.annotated.bed",sep = '\t', header = T, quote = "", stringsAsFactors = F)
DMAR_new$Type <- "None"
DMAR_new$browser <- 0

DMAR_new[DMAR_new$DMRsize >0 & DMAR_new$DARsize == 0,]$Type <- "soloDMR"
DMAR_new[DMAR_new$DMRsize == 0 & DMAR_new$DARsize > 0,]$Type <- "soloDAR"
DMAR_new[DMAR_new$DMRsize >0 & DMAR_new$DARsize > 0,]$Type <- "DMAR"
DMAR_new[DMAR_new$Type == "soloDMR",]$browser <- 1
DMAR_new[DMAR_new$Type == "soloDAR",]$browser <- 2
DMAR_new[DMAR_new$Type == "DMAR",]$browser <- 3

DMAR_Iri <- DMAR_new[DMAR_new$DMRs24vIri  !=0 | DMAR_new$DARs24vIri != 0,]
DMAR_Mel <- DMAR_new[DMAR_new$DMRs24vMel  !=0 | DMAR_new$DARs24vMel != 0,]
DMAR_24 <- DMAR_new[DMAR_new$DMRs15vs24  !=0 | DMAR_new$DARs15vs24 != 0,]
DMAR_MvI <- DMAR_new[DMAR_new$DMRMelvIri  !=0 | DMAR_new$DARMelvIri != 0,]

setwd("/bar/jjang/Pigment_project/DMAR")
write.table(DMAR_Iri[,c(1,2,3,33)], "DMAR_Iri.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(DMAR_Mel[,c(1,2,3,33)], "DMAR_Mel.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(DMAR_24[,c(1,2,3,33)], "DMAR_24hpf.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(DMAR_MvI[,c(1,2,3,33)], "DMAR_MelvsIri.bed", sep = "\t", quote = F, col.names = F, row.names =F)



####### USE FIMO TO scan FOR IRI-SPECIFIC TF motifs ###########
DMAR_Iri_accessible <- DMAR_new[DMAR_new$DMRs24vIri  >0 | DMAR_new$DARs24vIri >0 | DMAR_new$DMRMelvIri  >0 | DMAR_new$DARMelvIri <0,]
write.table(DMAR_Iri_accessible[,c(1,2,3,4)], "DMAR_Iri_Accesible_for_FIMO.bed", sep = "\t", quote = F, col.names = F, row.names =F)

DMAR_Mel_accessible <- DMAR_new[DMAR_new$DMRs24vMel  >0 | DMAR_new$DARs24vMel >0 | DMAR_new$DMRMelvIri  <0 | DMAR_new$DARMelvIri >0,]
write.table(DMAR_Mel_accessible[,c(1,2,3,4)], "DMAR_Mel_Accesible_for_FIMO.bed", sep = "\t", quote = F, col.names = F, row.names =F)


twoBitToFa -bed=FIMO_motif_scan/DMAR_Iri_Accesible_for_FIMO.bed danRer10.2bit DMAR_Iri_Accesible_for_FIMO.fasta

# module load meme
fimo Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta
mast Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta
mcast Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta

python3 Parse_fimo.py fimo.txt

FIMO <- read.table("fimo_parsed.txt",sep = "\t", header = T, quote = "", stringsAsFactors =F)
DMAR_Iri_FIMO <- merge(DMAR_Iri_accessible[,c(1:32)],FIMO, by = "chrompos", all.x = T)
DMAR_Iri_FIMO[is.na(DMAR_Iri_FIMO)] <- 0

write.table(DMAR_Iri_FIMO[,c(2:4,1,5:52)], "accessibleDMAR_Iri_wFIMOmotif.txt", sep = "\t", col.names = T, row.names = F, quote =F)
write.table(DMAR_Iri_FIMO[,c(2:4,1,5:52)], "saccessibleDMAR_Iri_wFIMOmotif_noheader.txt", sep = "\t", col.names = F, row.names = F, quote =F)

name <- DMAR_Iri_FIMO[,1]
df.OG2 <- DMAR_Iri_FIMO[,c(33:52)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(DMAR_Iri_FIMO[,c(33:52)])
kclus2 <- kmeans(df.OG2,20)
split2 <-  kclus2$cluster
ht2 = Heatmap(df.OG2, column_title = "TFs",name= "Motif\nOccurrence",col = colorRamp2(c(0, 1, 3), c("#fffef2","#ffb0b0","#800000")), 
    cluster_rows = F, cluster_columns = T,show_row_names = FALSE,split = split2,clustering_method_columns= "ward.D2")
#ht2 = Heatmap(df.OG2, column_title = "TFs",name= "Motif\nOccurrence",col = colorRamp2(c(0, 1, 3), c("#fffef2","#ffb0b0","#800000")), 
#    cluster_rows = T, cluster_columns = F,show_row_names = FALSE)
ht2

DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr18" & DMAR_Iri_FIMO$start > 16227000 & DMAR_Iri_FIMO$end < 16260000,] #ALX1 
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr8" & DMAR_Iri_FIMO$start > 24785000 & DMAR_Iri_FIMO$end < 24815000,] #ALX3 
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr7" & DMAR_Iri_FIMO$start > 26632000 & DMAR_Iri_FIMO$end < 26704000,] #ALX4b 
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr18" & DMAR_Iri_FIMO$start > 38323000 & DMAR_Iri_FIMO$end < 38383000,] #ALX4b 
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr6" & DMAR_Iri_FIMO$start > 15861000 & DMAR_Iri_FIMO$end < 15909000,] #gbx2 
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr18" & DMAR_Iri_FIMO$start > 48281000 & DMAR_Iri_FIMO$end < 48373000,] #CABZ01069595.1 aka ETS1
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr4" & DMAR_Iri_FIMO$start > 5975000 & DMAR_Iri_FIMO$end < 6041000,] #tfec
DMAR_Iri_FIMO[DMAR_Iri_FIMO$chr == "chr11" & DMAR_Iri_FIMO$start > 40620000 & DMAR_Iri_FIMO$end < 40790000,] #pax7a


ALX1 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$ALX1 > 0,] #2167
ALX3 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$ALX3 > 0,] #2476
ALX4 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$ALX4 > 0,] #1988
GBX2 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$GBX2 > 0,] #1721
SOX10 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$SOX10 > 0,] #6315
SOX9 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$SOX9 > 0,] #2773
ETS1 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$ETS1 > 0,] #2271
PAX7 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$PAX7 > 0,] #1600
TFEC <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$TFEC > 0,] #4848
SP3 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$SP3 > 0,] #3315
TFAP2A <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$TFAP2A > 0,] #1846
JUNB <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$JUNB > 0,] #3132
NR4A1 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$NR4A1 > 0,] #3419
BHLHE41 <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$BHLHE41 > 0,] #1689


write.table(ALX1[,c(2,3,4)], "DMAR_Iri_Accesible_wALX1_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ALX3[,c(2,3,4)], "DMAR_Iri_Accesible_wALX3_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ALX4[,c(2,3,4)], "DMAR_Iri_Accesible_wALX4_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(GBX2[,c(2,3,4)], "DMAR_Iri_Accesible_wGBX2_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SOX10[,c(2,3,4)], "DMAR_Iri_Accesible_wSOX10_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SOX9[,c(2,3,4)], "DMAR_Iri_Accesible_wSOX9_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ETS1[,c(2,3,4)], "DMAR_Iri_Accesible_wETS1_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(PAX7[,c(2,3,4)], "DMAR_Iri_Accesible_wPAX7_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(TFEC[,c(2,3,4)], "DMAR_Iri_Accesible_wTFEC_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SP3[,c(2,3,4)], "DMAR_Iri_Accesible_wSP3motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(TFAP2A[,c(2,3,4)], "DMAR_Iri_Accesible_wTFAP2_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(JUNB[,c(2,3,4)], "DMAR_Iri_Accesible_wJUNB_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(NR4A1[,c(2,3,4)], "DMAR_Iri_Accesible_wNR4A1_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(BHLHE41[,c(2,3,4)], "DMAR_Iri_Accesible_wBHLHE41_motif.bed", sep = "\t", quote = F, col.names = F, row.names =F)

### pull out regions with multiple alx motifs ##
DMAR_Iri_FIMO_alx <- DMAR_Iri_FIMO[DMAR_Iri_FIMO$ALX1 > 0 | DMAR_Iri_FIMO$ALX3 > 0 |DMAR_Iri_FIMO$ALX4 > 0,c(2:4,1,5:36,39,44)]
DMAR_Iri_FIMO_alx$alxOccurancy <- apply(DMAR_Iri_FIMO_alx, 1, function(x) {as.numeric(max(x[33],x[34],x[35]))})


#number of regions with motifs
count <-data.frame(c("ALX1", "ALX3","ALX4","GBX2","SOX10","SOX9","ETS1","PAX7","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"),c(2167,2476,1988,1721,6315,2773,2271,1600,4848,3315,1846,3132,3419,1689))
colnames(count) <- c("Motif","DMAR")
count$Motif <- factor(count$Motif, levels = c("ALX1","ALX3","ALX4","GBX2","PAX7","ETS1","SOX10","SOX9","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"))
mypalette <-c(brewer.pal(9,"RdYlGn"),brewer.pal(9,"RdYlBu"),brewer.pal(10,"PuOr"))
p <- ggplot(count, aes(x=Motif, y=DMAR, fill=Motif)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Number of Iri-specific Dynamic Regions with Motifs (26,486 total)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Motif", y = "# of regions")+scale_fill_manual(values = mypalette[c(1:9,16,17,18,25,26,27,28)])+scale_color_manual(values=mypalette[c(1:9,16,17,18,25,26,27,28)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p


#annotation of regions with motifs
Annotation_merged <-Reduce(function(x, y) merge(x, y,by = "Var1", all.x = T), list(as.data.frame(table(ALX1$Annotation)),
as.data.frame(table(ALX3$Annotation)),
as.data.frame(table(ALX4$Annotation)),
as.data.frame(table(GBX2$Annotation)),
as.data.frame(table(SOX10$Annotation)),
as.data.frame(table(SOX9$Annotation)),
as.data.frame(table(ETS1$Annotation)),
as.data.frame(table(PAX7$Annotation)),
as.data.frame(table(TFEC$Annotation)),
as.data.frame(table(SP3$Annotation)),
as.data.frame(table(TFAP2A$Annotation)),
as.data.frame(table(JUNB$Annotation)),
as.data.frame(table(NR4A1$Annotation)),
as.data.frame(table(BHLHE41$Annotation))))

Annotation_merged[is.na(Annotation_merged)]<-0
colnames(Annotation_merged) <- c("Annotation","ALX1","ALX3","ALX4","GBX2","SOX10","SOX9","ETS1","PAX7","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41")
str(Annotation_merged)

Annotation_merged$ALX1 <-  Annotation_merged$ALX1/sum(Annotation_merged$ALX1)*100
Annotation_merged$ALX3 <-  Annotation_merged$ALX3/sum(Annotation_merged$ALX3)*100
Annotation_merged$ALX4 <-  Annotation_merged$ALX4/sum(Annotation_merged$ALX4)*100
Annotation_merged$GBX2 <-  Annotation_merged$GBX2/sum(Annotation_merged$GBX2)*100
Annotation_merged$ETS1 <-  Annotation_merged$ETS1/sum(Annotation_merged$ETS1)*100
Annotation_merged$SOX10 <-  Annotation_merged$SOX10/sum(Annotation_merged$SOX10)*100
Annotation_merged$SOX9 <-  Annotation_merged$SOX9/sum(Annotation_merged$SOX9)*100
Annotation_merged$PAX7 <-  Annotation_merged$PAX7/sum(Annotation_merged$PAX7)*100
Annotation_merged$TFEC <-  Annotation_merged$TFEC/sum(Annotation_merged$TFEC)*100
Annotation_merged$SP3 <-  Annotation_merged$SP3/sum(Annotation_merged$SP3)*100
Annotation_merged$TFAP2 <-  Annotation_merged$TFAP2/sum(Annotation_merged$TFAP2)*100
Annotation_merged$JUNB <-  Annotation_merged$JUNB/sum(Annotation_merged$JUNB)*100
Annotation_merged$NR4A1 <-  Annotation_merged$NR4A1/sum(Annotation_merged$NR4A1)*100
Annotation_merged$BHLHE41 <-  Annotation_merged$BHLHE41/sum(Annotation_merged$BHLHE41)*100

Amelt <- melt(Annotation_merged,id.vars = "Annotation")
colnames(Amelt)[2]<- "Type"
Amelt$Type <- factor(Amelt$Type, levels = c("ALX1","ALX3","ALX4","GBX2","PAX7","ETS1","SOX10","SOX9","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"))
Amelt$Annotation <- factor(Amelt$Annotation, levels = c("non-coding","Intergenic","promoter-TSS","5' UTR","exon","intron","3' UTR","TTS"))

mypalette <-c(brewer.pal(9,"RdYlGn"),brewer.pal(9,"RdYlBu"),brewer.pal(10,"PuOr"))
p <- ggplot(Amelt, aes(x=Annotation, y=value, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Annotation of Regions with Motifs")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Annotation", y = "Percent")+scale_fill_manual(values = mypalette[c(1:9,16,17,18,25,26,27,28)])+scale_color_manual(values=mypalette[c(1:9,16,17,18,25,26,27,28)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p



#Pull out genes turned on in Iridophores (>5rpkm) and find closest DEG to regions w/ motif#
Iri_on_genes_all <- DEG2[DEG2$s24vIri_log2_change <0 & DEG2$Iri > 5,]
Iri_on_genes_all[order(-Iri_on_genes_all$Iri),] #look at top differentially expressed genes

colnames(prom) <- c("chr","start","end","gene","chrompos")
Iri_on_genes_all_prom <- merge(Iri_on_genes_all,prom, by = "gene", all.x = T)
Iri_on_genes_all_prom[is.na(Iri_on_genes_all_prom)] <-0

write.table(Iri_on_genes_all_prom[,c(16:18,1)],"Iri_allDEGs_on_promoter_locations.bed", sep = "\t", col.names = F, row.names = F, quote =F)

for i in *motif.bed; do sort -k1,1 -k2,2n $i > ${i/motif.bed/motif.sorted.bed};done
rm *motif.bed
sort -k1,1 -k2,2n Iri_allDEGs_on_promoter_locations.bed > Iri_allDEGs_on_promoter_locations.sorted.bed
rm Iri_allDEGs_on_promoter_locations.bed
for i in *motif.sorted.bed; do bedtools closest -d -a $i -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/Iri_allDEGs_on_promoter_locations.sorted.bed > ${i/motif.sorted.bed/motif.closestDEG.bed};done


ALX1_closest <- read.table("DMAR_Iri_Accesible_wALX1_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
ALX3_closest <- read.table("DMAR_Iri_Accesible_wALX3_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
ALX4_closest <- read.table("DMAR_Iri_Accesible_wALX4_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
GBX2_closest <- read.table("DMAR_Iri_Accesible_wGBX2_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
SOX10_closest <- read.table("DMAR_Iri_Accesible_wSOX10_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
SOX9_closest <- read.table("DMAR_Iri_Accesible_wSOX9_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
ETS1_closest <- read.table("DMAR_Iri_Accesible_wETS1_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
PAX7_closest <- read.table("DMAR_Iri_Accesible_wPAX7_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
TFEC_closest <- read.table("DMAR_Iri_Accesible_wTFEC_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
SP3_closest <- read.table("DMAR_Iri_Accesible_wSP3motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
TFAP2_closest <- read.table("DMAR_Iri_Accesible_wTFAP2_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
JUNB_closest <- read.table("DMAR_Iri_Accesible_wJUNB_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
NR4A1_closest <- read.table("DMAR_Iri_Accesible_wNR4A1_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)
BHLHE41_closest <- read.table("DMAR_Iri_Accesible_wBHLHE41_motif.closestDEG.bed", sep = "\t", quote = "", header = F, stringsAsFactors =F)



#distance between DEG and region
dist <-data.frame(c(rep("ALX1",nrow(ALX1_closest)), rep("ALX3",nrow(ALX3_closest)),rep("ALX4",nrow(ALX4_closest)),rep("GBX2",nrow(GBX2_closest)),rep("SOX10",nrow(SOX10_closest)),rep("SOX9",nrow(SOX9_closest)),rep("ETS1",nrow(ETS1_closest)),rep("PAX7",nrow(PAX7_closest)),rep("TFEC",nrow(TFEC_closest)),rep("SP3",nrow(SP3_closest)),rep("TFAP2",nrow(TFAP2_closest)),rep("JUNB",nrow(JUNB_closest)),rep("NR4A1",nrow(NR4A1_closest)),rep("BHLHE41",nrow(BHLHE41_closest))),c(ALX1_closest$V8,ALX3_closest$V8,ALX4_closest$V8,GBX2_closest$V8,SOX10_closest$V8,SOX9_closest$V8,ETS1_closest$V8,PAX7_closest$V8,TFEC_closest$V8,SP3_closest$V8,TFAP2_closest$V8,JUNB_closest$V8,NR4A1_closest$V8,BHLHE41_closest$V8))
colnames(dist) <- c("Motif","Distance")
dist$Motif <- factor(dist$Motif, levels = c("ALX1","ALX3","ALX4","GBX2","PAX7","ETS1","SOX10","SOX9","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"))
dist$Distance_kb <- dist$Distance/1000
mypalette <-c(brewer.pal(9,"RdYlGn"),brewer.pal(9,"RdYlBu"),brewer.pal(10,"PuOr"))

p <- ggplot(dist[dist$Distance_kb <= 5000,], aes(x=Distance_kb, color = Motif)) + geom_histogram(aes(y=..density..),alpha = 0.00,position="identity",binwidth = 100,size=1)+ggtitle("Distribution of Distance between\nDynamic Region and closest DEG")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Distance (kb)", y = "Density")+scale_fill_manual(values = c(mypalette[c(1:4)],"#FFE500",mypalette[c(6:9,16,17,18,25,26,27,28)]))+scale_color_manual(values=c(mypalette[c(1:4)],"#FFE500",mypalette[c(6:9,16,17,18,25,26,27,28)]))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
pp <- ggplot(dist[dist$Distance_kb <= 1000,], aes(x=Distance_kb, color = Motif)) + geom_histogram(aes(y=..density..),alpha = 0.00,position="identity",binwidth = 20,size=1)+ggtitle("Distribution of Distance between\nDynamic Region and closest DEG")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "Distance (kb)", y = "Density")+scale_fill_manual(values = c(mypalette[c(1:4)],"#FFE500",mypalette[c(6:9,16,17,18,25,26,27,28)]))+scale_color_manual(values=c(mypalette[c(1:4)],"#FFE500",mypalette[c(6:9,16,17,18,25,26,27,28)]))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
multiplot(p,pp,cols = 1)

#number of DEGs withint 100kb of region with motif
count <-data.frame(c("ALX1", "ALX3","ALX4","GBX2","SOX10","SOX9","ETS1","PAX7","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"),c(length(unique(ALX1_closest[ALX1_closest$V8 <=100000,]$V7)),length(unique(ALX3_closest[ALX3_closest$V8 <=100000,]$V7)),length(unique(ALX4_closest[ALX4_closest$V8 <=100000,]$V7)),length(unique(GBX2_closest[GBX2_closest$V8 <=100000,]$V7)),length(unique(SOX10_closest[SOX10_closest$V8 <=100000,]$V7)),length(unique(SOX9_closest[SOX9_closest$V8 <=100000,]$V7)),length(unique(ETS1_closest[ETS1_closest$V8 <=100000,]$V7)),length(unique(PAX7_closest[PAX7_closest$V8 <=100000,]$V7)),length(unique(TFEC_closest[TFEC_closest$V8 <=100000,]$V7)),length(unique(SP3_closest[SP3_closest$V8 <=100000,]$V7)),length(unique(TFAP2_closest[TFAP2_closest$V8 <=100000,]$V7)),length(unique(JUNB_closest[JUNB_closest$V8 <=100000,]$V7)),length(unique(NR4A1_closest[NR4A1_closest$V8 <=100000,]$V7)),length(unique(BHLHE41_closest[BHLHE41_closest$V8 <=100000,]$V7))))
colnames(count) <- c("Motif","DMAR")
count$Motif <- factor(count$Motif, levels = c("ALX1","ALX3","ALX4","GBX2","PAX7","ETS1","SOX10","SOX9","TFEC","SP3","TFAP2","JUNB","NR4A1","BHLHE41"))
mypalette <-c(brewer.pal(9,"RdYlGn"),brewer.pal(9,"RdYlBu"),brewer.pal(10,"PuOr"))
p <- ggplot(count, aes(x=Motif, y=DMAR, fill=Motif)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Number of DEGS within 100kb of Iri-specific Dynamic Regions with Motifs")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Motif", y = "# of DEGS")+scale_fill_manual(values = mypalette[c(1:9,16,17,18,25,26,27,28)])+scale_color_manual(values=mypalette[c(1:9,16,17,18,25,26,27,28)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p

###Pull out DEGs that are within 100kb of motif region####
write.table(ALX1_closest[ALX1_closest$V8 <=100000,7], "DEGs_within100kb_RegionswALX1_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ALX3_closest[ALX3_closest$V8 <=100000,7], "DEGs_within100kb_RegionswALX3_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ALX4_closest[ALX4_closest$V8 <=100000,7], "DEGs_within100kb_RegionswALX4_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(GBX2_closest[GBX2_closest$V8 <=100000,7], "DEGs_within100kb_RegionswGBX2_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SOX10_closest[SOX10_closest$V8 <=100000,7], "DEGs_within100kb_RegionswSOX10_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SOX9_closest[SOX9_closest$V8 <=100000,7], "DEGs_within100kb_RegionswSOX9_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(ETS1_closest[ETS1_closest$V8 <=100000,7], "DEGs_within100kb_RegionswETS1_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(PAX7_closest[PAX7_closest$V8 <=100000,7], "DEGs_within100kb_RegionswPAX7_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(TFEC_closest[TFEC_closest$V8 <=100000,7], "DEGs_within100kb_RegionswTFEC_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(SP3_closest[SP3_closest$V8 <=100000,7], "DEGs_within100kb_RegionswSP3motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(TFAP2_closest[TFAP2_closest$V8 <=100000,7], "DEGs_within100kb_RegionswTFAP2_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(JUNB_closest[JUNB_closest$V8 <=100000,7], "DEGs_within100kb_RegionswJUNB_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(NR4A1_closest[NR4A1_closest$V8 <=100000,7], "DEGs_within100kb_RegionswNR4A1_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(BHLHE41_closest[BHLHE41_closest$V8 <=100000,7], "DEGs_within100kb_RegionswBHLHE41_motif.txt", sep = "\t", quote = F, col.names = F, row.names =F)

write.table(Random1_closest[Random1_closest$V8 <=100000,7], "DEGs_within100kb_RandomRegions1.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(Random2_closest[Random2_closest$V8 <=100000,7], "DEGs_within100kb_RandomRegions2.txt", sep = "\t", quote = F, col.names = F, row.names =F)
write.table(Random3_closest[Random3_closest$V8 <=100000,7], "DEGs_within100kb_RandomRegions3.txt", sep = "\t", quote = F, col.names = F, row.names =F)





##Identify accessible DMARs near Iri-specific up-regulated DEGs (50kb distance from promoter)
setwd("Combined_analysis/DMAR_DEG_ALX")
Iri_on_genes_all_prom$start50kb <- Iri_on_genes_all_prom$start-50000
Iri_on_genes_all_prom$end50kb <- Iri_on_genes_all_prom$end+50000
Iri_on_genes_all_prom[Iri_on_genes_all_prom$start50kb < 0,]$start50kb <-0

write.table(Iri_on_genes_all_prom[,c(16,22,23,1,10)],"Iri_DEGs_on_promoter_50kb_extended.bed", sep = "\t", col.names = F, row.names = F, quote =F)

bedtools intersect -wao -a Iri_DEGs_on_promoter_50kb_extended.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/accessibleDMAR_Iri_wFIMOmotif_noheader.txt > Iri_DEGs_on_promoter_50kb_extended_wDMAR_FIMO.txt

Iri_DEG_50kb_DMAR_FIMO <- read.table("Combined_analysis/DMAR_DEG_ALX/Iri_DEGs_on_promoter_50kb_extended_wDMAR_FIMO.txt", header = F, stringsAsFactors = F, sep = "\t", quote = "")
Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO == "."] <- 0

colnames(Iri_DEG_50kb_DMAR_FIMO) <- c("chr","start50kb","end50kb","gene","genename",colnames(DMAR_Iri_FIMO),"overlap")
colnames(Iri_DEG_50kb_DMAR_FIMO)[6:9] <- c("chr","start","end","chrompos")

Iri_DEG_50kb_DMAR_FIMO[,c(38:57)] <- sapply(Iri_DEG_50kb_DMAR_FIMO[,c(38:57)],as.numeric)


Small_mlc_biosyn <- c("ENSDARG00000001873","ENSDARG00000004979","ENSDARG00000005897","ENSDARG00000006031","ENSDARG00000007366","ENSDARG00000008310","ENSDARG00000010267","ENSDARG00000010862","ENSDARG00000012076","ENSDARG00000013721","ENSDARG00000014866","ENSDARG00000015623","ENSDARG00000016375","ENSDARG00000016733","ENSDARG00000017444","ENSDARG00000019613","ENSDARG00000020956","ENSDARG00000021154","ENSDARG00000023176","ENSDARG00000023287","ENSDARG00000025012","ENSDARG00000031616","ENSDARG00000032206","ENSDARG00000033361","ENSDARG00000033413","ENSDARG00000033594","ENSDARG00000035872","ENSDARG00000036239","ENSDARG00000036833","ENSDARG00000037551","ENSDARG00000038297","ENSDARG00000038865","ENSDARG00000039429","ENSDARG00000039914","ENSDARG00000040314","ENSDARG00000040492","ENSDARG00000040869","ENSDARG00000042221","ENSDARG00000042856","ENSDARG00000043457")
nitrogen_compound_cat <- c("ENSDARG00000004979","ENSDARG00000005897","ENSDARG00000007366","ENSDARG00000008310","ENSDARG00000010267","ENSDARG00000012076","ENSDARG00000013721","ENSDARG00000014866","ENSDARG00000015623","ENSDARG00000017444","ENSDARG00000019613","ENSDARG00000021154","ENSDARG00000023176","ENSDARG00000025012","ENSDARG00000031616","ENSDARG00000033413","ENSDARG00000033594","ENSDARG00000036833","ENSDARG00000037551","ENSDARG00000038297","ENSDARG00000039429","ENSDARG00000039914","ENSDARG00000040492","ENSDARG00000040869","ENSDARG00000042221","ENSDARG00000042856","ENSDARG00000043457","ENSDARG00000001953","ENSDARG00000002986","ENSDARG00000003311","ENSDARG00000003961","ENSDARG00000005034","ENSDARG00000005464","ENSDARG00000006491","ENSDARG00000007494","ENSDARG00000008153","ENSDARG00000008931","ENSDARG00000009612","ENSDARG00000010130","ENSDARG00000011506","ENSDARG00000012355","ENSDARG00000014572","ENSDARG00000015524","ENSDARG00000016835","ENSDARG00000017049","ENSDARG00000018329","ENSDARG00000019207","ENSDARG00000023448","ENSDARG00000025094","ENSDARG00000028628","ENSDARG00000029075","ENSDARG00000029621","ENSDARG00000030340","ENSDARG00000030598","ENSDARG00000032116","ENSDARG00000033539","ENSDARG00000033666","ENSDARG00000034262","ENSDARG00000035555","ENSDARG00000036282","ENSDARG00000037064","ENSDARG00000037101","ENSDARG00000037191","ENSDARG00000038293","ENSDARG00000038524","ENSDARG00000038533","ENSDARG00000042124","ENSDARG00000042562")
lipid_metabolic <- c("ENSDARG00000004979","ENSDARG00000012076","ENSDARG00000033413","ENSDARG00000037551","ENSDARG00000005034","ENSDARG00000006491","ENSDARG00000008153","ENSDARG00000010130","ENSDARG00000011506","ENSDARG00000014572","ENSDARG00000019207","ENSDARG00000029621","ENSDARG00000033666","ENSDARG00000034262","ENSDARG00000035555","ENSDARG00000037064","ENSDARG00000037101","ENSDARG00000037191","ENSDARG00000038524","ENSDARG00000023287","ENSDARG00000035872","ENSDARG00000038865","ENSDARG00000003584","ENSDARG00000003635","ENSDARG00000007108","ENSDARG00000012194","ENSDARG00000012829","ENSDARG00000013711","ENSDARG00000017811","ENSDARG00000017882","ENSDARG00000018361","ENSDARG00000019228","ENSDARG00000020509","ENSDARG00000023820","ENSDARG00000025555","ENSDARG00000027469","ENSDARG00000032816","ENSDARG00000035544","ENSDARG00000036415","ENSDARG00000036636","ENSDARG00000040248","ENSDARG00000040469","ENSDARG00000042780")
purine_nucle_meta <- c("ENSDARG00000037551","ENSDARG00000037191","ENSDARG00000010267","ENSDARG00000015623","ENSDARG00000039429","ENSDARG00000040492","ENSDARG00000042221","ENSDARG00000002986","ENSDARG00000009612","ENSDARG00000015524","ENSDARG00000017049","ENSDARG00000018329","ENSDARG00000030340","ENSDARG00000033539","ENSDARG00000038293","ENSDARG00000042124","ENSDARG00000042562","ENSDARG00000010862","ENSDARG00000013561","ENSDARG00000032868")
guanine_synthesis <- c("ENSDARG00000011683","ENSDARG00000014866","ENSDARG00000079848","ENSDARG00000029524","ENSDARG00000016706","ENSDARG00000033539","ENSDARG00000101089","ENSDARG00000003750","ENSDARG00000004517","ENSDARG00000099222","ENSDARG00000056640","ENSDARG00000039452","ENSDARG00000012987","ENSDARG00000013561","ENSDARG00000018178","ENSDARG00000028000","ENSDARG00000012801","ENSDARG00000057661","ENSDARG00000025012","ENSDARG00000040988","ENSDARG00000005897","ENSDARG00000043457","ENSDARG00000054191","ENSDARG00000001873","ENSDARG00000104414","ENSDARG00000052816","ENSDARG00000040314","ENSDARG00000018266","ENSDARG00000040492","ENSDARG00000098646","ENSDARG00000016733","ENSDARG00000100003","ENSDARG00000099776","ENSDARG00000042336","ENSDARG00000019702","ENSDARG00000103826","ENSDARG00000005423","ENSDARG00000039007","ENSDARG00000020956","ENSDARG00000099730","ENSDARG00000099860","ENSDARG00000019644","ENSDARG00000071076","ENSDARG00000099079","ENSDARG00000055652","ENSDARG00000017772","ENSDARG00000103849","ENSDARG00000017049","ENSDARG00000104487","ENSDARG00000075132")


Iri_DEG_50kb_DMAR_FIMO_Small_mlc <- Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene %in% Small_mlc_biosyn,]

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Small Molecule Biosynthesis GO DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#ffb0b0","#800000")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = T, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2


Iri_DEG_50kb_DMAR_FIMO_nitrogen <- Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene %in% nitrogen_compound_cat,]

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Cellular Nitrogen Catabolic GO DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#ffb0b0","#800000")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = T, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2


Iri_DEG_50kb_DMAR_FIMO_lipid <- Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene %in% lipid_metabolic,]

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Lipid Metabolic GO DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#ffb0b0","#800000")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = T, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2


Iri_DEG_50kb_DMAR_FIMO_purine <- Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene %in% purine_nucle_meta,]

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Purine Metabolic GO DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#ffb0b0","#800000")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = T, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2



Iri_DEG_50kb_DMAR_FIMO_guanine <- Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene %in% guanine_synthesis,]

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Guanine Synthesis Cycle DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#ffb0b0","#800000")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = T, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2


##number of accessible DMARs near DEGs vs non-DEGS
Mel_on_genes_all <- DEG2[DEG2$s24vMel_log2_change <0 & DEG2$Mel > 5,]
nonDEGs <- DEG2[DEG2$s24vMel_log2_change == 0 & DEG2$s24vIri_log2_change == 0,]

Mel_on_genes_all_prom <- merge(Mel_on_genes_all,prom, by = "gene", all.x = T)
Mel_on_genes_all_prom[is.na(Mel_on_genes_all_prom)] <-0
Mel_on_genes_all_prom$start50kb <- Mel_on_genes_all_prom$start-50000
Mel_on_genes_all_prom$end50kb <- Mel_on_genes_all_prom$end+50000
Mel_on_genes_all_prom[Mel_on_genes_all_prom$start50kb < 0,]$start50kb <-0

write.table(Mel_on_genes_all_prom[,c(16,22,23,1,10)],"Mel_DEGs_on_promoter_50kb_extended.bed", sep = "\t", col.names = F, row.names = F, quote =F)

bedtools intersect -wao -a Mel_DEGs_on_promoter_50kb_extended.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/DMAR_Mel_Accesible_for_FIMO.bed > Mel_DEGs_on_promoter_50kb_extended_wDMAR_FIMO.txt

Mel_DEG_50kb_DMAR_FIMO <- read.delim("Combined_analysis/DMAR_DEG_ALX/Mel_DEGs_on_promoter_50kb_extended_wDMAR_FIMO.txt", header = F, stringsAsFactors = F, sep = "\t")
Mel_DEG_50kb_DMAR_FIMO[Mel_DEG_50kb_DMAR_FIMO == "."] <- 0

Mel_DEG_50kb_DMAR_FIMO$count <- 1


nonDEGs_prom <- merge(nonDEGs,prom, by = "gene", all.x = T)
nonDEGs_prom[is.na(nonDEGs_prom)] <-0
nonDEGs_prom <- merge(nonDEGs,prom, by = "gene", all.x = T)
nonDEGs_prom[is.na(nonDEGs_prom)] <-0
nonDEGs_prom$start50kb <- nonDEGs_prom$start-50000
nonDEGs_prom$end50kb <- nonDEGs_prom$end+50000
nonDEGs_prom[nonDEGs_prom$start50kb < 0,]$start50kb <-0

write.table(nonDEGs_prom[,c(16,22,23,1,10)],"nonDEGs_promoter_50kb_extended.bed", sep = "\t", col.names = F, row.names = F, quote =F)

bedtools intersect -wao -a nonDEGs_promoter_50kb_extended.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/DMAR_Mel_Accesible_for_FIMO.bed >nonDEGs_promoter_50kb_extended_wDMAR_MEL_FIMO.txt

bedtools intersect -wao -a nonDEGs_promoter_50kb_extended.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/DMAR_Iri_Accesible_for_FIMO.bed >nonDEGs_promoter_50kb_extended_wDMAR_Iri_FIMO.txt

Mel_nonDEG_50kb_DMAR_FIMO <- read.delim("/scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/DMAR_DEG_ALX/nonDEGs_promoter_50kb_extended_wDMAR_MEL_FIMO.txt", header = F, stringsAsFactors = F, sep = "\t")
Mel_nonDEG_50kb_DMAR_FIMO[Mel_nonDEG_50kb_DMAR_FIMO == "."] <- 0

Mel_nonDEG_50kb_DMAR_FIMO$count <- 1

library(plyr)
Mel_DEG_50kb_DMAR_FIMO_agg <- ddply(Mel_DEG_50kb_DMAR_FIMO[Mel_DEG_50kb_DMAR_FIMO$V4 != "",],.(V4),summarise,DMARs=sum(count))
Mel_nonDEG_50kb_DMAR_FIMO_agg <- ddply(Mel_nonDEG_50kb_DMAR_FIMO[Mel_nonDEG_50kb_DMAR_FIMO$V4 != "",],.(V4),summarise,DMARs=sum(count))

df_mel <- data.frame(c(rep("50kb",nrow(Mel_DEG_50kb_DMAR_FIMO_agg)),rep("non50kb",nrow(Mel_nonDEG_50kb_DMAR_FIMO_agg)),c(Mel_DEG_50kb_DMAR_FIMO_agg$DMARs,Mel_nonDEG_50kb_DMAR_FIMO_agg$DMARs))
colnames(df_mel) <- c("Type","DMAR")
df_mel[is.na(df_mel)] <- "."
mypalette2 <- brewer.pal(12,"Paired")
b <-ggplot(df_mel[df_mel$DMAR != ".",], aes(x=Type,y=as.numeric(DMAR)))+ geom_violin(aes(fill = Type, scale = 2), trim = F) +geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE,width=0.1)+ggtitle(paste("Distribution of DMAR count around DEGs in Melanophores"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Type", y = "# of DMARs")+scale_fill_manual(values = mypalette2[c(2,4,1,3)])+scale_color_manual(values=mypalette2[c(2,4,1,3)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b


#Pull out DMARs of each GO term and run homer enrichment
write.table(Iri_DEG_50kb_DMAR_FIMO_Small_mlc[,c(6:8)],"Iri_DEG_50kb_DMAR_FIMO_Small_mlc_forHomer.bed", sep = '\t', col.names = F, row.names = F, quote = F)
write.table(Iri_DEG_50kb_DMAR_FIMO_nitrogen[,c(6:8)],"Iri_DEG_50kb_DMAR_FIMO_nitrogen_forHomer.bed", sep = '\t', col.names = F, row.names = F, quote = F)
write.table(Iri_DEG_50kb_DMAR_FIMO_lipid[,c(6:8)],"Iri_DEG_50kb_DMAR_FIMO_lipid_forHomer.bed", sep = '\t', col.names = F, row.names = F, quote = F)
write.table(Iri_DEG_50kb_DMAR_FIMO_purine[,c(6:8)],"Iri_DEG_50kb_DMAR_FIMO_purine_forHomer.bed", sep = '\t', col.names = F, row.names = F, quote = F)

for i in *_forHomer.bed ; do findMotifsGenome.pl $i danRer10 ./HOMER/${i/.bed/.HOMER} -size given;done ###mostly tfec comes out as enrichment

####Add annotation, type and DEG and Gene for GO terms 50kb
Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif <- merge(Iri_DEG_50kb_DMAR_FIMO_guanine[Iri_DEG_50kb_DMAR_FIMO_guanine$ALX1 > 0 | Iri_DEG_50kb_DMAR_FIMO_guanine$ALX3 > 0 |Iri_DEG_50kb_DMAR_FIMO_guanine$ALX4 > 0 |Iri_DEG_50kb_DMAR_FIMO_guanine$GBX2 > 0 |Iri_DEG_50kb_DMAR_FIMO_guanine$SOX10 > 0 |Iri_DEG_50kb_DMAR_FIMO_guanine$ETS1 > 0 |Iri_DEG_50kb_DMAR_FIMO_guanine$TFEC > 0,],DEG3[,c(1:15)], by ="gene",all.x=T)

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Guanine Synthesis Cycle DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#ffb3ba","#ffdfba","#ffffba","#baffc9","#bae1ff"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_guanine_wMotif[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4

#purine
Iri_DEG_50kb_DMAR_FIMO_purine_wMotif <- merge(Iri_DEG_50kb_DMAR_FIMO_purine[Iri_DEG_50kb_DMAR_FIMO_purine$ALX1 > 0 | Iri_DEG_50kb_DMAR_FIMO_purine$ALX3 > 0 |Iri_DEG_50kb_DMAR_FIMO_purine$ALX4 > 0 |Iri_DEG_50kb_DMAR_FIMO_purine$GBX2 > 0 |Iri_DEG_50kb_DMAR_FIMO_purine$SOX10 > 0 |Iri_DEG_50kb_DMAR_FIMO_purine$ETS1 > 0 |Iri_DEG_50kb_DMAR_FIMO_purine$TFEC > 0,],DEG3[,c(1:15)], by ="gene",all.x=T)

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Purine Metabolism DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#ffd4e5","#ffb3ba","#ffdfba","#ffffba","#baffc9"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_purine_wMotif[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4

#lipid
Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif <- merge(Iri_DEG_50kb_DMAR_FIMO_lipid[Iri_DEG_50kb_DMAR_FIMO_lipid$ALX1 > 0 | Iri_DEG_50kb_DMAR_FIMO_lipid$ALX3 > 0 |Iri_DEG_50kb_DMAR_FIMO_lipid$ALX4 > 0 |Iri_DEG_50kb_DMAR_FIMO_lipid$GBX2 > 0 |Iri_DEG_50kb_DMAR_FIMO_lipid$SOX10 > 0 |Iri_DEG_50kb_DMAR_FIMO_lipid$ETS1 > 0 |Iri_DEG_50kb_DMAR_FIMO_lipid$TFEC > 0,],DEG3[,c(1:15)], by ="gene",all.x=T)

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Lipid Metabolism DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#ffb3ba","#ffdfba","#ffffba","#baffc9","#bae1ff"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_lipid_wMotif[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,100,150,200), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4


#Nitrogen
Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif <- merge(Iri_DEG_50kb_DMAR_FIMO_nitrogen[Iri_DEG_50kb_DMAR_FIMO_nitrogen$ALX1 > 0 | Iri_DEG_50kb_DMAR_FIMO_nitrogen$ALX3 > 0 |Iri_DEG_50kb_DMAR_FIMO_nitrogen$ALX4 > 0 |Iri_DEG_50kb_DMAR_FIMO_nitrogen$GBX2 > 0 |Iri_DEG_50kb_DMAR_FIMO_nitrogen$SOX10 > 0 |Iri_DEG_50kb_DMAR_FIMO_nitrogen$ETS1 > 0 |Iri_DEG_50kb_DMAR_FIMO_nitrogen$TFEC > 0,],DEG3[,c(1:15)], by ="gene",all.x=T)

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Cellular Nitrogen compound catabolic DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#c9c9ff","#ffd4e5","#ffb3ba","#ffdfba","#ffffba","#baffc9","#bae1ff"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_nitrogen_wMotif[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4

#Small molecule metabolism
Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif <- merge(Iri_DEG_50kb_DMAR_FIMO_Small_mlc[Iri_DEG_50kb_DMAR_FIMO_Small_mlc$ALX1 > 0 | Iri_DEG_50kb_DMAR_FIMO_Small_mlc$ALX3 > 0 |Iri_DEG_50kb_DMAR_FIMO_Small_mlc$ALX4 > 0 |Iri_DEG_50kb_DMAR_FIMO_Small_mlc$GBX2 > 0 |Iri_DEG_50kb_DMAR_FIMO_Small_mlc$SOX10 > 0 |Iri_DEG_50kb_DMAR_FIMO_Small_mlc$ETS1 > 0 |Iri_DEG_50kb_DMAR_FIMO_Small_mlc$TFEC > 0,],DEG3[,c(1:15)], by ="gene",all.x=T)

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,9])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Small Molecule Biosynthesis DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#c9c9ff","#ffd4e5","#ffb3ba","#ffdfba","#ffffba","#baffc9","#bae1ff"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,5])
df.OG2 <- Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_DEG_50kb_DMAR_FIMO_Small_mlc_wMotif[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4


####Identify motifs near known iridophore genes and TFs #####
Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$gene == "ENSDARG00000091086",]$genename <- "alx3"
#combined  gbx2, atic (ENSDARG00000016706), sox10(ENSDARG00000077467),tfec(ENSDARG00000098745),mitfa(ENSDARG00000003732),tfap2a(ENSDARG00000059279),tfap2e(ENSDARG00000008861),foxd3(ENSDARG00000021032), CABZ01069595.1 (ENSDARG00000024431) since not DEG
Iri_MarkerGenes <- prom[prom$gene == "ENSDARG00000003732"|prom$gene =="ENSDARG00000002933"|prom$gene == "ENSDARG00000016706" |prom$gene == "ENSDARG00000077467" | prom$gene == "ENSDARG00000098745" |prom$gene == "ENSDARG00000003732" |prom$gene == "ENSDARG00000059279" |prom$gene == "ENSDARG00000008861" | prom$gene == "ENSDARG00000021032"| prom$gene == "ENSDARG00000024431",]
Iri_MarkerGenes2 <- merge(Iri_MarkerGenes,DEG3[,c(1:15)],by = "gene", all.x = T)
Iri_MarkerGenes2[is.na(Iri_MarkerGenes2)] <- 0
Iri_MarkerGenes2[Iri_MarkerGenes2$gene == "ENSDARG00000021032",]$genename.x <- "foxd3"

Iri_MarkerGenes2$start25kb <- Iri_MarkerGenes2$start-25000
Iri_MarkerGenes2$end25kb <- Iri_MarkerGenes2$end+25000
Iri_MarkerGenes2$start50kb <- Iri_MarkerGenes2$start-50000
Iri_MarkerGenes2$end50kb <- Iri_MarkerGenes2$end+50000
Iri_MarkerGenes2[Iri_MarkerGenes2$start50kb < 0,]$start50kb <-0
Iri_MarkerGenes2[Iri_MarkerGenes2$start25kb < 0,]$start25kb <-0

write.table(Iri_MarkerGenes2[,c(2,22,23,1,14)],"Iri_MarkerGenes_promoter_50kb_extended.bed", sep = "\t", col.names = F, row.names = F, quote =F)

bedtools intersect -wao -a Iri_MarkerGenes_promoter_50kb_extended.bed -b /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/FIMO_motif_scan/accessibleDMAR_Iri_wFIMOmotif_noheader.txt > Iri_MarkerGenes_promoter_50kb_extended_wDMAR_FIMO.txt

Iri_Marker_50kb_DMAR_FIMO <- read.table("/scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/DMAR_DEG_ALX/Iri_MarkerGenes_promoter_50kb_extended_wDMAR_FIMO.txt", header = F, stringsAsFactors = F, sep = "\t", quote = "")

colnames(Iri_Marker_50kb_DMAR_FIMO) <- c("chr","start50kb","end50kb","gene","genename",colnames(DMAR_Iri_FIMO),"overlap")
colnames(Iri_Marker_50kb_DMAR_FIMO)[6:9] <- c("chr","start","end","chrompos")

Iri_Marker_50kb_DMAR_FIMO[Iri_Marker_50kb_DMAR_FIMO == "."]<-0

Iri_Marker_50kb_DMAR_FIMO[,c(38:57)] <- sapply(Iri_Marker_50kb_DMAR_FIMO[,c(38:57)],as.numeric)

Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$genename == "alx1" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx3"| Iri_DEG_50kb_DMAR_FIMO$genename == "alx4a" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx4b" | Iri_DEG_50kb_DMAR_FIMO$genename == "ltk"| Iri_DEG_50kb_DMAR_FIMO$genename == "ednrba" | Iri_DEG_50kb_DMAR_FIMO$genename == "pnp4a",]

Iri_Marker_50kb_DMAR_FIMO_ALL <- rbind(Iri_Marker_50kb_DMAR_FIMO,Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$genename == "gbx2" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx1" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx3"| Iri_DEG_50kb_DMAR_FIMO$genename == "alx4a" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx4b" | Iri_DEG_50kb_DMAR_FIMO$genename == "ltk"| Iri_DEG_50kb_DMAR_FIMO$genename == "ednrba"| Iri_DEG_50kb_DMAR_FIMO$genename == "pnp4a",]
)

Iri_Marker_50kb_DMAR_FIMO_ALL <- merge(Iri_Marker_50kb_DMAR_FIMO_ALL[Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "tfap2a" | Iri_Marker_50kb_DMAR_FIMO_ALL$GBX2 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ALX1 > 0 | Iri_Marker_50kb_DMAR_FIMO_ALL$ALX3 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ALX4 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$SOX10 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ETS1 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$TFEC > 0,],DEG3[,c(0:15)], by = "gene", all.x = T)
Iri_Marker_50kb_DMAR_FIMO_ALL[is.na(Iri_Marker_50kb_DMAR_FIMO_ALL)] <- 0


#Iri Marker gene specific#
Iri_Marker_50kb_DMAR_FIMO_Marker <- Iri_Marker_50kb_DMAR_FIMO_ALL[Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "ltk" | Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "atic" | Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "ednrba" | Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "pnp4a",]
Iri_Marker_50kb_DMAR_FIMO_Marker<- Iri_Marker_50kb_DMAR_FIMO_Marker[order(Iri_Marker_50kb_DMAR_FIMO_Marker$genename),]

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_Marker[,9])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_Marker[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_Marker[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Iri Marker Genes DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_Marker[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_Marker[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_Marker[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#ffb3ba","#ffdfba","#ffffba","#baffc9"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_Marker[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_Marker[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_Marker[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_Marker[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_Marker[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_Marker[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_Marker[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_Marker[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_Marker[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4



#Iri TF specific#
Iri_Marker_50kb_DMAR_FIMO_TF <- Iri_Marker_50kb_DMAR_FIMO_ALL[Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "ltk" & Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "atic" & Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "ednrba",]
Iri_Marker_50kb_DMAR_FIMO_TF<- Iri_Marker_50kb_DMAR_FIMO_TF[order(Iri_Marker_50kb_DMAR_FIMO_TF$genename),]

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,9])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(38,39,40,41,44,46,49)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(38,39,40,41,44,46,49)])
ht2 = Heatmap(df.OG2, column_title = "Iri TFs DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#fffef2","#17c0f2","#0084a1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(35)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("white","#ffd4e5","#ffb3ba","#ffdfba","#ffffba","#baffc9"),
    cluster_rows = F, cluster_columns = F)
ht

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(37)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("white","#d2a8dd","#a6fdfc","#ffe3d6"),
    cluster_rows = F, cluster_columns = F)
ht1

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(64)] #Motif occurrence
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3

name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4

ht2+ht+ht1+ht3+ht4
```




``` {bash}

############# Centrimo motif TF footprinting ##############
################## {Bash script} ##################
#!/bin/bash
# Author: Hyung Joo Lee

## SOFTWARE
ml meme/4.11.2

## INPUT
motif=$( cat jaspar2018_motif_run1.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
genome=/scratch/twlab/hlee/genomes/danRer10/danRer10.fa
db=JASPAR2018_CORE_vertebrates_non-redundant_run1.meme

## OUTPUT
txt=${motif}/${motif}.danRer10.fimo.txt
log=${motif}/${motif}.danRer10.fimo.log

## COMMANDS
# FIMO (MEME)
mkdir -p $motif
fimo --max-stored-scores 10000000 --motif $motif --text --thresh 1e-5 --verbosity 4 $db $genome >$txt 2>$log


#make FIMO for Alx3,
fimo --max-stored-scores 10000000 --text --motif PH0001.1 --thresh 1e-5 --verbosity 4 Alx3_FIMO.meme /bar/genomes/danRer10/danRer10.fa >Alx3.danRer10.fimo.txt 2>Alx3.danRer10.fimo.log

awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$9}' Alx3.danRer10.fimo.txt > Alx3.danRer10.fimo2.txt #make sure to remove the first column!

#2)reprocessed the bed files to make tagmented end location files

## SOFTWARE
ml kentsrc/20170117

## INPUT
tagAlign=$( ls WangT_ATAC_sp7*[04]dpa_rep0.filt.nodup.PE2SE.tn5.tagAlign.gz | sed "${SLURM_ARRAY_TASK_ID}q;d" )
chromSize=danRer10.chrom.sizes

## OUTPUT
base=${tagAlign##*/}
base=${base%.filt.nodup.PE2SE.tn5.tagAlign.gz}
bed_tmp_f=${base}.for.tmp.bed.gz
bed_tmp_r=${base}.rev.tmp.bed.gz
bg_f=${base}.for.bg.gz
bg_r=${base}.rev.bg.gz

## COMMANDS
zcat $tagAlign | awk 'BEGIN{OFS="\t"} $6=="+" {$3=$2+1; print}' | sort -k1,1 -k2,2n | gzip > $bed_tmp_f
bedItemOverlapCount chromSize=$chromSize null $bed_tmp_f | gzip > $bg_f

zcat $tagAlign | awk 'BEGIN{OFS="\t"} $6=="-" {$2=$3-1; print}' | sort -k1,1 -k2,2n | gzip > $bed_tmp_r
bedItemOverlapCount chromSize=$chromSize null $bed_tmp_r | gzip > $bg_r

rm $bed_tmp_f $bed_tmp_r

for i in *bed; do sort -k1,1 -k2,2n $i > ${i/.bed/.sorted.bed};done 
awk '{$3=sprintf("%.0f",$3)}7' Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations
awk '{$3=sprintf("%.0f",$3)}7' Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations
awk '{$3=sprintf("%.0f",$3)}7' Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations
awk '{$3=sprintf("%.0f",$3)}7' Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations
awk '{$3=sprintf("%.0f",$3)}7' Mel_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Mel_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations
awk '{$3=sprintf("%.0f",$3)}7' Iri_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.bed > Iri_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed #correct scientific notations

awk '{$2=sprintf("%.0f",$2)}7' Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations
awk '{$2=sprintf("%.0f",$2)}7' Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations
awk '{$2=sprintf("%.0f",$2)}7' Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations
awk '{$2=sprintf("%.0f",$2)}7' Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations
awk '{$2=sprintf("%.0f",$2)}7' Mel_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Mel_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations
awk '{$2=sprintf("%.0f",$2)}7' Iri_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected.bed > Iri_combined_downsampled.Tshift.fixed.tagmentpositions.sorted.corrected2.bed #correct scientific notations

for i in *combined*sorted.corrected2.bed; do sort -k1,1 -k2,2n $i > ${i/sorted.corrected2.bed/sorted.bed};done 
rm *combined*sorted.corrected.bed
rm *combined*sorted.corrected2.bed
for i in *combined*sorted.bed; do bedGraphToBigWig $i danRer10.chrom.sizes ${i/sorted.bed/bw};done



#3)ran centipede on atac bigwig files 

#SBATCH --mem=20G
#SBATCH --array=1-1638%100
#SBATCH --job-name=centipede

## SOFTWARE
ml R
ml r-modules
#bed_bedGraph.pl=~/bin/bed_bedGraph.pl
#bwtool=~/bin/bwtool
#fitCentipede.R=~/bin/fitCentipede.R

## INPUT
fimo_dir=Combined_analysis/centipede_motif_footprint/8_fimo
fimo=$( ls ${fimo_dir}/*.fimo.txt.gz | sed "${SLURM_ARRAY_TASK_ID}q;d" )
phastCons=danRer10.vertebrate.phastCons8way.bg.gz
atac_sp7ne0dpa_f=Combined_analysis/centipede_motif_footprint/Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw
atac_sp7ne0dpa_r=Combined_analysis/centipede_motif_footprint/Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw

## OUTPUT
motif=${fimo##*/}
motif=${motif%.danRer10.fimo.txt.gz}
motif_site=${motif}/${motif}_sites_danRer10.bed
X_sp7ne0dpa=${motif}/${motif}_sites_centipedeX_ATAC_Iri.txt.gz
Y=${motif}/${motif}_sites_centipedeY_phastCons8way.txt.gz

## COMMANDS
#1 make motif site bed file
mkdir -p $motif
zcat $fimo | awk 'BEGIN{OFS="\t"} NR>1 {print $2,$3-1,$4,NR-1,$6,$5}' > $motif_site

#2 make phastCons file (Y)
Combined_analysis/centipede_motif_footprint/bed_bedGraph.pl $motif_site $phastCons > ${motif}/tmp.bg
paste $motif_site ${motif}/tmp.bg | cut -f1-3,5-6,10 - | gzip >$Y #chr start end lor strand phastCons
rm ${motif}/tmp.bg

#3 and #4 make matrix ATAC signal around motif sites (X)
grep "+" $motif_site > ${motif}/${motif}_sites_danRer10.for.bed
grep "-" $motif_site > ${motif}/${motif}_sites_danRer10.rev.bed

# Determine the length of motif
width=$( head -1 $motif_site | awk '{print $3-$2}' )
gzip $motif_site
left=$(( $width/2 + 101 ))
right=$(( $width/2 + 100 - ($width+1)%2 ))

# 1st data (Iri)
# motifs found on forward strand
Combined_analysis/centipede_motif_footprint/bwtool matrix ${left}:$right <(cut -f1-3 ${motif}/${motif}_sites_danRer10.for.bed) $atac_sp7ne0dpa_f,$atac_sp7ne0dpa_r /dev/stdout |
      sed "s/NA/0/g; s/\.00\t/\t/g; s/\.00$//g" |
      paste ${motif}/${motif}_sites_danRer10.for.bed - > ${motif}/${motif}_sites_centipedeX_ATAC_Iri.for.txt
# motifs found on reverse strand: reverse the order of data
Combined_analysis/centipede_motif_footprint/bwtool matrix ${left}:$right <(cut -f1-3 ${motif}/${motif}_sites_danRer10.rev.bed) $atac_sp7ne0dpa_f,$atac_sp7ne0dpa_r /dev/stdout |
      sed "s/NA/0/g; s/\.00\t/\t/g; s/\.00$//g" > ${motif}/${motif}_sites_centipedeX_ATAC_Iri.rev.txt
awk 'BEGIN{ORS=""} {for(i=NF;i>0;i--) {print $i; if (i==1) print "\n"; else print "\t"}}' ${motif}/${motif}_sites_centipedeX_ATAC_Iri.rev.txt |
      paste ${motif}/${motif}_sites_danRer10.rev.bed - |
      cat ${motif}/${motif}_sites_centipedeX_ATAC_Iri.for.txt - |
      sort -k1,1 -k2,2n | cut -f 7- | gzip > $X_sp7ne0dpa
rm ${motif}/${motif}_sites_centipedeX_ATAC_Iri.for.txt ${motif}/${motif}_sites_centipedeX_ATAC_Iri.rev.txt
./fitCentipede2.R $X_sp7ne0dpa $Y ${motif}/${motif}_sites_ATAC_Iri &>${motif}/fitCentipede.R.${motif}_Iri.log


ht2+ht+ht1+ht3+ht4
```

```{R」
################# DMR, DAR and DMAR correlate with gene expression #####################
################## {R script} ##################
DMAR <- read.table("All_DMAR_Combined_wINFO_NEW_042620.bed",sep = "\t",header =T, stringsAsFactors =F)
DEG <- read.table("/scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/RNA_only/DESeq2_051718/DEGs_combined_samples_p0.01.txt",sep = "\t",header =T, stringsAsFactors =F)

DEG_TPM <- merge(DEG,geneTPM,by="gene")
DEG_TPM <- DEG_TPM[DEG_TPM$s15 >5 | DEG_TPM$s24 >5 | DEG_TPM$Mel >5 | DEG_TPM$Iri >5,]
DEG_TPM<-DEG_TPM[complete.cases(DEG_TPM),]

hypoDMR_s15vs24 <- DMAR[DMAR$DMRs15vs24 >0,]
hyperDMR_s15vs24 <- DMAR[DMAR$DMRs15vs24 <0,]
hypoDMR_s24vsMel <- DMAR[DMAR$DMRs24vMel >0,]
hyperDMR_s24vsMel <- DMAR[DMAR$DMRs24vMel <0,]
hypoDMR_s24vsIri <- DMAR[DMAR$DMRs24vIri >0,]
hyperDMR_s24vsIri <- DMAR[DMAR$DMRs24vIri <0,]

openingDAR_s15vs24 <- DMAR[DMAR$DARs15vs24 >0,]
closingDAR_s15vs24 <- DMAR[DMAR$DARs15vs24 <0,]
openingDAR_s24vsMel <- DMAR[DMAR$DARs24vMel >0,]
closingDAR_s24vsMel <- DMAR[DMAR$DARs24vMel <0,]
openingDAR_s24vsIri <- DMAR[DMAR$DARs24vIri >0,]
closingDAR_s24vsIri <- DMAR[DMAR$DARs24vIri <0,]
openingDAR_s24vsMI <- DMAR[DMAR$DARs24vMel >0 & DMAR$DARs24vIri >0,]
closingDAR_s24vsMI <- DMAR[DMAR$DARs24vMel <0 & DMAR$DARs24vIri <0,]
openingDAR_s24vsMel_specific <- DMAR[DMAR$DARs24vMel >0 & DMAR$DARs24vIri <=0,]
closingDAR_s24vsMel_specific <- DMAR[DMAR$DARs24vMel <0 & DMAR$DARs24vIri >=0,]
openingDAR_s24vsIri_specific <- DMAR[DMAR$DARs24vIri >0 & DMAR$DARs24vMel <=0,]
closingDAR_s24vsIri_specific <- DMAR[DMAR$DARs24vIri <0 & DMAR$DARs24vMel >=0,]


openingDMAR_s15vs24 <- DMAR[DMAR$DMRs15vs24 >0 & DMAR$DARs15vs24 >0,]
closingDMAR_s15vs24 <- DMAR[DMAR$DMRs15vs24 <0 & DMAR$DARs15vs24 <0,]
openingDMAR_s24vsMel <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DARs24vMel >0,]# hypo opening
closingDMAR_s24vsMel <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DARs24vMel <0,]# hypo closing
openingDMAR_s24vsIri <- DMAR[DMAR$DMRs24vIri >0 & DMAR$DARs24vIri >0,]# hypo opening
closingDMAR_s24vsIri <- DMAR[DMAR$DMRs24vIri >0 & DMAR$DARs24vIri <0,]# hypo closing

openingDMAR_s24vsMI <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DARs24vMel >0 & DMAR$DMRs24vIri >0 & DMAR$DARs24vIri >0,]# hypo opening
closingDMAR_s24vsMI <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DARs24vMel <0 & DMAR$DMRs24vIri <0 & DMAR$DARs24vIri <0,]# hypo closing
openingDMAR_s24vsMel_specific <- openingDMAR_s24vsMel[!(openingDMAR_s24vsMel$chrompos %in% openingDMAR_s24vsIri$chrompos),]# hypo opening
closingDMAR_s24vsMel_specific <- closingDMAR_s24vsMel[!(closingDMAR_s24vsMel$chrompos %in% closingDMAR_s24vsIri$chrompos),]# hypo closing
openingDMAR_s24vsIri_specific <- openingDMAR_s24vsIri[!(openingDMAR_s24vsIri$chrompos %in% openingDMAR_s24vsMel$chrompos),]# hypo opening
closingDMAR_s24vsIri_specific <- closingDMAR_s24vsIri[!(closingDMAR_s24vsIri$chrompos %in% closingDMAR_s24vsMel$chrompos),]# hypo closing




Gene.Prom <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t", header =F, stringsAsFactors =F)

Gene.Prom.DEG15v24 <- Gene.Prom[Gene.Prom$V4 %in% DEG_TPM[DEG_TPM$s15v24hpf_log2_change != 0,]$gene,]
Gene.Prom.DEG24vMel <- Gene.Prom[Gene.Prom$V4 %in% DEG_TPM[DEG_TPM$s24vMel_log2_change != 0,]$gene,]
Gene.Prom.DEG24vIri <- Gene.Prom[Gene.Prom$V4 %in% DEG_TPM[DEG_TPM$s24vIri_log2_change != 0,]$gene,]

setwd("DMR_DAR_DMAR_GeneExpression")
write.table(Gene.Prom.DEG15v24,"GRCz10.85.GENE.PROMOTER.DEGs15v24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(Gene.Prom.DEG24vMel,"GRCz10.85.GENE.PROMOTER.DEGs24vMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(Gene.Prom.DEG24vIri,"GRCz10.85.GENE.PROMOTER.DEGs24vIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)

write.table(hypoDMR_s15vs24,"hypoDMR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hyperDMR_s15vs24,"hyperDMR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hypoDMR_s24vsMel,"hypoDMR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hyperDMR_s24vsMel,"hyperDMR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hypoDMR_s24vsIri,"hypoDMR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hyperDMR_s24vsIri,"hyperDMR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)

write.table(openingDAR_s15vs24,"openingDAR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDAR_s15vs24,"closingDAR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDAR_s24vsMel,"openingDAR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDAR_s24vsMel,"closingDAR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDAR_s24vsIri,"openingDAR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDAR_s24vsIri,"closingDAR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)

write.table(openingDMAR_s15vs24,"openingDMAR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDMAR_s15vs24,"closingDMAR_s15vs24.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDMAR_s24vsMel,"openingDMAR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDMAR_s24vsMel,"closingDMAR_s24vsMel.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDMAR_s24vsIri,"openingDMAR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDMAR_s24vsIri,"closingDMAR_s24vsIri.bed", sep = "\t", col.names = F, row.names =F,quote = F)


################## {Bash script} ##################
bedtools closest -d -a hypoDMR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a hyperDMR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a hypoDMR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a hyperDMR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a hypoDMR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a hyperDMR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed

bedtools closest -d -a openingDAR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a closingDAR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a openingDAR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a closingDAR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a openingDAR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a closingDAR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed

bedtools closest -d -a openingDMAR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a closingDMAR_s15vs24.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs15v24.sorted.bed > closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed
bedtools closest -d -a openingDMAR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a closingDMAR_s24vsMel.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a openingDMAR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a closingDMAR_s24vsIri.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed

################## {R script} ##################
hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)

openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)

openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24 <- read.table("closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel <- read.table("closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri <- read.table("closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)

colnames(hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri ) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance") 
hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")

a <- data.frame(hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
a$type <- "hypoDMR"
a$comp <- "s15v24"
b <- data.frame(hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
b$type <- "hyperDMR"
b$comp <- "s15v24"
c<- data.frame(hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
c$type <- "hypoDMR"
c$comp <- "s24vMel"
d<- data.frame(hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
d$type <- "hyperDMR"
d$comp <- "s24vMel"
e<- data.frame(hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
e$type <- "hypoDMR"
e$comp <- "s24vIri"
f<- data.frame(hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
f$type <- "hyperDMR"
f$comp <- "s24vIri"

colnames(openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri ) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")

aa <- data.frame(openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
aa$type <- "openingDAR"
aa$comp <- "s15v24"
bb <- data.frame(closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
bb$type <- "closingDAR"
bb$comp <- "s15v24"
cc<- data.frame(openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
cc$type <- "openingDAR"
cc$comp <- "s24vMel"
dd<- data.frame(closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
dd$type <- "closingDAR"
dd$comp <- "s24vMel"
ee<- data.frame(openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
ee$type <- "openingDAR"
ee$comp <- "s24vIri"
ff<- data.frame(closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
ff$type <- "closingDAR"
ff$comp <- "s24vIri"

colnames(openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp <- merge(closingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24,DEG_TPM,by = "gene")
colnames(openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")#4773
openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")#1635
closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri ) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")#5862
openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")#1191
closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")

aaa <- data.frame(openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
aaa$type <- "openingDMAR"
aaa$comp <- "s15v24"
bbb <- data.frame("0",0,0)
bbb$type <- "closingDMAR"
bbb$comp <- "s15v24"
ccc<- data.frame(openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
ccc$type <- "openingDMAR"
ccc$comp <- "s24vMel"
ddd<- data.frame(closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
ddd$type <- "closingDMAR"
ddd$comp <- "s24vMel"
eee<- data.frame(openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene, openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
eee$type <- "openingDMAR"
eee$comp <- "s24vIri"
fff<- data.frame(closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
fff$type <- "closingDMAR"
fff$comp <- "s24vIri"

colnames(a)[1:3] <-c("gene","geneFC","distance")
colnames(aa)[1:3] <-c("gene","geneFC","distance")
colnames(aaa)[1:3] <-c("gene","geneFC","distance")
colnames(b)[1:3] <-c("gene","geneFC","distance")
colnames(bb)[1:3] <-c("gene","geneFC","distance")
colnames(bbb)[1:3] <-c("gene","geneFC","distance")
colnames(c)[1:3] <-c("gene","geneFC","distance")
colnames(cc)[1:3] <-c("gene","geneFC","distance")
colnames(ccc)[1:3] <-c("gene","geneFC","distance")
colnames(d)[1:3] <-c("gene","geneFC","distance")
colnames(dd)[1:3] <-c("gene","geneFC","distance")
colnames(ddd)[1:3] <-c("gene","geneFC","distance")
colnames(e)[1:3] <-c("gene","geneFC","distance")
colnames(ee)[1:3] <-c("gene","geneFC","distance")
colnames(eee)[1:3] <-c("gene","geneFC","distance")
colnames(f)[1:3] <-c("gene","geneFC","distance")
colnames(ff)[1:3] <-c("gene","geneFC","distance")
colnames(fff)[1:3] <-c("gene","geneFC","distance")


geneFC_DMR_DAR_DMAR<-rbind(a,b,c,d,e,f,aa,bb,cc,dd,ee,ff,aaa,bbb,ccc,ddd,eee,fff)
geneFC_DMR_DAR_DMAR$type <- factor(geneFC_DMR_DAR_DMAR$type, levels = c("hyperDMR","hypoDMR","closingDAR","openingDAR","closingDMAR","openingDMAR"))


pppp<-ggplot(geneFC_DMR_DAR_DMAR[geneFC_DMR_DAR_DMAR$distance < 25000,], aes(x=type, y= -geneFC)) +
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 25kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1:4,9,10)])+scale_color_manual(values=mypalette[c(1:4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
pppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")



ppppp<-ggplot(geneFC_DMR_DAR_DMAR[geneFC_DMR_DAR_DMAR$distance < 50000,], aes(x=type, y=-geneFC)) + #USE THIS FOR PAPER
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 50kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1:4,9,10)])+scale_color_manual(values=mypalette[c(1:4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
ppppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")

##plot solo DMR vs soloDAR vs DMAR ##
hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO <-hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp[!(hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$chrompos %in% openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$chrompos),]
hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO <-hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp #no closingDMAR_s15vs24
hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos %in% openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos %in% closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos %in% openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]
hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos %in% closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]

openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO <-openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp[!(openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$chrompos %in% openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$chrompos),]
closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO <-closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp#no closingDMAR_s15vs24
openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos %in% openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos %in% closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos %in% openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]
closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos %in% closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]

a1 <- data.frame(hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$gene,hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$s15v24hpf_log2_change,hypoDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$distance)
a1$type <- "hypoDMR"
a1$comp <- "s15v24"
b1 <- data.frame(hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$gene,hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$s15v24hpf_log2_change,hyperDMR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$distance)
b1$type <- "hyperDMR"
b1$comp <- "s15v24"
c1<- data.frame(hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$gene,hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$s24vMel_log2_change ,hypoDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$distance)
c1$type <- "hypoDMR"
c1$comp <- "s24vMel"
d1<- data.frame(hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$gene,hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$s24vMel_log2_change,hyperDMR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$distance)
d1$type <- "hyperDMR"
d1$comp <- "s24vMel"
e1<- data.frame(hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$gene,hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$s24vIri_log2_change,hypoDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$distance)
e1$type <- "hypoDMR"
e1$comp <- "s24vIri"
f1<- data.frame(hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$gene,hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$s24vIri_log2_change,hyperDMR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$distance)
f1$type <- "hyperDMR"
f1$comp <- "s24vIri"

aa1 <- data.frame(openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$gene,openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$s15v24hpf_log2_change,openingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$distance)
aa1$type <- "openingDAR"
aa1$comp <- "s15v24"
bb1 <- data.frame(closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$gene,closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$s15v24hpf_log2_change,closingDAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp.SOLO$distance)
bb1$type <- "closingDAR"
bb1$comp <- "s15v24"
cc1<- data.frame(openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$gene,openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$s24vMel_log2_change ,openingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$distance)
cc1$type <- "openingDAR"
cc1$comp <- "s24vMel"
dd1<- data.frame(closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$gene,closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$s24vMel_log2_change ,closingDAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO$distance)
dd1$type <- "closingDAR"
dd1$comp <- "s24vMel"
ee1<- data.frame(openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$gene,openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$s24vIri_log2_change,openingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$distance)
ee1$type <- "openingDAR"
ee1$comp <- "s24vIri"
ff1<- data.frame(closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$gene,closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$s24vIri_log2_change,closingDAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO$distance)
ff1$type <- "closingDAR"
ff1$comp <- "s24vIri"

aaa <- data.frame(openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$gene,openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$s15v24hpf_log2_change,openingDMAR_s15vs24.GENE.PROMOTER.DEGs15v24.GeneExp$distance)
aaa$type <- "openingDMAR"
aaa$comp <- "s15v24"
bbb <- data.frame("0",0,0)
bbb$type <- "closingDMAR"
bbb$comp <- "s15v24"
ccc<- data.frame(openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,openingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
ccc$type <- "openingDMAR"
ccc$comp <- "s24vMel"
ddd<- data.frame(closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change ,closingDMAR_s24vsMel.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
ddd$type <- "closingDMAR"
ddd$comp <- "s24vMel"
eee<- data.frame(openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene, openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,openingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
eee$type <- "openingDMAR"
eee$comp <- "s24vIri"
fff<- data.frame(closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,closingDMAR_s24vsIri.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
fff$type <- "closingDMAR"
fff$comp <- "s24vIri"

colnames(a1)[1:3] <-c("gene","geneFC","distance")
colnames(aa1)[1:3] <-c("gene","geneFC","distance")
colnames(aaa)[1:3] <-c("gene","geneFC","distance")
colnames(b1)[1:3] <-c("gene","geneFC","distance")
colnames(bb1)[1:3] <-c("gene","geneFC","distance")
colnames(bbb)[1:3] <-c("gene","geneFC","distance")
colnames(c1)[1:3] <-c("gene","geneFC","distance")
colnames(cc1)[1:3] <-c("gene","geneFC","distance")
colnames(ccc)[1:3] <-c("gene","geneFC","distance")
colnames(d1)[1:3] <-c("gene","geneFC","distance")
colnames(dd1)[1:3] <-c("gene","geneFC","distance")
colnames(ddd)[1:3] <-c("gene","geneFC","distance")
colnames(e1)[1:3] <-c("gene","geneFC","distance")
colnames(ee1)[1:3] <-c("gene","geneFC","distance")
colnames(eee)[1:3] <-c("gene","geneFC","distance")
colnames(f1)[1:3] <-c("gene","geneFC","distance")
colnames(ff1)[1:3] <-c("gene","geneFC","distance")
colnames(fff)[1:3] <-c("gene","geneFC","distance")


geneFC_DMR_DAR_DMAR2<-rbind(a1,b1,c1,d1,e1,f1,aa1,bb1,cc1,dd1,ee1,ff1,aaa,bbb,ccc,ddd,eee,fff)
geneFC_DMR_DAR_DMAR2$type <- factor(geneFC_DMR_DAR_DMAR2$type, levels = c("hyperDMR","hypoDMR","closingDAR","openingDAR","closingDMAR","openingDMAR"))


pppp<-ggplot(geneFC_DMR_DAR_DMAR2[geneFC_DMR_DAR_DMAR2$distance < 25000,], aes(x=type, y= -geneFC)) +
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 25kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1:4,9,10)])+scale_color_manual(values=mypalette[c(1:4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
pppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")


ppppp<-ggplot(geneFC_DMR_DAR_DMAR2[geneFC_DMR_DAR_DMAR2$distance < 50000,], aes(x=type, y=-geneFC)) +
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 50kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1:4,9,10)])+scale_color_manual(values=mypalette[c(1:4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
ppppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")

##### pigment specific #########
hypoDMR_s24vsMel_specific <- DMAR[DMAR$DMRs24vMel >0 & DMAR$DMRs24vIri <=0,]
hyperDMR_s24vsMel_specific <- DMAR[DMAR$DMRs24vMel <0 & DMAR$DMRs24vIri >=0,]
hypoDMR_s24vsIri_specific <- DMAR[DMAR$DMRs24vIri >0 & DMAR$DMRs24vMel <=0,]
hyperDMR_s24vsIri_specific <- DMAR[DMAR$DMRs24vIri <0 & DMAR$DMRs24vMel >=0,]


openingDAR_s24vsMel_specific <- DMAR[DMAR$DARs24vMel >0 & DMAR$DARs24vIri <=0,]
closingDAR_s24vsMel_specific <- DMAR[DMAR$DARs24vMel <0 & DMAR$DARs24vIri >=0,]
openingDAR_s24vsIri_specific <- DMAR[DMAR$DARs24vIri >0 & DMAR$DARs24vMel <=0,]
closingDAR_s24vsIri_specific <- DMAR[DMAR$DARs24vIri <0 & DMAR$DARs24vMel >=0,]


openingDMAR_s24vsMel_specific <- openingDMAR_s24vsMel[!(openingDMAR_s24vsMel$chrompos %in% openingDMAR_s24vsIri$chrompos),]# hypo opening
closingDMAR_s24vsMel_specific <- closingDMAR_s24vsMel[!(closingDMAR_s24vsMel$chrompos %in% closingDMAR_s24vsIri$chrompos),]# hypo closing
openingDMAR_s24vsIri_specific <- openingDMAR_s24vsIri[!(openingDMAR_s24vsIri$chrompos %in% openingDMAR_s24vsMel$chrompos),]# hypo opening
closingDMAR_s24vsIri_specific <- closingDMAR_s24vsIri[!(closingDMAR_s24vsIri$chrompos %in% closingDMAR_s24vsMel$chrompos),]# hypo closing


write.table(hypoDMR_s24vsMel_specific,"hypoDMR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hyperDMR_s24vsMel_specific,"hyperDMR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hypoDMR_s24vsIri_specific,"hypoDMR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(hyperDMR_s24vsIri_specific,"hyperDMR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)


write.table(openingDAR_s24vsMel_specific,"openingDAR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDAR_s24vsMel_specific,"closingDAR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDAR_s24vsIri_specific,"openingDAR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDAR_s24vsIri_specific,"closingDAR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)

write.table(openingDMAR_s24vsMel_specific,"openingDMAR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDMAR_s24vsMel_specific,"closingDMAR_s24vsMel_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(openingDMAR_s24vsIri_specific,"openingDMAR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)
write.table(closingDMAR_s24vsIri_specific,"closingDMAR_s24vsIri_specific.bed", sep = "\t", col.names = F, row.names =F,quote = F)

################## {Bash script} ##################
bedtools closest -d -a hypoDMR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a hyperDMR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed >hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a hypoDMR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a hyperDMR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed >hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed


bedtools closest -d -a openingDAR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a closingDAR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed >closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a openingDAR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a closingDAR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed >closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed

bedtools closest -d -a openingDMAR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed > openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a closingDMAR_s24vsMel_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vMel.sorted.bed >closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed
bedtools closest -d -a openingDMAR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed > openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed
bedtools closest -d -a closingDMAR_s24vsIri_specific.sorted.bed -b GRCz10.85.GENE.PROMOTER.DEGs24vIri.sorted.bed >closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed

################## {R script} ##################
hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)


openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)

openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel <- read.table("closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)
closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri <- read.table("closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.sorted.bed", sep = "\t", header =F, stringsAsFactors =F)


colnames(hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp <- merge(closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel,DEG_TPM,by = "gene")
colnames(openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")
colnames(closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri) <- c(colnames(DMAR),"chr_prom","start_prom","end_prom","gene","prom_label","distance")
closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp <- merge(closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri,DEG_TPM,by = "gene")


hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel$chrompos %in% openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel$chrompos %in% closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri$chrompos %in% openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]
hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri$chrompos %in% closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]
openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel$chrompos %in% openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp.SOLO <-closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp[!(closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel$chrompos %in% closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$chrompos),]
openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri$chrompos %in% openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]
closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp.SOLO <-closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp[!(closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri$chrompos %in% closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$chrompos),]


a5 <- data.frame(hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,hypoDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
a5$type <- "hypoDMRspecific"
a5$comp <- "s24vMel"
b5 <- data.frame(hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,hyperDMR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
b5$type <- "hyperDMRspecific"
b5$comp <- "s24vMel"
c5 <- data.frame(hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,hypoDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
c5$type <- "hypoDMRspecific"
c5$comp <- "s24vIri"
d5 <- data.frame(hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,hyperDMR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
d5$type <- "hyperDMRspecific"
d5$comp <- "s24vIri"
a3 <- data.frame(openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,openingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
a3$type <- "openingDARspecific"
a3$comp <- "s24vMel"
b3 <- data.frame(closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,closingDAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
b3$type <- "closingDARspecific"
b3$comp <- "s24vMel"
c3 <- data.frame(openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,openingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
c3$type <- "openingDARspecific"
c3$comp <- "s24vIri"
d3 <- data.frame(closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,closingDAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
d3$type <- "closingDARspecific"
d3$comp <- "s24vIri"
a4 <- data.frame(openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,openingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
a4$type <- "openingDMARspecific"
a4$comp <- "s24vMel"
b4 <- data.frame(closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$gene,closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$s24vMel_log2_change,closingDMAR_s24vsMel_specific.GENE.PROMOTER.DEGs24vsMel.GeneExp$distance)
b4$type <- "closingDMARspecific"
b4$comp <- "s24vMel"
c4 <- data.frame(openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,openingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
c4$type <- "openingDMARspecific"
c4$comp <- "s24vIri"
d4 <- data.frame(closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$gene,closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$s24vIri_log2_change,closingDMAR_s24vsIri_specific.GENE.PROMOTER.DEGs24vsIri.GeneExp$distance)
d4$type <- "closingDMARspecific"
d4$comp <- "s24vIri"



colnames(a3)[1:3] <-c("gene","geneFC","distance")
colnames(a4)[1:3] <-c("gene","geneFC","distance")
colnames(b3)[1:3] <-c("gene","geneFC","distance")
colnames(b4)[1:3] <-c("gene","geneFC","distance")
colnames(c3)[1:3] <-c("gene","geneFC","distance")
colnames(c4)[1:3] <-c("gene","geneFC","distance")
colnames(d3)[1:3] <-c("gene","geneFC","distance")
colnames(d4)[1:3] <-c("gene","geneFC","distance")
colnames(a5)[1:3] <-c("gene","geneFC","distance")
colnames(b5)[1:3] <-c("gene","geneFC","distance")
colnames(c5)[1:3] <-c("gene","geneFC","distance")
colnames(d5)[1:3] <-c("gene","geneFC","distance")


geneFC_DMR_DAR_DMAR3<-rbind(geneFC_DMR_DAR_DMAR,a3,b3,c3,d3,a4,b4,c4,d4,a5,b5,c5,d5)
geneFC_DMR_DAR_DMAR3$type <- factor(geneFC_DMR_DAR_DMAR3$type, levels = c("hypoDMR","hyperDMR","closingDAR","openingDAR","closingDMAR","openingDMAR","hypoDMRspecific","hyperDMRspecific","closingDARspecific","openingDARspecific","closingDMARspecific","openingDMARspecific"))


pppp<-ggplot(geneFC_DMR_DAR_DMAR3[geneFC_DMR_DAR_DMAR3$distance < 25000,], aes(x=type, y= -geneFC)) +
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 25kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10,1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,9,10,1,2,3,4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
pppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")


ppppp<-ggplot(geneFC_DMR_DAR_DMAR3[geneFC_DMR_DAR_DMAR3$distance < 50000,], aes(x=type, y=-geneFC)) +
  geom_boxplot(aes(fill=type),position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA)+ggtitle("Closest DEGs FC (within 50kb)")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "gene log2(FC)")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10,1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,9,10,1,2,3,4,9,10)])+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
ppppp+facet_grid(.~comp)+geom_hline(yintercept=0,linetype="dashed",color="red")




### Generate bar plots representing frequency and distribution of iridophore-associated DM/ARs with a particular TF motif ###
# ALX1_ALX3_ALX4_GBX2_TFEC_SOX10_MOTIF_occurrence_in_Iri_DMR_DAR_DMAR.pdf
data <- read.table("MOTIF_occurrence_in_24vsIri_DMAR_table.txt", sep = "\t", header = T, stringsAsFactors=F)
data<-as.data.frame(data)
data$MOTIF<-as.factor(data$MOTIF)
data$MOTIF <- factor(data$MOTIF, levels = c("ALX1","ALX3","ALX4","GBX2","TFEC","SOX10"))
data$DMAR_type <- factor(data$DMAR_type, levels = c("solo_closeDAR","solo_openDAR","solo_hyperDMR","solo_hypoDMR","hypo_closing_DMAR","hypo_opening_DMAR"))

mypalette <- brewer.pal(n = 12, name = "Paired")
#mypalette <-c("#740001","#ae0001","#e35d6a","#ffdfba","#ffffba","#baffc9","#bae1ff","#428bca","#d896ff")
p1 <- ggplot(data, aes(x=MOTIF, y=DMAR_wMOTIF,fill=DMAR_type,label = DMAR_wMOTIF)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("Motifs occurrence")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Motif", y = "# of DM/ARs with motif")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p1

p2<- ggplot(data, aes(x=MOTIF, y=MOTIF_pct, label=MOTIF_pct, fill=DMAR_type)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("Motifs occurrence")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Motif", y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p2

pdf("ALX1_ALX3_ALX4_GBX2_TFEC_SOX10_MOTIF_occurrence_in_Iri_DMR_DAR_DMAR.pdf", width = 10, height = 8)
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

# DMR_DAR_DMAR_occurence.pdf
data <- read.table("total_DMAR_distribution.txt", sep = "\t", header = T, stringsAsFactors=F)
data<-as.data.frame(data)
data$DMAR_type <- factor(data$DMAR_type, levels = c("solo_closeDAR","solo_openDAR","solo_hyperDMR","solo_hypoDMR","hypo_closing_DMAR","hypo_opening_DMAR"))
data$DMAR <-"DMAR"
p3 <- ggplot(data, aes(x=DMAR, y=DMAR_number, fill=DMAR_type,label=DMAR_number)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("DM/AR number")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()
p4 <- ggplot(data, aes(x=DMAR, y=pct, fill=DMAR_type,label=pct)) +
  geom_bar(position="stack",stat="identity",colour = "black")+
  ggtitle("DM/AR number")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs( y = "pct")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_color_manual(values=mypalette[c(1,2,3,4,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ 
  theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))+geom_text(size = 3, fontface = "bold",position = position_stack(vjust = 0.5))+coord_flip()

pdf("DMR_DAR_DMAR_occurence.pdf", width = 12, height = 4)
grid.draw(rbind(ggplotGrob(p3), ggplotGrob(p4), size = "last"))
dev.off()
```
