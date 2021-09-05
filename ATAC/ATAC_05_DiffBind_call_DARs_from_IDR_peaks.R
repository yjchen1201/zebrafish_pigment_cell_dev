#Tool version
#- DiffBind v 2.2.12: Used to call differentially accessible regions/peaks

# DiffBIND to call DARs (with IDR peaks and FDR 0.001) 
library(DiffBind)
setwd("DiffBind_IDRpeaks")
samples <- read.csv("Pigment_Diffbind_ATAC.csv") #This file is provided on github
data<- dba(sampleSheet="Pigment_Diffbind_ATAC.csv")
data<- dba.count(data)
plot(data)
data$config$th = 0.01
data <- dba.contrast(data, categories=DBA_TISSUE, minMembers = 2)
data<- dba.analyze(data,method = DBA_DESEQ2)

#8 Samples, 262704 sites in Pigment_Diffbind_ATAC.csv matrix:
#     ID Tissue Factor  Condition  Treatment Replicate Caller Intervals FRiP
#1  Mel1    Mel     ER  Resistant Full-Media         1 counts    262704 0.36#
#2  Mel2    Mel     ER  Resistant Full-Media         2 counts    262704 0.33
#3  Iri1    Iri     ER Responsive Full-Media         1 counts    262704 0.36
#4  Iri2    Iri     ER Responsive Full-Media         2 counts    262704 0.42
#5 s24-1    s24     ER Responsive Full-Media         1 counts    262704 0.45
#6 s24-2    s24     ER Responsive Full-Media         2 counts    262704 0.48
#7 s15-1    s15     ER Responsive Full-Media         1 counts    262704 0.42
#8 s15-2    s15     ER  Resistant Full-Media         2 counts    262704 0.36

#6 Contrasts:
#  Group1 Members1 Group2 Members2 DB.DESeq2
#1    Mel        2    Iri        2     22098
#2    Mel        2    s24        2     55958
#3    Mel        2    s15        2     53249
#4    Iri        2    s24        2     52791
#5    Iri        2    s15        2     55948
#6    s24        2    s15        2      2783

#120350 15somite_neg.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#172785 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#140323 24hpf_neg.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#167412 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#114207 Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#103868 Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt


mypalette <- brewer.pal(12,"Paired")
me <- data.frame(c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"),c(0.42,0.36,0.45,0.48,0.36,0.33,0.36,0.42),c("172,785\nIDR peaks","172,785\nIDR peaks","167,412\nIDR peaks","167,412\nIDR peaks","103,868\nIDR peaks","103,868\nIDR peaks","114207\nIDR peaks","114,207\nIDR peaks"))
colnames(me) <- c("Type", "FRiP","IDR")
me$Type <- factor(me$Type, levels = c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"))
p <- ggplot(data=me, aes(x=Type, y=FRiP, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("ATAC-seq QC (FRiP)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "Fraction of Reads in Peaks (FRiP)")+scale_y_continuous(lim = c(0,0.54))+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+geom_text(data=me,aes(x=Type,y=0.52,label=IDR),size=7)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
p<-p+geom_hline(aes(yintercept = 0.3), linetype = "dashed", lwd = 1.5, color = "#333333")
p+geom_hline(aes(yintercept = 0.2), linetype = "dashed", lwd = 1, color = "#777777")


plot(data,contrast =1)
dba.plotPCA(data,DBA_TISSUE,label=DBA_TISSUE)
data.DB_15v24 <- dba.report(data, th = 0.001,contrast = 6)
data.DB_24vIri <- dba.report(data, th = 0.001, contrast = 4)
data.DB_24vMel <- dba.report(data, th = 0.001, contrast = 2)
data.DB_MelvIri <- dba.report(data, th = 0.001, contrast = 1)

data.DB_15v24 <- data.frame(data.DB_15v24)
data.DB_24vMel <- data.frame(data.DB_24vMel)
data.DB_24vIri <- data.frame(data.DB_24vIri)
data.DB_MelvIri <- data.frame(data.DB_MelvIri)

data.DB_15v24$seqnames <-as.character(data.DB_15v24$seqnames)
data.DB_24vMel$seqnames <-as.character(data.DB_24vMel$seqnames)
data.DB_MelvIri$seqnames <-as.character(data.DB_MelvIri$seqnames)
data.DB_24vIri$seqnames <-as.character(data.DB_24vIri$seqnames)

write.table(data.DB_15v24[,c(1,2,3,9)],"Diffbind_15v24_ALL_peaks.bed",sep = "\t", quote = F, col.names =F,row.names =F)
write.table(data.DB_24vMel[,c(1,2,3,9)],"Diffbind_24vMel_ALL_peaks.bed",sep = "\t", quote = F, col.names =F,row.names =F)
write.table(data.DB_24vIri[,c(1,2,3,9)],"Diffbind_24vIri_ALL_peaks.bed",sep = "\t", quote = F, col.names =F,row.names =F)
write.table(data.DB_MelvIri[,c(1,2,3,9)],"Diffbind_MelvIri_ALL_peaks.bed",sep = "\t", quote = F, col.names =F,row.names =F)

# count diff bind open/closing peaks in each comparison
nrow(data.DB_15v24[data.DB_15v24$Fold <0,]) #closing in 24hpf 8

nrow(data.DB_15v24[data.DB_15v24$Fold >0,]) #opening in 24hpf 886

nrow(data.DB_24vMel[data.DB_24vMel$Fold <0,]) #Closing in Melanophore 26947

nrow(data.DB_24vMel[data.DB_24vMel$Fold >0,]) #opening in Melanophore 9390

nrow(data.DB_24vIri[data.DB_24vIri$Fold <0,]) #closing in Iridophore 23717

nrow(data.DB_24vIri[data.DB_24vIri$Fold >0,]) #opening in Iridophore 10712

nrow(data.DB_MelvIri[data.DB_MelvIri$Fold <0,]) #closed in Mel and open in Iri 8110

nrow(data.DB_MelvIri[data.DB_MelvIri$Fold >0,]) #open in Mel and closed in Iri 6101

d_consensus <- dba.peakset(data, consensus=c(DBA_TISSUE), minOverlap=0.5)
par(mfrow=c(2,2))
dba.plotHeatmap(data,th = 0.001, contrast=6, correlations=F)
dba.plotHeatmap(data,th = 0.001, contrast=4, correlations=FALSE)
dba.plotHeatmap(data,th = 0.001, contrast=2, correlations=FALSE)
dba.plotHeatmap(data,th = 0.001, contrast=1, correlations=FALSE)

#dba.plotMA(data, contrast =6,th = 0.001,bXY = T)
dba.plotMA(data, contrast =6,th = 0.001,bXY = F)
#dba.plotMA(data, contrast =4,th = 0.001,bXY = T)
dba.plotMA(data, contrast =4,th = 0.001,bXY = F)
#dba.plotMA(data, contrast =2,th = 0.001,bXY = T)
dba.plotMA(data, contrast =2,th = 0.001,bXY = F)
#dba.plotMA(data, contrast =1,th = 0.001,bXY = T)
dba.plotMA(data, contrast =1,th = 0.001,bXY = F)

data.block <- dba.report(data,th = 0.001, method=DBA_ALL_METHODS_BLOCK,bDB=TRUE,bAll=TRUE)
dba.plotVenn(data.block,c(6,4,2),label1="15NCC vs 24NCC",label2="24NCC vs Mel",label3="24NCC vs Iri")

############################################ {Bash script} ##################################################
# Make master diff peak
cat Diffbind_15v24_ALL_peaks.bed  Diffbind_24vIri_ALL_peaks.bed  Diffbind_24vMel_ALL_peaks.bed  Diffbind_MelvIri_ALL_peaks.bed > All_DiffBind_peaks.txt
sort -k1,1 -k2,2n All_DiffBind_peaks.txt > All_DiffBind_peaks.sorted.txt
bedtools merge -i All_DiffBind_peaks.sorted.txt > All_DiffBind_peaks_merged.txt
rm All_DiffBind_peaks.txt
rm All_DiffBind_peaks.sorted.txt

# Identify if DAR peaks is called in ALL DiffBind merged peaks#
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_15v24_ALL_peaks.bed > All_DiffBind_peaks_merged_15v24.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_24vMel_ALL_peaks.bed > All_DiffBind_peaks_merged_24vMel.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_24vIri_ALL_peaks.bed > All_DiffBind_peaks_merged_24vIri.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_MelvIri_ALL_peaks.bed > All_DiffBind_peaks_merged_MelvIri.txt
##############################################################################################################

# Import files
s15_s24<- read.table("All_DiffBind_peaks_merged_15v24.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_M <- read.table("All_DiffBind_peaks_merged_24vMel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_I <- read.table("All_DiffBind_peaks_merged_24vIri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
M_I <- read.table("All_DiffBind_peaks_merged_MelvIri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

# Add chrom position column, eg: chr1:234-567
s15_s24$chrompos<- paste(s15_s24$V1,":",s15_s24$V2,"-",s15_s24$V3,sep = "")
s24_M$chrompos<- paste(s24_M$V1,":",s24_M$V2,"-",s24_M$V3,sep = "")
s24_I$chrompos<- paste(s24_I$V1,":",s24_I$V2,"-",s24_I$V3,sep = "")
M_I$chrompos<- paste(M_I$V1,":",M_I$V2,"-",M_I$V3,sep = "")

# Merge all the list
DAR <-Reduce(function(x, y) merge(x, y,by = "chrompos", all = T), list(s15_s24[,c(9,1:3,7)],s24_M[,c(9,7)],s24_I[,c(9,7)],M_I[,c(9,7)]))
# Add column name
colnames(DAR) <- c("chrompos","chr","start","end","s15vs24","s24vMel","s24vIri","MelvIri")
# Add DAR size information
DAR$size <- DAR$end-DAR$start
DAR[DAR == "."] <- 0
DAR$s15vs24 <- as.numeric(DAR$s15vs24)
DAR$s24vMel <- as.numeric(DAR$s24vMel)
DAR$s24vIri <- as.numeric(DAR$s24vIri)
DAR$MelvIri <- as.numeric(DAR$MelvIri)

DAR <- DAR[,c(2,3,4,1,9,5,6,7,8)]
DAR <- DAR[!duplicated(DAR$chrompos),]
# Create the DAR table
write.table(DAR, "All_DiffBind_peaks_merged_wINFO.bed", row.names = F, col.names = F, sep = "\t",quote =F)
# Re-import DAR info
DAR <- read.table("All_DiffBind_peaks_merged_wINFO.bed",header =F, sep = "\t", stringsAsFactors = F, quote = "")
# Add column names
colnames(DAR) <- c("chr","start","end","chrompos","size","s15vs24","s24vMel","s24vIri","MelvIri")

# Count DAR numbers
str(DAR[DAR$s15vs24 !=0,]) #894 
str(DAR[DAR$s24vMel !=0,]) #36336
str(DAR[DAR$s24vIri !=0,]) #34428
str(DAR[DAR$MelvIri !=0,]) #14211

# Extract cell type-specific/shared opening/closing DARs
Mel_specific_closing <- DAR[DAR$s24vMel <0 & DAR$s24vIri >=0,] #10289
Mel_specific_opening <- DAR[DAR$s24vMel >0 & DAR$s24vIri <=0,] #7012
Mel_specific_opening_stringent <- DAR[DAR$s24vMel >0 & DAR$s24vIri <=0 & DAR$MelvIri >0,] #3780
Mel_specific_closing_stringent <- DAR[DAR$s24vMel <0 & DAR$s24vIri >=0 & DAR$MelvIri <0,] #1304

Iri_specific_closing <- DAR[DAR$s24vMel >=0 & DAR$s24vIri <0,] #7059
Iri_specific_opening <- DAR[DAR$s24vMel <=0 & DAR$s24vIri >0,] #8334
Iri_specific_opening_stringent <- DAR[DAR$s24vMel <=0 & DAR$s24vIri >0 & DAR$MelvIri <0,] #4589
Iri_specific_closing_stringent <- DAR[DAR$s24vMel >=0 & DAR$s24vIri <0 & DAR$MelvIri >0,] #596

Mel_Iri_shared_closing <- DAR[DAR$s24vMel <0 & DAR$s24vIri <0,] #16657
Mel_Iri_shared_opening <- DAR[DAR$s24vMel >0 & DAR$s24vIri >0,] #2378

# Save them as table
write.table(Mel_specific_closing, "Mel_specific_closing_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_opening, "Mel_specific_opening_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_opening_stringent, "Mel_specific_opening_DAR_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_closing_stringent, "Mel_specific_closing_DAR_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Iri_specific_closing, "Iri_specific_closing_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_opening, "Iri_specific_opening_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_opening_stringent, "Iri_specific_opening_DAR_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_closing_stringent, "Iri_specific_closing_DAR_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Mel_Iri_shared_closing , "Mel_Iri_shared_closing_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_Iri_shared_opening, "Mel_Iri_shared_opening_DAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

# DAR size distribution plot
NCC_DAR <- DAR[DAR$s15vs24 != 0,]
M_DAR <- DAR[DAR$s24vMel != 0,]
I_DAR <- DAR[DAR$s24vIri != 0,]
MI_DAR <- DAR[DAR$MelvIri != 0,]

ps24 <-ggplot(NCC_DAR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[4],fill = mypalette[4])+ggtitle(paste("DAR size distribution\n(15somite vs 24hpf)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "DAR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM <-ggplot(M_DAR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[6],fill = mypalette[6])+ggtitle(paste("DAR size distribution\n(24hpf vs Melanophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DAR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI <-ggplot(I_DAR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[8],fill = mypalette[8])+ggtitle(paste("DAR size distribution\n(24hpf vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DAR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI <-ggplot(MI_DAR, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[10],fill = mypalette[10])+ggtitle(paste("DAR size distribution\n(Melanophore vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DAR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
multiplot(ps24,pI,pM,pMI, cols = 2)
