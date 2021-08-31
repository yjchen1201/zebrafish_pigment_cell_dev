### Tool version
- cutadapt v1.10: Used to trim adapter reads 
- samtools version 1.3.1: Used to sort and downsample bam files for downstream processing
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


### Example used for Analysis
```{bash}
############################ {Bash script} ######################################
# adapter trimming
cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 25 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1_fastq.gz" "SAMPLE_R2_fastq.gz"; done
# Align to danRer10 using bwa mem
bwa mem "danRer10.fa" "SAMPLE_R1_trimmed.fastq.gz" "SAMPLE_R2_trimmed.fastq.gz" > "SAMPLE.sam"
# convert to bam file and only MAPQ>=30
samtools view -b -q 30 -o "SAMPLE_MAPQ30.bam" "SAMPLE_R1.sam" #filter quality based on 30
# sort bam
samtools sort -o "SAMPLE_MAPQ30.sorted.bam" "SAMPLE_MAPQ30.bam"
# Remove duplicates
java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=True REMOVE_DUPLICATES=True "SAMPLE_MAPQ30.sorted.bam" "SAMPLE_MAPQ30.sorted.rmdup" "SAMPLE_MAPQ30.sorteddup_removed.txt" ASSUME_SORTED=true

#Downsample reads to 35million (Example using samtool view -s subsample command)
samtools view -h -b -s 0.3125 "SAMPLE_MAPQ30.sorted.rmdup.bam" > "SAMPLE.downsampled.bam"

#Calculation Reads fraction in peaks (RFIP)#
for i in *downsampled.bam; do samtools view $i | wc -l ;done
#35020476 15somite_pos_Rep1.downsampled.bam
#34969034 15somite_pos_Rep2.downsampled.bam
#35007919 Iri_Rep1.downsampled.bam
#35016630 Iri_Rep2.downsampled.bam
#34991229 Mel_Rep1.downsampled.bam
#35016721 Mel_Rep4.downsampled.bam
#35040232 24hpf_pos_Rep1.downsampled.bam
#34975335 24hpf_pos_Rep2.downsampled.bam

# using MethylQA to covert bam to single-end read bed file
methylQA density -Q 10 -r -T -E 0 -I 2000 -o "SAMPLE_MAPQ30.sorted.rmdup.downsampled" "danRer10_lite.size" "SAMPLE_MAPQ30.sorted.rmdup.downsampled.bam"

# shift transposome inserted site
awk -v OFS="\t" '{if ($6=="+") print $1,$2+4,$3+4,$4,$5,$6; else if ($6=="-") print $1,$2-5,$3-5,$4,$5,$6}' "SAMPLE_MAPQ30.sorted.rmdup.downsampled.extend.bed" > "SAMPLE_MAPQ30.sorted.rmdup.downsampled.Tshift.bed}"
#NOTE THAT THIS CAN LEAD TO NEGATIVE VALUES FOR READS AT ENDS OF CHROMOSOMES, so make sure the start site is >0
```

```{R}
############################ {R script} ######################################
# make Genome browser files on just the ends of ATAC reads 
# import shifted bed file, Iri and Mel samples
Iri1 <- read.table("Iri_Rep1.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
Iri2 <- read.table("Iri_Rep2.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
Mel1 <- read.table("Mel_Rep1.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
Mel2 <- read.table("Mel_Rep2.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)

# fix the start site to 1 based
Iri1$start <- apply(Iri1, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
Iri1$end <- Iri1$start+1
Iri1$count <- 1
Iri1_aggregate <- ddply(Iri1, .(V1,start,end), summarise, coverage = sum(count))

Iri2$start <- apply(Iri2, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
Iri2$end <- Iri2$start+1
Iri2$count <- 1
Iri2_aggregate <- ddply(Iri2, .(V1,start,end), summarise, coverage = sum(count))

Mel1$start <- apply(Mel1, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
Mel1$end <- Mel1$start+1
Mel1$count <- 1
Mel1_aggregate <- ddply(Mel1, .(V1,start,end), summarise, coverage = sum(count))

Mel2$start <- apply(Mel2, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
Mel2$end <- Mel2$start+1
Mel2$count <- 1
Mel2_aggregate <- ddply(Mel2, .(V1,start,end), summarise, coverage = sum(count))

#Create shifted bed files
write.table(Iri1_aggregate, "Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Iri2_aggregate, "Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel1_aggregate, "Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel2_aggregate, "Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)

# Import shifted bed files, 24hpf and 15 somite samples
s24_1 <- read.table("24hpf_pos_Rep1.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
s24_2 <- read.table("24hpf_pos_Rep2.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
s15_1 <- read.table("15somite_pos_Rep1.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
s15_2 <- read.table("15somite_pos_Rep2.downsampled.Tshift.fixed.bed", sep = "\t", header = F, stringsAsFactors = F)
#fix the start site to 1 based
s24_1$start <- apply(s24_1, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
s24_1$end <- s24_1$start+1
s24_1$count <- 1
s24_1_aggregate <- ddply(s24_1, .(V1,start,end), summarise, coverage = sum(count))

s24_2$start <- apply(s24_2, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
s24_2$end <- s24_2$start+1
s24_2$count <- 1
s24_2_aggregate <- ddply(s24_2, .(V1,start,end), summarise, coverage = sum(count))

s15_1$start <- apply(s15_1, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
s15_1$end <- s15_1$start+1
s15_1$count <- 1
s15_1_aggregate <- ddply(s15_1, .(V1,start,end), summarise, coverage = sum(count))

s15_2$start <- apply(s15_2, 1, function(x) {if (x[6] == "+") {as.numeric(x[2])} else {as.numeric(x[3])}})
s15_2$end <- s15_2$start+1
s15_2$count <- 1
s15_2_aggregate <- ddply(s15_2, .(V1,start,end), summarise, coverage = sum(count))

#Create shifted bed files
write.table(s15_1_aggregate, "15somite_pos_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(s15_2_aggregate, "15somite_pos_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(s24_1_aggregate, "24hpf_pos_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(s24_2_aggregate, "24hpf_pos_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)

#combine Replicates together for Centipede
Iri1_aggregate <- read.table("Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
Iri2_aggregate <- read.table("Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
Iricombined_aggregate <- merge(Iri1_aggregate,Iri2_aggregate, by=c("V1","V2","V3"), all = T)
Iricombined_aggregate[is.na(Iricombined_aggregate)] <-0
Iricombined_aggregate$combined <- Iricombined_aggregate$V4.x+Iricombined_aggregate$V4.y 
write.table(Iricombined_aggregate[,c(1,2,3,6)], "Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)

Mel1_aggregate <- read.table("Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
Mel2_aggregate <- read.table("Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
Melcombined_aggregate <- merge(Mel1_aggregate,Mel2_aggregate, by=c("V1","V2","V3"), all = T)
Melcombined_aggregate[is.na(Melcombined_aggregate)] <-0
Melcombined_aggregate$combined <- Melcombined_aggregate$V4.x+Melcombined_aggregate$V4.y 
write.table(Melcombined_aggregate[,c(1,2,3,6)], "Mel_combined_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)

s24_1_aggregate <- read.table("24hpf_pos_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
s24_2_aggregate <- read.table("24hpf_pos_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", header = F, stringsAsFactors = F)
s24_combined_aggregate <- merge(s24_1_aggregate,s24_2_aggregate, by=c("V1","V2","V3"), all = T)
s24_combined_aggregate[is.na(s24_combined_aggregate)] <-0
s24_combined_aggregate$combined <- s24_combined_aggregate$V4.x+s24_combined_aggregate$V4.y 
write.table(s24_combined_aggregate[,c(1,2,3,6)], "24hpf_pos_combined_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
```
```{bash}
################################ {Bash script ###################################
#MACS2 call peak on downsampled (35million) files
macs2 callpeak -t "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed" -f BED -g 1.4e+9 -n "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01" -B --verbose 3 --SPMR --keep-dup all --nomodel -s 75 --extsize 73 --shift -37 -p .01 --outdir ../macs2

# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
sort -k 8gr,8gr "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak" | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc >"SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz"

#Find IDR peaks using 0.05 as the threshold 
idr --samples "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz" "SAMPLE_Rep2_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz" --input-file-type "narrowPeak" --rank "p.value" -o "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01.IDR_peaks.txt" --log-output-file "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01.IDR_peaks.log" --soft-idr-threshold 0.05 --plot --verbose

# count IDR peaks
for i in *peaks.txt; do wc -l $i;done
#172785 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#167412 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#114207 Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#103868 Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
```

```{R}
############################################ {R script } ##################################################
library(reshape)
library(RColorBrewer)
library(ggplot2)
mypalette <- c(brewer.pal(6,"Pastel1"),brewer.pal(6,"Pastel2"))
setwd("IDR_peaks_0.05")
Mel<- read.table("IDR_peaks_0.05/Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
Iri<- read.table("IDR_peaks_0.05/Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
hpf24<- read.table("IDR_peaks_0.05/24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
s15<- read.table("IDR_peaks_0.05/15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")

Mel$size <- Mel$V3-Mel$V2
Iri$size <- Iri$V3-Iri$V2
hpf24$size <- hpf24$V3-hpf24$V2
s15$size <- s15$V3-s15$V2

blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

d<-data.frame(c(rep("Mel",nrow(Mel)),rep("Iri",nrow(Iri)),rep("s24",nrow(hpf24)),rep("s15",nrow(s15))),c(Mel$size,Iri$size,hpf24$size,s15$size))
colnames(d) <- c("NCC","Size")
a<-ggplot(d, aes(Size, fill = NCC, colour = NCC)) + geom_density(alpha = 0.5, adjust = 0.5)+ggtitle("test")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12))+
      labs(x = "Size", y = "Density")+scale_x_continuous(lim = c(0,1000))+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
```
``` {bash}
############################################ {Bash script } ##################################################
# Merge all Peaks into master peak#
cat Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/ATAC_only/IDR_peaks_0.05/15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks.txt
# sort based on Choord
sort -k1,1 -k2,2n All_IDR_peaks.txt > All_IDR_peaks.sorted.txt
# Merge peaks
bedtools merge -i All_IDR_peaks.sorted.txt > All_IDR_peaks_merged.txt
# count merged peak numbers
wc -l  All_IDR_peaks_merged.txt # 262676

#Identify if sample peak is called in ALL IDR merged peaks#
bedtools intersect -c -a All_IDR_peaks_merged.txt -b 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_15s.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_24hpf.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_Mel.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_Iri.txt
```

```{R}
############################################ {R script } ##################################################
#Then merge the IDR peaks in R
Mel<- read.table("All_IDR_peaks_merged_Mel.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
Iri<- read.table("All_IDR_peaks_merged_Iri.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
s24<- read.table("All_IDR_peaks_merged_24hpf.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
s15<- read.table("All_IDR_peaks_merged_15s.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")

Mel$peakname <- paste(Mel$V1,":",Mel$V2,"-",Mel$V3,sep = "")
Iri$peakname <- paste(Iri$V1,":",Iri$V2,"-",Iri$V3,sep = "")
s24$peakname <- paste(s24$V1,":",s24$V2,"-",s24$V3,sep = "")
s15$peakname <- paste(s15$V1,":",s15$V2,"-",s15$V3,sep = "")
s15$size <- s15$V3-s15$V2
t <- merge(s15[,c("peakname","V1","V2","V3","size","V4")], s24[,c("peakname","V4")], by = "peakname", all = T)
t <- merge(t, Mel[,c("peakname","V4")], by = "peakname", all = T)
t <- merge(t, Iri[,c("peakname","V4")], by = "peakname", all = T)
write.table(t[c(2,3,4,1,5,6,7,8,9)], "All_IDR_peaks_Merged.txt", sep = "\t", col.names = F, row.names = F, quote = F)

colnames(t) <- c("peakname","chr","start","end","size","s15","s24","Mel","Iri")

allpeaks # 262676
all15s_peaks <- t[t$s15 > 0,] #165643
all24s_peaks <- t[t$s24 > 0,] #159141
allMel_peaks <- t[t$Mel > 0,] #98716
allIri_peaks <- t[t$Iri > 0,] #107012

only15s_peaks <- t[t$s15 > 0 & t$s24 == 0 & t$Mel == 0 & t$Iri == 0,] #39812
s15_24s_peaks <- t[t$s15 > 0 & t$s24 > 0 & t$Mel == 0 & t$Iri == 0,] #49871
s15_Mel_peaks <- t[t$s15 > 0 & t$s24 == 0 & t$Mel > 0 & t$Iri == 0,] #2528
s15_Iri_peaks <- t[t$s15 > 0 & t$s24 == 0 & t$Mel == 0 & t$Iri > 0,] #3000
s15_24s_Mel_peaks <- t[t$s15 > 0 & t$s24 > 0 & t$Mel > 0 & t$Iri == 0,] #9561
s15_24s_Iri_peaks <- t[t$s15 > 0 & t$s24 > 0 & t$Mel == 0 & t$Iri > 0,] #13893
s15_24s_Mel_Iri_peaks <- t[t$s15 > 0 & t$s24 > 0 & t$Mel > 0 & t$Iri > 0,] #45206
s15_Mel_Iri_peaks <- t[t$s15 > 0 & t$s24== 0 & t$Mel > 0 & t$Iri > 0,] #1772
only24s_peaks <- t[t$s15 == 0 & t$s24 > 0 & t$Mel == 0 & t$Iri == 0,] #28006
Mel_24s_peaks <- t[t$s15 == 0 & t$s24 > 0 & t$Mel > 0 & t$Iri == 0,] #3075
Iri_24s_peaks <- t[t$s15 == 0 & t$s24 > 0 & t$Mel == 0 & t$Iri > 0,] #4139
Mel_24s_Iri_peaks  <- t[t$s15 == 0 & t$s24 > 0 & t$Mel > 0 & t$Iri > 0,] #5390
onlyMel_peaks <- t[t$s15 == 0 & t$s24 == 0 & t$Mel > 0 & t$Iri == 0,] #22811
onlyIri_peaks <- t[t$s15 == 0 & t$s24 == 0 & t$Mel == 0 & t$Iri > 0,] #25239
Mel_Iri_peaks <- t[t$s15 == 0 & t$s24 == 0 & t$Mel > 0 & t$Iri > 0,] #8373

Early <- rbind(only15s_peaks,s15_24s_peaks,only24s_peaks)#117689
Intermediate <- rbind(s15_24s_Mel_peaks,Mel_24s_peaks,s15_24s_Iri_peaks,Iri_24s_peaks,s15_Mel_Iri_peaks,Mel_24s_Iri_peaks,Mel_Iri_peaks)#46203
Mel_Iri_shared_dynamce <- rbind(s15_Mel_Iri_peaks,Mel_Iri_peaks) #10145
Mel_specific <- rbind(s15_Mel_peaks,onlyMel_peaks)#25339
Iri_specific <- rbind(s15_Iri_peaks,onlyIri_peaks)#28239
allpeaks <-rbind(s15_24s_Mel_Iri_peaks,Early,Intermediate,Mel_specific, Iri_specific)

write.table(Mel_specific[,c(2,3,4,1)],"Mel_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Iri_specific[,c(2,3,4,1)],"Iri_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel_Iri_shared_dynamce[,c(2,3,4,1)],"Mel_Iri_shared_dynamic_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(rbind(s15_Mel_Iri_peaks[,c(2,3,4,1)],Mel_Iri_peaks[,c(2,3,4,1)]),"Mel_Iri_shared_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(rbind(s15_24s_peaks[,c(2,3,4,1)],only24s_peaks[,c(2,3,4,1)]),"24hpf_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)

#heatmap#
name <- allpeaks[,1]
df.OG <- allpeaks[,c(6:9)]
row.names(df.OG) <- name

library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
ht = Heatmap(df.OG, col = colorRamp2(c(0, 1, 2,3), c("#f6f6f6", "#74d600","#7bc043","#028900")), 
    cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = FALSE)
ht
```

```{bash}
############################################ {Bash script } ##################################################v
#Identify IDR differentially accessible regions(DARs) using bedtools intersect#
bedtools intersect -v -a Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Mel_only_IDR_peak_MelvsIri.txt
bedtools intersect -v -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -a Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Iri_only_IDR_peak_MelvsIri.txt

bedtools intersect -v -a 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > 24hpf_only_IDR_peak_24hpfvsMel.txt
bedtools intersect -v -a 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > 24hpf_only_IDR_peak_24hpfvsIri.txt

bedtools intersect -v -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -a Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Mel_only_IDR_peak_24hpfvsMel.txt
bedtools intersect -v -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -a Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Iri_only_IDR_peak_24hpfvsIri.txt

bedtools intersect -wo -a Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Both_IDR_peak_MelvsIri.txt
bedtools intersect -wo -a 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Both_IDR_peak_24hpfvsMel.txt
bedtools intersect -wo -a 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Both_IDR_peak_24hpfvsIri.txt
```

```{R}
############################################ {R script } ##################################################
### DiffBIND to call DARs (with IDR peaks and FDR 0.001) ("R analysis")###
library(DiffBind)
setwd("DiffBind/DiffBind_IDRpeaks")
samples <- read.csv("DiffBind/Pigment_Diffbind_ATAC.csv")
data<- dba(sampleSheet="DiffBind/Pigment_Diffbind_ATAC.csv")
data<- dba.count(data)
plot(data)
data$config$th = 0.01
data <- dba.contrast(data, categories=DBA_TISSUE, minMembers = 2)
data<- dba.analyze(data,method = DBA_DESEQ2)

#8 Samples, 262704 sites in matrix:
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
```
```{bash}
############################################ {Bash script } ##################################################
#make master diff peak
cat Diffbind_15v24_ALL_peaks.bed  Diffbind_24vIri_ALL_peaks.bed  Diffbind_24vMel_ALL_peaks.bed  Diffbind_MelvIri_ALL_peaks.bed > All_DiffBind_peaks.txt
sort -k1,1 -k2,2n All_DiffBind_peaks.txt > All_DiffBind_peaks.sorted.txt
bedtools merge -i All_DiffBind_peaks.sorted.txt > All_DiffBind_peaks_merged.txt
wc -l  All_DiffBind_peaks_merged.txt 
rm All_DiffBind_peaks.txt
rm All_DiffBind_peaks.sorted.txt

bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_15v24_ALL_peaks.bed > All_DiffBind_peaks_merged_15v24.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_24vMel_ALL_peaks.bed > All_DiffBind_peaks_merged_24vMel.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_24vIri_ALL_peaks.bed > All_DiffBind_peaks_merged_24vIri.txt
bedtools intersect -wao -a All_DiffBind_peaks_merged.txt -b Diffbind_MelvIri_ALL_peaks.bed > All_DiffBind_peaks_merged_MelvIri.txt
```

```{R}
############################################ {R script } ##################################################
s15_s24<- read.table("All_DiffBind_peaks_merged_15v24.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_M <- read.table("All_DiffBind_peaks_merged_24vMel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
s24_I <- read.table("All_DiffBind_peaks_merged_24vIri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
M_I <- read.table("All_DiffBind_peaks_merged_MelvIri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

s15_s24$chrompos<- paste(s15_s24$V1,":",s15_s24$V2,"-",s15_s24$V3,sep = "")
s24_M$chrompos<- paste(s24_M$V1,":",s24_M$V2,"-",s24_M$V3,sep = "")
s24_I$chrompos<- paste(s24_I$V1,":",s24_I$V2,"-",s24_I$V3,sep = "")
M_I$chrompos<- paste(M_I$V1,":",M_I$V2,"-",M_I$V3,sep = "")


DAR <-Reduce(function(x, y) merge(x, y,by = "chrompos", all = T), list(s15_s24[,c(9,1:3,7)],s24_M[,c(9,7)],s24_I[,c(9,7)],M_I[,c(9,7)]))
colnames(DAR) <- c("chrompos","chr","start","end","s15vs24","s24vMel","s24vIri","MelvIri")
DAR$size <- DAR$end-DAR$start
DAR[DAR == "."] <- 0
DAR$s15vs24 <- as.numeric(DAR$s15vs24)
DAR$s24vMel <- as.numeric(DAR$s24vMel)
DAR$s24vIri <- as.numeric(DAR$s24vIri)
DAR$MelvIri <- as.numeric(DAR$MelvIri)

DAR <- DAR[,c(2,3,4,1,9,5,6,7,8)]
DAR <- DAR[!duplicated(DAR$chrompos),]
write.table(DAR, "All_DiffBind_peaks_merged_wINFO.bed", row.names = F, col.names = F, sep = "\t",quote =F)
DAR <- read.table("All_DiffBind_peaks_merged_wINFO.bed",header =F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(DAR) <- c("chr","start","end","chrompos","size","s15vs24","s24vMel","s24vIri","MelvIri")


str(DAR[DAR$s15vs24 !=0,]) #894 
str(DAR[DAR$s24vMel !=0,]) #36336
str(DAR[DAR$s24vIri !=0,]) #34428
str(DAR[DAR$MelvIri !=0,]) #14211

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

#DAR size distribution plot
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
```
