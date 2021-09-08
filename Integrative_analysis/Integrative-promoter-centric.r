#### Integrative_promoter_centric.r  ####
# For generating Figure 2c
# Obtain promoter region ranges
## Using danRer10 1kb upstream TSS as promoter from UCSC genome browser
## Pull out gene info from GTF file
system("awk '($3==\"gene\"){OFS=\"\\t\"; if ($7~/+/){print $1,$4-1000,$4+500,$10}; if ($7~/-/){print $1,$5-500,$5+1000,$10}}' Danio_rerio.GRCz10.85.gtf | sed 's/[\";]//g;' > Danio_rerio.GRCz10.85.GENE.PROMOTER.bed")
system("awk '($3==\"gene\"){OFS=\"\\t\"; if ($7~/+/){print $1,$4,$5,$10}; if ($7~/-/){print $1,$4,$5,$10}}' Danio_rerio.GRCz10.85.gtf| sed 's/[\";]//g;' > Danio_rerio.GRCz10.85.GENE.LOCATION.bed")
# Load promoter info 
prom <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t", stringsAsFactors =F, header = F, quote ="")
prom[prom$V2 <0,]$V2 <- 0
prom$chrompos <- paste(prom$V1,":",prom$V2,"-",prom$V3,sep = "")

write.table(prom,"Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t",quote = F, col.names = F, row.names = F)

# Overlap with CpG methylation
## Convert DSS.txt file from WGBS/WGBS_02_DML-DMR.r to BED format. Using Mel_Combined_DSS.txt as an example
system("awk 'NR>1{OFS=\"\\t\";print $1,$2,$2,$3,$4,1}' Mel_Combined_DSS.txt > Mel_Combined_DSS.bed")
## Intersect promoter bed file with the DSS bed files 
system("bedtools intersect -wo -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b 15somite_NCC_Combined_DSS.bed > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_15s.txt")
system("bedtools intersect -wo -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b 24hpf_NCC_Combined_DSS.bed > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_24s.txt")
system("bedtools intersect -wo -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b Mel_Combined_DSS.bed > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Mel.txt")
system("bedtools intersect -wo -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b Iri_Combined_DSS.bed > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wMETH_Iri.txt")

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
## Save as a single matrix for all promoter methylation 
m<-Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3"), all.x = T), list(prom,prom_meth15_aggregate[,c(1,2,3,6,7)],prom_meth24_aggregate[,c(1,2,3,6,7)],prom_methM_aggregate[,c(1,2,3,6,7)],prom_methI_aggregate[,c(1,2,3,6,7)]))

# Overlap IDR peaks
## narrowPeak.gz.IDR_peaks.txt files are from ATAC/ATAC_03_MACS2_Callpeak_IDR_peaks.sh
system("bedtools intersect -c -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks15s.bed")
system("bedtools intersect -c -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks24s.bed")
system("bedtools intersect -c -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksMel.bed")
system("bedtools intersect -c -a Danio_rerio.GRCz10.85.GENE.PROMOTER.bed -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksIri.bed")

prom_IDR15<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks15s.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDR24<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaks24s.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDRM<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksMel.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
prom_IDRI<- read.table("ENSEMBL_danRer10_gene_UNIQUE_promoters_1kb_upstream_wIDRpeaksIri.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
## Merge IDR peak with methylation matrix
m_IDR<-Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3"), all.x = T), list(m,prom_IDR15[,c(1,2,3,6)],prom_IDR24[,c(1,2,3,6)],prom_IDRM[,c(1,2,3,6)],prom_IDRI[,c(1,2,3,6)]))
names(m_IDR) <- c("chr","start","end","gene","chrompos","s15_CpGcount","s15_AveMeth","s24_CpGcount","s24_AveMeth","Mel_CpGcount","Mel_AveMeth","Iri_CpGcount","Iri_AveMeth","s15_IDRpeak","s24_IDRpeak","Mel_IDRpeak","Mel_IDRpeak")

m_IDR[is.na(m_IDR)]<-0

write.table(m_IDR, "Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo.bed", sep = "\t", quote =F, col.names =T, row.names = F)
write.table(m_IDR, "Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo2.bed", sep = "\t", quote =F, col.names =F, row.names = F)

# Find promoters intersect with DMARs
system("bedtools intersect -wao -a Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo2.bed -b All_DMAR_Combined_wINFO.annotated2.bed > Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation.bed")
prom_DMAR <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
names(prom_DMAR) <- c("chr","start","end","gene","chrompos","s15_CpGcount","s15_AveMeth","s24_CpGcount","s24_AveMeth","Mel_CpGcount","Mel_AveMeth","Iri_CpGcount","Iri_AveMeth","s15_IDRpeak","s24_IDRpeak","Mel_IDRpeak","Iri_IDRpeak", "chr_DMAR","start_DMAR" , "end_DMAR"  ,    "chrompos_DMAR",  "DMARsize", "DMRsize" , "DMRs15vs24" , "DMRs24vMel", "DMRs24vIri",  "DMRMelvIri" ,"DARsize" ,   "DARs15vs24",  "DARs24vMel", "DARs24vIri","DARMelvIri", "IDR_s15" ,"IDR_s24"  , "IDR_M"  , "IDR_I" ,    "Meth_s15"   ,"Meth_s24" , "Meth_Mel" , "Meth_Iri"   , "CpG_s15",  "CpG_s24",  "CpG_Mel" ,"CpG_Iri", "aveCpG"   , "CpGdensity100bp","Annotation", "Annotation2","prom_DMAR_overlap")
prom_DMAR[is.na(prom_DMAR)] <- 0
prom_DMAR[prom_DMAR =="."] <- 0

# Add gene expression info 
DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) 
DEG2 <- read.table("DEGs_combined_samples_p0.01_wTFinfo.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) 
colnames(DEG2_genenames) <- c("gene","genename")
DEG3 <- merge(DEG2,DEG2_genenames, by = "gene", all.x = T)
DEG4<-merge(geneTPM, DEG3,by = "gene", all.x =T) ##ALL GENES
DEG4[is.na(DEG4)] <-0
DEG4<-DEG4[DEG4$s15 > 5 | DEG4$s24 > 5 |DEG4$Mel > 5 |DEG4$Iri > 5,] #filter genes with at least 5 rpkm

prom_DMAR2 <- merge(prom_DMAR,DEG4[,c(1:13)],by="gene",all.x = T)
prom_DMAR2 <- prom_DMAR2[,c(2,3,4,1,5:61)]
prom_DMAR2[is.na(prom_DMAR2)] <-0
write.table(prom_DMAR2,"Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation_wDEGs.bed", sep = "\t", col.names = T, row.names =F, quote =F)

prom_DMAR2 <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER_wMETHinfo_wIDRinfo_wDMAR_wAnnotation_wDEGs.bed", sep = "\t", header = T, stringsAsFactors =F, quote = "")
prom_DMAR2[,c(22:46)] <- sapply(prom_DMAR2[,c(22:46)],as.numeric)

# DEG promoter centric analysis
prom_DMAR2_Mel_On <- prom_DMAR2[prom_DMAR2$s24vMel_log2_change < 0,]
prom_DMAR2_Iri_On <- prom_DMAR2[prom_DMAR2$s24vIri_log2_change < 0,]

prom_DMAR2_Mel_Off <- prom_DMAR2[prom_DMAR2$s24vMel_log2_change > 0 ,]
prom_DMAR2_Iri_Off <- prom_DMAR2[prom_DMAR2$s24vIri_log2_change > 0,]

## Mel On
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

## Iri on
prom_DMAR2_Iri_On_nochange <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri == 0 & prom_DMAR2_Iri_On$DARs24vIri ==0,]
prom_DMAR2_Iri_On_only_hypo<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri > 0 & prom_DMAR2_Iri_On$DARs24vIri ==0,]
prom_DMAR2_Iri_On_only_opening<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri == 0 & prom_DMAR2_Iri_On$DARs24vIri >0,]
prom_DMAR2_Iri_On_both_hypo_open<- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$DMRs24vIri > 0 & prom_DMAR2_Iri_On$DARs24vIri >0,]

gl <- c(prom_DMAR2_Iri_On_nochange$chrompos,prom_DMAR2_Iri_On_only_hypo$chrompos,prom_DMAR2_Iri_On_only_opening$chrompos,prom_DMAR2_Iri_On_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Iri_On_others <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_On$chrompos %out% gl,]

vals <- c(nrow(prom_DMAR2_Iri_On_nochange),nrow(prom_DMAR2_Iri_On_only_hypo),nrow(prom_DMAR2_Iri_On_only_opening),nrow(prom_DMAR2_Iri_On_both_hypo_open),nrow(prom_DMAR2_Iri_On_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hypoDMR", "Only openDAR", "Opening DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

bb<-waffle::waffle(vals,colors = mypalette1,rows = 18,size = 0.5,title = "Iri DEG's promoter dynamics (On)")


## Mel off
prom_DMAR2_Mel_Off_nochange <- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel == 0 & prom_DMAR2_Mel_Off$DARs24vMel ==0,]
prom_DMAR2_Mel_Off_only_hypo<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel > 0 & prom_DMAR2_Mel_Off$DARs24vMel ==0,]
prom_DMAR2_Mel_Off_only_opening<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel == 0 & prom_DMAR2_Mel_Off$DARs24vMel <0,]
prom_DMAR2_Mel_Off_both_hypo_open<- prom_DMAR2_Mel_Off[prom_DMAR2_Mel_Off$DMRs24vMel > 0 & prom_DMAR2_Mel_Off$DARs24vMel <0,]

gl <- c(prom_DMAR2_Mel_Off_nochange$chrompos,prom_DMAR2_Mel_Off_only_hypo$chrompos,prom_DMAR2_Mel_Off_only_opening$chrompos,prom_DMAR2_Mel_Off_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Mel_Off_others <- prom_DMAR2_Mel_On[prom_DMAR2_Mel_Off$chrompos %out% gl,]

vals <- c(nrow(prom_DMAR2_Mel_Off_nochange),nrow(prom_DMAR2_Mel_Off_only_hypo),nrow(prom_DMAR2_Mel_Off_only_opening),nrow(prom_DMAR2_Mel_Off_both_hypo_open),nrow(prom_DMAR2_Mel_Off_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hypoDMR", "Only closingDAR", "Closing DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

cc<-waffle::waffle(vals,colors = mypalette1,rows = 18,size = 0.5,title = "Mel DEG's promoter dynamics (Off)")


## Iri off
prom_DMAR2_Iri_Off_nochange <- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri == 0 & prom_DMAR2_Iri_Off$DARs24vIri ==0,]
prom_DMAR2_Iri_Off_only_hypo<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri < 0 & prom_DMAR2_Iri_Off$DARs24vIri ==0,]
prom_DMAR2_Iri_Off_only_opening<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri == 0 & prom_DMAR2_Iri_Off$DARs24vIri <0,]
prom_DMAR2_Iri_Off_both_hypo_open<- prom_DMAR2_Iri_Off[prom_DMAR2_Iri_Off$DMRs24vIri < 0 & prom_DMAR2_Iri_Off$DARs24vIri <0,]

gl <- c(prom_DMAR2_Iri_Off_nochange$chrompos,prom_DMAR2_Iri_Off_only_hypo$chrompos,prom_DMAR2_Iri_Off_only_opening$chrompos,prom_DMAR2_Iri_Off_both_hypo_open$chrompos)

'%out%' <- function(x,y)!('%in%'(x,y))
prom_DMAR2_Iri_Off_others <- prom_DMAR2_Iri_On[prom_DMAR2_Iri_Off$chrompos %out% gl,]

vals <- c(nrow(prom_DMAR2_Iri_Off_nochange),nrow(prom_DMAR2_Iri_Off_only_hypo),nrow(prom_DMAR2_Iri_Off_only_opening),nrow(prom_DMAR2_Iri_Off_both_hypo_open),nrow(prom_DMAR2_Iri_Off_others))
val_names <- sprintf("%s (%s)", c("No change", "Only hyperDMR", "Only closingDAR", "Closing DMAR","Other"), scales::percent(round(vals/sum(vals), 2)))
names(vals) <- val_names

dd<-waffle::waffle(vals,colors = mypalette1,rows =18,size = 0.5,title = "Iri DEG's promoter dynamics (Off)")

multiplot(aa,bb,cc,dd,cols=1)