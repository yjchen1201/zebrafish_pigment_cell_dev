# DMAR master list 
## Prepare All_DMAR_Combined_wINFO with annotation. Code matching with 042620 version 
DMAR <- read.table("All_DMAR_Combined_wINFO.bed", sep = "\t", header = T, quote = "",stringsAsFactors = F)
anno_DMAR <- read.delim("All_DMAR_Combined_wINFO.annotated.txt",skip =1, header = F, quote = "",stringsAsFactors = F)
colnames(anno_DMAR)[1] <- "chrompos"
DMAR <- merge(DMAR,anno_DMAR[,c(1,8,9)],by = "chrompos",all.x = T)
colnames(DMAR)[34] <- "Annotation"
DMAR$Annotation <- apply(DMAR, 1, function(x) {unlist(strsplit(x[34], " \\("))[1]})
DMAR$V9 <- apply(DMAR, 1, function(x) {unlist(strsplit(x[35], " \\("))[1]})
colnames(DMAR)[35] <- "Annotation2"
## Write to table
DMAR <- DMAR[,c(2,3,4,1,5:35)]
write.table(DMAR, "All_DMAR_Combined_wINFO.annotated.bed", row.names = F, col.names = T, sep = "\t",quote =F)
## without header
write.table(DMAR, "All_DMAR_Combined_wINFO.annotated2.bed", row.names = F, col.names = F, sep = "\t",quote =F)

## Plot size distribution of DAR, DMR, DMAR [Figure S3b]
d <- data.frame(c(rep("DMR",nrow(DMAR[DMAR$DMRsize > 0,])),rep("DAR",nrow(DMAR[DMAR$DARsize > 0,])),rep("DMAR",nrow(DMAR))),c(DMAR[DMAR$DMRsize > 0,]$DMRsize,DMAR[DMAR$DARsize > 0,]$DARsize,DMAR$DMARsize))
colnames(d) <- c("Type","Size")
a <-ggplot(d, aes(Size, fill = Type, colour = Type)) + geom_density(alpha = 0, adjust = 1)+ggtitle(paste("DMAR size distribution"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Size (bp)", y = "Density")+scale_x_continuous(limits = c(0,2000))+scale_fill_manual(values = mypalette2[c(1,2,3)])+scale_color_manual(values=mypalette2[c(1,2,3)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a

## Plot CpG density and region size for solo DMR, solo DAR and DMARs [Figure S3c]
### hypothesize that regions with no change in methylation = CGIs and dynamic change = low CpG density
All_DMRs <- DMAR[DMAR$DMRsize != 0,] #32715
All_DARs <- DMAR[DMAR$DARsize !=0,] #54921
Only_DARs <- DMAR[DMAR$DMRsize == 0 & DMAR$DARsize >0,] #39669
Only_DMRs <- DMAR[DMAR$DARsize == 0 & DMAR$DMRsize >0,] #17463
dynamic_DMARs <- DMAR[DMAR$DMRsize != 0 & DMAR$DARsize !=0,] #15252

CD<-data.frame(c(rep("DMR",nrow(All_DMRs)), rep("DAR",nrow(All_DARs)), rep("soloDMR",nrow(Only_DMRs)), rep("soloDAR",nrow(Only_DARs)), rep("dynamicDMAR",nrow(dynamic_DMARs))),  c(All_DMRs$CpGdensity100bp,All_DARs$CpGdensity100bp,Only_DMRs$CpGdensity100bp,Only_DARs$CpGdensity100bp,dynamic_DMARs$CpGdensity100bp))
colnames(CD) <- c("Type","CpG_density")
CD$Type <- factor(CD$Type, levels = c("DMR","soloDMR","DAR","soloDAR","dynamicDMAR"))

b <-ggplot(CD, aes(x= Type,y=CpG_density))+ geom_violin(aes(fill = Type), trim = FALSE) +geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE,width=0.3)+ggtitle(paste("CpGs density distribution in Epigenetically Dynamic Regions"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "CpG density (# of CpGs in 100bp)", y = "Density")+scale_fill_manual(values = mypalette4)+scale_color_manual(values=mypalette4)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
b
# Identify different classes of DMARs 
## Calculate proportion of solitary hyper DMRs in hyper DMRs.
str(Only_DMRs[Only_DMRs$DMRs24vMel <0 | Only_DMRs$DMRs24vIri <0,]) #650
str(DMAR[DMAR$DMRs24vMel <0 | DMAR$DMRs24vIri <0,]) #695    ~93.5% hyperDMRs are solitary hyper DMRs and occur in regions with no peaks (650-7-9-15)/650*100 = >95%

## Obtain DMR specific subsets. Count number of regions in each category [Figure S4a]
Only_DMRs_15v24_hypo <-  Only_DMRs[Only_DMRs$DMRs15vs24>0,] #287
Only_DMRs_15v24_hyper <-  Only_DMRs[Only_DMRs$DMRs15vs24<0,] #155
Only_DMRs_Mel_specific_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel >0 & Only_DMRs$DMRs24vIri <= 0,]  #6701 
Only_DMRs_Mel_specific_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel <0 & Only_DMRs$DMRs24vIri >= 0,] #210
Only_DMRs_Iri_specific_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel <=0 & Only_DMRs$DMRs24vIri > 0,] #5672
Only_DMRs_Iri_specific_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel >= 0 & Only_DMRs$DMRs24vIri < 0,] #289
Only_DMRs_M_I_shared_hypo <- Only_DMRs[Only_DMRs$DMRs24vMel >0 & Only_DMRs$DMRs24vIri > 0,] #3886
Only_DMRs_M_I_shared_hyper <- Only_DMRs[Only_DMRs$DMRs24vMel <0 & Only_DMRs$DMRs24vIri < 0,] #151

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

write.table(Only_DMRs_M_I_shared_hypo_111, "Solo_DMRs_Mel_Iri_shared_hypo_alwaysOPEN.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F) 

Only_DMRs_M_I_shared_hyper_000 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M ==0,] #110
Only_DMRs_M_I_shared_hyper_100 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 >0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M ==0,] #17
Only_DMRs_M_I_shared_hyper_001 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I ==0 & Only_DMRs_M_I_shared_hyper$IDR_M > 0,] #0
Only_DMRs_M_I_shared_hyper_010 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I >0 & Only_DMRs_M_I_shared_hyper$IDR_M == 0,] #1
Only_DMRs_M_I_shared_hyper_011 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 ==0 & Only_DMRs_M_I_shared_hyper$IDR_I >0 & Only_DMRs_M_I_shared_hyper$IDR_M > 0,] #2
Only_DMRs_M_I_shared_hyper_111 <- Only_DMRs_M_I_shared_hyper[Only_DMRs_M_I_shared_hyper$IDR_s24 >0 & Only_DMRs_M_I_shared_hyper$IDR_I >0& Only_DMRs_M_I_shared_hyper$IDR_M >0,] #15

## Obtain DAR specific subsets. Count number of regions in each category [Figure S4a]
Only_DARs_15v24_opening <-  Only_DARs[Only_DARs$DARs15vs24>0,] #434
Only_DARs_15v24_closing <-  Only_DARs[Only_DARs$DARs15vs24<0,] #8

Only_DARs_Mel_specific_open <- Only_DARs[Only_DARs$DARs24vMel >0 & Only_DARs$DARs24vIri <= 0,]  #3605
Only_DARs_Mel_specific_closed <- Only_DARs[Only_DARs$DARs24vMel <0 & Only_DARs$DARs24vIri >= 0,] #8807

str(Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$Meth_Mel <30,]) #754
str(Only_DARs_Mel_specific_closed[Only_DARs_Mel_specific_closed$Meth_Mel <30,]) #3739

Only_DARs_Iri_specific_open <- Only_DARs[Only_DARs$DARs24vMel <=0 & Only_DARs$DARs24vIri > 0,] #3645
Only_DARs_Iri_specific_closed <- Only_DARs[Only_DARs$DARs24vMel >= 0 & Only_DARs$DARs24vIri < 0,] #6043

str(Only_DARs_Iri_specific_open[Only_DARs_Iri_specific_open$Meth_Iri <30,]) #1398
str(Only_DARs_Iri_specific_closed[Only_DARs_Iri_specific_closed$Meth_Iri <30,]) #902

Only_DARs_M_I_shared_open <- Only_DARs[Only_DARs$DARs24vMel >0 & Only_DARs$DARs24vIri > 0,] #589
Only_DARs_M_I_shared_closed <- Only_DARs[Only_DARs$DARs24vMel <0 & Only_DARs$DARs24vIri < 0,] #15058

write.table(Only_DARs_Mel_specific_open, "Solo_DARs_Mel_specific_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Mel_specific_closed, "Solo_DARs_Mel_specific_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Iri_specific_open, "Solo_DARs_Iri_specific_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_Iri_specific_closed, "Solo_DARs_Iri_specific_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_M_I_shared_open, "Solo_DARs_Mel_Iri_shared_opening.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Only_DARs_M_I_shared_closed, "Solo_DARs_Mel_Iri_shared_closing.DMAR.bed", row.names = F, col.names = F, sep = "\t",quote =F)

Functional_Only_DARs_Mel_specific_open <- Only_DARs_Mel_specific_open[Only_DARs_Mel_specific_open$Meth_Mel >0 & Only_DARs_Mel_specific_open$Meth_Mel <= 30 & Only_DARs_Mel_specific_open$Meth_s24 >0 & Only_DARs_Mel_specific_open$Meth_s24 <=30,]
as.data.frame(table(Functional_Only_DARs_Mel_specific_open$Annotation))
# Var1 Freq
# 1       5' UTR    1
# 2         exon    2
# 3   Intergenic   54
# 4       intron   19
# 5 promoter-TSS   23
# 6          TTS    3

# Distribution of CpG count and DMR/DMAR size of solo hypoDMRs and DMARs [Figure S4b]
## Extract specific DMR information out
soloMel_hypoDMR<-(DMAR[DMAR$DMRs24vMel > 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri <=0 & DMAR$DARs24vIri == 0 & DMAR$DARsize == 0 & DMAR$DMRsize >0,]) #6701 solo Mel hypoDMR 
soloMel_hyperDMR<-(DMAR[DMAR$DMRs24vMel < 0 & DMAR$DARs24vMel == 0 & DMAR$DMRs24vIri >= 0 & DMAR$DARs24vIri == 0 & DMAR$DARsize == 0 & DMAR$DMRsize >0,]) #210 solo Mel hyperDMR
soloMel_openDAR<-(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel > 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri <= 0 & DMAR$DMRsize == 0 & DMAR$DARsize >0,]) #3605 solo Mel openDAR
soloMel_closeDAR<-(DMAR[DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel < 0 & DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri >= 0 & DMAR$DMRsize == 0 & DMAR$DARsize >0,]) #8807 solo Mel closeDAR

soloIri_hypoDMR<-(DMAR[DMAR$DMRs24vIri > 0 & DMAR$DARs24vIri == 0 & DMAR$DMRs24vMel <=0 & DMAR$DARs24vMel == 0 &DMAR$DARsize == 0 & DMAR$DMRsize >0,]) #5672 solo Iri hypoDMR 
soloIri_hyperDMR<-(DMAR[DMAR$DMRs24vIri < 0 & DMAR$DARs24vIri == 0 & DMAR$DMRs24vMel >= 0 & DMAR$DARs24vMel == 0 & DMAR$DARsize == 0 & DMAR$DMRsize >0,]) #289 solo Iri hyperDMR
soloIri_openDAR<-(DMAR[DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri > 0 & DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel <= 0 & DMAR$DMRsize == 0 & DMAR$DARsize >0,]) #3645 solo Iri openDAR
soloIri_closeDAR<-(DMAR[DMAR$DMRs24vIri == 0 & DMAR$DARs24vIri < 0 & DMAR$DMRs24vMel == 0 & DMAR$DARs24vMel >= 0 & DMAR$DMRsize == 0 & DMAR$DARsize >0,]) #6043 solo Iri closeDAR

soloMel_hypoDMR_open <- soloMel_hypoDMR[soloMel_hypoDMR$IDR_s24 > 0 & soloMel_hypoDMR$IDR_M > 0,] #2370
soloMel_hypoDMR_other <- soloMel_hypoDMR[(soloMel_hypoDMR$IDR_s24 > 0 & soloMel_hypoDMR$IDR_M == 0) | (soloMel_hypoDMR$IDR_s24 == 0 & soloMel_hypoDMR$IDR_M > 0),] #2338
soloMel_hypoDMR_close <- soloMel_hypoDMR[soloMel_hypoDMR$IDR_s24 == 0 & soloMel_hypoDMR$IDR_M == 0,] #1993
soloIri_hypoDMR_open <- soloIri_hypoDMR[soloIri_hypoDMR$IDR_s24 > 0 & soloIri_hypoDMR$IDR_I > 0,] #1643
soloIri_hypoDMR_other <- soloIri_hypoDMR[(soloIri_hypoDMR$IDR_s24 > 0 & soloIri_hypoDMR$IDR_I == 0) | (soloIri_hypoDMR$IDR_s24 == 0 & soloIri_hypoDMR$IDR_I > 0),] #2046
soloIri_hypoDMR_close <- soloIri_hypoDMR[soloIri_hypoDMR$IDR_s24 == 0 & soloIri_hypoDMR$IDR_I == 0,] #1983

## Generate violin plot
t<-data.frame(c(soloMel_hypoDMR_open$DMRsize,soloMel_hypoDMR_other$DMRsize,soloMel_hypoDMR_close$DMRsize,soloIri_hypoDMR_open$DMRsize,soloIri_hypoDMR_other$DMRsize,soloIri_hypoDMR_close$DMRsize,openingDMAR_s24vsMel$DMRsize,closingDMAR_s24vsMel$DMRsize,openingDMAR_s24vsIri$DMRsize, closingDMAR_s24vsIri$DMRsize),c(soloMel_hypoDMR_open$aveCpGcount,soloMel_hypoDMR_other$aveCpGcount,soloMel_hypoDMR_close$aveCpGcount,soloIri_hypoDMR_open$aveCpGcount,soloIri_hypoDMR_other$aveCpGcount,soloIri_hypoDMR_close$aveCpGcount,openingDMAR_s24vsMel$aveCpGcount,closingDMAR_s24vsMel$aveCpGcount,openingDMAR_s24vsIri$aveCpGcount, closingDMAR_s24vsIri$aveCpGcount),c(soloMel_hypoDMR_open$CpGdensity100bp,soloMel_hypoDMR_other$CpGdensity100bp,soloMel_hypoDMR_close$CpGdensity100bp,soloIri_hypoDMR_open$CpGdensity100bp,soloIri_hypoDMR_other$CpGdensity100bp,soloIri_hypoDMR_close$CpGdensity100bp,openingDMAR_s24vsMel$CpGdensity100bp,closingDMAR_s24vsMel$CpGdensity100bp,openingDMAR_s24vsIri$CpGdensity100bp, closingDMAR_s24vsIri$CpGdensity100bp),c(rep("hypoDMR_open",nrow(soloMel_hypoDMR_open)),rep("hypoDMR_other",nrow(soloMel_hypoDMR_other)),rep("hypoDMR_other",nrow(soloMel_hypoDMR_close)),rep("hypoDMR_open",nrow(soloIri_hypoDMR_open)),rep("hypoDMR_other",nrow(soloIri_hypoDMR_other)),rep("hypoDMR_other",nrow(soloIri_hypoDMR_close)),rep("openingDMAR",nrow(openingDMAR_s24vsMel)),rep("closingDMAR",nrow(closingDMAR_s24vsMel)),rep("openingDMAR",nrow(openingDMAR_s24vsIri)),rep("closingDMAR",nrow(closingDMAR_s24vsIri))),c(rep("s24vMel",nrow(soloMel_hypoDMR_open)),rep("s24vMel",nrow(soloMel_hypoDMR_other)),rep("s24vMel",nrow(soloMel_hypoDMR_close)),rep("s24vIri",nrow(soloIri_hypoDMR_open)),rep("s24vIri",nrow(soloIri_hypoDMR_other)),rep("s24vIri",nrow(soloIri_hypoDMR_close)),rep("s24vMel",nrow(openingDMAR_s24vsMel)),rep("s24vMel",nrow(closingDMAR_s24vsMel)),rep("s24vIri",nrow(openingDMAR_s24vsIri)),rep("s24vIri",nrow(closingDMAR_s24vsIri))))
colnames(t) <-c("DMRsize","aveCpG","CpGdensity","type","comp")
t$type <- factor(t$type, levels = c("hypoDMR_open","hypoDMR_other","closingDMAR","openingDMAR"))

### Plot average CpG count in DM/ARs
p1<-ggplot(t, aes(x=type, y=aveCpG)) +
  geom_violin(aes(fill=type),trim=T,adjust =2)+geom_boxplot(fill="white",position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA,width = 0.2)+ggtitle("Ave CpG Count")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "# of CpG")+scale_fill_manual(values = c("#c8c7c7","#696a6a",mypalette[c(9,10)]))+scale_color_manual(values=c("#c8c7c7","#696a6a",mypalette[c(9,10)]))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p1+facet_grid(.~comp)

### Plot DMR size in DM/ARs
p2<-ggplot(t, aes(x=type, y=DMRsize)) +
  geom_violin(aes(fill=type),trim=T,adjust =2)+geom_boxplot(fill="white",position = position_dodge(0.9),outlier.shape=NA, outlier.size=NA,width = 0.2)+ggtitle("DMR size")+ 
  theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 14), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
  labs(x = "Sample", y = "DMR size")+scale_fill_manual(values = c("#c8c7c7","#696a6a",mypalette[c(9,10)]))+scale_color_manual(values=c("#c8c7c7","#696a6a",mypalette[c(9,10)]))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p2+facet_grid(.~comp)

# Make heatmap Mel specific DAR vs meth [Figure S4c-d]
library(ComplexHeatmap) ## For heatmap
library(circlize) ## For color options
name <- Only_DARs_Mel_specific_open[,4]
df.OG <- Only_DARs_Mel_specific_open[,c(21:23)]
row.names(df.OG) <- name

kclus <- kmeans(df.OG,8)
split <-  kclus$cluster
ht = Heatmap(df.OG, column_title = "Mel-specific opening soloDAR)",name= "Methylation",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = FALSE,split = split)
ht
## Use the same approach to generate methylation heatmap for Mel-specific closing soloDAR, Iri-specific opening soloDAR and Iri-specific closing soloDAR.

# Calculate the counts of differentially methylated and accessible regions [Numbers used in Figure 2a]
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

# Heatmap for DNA methylation levels of opening DARs identified in early NCC to late NCC transition [Figure 2b]
## Manual order of rows in heatmap 
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

All <- DMAR[DMAR$DARs15vs24>0 & DMAR$aveCpGcount > 0,]
name <- All[,4]
df.OG2 <- All[,c("Meth_s15","Meth_s24","Meth_Mel","Meth_Iri")]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("15somite", "24hpf","Mel", "Iri")

ht4 = Heatmap(df.OG2[roworder2,], column_title = "Methylation of 15somite to 24hpf DARs",name= "Methylaton",col = colorRamp2(c(0, 50, 100), c("#0392cf", "#fdf498","#ee4035")), cluster_rows=FALSE, cluster_columns = FALSE,show_row_names = FALSE)
ht4

# Bar chart illustrating number of DMRs, DARs and DEGs [Figure 1c-e]
mypalette <- brewer.pal(12,"Paired")
## Add countlist generating code
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
# Bar chart illustrating number of DMARs [Figure 2a]
## Left four pairs of comparisons
dmarplot <- DMAR_counts[DMAR_counts$Type == "DMAR" | DMAR_counts$Type == "DMAR2",]
dmarplot <- dmarplot[dmarplot$Name != "hypercloseDMAR" & dmarplot$Name != "hyperopenDMAR",]
dmarplot$Cell <- factor(dmarplot$Cell, levels = c("s15vs24","shared","s24vMel","s24vIri","MelvIri"))
dmarplot$Name <- factor(dmarplot$Name, levels = c("solohyperDMR","solohypoDMR","solocloseDAR","soloopenDAR","hypocloseDMAR","hypoopenDMAR"))

p <- ggplot(data=dmarplot, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DMAR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DMARs")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10)])+scale_y_continuous(lim = c(0,16000))+scale_color_manual(values=mypalette[c(1,2,3,4,9,10)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=dmarplot,position = position_dodge(width =0.9),aes(x=Cell,y=Count+400,label=Count),fontface='bold',hjust=0.5,size=5)
p
## Right one pair of comparison
dmarplot2 <- DMAR_counts[DMAR_counts$Type == "DMAR" | DMAR_counts$Type == "DMAR2",]
dmarplot2 <- dmarplot2[dmarplot2$Cell == "MelvIri",]
dmarplot2$Name <- factor(dmarplot2$Name, levels = c("solohyperDMR","solohypoDMR","solocloseDAR","soloopenDAR","hypocloseDMAR","hypoopenDMAR","hypercloseDMAR","hyperopenDMAR"))

p <- ggplot(data=dmarplot2, aes(x=Cell, y=Count, fill=Name)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("DMAR count")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Cell state", y = "Number of DMARs")+scale_fill_manual(values = mypalette[c(1,2,3,4,9,10,7,8)])+scale_y_continuous(lim = c(0,16000))+scale_color_manual(values=mypalette[c(1,2,3,4,9,10,7,8)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=dmarplot2,position = position_dodge(width =0.9),aes(x=Cell,y=Count+400,label=Count),fontface='bold',hjust=0.5,size=5)
p

# Annotations of each type of DMAR [Figure 2d]
Annotation_merged <-Reduce(function(x, y) merge(x, y,by = "Var1", all.x = T), list(as.data.frame(table(All_DMRs$Annotation)),
as.data.frame(table(All_DARs$Annotation)),
as.data.frame(table(Only_DMRs$Annotation)),
as.data.frame(table(Only_DARs$Annotation)),
as.data.frame(table(dynamic_DMARs$Annotation))))
Annotation_merged[is.na(Annotation_merged)]<-0
colnames(Annotation_merged) <- c("Annotation","DMRs","DARs","soloDMRs","soloDARs","dynamicDMARs")
str(Annotation_merged)

Annotation_merged$DMR <-  Annotation_merged$DMRs/sum(Annotation_merged$DMRs)*100
Annotation_merged$DAR <-  Annotation_merged$DARs/sum(Annotation_merged$DARs)*100
Annotation_merged$soloDMR <-  Annotation_merged$soloDMRs/sum(Annotation_merged$soloDMRs)*100
Annotation_merged$soloDAR <-  Annotation_merged$soloDARs/sum(Annotation_merged$soloDARs)*100
Annotation_merged$dynamicDMAR <-  Annotation_merged$dynamicDMARs/sum(Annotation_merged$dynamicDMARs)*100

Amelt <- melt(Annotation_merged[,c(1,7,8,9,10,11)],id.vars = "Annotation")
colnames(Amelt)[2]<- "Type"
Amelt$Type <- factor(Amelt$Type, levels = c("DMR","soloDMR","DAR","soloDAR","dynamicDMAR"))
Amelt$Annotation <- factor(Amelt$Annotation, levels = c("non-coding","Intergenic","promoter-TSS","5' UTR","exon","intron","3' UTR","TTS"))

mypalette4 <- c("#c67363","#de937c","#8ba6c8","#97c2d5","#edaa1e")
p <- ggplot(Amelt, aes(x=Annotation, y=value, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Annotation of Epigenetically Dynamic Regions")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Annotation", y = "Percent of peaks")+scale_fill_manual(values = mypalette4)+scale_color_manual(values=mypalette4)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title=element_text(size=16,face = "bold"),legend.text=element_text(size=12,face = "bold"))
p