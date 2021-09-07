# Import DAR DMR info file
DMAR_DAR<- read.table("Combined_DMR_DAR_wDARinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_DMR<- read.table("Combined_DMR_DAR_wDMRinfo.bed",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR15<- read.table("Combined_DMR_DAR_wIDRinfo_15s.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDR24<- read.table("Combined_DMR_DAR_wIDRinfo_24hpf.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRM<- read.table("Combined_DMR_DAR_wIDRinfo_Mel.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)
DMAR_IDRI<- read.table("Combined_DMR_DAR_wIDRinfo_Iri.txt",sep = "\t", header = F,quote = "", stringsAsFactors = F)

# Add chrom position column
DMAR_DAR$chrompos<- paste(DMAR_DAR$V1,":",DMAR_DAR$V2,"-",DMAR_DAR$V3,sep = "")
DMAR_DMR$chrompos<- paste(DMAR_DMR$V1,":",DMAR_DMR$V2,"-",DMAR_DMR$V3,sep = "")
DMAR_IDR15$chrompos<- paste(DMAR_IDR15$V1,":",DMAR_IDR15$V2,"-",DMAR_IDR15$V3,sep = "")
DMAR_IDR24$chrompos<- paste(DMAR_IDR24$V1,":",DMAR_IDR24$V2,"-",DMAR_IDR24$V3,sep = "")
DMAR_IDRM$chrompos<- paste(DMAR_IDRM$V1,":",DMAR_IDRM$V2,"-",DMAR_IDRM$V3,sep = "")
DMAR_IDRI$chrompos<- paste(DMAR_IDRI$V1,":",DMAR_IDRI$V2,"-",DMAR_IDRI$V3,sep = "")

# Methylation value, filter out DMRs with less than ave 5 cov
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
DMAR$DMRs15vs24 <- as.numeric(DMAR$DMRs15vs24)
DMAR$DMRs24vMel <- as.numeric(DMAR$DMRs24vMel)
DMAR$DMRs24vIri <- as.numeric(DMAR$DMRs24vIri)
DMAR$DMRMelvIri <- as.numeric(DMAR$DMRMelvIri)
DMAR$DARs15vs24 <- as.numeric(DMAR$DARs15vs24)
DMAR$DARs24vMel <- as.numeric(DMAR$DARs24vMel)
DMAR$DARs24vIri <- as.numeric(DMAR$DARs24vIri)
DMAR$DARMelvIri <- as.numeric(DMAR$DARMelvIri)

DMAR <- DMAR[!duplicated(DMAR$chrompos),]

DMAR <- DMAR[,c(2,3,4,1,31,5:30)]
DMAR$aveCpGcount <- apply(DMAR,1,function(x) {(as.numeric(x[24])+as.numeric(x[25])+as.numeric(x[26])+as.numeric(x[27]))/4})
DMAR$CpGdensity100bp <- apply(DMAR,1,function(x) {(as.numeric(x[32])/as.numeric(x[5]))*100})
DMAR$DARMelvIri <- DMAR$DARMelvIri*-1   ######### the directionality of DAR is switch for Mel v Iri for DARs!!!!###################

write.table(DMAR, "/scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/Combined_analysis/All_DMAR_Combined_wINFO_NEW_042620.bed", row.names = F, col.names = T, sep = "\t",quote =F)
