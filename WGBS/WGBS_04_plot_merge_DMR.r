# Set up plotting basics
source("PlotSetUp.r")
# Create a dataframe for a union set of DMR collected from all pairs of collections. 
## Merge all samples d30
combined_dmr <- read.table("Combined_DMRs_d30_p0.01.txt",sep = "\t", header = F, stringsAsFactors = F)
## Plot DMR size distribution
combined_dmr$size <- combined_dmr$V3-combined_dmr$V2
hist(combined_dmr$size,breaks = 1000, xlim = c(0,2000),xlab = "DMR Size", main = "Combined DMR size distribution (d30)", col = "grey")
## Read DMR files
cDMR15v24_d30_p0.01 <- read.table("Combined_DMR15v24_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vMel_d30_p0.01 <- read.table("Combined_DMR24vMel_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMR24vIri_d30_p0.01 <- read.table("Combined_DMR24vIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
cDMRMelvIri_d30_p0.01 <- read.table("Combined_DMRMelvIri_d30_p0.01_wSMOOTHING.txt",sep = "\t", header = F, stringsAsFactors = F)
## Add a chromosome position column, like chr1:12345-67890
cDMR15v24_d30_p0.01$chrompos<- paste(cDMR15v24_d30_p0.01$V1,":",cDMR15v24_d30_p0.01$V2,"-",cDMR15v24_d30_p0.01$V3,sep = "")
cDMR24vMel_d30_p0.01$chrompos<- paste(cDMR24vMel_d30_p0.01$V1,":",cDMR24vMel_d30_p0.01$V2,"-",cDMR24vMel_d30_p0.01$V3,sep = "")
cDMR24vIri_d30_p0.01$chrompos<- paste(cDMR24vIri_d30_p0.01$V1,":",cDMR24vIri_d30_p0.01$V2,"-",cDMR24vIri_d30_p0.01$V3,sep = "")
cDMRMelvIri_d30_p0.01$chrompos<- paste(cDMRMelvIri_d30_p0.01$V1,":",cDMRMelvIri_d30_p0.01$V2,"-",cDMRMelvIri_d30_p0.01$V3,sep = "")
## Merge DMR files
combined_dmr_wInfo <-Reduce(function(x, y) merge(x, y, by=c("chrompos"), all = T), list(cDMR15v24_d30_p0.01[,c("V1","V2","V3","chrompos","V11","V12")],cDMR24vMel_d30_p0.01[,c("chrompos","V11","V12")],cDMR24vIri_d30_p0.01[,c("chrompos","V11","V12")],cDMRMelvIri_d30_p0.01[,c("chrompos","V11","V12")]))
### Change column names to match the data
colnames(combined_dmr_wInfo) <- c("chrompos","chr","start","end","s15vs24","s24vMel","s24vIri","MelvsIri")
### Calculate region sizes
combined_dmr_wInfo$size <- combined_dmr_wInfo$end-combined_dmr_wInfo$start
### Replace "." which stands for no data with 0
combined_dmr_wInfo[combined_dmr_wInfo == "."] <- 0
### Change data type to numeric
combined_dmr_wInfo$s15vs24 <- as.numeric(combined_dmr_wInfo$s15vs24)
combined_dmr_wInfo$s24vMel <- as.numeric(combined_dmr_wInfo$s24vMel)
combined_dmr_wInfo$s24vIri <- as.numeric(combined_dmr_wInfo$s24vIri)
combined_dmr_wInfo$MelvsIri <- as.numeric(combined_dmr_wInfo$MelvsIri)

## Remove duplicated regions from the data frame. 
combined_dmr_wInfo <- combined_dmr_wInfo[!duplicated(combined_dmr_wInfo$chrompos),]
write.table(combined_dmr_wInfo, "Combined_DMRs_d30_p0.01_wINFO.bed", row.names = F, col.names = F, sep = "\t",quote =F)

# Plot methylation differences and also correlation in shared DMRs
## Extract sample-associated DMR sets
s24 <-combined_dmr_wInfo[combined_dmr_wInfo$s15vs24 !=0,] 
M<-combined_dmr_wInfo[combined_dmr_wInfo$s24vMel !=0,]
I<- combined_dmr_wInfo[combined_dmr_wInfo$s24vIri !=0,]
MI <-combined_dmr_wInfo[combined_dmr_wInfo$MelvsIri !=0,]
## Density plot for methylation change distribution
ps24 <-ggplot(s24, aes(s15vs24)) + geom_density(alpha = 0.5,adjust =0.15,color = mypalette[3],fill = mypalette[3])+ggtitle(paste("DMR methylation change from 15somite to 24hpf"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
ps24
pM <-ggplot(M, aes(s24vMel)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[5],fill = mypalette[5])+ggtitle(paste("DMR methylation change from 24hpf to Melanophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM
pI <-ggplot(I, aes(s24vIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[7],fill = mypalette[7])+ggtitle(paste("DMR methylation change from 24hpf to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI
pMI <-ggplot(MI, aes(MelvsIri)) + geom_density(alpha = 0.5,adjust =1,color = mypalette[9],fill = mypalette[9])+ggtitle(paste("DMR methylation change from Melanophore to Iridophore"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Methylation change", y = "Density")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI
multiplot(ps24,pI,pM,pMI, cols = 2)

# DMR size distribution
ps24 <-ggplot(s24, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[4],fill = mypalette[4])+ggtitle(paste("DMR size distribution\n(15somite vs 24hpf)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pM <-ggplot(M, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[6],fill = mypalette[6])+ggtitle(paste("DMR size distribution\n(24hpf vs Melanophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pI <-ggplot(I, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[8],fill = mypalette[8])+ggtitle(paste("DMR size distribution\n(24hpf vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
pMI <-ggplot(MI, aes(size)) + geom_density(alpha = 0.5,adjust =0.4,color = mypalette[10],fill = mypalette[10])+ggtitle(paste("DMR size distribution\n(Melanophore vs Iridophore)"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
    labs(x = "DMR size", y = "Density")+scale_x_continuous(lim = c(0,2000))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.position="none")
multiplot(ps24,pI,pM,pMI, cols = 2)
# Features of shared and cell type-specific DMRs 
## Pull out melanophore-specific
Mel_specific_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <0 & combined_dmr_wInfo$s24vIri >=0,] 
Mel_specific_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri <=0,] 
## Using more strigent way to call Mel hypoDMR
Mel_specific_hypoDMR_stringent <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri <=0 & combined_dmr_wInfo$MelvsIri <0,] 
## Pull out iridophore-specific
Iri_specific_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >=0 & combined_dmr_wInfo$s24vIri <0,] 
Iri_specific_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <=0 & combined_dmr_wInfo$s24vIri >0,] 
# Using more strigent way to call Iri hypoDMR
Iri_specific_hypoDMR_stringent <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <=0 & combined_dmr_wInfo$s24vIri >0 & combined_dmr_wInfo$MelvsIri >0,] 
## Pull out shared DMRs
Mel_Iri_shared_hyperDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel <0 & combined_dmr_wInfo$s24vIri <0,] 
Mel_Iri_shared_hypoDMR <- combined_dmr_wInfo[combined_dmr_wInfo$s24vMel >0 & combined_dmr_wInfo$s24vIri >0,] 

## Write hypo DMRs to files
write.table(Mel_specific_hyperDMR, "Mel_specific_hyperDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_hypoDMR, "Mel_specific_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_specific_hypoDMR_stringent, "Mel_specific_hypoDMR_d30_p0.01_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Iri_specific_hyperDMR, "Iri_specific_hyperDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_hypoDMR, "Iri_specific_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Iri_specific_hypoDMR_stringent, "Iri_specific_hypoDMR_d30_p0.01_stringent.bed", row.names = F, col.names = F, sep = "\t",quote =F)

write.table(Mel_Iri_shared_hyperDMR , "Mel_Iri_shared_hyperDMR _d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)
write.table(Mel_Iri_shared_hypoDMR, "Mel_Iri_shared_hypoDMR_d30_p0.01.bed", row.names = F, col.names = F, sep = "\t",quote =F)

## Plot DMR methylation change comparison
Mel_Iri_shared_hyperDMR$Type <- "Shared HyperDMR"
Mel_Iri_shared_hypoDMR$Type <- "Shared HypoDMR"
Mel_specific_hyperDMR$Type <- "Mel-specific HyperDMR"
Mel_specific_hypoDMR$Type <- "Mel-specific HypoDMR"
Iri_specific_hyperDMR$Type <- "Iri-specific HyperDMR"
Iri_specific_hypoDMR$Type <- "Iri-specific HypoDMR"

shared <- rbind(Mel_Iri_shared_hyperDMR,Mel_Iri_shared_hypoDMR,Mel_specific_hyperDMR,Mel_specific_hypoDMR[Mel_specific_hypoDMR$s24vMel > 0.2 | Mel_specific_hypoDMR$s24vMel < -0.2,],Iri_specific_hyperDMR,Iri_specific_hypoDMR)
ggplot(shared[shared$type == "Shared hyperDMR" | shared$type == "Shared hypoDMR",], aes(x=s24vMel, y=s24vIri,fill=Type,color=Type)) +
  geom_point(size = 0.2) + geom_jitter(data = shared[shared$type != "Shared hyperDMR" | shared$type != "Shared hypoDMR",],aes(x=s24vMel, y=s24vIri,fill=Type,color=Type),width = 0.01, height = 0.01,size = 0.2)+
  ggtitle(paste("DMR methylation change comparison from 24hpf to Pigment cells"))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 20), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      scale_fill_manual(values = mypalette3)+scale_color_manual(values=mypalette3)+labs(x = "Methylaton change from 24hpf to Melanophore", y = "Methylaton change from 24hpf to Iridophore")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+ theme(legend.title = element_text(size = 12, face = "bold"),legend.text = element_text(size = 10, face = "bold"))+ theme(legend.position="top")+ guides(color=guide_legend(override.aes=list(fill=NA,size = 6)))+ theme(legend.key=element_blank())


