### Tool version
#- R v3.3.0: Used to analyze sequencing library results
#- DSS v2.14.0 :Used to call differentially methylated regions
#- DESeq2 v1.12.4: Used to call differentially expressed genes
#- DiffBind v 2.2.12: Used to call differentially accessible regions/peaks

####QC for individual WGBS Replicates####
setwd(<WGBS_folder>) # use the directory where outputs from convert_DSS.py are saved as working directory
# Setting up plotting basics
source("PlotSetUp.r")
## Load methylation profile 
s15_R1 <- read.table("15somite_NCC_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s15_R2 <- read.table("15somite_NCC_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s24_R1 <- read.table("24hpf_NCC_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
s24_R2 <- read.table("24hpf_NCC_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Mel_R1 <- read.table("Mel_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Mel_R2 <- read.table("Mel_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Iri_R1 <- read.table("Iri_Rep1_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
Iri_R2 <- read.table("Iri_Rep2_DSS.txt",sep = "\t", header = T, stringsAsFactors = F)
## Calculate methylation level at each CpG site, and round the resulting ratio to 2 digits to the right of decimal.
s15_R1$M <- round(s15_R1$X/s15_R1$N, digits = 2)
s15_R2$M <- round(s15_R2$X/s15_R2$N, digits = 2)
s24_R1$M <- round(s24_R1$X/s24_R1$N, digits = 2)
s24_R2$M <- round(s24_R2$X/s24_R2$N, digits = 2)
Mel_R1$M <- round(Mel_R1$X/Mel_R1$N, digits = 2)
Mel_R2$M <- round(Mel_R2$X/Mel_R2$N, digits = 2)
Iri_R1$M <- round(Iri_R1$X/Iri_R1$N, digits = 2)
Iri_R2$M <- round(Iri_R2$X/Iri_R2$N, digits = 2)

# Averge methylation
filter = 5 # Using 5 as cutoff. Only include CpGs with a coverage no less than 5.
## Calculate average methylation for each dataset
me <- data.frame(c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"),c(mean(s15_R1[s15_R1$N >= filter,]$M),mean(s15_R2[s15_R2$N >= filter,]$M),mean(s24_R1[s24_R1$N >= filter,]$M),mean(s24_R2[s24_R2$N >= filter,]$M),mean(Mel_R1[Mel_R1$N >= filter,]$M),mean(Mel_R2[Mel_R2$N >= filter,]$M),mean(Iri_R1[Iri_R1$N >= filter,]$M),mean(Iri_R2[Iri_R2$N >= filter,]$M)))
colnames(me) <- c("Type", "Methylation")
me$Type <- factor(me$Type, levels = c("15somite Rep1","15somite Rep2","24hpf Rep1","24hpf Rep2","Mel Rep1","Mel Rep2","Iri Rep1","Iri Rep2"))
p <- ggplot(data=me, aes(x=Type, y=Methylation, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Average Methylation across samples (coverage filter >=5)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "Average Methylation")+scale_y_continuous(lim = c(0,1))+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=me,aes(x=Type,y=Methylation+0.03,label=paste(round(Methylation,3)*100,"%",sep ="")),size=8)+ theme(legend.position="none")
p

# Number of CpG covered 
## Calculate number of CpGs covered in each sample with coverage cutoff at 0,5,10
co <- data.frame(c(rep("Mel-R1", 3),rep("Mel-R2", 3),rep("Iri-R1", 3),rep("Iri-R2", 3),rep("24hpf-R1",3),rep("24hpf-R2", 3),rep("s15-R1", 3),rep("s15-R2", 3)),c(nrow(Mel_R1[Mel_R1$N>0,]), nrow(Mel_R1[Mel_R1$N>=5,]), nrow(Mel_R1[Mel_R1$N>=10,]),nrow(Mel_R2[Mel_R2$N>0,]), nrow(Mel_R2[Mel_R2$N>=5,]), nrow(Mel_R2[Mel_R2$N>=10,]),nrow(Iri_R1[Iri_R1$N>0,]), nrow(Iri_R1[Iri_R1$N>=5,]), nrow(Iri_R1[Iri_R1$N>=10,]),nrow(Iri_R2[Iri_R2$N>0,]), nrow(Iri_R2[Iri_R2$N>=5,]), nrow(Iri_R2[Iri_R2$N>=10,]),nrow(s24_R1[s24_R1$N>0,]), nrow(s24_R1[s24_R1$N>=5,]), nrow(s24_R1[s24_R1$N>=10,]),nrow(s24_R2[s24_R2$N>0,]), nrow(s24_R2[s24_R2$N>=5,]), nrow(s24_R2[s24_R2$N>=10,]),nrow(s15_R1[s15_R1$N>0,]), nrow(s15_R1[s15_R1$N>=5,]), nrow(s15_R1[s15_R1$N>=10,]),nrow(s15_R2[s15_R2$N>0,]), nrow(s15_R2[s15_R2$N>=5,]), nrow(s15_R2[s15_R2$N>=10,])),c(rep(c(">0",">=5",">=10"),8)))
colnames(co) <- c("Type", "count","Filter")
co$Type <- factor(co$Type, levels = c("s15-R1","s15-R2","24hpf-R1","24hpf-R2","Mel-R1","Mel-R2","Iri-R1","Iri-R2"))
co$Filter <- factor(co$Filter, levels = c(">0",">=5",">=10"))
p <- ggplot(data=co, aes(x=Type, y=count, fill=Filter)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle(paste("CpG coverage cutoff",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "CpG count")+scale_fill_manual(values = mypalette2)+scale_color_manual(values=mypalette2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p

# Methylation distribution
## Only include CpGs with a coverage no less than 5.
filter = 5
meth <- data.frame(c(rep("Mel-R1", nrow(Mel_R1[Mel_R1$N >= filter,])),rep("Mel-R2", nrow(Mel_R2[Mel_R2$N >= filter,])),rep("Iri-R1", nrow(Iri_R1[Iri_R1$N >= filter,])),rep("Iri-R2", nrow(Iri_R2[Iri_R2$N >= filter,])),rep("24hpf-R1", nrow(s24_R1[s24_R1$N >= filter,])),rep("24hpf-R2", nrow(s24_R2[s24_R2$N >= filter,])),rep("s15-R1", nrow(s15_R1[s15_R1$N >= filter,])),rep("s15-R2", nrow(s15_R2[s15_R2$N >= filter,]))), c(Mel_R1[Mel_R1$N >= filter,]$M,Mel_R2[Mel_R2$N >= filter,]$M,Iri_R1[Iri_R1$N >= filter,]$M,Iri_R2[Iri_R2$N >= filter,]$M,s24_R1[s24_R1$N >= filter,]$M,s24_R2[s24_R2$N >= filter,]$M,s15_R1[s15_R1$N >= filter,]$M,s15_R2[s15_R2$N >= filter,]$M))
colnames(meth) <- c("Type", "Meth")
a <-ggplot(meth, aes(Meth, fill = Type, colour = Type)) + geom_density(alpha = 0.05, adjust = 3)+ggtitle(paste("Pigment CpG methylation distribution (coverage filter >=",filter,")",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Coverage", y = "Density")+scale_fill_manual(values = mypalette)+scale_color_manual(values=mypalette)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
a

# Principal component analysis 
## Only include CpGs with a coverage no less than 5.
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

# Plot samples on PC axes 
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

# Plot pair-wise correlation between samples
gsc_cor <- cor(m[,3:10])
corrplot(gsc_cor, order="hclust", addrect=4, method ="number")
corrplot.mixed(gsc_cor, order="hclust", addrect=4)
corrplot(gsc_cor, order="hclust", addrect=4)


#### QC for replicate-combined WGBS ####
# Combine replicates
## Combine 15somite NCC Rep1 with Rep2
s15pos_combined <- merge(s15_R1,s15_R2,by=c("chr","pos"), all = TRUE)
s15pos_combined[is.na(s15pos_combined)]<-0
s15pos_combined$N <- s15pos_combined$N.x+s15pos_combined$N.y
s15pos_combined$X <- s15pos_combined$X.x+s15pos_combined$X.y
s15pos_combined <-s15pos_combined[,c(1,2,9,10)]
write.table(s15pos_combined, "15somite_NCC_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

## Combine 24hpf NCC Rep1 with Rep2
s24pos_combined <- merge(s24_R1,s24_R2,by=c("chr","pos"), all = TRUE)
s24pos_combined[is.na(s24pos_combined)]<-0
s24pos_combined$N <- s24pos_combined$N.x+s24pos_combined$N.y
s24pos_combined$X <- s24pos_combined$X.x+s24pos_combined$X.y
s24pos_combined <-s24pos_combined[,c(1,2,9,10)]
write.table(s24pos_combined, "24hpf_NCC_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

## Combine Mel Rep1 with Rep2
Mel_combined <- merge(Mel_R1,Mel_R2,by=c("chr","pos"), all = TRUE)
Mel_combined[is.na(Mel_combined)]<-0
Mel_combined$N <- Mel_combined$N.x+Mel_combined$N.y
Mel_combined$X <- Mel_combined$X.x+Mel_combined$X.y
Mel_combined <-Mel_combined[,c(1,2,9,10)]
write.table(Mel_combined, "Mel_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

## Combine Iri Rep1 with Rep2
Iri_combined <- merge(Iri_R1,Iri_R2,by=c("chr","pos"), all = TRUE)
Iri_combined[is.na(Iri_combined)]<-0
Iri_combined$N <- Iri_combined$N.x+Iri_combined$N.y
Iri_combined$X <- Iri_combined$X.x+Iri_combined$X.y
Iri_combined <-Iri_combined[,c(1,2,9,10)]
write.table(Iri_combined, "Iri_Combined_DSS.txt",col.names = T,row.names = F, sep = "\t", quote = F)

# Calculate methylation level and round that to two digits to the right of decimal
s15pos_combined$M <- round(s15pos_combined$X/s15pos_combined$N, digits = 2)
s24pos_combined$M <- round(s24pos_combined$X/s24pos_combined$N, digits = 2)
Mel_combined$M <- round(Mel_combined$X/Mel_combined$N, digits = 2)
Iri_combined$M <- round(Iri_combined$X/Iri_combined$N, digits = 2)

# Average methylation
## Only include CpGs with a coverage no less than 5.
filter = 5 
me <- data.frame(c("15somite","24hpf","Mel","Iri"),c(mean(s15pos_combined[s15pos_combined$N >= filter,]$M),mean(s24pos_combined[s24pos_combined$N >= filter,]$M),mean(Mel_combined[Mel_combined$N >= filter,]$M),mean(Iri_combined[Iri_combined$N >= filter,]$M)))
colnames(me) <- c("Type", "Methylation")
me$Type <- factor(me$Type, levels = c("15somite","24hpf","Mel","Iri"))
# Plot average methylation
p1 <- ggplot(data=me, aes(x=Type, y=Methylation, fill=Type)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle("Average Methylation across samples (coverage filter >=5)")+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "Average Methylation")+scale_y_continuous(lim = c(0,1))+scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))+geom_text(data=me,aes(x=Type,y=Methylation+0.03,label=paste(round(Methylation,3)*100,"%",sep ="")),size=8)+ theme(legend.position="none")
p1

# Number of CpG covered 
## Coverage cutoff at 0, 5, 10.
co <- data.frame(c(rep("Mel", 3),rep("Iri", 3),rep("24hpf", 3),rep("s15", 3)),c(nrow(Mel_combined[Mel_combined$N>0,]), nrow(Mel_combined[Mel_combined$N>=5,]), nrow(Mel_combined[Mel_combined$N>=10,]),nrow(Iri_combined[Iri_combined$N>0,]), nrow(Iri_combined[Iri_combined$N>=5,]), nrow(Iri_combined[Iri_combined$N>=10,]),nrow(s24pos_combined[s24pos_combined$N>0,]), nrow(s24pos_combined[s24pos_combined$N>=5,]), nrow(s24pos_combined[s24pos_combined$N>=10,]),nrow(s15pos_combined[s15pos_combined$N>0,]), nrow(s15pos_combined[s15pos_combined$N>=5,]), nrow(s15pos_combined[s15pos_combined$N>=10,])),c(rep(c(">0",">=5",">=10"),4)))
colnames(co) <- c("Type", "count","Filter")
co$Type <- factor(co$Type, levels = c("s15","24hpf","Mel","Iri"))
co$Filter <- factor(co$Filter, levels = c(">0",">=5",">=10"))
# Plot CpG coverage cutoff
p2 <- ggplot(data=co, aes(x=Type, y=count, fill=Filter)) +geom_bar(stat="identity",colour = "black", position=position_dodge())+ggtitle(paste("CpG coverage cutoff",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Library", y = "CpG count")+scale_fill_manual(values = mypalette2)+scale_color_manual(values=mypalette2)+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p2

# Methylation distribution
## Only include CpGs with a coverage no less than 5.
filter = 5
meth <- data.frame(c(rep("Mel", nrow(Mel_combined[Mel_combined$N >= filter,])),rep("Iri", nrow(Iri_combined[Iri_combined$N >= filter,])),rep("24hpf", nrow(s24pos_combined[s24pos_combined$N >= filter,])),rep("s15", nrow(s15pos_combined[s15pos_combined$N >= filter,]))), c(Mel_combined[Mel_combined$N >= filter,]$M,Iri_combined[Iri_combined$N >= filter,]$M,s24pos_combined[s24pos_combined$N >= filter,]$M,s15pos_combined[s15pos_combined$N >= filter,]$M))
colnames(meth) <- c("Type", "Meth")
meth$Type <- factor(meth$Type, levels = c("s15","24hpf","Mel","Iri"))
# Plot CpG methylation distribution
p3 <-ggplot(meth, aes(Meth, fill = Type, colour = Type)) + geom_density(alpha = 0.1, adjust = 4)+ggtitle(paste("Pigment CpG methylation distribution (coverage filter >=",filter,")",sep =""))+theme(plot.title = element_text(hjust = 0.5, face = "bold",size = 16), axis.title=element_text(size=12,face = "bold"),axis.text.x = element_text(face = "bold",size = 12),axis.text.y = element_text(face = "bold",size = 12))+
      labs(x = "Coverage", y = "Density")+scale_fill_manual(values = mypalette[c(1,3,5,7)])+scale_color_manual(values=mypalette[c(1,3,5,7)])+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(panel.background = element_rect(fill = 'white',colour = 'black'))
p3

# PCA
## Only include CpGs with a coverage no less than 5.
filter = 5
m<-Reduce(function(x, y) merge(x, y, by=c("chr","pos")), list(Mel_combined[Mel_combined$N >= filter,c(1,2,5)],Iri_combined[Iri_combined$N >= filter,c(1,2,5)],s24pos_combined[s24pos_combined$N >= filter,c(1,2,5)],s15pos_combined[s15pos_combined$N >= filter,c(1,2,5)]))
names(m) <- c("chr","start","Mel","Iri","24hpf","s15")
m.pca <-prcomp(m[,3:6])
summary(m.pca)

#Importance of components:
#                          PC1     PC2     PC3     PC4
#Standard deviation     0.5325 0.12170 0.07764 0.06329
#Proportion of Variance 0.9194 0.04802 0.01954 0.01299
#Cumulative Proportion  0.9194 0.96747 0.98701 1.00000

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
# Sample correlation plot
gsc_cor <- cor(m[,3:6])
corrplot(gsc_cor, order="hclust", addrect=4, method ="number")
corrplot.mixed(gsc_cor, order="hclust", addrect=4)
corrplot(gsc_cor, order="hclust", addrect=4)


#### Differential DNA methylation analysis with DSS####
# Load package for differential analysis
library(DSS)
# Create an object of BSseq class, 
BSobj <-makeBSseqData(list(s15pos_combined,s24pos_combined,Mel_combined,Iri_combined), c("s15","24hpf","Mel","Iri"))

# Conduct differential methylation loci (DML) tests between each two samples
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
