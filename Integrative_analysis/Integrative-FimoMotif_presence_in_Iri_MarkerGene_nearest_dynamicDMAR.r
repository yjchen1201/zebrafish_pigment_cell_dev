####### USE FIMO TO scan FOR IRI-SPECIFIC TF motifs ###########
# Read DMAR list generated from "Integrative-generate_DMAR_master_list_02.r"
DMAR <- read.table("All_DMAR_Combined_wINFO.bed",sep = "\t",header =T, stringsAsFactors =F)
DMAR_Iri_accessible <- DMAR[DMAR_new$DMRs24vIri  >0 | DMAR_new$DARs24vIri >0 | DMAR_new$DMRMelvIri  >0 | DMAR_new$DARMelvIri <0,]
write.table(DMAR_Iri_accessible[,c(1,2,3,4)], "DMAR_Iri_Accesible_for_FIMO.bed", sep = "\t", quote = F, col.names = F, row.names =F)

#DMAR_Mel_accessible <- DMAR_new[DMAR_new$DMRs24vMel  >0 | DMAR_new$DARs24vMel >0 | DMAR_new$DMRMelvIri  <0 | DMAR_new$DARMelvIri >0,]
#write.table(DMAR_Mel_accessible[,c(1,2,3,4)], "DMAR_Mel_Accesible_for_FIMO.bed", sep = "\t", quote = F, col.names = F, row.names =F)

# Use UCSC twoBitToFa to convert regions in bed file to fasta
system("module load kentUCSC")
system("twoBitToFa -bed=DMAR_Iri_Accesible_for_FIMO.bed danRer10.2bit DMAR_Iri_Accesible_for_FIMO.fasta")

# Use fimo to scan for Iri TFs
system("module meme")
system("fimo Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta")
#system("mast Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta")
#system("mcast Iri_TFs_FIMO.meme DMAR_Iri_Accesible_for_FIMO.fasta")

# Use Parse_fimo.py code to parse out motif presence in bash
system("module load python3")
system("python3 Parse_fimo.py fimo.txt") #fimo.txt is the output from fimo step

# Read the output generated using Parsed_fimo.py
FIMO <- read.table("fimo_parsed.txt",sep = "\t", header = T, quote = "", stringsAsFactors =F)

# merge DMAR_Iri_accessible table with FIMO table by chrompos column
DMAR_Iri_FIMO <- merge(DMAR_Iri_accessible[,c(1:32)],FIMO, by = "chrompos", all.x = T)
DMAR_Iri_FIMO[is.na(DMAR_Iri_FIMO)] <- 0

# Save Iri accessibleDMAR with FIMO scanned motifs
#write.table(DMAR_Iri_FIMO[,c(2:4,1,5:52)], "accessibleDMAR_Iri_wFIMOmotif.txt", sep = "\t", col.names = T, row.names = F, quote =F)
write.table(DMAR_Iri_FIMO[,c(2:4,1,5:52)], "accessibleDMAR_Iri_wFIMOmotif_noheader.txt", sep = "\t", col.names = F, row.names = F, quote =F)


# Input file "Danio_rerio.GRCz10.85.GENE.PROMOTER.bed" is generated from "Integrative-promoter-centric.r"
# import prmoter information
prom <- read.table("Danio_rerio.GRCz10.85.GENE.PROMOTER.bed", sep = "\t", stringsAsFactors =F, header = F, quote ="")
# extract Iri marker genee promoter info out
Iri_MarkerGenes <- prom[prom$gene == "ENSDARG00000003732"|prom$gene =="ENSDARG00000002933"|prom$gene == "ENSDARG00000016706" |prom$gene == "ENSDARG00000077467" | prom$gene == "ENSDARG00000098745" |prom$gene == "ENSDARG00000003732" |prom$gene == "ENSDARG00000059279" |prom$gene == "ENSDARG00000008861" | prom$gene == "ENSDARG00000021032"| prom$gene == "ENSDARG00000024431",]

# Read DEG info
## "DEG_p0.01_ENSEMBL_to_Gene_Name.txt" and "DEGs_combined_samples_p0.01_wTFinfo.txt" are from RNA/RNA_02_DEGanalysis.r
DEG2_genenames<- read.table("DEG_p0.01_ENSEMBL_to_Gene_Name.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) 
DEG2 <- read.table("DEGs_combined_samples_p0.01_wTFinfo.txt",header = T, ,quote = "", sep = "\t",stringsAsFactors = F) 
colnames(DEG2_genenames) <- c("gene","genename")
DEG3 <- merge(DEG2,DEG2_genenames, by = "gene", all.x = T)

# Add DEG info to Iri Marker gene table
Iri_MarkerGenes2 <- merge(Iri_MarkerGenes,DEG3[,c(1:15)],by = "gene", all.x = T)
Iri_MarkerGenes2[is.na(Iri_MarkerGenes2)] <- 0
Iri_MarkerGenes2[Iri_MarkerGenes2$gene == "ENSDARG00000021032",]$genename.x <- "foxd3"

# Extend promoter region to 50kb on both side
Iri_MarkerGenes2$start50kb <- Iri_MarkerGenes2$start-50000
Iri_MarkerGenes2$end50kb <- Iri_MarkerGenes2$end+50000
# Correct the start to 0 if it is <0
Iri_MarkerGenes2[Iri_MarkerGenes2$start50kb < 0,]$start50kb <-0

# Save Iri marker gene promoter extended information table
write.table(Iri_MarkerGenes2[,c(2,22,23,1,14)],"Iri_MarkerGenes_promoter_50kb_extended.bed", sep = "\t", col.names = F, row.names = F, quote =F)


# Use bedtools to intersect Iri markergene promoter extended regions with Iri accessible DMARs with FIMO motifs
system("module load bedtools")
system("bedtools intersect -wao -a Iri_MarkerGenes_promoter_50kb_extended.bed -b accessibleDMAR_Iri_wFIMOmotif_noheader.txt > Iri_MarkerGenes_promoter_50kb_extended_wDMAR_FIMO.txt")

# Import the intersected file
Iri_Marker_50kb_DMAR_FIMO <- read.table("Iri_MarkerGenes_promoter_50kb_extended_wDMAR_FIMOMotif.txt", header = F, stringsAsFactors = F, sep = "\t", quote = "")
colnames(Iri_Marker_50kb_DMAR_FIMO) <- c("chr","start50kb","end50kb","gene","genename",colnames(DMAR_Iri_FIMO),"overlap")
colnames(Iri_Marker_50kb_DMAR_FIMO)[6:9] <- c("chr","start","end","chrompos")
Iri_Marker_50kb_DMAR_FIMO[Iri_Marker_50kb_DMAR_FIMO == "."]<-0
Iri_Marker_50kb_DMAR_FIMO[,c(38:45)] <- sapply(Iri_Marker_50kb_DMAR_FIMO[,c(38:45)],as.numeric)
Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$genename == "alx1" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx3"| Iri_DEG_50kb_DMAR_FIMO$genename == "alx4a" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx4b" | Iri_DEG_50kb_DMAR_FIMO$genename == "ltk"| Iri_DEG_50kb_DMAR_FIMO$genename == "ednrba" | Iri_DEG_50kb_DMAR_FIMO$genename == "pnp4a",]
Iri_Marker_50kb_DMAR_FIMO_ALL <- rbind(Iri_Marker_50kb_DMAR_FIMO,Iri_DEG_50kb_DMAR_FIMO[Iri_DEG_50kb_DMAR_FIMO$genename == "gbx2" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx1" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx3"| Iri_DEG_50kb_DMAR_FIMO$genename == "alx4a" | Iri_DEG_50kb_DMAR_FIMO$genename == "alx4b" | Iri_DEG_50kb_DMAR_FIMO$genename == "ltk"| Iri_DEG_50kb_DMAR_FIMO$genename == "ednrba"| Iri_DEG_50kb_DMAR_FIMO$genename == "pnp4a",])
Iri_Marker_50kb_DMAR_FIMO_ALL <- merge(Iri_Marker_50kb_DMAR_FIMO_ALL[Iri_Marker_50kb_DMAR_FIMO_ALL$genename == "tfap2a" | Iri_Marker_50kb_DMAR_FIMO_ALL$GBX2 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ALX1 > 0 | Iri_Marker_50kb_DMAR_FIMO_ALL$ALX3 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ALX4 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$SOX10 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$ETS1 > 0 |Iri_Marker_50kb_DMAR_FIMO_ALL$TFEC > 0,],DEG3[,c(0:15)], by = "gene", all.x = T)
Iri_Marker_50kb_DMAR_FIMO_ALL[is.na(Iri_Marker_50kb_DMAR_FIMO_ALL)] <- 0

#Iri TF specific#
Iri_Marker_50kb_DMAR_FIMO_TF <- Iri_Marker_50kb_DMAR_FIMO_ALL[Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "ltk" & Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "atic" & Iri_Marker_50kb_DMAR_FIMO_ALL$genename != "ednrba",]
Iri_Marker_50kb_DMAR_FIMO_TF<- Iri_Marker_50kb_DMAR_FIMO_TF[order(Iri_Marker_50kb_DMAR_FIMO_TF$genename),]
Iri_Marker_50kb_DMAR_FIMO_TF <- Iri_Marker_50kb_DMAR_FIMO_TF[Iri_Marker_50kb_DMAR_FIMO_TF$ALX1 >0 |Iri_Marker_50kb_DMAR_FIMO_TF$ALX3 >0 |Iri_Marker_50kb_DMAR_FIMO_TF$ALX4 >0 |Iri_Marker_50kb_DMAR_FIMO_TF$GBX2>0 |Iri_Marker_50kb_DMAR_FIMO_TF$SOX10 >0 |Iri_Marker_50kb_DMAR_FIMO_TF$ETS1 >0 |Iri_Marker_50kb_DMAR_FIMO_TF$TFEC >0,]

######## Generate supFig7c heatmap plot ###########
name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,9])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(38:44)]
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(38:45)])
ht2 = Heatmap(df.OG2, column_title = "Iri TFs DMARs\n(50kb from DEG promoter)",name= "Motif\nOccurrence",col = colorRamp2(c(0, 2, 4), c("#FFFEF2","#17C0F2","#0084A1")),row_names_side = "left",row_names_max_width = unit(20, "cm"),
    cluster_rows = F, cluster_columns = F,show_row_names = T,clustering_method_rows= "ward.D",row_dend_side = c("right"),)
ht2
name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- (Iri_Marker_50kb_DMAR_FIMO_TF[,c(35)])
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(35)])
ht = Heatmap(df.OG2,name= "Annotation",col = c("#FFD4E5","#FFB3BA","#FFDFBA","#FFFFBA","#BAFFC9"),
    cluster_rows = F, cluster_columns = F)
ht
name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(37)]
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(37)])
ht1 = Heatmap(df.OG2,name= "Type",col = c("#D2A8DD","#A6FDFC","#FFE3D6"),
    cluster_rows = F, cluster_columns = F)
ht1
name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(64)]
df.OG2 <- df.OG2*-1
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(64)])
ht3=  Heatmap(df.OG2,name= "Fold Change",col = colorRamp2(c(-5, 0, 5), c("#BAFFC9", "#FFFEF2","#C9C9FF")),
    cluster_rows = F, cluster_columns = FALSE)+scale_x_discrete(labels=c("s24vIri_log2_change" = "24hpf vs Iri"))
ht3
name <- make.unique(Iri_Marker_50kb_DMAR_FIMO_TF[,5])
df.OG2 <- Iri_Marker_50kb_DMAR_FIMO_TF[,c(70,72)] #Motif occurrence
row.names(df.OG2) <- name
colnames(df.OG2) <- colnames(Iri_Marker_50kb_DMAR_FIMO_TF[,c(70,72)])
ht4=  Heatmap(df.OG2,name= "RPKM",col = colorRamp2(c(0,50,200,300,400), c("#BAE1FF","#FFFFBA","#FFB3BA","#E0301E","#800000")), cluster_rows = F, cluster_columns = FALSE)
ht4
ht2+ht+ht1+ht3+ht4

######## Generate supFig8 stacked bar plot ###########
#Import table (can be found in "Input_files" folder)
data<-read.table(file="supFig8_table_Motif_presence_in_IriDMAR_nearest_IriDEGs.txt", header=TRUE,sep="\t")
#Reorder the genes, and motif combination
data$genename <- reorder(data$genename, data$Order)
data$Motif_presence <- factor(data$Motif_presence, levels=c("ALX only","SOX10 only","TFEC only", "ALX & SOX10","ALX & TFEC","SOX10 & TFEC","ALX & SOX10 & TFEC","None"))
#choose color
mypalette <-c(brewer.pal(8,"Set3"))
#plot
p3<-ggplot(data, aes(fill=Motif_presence, y=DMAR_count, x=genename))+
  geom_bar(position="stack", stat="identity")+
  labs(y= "DM/AR counts", x="Gene name")+
  coord_flip()+
  theme_bw() + 
  theme(
    panel.border = element_rect(colour = "black",size = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black",size = 0.25),
    axis.title=element_text(size=12,face = "bold"),
    axis.text.x = element_text(face = "bold",size = 12, color="black"), 
    axis.text.y = element_text(face = "bold.italic",size = 8, color="black")
)+
  theme(panel.background = element_rect(fill = 'white',colour = 'black'))+
  theme(legend.title=element_text(size=12,face = "bold"),legend.text=element_text(size=10,face = "bold"))+
  scale_fill_discrete(name = "Motif presence", breaks=c("ALX only","SOX10 only","TFEC only", "ALX & SOX10","ALX & TFEC","SOX10 & TFEC","ALX & SOX10 & TFEC","None"))+
  theme(legend.justification=c(0.95,0.95), legend.position=c(0.95,0.95))+
  scale_fill_manual(values = mypalette)+ scale_color_manual(values=mypalette)
p3
