# Use https://biit.cs.ut.ee/gprofiler/orth for gene conversion for AME ##
#Some gene names still missing. Manually curate using http://www.zebrafishmine.org/template.do?name=Human_Gene_Curated_Zeb_Ortho. Use https://resources.altius.org/~jvierstra/projects/motif-clustering/ to identify cluster and DBD

#Load AME findings on R
masterTF_orth_list <- read.csv("MASTER-LIST-TF-orthologs_converted_ensemble_withMotifClusterInfo.csv", header = T, sep = ",", stringsAsFactors = F)
colnames(masterTF_orth_list) <- c("motif_alt_ID","gene","genename","Cluster","DBD")

# Import the output from ame motif finding
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

#Import DEG information 
DEG2_Z <- read.table("DEG_TPM_Zscore_p0.01_TPMfilter5.txt",header = T, sep = "\t",stringsAsFactors =F) 

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
Iri_hypoDMR_AME_filtered <-Iri_hypoDMR_AME_genename_DEG[Iri_hypoDMR_AME_genename_DEG$s24vIri_log2_change <0 | Iri_hypoDMR_AME_genename_DEG$MelvIri_log2_change <0 & Iri_hypoDMR_AME_genename_DEG$Iri > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri")] #upreg in iri **Use this
Iri_openDAR_AME_filtered <-Iri_openDAR_AME_genename_DEG[Iri_openDAR_AME_genename_DEG$s24vIri_log2_change <0 | Iri_openDAR_AME_genename_DEG$MelvIri_log2_change <0 & Iri_openDAR_AME_genename_DEG$Iri > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"),] #upreg in iri 
Iri_openDMAR_AME_filtered <-Iri_openDMAR_AME_genename_DEG[Iri_openDMAR_AME_genename_DEG$s24vIri_log2_change <0 | Iri_openDMAR_AME_genename_DEG$MelvIri_log2_change <0 & Iri_openDMAR_AME_genename_DEG$Iri > 5,c("genename.x","pvalue_log","DBD","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri"),] #upreg in iri 

Iri_AME_masterlist <- unique(rbind(Iri_hypoDMR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Iri_openDAR_AME_filtered[,c(1,3,4,5,6,7,8,9)],Iri_openDMAR_AME_filtered[,c(1,3,4,5,6,7,8,9)]))

Iri_AME_masterlist_withLogP <- Reduce(function(x, y) merge(x, y, by=c("genename.x","s24vMel_log2_change","s24vIri_log2_change","MelvIri_log2_change","s24","Mel","Iri","DBD"), all = T), list(Iri_AME_masterlist,Iri_hypoDMR_AME_filtered,Iri_openDAR_AME_filtered,Iri_openDMAR_AME_filtered))
Iri_AME_masterlist_withLogP[is.na(Iri_AME_masterlist_withLogP)] <- 0

Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),]

Iri_AME_masterlist_withLogP_all3 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x >0 & Iri_AME_masterlist_withLogP$pvalue_log.y >0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]
Iri_AME_masterlist_withLogP_all2 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x ==0 & Iri_AME_masterlist_withLogP$pvalue_log.y >0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]
Iri_AME_masterlist_withLogP_all1 <- Iri_AME_masterlist_withLogP[Iri_AME_masterlist_withLogP$pvalue_log.x ==0 & Iri_AME_masterlist_withLogP$pvalue_log.y ==0 & Iri_AME_masterlist_withLogP$pvalue_log >0,]

Iri_AME_masterlist_withLogP_reordered <- rbind(Iri_AME_masterlist_withLogP_all3,Iri_AME_masterlist_withLogP_all2,Iri_AME_masterlist_withLogP_all1)
Iri_AME_masterlist_withLogP_reordered<-unique(Iri_AME_masterlist_withLogP_reordered)

### Generate Fig3b heatmap plot ###
library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options  
# Re-order motif based on DBD and plot #
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
ht = Heatmap(df.OG, column_title = "TFs",name= "TPM",col = colorRamp2(c(0,20,40,60,80), c("#bae1ff","#ffffba","#ffb3ba","#e0301e","#800000")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht

df.OG <- Iri_AME_masterlist_withLogP[order(Iri_AME_masterlist_withLogP[,"DBD"],Iri_AME_masterlist_withLogP[,"genename.x"]),c(9,10,11)]
row.names(df.OG) <- name
#kclus <- kmeans(df.OG,10)
#split <-  kclus$cluster
ht3 = Heatmap(df.OG, column_title = "TFs",name= "-log pval",col = colorRamp2(c(0,10,200,300,400), c("#ffffff","#f1ebe1","#c0cfb2","#8ba888","#44624a")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = F)
ht3 #12 x 3

### Generate Fig3a heatmap plot ###
# Mel open AME motif figure #
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

# Re-order motif based on DBD and plot #
name <- Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),1]
df.OG2 <- -Mel_AME_masterlist_withLogP[order(Mel_AME_masterlist_withLogP[,"DBD"],Mel_AME_masterlist_withLogP[,"genename.x"]),c(2,3,4)]
row.names(df.OG2) <- name
colnames(df.OG2) <- c("24hpf vs Mel", "24hpf vs Iri", "Mel vs Iri")
htname =  Heatmap(df.OG2, column_title = "DEG Transcription Factors",name= "log2 Fold Change",col = colorRamp2(c(-5, 0, 5), c("#baffc9", "#fffef2","#c9c9ff")), 
    cluster_rows = F, cluster_columns = FALSE,show_row_names = T)
htname

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
     
                                      
                                      
### Generate Fig 3c bar plots representing frequency and distribution of iridophore-associated DM/ARs with a particular TF motif ###
# ALX1_ALX3_ALX4_GBX2_TFEC_SOX10_MOTIF_occurrence_in_Iri_DMR_DAR_DMAR plot
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

# DMR_DAR_DMAR_occurence plot
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
