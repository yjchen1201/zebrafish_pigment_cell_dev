
# Merge the IDR peaks in R
## Input files are from ATAC_03_MACS2_Callpeak_IDR_peaks.sh output
Mel<- read.table("All_IDR_peaks_merged_Mel.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
Iri<- read.table("All_IDR_peaks_merged_Iri.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
s24<- read.table("All_IDR_peaks_merged_24hpf.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")
s15<- read.table("All_IDR_peaks_merged_15s.txt",, header = F, sep = "\t", stringsAsFactors = F, quote = "")

# Add "peakname" column, eg: chr1:234-567
Mel$peakname <- paste(Mel$V1,":",Mel$V2,"-",Mel$V3,sep = "")
Iri$peakname <- paste(Iri$V1,":",Iri$V2,"-",Iri$V3,sep = "")
s24$peakname <- paste(s24$V1,":",s24$V2,"-",s24$V3,sep = "")
s15$peakname <- paste(s15$V1,":",s15$V2,"-",s15$V3,sep = "")

# Get IDR peak size
s15$size <- s15$V3-s15$V2

# Merge IDR peak information based on peakname
t <- merge(s15[,c("peakname","V1","V2","V3","size","V4")], s24[,c("peakname","V4")], by = "peakname", all = T)
t <- merge(t, Mel[,c("peakname","V4")], by = "peakname", all = T)
t <- merge(t, Iri[,c("peakname","V4")], by = "peakname", all = T)

# Save the merged peak info table
write.table(t[c(2,3,4,1,5,6,7,8,9)], "All_IDR_peaks_Merged.txt", sep = "\t", col.names = F, row.names = F, quote = F)
# Add column names
colnames(t) <- c("peakname","chr","start","end","size","s15","s24","Mel","Iri")

# Extract IDR peaks based on different criteria
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

# Save Cell type-specific and shared IDR peaks table
write.table(Mel_specific[,c(2,3,4,1)],"Mel_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Iri_specific[,c(2,3,4,1)],"Iri_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel_Iri_shared_dynamce[,c(2,3,4,1)],"Mel_Iri_shared_dynamic_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(rbind(s15_Mel_Iri_peaks[,c(2,3,4,1)],Mel_Iri_peaks[,c(2,3,4,1)]),"Mel_Iri_shared_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(rbind(s15_24s_peaks[,c(2,3,4,1)],only24s_peaks[,c(2,3,4,1)]),"24hpf_specific_IDRpeaks.bed", sep = "\t", col.names = F, row.names = F, quote = F)

#Plot heatmap
name <- allpeaks[,1]
df.OG <- allpeaks[,c(6:9)]
row.names(df.OG) <- name

library("ComplexHeatmap") ## For heatmap
library("circlize") ## For color options
ht = Heatmap(df.OG, col = colorRamp2(c(0, 1, 2,3), c("#f6f6f6", "#74d600","#7bc043","#028900")), 
    cluster_rows = FALSE, cluster_columns = FALSE,show_row_names = FALSE)
ht
