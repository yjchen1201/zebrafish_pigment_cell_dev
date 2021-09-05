#R/3.4.1
#library(dplyr)
# Make genome browser files on just the ends of ATAC reads 

# Import shifted bed file (Iridophore and Melanophore samples)
## Input files: from ATAC_01_preprocessing: SAMPLE_MAPQ30.sorted.rmdup.downsampled.Tshift.bed
Iri1 <- read.table("Iri_Rep1_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
Iri2 <- read.table("Iri_Rep2_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
Mel1 <- read.table("Mel_Rep1_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
Mel2 <- read.table("Mel_Rep2_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)

# Fix the start site to 1 based
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

# Create shifted bed files
write.table(Iri1_aggregate, "Iri_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Iri2_aggregate, "Iri_Rep2_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel1_aggregate, "Mel_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(Mel2_aggregate, "Mel_Rep4_downsampled.Tshift.fixed.tagmentpositions.bed", sep = "\t", col.names = F, row.names = F, quote = F)



# Import shifted bed file (24hpf and 15 somite samples)
## Input files: from ATAC_01_preprocessing: SAMPLE_MAPQ30.sorted.rmdup.downsampled.Tshift.bed
s24_1 <- read.table("24hpf_pos_Rep1_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
s24_2 <- read.table("24hpf_pos_Rep2_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
s15_1 <- read.table("15somite_pos_Rep1_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
s15_2 <- read.table("15somite_pos_Rep2_MAPQ30.sorted.rmdup.downsampled.Tshift.bed", sep = "\t", header = F, stringsAsFactors = F)
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



#Combine replicates together for Centipede
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