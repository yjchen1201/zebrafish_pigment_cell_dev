# Tools:
#- Bedtools: Use intersect function
#- Deeptools: Use to calculate ATAC/WGBS signal matrix and plot signal heatmap


# Extract only chr, start, end, chrompo columns from Iri DM/AR files and sort 
## Input files listed below are generated using script in "Integrative-DMAR-analysis.r":
### Iri_solo_hypoDMR_s24vsIri.bed
### Iri_solo_hyperDMR_s24vsIri.bed
### Iri_solo_openDAR_s24vsIri.bed
### Iri_solo_closeDAR_s24vsIri.bed
### Iri_hypoclosing_DMAR_s24vsIri.bed
### Iri_hypoopening_DMAR_s24vsIri.bed

for i in *bed; do awk -v OFS="\t" '{print $1,$2,$3,$4}' $i | sort -k1,1 -k2,2n > ${i/.bed/.onlyCoord.srt.bed};done 

# Output files:
#Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.srt.bed
#Iri_hypoopening_DMAR_s24vsIri.onlyCoord.srt.bed
#Iri_solo_hyperDMR_s24vsIri.onlyCoord.srt.bed
#Iri_solo_hypoDMR_s24vsIri.onlyCoord.srt.bed
#Iri_solo_openDAR_s24vsIri.onlyCoord.srt.bed
#Iri_solo_closeDAR_s24vsIri.onlyCoord.srt.bed

# Find out motif presence in Iri-associated DM/ARs using bedtools intersect
## MOTIF DBD location files are generated using code ..., you can also find these files in "Input_files" folder
#ALX4_DBD_sites_danRer10.bed
#GBX2_DBD_sites_danRer10.bed
#SOX10_DBD_sites_danRer10.bed
#TFEC_DBD_sites_danRer10.bed
for i in *onlyCoord.srt.bed; do bedtools intersect -wa -a $i -b ALX4_DBD_sites_danRer10.bed > ${file/.srt.bed/.with_ALX4_MOTIF.bed};done &
for i in *onlyCoord.srt.bed; do bedtools intersect -wa -a $i -b GBX2_DBD_sites_danRer10.bed > ${file/.srt.bed/.with_GBX2_MOTIF.bed};done &
for i in *onlyCoord.srt.bed; do bedtools intersect -wa -a $i -b TFEC_DBD_sites_danRer10.bed > ${file/.srt.bed/.with_TFEC_MOTIF.bed};done &
for i in *onlyCoord.srt.bed; do bedtools intersect -wa -a $i -b SOX10_DBD_sites_danRer10.bed > ${file/.srt.bed/.with_SOX10_MOTIF.bed};done &

# Use deeptools to generate heatmap plot for DNA methylation and ATAC signal across iridophore-associated DM/ARs with a particular TF motif # [Fig.3d and SupFig.7a]
# The ATAC signal bigWig files are generated using code in "ATAC/ATAC_01_preprocessing.sh"
# The WGBS singal bigWig files are generated using code in "WGBS/WGBS_02_DML-DMR.r"
## 24NCC and Iri epi signals at those regions ##  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Iri.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Iri_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Iri_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Iri_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Iri.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Iri_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Iri_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Iri_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Iri.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Iri_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Iri_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Iri_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Iri.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Iri_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Iri_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Iri_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
# Plot heatmap
plotHeatmap -m Iri_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Iri_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "s24vIri soloDAR soloDMAR DMAR with ALX4 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
plotHeatmap -m Iri_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Iri_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "s24vIri soloDAR soloDMAR DMAR with GBX2 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
plotHeatmap -m Iri_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Iri_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "s24vIri soloDAR soloDMAR DMAR with TFEC motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
plotHeatmap -m Iri_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Iri_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "s24vIri soloDAR soloDMAR DMAR with SOX10 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
  
  
## 24NCC and Mel epi signals at those regions ##  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Mel.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Mel_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_ALX4_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Mel_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Mel_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Mel.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Mel_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_GBX2_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Mel_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Mel_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Mel.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Mel_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_TFEC_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Mel_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Mel_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &
  
computeMatrix scale-regions -bs 100 -p max \
-S /scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_s24.bw \
/scratch/jjang/PIGMENT_PROJECT/Rebuttal/DMR_BioRep/SmoothedCpG_Methylation_Mel.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_s24pos_merged_forDeeptools.bw \
/scratch/jjang/PIGMENT_PROJECT/Pigment_ATAC_ALL/ATAC_Mel_merged_forDeeptools.bw \
-R Iri_solo_openDAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_solo_closeDAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_solo_hypoDMR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_hypoopening_DMAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
Iri_hypoclosing_DMAR_s24vsIri.onlyCoord.with_SOX10_MOTIF.bed \
--beforeRegionStartLength 1500 --regionBodyLength 1000 --afterRegionStartLength 1500 \
-o Mel_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz > Mel_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz.log 2>&1 &

# Plot heatmap
plotHeatmap -m Mel_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Mel_ALX4_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "Mel in s24vIri soloDAR soloDMAR DMAR with ALX4 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30 
plotHeatmap -m Mel_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Mel_GBX2_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "Mel in s24vIri soloDAR soloDMAR DMAR with GBX2 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
plotHeatmap -m Mel_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Mel_TFEC_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "Mel in s24vIri soloDAR soloDMAR DMAR with TFEC motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
plotHeatmap -m Mel_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.gz -out Mel_SOX10_in_soloDAR_soloDMR_DMAR_s24vIri_bs100_1500UpAndDown.mat.pdf -T "Mel in s24vIri soloDAR soloDMAR DMAR with SOX10 motif" --colorList 'blue,yellow,red' 'blue,yellow,red' 'white,green' 'white,green' --missingDataColor "white" --zMin 0 0 0 0 --zMax 1 1 120 120 --yMin 0 0 0 0 --yMax 1 1 110 110  --heatmapWidth 12 --refPointLabel "Center" --heatmapHeight 30
