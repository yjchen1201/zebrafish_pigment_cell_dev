# Combine DMR and DAR regions together to one master list
## Input file "Combined_DMRs_d30_p0.01_wINFO.bed" is generated from WGBS/WGBS_04_plot_merge_DMR.r
## Input file "All_DiffBind_peaks_merged_wINFO.bed" is generated from ATAC/ATAC_05_DiffBind_call_DARs_from_IDR_peaks.r
cat Combined_DMRs_d30_p0.01_wINFO.bed All_DiffBind_peaks_merged_wINFO.bed >Combined_DMR_DAR.bed
# Sort regions based on choord
sort -k1,1 -k2,2n Combined_DMR_DAR.bed > Combined_DMR_DAR.sorted.bed
# Merge overlapped regions
bedtools merge -i Combined_DMR_DAR.sorted.bed > Combined_DMR_DAR.bed

# Add DAR and DMR information to merged region list
## Input file "Combined_DMRs_d30_p0.01_wINFO.bed" is generated from WGBS/WGBS_04_plot_merge_DMR.r
## Input file "All_DiffBind_peaks_merged_wINFO.bed" is generated from ATAC/ATAC_05_DiffBind_call_DARs_from_IDR_peaks.r
bedtools intersect -wao -a Combined_DMR_DAR.bed -b Combined_DMRs_d30_p0.01_wINFO.bed > Combined_DMR_DAR_wDMRinfo.bed
bedtools intersect -wao -a Combined_DMR_DAR.bed -b All_DiffBind_peaks_merged_wINFO.bed > Combined_DMR_DAR_wDARinfo.bed

# Add IDR peak information
## Input *downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt files are generated from ATAC/ATAC_03_MACS2_Callpeak_IDR_peaks.sh
bedtools intersect -c -a Combined_DMR_DAR.bed -b 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_15s.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_24hpf.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_Mel.txt
bedtools intersect -c -a Combined_DMR_DAR.bed -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > Combined_DMR_DAR_wIDRinfo_Iri.txt

# Calculate average methylation of DMARs
## Input files generated from WGBS/WGBS_02_DML-DMR.r: 
### 15somite_NCC_Combined_DSS.txt 
### 24hpf_NCC_Combined_DSS.txt 
### Mel_NCC_Combined_DSS.txt 
### Iri_NCC_Combined_DSS.txt 
for i in *Combined_DSS.txt; do sed '1d' $i | awk 'BEGIN{OFS="\t"}{print $1,$2,$2,$3,$4,1}' > ${i/.txt/.bed};done
bedtools intersect -wo -a Combined_DMR_DAR.bed -b 15somite_NCC_Combined_DSS.bed > Combined_DMR_DAR_wMETH_15s.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b 24hpf_NCC_Combined_DSS.bed > Combined_DMR_DAR_wMETH_24hpf.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b Mel_Combined_DSS.bed > Combined_DMR_DAR_wMETH_Mel.txt
bedtools intersect -wo -a Combined_DMR_DAR.bed -b Iri_Combined_DSS.bed > Combined_DMR_DAR_wMETH_Iri.txt
