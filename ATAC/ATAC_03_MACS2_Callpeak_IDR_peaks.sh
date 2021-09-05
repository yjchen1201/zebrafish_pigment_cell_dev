# Tool version
#- macs2 v2.1.1: Used to call peaks for ATAC-seq
#- python3: To use idr
#- idr/2.0.2: Use to identify reproducible peaks

# MACS2 call peak on downsampled (35million) files
macs2 callpeak -t "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.bed" -f BED -g 1.4e+9 -n "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01" -B --verbose 3 --SPMR --keep-dup all --nomodel -s 75 --extsize 73 --shift -37 -p .01 --outdir ../macs2

# Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
sort -k 8gr,8gr "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak" | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}' | gzip -nc >"SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz"

# Find IDR peaks using 0.05 as the threshold 
idr --samples "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz" "SAMPLE_Rep2_downsampled.Tshift.fixed.tagmentpositions.p01_peaks.narrowPeak.gz" --input-file-type "narrowPeak" --rank "p.value" -o "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01.IDR_peaks.txt" --log-output-file "SAMPLE_Rep1_downsampled.Tshift.fixed.tagmentpositions.p01.IDR_peaks.log" --soft-idr-threshold 0.05 --plot --verbose

# Count IDR peaks
for i in *peaks.txt; do wc -l $i;done
#172785 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#167412 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#114207 Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt
#103868 Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt


# Merge all Peaks into master peak#
cat Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks.txt
# sort based on Choord
sort -k1,1 -k2,2n All_IDR_peaks.txt > All_IDR_peaks.sorted.txt
# Merge peaks
bedtools merge -i All_IDR_peaks.sorted.txt > All_IDR_peaks_merged.txt
# count merged peak numbers
wc -l  All_IDR_peaks_merged.txt # 262676

#Identify if sample peak is called in ALL IDR merged peaks#
bedtools intersect -c -a All_IDR_peaks_merged.txt -b 15somite_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_15s.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b 24hpf_pos.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_24hpf.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b Mel.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_Mel.txt
bedtools intersect -c -a All_IDR_peaks_merged.txt -b Iri.downsampled.Tshift.p01_peaks.narrowPeak.gz.IDR_peaks.txt > All_IDR_peaks_merged_Iri.txt
