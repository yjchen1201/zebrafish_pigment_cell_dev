# Tool version
#- bedtools v2.27.1: closest, intersect, shuffle command used as specified in manuscript 
#- meme v5.0.3: Used AME command to find motif enrichment and FIMO command to scan for motif presence 

######### Use AME to do motif enrichment analysis #########

#make background file of all peaks 
## Input file is generated from ATAC/ATAC_04_Merge_IDR_peaks_identify_specific_shared_IDR_peaks.r
# Extract chr, start, end, and chrompos information
awk -v OFS="\t" '{print $1, $2, $3, $4}' All_IDR_peaks_Merged.txt > All_IDR_peaks_Merged_forfasta.txt
# Convert to fasta using KentUCSC tool twoBitToFa
## the danRer10.2bit was generated from danRer10.fa using UCSC faToTwoBit https://genome.ucsc.edu/goldenPath/help/twoBit.html
twoBitToFa -bed=All_IDR_peaks_Merged_forfasta.txt danRer10.2bit All_IDR_peaks_Merged.fa


# Make fasta files of DARs
## Input *hypoDMR_d30_p0.01.bed files are generated from /WGBS/WGBS_04_plot_merge_DMR.r
## Input *opening_DAR.bedATAC files are generated from /ATAC_05_DiffBind_call_DARs_from_IDR_peaks.r

#hypoDMRs#
awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_hypoDMR_d30_p0.01.bed > Iri_specific_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Iri_specific_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Iri_specific_hypoDMR_d30_p0.01.fasta

ame --verbose 1 --oc Iri_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_hypoDMR_d30_p0.01.bed > Mel_specific_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_specific_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Mel_specific_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_hypoDMR_d30_p0.01.bed > Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Mel_Iri_shared_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_Iri_shared_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_hypoDMR_d30_p0.01.bed > Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypoDMR_d30_p0.01.onlyCoord.bed danRer10.2bit Mel_Iri_shared_hypoDMR_d30_p0.01.fasta
ame --verbose 1 --oc Mel_Iri_shared_hypoDMR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_hypoDMR_d30_p0.01.fasta JASPAR_2016_Vertebrate.meme

#DARs#
awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_opening_DAR.bed > Iri_specific_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_opening_DAR.onlyCoord.bed danRer10.2bit Iri_specific_opening_DAR.fasta
ame --verbose 1 --oc Iri_specific_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_opening_DAR.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Iri_specific_closing_DAR.bed > Iri_specific_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_closing_DAR.onlyCoord.bed danRer10.2bit Iri_specific_closing_DAR.fasta
ame --verbose 1 --oc Iri_specific_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_closing_DAR.fasta JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_opening_DAR.bed > Mel_specific_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_opening_DAR.onlyCoord.bed danRer10.2bit Mel_specific_opening_DAR.fasta
ame --verbose 1 --oc Mel_specific_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_opening_DAR.fasta JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_specific_closing_DAR.bed > Mel_specific_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_closing_DAR.onlyCoord.bed danRer10.2bit Mel_specific_closing_DAR.fasta
ame --verbose 1 --oc Mel_specific_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_closing_DAR.fasta JASPAR_2016_Vertebrate.meme


awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_opening_DAR.bed > Mel_Iri_shared_opening_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_opening_DAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_opening_DAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_openDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_opening_DAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $1,$2,$3,$1":"$2"-"$3}' Mel_Iri_shared_closing_DAR.bed > Mel_Iri_shared_closing_DAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_closing_DAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_closing_DAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_closeDAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_closing_DAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme


#DMARs#
awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Iri_specific_hypo_opening.DMAR.bed > Iri_specific_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Iri_specific_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Iri_specific_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Iri_specific_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Iri_specific_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Mel_specific_hypo_opening.DMAR.bed > Mel_specific_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Mel_specific_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Mel_specific_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Mel_specific_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_specific_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme

awk 'OFS="\t" {print $4,$1,$2,$4":"$1"-"$2}' Mel_Iri_shared_hypo_opening.DMAR.bed > Mel_Iri_shared_hypo_opening.DMAR.onlyCoord.bed
twoBitToFa -bed=Mel_Iri_shared_hypo_opening.DMAR.onlyCoord.bed danRer10.2bit Mel_Iri_shared_hypo_opening.DMAR.fasta
ame --verbose 1 --oc Mel_Iri_shared_openDMAR -scoring totalhits --method fisher --control All_IDR_peaks_Merged.fa Mel_Iri_shared_hypo_opening.DMAR.fasta /scratch/jjang/PIGMENT_PROJECT/Pigment_integrative_analysis_051418/JASPAR_2016_Vertebrate.meme
