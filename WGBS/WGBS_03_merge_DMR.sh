# Merge all DMRs into master list
## In the same working directory with outputs from R
cat DMR15v24_d30_p0.01_wSMOOTHING.txt DMR24vIri_d30_p0.01_wSMOOTHING.txt DMR24vMel_d30_p0.01_wSMOOTHING.txt DMRMelvIri_d30_p0.01_wSMOOTHING.txt> Combined_DMRs_d30_p0.01.txt
sort -k1,1 -k2,2n Combined_DMRs_d30_p0.01.txt > Combined_DMRs_d30_p0.01.sorted.txt #delete headers manually
awk -F',' '{gsub(/"/, "", $1); print $1}' Combined_DMRs_d30_p0.01.sorted.txt >  Combined_DMRs_d30_p0.01.sorted.fixed.txt #remove quotes from column 1

# merge combined DMRs
bedtools merge -i Combined_DMRs_d30_p0.01.sorted.fixed.txt > Combined_DMRs_d30_p0.01.txt 
wc -l  Combined_DMRs_d30_p0.01.txt 
rm Combined_DMRs_d30_p0.01.sorted.txt Combined_DMRs_d30_p0.01.sorted.fixed.txt

#pull out DMRs associated with combined_DMR
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR15v24_d30_p0.01_wSMOOTHING.txt > Combined_DMR15v24_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR24vMel_d30_p0.01_wSMOOTHING.txt > Combined_DMR24vMel_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMR24vIri_d30_p0.01_wSMOOTHING.txt > Combined_DMR24vIri_d30_p0.01_wSMOOTHING.txt
bedtools intersect -wao -a Combined_DMRs_d30_p0.01.txt -b  DMRMelvIri_d30_p0.01_wSMOOTHING.txt > Combined_DMRMelvIri_d30_p0.01_wSMOOTHING.txt
