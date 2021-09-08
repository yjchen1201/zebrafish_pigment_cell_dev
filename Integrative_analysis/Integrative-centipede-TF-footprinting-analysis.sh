#### TF footprinting ####
# Requires meme. Using version 4.11.2
# Perform fimo for Alx3, Alx3_meth, Alx4_meth
fimo --max-stored-scores 10000000 --text --motif JJ0001.1 --thresh 1e-5 --verbosity 4 Alx3_meth_FIMO.meme danRer10.fa > Alx3_meth.danRer10.fimo.txt 
fimo --max-stored-scores 10000000 --text --motif JJ0002.1 --thresh 1e-5 --verbosity 4 Alx4_meth_FIMO.meme danRer10.fa > Alx4_meth.danRer10.fimo.txt 
fimo --max-stored-scores 10000000 --text --motif PH0001.1 --thresh 1e-5 --verbosity 4 Alx3_FIMO.meme danRer10.fa > Alx3.danRer10.fimo.txt 
## First column of the .fimo.txt file is the pattern name. Remove for downstream analysis. 
awk -v OFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$9}' Alx3.danRer10.fimo.txt > Alx3.danRer10.fimo2.txt #make sure to remove the first column
## Repeat the same process for Alx3_meth.danRer10.fimo.txt and Alx4_meth.danRer10.fimo.txt

# Transform ATAC-seq data into tagment end format 
## Files are generated in ATAC/ATAC_02_generate_input_files_for_peakcalling_and_centipede.r
## Change directory to where tagmentpositions.bed are located, correct annotations and sort
for i in *combined*bed; do awk '{$3=sprintf("%.0f",$3)}7' $i | '{$2=sprintf("%.0f",$2)}7' | sort -k1,1 -k2,2n > ${i/.bed/.sorted.bed};done 
## Convert sorted bed into bigwig format
for i in *combined*sorted.bed; do bedGraphToBigWig $i danRer10.chrom.sizes ${i/sorted.bed/bw};done

# Preparation for centipede
## 1 make motif site bed file
awk 'BEGIN{OFS="\t"} NR>1 {print $2,$3-1,$4,NR-1,$6,$5}' Alx3.danRer10.fimo.txt > Alx3_sites_danRer10.bed

## 2 make phastCons file (Y)
bed_bedGraph.pl Alx3_sites_danRer10.bed danRer10.vertebrate.phastCons8way.bg.gz > tmp.bg
paste Alx3_sites_danRer10.bed tmp.bg | cut -f1-3,5-6,10 - | gzip > Alx3_sites_centipedeY_phastCons8way.txt.gz #chr start end lor strand phastCons
rm tmp.bg

## 3 and #4 make matrix ATAC signal around motif sites (X)
### Using Alx3 as an example
grep "+" Alx3_sites_danRer10.bed > Alx3_sites_danRer10.for.bed
grep "-" Alx3_sites_danRer10.bed > Alx3_sites_danRer10.rev.bed

## Determine the length of motif
width=$( head -1 Alx3_sites_danRer10.bed | awk '{print $3-$2}' )
gzip Alx3_sites_danRer10.bed
left=$(( $width/2 + 101 ))
right=$(( $width/2 + 100 - ($width+1)%2 ))

# Focus on Iridophore
## Motifs found on forward strand
bwtool matrix ${left}:$right <(cut -f1-3 Alx3_sites_danRer10.for.bed) Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw,Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw /dev/stdout |
	sed "s/NA/0/g; s/\.00\t/\t/g; s/\.00$//g" |
	paste Alx3_sites_danRer10.for.bed - > Alx3_sites_centipedeX_ATAC_Iri.for.txt
## Motifs found on reverse strand: reverse the order of data
bwtool matrix ${left}:$right <(cut -f1-3 Alx3_sites_danRer10.rev.bed) Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw,Iri_combined_downsampled.Tshift.fixed.tagmentpositions.bw /dev/stdout |
	sed "s/NA/0/g; s/\.00\t/\t/g; s/\.00$//g" > Alx3_sites_centipedeX_ATAC_Iri.rev.txt
awk 'BEGIN{ORS=""} {for(i=NF;i>0;i--) {print $i; if (i==1) print "\n"; else print "\t"}}' Alx3_sites_centipedeX_ATAC_Iri.rev.txt |
	paste Alx3_sites_danRer10.rev.bed - |
	cat Alx3_sites_centipedeX_ATAC_Iri.for.txt - |
	sort -k1,1 -k2,2n | cut -f 7- | gzip > Alx3_sites_centipedeX_ATAC_Iri.txt.gz
rm Alx3_sites_centipedeX_ATAC_Iri.for.txt Alx3_sites_centipedeX_ATAC_Iri.rev.txt
## Plot footprints
./fitCentipede2.R Alx3_sites_centipedeX_ATAC_Iri.txt.gz Alx3_sites_centipedeY_phastCons8way.txt.gz Alx3_sites_ATAC_Iri
