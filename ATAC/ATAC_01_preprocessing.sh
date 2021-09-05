# Tool version
#- cutadapt v1.10: Used to trim adapter reads 
#- samtools version 1.3.1: Used to sort and downsample bam files for downstream processing
#- picard v2.8.1: Used to mark and remove duplicate reads
#- methylQA v0.1.6: Use to covert bam to single-end read bed file

# Adapter trimming
## Input SAMPLE_R1.fastq.gz, SAMPLE_R2.fastq.gz. 
## Output SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 25 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1_fastq.gz" "SAMPLE_R2_fastq.gz"; done

# Reads alignment
## Input SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
## Output SAMPLE.sam
bwa mem "danRer10.fa" "SAMPLE_R1_trimmed.fastq.gz" "SAMPLE_R2_trimmed.fastq.gz" > "SAMPLE.sam"

#SAM->BAM, and only include reads with MAPQ>=30
## Input SAMPLE.sam
- samtools version 1.3.1: Used to sort and downsample bam files for downstream processing
## Output SAMPLE_MAPQ30.bam
samtools view -b -q 30 -o "SAMPLE_MAPQ30.bam" "SAMPLE.sam" 


# sort bam
## Input SAMPLE_MAPQ30.bam
## Output SAMPLE_MAPQ30.sorted.bam
samtools sort -o "SAMPLE_MAPQ30.sorted.bam" "SAMPLE_MAPQ30.bam"


# Remove duplicates
## Input SAMPLE_MAPQ30.sorted.bam
## Output SAMPLE_MAPQ30.sorted.rmdup.bam
java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=True REMOVE_DUPLICATES=True "SAMPLE_MAPQ30.sorted.bam" "SAMPLE_MAPQ30.sorted.rmdup" "SAMPLE_MAPQ30.sorteddup_removed.txt" ASSUME_SORTED=true


#Downsample reads to 35million (Example using samtool view -s subsample command)
## Input SAMPLE_MAPQ30.sorted.rmdup.bam
## Output SAMPLE.downsampled.bam
samtools view -h -b -s 0.3125 "SAMPLE_MAPQ30.sorted.rmdup.bam" > "SAMPLE_MAPQ30.sorted.rmdup.downsampled.bam"

#Calculate reads fraction in peaks (RFIP)#
## Input SAMPLE_MAPQ30.sorted.rmdup.downsampled.bam
for i in *downsampled.bam; do samtools view $i | wc -l ;done
#35020476 15somite_pos_Rep1.downsampled.bam
#34969034 15somite_pos_Rep2.downsampled.bam
#35007919 Iri_Rep1.downsampled.bam
#35016630 Iri_Rep2.downsampled.bam
#34991229 Mel_Rep1.downsampled.bam
#35016721 Mel_Rep4.downsampled.bam
#35040232 24hpf_pos_Rep1.downsampled.bam
#34975335 24hpf_pos_Rep2.downsampled.bam

# Using MethylQA to covert bam to single-end read bed file
## Input SAMPLE_MAPQ30.sorted.rmdup.downsampled.bam
## Output SAMPLE_MAPQ30.sorted.rmdup.downsampled.extend.bed
methylQA density -Q 10 -r -T -E 0 -I 2000 -o "SAMPLE_MAPQ30.sorted.rmdup.downsampled" "danRer10_lite.size" "SAMPLE_MAPQ30.sorted.rmdup.downsampled.bam"


# Shift transposome inserted site
## Input SAMPLE_MAPQ30.sorted.rmdup.downsampled.extend.bed
## Output SAMPLE_MAPQ30.sorted.rmdup.downsampled.Tshift.bed
awk -v OFS="\t" '{if ($6=="+") print $1,$2+4,$3+4,$4,$5,$6; else if ($6=="-") print $1,$2-5,$3-5,$4,$5,$6}' "SAMPLE_MAPQ30.sorted.rmdup.downsampled.extend.bed" > "SAMPLE_MAPQ30.sorted.rmdup.downsampled.Tshift.bed"

#NOTE THAT THIS CAN LEAD TO NEGATIVE VALUES FOR READS AT ENDS OF CHROMOSOMES, so make sure the start site is >0
```
