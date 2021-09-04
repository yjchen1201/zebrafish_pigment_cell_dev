# Adapter trimming
## Input SAMPLE_R1.fastq.gz, SAMPLE_R2.fastq.gz. 
## Output SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1.fastq.gz" "SAMPLE_R2.fastq.gz"

# Reads alignment
## Input SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
## Output SAMPLE_R1_trimmed_bismark_bt2_pe.bam
bismark --bowtie2 --skip 9 -N 1 -L 28 --score_min L,0,-0.6 "reference genome path" -1 "SAMPLE_R1_trimmed.fastq.gz" -2 "SAMPLE_R2_trimmed.fastq.gz"

# SAM->BAM, and only include reads mapped to the input bed region 
## Input SAMPLE_R1_trimmed_bismark_bt2_pe.bam
## Output SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam
samtools view -b -L "danRer10_genome_lite.size.bed" "SAMPLE_R1_trimmed_bismark_bt2_pe.bam" > "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam"

# Sort 
## Input SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam
## Output SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam
samtools sort -o "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam" "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.bam"  

# Remove duplicate
## Input SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam
## Output SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam
java -jar $PICARD MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=True REMOVE_DUPLICATES=True I="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.bam" O="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam" M="SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.dup_removed.txt" ASSUME_SORTED=true

# Sort bam file by read name
## Input SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam
## Output SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam
samtools sort -n -o "SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam" "SAMPLE_R1_trimmed_bismark_bt2_pe_chromlite.sorted.rmdup.bam"

# Extract the methylation information from the Bismark alignment output
## Input SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam
## Output ?
bismark_methylation_extractor -p --no_overlap --comprehensive --gzip --no_header --bedgraph --counts --ignore 10 --ignore_r2 10 --report -o ./ --genome_folder <path to genome folder> "SAMPLE_R1_trimmed_bismark_chromlite.resorted.rmdup.bam"

##Bismark ON LAMBDA DNA###
bismark --quiet <path to lambda genome> -1 "SAMPLE_R1.fastq.gz" -2 "SAMPLE_R2.fastq.gz"


###convert bismark coverage output into DSS input file format
## Input ?
## Output ?
python3 convert_DSS.py Merged_ALLRUNS_Iri_Rep1_trimmed_bismark_resorted.rmdup.bismark.cov