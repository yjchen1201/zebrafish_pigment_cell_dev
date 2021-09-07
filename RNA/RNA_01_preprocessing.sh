# RNA-seq data preprocessing
# Adapter trimming
## Input SAMPLE_R1.fastq.gz, SAMPLE_R2.fastq.gz
## Output SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --minimum-length 50 -o "SAMPLE_R1_trimmed.fastq.gz" -p "SAMPLE_R2_trimmed.fastq.gz" "SAMPLE_R1.fastq.gz" "SAMPLE_R2.fastq.gz"

# Alignment using STAR, output is BAM file sorted by coordinate
## Input SAMPLE_R1_trimmed.fastq.gz, SAMPLE_R2_trimmed.fastq.gz
## Output SAMPLE_sortedByCoord.bam
## We used danRer10 as genome reference, and Danio_rerio.GRCz10.85.gtf as annotation
STAR --readFilesCommand zcat --outTmpDir <out_temp_dirctory> --clip5pNbases 11 --outSAMtype BAM SortedByCoordinate --genomeDir <genome_directory> --sjdbGTFfile <annotation_gtf_path> --outFileNamePrefix <output_prefix> --readFilesIn "SAMPLE_R1_trimmed.fastq.gz" "SAMPLE_R2_trimmed.fastq.gz"

# Normalziation to RPM
## Input SAMPLE_sortedByCoord.bam
## Output SAMPLE_sortedByCoord.out.bam
STAR --runMode inputAlignmentsFromBAM --inputBAMfile "SAMPLE_sortedByCoord.bam" --outWigType bedGraph --outWigStrand Unstranded --outWigNorm RPM --outFileNamePrefix "SAMPLE_sortedByCoord.out.bam"

# Count alignment by gene
## Input SAMPLE_sortedByCoord.out.bam
## Output SAMPLE_sortedByCoord.counts
## Danio_rerio.GRCz10.85.gtf as annotation
htseq-count -f bam -r pos -s no -t transcript "SAMPLE_sortedByCoord.out.bam" <annotation_gtf_path> > "SAMPLE_sortedByCoord.counts"
# Count alignment by transcript
stringtie "SAMPLE_sortedByCoord.out.bam" -o "RefSeq.gtf" -A "RefSeq.gene.abundance.txt" -G "danRer10_RefSeq_forstringtie.gtf"

#Filter NM genes only
python3 Filter_only_NM.py "SAMPLE_abundance.txt"
#Make list of all genes expressed in NCC and pigment cells (RefSeq)
## Output Combined_gene_list_RefSeq.txt
python3 All_expressed_genelist_RefSeq.py RNA_15somiteNCC_Rep1.gene.abundance.filteredNM.txt RNA_15somiteNCC_Rep2.gene.abundance.filteredNM.txt RNA_24hpfNCC_Rep1.gene.abundance.filteredNM.txt RNA_24hpfNCC_Rep2.gene.abundance.filteredNM.txt RNA_Mel_Rep1.gene.abundance.filteredNM.txt RNA_Mel_Rep2.gene.abundance.filteredNM.txt RNA_Iri_Rep1.gene.abundance.filteredNM.txt RNA_Iri_Rep2.gene.abundance.filteredNM.txt 
