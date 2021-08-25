data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_15somite_pos_Rep2_RNA_AR008_ACTTGAATC_S5_R1_trimmed..gzAligned.sortedByCoord.out.sorted.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_15somite_pos_Rep2_RNA_AR008_ACTTGAATC_S5_R1_trimmed..gzAligned.sortedByCoord.out.sorted.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_15s_pos_Rep4_RNA_AR022_CGTACGTAA_S6_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_15s_pos_Rep4_RNA_AR022_CGTACGTAA_S6_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_24GFP_pos_Rep1_AR001_ATCACG_S1_R1_trimmed_001.Aligned.sortedByCoord.out.sorted.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_24GFP_pos_Rep1_AR001_ATCACG_S1_R1_trimmed_001.Aligned.sortedByCoord.out.sorted.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_24GFP_pos_Rep2_AR008_ACTTGA_S3_R1_trimmed_001.Aligned.sortedByCoord.out.sorted.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_24GFP_pos_Rep2_AR008_ACTTGA_S3_R1_trimmed_001.Aligned.sortedByCoord.out.sorted.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Mel_R1_RNA_AR009_GATCAGATC_S2_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Mel_R1_RNA_AR009_GATCAGATC_S2_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Mel_R2_RNA_AR010_TAGCTTATC_S3_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Mel_R2_RNA_AR010_TAGCTTATC_S3_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_RNA_IRI_rep3_R1.trimmed..gzAligned.sortedByCoord.out.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_RNA_IRI_rep3_R1.trimmed..gzAligned.sortedByCoord.out.gene.abundance.txt",'w')
out.write(outdata)

data = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Iri_R2_RNA_AR021_GTTTCGGAA_S5_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt", 'r')
outdata = ""
for line in data:
	if line[0] == "G":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[0][0:2] == "EN":
			 outdata += line[2]+"\t"+line[4]+"\t"+line[5]+"\t"+line[3]+"\t"+line[0]+"\t"+line[6]+"\t"+line[7]+"\t"+line[8]+"\n"
out = open("/scratch/jjang/PIGMENT_PROJECT/Pigment_RNA_ALL/stringtie/Merged_WangT_Iri_R2_RNA_AR021_GTTTCGGAA_S5_R1_trimmedAligned.sortedByCoord.out.gene.abundance.txt",'w')
out.write(outdata)