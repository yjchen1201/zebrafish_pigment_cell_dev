#example: python3 convert_DSS.py Merged_ALLRUNS_Iri_Rep1_trimmed_bismark_resorted.rmdup.bismark.cov
import sys
if len(sys.argv)!=2:
	sys.exit(__doc__)
seq = sys.argv[1]
seqname = seq.split("_")[2]
rep = seq.split("_")[3]
print(seqname,rep)
cpg_data = open("CpG_sites_zfish.bed",'r')
cpg_dic = {}
for line in cpg_data:
	line = line.split("\t")
	chrbase = str(line[0])+"."+str(line[1])
	if str(line[0]) not in cpg_dic.keys():
		cpg_dic[str(line[0])] = {}
		cpg_dic[str(line[0])][line[1]]= [line[0],line[1],0,0]
	else:
		cpg_dic[str(line[0])][line[1]]= [line[0],line[1],0,0]
data = open(seq, 'r')
out = ""
out += "chr"+"\t"+"pos"+"\t"+"N"+"\t"+"X"+"\n"
for line in data:
	line = line.rstrip().split("\t")
	if len(line[0]) <6: 
		if str(line[0]) != "chrM":
			if str(int(line[1])-2) in cpg_dic[line[0]].keys():
				totalcov = int(line[4])+int(line[5])
				Ccount = int(line[4])
				cpg_dic[line[0]][str(int(line[1])-2)][2]+=totalcov
				cpg_dic[line[0]][str(int(line[1])-2)][3]+=Ccount
			else:
				if str(int(line[1])-1) in cpg_dic[line[0]].keys():
					totalcov = int(line[4])+int(line[5])
					Ccount = int(line[4])
					cpg_dic[line[0]][str(int(line[1])-1)][2]+=totalcov
					cpg_dic[line[0]][str(int(line[1])-1)][3]+=Ccount
for i in sorted(cpg_dic.keys()):
	for pos in sorted(cpg_dic[i].keys()):
		if int(cpg_dic[i][pos][2]) >0:
			out += cpg_dic[i][pos][0]+"\t"+ cpg_dic[i][pos][1]+"\t"+str(cpg_dic[i][pos][2])+"\t"+str(cpg_dic[i][pos][3])+"\n"
outfile = open(seqname+"_"+rep+"_DSS.txt",'w')
outfile.write(out)
