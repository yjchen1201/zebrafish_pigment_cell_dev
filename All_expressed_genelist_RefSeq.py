#python3 All_expressed_genelist_RefSeq.py Input1.gene.abundance.filteredNM.txt Input2.gene.abundance.filteredNM.txt etc
import sys
if len(sys.argv)!=9:
	sys.exit(__doc__)
outdata = ""
glist = []
for i in range(1,9):
	data = open(sys.argv[i],'r')
	print(data)
	for line in data:
		line = line.split("\t")
		if line[0] not in glist:
			glist.append(line[0])
			outdata += line[0]+"\n"
out = open("Combined_gene_list_RefSeq.txt",'w')
out.write(outdata)