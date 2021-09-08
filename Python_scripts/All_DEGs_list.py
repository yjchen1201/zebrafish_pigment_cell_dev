#python3 All_DEGs_list.py Input1_DEG_files.txt Input2_DEG_files.txt etc.
import sys
if len(sys.argv)!=5:
	sys.exit(__doc__)
outdata = ""
glist = []
for i in range(1,5):
	data = open(sys.argv[i],'r')
	print(data)
	for line in data:
		line = line.split("\t")
		if line[0] not in glist:
			glist.append(line[0])
			outdata += line[0]+"\n"
out = open("Combined_DEGs_list.txt",'w')
out.write(outdata)
