import sys
if len(sys.argv)!=2:
	sys.exit(__doc__)
infile = sys.argv[1]
data = open(infile,'r')
outdata = ""
for line in data:
	if line[0] == "N":
		outdata += line
outfile = infile.split(".txt")[0]+".filteredNM.txt"
out = open(outfile,'w')
out.write(outdata)
