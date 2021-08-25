import sys
if len(sys.argv)!=2:
	sys.exit(__doc__)
seq = sys.argv[1]
data = open(seq, "r")
outdata = ""
outdata+= "chrompos"+"\t"+"ALX1"+"\t"+"ALX3"+"\t"+"ALX4"+"\t"+"GBX2"+"\t"+"SOX10"+"\t"+"ETS1"+"\t"+"TFEC"+"\n"
dic = {}
for line in data:
	print(line)
	if line[0] == "m":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[2] not in dic.keys():
			dic[line[2]] = [0,0,0,0,0,0,0]
			if line[1] == "Alx1":
				dic[line[2]][0]+=1
			if line[1] == "ALX3":
				dic[line[2]][1]+=1
			if line[1] == "ALX4":
				dic[line[2]][2]+=1
			if line[1] == "GBX2":
				dic[line[2]][3]+=1
			if line[1] == "SOX10":
				dic[line[2]][4]+=1
			if line[1] == "ETS1":
				dic[line[2]][5]+=1
			if line[1] == "TFEC":
				dic[line[2]][6]+=1
		else:
			if line[1] == "Alx1":
				dic[line[2]][0]+=1
			if line[1] == "ALX3":
				dic[line[2]][1]+=1
			if line[1] == "ALX4":
				dic[line[2]][2]+=1
			if line[1] == "GBX2":
				dic[line[2]][3]+=1
			if line[1] == "SOX10":
				dic[line[2]][4]+=1
			if line[1] == "ETS1":
				dic[line[2]][5]+=1
			if line[1] == "TFEC":
				dic[line[2]][6]+=1
for i in dic.keys():
	outdata += i+"\t"+str(dic[i][0])+"\t"+str(dic[i][1])+"\t"+str(dic[i][2])+"\t"+str(dic[i][3])+"\t"+str(dic[i][4])+"\t"+str(dic[i][5])+"\t"+str(dic[i][6])+"\n"
out = open(seq.strip(".txt")+"_parsed.txt",'w')
out.write(outdata)
		