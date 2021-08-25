#prints out 1)ID   2)Alx1 motif   3)ALX3 motif   4)Alx4 motif   5)GBX2 motif   6)PAX7 motif   7)SOX9 motif   8)SOX10 motif     9)Sox5 motif   10)ELK motif   11)ETS1  12)TFEC  13)HEY1 14)Mycn 15)SP3 16)Nfe2l2 17)Gabpa 18)TFAP2A 19)BHLHE41 20)JUNB 21)NR4A2 22)NR4A1
import sys
if len(sys.argv)!=2:
	sys.exit(__doc__)
seq = sys.argv[1]
data = open(seq, "r")
outdata = ""
outdata+= "chrompos"+"\t"+"ALX1"+"\t"+"ALX3"+"\t"+"ALX4"+"\t"+"GBX2"+"\t"+"PAX7"+"\t"+"SOX9"+"\t"+"SOX10"+"\t"+"ELK4"+"\t"+"ETS1"+"\t"+"GABPA"+"\t"+"NFE2l2"+"\t"+"TFEC"+"\t"+"HEY1"+"\t"+"MYCN"+"\t"+"BHLHE41"+"\t"+"TFAP2A"+"\t"+"SP3"+"\t"+"JUNB"+"\t"+"NR4A2"+"\t"+"NR4A1"+"\n"
dic = {}
for line in data:
	if line[0] == "#":
		pass
	else:
		line = line.rstrip().split("\t")
		if line[2] not in dic.keys():
			dic[line[2]] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			if line[1] == "Alx1":
				dic[line[2]][0]+=1
			if line[1] == "Alx3":
				dic[line[2]][1]+=1
			if line[1] == "Alx4":
				dic[line[2]][2]+=1
			if line[1] == "Gbx2":
				dic[line[2]][3]+=1
			if line[1] == "Pax7":
				dic[line[2]][4]+=1
			if line[1] == "SOX9":
				dic[line[2]][5]+=1
			if line[1] == "SOX10":
				dic[line[2]][6]+=1
			if line[1] == "ELK4":
				dic[line[2]][7]+=1
			if line[1] == "ETS1":
				dic[line[2]][8]+=1
			if line[1] == "TFEC":
				dic[line[2]][11]+=1
			if line[1] == "HEY1":
				dic[line[2]][12]+=1
			if line[1] == "Mycn":
				dic[line[2]][13]+=1
			if line[1] == "SP3":
				dic[line[2]][16]+=1
			if line[1] == "Nfe2l2":
				dic[line[2]][10]+=1
			if line[1] == "Gabpa":
				dic[line[2]][9]+=1
			if line[1] == "TFAP2A":
				dic[line[2]][15]+=1
			if line[1] == "BHLHE41":
				dic[line[2]][14]+=1
			if line[1] == "JUNB":
				dic[line[2]][17]+=1
			if line[1] == "NR4A2":
				dic[line[2]][18]+=1
			if line[1] == "NR4A1":
				dic[line[2]][19]+=1
		else:
			if line[1] == "Alx1":
				dic[line[2]][0]+=1
			if line[1] == "Alx3":
				dic[line[2]][1]+=1
			if line[1] == "Alx4":
				dic[line[2]][2]+=1
			if line[1] == "Gbx2":
				dic[line[2]][3]+=1
			if line[1] == "Pax7":
				dic[line[2]][4]+=1
			if line[1] == "SOX9":
				dic[line[2]][5]+=1
			if line[1] == "SOX10":
				dic[line[2]][6]+=1
			if line[1] == "ELK4":
				dic[line[2]][7]+=1
			if line[1] == "ETS1":
				dic[line[2]][8]+=1
			if line[1] == "TFEC":
				dic[line[2]][11]+=1
			if line[1] == "HEY1":
				dic[line[2]][12]+=1
			if line[1] == "Mycn":
				dic[line[2]][13]+=1
			if line[1] == "SP3":
				dic[line[2]][16]+=1
			if line[1] == "Nfe2l2":
				dic[line[2]][10]+=1
			if line[1] == "Gabpa":
				dic[line[2]][9]+=1
			if line[1] == "TFAP2A":
				dic[line[2]][15]+=1
			if line[1] == "BHLHE41":
				dic[line[2]][14]+=1
			if line[1] == "JUNB":
				dic[line[2]][17]+=1
			if line[1] == "NR4A2":
				dic[line[2]][18]+=1
			if line[1] == "NR4A1":
				dic[line[2]][19]+=1
for i in dic.keys():
	outdata += i+"\t"+str(dic[i][0])+"\t"+str(dic[i][1])+"\t"+str(dic[i][2])+"\t"+str(dic[i][3])+"\t"+str(dic[i][4])+"\t"+str(dic[i][5])+"\t"+str(dic[i][6])+"\t"+str(dic[i][7])+"\t"+str(dic[i][8])+"\t"+str(dic[i][9])+"\t"+str(dic[i][10])+"\t"+str(dic[i][11])+"\t"+str(dic[i][12])+"\t"+str(dic[i][13])+"\t"+str(dic[i][14])+"\t"+str(dic[i][15])+"\t"+str(dic[i][16])+"\t"+str(dic[i][17])+"\t"+str(dic[i][18])+"\t"+str(dic[i][19])+"\n"
out = open(seq.strip(".txt")+"_parsed.txt",'w')
out.write(outdata)
		