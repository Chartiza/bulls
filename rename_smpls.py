#!/usr/bin/python

sootv = {}

#Read file sootvetstviya
for l in open ("filesootv"):
	data = l.strip().split("\t")
	if data[0] not in sootv:
		sootv[data[0]] = data[1]

#Read FinalReport file
for l in open('Ire30_GP'):
	data = l.strip().split("\t")
	if data[1] in sootv:
		print(data[0]+"\t"+sootv[data[1]]+"\t"+data[2]+"\t"+data[3]+"\t"+"\t"+data[4]+"\t"+data[5])
		