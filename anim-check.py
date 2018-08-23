##############
#check new animals
#############

#python anim-check.py MGFanimals.txt z373.txt

import sys
import csv

animals={}

#Read MGFanimals.txt
for l in open(sys.argv[1]):
	data=l.rstrip().split('\t')
	if data[1] not in animals:
		animals[data[4]]=[data[1]+"\t"+data[6]]
		#print[data[4]]

#Read check animals list
for l in open(sys.argv[2]):
	data=l.rstrip().split('\t')
	if data[1] in animals:
		print ("For sample {} {}:".format(data[1],data[3]))
		print(animals[data[1]])
		


