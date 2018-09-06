#-*- coding: utf-8 -*-
import pandas as pd
import sys
import csv
import os
from itertools import islice
from subprocess import call

print 'Укажите путь до FR-filename:'
finalreport = '/5500disk/home/chara/DBs/DB_genotypes/FinalReports/FRgood_20180718'
#finalreport = raw_input()
frfilename = finalreport.rstrip().split('/')[-1] #comFR_FRgood_20180718

try:
	#os.makedirs('./' + frfilename + '/data')
	os.makedirs('./'+frfilename+'/plink_'+frfilename)
except OSError:
	pass

#_______________commonSNP___________________________
def common():
    print ('Choose common SNP...')
    comsnp=set()
    dirname='./'+frfilename+'/'
    #Read needed SNP list
    for l in open('/5500disk/home/chara/DBs/DB_genotypes/SampleSheets/SNP_groups/CommonSNP_all_chip'):
	data=l.rstrip().split('\t')
	if data[0] not in comsnp:
		comsnp.add(data[0])
    #Read FinalReport
    f = open('./'+frfilename + '/' + frfilename, 'w')
    for l in open(finalreport):
	data=l.rstrip().split('\t')
	if data[0] in comsnp:
		f.write(l)
#____________________KING___________________________
def king():
    print ('Запуск KING...')
    #make_tped_tfam
    command = ('python make_plink.py ' + finalreport+' /5500disk/home/chara/DBs/DB_genotypes/SampleSheets/SS_DB').split(' ')
    print (' '.join(command))
    call(command)
    #plink
    dirname='./'+frfilename+'/plink_'+frfilename+'/'
    tped=dirname+frfilename+'.tped'
    tfam=dirname+frfilename+'.tfam'
    #command = ('plink --tped '+tped+' --tfam '+tfam+' --missing --chr-set 30 --out '+dirname+frfilename).split(' ')
    command = ('plink --tped '+tped+' --tfam '+tfam+' --chr-set 30 --out '+dirname+frfilename+' --missing-genotype -').split(' ')
    print (' '.join(command))
    call(command)
    #check_CallRate, MAF and HW
    bed=dirname+frfilename
    command = ('plink --bfile '+bed+' --missing --hardy --freq --chr-set 30 --make-bed --out '+dirname+frfilename).split(' ')
    print (' '.join(command))
    call(command)
    #king
    bed=dirname+frfilename+'.bed'
    command = ('king -b '+bed +' --kinship --prefix '+dirname+frfilename+'--related --degree 3').split(' ')
    #print (' '.join(command))
    #call(command)

#____________________stat________________________
def stat():
	print ('Запуск расчета статистики...')
	dirname='./'+frfilename+'/plink_'+frfilename+'/'
	cr = pd.read_csv(dirname + frfilename + '.imiss', sep = '\s+') #callrate smpls
	sr = pd.read_csv(dirname + frfilename + '.lmiss', sep = '\s+') #callrate SNPs
	mf = pd.read_csv(dirname + frfilename + '.frq', sep = '\s+') #MAF
	hw = pd.read_csv(dirname + frfilename + '.hwe', sep = '\s+') #HWE
	bcr = cr.loc[cr['F_MISS']>0.09]
	#calculate bad values
	bsr = sr.loc[sr['F_MISS']>0.091]
	bmf = mf.loc[mf['MAF']<0.05]
	hw.columns = ['CHR', 'SNP', 'TEST', 'A1', 'A2', 'GENO', 'OHET', 'EHET', 'P']
	hw['dif'] = hw.OHET-hw.EHET
	bhw = hw.loc[hw['dif']>0.15]
	#print stat to file
	f = open('./'+frfilename + '/' + frfilename + '.log', 'w')
	p = 'Genotypes batch: ' + frfilename + '\n'
	f.write("".join(p))
	p = 'Bad qual Sampls, CR<90: ' + str(len(bcr.index)) + '\n'
	f.write("".join(p))
	p = 'Bad qual SNP: ' + str(len(bsr.index)) + '\n'
	f.write("".join(p))
	p = 'SNP with MAF < 0.05: ' + str(len(bmf.index)) + '\n'
	f.write("".join(p))
	p = 'SNP with HWEdiff > 0.15: ' + str(len(bhw.index)) + '\n'
	f.write("".join(p))
	#print bad samples to file
	bcr.to_csv('./'+frfilename + '/BadSmpl.' + frfilename + '.csv', sep = '\t', index = None)

def newFR():
	smpl_to_cut=set()
	#dirname='./'+frfilename+'/plink_'+frfilename+'/'
	print ('Eeclude bad smpls >>> ')
	#Read BedSMPLS.csv
	for l in open('./'+frfilename + '/BadSmpl.' + frfilename + '.csv'):
		data=l.rstrip().split('\t')
		smpl_to_cut.add(data[1])
		print(data[1])
	#Read FinalReport
	f = open(finalreport, 'r')
	lines = f.readlines()
	f.close()
	f = open(finalreport, 'w')
	for l in lines:
		data=l.rstrip().split('\t')
		if data[1] not in smpl_to_cut:
			#print (data)
			f.write(l)
	print ('Eeclude bad SNPs >>> ')
	for l in open('./'+frfilename + '/BadSmpl.' + frfilename + '.csv'):
		data=l.rstrip().split('\t')
		smpl_to_cut.add(data[1])
		print(data[1])



common()
finalreport = './'+frfilename + '/' + frfilename
king()
stat()
newFR()