#-*- coding: utf-8 -*-
import pandas as pd
import sys
import os
from itertools import islice
from subprocess import call
print 'Укажите путь до FR-filename:'
#finalreport = ''
finalreport = raw_input()
frfilename = finalreport.rstrip().split('/')[-1]

try:
	os.makedirs('./'+frfilename+'/plink_'+frfilename)
except OSError:
	pass
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
    print (' '.join(command))
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
	sys.stdout = open('./'+frfilename + '/' + frfilename + '.log', 'w')
	p = 'Genotypes batch: ' + frfilename
	print("".join(p))
	p = 'Bad qual Sampls, CR<90: ' + str(len(bcr.index))
	print("".join(p))
	p = 'Bad qual SNP: ' + str(len(bsr.index))
	print("".join(p))
	p = 'SNP with MAF < 0.05: ' + str(len(bmf.index))
	print("".join(p))
	p = 'SNP with HWEdiff > 0.15: ' + str(len(bhw.index))
	print("".join(p))
	#print bad samples to file
	sys.stdout = open('./'+frfilename + '/BadSmpl.' + frfilename + '.csv', 'w')
	print(bcr)

king()
stat()
