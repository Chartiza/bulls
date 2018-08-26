#-*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sys
import os
from itertools import islice
from subprocess import call
print 'Укажите путь до FR-filename:'
finalreport = '/5500disk/home/chara/DBs/DB_genotypes/FinalReports/FRgood_20180718'
#finalreport = raw_input()
frfilename = finalreport.rstrip().split('/')[-1]
partia = sys.argv[1]

try:
	os.makedirs('./'+frfilename+'/plink_'+frfilename)
except OSError:
	pass
#_____________________GrepSmpls____________________
def makeFR(list1, frfilename):
    print ('Make newFR...')
    smpl=set()
    dirname='./'+frfilename+'/plink_'+frfilename+'/'
    f = open(dirname+frfilename, 'w')
    for l in open('./'+list1):
        data=l.strip().split('\t')
        if data[1] == 'frstlist':
            smpl.add(data[0])
    for l in open(finalreport):
        data = l.strip().split('\t')
        if data[1] in smpl:
            f.write(l)
    f.close()

    print ('Запуск KING...')
    #make_tped_tfam
    command = ('python make_plink.py ' + dirname+frfilename + ' /5500disk/home/chara/DBs/DB_genotypes/SampleSheets/SS_DB').split(' ')
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
    #king
    bed=dirname+frfilename+'.bed'
    command = ('king -b '+bed +' --kinship --prefix '+dirname+frfilename+'--related --degree 3').split(' ')
    print (' '.join(command))
    call(command)
	#сравнить КК
    print ('Make compare..')
    #make res file
    f = open('./'+frfilename+'/KKcompare.csv', 'w')
	#read kinsip value
    kval={}
    for l in open('./pattern/Kinship.csv'):
        data = l.strip().split('\t')
        b = ''.join([data[0]]) + ''.join([data[1]])
    	if b not in kval:
        	kval[b]=''.join([data[2]])
	#read king value
    for l in open(dirname+frfilename+'--related.kin0'):
        data = l.strip().split('\t')
        if data[0] != 'FID1':
			a = ''.join([data[1]]) + ''.join([data[3]])
			b = ''.join([data[3]]) + ''.join([data[1]])
			if b in kval:
				f.write(data[1] + '\t' + data[3] + '\t' + ''.join(kval[b]) + '\t' + data[7] +'\n')
			elif a in kval:
				f.write(data[1] + '\t' + data[3] + '\t' + ''.join(kval[a]) + '\t' + data[7] +'\n')
			else:
				f.write(data[1] + '\t' + data[3] + '\t' + '0' + '\t' + data[7]+'\n')
    f.close()
    #Print diff between king and kinship
    fr = pd.read_csv('./'+frfilename+'/KKcompare.csv', sep = '\s+', index_col = False, header=None)
    k2 = []
    k3 = []

    for row in fr[2]:
    	if row < 0.0442:
        	k2.append('0')
    	elif row > 0.0442 and row < 0.0884:
        	k2.append('3rd')
    	elif row > 0.0884 and row < 0.177:
        	k2.append('2nd')
    	elif row > 0.177 and row < 0.354:
        	k2.append('1st')
    	elif row > 0.354:
        	k2.append('rep')
    	else:
        	k2.append('else')

    for row in fr[3]:
    	if row < 0.0442:
        	k3.append('0')
    	elif row > 0.0442 and row < 0.0884:
        	k3.append('3rd')
    	elif row > 0.0884 and row < 0.177:
        	k3.append('2nd')
    	elif row > 0.177 and row < 0.354:
        	k3.append('1st')
    	elif row > 0.354:
        	k3.append('rep')
    	else:
        	k3.append('else')

    fr['kinship'] = k2
    fr['king'] = k3
    
    #calculate number good/bad pairs
    fr['kinship'].equals(fr['king'])
    fr['res'] = np.where(fr['kinship']!=fr['king'], 'False', 'True')
    fr.columns = ['anim1', 'anim2', 'kinshipVal', 'kingVal', 'kinship', 'king', 'result']
    fr.to_csv('./'+frfilename+'/DEGREE.csv', sep='\t', encoding='utf-8', index = False)
    
    bd = fr.loc[fr['result'] == 'False']
    bd1 = bd[['anim1', 'anim2']].melt()
    bad = pd.DataFrame(bd1.value.value_counts())
    
    gd = fr.loc[fr['result'] == 'True']
    gd1 = gd[['anim1', 'anim2']].melt()
    good = pd.DataFrame(gd1.value.value_counts())

    ls = pd.read_table('pattern/IlyaGtList-1.txt', header=None).set_index([0])
    r = pd.concat([bad,good,ls], axis = 1).fillna(0)
    r.columns = ['bad', 'good', 'group']
    r['perc']=r['bad']/(r['bad']+r['good'])
    r['tgroup']=np.where(r['group'] == 'notinfirst', 'NA', 'test')
    r['tgroup']=np.where((r['group'] != 'notinfirst') & (r['good'] == 0) & (r['bad'] == 0), 'NoData', r['tgroup'])
    r.to_csv('./'+frfilename+'_COMPARE_'+partia+'.csv', sep='\t', encoding='utf-8')
    
    a = r['tgroup'].value_counts()
    print(a)

    #def the worst
    pl = pd.read_csv('./'+frfilename+'_COMPARE_'+partia+'.csv', sep='\s+')
    wr = pl.loc[pl['perc'] != 'NoData']
    wr1 = wr.sort_values(['perc', 'bad'], ascending=False)
    w = wr1.index[0]
    print(w)
    f = open('./'+frfilename+'_worst.csv', 'w')
    f.write(str(w))
    f.close()   

makeFR('pattern/IlyaGtList-1.txt', frfilename)
#king()
