#-*- coding: utf-8 -*-
import pandas as pd
from subprocess import call

for i in range(2,101):
    pl = pd.read_csv('FRgood_20180718_worst.csv', sep='\s+', header = None)
    ls = pd.read_table('temp/IlyaGtList-1.txt', header = None).set_index([0])
    ls.columns=['group']
    a = int(pl.iloc[0])
    #print(a)
    ls.loc[a, 'group']='exclude'
    ls.to_csv('temp/IlyaGtList-1.txt', sep='\t', encoding='utf-8')
    #run main script
    command = ('python KKcompare.py '+str(i)).split(' ')
    print (' '.join(command))
    call(command)
