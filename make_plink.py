#python make_plink.py  FR SS
import sys
import os
file_fr = sys.argv[1]
newfilename=file_fr.rstrip().split('/')[-1]
#print (newfilename)
file_ss = sys.argv[2]
dirname='./'+newfilename+'/plink_'+newfilename+'/'
try:
	os.mkdir(dirname)
except OSError:
	pass
snp={}
samples_list=set()
snps_list=set()
chrom={}
poz={}
for_tped={}
with open(dirname+newfilename+'.tped','w') as tped:
	with open(dirname+newfilename+'.tfam','w') as tfam:
		for l in open(file_fr): #Read FinalReport
			data=l.rstrip().split('\t')

			samples_list.add(data[1])
			snps_list.add(data[0])

			snp[data[0]+'\t'+data[1]]=data[2]+'\t'+data[3]

		samples_list=tuple(samples_list)
		snps_list=tuple(snps_list)
		#print (len(samples_list))
		#print (len(snp))

		for i in snps_list:
			for j in samples_list:
				if i+'\t'+j not in snp:
					snp[i+'\t'+j]='-'+'\t'+'-'
				if i not in for_tped:
					for_tped[i]=[snp[i+'\t'+j]]
				else:
					for_tped[i].append(snp[i+'\t'+j])

		for l in open(file_ss):#Find Chr and pos for snps
				data=l.rstrip().split('\t')
				chrom[data[0]]=data[1]
				poz[data[0]]=data[2]
		#print (len(chrom))

		for l in for_tped:#write .tped for all samples
			if l in chrom:
				wr_tped=[chrom[l],l,'0',poz[l],'\t'.join(for_tped[l])]
				tped.write('\t'.join(wr_tped)+'\n')
		c=0
		for l in samples_list:#write .tfam for all samples
			c+=1
			wr_tfam=[str(c),l,'0','0','1','0']
			tfam.write('\t'.join(wr_tfam)+'\n')