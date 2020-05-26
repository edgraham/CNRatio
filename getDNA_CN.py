#! /usr/bin/env python
from Bio import SeqIO
import os
import sys
import numpy as np
import collections
conversionDictN = {'C':3,'G':5,'A':5,'T':2}
conversionDictC = {'C':4,'G':4,'A':4,'T':4}
AAconversionDictN = {'G':1,'A':1,'V':1,'C':1,'P':1,'L':1,'I':1,'M':1,'W':2,'F':1,'K':2,'R':4,'H':3,'S':1,'T':1,'Y':1,'N':2,'Q':2,'D':1,'E':1}
AAconversionDictC = {'G':2,'A':3,'V':5,'C':3,'P':5,'L':6,'I':6,'M':5,'W':11,'F':9,'K':6,'R':6,'H':5,'S':3,'T':4,'Y':9,'N':4,'Q':5,'D':4,'E':5}
AAList = ['G','A','V','C','P','L','I','M','W','F','K','R','H','S','T','Y','N','Q','D','E']
def getCarbonNitrogenCount(g):
	C = 0
	N = 0
	Total = 0
	CN_ALL_CONTIGS =[]
	for record in SeqIO.parse(g,"fasta"):
		tmp = [y for y in str(record.seq)]
		counttmp = dict(collections.Counter(tmp))
		C = C+(counttmp['C']*4)+(counttmp['G']*4)+(counttmp['A']*4)+(counttmp['T']*4)
		N = N+(counttmp['C']*3)+(counttmp['G']*5)+(counttmp['A']*5)+(counttmp['T']*2)
		Ctmp = (counttmp['C']*4)+(counttmp['G']*4)+(counttmp['A']*4)+(counttmp['T']*4)
		Ntmp = (counttmp['C']*3)+(counttmp['G']*5)+(counttmp['A']*5)+(counttmp['T']*2)
		CN_tmp = float(Ctmp)/float(Ntmp)
		CN_ALL_CONTIGS.append(CN_tmp)
	CN_Total = float(C)/float(N)
	CN_Mean = sum(CN_ALL_CONTIGS)/len(CN_ALL_CONTIGS)
	Min =min(CN_ALL_CONTIGS)
	Max = max(CN_ALL_CONTIGS)
	return CN_Total, CN_Mean, Min, Max
def getCarbonNitrogenGeneCalls(g):
	C = 0
	N = 0
	perGeneC=[]
	perGeneN=[]
	PerGeneCN = []
	for record in SeqIO.parse(g,"fasta"):
                tmp = [y for y in str(record.seq)]
		counttmp = dict(collections.Counter(tmp))
		tmpDictC ={}
		tmpDictN = {}
		for k in AAList:
			try:
				tmpDictC[k] = counttmp[k]
				tmpDictN[k]=counttmp[k]
			except KeyError:
				continue
		tmpDictC.update((x,y*AAconversionDictC[x]) for x,y in tmpDictC.iteritems())
		C =C + sum(tmpDictC.values())
		N = N + sum(tmpDictN.values())
		Ctmp = sum(tmpDictC.values())
		Ntmp = sum(tmpDictN.values())
		perGeneC.append(float(Ctmp))
		perGeneN.append(float(Ntmp))
		PerGeneCN.append(float(Ctmp)/float(Ntmp))
	CN_Total = float(C)/float(N)
	CN_Mean = sum(PerGeneCN)/len(PerGeneCN)
	MinCN = min(PerGeneCN)
	MaxCN = max(PerGeneCN)
	N_total=N
	C_total=C
	N_Min=min(perGeneN)
	C_Min=min(perGeneC)
	N_Max=max(perGeneN)
	C_Max=max(perGeneC)
	N_Mean=sum(perGeneN)/len(perGeneN)
	C_Mean=sum(perGeneC)/len(perGeneC)
	return CN_Total, CN_Mean, MinCN, MaxCN, N_total, C_total, N_Mean, C_Mean, C_Min, N_Min, N_Max, C_Max

geneCallList = [f for f in os.listdir(sys.argv[1]) if f.endswith('.faa')]

output2 = open('CN_AA.txt','a')
output2.write('Genome\tCN_Total\tCN_Mean\tMinCN\tMaxCN\tN_Total\tC_Total\tN_Mean\tC_Mean\tC_Min\tN_Min\tN_Max,\tC_Max\n')
for genome in geneCallList:
        genomeName = genome.rsplit('.',1)[0]
	print genomeName
        CN_Total, CN_Mean, MinCN, MaxCN, N_total, C_total, N_Mean, C_Mean, C_Min, N_Min, N_Max, C_Max = getCarbonNitrogenGeneCalls(genome)
        output2.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(genomeName,CN_Total, CN_Mean, MinCN, MaxCN, N_total, C_total, N_Mean, C_Mean, C_Min, N_Min, N_Max, C_Max))
