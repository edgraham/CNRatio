#! /usr/bin/env python
from Bio import SeqIO
import argparse
import numpy as np
import textwrap, glob
import pandas as pd

SynonymousCodons = { 
'CYS': ['TGT', 'TGC'], 
'ASP': ['GAT', 'GAC'], 
'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
'GLN': ['CAA', 'CAG'], 
'MET': ['ATG'], 
'ASN': ['AAC', 'AAT'], 
'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
'LYS': ['AAG', 'AAA'], 
 'STOP': ['TAG', 'TGA', 'TAA'], 
'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
'PHE': ['TTT', 'TTC'], 
'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
'ILE': ['ATC', 'ATA', 'ATT'], 
'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
'HIS': ['CAT', 'CAC'], 
'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
'TRP': ['TGG'], 
'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
'GLU': ['GAG', 'GAA'], 
'TYR': ['TAT', 'TAC'] 
} 

CodonDict = {'CYS':'TGT', 'CYS':'TGC', 'ASP':'GAT','ASP': 'GAC', 'SER':'TCT', 'SER':'TCG', 'SER':'TCA', 'SER':'TCC', 'SER':'AGC', 'SER':'AGT','GLN':'CAA', 'GLN':'CAG', 'MET':'ATG', 'ASN':'AAC', 'ASN':'AAT','PRO':'CCT', 'PRO':'CCG', 'PRO':'CCA', 'PRO':'CCC', 'LYS':'AAG', 'LYS':'AAA', 'STOP':'TAG','STOP':'TGA','STOP':'TAA', 'THR':'ACC', 'THR':'ACA', 'THR':'ACG', 'THR':'ACT', 'PHE':'TTT', 'PHE':'TTC', 'ALA':'GCA', 'ALA':'GCC', 'ALA':'GCG', 'ALA':'GCT','GLY':'GGT', 'GLY':'GGG', 'GLY':'GGA', 'GLY':'GGC','ILE':'ATC', 'ILE':'ATA', 'ILE':'ATT','LEU':'TTA', 'LEU':'TTG', 'LEU':'CTC', 'LEU':'CTT', 'LEU':'CTG', 'LEU':'CTA','HIS':'CAT', 'HIS':'CAC','ARG':'CGA', 'ARG':'CGC', 'ARG':'CGG', 'ARG':'CGT', 'ARG':'AGG', 'ARG':'AGA', 'TRP': 'TGG', 'VAL':'GTA', 'VAL':'GTC', 'VAL':'GTG', 'VAL':'GTT','GLU':'GAG', 'GLU':'GAA','TYR':'TAT', 'TYR':'TAC'} 
#ADD IN EFFECTIVE NUMBER OF CODONS Nc (Wright 1990)
#GC Percent
#GC CONTENT OF 3RD POSITION OF SYNONYMOUS CODONS (GC3s)
#CALCULATE USAGE OF EACH BASE AT WOBBLE
CodonDict ={'TGT':0, 'TGC':0,'GAT':0, 'GAC':0,'TCT':0, 'TCG':0, 'TCA':0, 'TCC':0, 'AGC':0, 'AGT':0,'CAA':0, 'CAG':0,'ATG':0,'AAC':0, 'AAT':0,'CCT':0, 'CCG':0, 'CCA':0, 'CCC':0,'AAG':0, 'AAA':0,'TAG':0, 'TGA':0, 'TAA':0,'ACC':0, 'ACA':0, 'ACG':0, 'ACT':0,'TTT':0, 'TTC':0,'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0,'GGT':0, 'GGG':0, 'GGA':0, 'GGC':0,'ATC':0, 'ATA':0, 'ATT':0,'TTA':0, 'TTG':0, 'CTC':0, 'CTT':0, 'CTG':0, 'CTA':0,'CAT':0, 'CAC':0,'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0, 'AGG':0, 'AGA':0,'TGG':0,'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0,'GAG':0, 'GAA':0,'TAT':0, 'TAC':0}
def find_len_syn(term):
    x = list(SynonymousCodons.values())
    for y in x:
        if term in y:
            return len(y), y

def BuildCodonFrequencies_all(AA,con):
    codon_table = ['TGT','TGC','GAT','GAC','TCT','TCG','TCA','TCC','AGC','AGT','CAA','CAG','ATG','AAC','AAT','CCT','CCG','CCA','CCC','AAG','AAA','TAG','TGA','TAA','ACC','ACA','ACG','ACT','TTC','TTT','GCA','GCC','GCG','GCT','GGT','GGG','GGA','GGC','ATC','ATA','ATT','TTA','TTG','CTC','CTT','CTG','CTA','CAT','CAC','CGA','CGC','CGG','CGT','AGG','AGA','TGG','GTA','GTC','GTG','GTT','GAG','GAA','TAT','TAC']
    dataframe_freq = pd.DataFrame()
    dataframe_rscu = pd.DataFrame()
    dataframe_raw = pd.DataFrame()
    dataframe_start=pd.DataFrame()
    for handle in glob.glob('*'+str(con)):
        organism = handle.rsplit(".",1)[0]
        amino=dict(AA)
        keys = amino.keys()
        count = 0
        amino_count = 0
        start = dict()
        for seq_record in SeqIO.parse(handle, "fasta") :
                count += 1
                seq = str(seq_record.seq)
                orfcodons = textwrap.wrap(seq, 3)
                start_codon = orfcodons[0]
                if set(codon_table)> set(orfcodons):
                    for item in orfcodons:
                        if str(item) in keys:
                            amino[item]+=1
                            amino_count +=1
                    if start_codon in start.keys():
                        start[start_codon] += 1
                    else:
                        start.setdefault(start_codon,0)
                        start[start_codon]+=1    
        rscu = dict()
        if all(value == 0 for value in amino.values()) == True:
            for y in amino.keys():
                if amino[y] == 0:
                    rscu[y] = 0
        else:    
            for y in amino.keys():
                if amino[y] == 0:
                    rscu[y] = 0
                else:
                    length, synonymous_list = find_len_syn(str(y))
                    rscu[y] = float(amino[y])/((1.0/float(length))*((sum(amino[a]for a in synonymous_list))))
        amino_freq = dict()
        if all(value == 0 for value in amino.values()) == True:
            for y in amino.keys():
                if amino[y] == 0:
                    amino_freq[y] = 0
        else:    
            for y in amino.keys():
                if amino[y] == 0:
                    amino_freq[y] = 0
                else:
                    length, synonymous_list = find_len_syn(str(y))
                    amino_freq[y] = float(amino[y])/((sum(amino[a]for a in synonymous_list)))
        start_freq = dict()
        for y in start.keys():
            start_freq[y]=float(float(start[y])/float(sum(start.values())))
        amino_freq_mapped = map(list,amino_freq.items())                    
        rscu_mapped = map(list,rscu.items())
        amino_mapped = map(list,amino.items())
        start_mapped = map(list,start_freq.items())
        appended_start = []
        appended_rscu = []
        appended_raw_amino = []
        appended_freq=[]
        for codon in start_mapped:
            codon.insert(0,str(organism))
            appended_start.append(codon)
        for codon in amino_mapped:
            codon.insert(0,str(organism))
            appended_raw_amino.append(codon)
        for codon in rscu_mapped:
            codon.insert(0,str(organism))
            appended_rscu.append(codon)
        for codon in amino_freq_mapped:
            codon.insert(0,str(organism))
            appended_freq.append(codon)       
        dictionary_for_rscu = {}
        for org, aminoacid,freq in appended_rscu:
            dictionary_for_rscu.setdefault(org,{})[aminoacid]=freq
        dataframe1 = pd.DataFrame.from_dict(dictionary_for_rscu, orient='index')
        dataframe_rscu = dataframe_rscu.append(dataframe1)
        dictionary_for_raw = {}
        for org, aminoacid,freq in appended_raw_amino:
            dictionary_for_raw.setdefault(org,{})[aminoacid]=freq
        dataframe2 = pd.DataFrame.from_dict(dictionary_for_raw, orient='index')
        dataframe_raw = dataframe_raw.append(dataframe2)
        dictionary_for_freq = {}
        for org, aminoacid,freq in appended_freq:
            dictionary_for_freq.setdefault(org,{})[aminoacid]=freq
        dataframe3 = pd.DataFrame.from_dict(dictionary_for_freq, orient='index')
        dataframe_freq = dataframe_freq.append(dataframe3)
        dictionary_for_start = {}
        for org, aminoacid,freq in appended_start:
            dictionary_for_start.setdefault(org,{})[aminoacid]=freq
        dataframe4 = pd.DataFrame.from_dict(dictionary_for_start, orient='index')
        dataframe_start = dataframe_start.append(dataframe4)     
        print "There are %s proteins in organism " % count + str(organism) + " with %s total codons" % amino_count
        count = 0  
        amino_count = 0    
    y = list(dataframe_rscu)
    dataframe_start =dataframe_start.replace(r'\s+', np.nan, regex=True)
    dataframe_start = dataframe_start.fillna(0)  
    reorder_rscu = dataframe_rscu.reindex(columns=['TGT','TGC','GAT','GAC','TCT','TCG','TCA','TCC','AGC','AGT','CAA','CAG','ATG','AAC','AAT','CCT','CCG','CCA','CCC','AAG','AAA','TAG','TGA','TAA','ACC','ACA','ACG','ACT','TTC','TTT','GCA','GCC','GCG','GCT','GGT','GGG','GGA','GGC','ATC','ATA','ATT','TTA','TTG','CTC','CTT','CTG','CTA','CAT','CAC','CGA','CGC','CGG','CGT','AGG','AGA','TGG','GTA','GTC','GTG','GTT','GAG','GAA','TAT','TAC']).transpose() 
    reorder_raw = dataframe_raw.reindex(columns=['TGT','TGC','GAT','GAC','TCT','TCG','TCA','TCC','AGC','AGT','CAA','CAG','ATG','AAC','AAT','CCT','CCG','CCA','CCC','AAG','AAA','TAG','TGA','TAA','ACC','ACA','ACG','ACT','TTC','TTT','GCA','GCC','GCG','GCT','GGT','GGG','GGA','GGC','ATC','ATA','ATT','TTA','TTG','CTC','CTT','CTG','CTA','CAT','CAC','CGA','CGC','CGG','CGT','AGG','AGA','TGG','GTA','GTC','GTG','GTT','GAG','GAA','TAT','TAC']).transpose() 
    reorder_freq = dataframe_freq.reindex(columns=['TGT','TGC','GAT','GAC','TCT','TCG','TCA','TCC','AGC','AGT','CAA','CAG','ATG','AAC','AAT','CCT','CCG','CCA','CCC','AAG','AAA','TAG','TGA','TAA','ACC','ACA','ACG','ACT','TTC','TTT','GCA','GCC','GCG','GCT','GGT','GGG','GGA','GGC','ATC','ATA','ATT','TTA','TTG','CTC','CTT','CTG','CTA','CAT','CAC','CGA','CGC','CGG','CGT','AGG','AGA','TGG','GTA','GTC','GTG','GTT','GAG','GAA','TAT','TAC']).transpose() 
    reorder_rscu["AA"] = pd.Series(CodonDict,dtype=np.unicode_)
    reorder_raw["AA"] = pd.Series(CodonDict,dtype=np.unicode_)
    reorder_freq["AA"] = pd.Series(CodonDict,dtype=np.unicode_)
    start_tab = dataframe_start.to_csv("start_codon_freq.txt",sep="\t")
    rscu_tab = reorder_rscu.to_csv("rscu.txt",sep="\t")
    raw_tab = reorder_raw.to_csv("raw_codon_count.txt",sep="\t")
    freq_tab = reorder_freq.to_csv("codon_freq.txt",sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Codon_Bias', usage='%(prog)s -f [fa,fasta,fna]')
    parser.add_argument("-f", dest="inputFile", help="Specify suffix linking fasta files (e.g fa,fna,faa)")
    args = parser.parse_args()
    print('Calculating Overall Codon Usage Frequencies')
    print('********************************************')
    BuildCodonFrequencies_all(CodonDict,args.inputFile)
