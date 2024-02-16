import pickle, os, sys, threading, Bio
import numpy as np
from Bio import Entrez, SeqIO, pairwise2, AlignIO, Phylo
with open('Hwk4Data.pkl','rb')as f:dd=pickle.load(f)
with open('Hwk4Answers.pkl','rb')as f:ad=pickle.load(f)
Entrez.email = "brian.munsky@colostate.edu"  # Insert your email here
def get_genbank(accession_number):handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text");record = SeqIO.read(handle, "genbank");handle.close();return record
se = get_genbank(dd['accession']);q1 = len(se.seq); q2 = se.seq[:100] 
nt = 'ATGC';ct = [];lb = [];q3 = {}
sequences = {'seq1':'NM_004417.4', 'seq2':'NM_013642.3', 'seq3':'NM_001085359.3', 'seq4':'XM_038494850.1', 'seq5':'NM_001046452.2', 'seq6':'NM_001257450.2', 'seq7':'XM_030283971.3', 'seq8':'XM_002916919.4', 'seq9':'NM_001256075.1'}
for i in range(4):
    for j in range(4):
        lb.append(nt[i]+nt[j]);q3[nt[i]+nt[j]] = 0
for i in range(len(se.seq)-1): q3[se.seq[i:i+2]] += 1
tc = [];tl = [];q4 = {}
for i in range(4):
    for j in range(4):
        for k in range(4):
            tl.append(nt[i]+nt[j]+nt[k]);q4[nt[i]+nt[j]+nt[k]] = 0
for i in range(len(se.seq)-2):q4[se.seq[i:i+3]] += 1
def fao(se):
    o = [];ol = []      
    for f in range(3): 
        ls = f  
        for oi in se[f:].translate(to_stop=False).split('*'): 
            st = ls + oi.find('M') * 3 if 'M' in oi else None 
            sp = ls + (len(oi))*3 if 'M' in oi else None      
            if st is not None and sp is not None:            
                o.append(se[st:sp])               
                ol.append(len(oi))                  
            ls += len(oi)*3+3
    return o,ol 
o, ol = fao(se.seq);i = np.argmax(ol);
q5 = o[i];q6 = q5.translate() 
sd = {}; q8 = []; 
for k in sequences: sd[k] = get_genbank(sequences[k]); q8.append(sd[k].annotations['organism'])
def test1():assert ad['q1'] == q1
def test2():assert ad['q2'] == q2
def test3():assert ad['q3'] == q3
def test4():assert ad['q4'] == q4
def test5():assert ad['q5'] == q5
def test6():assert ad['q6'] == q6
def test8():assert ad['q8'] == q8