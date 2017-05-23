#!/usr/bin/python
# -*- coding: utf-8 -*-

import os.path, fileinput
from Bio import SeqIO

# Setting path and files
abspath = os.path.split(os.path.dirname(os.path.realpath(__file__)))
seqfile = "%s/Slavoff_cDNAs.fasta" %abspath[0]
outfile = open("%s/Primers.txt" %abspath[0],"w")
outfile.write("Name\tForward(5\'-3\')\tTm\tGC%\tLength\tReverse(5\'-3\')\tTm\tGC%\tLength\n")

# Input param values
Tminf = int(raw_input("Please enter minimum Tm (°C) : "))
lenmin = int(raw_input("Please enter minimum primer length (bp) : "))
lenmax = int(raw_input("Please enter maximum primer length (bp) : "))
length = range(lenmin, lenmax+1)
# Default param values
replace = {'A':'T','T':'A','C':'G','G':'C'}

# Primer search
for record in SeqIO.parse(seqfile,'fasta'):
    cdna = 'ATG' + str(record.seq[3:-3]) + 'TAG'
    rcdna = ''
    for n in cdna[::-1]:
        rcdna += replace[n]
    l = len(rcdna)
    sense, anti = [], []
    sd, ad = {'A':0,'T':0,'C':0,'G':0}, {}
    for nucl in sd.keys():
        sd[nucl] = cdna[:9].count(nucl)
        ad[nucl] = rcdna[:9].count(nucl)
    for i in length:
        sd[cdna[i-1]] += 1
        ad[rcdna[i-1]] += 1
        Tm = (sd['A']+sd['T'])*2 + (sd['C']+sd['G'])*4
        rTm = (ad['A']+ad['T'])*2 + (ad['C']+ad['G'])*4
        GC = (sd['C']+sd['G'])/float(i)*100
        rGC = (ad['C']+ad['G'])/float(i)*100
        if Tm>=Tminf:
            sense.append((cdna[:i],Tm,GC))
        if rTm>=Tminf:
            anti.append((rcdna[:i],rTm,rGC))
    outfile.write("%s\t"%record.id)
    if len(sense)>0:
        outfile.write("%s\t%d\t%0.1f\t%d\t" %(sense[0][0],sense[0][1],sense[0][2],len(sense[0][0])))
    if len(anti)>0:
        outfile.write("%s\t%d\t%0.1f\t%d\n" %(anti[0][0],anti[0][1],anti[0][2],len(anti[0][0])))
    else: outfile.write("\n")
outfile.close()
