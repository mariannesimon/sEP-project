#!/usr/bin/python

import os.path, sys, re, random, numpy as np
from Bio import SeqIO

# Set variables
try:
    elm_fname, seq_fname, out_name = sys.argv[1], sys.argv[2], sys.argv[3]
except:
    print "Problem with file names."
N = 500

# Shuffle sequence
def shuffle(seq):
    s = ""
    return s.join(random.sample(seq,len(seq)))

# Motif search
with open(elm_fname,"r") as elm_file, open(out_name,"w") as outfile:
    raw_elm = elm_file.readlines()[6:]
    elmlist = [[line.split("\t")[i][1:-1] for i in range(8)] for line in raw_elm]
    outfile.write("sPEP id\tsPEP length\tStart\tEnd\tMatch\tELM_ID\tPvalue\n")
    # Read fasta file
    for spep in SeqIO.parse(seq_fname,"fasta"):
        pval, allmatch = [], []
        for elm in elmlist:
            matches = [(m.start(),m.end(),m.group(),elm[1]) for m in re.finditer(r"%s" %elm[4], str(spep.seq))]
            # Pvalue estimation for the motif if found
            if len(matches)>0:
                allmatch.append(matches)
                p = 0
                i = 0
                # Shuffle N times target sequence
                while i<N:
                    num = 0
                    for m in re.finditer(r"%s" %elm[4], shuffle(str(spep.seq))):
                        num += 1
                    p += int(num>=len(matches))
                    i += 1
                pval.append(p/float(N))
        # BH False discovery rate control
        for i in np.argsort(pval):
            if pval[i]<=(i+1)/float(len(pval))*0.05:
                for j in allmatch[i]: # one line for each elm
                    outfile.write("%s\t%d\t%d\t%d\t%s\t%s\t%0.4f\n" \
                    %(spep.id,len(str(spep.seq)),j[0],j[1],j[2],j[3],pval[i]))
