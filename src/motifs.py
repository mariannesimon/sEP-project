#!/usr/bin/python

import os.path, sys, re, random
from Bio import SeqIO

path = os.path.split(os.path.dirname(os.path.realpath(__file__)))
elm_fname = sys.argv[1]
seq_fname = sys.argv[2]
out_name = sys.argv[3]

def shuffle(seq):
    s = ""
    return s.join(random.sample(seq,len(seq)))

elm_original = False
# Modify original file
if elm_original:
    outfile = open("%s/motifs_elm.txt" %path[0],"w")
    with open("%s/data/elm_classes.tsv" %path[0],"r") as elmfile:
        raw_data = elmfile.readlines()[6:]
        for line in raw_data:
            line = line.split("\t")
            outfile.write("%s\t%s\t%s\n"%(line[4][1:-1],line[1][1:-1],line[5][1:-1]))
    outfile.close()

# Motif search
outfile = open("%s/%s" %(path[0],out_name),"w")
with open(elm_fname,"r") as elm_file:
    elmlist = [elmline.split("\t") for elmline in elm_file.readlines()]
    for spep in SeqIO.parse(seq_fname,"fasta"):
        outfile.write("# %s\n# Start\tEnd\tMatch\tID\tPvalue\tProbability\n"%spep.id)
        for elm in elmlist:
            matches = []
            score = 0
            for m in re.finditer(r"%s" %elm[0], str(spep.seq)):
                score += 1
                matches.append((m.start(),m.end(),m.group()))
            if score>0:
                pval = 0
                for i in range(100):
                    num = 0
                    seq = shuffle(str(spep.seq))
                    for m in re.finditer(r"%s" %elm[0], seq):
                        num += 1
                    pval += int(num>=score)
                pval = pval/100.0
                if pval<=0.1:
                    for m in matches:
                        outfile.write("\t\t%d\t%d\t%s\t%s\t%0.4f\t%s" \
                        %(m[0],m[1],m[2],elm[1],pval,elm[2]))
outfile.close()
