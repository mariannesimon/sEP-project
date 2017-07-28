#!/usr/bin/python
import os.path, sys, re, random
from argparse import ArgumentParser
from Bio import SeqIO

# Set arguments
parser = ArgumentParser(description='ELM motif search')
parser.add_argument('--motifs', '-m', required=True,
                    help='file containing the ELM motif classes')
parser.add_argument('-i',help='input file in fasta format',required=True,metavar='INPUT_NAME')
parser.add_argument('-o',help='output file',required=True,metavar='OUTPUT_NAME')
parser.add_argument('-N',help='number of iterations for p-value estimation (default=500)',
                    default=500, type=int, metavar='ITERATIONS')
args = parser.parse_args()

# Set variables
elm_fname, seq_fname, out_name = args.motifs, args.i, args.o
N = args.N

# Shuffle sequence
def shuffle(seq):
    return ''.join(random.sample(seq,len(seq)))

# Motif search
with open(elm_fname,"r") as elm_file, open(out_name,"w") as outfile:
    raw_elm = elm_file.readlines()[6:]
    elmlist = [[line.split("\t")[i][1:-1] for i in range(8)] for line in raw_elm]
    outfile.write("sPEP id\tsPEP length\tStart\tEnd\tMatch\tELM_ID\tPvalue\n")
    # Read fasta file
    for spep in SeqIO.parse(seq_fname,"fasta"):
        # Motif search using regex
        matches = {(elm[1],elm[4]): [(m.start(),m.end(),m.group()) for m in re.finditer(r"%s" %elm[4], str(spep.seq)) ] \
                    for elm in elmlist if len(re.findall(r"%s" %elm[4], str(spep.seq)))>0}
        # Pvalue estimation
        pvalues = {elm[0]:0 for elm in matches.keys()}
        for i in xrange(N):
            pvalues = {elm[0]: pvalues[elm[0]]+int(len(re.findall(r"%s"%elm[1],shuffle(str(spep.seq))))>=len(matches[elm])) \
                        for elm in matches.keys()}
        pvalues = {elm[0]: pvalues[elm[0]]/float(N) for elm in matches.keys()}
        # print (sorted(pvalues.values()).index(pvalues[matches.keys()[0][0]])+1)/float(len(pvalues.values()))*0.05
        lines = [''.join( ['{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(spep.id, len(str(spep.seq)), m[0], m[1], m[2], elm[0], pvalues[elm[0]]) \
                    for m in matches[elm]] ) for elm in matches.keys() if pvalues[elm[0]]<0.01]
        outfile.writelines(lines)
