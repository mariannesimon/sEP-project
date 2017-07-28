#!/usr/bin/python
# -*- coding: utf-8 -*-
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData
from Bio import SeqIO
import sys

# Input file
folder, seq_file = sys.argv[1], sys.argv[2]
if not seq_file:
    sys.exit(1)

# Output file
out_file = open('{}/protParams.tsv'.format(folder),'w')
out_file.write('id\tlen\tmw\tpI\tinstability')
for letter in IUPACData.protein_letters_1to3.values():
    out_file.write('\t{}'.format(letter))
out_file.write('\n')

# Read input and get params
for record in SeqIO.parse(seq_file,'fasta'):
    pA = ProteinAnalysis(str(record.seq))
    out_file.write('{}\t{}\t{:0.1f}\t{:0.2f}\t{}'.format(record.id.replace(':','_'),\
        len(str(record.seq)), pA.molecular_weight(),pA.isoelectric_point(), pA.instability_index()))
    for val in pA.count_amino_acids().values():
        out_file.write('\t{}'.format(val))
    out_file.write('\n')

out_file.close()
