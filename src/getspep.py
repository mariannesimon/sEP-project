#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
import os, re

# Set path
path = os.path.split(os.path.dirname(os.path.realpath(__file__)))
# Slice the identifiers of Biomart file
make_id = lambda identifier: identifier.split("|")[0]
# Write to stdout
out = lambda s, m: ">%s\n%s\n" %(s[0], m)
# Find regexp
def findRE(indexname,fastaID,spep,frame):
    """ Find Regexp in fasta file from sPEP short sequence.
        The codons are translated according to NCBI standart table,
        alternative start codons are not considered. """
    if fastaID!="None":
        le = int(spep[-1])
        seq = str(indexname[fastaID].seq[frame:].translate())
        if spep[6]=="Other":
            m = re.search(r"[A-Z]*%s[A-Z]*\*" %spep[5],seq)
        else:
            m = re.search(r"%s[A-Z]*%s[A-Z]*\*" \
            %(str(Seq(spep[6]).translate()),spep[5]),seq)
        if m:
            if m.group()[0]!="M":
                st = max(m.start(),max(0,m.end()-le-1))
            elif m.group()[0]=="M" and m.end()-m.start()>le:
                dist = m.end() - m.start() - le
                st = m.start()
                for p in [pos.start() for pos in re.finditer(r"M", m.group())]:
                    if abs(m.end()-m.start()-p-le) <= dist:
                        dist = abs(m.end()-m.start()-p-le)
                        st = m.start() + p
            else: st = m.start()
            cdna = str(indexname[fastaID].seq[i+3*st:i+3*m.end()]).upper()
            prot = seq[st:m.end()-1]
            return (cdna,prot)

# __main__
with open("%s/data/sPEP_Slavoff_Annotated.txt" %path[0],"r") as slavoff_table, \
    open("%s/Slavoff_cDNAs.fasta" %path[0],"w") as output_cdna, \
    open("%s/Slavoff_sPEPs.fasta" %path[0],"w") as output_spep:
    colnames = slavoff_table.readline().split("\t")
    # Galaxy file
    galaxy = SeqIO.index("%s/data/Galaxy_DNA.fasta" %path[0],"fasta")
    # Biomart file
    biomart = SeqIO.index("%s/data/Biomart_cDNA.fasta" %path[0],"fasta",key_function=make_id)
    # Search in fasta files
    notFound = []
    for spep in slavoff_table:
        spep = spep.split("\t")
        for i in range(0,3):
            match = findRE(galaxy,spep[0],spep,i)
            if match:
                output_cdna.write(out(spep,match[0]))
                output_spep.write(out(spep,match[1]))
                break
        if not match:
            match = []
            for transcript in spep[2].split(";"):
                for i in range(0,3):
                    tmatch = findRE(biomart,transcript,spep,i)
                    if tmatch and tmatch not in match:
                        match.append(tmatch)
            if len(match)==1:
                output_cdna.write(out(spep,match[0][0]))
                output_spep.write(out(spep,match[0][1]))
            elif len(match)>1:
                for s in match:
                    if len(s[1])==int(spep[-1]):
                        output_cdna.write(out(spep,s[0]))
                        output_spep.write(out(spep,s[1]))
                        break
        if not match or len(match)==0:
            notFound.append(spep[0])
    print notFound, len(notFound)
