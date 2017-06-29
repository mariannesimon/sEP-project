#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from html import HTML
import os, re, sys, sqlite3

# Set path
path = os.path.split(os.path.dirname(os.path.realpath(__file__)))

# Slice the identifiers of Biomart file
make_id = lambda identifier: identifier.split("|")[0]

# Find regexp according to standard db translation codons
def findRE(index,fastaid,spep,frame,start,le):
    if fastaid!="None":
        seq = str(index[fastaid].seq[frame:].translate())
        if start=="Other":
            m = re.search(r"[A-Z]*%s[A-Z]*\*" %spep,seq)
        else:
            m = re.search(r"%s[A-Z]*%s[A-Z]*\*" %(str(Seq(start).translate()),spep),seq)
        if m:
            st = m.start()
            if m.group()[0]!="M":
                st = max(m.start(),max(0,m.end()-le-1))
            elif m.group()[0]=="M" and m.end()-m.start()>le:
                dist = m.end() - m.start() - le
                for p in [pos.start() for pos in re.finditer(r"M", m.group())]:
                    if abs(m.end()-m.start()-p-le) <= dist:
                        dist = abs(m.end()-m.start()-p-le)
                        st = m.start() + p
            return (seq[st:m.end()-1],str(index[fastaid].seq[i+3*st:i+3*m.end()]).upper())
    else: return None

# Write in fasta output files
def wFasta(spep,match):
    global cdna_file, spep_file
    cdna_file.write(">%s\n%s\n" %(spep[0],match[1]))
    spep_file.write(">%s\n%s\n" %(spep[0],match[0]))

# Write in database
def wdb(spep,transcript,match):
    global c
    c.execute('INSERT INTO spep VALUES (?,?,?,?,?,?,?)',\
    (spep[0],int(spep[-1]),spep[1],transcript,spep[3],match[0],match[1]))

# Output files
cdna_file = open("%s/Slavoff_cDNA.fasta" %path[0],"w")
spep_file = open("%s/Slavoff_sPEP.fasta" %path[0],"w")

# Database
conn = sqlite3.connect('%s/spepDB.db'%path[0])
c = conn.cursor()
try:
    c.execute('DROP TABLE spep')
except sqlite3.OperationalError:
    pass

c.execute('''CREATE TABLE spep (id varchar(30) primary key not null, length int,
gene varchar(30), transcript varchar(30), location text, aa text, cdna text)''')

with open("%s/data/sPEP_Slavoff_Annotated.tsv" %path[0],"r") as slavoff_db:
    colnames = slavoff_db.readline().split("\t")
    # Galaxy file
    galaxy = SeqIO.index("%s/data/Galaxy_DNA.fasta" %path[0],"fasta")
    # Biomart file
    biomart = SeqIO.index("%s/data/Biomart_cDNA.fasta" %path[0],"fasta",key_function=make_id)
    # Search in fasta files
    notFound = []
    for spep in slavoff_db:
        spep = spep.split("\t")
        match = []
        trlist = []
        for transcript in spep[2].split(";"):
            for i in range(0,3):
                tmatch = findRE(biomart,transcript,spep[5],i,spep[6],int(spep[-1]))
                if tmatch:
                    match.append(tmatch)
                    trlist.append(transcript)
        if len(match)==0:
            for i in range(0,3):
                match = findRE(galaxy,spep[0],spep[5],i,spep[6],int(spep[-1]))
                if match:
                    wdb(spep,"None",match)
                    wFasta(spep,match)
                    break
        elif len(match)==1:
            wdb(spep,trlist[0],match[0])
            wFasta(spep,match[0])
        else:
            dist = int(spep[-1])
            m = None
            for s in match:
                if abs(len(s[0])-int(spep[-1]))<dist:
                    dist = abs(len(s[0])-int(spep[-1]))
                    m = s
            wdb(spep,trlist[match.index(m)],m)
            wFasta(spep,m)
        if not match or len(match)==0:
            notFound.append(spep[0])

conn.commit()
conn.close()
cdna_file.close()
spep_file.close()
