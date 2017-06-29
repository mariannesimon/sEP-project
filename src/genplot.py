#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os

# Files
path = os.path.split(os.path.dirname(os.path.realpath(__file__)))
spep = sys.argv[1]
anchor_file = open('{}/{}/anchor_tmp/{}'.format(path[0],sys.argv[2],spep),'r')
surfp_file = open("{}/{}/NetSurfP.txt".format(path[0],sys.argv[2]),"r")
tmh_file = open("{}/{}/TMHMM.txt".format(path[0],sys.argv[2]),"r")
drnap_file = open("{}/{}/DRNApred.txt".format(path[0],sys.argv[2]),"r")
try:
    jpred = open("{}/{}/jpred_tmp/{}".format(path[0],sys.argv[2],spep),'r')
except:
    jpred = None

# Color palettes
# pal2 = [u'#ff0000', u'#ff5a00', u'#ff8100', u'#ff8f03', u'#ffb608']
# pal = [u'#73CCB2', u'#4C9999', u'#2E7085', u'#175275',u'#003366']
# pal3 = [u'#FF3300',u'#FF5180',u'#FF7000',u'#FF9600',u'#D6CC0A',u'#85CC1F',u'#70CC24', u'#5CCC29',u'#47CC2E',u'#33CC33']
# pal4 = [u'#FF0000',u'#FF1400',u'#FF2900',u'#FF3D00',u'#FF5180',u'#FF6600',u'#FF7A00', u'#FF8F00',u'#FFAD00',u'#FFCC00']
rb = [u'#EB5757',u'#DF5353',u'#D44E4E',u'#BC4646',u'#A43D3D',u'#8D3434',u'#762C2C',u'#5E2323',u'#3B1616',u'#170909' ]
rb2 = [u'#57D9D4',u'#4CBFD9',u'#3D99E0',u'#2966EB',u'#1433F5']
bb=[u'#88C0CA',u'#7CBAC5',u'#70B4BF',u'#64AEBA',u' #4CA1AF',u'#44919E',u'#397983',u'#2E6169',u'#1B383D',u'#000000']

# Variables
seq, struct, acc, sp, rsa = '', '', '', [], []

# NetSurfP
surflines = surfp_file.readlines()[15:]
for line in surflines:
    line = line.split()
    if len(line)==10 and line[2]==spep:
        seq += line[1]
        acc += line[0]*(line[0]=='B') + '-'*(line[0]=='E')
        rsa.append(int(10*float(line[4])))
        if float(line[7])>=0.5:
            struct+="H"
            sp.append(int(10*float(line[7]))-5)
        elif float(line[8])>=0.5:
            struct+="E"
            sp.append(int(10*float(line[8]))-5)
        elif float(line[9])>=0.5:
            struct+="-"
            sp.append(int(10*float(line[9]))-5)
        else:
            struct+="-"
            sp.append(0)
# Jpred
jnet,sol,conf='','',[]
if jpred:
    res = jpred.readlines()
    jnet = ''.join(res[4].split(':')[1].split(','))
    conf = res[5].split(':')[1].split(',')
    sol = ''.join(res[6].split(':')[1].split(','))
    hmm = ''.join(res[9].split(':')[1].split(','))

# Anchor
bind,disor,prob='',[],[]
for line in anchor_file:
    if line[0]!="#":
        bind += line.split()[1]*(line.split()[3]=='1') + '-'*(line.split()[3]=='0')
        prob.append(int(10*float(line.split()[2])))
        disor.append(int(10*float(line.split()[4])))
        # disor.append(int(10*float(line.split()[4])))

# TMHMM
topo=''
for line in tmh_file:
    if line.split()[0]==spep:
        topo = line.split()[-1].split('=')[-1]
        break

# DRNApred
drna=[]
rna,dna,pdna,prna='','',[],[]
drnalines = drnap_file.readlines()
for i in xrange(len(drnalines)):
    if drnalines[i][0]==">":
        if drnalines[i].split(">")[-1]=="{}\n".format(spep):
            drna = drnalines[i+2:i+2+len(struct)]
            break
for line in drna:
    dna += (line.split()[2]=='1')*line.split()[0] + (line.split()[2]=='0')*'-'
    pdna.append(int(float(line.split()[1])*10))
    prna.append(int(float(line.split()[3])*10))
    rna += (line.split()[4]=='1')*line.split()[0] + (line.split()[4]=='0')*'-'

anchor_file.close()
surfp_file.close()
tmh_file.close()
drnap_file.close()
try:
    jpred.close()
except:
    pass

# Print row
row,bar='',''
i = 1
while i<len(seq):
    row+=str(i)
    bar+='&#124;'
    nchar=8
    row+=' '*min(nchar,len(seq)-i-len(str(i))+1)
    bar+=' '*min(nchar+len(str(i))-1,len(seq)-i-len(str(i))+1)
    i += len(str(i))+nchar

# HTML output
cjnet, cs, crsa, cdis, crna, cdna, cbind = '', '', '', '', '', '',''
for i in range(len(struct)):
    crsa+='<font color={}>{}</font>'.format(rb[rsa[i]],acc[i])
    cs+='<font color={}>{}</font>'.format(bb[sp[i]+5],struct[i])
    cdis+='<font color={}>{}</font>'.format(bb[disor[i]],seq[i])
    if len(prna)>0:
        crna+='<font color={}>{}</font>'.format(bb[prna[i]],rna[i])
    if len(pdna)>0:
        cdna+='<font color={}>{}</font>'.format(bb[pdna[i]],dna[i])
    cbind+='<font color={}>{}</font>'.format(bb[prob[i]],bind[i])
    if len(conf)>0:
        cjnet+='<font color={}>{}</font>'.format(rb[int(conf[i])],jnet[i])
if sol=='':
    sol='\n'

h = "<pre>"
h+="""{:22}{}\n{:22}{}\n{:22}{}\n\n{:22}{}\n{:22}{}\n\n{:22}{}\n{:22}{}\n{:22}{}\n{:22}{}\n\n{:22}{}{:22}{}
""".format("",row,"",bar,"Sequence : ",seq,"NetsurfP SS : ",cs,"Jnet SS : ",cjnet, "Binding regions : ",\
cbind, "DNA binding : ",cdna,"RNA binding : ",crna,"Disordered regions : ",cdis, "25% Solvant Access : ",sol,"25% Exposure : ",acc)
h += """\n<b>Colour codes</b> : <font color=#88C0CA>low probability</font> to <font color=#000000>high probability</font>
               <font color=#EB5757>low confidence</font> to <font color=#000000>high confidence</font>"""
h += "</p></pre>"
print h
