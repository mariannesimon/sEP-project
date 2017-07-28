#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os

# Open files
names_file = open('{}/names'.format(sys.argv[1]),'r')
anchor_file = open('{}/Anchor.txt'.format(sys.argv[1]),'r')
surfp_file = open("{}/NetSurfP.txt".format(sys.argv[1]),"r")
tmh_file = open("{}/TMHMM.txt".format(sys.argv[1]),"r")
drnap_file = open("{}/DRNApred.txt".format(sys.argv[1]),"r")

# Read files
surfp_lines = surfp_file.readlines()[15:]
anchor_lines = anchor_file.readlines()
tmh_lines = tmh_file.readlines()
drna_lines = drnap_file.readlines()
names = [line.split()[0] for line in names_file.readlines()]

# Color palettes
rb = [u'#EB5757',u'#DF5353',u'#D44E4E',u'#BC4646',u'#A43D3D',u'#8D3434',u'#762C2C',u'#5E2323',u'#3B1616',u'#170909' ]
bb=[u'#A5E0EB',u'#88C0CA',u'#7CBAC5',u'#70B4BF',u'#64AEBA',u' #4CA1AF',u'#44919E',u'#397983',u'#2E6169',u'#1B383D',u'#000000']

# Write summary for every name
for name in names:
    # Open and read Jpred file
    jnet,sol,hmm,pssm,conf='','','','',[]
    try:
        jpred = open("{}/jpred_tmp/{}".format(sys.argv[1],name),'r')
        res = jpred.readlines()
        jnet = ''.join(res[4].split(':')[1].split(','))
        conf = res[5].split(':')[1].split(',')
        sol = ''.join(res[6].split(':')[1].split(','))
        hmm = ''.join(res[9].split(':')[1].split(','))
        if res[10].split(':')[0]=='JNETPSSM':
            pssm = ''.join(res[10].split(':')[1].split(','))
        jpred.close()
    except:
        jpred = None
        pass
    # NetsurfP output
    seq, struct, acc, sp, rsa = '', '', '', [], []
    for line in surfp_lines:
        line = line.split()
        if len(line)==10 and line[2]==name:
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
    # Anchor
    bind,disor,prob='',[],[]
    for line in anchor_lines:
        if ''.join(line.split()[-1::])==name:
            bind += line.split()[1]*(line.split()[3]=='1') + '-'*(line.split()[3]=='0')
            prob.append(int(10*float(line.split()[2])))
            disor.append(int(10*float(line.split()[4])))
    # TMHMM
    topo=''
    for line in tmh_lines:
        if line.split()[0]==name:
            topo = line.split()[-1].split('=')[-1]
            break

    # DRNApred
    rna, dna, pdna, prna, drna = '', '', [], [], []
    for i in xrange(len(drna_lines)):
        if drna_lines[i][0]==">":
            if drna_lines[i].split(">")[-1]=="{}\n".format(name):
                drna = drna_lines[i+2:i+2+len(struct)]
                break
    for line in drna:
        dna += (line.split()[2]=='1')*line.split()[0] + (line.split()[2]=='0')*'-'
        pdna.append(int(float(line.split()[1])*10))
        prna.append(int(float(line.split()[3])*10))
        rna += (line.split()[4]=='1')*line.split()[0] + (line.split()[4]=='0')*'-'

    # Print format
    row,bar='',''
    i = 1
    while i<len(seq):
        row+=str(i)
        bar+='&#124;'
        nchar=8
        row+=' '*min(nchar,len(seq)-i-len(str(i))+1)
        bar+=' '*min(nchar+len(str(i))-1,len(seq)-i-len(str(i))+1)
        i += len(str(i))+nchar

    # Colour codes
    cjnet, cs, crsa, cdis, crna, cdna, cbind = '', '', '', '', '', '',''
    for i in range(len(struct)):
        crsa+='<font color={}>{}</font>'.format(rb[rsa[i]],acc[i])
        cs+='<font color={}>{}</font>'.format(bb[sp[i]+5],struct[i])
        cdis+='<font color={}>{}</font>'.format(bb[disor[i]],seq[i])
        cbind+='<font color={}>{}</font>'.format(bb[prob[i]],bind[i])
        if len(prna)>0:
            crna+='<font color={}>{}</font>'.format(bb[prna[i]],rna[i])
        if len(pdna)>0:
            cdna+='<font color={}>{}</font>'.format(bb[pdna[i]],dna[i])
        if len(conf)>0:
            cjnet+='<font color={}>{}</font>'.format(rb[int(conf[i])],jnet[i])
    if sol=='':
        sol='\n'

    c1, c2 = ['0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1'], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    cc1, cc2 = '', ''
    for i in range(len(c1)):
        cc1 += '<font color={}>{}</font> '.format(bb[i],c1[i])
    for i in range(len(rb)):
        cc2 += '<font color={}>{}</font> '.format(rb[i],c2[i])

    # HTML output
    out = open('{}/summary/{}'.format(sys.argv[1],name), 'w')
    h = "<pre>"
    h+="""{:30}{}\n{:30}{}\n{:30}{}\n\n{:30}{}\n{:30}{}\n{:30}{}\n\n{:30}{}\n{:30}{}\n{:30}{}\n{:30}{}\n\n{:30}{}{:30}{}
    """.format("",row,"",bar,"Sequence : ",seq,"NetsurfP SS : ",cs,"Jpred SS pred : ",cjnet, "Jpred SS from alignment : ", pssm, "Binding regions : ",\
    cbind, "DNA binding : ",cdna,"RNA binding : ",crna,"Disordered regions : ",cdis, "Jpred 25% Solvant Access : ",sol,"NetSurfP 25% Exposure : ",acc)
    h += """\n<b>Colour codes</b>  Probability : {}\n{:14}Confidence : {}""".format(cc1,'',cc2)
    h += "</pre>"
    out.writelines(h)

anchor_file.close()
surfp_file.close()
tmh_file.close()
drnap_file.close()
