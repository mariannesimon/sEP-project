#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3, os, sys
from argparse import ArgumentParser
from SQLite import SQLite
from collections import OrderedDict
import Tkinter as Tk

def wHTML(title, data, colnames):
    h = "<!DOCTYPE html>\n<html><head>\n"
    h+="""<style>body{font-family:"Helvetica Neue",Helvetica,Arial;text-align:center;}
        table{overflow-x:auto;border-collapse:collapse;width:100%;}
        th,td{padding:8px 10px;text-align:left;max-width:250px;overflow-x:auto;border-bottom:1px solid #ddd;}
        th{background-color:#f5f5f5;}tr:hover{background-color:#f5f5f5;}.toggle{display:none;}
        .toggle:hover{background-color:white;}.toggle td{text-align:left;}a{text-decoration:none;}
        pre{white-space:pre;white-space:pre-wrap;word-wrap:break-word;}</style>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
        <script>$(document).ready(function() { $('tr').click( function() { $(this).next('.toggle').toggle() ; });
        $('.show').click(function(){ $('.toggle').show(); }); $('.hide').click(function(){ $('.toggle').hide(); }); });
        </script>\n"""
    h += """</head>\n<h1>{}</h1>\n
    <p>Click a row to show/hide the spep structure</p><p><a href='#' class='show'>Show all</a> /
    <a href='#' class='hide'> Hide all</a></p><table>\n<tr>\n""".format(title)
    le = len(colnames)
    for name in colnames:
        h += "<th><b>{}</b></th>".format(name)
    h += "</tr>\n"
    for row in data:
        h += "<tr>\n"
        for elt in row:
            if elt and elt!=row['Structure']:
                h += "<td>{}</td>".format('<br>'.join(str(elt).split("\n")))
            elif not elt:
                h += "<td></td>"
        try:
            htfile = open(row['Structure'],'r')
            h += "</tr>\n<tr class='toggle'><td colspan=" + str(le) + ">"
            for line in htfile:
                h += line
            h += "</td></tr>"
        except: pass
    h += "</table>\n<p><a href='#' class='show'>Show all</a><br><a href='#'' class='hide'>Hide all</a></p>\n</html>"
    return h

def wCSV(data, colnames):
    h = ','.join(colnames)
    h += "\n"
    for row in data:
        elts = [' '.join(str(elt).replace('&#8211','-').split("\n")) for elt in row if elt != row['Structure']]
        if len(elts)<len(row)-1:
            for i in range(len(row)-1-len(elts)):
                elts.append('')
        h += ','.join(elts)
        h += "\n"
    return h

def wFasta(data, type='aa'):
    h = ''
    for row in data:
        if type=='aa':
            h += '>{}\n{}\n'.format(row['Id'],row['AA'])
        else:
            h += '>{}\n{}\n'.format(row['Id'],row['cDNA'])
    return h


if __name__=='__main__':

    parser = ArgumentParser(description='Script to edit and access SQLite database')
    parser.add_argument('--directory', '-d',
                        help='''Folder containing the pipeline output files. If provided, the database
                        will be updated according to the files.''')
    parser.add_argument('--annotation','-a', help='''File containing annotation for sEPs.
                        Only used when creating a new database.''')
    parser.add_argument('--sorfile','-s', help='''File containing annotation for sORFs.
                        Only used when creating a new database.''')
    parser.add_argument('--outfile','-o', help='Output file name.')
    parser.add_argument('--format', '-f', help='Output format (tsv, html or fasta). Only used if output file provided.')
    parser.add_argument('-db', help='Database name')
    parser.add_argument('--table', '-t', help='Table name.')
    args = parser.parse_args()

    if args.directory:
        res = raw_input('Do you want to update all rows in database ? y/n : ')
    else:
        res='n'
    # Database connection
    dbconn = SQLite(args.db, args.table)

    # Create database
    if not dbconn.select() or len(dbconn.select())==0:
        columns = OrderedDict([('Id','varchar(30) primary key not null'), ('Organism','varchar(100)'),('Length', 'int'), ('Start_codon','varchar(10)'),\
            ('Coordinates','varchar(100)'),('Transcript','varchar(50)'),('Gene','varchar(30)'), ('Annotation','text'), \
            ('Exons','varchar(50)'), ('Origin','varchar(100)'),('AA','text'), ('cDNA','text'), ('pI','float'), ('Molecular_weight','float'),
            ('Instability','float'),('Target','varchar(30)'), ('Signal_peptide','varchar(30)'), \
            ('Topology','varchar(50)'), ('Motifs','text'), ('Structure','varchar(255)')])
        dbconn.create(args.table, columns)
        insert_columns = tuple(columns.keys()[:11])
        try:
            sep_file = open(args.annotation,'r')
            sep_file.readline()
            for line in sep_file:
                l = line.split(',')
                values = tuple([l[0], l[8], l[7], l[1], l[4], l[3], l[6], l[5], 'Mass Spec', l[9][:-2], ''])
                dbconn.insert(data=values,columns=insert_columns)
        except:
            print 'Could not open %s' %args.annotation
        try:
            sorf_file = open(args.sorfile,'r')
            sorf_file.readline()
            for line in sorf_file:
                l = line.split('\t')
                values = tuple([l[0].replace(':','_'), l[9], l[10], 'chr{}:{}-{}'.format(l[2],l[4],l[5]),l[25], '', l[16], '', 'Ribo Prof', l[14], l[13]])
                dbconn.insert(data=values,columns=insert_columns)
        except:
            print 'Could not open %s' %args.sorfile

    # Open files and edit database
    if res=='y' or res=='Y' or res=='yes':
        try:
            count=0
            names = open("{}/names".format(args.directory),'r')
            for name in names:
                name = name[:-1]
                count+=1
                if not dbconn.selectWhere('*','Id',name):
                    dbconn.insert(data=tuple([name,'{}/summary/{}'.format(args.directory,name)]), columns=tuple(['Id','Structure']))
                else:
                    dbconn.update('Structure', '{}/summary/{}'.format(args.directory,name), 'Id', name)
            names.close()
        except:
            pass
        try:
            pa_file = open("{}/protParams.tsv".format(args.directory),'r')
            pa_file.readline()
            for line in pa_file:
                line = line.split('\t')
                dbconn.updateMany(('pI','Molecular_weight','Instability'),(line[3],line[2],line[4]),'Id',line[0])
            pa_file.close()
        except:
            pass
        try:
            target_file = open("{}/TargetP.txt".format(args.directory),'r')
            for line in target_file.readlines()[8:8+count]:
                line = line.split()
                if int(line[-1]) <= 3:
                    dbconn.update("Length",int(line[1]),'id',line[0])
                    if line[-2]=="M":
                        dbconn.update("Target","Mitochondrial",'Id',line[0])
                    elif line[-2]=="S":
                        dbconn.update("Target","Secretory",'Id',line[0])
                    else:
                        dbconn.update("Target",'Other','Id',line[0])
                else:
                    dbconn.update("Target","&#8211",'Id',line[0])
            target_file.close()
        except:
            pass
        try:
            signal_file = open("{}/SignalP.txt".format(args.directory),'r')
            for line in signal_file.readlines()[2:]:
                line = line.split()
                if line[9]=="N":
                    dbconn.update("Signal_peptide","No",'Id',line[0].replace(':','_'))
                else:
                    dbconn.update("Signal_peptide","Yes",'Id',line[0].replace(':','_'))
            signal_file.close()
        except:
            pass
        try:
            tmh_file = open("{}/TMHMM.txt".format(args.directory),'r')
            for line in tmh_file:
                dbconn.update('Topology',line.split()[-1].split('=')[-1],'Id',line.split()[0].replace(':','_'))
            tmh_file.close()
        except:
            pass
        try:
            motif_file = open("{}/Motifs.tsv".format(args.directory),'r')
            name, motif = '', ''
            for line in motif_file.readlines()[1:]:
                line = line.split("\t")
                if line[0].replace(':','_')==name:
                    motif += line[2] + "-" + line[3] + ": " + line[5] + "\n"
                else:
                    dbconn.update("Motifs",motif,'Id',name)
                    name = line[0].replace(':','_')
                    motif = line[2] + "-" + line[3] + ": " + line[5] + "\n"
            motif_file.close()
        except:
            pass

    # Output
    if args.outfile and args.format:
        data = dbconn.select('''Id , Organism, Length, Start, Transcript, Gene, Annotation,
            Exons, Origin, pI, Molecular_weight, Target, Signal_peptide, Topology, Motifs, AA, Structure''')
        head = ['sPEP id','Organism','Length','Start codon','Transcript','Gene','Annotation',\
            'Number of exons','Origin','pI','Molecular weight','Subcellular loc.','Signal Peptide', 'Topology', 'Motifs', 'AA']
        out = open(args.outfile, 'w')
        if args.format.lower()=='csv':
            out.writelines(wCSV(data, head))
        elif args.format.lower()=='html':
            out.writelines(wHTML('sPEPs table', data, head))
        elif args.format.lower()=='fasta' or args.format.lower()=='fa':
            out.writelines(wFasta(data))
    dbconn.close()
