#!/usr/bin/python
# -*- coding: utf-8 -*-

import sqlite3, os, sys
from SQLite import SQLite

def wHTML(connexion, columns, title):
    global path, folder
    h = "<!DOCTYPE html>\n<html><head>\n"
    h+="""<link rel='stylesheet' href='style.css'>
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
        <script>$(document).ready(function() {
          $('tr').click( function() { $(this).next('.toggle').toggle() ; });
          $('.show').click(function(){ $('.toggle').show(); });
          $('.hide').click(function(){ $('.toggle').hide(); });
        });</script>\n"""
    h += """</head>\n<h1>{}</h1>\n<p>Click a row to show/hide the spep structure<br>
    <a href='#' class='show'>Show all</a><br><a href='#' class='hide'>
    Hide all</a></p><table>\n<tr>\n""".format(title)
    le = len(columns)
    for name in columns:
        h += "<th><b>{}</b></th>".format(name)
    try:
        speps=connexion.select('id, length, gene, location, target, signal, motifs')
    except:
        pass
    h += "</tr>\n"
    for row in speps:
        h += "<tr>\n"
        spep = row[0]
        for elt in row:
            h += "<td>{}</td>\n".format('<br>'.join(str(elt).split("\n")))
        try:
            htfile = open('{}/{}/cbs-anchor/{}'.format(path[0],folder,spep))
            h += "</tr>\n"
            h += "<tr class='toggle'><td colspan=" + str(le) + ">"
            for line in htfile:
                h += line
            h += "</td></tr>"
        except: pass
    h += "</table>\n<p><a href='#' class='show'>Show all</a><br><a href='#'' class='hide'>Hide all</a></p>\n</html>"
    return h

def wTSV(cursor,table,columns,title):
    h = "{}\n".format(title)
    for name in columns:
        h += "{}\t".format(name)
    h += "\n"
    try:
        cursor.execute('select * from {}'.format(table))
    except:
        pass
    for row in cursor.fetchall():
        for elt in row:
            h += "{}\t".format(str(elt))
        h += "\n"
    return h

if __name__=='__main__':
    # Set path and variables
    path = os.path.split(os.path.dirname(os.path.realpath(__file__)))
    folder = sys.argv[1]
    le = int(sys.argv[2])
    database = '%s/spepDB.db' %path[0]
    table = 'spep'
    dbconn = SQLite(database, table)
    columns = ('id', 'length', 'gene', 'location', 'target', 'signal', 'motifs')
    # Files
    sep_annot = open("{}/data/SEP_Annotated.csv".format(path[0]),'r')
    targf = open("{}/{}/TargetP.txt".format(path[0],folder),'r')
    signf = open("{}/{}/SignalP.txt".format(path[0],folder),'r')
    motf = open("{}/{}/Motifs.tsv".format(path[0],folder),'r')
    names = open("{}/{}/names".format(path[0],folder),'r')
    sep_annot.readline()
    dbconn.create('spep', '''(id varchar(30) primary key not null, length int, coordinates varchar(100),
    gene varchar(30), location text, exons varchar(50), aa text, cdna text)''')
    for line in sep_annot:
        line = line.split(',')
        tup = tuple([line[0],line[7],line[1],line[3],line[5],line[4],line[8],''])
        dbconn.insert(data=tup)
    # Update columns
    dbconn.addColumn("structure","varchar(30)")
    dbconn.addColumn("target","varchar(30)")
    dbconn.addColumn("signal","varchar(30)")
    dbconn.addColumn("motifs","text")
    for name in names:
        if not dbconn.selectWhere('*','id',name[:-1]):
            dbconn.insert(data="('{}')".format(name[:-1]),columns='(id)')
    for line in targf.readlines()[8:8+le]:
        line = line.split()
        if int(line[-1]) <= 3:
            dbconn.update("length",line[1],'id',line[0])
            if line[-2]=="M":
                dbconn.update("target","Mitochondrial",'id',line[0])
            elif line[-2]=="S":
                dbconn.update("target","Secretory",'id',line[0])
            else:
                dbconn.update("target",'Other','id',line[0])
        else:
            dbconn.update("target","&#8211",'id',line[0])
    for line in signf.readlines()[2:]:
        line = line.split()
        if line[9]=="N":
            dbconn.update("signal","No",'id',line[0])
        else:
            dbconn.update("signal","Yes",'id',line[0])
    # for line in motf.readlines()[1:]:
    #     line = line.split("\t")
    #     name = line[0].replace(':','_')
    #     if dbconn.selectWhere("motifs","id",name):
    #         d = dbconn.selectWhere("motifs","id",name)[0][0] + line[2] + "-" + line[3] + ": " + line[5] + "\n"
    #     else:
    #         d = line[2] + "-" + line[3] + ": " + line[5] + "\n"
    #     dbconn.update("motifs",d,'id',name)
    # HTML output
    print wHTML(dbconn, columns, 'sPEPs table')
    dbconn.close()
