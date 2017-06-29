#!/usr/bin/python
# -*- coding: utf-8 -*-
import sqlite3, sys

class SQLite:
    def __init__(self, dbname, table=None):
        try:
            self.conn = sqlite3.connect(dbname)
            self.cursor = self.conn.cursor()
            if table:
                self.table = table
        except:
            print "Error when connecting to database"
            sys.exit()

    def create(self, tablename, values):
        try:
            self.cursor.execute("drop table {}".format(tablename))
            self.cursor.execute("create table {} {}".format(tablename,values))
            self.table = tablename
        except:
            pass

    def close(self):
        self.conn.commit()
        self.conn.close()

    def columns(self):
        try:
            self.cursor.execute('pragma table_info(%s)' %self.table)
            return [tup[1] for tup in self.cursor.fetchall()]
        except:
            print "ERROR: Problem when fetching column names"

    def insert(self, data, columns=()):
        try:
            if len(columns)>0:
                self.cursor.execute('insert into {} {} values {}'.format(self.table,columns,data))
            else:
                self.cursor.execute('insert into {} values {}'.format(self.table,data))
        except:
            print "ERROR : Could not insert {} into table".format(data)
            sys.exit()

    def update(self, col, data, cond, row):
        try:
            self.cursor.execute("""update {} set {}=? where
                {}=?""".format(self.table,col,cond,),(data,row))
        except:
            pass

    def addColumn(self, col, coltype):
        try:
            self.cursor.execute('alter table {} add column {} {}'.format(self.table,col,coltype))
        except:
            pass
            # print "WARNING: Could not insert column {}".format(col)

    def select(self, colnames=''):
        try:
            if colnames=='' or colnames=='*':
                self.cursor.execute('select * from {}'.format(self.table))
            else:
                self.cursor.execute("select {} from {}".format(colnames, self.table))
            return self.cursor.fetchall()
        except:
            print "WARNING: Could not fetch data"

    def selectWhere(self, colnames, column, value):
        try:
            if colnames=='' or colnames=='*':
                self.cursor.execute("select * from {} where {}='{}'".format(self.table, column, value))
            else:
                self.cursor.execute("select {} from {} where {}='{}'".format(colnames, self.table, column, value))
            return self.cursor.fetchall()
        except:
            print "WARNING : Could not fetch {}={}".format(column,value)
