#!/usr/bin/python
# -*- coding: utf-8 -*-
import sqlite3 as sqlite, sys

class SQLite:
    def __init__(self, dbname, table=None):
        try:
            self.conn = sqlite.connect(dbname)
            self.conn.row_factory = sqlite.Row
            self.cursor = self.conn.cursor()
            if table:
                self.table = table
        except sqlite.Error, e :
            print e
            sys.exit(1)

    def create(self, tablename, dict):
        # Drops existing table if exists and creates a new one
        # tablename : string
        # dict : dictionary of column names (keys) and column attributes (values)
        #   ex of dict : {'id', 'int primary key not null'}
        values = ', '.join([' '.join([key, dict[key]]) for key in dict.keys()])
        try:
            self.cursor.execute("drop table {}".format(tablename))
        except sqlite.Error, e:
            pass
        finally:
            self.cursor.execute("create table {} ({})".format(tablename,values))
            self.table = tablename

    def close(self):
        # Saves changes and closes connection with database
        self.conn.commit()
        self.conn.close()

    def getColumns(self):
        # Get all column information of table
        try:
            self.cursor.execute('pragma table_info(%s)' %self.table)
            return [tup[1] for tup in self.cursor.fetchall()]
        except sqlite.Error, e:
            print e

    def insert(self, data, columns=()):
        # Inserts new row in table in specified columns
        # data : tuple of values
        # columns : tuple of column names
        try:
            if len(columns)>0:
                self.cursor.execute('insert into {} {} values {}'.format(self.table,columns,data))
            elif len(columns)==0 and len(data)==1:
                self.cursor.execute('insert into {} (id) values (?)'.format(self.table,columns),(data[0],))
            else:
                self.cursor.execute('insert into {} values {}'.format(self.table,data))
        except sqlite.Error, e:
            print e
            sys.exit(1)

    def update(self, col, data, cond, row):
        # Updates one column (col) with value (data) where column value (cond) is equal to specified value (row)
        # col, data, cond, row : strings
        try:
            self.cursor.execute("update {} set {}=? where {}=?".format(self.table,col,cond,),(data,row))
        except sqlite.Error, e:
            print e

    def updateMany(self, cols, data, cond, row):
        # Updates several columns
        # cols : tuple of column names
        # data, cond, row : strings
        st = ', '.join(['='.join(s) for s in zip(cols,data)])
        try:
            self.cursor.execute("update {} set {} where {}=?".format(self.table,st,cond),(row,))
        except sqlite.Error, e:
            print e

    def addColumn(self, col, coltype):
        # col : column name
        # coltype : string with column attributes (ex : 'int primary key not null')
        try:
            self.cursor.execute('alter table {} add column {} {}'.format(self.table,col,coltype))
        except:
            pass

    def select(self, colnames=''):
        # selects all rows for specified columns
        # colnames : string of column names separated by a comma
        try:
            if colnames=='' or colnames=='*':
                self.cursor.execute('select * from {}'.format(self.table))
            else:
                self.cursor.execute("select {} from {}".format(colnames, self.table))
            return self.cursor.fetchall()
        except sqlite.Error, e:
            print e

    def selectWhere(self, colnames, column, value):
        # selects rows which meet the specified condition for specified columns
        # colnames : string of column names separated by comma
        # column, value : strings of values
        try:
            if colnames=='' or colnames=='*':
                self.cursor.execute("select * from {} where {}=?".format(self.table, column),(value,))
            else:
                self.cursor.execute("select {} from {} where {}=?".format(colnames, self.table, column),(value,))
            return self.cursor.fetchall()
        except sqlite.Error, e:
            print e
