#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#  
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not, write to 
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_mysql
Description:  RJE Module for converting delimited text files to MySQL tables
Version:      1.2
Last Edit:    05/09/11
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to take a bunch of delimited text files and create appropriate build statements for making
    MySQL tables, checking the content if desired. The main MySQL class has the global parameters for the reading and
    handling of input files, which are converted into Table objects. Each Table object itself has a number of Class
    objects. Unless memsaver=T, rje_mysql will check to see if contents for each field are unique.

    The script processes files and generates builds statements in the following way:
    
    1.  The next file to be processed is displayed along with the default table name, its first line (split into
        fields), the delimiter used, and the number of fields. The list of headers (set by header=T/F) are displayed
        along with the second line of the file, also split into fields. The delimiter is a tab for *.tdt files, a
        comma for *.csv files, else is set by delimit=X.

    2.  There is then a user menu with the following options:
        < S >kip file = Skip this file and move on to the next one.
        < T >able Name/Descriptions = option to change the table name (default = file name without extension) and the
            description for the table that is printed in the build file's comment lines
        < H >eader Names/Descriptions = option to manually change the names and descriptions for field headers
        < A >utoname Headers = if the file does not meet the given header=T/F option, this will automatically rename
            all field headers using the first line if desired, or simply "field1" to "fieldn".
        Change < D >elimiter = change the delimiter selected for the file and return to step 1 above.
        < P >roceed = keep current settings and process file

    3.  The file is then read through to assign Types to each field:
        (a) Field contents are scanned to see if they contain pure numeric values and, if so, pure integer values.
            Otherwise, string contents are assumed.
        (b) The presence of Null (empty) values are recorded. Also, unless memsaver=T, the field contents is analysed
            to determine whether it is unique for each row.
        (c) The minimum and maximum lengths are recorded for each field. Minimum and maximum values are recorded for
            numeric values

    4.  Field types are assigned according to the field content characteristics:
        (a) Numeric fields are set as "Unsigned" if the minimum value is >= 0. This is used with the minimum and
            maximum values to assign the integer type if an integer or FLOAT if non-integer.
        (b) String fields are assigned a type based on the minimum and maximum lengths of the contents and whether
            these two values are different. This way a CHAR, VARCHAR or BLOB is assigned of the correct length.
        If i=1+, an option is given to change the read characteristics of each field.

    5.  Unique fields without null entries are identified as potential Primary Keys. If i=-1 then the first is used,
        else the user can choose. If there are no potential fields, are the user chooses none, then an
        auto-incrementing INT(10) field is added as a primary key.

    6.  Each field can be selected to be indexed. By default, any non-numeric fields which have a maximum length in
        the range set by indexlen=X,Y are indexed.

    7.  Based on all the information given, the build statement is generated. This consists of:
        - Comment lines first identify the name, description, file and number of lines.
        - DROP TABLE IF EXISTS command to clear existing data
        - CREATE TABLE using the Types, Primary Key and Indexed fields determined above
        - LOAD DATA statement, including the list of fields in the file if an auto-incremented key was added
        - DELETE FROM statement to remove header lines, allowing for data truncation upon loading
        Upon hitting <ENTER>, this will be written to the given output file (buildfile=FILE & combine=T/F)

Commandline:
    filelist=FILE(s): Input file or files. May be comma-separated (FILE1,FILE2) and include wildcards. [*.tdt,*.csv]
    subfolders=T/F  : Whether to look in subfolders [False]
    mysql=T/F       : Whether to assing data types and check data (else just report lengths of field contents) [True]
    checktypes=T/F  : Whether to check Data Types for given file [True]
    buildfile=FILE  : Output file for MySQL Build statements [mysql_build.txt]
    combine=T/F     : Whether to combine build statements in one file (True) or have separate file per table (False) [True]
    append=T/F      : Whether to append output files or generate new [True]
    header=T/F      : By default, the first line will be read as a header [True]
    memsaver=T/F    : Will not check for Unique fields if True. [False]
    indexlen=X,Y    : Will index all fields with non-numerics between X and Y letters by default [4,30]
    sqldump=FILE    : Read in an SQL dump file and convert into CSV [None]
    
Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy
import os
#import random
import re
import string
import sys
import time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
#########################################################################################################################
### History
# 0.0 - Initial Compilation.
# 1.0 - Initial working version with basic functionality complete
# 1.1 - Added checking of field names against reserved words
# 1.2 - Added reading of MySQL dump file and converting into delimited (CSV) files.
#########################################################################################################################
### Major Functionality to Add
# [ ] : Add default values to Fields (info['Default'])
# [ ] : Add more log printing
# [ ] : Major tidy and overhaul!
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    try:
        start_time = time.time()
        program = 'RJE_MYSQL'
        version = '1.2'
        last_edit = 'September 11'  
        description = 'MySQL Build Module'
        author = 'Dr Richard J. Edwards.'
        info = rje.Info(program,version,last_edit,description,author,start_time)
        return info
    except:
        print 'Problem making Info object.'
        raise
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if info == None:
            info = makeInfo()
        if out == None:
            out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print "\n\nHelp for " + info.program + " V:" + info.version + ": " + time.asctime(time.localtime(info.start_time)) + "\n"
            out.verbose(-1,-1,text=__doc__)
            sys.exit()
        elif out.stat['Interactive'] > 1:    # Ask for more commands
            cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit:
        sys.exit()
    except KeyboardInterrupt:
        sys.exit()
    except:
        print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    ### <0> ### Objects setup
    try:
        ## <a> ## Initial Command Setup & Info
        cmd_list = sys.argv[1:]
        info = makeInfo()
        cmd_list = rje.getCmdList(cmd_list,info=info)      ### Load defaults from program.ini
        ## <b> ## Out object
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(0,2,cmd_list,1)
        out.printIntro(info)
        ## <c> ## Additional commands
        cmd_list = cmdHelp(info,out,cmd_list)
        ## <d> ## Log
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit:
        sys.exit()
    except:
        print "Problem during initial setup."
        raise
#########################################################################################################################
### MySQL reserved words
mysql_reserved_words = ['ACTION', 'ADD', 'AFTER', 'AGAINST', 'AGGREGATE', 'ALGORITHM', 'ALL', 'ALTER', 'ANALYZE',
                        'AND', 'ANY', 'AS', 'ASC', 'ASCII', 'ASENSITIVE', 'AUTO_INCREMENT', 'AVG', 'AVG_ROW_LENGTH',
                        'BACKUP', 'BDB', 'BEFORE', 'BEGIN', 'BERKELEYDB', 'BETWEEN', 'BIGINT', 'BINARY', 'BINLOG',
                        'BIT', 'BLOB', 'BOOL', 'BOOLEAN', 'BOTH', 'BTREE', 'BY', 'BYTE', 'CACHE', 'CALL', 'CASCADE',
                        'CASCADED', 'CASE', 'CHAIN', 'CHANGE', 'CHANGED', 'CHAR', 'CHARACTER', 'CHARSET', 'CHECK',
                        'CHECKSUM', 'CIPHER', 'CLIENT', 'CLOSE', 'COLLATE', 'COLLATION', 'COLUMN', 'COLUMNS', 'COMMENT',
                        'COMMIT', 'COMMITTED', 'COMPACT', 'COMPRESSED', 'CONCURRENT', 'CONDITION', 'CONNECTION',
                        'CONSISTENT', 'CONSTRAINT', 'CONTAINS', 'CONTINUE', 'CONVERT', 'CREATE', 'CROSS', 'CUBE',
                        'CURRENT_DATE', 'CURRENT_TIME', 'CURRENT_TIMESTAMP', 'CURRENT_USER', 'CURSOR', 'DATA',
                        'DATABASE', 'DATABASES', 'DATE', 'DATETIME', 'DAY', 'DAY_HOUR', 'DAY_MICROSECOND', 'DAY_MINUTE',
                        'DAY_SECOND', 'DEALLOCATE', 'DEC', 'DECIMAL', 'DECLARE', 'DEFAULT', 'DEFINER', 'DELAYED',
                        'DELAY_KEY_WRITE', 'DELETE', 'DESC', 'DESCRIBE', 'DES_KEY_FILE', 'DETERMINISTIC', 'DIRECTORY',
                        'DISABLE', 'DISCARD', 'DISTINCT', 'DISTINCTROW', 'DIV', 'DO', 'DOUBLE', 'DROP', 'DUAL',
                        'DUMPFILE', 'DUPLICATE', 'DYNAMIC', 'EACH', 'ELSE', 'ELSEIF', 'ENABLE', 'ENCLOSED', 'END',
                        'ENGINE', 'ENGINES', 'ENUM', 'ERRORS', 'ESCAPE', 'ESCAPED', 'EVENTS', 'EXECUTE', 'EXISTS',
                        'EXIT', 'EXPANSION', 'EXPLAIN', 'EXTENDED', 'FALSE', 'FAST', 'FETCH', 'FIELDS', 'FILE', 'FIRST',
                        'FIXED', 'FLOAT', 'FLOAT4', 'FLOAT8', 'FLUSH', 'FOR', 'FORCE', 'FOREIGN', 'FOUND', 'FRAC_SECOND',
                        'FROM', 'FULL', 'FULLTEXT', 'FUNCTION', 'GEOMETRY', 'GEOMETRYCOLLECTION', 'GET_FORMAT', 'GLOBAL',
                        'GOTO', 'GRANT', 'GRANTS', 'GROUP', 'HANDLER', 'HASH', 'HAVING', 'HELP', 'HIGH_PRIORITY',
                        'HOSTS', 'HOUR', 'HOUR_MICROSECOND', 'HOUR_MINUTE', 'HOUR_SECOND', 'IDENTIFIED', 'IF', 'IGNORE',
                        'IMPORT', 'IN', 'INDEX', 'INDEXES', 'INFILE', 'INNER', 'INNOBASE', 'INNODB', 'INOUT',
                        'INSENSITIVE', 'INSERT', 'INSERT_METHOD', 'INT', 'INT1', 'INT2', 'INT3', 'INT4', 'INT8',
                        'INTEGER', 'INTERVAL', 'INTO', 'INVOKER', 'IO_THREAD', 'IS', 'ISOLATION', 'ISSUER', 'ITERATE',
                        'JOIN', 'KEY', 'KEYS', 'KILL', 'LABEL', 'LANGUAGE', 'LAST', 'LEADING', 'LEAVE', 'LEAVES', 'LEFT',
                        'LEVEL', 'LIKE', 'LIMIT', 'LINES', 'LINESTRING', 'LOAD', 'LOCAL', 'LOCALTIME', 'LOCALTIMESTAMP',
                        'LOCK', 'LOCKS', 'LOGS', 'LONG', 'LONGBLOB', 'LONGTEXT', 'LOOP', 'LOW_PRIORITY', 'MASTER',
                        'MASTER_CONNECT_RETRY', 'MASTER_HOST', 'MASTER_LOG_FILE', 'MASTER_LOG_POS', 'MASTER_PASSWORD',
                        'MASTER_PORT', 'MASTER_SERVER_ID', 'MASTER_SSL', 'MASTER_SSL_CA', 'MASTER_SSL_CAPATH',
                        'MASTER_SSL_CERT', 'MASTER_SSL_CIPHER', 'MASTER_SSL_KEY', 'MASTER_USER', 'MATCH',
                        'MAX_CONNECTIONS_PER_HOUR', 'MAX_QUERIES_PER_HOUR', 'MAX_ROWS', 'MAX_UPDATES_PER_HOUR',
                        'MAX_USER_CONNECTIONS', 'MEDIUM', 'MEDIUMBLOB', 'MEDIUMINT', 'MEDIUMTEXT', 'MERGE',
                        'MICROSECOND', 'MIDDLEINT', 'MIGRATE', 'MINUTE', 'MINUTE_MICROSECOND', 'MINUTE_SECOND',
                        'MIN_ROWS', 'MOD', 'MODE', 'MODIFIES', 'MODIFY', 'MONTH', 'MULTILINESTRING', 'MULTIPOINT',
                        'MULTIPOLYGON', 'MUTEX', 'NAME', 'NAMES', 'NATIONAL', 'NATURAL', 'NCHAR', 'NDB', 'NDBCLUSTER',
                        'NEW', 'NEXT', 'NO', 'NONE', 'NOT', 'NO_WRITE_TO_BINLOG', 'NULL', 'NUMERIC', 'NVARCHAR',
                        'OFFSET', 'OLD_PASSWORD', 'ON', 'ONE', 'ONE_SHOT', 'OPEN', 'OPTIMIZE', 'OPTION', 'OPTIONALLY',
                        'OR', 'ORDER', 'OUT', 'OUTER', 'OUTFILE', 'PACK_KEYS', 'PARTIAL', 'PASSWORD', 'PHASE', 'POINT',
                        'POLYGON', 'PRECISION', 'PREPARE', 'PREV', 'PRIMARY', 'PRIVILEGES', 'PROCEDURE', 'PROCESSLIST',
                        'PURGE', 'QUARTER', 'QUERY', 'QUICK', 'RAID0', 'RAID_CHUNKS', 'RAID_CHUNKSIZE', 'RAID_TYPE',
                        'READ', 'READS', 'REAL', 'RECOVER', 'REDUNDANT', 'REFERENCES', 'REGEXP', 'RELAY_LOG_FILE',
                        'RELAY_LOG_POS', 'RELAY_THREAD', 'RELEASE', 'RELOAD', 'RENAME', 'REPAIR', 'REPEAT', 'REPEATABLE',
                        'REPLACE', 'REPLICATION', 'REQUIRE', 'RESET', 'RESTORE', 'RESTRICT', 'RESUME', 'RETURN',
                        'RETURNS', 'REVOKE', 'RIGHT', 'RLIKE', 'ROLLBACK', 'ROLLUP', 'ROUTINE', 'ROW', 'ROWS',
                        'ROW_FORMAT', 'RTREE', 'SAVEPOINT', 'SCHEMA', 'SCHEMAS', 'SECOND', 'SECOND_MICROSECOND',
                        'SECURITY', 'SELECT', 'SENSITIVE', 'SEPARATOR', 'SERIAL', 'SERIALIZABLE', 'SESSION', 'SET',
                        'SHARE', 'SHOW', 'SHUTDOWN', 'SIGNED', 'SIMPLE', 'SLAVE', 'SMALLINT', 'SNAPSHOT', 'SOME',
                        'SONAME', 'SOUNDS', 'SPATIAL', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE', 'SQLWARNING',
                        'SQL_BIG_RESULT', 'SQL_BUFFER_RESULT', 'SQL_CACHE', 'SQL_CALC_FOUND_ROWS', 'SQL_NO_CACHE',
                        'SQL_SMALL_RESULT', 'SQL_THREAD', 'SQL_TSI_DAY', 'SQL_TSI_FRAC_SECOND', 'SQL_TSI_HOUR',
                        'SQL_TSI_MINUTE', 'SQL_TSI_MONTH', 'SQL_TSI_QUARTER', 'SQL_TSI_SECOND', 'SQL_TSI_WEEK',
                        'SQL_TSI_YEAR', 'SSL', 'START', 'STARTING', 'STATUS', 'STOP', 'STORAGE', 'STRAIGHT_JOIN',
                        'STRING', 'STRIPED', 'SUBJECT', 'SUPER', 'SUSPEND', 'TABLE', 'TABLES', 'TABLESPACE', 'TEMPORARY',
                        'TEMPTABLE', 'TERMINATED', 'TEXT', 'THEN', 'TIME', 'TIMESTAMP', 'TIMESTAMPADD', 'TIMESTAMPDIFF',
                        'TINYBLOB', 'TINYINT', 'TINYTEXT', 'TO', 'TRAILING', 'TRANSACTION', 'TRIGGER', 'TRIGGERS',
                        'TRUE', 'TRUNCATE', 'TYPE', 'TYPES', 'UNCOMMITTED', 'UNDEFINED', 'UNDO', 'UNICODE', 'UNION',
                        'UNIQUE', 'UNKNOWN', 'UNLOCK', 'UNSIGNED', 'UNTIL', 'UPDATE', 'USAGE', 'USE', 'USER',
                        'USER_RESOURCES', 'USE_FRM', 'USING', 'UTC_DATE', 'UTC_TIME', 'UTC_TIMESTAMP', 'VALUE', 'VALUES',
                        'VARBINARY', 'VARCHAR', 'VARCHARACTER', 'VARIABLES', 'VARYING', 'VIEW', 'WARNINGS', 'WEEK',
                        'WHEN', 'WHERE', 'WHILE', 'WITH', 'WORK', 'WRITE', 'X509', 'XA', 'XOR', 'YEAR', 'YEAR_MONTH',
                        'ZEROFILL']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CLASSES                                                                                                 #
#########################################################################################################################

#########################################################################################################################
### Main MySQL Class:                                                                                                   #
#########################################################################################################################
class MySQL(rje.RJE_Object):     
    '''
    MySQL Class. Author: Rich Edwards (2005).

    Info:str
    - Name = 'Project' name - output build file name
    - FileList = Input file or files. May be comma-separated (FILE1,FILE2) and include wildcards. [*.tdt,*.csv]
    - SQLDump = Read in an SQL dump file and convert into CSV [None]
    
    Opt:boolean
    - SubFolders = Whether to look in subfolders [False]
    - MySQL = Whether to assing data types and check data (else just report lengths of field contents) [True]
    - CheckTypes = Whether to check Data Types for given file [True]
    - Combine = Whether to combine build statements in one file (True) or have separate file per table (False) [True]
    - Append = Whether to append output files or generate new [False]
    - Header = By default, the first line will be read as a header [True]
    - MemSaver = Will not check for Unique fields if False. [True]

    Stat:numeric
    - MinIndexLen,MaxIndexLen = Will index all fields with non-numerics between X and Y letters by default [4,10]

    Obj:RJE_Objects
    '''
    ### Attributes
    tablelist = []  ### List of tables
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['FileList']
        - Opt:boolean ['SubFolders','MySQL','CheckTypes','Combine','Header','MemSaver']
        - Stats:float ['MinIndexLen','MaxIndexLen']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['FileList','SQLDump']
        self.optlist = ['SubFolders','MySQL','CheckTypes','Combine','Header','MemSaver']
        self.statlist = ['MinIndexLen','MaxIndexLen']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None)
        self.info['Name'] = 'mysql_build.txt'
        self.info['FileList'] = '*.csv,*.tdt'
        self.opt['Append'] = True
        self.opt['SubFolders'] = False
        self.opt['MemSaver'] = False
        self.setStat({'MinIndexLen':4,'MaxIndexLen':30})
        ### <c> ### Other Attributes
        self.tablelist = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='Name',arg='buildfile')
                self._cmdRead(cmd,type='info',att='Name')
                self._cmdRead(cmd,type='info',att='FileList')
                self._cmdRead(cmd,type='file',att='SQLDump')
                self._cmdRead(cmd,type='opt',att='MySQL')
                self._cmdRead(cmd,type='opt',att='Combine')
                self._cmdRead(cmd,type='opt',att='Header')
                self._cmdRead(cmd,type='opt',att='MemSaver')
                self._cmdRead(cmd,type='min',att='MinIndexLen',arg='indexlen')
                self._cmdRead(cmd,type='max',att='MaxIndexLen',arg='indexlen')
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Class Methods                                                                                           #
#########################################################################################################################
    def run(self):      ### Overall Run Method
        '''
        Overall Run Method.
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            if rje.exists(self.info['SQLDump']): return self.sqlConvert()
            if self.opt['Combine'] and not self.opt['Append'] and os.path.exists(self.info['Name']) and self.stat['Interactive'] >= 0:
                if rje.yesNo('Output file %s already exists. Overwrite?' % self.info['Name']):
                    os.unlink(self.info['Name'])
                else:
                    return
            filelist = rje.getFileList(callobj=self,filelist=rje.split(self.info['FileList'],','),subfolders=self.opt['SubFolders'])
            self.verbose(1,1,'%d Putative files to convert to MySQL.' % len(filelist),2)

            ### <1> ### Work Through Each File
            _stage = '<1> Work Through Files'
            processed = 0
            for file in filelist:
                processed += self._processFile(file)
            self.verbose(1,1,'%d files processed successfully.' % processed,2)
                
        except:
            self.log.errorLog('Error in run(%s)' % _stage,printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _processFile(self,file,delimit=None):    ### Processes given File
        '''
        Processes given File.
        >> file:str = file name
        >> delimit:str = delimit to over-ride automatically determined delimiter
        << ok:int = 1 if processed, 0 if skipped/failed
        '''
        try:
            ### Get Table Name and Fields ###
            tablename = rje.baseFile(file,True)
            tabledesc = 'Generated from %s' % file
            if not delimit:
                delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(os.path.splitext(file)[1]))
            TABLE = open(file,'r')
            firstline = rje.readDelimit(TABLE.readline(),delimit)
            secondline = rje.readDelimit(TABLE.readline(),delimit)
            TABLE.close()
            self.verbose(0,3,'\n\nNext file: %s\nName: %s' % (file,tablename),1)
            self.verbose(0,3,'\nDelimit: "%s"\nFirst line: %s' % (rje.replace(delimit,'\\','\\\\'),firstline),1)
            self.verbose(0,2,'=> %d Fields.' % len(firstline),1)
            headers = firstline[0:]
            desclist = firstline[0:]
            usehead = self.opt['Header']
            if not usehead:
                for f in range(len(firstline)):
                    headers[f] = 'field%d' % (f+1)
            for f in range(len(firstline)):
                desclist[f] = '%s %s' % (os.path.basename(file),headers[f])
            self.verbose(0,2,'\n=> Headers: %s' % headers,1)
            if secondline:
                self.verbose(0,3,'Second line: %s' % secondline,1)

            ### Options Menu ###
            while self.stat['Interactive'] >= 0:     # Cycle Through Menu
                choice = rje.choice('\n<S>kip file, <T>able Name/Description, <H>eader Names/Descriptions, <A>utoname Headers, Change <D>elimiter, <P>roceed?',default='P').upper()
                if choice == 'S': # Skip this file
                    return 0
                elif choice == 'T': # Rename table
                    tablename = rje.choice('New Table Name?',tablename,confirm=True)
                    tabledesc = rje.choice('New Table Description?',tabledesc,confirm=True)
                elif choice == 'H': # Rename headers
                    for f in range(len(firstline)):
                        headers[f] = rje.choice('New header for field: %s' % firstline[f],headers[f],True)
                        desclist[f] = rje.choice('Description for field: %s' % headers[f],desclist[f],True)
                elif choice == 'A': # Autoname headers
                    if usehead:
                        usehead = rje.yesNo('First line is field headers?')
                    else:
                        usehead = rje.yesNo('First line is field headers?','N')
                    if usehead:
                        headers = firstline[0:]
                    else:
                        for f in range(len(firstline)):
                            headers[f] = 'field%d' % (f+1)
                    for f in range(len(firstline)):
                        desclist[f] = '%s %s' % (os.path.basename(file),headers[f])                            
                elif choice == 'D' and rje.yesNo('Change delimiter from "%s"?' % rje.replace(delimit,'\\','\\\\')):
                    delimit = rje.choice('New text delimiter?',delimit,True)
                    return self._processFile(file,delimit)
                elif choice == 'P':   # Proceed
                    forbidden_fields = []
                    for h in headers:
                        if h.upper() in mysql_reserved_words:
                            forbidden_fields.append(h)
                    if forbidden_fields and not rje.yesNo('\nWARNING! The following field names are reserved MySQL words and may lead to errors:\n - %s.\nProceed?' % forbidden_fields,default='N'):
                        continue
                    break

            ### Check and Log reserved words
            for h in headers:
                if h in mysql_reserved_words:
                    self.log.printLog('#FIELD','WARNING: %s is a MySQL reserved word and may cause an error as a field name.' % h)

            ### Create New Table Object ###
            table = Table(self.log,self.cmd_list)
            table.setInfo({'Name':tablename,'Description':tabledesc,'File':file})
            table.opt['MemSaver'] = self.opt['MemSaver']
            #print self.opt['MemSaver'],table.opt['Unique']
            if usehead:
                table.headers = firstline
            else:
                table.headers = ['']
            ## Fields ##
            for f in range(len(firstline)):
                table._addField(headers[f])
                table.fieldlist[f].info['Description'] = desclist[f]
            table._readFile()

            ### Make MySQL Field Types ###
            if not self.opt['MySQL']:
                table.summarise()
                return 1
            table._fieldTypes()
            table._makePrimaryKey()
            table._setIndices(self.stat['MinIndexLen'],self.stat['MaxIndexLen'])

            ### Check data integrity ###
            #if self.opt['CheckTypes']:
            #!# Add check data stuff    

            ### Add to Build File ###
            buildtext = table._buildText(delimit)
            self.verbose(0,0,'\n\n%s\n\n' % buildtext,2)
            if self.info['Name'] == 'None':     ### Make build file 'mysql_build.txt'
                return 1
            if self.opt['Combine']:
                buildfile = self.info['Name'] 
                BUILDFILE = open(buildfile,'a') # This was deleted at start if append=F
            else:
                buildfile = '%s.%s' % (table.info['Name'],self.info['Name'])
                BUILDFILE = open(buildfile,'w')
            BUILDFILE.write('%s\n\n\n\n' % buildtext)
            BUILDFILE.close()
            self.log.printLog('#TAB','%s build statement output to %s.' % (table.info['Name'],buildfile))

            ### Finish ###                
            return 1
        except:
            self.log.errorLog('Error in processFile(%s)' % file,printerror=True,quitchoice=True)
            return 0
#########################################################################################################################
    def sqlConvert(self):   ### Crude converter for SQL dump file into CSV files
        '''Crude converter for SQL dump file into CSV files'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sql = self.info['SQLDump']
            base = rje.baseFile(sql)
            DUMP = open(sql,'r')
            CSV = None
            table = ''
            line = DUMP.readline()
            create = False
            tx = 0
            ### ~ [1] ~ Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while line:
                if rje.matchExp('CREATE TABLE `(\S+)`',line):
                    create = True
                    table = rje.matchExp('CREATE TABLE `(\S+)`',line)[0]
                    csv = '%s.%s.csv' % (base,table)
                    self.printLog('#CSV','Creating %s: %s' % (table,csv))
                    CSV = open(csv,'w')
                    chead = []
                elif create and rje.matchExp('^\s+`(\S+)`',line): chead.append(rje.matchExp('^\s+`(\S+)`',line)[0])
                elif create and rje.matchExp('^\s+(PRIMARY)',line):
                    CSV.write('%s\n' % rje.join(chead,','))
                    self.printLog('#HEAD','%s: %s' % (table,rje.join(chead,',')))
                    create = False
                elif line.find('INSERT INTO `%s`' % table) == 0:
                    dx = 0
                    for entry in rje.split(line,'(')[1:]:
                        output = rje.readDelimit(rje.replace(entry[:entry.find(')')],"'",'"'),',')
                        rje.writeDelimit(CSV,output,',')
                        dx += 1
                elif line.find('UNLOCK TABLES') == 0:
                    self.printLog('#TAB','%s: %s entries' % (table,rje.iStr(dx)))
                    CSV.close(); tx += 1
                line = DUMP.readline()
            self.printLog('#SQL','%d tables read from %s' % (tx,sql))
        except: self.errorLog('sqlConvert stuffed')            
#########################################################################################################################
### End of MySQL Class                                                                                                  #
#########################################################################################################################
       
                                                    ### ~ ### ~ ###

#########################################################################################################################
### Table Class:                                                                                                        #
#########################################################################################################################
class Table(rje.RJE_Object):     
    '''
    MySQL Table Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Table name
    - File = Name of input file
    - Description = Table description
    
    Opt:boolean

    Stat:numeric
    - Lines = Number of lines read in
    - Screened = Number of header lines screened out
    - BadLines = Number of bad lines (wrong number of fields) screened out

    Obj:RJE_Objects
    - PrimaryKey = Field Object
    '''
    ### Attributes
    headers = []    # List of headers to use to screen lines when reading table
    fieldlist = []  # List of fields
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['File','Description']
        - Opt:boolean []
        - Stats:float ['Lines','Screened','BadLines']
        - Obj:RJE_Object ['PrimaryKey']
        '''
        ### <a> ### Basics 
        self.infolist = ['File','Description']
        self.optlist = []
        self.statlist = ['Lines','Screened','BadLines']
        self.objlist = ['PrimaryKey']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        ### <c> ### Other Attributes
        self.fieldlist = []
#########################################################################################################################
    ### <2> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _addField(self,fieldname):      ### Adds a field to the table
        '''
        Adds a field to the table.
        >> fieldname:str = Field Name
        '''
        try:
            self.fieldlist.append(Field(self.log,self.cmd_list))
            self.fieldlist[-1].info['Name'] = fieldname
            self.fieldlist[-1].opt['Unique'] = not self.opt['MemSaver']
        except:
            self.log.errorLog('Error in _addField(%s)' % fieldname,printerror=True,quitchoice=True)
#########################################################################################################################
    def _readFile(self):      ### Reads data from file and sets field attributes accordingly
        '''
        Reads data from file and sets field attributes accordingly.
        '''
        try:
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(os.path.splitext(self.info['File'])[1]))
            #self.verbose(0,3,'\nReading %s data from %s...' % (self.info['Name'],self.info['File']),0)
            self.verbose(0,4,'',1)
            self.log.printLog('#TAB','Reading %s data from %s...' % (self.info['Name'],self.info['File']))
            TABLE = open(self.info['File'],'r')
            line = TABLE.readline()
            rx = 0
            while line not in [None,'']:
                nextline = rje.readDelimit(line,delimit)
                if nextline == ['']:
                    line = TABLE.readline()
                    continue
                if '%s' % nextline == '%s' % self.headers:  # Skip line
                    self.verbose(0,3,'Skipping header line...',0)
                    line = TABLE.readline()
                    self.stat['Screened'] += 1
                    continue
                if len(nextline) != len(self.fieldlist):
                    self.log.errorLog('Only %d fields read. Should be %d: %s' % (len(nextline),len(self.fieldlist),nextline),False,False)
                    line = TABLE.readline()
                    self.stat['BadLines'] += 1
                    continue
                for f in range(len(nextline)):
                    self.fieldlist[f]._readData(nextline[f])
                rx += 1
                rje.progressPrint(self,rx)
                line = TABLE.readline()
            self.stat['Lines'] = rx + self.stat['Screened'] + self.stat['BadLines']
            self.verbose(0,4,'',1)
            self.log.printLog('#TAB','%s lines of %d fields successfully read from %s.' % (rje.integerString(rx),len(self.fieldlist),self.info['File']))
            #self.verbose(0,1,'%s lines of %d fields successfully read from %s.' % (rje.integerString(rx),len(self.fieldlist),self.info['File']),2)
        except:
            self.log.errorLog('Error in _readFile',printerror=True,quitchoice=True)
#########################################################################################################################
    def summarise(self):  ### Print details to screen
        '''
        Print details to screen.
        '''
        self.verbose(0,3,self.details(),0)
        for field in self.fieldlist:
            self.verbose(0,3,'*** Field %d ***' % (self.fieldlist.index(field)+1),1)
            self.verbose(0,2,field.details(),0)
        self.verbose(0,1,'%d Fields in total.' % len(self.fieldlist),1)
#########################################################################################################################
    def _fieldTypes(self):  ### Determines field types from attributes
        '''
        Determines field types from attributes.
        '''
        try:
            self.verbose(0,3,self.details(),0)
            self.verbose(0,1,'%d Fields in total.' % len(self.fieldlist),2)
            for field in self.fieldlist:
                self.verbose(0,3,'*** Field %d ***' % (self.fieldlist.index(field)+1),1)
                field._autoType()
                self.verbose(0,2,field.details(),0)
                while self.stat['Interactive'] > 0 and not rje.yesNo('Keep %s details?' % field.info['Name']):
                    field.edit()
        except:
            self.log.errorLog('Error in _readFile',printerror=True,quitchoice=True)
#########################################################################################################################
    def _makePrimaryKey(self):   ### Sets Primary Key for table
        '''
        Sets Primary Key for table.
        '''
        try:
            possible_keys = []
            field_dic = {}
            for field in self.fieldlist:
                if field.opt['Unique'] and field.opt['NotNull']:
                    possible_keys.append(field.info['Name'])
                    field_dic[field.info['Name']] = field
            self.verbose(0,2,'\n\n%d possible Primary Key fields: %s' % (len(possible_keys),possible_keys),1)
            if possible_keys == []:
                self._addAutoIncrementKey()
                return
            #elif len(possible_keys) == 1 and self.stat['Interactive'] < 1:
            #    self.setPrimaryKey(field_dic[possible_keys[0]])
            #    return
            for option in possible_keys:
                if self.stat['Interactive'] < 0 or rje.yesNo('Make %s Primary Key?' % option):
                    self.setPrimaryKey(field_dic[option])
                    return
            if rje.yesNo('Add auto-incrementing ID as Primary Key?'):
                self._addAutoIncrementKey()
                return
            self._makePrimaryKey()
        except:
            self.log.errorLog('Error in _readFile',printerror=True,quitchoice=True)
#########################################################################################################################
    def setPrimaryKey(self,keyfield):  ### Sets keyfield as primary key
        '''
        Sets keyfield as primary key.
        >> keyfield:Field object = Primary Key
        '''
        for field in self.fieldlist:
            field.opt['PrimaryKey'] = False
        keyfield.opt['PrimaryKey'] = True
        self.log.printLog('#KEY','Primary Key for Table %s = %s.' % (self.info['Name'],keyfield.info['Name']))
#########################################################################################################################
    def _addAutoIncrementKey(self): ### Adds auto incrementing field
        '''
        Adds auto incrementing field.
        '''
        self.log.printLog('#KEY','Adding auto-incrementing ID field for Primary Key for Table %s.' % self.info['Name'])
        newfield = Field(self.log,self.cmd_list)
        newfield.info['Name'] = '%s_id' % self.info['Name'].lower()
        if self.stat['Interactive'] >= 0:
            if rje.yesNo('Name auto-incrementing ID field "auto_id"? (<N>o for custom name.)'):
                newfield.info['Name'] = 'auto_id'
            else:
                newfield.info['Name'] = rje.choice('Name for auto-incrementing ID field?',default=newfield.info['Name'],confirm=True)
        newfield.info['Type'] = 'INT(10)'
        newfield.info['Description'] = 'Auto-incrementing Primary Key ID'
        newfield.setOpt({'Unsigned':True,'NotNull':True,'AutoIncrement':True,'Numeric':True})
        self.fieldlist = [newfield] + self.fieldlist
        self.setPrimaryKey(newfield)
#########################################################################################################################
    def _setIndices(self,min,max):  ### Sets index on fields according to min and max lengths
        '''
        Sets index on fields according to min and max lengths.
        >> min & max:int = inclusive values between which the max length of a text field must fall to be indexed.
        '''
        if self.stat['Interactive'] >= 0 and not rje.yesNo('Add additional indexed fields?'):
            return
        for field in self.fieldlist:
            if field.opt['PrimaryKey']:
                continue
            field.opt['Index'] = not field.opt['Numeric']
            if field.stat['MaxLen'] < min or field.stat['MaxLen'] > max:
                field.opt['Index'] = False
            if self.stat['Interactive'] >= 0 and field.opt['Index']:
                field.opt['Index'] = rje.yesNo('Index %s?' % field.info['Name'])
            elif self.stat['Interactive'] >= 0:
                field.opt['Index'] = rje.yesNo('Index %s?' % field.info['Name'],default='N')
#########################################################################################################################
    def _buildText(self,delimit='\t'):   ### Returns the auto-generated build statement for making a MySQL table
        '''
        Returns the auto-generated build statement for making a MySQL table.
        << buildtext:str = multiline build statement
        '''
        ### Setup ###
        desc_add = 35
        fieldlen = 0
        indexlist = []
        fieldlist = []
        for field in self.fieldlist:
            fieldlist.append(field.info['Name'])
            if len(field.info['Name']) > fieldlen:
                fieldlen = len(field.info['Name'])
            if field.opt['Index'] and not field.opt['PrimaryKey']:
                indexlist.append(field.info['Name'])
            elif field.opt['PrimaryKey']:
                indexlist = [field.info['Name']] + indexlist
        ### Header ###
        buildlist = ['### %s Table - %s ###' % (self.info['Name'],self.info['Description'])]
        buildlist.append('# - Build Statement generated by rje_mysql from %s' % self.info['File'])
        buildlist.append('# - Generated %s' % time.asctime(time.localtime(time.time())))
        buildlist.append('# - %s total lines: %d header lines & %s bad lines' % (rje.integerString(self.stat['Lines']),self.stat['Screened'],rje.integerString(self.stat['BadLines'])))
        buildlist.append('')
        ### Create Table ###
        buildlist.append('DROP TABLE IF EXISTS %s;' % self.info['Name'])
        buildlist.append('CREATE TABLE %s (' % self.info['Name'])
        ## Fields ##
        for field in self.fieldlist:
            fieldname = field.info['Name']
            while len(fieldname) < fieldlen:
                fieldname += ' '
                #self.deBug('<%s> (%d) vs %d letters' % (fieldname,len(fieldname),fieldlen))
            fieldbuild = '    %s %s' % (fieldname,field.info['Type'])
            if field.opt['Numeric'] and field.opt['Unsigned']:
                fieldbuild += ' UNSIGNED'
            if field.opt['NotNull']:
                fieldbuild += ' NOT NULL'
            else:
                fieldbuild += ' DEFAULT NULL'
            if field.opt['AutoIncrement']:
                fieldbuild += ' AUTO_INCREMENT'
                desc_add = 50
            fieldbuild += ','
            while len(fieldbuild) < (fieldlen+desc_add):
                fieldbuild += ' '
            fieldbuild += '# %s' % field.info['Description']
            buildlist.append(fieldbuild)
        ## Primary Key and Indices ##
        buildlist.append('    PRIMARY KEY (%s)' % indexlist.pop(0))
        while indexlist:
            buildlist[-1] = buildlist[-1] + ','
            buildlist.append('    INDEX (%s)' % indexlist.pop(0))
        buildlist.append('    );')
        ### Load Data ###
        loadcmd = "LOAD DATA LOCAL INFILE '%s' INTO TABLE %s" % (rje.replace(self.info['File'],'\\','\\\\'),self.info['Name'])
        if delimit != '\t':
            loadcmd += " FIELDS TERMINATED BY ',' OPTIONALLY ENCLOSED BY '\"'"
        if self.fieldlist[0].opt['AutoIncrement']:
            loadcmd += ' (%s)' % rje.join(fieldlist[1:],', ')
        loadcmd += ';'
        buildlist.append(loadcmd)
        ### Remove headers ###
        if self.headers != ['']:
            fd = 0
            if self.fieldlist[0].opt['AutoIncrement']:
                fd = 1
            buildlist.append("DELETE FROM %s WHERE %s = '%s';" % (self.info['Name'],self.fieldlist[fd].info['Name'],self.headers[0][:self.fieldlist[fd].stat['MaxLen']]))
        return rje.join(buildlist,sep='\n')
#########################################################################################################################
### End of Table Class                                                                                                  #
#########################################################################################################################
       
                                                    ### ~ ### ~ ###

#########################################################################################################################
### Field Class:                                                                                                        #
#########################################################################################################################
class Field(rje.RJE_Object):     
    '''
    MySQL Field Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Field header name
    - Type = Type of MySQL field
    - Description = Description of Field
    
    Opt:boolean
    - Numeric = Whether it is a numeric field (define subtypes later)
    - Integer = Whether field is an integer
    - Unique = Whether entries should be unique
    - PrimaryKey = Whether it is a primary key
    - Index = Whether to index table on field
    - Unsigned = MySQL option
    - NotNull = MySQL option. No missing data.
    - AutoIncrement = AutoIncrementing field

    Stat:numeric
    - MaxLen = Maximum length in characters of field content
    - MinLen = Minimum length in characters of field content
    - MaxNum = Maximum number of field content
    - MinNum = Minimum number in characters of field content

    Obj:RJE_Objects
    '''
    ### Attributes
    datalist = []   # Used to determine Unique setting if memsaver=F
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','Type','Description']
        - Opt:boolean ['Numeric','Integer','Unique','PrimaryKey','Index','Unsigned','NotNull','AutoIncrement']
        - Stats:float ['MaxLen','MinLen','MaxNum','MinNum']
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','Description']
        self.optlist = ['Numeric','Integer','Unique','PrimaryKey','Index','Unsigned','NotNull','AutoIncrement']
        self.statlist = ['MaxLen','MinLen','MaxNum','MinNum']
        self.objlist = []
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None)
        self.info['Description'] = 'No description given.'
        self.stat['MinLen'] = 1e10
        self.opt['Numeric'] = True
        self.opt['Integer'] = True
        self.opt['Unsigned'] = True
        self.opt['NotNull'] = True
        ### <c> ### Other Attributes
        self.datalist = []  
#########################################################################################################################
    ### <2> ### Class Methods                                                                                           #
#########################################################################################################################
    def _readData(self,data):   ### Updates attributes based on data
        '''
        Updates attributes based on data.
        >> data:str = Read data
        '''
        try:
            ### Unique ###
            if self.opt['Unique']:
                if data in self.datalist:   # Not unique!
                    self.verbose(2,4,'%s repeat of %s: Not Unique...' % (self.info['Name'],data),0)
                    self.opt['Unique'] = False
                    self.datalist = []  # Save memory
                else:
                    self.datalist.append(data)
            ### Null ###
            if data == '':
                self.opt['NotNull'] = False
                return
            ### Number ###
            if self.opt['Numeric']:
                try:
                    number = string.atof(data)
                    if number < self.stat['MinNum']:
                        self.stat['MinNum'] = number
                    if number > self.stat['MaxNum']:
                        self.stat['MaxNum'] = number
                    if self.opt['Integer'] and int(number) != number:
                        self.opt['Integer'] = False
                    if number < 0:
                        self.opt['Unsigned'] = False
                except: # Not a number!
                    self.opt['Numeric'] = False
                    self.opt['Integer'] = False
            ### Length ###
            if len(data) < self.stat['MinLen']:
                self.stat['MinLen'] = len(data)
            if len(data) > self.stat['MaxLen']:
                self.stat['MaxLen'] = len(data)
        except:
            self.log.errorLog('Error in _method(%s)' % _stage,printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _autoType(self):    ### Suggests best data type based on statistics
        '''
        Sets self.info['Type'] best data type based on statistics.
        : Numeric Types
        : - TINYINT -128 to 127
        : - SMALLINT -32,768 to 32,767
        : - MEDIUMINT -8,388,608 to 8,388,607
        : - INT -2,147,483,648 to 2,147,483,647
        : - BIGINT -9,223,372,036,054,775,808 to 9,223,372,036,054,775,807
        : - FLOAT (single 4b)
        : - DOUBLE
        : - DECIMAL(digits,dp)
        : Text Types 
        : - CHAR(n) = fixed length <=255
        : - VARCHAR(n) = variable length <=255
        : - BLOB = case-insensitive >255
        : - TEXT = case-sensitive >255
        : - For TEXT and BLOB:
        : - n < 65,536
        : - TINY... <256
        : - MEDIUM... < 16,777,216
        : - LONG...  <4,294,967,295
        '''
        int_types = ['TINYINT','SMALLINT','MEDIUMINT','INT','BIGINT']
        int_max = {'TINYINT':127,'SMALLINT':32767,'MEDIUMINT':8388607,'INT':2147483647,'BIGINT':9223372036054775807}
        if self.opt['Numeric'] and self.opt['Integer']:
            for type in int_types:
                if self.opt['Unsigned']:
                    (min,max) = (-1-int_max[type],int_max[type])
                else:
                    (min,max) = (0,2*int_max[type]+1)
                if self.stat['MinNum'] >= min and self.stat['MaxNum'] <= max:
                    self.info['Type'] = type
                    #self.info['Type'] = '%s(%d)' % (type,self.stat['MaxLen'])
                    return
            self.info['Type'] = 'DOUBLE'
        elif self.opt['Numeric']:
            self.info['Type'] = 'FLOAT'
        elif self.stat['MaxLen'] <= 255:
            if self.stat['MaxLen'] == self.stat['MinLen']:
                self.info['Type'] = 'CHAR(%d)' % self.stat['MaxLen']
            else:
                self.info['Type'] = 'VARCHAR(%d)' % self.stat['MaxLen']
        else:   # Blob!
            if self.stat['MaxLen'] < 65536:
                self.info['Type'] = 'BLOB'
            elif self.stat['MaxLen'] < 16777216:
                self.info['Type'] = 'MEDIUMBLOB'
            else:
                self.info['Type'] = 'LONGBLOB'
#########################################################################################################################
### End of Field Class                                                                                                  #
#########################################################################################################################
       
#########################################################################################################################
## End of SECTION II                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    try:
        ### <0>  ### Basic Setup of Program
        [info,out,mainlog,cmd_list] = setupProgram()
        
        ### <1> ### Rest...
        mysql = MySQL(mainlog,['append=T']+cmd_list)
        mysql.run()
        

        ### <X> ### End
    except SystemExit:
        return  # Fork exit etc.
    except KeyboardInterrupt:
        mainlog.errorLog("User terminated.\n")
    except:
        print "Unexpected error:", sys.exc_info()[0]
    mainlog.printLog('#LOG', "%s V:%s End: %s\n" % (info.program, info.version, time.asctime(time.localtime(time.time()))), 1)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try:
        runMain()
    except:
        print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
