#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_db
Description:  R Edwards Relational Database module
Version:      1.9.1
Last Edit:    09/03/19
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to read in and store data as a series of tables in a similar fashion to a database, e.g.
    MySQL. Although undoubtedly slower than such dedicated software for querying etc., this module is primarily designed
    for used to make links and manipulate data and tables on the fly in other python modules.

    The main Database Class will control the linking etc. of tables, which are in turn stored in the Table Class.

    NB. This module should not be confused with rje_dbase, which is for downloading and processing public databases.

Commandline:
    dbindex=T/F : Whether to run in "index" mode, storing a file position rather than all data (read only) !Not yet implemented! [False]
    Additional Commandline functionality will be added with time.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_scoring, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added merge tables option.
    # 0.2 - Miscellaneous updates to various methods.
    # 0.3 - Minor doc tweaks and added keepFields().
    # 0.4 - Improved use of AutoID and added Table.autoID() method.
    # 0.5 - Initial coding of index mode. (Not yet fully functional.)
    # 1.0 - Working, so upgraded to version 1.0!
    # 1.1 - Added sortedEntries() function.
    # 1.2 - Added Table.hasField(field). Add openTable(), readEntry() and readSet() methods.
    # 1.3 - Minor modifications for SLiMCore FUPC development.
    # 1.4 - Added list checking with addEmptyTable.
    # 1.5 - Fixed occasional key error following addField. Added indexReport() method.
    # 1.6 - Added option to save a subset of entries using saveToFile(savekeys=LIST).
    # 1.7.0 - Added splitchar to table splitting.
    # 1.7.1 - Reinstated raise error if expected table missing.
    # 1.7.2 - Fixed numerical join issue during Table.compress().
    # 1.7.3 - Added lower case enforcement of headers for reading tables from file.
    # 1.7.4 - Added optional restricted Field set for output.
    # 1.7.5 - Added more error messages and tableNames() method.
    # 1.7.6 - Added table.opt['Formatted'] = Whether table data has been successfully formatted using self.dataFormat()
    # 1.7.7 - Added option to constrain table splitting to certain field values.
    # 1.8.0 - Added option to store keys as tuples for correct sorting. (Make default at some point.)
    # 1.8.1 - Added sfdict to saveTable output.
    # 1.8.2 - Fixed minor readSet bug.
    # 1.8.3 - Minor debugging message changes.
    # 1.8.4 - Cosmetic log message changes.
    # 1.8.5 - Added saveToFileName() function.
    # 1.8.6 - Minor IndexReport tweak.
    # 1.9.0 - Added comment output to saveToFile().
    # 1.9.1 - Updated logging of adding/removing fields: default is now when debugging only.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Create general module class structure of Databases and Tables
    # [Y] : Add methods for loading and cross-linking tables
    # [Y] : Use rje_scoring to add filters etc. to Tables. Also new fields as formulas or strings.
    # [Y] : Add option to add a table from a file, adding a new Primary Key identifer field "AutoID"
    # [ ] : Add methods for generating new tables with custom capabilities, such as the max/min score etc.
    # [ ] : Add "index" mode, which will point to (and read) an entry from a file rather than loading all data.
    # [ ] : Add use of mySQL in place of this object (but using same methods calls).
    # [ ] : Check/add capacity to ignore # comment lines in file header before getting to actual field headers.
    # [ ] : Improve consistency and handling of AutoID fields in joins and merges etc.
    # [ ] : Add ability to have field tuples as indexes.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copy_right) = ('RJE_DB', '1.9.1', 'March 2019', '2008')
    description = 'R Edwards Relational Database module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please report bugs to Richard.Edwards@UNSW.edu.au']
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?','N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Database Class                                                                                          #
#########################################################################################################################
class Database(rje.RJE_Object):
    '''
    Database Class. Author: Rich Edwards (2007).

    Info:str
    - Name = Name for Database
    
    Opt:boolean
    - DBIndex = Whether to run in "index" mode, storing a file position rather than all data (experimental) [False]

    Stat:numeric

    List:list
    - Tables = List of Table Objects

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = []
        self.optlist = ['DBIndex']
        self.statlist = []
        self.listlist = ['Tables']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        ### Other Attributes ###
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdRead(cmd,type='opt',att='DBIndex')  # No need for arg if arg = att.lower()
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def tables(self): return self.list['Tables'][0:]
    def deleteTable(self,table): self.list['Tables'].remove(self.getTable(table))
    def tableNames(self):
        tabnames = []
        for table in self.tables(): tabnames.append(table.name())
        return tabnames
#########################################################################################################################
    def saveDB(self,tables=[],splitfield=None,outdir='',replace=True,log=True,logskip=True):   ### Outputs all tables to files.
        '''
        Outputs all tables to files.
        >> tables:list [] = Optional list of tables to output. If blank will save all tables.
        >> splitfield:str [None] = Optional field to split tables on prior to output. Use split field as basefile.
        >> outdir:str [None] = Optional output directory.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            savex = 0; tabx = 0
            if not tables: tables = self.list['Tables'][0:]
            logskip = logskip and log
            delimit = self.info['Delimit']
            runpath = self.getStr('RunPath') != rje.makePath(os.path.abspath(os.curdir)) or self.debugging()
            ### ~ [1] Regular saving ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not splitfield:
                for table in self.list['Tables']:
                    if outdir: filename = '%s%s.%s.%s' % (outdir,rje.baseFile(self.basefile(),strip_path=True),table.info['Name'],rje.delimitExt(delimit))
                    else: filename = '%s.%s.%s' % (self.basefile(runpath=runpath),table.info['Name'],rje.delimitExt(delimit))
                    if replace or not os.path.exists(filename): table.saveToFile(filename); savex += 1
                    elif logskip and os.path.exists(filename): self.printLog('#SKIP','Skipping save of existing %s' % filename)
                if log: self.printLog('#DBOUT','%d of %d tables saved.' % (savex,len(self.tables())))
                return savex
            ### ~ [2] Special split Table output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for table in tables[0:]:
                if splitfield not in table.fields():
                    self.warnLog('No splitfield=%s output for table %s' % (splitfield,table.name()))
                    continue
                splitdict = self.splitTable(table,splitfield,asdict=True,keepfield=True,splitchar=None)
                for skey in rje.sortKeys(splitdict):
                    stable = splitdict[skey]; tabx += 1
                    filename = '%s%s.%s.%s' % (outdir,skey,table.info['Name'],rje.delimitExt(delimit))
                    if replace or not os.path.exists(filename): stable.saveToFile(filename); savex += 1
                    elif logskip and os.path.exists(filename): self.printLog('#SKIP','Skipping save of existing %s' % filename)
                    self.deleteTable(stable)
            if log: self.printLog('#DBOUT','%s of %s tables saved.' % (rje.iStr(savex),rje.iStr(tabx)))
        except: self.errorLog('Major disaster during db.saveDB() - %s of %s tables saved.' % (rje.iStr(savex),rje.iStr(tabx)))
        return savex
#########################################################################################################################
    def setBasefile(self,basefile=None,cascade=True): ### Sets basefile and cascades to daughter objects
        '''Sets basefile and cascades to daughter objects.'''
        if basefile: self.setStr({'Basefile':basefile})
        else: basefile = self.basefile()
        for obj in self.obj.values() + self.list['Tables']:
            try:
                if obj and cascade and obj.basefile() != basefile: obj.setBasefile(basefile)
            except: pass
#########################################################################################################################
    ### <2> ### Main Table Generation/Manipulation Methods                                                              #
#########################################################################################################################
    def getTable(self,name,errors=False):    ### Returns table object with given name
        '''Returns table object with given name.'''
        try:### ~ [1] Check each table for given name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not name: raise ValueError
            if name in self.list['Tables']: return name     # Actual table given
            for table in self.list['Tables']:
                if table.info['Name'] == name: return table
            if name and type(name) != str: return name      # Assume Table object not in self.list['Tables']
            raise ValueError
        except:
            if errors: self.log.errorLog('Cannot find table "%s"' % name)
        return None
#########################################################################################################################
    def openTable(self,filename=None,mainkeys=[],datakeys='All',delimit=None,headers=[],name=None,expect=True,replace=False):   ### Adds but does not load
        '''
        Adds a table using the settings given but stops after header establishment - does not read data.
        >> name:str [None] = Name for table. If None, will try to extract from filename.
        >> expect:bool [True] = Expect table's existence and raise error if missing (True) or just return None (False)
        >> replace:bool [False] = Whether to replace another table with that name already in self.list['Tables']
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if datakeys != 'All': self.warnLog('Cannot use datakeys with db.openTable()!',dev=True)
            if not delimit and filename: delimit = rje.delimitFromExt(filename=filename)
            elif not delimit: delimit = self.getStr('Delimit')
            if not filename and name: filename = '%s.%s.%s' % (self.basefile(),name,rje.delimitExt(delimit))
            if not rje.exists(filename):
                if expect: raise IOError
                else: return None
            table = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
            if name: table.info['Name'] = name
            else: name = table.info['Name'] = filename

            ### ~ [1] Open file and sort out fields and keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            FILE = table.obj['File'] = open(filename,'r')
            table.list['Fields'] = rje.readDelimit(FILE.readline(),delimit)
            if headers:
                if headers != table.list['Fields']: FILE.seek(0)
                table.list['Fields'] = headers
            else: headers = table.list['Fields']
            if not mainkeys: mainkeys = headers[:1]
            table.list['Keys'] = mainkeys[0:]

            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getTable(name):
                if replace: self.deleteTable(name); self.warnLog('Warning: existing table "%s" replaced.' % name)
                else: self.warnLog('Warning: table "%s" already found - might cause issues.' % name)
            self.list['Tables'].append(table)
            self.printLog('#TABLE','Table "%s" opened: %d fields.' % (table.info['Name'],table.fieldNum()))
            return table
        except:
            if expect: self.errorLog('Error with openTable: %s not added' % filename)
            return None
#########################################################################################################################
    def dbFileName(self,name): return '%s.%s.%s' % (self.basefile(),name,rje.delimitExt(self.getStr('Delimit')))
#########################################################################################################################
    def addTable(self,filename=None,mainkeys=[],datakeys='All',delimit=None,headers=[],ignore=['#'],lists=False,name=None,expect=True,replace=False,uselower=False):   ### Adds
        '''
        Adds a table using the Table.loadDataDict() method (see this for further parameter explanations).
        >> filename:str = Name of file to load from
        >> mainkeys:list = List of headers to be used as key for returned dictionary. If None, will use first header.
        >> datakeys = List of headers to be used as keys for data returned for each mainkey (all headers if [])
        >> delimit = string delimiter. If None, will identify from filename
        >> headers = List of headers to use instead of reading from first line or use datakeys
        >> ignore:list = Leading strings for lines to ignore (e.g. #)
        >> lists:bool [False] = whether to return values as lists. (Otherwise, later entries will overwrite earlier ones)
        >> name:str [None] = Name for table. If None, will try to extract from filename.
        >> expect:bool [True] = Expect table's existence and raise error if missing (True) or just return None (False)
        >> replace:bool [False] = Whether to replace another table with that name already in self.list['Tables']
        >> uselower:bool [False] = Whether to convert headers into lower case.
        << returns new Table object
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not delimit and filename: delimit = rje.delimitFromExt(filename=filename)
            elif not delimit: delimit = self.info['Delimit']
            if not filename and name: filename = '%s.%s.%s' % (self.basefile(),name,rje.delimitExt(delimit))
            #self.deBug(filename)
            if not rje.exists(filename): raise IOError
            try:
                if not mainkeys:
                    fline = open(filename,'r').readline()
                    fline = rje.chomp(fline)
                    mainkeys = string.split(fline,delimit)[:1]
                #string.join(mainkeys)  #?# Why was this here?!
            except: self.errorLog('Problem with mainkeys list given to addTable',printerror=False); raise ValueError
            table = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
            #x#if name =='TP': table.opt['DeBug'] = True
            if len(mainkeys) == 1 and mainkeys[0].lower() in ['all','auto','#']: mainkeys = mainkeys[0].lower()
            if mainkeys not in ['All','all','auto','#'] and datakeys not in ['All','all']:
                try:
                    for key in mainkeys:
                        if key not in datakeys: datakeys.insert(0,key)
                except: pass
            if uselower:
                for flist in [mainkeys,datakeys,headers]:
                    if flist in ['All','all','auto','#']: continue
                    for i in range(len(flist)): flist[i] = flist[i].lower()
            ### ~ [2] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DBIndex'): table.loadDataIndex(filename,mainkeys,datakeys,delimit,headers,ignore,uselower=uselower)
            else: table.loadDataDict(filename,mainkeys,datakeys,delimit,headers,ignore,lists,uselower=uselower)
            table.list['Keys'] = table.list['Keys'][0:]
            if name and self.getTable(name):
                if replace: self.deleteTable(name); self.warnLog('Warning: existing table "%s" replaced.' % name)
                else: self.warnLog('Warning: table "%s" already found - might cause issues.' % name)
            self.list['Tables'].append(table)
            if name and table.getStr('Name') != name:
                self.printLog('#LOAD','Loaded table "%s" renamed "%s".' % (table.getStr('Name'),name))
                table.setStr({'Name':name})
            self.printLog('#TABLE','Table "%s" added: %d fields; %s entries.' % (table.info['Name'],table.fieldNum(),rje.integerString(table.entryNum())))
            return table
        except:
            if expect: self.errorLog('Error with addTable: %s not added' % filename); raise
            return None
#########################################################################################################################
    def addEmptyTable(self,name,fields,keys,log=True): ### Adds empty table
        '''
        Adds an empty table to the self.list['Tables'].
        >> name:str = Name for new table. Should be unique
        >> fields:list = List of fields (strings) for new table.
        >> keys:list = Subset of fields that define unique entries for Table. (Use a counter field if no such combo.)
        << returns new Table object
        '''
        table = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
        table.info['Name'] = name
        table.list['Fields'] = fields[0:]
        table.list['Keys'] = rje.asList(keys)
        self.list['Tables'].append(table)
        self.bugLog('#TABLE','Empty Table "%s" added: %d fields; %s entries.' % (table.info['Name'],table.fieldNum(),rje.integerString(table.entryNum())),log=log,screen=log)
        return table
#########################################################################################################################
    def mergeTables(self,table1,table2,overwrite=True,matchfields=True):    ### Merges table2 into table1 and removes
        '''Merges table2 into table1 and removes.'''
        try:### ~ [1] Perform merge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fields1 = table1.fields(); fields2 = table2.fields()
            if table1.keys() != table2.keys():
                if table1.keys() in [['AutoID'],['#']] and table2.keys() in [['AutoID'],['#']]:
                    table2.renameField(table2.keys()[0],table1.keys()[0])
                else: raise ValueError('Failed to merge %s and %s: different keys!' % (table1.name(),table2.name()))
            if matchfields and fields1 != fields2:
                self.errorLog('%s: %s' % (table1.name(),string.join(fields1,', ')))
                self.errorLog('%s: %s' % (table2.name(),string.join(fields2,', ')))
                raise ValueError
            if table1.keys() in [['AutoID'],['#']] and table1.entryNum():
                tkey = table1.keys()[0]
                table1.dataFormat({'AutoID':'int','#':'int'})
                table2.dataFormat({'AutoID':'int','#':'int'})
                ex = max(table1.dataKeys())
                newdata = {}
                for k2 in table2.dataKeys():
                    ex += 1
                    newdata[ex] = table2.data().pop(k2)
                    newdata[ex][tkey] = ex
                table2.dict['Data'] = newdata
                #!# Renumber IDs for second table
            rje.combineDict(table1.dict['Data'],table2.dict['Data'],overwrite=overwrite)
            if table2 in self.list['Tables']: self.list['Tables'].remove(table2)
            self.printLog('#MERGE','Merged table %s into %s -> %s entries.' % (table2.info['Name'],table1.info['Name'],rje.iLen(table1.entries())))
            ### ~ [2] Check fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not matchfields: return True
            for field in fields2:
                if field not in fields1: table1.list['Fields'].append(field)
            (mx,mtot) = (0.0,table1.entryNum())
            for entry in table1.entries():
                self.progLog('\r#CHECK','Matching Field Entries for merged tables: %.2f%%' % (mx/mtot)); mx += 100.0
                for field in table1.list['Fields']:
                    if field not in entry: entry[field] = ''
            self.printLog('\r#CHECK','Matching Field Entries for merged tables complete.',log=False)
            return table1
        except ValueError: self.errorLog('Field mismatch! Cannot merge %s and %s.' % (table1.name(),table2.name()))
        except: self.errorLog('MergeTables(%s,%s) error' % (table1,table2))
        return False
#########################################################################################################################
    def joinTables(self,name='',join=[],newkey=[],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True,warnings=True):   ### Makes a new table using join of [(Table,Field[,Fieldlist])]
        '''
        Makes a new table by joining existing tables, using the given fields, in priority order listed. If a Field from a
        latter table already exists, it will be added as "Table_Field" instead. If multiple combinations arise for the
        same new key, only the first will be kept.
        >> name:str [''] = Name for new table. If not given will become "TableX"
        >> join:list of (Table,Field[,Fieldlist]) tuples, where Table is a table name and Field is a Field name or
            formula to be used for the join. Fieldlist is an optional list of Fields from that Table to include in the
            new table. If Field does not exist, it will be added. (Field may be a Formula.)
        >> newkey:list [] = If None, will make a new "AutoID" key field. Else, will use the given Fields.
        >> cleanup:bool [True] = If True will delete any Fields generated just to make the join
        >> delimit:str ['\t'] = Delimiter to be used to join the key fields
        >> delimit:str ['\t'] = Delimiter to be used to join the key fields
        >> empties:bool [True] = Whether to keep entries that do not link to 1+ tables with empty values or delete them.
        >> check:bool [False] = Whether to check for entries that are not being joined.
        >> keeptable:bool [True] = Whether to add new table to self.list['Tables']
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(join) < 2: self.errorLog('Cannot join %d table!' % len(join),printerror=False); raise ValueError
            ## ~ [1a] Setup lists etc. for use in method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            deletefields = []   # List of (Table,NewField) to cleanup
            newfields = []      # List of fields for the new table to have
            fieldmap = {}       # Mapping of (newfield:(table,field))
            autoid = False
            if not newkey:
                newkey = ['AutoID']
                newfields = ['AutoID']
                autoid = True
            newtable = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
            newtable.info['Delimit'] = delimit
            ## ~ [1b] Setup new table object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getTable(name):
                self.errorLog('Table "%s" already exists! Will make new name.' % name,printerror=False)
                name = ''
            if not name:
                x = 1
                while self.getTable('Table%d' % x): x += 1
                name = 'Table%d' % x
            newtable.info['Name'] = name
            ## ~ [1c] Add new fields and prepare tables for join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('\r#JOIN','Joining %d tables...' % len(join))
            for jointup in join:
                fieldlist = []              # List of fields to include in new table
                if len(jointup) == 3: (table,field,fieldlist) = jointup
                else: (table,field) = jointup
                table = self.getTable(table,errors=True)
                if not table: raise ValueError
                if not fieldlist: fieldlist = table.list['Fields'][0:]
                if field not in table.list['Fields']:
                    table.makeField(field)
                    if cleanup: deletefields.append((table,field))
                if field in fieldlist and field in newfields: fieldlist.remove(field)   # No need to duplicate
                for tfield in fieldlist:
                    if tfield in table.list['Fields']:
                        if tfield in newfields:
                            if tfield not in ['#','AutoID']: newfields.append('%s_%s' % (table.info['Name'],tfield))
                        else: newfields.append(tfield)
                        fieldmap[newfields[-1]] = (table,tfield)
                    else: self.errorLog('Field "%s" not found in table "%s"' % (tfield,table.info['Name']),printerror=False)
                #?# Should there be some form of indexing here = making a new dictionary of {value:[keys]} #?#
                table.index(field)
            ## ~ [1d] Neaten up fields and keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newtable.list['Fields'] = newfields[0:]
            for nkey in newkey:
                if nkey not in newfields:
                    self.errorLog('Cannot use field "%s" in new key - not found in table fields!' % nkey,printerror=False)
                    raise ValueError

            ### ~ [2] Generate new table and make join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Make shell of new table from first join table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpx = 0    # Temporary key = unique counter. Will replace with real key after join
            newdata = {}
            itable = self.getTable(join[0][0],errors=True)
            ifield = join[0][1]
            tx = 0.0; txx = len(itable.dict['Index'][ifield]); jx = 1
            for ijoin in rje.sortKeys(itable.dict['Index'][ifield]):
                self.progLog('\r#JOIN','Joining table 1 of %d: %.2f%%  ' % (len(join),tx/txx)); tx += 100.0
                for ikey in itable.dict['Index'][ifield][ijoin]:
                    newdata[tmpx] = {}
                    for field in newfields:
                        if autoid and field == 'AutoID': newdata[tmpx][field] = tmpx
                        elif fieldmap[field][0] == itable: newdata[tmpx][field] = itable.dict['Data'][ikey][fieldmap[field][1]]
                    newdata[tmpx][ifield] = ijoin       # May replace if weird renaming happening
                    tmpx += 1
            ## ~ [2b] Add one more table at a time to the new table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for joins in join[1:]:
                jtable = self.getTable(joins[0],errors=True)
                jfield = joins[1]
                (olddata,newdata,tmpx) = (newdata,{},0)
                icheck = jtable.dict['Index'][jfield].keys()[0:]
                jx += 1; tx = 0.0; txx = len(olddata)
                for tkey in rje.sortKeys(olddata):
                    self.progLog('\r#JOIN','Joining table %d of %d: %.2f%%  ' % (jx,len(join),tx/txx)); tx += 100.0
                    ijoin = olddata[tkey][ifield]
                    if ijoin in jtable.dict['Index'][jfield]:   # Add each
                        for jkey in jtable.dict['Index'][jfield][ijoin]:
                            newdata[tmpx] = {}
                            for field in newfields:
                                if field in olddata[tkey]: newdata[tmpx][field] = olddata[tkey][field]
                                elif field in fieldmap and fieldmap[field][0] == jtable: newdata[tmpx][field] = jtable.dict['Data'][jkey][fieldmap[field][1]]
                            newdata[tmpx][ifield] = ijoin       # May replace if weird renaming happening
                            if autoid: newdata[tmpx]['AutoID'] = tmpx
                            tmpx += 1
                        if ijoin in icheck: icheck.remove(ijoin)
                    elif empties:   # Add empty
                        newdata[tmpx] = olddata.pop(tkey)
                        if autoid: newdata[tmpx]['AutoID'] = tmpx
                        tmpx += 1
                    else: continue  # Do nothing
                if check and icheck: self.printLog('\r#JOIN','Warning! %s "%s" values not joined from "%s"' % (rje.integerString(len(icheck)),jfield,jtable.name()))
            ## ~ [2c] Put together using actual keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug(newdata)
            del olddata
            tx = 0.0; txx = len(newdata); ax = 0
            #self.debug(newkey)
            for xkey in rje.sortKeys(newdata):
                self.progLog('\r#JOIN','Joined %d tables. Generating new keys: %.2f%%' % (len(join),tx/txx)); tx += 100.0
                if autoid:
                    ax += 1
                    tkey = rje.preZero(ax,max(len(newdata),9999))
                else:
                    tkey = []
                    for nkey in newkey:
                        if nkey not in newdata[xkey]: newdata[xkey][nkey] = ''  # Must have something for key
                        tkey.append('%s' % newdata[xkey][nkey])
                    if newtable.getBool('TupleKeys'): tkey = tuple(tkey)
                    else: tkey = string.join(tkey,delimit)
                if tkey in newtable.dict['Data']:
                    if warnings: self.warnLog('Duplicate values for join keys "%s" dropped from table "%s"' % (tkey,newtable.info['Name']))
                    self.debug(newdata[xkey])
                else: newtable.dict['Data'][tkey] = newdata.pop(xkey)
            del newdata

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#JOIN','New join table "%s" created: %d fields; %s entries' % (newtable.info['Name'],newtable.fieldNum(),rje.integerString(newtable.entryNum())))
            for (table,field) in deletefields: table.deleteField(field)     # List of (Table,NewField) to cleanup
            newtable.fillBlanks()
            newtable.list['Keys'] = newkey
            if keeptable: self.list['Tables'].append(newtable)
            return newtable
                
        except: self.log.errorLog('Problem during joinTables()')
        return None
#########################################################################################################################
    def splitTable(self,table,field,asdict=False,keepfield=False,splitchar=None,values=[],add=True):   ### Splits table based on unique values of given field
        '''
        Splits table based on unique values of given field. To split on multiple fields, first combine these into a
        single field using table.makeField(). New tables are named X_Y, where X is the name of the first table, and Y is
        the content of the field used for the split. If these tables already exist, a warning will be given.
        >> table:Table to split
        >> field:str = Field to split table on
        >> asdict:bool [False] = whether to return new tables dictionary of {field value:table}
        >> keepfield:bool [False] = whether to retain the field used for the split.
        >> splitchar:str [None] = character to split entries on when making index keys. (No split if no str given)
        >> values:list [] = optional list of field values to restrict split to.
        >> add:bool [True] = Whether to add to self.list['Tables']
        << returns list of new tables
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table = self.getTable(table)
            newtables = []
            ## ~ [1a] Check field name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if field not in table.list['Fields']:
                self.log.errorLog('Field "%s" missing from table "%s"' % (field,table.name()),printerror=False)
                return []
            ## ~ [1b] Index table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            del_index = field not in table.dict['Index']
            index = table.index(field,force=True,splitchar=splitchar)
            if values:
                for ikey in rje.sortKeys(index):    # Add each table in alphabetical order
                    if ikey not in values: index.pop(ikey)
                if len(values) != len(index): self.warnLog('%s %s values given for table "%s" split but only %s found.' % (rje.iLen(values),field,table.name(),rje.iLen(index)))
            self.printLog('#SPLIT','Splitting table "%s" on "%s": %d new tables' % (table.name(),field,len(index)))
            if self.interactive() > 0 and not rje.yesNo('Proceed?'):
                self.printLog('#SPLIT','Split cancelled by user.')
                return []

            ### ~ [2] Split table on field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            idict = {}
            for ikey in rje.sortKeys(index):    # Add each table in alphabetical order
                ## ~ [2a] Create and name table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                itable = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
                newname = '%s_%s' % (table.name(),ikey)
                if self.getTable(newname): self.printLog('#WARNING','Table "%s" already exists!' % newname)
                itable.info['Name'] = newname
                itable.info['Delimit'] = table.info['Delimit']
                for datatype in ['Keys','Fields']: itable.list[datatype] = table.list[datatype][0:]
                ## ~ [2b] Transfer data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for tkey in index[ikey]: itable.dict['Data'][tkey] = rje.combineDict({},table.dict['Data'][tkey])
                if not keepfield: itable.dropField(field)
                ## ~ [2c] Add table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newtables.append(itable)
                if add: self.list['Tables'].append(itable)
                idict[ikey] = itable
                self.printLog('#TABLE','Table "%s" added: %d fields; %s entries.' % (itable.name(),itable.fieldNum(),rje.integerString(itable.entryNum())))

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SPLIT','%s new tables added' % len(newtables))
            if asdict: return idict
            return newtables
        except: self.errorLog('Problem during Database.splitTable()')
        for ntable in newtables:
            if ntable in self.list['Tables']: self.list['Tables'].remove(ntable)
        return []
#########################################################################################################################
    def copyTable(self,table,newname,replace=True,add=True):  ### Makes a copy of the table
        '''Makes a copy of the table.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            table = self.getTable(table)
            itable = Table(self.log,self.cmd_list+['basefile=%s' % self.info['Basefile']])
            if self.getTable(newname):
                if replace:
                    self.list['Tables'].remove(self.getTable(newname))
                    self.printLog('#REPTAB','Table "%s" replaced.' % newname)
                else: self.printLog('#WARNING','Table "%s" already exists!' % newname)
            ### ~ [2] ~ Make new table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            itable.info['Name'] = newname
            itable.info['Basefile'] = table.info['Basefile']
            itable.info['Delimit'] = table.info['Delimit']
            for datatype in ['Keys','Fields']: itable.list[datatype] = table.list[datatype][0:]
            for tkey in table.datakeys(): itable.dict['Data'][tkey] = rje.combineDict({},table.dict['Data'][tkey])
            ## ~ [2a] Add table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if add: self.list['Tables'].append(itable)
            self.bugLog('#TABLE','Table "%s" copied from "%s": %d fields; %s entries.' % (itable.name(),table.name(),itable.fieldNum(),rje.integerString(itable.entryNum())))
            return itable
        except: self.errorLog('Problem during Database.copyTable()')
        return None
#########################################################################################################################
    def fieldConvertPrefix(self,convert,prefix,log=True):    ### Converts all matching field names in all tables
        '''
        Converts all matching field names in all tables. Good for standardising fields from multiple sources.
        >> convert:list (or str) = Prefixes to convert (e.g. Query,Qry)
        >> prefix:str = Prefix to use in place of all convert prefixes
        >> log:bool = Whether to log each renaming
        << (fx,tx) = Return tuple of numbers of (fields,tables) in which conversion took place.
        '''
        fx = tx = 0
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if type(convert) == str: convert = [convert]
            if prefix in convert: convert.remove(prefix)
            ### ~ [2] ~ Convert ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for table in self.tables():
                cfields = []
                for field in table.fields():
                    for conpref in convert:
                        if field.startswith(conpref):
                            newname = prefix + field[len(conpref):]
                            table.renameField(field,newname,log)
                            cfields.append(newname)
                if cfields:
                    self.printLog('#FIELD','Converted %d "%s" field names -> %s' % (len(cfields),table.name(),string.join(cfields,'; ')))
                    fx += len(cfields)
                    tx += 1
        except: self.errorLog('Problem during Database.fieldConvertPrefix()')
        return (fx,tx)
#########################################################################################################################
### End of SECTION II: Database Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: Table Class                                                                                            #
#########################################################################################################################
class Table(rje.RJE_Object):     
    '''
    Table Class. Author: Rich Edwards (2007).

    Info:str
    - Name = Name for Database Table
    - Delimit = Delimiter splitting up keys (if appropriate - see rje.dataDict())
    - Source = Name of source file
    
    Opt:boolean
    - DBIndex = Whether to run in "index" mode, storing a file position rather than all data (experimental) [False]
    - Formatted = Whether table data has been successfully formatted using self.dataFormat()
    - Lists = Whether individual data entries may be list (see rje.dataDict())
    - TupleKeys = Whether to store keys as tuples for correct sorting.

    Stat:numeric

    List:list
    - Keys = List of Data columns used to make unique keys (see rje.dataDict())
    - Fields = List of column headers in data dictionary
    - MatchData = List of data populated during self.readSet()

    Dict:dictionary
    - Data = Dictionary of actual table data dictionaries {key:data}
    - DataTypes = Dictionary of {column:data type}
    - Index = Dictionary of {Index:{Index value:[keylist]}}

    Obj:RJE_Objects
    - File = Open file handle for reading (Index mode)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Delimit','Source']
        self.optlist = ['DBIndex','Formatted','Lists','TupleKeys']
        self.statlist = []
        self.listlist = ['Keys','Fields','MatchData']
        self.dictlist = ['Data','DataTypes','Index']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Delimit':'\t'})
        self.setBool({'Formatted':False,'TupleKeys':False})
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options ### 
                self._cmdRead(cmd,type='opt',att='DBIndex')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='opt',att='TupleKeys')  # No need for arg if arg = att.lower()
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Table attribute methods                                                                                 #
#########################################################################################################################
    def formatted(self): return self.getBool('Formatted')
    def entryNum(self): return len(self.dict['Data'])
    def iNum(self): return rje.iStr(self.entryNum())
    def fieldNum(self): return len(self.list['Fields'])
    def name(self): return self.info['Name']
    def rename(self,newname): self.info['Name'] = newname
    def keys(self): return self.list['Keys'][0:]
    def hasField(self,field): return field in self.fields()
    def fields(self): return self.list['Fields'][0:]
    def entry(self,key): return self.data(key,expect=True)
    def data(self,key=None,makekey={},keylist=[],expect=False):
        '''
        Returns data dictionary. If key given, will just return single entry for key (or None if key not found).
        >> key:str [None] = DataKey for data dictionary
        >> makekey:dict [{}] = Make the key from an example entry dictionary
        >> keylist:list [] = List of key values to join with delimiter or convert to tuple
        >> expect:bool [False] = Whether to raise exception if missing, else return None.
        '''
        if not (key or makekey or keylist): return self.dict['Data']
        try:
            if key: return self.dict['Data'][key]
            elif makekey: return self.dict['Data'][self.makeKey(makekey)]
            elif self.getBool('TupleKeys'): return self.dict['Data'][tuple(keylist)]
            else: return self.dict['Data'][string.join(keylist,self.info['Delimit'])]
        except:
            if expect: raise
            else: return None
    def datakeys(self): return rje.sortKeys(self.dict['Data'])
    def dataKeys(self): return self.datakeys()
    def indexKeys(self,index,force=False,log=True):
        '''Returns sorted keys for given index.'''
        try: return rje.sortKeys(self.index(index,force,log=log))
        except: self.errorLog('%s indexKeys() error!' % self.name()); raise
    def indexReport(self,index,logstr=None,force=False):    ### Summarises index numbers to printLog
        ikeys = self.indexKeys(index,force)
        for ikey in ikeys:
            if logstr and not logstr.startswith('#'): logstr = '#%s' % logstr
            if logstr: self.printLog(logstr,'%s "%s" %s:%s entries.' % (rje.iLen(self.indexEntries(index,ikey)),ikey,self.name(),index))
            else: self.printLog('#%s' % index.upper()[:6],'%s %s %s:%s entries.' % (rje.iLen(self.indexEntries(index,ikey)),ikey,self.name(),index))
        return ikeys
    def entries(self,keys=None,sorted=False):
        '''Returns all entries as a list.'''
        if not keys and not sorted: return self.dict['Data'].values()
        if sorted and not keys: return self.entryList(self.datakeys())
        if type(keys) != list and keys in self.dict['Data']: return self.dict['Data'][keys]
        return self.entryList(keys)
    def clear(self):    ### Deletes data and indices but leaves rest of structure in place
        self.dict['Data'] = {}
        self.dict['Index'] = {}
#########################################################################################################################
    def entrySummary(self,entry,fields=[],invert=False):   ### Returns a string summary of an entry - useful for debugging etc.
        '''Returns a string summary of an entry - useful for debugging etc.'''
        if not fields: fields = self.fields()
        estr = 'Entry "%s"\n' % self.makeKey(entry)
        for field in self.fields():
            if field in fields != invert: estr += '|-- %s:\t%s\n' % (field,entry[field])
        return estr
#########################################################################################################################
    def indexDataKeys(self,index,value=None):     ### Return list of datakeys from index & value
        '''Return list of datakeys from index & value.'''
        try: return self.index(index)[value]
        except: return []
#########################################################################################################################
    def indexEntries(self,index,value=None,asdict=False):     ### Return list of entries from index & value
        '''Return list of entries from index & value.'''
        if not asdict:
            if type(value) == list:
                ilist = []
                for val in value:
                    try: ilist += self.entryList(self.index(index)[val])
                    except: pass
                return ilist
            else:
                try: return self.entryList(self.index(index)[value])
                except: return []
        try:
            idict = {}
            for ikey in self.index(index): idict[ikey] = self.indexEntries(index,ikey)
            return idict
        except: return {}
#########################################################################################################################
    def sortedEntries(self,field,reverse=False):    ### Returns list of entries, sorted by field.                   #V1.1
        '''
        Returns list of entries, sorted by field.
        >> field:str = Field on which to sort entries.append
        >> reverse:bool [False] = Whether to reverse sort (big -> small)
        << list of entries
        '''
        entries = []
        sortkeys = rje.sortKeys(self.index(field),revsort=reverse)
        for ikey in sortkeys: entries += self.indexEntries(field,ikey)
        return entries
#########################################################################################################################
    def entryList(self,keylist):    ### Return list of entries from list of keys
        '''Return list of entries from list of keys.'''
        entries = []
        for key in keylist: entries.append(self.data()[key])
        return entries
#########################################################################################################################
    def remakeKeys(self,warnings=True):   ### Remakes all the keys (following changes to entries)
        '''Remakes all the keys (following changes to entries)'''
        newdict = {}
        for okey in rje.sortKeys(self.dict['Data']):
            entry = self.dict['Data'][okey]
            newkey = self.makeKey(entry)
            if warnings and newkey in newdict: self.printLog('#KEY','Warning: %s over-written!' % str(newkey))
            newdict[newkey] = entry
        self.dict['Data'] = newdict
        self.dict['Index'] = {}
        self.bugLog('#KEY','Keys remade for %s %s table entries.' % (rje.integerString(len(newdict)),self.name()))
#########################################################################################################################
    def addEntry(self,entry,warn=True,overwrite=True,splitchar=None,remake=True): ### Adds entry to self.dict['Data']
        '''Adds entry to self.dict['Data']. Makes key using self.makeKey().'''
        for f in self.fields():
            if f not in entry: entry[f] = ''
        ekey = self.makeKey(entry)
        if ekey in self.dict['Data']:
            if not overwrite: return None
            if warn: self.warnLog('%s entry "%s" being overwritten' % (self.name(),ekey),'entry_overwrite',suppress=True)
        #i# Remake will remove excess fields and break linkage with original dictionary. Otherwise the SAME object is used!
        if remake:
            newentry = {}
            for field in self.fields():
                try: newentry[field] = entry[field]
                except: newentry[field] = ''
            entry = newentry
            self.dict['Data'][ekey] = entry
        else: self.dict['Data'][ekey] = entry
        for ikey in self.dict['Index'].keys():
            try:
                idata = entry[ikey]
                if splitchar:
                    for i in string.split(entry[ikey],splitchar):
                        if i not in self.dict['Index'][ikey]: self.dict['Index'][ikey][i] = [ekey]
                        elif ekey not in self.dict['Index'][ikey][i]:
                            self.dict['Index'][ikey][i].append(ekey)
                            self.dict['Index'][ikey][i].sort()
                else:
                    if idata not in self.dict['Index'][ikey]: self.dict['Index'][ikey][idata] = [ekey]
                    elif ekey not in self.dict['Index'][ikey][idata]:
                        self.dict['Index'][ikey][idata].append(ekey)
                        self.dict['Index'][ikey][idata].sort()
            except: self.dict['Index'].pop(ikey); self.debug('%s not in %s' % (ikey,entry))
        return entry
#########################################################################################################################
    def makeKey(self,entry):    ### Returns what the key should be for a given entry dictionary
        '''
        Returns what the key should be for a given entry dictionary. Uses self.list['Keys'] and self.info['Delimit']
        >> entry:dict = Data dictionary from self.dict['Data'].
        '''
        ### ~ [1] Generate new key from keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        newkey = []
        if len(self.keys()) == 1:
            kfield = self.keys()[0]
            if kfield in ['#','AutoID'] and not entry.get(kfield):
                if self.entryNum(): entry[kfield] = max(self.dataKeys()) + 1
                else: entry[kfield] = 1
            return entry[kfield]
        if self.getBool('TupleKeys'):
            for kfield in self.keys(): newkey.append(entry[kfield])
            return tuple(newkey)
        else:
            for kfield in self.keys(): newkey.append(str(entry[kfield]))
            return string.join(newkey,self.info['Delimit'])
#########################################################################################################################
    def autoID(self,log=True,startx=0):   ### Replaces existing keys with AutoID
        '''Replaces existing keys with AutoID.'''
        try:### ~ [1] Replace existing keys with AutoID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if '#' in self.fields(): self.newKey(['#'],True); return True
            self.list['Fields'].insert(0,'#')
            ax = startx
            for oldkey in self.dataKeys():
                self.dict['Data'][ax] = self.dict['Data'].pop(oldkey)
                self.dict['Data'][ax]['#'] = ax
                ax += 1
            if log: self.printLog('#AUTO','AutoID # field added for %s entries of Table "%s".' % (rje.iStr(self.entryNum()),self.name()))
        except:
            self.errorLog('Problem making AutoID Key for Table "%s"' % self.info['Name']); raise
#########################################################################################################################
    def newKey(self,fieldlist,startfields=False,strict=False):       ### Makes fieldlist new key
        '''
        Makes fieldlist new key.
        >> strict:bool [False] = whether to kill if data is overwritten. Otherwise will warn.
        '''
        oldkeys = self.list['Keys'][0:]
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: string.split(fieldlist); fieldlist = [fieldlist]   # Only works if string, not list! (One field!)
            except: pass
            for field in fieldlist:
                if field not in self.list['Fields']:
                    self.errorLog('New Key "%s" missing from "%s" Fields!' % (field,self.info['Name']),printerror=False)
                    #self.deBug('%s' % self.list['Fields'])
                    raise ValueError
            newdata = {}
            self.list['Keys'] = fieldlist[0:]
            ptxt = 'Making new key (%s) for Table "%s"' % (string.join(fieldlist,'|'),self.info['Name'])
            (ex,etot) = (0.0,self.entryNum())
            ### ~ [2] Make new dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in self.entries():
                self.bugProg('\r#KEY','%s: %.2f%%' % (ptxt,ex/etot)); ex += 100.0
                newkey = self.makeKey(entry)
                if newkey in newdata:
                    self.debug('%s:\n%s\n%s' % (newkey,newdata[newkey],entry))
                    if strict:
                        self.errorLog('New %s Key "%s" is not unique!' % (self.name(),str(newkey)),printerror=False); raise ValueError
                    else:
                        self.warnLog('New %s Key "%s" is not unique!' % (self.name(),str(newkey)))
                newdata[newkey] = entry
            self.dict['Data'] = newdata
            self.bugLog('\r#KEY','%s complete.' % (ptxt))
            ### ~ [3] Field List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if startfields:
                for field in fieldlist: self.list['Fields'].remove(field)
                self.list['Fields'] = fieldlist + self.list['Fields']
            self.dict['Index'] = {}
        except:
            self.errorLog('Problem making new Key for Table "%s"' % self.info['Name'])
            self.list['Keys'] = oldkeys
#########################################################################################################################
    def orderedDataList(self,field,empties=True): return self.dataList(self.entries(sorted=True),field,sortunique=False,empties=empties)
#########################################################################################################################
    def dataList(self,entries,field,sortunique=True,empties=True):   ### Returns list of field values for entries
        '''
        Returns list of field values for entries.
        >> entries:list = list of database table entries
        >> field:str = Field for which to return data from
        >> sortunique:bool = whether data should be sorted and unique
        >> empties:bool = whether to include empty ('') values in the list
        '''
        dlist = []
        for entry in entries: dlist.append(entry[field])
        if sortunique: dlist = rje.sortUnique(dlist)
        while '' in dlist and not empties: dlist.remove('')
        return dlist
#########################################################################################################################
    def indexDataList(self,index,key,field,sortunique=True):
        '''Returns the contents of field for all those entries for which index matches key.'''
        return self.dataList(self.indexEntries(index,key),field,sortunique)
#########################################################################################################################
    ### <3> ### Data loading/saving methods                                                                             #
#########################################################################################################################
    def loadFromFile(self,uselower=False):      ### Loads data from file using own attributes
        '''Loads data from file using own attributes.'''
        try:
            if self.getBool('DBIndex'): self.loadDataIndex(self.info['Source'],mainkeys=self.list['Keys'],datakeys=self.list['Fields'],uselower=uselower)
            else: self.dict['Data'] = rje.dataDict(self,self.info['Source'],mainkeys=self.list['Keys'],datakeys=self.list['Fields'],lists=self.opt['Lists'],uselower=uselower)
            if self.getBool('TupleKeys'): self.remakeKeys()
        except: self.log.errorLog(rje_zen.Zen().wisdom(),quitchoice=True)
#########################################################################################################################
    def loadDataDict(self,filename,mainkeys=[],datakeys='All',delimit=None,headers=[],ignore=['#'],lists=False,add=False,screen=True,uselower=False):   ### Loads
        '''
        Loads data from file and sets own attributes to match.
        >> filename:str = Name of file to load from
        >> mainkeys:list = List of headers to be used as key for returned dictionary. If None, will use first header.
        >> datakeys = List of headers to be used as keys for data returned for each mainkey (all headers if [])
        >> delimit = string delimiter. If None, will identify from filename
        >> headers = List of headers to use instead of reading from first line or use datakeys
        >> ignore:list = Leading strings for lines to ignore (e.g. #)
        >> lists:bool [False] = whether to return values as lists. (Otherwise, later entries will overwrite earlier ones)
        >> add:bool [False] = whether to add to existing data dictionary rather than replacing it
        >> uselower:bool [False] = Whether to convert headers into lower case.
        '''
        try:### ~ [1] Setup own attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not delimit: delimit = rje.delimitFromExt(filename=filename)
            self.progLog('\r#LOAD','Loading data from "%s"...' % (filename),screen=screen)
            if len(mainkeys) == 1 and mainkeys[0].lower() in ['all','auto','#']: mainkeys = mainkeys[0].lower()
            autoid = mainkeys in ['All','all','auto','#']
            if not add:
                self.info['Source'] = filename
                self.info['Name'] = rje.baseFile(filename,True)
                self.list['Keys'] = mainkeys
                self.opt['Lists'] = lists
                self.info['Delimit'] = delimit
                if headers: self.list['Fields'] = headers[0:]
                elif datakeys and datakeys != 'All' and not autoid:
                    for key in mainkeys:
                        if key not in datakeys: datakeys = [key] + datakeys[0:]
                    self.list['Fields'] = datakeys
            getheaders = not self.list['Fields']
            if not mainkeys: mainkeys = datakeys[0:]
            ### ~ [2] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if add: rje.combineDict(self.dict['Data'],rje.dataDict(self,filename,mainkeys,datakeys,delimit,headers,getheaders,ignore,lists,uselower=uselower))
            else: self.dict['Data'] = rje.dataDict(self,filename,mainkeys,datakeys,delimit,headers,getheaders,ignore,lists,uselower=uselower,debug=False) #self.debugging())
            if getheaders: self.list['Fields'] = self.dict['Data'].pop('Headers')[0:]
            for field in self.list['Fields'][0:]:
                while self.list['Fields'].count(field) > 1:
                    self.list['Fields'].reverse()
                    self.list['Fields'].remove(field)
                    self.list['Fields'].reverse()                    
            self.printLog('\r#LOAD','Loading data from "%s" complete.' % (filename),log=False,screen=screen)
            ### ~ [3] Add Auto ID if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if autoid:
                self.list['Fields'] = ['#']+self.list['Fields']
                self.list['Keys'] = ['#']
                ax = 0
                for dkey in rje.sortKeys(self.dict['Data']):
                    entry = self.dict['Data'].pop(dkey); ax += 1
                    #?# Why? entry['#'] = rje.preZero(ax,max(len(self.dict['Data']),9999))
                    entry['#'] = ax
                    self.dict['Data'][entry['#']] = entry
                self.bugLog('#KEY','Added AutoID key for %s entries' % rje.integerString(ax))
            elif self.getBool('TupleKeys'): self.remakeKeys()
        except: self.errorLog('Problem loading data from %s! Check keys: %s' % (filename,mainkeys)); raise
#########################################################################################################################
    def loadDataIndex(self,filename,mainkeys=[],datakeys='All',delimit=None,headers=[],ignore=[],screen=True,uselower=False):   ### Loads
        '''
        Loads data from file into dictionary of file positions and sets own attributes to match.
        >> filename:str = Name of file to load from
        >> mainkeys:list = List of headers to be used as key for returned dictionary. If None, will use first header.
        >> datakeys = List of headers to be used as keys for data returned for each mainkey (all headers if [])
        >> delimit = string delimiter. If None, will identify from filename
        >> headers = List of headers to use instead of reading from first line or use datakeys
        >> ignore:list = Leading strings for lines to ignore (e.g. #)
        >> lists:bool [False] = whether to return values as lists. (Otherwise, later entries will overwrite earlier ones)
        >> add:bool [False] = whether to add to existing data dictionary rather than replacing it
        >> uselower:bool [False] = Whether to convert headers into lower case.
        '''
        try:### ~ [1] Setup own attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            raise ValueError #!# This is not ready! #!#
            if not delimit: delimit = rje.delimitFromExt(filename=filename)
            self.printLog('#LOAD','Loading data from "%s"...' % (filename),newline=False,log=False,screen=screen)
            autoid = mainkeys in ['All','all','auto','#']
            self.info['Source'] = filename
            self.info['Name'] = rje.baseFile(filename,True)
            self.list['Keys'] = mainkeys
            self.opt['Lists'] = lists
            self.info['Delimit'] = delimit
            if headers: self.list['Fields'] = headers
            elif datakeys and datakeys != 'All' and not autoid:
                for key in mainkeys:
                    if key not in datakeys: datakeys = [key] + datakeys[0:]
                self.list['Fields'] = datakeys
            getheaders = not self.list['Fields']
            if not mainkeys: mainkeys = datakeys
            ### ~ [2] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            warnx = 0; warn10 = []
            self.obj['File'] = FILE = open(filename,'r'); fpos = 0
            ## ~ [2a] Get headers, if not given ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while not headers:
                fline = FILE.readline()
                if rje.ignoreLine(fline,ignore): continue
                if uselower: fline = fline.lower()
                headers = rje.readDelimit(fline,delimit)
            ## ~ [2b] Establish keys for data dictionary. Mainkeys form dictionary keys. Datakeys determine which data is returned. ~ ##
            autoid = mainkeys in ['All','all','auto','#']
            if not mainkeys: mainkeys = headers[:1]
            if not datakeys:
                datakeys = headers[0:]
                for key in mainkeys: datakeys.remove(key)
            if datakeys in ['All','all']: datakeys = headers[0:]
            keylen = 0
            if autoid: keylen = len(headers) - 1
            else: 
                for key in mainkeys:
                    try: keylen = max(keylen,headers.index(key))
                    except:
                        self.errorLog('Key "%s" not found in %s headers: %s' % (key,filename,headers))
                        raise ValueError
            ## ~ [2c] Read data from file into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            datadict = {}
            fpos = FILE.tell()
            FILE.seek(0,2); fend = FILE.tell(); FILE.seek(fpos); ix = 0       
            fline = FILE.readline(); 
            while fline:
                fprev = fpos
                fpos = FILE.tell()
                try:
                    ## Check whether line is to be ignored ##
                    if rje.ignoreLine(fline,ignore):
                        fline = FILE.readline()
                        continue
                    ## Convert to data list and check for headers ##
                    data = rje.readDelimit(fline,delimit)
                    if len(data) < keylen:
                        fline = FILE.readline()
                        continue
                    linedata = {}
                    for h in range(len(headers)):
                        try: linedata[headers[h]] = data[h]
                        except: linedata[headers[h]] = ''
                    ## Main Key ##
                    ix += 1
                    if autoid: mainkey = ix
                    else:
                        mainkey = []
                        for key in mainkeys: mainkey.append(linedata[key])
                        if not string.join(mainkey,''): fline = FILE.readline(); continue
                        mainkey = string.join(mainkey,delimit)
                        #x#callobj.deBug('%s:%s' % (mainkey,datadict.has_key(mainkey)))
                    if not datadict.has_key(mainkey): datadict[mainkey] = fprev
                    else:
                        warnx += 1
                        if warnx <= 10: warn10.append(string.replace(mainkey,'\t',','))
                        if debug: self.deBug('Dup: %s' % mainkey)
                    fline = FILE.readline()
                except: self.deBug(fline); raise
            ## ~ [2d] Finish and return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if warnx: 
                self.printLog('\r#WARN','Warning: %s entries overwritten due to common key' % rje.iStr(warnx))
                if warnx > 10: self.printLog('\r#WARN','Dups: %s ...' % string.join(warn10,' | '))
                else: self.printLog('\r#WARN','Dups: %s.' % string.join(warn10,' | '))
            ### ~ [3] Update Table Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Data'] = datadict
            if getheaders: self.list['Fields'] = headers[0:]
            for field in self.list['Fields'][0:]:
                while self.list['Fields'].count(field) > 1:
                    self.list['Fields'].reverse()
                    self.list['Fields'].remove(field)
                    self.list['Fields'].reverse()                    
            self.printLog('\r#LOAD','Loading data from "%s" complete.' % (filename),log=False,screen=screen)
            ### ~ [3] Add Auto ID if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if autoid:
                self.list['Fields'] = ['#']+self.list['Fields']
                self.list['Keys'] = ['#']
                ax = 0
                for dkey in rje.sortKeys(self.dict['Data']):
                    entry = self.dict['Data'].pop(dkey); ax += 1
                    entry['#'] = rje.preZero(ax,max(len(self.dict['Data']),9999))
                    self.dict['Data'][entry['#']] = entry
                self.bugLog('#KEY','Added AutoID key for %s entries' % rje.integerString(ax))
            elif self.getBool('TupleKeys'): self.remakeKeys()
        except: self.errorLog('Problem loading data from %s! Check keys: %s' % (filename,mainkeys)); raise
#########################################################################################################################
    def saveToFileName(self,delimit=None):  ### Returns delimited file name
        try:### ~ [1] Setup parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not delimit: delimit = self.info['Delimit']
            runpath = self.getStr('RunPath') != rje.makePath(os.path.abspath(os.curdir)) or self.debugging()
            filename = '%s.%s.%s' % (self.basefile(runpath=runpath),self.info['Name'],rje.delimitExt(delimit))
            return filename
        except: self.errorLog('Problem generating table "%s" filename' % (self.info['Name']))
#########################################################################################################################
    def saveToFile(self,filename=None,delimit=None,backup=True,append=False,savekeys=[],savefields=[],sfdict={},log=True,headers=True,comments=[],buglog=False):    ### Saves data to delimited file
        '''
        Saves data to delimited file.
        >> filename:str [None] = Output file name (will use self.info['Name'] if None)
        >> delimit:str [None] = Delimiter. Will glean from filename or use self.info['Delimit'] if no name
        >> backup:bool [True] = Whether to check for, and backup, existing file
        >> append:bool [True] = Whether to append existing file
        >> savekeys:list [] = Optional list of keys to save subset of Table. Can be used for sorted output.
        >> savefields:list [] = Optional list of fields to save subset of Table.
        >> sfdict:dict {} = Optional dictionary of {field:sigfig} for formatting output.
        >> log:bool [True] = Whether to log output.
        >> headers:bool [True] = Whether to output the field headers.
        >> comments:list [] = Add a list of comment lines to the start of the file if not appending. Should usually start with #.
        >> buglog:bool [False] = Special debugging output to log.
        '''
        try:### ~ [1] Setup parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename and not delimit: delimit = rje.delimitFromExt(filename=filename,write=True)
            elif not delimit: delimit = self.info['Delimit']
            if not filename: filename = self.saveToFileName(delimit)
            rje.mkDir(self,filename)
            if backup and not append: rje.backup(self,filename,appendable=False)
            elif os.path.exists(filename) and not append:
                try: os.unlink(filename)
                except:
                    if self.debugging(): self.errorLog('Odd behaviour trying to delete %s' % filename)
            if savekeys: outkeys = savekeys[0:]
            else: outkeys = rje.sortKeys(self.dict['Data'])
            if not savefields: savefields = self.list['Fields'][0:]
            for field in savefields:
                if '*' in sfdict and '*' not in self.fields() and field not in sfdict: sfdict[field] = sfdict['*']
            ### ~ [2] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log: self.progLog('\r#SAVE','Saving table "%s"...' % (self.info['Name']))
            if append and rje.exists(filename): OUT = open(filename,'a')
            else:
                OUT = open(filename,'w')
                if comments:
                    hashwarn = 0
                    for comment in comments:
                        if comment.startswith('@') and filename.endswith('.sam'): pass
                        elif not comment.startswith('#'): hashwarn += 1
                        OUT.write('%s\n' % comment)
                    if hashwarn and log:
                        self.warnLog('%d of %d comments in "%s" did not start with "#"!' % (hashwarn,len(comments),filename),warntype="hashwarn",suppress=True)
                    if log: self.printLog('#SAVE','%s: %d leading comments' % (filename,len(comments)))
                outlist = []
                for field in savefields: outlist.append('%s' % field)
                if headers: OUT.write('%s\n' % string.join(outlist,delimit))
            (p,px) = (0.0,len(outkeys))
            sx = 0
            for key in outkeys:
                if log: self.progLog('\r#SAVE','Saving table "%s": %.1f%%' % (self.info['Name'],p/px)); p += 100.0
                entry = self.dict['Data'][key]
                outlist = []
                for field in savefields:
                    if field in entry:
                        if field in sfdict:
                            try: outlist.append('%s' % rje.sf(entry[field],sfdict[field]))
                            except: outlist.append('%s' % entry[field])
                        else: outlist.append('%s' % entry[field])
                    else: outlist.append('')
                    if outlist[-1].find(delimit) >= 0:
                        outlist[-1].replace('"','\"')
                        outlist[-1] = '"%s"' % outlist[-1]
                if self.debugging() and buglog: self.printLog('#BUGOUT',string.join(outlist,delimit))
                OUT.write('%s\n' % string.join(outlist,delimit)); sx += 1
            OUT.close()
            if log:
                if sx: self.printLog('\r#SAVE','Table "%s" saved to "%s": %s entries.' % (self.info['Name'],filename,rje.iStr(sx)))
                else: self.printLog('\r#SAVE','Table "%s" saved to "%s": headers only.' % (self.info['Name'],filename))
            return filename
        except: self.errorLog('Problem saving table "%s" to "%s"' % (self.info['Name'],filename))
#########################################################################################################################
    def readEntry(self,add=True,close=True):    ### Reads next entry from open file                                 #V2.0
        '''Reads next entry from open file.'''
        try:### ~ [1] Read line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = rje.readDelimit(self.obj['File'].readline(),self.getStr('Delimit'))
            if not data:
                if close: self.obj['File'].close(); self.obj['File'] = None
                return None
            ### ~ [2] Convert to dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            entry = {}
            for field in self.fields(): entry[field] = data.pop(0)
            if data: self.warnLog('Excess %s fields during readEntry(%s)' % (self.name(),data),'excessdata',suppress=True)
            if add: self.addEntry(entry)
            return entry
        except: self.errorLog('Problem during %s readEntry()' % self.name()); return None
#########################################################################################################################
    def readSet(self,fields,entry=None,clear=True):    ### Reads in a set of entries with matching fields data     #V2.0
        '''
        Reads in a set of entries with matching fields data.
        >> fields:list = Fields to match (starting with next entry and assuming sorted file if entry not given)
        >> entry:dict [None] = Entry to use for data in place of next read entry. [NOT WORKING!]
        >> clear:bool [True] = Whether to clear existing entries before reading.
        << returns list of field values matching the set read in
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add check of fields? #!#
            #!# NOTE: Giving an entry does not seem to work!! #!#
            FILE = self.obj['File']
            if not FILE: return []  # Closed or failed to open
            fstart = FILE.tell()
            fend = rje.endPos(FILE)
            if clear: self.dict['Data'] = {}
            self.list['MatchData'] = matchdata = []; prex = self.entryNum()
            esorted = not entry  # Whether to scan through until match broken
            if not entry: entry = self.readEntry(add=True)
            if not entry: return []
            for field in fields: matchdata.append(entry[field])
            if not esorted: FILE.seek(0)
            ### ~ [1] Read entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('\r#READ','Reading %s entries (%s)...' % (self.name(),string.join(matchdata,'|')))
            while FILE:
                fpos = FILE.tell()
                nextentry = self.readEntry(add=False,close=esorted)
                for f in range(len(fields)):
                    if not nextentry or nextentry[fields[f]] != matchdata[f]: nextentry = None; break
                if not nextentry and esorted:
                    try: FILE.seek(fpos)    # Reset for next batch
                    except: pass            # Should mean that file object was closed
                    break                   # Stop looping
                elif nextentry: self.addEntry(nextentry)
                elif FILE.tell() >= fend: break
            if prex: self.printLog('\r#READ','%s %s entries read (%s) -> %s entries' % (rje.iStr(self.entryNum()-prex),self.name(),string.join(matchdata,'|'),rje.iStr(self.entryNum())))
            else: self.printLog('\r#READ','%s %s entries read (%s).   ' % (rje.iStr(self.entryNum()),self.name(),string.join(matchdata,'|')))
            if not esorted: FILE.seek(fstart)
            return matchdata
        except: self.errorLog('Something went wrong during %s readSet()' % self.name())
#########################################################################################################################
    ### <4> ### Make/customise field methods                                                                            #
#########################################################################################################################
    def newField(self,fieldname,after='',evalue=None,log=None): return self.makeField(fieldname=fieldname,after=after,evalue=evalue,log=log)
    def addField(self,fieldname,after='',evalue=None,log=None): return self.makeField(fieldname=fieldname,after=after,evalue=evalue,log=log)
    def addFields(self,fields,after='',evalue=None,log=None):
        for fieldname in fields: self.makeField(fieldname=fieldname,after=after,evalue=evalue,log=log)
    def makeField(self,formula='',fieldname='',after='',evalue=None,log=None,warn=True): ### Adds a field using data, If cannot calculate, will concatenate
        '''
        Adds a field using data, If cannot calculate, will concatenate using #Field# replacements (i.e. the formula
        should contain #Field# and this will be replaced (w/o #) with the contents of that field.
        >> formula:str [''] = formula for generating new field data.
        >> fieldname:str [''] = Name of new field
        >> after:str [''] = Name of field for new field to follow (at end if blank)
        >> log:bool [None] = Whether to output to log. If None, will set to self.debugging()
        >> warn:bool [True] = Whether to warn if field already exists
        << Returns True/False whether field added or not
        '''
        try:### ~ [1] Setup and test formula ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log == None: log = self.debugging()
            ## ~ [1a] Check for existence of field already ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not fieldname: fieldname = formula
            if fieldname in self.list['Fields']:
                if warn: self.warnLog('Cannot add field "%s" to table "%s": already exists!' % (fieldname,self.info['Name']))
                return False
            ## ~ [1b] Check whether field can be calculated using rje.formula() ~~~~~~~~~~~~~~~~~~~ ##
            if log: self.progLog('\r#FIELD','Adding "%s" to "%s"' % (fieldname,self.info['Name']))
            if formula: calculate = rje.formula(None,formula,varlist=self.list['Fields'],check=True,calculate=False)
            else: calculate = False
            
            ### ~ [2] Add field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (p,px) = (0.0,self.entryNum())
            for key in self.dict['Data']:
                if log: self.progLog('\r#FIELD','Adding "%s" to "%s": %.1f%%' % (fieldname,self.info['Name'],p/px))
                p += 100.0
                try:
                    data = self.dict['Data'][key]
                    if calculate: data[fieldname] = rje.formula(self,formula,data,varlist=self.list['Fields'],check=False,calculate=True)
                    elif formula:
                        value = formula[0:]
                        for field in self.list['Fields']: value = value.replace('#%s#' % field,'%s' % data[field])
                        data[fieldname] = value
                    elif evalue != None: data[fieldname] = evalue
                    else: data[fieldname] = ''
                except:
                    print self.dict['Data']
                    print data, fieldname; raise
            if after and after in self.list['Fields']: self.list['Fields'].insert(self.list['Fields'].index(after)+1,fieldname)
            else: self.list['Fields'].append(fieldname)
            if log: self.printLog('\r#FIELD','Added field "%s" to table "%s"' % (fieldname,self.info['Name']))
            return True
        except: return self.log.errorLog(rje_zen.Zen().wisdom(),quitchoice=True)
#########################################################################################################################
    def dropField(self,field,log=None): return self.deleteField(field,log)
    def deleteField(self,field,log=None):    ### Deletes given field from table
        '''Deletes given field from table.'''
        try:### ~ [1] Delete field from table attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log == None: log = self.debugging()
            newkey = False
            if field in self.dict['Index']: self.dict['Index'].pop(field)
            if field in self.list['Keys']: self.list['Keys'].remove(field); newkey = True
            if field in self.list['Fields']: self.list['Fields'].remove(field)
            elif not newkey: return    # Nothing to delete
            ### ~ [2] Delete field from table data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for data in self.dict['Data'].values():
                try: data.pop(field)
                except: pass
            if log: self.bugLog('#FIELD','Field "%s" removed from table "%s"' % (field,self.info['Name']))
            if newkey:
                newdata = {}
                for oldkey in self.dataKeys():
                    entry = self.dict['Data'].pop(oldkey)
                    newdata[self.makeKey(entry)] = entry
                self.dict['Data'] = newdata
                if log: self.bugLog('#Key','Key updated ("%s" deleted) in table "%s".' % (field,self.info['Name']))
        except: return self.log.errorLog(rje_zen.Zen().wisdom(),quitchoice=True)
#########################################################################################################################
    def renameField(self,field,newname,log=True):    ### Renames field in table
        '''Renames field in table.'''
        try:### ~ [1] Change field in table attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newkey = False
            if field == newname: return
            if field in self.list['Fields']: self.list['Fields'][self.list['Fields'].index(field)] = newname
            else: self.errorLog('Cannot find %s in %s fields' % (field,self.name())); raise ValueError
            if field in self.dict['Index']: self.dict['Index'][newname] = self.dict['Index'].pop(field)
            if field in self.list['Keys']: self.list['Keys'][self.list['Keys'].index(field)] = newname; newkey = True
            ### ~ [2] Change field in  table data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mx = 0
            for data in self.entries():
                try: data[newname] = data.pop(field)
                except: mx += 1; pass
            if log: self.bugLog('#FIELD','Field "%s" renamed "%s" in table "%s" (%s entries w/o data)' % (field,newname,self.info['Name'],rje.integerString(mx)))
            if newkey:
                newdata = {}
                for oldkey in self.dataKeys():
                    entry = self.dict['Data'].pop(oldkey)
                    newdata[self.makeKey(entry)] = entry
                self.dict['Data'] = newdata
                if log: self.bugLog('#Key','Key updated ("%s" renamed "%s") in table "%s".' % (field,newname,self.info['Name']))
        except:
            #try: self.deBug(entry)
            #except: pass
            return self.log.errorLog('Problem renaming "%s" field "%s"' % (self.info['Name'],field))
#########################################################################################################################
    def splitField(self,field,splitlist,split='|',replace=False):   ### Splits field into splitlist fields using split character
        '''
        Splits field into splitlist fields using split character.
        >> field:str = Field to split
        >> splitlist:list = Names of new Fields to create
        >> split:str = character for string.split()
        >> replace:bool [False] = Whether to delete field when finished
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for newf in splitlist:
                if newf in self.fields(): return self.printLog('Cannot split field to overwrite exisiting field',printerror=False)
            self.list['Fields'] += splitlist
            ### ~ [2] Split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dkey in self.data():
                newdata = string.split(self.data()[dkey][field],split)
                if len(newdata) != len(splitlist): self.log.errorLog('Split field for "%s" gives wrong number of data' % dkey,printerror=False)
                for i in range(len(splitlist)): self.dict['Data'][dkey][splitlist[i]] = newdata[i]
            self.printLog('#SPLIT','%d new fields generated for "%s" from "%s"' % (len(splitlist),self.name(),field))
            if replace: self.deleteField(field)
            return True
        except: self.log.errorLog(rje_zen.Zen().wisdom())
        for newf in splitlist: self.deleteField(newf)
#########################################################################################################################
    def joinFields(self,field,joinlist,join='|',replace=False):   ### Joins joinlist into field using join character
        '''
        Makes a new field be joining fields in joinlist
        >> field:str = Name of new Field to create
        >> joinlist:list = Names of Fields to join
        >> join:str = character for string.join()
        >> replace:bool [False] = Whether to delete old fields when finished
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if field in self.fields(): return self.errorLog('Cannot join fields to overwrite existing field',printerror=False)
            for oldf in joinlist:
                if oldf not in self.fields(): return self.errorLog('"%s" Field "%s" does not exist' % (self.info['Name'],oldf),printerror=False)
            ### ~ [2] Join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for dkey in self.data():
                newdata = []
                for jfield in joinlist: newdata.append(str(self.data()[dkey][jfield]))
                self.dict['Data'][dkey][field] = string.join(newdata,join)
            self.list['Fields'].append(field)
            self.bugLog('#JOIN','New field "%s" (%s) generated for "%s"' % (field,string.join(joinlist,join),self.name()))
            if replace:
                for oldf in joinlist: self.deleteField(oldf)
            return True
        except: self.log.errorLog(rje_zen.Zen().wisdom())
        self.deleteField(field)        
#########################################################################################################################
    def rankField(self,field,newfield='',rev=False,absolute=True,lowest=False,unique=False,warn=True): ### Add ranks of field as new field
        '''
        Add ranks of field as new field.
        >> field:str = Field to rank on
        >> newfield:str [field.Rank] = New field name 
        >> rev:Boolean [False] = if True will return 0 for Highest
        >> absolute:boolean [True] = return 1 to n, rather than 0 to 1
        >> lowest:boolean [False] = returns lowest rank rather mean rank in case of ties
        >> unique:boolean [False] = give each element a unique rank (ties rank in random order)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not newfield: newfield = '%s.Rank' % field
            if newfield in self.fields() and warn: self.warnLog('Field "%s" will be over-written with new ranking' % newfield)
            scorelist = []
            entries = self.entries()    # Fix order
            for entry in entries: scorelist.append(entry[field])
            ### ~ [1] Rank and add newfield ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ranks = rje.rankList(scorelist,rev,absolute,lowest,unique)
            for entry in entries: entry[newfield] = ranks.pop(0)
            if newfield not in self.fields(): self.list['Fields'].append(newfield)
            self.printLog('#RANK','Added ranked "%s" as field "%s" to Table "%s"' % (field,newfield,self.info['Name']))
        except: self.errorLog('Problem with Table.rankField()')
#########################################################################################################################
    def rankFieldByIndex(self,index,rankfield,newfield='',rev=False,absolute=True,lowest=False,unique=False,warn=True):    ### Add ranks of field with index group as new field
        '''
        Add ranks of field with index group as new field.
        >> index:str = field upon which to separate rank lists.
        >> rankfield:str = field with data to rank.
        >> newfield:str [rankfield.Rank] = new rank field name.
        >> rev:Boolean [False] = if True will return 0 for Highest
        >> absolute:boolean [True] = return 1 to n, rather than 0 to 1
        >> lowest:boolean [False] = returns lowest rank rather mean rank in case of ties
        >> unique:boolean [False] = give each element a unique rank (ties rank in random order)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.index(index)
            if not newfield: newfield = '%s.Rank' % rankfield
            if newfield in self.fields() and warn: self.warnLog('Field "%s" will be over-written with new ranking' % newfield)
            ### ~ [1] Rank and add newfield ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for group in self.dict['Index'][index]:
                scorelist = []
                entries = self.entryList(self.dict['Index'][index][group])    # Fix order
                for entry in entries: scorelist.append(entry[rankfield])
                ranks = rje.rankList(scorelist,rev,absolute,lowest,unique)
                for entry in entries: entry[newfield] = ranks.pop(0)
            self.list['Fields'].append(newfield)
        except: self.errorLog('Problem with Table.rankFieldByIndex()')
#########################################################################################################################
    ### <5> ### Indexing methods                                                                                        #
#########################################################################################################################
    def index(self,index,force=False,fullreport=False,log=True,make=False,splitchar=None):   ### Indexes table on field given and returns dictionary
        '''
        Indexes table on field given.
        >> index:str = Field on which to index table. 
        >> force:bool [False] = whether to regenerate existing index
        >> fullreport:bool [False] = whether to report if index already exists
        >> log:bool [True] = whether to report to log.errorLog
        >> make:bool [False] = whether to make create new using formula if missing.
        >> splitchar:str [None] = character to split entries on when making index keys. (No split if no str given)
        << returns index dictionary, so can be used access index, creating if necessary
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            itype = type(index)     # Whether index keys are a string (single field) or list (2+ fields)
            if itype not in [str,list,tuple]: raise TypeError('Index keys must be str, tuple or list, not %s' % itype)
            if itype == list:     # Given a list of fields. NOTE: May cause issues outside self.index()
                index = tuple(index); itype = tuple
            if index in self.dict['Index'] and not force:
                if fullreport and log: self.bugLog('#INDEX','Table "%s" indexed on "%s" already' % (self.info['Name'],index))
                return self.dict['Index'][index]
            self.dict['Index'][index] = {}
            if itype == tuple:
                for field in index:
                    if field not in self.list['Fields']: raise ValueError('Table "%s" has no Field "%s" for Index' % (self.info['Name'],field))
            elif index not in self.list['Fields']:
                if make or (self.i() > 0 and rje.yesNo('Table %s missing field "%s" for index. Make?' % (self.info['Name'],index),'N')): self.makeField(index,log=log)
                else:
                    raise ValueError('Table "%s" has no Field "%s" to Index' % (self.info['Name'],index))
            missing = 0
            ### ~ [2] Generate Index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (p,px) = (0.0,self.entryNum())
            for key in rje.sortKeys(self.dict['Data']):
                if log: self.bugLog('\r#INDEX','Indexing "%s" on "%s": %.1f%%' % (self.info['Name'],index,p/px),newline=False,log=False)
                p += 100.0
                try:
                    if itype == tuple:
                        tkeys = {}  # Field: [Values]
                        for field in index:
                            tkeys[field] = [self.dict['Data'][key][field]]
                            if splitchar: tkeys[field] = string.split(self.dict['Data'][key][field],splitchar)
                        # Now need to make all possible combos
                        ikeys = [()]
                        ifields = list(index)[0:]
                        while ifields:
                            prekeys = ikeys
                            ikeys = []
                            field = ifields.pop(0)
                            for pkey in prekeys:
                                for ival in tkeys[field]:
                                    combo = pkey + (ival,)
                                    if combo not in ikeys: ikeys.append(combo)
                        for i in ikeys:
                            if i not in self.dict['Index'][index]: self.dict['Index'][index][i] = []
                            self.dict['Index'][index][i].append(key)
                    else:
                        if splitchar:
                            for i in string.split(self.dict['Data'][key][index],splitchar):
                                if i not in self.dict['Index'][index]: self.dict['Index'][index][i] = []
                                if key not in self.dict['Index'][index][i]: self.dict['Index'][index][i].append(key)
                        else:
                            i = self.dict['Data'][key][index]
                            if i not in self.dict['Index'][index]: self.dict['Index'][index][i] = []
                            self.dict['Index'][index][i].append(key)
                except KeyboardInterrupt: raise
                except:
                    if self.debugging(): self.errorLog('Oops')
                    missing += 1
            if log: self.bugLog('\r#INDEX','Indexed "%s" on "%s": %s unique; %s missing entries' % (self.info['Name'],index,rje.integerString(len(self.dict['Index'][index])),rje.integerString(missing)),log=log)
            return self.dict['Index'][index]
        except: self.log.errorLog('Problem with Table.index()')
        try:
            self.dict['Index'].pop(index)
            self.log.errorLog('Index "%s" may be damaged: deleted.' % index,printerror=False)
        except: pass
        raise ValueError('Failed to generate "%s" index for %s.' % (index, self.name()))
#########################################################################################################################
    def subset(self,index,value,copy=False):   ### Returns a subset of the data as a dictionary
        '''Returns a subset of the data as a dictionary.'''
        try:### ~ [1] Make a new dictionary from index keys & entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            subdata = {}
            for ekey in self.index(index)[value]:
                if copy: subdata[ekey] = rje.combineDict({},self.dict['Data'][ekey],copyblanks=True)
                else: subdata[ekey] = self.dict['Data'][ekey]
            return subdata
        except: return {}
#########################################################################################################################
    ### <6> ### Major Table Manipulation Methods                                                                        #
#########################################################################################################################
    def reshapeWide(self,key,reshape=[],evalue=None):     ### Reshapes the table to add Fields and reduce entries
        '''
        Reshapes the table to add Fields and reduce entries. Fields listed in reshape will be duplicated. The chosen
        key must be part of the Key. All remaining fields should be constant for any given combination of the remaining
        Keys. New fields are named X|Y, where X is the original field name and Y is the value of key
        >> key:str = Key field containing data to use as suffix for reshaped fields
        >> reshape:list = List of fields to be duplicated based on field
        >> evalue:str = Value for missing combinations.
        '''
        try:### ~ [1] Setup and check suitability of given fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if key not in self.keys(): return self.log.errorLog('Cannot reshapeWide(): key "%s" not in keys' % key, printerror=False)
            for field in reshape:
                if field in self.keys(): return self.log.errorLog('Cannot reshapeWide(): reshape fields cannot be in keys', printerror=False)
            ### ~ [2] Remove key and add new fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            etot = self.entryNum(); ftot = self.fieldNum(); ex = 0.0
            self.list['Keys'].remove(key)
            newfields = []
            for keyval in rje.sortKeys(self.index(key)):
                for field in reshape: self.addField(fieldname='%s|%s' % (field,keyval)); newfields.append('%s|%s' % (field,keyval))
                for entry in self.indexEntries(key,keyval):
                    self.progLog('\r#WIDE','Adding wide field data: %.2f%%' % (ex/etot)); ex += 100.0
                    for field in reshape: entry['%s|%s' % (field,keyval)] = entry[field]
            self.printLog('\r#WIDE','Adding wide field data: %d -> %d Fields.' % (ftot,self.fieldNum()))
            for field in [key] + reshape: self.deleteField(field)
            ### ~ [3] Compress ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.compress(self.keys(),default='text')
            if evalue != None:
                for field in newfields:
                    for entry in self.entries():
                        if not entry[field]: entry[field] = evalue

            self.printLog('\r#WIDE','Reshaped wide on %s field: %s -> %s Entries.' % (key,rje.integerString(etot),rje.integerString(self.entryNum())))
            return True
        except:
            self.deBug(self.info); self.deBug(self.list)
            return self.errorLog('Major problem during Table.reshapeWide()')            
#########################################################################################################################
    def reshapeLong(self,newfield,reshape=[],prog=True):  ### Reshapes the table to add entries and compress fields
        '''
        Reshapes the table to add entries and compress fields. This is the reverse of reshapeWide().
        >> newfield:str = This is the name of the field to be created. It will be added to the keys
        >> reshape:list = List of fields to create from X|Y, where X is the field and Y is the value to go into newfield
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if newfield in self.fields(): return self.errorLog('Cannot reshapeLong(): newfield already exists!', printerror=False)
            for field in reshape:
                if field in self.fields(): return self.errorLog('Cannot reshapeLong():  reshape field already exists!', printerror=False)
            ## ~ [1a] Evaluate new fields and data values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            widefields = {}     # Dictionary of {reshape:[widefieldlist]}
            widevalues = []     # List of the different values that should end up in newfield
            for pre in reshape:
                widefields[pre] = []
                for field in self.fields():
                    if field.find('%s|' % pre) == 0:
                        widefields[pre].append(field)
                        val = field[len('%s|' % pre):]
                        if val not in widevalues: widevalues.append(val)
            if prog: self.printLog('#LONG','Found %s reshape values from %s fields' % (len(widevalues),len(widefields)))
            ## ~ [1b] Update key and field lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['Keys'].append(newfield)
            self.list['Fields'].append(newfield)
            for pre in reshape:
                self.list['Fields'].append(pre)
                for old in widefields[pre]: self.list['Fields'].remove(old)

            ### ~ [2] Reshape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newdata = {}
            (dx,dtot) = (0.0,self.entryNum())
            old_data = self.data()
            for dkey in rje.sortKeys(old_data):        # Work through each existing key
                if prog: self.progLog('\r#LONG','Reshaping long: %.2f%%' % (dx/dtot)); dx += 100.0
                oldwide = old_data.pop(dkey)
                for val in widevalues:          # Work through each value for the new field
                    valdata = {}
                    for oldfield in oldwide:
                        if oldfield in self.list['Fields']: valdata[oldfield] = oldwide[oldfield]
                    valdata[newfield] = val     # This will now form part of the key
                    for pre in reshape:         # Pull out relevant variant of each reshape field
                        try: valdata[pre] = oldwide['%s|%s' % (pre,val)]
                        except: self.errorLog('No field "%s|%s" found!' % (pre,val),printerror=False)
                    newdata[self.makeKey(valdata)] = valdata
            if prog: self.printLog('\r#LONG','Reshaped %d fields long.' % len(widefields))

            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Data'] = newdata
            #for pre in reshape:
            #    for old in widefields[pre]: self.deleteField(old)
            self.printLog('\r#LONG','Reshaped table "%s": %d fields; %s entries' % (self.name(),self.fieldNum(),rje.integerString(self.entryNum())))
            return True
        except: return self.log.errorLog('Major problem during Table.reshapeLong()')            
#########################################################################################################################
    def compress(self,newkeys,rules={},default='mean',best=[],joinchar='|'): ### Compresses table to newkeys, keeping only one entry per key
        '''
        Compresses table to newkeys, keeping only one entry per key. Multiple entries are compressed according to rules:
        - min/max/mean/median/geomean/sum = ways to treat numeric data
        - str/text will treat field as text and just use the first non-empty value for it
        - min/max will use the first/last sorted (non-empty) string datum. All others will use the first unsorted datum.
        - list will concatenate entries, separated by joinchar (sorted & unique)
        >> newkeys:list = Fields to be used for new unique keys
        >> rules:dict = Dictionary of {Field:compression rule}
        >> default:str = Default compression rule for fields not in rules.
        >> best:list = List of fields to work through in order, keeping only the best. The first field is looked at first
            and subsequent fields only looked at in the case of a tie. If a field is not present in rules with 'min', the
            maximum value is used in each case.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: string.join(newkeys)
            except: newkeys = [newkeys]
            for field in newkeys:
                if field not in self.fields(): return self.log.errorLog('Cannot compress(): new key field "%s" missing!' % field, printerror=False)
            ## ~ [1a] Check/update rules ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            okrules = ['str','text','min','max','mean','median','geomean','list','sum']
            if default not in okrules: return self.log.errorLog('Cannot compress(): default rule "%s" not recognised!' % default, printerror=False)
            for field in self.fields():
                if field in newkeys: continue
                if field not in rules: rules[field] = default
                elif rules[field] not in okrules: return self.log.errorLog('Cannot compress(): rule "%s" not recognised!' % rules[field], printerror=False)
            ## ~ [1b] Create temporary key field and index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.deleteField('TempCompressFieldIndex',log=False)          # 'TempCompressFieldIndex' cannot already exist
            self.joinFields('TempCompressFieldIndex',newkeys,join=self.info['Delimit'])   
            index = self.index('TempCompressFieldIndex',force=True,fullreport=True)
            #self.deBug(rje.sortKeys(index))

            ### ~ [2] Compress Data into new data dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newdata = {}
            for nkey in index:      ### Each nkey should have a list of oldkeys to compress
                oldkeys = index[nkey]
                newdata[nkey] = {}
                if len(oldkeys) == 1: newdata[nkey] = self.data()[oldkeys[0]]   # No need to compress!
                elif best:  # Special compression using certain fields to determine whole entry
                    ok = 0
                    for i in range(1,len(oldkeys)):
                        for bfield in best:
                            try: usemin = rules[bfield].lower() == 'min'
                            except: usemin = False
                            try: (vi,vok) = (float(self.data()[oldkeys[i]][bfield]), float(self.data()[oldkeys[ok]][bfield]))
                            except: self.errorLog('%s: %s (%s) v %s (%s)' % (bfield,oldkeys[i],self.data()[oldkeys[i]][bfield], oldkeys[ok], self.data()[oldkeys[ok]][bfield])); raise
                            #self.deBug('%s %d %s: %s v %s (%d)' % (nkey,i,bfield,vi,vok,ok))
                            if vi == vok: continue                                      # Need to look at next field
                            elif (usemin and vi < vok) or (not usemin and vi > vok):     # Best so far
                                ok = i
                                break
                            else: break                                                 # Not good enough!
                    #self.deBug('=> %d' % ok)
                    for field in self.fields(): newdata[nkey][field] = self.data()[oldkeys[ok]][field]
                else:       # Normal compression, treating each field separately
                    ## ~ [2b] Deal with each field in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for key in newkeys: newdata[nkey][key] = self.data()[oldkeys[0]][key]
                    for field in rules:
                        if field in newkeys: continue   # No need to process this as will be unique
                        if field not in self.fields(): continue # Excess rules over fields!
                        values = []
                        numeric = not rules[field] in ['text','list','str']
                        for okey in oldkeys:
                            if not field in self.data()[okey] or self.data()[okey][field] == '': continue
                            if not numeric:  break
                            try:
                                self.data()[okey][field] / 3
                                values.append(self.data()[okey][field])
                            except:
                                try: values.append(string.atof(self.data()[okey][field]))
                                except: numeric = False
                        if not numeric:
                            values = []
                            for okey in oldkeys:
                                if not field in self.data()[okey] or self.data()[okey][field] == '': continue
                                if rules[field] == 'list': values += string.split(str(self.data()[okey][field]),joinchar)
                                else: values.append(self.data()[okey][field])
                        if not values: newdata[nkey][field] = ''
                        elif not numeric:
                            if rules[field] in ['min','max','list']: values.sort()
                            if rules[field] == 'max': newdata[nkey][field] = values[-1]
                            elif rules[field] == 'list':
                                #self.bugPrint('%s' % values)
                                sortdata = rje.sortUnique(values)
                                #self.debug('%s' % sortdata)
                                try: newdata[nkey][field] = string.join(sortdata,joinchar)
                                except:
                                    joindata = []
                                    for x in sortdata: joindata.append('%s' % x)
                                    newdata[nkey][field] = string.join(joindata,joinchar)
                            else: newdata[nkey][field] = values[0]
                        else:
                            if rules[field] == 'min': newdata[nkey][field] = min(values)
                            elif rules[field] == 'max': newdata[nkey][field] = max(values)
                            elif rules[field] == 'mean': newdata[nkey][field] = float(sum(values)) / len(values)
                            elif rules[field] == 'median' and rje.isOdd(len(values)): newdata[nkey][field] = values[len(values)/2]
                            elif rules[field] == 'median': newdata[nkey][field] = (values[len(values)/2]+values[(len(values)/2)-1])/2.0
                            elif rules[field] == 'geomean': newdata[nkey][field] = rje.geoMean(values)
                            elif rules[field] == 'sum': newdata[nkey][field] = sum(values)
            self.dict['Data'] = newdata
            self.list['Keys'] = newkeys
            self.dict['Index'] = {}
            self.deleteField('TempCompressFieldIndex')          # Cleanup after index
            if self.getBool('TupleKeys'): self.remakeKeys()     # Replace delimit join with tuples
            self.printLog('#CMPRSS','Compressed table "%s": %d fields; %s entries' % (self.name(),self.fieldNum(),rje.integerString(self.entryNum())))                        
        except: return self.errorLog('Major problem during Table.compress()')            
#########################################################################################################################
    def dropEntry(self,entry):    ### Drops specific entry from Table
        ekey = self.makeKey(entry)
        for ikey in self.dict['Index'].keys():
            try:
                self.dict['Index'][ikey][entry[ikey]].remove(ekey)
                if not self.dict['Index'][ikey][entry[ikey]]: self.dict['Index'][ikey].pop(entry[ikey])
            except: self.dict['Index'].pop(ikey)
        self.dict['Data'].pop(ekey)
#########################################################################################################################
    def dropEntries(self,filters,inverse=False,log=True,logtxt='',purelist=False,keylist=False):    ### Drops certain entries from Table
        '''
        Drops certain entries from Table.
        >> filters:list of str = Criteria that entries must meet
        >> inverse:bool [False] = Whether entries matching criteria should be exclusively kept rather than dropped.
        >> log:bool [True] = Whether to report reduction in log
        >> logtxt:str [''] = Text to add to log file line prior to count of rejected entries.
        >> purelist:bool [False] = Whether filters is actually a pure list of entries.
        >> keylist:bool [False] = Whether filters is actually a list of entry keys.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if type(filters) == str: filters = [filters]
            prex = self.entryNum()
            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if purelist:
                for entry in filters: self.dropEntry(entry)
            elif keylist: return self.dropEntries(self.entryList(filters),inverse,log,logtxt,purelist=True)
            else:
                statfilter = rje_scoring.setupStatFilter(self,self.list['Fields'],filters)
                #self.debug(statfilter)
                self.dict['Data'] = rje_scoring.statFilter(self,self.dict['Data'],statfilter,inverse,filtermissing=True)
            ### ~ [2] Report ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if prex != self.entryNum():
                if log and logtxt: self.printLog('#DROP','%s: %s %s entries reduced to %s entries' % (logtxt,rje.integerString(prex),self.info['Name'],rje.integerString(self.entryNum())))
                elif log: self.printLog('#DROP','%s %s entries reduced to %s entries' % (rje.integerString(prex),self.info['Name'],rje.integerString(self.entryNum())))
                if not purelist and not keylist: self.dict['Index'] = {}
        except: return self.log.errorLog('Major problem during Table.dropEntries()')            
#########################################################################################################################
    def dropEntriesDirect(self,field,values,inverse=False,log=True,force=False):    ### Drops certain entries from Table
        '''
        Drops certain entries from Table.
        >> filters:list of str = Criteria that entries must meet
        >> inverse:bool [False] = Whether entries matching criteria should be exclusively kept rather than dropped.
        >> log:bool [True] = Whether to report reduction in log
        >> force:bool [False] = Whether to force regeneration of index
        '''
        try:### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = self.entryNum()
            for ikey in rje.sortKeys(self.index(field,force=force)):
                if ikey in values and not inverse:
                    for dkey in self.index(field)[ikey][0:]: self.dropEntry(self.data(dkey)) #self.dict['Data'].pop(dkey)
                elif inverse and ikey not in values:
                    for dkey in self.index(field)[ikey][0:]: self.dropEntry(self.data(dkey)) #self.dict['Data'].pop(dkey)
            if prex != self.entryNum() and log: self.printLog('#DROP','%s %s entries reduced to %s entries on %s.' % (rje.integerString(prex),self.info['Name'],rje.integerString(self.entryNum()),field))
            #if prex != self.entryNum(): self.dict['Index'] = {}
        except TypeError:
            try:
                check = values[0:]
                check.sort()
                self.errorLog('Major problem during Table.dropEntriesDirect()'); raise
            except: self.dropEntriesDirect(field,[values],inverse,log)
        except: self.log.errorLog('Major problem during Table.dropEntriesDirect()'); raise
#########################################################################################################################
    def dropIndexEntries(self,index,values,inverse=False,log=True,force=False):    ### Drops certain entries from Table
        '''
        Drops certain entries from Table.
        >> index:str = Table index field
        >> values:list/str = list of values for index
        >> inverse:bool [False] = Whether entries matching criteria should be exclusively kept rather than dropped.
        >> log:bool [True] = Whether to report reduction in log
        >> force:bool [False] = Whether to force regeneration of index
        '''
        try:### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = self.entryNum()
            if type(values) == str: values = [values]
            entries = []
            for value in values: entries += self.indexEntries(index,value)  # Returns list of entries from index & value
            for entry in self.entries():
                if entry in entries:
                    if not inverse: self.dropEntry(entry)
                elif inverse: self.dropEntry(entry)
            if prex != self.entryNum() and log: self.printLog('#DROP','%s %s entries reduced to %s entries on %s.' % (rje.integerString(prex),self.info['Name'],rje.integerString(self.entryNum()),index))
        except: self.log.errorLog('Major problem during Table.dropIndexEntries()'); raise
#########################################################################################################################
    def dropFields(self,fields,inverse=False,log=None): ### Drops certain fields from Table
        '''
        Drops certain entries from Table.
        >> fields:list of str = Fields to delete
        >> inverse:bool [False] = Whether fields should be exclusively kept rather than dropped.
        >> log:bool [True] = Whether to report reduction in log
        '''
        try:### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = self.fieldNum()
            if inverse:
                for field in self.fields():
                    if field not in fields: self.deleteField(field,log)
            else:
                for field in fields: self.deleteField(field,log)
            if log == None: log = True
            if log: self.printLog('#DROP','%s %s fields reduced to %s fields' % (rje.integerString(prex),self.info['Name'],rje.integerString(self.fieldNum())))
        except: return self.log.errorLog('Major problem during Table.dropEntries()')            
#########################################################################################################################
    def keepFields(self,fields,log=None):   ### Inverse dropFields() - only keeps certain fields in Table
        '''Inverse dropFields() - only keeps certain fields in Table.'''
        return self.dropFields(fields,True,log)
#########################################################################################################################
    def setFields(self,fields,log=None):   ### Inverse dropFields() - only keeps certain fields in Table
        '''Inverse dropFields() - only keeps certain fields in Table. Adds missing fields and reorders.'''
        self.dropFields(fields,True,log)
        for field in fields:
            if field not in self.fields(): self.addField(field,log=log)
        self.list['Fields'] = fields[0:]
        self.printLog('#FIELD','Set %s fields: %s' % (self.name(),string.join(fields[0:],', ')))
#########################################################################################################################
    def fillBlanks(self,blank='',fields=[],fillempty=False,prog=True,log=True):      ### Fills in missing fields in entries with blanks
        '''Fills in missing fields in entries with blanks.'''
        if not fields: fields = self.fields()
        (ex,etot) = (0.0,self.entryNum())
        bx = 0
        for entry in self.entries():
            if prog: self.bugProg('\r#FILL','Filling %s blanks: %.2f%%' % (self.info['Name'],ex/etot)); ex += 100.0
            for field in fields:
                if field not in entry: entry[field] = blank; bx += 1
                elif fillempty and not string.join(string.split('%s' % entry[field])): entry[field] = blank; bx += 1
        if prog or log: self.bugLog('\r#FILL','Filled %s %s blank values with "%s"' % (rje.iStr(bx),self.info['Name'],blank),log=log)
#########################################################################################################################
    ### <7> ### Data Reformatting Methods                                                                               #
#########################################################################################################################
    def dataFormat(self,reformat={},skipblank=True):  ### Reformats fields and update self.dict['DataTypes']
        '''
        Reformats fields using (and updating) self.dict['DataTypes'].
        >> reformat:dict = Dictionary of {field:data type} to populate self.dict['DataTypes'].
        - str/int/num/bool/estr
        >> skipblank:bool [True] = skip empty ('') entries
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setBool({'Formatted':False})
            intwarn = []; rekey = False
            for field in reformat:
                self.dict['DataTypes'][field] = reformat[field].lower()
                if field in self.keys(): rekey = True
            ### ~ [2] Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (ex,etot,fx) = (0.0,self.entryNum(),0)
            for oldkey in self.dataKeys():
                entry = self.data()[oldkey]
                self.bugProg('\r#FORMAT','Reformatting %s: %.2f%%' % (self.info['Name'],ex/etot)); ex += 100.0
                for field in self.dict['DataTypes']:
                    if field not in entry: continue
                    try:
                        if entry[field] == '':
                            if not skipblank: entry[field] = {'str':'','int':0,'num':0.0,'flo':0.0,'boo':False}[self.dict['DataTypes'][field][:3]]
                            continue                            
                        if self.dict['DataTypes'][field][:3] == 'str': entry[field] = str(entry[field])
                        if self.dict['DataTypes'][field][:3] in ['int']:
                            try: entry[field] = int(entry[field])
                            except:
                                entry[field] = int(float(entry[field]))
                                if field not in intwarn:
                                    self.printLog('\r#FWARN','Integer field "%s" might have contained float values. ' % field)
                                    intwarn.append(field)
                        if self.dict['DataTypes'][field][:3] in ['num','flo','est']: entry[field] = float(entry[field])
                        if self.dict['DataTypes'][field][:4] == 'estr': entry[field] = rje.expectString(entry[field])
                        if self.dict['DataTypes'][field][:4] == 'bool':
                            if entry[field]:
                                if str(entry[field]).lower() in ['0','false','f','no','n']: entry[field] = False
                                else: entry[field] = True
                            else: entry[field] = False
                    except:
                        fx += 1
                        self.deBug('%s %s - %s?' % (field,entry[field],self.dict['DataTypes'][field]))
                if rekey:
                    newkey = self.makeKey(entry)
                    self.dict['Data'][newkey] = self.dict['Data'].pop(oldkey)
            self.bugLog('\r#FORMAT','Reformatting %s complete: %s format errors' % (self.info['Name'],fx))
            if not fx: self.setBool({'Formatted':True})
        except: return self.log.errorLog('Major problem during Table.dataFormat()')            
#########################################################################################################################
### End of SECTION III: Table Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: print '\n\n *** No standalone functionality! *** \n\n'
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
