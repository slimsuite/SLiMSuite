#!/usr/bin/python

# See below for name and description
# Copyright (C) 2012 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_xref
Description:  Generic identifier cross-referencing 
Version:      1.8.0
Last Edit:    28/06/15
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is primarily for use with other programs/modules to handle database cross-referencing. Initially, it is 
    designed to be able to use the output from rje_genemap but might eventually replace this module by generating the 
    xref data from raw data in the first place. Although it is designed and documented with aliases in mind, e.g. mapping
    UniProt IDs to HGNC, it can also be used for more general mapping of, for example, Pfam domains and genes.
    
    This module is designed to work with a main xrefdata table (db table XRef) that contains 1:1:1 values for each
    identifier. Where there are multiple values for a given field, these will be joined with '|' characters (or
    splitchar=X) and likewise split on '|' to generate aliases from the xrefdata table. (These will be sorted if
    sortxref=T. This may slow the program down a little.) Split fields will have whitespace removed. To avoid a field
    being split, use splitskip=LIST.

    If multiple xrefdata files are provided, these will be combined into a single table. The first file becomes the
    master XRef table and any additional input will be added to this table provided that (a) it has a KeyID or AltKeys
    field, or (b) it can be mapped onto an existing KeyID via a given MapField.

    NOTE. If newheaders=LIST is used, it will apply to the FIRST xrefdata table only and the new fields will be used for
    all subsequent processing. All other options should therefore use these new headers. If multiple files need new
    headers, the program may need to be run several times before combining them. If newheaders replaces KeyID then KeyID
    will be updated.

    If fullmap=T then ALL MapFields will be used for mapping, even if this produces multiple mappings. Otherwise, the
    first successful mapping will be used. The KeyID and AltKeys fields are automatically added to the
    front of any MapField list (and thus will be used first for mapping). If an ID List is provided (idlist=LIST) then
    these final XRef data are restricted to the KeyIDs in this list, once all mapping has been done.

    If FileXRef file is given (filexref=FILE) then the XRef data will be mapped onto this file (via MapFields) and
    subsequent output will correspond to the mappings and IDs in this combined file. XRef fields to be added to the file
    can be limited with xrefs=LIST. Sorted (unique) lists of mapped files can be produced by xreflist=LIST. If no
    filexref file is given then output is for the entire XRef data table. Note that the default is to produce this table.

    Some database identifiers are prefixed with the database and a colon. e.g. Entrez gene 10840 may be written
    ENTREZ:10840. This will be recognised and standardised by either removing the prefix or, if the field is in
    dbprefix=LIST, by enforcing the DB:ID format. The database is always taken from the XRef field, so this must match,
    although the case will be switched to uppercase unless in keepcase=LIST.

Commandline:
    ### ~ Input/Field Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    xrefdata=FILES  : List of files with delimited data of identifier cross-referencing (wildcards allowed) []
    newheaders=LIST : List of new Field headers for XRefData (will replace old - must be complete) []
    keepcase=LIST   : Any fields matching keepcase will retain mixed case, otherwise be converted to upper case ['Desc','Description','Name']
    dbprefix=LIST   : List of fields that should have the field added to the ID as a prefix, e.g. HGNC:0001 []
    stripvar=CDICT  : Remove variants using Field:Char list, e.g. Uniprot:-,GenPept:. []
    compress=LIST   : Compress listed fields into lists (using splitchar) to allow 1:many mapping in xrefdata. []
    splitchar=X     : Character on which to split fields for multiple alias processing ['|']
    splitcsv=T/F    : Whether to also split fields based on comma separation [True]
    splitskip=LIST  : List of fields to bypass for field splitting ['Desc','Description','Name']
    sortxref=T/F    : Whether to sort multiple xref data alphabetically [True]
    keyid=X         : Key field header to be used in main Data dictionary - aliases map to this ['Gene']
    comments=LIST   : List of comment line prefixes marking lines to ignore (throughout file) ['//','%']
    xreformat=T/F   : Whether to apply field reformatting to input xrefdata (True) or just xrefs to map (False) [False]
    yeastxref=T/F   : First xrefdata file is a yeast.txt file to convert. (http://www.uniprot.org/docs/yeast.txt)

    ### ~ XRef/Processing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    altkeys=LIST    : Alternative fields to look for in Alias Data ['Symbol','HGNC symbol']
    onetomany=T/F   : Whether to keep potential one-to-many altkeys IDs [False]
    mapfields=X     : Fields to be used for Alias mapping plus KeyID. (Must be in XRef). []
    maptomany=T/F   : Whether to keep potential one-to-many mapfield IDs [True]
    fullmap=T/F     : Whether to map onto ALL map fields or stop at first hit [False]
    uniquexref=T/F  : Whether to restrict analysis to unique XRef IDs [False]
    mapxref=LIST    : List of identifiers to map to KeyIDs using mapfields []
    filexref=FILE   : File to XRef and expand with xrefs before re-saving []
    badid=LIST      : List of XRef IDs to ignore ['!FAILED!','None','N/A','-']
    aliases=LIST    : Combine XRef fields into single 'Aliases' field (and remove KeyID if found)

    ### ~ Join Method Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    join=LIST       : Run in join mode for list of FILE:key1|...|keyN:JoinField []
    naturaljoin=T/F : Whether to only output entries that join to all tables [False]

    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : Basefile for output files [Default: filexref or first xrefdata input file w/o path]
    savexref=T/F    : Save the xrefdata table (*.xref.tdt) following compilation of data [True]
    idlist=LIST     : Subset of key IDs to map onto. (All if blank) []
    xrefs=LIST      : List of XRef (or join) fields to keep (blank/* for all) []
    xreflist=LIST   : List of XRef fields to output as sorted (unique) lists (*.*.txt) (* for all) []
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Added xfrom and xto fields and xMap() function for mapping from one ID set to another.
    # 1.1 - Added output of ID lists to text files. Major reworking. Tested with HPRD and HGNC.
    # 1.2 - Added join=LIST Run in join mode for list of FILE:key1|...|keyN:JoinField [] and naturaljoin=T/F.
    # 1.3.0 - Added compress=LIST to handle 1:many input data. []
    # 1.3.1 - Fixed xref list bug.
    # 1.4.0 - Added optional Mapping dictionary for speeding up recurring mapping (should avoid if memsaver=F).
    # 1.5.0 - Added stripvar=CDICT removal of variants using Field:Char list, e.g. Uniprot:-,GenPept:. []
    # 1.6.0 - Added mapxref=LIST List of identifiers to map to KeyIDs using mapfields []
    # 1.7.0 - Added comments=LIST ist of comment line prefixes marking lines to ignore (throughout file) ['//','%']
    # 1.7.1 - Added xreformat=T/F : Whether to apply field reformatting to input xrefdata (True) or just xrefs to map (False) [False]
    # 1.8.0 - Added recognition and parsing of yeast.txt XRef file from Uniprot (http://www.uniprot.org/docs/yeast.txt)
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add methods for returning aliases or mapping.
    # [X] : Add default xto and xfrom fields for mapping.
    # [Y] : Replace attributes
    # [Y] : Revise the module loading method
    # [Y] : Add the addAlias method. (Expand/revise db index as it goes along.)
    # [Y] : Add revised mapping method using expanded rje_db index method with splitchar=X and splitfield=True
    # [Y] : Add output of ID lists
    # [X] : Add parsing of HGNC data using rje_genemap prior to run?
    # [?] : savetab=T/F     : Save individual alias/mapping tables (*.*.tdt) following compilation of data [False]
    # [Y] : Add some kind of integrity check for multiple/cyclical mapping?
    # [+] : Deal with db:ID mappings in both add Alias and mapping methods. Check both but make a switch for output.
    # [Y] : |-- Remember to check whether field is in Upper Case. (Should use field header for db)
    # [+] : Add uniquexref=T/F to restrict to 1:1
    # [Y] : Add some checks of synonyms etc. (1:1 vs 1:many)
    # [+] : Add options to combine fields, e.g. combine altkeys into "Aliases"
    # [ ] : Still need to do a thorough testing of all the [+] items.
    # [ ] : Add idlist and keyid function to join easily enough: filter on index(keyid) following join.
    # [ ] : Add newheaders=LIST function to join: initially map to xrefs and then reorder ['AutoID']+newheaders.
    # [ ] : Add reading of single join element and only using newheaders and idlist/keyid filter. (KeyID from jfield?)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('RJE_XRef', '1.8.0', 'June 2015', '2014')
    description = 'Generic identifier cross-referencing'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmdhelp = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmdhelp > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)   # Reads arguments and load defaults from program.ini
        out = rje.Out(cmd_list=cmd_list)                    # Sets up Out object for controlling output to screen
        out.verbose(2,2,cmd_list,1)                         # Prints full commandlist if verbosity >= 2 
        out.printIntro(info)                                # Prints intro text using details from Info object
        cmd_list = cmdHelp(info,out,cmd_list)               # Shows commands (help) and/or adds commands from user
        log = rje.setLog(info,out,cmd_list)                 # Sets up Log object for controlling log file output
        return (info,out,log,cmd_list)                      # Returns objects for use in program
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Problem during initial setup.'; raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: XRef Class                                                                                              #
#########################################################################################################################
class XRef(rje_obj.RJE_Object):     
    '''
    XRef Class. Author: Rich Edwards (2013).

    Str:str
    - FileXRef = File to XRef and expand with xrefs before re-saving []
    - KeyID = Key field header to be used in main Data dictionary - aliases map to this ['Gene']
    - SplitChar = Character on which to split fields for multiple alias processing ['|']

    Bool:boolean
    - FullMap = Whether to map onto ALL map fields or stop at first hit [False]
    - MapToMany = Whether to keep potential one-to-many mapfield IDs [True]
    - NaturalJoin=T/F : Whether to only output entries that join to all tables [False]
    - OneToMany = Whether to keep potential one-to-many altkeys IDs [False]
    - SaveXRef = Save the xrefdata table (*.xref.tdt) following compilation of data [True]
    - SortXRef = Whether to sort multiple xref data alphabetically [True]
    - SplitCSV = Whether to also split fields based on comma separation [True]
    - UniqueXRef = Whether to restrict analysis to unique XRef IDs [False]
    - XReformat=T/F   : Whether to apply field reformatting to input xrefdata (True) or just xrefs to map (False) [False]
    - YeastXRef=T/F   : Whether the first xrefdata is a yeast.txt file. (http://www.uniprot.org/docs/yeast.txt)

    Int:integer

    Num:float

    List:list
    - Aliases = Combine XRef fields into single 'Aliases' field (and remove KeyID if found)
    - AltKeys = Alternative fields to look for in Alias Data ['Symbol','HGNC symbol']
    - BadID = List of XRef IDs to ignore ['!FAILED!','None','NA','N/A']
    - Comments = List of comment line prefixes marking lines to ignore (throughout file) ['//','%']
    - Compress = Compress listed fields into lists (using splitchar) to allow 1:many mapping in xrefdata. []
    - DBPrefix = List of fields that should have the field added to the ID as a prefix, e.g. HGNC:0001 []
    - IDList = Subset of key IDs to map onto. (All if blank) []
    - Join=LIST       : Run in join mode for list of FILE:key1|...|keyN;JoinField []
    - KeepCase = Any fields matching keepcase will retain mixed case, otherwise be converted to upper case ['Desc','Description']
    - MapFields = Fields to be used for Alias mapping plus KeyID. (Must be in XRef). []
    - MapXRef=LIST    : List of identifiers to map to KeyIDs using mapfields []
    - NewHeaders = List of new Field headers for XRefData (will replace old - must be complete) []
    - SplitSkip = List of fields to bypass for field splitting ['Desc','Description']
    - XRefData = Delimited data of identifier cross-referencing []
    - XRefs = List of XRef fields to output (* for all) [*]
    - XRefList = List of XRef fields to output as sorted (unique) lists (*.*.txt) (* for all) []

    Dict:dictionary
    - Mapping = Optional dictionary of {XField:{xref:mapped id}}
    - StripVar = Remove variants using Field:Char list, e.g. Uniprot:-,GenPept:. []

    Obj:RJE_Objects
    - DB = rje_db.Database object - containing XRef and alias tables. (Alias tables named by aliasdata fields)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['FileXRef','KeyID','SplitChar']
        self.boollist = ['FullMap','MapToMany','NaturalJoin','OneToMany','SaveXRef','SortXRef','SplitCSV','UniqueXRef','XReformat','YeastXRef']
        self.intlist = []
        self.numlist = []
        self.listlist = ['Aliases','AltKeys','BadID','Comments','Compress','DBPrefix','IDList','Join','KeepCase','MapFields','MapXRef','NewHeaders','SplitSkip','XRefs','XRefList','XRefData']
        self.dictlist = ['Mapping','StripVar']
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'KeyID':'Gene','SplitChar':'|'})
        self.setBool({'FullMap':False,'MapToMany':True,'OneToMany':False,'SortXRef':True,'SplitCSV':True,'SaveXRef':True,'UniqueXRef':False,'XReformat':False,'YeastXRef':False})
        self.setInt({})
        self.setNum({})
        self.list['AltKeys'] = ['Symbol','HGNC symbol']
        self.list['SplitSkip'] = ['Desc','Description','Name']
        self.list['KeepCase'] = ['Desc','Description','Name']
        self.list['XRefList'] = []
        self.list['BadID'] = ['!FAILED!','None','N/A','-']
        self.list['Comments'] = ['//','%']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['KeyID','SplitChar'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['FileXRef'])  # String representing file path
                self._cmdReadList(cmd,'bool',['FullMap','MapToMany','NaturalJoin','OneToMany','SaveXRef','SortXRef','SplitCSV','UniqueXRef','XReformat','YeastXRef'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Aliases','AltKeys','BadID','Comments','Compress','DBPrefix','IDList','Join','KeepCase','MapFields','MapXRef','NewHeaders','SplitSkip','XRefs','XRefList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['XRefData']) # List of files using wildcards and glob
                self._cmdReadList(cmd,'cdict',['StripVar']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['Join']: return self.join()
            if not self.setup(): return False
            keyid = self.getStr('KeyID')
            xdb = self.db('xref')

            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Add XRef data to other file and/or reduce XRef fields ~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.exists(self.getStr('FileXRef')): xdb = self.fileXRef(self.getStr('FileXRef')) # Combine
            elif self.list['MapXRef']: xdb = self.fileXRef(mapxref=self.list['MapXRef'])
            elif self.list['XRefs']:
                xdb.keepFields([keyid]+self.list['XRefs'])
                newhead = []
                for field in [keyid]+self.list['XRefs']:
                    if field in newhead: continue
                    elif field in xdb.fields(): newhead.append(field)
                xdb.list['Fields'] = newhead[0:]
                self.printLog('#XLIST','Restrict output to %s.' % string.join(newhead,', '))
            ## ~ [2b] ~ Reduce to IDList will happen after all other mapping done ~~~~~~~~~~~~~~~~~ ##
            if keyid not in self.list['KeepCase']:
                self.list['IDList'] = rje.listUpper(self.list['IDList'])
                self.printLog('#CASE','%s "%s" identifiers converted to uppercase.' % (rje.iLen(self.list['IDList']),keyid))
            if self.list['IDList']: xdb.dropEntriesDirect(keyid,self.list['IDList'],inverse=True,log=True)
            if not xdb.entryNum():
                self.warnLog('No XRef entries retained.')
                return False

            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SaveXRef'): self.dict['Output']['tab'] = xdb.saveToFile()
            if '*' in self.list['XRefList']: self.list['XRefList'] = xdb.fields()[0:]
            for field in self.list['XRefList']:
                xfile = self.xrefList(field)
                if field.lower() in self.dict['Output']: self.dict['Output'][field.lower()] = xfile
            ## ~ [3a] Additional REST Outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('Rest'):
                self.dict['Output']['failed'] = []
                mdb = self.db().addEmptyTable('mapped',['MapID','XRef'],['MapID'],log=True)
                xfield = self.getStr('KeyID')
                for mkey in self.dict['Mapping'][xfield]:
                    mxref = self.dict['Mapping'][xfield][mkey]
                    if mxref: mdb.addEntry({'MapID':mkey,'XRef':mxref})
                    else: self.dict['Output']['failed'].append(mkey)
                self.dict['Output']['mapped'] = mdb.saveToFile()
                self.dict['Output']['failed'].sort()
                if self.dict['Output']['failed']: self.dict['Output']['failed'] = string.join(self.dict['Output']['failed'],'\n')
                else: self.dict['Output']['failed'] = 'None'
            return True
        except:
            self.errorLog(rje_obj.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if '*' in self.list['XRefs']: self.list['XRefs'] = []
            db = self.db()
            keyid = self.getStr('KeyID')
            ## ~ [0a] FileXRef file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.exists(self.getStr('FileXRef'),True):
                if self.baseFile().lower() in ['','none']: self.setBaseFile(rje.baseFile(self.getStr('FileXRef'),strip_path=True))
            else:
                self.printLog('#ERR','IOError: FileXRef file "%s" not found!' % self.getStr('FileXRef'))
                raise IOError
            ## ~ [0b] XRefData files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xfiles = []
            for xfile in self.list['XRefData']:
                if not rje.checkForFile(xfile):
                    if rje.yesNo('Cannot find %s. Quit %s?' % (xfile,self.prog())):
                        self.printLog('#ERR','IOError: XRefData file "%s" not found!' % xfile)
                        raise IOError
                else: xfiles.append(xfile)
            if not xfiles: return False     # Add functions for filexref alone?
            xrefdata = xfiles.pop(0)
            ## ~ [0b] BaseFile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.baseFile().lower() in ['','none']: self.setBaseFile(rje.baseFile(xrefdata,strip_path=True))
            #self.debug(db.baseFile())
            ## ~ [0c] Yeast.txt file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('YeastXRef'): xrefdata = self.yeastXRef(xrefdata)

            ### ~ [1] Load XRefData ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Apply NewHeader and compressed field lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['NewHeaders'] and keyid not in self.list['NewHeaders']:
                raise ValueError('KeyID "%s" not found in NewHeaders: %s' % (keyid,self.list['NewHeaders']))
            if self.list['Compress']: loadkey = '#'
            else: loadkey = [keyid]
            if self.list['NewHeaders']:
                xdb = db.addTable(filename=xrefdata,mainkeys=loadkey,name='xref',headers=self.list['NewHeaders'],ignore=self.list['Comments'])
            else:
                xdb = db.addTable(filename=xrefdata,mainkeys=loadkey,name='xref',ignore=self.list['Comments'])
            if self.list['Compress']:
                if '*' in self.list['Compress']: self.list['Compress'] = xdb.fields(); self.list['Compress'].remove(keyid)
                rules = {}
                for field in self.list['Compress']: rules[field] = 'list'
                xdb.compress([keyid],rules,default='str',joinchar=self.getStr('SplitChar'))
                xdb.dropField('#')
            if not xdb: raise ValueError
            if self.list['NewHeaders'] and self.list['NewHeaders'] != xdb.fields():
                if len(self.list['NewHeaders']) != xdb.fieldNum():
                    self.errorLog('%d New headers given (newheaders=LIST) but %d fields.' % (len(self.list['NewHeaders']),xdb.fieldNum()),printerror=False)
                    raise ValueError
                for i in range(xdb.fieldNum()):
                    self.printLog('#FIELD','Loaded field "%s" renamed "%s".' % (xdb.fields()[i],self.list['NewHeaders'][i]))
                    xdb.renameField(xdb.fields()[i],self.list['NewHeaders'][i])
                keyid = xdb.keys()[0]
                if self.str['KeyID'] != keyid:
                    self.printLog('#KEYID','KeyID changed from "%s" to "%s".' % (self.str['KeyID'],keyid))
                    self.str['KeyID'] = keyid
                else: self.printLog('#KEYID','KeyID = "%s".' % (keyid))
            self.list['SplitSkip'].append(keyid)
            ## ~ [1b] Reformat entry fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for field in xdb.fields():
                if self.getBool('XReformat'):
                    for entry in xdb.entries(): self.reformatEntryField(entry,field)
                    rftxt = []
                    if field in self.list['DBPrefix']: rftxt.append('DBPrefix')
                    if field in self.dict['StripVar']: rftxt.append('StripVar')
                    if field not in self.list['KeepCase']: rftxt.append('Case')
                    if rftxt: rftxt = '%s ' % string.join(rftxt,'/')
                    else: rftxt = ''
                    self.printLog('#FORMAT','Applied %sreformatting to "%s".' % (rftxt,field))
                    if field == keyid: xdb.remakeKeys()
                if field in self.list['AltKeys'] + self.list['MapFields']:
                    if field in self.list['SplitSkip']: xdb.index(field)
                    else: xdb.index(field,splitchar=self.getStr('SplitChar'))
            ## ~ [1c] Check synonym mismatch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.checkOneToMany(xrefdata): return False

            ### ~ [2] Add additional xref data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for xfile in xfiles:
                tdb = db.addTable(xfile,mainkeys='#',lists=False,name='AliasTemp',expect=True,replace=True)
                for field in self.list['AltKeys']:
                    if field in tdb.fields(): tdb.renameField(field,keyid); break
                for field in tdb.fields()[0:]:
                    if field == keyid: continue
                    if self.list['XRefs'] and field not in xdb.fields() + self.list['XRefs'] + tdb.keys() + self.list['Aliases']:
                        tdb.dropField(field)
                        continue
                    if self.getBool('XReformat'):
                        for entry in tdb.entries(): self.reformatEntryField(entry,field)
                        rftxt = []
                        if field in self.list['DBPrefix']: rftxt.append('DBPrefix')
                        if field in self.dict['StripVar']: rftxt.append('StripVar')
                        if field not in self.list['KeepCase']: rftxt.append('Case')
                        if rftxt: rftxt = '%s ' % string.join(rftxt,'/')
                        else: rftxt = ''
                        self.printLog('#FORMAT','Applied %sreformatting to "%s".' % (rftxt,field))
                #for field in self.list['MapFields']:
                #    if field in self.list['SplitSkip']: tdb.index(field)
                #    else: tdb.index(field,splitchar=self.getStr('SplitChar'))
                for entry in tdb.entries(): self.addXRef(entry)
                db.deleteTable('AliasTemp')
            ## ~ [2a] Check synonym mismatch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.checkOneToMany(xfile): return False  #!# Should not be necessary so recode if it happens!
            ## ~ [2b] Check XRefs against NewHeaders ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for xref in self.list['XRefs'][0:]:
                if xref not in xdb.fields():
                    self.warnLog('XRef field "%s" not found in headers: removed from XRef list' % xref,quitchoice=True)
                    self.list['XRefs'].remove(xref)

            ### ~ [3] Combine Aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            splitchar = self.getStr('SplitChar')
            if self.list['Aliases']:
                for field in self.list['Aliases'][0:]:
                    if field not in xdb.fields():
                        self.printLog('#ALIAS','Aliases field "%s" not found!' % field)
                        self.list['Aliases'].remove(field)
                for entry in xdb.entries():
                    if 'Aliases' in entry and entry['Aliases']: entry['Aliases'] = string.split(entry['Aliases'],splitchar)
                    else: entry['Aliases'] = []
                    for field in self.list['Aliases']:
                        for alias in string.split(entry[field],splitchar):
                            if alias and alias not in entry['Aliases']: entry['Aliases'].append(alias)
                    while '' in entry['Aliases']: entry['Aliases'].remove('')
                    entry['Aliases'].sort()
                    entry['Aliases'] = string.join(entry['Aliases'],splitchar)
                if 'Aliases' not in xdb.fields(): xdb.list['Fields'].append('Aliases')
                if 'Aliases' not in self.list['XRefs']: self.list['XRefs'].append('Aliases')
                for field in self.list['Aliases'][0:]:
                    if self.list['XRefs'] and field not in self.list['XRefs'] + xdb.keys(): xdb.dropField(field)

            self.restSetup()
            self.printLog('#SETUP','XRef setup complete.')
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def yeastXRef(self,xrefdata): ### Convert a yeast.txt file to *.tdt file
        '''Convert a yeast.txt file to *.tdt file.'''
        try:### ~ [0]  Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#YEAST',xrefdata)
            YTXT = open(xrefdata,'r')
            ytabfile = '%s.tdt' % self.basefile()
            yfields = ['Gene','Synonyms','OLN','Uniprot','UniprotID','SGD','Size','3D','CH']
            ydb = self.db().addEmptyTable('yeast',yfields,['Gene'],log=True)
            ### ~ [1]  Process Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ystart = [0]; yi = 1
            yline = YTXT.readline(); self.debug(yline)
            while not yline.startswith('_'): yline = YTXT.readline()
            while yi < len(yline):
                if yline[yi] == '_' and yline[yi-1] != '_': ystart.append(yi)
                yi += 1
            ystart.append(len(yline))   # List of field position starts
            yline = YTXT.readline()
            while not yline.startswith('_'): yline = YTXT.readline()
            ### ~ [2]  Process Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            yline = YTXT.readline(); self.debug(yline)
            while yline and not yline.startswith('-'):
                yline = string.replace(yline,';',' ')
                yentry = {}
                for yx in range(len(ystart)-1):
                    data = string.split(yline[ystart[yx]:ystart[yx+1]])
                    self.bugPrint(data)
                    if yx:
                        ykey = yfields[yx+1]
                        if data: yentry[ykey] = data[0]
                        else: yentry[ykey] = ''
                    else:
                        if data:
                            yentry['Gene'] = data[0]
                            yentry['Synonyms'] = string.join(data[1:],self.getStr('SplitChar'))
                        else: break
                self.debug(yentry)
                if yentry: ydb.addEntry(yentry)
                yline = YTXT.readline()
            ydb.saveToFile(ytabfile)
            return ytabfile
        except: self.errorLog('Problem during %s yeastXRef.' % self); return False  # Setup failed
#########################################################################################################################
    def checkOneToMany(self,filestr):   ### Checks and optionally purges one-to-many mappings
        '''Checks and optionally purges one-to-many mappings.'''
        try:### ~ [0]  Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('OneToMany') and self.getBool('MapToMany') and (self.i() < 1 or not rje.yesNo('Check alias/mapping amibiguities?',default='N')):
                return True
            xdb = self.db('xref')
            ambig = {}  # Dictionary of IDs with ambiguity
            mwx = 0     # Total count of ambiguous IDs
            prevfields = []

            ### ~ [1] Check AltKeys for ambiguity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for field in self.list['AltKeys']:
                if field not in xdb.fields(): continue
                ambig[field] = []
                ix = 0.0; itot = len(xdb.index(field))
                for fid in rje.sortKeys(xdb.index(field)):
                    self.progLog('\r#FID','Checking "%s" for ambiguity: %.2f%%' % (field,ix/itot)); ix += 100.0
                    if not fid or fid in self.list['BadID']: continue
                    if fid in xdb.dataKeys() and xdb.index(field)[fid] != [fid]:
                        self.warnLog('"%s" alias %s is KeyID but points to %s.' % (field,fid,string.join(xdb.index(field)[fid])),'wrongkeyid',suppress=self.i()>0)
                        mwx += 1; ambig[field].append(fid)
                    elif len(xdb.index(field)[fid]) > 1:
                        self.warnLog('"%s" alias %s points to multiple KeyID: %s.' % (field,fid,string.join(xdb.index(field)[fid])),'ambiguity',suppress=self.i()>0)
                        mwx += 1; ambig[field].append(fid)
                    else:
                        for afield in prevfields:
                            if fid in xdb.index(afield) and xdb.index(afield)[fid] != xdb.index(field)[fid]:
                                self.warnLog('"%s" alias %s found in "%s" but points to different Key ID: %s vs %s.' % (field,fid,afield,xdb.index(afield)[fid][0],xdb.index(field)[fid][0]),'ambiguity',suppress=self.i()>0)
                                mwx += 1; ambig[field].append(fid)
                                if fid not in ambig[field]: ambig[afield].append(fid)
                self.printLog('\r#FID','%s "%s" alias ambiguities found (from %s).' % (rje.iLen(ambig[field]),field,rje.iStr(itot)))
                prevfields.append(field)
            ## ~ [1b] Purge/Report ambiguities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if mwx and (not self.getBool('OneToMany') or (self.i() > 0 and rje.yesNo('Remove %s alias amibiguities?' % rje.iStr(mwx),default='N'))):
                for field in ambig:
                    for fid in ambig[field]:
                        for kid in xdb.index(field).pop(fid):
                            if field in self.list['SplitSkip']: xdb.data(kid)[field] = ''
                            else:
                                fdata = string.split(xdb.data(kid)[field],self.getStr('SplitChar'))
                                fdata.remove(fid)
                                xdb.data(kid)[field] = string.join(fdata,self.getStr('SplitChar'))
                    if ambig[field]: self.printLog('#PURGE','%s ambiguous "%s" aliases purged' % (rje.iLen(ambig[field]),field))
            elif mwx: self.warnLog('%s alias ambiguities found!' % rje.iStr(mwx),quitchoice=True)
            else: self.printLog('#CHECK','No alias ambiguities found in %s' % filestr)
            #self.deBug('...')

            ### ~ [2] Check MapFields for ambiguity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mwx = 0
            if not self.getBool('MapToMany') or (self.i() > 0 and rje.yesNo('Check mapping amibiguities?',default='N')):
                for field in self.list['MapFields']:
                    if field not in xdb.fields() or field in ambig: continue
                    ambig[field] = []
                    ix = 0.0; itot = len(xdb.index(field))
                    for fid in rje.sortKeys(xdb.index(field)):
                        self.progLog('\r#FID','Checking "%s" for ambiguity: %.2f%%' % (field,ix/itot)); ix += 100.0
                        if not fid or fid in self.list['BadID']: continue
                        if len(xdb.index(field)[fid]) > 1:
                            self.warnLog('"%s" mapping ID %s points to multiple KeyID: %s.' % (field,fid,string.join(xdb.index(field)[fid])),'multimap',suppress=self.i()>0)
                            mwx += 1; ambig[field].append(fid)
                    self.printLog('\r#FID','%s "%s" mapping ambiguities found (from %s).' % (rje.iLen(ambig[field]),field,rje.iStr(itot)))
            ## ~ [1b] Purge/Report ambiguities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if mwx and (not self.getBool('MapToMany') or (self.i() > 0 and rje.yesNo('Remove %s mapping amibiguities?' % rje.iStr(mwx),default='N'))):
                for field in self.list['MapFields']:
                    if field not in xdb.fields() or field in ambig: continue
                    for fid in ambig[field]:
                        for kid in xdb.index(field).pop(fid):
                            if field in self.list['SplitSkip']: xdb.data(kid)[field] = ''
                            else:
                                fdata = string.split(xdb.data(kid)[field],self.getStr('SplitChar'))
                                fdata.remove(fid)
                                xdb.data(kid)[field] = string.join(fdata,self.getStr('SplitChar'))
                    if ambig[field]: self.printLog('#PURGE','%s ambiguous "%s" mapping IDs purged' % (rje.iLen(ambig[field]),field))
            elif mwx: self.warnLog('%s mapping ambiguities found!' % rje.iStr(mwx),quitchoice=True)
            else: self.printLog('#CHECK','No mapping ambiguities found in %s' % filestr)

            return True
        except: self.errorLog('Problem during %s.checkOneToMany().' % self); return False  # Setup failed
#########################################################################################################################
    def reformatEntryField(self,entry,field,inplace=True):  ### Reformats entry field *in place* using object parameter settings.
        '''Reformats entry using object field parameter settings.'''
        try:###  ~ [1] Split data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fupper = field.upper()
            if field not in self.list['KeepCase']: edata = entry[field].upper()
            else: edata = entry[field]
            if field in self.list['SplitSkip']: edata = [edata]
            else:
                if self.getBool('SplitCSV'): edata = string.replace(edata,',',self.getStr('SplitChar'))
                edata = string.split(edata,self.getStr('SplitChar'))
                edata = string.split(string.join(edata))    # Removes additional whitespace (even if splitchar = " ")
            ### ~ [2] DBPrefix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdata = []
            for estr in edata:
                if not estr: continue
                # Reformat. Note. Prefix without : is no longer recognised.
                if rje.matchExp('^%s:(\S.*)$' % field,estr): estr = rje.matchExp('^%s:(\S.*)$' % field,estr)[0]
                elif rje.matchExp('^%s:(\S.*)$' % fupper,estr): estr = rje.matchExp('^%s:(\S.*)$' % fupper,estr)[0]
                if field in self.dict['StripVar']: estr = string.split(estr,self.dict['StripVar'][field])[0]
                if field in self.list['DBPrefix']:    # Automatically removed above: add back for given fields
                    if field in self.list['KeepCase']: estr = '%s:%s' % (field,estr)
                    else: estr = '%s:%s' % (fupper,estr)
                # Add reformatted ID unless Bad or already in list. (Want to merge reformatted variants.)
                if estr in xdata: continue
                elif estr in self.list['BadID']: continue
                else: xdata.append(estr)
            ### ~ [3] Rejoin data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SortXRef'): xdata.sort()
            if inplace: entry[field] = string.join(xdata,self.getStr('SplitChar'))
            else: return string.join(xdata,self.getStr('SplitChar'))
        except: self.errorLog('Problem with reformatEntryField(%s)' % field); raise
#########################################################################################################################
    def addXRef(self,entry):   ### Add entry data to main XRef table
        '''Add entry data to main XRef table.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            keyid = self.getStr('KeyID')
            ### ~ [1] Get list of XRef keys for mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xkeys = self.xref(entry)
            if not xkeys: return
            if self.getBool('UniqueXRef'): xkeys = [xkeys]
            ### ~ [2] Update data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adata = {}
            for afield in entry:
                if afield == keyid or not entry[afield]: continue
                if afield in self.list['SplitSkip']: adata[afield] = entry[afield]
                else: adata[afield] = string.split(entry[afield],self.getStr('SplitChar'))
            #self.debug('%s -> %s' % (entry,adata))
            for xid in xkeys:
                xentry = xdb.data(xid)
                #self.debug('%s -> %s' % (xid,xentry))
                for afield in adata:
                    if afield == '#': continue
                    if afield not in xdb.fields(): xdb.addField(afield,evalue=''); xdb.index(afield)
                    if afield in self.list['SplitSkip']:
                        if not xentry[afield]: xentry[afield] = adata[afield]
                    else:
                        if xentry[afield]: xdata = string.split(xentry[afield],self.getStr('SplitChar'))
                        else: xdata = []
                        for a in adata[afield]:
                            if a and a not in xdata + self.list['BadID']:
                                xdata.append(a)
                                if a not in xdb.index(afield): xdb.dict['Index'][afield][a] = []
                                xdb.dict['Index'][afield][a].append(xid)
                        if self.getBool('SortXRef'): xdata.sort()
                        xentry[afield] = string.join(xdata,self.getStr('SplitChar'))
                    #self.bugPrint('%s ... %s' % (afield,xentry))
                #self.debug('%s -> %s' % (xid,xentry))
        except: self.errorLog('Problem with addXRef()'); raise
#########################################################################################################################
    def xref(self,mapdict,xfield=None,mapfields=[],altfields=[],fullmap=None,unique=None,usedict=False,strictmap=False):    ### Cross-reference xentry to keyid
        '''
        Cross-reference xid to keyid (or idfield if given). Identifiers are always mapped first onto the keyid field, and
        then subsequently return the xfield mapping if this is not the keyid. If mapdict is a dictionary then any
        mapfields fields will try to be mapped onto XRef fields. If a mapfields field is in altfields then it will also
        try to be mapped directly onto KeyID, else against an altfield (if present), i.e. all xfield + altfields values
        will try to be mapped against xfield + altfields. If mapdict is a string then the value will be used for all
        mapfields. If strictmap=False, it will also be used for xfield and all altfields.
        >> mapdict:dict = Dictionary of {field:value} to map onto KeyID(s). Can be str, assigned to ALL mapfields.
        >> xfield:str [None] = Field to map IDs onto. Will use self.str['KeyID'] if None
        >> mapfields:list [] = List of fields to use for mapping from. (self.list['MapFields'] if None)
        >> altfields:list [] = List of alternative fields for mapping xfield identifiers onto. (self.list['AltKeys'] if None and xfield==keyid)
        >> fullmap:bool [None] = Whether to use all mapping fields or stop at first hit.
        >> unique:bool [None] = Whether to only return unique mappings, None if missing and False if 2+ mapped IDs.
        >> usedict:bool [False] = Whether to use Mapping dictionary to speed up recurrent mappings. (mapdict=str only)
        >> strictmap:bool [False] = Whether to strictly stick to given mapfields list, or add xfield and altfields.
        << xlist:list = List of mapped IDs, or single ID/None/False if unique=True.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# This method needs reworking, tidying and documenting. #!#
            xdb = self.db('xref')
            if not xfield: xfield = self.getStr('KeyID')    # By default, the KeyID is returned.
            if xfield not in xdb.fields(): raise ValueError('XRef field "%s" not found!' % xfield)
            if usedict and xfield not in self.dict['Mapping']: self.dict['Mapping'][xfield] = {}
            ## ~ [0a] Setup mapping field lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not mapfields: mapfields = self.list['MapFields']    # This is the subset of mapdict fields to use
            if not altfields and xfield == self.getStr('KeyID'):
                altfields = self.list['AltKeys']    # This is the set of additional fields to map KeyID fields onto
            # Note that xfield + altfields in the mapdict and xref will all try to be mapped onto each other.
            # xfield and altfields should therefore be of the same type
            if fullmap == None: fullmap = self.getBool('FullMap')
            if unique == None: unique = self.getBool('UniqueXRef')
            ## ~ [0b] Setup the dictionary of mapping fields and values ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            entry = mapdict
            try:
                entry.keys()
                usedict = False     # Only usedict if mapdict is a single str
            except:
                entry = {}
                for mfield in mapfields: entry[mfield] = mapdict    # Enables mapping via any mapfield.
                if not strictmap:                                   # Try mapping to keyid identifiers directly first.
                    for mfield in [xfield] + altfields: entry[mfield] = mapdict
            if usedict and mapdict in self.dict['Mapping'][xfield]: return self.dict['Mapping'][xfield][mapdict]
            ### ~ [1] Get list of XRef keys for mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xkeys = []      # First get the set of xref keys that mapdict maps onto
            ## ~ [1a] KeyID and AltKeys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for tfield in [xfield] + altfields:     # Try mapping ONTO these fields, starting with KeyID
                if tfield not in xdb.fields(): continue
                elif strictmap and tfield not in mapfields: continue
                for mfield in [xfield] + altfields: # Try mapping FROM these fields, starting with KeyID
                    if mfield not in entry: continue
                    elif strictmap and tfield not in mapfields: continue
                    #mupper = xfield.upper()
                    mdata = self.reformatEntryField(entry,mfield,inplace=False) # Should now match target dbprefix
                    if mfield == self.getStr('KeyID') or mfield in self.list['SplitSkip']: mids = [mdata]
                    else: mids = string.split(mdata,self.getStr('SplitChar'))
                    for xid in mids[0:]:
                        if xid in self.list['BadID'] or not xid: continue
                        if xid in xdb.index(tfield): xkeys += xdb.index(tfield)[xid]
                        #elif '%s:%s' % (xfield,id) in xdb.index(xfield): xkeys += xdb.index(xfield)['%s:%s' % (xfield,id)]
                        #elif '%s:%s' % (mupper,id) in xdb.index(xfield): xkeys += xdb.index(xfield)['%s:%s' % (mupper,id)]
                    self.bugPrint('%s -> %s: %s? => %s' % (mfield,tfield,mids,xkeys))
                    if xkeys and not fullmap: break
            ## ~ [1b] MapFields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug(mapfields)
            for mfield in mapfields:    # This is a direct mapping of mfield content to mfield content
                if xkeys and not fullmap: break
                if mfield not in entry: continue
                if mfield not in xdb.fields(): continue
                #mupper = mfield.upper()
                mdata = self.reformatEntryField(entry,mfield,inplace=False) # Should now match target dbprefix
                if mfield in self.list['SplitSkip']: mids = [mdata]
                else: mids = string.split(mdata,self.getStr('SplitChar'))
                for xid in mids[0:]:
                    if xid and xid not in self.list['BadID'] and xid in xdb.index(mfield): xkeys += xdb.index(mfield)[xid]
                    #elif '%s:%s' % (mfield,id) in xdb.index(mfield): xkeys += xdb.index(mfield)['%s:%s' % (mfield,id)]
                    #elif '%s:%s' % (mupper,id) in xdb.index(mfield): xkeys += xdb.index(mfield)['%s:%s' % (mupper,id)]
                #self.bugPrint('%s -> %s: %s? => %s' % (mfield,mfield,mids,xkeys))
                if xkeys and not fullmap: continue
            xlist = rje.sortUnique(xkeys)
            if len(xlist) > 1: self.debug('Final %s => %s' % (xfield,xlist))
            else: self.bugPrint('Final %s => %s' % (xfield,xlist))
            ### ~ [2] Map mapped KeyIDs onto idfield ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xfield != self.getStr('KeyID'):
                klist = xlist[0:]
                xlist = []
                for kid in klist:
                    kdata = xdb.data(kid)[xfield]
                    if kdata in self.list['BadID'] or not kdata: continue
                    if xfield in self.list['SplitSkip']: xlist += [kdata]
                    else: xlist += string.split(kdata,xdb.getStr('SplitChar'))
                xlist = rje.sortUnique(xlist)
            while '' in xlist: xlist.remove('')
            ### ~ [3] Return as determined by Unique ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if unique:
                if not xlist: xref = None
                elif len(xlist) == 1: xref = xlist[0]
                else: xref = False
            else:
                if len(xlist) > 1: self.warnLog('Entry mapped to multiple %s IDs: %s' % (xfield,string.join(xlist,'|')),'one-to-many',suppress=self.i()>0)
                xref = xlist[0:]
            if usedict: self.dict['Mapping'][xfield][mapdict] = xref
            return xref
        except: self.errorLog('%s.xref(%s -> %s) failure' % (self.prog(),mapdict,xfield)); return False
#########################################################################################################################
    def fileXRef(self,filexref=None,unique=False,mapxref=[]): ### Map appropriate XRef field data onto another data file
        '''
        Map appropriate XRef field data onto another data file.
        >> filexref:str = File to map XRefs onto.
        >> unique:bool [False] = Whether to only return unique mappings.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            xdb = self.db('xref')
            keyid = self.getStr('KeyID')
            if self.list['XRefs']: xfields = self.list['XRefs']
            else: xfields = xdb.fields()
            if filexref: fdb = db.addTable(filexref,mainkeys='#',lists=False,name='filexref',expect=True,replace=True)
            else:
                fdb = db.addEmptyTable('xmap',['#','MapXRef',keyid],['#'],log=True)
                if 'MapXRef' not in self.list['KeepCase']:
                    self.list['MapXRef'] = rje.listUpper(self.list['MapXRef'])
                    self.printLog('#CASE','%s "MapXRef" identifiers converted to uppercase.' % (rje.iLen(self.list['MapXRef'])))
                for mid in self.list['MapXRef']:
                    kid = self.xref(mid,unique=True,usedict=True)
                    if not kid: kid = ''
                    fdb.addEntry({'#':fdb.entryNum(),'MapXRef':mid,keyid:kid})
            jfields = fdb.fields()
            for xfield in xfields[0:]:
                if xfield not in jfields:
                    jfields.append(xfield)
                    if not mapxref: self.printLog('#ADD','Field "%s" added to %s for output.' % (xfield,filexref))
                if xfield not in xdb.fields():
                    self.warnLog('Xref field "%s" not found in xrefdata.' % xfield)
                    xfields.remove(xfield)
            ndb = db.addEmptyTable('join',jfields,['#'],log=True)
            ### ~ [1] Merge data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nx = 0  # Entries with no mapping
            for entry in fdb.entries()[0:]:
                #X#for field in fdb.fields(): self.reformatEntryField(entry,field)
                xids = self.xref(entry,unique=unique,usedict=True)
                if not xids:
                    self.warnLog('No mapping for %s' % entry,dev=True,suppress=True)
                    nx +=1; continue
                elif unique: xids = [xids]
                for xid in xids:
                    for xfield in xfields: entry[xfield] = xdb.data(xid)[xfield]
                    ndb.addEntry(entry)
            ### ~ [2] Clean up tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db.deleteTable('xref')
            if filexref: db.deleteTable('filexref')
            else: db.deleteTable('xmap')
            ndb.info['Name'] = 'xref'
            return ndb
        except: self.errorLog('Problem with fileXRef()'); raise
#########################################################################################################################
    ### <3> ### XRef Output Methods                                                                                     #
#########################################################################################################################
    def xrefList(self,field):   ### Output list of field values to *.*.txt
        '''Output list of field values to *.*.txt.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if field in ['#','AutoID']: return self.printLog('#SKIP','Skipping XRefList output for field "%s".' % field)
            xdb = self.db('xref')
            xfile = '%s.%s.txt' % (self.basefile(),field)
            rje.backup(self,xfile,appendable=False)
            ### ~ [1] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if field in self.list['SplitSkip']: xlist = rje.sortKeys(xdb.index(field))
            else: xlist = rje.sortKeys(xdb.index(field,splitchar=self.getStr('SplitChar')))
            if '' in xlist: xlist.remove('')
            if ' ' in xlist: xlist.remove(' ')
            open(xfile,'w').write(string.join(xlist,'\n'))
            if not self.getStrLC('Rest'): self.printLog('#OUT','Saved %s %s IDs to %s.' % (rje.iLen(xlist),field,xfile))
            return xfile
        except: self.errorLog('Problem with xrefList(%s)' % field); raise
#########################################################################################################################
    ### <4> ### Special Join Method                                                                                     #
#########################################################################################################################
    def join(self): ### Special standalone table join method
        '''Special standalone table join method for list of FILE:key1|...|keyN:JoinField [].'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            joinlist = []
            for joiner in self.list['Join']:
                try:
                    [jfile,jkeys,jjoin] = string.split(joiner,':')
                    jkeys = string.split(jkeys,'|')
                    jdb = db.addTable(jfile,jkeys,ignore=self.list['Comments'])
                    if self.list['XRefs']: joinlist.append((jdb,jjoin,rje.listIntersect(jdb.fields(),self.list['XRefs'])))
                    else: joinlist.append((jdb,jjoin))
                    if not self.getStrLC('Basefile'): self.basefile(rje.baseFile(jfile))
                except: self.errorLog('Failed to process joiner "%s"' % joiner); raise
            ### ~ [1] Make join table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(joinlist) > 1: jdb = db.joinTables('join',joinlist,empties=not self.getBool('NaturalJoin'))
            else: jdb = joinlist[0][0]
            ### ~ [2] Rename Fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dev():  # Rename fields
                if self.list['XRefs'] and self.list['NewHeaders']:
                    if len(self.list['XRefs']) != len(self.list['NewHeaders']):
                        self.errorLog('Cannot rename joined fields if xrefs=LIST and newheaders=LIST are of different lengths.',printerror=False)
                    else:
                        for i in range(len(self.list['NewHeaders'])):
                            jdb.renameField(self.list['XRefs'][i],self.list['NewHeaders'][i])
            ### ~ [3] Filter on IDList~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                if self.list['IDList']:
                    if self.getStr('KeyID') in jdb.fields():
                        jdb.dropEntriesDirect(self.getStr('KeyID'),self.list['IDList'],inverse=True)
                    else: self.errorLog('Cannot filter on IDList: %s not found in fields!' % self.getStr('KeyID'))
            ### ~ [4] Save to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdb.saveToFile()                #!# Add outfile parameter? #!#
            jdb.setStr({'Name':'xref'})
            for field in self.list['XRefList']:
                if field in jdb.fields(): self.xrefList(field)
            return jdb
        except: self.errorLog('Problem with XRef.join()'); return False
#########################################################################################################################
    ### <5> ### XRef REST Methods                                                                                       #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT:

        tab = main table of identified elements. [tdt]
        mapped = pairs of provided identifiers and the primary ID mapped onto. [tdt]
        failed = list of identifiers that failed to map. [list]

        In addition, there will be a tab per field of the XRef file listing the sorted unique identifiers mapped.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb = self.db('xref')
            if self.getStrLC('Rest'): self.list['XRefList'] = xdb.fields()[0:]
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        xdb = self.db('xref')
        outlist = ['tab','mapped','failed']
        for field in xdb.fields():
            if field.lower() in ['log','xref','status']: continue
            if field not in outlist: outlist.append(field.lower())
        return outlist
#########################################################################################################################
### End of SECTION II: XRef Class                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################

#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: XRef(mainlog,cmd_list).run()

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
