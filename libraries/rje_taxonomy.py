#!/usr/bin/python

# See below for name and description
# Copyright (C) 2014 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_taxonomy
Description:  Downloads, reads and converts Uniprot species codes and NCBI Taxa IDs
Version:      1.3.0
Last Edit:    22/06/18
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to download and interconvert between NCBI Taxa IDs and Uniprot species codes and species
    names. It uses two main files: speclist.txt from Uniprot and node.dmp from NCBI taxonomy. These will be downloaded
    if missing and download=T, else can be manually downloaded from:

    - http://www.uniprot.org/docs/speclist.txt
    - ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

    These will be saved in the directory given by taxdir=PATH/ (default ./SourceData/) and will have the download
    date inserted. This enables several versions to be stored together if desired and selected using the sourcedate=DATE
    option, where DATE is in the form YYYY-MM-DD.

    Alternatively, specfile=FILE, taxmap=FILE and namemap=FILE can be used to give other files of the same format. The
    SpecFile should follow key format features of the organism codes but do not need the headers. The TaxMap
    file simply uses the first three columns of the nodes.dump file (separated by "\t|\t"), which correspond to (1) the
    Taxa ID, (2) the parent Taxa ID, and (3) the rank of the entry. These ranks are used to determine which taxa are
    output if rankonly=T. By default, this is "species" and "subspecies". NameMap uses all four fields.

    To extract and/or convert a set of Taxa IDs, a list of taxa should be given using taxin=LIST, where LIST is either a
    comma separated list of taxa, or a file containing one taxon per line. Taxa may be a mix of NCBI Taxa IDs, Uniprot
    species codes and case-insensitive (but exact match) species names. By default, all taxa will be combined but if
    batchmode=T then each TaxIn element will be processed individually. When batchmode=T, individual list elements can be
    files containing taxa. (This will not work in batchmode=F, unless only a single file is given.)

    Taxa are first mapped on to NCBI Taxa IDs. Unless nodeonly=T, taxa will also be mapped on to all of their child taxa,
    as defined in nodes.dmp. If rankonly=T then only those taxa with a rank matching ranktypes=LIST will be retained. IDs
    can be further restricted by supplying a list with restrictid=LIST, which will limit mapped IDs to those within the
    list given. (Note that this list could itself be created by a previous file of rje_taxonomy and given as a file.)
    Taxa IDs are then mapped on the Uniprot species codes and species names using the SpecFile data. If missing from this
    file, scientific names will be pulled out of the NCBI NameMap file instead. NOTE: Uniprot is used first because NCBI
    has more redundant taxonomy assignments.

    Output is determined by the taxout=LIST option, which is set by default to 'taxid'. Four possible output types are
    permitted:
    - taxid = NCBI Taxa IDs (e.g. 9606 or 7227)
    - spcode = Uniprot species codes (e.g. HUMAN or DROME)
    - name = Scientific name (e.g. Homo sapiens or Drosophila melanogaster)
    - common = Common name (e.g. Human or Fruit fly)

    These will be output as lists to BASE.TYPE.txt files, where BASE is set by basefile=X (using the first taxin=LIST
    element if missing) and TYPE is the taxout type. If batchmode=T then a separate set of files will be made for each
    element of the TaxIn list, using BASE.TAXIN.TYPE.txt file naming.

Commandline:
    ### ~ SOURCE DATA OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    specfile=FILE       : Uniprot species code download. [speclist.txt]
    taxmap=FILE         : NCBI Node Dump File [nodes.dmp]
    namemap=FILE        : NCBI Name mapping file [names.dmp]
    taxdir=PATH/        : Will look in this directory for input files if not found ['./SourceData/']
    sourcedate=DATE     : Source file date (YYYY-MM-DD) to preferentially use [None]
    download=T/F        : Whether to download files directly from websites where possible if missing [True]
    ### ~ TAXONOMY CONVERSION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    taxin=LIST          : List of Taxa IDs, Uniprot species codes (upper case) and/or common/scientific names []
    batchmode=T/F       : Treat each element of taxin as a separate run (will be used for output basefile) [False]
    taxout=LIST         : List of output formats (taxid/spcode/name/common/all) [taxid]
    nodeonly=T/F        : Whether to limit output to the matched nodes (i.e. no children) [False]
    rankonly=T/F        : Whether to limit output to species-level taxonomic codes [False]
    ranktypes=LIST      : List of Taxon types to include if rankonly=True [species,subspecies,no rank]
    restrictid=LIST     : List of Taxa IDs to restrict output to (i.e. output overlaps with taxin) []
    basefile=X          : Results file prefix. Will use first taxin=LIST term if missing [None]
    taxtable=T/F        : Whether to output results in a table rather than text lists [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Initial working version with rje_ensembl.
    # 1.0 - Fully functional version with modified viral species code creation.
    # 1.1.0 - Added parsing of yeast strains.
    # 1.2.0 - Added storage of Parents.
    # 1.3.0 - taxtable=T/F        : Whether to output results in a table rather than text lists [False]
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Add a batch mode
    # [Y] : Add a restricted ID list.
    # [ ] : Add a special option for generating a cut-down nodes.dmp restricted to (real/virtual) Uniprot species?
    # [Y] : Add getSpCode() method and replace in rje_ensembl etc.
    # [ ] : Add option [TRUE] for loading and saving of an spcode.DATE.tdt file for invented species codes.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_Taxonomy', '1.3.0', 'June 2018', '2014')
    description = 'Downloads, reads and converts Uniprot species codes and NCBI Taxa IDs'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_obj.zen()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        cmd_help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if cmd_help > 0:
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class Taxonomy(rje_obj.RJE_Object):
    '''
    Taxonomy Class. Author: Rich Edwards (2014).

    Str:str
    - NameMap = NCBI Name mapping file [names.dmp]
    - SourceDate = Source file date (YYYY-MM-DD) to preferentially use [None]
    - SourcePath = Set by taxdir=PATH. Will look in this directory for input files if not found ['SourceData/']
    - SpecFile = Uniprot species code download
    - TaxMap = NCBI Node Dump File

    Bool:boolean
    - BatchMode = Treat each element of taxin as a separate run (will be used for output basefile) [False]
    - Download = Whether to download files directly from websites where possible if missing [True]
    - NodeOnly = Whether to limit output to the matched nodes (i.e. no children) [False]
    - RankOnly = Whether to limit output to species-level taxonomic codes [False]
    - Setup = Whether the setup() has been run.
    - TaxTable=T/F        : Whether to output results in a table rather than text lists [False]

    Int:integer

    Num:float
    
    List:list
    - RestrictID = List of Taxa IDs to restrict output to (i.e. output overlaps with taxin) []
    - RankID = List of Taxa IDs corresponding to Species []
    - RankTypes = List of Taxon types to include if rankonly=True [species,subspecies]
    - TaxIn = List of Taxa IDs, Uniprot species codes (upper case) and/or common/scientific names []
    - TaxOut = List of output formats (taxid/spcode/name/common) [taxid]

    Dict:dictionary
    - NameMap = Dictionary of {TaxID:Scientific Name}       ??? Is this used ???
    - Parent = Dictionary of {TaxID:Parent TaxID}
    - Rank = Dictionart of {TaxID:Taxonomic ranking}
    - SourceURL = dictionary of {Data Source: download URL}
    - TaxDict = Dictionary of {'TaxID':{'SpCode':X,'Name':X,'Common':X}}
    - TaxMap = Dictionary of {'TaxID':[Child TaxIDs]}

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['NameMap','SourceDate','SourcePath','SpecFile','TaxMap']
        self.boollist = ['BatchMode','Download','NodeOnly','RankOnly','Setup','TaxTable']
        self.intlist = []
        self.numlist = []
        self.listlist = ['RankTypes','RestrictID','RankID','TaxIn','TaxOut']
        self.dictlist = ['NameMap','Parent','Rank','SourceURL','TaxDict','TaxMap']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'TaxOut':'taxid','SourcePath':'SourceData/','SpecFile':'speclist.txt','TaxMap':'nodes.dmp',
                     'NameMap':'names.dmp'})
        self.setBool({'Download':True,'Setup':False,'TaxTable':False})
        self.setInt({})
        self.setNum({})
        self.list['TaxOut'] = ['taxid']
        self.list['RankTypes'] = ['species','subspecies','no rank']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.dict['SourceURL'] = {'SpecFile':'http://www.uniprot.org/docs/speclist.txt',
             'TaxMap':'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz',
             'NameMap':'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'}
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self._setForkAttributes()   # Delete if no forking
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
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                self._cmdReadList(cmd,'date',['SourceDate'])   # String representing date YYYY-MM-DD
                self._cmdRead(cmd,'path','SourcePath','taxdir')  # String representing directory path
                self._cmdReadList(cmd,'file',['NameMap','SpecFile','TaxMap'])  # String representing file path
                self._cmdReadList(cmd,'bool',['BatchMode','Download','NodeOnly','RankOnly','TaxTable'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['RankTypes','RestrictID','TaxIn','TaxOut'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.getBool('TaxTable'): self.setBool({'BatchMode':True})
            ### ~ [2] ~ Single Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('BatchMode'): return self.mapTaxa(self.list['TaxIn'],self.list['TaxOut'],self.getBool('NodeOnly'),self.getBool('RankOnly'))
            ### ~ [3] ~ Batch Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('TaxTable'):
                tdb = self.db().addEmptyTable('taxtable',['TaxIn']+self.list['TaxOut'],['TaxIn'])
            basefile = self.baseFile()
            for taxa in self.list['TaxIn'][0:]:
                self._cmdReadList('taxin=%s' % taxa,'list',['TaxIn'])  # List of strings (split on commas or file lines)
                self.setBaseFile('%s.%s' % (basefile,rje.baseFile(taxa,strip_path=True)))
                taxdict = self.mapTaxa(self.list['TaxIn'],self.list['TaxOut'],self.getBool('NodeOnly'),self.getBool('RankOnly'),savetaxout=not self.getBool('TaxTable'))
                if self.getBool('TaxTable'):
                    tentry =  {'TaxIn':taxa}
                    for tfield in taxdict: tentry[tfield] = string.join(taxdict[tfield],'|')
                    tdb.addEntry(tentry)
            self.baseFile(basefile)
            if self.getBool('TaxTable'): tdb.saveToFile()
            return True
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self,force=False,parents=True):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Setup') and not force: self.printLog('#SETUP','Taxonomy setup already complete.'); return True
            if not self.setupSourceData(): raise IOError
            if not self.getStrLC('Basefile'):
                if self.getBool('BatchMode'): self.setBaseFile('batch')
                elif self.list['TaxIn']: self.setBaseFile(rje.baseFile(self.list['TaxIn'][0],strip_path=True))
            self.list['TaxOut'] = string.join(self.list['TaxOut']).lower().split()
            if 'all' in self.list['TaxOut']: self.list['TaxOut'] = ['taxid','spcode','name','common']
            self.list['RankID'] = []
            ### ~ [2] TaxMap Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxmap = self.dict['TaxMap'] = {}
            tx = 0; px = 0; fx = 0
            for tline in open(self.getStr('TaxMap'),'r').readlines():
                self.progLog('\r#TAXID','Reading %s: %s TaxID' % (self.getStr('TaxMap'),rje.iStr(tx)))
                #try: (child,parent,taxtype) = rje.matchExp('^(\d+)\s+\|\s+(\d+)\s+\|\s+(\S+)\s+',tline)
                try: (child,parent,taxtype) = string.split(tline,'\t|\t')[:3]
                except: fx += 1; self.debug(tline); continue
                self.dict['Rank'][child] = taxtype
                if parent not in taxmap: taxmap[parent] = []
                if not taxmap[parent]: px += 1
                if taxtype in self.list['RankTypes']: self.list['RankID'].append(child)
                if child not in taxmap: taxmap[child] = []
                taxmap[parent].append(child); tx += 1
                if child in self.dict['Parent']: self.warnLog('Child TaxID "%s" already has parent!' % child)
                if parents and child != parent: self.dict['Parent'][child] = parent
            self.printLog('\r#TAXID','%s TaxID (%s parent taxa) read from %s; %s failed.' % (rje.iStr(tx),rje.iStr(px),self.getStr('TaxMap'),rje.iStr(fx)))
            self.printLog('#SPEC','%s TaxID mapped to %s RankTypes' % (rje.iLen(self.list['RankID']),string.join(self.list['RankTypes'],'/')))
            if self.test():
                pcheck = 0
                for tax in taxmap:
                    if taxmap[tax]: pcheck += 1
                self.printLog('#TEST','%s parent taxa with listed children' % rje.iStr(pcheck))
                if px != pcheck: raise ValueError
            ### ~ [3] NameMap Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('MemSaver'):
                taxdict = self.dict['TaxDict']
                ## ~ [3a] SpecFile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                tx = 0; cx = 0; taxid = None
                for tline in open(self.getStr('SpecFile'),'r').readlines():
                    self.progLog('\r#SPEC','Reading %s species data: %s TaxID' % (self.getStr('SpecFile'),rje.iStr(tx)))
                    nmatch = rje.matchExp('^(\S+)\s+\S+\s+(\d+):\s+N=(\S.+)\s*$',tline)
                    if nmatch:
                        taxid = nmatch[1]; tx += 1
                        taxdict[taxid] = {'spcode': nmatch[0], 'name': nmatch[2]}
                    elif taxid and rje.matchExp('C=(\S.+)\s*$',tline): taxdict[taxid]['common'] = rje.matchExp('C=(\S.+)\s*$',tline)[0]; cx += 1
                self.printLog('\r#SPEC','%s species codes/names and %s common names read from %s.' % (rje.iStr(tx),rje.iStr(cx),self.getStr('SpecFile')))
                ## ~ [3b] NCBI names.dmp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                tx = 0
                for tline in open(self.getStr('NameMap'),'r').readlines():
                    self.progLog('\r#SPEC','Reading %s species names: %s TaxID' % (self.getStr('NameMap'),rje.iStr(tx)))
                    tdata = string.split(tline,'\t|\t')
                    if not tdata[3].startswith('scientific name'): continue
                    taxid = tdata[0]
                    if taxid not in taxdict: taxdict[taxid] = {'name': tdata[1]}; tx += 1
                self.printLog('\r#SPEC','%s extra species names read from %s.' % (rje.iStr(tx),self.getStr('NameMap')))
            ### ~ [4] Species code table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            spfile = '%sspcode.%s.tdt' % (self.getStr('SourcePath'),self.getStr('SourceDate'))
            self.db().addTable(spfile,['Species'],name='SpCode',expect=False)
            ### ~ [5] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setBool({'Setup':True})
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def setupSourceData(self):    ### Main class setup method.                                                      #V0.0
        '''Setup and optionally download source data.'''
        try:### ~ [0] Setup Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            sdb = db.addEmptyTable('Source',['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            rje.mkDir(self,self.getStr('SourcePath'),True)
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~ SETUP TAXONOMIC DATA ~~~~~~~~~~~~~~~~~~~~~ ##',timeout=False)
            ### ~ [1] Uniprot species codes file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.sourceDataFile('SpecFile'): raise IOError
            ### ~ [2] NCBI TaxID mapping file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.sourceDataFile('TaxMap'): raise IOError
            if not self.sourceDataFile('NameMap'): raise IOError
            return True
        except: self.errorLog('%s.setupSourceData failure' % self); return False
#########################################################################################################################
    ### <3> ### Main Taxonomy Mapping Methods                                                                           #
#########################################################################################################################
    def mapTaxa(self,taxin,taxout=['spcode'],nodeonly=False,rankonly=False,savetaxout=True):    ### Takes a list of Taxa and returns mapped Taxa data
        '''
        Takes a list of Taxa and returns mapped Taxa data.
        >> taxin:str or list of taxon identifiers to map from.
        >> taxout:str or list of taxa output formats
        >> nodeonly:bool = whether to limit TaxID mapping to the precise matching nodes (else include children)
        >> rankonly:bool = whether to limit TaxID to those matching self.list['RankTypes'] taxon types.
        >> savetaxout:bool [True] = Whether to save the TaxOut list to a text file
        << taxoutlist:list of mapped taxa if taxout is a string, OR
        << taxoutdict:dict of mapped taxa if taxout is a list
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tlist = True
            try: taxout.sort()
            except: tlist = False
            if tlist:
                if not taxout: return {}
                taxout = [taxout]
            elif not taxout: return []
            ### ~ [2] ~ Map to TaxID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxid = self.mapToTaxID(self.list['TaxIn'],nodeonly,rankonly)
            if self.list['RestrictID']:
                tx = len(taxid)
                taxid = rje.listIntersect(taxid,self.list['RestrictID'])
                self.printLog('#TAXID','%s of %s TaxID in %s Restricted IDs.' % (rje.iLen(taxid),rje.iStr(tx),rje.iLen(self.list['RestrictID'])))
            ### ~ [3] ~ Map TaxID and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxdict = {}; taxoutdict = {}
            for taxout in self.list['TaxOut']:
                taxout = taxout.lower()
                if taxout == 'taxid':
                    taxoutlist = taxid
                elif taxout in ['spcode','name','common']:
                    if not taxdict: taxdict = self.taxDict(taxid)
                    taxoutlist = []
                    for t in taxid:
                        try: taxoutlist.append(taxdict[t][taxout])
                        except: self.warnLog('No "%s" data for TaxID %s' % (taxout, t),'Missing_%s' % taxout,suppress=True)
                    taxoutlist.sort()
                else: self.errorLog('TaxOut format "%s" not recognised' % taxout,printerror=False); continue
                taxoutdict[taxout] = taxoutlist
                if savetaxout:
                    if not taxoutlist: self.printLog('#OUT','No %s IDs to output' % taxout); continue
                    tfile = '%s.%s.txt' % (self.baseFile(),taxout)
                    rje.backup(self,tfile)
                    open(tfile,'w').write(string.join(taxoutlist,'\n'))
                    self.printLog('#OUT','%s %s IDs output to %s.' % (rje.iLen(taxoutlist), taxout, tfile))
            if tlist: return taxoutdict
            return taxoutlist
        except: self.errorLog('Problem during %s mapTaxa.' % self); raise
#########################################################################################################################
    def mapToTaxID(self,taxa,nodeonly=False,rankonly=False,log=True,warn=True):  ### Maps taxa onto TaxID. If taxa is a list, will process each element.
        '''Maps taxa onto TaxID. If taxa is a list, will process each element. Returns a list.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not taxa: return []
            taxid = []
            ### ~ [1] Taxa List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tlist = True
            try: taxa.sort()
            except: tlist = False
            if tlist:
                tx = 0.0; ttot = len(taxa)
                if ttot > 1:
                    for t in taxa:
                        if log: self.progLog('\r#TAXID','Mapping to TaxID: %.1f%%' % (tx/ttot)); tx += 100.0
                        taxid += self.mapToTaxID(t,nodeonly,rankonly,log=False)
                    taxid = rje.sortUnique(taxid)
                    if log:
                        if ttot > 1: self.printLog('\r#TAXID','Mapped %s taxa to %s TaxID' % (rje.iStr(ttot),rje.iLen(taxid)))
                else:
                    t = taxa[0]
                    if log: self.progLog('\r#TAXID','Mapping %s to TaxID...' % t)
                    taxid = rje.sortUnique(self.mapToTaxID(t,nodeonly,rankonly,log=False))
                    if log: self.printLog('\r#TAXID','Mapped %s to %s TaxID' % (t,rje.iLen(taxid)))
                return taxid
            ### ~ [2] Individual taxa ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxmap = self.dict['TaxMap']; rankid = self.list['RankID']
            taxa = '%s' % taxa
            ## ~ [2a] Taxa ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.matchExp('^(\d+)$', taxa):
                #if taxa not in taxmap: self.taxaChildren(taxa)
                #if taxa in rankid: return [taxa]
                if nodeonly:
                    if taxa in rankid or not rankonly: return [taxa]
                    else: return []
                if taxa not in taxmap:
                    if warn: self.warnLog('Cannot find TaxID %s!' % taxa,'Missing_TaxID',suppress=True)
                    return []
                parents = [taxa]
                while parents:
                    taxa = parents.pop(0)
                    #if taxa not in taxmap: self.taxaChildren(taxa)
                    if not rankonly or taxa in rankid: taxid.append(taxa)
                    parents += taxmap[taxa]
                return taxid
            ## ~ [2b] Species Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if taxa == string.replace(taxa.upper(),' ',''):
                greplines = os.popen('grep "%s" %s' % (taxa, self.getStr('SpecFile'))).readlines()
                for entry in greplines:
                    try: taxid.append(rje.matchExp('^%s\s+\S+\s+(\d+):' % taxa,entry)[0])
                    except: pass
                if not taxid and warn: self.warnLog('Cannot find Species Code "%s"!' % taxa,'Missing_SpCode',suppress=True)
                if len(taxid) > 1: self.warnLog('Species Code "%s" hits %d Taxa ID (%s)' % (taxa, len(taxid), string.join(taxid,'|')))
                return self.mapToTaxID(taxid,nodeonly,rankonly,log=False) #taxid
            ### ~ [3] Species name etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxa = taxa.replace('_',' ')
            ## ~ [3a] Grep from Uniprot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            greplines = os.popen('grep -B 2 -i "%s" %s' % (taxa, self.getStr('SpecFile'))).readlines()
            gtaxid = None; comid = []; synid = []
            for entry in greplines:
                try: gtaxid = rje.matchExp('^\S+\s+\S+\s+(\d+):',entry)[0]
                except: pass
                if rje.matchExp('s=(%s)\s*$' % taxa.lower(),entry.lower()): synid.append(gtaxid)
                elif rje.matchExp('c=(%s)\s*$' % taxa.lower(),entry.lower()): comid.append(gtaxid)
                elif rje.matchExp('=(%s)\s*$' % taxa.lower(),entry.lower()): taxid.append(gtaxid)
            if not taxid: taxid = comid
            if not taxid: taxid = synid
            if not taxid and warn: self.warnLog('Cannot find Taxon name "%s" in Uniprot!' % taxa,'Missing Taxon',suppress=True)
            if len(taxid) > 1:
                #self.bugPrint(string.join(greplines))
                #self.debug('%s %s %s' % (taxid,comid,synid))
                if warn: self.warnLog('Species Code "%s" hits %d Taxa ID (%s)' % (taxa, len(taxid), string.join(taxid,'|')))
            if taxid: return self.mapToTaxID(taxid,nodeonly,rankonly,log=False) #taxid
            #self.debug(taxid)
            ## ~ [3b] Grep from NCBI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            greplines = os.popen('grep -i -e "\t%s\t" %s' % (taxa, self.getStr('NameMap'))).readlines()
            for entry in greplines:
                try:
                    #gtaxid = rje.matchExp('^(\d+)\s+\S\s+(\S.+)$',entry)
                    gtaxid = string.split(entry,'\t|\t')
                    if gtaxid[1].lower() == taxa.lower(): taxid.append(gtaxid[0])
                    elif gtaxid[2] and gtaxid[2].lower() == taxa.lower(): taxid.append(gtaxid[0])
                except: pass
            if len(taxid) > 1 and warn: self.warnLog('Species Code "%s" hits %d Taxa ID (%s)' % (taxa, len(taxid), string.join(taxid,'|')))
            return self.mapToTaxID(taxid,nodeonly,rankonly,log=False) #taxid
        except: self.errorLog('%s.mapToTaxID() error' % (self)); raise
#########################################################################################################################
    def taxDict(self,taxid,store=False,skipuni=False):    ### Extracts taxonomy details from SpecFile for taxid
        '''Extracts taxonomy details from SpecFile for taxid. If taxid is a list, will process each element.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxdict = {}
            ### ~ [1] Taxa List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tlist = True
            try: taxid.sort()
            except: tlist = False
            if tlist:
                tx = 0.0; ttot = len(taxid); mx = 0
                for t in taxid:
                    self.progLog('\r#SPEC','Extracting Uniprot species details: %.1f%%' % (tx/ttot)); tx += 100.0
                    taxdict[t] = self.taxDict(t,store)
                    if not taxdict[t]: mx += 1
                self.printLog('\r#SPEC','Extracted Uniprot/NCBI species details for %s TaxID: %s missing' % (rje.iStr(ttot),rje.iStr(mx)))
                return taxdict
            ### ~ [2] Individual taxa ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxid = '%s' % taxid
            if taxid in self.dict['TaxDict']: return self.dict['TaxDict'][taxid]
            if not skipuni:
                greplines = os.popen('grep -A 1 " %s:" %s' % (taxid, self.getStr('SpecFile'))).readlines()
                for entry in greplines:
                    nmatch = rje.matchExp('^(\S+)\s+\S+\s+(\d+):\s+N=(\S.+)\s*$',entry)
                    if nmatch and nmatch[1] != taxid: break # Next taxon
                    if nmatch: taxdict['spcode'] = nmatch[0]; taxdict['name'] = nmatch[2]
                    elif rje.matchExp('C=(\S.+)\s*$',entry): taxdict['common'] = rje.matchExp('C=(\S.+)\s*$',entry)[0]
            #if not taxdict and taxid in self.list['RankID']: self.warnLog('Cannot find TaxID "%s" in %s!' % (taxid,self.getStr('SpecFile')),'Missing_TaxID',suppress=True)
            ## ~ [2b] ~ Adding missing scientific names from NameMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not taxdict:
                for entry in os.popen('grep -i -e "^%s\t" %s' % (taxid, self.getStr('NameMap'))).readlines():
                    tdata = string.split(entry,'\t|\t')
                    if not tdata[3].startswith('scientific name'): continue
                    tname = tdata[1]
                    if 'name' in taxdict: self.warnLog('TaxID %d hits "%s" and "%s"!' % (taxid, taxdict[name],tname))
                    else: taxdict['name'] = tname
            return taxdict
        except: self.errorLog('%s.taxDict() error' % (self)); raise
#########################################################################################################################
    def taxaChildren(self,taxid):   ### Extracts TaxID children from TaxMap file and updates RankID and TaxMap dicts.
        '''Extracts TaxID children from TaxMap file and updates RankID and TaxMap dicts.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NB. This is very slow and so reading the while.
            self.debug(taxid)
            taxmap = self.dict['TaxMap']
            if taxid in taxmap: return taxmap[taxid]
            ### ~ [1] Parse from TaxMap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            taxmap[taxid] = []
            for tline in os.popen('grep -e "\s%s\s" %s' % (taxid,self.getStr('TaxMap'))).readlines():
                try: (child,parent,taxtype) = rje.matchExp('^(\d+)\s+\|\s+(\d+)\s+\|\s+(\S+)\s+',tline)
                except: continue
                if parent not in taxmap: taxmap[parent] = []
                taxmap[parent].append(child)
                if taxtype in ['species','subspecies']: self.list['RankID'].append(child)
                self.progLog('\r#TAXID','Reading %s: %s TaxID' % (self.getStr('TaxMap'),rje.iLen(taxmap)))
            return taxmap[taxid]
        except: self.errorLog('%s.taxaChildren(%s) error' % (self,taxid)); raise
#########################################################################################################################
    def getSpCode(self,taxa,invent=True,warn=True):  # Returns a species code for input taxon. Can invent a six-letter code.
        '''
        Returns a species code for input taxon. Will create six-letter code if invent=T and spcode not found.
        >> taxon:str or list = Taxon to map. Can be TaxID, SpCode or species name (scientific or Uniprot common)
        >> invent:bool [True] = Whether to invent spcode from species name if real spcode missing.
        >> warn:bool [True] = Whether to warn if invented spcode from species name if real spcode missing.
        << spcode:str or dict = Upper case Uniprot species code (or invented equivalent). Dict if input is list.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tlist = True
            try: taxid.sort()
            except: tlist = False
            ### ~ [1] Setup dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if tlist:
                spdict = {}
                for taxon in taxa: spdict[taxon] = self.getSpCode(taxon,invent)
                return spdict
            ### ~ [2] Get SpCode for one taxon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not taxa: return ''
            if rje.matchExp('^(\d+)$', taxa): taxid = [taxa]; taxa = self.getSpecies(taxa)
            else:
                if taxa == string.replace(taxa.upper(),' ',''): return taxa
                taxid = self.mapToTaxID(taxa,nodeonly=True,rankonly=False,log=False,warn=False)
            if len(taxid) > 1 and warn: self.warnLog('"%s" maps to %d TaxIDs! (Will use %s)' % (taxa,len(taxid),taxid[0]),'SpCode Multiple TaxID',suppress=True)
            if not taxid and warn: self.warnLog('Could not find TaxID for "%s"' % taxa,'SpCode Missing TaxID',suppress=True)
            else:
                taxid = taxid[0]
                taxdict = self.taxDict(taxid,store=not self.getBool('MemSaver'))
                try:
                    if taxdict['spcode']: return taxdict['spcode']
                except: pass
            ### ~ [3] Invent SpCode for one taxon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not invent: return ''
            sdb = self.db('SpCode')
            if not sdb: sdb = self.db().addEmptyTable('SpCode',['Species','SpCode'],['Species'])
            if sdb.data(taxa): return sdb.data(taxa)['SpCode']
            twords = string.split(taxa.replace('_',' ').upper())
            if 'VIRUS' in taxa.upper():
                twords = string.split(taxa.upper().replace('_',' ').replace('VIRUS',' VIRUS').replace('-',' '))
                spcode = ''; virus = False
                for t in twords[1:]:
                    if t == 'VIRUS': virus = True
                    if virus and t[0] not in string.ascii_uppercase + string.digits: break
                    if virus and len(t) < 4 and t[0] in string.digits: spcode += t
                    else: spcode += t[0]
                i = 1
                while len(twords[0][:i]+spcode) < 6 and i <= len(twords[0]): i += 1
                spcode = twords[0][:i]+spcode
            elif taxa.startswith('Saccharomyces cerevisiae') and len(twords) == 3: spcode = twords[2]
            elif len(twords) > 1: spcode = twords[0][:3] + twords[1][:3]
            else: spcode = twords[0][:6]
            if not taxa: spcode = taxid
            ## ~ [3a] Clean up SPCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmp = string.join(spcode); spcode = ''
            for x in tmp:
                if x in string.ascii_uppercase + string.digits: spcode += x
            ## ~ [3b] Return SPCODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if warn: self.warnLog('Species code %s invented for "%s"' % (spcode,taxa),'Invented SpCode',suppress=True)
            sdb.addEntry({'Species':taxa,'SpCode':spcode})
            sdb.saveToFile('%sspcode.%s.tdt' % (self.getStr('SourcePath'),self.getStr('SourceDate')),backup=False,append=True,savekeys=[taxa],log=False)
            return spcode
        except:
            if tlist: self.errorLog('%s.getSpCode(%s) error' % (self,taxa)); raise
            else: self.errorLog('%s.getSpCode() error' % (self)); raise
#########################################################################################################################
    def getSpecies(self,taxon):
        try: return self.taxDict(self.mapToTaxID(taxon,nodeonly=True)[0],store=not self.getBool('MemSaver'))['name']
        except:
            if self.dev() or self.test(): self.errorLog('getSpecies(%s) error' % taxon)
            return ''
#########################################################################################################################
    def getCommon(self,taxon):
        try: return self.taxDict(self.mapToTaxID(taxon,nodeonly=True)[0],store=not self.getBool('MemSaver'))['common']
        except:
            if self.dev() or self.test(): self.errorLog('getCommon(%s) error' % taxon)
            return ''
#########################################################################################################################
### End of SECTION II: Taxonomy Class                                                                                   #
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
    try: Taxonomy(mainlog,cmd_list).run()

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
