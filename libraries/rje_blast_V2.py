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
Module:       rje_blast
Description:  BLAST+ Control Module
Version:      2.18.0
Last Edit:    30/10/17
Copyright (C) 2013  Richard J. Edwards - See source code for GNU License Notice

Function:
    This is an updated BLAST module to utilise the improved BLAST+ library rather than the old Legacy BLAST. During the
    upgrade, other improvements are also being made to the module organisation in line with more recent tools in the
    SeqSuite package. In particular, BLAST Search, Hit and PWAln objects are being replaced by rje_Database tables and
    entries. This will allow greater flexibility in summary outputs for future. GABLAM statistics will also be entered
    directly into a Database table. To minimise memory requirements, these tables can be cleared as each BLAST result is
    read in if the data is not needed. This revised structure will also enable reading of tabular results from other
    searches as required in future.

    In Version 2.0, the old BLAST options are included but these will be upgraded to newer BLAST+ options. To run the
    old version, use BLASTRun.run(oldblast=True).

Commandline:
    ## Search Options ##    
    blastprog=X     : BLAST program to use. blastp=X also recognised. (BLAST -p X) [blastp]
    blasti=FILE     : Input file (BLAST -i FILE) [None]
    blastd=FILE     : BLAST database (BLAST -d FILE) [None]
    formatdb=T/F    : Whether to (re)format BLAST database [False]

    blaste=X        : E-Value cut-off for BLAST searches (BLAST -e X) [1e-4]
    blastv=X        : Number of one-line hits per query (BLAST -v X) [500]
    blastb=X        : Number of hit alignments per query (BLAST -b X) [500]
    tophits=X       : Sets max number of BLAST hits returned (blastb and blastv) [500]

    blastf=T/F      : Complexity Filter (BLAST -F X) [True]
    blastcf=T/F     : Use BLAST Composition-based statistics (BLAST -C X) [False]
    blastg=T/F      : Gapped BLAST (BLAST -g X) [True]
    softmask=T/F    : Whether to use soft masking for searches [True]

    blastopt=FILE   : File containing raw BLAST options (applied after all others) []

    ## Standalone Run Options ##
    savelocal=LIST  : Whether to generate extra output for the local BLAST hits table (GFF3/SAM/TDT/TDTSEQ) []
    reftype=X       : Whether to map SAM/GFF3 hits onto the Qry, Hit, Both or Combined [Hit]
    qassemblefas=T/F: Special mode for running with outfmt=4 and then converting to fasta file [False]
    qcomplete=T/F   : Whether the query sequence should be full-length in qassemblefas output [False]
    qconsensus=X    : Whether to convert QAssemble alignments to consensus sequences (None/Hit/Full) [None]
    qfasdir=PATH    : Output directory for QAssemble alignments [./QFAS/]

    ## GABLAM Parameters ##
    gablamfrag=X    : Length of gaps between mapped residue for fragmenting local hits [100]
    fragmerge=X     : Max Length of gaps between fragmented local hits to merge [0]
    localcut=X      : Cut-off length for local alignments contributing to global GABLAM stats) [0]
    localidcut=PERC : Cut-off local %identity for local alignments contributing to global GABLAM stats [0.0]
    qassemble=T/F   : Whether to fully assemble query stats from all hits [False]
    selfsum=T/F     : Whether to also include self hits in qassemble output [False] * qassemble must also be T *

    ## Output options ##
    blasto=FILE     : Output file (BLAST -o FILE) [*.blast]
    restab=LIST     : Whether to output summary results tables (Run/Search/Hit/Local/GABLAM) [Search,Hit]
    runfield=T/F    : Whether to include Run Field in summary tables. (Useful if appending.) [False]

    ## System Parameters ##
    blastpath=PATH  : Path to BLAST programs ['']
    blast+path=PATH : Path to BLAST+ programs (will use blastpath if not given) ['']
    legacy=T/F      : Whether to run in "legacy" mode using old BLAST commands etc. (Currently uses BLAST) [False]
    oldblast=T/F    : Whether to run with old BLAST programs rather than new BLAST+ ones [False]
    blasta=X        : Number of processors to use (BLAST -a X) [1]
    blastforce=T/F  : Whether to force regeneration of new BLAST results if already existing [False]
    ignoredate=T/F  : Ignore date stamps when deciding whether to regenerate files [False]
    gzip=T/F        : Whether to gzip (and gunzip) BLAST results files if keeping (not Windows) [True]

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, re, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_blast_V1, rje_db, rje_obj, rje_seq, rje_seqlist, rje_sequence, rje_zen
import rje_dismatrix_V2 as rje_dismatrix
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 2.0 - Initial Compilation from rje_blast_V1 V1.14.
    # 2.1 - Tweaking code to work with GOPHER 3.x - removing self.info etc. Added blastObj() method.
    # 2.2 - Added gablamData() to return old-style GABLAM dictionary from table.
    # 2.3 - Added blastCluster() method to return UPC clustering and GABLAM distance matrix from a file.
    # 2.4 - Scrapped BLAST "Run" field to simplify code - keep a single run per BLASTRun object.
    # 2.5 - Minor modifications for SLiMCore UPC generation.
    # 2.6 - Minor bug fixes.
    # 2.7 - Fixed occasional oneline versus description mismatch error. Fixed some localhits bugs.
    # 2.7.1 - Added capacity to keep alignments following GABLAM calculations.
    # 2.7.2 - Fixed bug with hitToSeq fasta output for rje_seqlist.SeqList objects.
    # 2.8.0 - A more significant BLAST e-value setting will filter read results.
    # 2.9.0 - Added     qassemble=T/F   : Whether to fully assemble query stats from all hits [False].
    # 2.9.1 - Updated default BLAST and BLAST+ paths to '' for added modules.
    # 2.10.0 - Added nocoverage calculation based on local alignment table.
    # 2.11.0 - Added localFragFas output method.
    # 2.11.1 - Fixed snp local table revcomp bug. [Check this!]
    # 2.11.2 - Fixed GABLAM calculation bug when '*' in protein sequences.
    # 2.12.0 - Added localidcut %identity filter for GABLAM calculations.
    # 2.13.0 - Added GFF and SAM output for BLAST local tables for GABLAM, PAGSAT etc.
    # 2.14.0 - Updated gablamfrag=X and fragmerge=X usage. Fixed localFragFas position output.
    # 2.15.0 - Fragmerge no longer removes flanks and can be negative for enforced overlap!
    # 2.16.0 - Added qassemblefas mode for generating fasta file from outfmt 4 run.
    # 2.16.1 - Improved error messages for BLAST QAssembly.
    # 2.17.0 - qconsensus=X    : Whether to convert QAssemble alignments to consensus sequences (None/Hit/Full) [None]
    # 2.17.1 - Modified QAssembleFas output sequence names for better combining of hits. Added QFasDir.
    # 2.17.2 - Modified QAssembleFas output file names for better re-running. Fixed major QConsensus Bug.
    # 2.18.0 - Added REST output. Fixed QConsensus=Full bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [ ] : Check/Implement old V1.14 wishlist:
    # |- [ ] : Output of alignment with options for line lengths and numbers +-
    # |- [ ] : OrderAln by any stat ['BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd']
    # |- [ ] : Add standalone running for blast searching?
    # |- [ ] : Replace blast.search with list['Search'] etc?
    # |- [ ] : Add documentation for implementation details - each method?
    # |- [ ] : Locate which classes/methods call BLASTRun.hitToSeq and look to improve reporting etc.
    # |- [ ] : Fix DNA implementation of GABLAM to allow Ordered GABLAM in either direction.
    # |- [ ] : Add positional information to GABLAM dictionary - start and end of aligned portions
    # |- [ ] : Add oritentation of Query and Hit for DNA GABLAM
    # |- [ ] : Check/fix the database format checking of DNA databases
    # |- [Y] : Update to new Module Structre (V2.0)
    # |- [ ] : Check and Tidy this To Do list!
    # [ ] : Update checkDB to compare number of query sequences with number of searches.
    # [Y] : Replace formatDB with makeBlastDB
    # |- [ ] : Check out makeblastdb options: title, parse_seqids, hash_index, mask_data, out, max_file_sz, logfile, taxid, taxid_map
    # [Y] : Replace formatcmd with blastdbcmd
    # [Y] : Scrap the "Run" field. Just use multiple objects for multiple runs! (New Class that manages these?)
    # [ ] : Consider adding extra info to Query table, such as top hit e-value (e.g. GABLAM HitSum).
    # [ ] : Make sure that BLAST+ does not receive duplicate commands if blastopt is used. (blastopt over-rules)
    # [Y] : Add selfhit removal option for qassemble
    # [ ] : Add Description to Search table.
    # [ ] : Add GFF3 and SAM Output.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, cyear) = ('RJE_BLAST', '2.18.0', 'October 2017', '2013')
    description = 'BLAST+ Control Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
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
formats = {'Complexity Filter':'bool','Composition Statistics':'bool','GappedBLAST':'bool','Rank':'int',
           'OneLine':'int','HitAln':'int','DBLen':'int','DBNum':'int','Length':'int','Hits':'int','FragMerge':'int',
           'BitScore':'int','E-Value':'num','GablamFrag':'int','LocalIDCut':'num','LocalCut':'int','GABLAM':'bool',
           'AlnID':'int','Expect':'num','Identity':'int','Positives':'int','QryStart':'int','QryEnd':'int',
           'SbjStart':'int','SbjEnd':'int','Ordered':'bool','Start':'int','End':'int','Len':'num','Sim':'num','ID':'num'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: BLASTRun Class                                                                                          #
#########################################################################################################################
class BLASTRun(rje_obj.RJE_Object):     
    '''
    BLASTRun Class. Author: Rich Edwards (2013). Note that these object parameters store details of the latest/current
    BLAST being performed but multiple BLAST runs can now be stored in the Database table so long as different Names
    (corresponding to the output file) are used.

    Str:str
    - Name = Output File Name (BLAST+ -out X)
    - Type = Program (blastp etc.) (BLAST+ program)
    - DBase = Database to search (BLAST+ -db X)
    - InFile = Input file name (BLAST -query X)
    - OptionFile = file containing a string of BLAST options to append to commandline
    - QConsensus = Whether to convert QAssemble alignments to consensus sequences (None/Hit/Full) [None]
    - QFasDir=PATH    : Output directory for QAssemble alignments [./QFAS/]
    - BLAST Path = path to blast programs
    - BLAST+ Path = path to blast+ programs
    - BLASTCmd = system command used to generate BLAST in self.blast()
    - BLASTOpt = Additional BLAST options
    - RefType=X       : Whether to map SAM/GFF3 hits onto the Qry, Hit or Both [Hit]

    Bool:boolean
    - Composition Statistics
    - Complexity Filter = (BLAST -F) [True]
    - BlastForce = Whether to force regeneration of new BLAST results if already existing [False]
    - FormatDB = whether to (re)format database before blasting
    - GappedBLAST = Gapped BLAST (BLAST -g X) [True]
    - GZip = Whether to (un)zip main BLAST results file if keeping [True]
    - IgnoreDate = Ignore date stamps when deciding whether to regenerate files [False]
    - OldBLAST = Whether to run in "oldblast" mode using old BLAST commands etc. [False]
    - QAssemble = Whether to fully assemble query stats from all hits [False]
    - QAssembleFas=T/F: Special mode for running with outfmt=4 and then converting to fasta file [False]
    - QComplete=T/F   : Whether the query sequence should be full-length in qassemblefas output [False]
    - RunField = Whether to inlcude Run Field in summary tables. (Useful if appending.) [False]
    - SelfSum = Whether to also include self hits in qassemble output [False] * qassemble must also be T *
    - SoftMask = Whether to use soft masking for searches [True]
    
    Int:integer
    - OneLine = Number of one-line hits per query (BLAST -v X) [500]
    - HitAln  = Number of hit alignments per query (BLAST -b X) [250]
    - DBLen = Length of Database (letters)
    - DBNum = Number of Sequences in Database
    - BLASTa = Number of processors to use (BLAST -a X) [1]
    - FragMerge=X     : Max Length of gaps between fragmented local hits to merge [0]
    - GablamFrag = Length of gaps between mapped residue for fragmenting local hits [100]
    - LocalCut = Cut-off length for local alignments contributing to global GABLAM stats) [0]

    Num:float
    - E-Value = e-value cut-off (BLAST -e X) [1e-4]
    - LocalIDCut=PERC : Cut-off local %identity for local alignments contributing to global GABLAM stats [0.0]

    List:list
    - ResTab = Whether to output summary results tables (Run/Search/Hit/Local/GABLAM) [Search,Hit]
    - SaveLocal=LIST  : Whether to generate extra output for the local BLAST hits table (GFF3/SAM) []

    Dict:dictionary
    - QAssemble = {Query:Combined GABLAM string} for QAssemble=T

    Obj:RJE_Objects
    - DB = rje_db.Database object. Contains Run, Search, Hit, PWAln and GABLAM tables. (Potentially.)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Name','Type','DBase','InFile','OptionFile','QConsensus','QFasDir','BLAST Path','BLAST+ Path','BLASTCmd','BLASTOpt','RefType','REST']
        self.boollist = ['Complexity Filter','FormatDB','Composition Statistics','GappedBLAST','IgnoreDate','GZip',
                         'OldBLAST','QAssemble','SoftMask','SelfSum','BlastForce','QAssembleFas','QComplete']
        self.intlist = ['OneLine','HitAln','DBLen','DBNum','BLASTa','FragMerge','GablamFrag','LocalCut']
        self.numlist = ['E-Value','LocalIDCut']
        self.listlist = ['ResTab','SaveLocal']
        self.dictlist = ['QAssemble']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'BLAST Path':'','BLAST+ Path':'','Type':'blastp','RefType':'Hit','QFasDir':rje.makePath('./QFAS/')})
        self.setBool({'FormatDB':False,'Composition Statistics':False,'Complexity Filter':True,'IgnoreDate':False,'SoftMask':True,
                      'GappedBLAST':True,'RunField':False})  #!# Check defaults!!
        self.setInt({'OneLine':500,'HitAln':500,'GablamFrag':100,'FragMerge':0})
        self.setNum({'E-Value':0.0001})
        self.list['ResTab'] = ['Search','Hit']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] ~ Read Commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                self._cmdRead(type='path',att='BLAST Path',arg='blastpath',cmd=cmd)
                self._cmdRead(type='path',att='BLAST+ Path',arg='blast+path',cmd=cmd)
                self._cmdRead(type='file',att='Name',arg='blasto',cmd=cmd)
                self._cmdRead(type='str',att='Type',arg='blastp',cmd=cmd.lower())
                self._cmdRead(type='str',att='Type',arg='blastprog',cmd=cmd.lower())
                self._cmdRead(type='file',att='DBase',arg='blastd',cmd=cmd)
                self._cmdRead(type='file',att='InFile',arg='blasti',cmd=cmd)
                self._cmdRead(type='bool',att='Complexity Filter',arg='blastf',cmd=cmd)
                self._cmdRead(type='bool',att='Composition Statistics',arg='blastcf',cmd=cmd)
                self._cmdReadList(cmd,'bool',['BlastForce','FormatDB','IgnoreDate','GZip','OldBLAST','QAssemble','SelfSum','SoftMask','RunField','QAssembleFas','QComplete'])
                self._cmdReadList(cmd,'int',['BLASTa','FragMerge','GablamFrag','LocalCut'])
                self._cmdRead(type='bool',att='GappedBLAST',arg='blastg',cmd=cmd)
                self._cmdRead(type='bool',att='OldBLAST',arg='legacy',cmd=cmd)
                self._cmdRead(type='float',att='E-Value',arg='blaste',cmd=cmd)
                self._cmdRead(type='int',att='OneLine',arg='blastv',cmd=cmd)
                self._cmdRead(type='int',att='HitAln',arg='blastb',cmd=cmd)
                self._cmdRead(type='int',att='OneLine',arg='tophits',cmd=cmd)
                self._cmdRead(type='int',att='HitAln',arg='tophits',cmd=cmd)
                self._cmdRead(type='file',att='OptionFile',arg='blastopt',cmd=cmd)
                self._cmdReadList(cmd,'str',['RefType','QConsensus','REST'])   # Normal strings
                self._cmdReadList(cmd,'path',['QFasDir'])  # String representing directory path
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                self._cmdReadList(cmd,'float',['LocalIDCut']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ResTab','SaveLocal'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [2] ~ Check BLAST(+) Path settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        pathfail = False
        if self.oldBLAST(): blast = 'BLAST'
        else: blast = 'BLAST+'
        if not self.getStrLC('%s Path' % blast): self.setStr({'%s Path' % blast:''})
        pathfound = os.path.exists(self.getStr('%s Path' % blast))
        oldfound = os.path.exists(self.getStr('%s Path' % blast)+'formatdb')
        newfound = os.path.exists(self.getStr('%s Path' % blast)+'makeblastdb')
        #?# Should BLAST installation detection over-ride blastpath setting #?#
        if not self.getStr('%s Path' % blast):
            if os.popen({'BLAST':'formatdb','BLAST+':'makeblastdb -help'}[blast]).read():
                self.printLog('#NCBI','Installation of %s detected. Path not required.' % blast)
            else:
                self.printLog('#NCBI','Installation of %s not detected and no path given.' % blast)
                pathfail = True
        elif not pathfound:
            self.errorLog('%s Path not found: "%s"' % (blast,self.getStr('%s Path' % blast)),printerror=False)
            pathfail = True
        elif self.oldBLAST() and not oldfound:
            if newfound: self.errorLog('BLAST path seems to point to BLAST+ programs! Use oldblast=F for BLAST+.',printerror=False)
            else: self.errorLog('%s programs not found in %s directory!' % (blast,blast),printerror=False)
            pathfail = True
        elif not self.oldBLAST() and not newfound:
            if oldfound: self.errorLog('BLAST+ path seems to point to BLAST programs! Use oldblast=T for old BLAST.',printerror=False)
            else: self.errorLog('%s programs not found in %s directory!' % (blast,blast),printerror=False)
            pathfail = True
        if pathfail and not self.getStr('%s Path' % blast): self.warnLog('Installation of %s not detected and no path given. %s functions may fail. (Check %spath=PATH/).' % (blast,blast,blast.lower()),'blastpath',quitchoice=False)
        elif pathfail: self.warnLog('Cannot execute %s functions without correct %s Path (%spath=PATH/). Bad things might happen if you proceed.' % (blast,blast,blast.lower()),'fatal',quitchoice=True)
        ### ~ [3] ~ Special ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not rje.exists(self.getStr('OptionFile')): self.setStr({'BLASTOpt':self.getStr('OptionFile')})
        if self.getBool('FormatDB'):
            if self.getStr('Type') in ['blastn','tblastn','tblastx']: self.formatDB(protein=False,force=self.force())
            else: self.formatDB(force=self.force())
        restab = string.split(string.join(self.list['ResTab']).lower())
        self.list['ResTab'] = []
        for tabtype in ['Run','Search','Hit','Local','GABLAM']:
            if tabtype.lower() in restab: self.list['ResTab'].append(tabtype)
        if self.getBool('Force') and not self.force(): self.warnLog('force=T blastforce=F. BLAST results will not be regenerated.',warntype='blastforce',suppress=True)
        if self.getNum('LocalIDCut') < 1.0: self.num['LocalIDCut'] *= 100.0
#########################################################################################################################
    def nt(self):   ### Returns whether BLAST result alignments are nucleotide
        '''Returns whether BLAST result alignments are nucleotide.'''
        return self.getStrLC('Type') == 'blastn'
#########################################################################################################################
    def force(self): return self.getBool('BlastForce')
#########################################################################################################################
    def blastPath(self):    ### Returns appropriate BLAST Path
        '''Returns appropriate BLAST Path.'''
        if self.getBool('OldBLAST') or self.getStr('BLAST+ Path').lower() in ['','none']: return self.getStr('BLAST Path')
        return self.getStr('BLAST+ Path')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self,oldblast=False,save=True,format=False,clear=True,gablam=True):  ### Main run method                            #V2.0
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if oldblast: self.setBool({'OldBLAST':True})
            if self.getBool('OldBLAST'): return rje_blast_V1.BLASTRun(self.log,self.cmd_list).blast(use_existing=not self.force(),log=True)
            if self.getStrLC('Rest'): self.restSetup()
            ## Special QAssembleFas mode
            if self.getBool('QAssembleFas'): return self.qAssembleFas()
            ## Setup tables etc.
            if not self.setup(load=not self.getStrLC('Rest'),run=True): raise ValueError
            self.list['SaveLocal'] = rje.sortUnique(rje.listLower(self.list['SaveLocal']))
            if rje.listIntersect(['sam','tdtseq','tdt'],self.list['SaveLocal']) and 'Local' not in self.list['ResTab']:
                self.list['ResTab'].append('Local')
            needtoblast = self.force()
            for table in self.list['ResTab']: needtoblast = needtoblast or not self.db(table).entryNum()
            #self.deBug('needtoblast: %s' % needtoblast)
            if needtoblast and not self.getBool('FormatDB'):
                if self.getStr('Type') in ['blastn','tblastn','tblastx']: self.formatDB(protein=False,force=self.force())
                else: self.formatDB(force=self.force())
            if needtoblast:
                if clear: self.clear()
                if self.db('Search').entryNum(): self.warnLog('Previous BLAST results not cleared. May be overwritten if same queries used.')
            ### ~ [2] ~ Main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if needtoblast and rje.checkForFile(self.getStr('DBase')) and rje.checkForFile(self.getStr('InFile')):
                self.blast(use_existing=not self.force())
            if needtoblast and rje.checkForFile(self.getStr('Name')):
                self.readBLAST(gablam=gablam,local=True,keepaln='sam' in self.list['SaveLocal'])
                for table in self.db().list['Tables']:
                    if save and table.name() in self.list['ResTab']:
                        if self.getBool('RunField') and table.name() != 'Run': table.addField('Run',evalue=self.getStr('Name')); table.list['Fields'].insert(0,table.list['Fields'].pop(-1))
                        table.dataFormat({'Rank':'int','AlnID':'int'})
                        if table.name() == 'Hit':
                            table.newKey(['Query','Rank','Hit'])
                            table.dropFields(['GablamFrag','LocalCut','GABLAM'])
                        if table.name() == 'Local':
                            if 'tdtseq' in self.list['SaveLocal']: table.saveToFile(savefields=table.fields()[:14])
                            else: table.saveToFile(savefields=table.fields()[:12])
                        else: table.saveToFile()
                        self.dict['Output'][table.name().lower()] = table.saveToFileName()
                        if self.getBool('RunField') and table.name() != 'Run': table.dropField('Run')
            if format: self.formatTables()
            ### ~ [3] ~ Add GFF/SAM output etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Save additional formats. RefType=Qry/Hit/Both to control output.
            if self.list['SaveLocal']:
                for ftype in self.list['SaveLocal']:
                    if ftype.startswith('gff'): self.saveGFF(ftype)
                    if ftype == 'sam': self.saveSAM()
        except: self.errorLog('Error during main BLASTRun.run()',quitchoice=True)
#########################################################################################################################
    def setup(self,load=True,run=False):    ### Main class setup method.                                                      #V2.0
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Name').lower() in ['','none']:
                self.setStr({'Name':'%s.%s.blast' % (rje.baseFile(self.getStr('InFile')),rje.baseFile(self.getStr('DBase'),strip_path=True))})
                self.printLog('#NAME',self.getStr('Name'))
            if self.baseFile().lower() in ['','none']: self.baseFile(rje.baseFile(self.getStr('Name')))
            db = self.db()
            if not db:
                if run: db = self.obj['DB'] = rje_db.Database(self.log,['tuplekeys=T']+self.cmd_list)
                else: db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db.baseFile(self.baseFile())
            load = load and not self.force()
            if load: self.warnLog('WARNING: Correct number handling not implemented - tables might sort oddly if reloaded')
            ## ~ [1a] BLAST Runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'Run' not in self.list['ResTab'] or not load or not db.addTable(mainkeys=['Run'],name='Run',expect=False):
                #self.deBug(self.db('Run',add=False))
                if self.db('Run',add=False): self.db('Run').clear()
                else: db.addEmptyTable('Run',['Run','Type','E-Value','DBase','InFile','BLASTCmd','Complexity Filter','Composition Statistics','SoftMask','GappedBLAST','OneLine','HitAln','DBLen','DBNum'],['Run'],log=False)
            ## ~ [1b] BLAST Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'Search' not in self.list['ResTab'] or not load or not db.addTable(mainkeys=['Query'],name='Search',expect=False):
                if self.db('Search',add=False): self.db('Search').clear()
                elif self.getBool('QAssemble'): db.addEmptyTable('Search',['Query','Length','Hits','MaxScore','TopE','Coverage','Identity','Positives'],['Query'],log=False)
                else: db.addEmptyTable('Search',['Query','Length','Hits','MaxScore','TopE'],['Query'],log=False)
            ## ~ [1c] BLAST Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'Hit' not in self.list['ResTab'] or not load or not db.addTable(mainkeys=['Query','Hit'],name='Hit',expect=False):
                if self.db('Hit',add=False): self.db('Hit').clear()
                else: db.addEmptyTable('Hit',['Query','Rank','Hit','Description','BitScore','E-Value','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Query','Hit'],log=False)
            ## ~ [1d] Local Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'Local' not in self.list['ResTab'] or not load or not db.addTable(mainkeys=['Query','Hit','AlnID'],name='Local',expect=False):
                if self.db('Local',add=False): self.db('Local').clear()
                else: db.addEmptyTable('Local',['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
                                       ['Query','Hit','AlnID'],log=False)
            ## ~ [1e] GABLAM Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'GABLAM' not in self.list['ResTab'] or not load or not db.addTable(mainkeys=['Query','Hit','QryHit'],name='GABLAM',expect=False):
                if self.db('GABLAM',add=False): self.db('GABLAM').clear()
                else: db.addEmptyTable('GABLAM',['Query','Hit','QryHit','GABLAM Start','GABLAM End','GABLAM Dirn','GABLAM Len','GABLAM Sim','GABLAM ID','GABLAM Frag',
                                                 'GABLAMO Start','GABLAMO End','GABLAMO Dirn','GABLAMO Len','GABLAMO Sim','GABLAMO ID','GABLAMO Frag'],
                                       ['Query','Hit','QryHit'],log=False)
            ### ~ [2] Clean up Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for table in self.db().tables():
                if 'Run' in table.fields() and table.name() != 'Run': table.dropField('Run')
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        The standard REST call is in the form: `blast&blasti=FASFILE&blastd=FASFILE&blastprog=X`. By default, results
        will be parsed into summary tables. If `qassemblefas=T`, the summary tables will be replaced with alignments of
        each query with its hits. Each local alignment is a separate sequence in the alignment unless `qconsensus=X` is
        used to convert QAssemble alignments to consensus sequences.

        * `qconsensus=Hit` converts QAssemble alignments to consensus sequences by Hit sequence.
        * `qconsensus=Full` converts all QAssemble alignments to a single consensus sequence.

        When `qconsensus=X` is used, the most abundant amino acid or nucleotide at each position is used. In the case of
        a tie, the query sequence is used if it's one of the options, else the highest ranked one is used.

        Run with &rest=docs for program documentation and more options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        After running, click on the `blastres` tab to see overall BLAST search results. The following tables and sequence
        file will usually also be output:

        ### Available REST Outputs:
        blastres = Main BLAST results output file
        search = Delimited summary of each search [tab]
        hit = Delimited summary of all hit [tab]
        local = Delimited table of all local BLAST alignments [tab]
        blasti = BLAST search input (queries) [fas]
        blastd = BLAST search database [fas]
        qassembly = Assembled QAssembly hits as aligned fasta for first query [fas]

        NOTE: If run in `&qassemblefas=T` mode, the `search`, `hit` and `local` tables will not be output. The first
        query will have QAssembly output in the `qassembly` tab. All queries will have also their own tab with their
        qassembly alignment, named after their accession numbers.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder():
                if outfmt not in self.dict['Output'] or not self.dict['Output'][outfmt]:
                    self.dict['Output'][outfmt] = 'No output generated.'
            self.dict['Output']['blastres'] = 'Name'
            self.dict['Output']['blasti'] = 'InFile'
            self.dict['Output']['blastd'] = 'DBase'
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        outlist = ['blastres']
        if self.getBool('QAssembleFas'): outlist.append('qassembly')
        else:
            for table in self.list['ResTab']: outlist.append(table.lower())
        outlist += ['blasti','blastd']
        for outfmt in rje.sortKeys(self.dict['Output']):
            if outfmt not in outlist: outlist.append(outfmt)
        if self.dev(): self.printLog('#REST',string.join(outlist,'; '))
        return outlist
#########################################################################################################################
    ### <3> ### General Class Data Methods                                                                              #
#########################################################################################################################
    def formatTables(self): ### Formats search results stored in tables for other programs to use properly          #V2.0
        '''Formats search results stored in tables for other programs to use properly.'''
        for table in self.db().list['Tables']: table.dataFormat(formats)
#########################################################################################################################
    def clear(self): self.setup(load=False)   # Currently, this will generate a new Database object full of empty tables.     #V2.0
#########################################################################################################################
    def addDBRun(self):     ## Adds DBRun using current option settings                                             #V2.0
        '''Adds DBRun using current option settings.'''
        dbrun = {}
        for add_dict in [self.str,self.bool,self.int,self.num]: rje.combineDict(dbrun,add_dict)
        self.db('Run').addEntry(dbrun)
#########################################################################################################################
    def name(self): return self.getStr('Name')
    def legacy(self): return self.getBool('OldBLAST')
    def oldBLAST(self): return self.getBool('OldBLAST')
#########################################################################################################################
    def queries(self): return self.db('Search').dataKeys()[0:]
    def searchNum(self):   ### Returns number of Searches                                                           #V2.4
        '''Returns number of Searches for specific run or all if no run given.'''
        return self.db('Search').entryNum()
#########################################################################################################################
    def hitNum(self,query=None):    ### Returns number of Hits for specific query or all if no query given          #V2.4
        '''Returns number of Hits for specific run or all if no run given.'''
        if query: return self.db('Search').data(query)['Hits']
        else: return self.db('Hit').entryNum()
#########################################################################################################################
    def queryHits(self,query=None): ### Returns a list of hits for a given Query                                    #V2.4
        '''
        Returns a list of hits for a given Query.
        >> query:str [None] = BLAST Search Query. If None will return all hits.
        << Returns a list of hits sorted by BLAST rank.
        '''
        try:
            qhits = []
            for rank in rje.sortKeys(self.db('Hit').index('Rank',log=False)):
                for hentry in self.db('Hit').indexEntries('Rank',rank):
                    if (hentry['Query'] == query or not query) and hentry['Hit'] not in qhits: qhits.append(hentry['Hit'])
            return qhits
        except: return []
#########################################################################################################################
    def hitData(self,hit,query): ### Returns entry from Hit dictionary                                              #V2.4
        '''
        Returns entry from Hit dictionary.
        >> hit:str = BLAST Hit
        >> query:str = BLAST Search Query
        '''
        try:
            hkey = self.db('Hit').makeKey({'Hit':hit,'Query':query})
            return self.db('Hit').data(hkey)
        except: return {}
#########################################################################################################################
    def gablamData(self,hit,query):    ### Returns ['Query','Hit']:{GABLAM} dictionary                              #V2.4
        '''
        Returns ['Query','Hit']:{GABLAM} dictionary.
        >> hit:str = BLAST Hit
        >> query:str = BLAST Search Query
        '''
        gdb = self.db('GABLAM')
        qkey = gdb.makeKey({'Hit':hit,'Query':query,'QryHit':'Query'})
        hkey = gdb.makeKey({'Hit':hit,'Query':query,'QryHit':'Hit'})
        return {'Query':rje.combineDict({},gdb.data()[qkey]),'Hit':rje.combineDict({},gdb.data()[hkey])}
#########################################################################################################################
    def localData(self,hit,query): ### Returns Local data entries as list in Alignment order                        #V2.4
        '''
        Returns Local data entries as list in Alignment order.
        >> hit:str = BLAST Hit
        >> query:str = BLAST Search Query
        '''
        try:
            localdata = []; sortdata = {}
            for lkey in rje.listIntersect(self.db('Local').index('Query',log=False)[query],self.db('Local').index('Hit',log=False)[hit]):
                entry = self.db('Local').data(lkey)
                sortdata[entry['AlnID']] = entry
            for lkey in rje.sortKeys(sortdata): localdata.append(sortdata[lkey])
            return localdata
        except:
            self.errorLog('Problem with localData(%s vs %s)' % (query,hit))
            return []
#########################################################################################################################
    def saveHitIDs(self,outfile=None,search=None,noempty=True,backups=True):   ### Saves IDs of hits in file (e.g. for later fastacmd extraction)
        '''
        Saves IDs of hits in file (e.g. for later fastacmd extraction).
        >> outfile:str = outfile to use.
        >> search:str = optional Search ID
        >> noempty:bool [True] = Whether to create empty files (False) or skip (True)
        '''
        try:### ~ [1] ~ Save simple list of IDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfile == None: outfile = '%s.id' % self.baseFile()
            hdb = self.db('Hit')
            idlist = self.queryHits(search)
            if not idlist and noempty: return -1
            if backups: rje.backup(self,outfile)
            elif os.path.exists(outfile): os.unlink(outfile)
            open(outfile, 'a').write('%s\n' % string.join(idlist,'\n'))
            self.dict['Output']['id'] = outfile
            return len(idlist)
        except: self.errorLog('Major error during BLASTRun.saveHitIDs().'); return -1
#########################################################################################################################
    ### <4> ### BLAST Search Methods                                                                                    #
#########################################################################################################################
    def blast(self,wait=True,type=None,cleandb=False,use_existing=False,log=True):    ### Performs BLAST using object attributes
        '''
        Performs BLAST using object attributes.
        >> wait:boolean  = whether to wait for BLAST. [True]
        >> type:str = type of BLAST search [None]
        >> cleandb:bool = whether to cleanup (delete) searchDB files after search [False]
        >> use_existing:bool = if True, will check for existing result and use if newer than files
        >> log:bool = Whether to log BLAST run
        '''
        try:### ~ [1] ~ Setup BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Check for Existing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if use_existing and self.checkBLAST(logcheck=log):                ### BLAST Results exist
                if self.getBool('IgnoreDate'): return True       ### Don't check age!
                if rje.isYounger(self.getStr('DBase'),self.getStr('Name')) == self.getStr('Name') and rje.isYounger(self.getStr('InFile'),self.getStr('Name')) == self.getStr('Name'): return True
                #self.debug(rje.isYounger(self.getStr('DBase'),self.getStr('Name')) == self.getStr('Name'))
                #self.debug(rje.isYounger(self.getStr('InFile'),self.getStr('Name')) == self.getStr('Name'))
                #self.debug('Checked but bad date')
            ## ~ [1b] ~ Setup BLAST Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            blastpath = self.blastPath()
            if self.getInt('OneLine') < self.getInt('HitAln'):
                self.printLog('#CMD','Reduced reported alignments from %s to match one-line reporting number of %s.' % (rje.integerString(self.getInt('HitAln')),rje.integerString(self.getInt('OneLine'))))
                self.int['HitAln'] = self.int['OneLine']
            ### ~ [2] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.getStr('Name'),log=True)
            if type: self.str['Type'] = type
            ## ~ [2a] ~ Make system command for BLAST call ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('OldBLAST'):
                command = blastpath + 'blastall'
                command += ' -p %s' % self.getStr('Type')
                command += ' -i %s' % self.getStr('InFile')
                command += ' -d %s' % self.getStr('DBase')
                command += ' -o %s' % self.getStr('Name')
                command += ' -F %s' % str(self.getBool('Complexity Filter'))[0] 
                try: command += ' -C %s' % str(self.getBool('Composition Statistics'))[0]
                except: pass    # For pickled BLASTs
                command += ' -e %e' % self.getNum('E-Value')
                command += ' -v %d' % self.getInt('OneLine')
                command += ' -b %d' % min(self.getInt('HitAln'),self.getInt('OneLine'))
                if not self.getBool('GappedBLAST'): command += ' -g F'
                if self.getInt('BLASTa') > 1: command += ' -a %d' % self.getInt('BLASTa')
            else:
                command = blastpath + self.getStr('Type')
                command += ' -query %s' % self.getStr('InFile')
                command += ' -db %s' % self.getStr('DBase')
                command += ' -out %s' % self.getStr('Name')
                if self.getStr('Type') in ['blastn']:
                    if self.getBool('Complexity Filter'): self.warnLog('Cannot use Complexity Filter "seg" option in %s' % self.getStr('Type'),'segwarn',suppress=True)
                    if self.getBool('Composition Statistics'): self.warnLog('Cannot use Composition Statistics "comp_based_stats" option in %s' % self.getStr('Type'),'comp_based_statswarn',suppress=True)
                else:
                    if self.getBool('Complexity Filter'): command += ' -seg yes'
                    else: command += ' -seg no'
                    try: command += ' -comp_based_stats %s' % str(self.getBool('Composition Statistics'))[0]
                    except: pass    # For pickled BLASTs
                command += ' -soft_masking %s' % str(self.getBool('SoftMask')).lower()
                command += ' -evalue %e' % self.getNum('E-Value')
                command += ' -num_descriptions %d' % max(self.getInt('HitAln'),self.getInt('OneLine'))
                num_alignments = min(self.getInt('HitAln'),self.getInt('OneLine'))
                #if num_alignments: command += ' -num_alignments %d' % (max(num_alignments,2))   #!# WHY?! Was there a BLAST bug?
                if num_alignments > 0: command += ' -num_alignments %d' % num_alignments
                else: command += ' -num_alignments 0'
                if not self.getBool('GappedBLAST'): command += ' -ungapped'
                if self.getInt('BLASTa') > 1: command += ' -num_threads %d' % self.getInt('BLASTa')
            if self.getStr('BLASTOpt').lower() not in ['','none']: command = '%s %s' % (command,self.getStr('BLASTOpt'))
            if rje.exists(self.getStr('OptionFile')):
                for line in open(self.getStr('OptionFile'),'r').readlines(): command = '%s %s' % (command,rje.chomp(line))
            if not wait: command += ' &'
            ##~ [2b] ~ Perform BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.str['BLASTCmd'] = command
            if log: self.printLog('\r#SYS',command)
            elif self.v() > 1: self.printLog('\r#SYS',command,log=False)
            os.system(command)
            ## ~ [2c] ~ Cleanup database ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if cleandb: cleanupDB(self,self.getStr('DBase'))
        except: self.errorLog('Fatal Error during BLASTRun.blast()'); raise
#########################################################################################################################
    def checkBLAST(self,resfile=None,logcheck=True,expect=False,checkeval=True):  ### Checks that each BLAST started has an end                 #V2.0
        '''
        Checks that each BLAST started has an end, thus identifying BLAST runs that have been terminated prematurely and
        need to be re-run.
        >> resfile:str = Results File (set as self.info['Name'])
        >> logcheck:boolean = Whether to print findings to log [True]
        >> expect:bool [False] = Whether complete BLAST results file is expected to exist.
        >> checkeval:bool [True] = Whether to check that read E-value is <= current setting.
        << Returns number of Searches. (None if format etc. failed)
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            searchx = 0
            if resfile != None: self.str['Name'] = resfile
            if not os.path.exists(self.str['Name']):
                if logcheck and expect: self.printLog('#BLAST','Results file "%s" not found!' % self.str['Name'])
                return None
            RESFILE = open(self.str['Name'],'r');
            RESFILE.seek(0,2); fend = RESFILE.tell(); RESFILE.seek(0)
            line = RESFILE.readline()
            if not line:
                if logcheck: self.errorLog('%s found but no content!' % self.str['Name'],printerror=False)
                return None
            ### ~ [2] ~ Check Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtxt = 'Checking %s BLAST results' % self.str['Name']
            px = 0.0; self.progLog('\r#BLAST','%s: %.2f%%' % (logtxt,px/fend))
            last = 'Top'    # Need at least one search to be intact
            while line:
                self.progLog('\r#BLAST','%s: %.2f%%' % (logtxt,px/fend))
                if line.startswith('Query=') and last in [None,'Top']:
                    last = line   # New Search
                    searchx += 1
                    px = 100.0 * RESFILE.tell(); self.progLog('\r#BLAST','%s: %.2f%% (%s searches)' % (logtxt,px/fend,rje.iStr(searchx)))
                elif line.startswith('Query='): self.deBug('\n%s\n??\n%s\n' % (last,line))
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if last.find('Query=') != 0:    # Wrong place!
                        if logcheck: self.errorLog('%s check: Found "Query=" in wrong place!' % self.str['Name'],printerror=False)
                        return None
                    last = line
                elif line.find('***** No hits found *****') >= 0:  # No Hits
                    if last.find('Query=') != 0:    # Wrong place!
                        if logcheck: self.errorLog('%s check: Found "Query=" in wrong place!' % self.str['Name'],printerror=False)
                        return None
                    last = line
                elif line.find('Effective length') >= 0 or line.find('Effective search space') >= 0 or line.find('Gap Penalties:') >= 0 or line[:5] in ['BLAST','Matri'] or 'BLAST' in string.join(string.split(line)[:1]):
                    if not last: pass
                    elif last.find('Sequences producing significant alignments:') >= 0: last = None     # One-line hits
                    elif last.find('***** No hits found *****') >= 0: last = None                      # No Hits
                    #else: self.deBug('\n%s\n>>\n%s\n' % (last,line))
                elif checkeval and line.find('Number of sequences better than') >= 0:
                    search_eval = string.atof(rje.matchExp('Number of sequences better than\s+(\S+):', line)[0])
                    if self.num['E-Value'] > search_eval:   # Existing BLAST used a more stringent e-value
                        if self.i() < 0 or rje.yesNo('%s used blaste=%s. Current setting is less stringent blaste=%s. Reject BLAST results?' % (self.str['Name'],search_eval,self.num['E-Value'])):
                            last = 'Eval'; break
                        else: self.warnLog('%s used blaste=%s. Current setting is less stringent blaste=%s.' % (self.str['Name'],search_eval,self.num['E-Value']))
                line = RESFILE.readline()
            ### ~ [3] ~ Finish check and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            RESFILE.close()
            if last == 'Eval': self.printLog('\r#BLAST','%s: BLAST search too stringent (e=%s).' % (logtxt,search_eval),log=logcheck); return 0
            elif last: self.printLog('\r#BLAST','%s: BLAST incomplete?! (%s searches)' % (logtxt,rje.iStr(searchx)),log=logcheck); return 0
            else: self.printLog('\r#BLAST','%s: BLAST intact (%s searches).' % (logtxt,rje.iStr(searchx)),log=logcheck); return searchx
        except:
            try: RESFILE.close()
            except: pass
            self.errorLog('Error during BLASTRun.checkBLAST().'); return 0
#########################################################################################################################
    def checkProg(self,qtype='Unknown',stype='Unknown',log=True): ### Checks sequence types against blastp          #V2.0
        '''
        Checks sequence types against blastp.
        >> qtype:str [None] = Query sequence type
        >> stype:str [None] = SearchDB sequence type
        >> log:bool [True] = Whether to output result of check to log
        '''
        try:### ~ [1] ~ Check BLASTP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qtype = qtype.lower(); stype = stype.lower()
            problem = False
            if self.str['Type'] in ['blastp','blastx'] and stype[1:] == 'na': problem = True
            elif self.str['Type'] in ['blastp','tblastn'] and qtype[1:] == 'na': problem = True
            elif self.str['Type'] in ['blastn','tblastn','tblastx'] and stype[:4] in ['prot','pept']: problem = True
            elif self.str['Type'] in ['blastn','blastx','tblastx'] and qtype[:4] in ['prot','pept']: problem = True
            if not problem:
                if log: self.printLog('#PROG','BLAST Program OK (%s): %s vs %s' % (self.str['Type'],qtype,stype))
                return True
            ### ~ [2] ~ Predict BLASTP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blastp = ''
            if qtype[:4] in ['prot','pept']:
                if stype[:4] in ['prot','pept']: blastp = 'blastp'
                elif stype[1:] == 'na': blastp = 'tblastn'
            elif qtype[1:] == 'na':
                if stype[:4] in ['prot','pept']: blastp = 'blastx'
                elif stype[1:] == 'na': blastp = 'blastn'
            if blastp:
                if self.i() < 0 or rje.yesNo('%s vs %s: switch BLAST program to %s?' % (qtype,stype,blastp)):
                    self.str['Type'] = blastp
                    if log: self.printLog('#PROG','BLAST Program changed (%s): %s vs %s' % (self.str['Type'],qtype,stype))
                    return True
            else: blastp = 'Unknown'
            self.printLog('#ERR','WARNING: %s vs %s and blastp=%s (Recommended: %s)' % (qtype,stype,self.str['Type'],blastp))
            return False
        except: self.errorLog('Something went wrong with blast.checkProg()'); return False
#########################################################################################################################
    def readBLAST(self,resfile=None,clear=False,gablam=True,unlink=False,local=False,screen=True,log=False,keepaln=False):        #V2.0
        '''
        Reads BLAST Results into objects.
        >> resfile:str = Results File (set as self.info['Name'])
        >> clear:Boolean = whether to clear current searches (True) or just append (False) [False]
        >> gablam:Boolean = whether to calculate gablam statistics and clear alignments to save memory [True]
        >> unlink:Boolean = whether to delete BLAST results file after reading [False]
        >> local:Boolean = whether to store Local alignment dictionary with basic alignment data [False]
        >> screen:Bool [False] = whether to output reading details to screen
        >> log:Bool [False] = whether to output reading details to log
        >> keepaln:Bool [False] = whether to keep local alignment strings when calculating GABLAM.
        << returns True if (apparently) read to completion OK, else False
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if resfile != None: self.str['Name'] = resfile
            if clear or not self.db(): self.clear()
            if not resfile: resfile = self.str['Name']
            RESFILE = open(resfile,'r')
            RESFILE.seek(0,2) 
            fend = RESFILE.tell()
            presx = self.searchNum(); prehx = self.hitNum()
            RESFILE.seek(0); fpos = 0
            ## ~ [1a] Check BLAST Type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newtype = string.split(RESFILE.readline())[0].lower()
            if newtype != self.str['Type'].lower():
                if newtype in ['blastn','blastp','blastx','tblastn','tblastx','rpsblast','rpstblastn']:
                    self.str['Type'] = newtype
                    self.printLog('#BLAST','BLAST type changed to "%s" to match results file.' % newtype,screen=screen,log=log)
                else: self.errorLog('Do not recognise BLAST type "%s". Will try to read anyway.' % newtype)
            ## ~ [1b] Make sure Run in Run Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rdb = self.db('Run'); sdb = self.db('Search'); hdb = self.db('Hit')
            if resfile not in rdb.data(): rdb.addEntry({'Run':resfile,'Type':newtype})
            rentry = rdb.data(resfile)
            line = RESFILE.readline()
            while not rje.matchExp('^Database:\s+(\S.+)$',line): line = RESFILE.readline()
            rentry['DBase'] = rje.matchExp('^Database:\s+(\S.+)$',line)[0]
            line = RESFILE.readline()
            while not rje.matchExp('^\s+(\d+) sequences; (\d+) total',string.replace(line,',','')):
                rentry['DBase'] += rje.chomp(line)
                line = RESFILE.readline()
            (rentry['DBNum'],rentry['DBLen']) = rje.matchExp('^\s+(\d+) sequences; (\d+) total',string.replace(line,',',''))
            ## ~ [1c] Update Search Table for QAssemble=T ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('QAssemble'):
                for field in ['Coverage','Identity','Positives']:
                    if field not in sdb.fields(): sdb.addField(field)
            ### ~ [2] ~ Read in Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtxt = 'Reading %s BLAST results' % (newtype)
            RESFILE.seek(0)
            sentry = True
            while sentry:
                (sentry,fpos) = self.readNextBLASTSearch(fpos,resfile,gablam,local,RESFILE,keepaln)
                self.progLog('\r#BLAST','%s %.2f%%: %s searches; %s hits' % (logtxt,100.0*fpos/fend,rje.integerString(self.searchNum()-presx),rje.integerString(self.hitNum()-prehx)),screen=screen)
            logtxt = 'Reading %s %s BLAST results' % (resfile,newtype)
            ptxt = '%s complete: %s searches; %s hits' % (logtxt,rje.integerString(self.searchNum()-presx),rje.integerString(self.hitNum()-prehx))
            self.printLog('\r#BLAST','%s (vs %s Seq; %s letters)' % (ptxt,rje.iStr(rentry['DBNum']),rje.iStr(rentry['DBLen'])),log=log,screen=screen)
            ### ~ [5] ~ Unlink if necessary and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            RESFILE.close()
            #self.bugPrint(self.str['Name'])
            #self.debug('Exist:%s; Unlink:%s' % (os.path.exists(self.str['Name']),unlink))
            if unlink and os.path.exists(self.str['Name']):
                if self.i() < 1 or rje.yesNo('Delete results file: %s?' % self.str['Name']):
                    os.unlink(self.str['Name'])
                    self.printLog('#BLAST','Deleted BLAST results file: %s' % self.str['Name'])
            return True
        except:
            if self.checkBLAST():
                self.errorLog('Fatal Error during BLASTRun.readBLAST() despite intact BLAST results.')
                raise
            else: self.errorLog('BLASTRun.readBLAST() failed due to incomplete BLAST results file.')
        return False
#########################################################################################################################
    def readNextBLASTSearch(self,fpos=0,resfile=None,gablam=False,local=False,OPENRES=None,keepaln=False):  ### Reads BLAST Result into objects
        '''
        Reads BLAST Results into objects.
        >> fpos:Integer = position in results file to start reading from.
        >> resfile:str = Results File
        >> gablam:Boolean = whether to calculate gablam statistics and clear alignments to save memory [False]
        >> local:Boolean = whether to store Local alignment dictionary with basic alignment data [False]
        >> OPENRES:FILE = Open BLAST results file. [None]
        >> keepaln:Bool [False] = whether to keep local alignment strings when calculating GABLAM.
        << returns tuple: (Search table entry object or None if no more results to be read,fpos).
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if fpos < 0: return (None,-1)
            if not resfile: resfile = self.str['Name']
            if OPENRES: RESFILE = OPENRES
            else: RESFILE = open(resfile,'r')
            RESFILE.seek(fpos)
            hdb = self.db('Hit')
            ### ~ [2] ~ Read in Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            search_eval = 0.0; eval_hitnum = 0
            line = RESFILE.readline()
            search = None; readhits = False; hitalnx = 0; i = 0; hits = []; hitnames = [];
            hitaln = {} # Dictionary of {hit:[aln dictionaries]}
            while line:
                while i > 0 and line: fpos = RESFILE.tell(); line = RESFILE.readline(); i -= 1
                if not line: break
                ## ~ [2a] ~ Basic Search Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.find('Query=') == 0:    # New Search
                    if search: break            # Read one whole search, so stop and return
                    search = self.db('Search').addEntry({'Query':rje.matchExp('^Query=\s+(\S+)', line)[0],'Length':0,'Hits':0,'MaxScore':0.0,'TopE':self.getNum('E-Value')})
                elif rje.matchExp('^\s+\((\S+) letters\)',line):
                    len_match = string.replace(rje.matchExp('^\s+\((\S+) letters\)',line)[0],',','')
                    search['Length'] = string.atoi(len_match)
                elif rje.matchExp('^Length=(\S+)',line) and not eval_hitnum:
                    len_match = string.replace(rje.matchExp('^Length=(\d+)',line)[0],',','')
                    search['Length'] = string.atoi(len_match)
                elif line.find('Number of letters in database:') >= 0:
                    dblen = rje.matchExp('Number of letters in database:\s+(\d\S*)', line)[0]
                    dblen = re.sub('\D','',dblen)
                    self.int['DBLen'] = string.atoi(dblen)
                elif line.find('Number of sequences in database:') >= 0:
                    dbnum = rje.matchExp('Number of sequences in database:\s+(\d\S*)', line)[0]
                    dbnum = re.sub('\D','',dbnum)
                    self.int['DBNum'] = string.atoi(dbnum)
                elif rje.matchExp('(\S+) sequences; (\S+) total letters',line):
                    (dbnum,dblen) = rje.matchExp('(\S+) sequences; (\S+) total letters',line)
                    self.int['DBNum'] = string.atoi(re.sub('\D','',dbnum))
                    self.int['DBLen'] = string.atoi(re.sub('\D','',dblen))
                elif line.find('Number of sequences better than') >= 0:
                    search_eval = string.atof(rje.matchExp('Number of sequences better than\s+(\S+):', line)[0])
                    self.num['E-Value'] = min(self.num['E-Value'],search_eval)
                ## ~ [2b] ~ One-line hit data (BLASTHit) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if string.split(line)[-1] == 'N': (si,ei) = (-3,-2)
                    else: (si,ei) = (-2,-1)
                    rankx = 0
                    RESFILE.readline(); fpos = RESFILE.tell(); line = RESFILE.readline()   # Skip blank line
                    while rje.matchExp('^(\S+)\s.*\s(\S*\d)\s+(\S*\d)\s*$',line):
                        rankx += 1
                        match = string.split(string.replace(line,'lcl|',''))
                        #self.bugPrint(match[0])
                        evalue = match[ei]
                        if evalue.find('e') == 0: evalue = '1' + evalue
                        hite = float(evalue)
                        if hite <= self.num['E-Value']:     # Hit is OK
                            hit = {'Query':search['Query'],'Rank':rankx,'Hit':match[0],'Description':string.join(string.split(line)[1:si]),'BitScore':string.atof(match[si]),'E-Value':hite,'Length':0,'Aln':0}
                            if rankx == 1: search['TopE'] = hit['E-Value']; search['MaxScore'] = hit['BitScore']
                            hitkey = hdb.makeKey(hit)
                            if hitkey in hdb.data(): self.errorLog('Sequence "%s" hit by %s multiple times! Check that sequence names are unique in database.' % (match[0],search['Query']),printerror=False)
                            hit = hdb.addEntry(hit)
                            hits.append(hit); search['Hits'] = len(hits)
                            hitnames.append(hit['Hit'])
                            hitaln[hit['Hit']] = []
                        eval_hitnum += 1
                        fpos = RESFILE.tell(); line = RESFILE.readline()
                    hitalnx = 0; readhits = True
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    hitalnx = 0; readhits = True
                ## ~ [2c] ~ Aln Hit data (PWAln) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('>') == 0:
                    if not readhits: i += 1; continue   # '>' Character within name happened to occur at start of line
                    hitname = rje.matchExp('^>(\S+)\s*',string.replace(line,'lcl|',''))[0]
                    #try: self.bugPrint('%s: %d (%d)' % (hitname,hitnames.index(hitname),hitalnx))
                    #except: self.deBug('%s: %s' % (hitname,hitname in hitnames))
                    if hitalnx >= len(hits) and hitalnx >= eval_hitnum:
                        self.errorLog('Apparent hits for %s have exceeded the %d found. (Hit %s.) BLAST read-through or sequence name format error?' % (search['Query'],len(hits),hitname),False,False)                       
                    if hitalnx >= len(hits) or hitname != hits[hitalnx]['Hit']:      # Identify hit object
                        for hit in hits:
                            if hit['Hit'] == hitname: hitalnx = hits.index(hit)
                    if hitalnx >= len(hits) or hitname != hits[hitalnx]['Hit']:
                        if hitalnx >= eval_hitnum:
                            self.errorLog('Problem with BLAST results - %s single-line hits and alignments do not match' % search['Query'],printerror=False,quitchoice=True)
                        aln = None
                        i += 1; continue
                    hit = hits[hitalnx]
                    hit['Description'] = string.join(string.split(line)[1:])
                    hitalnx += 1
                    aln = None
                    while line:
                        fpos = RESFILE.tell(); line = RESFILE.readline()
                        if not line: break
                        ## Hit Length ##
                        #self.bugPrint('%s >> %s' % (rje.chomp(line),rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line)))
                        if rje.matchExp('^\s*Length\s*=\s*(\d+)',line):
                            hit['Length'] = string.atoi(rje.matchExp('^\s*Length\s*=\s*(\d+)',line)[0])
                        ## New Aln Block ##
                        elif hit['Length'] and rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line):
                            #self.deBug('%s' % aln)
                            if aln:
                                hitaln[hit['Hit']].append(aln); hit['Aln'] += 1
                                if local: self.db('Local').addEntry(aln)
                            scores = rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line)
                            evalue = scores[1]
                            if evalue.find('e') == 0: evalue = '1' + evalue
                            aln = {'Query':search['Query'],'Hit':hit['Hit'],'AlnID':hit['Aln']+1,'BitScore':string.atof(scores[0]),'Expect':string.atof(evalue),
                                   'Length':0,'Identity':0,'Positives':0,'QryStart':-1,'QryEnd':-1,'SbjStart':-1,'SbjEnd':-1,'QrySeq':'','SbjSeq':'','AlnSeq':''}
                            fpos = RESFILE.tell(); line = RESFILE.readline()
                            if re.search('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives = (\d+)/(\d+)\s?',line):
                                sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives\s+=\s+(\d+)/(\d+)\s?',line)
                                aln['Positives'] = string.atoi(sim[2])
                            else: sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?',line)
                            try:
                                aln['Length'] = string.atoi(sim[1])
                                aln['Identity'] = string.atoi(sim[0])
                            except: self.errorLog('Problem reading line "%s"' % line)
                            fpos = RESFILE.tell(); line = RESFILE.readline()
                            if line.find('Frame') < 0: RESFILE.seek(fpos)    # Not BLASTX or TBLASTN
                        ## Alignment ##
                        elif line.find('Query') == 0 and aln != None:
                            # Query Line
                            alnstructure = rje.matchExp('^(Query:*\s+\d+\s+)(\S+)\s+(\d+)',line)
                            if not alnstructure: alnstructure = rje.matchExp('^(Query:*\s+)(-+)',line) + ('%d' % aln['QryEnd'],)
                            if not alnstructure:    # Something's wrong!
                                self.errorLog('Problem reading aln query for %s hit %s.' % (search['Query'],hit['Hit']),printerror=False)
                                rje.combineDict(aln,{'QryStart':-1,'QryEnd':-1,'SbjStart':-1,'SbjEnd':-1,'QrySeq':'','SbjSeq':'','AlnSeq':''},copyblanks=True)
                                break
                            leader = len(alnstructure[0])
                            aln['QrySeq'] += alnstructure[1]
                            aln['QryEnd'] = string.atoi(alnstructure[2])
                            if aln['QryStart'] <= 0:
                                aln['QryStart'] = string.atoi(rje.matchExp('^Query:*\s+(\d+)\s',alnstructure[0])[0])
                            # Alignment Line
                            fpos = RESFILE.tell(); line = RESFILE.readline(); aln['AlnSeq'] += line[leader:(leader+len(alnstructure[1]))]
                            # Subject Line
                            fpos = RESFILE.tell(); line = RESFILE.readline(); subject = rje.matchExp('^Sbjct:*\s+(\d+)\s+(\S+)\s+(\d+)',line)
                            if not alnstructure: alnstructure = rje.matchExp('^(Sbjct:*\s+)(-+)',line) + ('%d' % aln['SbjEnd'],)
                            if not subject:
                                self.errorLog('Problem reading aln subject for %s hit %s.' % (search['Query'],hit['Hit']),printerror=False)
                                rje.combineDict(aln,{'QryStart':-1,'QryEnd':-1,'SbjStart':-1,'SbjEnd':-1,'QrySeq':'','SbjSeq':'','AlnSeq':''},copyblanks=True)
                                break
                            aln['SbjSeq'] += subject[1]
                            aln['SbjEnd'] = string.atoi(subject[2])
                            if aln['SbjStart'] <= 0: aln['SbjStart'] = string.atoi(subject[0])
                        ## End of Aln ##
                        elif line.find('>') == 0:
                            RESFILE.seek(fpos)  # Return to beginning of this line
                            break   # so '>' will still be caught as next Aln
                        ## End of All Alns ##
                        elif line.find('Database:') >= 0 and hit['Length']: break   # hit.alnNum() to avoid matching 'Database:' in hit description
                        elif line.find('Effective') >= 0 and hit['Length']: break   # hit.alnNum() to avoid matching 'Database:' in hit description
                        elif line.find('BLAST') == 0: break
                        elif line.find('TBLASTN') == 0: break
                        elif not hit['Length']: hit['Description'] = string.join([hit['Description']]+string.split(line)[0:])
                    if aln: 
                        hitaln[hit['Hit']].append(aln); hit['Aln'] += 1
                        if local: self.db('Local').addEntry(aln)
                ## ~ [2d] ~ Continue reading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##      
                i += 1
            if not OPENRES: RESFILE.close()
            ### ~ [3] ~ Filter on E-Value if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if search_eval > self.num['E-Value']:   # Filter any hits that are too weak
                pass    #!# Add filtering of hits here!
                #!# Renamed eval above and limited update of Search, Hit and Local tables
                #!# Possibly add additional BLAST filters here?
                #hdb = self.db('Hit')
                # Update search = self.db('Search').addEntry({'Query':rje.matchExp('^Query=\s+(\S+)', line)[0],'Length':0,'Hits':0,'MaxScore':0.0,'TopE':self.getNum('E-Value')})
                #self.db('Local')
            ### ~ [4] ~ GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gablam and search:
                for hit in hits:
                    gablam = self.globalFromLocal(search,hit,hitaln[hit['Hit']],keepaln=keepaln)     # This will also clear the alignment data.
                    for qh in ['Query','Hit']:
                        self.db('GABLAM').addEntry(rje.combineDict({'Query':search['Query'],'Hit':hit['Hit'],'QryHit':qh},gablam[qh]))
                ## ~ [4a] ~ QAssemble global stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getBool('QAssemble'):
                    search['Coverage'] = search['Identity'] = search['Positives'] = 0
                    if hits:
                        qrylen = search['Length']
                        try:
                            qgablam = self.dict['QAssemble'].pop(search['Query'])
                            # Update search (entry) dictionary
                            search['Coverage'] = qrylen - string.count(qgablam,'-')
                            search['Identity'] = string.count(qgablam,'|')
                            search['Positives'] = search['Identity'] + string.count(qgablam,'+')
                        except: self.debug(self.dict['QAssemble'])
            return (search,fpos)
        except:
            if not OPENRES: 
                try: RESFILE.close()
                except: pass
            self.errorLog('Error during BLASTRun.readNextBLASTSearch()'); raise
#########################################################################################################################
    def globalFromLocal(self,search,hit,localdata=None,keepaln=False,log=False):      ### Returns a dictionary of global alignment stats from Query-Hit local alignments
        '''
        Returns a dictionary of global alignment stats from Query-Hit local alignments. During this method, strings are
        made for Query and Hit sequences. These are "-" where there is no BLAST match, "X" over indels and mismatches,
        "|" for matches and "+" for positives (protein only).
        >> search:dict = entry in Search Table
        >> hit:dict = entry in Hit Table
        >> localdata:list = list of entries from Local table
        >> keepaln:bool [False] = Whether to store GABLAM alignments in Hit object
        << {Query:{},Hit:{}}, where each value is a dictionary of [GABLAM(O) ID, GABLAM(O) Sim, GABLAM(O) Len]
        '''
        try:
            ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qrylen = search['Length']
            query = search['Query']
            if not localdata:
                #self.debug('GABLAM %s vs hit %s' % (search['Query'],hit['Hit']))
                if search['Query'] not in self.db('Local').index('Query',log=False):
                    self.warnLog('Local hits missing for BLAST. AlnNum < OneLine?!','missing_local',suppress=True)
                    return {}   #
                localdata = self.localData(hit['Hit'],search['Query'])
            ## ~ [1a] Check for existing data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gdb = self.db('GABLAM')
            gq = gdb.data(gdb.makeKey({'Query':search['Query'],'Hit':hit['Hit'],'QryHit':'Query'}))
            gh = gdb.data(gdb.makeKey({'Query':search['Query'],'Hit':hit['Hit'],'QryHit':'Hit'}))
            if gq and gh:
                #!# Check QAssemble
                return {'Query':gq,'Hit':gh}
            ## ~ [1b] Setup stats for calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gablam = {} # Global Alignment from BLAST Local Alignment Matrix
            for aln in ['Qry','QryO']: gablam[aln] = ['-'] * int(qrylen)            # O = ordered
            for aln in ['Hit','HitO']: gablam[aln] = ['-'] * int(hit['Length'])     # O = ordered
            for qh in ['Qry','Hit']:
                gablam['%sDirn' % qh] = 'None'
                gablam['%sStart' % qh] = len(gablam[qh])
                gablam['%sEnd' % qh] = 0
                gablam['%sDirnO' % qh] = 'None'
                gablam['%sStartO' % qh] = len(gablam[qh])
                gablam['%sEndO' % qh] = 0
            lpairs = []     # List of paired query-hit residue tuples for ordered GABLAM
            orientation = {'None':{True:'Bwd',False:'Fwd'},'Fwd':{True:'Both',False:'Fwd'},
                           'Bwd':{True:'Bwd',False:'Both'},'Both':{True:'Both',False:'Both'}}
            ## ~ [1c] QAssemble stats setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if query not in self.dict['QAssemble']:
                qassemble = ['-'] * int(qrylen)
                if self.getBool('QAssemble'): self.dict['QAssemble'][query] = qassemble
            else: qassemble = self.dict['QAssemble'][query]
            if not self.getBool('QAssemble'): qassemble = []
            elif not self.getBool('SelfSum') and qrylen == int(hit['Length']):  # Check for self matches
                if query == hit['Hit']: qassemble = []
                elif string.split(hit['Hit'],'__')[-1] == query: qassemble = []
                elif string.split(query,'__')[-1] == hit['Hit']: qassemble = []
            ### ~ [2] Compile alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.deBug('GABLAMO for Hit %s: %d alignments' % (self.info['Name'],self.alnNum()))
            lcutx = 0; lidx = 0; coredirn = None
            for aln in localdata:
                if aln['Length'] < self.int['LocalCut']: lcutx += 1; continue
                if 100.0 * aln['Identity'] / aln['Length'] < self.getNum('LocalIDMin'): lidx += 1; continue
                ## ~ [2a] Assess for backwards hit as in nucleotide BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qbackwards = aln['QryEnd'] < aln['QryStart']       # Whether match reversed (e.g. BLASTX)
                sbackwards = aln['SbjEnd'] < aln['SbjStart']       # Whether match reversed (e.g. TBLASTN)
                gablam['QryDirn'] = orientation[gablam['QryDirn']][qbackwards]
                gablam['HitDirn'] = orientation[gablam['HitDirn']][sbackwards]
                if not coredirn: coredirn = gablam['HitDirn']
                gablam['QryStart'] = min(gablam['QryStart'],aln['QryStart'],aln['QryEnd'])
                gablam['HitStart'] = min(gablam['HitStart'],aln['SbjStart'],aln['SbjEnd'])
                gablam['QryEnd'] = max(gablam['QryEnd'],aln['QryStart'],aln['QryEnd'])
                gablam['HitEnd'] = max(gablam['HitEnd'],aln['SbjStart'],aln['SbjEnd'])
                #x#    self.log.printLog('#REV','Query has reverse orientation for alignment %s: Ignoring for GABLAM.' % self.aln.index(aln))
                #x#    continue
                ## ~ [2b] Assess aln for GABLAMO ordering with previous alns ~~~~~~~~~~~~~~~~~~~~~~ ##
                order_ok = True #x# not (qbackwards or sbackwards)      # GABLAMO is in forward direction only
                #!# Update method to allow GABLAMO in all directions #!#
                for pair in lpairs:
                    if qbackwards == sbackwards:
                        if aln['QryStart'] < pair[0] and aln['SbjStart'] > pair[1]: order_ok = False
                        elif aln['QryStart'] > pair[0] and aln['SbjStart'] < pair[1]: order_ok = False
                        elif aln['QryEnd'] < pair[0] and aln['SbjEnd'] > pair[1]: order_ok = False
                        elif aln['QryEnd'] > pair[0] and aln['SbjEnd'] < pair[1]: order_ok = False
                    else:
                        if aln['QryStart'] < pair[0] and aln['SbjStart'] < pair[1]: order_ok = False
                        elif aln['QryStart'] > pair[0] and aln['SbjStart'] > pair[1]: order_ok = False
                        elif aln['QryEnd'] < pair[0] and aln['SbjEnd'] < pair[1]: order_ok = False
                        elif aln['QryEnd'] > pair[0] and aln['SbjEnd'] > pair[1]: order_ok = False
                if order_ok:
                    lpairs.append((aln['QryStart'],aln['SbjStart']))
                    lpairs.append((aln['QryEnd'],aln['SbjEnd']))
                    gablam['QryDirnO'] = orientation[gablam['QryDirnO']][qbackwards]
                    gablam['HitDirnO'] = orientation[gablam['HitDirnO']][sbackwards]
                    gablam['QryStartO'] = min(gablam['QryStartO'],aln['QryStart'],aln['QryEnd'])
                    gablam['HitStartO'] = min(gablam['HitStartO'],aln['SbjStart'],aln['SbjEnd'])
                    gablam['QryEndO'] = max(gablam['QryEndO'],aln['QryStart'],aln['QryEnd'])
                    gablam['HitEndO'] = max(gablam['HitEndO'],aln['SbjStart'],aln['SbjEnd'])
                ## ~ [2c] Perform GABLAM calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qres = aln['QryStart'] - 1     # Count from zero in list
                hres = aln['SbjStart'] - 1
                shop = qhop = 1
                if qbackwards: qhop = -1
                if sbackwards: shop = -1
                iloop = 1       # Number of loops needed for each aln position
                if self.str['Type'] in ['blastx','tblastn']: iloop = 3
                try:
                    for r in range(aln['Length']): # Sim
                        for i in range(iloop):     # Break to exit this loop, depending on nt/protein GABLAM
                            if aln['AlnSeq'][r] == '+':
                                if gablam['Qry'][qres] != '|': gablam['Qry'][qres] = '+'
                                if qassemble and qassemble[qres] != '|': qassemble[qres] = '+'
                                if gablam['Hit'][hres] != '|': gablam['Hit'][hres] = '+'
                                if gablam['QryO'][qres] != '|' and order_ok: gablam['QryO'][qres] = '+'
                                if gablam['HitO'][hres] != '|' and order_ok: gablam['HitO'][hres] = '+'
                            elif aln['AlnSeq'][r] == '|' or re.search('[A-Za-z*]',aln['AlnSeq'][r]):   # ID
                                gablam['Qry'][qres] = '|'
                                if qassemble: qassemble[qres] = '|'
                                gablam['Hit'][hres] = '|'
                                if order_ok:
                                    gablam['QryO'][qres] = '|'
                                    gablam['HitO'][hres] = '|'
                            else:
                                if gablam['Qry'][qres] == '-': gablam['Qry'][qres] = 'X'
                                if qassemble and qassemble[qres] == '-': qassemble[qres] = 'X'
                                if order_ok and gablam['QryO'][qres] == '-': gablam['QryO'][qres] = 'X'
                                if gablam['Hit'][hres] == '-': gablam['Hit'][hres] = 'X'
                                if order_ok and gablam['HitO'][hres] == '-': gablam['HitO'][hres] = 'X'
                            if i < 2 and self.str['Type'] == 'blastx' and aln['QrySeq'][r] not in [' ','-']: qres += qhop    # DNA triplet
                            if i < 2 and self.str['Type'] == 'tblastn' and aln['SbjSeq'][r] not in [' ','-']: hres += shop    # DNA triplet
                        if re.search('[A-Za-z*]',aln['QrySeq'][r]): qres += qhop
                        if re.search('[A-Za-z*]',aln['SbjSeq'][r]): hres += shop
                    if qbackwards: qres += 2
                    if sbackwards: hres += 2
                    if qres != aln['QryEnd']:
                        #!# Could be DNA sequence #!#
                        dres = ((qres - aln['QryStart']) * 3) + aln['QryStart'] + 2
                        if dres != aln['QryEnd']:
                            self.errorLog('Hit %s: Query end position should be %d but reached %d (or %d) in Aln process!' % (hit['Hit'],aln['QryEnd'],qres,dres),False,False)
                            print aln
                            print gablam
                            raw_input('Continue?')
                            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                                gablam[aln] = ['X'] * qrylen    #!# Hit!! #!#
                            break
                    if hres != aln['SbjEnd']:
                        ## Check backwards match ##
                        bwd_hres = hres - aln['SbjStart']
                        #!# Could be DNA sequence #!#  dres = ((qres - aln['QryStart']) * 3) + aln['QryStart'] + 2
                        if aln['SbjStart'] != bwd_hres:
                            self.errorLog('Hit %s: Subject end position should be %d but reached %d in Aln process (Bwd: %d)!' % (hit['Hit'],aln['SbjEnd'],hres,bwd_hres),False,False)
                            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                                gablam[aln] = ['X'] * qrylen    #!# Hit!! #!#
                            break
                        else: self.errorLog('#ERR','Subject %s appears to have reverse orientation. GABLAMO Results will be wrong!' % hit['Hit'],printerror=False)
                except: self.errorLog('GABLAM problems: GABLAM stats will be wrong for %s Hit %s' % (search['Query'],hit['Hit']))
                if not keepaln:
                    for field in ['QrySeq','SbjSeq','AlnSeq']: aln[field] = ''

            ### ~ [3] Make GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.verbose(2,4,'\n*** Hit: %s ***' % hit['Hit'],1)
            if coredirn != gablam['HitDirnO']: gablam['HitDirnO'] = '%s.%s' % (coredirn,gablam['HitDirnO'])
            for aln in ['Qry','Hit','QryO','HitO']:  # O = ordered
                gablam[aln] = string.join(gablam[aln],'')
                #self.verbose(2,4,'%s: %s' % (aln, gablam[aln]),1)
                #self.bugPrint('%s: %s' % (aln, gablam[aln]))
                #self.deBug('%s: Len = %d; Unaln = %d' % (aln,len(gablam[aln]),string.count(gablam[aln],'-')))
            gdict = { 'Query':{}, 'Hit':{} }
            gdict['Query']['GABLAM Len'] = qrylen - string.count(gablam['Qry'],'-')
            gdict['Query']['GABLAM ID'] = string.count(gablam['Qry'],'|')
            gdict['Query']['GABLAM Sim'] = gdict['Query']['GABLAM ID'] + string.count(gablam['Qry'],'+')
            gdict['Query']['GABLAMO Len'] = qrylen - string.count(gablam['QryO'],'-')
            gdict['Query']['GABLAMO ID'] = string.count(gablam['QryO'],'|')
            gdict['Query']['GABLAMO Sim'] = gdict['Query']['GABLAMO ID'] + string.count(gablam['QryO'],'+')
            gdict['Hit']['GABLAM Len'] = hit['Length'] - string.count(gablam['Hit'],'-')
            gdict['Hit']['GABLAM ID'] = string.count(gablam['Hit'],'|')
            gdict['Hit']['GABLAM Sim'] = gdict['Hit']['GABLAM ID'] + string.count(gablam['Hit'],'+')
            gdict['Hit']['GABLAMO Len'] = hit['Length'] - string.count(gablam['HitO'],'-')
            gdict['Hit']['GABLAMO ID'] = string.count(gablam['HitO'],'|')
            gdict['Hit']['GABLAMO Sim'] = gdict['Hit']['GABLAMO ID'] + string.count(gablam['HitO'],'+')
            for gscore in ['Start','End','Dirn']:
                gdict['Query']['GABLAM %s' % gscore] = gablam['Qry%s' % gscore]
                gdict['Query']['GABLAMO %s' % gscore] = gablam['Qry%sO' % gscore]
                gdict['Hit']['GABLAM %s' % gscore] = gablam['Hit%s' % gscore]
                gdict['Hit']['GABLAMO %s' % gscore] = gablam['Hit%sO' % gscore]
            for gab in ['','O']:
                gdict['Hit']['GABLAM%s Frag' % gab] = []
                if self.getInt('GablamFrag') < 1: break     # No GABLAM Fragging
                i = 0
                #self.deBug('HIT SEQUENCE: %s' % gablam['Hit%s' % gab])
                #self.bugPrint('GABLAM%s: %s' % (gab,min(max(-1,gablam['Hit%s' % gab].find('|',i)),max(-1,gablam['Hit'].find('+',i)))))
                #self.deBug('Hit%s: Len = %d; Unaln = %d' % (gab,len(gablam['Hit%s' % gab]),string.count(gablam['Hit%s' % gab],'-')))
                while i > -1 and i < len(gablam['Hit%s' % gab]):
                    if log: self.progLog('\r#FRAG','GABLAM%s Fragging: %.2f%%' % (gab,i*100.0/len(gablam['Hit%s' % gab])))
                    i1 = max(-1,gablam['Hit%s' % gab].find('|',i))
                    i2 = max(-1,gablam['Hit%s' % gab].find('+',i))
                    i = max(i1,i2)
                    if i < 0: break
                    elif min(i1,i2) > -1: i = min(i1,i2)
                    fragend = fragstart = i
                    gx = 0
                    #i# Loop until the number of consecutive '-' characters reaches the gablamfrag limit.
                    #i# gablamfrag=1 will split at any gap.
                    #i# NOTE: It is impossible using gablamfrag here to divide separate local hits that are directly adjacent.
                    while gx < self.getInt('GablamFrag') and i < (len(gablam['Hit%s' % gab]) - 1):
                        if log: self.progLog('\r#FRAG','GABLAM%s Fragging: %.2f%%' % (gab,i*100.0/len(gablam['Hit%s' % gab])))
                        i += 1
                        if gablam['Hit%s' % gab][i] == '-': gx += 1
                        else: gx = 0; fragend = i
                    gdict['Hit']['GABLAM%s Frag' % gab].append((fragstart,fragend))
                    if gablam['Hit%s' % gab][i] != '-': i += 1
                    #self.deBug('%s' % (gdict['Hit']['GABLAM%s Frag' % gab]))
                if log: self.progLog('\r#FRAG','GABLAM%s Fragging done!     ' % gab)
            if lcutx and log: self.printLog('#CUT','%d local alignments < length %d ignored for GABLAM (localcut=X)' % (lcutx,self.int['LocalCut']))
            if lidx and log: self.printLog('#CUT','%d local alignments < %s%% local identity ignored for GABLAM (localidcut=X)' % (lidx,rje.sf(self.getNum('LocalIDCut'))))
            return gdict
        except:
            self.errorLog('Major problem during BLASTHit.globalFromLocal vs Hit %s (qylen=%d)' % (hit['Hit'],qrylen))
            gdict = { 'Query':{}, 'Hit':{} }
            for stat in ['GABLAM Len','GABLAM ID','GABLAM Sim','GABLAMO Len','GABLAMO ID','GABLAMO Sim']:
                gdict['Query'][stat] = -1
                gdict['Hit'][stat] = -1
            #self.dict['GABLAM'] = gdict
            #self.opt['GABLAM'] = True #!# Really?
            return gdict
#########################################################################################################################
    def formatDB(self,fasfile=None,protein=True,force=True,log=True,checkage=None,details=False): ### makeBLASTDB   #V2.1
        '''
        BLAST formats database given.
        >> fasfile:str = Name of file to form database [None]
        >> protein:bool = whether protein [True]
        >> force:bool = whether to overwrite an existing formatted DB [True]
        >> log:bool = whether to output format status to log [True]
        >> checkage:bool = whether to check age of DB files versus fasfile. If None uses IgnoreData option [None]
        >> details:bool = whether to output details of DB formatting to log file [False] 
        '''
        ### ~ [1] Check and update DBase as necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if checkage == None: checkage = not self.getAttribute('bool','IgnoreDate',default=False)
        if log: log = self.log
        if fasfile: self.setStr({'DBase':fasfile})
        elif self.getStr('DBase') == 'None': self.printLog('#ERR','Cannot format database "None"!'); raise ValueError
        dbase = self.getStr('DBase')
        ### ~ [2] Call formatdb if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if force or not checkForDB(dbfile=dbase,checkage=checkage,log=log,protein=protein):    #!# Check this #!#
            formatDB(dbase,self.blastPath(),protein,log,oldblast=self.getBool('OldBLAST'),details=details)
        else: self.printLog('#DB ','%s already formatted. (Force = False).' % self.getStr('DBase'),log=log)
#########################################################################################################################
    ### <5> ### BLAST Data manipulation                                                                                 #
#########################################################################################################################
    def hitToSeq(self,seqlist,hitlist=[],query=None,filename=None,appendfile=False,asdict=True):   ### Saves hits from given searches to sequence object/file
        '''
        Saves hits from given searches to sequence object and, if given, a file.
        >> seqlist:rje_seq.SeqList Object *Necessary!* 
        >> hitlist:list of Hit names [all if none]
        >> query:str [None] = Name of query to retrieve hits for. (Over-ruled by hitlist.)
        >> filename:str = Name of fasta output file - no save if None [None]
        >> appendfile:bool = Whether to append file
        >> asdict:bool [True] = Whether to return hit sequences as a dictionary of {Hit:Sequence} or just [Sequence] list
        << returns dictionary of {Hit:Sequence}
        '''
        try:### ~ [1] ~ Map hits to seqlist sequences and return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqobj = '%s' % seqlist
            newseqlist = seqobj.startswith('<rje_seqlist')
            if newseqlist: seqdic = seqlist.seqNameDic()
            else: seqdic = seqlist.seqNameDic(proglog=False)
            hitseq = {}; hx = 0; hitseqlist = []
            if not hitlist: hitlist = self.queryHits(query)
            for hit in hitlist:
                hitseq[hit] = None
                if seqdic.has_key(hit):
                    hitseq[hit] = seqdic[hit]; hitseqlist.append(seqdic[hit])
                    hx += 1
                else:
                    newseq = seqlist.seqFromBlastDBCmd(hit,self.getStr('DBase'))
                    if newseq:
                        hitseq[hit] = newseq; hx += 1
                        self.verbose(2,4,'Seq %d retrieved from %s: %s.' % (seqlist.seqNum(),self.getStr('DBase'),newseq.shortName()),1)
                    else: self.errorLog('No Seq retrieval for %s using BLASTRun.hitToSeq().' % hit,printerror=False)
            if hx and hitseq and filename and filename.lower() != 'none':
                if newseqlist: seqlist.saveSeq(seqs=hitseqlist,seqfile=filename,append=appendfile)
                else: seqlist.saveFasta(seqs=seqlist.seq[-hx:],seqfile=filename,append=appendfile)
                self.dict['Output']['hitseq'] = filename    # Hit sequences in Fasta format
            if asdict: return hitseq
            else: return rje.dictValueList(hitseq,hitlist)
        except: self.errorLog('Major error during BLASTRun.hitToSeq().');  raise
#########################################################################################################################
    def localFragFas(self,byquery=True,combined=True,fragdir='BLASTFAS/',outbase=None,minfrag=0,addflanks=0,revcomp=True,append=False):
        '''
        This method executes the local fragment output, merging where required and adding flanking regions if desired.
        It uses the Hit and Local tables to identify the fragments, so these need to be complete.
        NOTE: Merging is decided based on pre-flanked sequences. (Set fragmerge to >= addflanks if flanks should merge.)
        @param byquery:bool [True] = Whether to output individual query FASTA files.
        @param combined:bool [True] = Whether to generate combined output to OUTBASE.fas.
        @param fragdir:str ['BLASTFAS/'] = Output directory for individual query fasta files.
        @param outbase:str [self.baseFile()] = Basefile for combined output file (combined=True).
        @param minfrag:int [0] = Minimum length of a local hit fragment to be output.
        @param addflanks:int [0] = Length of flanking sequence to be added to outputs.
        @param revcomp:bool [True] = Whether to reverse complement reverse strand hits (DNA only).
        @param append:bool [False] = Whether fasta files should be appended rather than overwritten.
        @return:
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not outbase: outbase = self.baseFile()               # Basename for combined fasta file
            #i# Can call with addflanks = self.getInt('AddFlanks') to use read parameter
            #i# gabfrag is the (direction-specific) within-query local hit merge parameter
            #i# Added flanks are now indistinguishable from the core hits when merging
            #X#gabfrag = max(self.getInt('GablamFrag') - addflanks,0)  # Want to join fragments within gabfrag distance
            gabfrag = self.getInt('GablamFrag')      # Want to join qry-hit fragments within gabfrag distance
            #i# fragmerge is the (direction-specific) (within- and) between-query local hit merge parameter
            #fragmerge = max(self.getInt('FragMerge') - addflanks,0)  # Want to join fragments within gabfrag distance
            fragmerge = self.getInt('FragMerge')  # Want to join fragments within gabfrag distance - includes flanks
            ## ~ [0a] Setup source database data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hdb = self.db('Hit')                                    # Need this for Hit Lengths
            hqfield = 'Query'
            ldb = self.db('Local')                                  # Use this for fragment data
            lqfield = 'Query'
            if lqfield not in ldb.fields() and 'Qry' in ldb.fields(): lqfield = 'Qry'
            if not hdb.entryNum(): return self.printLog('#FAS','No Hit entries for %s localFragFas output.' % outbase)
            if not ldb.entryNum(): return self.printLog('#FAS','No Local entries for %s localFragFas output.' % outbase)
            ## ~ [0b] Setup objects for storing fragment data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hitlen = {}                                             # Dictionary of Hit:Length
            hitfrag = {}                                            # Dictionary of Hit:{Dirn:[tuple list]}
            qryfrag = {}                                            # Dictionary of Qry:{Hit:{Dirn:[tuple list]}}
            for qry in hdb.indexKeys(hqfield): qryfrag[qry] = {}
            for entry in hdb.entries():
                hitlen[entry['Hit']] = entry['Length']
                hitfrag[entry['Hit']] = {'Fwd':[],'Bwd':[]}
                qryfrag[entry[hqfield]][entry['Hit']] = {'Fwd':[],'Bwd':[]}
            ### ~ [1] ~ Generate Fragment Tuples for hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Tuples are always converted into (start,end) relative to the hit sequence
            #i# When in the Bwd direction, they will need to be reverse complemented and positions reversed
            #i# This happens if revcomp=T but not if revcomp=F. If revcomp=F, the direction is always set to Fwd.
            ex = 0.0; etot = ldb.entryNum()
            for entry in ldb.entries():
                self.progLog('\r#FRAG','Generating hit fragments: %.1f%%' % (ex/etot)); ex += 100.0
                if entry['Length'] < minfrag: continue
                hit = entry['Hit']
                qry = entry[lqfield]
                if entry['SbjStart'] > entry['SbjEnd']:
                    dirn = 'Bwd'
                    frag = (max(entry['SbjEnd']-addflanks,1),min(entry['SbjStart']+addflanks,hitlen[hit]))
                    if not revcomp: dirn = 'Fwd'
                else:
                    dirn = 'Fwd'
                    frag = (min(entry['SbjStart']-addflanks,hitlen[hit]),max(entry['SbjEnd']+addflanks,1))
                qryfrag[qry][hit][dirn].append(frag)
            self.progLog('\r#FRAG','Generating hit fragments: compressing.')
            ### ~ [2] ~ Collapse Fragment Tuples for hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] First, collapse within-queries using gabfrag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qfx = 0
            for qry in qryfrag:
                for hit in qryfrag[qry]:
                    for dirn in qryfrag[qry][hit]:
                        qryfrag[qry][hit][dirn] = rje.collapseTupleList(qryfrag[qry][hit][dirn],joindistance=gabfrag)
                        qfx += len(qryfrag[qry][hit][dirn])
                        if combined: hitfrag[hit][dirn] += qryfrag[qry][hit][dirn]
            self.printLog('\r#FRAG','Generated %s collapsed hit fragments.' % rje.iStr(qfx))
            ## ~ [2b] Second, collapse between-queries using fragmerge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if combined:
                hfx = 0
                for hit in hitfrag:
                    for dirn in hitfrag[hit]:
                        #hitfrag[hit][dirn] = rje.collapseTupleList(hitfrag[hit][dirn],joindistance=gabfrag)
                        hitfrag[hit][dirn] = rje.collapseTupleList(hitfrag[hit][dirn],joindistance=fragmerge,overlaps=True)
                        hfx += len(hitfrag[hit][dirn])
                self.printLog('\r#FRAG','%s -> %s combined hit fragments.' % (rje.iStr(qfx),rje.iStr(hfx)))
            ### ~ [3] ~ Output of fragments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug(self.getStr('DBase'))
            #self.debug(rje.exists(self.getStr('DBase')))
            hitseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('DBase'),'autoload=T','autofilter=F','mode=file'])
            hitseqdict = hitseqlist.makeSeqNameDic('short')
            ## ~ [3a] Output individual query files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if byquery:
                fasx = 0; qtot = len(qryfrag)
                for qry in rje.sortKeys(qryfrag):
                    qfile = '%s%s.fas' % (fragdir,qry)
                    qseq = []
                    for hit in rje.sortKeys(qryfrag[qry]):
                        (seqname,sequence) = hitseqlist.getSeq(hitseqdict[hit])
                        seqlen = len(sequence)
                        for (fragstart,fragend) in qryfrag[qry][hit]['Fwd']:
                            fragstart = max(1,fragstart)
                            fragend = min(fragend,seqlen)
                            fraqseq = sequence[fragstart-1:fragend]
                            fragname = string.split(seqname)
                            fragname[0] += '.%s-%s' % (rje.preZero(fragstart,seqlen),rje.preZero(fragend,seqlen))
                            if len(fragname) == 1: fragname.append('No description')
                            fragname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                            #fragname.insert(1,'%s GABLAM Hit:' % qry)
                            fragname = string.join(fragname)
                            qseq.append((fragname,fraqseq))
                        for (fragstart,fragend) in qryfrag[qry][hit]['Bwd']:
                            fragstart = max(1,fragstart)
                            fragend = min(fragend,seqlen)
                            fraqseq = rje_sequence.reverseComplement(sequence[fragstart-1:fragend])
                            fragname = string.split(seqname)
                            fragname[0] += '.%s-%s' % (rje.preZero(fragend,seqlen),rje.preZero(fragstart,seqlen))
                            fragname.insert(1,'RevComp')
                            fragname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                            #fragname.insert(1,'%s GABLAM Hit:' % qry)
                            fragname = string.join(fragname)
                            qseq.append((fragname,fraqseq))
                    if qseq:
                        hitseqlist.saveSeq(qseq,qfile,append=append,seqtuples=True)
                        fasx += 1
                self.printLog('\r#FAS','Generated %s fasta files for %s queries.' % (rje.iStr(fasx),rje.iStr(qtot)))
            ## ~ [3b] Output combined merged hit file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if combined:
                hfile = '%s.fas' % (outbase)
                hseq = []
                for hit in rje.sortKeys(hitfrag):
                    (seqname,sequence) = hitseqlist.getSeq(hitseqdict[hit])
                    seqlen = len(sequence)
                    for (fragstart,fragend) in hitfrag[hit]['Fwd']:
                        fragstart = max(1,fragstart)
                        fragend = min(fragend,seqlen)
                        fraqseq = sequence[fragstart-1:fragend]
                        fragname = string.split(seqname)
                        fragname[0] += '.%s-%s' % (rje.preZero(fragstart,seqlen),rje.preZero(fragend,seqlen))
                        if len(fragname) == 1: fragname.append('No description')
                        fragname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                        fragname = string.join(fragname)
                        hseq.append((fragname,fraqseq))
                    for (fragstart,fragend) in hitfrag[hit]['Bwd']:
                        fragstart = max(1,fragstart)
                        fragend = min(fragend,seqlen)
                        fraqseq = rje_sequence.reverseComplement(sequence[fragstart-1:fragend])
                        fragname = string.split(seqname)
                        fragname[0] += '.%s-%s' % (rje.preZero(fragend,seqlen),rje.preZero(fragstart,seqlen))
                        fragname.insert(1,'RevComp')
                        fragname[-1] += '|(Pos:%s..%s)' % (rje.iStr(fragstart),rje.iStr(fragend))
                        fragname = string.join(fragname)
                        hseq.append((fragname,fraqseq))
                hitseqlist.saveSeq(hseq,hfile,append=append,seqtuples=True)
                self.dict['Output']['fragfas'] = hfile    # Hit fragments in Fasta format
        except: self.errorLog('Major error during BLASTRun.localFragFas().');  raise
#########################################################################################################################
    def blastClusters(self,seqfile,seqdict={},dna=False,keepblast=False):   ### Performs BLAST and returns sequence clusters #V2.3
        '''
        Performs BLAST and returns dismatrix. Can then use dismatrix.cluster() to return sequence clusters.
        >> seqfile:str [None] = File to use for clustering.
        >> dna:bool [False] = Whether sequence is DNA or not
        >> seqdict:dict = Optional dictionary linking sequence names to alternative keys for dismatrix object
        >> keepblast:bool [False] = Whether to keep BLAST file or delete.
        << returns a dismatrix object with diskey distances
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.blastIDMatrix(seqfile,seqdict,dna,keepblast).cluster(maxdis=1.0,singletons=True)
        except: self.errorLog('Problem during BLASTRun.resultsMatrix()'); raise
#########################################################################################################################
    def blastIDMatrix(self,seqfile,seqdict={},dna=False,keepblast=False): ### Performs BLAST and returns dismatrix. #V2.3
        '''
        Performs BLAST and returns dismatrix. Can then use dismatrix.cluster() to return sequence clusters.
        >> seqfile:str [None] = File to use for clustering.
        >> dna:bool [False] = Whether sequence is DNA or not
        >> seqdict:dict = Optional dictionary linking sequence names to alternative keys for dismatrix object
        >> keepblast:bool [False] = Whether to keep BLAST file or delete.
        << returns a dismatrix object with diskey distances
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setStr({'InFile':seqfile,'DBase':seqfile,'Name':'%s.self.blast' % rje.baseFile(seqfile),'Type':'blastp'})
            seqx = rje_seq.SeqCount(self,seqfile)
            self.setInt({'OneLine':seqx,'HitAln':seqx})
            if dna: self.setStr({'Type':'blastn'})
            ## ~ [1a] FormatDB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.formatDB(fasfile=seqfile,force=self.force(),protein=not dna)
            if not checkForDB(dbfile=seqfile,checkage=False,log=self.log,protein=not dna):
                self.errorLog('FormatDB failed for unknown reasons.',printerror=False)
                raise IOError
            ### ~ [2] Perform and read BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.blast(cleandb=True,use_existing=not self.force())
            if not self.readBLAST(gablam=True,unlink=not keepblast):
                self.errorLog('Major problem with BLAST for unknown reasons.')
                raise IOError
            ### ~ [3] Make Distance Matrix from BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.resultsIDMatrix(seqdict)
        except: self.errorLog('Problem during BLASTRun.blastIDMatrix()'); raise            
#########################################################################################################################
    def resultsIDMatrix(self,seqdict={}): ### Generates dismatrix from BLAST results.              #V2.3
        '''
        Performs BLAST and returns dismatrix. Can then use dismatrix.cluster() to return sequence clusters.
        >> seqdict:dict = Optional dictionary linking sequence names to alternative keys for dismatrix object
        << returns a dismatrix object with GABLAM ID distances
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dismatrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
            if not seqdict:
                for query in self.db('Search').dataKeys(): seqdict[query] = query
            ### ~ [2] Make Distance Matrix from BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for query in self.db('Search').index('Query'):
                seq = seqdict[query]
                dismatrix.addDis(seq,seq,0.0)      # Always self-distance of zero, even without BLAST hits
            for entry in self.db('GABLAM').entries():
                if entry['QryHit'] != 'Query': continue
                qlen = self.db('Search').indexDataList('Query',entry['Query'],'Length')[0]
                seq = seqdict[entry['Query']]
                hit = entry['Hit']
                if seqdict.has_key(hit):
                    if seqdict[hit] != seq:     # Has not hit itself
                        dismatrix.addDis(seq,seqdict[hit],1 - (entry['GABLAM ID'] / float(qlen)))
                else: self.errorLog('No sequence for hit "%s"!' % hit,printerror=False)
            return dismatrix
        except: self.errorLog('Problem during BLASTRun.resultsIDMatrix()'); raise
#########################################################################################################################
    def qAssembleFas(self):   ### Special runtime for running with outfmt 4 and converting to fasta
        '''Special runtime for running with outfmt 4 and converting to fasta.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Name'):
                if self.getStr('InFile') == self.getStr('DBase'): self.setStr({'Name':'%s.self.blast' % (rje.baseFile(self.getStr('InFile')))})
                else: self.setStr({'Name':'%s.%s.blast' % (rje.baseFile(self.getStr('InFile')),rje.baseFile(self.getStr('DBase'),strip_path=True))})
                self.printLog('#NAME',self.getStr('Name'))
            if self.baseFile().lower() in ['','none']: self.baseFile(rje.baseFile(self.getStr('Name')))
            bfile = self.getStr('Name')
            ## Check formatting and BLASTing - might already have results
            needtoblast = self.force() or not rje.checkForFile(bfile)
            if needtoblast and not self.getBool('FormatDB'):
                if self.getStr('Type') in ['blastn','tblastn','tblastx']: self.formatDB(protein=False,force=self.force())
                else: self.formatDB(force=self.force())
            ### ~ [2] ~ BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if needtoblast:
                if not rje.checkForFile(self.getStr('InFile')): raise IOError('Cannot find BLAST query file: %s' % self.getStr('InFile'))
                if not rje.checkForFile(self.getStr('DBase')): raise IOError('Cannot find BLAST database file: %s' % self.getStr('DBase'))
                if self.getStrLC('BLASTOpt'): self.setStr({'BLASTOpt':'%s -outfmt 4' % self.getStr('BLASTOpt')})
                else: self.setStr({'BLASTOpt':'-outfmt 4'})
                self.blast(use_existing=not self.force())
            if not rje.checkForFile(bfile): raise IOError('Cannot find BLAST results file: %s' % bfile)
            ### ~ [3] ~ Make Fasta ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.assemblyToAlignments(bfile,basefile=self.baseFile(),save=True,returndict=False,fullquery=self.getBool('QComplete'))
        except: self.errorLog('Problem during BLASTRun.qAssembleFas()');
#########################################################################################################################
    def assemblyToAlignments(self,blastres='',basefile='',save=True,returndict=False,fullquery=False): ### Converts BLAST outfmt 4 into fasta alignment files.
        '''
        Converts BLAST outfmt 4 into fasta alignment files.
        >> blastres:str [''] = Name of BLAST results file in assembly format. Will use self.getStr('Name') if blank.
        >> basefile:str [''] = Basefile for output fasta files. Will use blastres if flank. (Will keep path.)
        >> save:bool [True] = Whether to save alignments to fasta files.
        >> returndict: bool [False] = Whether to return a dictionary of {query shortname:[(name,sequence)]}
        >> fullquery: bool [False] = Whether to extend query sequence to full length.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not save and not returndict: raise ValueError('Need to save or return alignment dictionary!')
            #!# These parameters could be set for direct input with the addition of, if not X:
            if not blastres: blastres = self.getStr('Name')
            if not rje.exists(blastres): raise IOError('Blast results file "%s" missing!' % blastres)
            if not basefile: basefile = blastres
            qcmd = ['seqin=%s' % self.getStr('InFile'),'mode=file','autoload=T','autofilter=F']
            qseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+qcmd)
            qseqdict = qseqlist.seqNameDic()
            hcmd = ['seqin=%s' % self.getStr('DBase'),'mode=file','autoload=T','autofilter=F']
            hseqlist = rje_seqlist.SeqList(self.log,self.cmd_list+hcmd)
            hseqdict = hseqlist.seqNameDic()
            alndict = {}    # Dictionary to return if returndict=True
            resx = 0        # Number of queries with BLAST results
            outx = 0        # Number of alignments output
            rje.mkDir(self,self.getStr('QFasDir'))
            ### ~ [1] ~ Process File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            query = None        # Identity of query (or None of not reading)
            hits = []           # List of hits in order of significance
            readhits = False    # Whether expecting BLAST hit information
            ASS = open(blastres,'r')
            bline = ASS.readline()
            blastp = string.split(bline)[0]
            if 'BLAST' not in blastp.upper(): raise ValueError('%s does not appear to be a BLAST results file!' % blastres)
            self.printLog('#ALN','Generating alignments for %s assembly %s' % (blastp,blastres))
            while 1:
                bline = ASS.readline()
                if not bline: break
                bdata = string.split(bline)
                if bline.startswith('Query='):                                      # New query
                    query = bdata[1]; alnseq = {'Query':''}; hits = []; resx +=1
                    if query not in qseqdict: raise ValueError('Query %s not found in query sequence file!' % query)
                    (qname,qseq) = qseqlist.getSeq(qseqdict[query],'tuple')
                elif bline.startswith('Sequences producing significant alignments'):# One-line hit info
                    while query:
                        bline = rje.chomp(ASS.readline())
                        bdata = string.split(bline)
                        if readhits and not bline: readhits = False; break
                        if bline:
                            hit = bdata[0]
                            if hit.startswith('lcl|'): hit = hit[4:]
                            alnseq[hit] = []    # Each hit might have multiple (overlapping!) local alignments
                            hits.append(hit)
                            readhits = True
                    self.debug(alnseq)
                elif bline.startswith('Query') and query:                           # Read alignment
                    self.debug(bdata)
                    alnstart = bline.find(bdata[2])
                    alnend = alnstart + len(bdata[2])
                    alnpos = {'Query':[int(bdata[1]),int(bdata[3])]}
                    alnseq['Query'] += bdata[2]
                    asshits = []    # List of hits being read in... [startpos,endpos,sequence]
                    assx = 0
                    while not bline.startswith('Lambda'):
                        bline = rje.chomp(ASS.readline())
                        if not bline or bline.startswith('Lambda'): continue
                        bdata = string.split(bline)
                        self.bugPrint(bdata)
                        self.bugPrint(alnseq)
                        self.debug(alnpos)
                        if bline.startswith('Query'):
                            alnpos['Query'][1] = int(bdata[3])
                            alnseq['Query'] += bdata[2]
                            alnstart = bline.find(bdata[2])
                            alnend = alnstart + len(bdata[2])
                            assx = 0    # Index for next expected asshit alignment
                            self.bugPrint('%s - %s' % (alnstart,alnend))
                        else:   # Hit alnseq data will be a list of [hit,startpos,endpos,sequence] lists
                            hit = bdata[0]
                            if len(bdata) == 2:     # No sequence data, just gaps
                                asshits[assx][3] += string.replace(bline[alnstart:alnend],' ','-')
                                assx += 1
                                continue
                            self.bugPrint(asshits)
                            if assx >= len(asshits) or hit != asshits[assx][0] or rje.modulus(int(bdata[1])-asshits[assx][2]) != 1: # New alignment
                                asshits.insert(assx,[hit,int(bdata[1]),int(bdata[3]),string.replace(bline[alnstart:alnend],' ','-')])
                                if len(asshits[assx][3]) < len(alnseq['Query']):
                                    asshits[assx][3] = '-' * (len(alnseq['Query']) - len(asshits[assx][3])) + asshits[assx][3]
                            else:   # Continue
                                asshits[assx][2] = int(bdata[3])
                                asshits[assx][3] += string.replace(bline[alnstart:alnend],' ','-')
                            if not string.split(bline[alnend-1:alnend]):    # End of sequence
                                alnseq[hit].append(asshits.pop(assx))
                            else: assx += 1
                            self.deBug(asshits)
                            continue
                    ## ~ Tidy hits and generate alnpos dictionary entries ~~ ##
                    while asshits: alnseq[asshits[0][0]].append(asshits.pop(0))

                    ## ~ Adjust length ~~~~~~~~~~~~~~~~~~~~~~ ##
                    if fullquery:
                        self.debug(alndict)
                        self.debug(alnpos)
                        alndict[query] = [(qname,qseq[:alnpos['Query'][0]-1]+alnseq['Query']+qseq[alnpos['Query'][1]:])]
                        alnext = ('-' * (alnpos['Query'][0]-1), '-' * (len(qseq) - alnpos['Query'][1]))
                    else:
                        qname = '%s (Region %s..%s) %s' % (query,rje.iStr(alnpos['Query'][0]),rje.iStr(alnpos['Query'][1]),string.join(string.split(qname)[1:]))
                        alndict[query] = [(qname,alnseq['Query'])]
                    ## ~ Make sequence tuples ~~~~~~~~~~~~~~~ ##
                    qacc = query
                    if rje.matchExp('^\S+_\S+__(\S+)$',query): qacc = rje.matchExp('^\S+_\S+__(\S+)$',query)[0]
                    fconsensus = []
                    for hit in hits:
                        hitname = hseqlist.getSeq(hseqdict[hit],'tuple')[0]
                        self.debug(alnseq[hit])
                        if len(alnseq[hit]) == 1 and self.getStrLC('QConsensus') != 'full':   # Old code
                            alnpos[hit] = alnseq[hit][0][1:3]
                            alnseq[hit] = alnseq[hit][0][3]
                            hitname = '%s.%s (QAssemble Region %s..%s) %s' % (hit,qacc,rje.iStr(alnpos[hit][0]),rje.iStr(alnpos[hit][1]),string.join(string.split(hitname)[1:]))
                            alnseq[hit] += '-' * (len(alnseq['Query']) - len(alnseq[hit]))
                            if fullquery: alndict[query].append((hitname,alnext[0]+alnseq[hit]+alnext[1]))
                            else: alndict[query].append((hitname,alnseq[hit]))
                        else:
                            hx = 0
                            hconsensus = []; hconsname = '%s.%s %d hits (QAssemble Regions )' % (hit,qacc,len(alnseq[hit]))
                            for [hcheck,spos,epos,hseq] in alnseq[hit]:
                                if hcheck != hit: raise ValueError
                                hx += 1
                                newname = '%s.%s.%d (QAssemble Region %s..%s) %s' % (hit,qacc,hx,rje.iStr(spos),rje.iStr(epos),string.join(string.split(hitname)[1:]))
                                hseq += '-' * (len(alnseq['Query']) - len(hseq))
                                # Add flanks if fullquery=T
                                if fullquery: hseq = alnext[0]+hseq+alnext[1]
                                if self.getStrLC('QConsensus') == 'hit':
                                    hconsname = '%s %s..%s|' % (hconsname[:-1], rje.iStr(spos), rje.iStr(epos))
                                    hconsensus.append(hseq)
                                elif self.getStrLC('QConsensus') == 'full': fconsensus.append(hseq)
                                else: alndict[query].append((newname,hseq))
                            hconsname = '%s) %s' % (hconsname[:-1], string.join(string.split(hitname)[1:]))
                            if hconsensus:
                                hseq = self.makeConsensus(alndict[query][0][1],hconsensus)
                                alndict[query].append((hconsname, hseq))

                        #>>> OLD >>>#
                        #hitname = '%s (Region %s..%s) %s' % (hit,rje.iStr(alnpos[hit][0]),rje.iStr(alnpos[hit][1]),string.join(string.split(hitname)[1:]))
                        #alnseq[hit] += '-' * (len(alnseq['Query']) - len(alnseq[hit]))
                        #if fullquery: alndict[query].append((hitname,alnext[0]+alnseq[hit]+alnext[1]))
                        #else: alndict[query].append((hitname,alnseq[hit]))
                    if self.getStrLC('QConsensus') == 'full':
                        hseq = self.makeConsensus(alndict[query][0][1], fconsensus)
                        hconsname = '%s.%s (%s)' % (qname,qacc,string.join(hits,'; '))
                        alndict[query].append((hconsname, hseq))
                elif bline.startswith('***** No hits found *****'): query = None    # No hits
                elif bline.startswith('Effective search space used') and query:     # Process
                    fasfile = '%s%s.%s.qaln.fas' % (self.getStr('QFasDir'),os.path.basename(basefile),query)
                    if self.getStrLC('QConsensus') == 'full': fasfile = fasfile[:-3] + 'fullcons.fas'
                    elif self.getStrLC('QConsensus') == 'hit': fasfile = fasfile[:-3] + 'hitcons.fas'
                    OUTFAS = open(fasfile,'w')
                    for (name,sequence) in alndict[query]: OUTFAS.write('>%s\n%s\n' % (name,sequence))
                    OUTFAS.close(); outx += 1
                    self.printLog('#OUT','%s sequence alignment out to %s' % (rje.iLen(alndict[query]),fasfile))
                    if outx < 2: self.dict['Output']['qassembly'] = fasfile
                    accnum = string.split(query,'__',maxsplit=1)[-1]
                    self.dict['Output'][accnum] = fasfile
                    if not returndict: alndict = {}
            ASS.close()
            self.printLog('#ALN','%s alignments output from %s assemblies of %s queries.' % (rje.iStr(outx),rje.iStr(resx),rje.iLen(qseqdict)))
            if not outx:
                self.warnLog('No alignments output. If BLAST hits are expected, check that correct BLAST program used.')
            if returndict: return alndict
            else: return outx
        except: self.errorLog('Problem during BLASTRun.assemblyToAlignments()'); raise
#########################################################################################################################
    def makeConsensus(self, qryseq, hitseq, mindepth=1):  ### Makes a consensus sequence from reference and subsequences
        '''
        Makes a consensus sequence from reference and subsequences.
        >> qryseq:str = Full length Reference sequence
        >> hitseq:list of strings making consensus (may be missing Cterm/3')
        '''
        try:  ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#CONS','Generating consensus from qryseq and %d hit sequeces' % len(hitseq))
            #x#if len(hitseq) == 1: return qryseq
            conseq = ''
            wild = {True: 'N', False: 'X'}[self.nt()]
            nocons = ['*', wild, '']
            if mindepth < 2: nocons.append('-')
            if mindepth > len(hitseq):
                self.warnLog('Reduced consensus mindepth=%d to sequence number (mindepth=%d)' % (mindepth, len(hitseq)))
                mindepth = len(hitseq)
            ### ~ [1] ~ Compile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(qryseq)):
                r = qryseq[i]
                ilist = []
                flist = []   # Full list including the nocons elements
                nongapx = 0
                for seq in hitseq:
                    flist.append(seq[i])
                    if seq[i:i + 1] not in nocons: ilist.append(seq[i])
                    if seq[i] != '-': nongapx += 1
                if ilist and r not in ilist: r = ilist[0]
                elif not ilist:
                    if wild in flist: r = wild
                    elif '*' in flist: r = '*'
                    else: r = '-'
                for x in ilist:
                    if x == '-' and nongapx >= mindepth: continue
                    if ilist.count(x) > ilist.count(r):
                        r = x
                    elif nongapx >= mindepth and r == '-':
                        r = x
                conseq = conseq + r
            return conseq
        except:
            self.errorLog(rje_zen.Zen().wisdom()); return qryseq
#########################################################################################################################
    def reduceLocal(self,locdb=None,queries=[],hits=[],sortfield='Identity',keepself=False,minloclen=0,minlocid=0):    ### Reduces local BLAST alignments to cover each hit region only once
        '''
        Reduces local BLAST alignments to cover each hit region only once. Local alignments are sorted by identity
        (unless sortfield=X changed) and processed in order. Other local alignments that overlap are truncated and
        updated. Any alignments completely overlapped are removed. Processes and returns a COPY of the table.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param sortfield:str ['Identity'] = LocalDB field used to sort local alignments.
        @param keepself:bool [False] = Whether to include self query-hit pairs in assessment.
        @param minloclen:int [0] = Minimum local length to keep.
        @param minlocid:pc [0] = Minimum local %identity (0-100) to keep.
        @return: copy of local table, filtered and reduced.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ BLAST Local Alignment Reduction ~~~~~ ##')
            if not locdb: locdb = self.db('Local')
            if sortfield not in locdb.fields():
                raise ValueError('"%s" is not a valid Local hit table field (%s)' % (sortfield,string.join(locdb.fields(),'|')))
            # Create a copy to protect initial data
            bdb = self.db().copyTable(locdb,'%s.reduced' % locdb.name())
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            # Filter if required
            if queries: bdb.dropEntriesDirect('Query',queries,inverse=True)
            if hits: bdb.dropEntriesDirect('Hit',hits,inverse=True)
            if not keepself: bdb.dropEntries(['Query=Hit'])
            bdb.dataFormat({'AlnID':'int','BitScore':'num','Expect':'num','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'})
            btot = bdb.entryNum()
            bdb.dropEntries(['Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing short local hits')

            if minlocid > 0.0:
                badidx = 0
                for entry in bdb.entries()[0:]:
                    if 100.0 * float(entry['Identity']) / int(entry['Length']) < minlocid:
                        badidx += 1
                        bdb.dropEntry(entry)
                self.printLog('#MINID','Dropped %s entries < LocalIDMin=%s%%' % (rje.iStr(badidx),rje.sf(minlocid)))

            mx = btot - bdb.entryNum()
            ### ~ [1] Cycle through local alignments, reducing as required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bentries = bdb.sortedEntries(sortfield,reverse=sortfield in ['Identity','Positives','Length','BitScore'])   # List of all entries (sorted) to process
            alignpos = {}; ax = 0   # Dictionary of {Hit:[(start,stop) list of positions included in local aln]}
            while bentries:
                ## ~ [1a] Grab next best remaining hit from bentries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = bentries.pop(0)     # This is best remaining hit
                #self.debug('QrySeq: %d; AlnSeq: %d; SbjSeq: %d' % (len(entry['QrySeq']),len(entry['AlnSeq']),len(entry['SbjSeq'])))
                ax += 1
                self.progLog('\r#LOCALN','Processing local alignments: %s -> %s' % (rje.iLen(bentries),rje.iStr(ax)))
                hit = entry['Hit']
                region = (min(entry['SbjStart'],entry['SbjEnd']),max(entry['SbjStart'],entry['SbjEnd']))
                ## ~ [1b] Update alignpos dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if hit not in alignpos: alignpos[hit] = []
                alignpos[hit].append(region)
                alignpos[hit] = rje.collapseTupleList(alignpos[hit])
                ## ~ [1c] Adjust/Filter remaining entries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ex = 0
                while ex < len(bentries):
                    xentry = bentries[ex]
                    yentry = None   # Will be created if splitting
                    # Skip if different Hit
                    if entry['Hit'] != xentry['Hit']: ex += 1; continue
                    # Check for overlapping regions and remove
                    xregion = (min(xentry['SbjStart'],xentry['SbjEnd']),max(xentry['SbjStart'],xentry['SbjEnd']))
                    # No overlap
                    if xregion[1] < region[0] or xregion[0] > region[1]: ex += 1; continue
                    # Completely overlapped: remove
                    elif xregion[0] >= region[0] and xregion[1] <= region[1]:
                        bdb.dropEntry(xentry)
                        bentries.pop(ex)
                        continue
                    # Middle covered: split
                    elif region[0] > xregion[0] and region[1] < xregion[1]:
                        self.bugPrint('\nEntry splitting: %s vs %s' % (xregion,region))
                        xalnx = max(bdb.indexDataList('Hit',hit,'AlnID'))
                        yentry = rje.combineDict({'AlnID':xalnx+1},xentry,overwrite=False)
                        self.printLog('#ALNID','%s vs %s Aln %d -> %d & %d' % (xentry['Query'],xentry['Hit'],xentry['AlnID'],xentry['AlnID'],yentry['AlnID']))
                        #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},xentry,overwrite=False,replaceblanks=False))
                        self.trimLocal(xentry,trimend='End',trimto=region[0],sortends=True)   # Trim the end back to region[0]
                        #self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},xentry,overwrite=False,replaceblanks=False))
                        #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},yentry,overwrite=False,replaceblanks=False))
                        self.trimLocal(yentry,trimend='Start',trimto=region[1],sortends=True)   # Trim the start back to region[1]
                        #self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},yentry,overwrite=False,replaceblanks=False))
                    # Overlap at one end
                    elif region[0] <= xregion[1] <= region[1]:  # End overlaps with focal entry
                        #self.bugPrint('\nEnd overlap: %s vs %s' % (xregion,region))
                        self.trimLocal(xentry,trimend='End',trimto=region[0],sortends=True)   # Trim the end back to region[0]
                    elif region[0] <= xregion[0] <= region[1]:  # Start overlaps with focal entry
                        #self.bugPrint('\nStart overlap: %s vs %s' % (xregion,region))
                        self.trimLocal(xentry,trimend='Start',trimto=region[1],sortends=True)   # Trim the start back to region[1]
                    else: raise ValueError('Entry filtering has gone wrong: %s vs %s' % (xregion,region))
                    ## Check lengths
                    if xentry['Length'] >= minloclen: ex += 1
                    else: bdb.dropEntry(xentry); bentries.pop(ex); mx += 1
                    if yentry and yentry['Length'] >= minloclen:
                        bdb.addEntry(yentry)
                        bentries.insert(ex,yentry)
                        ex += 1
                    elif yentry: mx += 1
            ### ~ [2] Check and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#LOCALN','Processing local alignments: %s -> %s (%s failed to meet minloclen=%d)' % (rje.iStr(btot),rje.iStr(ax),rje.iStr(mx),minloclen))
            if ax != bdb.entryNum(): raise ValueError('EntryNum mismatch following reduceLocal(): %s best entries but %s local alignments' % (ax,bdb.entryNum()))
            return bdb
        except: self.errorLog('Problem during BLASTRun.reduceLocal()'); return None
#########################################################################################################################
    def trimLocal(self,lentry,trimend,trimto,sortends=True,debug=False):  # Trims local alignment entry data to hit coordinates
        '''
        Trims local alignment entry data.
        @param lentry: local alignment entry
        @param trimend: which end to trim (Start/End)
        @param trimto: position to trim up to (and including)
        @param sortends: whether to switch Start/End if reversed match is
        @return: modified lentry. (NOTE: Modified in place.)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            for field in ['QrySeq','SbjSeq','AlnSeq']:
                if not lentry[field]:
                    self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
                    raise ValueError('Empty %s!' % field)
            fwd = lentry['SbjStart'] <= lentry['SbjEnd']
            if sortends and not fwd: trimend = {'Start':'End','End':'Start'}[trimend]
            if debug: self.bugPrint('Fwd:%s => %s' % (fwd,trimend))
            # Check need to trim
            if fwd and not lentry['SbjStart'] <= trimto <= lentry['SbjEnd']:
                self.warnLog('Fwd trimLocal() called but Start <= trimpos <= End not met.','trimlocal',suppress=True)
                return lentry
            if not fwd and not lentry['SbjStart'] >= trimto >= lentry['SbjEnd']:
                self.warnLog('Bwd trimLocal() called but End <= trimpos <= Start not met.','trimlocal',suppress=True)
                return lentry
            # Sort out starting positions in Query/Subject
            #qpos = lentry['QryStart'] #lentry['Qry%s' % trimend]
            spos = lentry['SbjStart']   #lentry['Sbj%s' % trimend]
            # Sort out starting position in alignment
            ai = 0
            if fwd: sdir = 1
            else: sdir = -1
            if trimend == 'Start': trimto += sdir  # For Start trim, want to go beyond trimto for alignment split
            ### ~ [1] Trim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Locate trimto ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            self.bugPrint('Fwd:%s Scan %d pos until %d reaches %s' % (fwd,len(lentry['SbjSeq']),spos,trimto))
            while (fwd and spos < trimto) or (not fwd and spos > trimto):
                ai += 1
                if lentry['SbjSeq'][ai] != '-': spos += sdir
            self.bugPrint('Reached AlnPos %d = SbjPos %d => Trim %s' % (ai,spos,trimend))
            ## ~ [1b] Split sequence and update stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if trimend == 'Start':  # Removing Start
                lentry['QryStart'] += len(lentry['QrySeq'][:ai]) - lentry['QrySeq'][:ai].count('-')
                lentry['QrySeq'] = lentry['QrySeq'][ai:]
                lentry['SbjStart'] += sdir * (len(lentry['SbjSeq'][:ai]) - lentry['SbjSeq'][:ai].count('-'))
                lentry['SbjSeq'] = lentry['SbjSeq'][ai:]
                if lentry['Positives']:
                    lentry['Positives'] -= lentry['AlnSeq'][:ai].count('+')
                    lentry['Positives'] -= lentry['AlnSeq'][:ai].count('|')
                lentry['Identity'] -= lentry['AlnSeq'][:ai].count('|')
                lentry['Length'] -= len(lentry['AlnSeq'][:ai])
                lentry['AlnSeq'] = lentry['AlnSeq'][ai:]
            else:   # Removing end
                lentry['QryEnd'] -= (len(lentry['QrySeq'][ai:]) - lentry['QrySeq'][ai:].count('-'))
                lentry['QrySeq'] = lentry['QrySeq'][:ai]
                lentry['SbjEnd'] -= sdir * (len(lentry['SbjSeq'][ai:]) - lentry['SbjSeq'][ai:].count('-'))
                lentry['SbjSeq'] = lentry['SbjSeq'][:ai]
                if lentry['Positives']:
                    lentry['Positives'] -= lentry['AlnSeq'][ai:].count('+')
                    lentry['Positives'] -= lentry['AlnSeq'][ai:].count('|')
                lentry['Identity'] -= lentry['AlnSeq'][ai:].count('|')
                lentry['Length'] -= len(lentry['AlnSeq'][ai:])
                lentry['AlnSeq'] = lentry['AlnSeq'][:ai]
            ### ~ [2] Return trimmed entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.bugPrint(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            return lentry
        except:
            self.debug(rje.combineDict({'QrySeq':'','SbjSeq':'','AlnSeq':''},lentry,overwrite=False,replaceblanks=False))
            self.errorLog('BLASTRun.trimLocal error'); raise
#########################################################################################################################
    def snpTableFromLocal(self,locdb=None,queries=[],hits=[],save=True):    ### Outputs a SNP table (similar format to nucmer) from local alignments
        '''
        Outputs a SNP table (similar format to nucmer) from local alignments.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param save:bool = Whether to save SNP table once generated.
        @return: SNP table
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            snpdb = db.addEmptyTable('SNP',['Locus','Pos','REF','ALT','AltLocus','AltPos'],keys=['Locus','Pos','AltLocus','AltPos'])
            snpdb.setBool({'TupleKeys':True})
            if not locdb: locdb = self.db('Local')
            if rje.listDifference(queries,locdb.index('Query')):
                self.warnLog('Queries for snpTableFromLocal() not found in local table','missing_queries',suppress=True)
            if rje.listDifference(hits,locdb.index('Hit')):
                self.warnLog('Hits for snpTableFromLocal() not found in local table','missing_hits',suppress=True)
            revcomp = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
            ### ~ [1] Generate SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for lentry in locdb.entries():
                if queries and lentry['Query'] not in queries: continue
                if hits and lentry['Hit'] not in hits: continue
                fwd = True
                sdir = 1
                if lentry['SbjStart'] > lentry['SbjEnd']: sdir = -1; fwd = False
                qpos = lentry['QryStart'] - 1
                spos = lentry['SbjStart'] - sdir
                for ai in range(len(lentry['AlnSeq'])):
                    if lentry['QrySeq'][ai] != '-': qpos += 1
                    if lentry['SbjSeq'][ai] != '-': spos += sdir
                    if lentry['QrySeq'][ai].upper() == lentry['SbjSeq'][ai].upper(): continue   # Identity
                    sentry = {'Locus':lentry['Query'],'Pos':qpos,'REF':lentry['QrySeq'][ai].upper(),
                              'ALT':lentry['SbjSeq'][ai].upper(),'AltLocus':lentry['Hit'],'AltPos':spos}
                    #!# Why reverse complement this?! ALready done in alignment sequences!!
                    #???# if not fwd: sentry['ALT'] = revcomp[sentry['ALT']]
                    snpdb.addEntry(sentry)
            ### ~ [2] Output SNP Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if save: snpdb.saveToFile()
            return snpdb
        except: self.errorLog('BLASTRun.snpTableFromLocal() error'); return None
#########################################################################################################################
    def noCoverage(self,locdb=None,queries=[],hits=[],save=True):    ### Outputs a table of regions not covered by local alignments
        '''
        Outputs a table of regions not covered by local alignments.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param save:bool = Whether to save SNP table once generated.
        @return: local region table
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not locdb: locdb = self.db('local')
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            locdb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            ## ~ [0a] CovDict and Lengths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            covdict = {'Qry':{},'Hit':{}}       # Tuples of coverage (to invert)
            lengths = {'Qry':{},'Hit':{}}
            for qry in locdb.index('Query'): covdict['Qry'][qry] = []; lengths['Qry'][qry] = None
            for hit in locdb.index('Hit'): covdict['Hit'][hit] = []; lengths['Hit'][hit] = None
            qdb = self.db('Search')
            if qdb:
                for entry in qdb.entries(): lengths['Qry'][entry['Query']] = int(entry['Length'])
            hdb = self.db('Hit')
            if hdb:
                for entry in hdb.entries(): lengths['Hit'][entry['Hit']] = int(entry['Length'])
            ### ~ [1] Make coverage lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in locdb.entries():
                covdict['Qry'][entry['Query']].append((min(entry['QryStart'],entry['QryEnd']),max(entry['QryStart'],entry['QryEnd'])))
                covdict['Hit'][entry['Hit']].append((min(entry['SbjStart'],entry['SbjEnd']),max(entry['SbjStart'],entry['SbjEnd'])))
            ### ~ [2] Make inverted no-coverage lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            nocovdb = db.addEmptyTable('nocoverage',['QryHit','Locus','Start','End'],keys=['QryHit','Locus','Start'])
            for Qry in covdict['Qry']:
                for (start,end) in rje.invertTupleList(covdict['Qry'][Qry],minx=1,maxx=lengths['Qry'][Qry]):
                    nocovdb.addEntry({'QryHit':'Qry','Locus':Qry,'Start':start,'End':end})
            for Hit in covdict['Hit']:
                for (start,end) in rje.invertTupleList(covdict['Hit'][Hit],minx=1,maxx=lengths['Hit'][Hit]):
                    nocovdb.addEntry({'QryHit':'Hit','Locus':Hit,'Start':start,'End':end})
            if save: nocovdb.saveToTable()
            return nocovdb
        except: self.errorLog('BLASTRun.noCoverage() error'); return None
#########################################################################################################################
    ### <6> ### Local Hit Reformatting                                                                                  #
#########################################################################################################################
    def saveGFF(self,filename=None,locdb=None,ftype='gff3',queries=[],hits=[],reftype='Hit',attributes=['*'],append=False):    ### Outputs local table to GFF format
        '''
        Outputs local table to GFF format.
        @param filename:str [self.basefile().locdb.ftype] = Output filename.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param ftype:str ['gff3'] = File type
        @param queries:list = Restrict analysis to search queries.
        @param hits:list = Restrict analysis to search hits.
        @param reftype:str = The field defining the genome that would be loaded in IGV.
        @param attributes:list ['*'] = List of table fields to populate attributes field, '*' = all "extra" fields.
        @param append:bool [False] = Whether to append rather than write to a new file.
        @return: True/False

        GFF3 format is defined as follows:

        Column 1: "seqid"
        The ID of the landmark used to establish the coordinate system for the current feature. IDs may contain any
        characters, but must escape any characters not in the set [a-zA-Z0-9.:^*$@!+_?-|]. In particular, IDs may not
        contain unescaped whitespace and must not begin with an unescaped ">".

        Column 2: "source"
        The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this
        feature. Typically this is the name of a piece of software, such as "Genescan" or a database name, such as
        "Genbank." In effect, the source is used to extend the feature ontology by adding a qualifier to the type
        creating a new composite type that is a subclass of the type in the type column.

        Column 3: "type"
        The type of the feature (previously called the "method"). This is constrained to be either a term from the
        Sequence Ontology or an SO accession number. The latter alternative is distinguished using the syntax SO:000000.
        In either case, it must be sequence_feature (SO:0000110) or an is_a child of it.
        - This will be set as "region".

        Columns 4 & 5: "start" and "end"
        The start and end coordinates of the feature are given in positive 1-based integer coordinates, relative to the
        landmark given in column one. Start is always less than or equal to end. For features that cross the origin of a
        circular feature (e.g. most bacterial genomes, plasmids, and some viral genomes), the requirement for start to be
        less than or equal to end is satisfied by making end = the position of the end + the length of the landmark
        feature.

        For zero-length features, such as insertion sites, start equals end and the implied site is to the right of the
        indicated base in the direction of the landmark.

        Column 6: "score"
        The score of the feature, a floating point number. As in earlier versions of the format, the semantics of the
        score are ill-defined. It is strongly recommended that E-values be used for sequence similarity features, and
        that P-values be used for ab initio gene prediction features.

        Column 7: "strand"
        The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for
        features that are not stranded. In addition, ? can be used for features whose strandedness is relevant, but
        unknown.

        Column 8: "phase"
        For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The
        phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning
        of this feature to reach the first base of the next codon. In other words, a phase of "0" indicates that the next
        codon begins at the first base of the region described by the current line, a phase of "1" indicates that the
        next codon begins at the second base of this region, and a phase of "2" indicates that the codon begins at the
        third base of this region. This is NOT to be confused with the frame, which is simply start modulo 3.

        For forward strand features, phase is counted from the start field. For reverse strand features, phase is counted
        from the end field.

        The phase is REQUIRED for all CDS features.

        Column 9: "attributes"
        A list of feature attributes in the format tag=value. Multiple tag=value pairs are separated by semicolons. URL
        escaping rules are used for tags or values containing the following characters: ",=;". Spaces are allowed in this
        field, but tabs must be replaced with the %09 URL escape. Attribute values do not need to be and should not be
        quoted. The quotes should be included as part of the value by parsers and not stripped.

        NOTE: Can add col=X for colouring features in IGV.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            raise ValueError('GFF output currently not available!')
            #if reftype.lower() == 'both':
            db = self.db()
            if not locdb: locdb = self.db('Local')
            locdb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
        except: self.errorLog('%s.saveGFF() error' % self.prog()); return None
#########################################################################################################################
    def saveSAM(self,filename=None,locdb=None,queries=[],hits=[],reftype=None,cigar=True,sequence=True,append=False):    ### Outputs local table to GFF format
        '''
        Outputs local table to GFF format.
        @param filename:str [self.basefile().locdb.sam] = Output filename.
        @param locdb:Table [self.db('local')] = Local hits database Table to modify.
        @param queries:list [] = Restrict analysis to search queries.
        @param hits:list [] = Restrict analysis to search hits.
        @param reftype:str ['Hit'] = The field defining the genome that would be loaded in IGV.
        @param cigar:bool [True] = Whether to output the CIGAR string if Aln fields are present.
        @param sequence:bool [True] = Whether to output the sequence if Aln fields are present.
        @param append:bool [False] = Whether to append rather than write to a new file.
        @return: True/False

        SAM format used for this output:
        Col Field Type Regexp/Range Brief description
        1 QNAME String [!-?A-~]{1,254} Query template NAME = The non-reftype sequence name
        2 FLAG Int [0,216-1] bitwise FLAG = The first alignment for a query should be 0; The rest 256. Add 16 for reverse hits?
        3 RNAME String \*|[!-()+-<>-~][!-~]* Reference sequence NAME = reftype sequence name
        4 POS Int [0,231-1] 1-based leftmost mapping POSition
        5 MAPQ Int [0,28-1] MAPping Quality = Set to 255 (unavailable)
        6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
        7 RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
        8 PNEXT Int [0,231-1] Position of the mate/next read
        9 TLEN Int [-231+1,231-1] observed Template LENgth
        10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
        11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not reftype: reftype = self.getStr('RefType')    # Default=Hit
            if not filename: filename = '%s.sam' % self.baseFile()
            if reftype.lower() == 'both':
                self.saveSAM('%s.qry.sam' % rje.baseFile(filename),locdb,queries,hits,'Qry',cigar,sequence,append=False)
                return self.saveSAM('%s.hit.sam' % rje.baseFile(filename),locdb,queries,hits,'Hit',cigar,sequence,append=False)
            elif reftype.lower() == 'combined':
                self.saveSAM(filename,locdb,queries,hits,'Qry',cigar,sequence,append=False)
                return self.saveSAM(filename,locdb,queries,hits,'Hit',cigar,sequence,append=True)
            db = self.db()
            ## ~ [0a] Generate SAM Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not locdb: locdb = self.db('Local')
            samdb = db.copyTable(locdb,'sam')
            seqfields = 'QrySeq' in samdb.fields() and 'SbjSeq' in samdb.fields()
            if 'Query' in samdb.fields(): samdb.renameField('Query','Qry')
            if 'AlnNum' in samdb.fields(): samdb.renameField('AlnNum','AlnID')
            if reftype.lower() in ['query','qry']: reftype = 'Qry'
            elif reftype.lower() in ['hit']: reftype = 'Hit'
            else: raise ValueError('reftype=%s not recognised' % reftype)
            qrytype = {'Hit':'Qry','Qry':'Hit'}[reftype]
            ## ~ [0b] Cleanup SAM Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if queries: samdb.dropEntriesDirect('Qry',queries,inverse=True,log=self.debugging(),force=False)
            if hits: samdb.dropEntriesDirect('Hit',hits,inverse=True,log=self.debugging(),force=False)
            samdb.dataFormat({'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Identity':'int'})
            samdb.renameField('SbjStart','HitStart')
            samdb.renameField('SbjEnd','HitEnd')
            samdb.renameField('SbjSeq','HitSeq')
            qlist = samdb.indexKeys(qrytype)
            # Generate SAM format fields
            #samdb.renameField(qrytype,'QNAME')
            samdb.addField('QNAME',evalue='')
            samdb.addField('FLAG',evalue=256)
            #samdb.renameField(reftype,'RNAME')
            #samdb.addField('POS',evalue=0)
            samdb.addField('MAPQ',evalue=255)
            samdb.addField('CIGAR',evalue='*')
            samdb.addField('RNEXT',evalue='*')
            samdb.addField('PNEXT',evalue=0)
            samdb.addField('TLEN',evalue=0) #!# At some point, learn how to make multi-segment mappings
            samdb.addField('SEQ',evalue='*')
            samdb.addField('QUAL',evalue='*')
            samdb.addField('NOTES',evalue='')
            qlist = []
            for entry in samdb.sortedEntries('Identity',reverse=True):
                self.bugPrint(entry)
                revhit = entry['HitStart'] > entry['HitEnd']
                qlen = max(entry['%sStart' % qrytype],entry['%sEnd' % qrytype])   #!# Replace with actual length
                #?# Having different sequences with the same name might be messing things up!
                #!# Replace with Name.Start-End and *add* to qlist to assign primary/secondary #!#
                #!# Check lengths from other table(s) and only add start/end if not full length
                if revhit and reftype == 'Hit': qname = '%s.%s-%s' % (entry[qrytype],rje.preZero(entry['%sEnd' % qrytype],qlen),rje.preZero(entry['%sStart' % qrytype],qlen))
                else: qname = '%s.%s-%s' % (entry[qrytype],rje.preZero(entry['%sStart' % qrytype],qlen),rje.preZero(entry['%sEnd' % qrytype],qlen))
                if qname not in qlist: qlist.append(qname); entry['FLAG'] = 0
                entry['QNAME'] = qname

                if revhit: entry['FLAG'] += 16     # Reverse complement sequence
                #i# If the Reference Genome is the Hit, rev hits need to be reverse complemented pre-CIGAR generation.
                if seqfields and revhit and reftype == 'Hit':
                    entry['QrySeq'] = rje_sequence.reverseComplement(entry['QrySeq'])
                    entry['HitSeq'] = rje_sequence.reverseComplement(entry['HitSeq'])
                    (entry['QryStart'],entry['QryEnd']) = (entry['QryEnd'],entry['QryStart'])
                    (entry['HitStart'],entry['HitEnd']) = (entry['HitEnd'],entry['HitStart'])
                # Generate CIGAR string
                if seqfields and cigar: entry['CIGAR'] = self.aln2CIGAR(entry['%sSeq' % qrytype],entry['%sSeq' % reftype])
                #i# If the Reference Genome is the Qry, rev hits need to be reverse complemented post-CIGAR generation.
                #if seqfields and revhit and reftype == 'Qry':
                #    entry['HitSeq'] = rje_sequence.reverseComplement(entry['HitSeq'])
                #!# Or do it for reftype=Hit because the SEQ *should* be the revcomp sequence versus the reference!
                if seqfields and revhit and reftype == 'Hit':
                    entry['HitSeq'] = rje_sequence.reverseComplement(entry['HitSeq'])
                if seqfields and sequence: entry['SEQ'] = string.replace(entry['%sSeq' % qrytype],'-','')
                # Add notes
                for field in ['BitScore','Expect','Length','Identity']:
                    entry['NOTES'] += '%s=%s;' % (field.lower(),entry[field])
                self.debug(entry)
            samdb.renameField('%sStart' % reftype,'POS')
            #samdb.renameField(qrytype,'QNAME')
            samdb.renameField(reftype,'RNAME')
            # Tidy fields
            samfields = ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','NOTES']
            samdb.newKey(['RNAME','POS','QNAME','AlnID'])
            samdb.keepFields(samfields+['AlnID'])
            ### ~ [2] Output SAM File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            samdb.saveToFile(filename,delimit='\t',backup=True,append=append,savefields=samfields,log=True,headers=False)
        except: self.errorLog('%s.saveSAM() error' % self.prog()); return None
#########################################################################################################################
    def aln2CIGAR(self,qryseq,trgseq):    ### Convert alignment into CIGAR string and returns.
        '''
        Convert alignment into CIGAR string and returns.
        >> qryseq:str = Alignment string of the query sequence.
        >> trgseq:str = Alignment string of the target sequence.
        << cigar:str = CIGAR string returned
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# These are CIGAR elements. Note that the "query" is the "read" and the "target" is the "genome".
            # M	 	Match (alignment column containing two letters). This could contain two different letters (mismatch) or two identical letters. USEARCH generates CIGAR strings containing Ms rather than X's and ='s (see below).
            # D	 	Deletion (gap in the target sequence).
            # I	 	Insertion (gap in the query sequence).
            # S	 	Segment of the query sequence that does not appear in the alignment. This is used with soft clipping, where the full-length query sequence is given (field 10 in the SAM record). In this case, S operations specify segments at the start and/or end of the query that do not appear in a local alignment.
            # H	 	Segment of the query sequence that does not appear in the alignment. This is used with hard clipping, where only the aligned segment of the query sequences is given (field 10 in the SAM record). In this case, H operations specify segments at the start and/or end of the query that do not appear in the SAM record.
            # =	 	Alignment column containing two identical letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
            # X	 	Alignment column containing a mismatch, i.e. two different letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
            #i# Using =/X rather than M - could add switch some time for this. S & H not used.
            if len(qryseq) != len(trgseq): raise ValueError('Alignments must be sequences of same length!')
            cigel = [0,'']  # [Element,count] of current CIGAR element to be considered.
            cigar = ''
            ### ~ [1] Generate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in xrange(len(qryseq)):
                qi = qryseq[i].upper()
                ti = trgseq[i].upper()
                if qi == ti: ci = '='
                elif qi == '-': ci = 'D'
                elif ti == '-': ci = 'I'
                else: ci = 'X'
                if ci == cigel[1]: cigel[0] += 1
                else:
                    if cigel[0]: cigar += '%d%s' % tuple(cigel)
                    cigel = [1,ci]
            if cigel[0]: cigar += '%d%s' % tuple(cigel)
            return cigar
        except: self.errorLog('%s.aln2CIGAR() error' % self.prog()); return None
#########################################################################################################################
### SECTION IIb: OLD BLASTRun Class                                                                                     #
#########################################################################################################################
    ### <2> ### BLAST Search                                                                                            #
#########################################################################################################################
    def OLDreadBLAST(self,resfile=None,clear=False,gablam=False,unlink=False,local=False,screen=True,log=False):       
        '''
        Reads BLAST Results into objects.
        >> resfile:str = Results File (set as self.info['Name'])
        >> clear:Boolean = whether to clear current searches (True) or just append (False) [False]
        >> gablam:Boolean = whether to calculate gablam statistics and clear alignments to save memory [False]
        >> unlink:Boolean = whether to delete BLAST results file after reading [False]
        >> local:Boolean = whether to store Local alignment dictionary with basic alignment data [False]
        >> screen:Bool [False] = whether to output reading details to screen
        >> log:Bool [False] = whether to output reading details to log
        << returns True if (apparently) read to completion OK, else False
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if resfile != None: self.info['Name'] = resfile
            resline = self.loadFromFile(filename=self.info['Name'],v=1,checkpath=False,chomplines=True)
            try:
                newtype = string.split(resline[0])[0].lower()
                if newtype in ['blastn','blastp','blastx','tblastn','tblastx','rpsblast','rpstblastn']: self.info['Type'] = newtype
            except:
                if not resline: self.errorLog('No lines read in from %s.' % self.info['Name'],printerror=False)
                else: self.errorLog('Error with lines read in from %s.' % self.info['Name'],printerror=False)
                raise IOError
            if clear: self.search = []
            ### ~ [2] ~ Read in Search results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            logtxt = 'Reading %s %s BLAST results' % (self.info['Name'],self.info['Type']); (sx,hx) = (0,0)
            self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
            search = None; readhits = False
            i = 0
            hitaln = 0
            while i < len(resline):
                line = resline[i]
                ## ~ [2a] ~ Basic Search Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.find('Query=') == 0:    # New Search
                    search = self._addSearch(); readhits = False
                    search.info['Name'] = rje.matchExp('^Query=\s+(\S+)', line)[0]
                    sx += 1
                elif re.search('^\s+\(\S+ letters\)',line):
                    len_match = string.replace(rje.matchExp('^\s+\((\S+) letters\)',line)[0],',','')
                    search.stat['Length'] = string.atoi(len_match)
                elif line.find('Number of letters in database:') >= 0:
                    dblen = rje.matchExp('Number of letters in database:\s+(\d\S*)', line)[0]
                    dblen = re.sub('\D','',dblen)
                    self.stat['DBLen'] = string.atoi(dblen)
                elif line.find('Number of sequences in database:') >= 0:
                    dbnum = rje.matchExp('Number of sequences in database:\s+(\d\S*)', line)[0]
                    dbnum = re.sub('\D','',dbnum)
                    self.stat['DBNum'] = string.atoi(dbnum)
                elif rje.matchExp('(\S+) sequences; (\S+) total letters',line):
                    (dbnum,dblen) = rje.matchExp('(\S+) sequences; (\S+) total letters',line)
                    self.stat['DBNum'] = string.atoi(re.sub('\D','',dbnum))
                    self.stat['DBLen'] = string.atoi(re.sub('\D','',dblen))
                elif line.find('Number of sequences better than') >= 0:
                    self.stat['E-Value'] = string.atof(rje.matchExp('Number of sequences better than\s+(\S+):', line)[0])
                ## ~ [2b] ~ One-line hit data (BLASTHit) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('Sequences producing significant alignments:') >= 0: # One-line hits
                    if string.split(line)[-1] == 'N': (si,ei) = (-3,-2)
                    else: (si,ei) = (-2,-1)
                    i += 2  # Skip blank line
                    while rje.matchExp('^(\S+)\s.*\s(\S*\d)\s+(\S*\d)\s*$',resline[i]):
                        #match = rje.matchExp('^(\S+)\s.*\s(\d\S*)\s+(\S+\d)\s*$',resline[i])
                        match = string.split(resline[i])
                        hit = search._addHit()
                        hit.setInfo({'Name':match[0],'Type':self.info['Type']})
                        hit.stat['BitScore'] = string.atof(match[si])
                        eval = match[ei]
                        if eval.find('e') == 0: eval = '1' + eval
                        hit.stat['E-Value'] = string.atof(eval)
                        i += 1
                        hx += 1
                    line = resline[i]   # End of one-lines (blank line)
                    self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
                    search.checkHitNames()
                    hitaln = 0; readhits = True
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    search.hit = []
                    self.progLog('\r#BLAST','%s: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx)),screen=screen)
                    hitaln = 0; readhits = True
                ## ~ [2c] ~ Aln Hit data (PWAln) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif line.find('>') == 0:
                    if not readhits: i += 1; continue   # '>' Character within name happened to occur at start of line
                    hitname = rje.matchExp('^>(\S+)\s*',line)[0]
                    if hitaln >= len(search.hit):
                        self.errorLog('Apparent hits for %s have exceeded the %d found. (Hit %s.) BLAST read-through or sequence name format error?' % (search.info['Name'],hitname,len(search.hit)),False,False)
                        raise ValueError
                    #if hitaln not in search.hit:
                    #    self.errorLog('Cannot find hitaln %s for %s hits. Sequence name format error?' % (hitaln,search.info['Name']))
                    #    raise ValueError
                    if hitname != search.hit[hitaln].info['Name']:      # Identify hit object
                        for hit in search.hit:
                            if hit.info['Name'] == hitname: hitaln = search.hit.index(hit)
                        if hitname != search.hit[hitaln].info['Name']:
                            self.errorLog('Problem with BLAST results - %s single-line hits and alignments do not match' % search.info['Name'],printerror=False,quitchoice=True)
                            i += 1; continue
                    hit = search.hit[hitaln]
                    hitaln += 1
                    aln = None
                    while i < (len(resline) - 1):
                        i += 1
                        line = resline[i]
                        ## Hit Length ##
                        if re.search('^\s+Length\s+=\s+(\d+)',line):
                            hit.stat['Length'] = string.atoi(rje.matchExp('^\s+Length = (\d+)',line)[0])
                        ## New Aln Block ##
                        elif hit.stat['Length'] and rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line):
                            aln = hit._addAln()
                            scores = rje.matchExp('^\s*Score\s+=\s+(\S+)\s.+Expect.+=\s+(\S+\d)\D*\s*',line)
                            aln.stat['BitScore'] = string.atof(scores[0])
                            eval = scores[1]
                            if eval.find('e') == 0: eval = '1' + eval
                            aln.stat['Expect'] = string.atof(eval)
                            i += 1
                            if re.search('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives = (\d+)/(\d+)\s?',resline[i]):
                                sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?.+Positives\s+=\s+(\d+)/(\d+)\s?',resline[i])
                                aln.stat['Positives'] = string.atoi(sim[2])
                            else: sim = rje.matchExp('Identities\s+=\s+(\d+)/(\d+)\s?',resline[i])
                            try:
                                aln.stat['Length'] = string.atoi(sim[1])
                                aln.stat['Identity'] = string.atoi(sim[0])
                            except: self.errorLog('Problem reading line "%s"' % resline[i])
                            if resline[i+1].find('Frame') >= 0: i += 1    # BLASTX or TBLASTN
                        ## Alignment ##
                        elif line.find('Query:') == 0 and aln != None:
                            # Query Line
                            alnstructure = rje.matchExp('^(Query:\s+\d+\s+)(\S+)\s+(\d+)',line)
                            if not alnstructure:    # Something's wrong!
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                aln._clear()
                                break
                            leader = len(alnstructure[0])
                            aln.info['QrySeq'] += alnstructure[1]
                            aln.stat['QryEnd'] = string.atoi(alnstructure[2])
                            if aln.stat['QryStart'] == 0:
                                aln.stat['QryStart'] = string.atoi(rje.matchExp('^Query:\s+(\d+)\s',alnstructure[0])[0])
                            # Alignment Line
                            i += 1; aln.info['AlnSeq'] += resline[i][leader:(leader+len(alnstructure[1]))]
                            # Subject Line
                            i += 1; subject = rje.matchExp('^Sbjct:\s+(\d+)\s*(\D+)\s+(\d+)',resline[i])
                            if not subject:
                                self.errorLog('Problem reading aln for %s hit %s.' % (search.info['Name'],hit.info['Name']),printerror=False)
                                aln._clear()
                                break
                            aln.info['SbjSeq'] += subject[1]
                            aln.stat['SbjEnd'] = string.atoi(subject[2])
                            if aln.stat['SbjStart'] == 0: aln.stat['SbjStart'] = string.atoi(subject[0])
                        ## End of Aln ##
                        elif line.find('>') == 0:
                            i -= 1  # i Advanced anyway, so subtract 1
                            break   # so '>' will still be caught as next Aln
                        ## End of All Alns ##
                        elif line.find('Database:') >= 0 and hit.alnNum() > 0: break   # hit.alnNum() to avoid matching 'Database:' in hit description
                        elif line.find('BLAST') == 0: break
                        elif line.find('TBLASTN') == 0: break
                ## ~ [2d] ~ Continue reading ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##      
                i += 1
            ptxt = '%s complete: %s searches; %s hits' % (logtxt,rje.integerString(sx),rje.integerString(hx))
            self.printLog('\r#BLAST','%s (vs %s Sequences of total %s letters)' % (ptxt,rje.integerString(self.stat['DBNum']),rje.integerString(self.stat['DBLen'])),log=log,screen=screen)
            ### ~ [3] ~ Local alignment options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if local:
                self.progLog('#LOCAL','Generating local aln dictionaries...',screen=screen)
                for search in self.search:
                    for hit in search.hit: hit.makeLocalDict()
                self.progLog('\r#LOCAL','Generation of local aln dictionaries complete!',screen=screen)
            ### ~ [4] ~ GABLAM calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if gablam:  
                logtxt = 'Calculating GABLAM statistics from %s BLAST results' % self.info['Name']
                _hitx = 0
                sx = 0.0
                for search in self.search:
                    self.progLog('\r#GABLAM','%s: %.1f%%' % (logtxt,sx/len(self.search)),screen=screen)
                    _hitx += search.hitNum()
                    search.gablam()
                    sx += 100.0
                logtxt = 'Calculation of GABLAM statistics from %s BLAST results (%s hits) complete.' % (self.info['Name'],rje.integerString(_hitx))
                self.printLog('\r#GABLAM',logtxt,log=log,screen=screen)
            ### ~ [5] ~ Unlink if necessary and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if unlink and os.path.exists(self.info['Name']): os.unlink(self.info['Name'])
            return True
        except:
            if self.checkBLAST():
                self.errorLog('Fatal Error during BLASTRun.readBLAST() despite intact BLAST results.')
                raise
            else:
                self.errorLog('BLASTRun.readBLAST() failed due to incomplete BLAST results file.')
                return False
#########################################################################################################################
    def _addSearch(self):   ### Adds and returns a new search object
        '''Adds and returns a new search object.'''
        newsearch = BLASTSearch(log=self.log,cmd_list=self.cmd_list)#; newsearch.setOpt({'DeBug':self.getBool('DeBug')})
        self.search.append(newsearch)
        return newsearch
#########################################################################################################################
    ### <3> ### Output                                                                                                  #
#########################################################################################################################
    def saveCutBLAST(self,outfile=None):   ### Saves a cutdown version of current BLAST (no alignments)
        '''
        Saves a cutdown version of current BLAST (no alignments), which has enough lines to be successfully read in by
        the BLASTRun class. (Obviously, no GABLAM statistics can be calculated!
        >> outfile:str = outfile to use if different from self.info['Name']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfile: self.info['Name'] = outfile
            CUT = open(self.info['Name'], 'w')
            ### ~ [2] ~ Output cut-down search data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for search in self.search:
                CUT.write('Query= %s\n\n' % search.info['Name'])
                CUT.write('         (%d letters)\n' % search.stat['Length'])
                CUT.write('Number of letters in database: %d\n' % self.stat['DBLen'])
                CUT.write('Number of sequences in database: %d\n' % self.stat['DBNum'])
                CUT.write('Number of sequences better than %e:\n' % self.stat['E-Value'])
                CUT.write('\nSequences producing significant alignments:\n\n')
                if search.hitNum() > 0:
                    for hit in search.hit:
                        CUT.write('%s  %.1f  %e\n' % (hit.info['Name'],hit.stat['BitScore'],hit.stat['E-Value']))
                else: CUT.write('***** No hits found ******\n')
                CUT.write('\n\n')
            CUT.close()
        except: self.errorLog('Major error during BLASTRun.saveCutBLAST().',quitchoice=True)
#########################################################################################################################
    def searchSeq(self,seqlist,proglog=True,inverse=False):   ### Returns dictionary of searches as sequences from seqlist (or None if missing)
        '''
        Returns dictionary of searchesas sequences from seqlist (or None if missing).
        >> seqlist:SeqList object
        >> proglog:bool [True] = whether to log progress 
        >> inverse:bool [False] = whether to reverse dictionary to sequence:search
        << {Search Object: Sequence Object}
        '''
        try:### ~ [1] ~ Setup (Clear objects) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            searchdic = {}
            for search in self.search: searchdic[search] = None
            ### ~ [2] ~ Map SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist:
                seqdic = seqlist.seqNameDic(proglog=False)  #proglog)
                (mx,hx) = (0,0.0)
                for search in self.search:
                    hx += 100.0
                    if proglog:
                        ltxt = 'Mapping %s searches onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
                        self.printLog('\r#SEARCH','%s: %.2f%%' % (ltxt,(hx/self.searchNum())),newline=False,log=False)
                    if seqdic.has_key(search.info['Name']):
                        searchdic[search] = seqdic[search.info['Name']]
                        mx += 1
                    else:
                        for key in seqdic.keys():
                            if key.find('|%s' % search.info['Name']) > 0:
                                searchdic[search] = seqdic[key]
                                mx += 1
                                break   #!# Add report_none=False?
                if proglog:
                    ltxt = 'Mapping %s searches onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
                    self.printLog('\r#SEARCH','%s: %s of %s mapped.' % (ltxt,rje.integerString(mx),rje.integerString(self.searchNum())))
            ### ~ [3] ~ Finish and return (empty if no seqlist) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if inverse:
                invdic = {}
                for search in searchdic:
                    if searchdic[search]: invdic[searchdic[search]] = search
                return invdic
            return searchdic                        
        except: self.errorLog('Major Problem with BLASTSearch.searchSeq().'); raise       
#########################################################################################################################
### END OF SECTION II: BLASTRun Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###
#########################################################################################################################
### SECTION III: BLASTSearch Class                                                                                      #
#########################################################################################################################
    def hitSeq(self,seqlist,proglog=True,inverse=False):   ### Returns dictionary of hits as sequences from seqlist (or None if missing)
        '''
        Returns dictionary of hits as sequences from seqlist (or None if missing).
        >> seqlist:SeqList object
        >> proglog:bool [True] = whether to log progress 
        >> inverse:bool [False] = whether to reverse dictionary to sequence:hit
        << {Hit Object: Sequence Object}
        '''
        try:### ~ [1] ~ Setup (Clear) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hitdic = {}
            for hit in self.hit: hitdic[hit] = None
            ltxt = 'Mapping %s hits onto %s sequences' % (self.info['Name'],seqlist.info['Name'])
            ### ~ [2] ~ Map SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqlist:
                seqdic = seqlist.seqNameDic(proglog=proglog)
                (mx,hx) = (0,0.0)
                for hit in self.hit:
                    hx += 100.0
                    if proglog: self.progLog('\r#HIT','%s: %.2f%%' % (ltxt,(hx/self.hitNum())))
                    if seqdic.has_key(hit.info['Name']):
                        hitdic[hit] = seqdic[hit.info['Name']]
                        mx += 1
                    else:
                        for key in seqdic.keys():
                            if key.find('|%s' % hit.info['Name']) > 0:
                                hitdic[hit] = seqdic[key]
                                mx += 1
                                break   #!# Add report_none=False?
                if proglog: self.printLog('\r#HIT','%s: %s of %s mapped.' % (ltxt,rje.integerString(mx),rje.integerString(self.hitNum())))
            ### Return (empty if no seqlist) ###
            if inverse:
                invdic = {}
                for hit in hitdic:
                    if hitdic[hit]: invdic[hitdic[hit]] = hit
                    return invdic
            return hitdic                        
        except: self.errorLog('Major Problem with BLASTSearch.hitSeq().'); raise       
#########################################################################################################################
### END OF SECTION II: BLASTRun Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION VI: GENERAL MODULE METHODS                                                                                  #
#########################################################################################################################
def formatDB(fasfile,blastpath,protein=True,log=None,oldblast=False,details=False):  ### Formats a blastDB             #V2.1
    '''
    Formats a blastDB.
    >> fasfile:str = file to format
    >> blastpath:str = path to BLAST programs
    >> protein:boolean = whether protein sequences
    >> log:Log Object
    >> oldblast:bool = whether to use old formatdb commands (True) or new makeblastdb command (False) [False]
    >> details:bool = whether to output details of DB formatting to log file [False]
    '''
    try:### ~ [1] Setup Command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        path = rje.makePath(blastpath)
        ## ~ [1a] OldBLAST formatdb command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if oldblast:
            if os.sep == '\\': command = blastpath + 'formatdb -i "%s" -o T' % fasfile
            else: command = blastpath + 'formatdb -i %s -o T -a F' % fasfile
            #!# -a F option is to remove SeqPortNew errors with -o T (??!!) I don't really understand what this means.
            if protein: command += ' -p T'
            else: command += ' -p F'
        ## ~ [1b] New makeblastdb command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        else:
            command = blastpath + 'makeblastdb -in "%s" -out "%s" -title "%s" -parse_seqids' % (fasfile,fasfile,os.path.abspath(fasfile))
            if protein: command += ' -dbtype prot'
            else: command += ' -dbtype nucl'
            #?# Other options: -mask_data refseq_seg.asnb             
        ### ~ [2] Execute command ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if log: log.printLog('#DB ',command,screen=log.v()>0)
        for oline in os.popen(command):  #?# Use os.popen() and catch/output STDOUT to log #?#
            if rje.chomp(oline) and log: log.printLog('#NCBI',oline,log=details,screen=details and log.v()>0)
            elif rje.chomp(oline) and details: print oline
    except:
        if log: log.errorLog('Major Problem during rje_blast.formatDB(%s).' % fasfile)
        else: print 'Major Problem during rje_blast.formatDB(%s).' % fasfile
        raise
#########################################################################################################################
def checkForDB(dbfile=None,checkage=True,log=None,protein=True,oldblast=False):     ### Checks for BLASTDB files and returns True or False as appropriate
    '''
    Checks for BLASTDB files and returns True or False as appropriate.
    >> dbfile:str = sequence file forming basis of database
    >> checkage:Boolean = also check that blastdb files are newer than dbfile
    >> log:Log Object
    >> protein:boolean = whether database is protein
    >> oldblast:bool = whether to use old formatdb commands (True) or new makeblastdb command (False) [False]
    '''
    try:### ~ [1] ~ Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(dbfile):
            if log: log.errorLog('%s missing' % dbfile,False,False)
            return False
        if protein: suffix = ['phr','pin','psd','psi','psq','pog']
        else: suffix = ['nhr','nin','nsd','nsi','nsq','nog']
        if oldblast: suffix.pop(-1)
        missing = []
        for suf in suffix:
            if os.path.exists('%s.%s' % (dbfile,suf)):
                if checkage and rje.isYounger('%s.%s' % (dbfile,suf),dbfile) == dbfile:
                    if log: log.errorLog('%s.%s too old' % (dbfile,suf),False,False)
                    return False
            elif os.path.exists('%s.00.%s' % (dbfile,suf)):
                if checkage and rje.isYounger('%s.00.%s' % (dbfile,suf),dbfile) == dbfile:
                    if log: log.errorLog('%s.%s too old' % (dbfile,suf),False,False)
                    return False
            else: missing.append(suf)
        if missing:
            if log and len(missing) < 5: log.errorLog('%s.%s missing' % (dbfile,string.join(missing,'/')),False,False)
            return False
        return True
    except:
        if log: log.errorLog('Major Problem during rje_blast.checkForDB().')
        else: print 'Major Problem during rje_blast.checkForDB().'
        raise
#########################################################################################################################
def cleanupDB(callobj=None,dbfile=None,deletesource=False):     ### Deletes files created by formatdb
    '''
    Deletes files created by formatdb.
    >> callobj:Object = object calling method
    >> dbfile:str = sequence file forming basis of database
    >> deletesource:boolean = whether to delete dbfile as well
    '''
    try:
        for suf in ['phr','pin','psd','psi','psq','pal','nhr','nin','nsd','nsi','nsq','pni','pnd','pog','nog']:
            if os.path.exists('%s.%s' % (dbfile,suf)): os.unlink('%s.%s' % (dbfile,suf))
            for x in range(100):
                xstr = rje.preZero(x,99)
                if os.path.exists('%s.%s.%s' % (dbfile,xstr,suf)): os.unlink('%s.%s.%s' % (dbfile,xstr,suf))
        if deletesource and os.path.exists(dbfile): os.unlink(dbfile)
    except:
        if callobj: callobj.log.errorLog('Major Problem during rje_blast.cleanupDB().')
        else: print 'Major Problem during rje_blast.cleanupDB().'
        raise
#########################################################################################################################
def expectString(_expect): return rje.expectString(_expect)
#########################################################################################################################
### END OF SECTION VI                                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def blastObj(log,cmd_list,type='New'):   ### Returns appropriate BLAST Object                                       #V2.1
    '''
    Returns appropriate BLAST Object.
    - New/V2/BLAST+ = New Object for BLAST+
    - Old/V1/BLAST = Old Object for BLAST
    - Dev = Will return V1 if oldblast=T or V2 otherwise.
    - OldBLAST = New Object running old BLAST legacy commands. (NB. Old BLAST, *not* BLAST+ legacy options!)
    '''
    try:
        if type.lower() in ['new','v2','blast+']: return BLASTRun(log,cmd_list)
        if type.lower() in ['dev']:
            blast = BLASTRun(log,cmd_list)
            #!# For now, return BLAST_V1 Object if oldblast=T. Switch at some point to use Legacy within BLAST #!#
            if blast.getBool('OldBLAST'): return rje_blast_V1.BLASTRun(log,cmd_list)
            else: return blast
        if type.lower() in ['legacy']: return BLASTRun(log,cmd_list+['legacy=T'])   # Currently OldBLAST
        if type.lower() in ['old','v1','blast']: return rje_blast_V1.BLASTRun(log,cmd_list)
        raise ValueError('BLAST type not recognised for blastObj()')
    except:
        log.errorLog('Problem with BLAST Object initiation (type "%s")' % type)
        raise
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
    try: BLASTRun(mainlog,['tuplekeys=T']+cmd_list).run()

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
