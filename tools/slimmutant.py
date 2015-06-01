#!/usr/bin/python

# See below for name and description
# Copyright (C) 2013 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMMutant
Description:  Short Linear Motif Mutation Analysis Tool
Version:      1.3
Last Edit:    16/09/14
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMMutant is a Short Linear Motif Mutation Analysis Tool, designed to identify and assess mutations that create
    and/or destroy SLiMs. There are three main run modes:

    - Mode 1. Generating mutant datasets for analysis [`generate=T`]

    The main input is: (1) a file of protein sequence mutations [`mutfile=FILE`] in a delimited text format with aa
    substitution [`mutfield=X`] and protein accession number [`protfield=X`] data; (2) a corresponding sequence file
    [`seqin=FILE`]; (3) a file of SLiMs to analyse [`motifs=FILE`]. This will process the data and generate two sequence
    files: `*.wildtype.fas` and `*.mutant.fas`. These files will be named after the input `mutfile` unless `basefile=X`
    is used.

    - Mode 2. Run SLiMProb on datasets [`slimprob=T`]

    This will run SLiMProb on the two datasets, once per `*.ini` file given by `slimini=LIST`. These runs should have
    distinct `runid=X` settings. If no `*.ini` files are given, as single run will be made using commandline settings.

    - Mode 3. Compile results of SLiMProb runs. [`analyse=T`]

    This will compare the `*.wildtype.fas` and `*.mutant.fas` results from the `*.occ.csv` file produced by SLiMProb.
    All mutations analysed will be identified from `*.mutant.fas`. SLiM occurrences are then matched up between wildtype
    and mutant versions of the same sequence. If none of the mutations have effected the SLiM prediction, then the
    wildtype and all mutant sequences will return the motif. If, on the other hand, mutations have created/destroyed
    motifs, occurrences will be missing from the wildtype and/or 1+ mutant sequences. All unaffected SLiM instances are
    first removed and altered SLiM instances  output to `*.MutOcc.csv`. Differences between mutants and wildtypes are
    calculated for each `RunID`-`Motif` combination and summary results output to `*.Mut_vs_WT.csv`. If `motlist=LIST` is
    given, analysis is restricted to a subset of motifs.

    Unless `basefile=FILE` is given, output files will be named after `mutfile=FILE` but output into the current run
    directory. If running in batch mode, basefile cannot be used.

    NOTE: SLiMMutant is still in development and has not been thoroughly tested or benchmarked.

Commandline:
    ### ~ SEQUENCE GENERATION METHODS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    generate=T/F	: Whether to run sequence generation pipeline [False]
    mutfile=FILE	: Delimited text file with sequence mutation info. Sets basefile. []
    mutfield=X		: Field in mutfile corresponding to AA subsitution data ['AAChange']
    protfield=X		: Field in mutfile corresponding to protein accession number ['Uniprot']
    splitfield=X    : Field in mutfile to split data on (saved as basefile.X.tdt) []
    seqin=FILE		: Input file with protein sequences []
    motifs=FILE		: Input file of SLiMs []
    motlist=LIST	: List of input SLiMs to restrict analysis to []
    mutflanks=X		: Generate for casemask=Upper of X aa flanking mutation (None if < 1) [0]
    minmutant=X     : Minimum number of mutants for output [100]
    maxmutant=X     : Maximum number of mutants for output [100000]

    ### ~ SLiMProb Run Methods ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    slimprob=T/F	: Whether to run SLiMProb on *.wildtype.fas and *.mutant.fas (*=basefile) [False]
    slimini=LIST	: Lists of INI file with settings for SLiMProb run. Should include runid=X and resdir=PATH. []
    resdir=PATH   	: Location of output files. SLiMProb resdir should be in slimini [SLiMMutant/ (and SLiMProb/)]

    ### ~ SLiMProb Results Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    analyse=T/F		: Whether to analyse the results of a SLiMProb run [False]
    resfile=FILE   	: Main SLiMProb results table (*.csv and *.occ.csv) [slimprob.csv]
    runid=X			: Limit analysis to SLiMProb RunID (blank = analyse all) []
    buildpath=PATH 	: Alternative path to look for existing intermediate files (e.g. *.upc) [SLiMProb/]

    ### ~ SLiMProb PPI Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    slimppi=T/F         : Whether to perform SLiMPPI analysis (will set analyse=T) [False]
    sourcepath=PATH/    : Will look in this directory for input files if not found ['SourceData/']
    sourcedate=DATE		: Source file date (YYYY-MM-DD) to preferentially use [None]
    ppisource=X			: Source of PPI data. (HINT/FILE) FILE needs 'Hub' and 'SpokeUni' fields. ['HINT']
    dmifile=FILE        : Delimited text file containing domain-motif interaction data ['elm_interaction_domains.tsv']

    ### ~ Batch running ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    batch=FILELIST  : List of mutfiles to run in batch mode. Wildcards allowed. []

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import math, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seq, rje_seqlist, rje_slim, rje_slimlist, rje_slimcore, rje_uniprot
import pingu_V4
import slimprob
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Working version with standalone functionality.
    # 1.1 - Minor tweaks to generate method to increase speed. (Make index in method.) Added splitfield=X.
    # 1.2 - Added a batch mode for mutfiles - all other options will be kept fixed. Added maxmutant and minmutant.
    # 1.3 - Added SLiMPPI analysis (will set analyse=T). Started basing on SLiMCore
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add option to treat combine all mutations for given sequence. (Currently each treated separately.)
    # [ ] : Add option to include indels. (Currently only single base substitutions.)
    # [ ] : Add option to extract UniProt AccNum and convert to sequences.
    # [Y] : Consider option to restrict analysis to window around mutations. (Much more efficient SLiMProb search!)
    # [ ] : Add option to use Ensembl transcript -> peptide mapping.
    # [Y] : Update to somehow use Megaslim wildtype RLC scores for mutants too. (Also disorder?)
    # [X] : Convert splitfield=X to splitfields=LIST? (But how to manage output directories?
    # [ ] : Remove flank analysis and dev only run once the new SLiMSearch method is up and running.
    # [ ] : Add a method that will check the mutations for altering the disorder state of the protein.
    # [ ] : Consider adding protein weighting to mutation frequency, where each protein gets equal weight (not each mut)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SLiMMutant', '1.3', 'September 2014', '2014')
    description = 'Short Linear Motif Mutation analysis'
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
### SECTION II: SLiMMutant Class                                                                                        #
#########################################################################################################################
class SLiMMutant(rje_slimcore.SLiMCore):
    '''
    SLiMMutant Class. Author: Rich Edwards (2014).

    Str:str
    - BuildPath = Alternative path to look for existing intermediate files (e.g. *.upc) [SLiMProb/]
    - DMIFile=FILE        : Delimited text file containing domain-motif interaction data ['elm_interaction_domains.tsv']
    - Motifs=FILE		: Input file of SLiMs []
    - MutFile=FILE	: Delimited text file with sequence mutation info []
    - MutField=X		: Field in mutfile corresponding to AA subsitution data ['AAChange']
    - ProtField=X		: Field in mutfile correspoding to protein accession number ['Uniprot']
    - ResDir = Location of output files. SLiMProb resdir should be in slimini [SLiMMutant/ (and SLiMProb/)]
    - ResFile = Main SLiMProb results table (*.csv and *.occ.csv) [slimprob.csv]
    - RunID = Limit analysis to SLiMProb RunID (blank = analyse all) [SLiMMutant]
    - Seqin=FILE		: Input file with protein sequences []
    - SplitField=X    : Field in mutfile to split data on (saved as basefile.X.tdt) []

    Bool:boolean
    - Analyse=T/F		: Whether to analyse the results of a SLiMProb run [False]
    - Generate=T/F	: Whether to run sequence generation pipeline [False]
    - SLiMPPI=T/F         : Whether to perform SLiMPPI analysis (will set analysis=T) [False]
    - SLiMProb=T/F	: Whether to run SLiMProb on *.wildtype.fas and *.mutant.fas (*=basefile) [False]

    Int:integer
    - MaxMutant = Maximum number of mutants for output [100000]
    - MinMutant = Minimum number of mutants for output [100]
    - MutFlanks = Generate for casemask=Upper of X aa flanking mutation (-1 = None) [0]

    Num:float
    
    List:list
    - Batch=FILELIST  : List of mutfiles to run in batch mode. Wildcards allowed. []
    - MotList = List of input SLiMs to restrict analysis to []
    - SlimINI = Lists of INI file with settings for SLiMProb run. Should include runid=X. []

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database Main Database object for loading in Mutation Data
    - PINGU = pingu_V4.PINGU object for managing/parsing PPI data.
    - SeqList = rje_seqlist.SeqList object for protein sequence data
    - SLiMCore = rje_slimcore.SLiMCore object for reading in UPC etc.
    - SlimList = rje_slimlist.SLiMList object for handling searched motifs
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['DMIFile','MutFile','MutField','ProtField','RunID','ResFile','SplitField']
        self.boollist = ['Analyse','Generate','SLiMPPI','SLiMProb']
        self.intlist = ['MaxMutant','MinMutant','MutFlanks']
        self.numlist = []
        self.listlist = ['Batch','MotList','SlimINI']
        self.dictlist = []
        self.objlist = ['PINGU']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
        self.setStr({'DMIFile':'elm_interaction_domains.tsv','MutField':'AAChange','ProtField':'Uniprot'})
        self.setBool({})
        self.setInt({'MaxMutant':100000,'MinMutant':100,'MutFlanks':0})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        #self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,['delimit=,']+self.cmd_list)
        self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F'])
        self.obj['SLiMCore'] = rje_slimcore.SLiMCore(self.log,self.cmd_list)
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list+['autoload=F'])
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self.coreCmd()  #!# Remove these core commands from repetition below #!#
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['MutField','ProtField','RunID','SplitField'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['DMIFile','MutFile','ResFile'])  # String representing file path
                self._cmdReadList(cmd,'bool',['Analyse','Generate','SLiMPPI','SLiMProb'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxMutant','MinMutant','MutFlanks'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['MotList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Batch','SlimINI']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getBool('SLiMPPI'): self.setBool({'Analyse':True})
        #if self.getInt('MutFlanks'): self.printLog('#FLANK','MutFlank Setting not currently available.')
        #self.setInt({'MutFlanks':0})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['Batch']: return self.batchRun()
            elif not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Generate') and not self.generate(): return False
            if self.getBool('SLiMProb') and not self.slimProb(): return False
            if self.getBool('Analyse') and not self.analyse(): return False
            return True
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def batchRun(self): ### Runs in batch mode
        '''Runs in batch mode.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if not self.basefile(return_none=''): self.basefile(self.getStr('ResDir').lower())  #!# Strip path!
            if not self.setup(): return False
            self.obj['SlimList'].loadMotifs()
            self.printLog('#BATCH','%s batch files recognised for run.' % rje.iLen(self.list['Batch']))
            self.printLog('#RUN','Generate:%s; SLiMSearch:%s; Analyse:%s; SLiMPPI:%s' % (self.getBool('Generate'),self.getBool('SLiMProb'),self.getBool('Analyse'),self.getBool('SLiMPPI')))
            ### ~ [2] ~ Individual SLiMMutant runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bx = 0
            mvwfields = []
            for batch in self.list['Batch']:
                bx += 1
                self.printLog('#BATCH','# ~~~~~~~~~~~~ Batch %s of %s: %s ~~~~~~~~~~~~ #' % (rje.iStr(bx),rje.iLen(self.list['Batch']),batch))
                slimmutant = SLiMMutant(self.log,self.cmd_list)
                slimmutant.setStr({'Basefile':'','MutFile':batch})
                slimmutant.list['Batch'] = []
                slimmutant.obj['PINGU'] = self.obj['PINGU']
                slimmutant.obj['SlimList'] = self.obj['SlimList']
                slimmutant.run()
                if self.getBool('SLiMPPI'): mvwfields = slimmutant.db('Mut_vs_WTDMI').fields()
                self.printLog('#BATCH','# ~~~~~~~~~~~~ End of %s batch run ~~~~~~~~~~~~ #' % (batch))
            #!# Add compilation: ResDir/ and basefile
            ### ~ [3] ~ Combine SLiMMutant runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('SLiMPPI'): return True
            if not mvwfields: return False
            cdb = self.db().addEmptyTable('bymutfile',['MutFile']+mvwfields,['MutFile','RunID','Motif','DomLink'])
            # RunID,Motif,Pattern,Mut,WT,DomLink,fMut,eMut,pMut,pWT
            for batch in self.list['Batch']:
                inbase = rje.baseFile(batch,strip_path=True)    #!# Will need to add resdir
                mvwfile = '%s.Mut_vs_WTDMI.csv' % inbase
                mdb = self.db().addTable(mvwfile,mainkeys=['RunID','Motif','DomLink'],name=mvwfile)
                mdb.dataFormat({'Mut':'int','WT':'int','fMut':'num','eMut':'num'})
                for entry in mdb.entries():
                    entry['eMut'] = entry['eMut'] * (entry['WT']+entry['Mut'])  # Now an expected number!
                    entry['MutFile'] = batch
                    cdb.addEntry(entry)
            mutdb = cdb
            motdb = self.db().copyTable(cdb,'bymotif',replace=True,add=True)
            mutdb.compress(['RunID','MutFile','DomLink'],default='sum')

            #!# Modify to use the post-masking mutfreq? (Or have as an option?)
            #!# Mask out mutations in globular domains at start? Report numbers?
            #!# Can the recurrence of mutations be incoporated as a filter/weighting?

            motdb.compress(['RunID','Motif','DomLink'],default='sum')
            #!# Calulate p* based on obs versus expected (Poisson)
            for cdb in [mutdb,motdb]:
                cdb.addField('GOF'); cdb.addField('GOFadj')
                for entry in cdb.entries():
                    try:
                        entry['GOF'] = math.log(max(0.5,entry['Mut'])/max(0.5,float(entry['WT'])),2)
                        entry['GOFadj'] = math.log(max(0.5,entry['Mut'])/max(0.5,float(entry['eMut'])),2)
                        entry['pMut'] = rje.logPoisson(entry['Mut'],entry['eMut'],callobj=self)
                        entry['pWT'] = rje.logPoisson(entry['WT'],(entry['WT']+entry['Mut'])-entry['eMut'],callobj=self)
                    except:
                        self.debug(entry)
                        raise
                cdb.saveToFile()


            return True
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.basefile(return_none='') and self.getStrLC('MutFile'):
                self.setBasefile(rje.baseFile(self.getStr('MutFile'),strip_path=True))
                self.printLog('#BASE','Basefile from MutFile: %s.' % self.basefile())
            ### ~ [2] Setup additional analysis data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['PINGU'] and self.getBool('SLiMPPI'):   # Might have been set up by batch run
                ## ~ [2a] Setup PINGU object and load PPI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pingu = self.obj['PINGU'] = pingu_V4.PINGU(self.log,self.cmd_list+['ppispec=HUMAN'])
                pingu.setupSourceData(setup_ppi=True,setup_dmi=False)
                self.obj['DB'] = pingu.db()
                ppidb = pingu.db('PPISource')
                if not ppidb: ppidb = pingu.db('PPI.HUMAN')
                if not ppidb: self.printLog('#FAIL','Cannot execute SLiMPPI run without PPI Data.'); return False
                ## ~ [2b] Load Uniprot-Pfam link data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                spec = 'HUMAN'
                ppidat = pingu.sourceDataFile('Uniprot.%s' % spec,expect=True,ask=True)  # If this exists, will use
                if not ppidat: self.warnLog('Something went wrong making/finding %s UniProt file.' % spec, quitchoice=pingu.getBool('Integrity'))
                pfamfile = rje.baseFile(ppidat) + 'uniprot-pfam.tdt'
                if rje.exists(pfamfile): pfamdb = pingu.db().addTable(pfamfile,mainkeys=['Uniprot'],name='Pfam')
                else:
                    ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=pfam'])
                    ppuni.setStr({'Name':ppidat})
                    #!# Replace this with generating/ready Pfam domain table? #!#
                    ppuni.readUniProt()     # Pfam domains stored in entry.dict['DB']['Pfam'] as as dictionary of {PfamID:Details}.
                    pfamdb = pingu.db().addEmptyTable('Pfam',['Uniprot','Pfam'],['Uniprot'])
                    for entry in ppuni.entries():
                        if 'Pfam' in entry.dict['DB'] and entry.dict['DB']['Pfam']:
                            pfamdb.addEntry({'Uniprot':entry.accNum(),'Pfam':string.join(rje.sortKeys(entry.dict['DB']['Pfam']),'|')})
                    pfamdb.saveToFile(pfamfile)
                self.printLog('#PFAM','%s %s Uniprot entries with Pfam domains.' % (rje.iStr(pfamdb.entryNum()),spec))
                self.printLog('#PFAM','%s different Pfam domains mapped to %s Uniprot.' % (rje.iLen(pfamdb.index('Pfam',splitchar='|')),spec))
                ## ~ [3c] Load DMI data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pingu.setStr({'DMIFile':self.getStr('DMIFile')})
                dmifile = pingu.sourceDataFile('DMIFile',expect=True,ask=True)
                dmidb = pingu.db().addTable(dmifile,mainkeys=["ELM identifier","Interaction Domain Id"],datakeys=["ELM identifier","Interaction Domain Id"],ignore=['#'],name='DMI',expect=False)
                if dmidb: dmidb.renameField("Interaction Domain Id",'Pfam')
                else: dmidb = pingu.db().addTable(dmifile,mainkeys=["ELM identifier",'Pfam'],datakeys=["ELM identifier",'Pfam'],ignore=[''],name='DMI',expect=False)
                if dmidb: dmidb.renameField("ELM identifier",'Motif')
                else: dmidb = pingu.db().addTable(dmifile,mainkeys=['Motif','Pfam'],datakeys=['Motif','Pfam'],ignore=[''],name='DMI',expect=True)
                if not dmidb: self.printLog('#FAIL','Cannot execute SLiMPPI run without DMIFile.'); return False
                dmidb.index('Motif')
                dmidb.index('Pfam')

            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <4> ### Generate Methods                                                                                        #
#########################################################################################################################
    def generate(self): ### Generates sequence data from mutations and wildtype sequences.
        '''
        Generates sequence data from mutations and wildtype sequences.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('MutFile')):
                self.errorLog('Cannot find MutFile "%s"!' % self.getStr('MutFile'),printerror=False)
                raise IOError
            db = self.db()
            seqlist = self.obj['SeqList']
            mfield = self.getStr('MutField')
            pfield = self.getStr('ProtField')
            if self.getStrLC('SplitField'): splitfield = self.getStr('SplitField')
            else: splitfield = None
            flankx = self.getInt('MutFlanks')
            ## ~ [1a] Load Mutation Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if splitfield:
                mdb = db.addTable(self.getStr('MutFile'),mainkeys=[self.getStr('ProtField'),self.getStr('MutField'),splitfield],lists=False,name='Mutations',expect=True,replace=True,ignore=[])
                mdb.dict['Index'][splitfield] = {}
            else: mdb = db.addTable(self.getStr('MutFile'),mainkeys=[self.getStr('ProtField'),self.getStr('MutField')],lists=False,name='Mutations',expect=True,replace=True,ignore=[])
            ## ~ [1b] Load Sequence Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist.loadSeq()
            if not seqlist.seqs(): raise ValueError
            seqdict = seqlist.makeSeqNameDic('accnum',clear=True,warnings=True)
            #self.debug(seqdict.keys()[:100])

            ### ~ [2] Process Mutation data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mdb.addField('SLiMMutantStatus'); mdb.addField('SLiMMutantAccNum',evalue=''); mdb.addField('SLiMMutantSequence',evalue='')
            mx = 0.0; mtot = mdb.entryNum()
            mindex = mdb.dict['Index']['SLiMMutantStatus'] = {}
            for status in ['NoMatch','NoSeq','WrongSeq','OK','Synonymous','Nonsense','NoStop','Unknown']: mindex[status] = []
            #nomatchx = 0; noseqx = 0; wrongseqx = 0
            flankseq = {}   # Dictionary of {'Uniprot':WTSeq with flanks}
            seqflanks = {}  # Dictionary of {'Uniprot':[(nflank,cflank)]}
            #for entry in mdb.entries():
            for mkey in mdb.dataKeys():
                entry = mdb.data(mkey)
                self.progLog('\r#MUT','Reading mutation data: %.1f%%' % (mx/mtot)); mx += 100
                if splitfield:
                    try: mdb.dict['Index'][splitfield][entry[splitfield]].append(mkey)
                    except: mdb.dict['Index'][splitfield][entry[splitfield]] = [mkey]
                mdata = rje.matchExp('^(\D)(\d+)(\D)$',entry[self.getStr('MutField')])
                if not mdata: mdata = rje.matchExp('^p.(\D)(\d+)(\D)$',entry[self.getStr('MutField')])
                if not mdata:
                    entry['SLiMMutantStatus'] = 'NoMatch'
                    mindex[entry['SLiMMutantStatus']].append(mkey)
                    #nomatchx += 1; mdb.data().pop(mkey)
                    continue
                if mdata[0] == mdata[2]: entry['SLiMMutantStatus'] = 'Synonymous'; mindex[entry['SLiMMutantStatus']].append(mkey); continue
                if '?' in mdata: entry['SLiMMutantStatus'] = 'Unknown'; mindex[entry['SLiMMutantStatus']].append(mkey); continue
                if mdata[0] == '*': entry['SLiMMutantStatus'] = 'NoStop'; mindex[entry['SLiMMutantStatus']].append(mkey); continue
                if mdata[2] == '*': entry['SLiMMutantStatus'] = 'Nonsense'; mindex[entry['SLiMMutantStatus']].append(mkey); continue
                if entry[self.getStr('ProtField')] not in seqdict: entry['SLiMMutantStatus'] = 'NoSeq'; mindex[entry['SLiMMutantStatus']].append(mkey); continue
                seq = seqlist.getSeq(seqdict[entry[self.getStr('ProtField')]])
                if not seq:
                    entry['SLiMMutantStatus'] = 'NoSeq'
                    mindex[entry['SLiMMutantStatus']].append(mkey)
                    #noseqx += 1; mdb.data().pop(mkey)
                    continue
                pos = int(mdata[1])
                if len(seq[1]) < pos or seq[1][pos-1] != mdata[0]:
                    if self.v() > 1:
                        if len(seq[1]) < pos: self.warnLog('Sequence too short! %s vs %s aa' % (entry[self.getStr('MutField')],len(seq[1])),'wrongseq',suppress=True)
                        else: self.warnLog('Sequence position mismatch %s = %s (%s)' % (entry[self.getStr('MutField')],seq[1][pos-1],seq[1][pos-3:][:5]),'wrongseq',suppress=True)
                    entry['SLiMMutantStatus'] = 'WrongSeq'
                    mindex[entry['SLiMMutantStatus']].append(mkey)
                    #wrongseqx += 1; mdb.data().pop(mkey)
                    continue
                entry['SLiMMutantStatus'] = 'OK'
                mindex[entry['SLiMMutantStatus']].append(mkey)
                mseq = entry['SLiMMutantSequence'] = seq[1][:pos-1] + mdata[2] + seq[1][pos:]
                if flankx > 0:
                    nflank = max(0,pos-1-flankx)
                    cflank = pos + flankx
                    entry['SLiMMutantSequence'] = mseq.lower()[:nflank] + mseq[nflank:cflank] + mseq.lower()[cflank:]
                    if entry[pfield] not in seqflanks: seqflanks[entry[pfield]] = []
                    seqflanks[entry[pfield]].append((nflank,cflank))
                    #!# This could also include the code to have all mutations in one sequence #!#
                    #self.debug(entry['SLiMMutantSequence'])
                    if entry[pfield] not in flankseq: flankseq[entry[pfield]] = seq[1].lower()
                    mseq = flankseq[entry[pfield]]
                    flankseq[entry[pfield]] = mseq[:nflank] + mseq.upper()[nflank:cflank] + mseq[cflank:]
                    #self.deBug(flankseq[entry[pfield]])
                if len(seq[1]) != len(entry['SLiMMutantSequence']): raise ValueError
                if mdata[2] == '*': entry['SLiMMutantSequence'] = entry['SLiMMutantSequence'][:pos-1]
            self.printLog('\r#MUT','Read mutation data for %s entries.' % rje.iStr(mtot))
            for status in rje.sortKeys(mdb.index('SLiMMutantStatus')):
                self.printLog('#MUT','%s: %s entries.' % (status,rje.iLen(mdb.index('SLiMMutantStatus')[status])))
            ox = 0.0; otot = len(mindex['OK'])
            if flankx > 0:  #?# Why is this needed in addition to above? (What does above do?!)
                for entry in mdb.indexEntries('SLiMMutantStatus','OK'):
                    self.progLog('\r#FLANK','Generating MutFlanks (%d) sequences: %.2f%%' % (flankx,ox/otot)); ox += 100.0
                    mseq = entry['SLiMMutantSequence']
                    for (nflank,cflank) in seqflanks[entry[pfield]]: mseq = mseq[:nflank] + mseq.upper()[nflank:cflank] + mseq[cflank:]
                    entry['SLiMMutantSequence'] = mseq
                self.printLog('\r#FLANK','Generated MutFlanks (%d) sequences for %s mutations.' % (flankx,rje.iStr(otot)))
            #self.printLog('#MUT','NoMatch: %s entries.' % (rje.iStr(nomatchx)))
            #self.printLog('#MUT','NoSeq: %s entries.' % (rje.iStr(noseqx)))
            #self.printLog('#MUT','WrongSeq: %s entries.' % (rje.iStr(wrongseqx)))
            #self.printLog('#MUT','OK: %s entries.' % (rje.iStr(mdb.entryNum())))
            ox = 0.0; otot = len(mindex['OK'])
            for entry in mdb.indexEntries('SLiMMutantStatus','OK'):
                self.progLog('\r#ACC','Creating mutant Accnum: %.1f%%' % (ox/otot)); ox += 100.0
                entry['SLiMMutantAccNum'] = '%s.%s' % (entry[self.getStr('ProtField')],entry[self.getStr('MutField')])
            self.printLog('\r#ACC','Created mutant Accnum for %s mutations.' % rje.iStr(otot))

            ### ~ [3] Output sequence data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.mkDir(self,self.basefile())
            ## ~ [3a] Split data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if splitfield:
                self.printLog('#SPLIT','%s different %s values.' % (rje.iLen(mdb.index(splitfield)),splitfield))
                for split in mdb.indexKeys(splitfield):
                    mutkeys = rje.listIntersect(mdb.index('SLiMMutantStatus')['OK'],mdb.index(splitfield)[split])
                    mutkeys.sort()
                    if len(mutkeys) < self.getInt('MinMutant'):
                        self.printLog('#MIN','%s %s mutants < %s (minmutant=X)' % (rje.iLen(mutkeys),split,rje.iStr(self.getInt('MinMutant'))))
                        continue
                    elif len(mutkeys) > self.getInt('MaxMutant'):
                        self.printLog('#MAX','%s %s mutants < %s (maxmutant=X)' % (rje.iLen(mutkeys),split,rje.iStr(self.getInt('MaxMutant'))))
                        continue
                    ## Save reduced mutant table ##
                    if self.dev(): mdb.dropField('SLiMMutantSequence')
                    mdb.saveToFile('%s.%s.tdt' % (self.basefile(),split),savekeys=mdb.index(splitfield)[split])     #!# Don't want to save Mutant Sequence! (Probably WRONG)
                    ## Save wildtype sequences ##
                    wtacc = rje.listIntersect(mdb.indexDataList('SLiMMutantStatus','OK',pfield),mdb.indexDataList(splitfield,split,pfield))
                    wtacc.sort()
                    wtfile = '%s.%s.wildtype.fas' % (self.basefile(),split); rje.backup(self,wtfile,appendable=False)
                    mutfile = '%s.%s.mutant.fas' % (self.basefile(),split); rje.backup(self,mutfile,appendable=False)
                    WT = open(wtfile,'w'); ox = 0.0
                    for acc in wtacc:
                        self.progLog('\r#WT','Saving wildtype sequences: %.1f%%' % (ox/otot)); ox += 100.0
                        seq = seqlist.getSeq(seqdict[acc])
                        if flankx > 0: WT.write('>%s\n%s\n' % (seq[0],flankseq[acc]))
                        else: WT.write('>%s\n%s\n' % seq)
                    self.printLog('\r#WT','Saved %s wildtype %s sequences to %s' % (rje.iLen(wtacc),split,wtfile))
                    WT.close()
                    if self.dev():
                        self.printLog('#DEV','No mutant sequence output in Dev mode.')
                        continue
                    ## Save mutant sequences ##
                    MUT = open(mutfile,'w')
                    muts = []; ox = 0.0
                    for mkey in mutkeys:
                        self.progLog('\r#MUT','Saving mutant sequences: %.2f%%' % (ox/otot)); ox += 100.0
                        entry = mdb.data(mkey)
                        seq = seqlist.getSeq(seqdict[entry[pfield]])
                        name = string.split(seq[0])
                        mut = '%s.%s' % (name[0],entry[mfield])
                        if mut in muts: continue
                        else: muts.append(mut)
                        MUT.write('>%s.%s %s\n%s\n' % (name[0],entry[mfield],string.join(name[1:]),entry['SLiMMutantSequence']))
                    self.printLog('\r#MUT','Saved %s mutant %s sequences to %s' % (rje.iLen(muts),split,mutfile))
                    MUT.close()
            ## ~ [3b] Full data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                wtfile = '%s.wildtype.fas' % self.basefile(); rje.backup(self,wtfile,appendable=False)
                mutfile = '%s.mutant.fas' % self.basefile(); rje.backup(self,mutfile,appendable=False)
                WT = open(wtfile,'w'); ox = 0.0
                wtacc = mdb.indexDataList('SLiMMutantStatus','OK',pfield)
                wtacc.sort()
                for acc in wtacc:
                    self.progLog('\r#WT','Saving wildtype sequences: %.1f%%' % (ox/otot)); ox += 100.0
                    seq = seqlist.getSeq(seqdict[acc])
                    if flankx > 0: WT.write('>%s\n%s\n' % (seq[0],flankseq[acc]))
                    else: WT.write('>%s\n%s\n' % seq)
                self.printLog('\r#WT','Saved wildtype sequences to %s' % (wtfile))
                WT.close()
                if self.dev(): mdb.dropField('SLiMMutantSequence')
                mdb.saveToFile()
                if self.dev(): return self.printLog('#DEV','No mutant sequence output in Dev mode.')
                MUT = open(mutfile,'w')
                muts = []; ox = 0.0
                mutkeys = mdb.index('SLiMMutantStatus')['OK']
                mutkeys.sort()
                for mkey in mutkeys:
                    self.progLog('\r#MUT','Saving mutant sequences: %.2f%%' % (ox/otot)); ox += 100.0
                    entry = mdb.data(mkey)
                    seq = seqlist.getSeq(seqdict[entry[pfield]])
                    name = string.split(seq[0])
                    mut = '%s.%s' % (name[0],entry[mfield])
                    if mut in muts: continue
                    else: muts.append(mut)
                    MUT.write('>%s.%s %s\n%s\n' % (name[0],entry[mfield],string.join(name[1:]),entry['SLiMMutantSequence']))
                self.printLog('\r#MUT','Saved mutant sequences to %s' % (mutfile))
                MUT.close()

            return True
        except: self.errorLog('%s.generate error' % self.prog()); return False
#########################################################################################################################
    ### <5> ### SLiMProb Methods                                                                                        #
#########################################################################################################################
    def slimProb(self): ### Runs SLiMProb on *.wildtype.fas and *.mutant.fas using the provided INI file list.
        '''
        Runs SLiMProb on *.wildtype.fas and *.mutant.fas using the provided INI file list.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dev(): return self.slimSearch()
            if self.getStrLC('RunID'): runid = self.getStr('RunID')
            else: runid = 'SLiMMutant'
            inbase = rje.baseFile(self.getStr('MutFile'))
            basefile = self.basefile()
            batch = 'batch='
            for stype in ['wildtype','mutant']:
                sfile = '%s.%s.fas' % (basefile,stype)
                if not rje.exists(sfile): sfile = '%s.%s.fas' % (inbase,stype)
                if not rje.exists(sfile):
                    self.errorLog('Cannot find "%s". Check basefile=X setting.' % sfile,printerror=False)
                    raise IOError
                batch += sfile + ','
            scmd = ['seqin=None',batch[:-1]]
            bcmd = ['slimmutant=T','efilter=F','maxsize=0','walltime=0','extras=1','runid=%s' % runid,'resfile=%s.slimmutant.occ.csv' % basefile]
            if self.getInt('MutFlanks') > 0: bcmd.append('casemask=Lower')

            ### ~ [2] Run SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['SlimINI']:
                append = False
                for ini_file in self.list['SlimINI']:
                    icmd = rje.longCmd(rje.iniCmds('',ini_file,iowarning=True,altpaths=False))
                    slimprob.SLiMProb(self.log,bcmd+self.cmd_list+icmd+scmd).run()
                    if not append: scmd.append('append=T'); append = True
            else:
                slimprob.SLiMProb(self.log,bcmd+self.cmd_list+scmd).run()

            return True
        except: self.errorLog('%s.slimProb error' % self.prog()); raise
#########################################################################################################################
    def slimSearch(self):   ### Replacement SLiMProb method that directly performs search and filter. (No stats.)
        '''Replacement SLiMProb method that directly performs search and filter. (No stats.).'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # NOTE: This is dev only for the moment. Will remove the need for flank analysis. (Why bother?)
            # NOTE: Mutations affecting the disorder status of a region should be checked independently from ELM +/-
            basefile = self.basefile()
            if not self.force() and rje.exists('%s.MutOcc.csv' % basefile) and (self.i() < 1 or rje.yesNo('MutOcc file found. Use this file and skip SLiMSearch?')): return True
            db = self.db()
            ## ~ [1a] Load Mutation Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mdb = self.db('Mutations',add=False)
            if not mdb: mdb = db.addTable(self.getStr('MutFile'),mainkeys=[self.getStr('ProtField'),self.getStr('MutField')],lists=False,name='Mutations',expect=True,replace=True,ignore=[])
            #mdb.dropEntriesDirect('SLiMMutantStatus',['OK'],inverse=True)   # Reduce to actual mutations
            for mkey in mdb.dataKeys():
                if mdb.data(mkey)['SLiMMutantStatus'] != 'OK': mdb.data().pop(mkey)
            self.printLog('#OK','Entries reduced to %s "OK" mutations.' % rje.iStr(mdb.entryNum()))
            mdb.index(self.getStr('ProtField'))
            ## ~ [1b] Load SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.slims(): self.obj['SlimList'].loadMotifs()
            slimdict = self.obj['SlimList'].slimDict(corelist=True)     # {corename:[SLiM Objects]}
            ## ~ [1c] SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            seqlist = self.obj['SeqList'] = rje_seqlist.SeqList(self.log,seqcmd)
            seqdict = seqlist.makeSeqNameDic('max')
            ## ~ [1d] Setup MutOcc table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ofields = string.split('RunID,Motif,Start_Pos,End_Pos,Pattern,Desc,UPC,Uniprot,Mut,Mutants,WT,Mutations',',')  # *.MutOcc.csv
            odb = self.db().addEmptyTable('MutOcc',ofields,['RunID','Motif','Uniprot','Start_Pos','End_Pos'])
            odb.setStr({'Delimit':','})
            ## ~ [1e] Setup MutExp table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # eGOF and eLOF are expected GOF and LOF given that number of mutations in THAT protein
            # Uses total RunID mutation frequency to establish expected numbers of mutations
            #!# Consider adding additional output using overall mutation numbers?
            efields = string.split('RunID,Motif,Uniprot,eGOF,eLOF,Mutations',',')  # *.MutOcc.csv
            odb = self.db().addEmptyTable('MutOcc',ofields,['RunID','Motif','Uniprot','Start_Pos','End_Pos'])
            odb.setStr({'Delimit':','})

            ### ~ [2] Run SLiMSearch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['SlimINI']: self.list['SlimINI'] = ['']
            for ini_file in self.list['SlimINI']:
                ## ~ [2a] Setup SLiMCore object for masking and RunID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if ini_file: icmd = rje.longCmd(rje.iniCmds('',ini_file,iowarning=True,altpaths=False))
                else:
                    runid = 'SLiMMutant'
                    if self.getStrLC('RunID'): runid = self.getStr('RunID')
                    icmd = ['runid=%s' % runid]
                mcmd = ['metmask=F','compmask=None','posmask=']
                slimcore = rje_slimcore.SLiMCore(self.log,mcmd+self.cmd_list+icmd)    # Loads masking and RunID options. Can use MegaSLiM.
                if ini_file:
                    if not slimcore.getStrLC('RunID'): slimcore.setStr({'RunID':ini_file})
                    self.printLog('#MASK','%s (%s): %s' % (ini_file,slimcore.getStr('RunID'),slimcore.maskText()))
                else: self.printLog('#MASK','%s: %s' % (slimcore.getStr('RunID'),slimcore.maskText()))
                runid = slimcore.getStr('RunID')
                ## ~ [2b] Mask sequences & mutation dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                maskedseq = {}
                mutations = {}  #!# Add option to load? #!#
                sx = 0.0; stot = len(mdb.indexKeys(self.getStr('ProtField')))
                for prot in mdb.indexKeys(self.getStr('ProtField')):
                    self.progLog('\r#MASK','%s Masking %s input sequences: %.2f%%' % (slimcore.maskText(),runid,sx/stot)); sx += 100.0
                    (name,sequence) = seqlist.getSeq(seqdict[prot])
                    sname = string.split(name)[0]
                    maskseq = slimcore.maskSequence(sequence,name,log=self.getBool('LogMask'),screen=False)
                    maskedseq[prot] = maskseq
                    mutations[prot] = []
                    for mutant in mdb.indexDataList(self.getStr('ProtField'),prot,self.getStr('MutField')):
                        (wt,pos,mut) = rje.matchExp('(\D)(\d+)(\D)',mutant)
                        pos = int(pos) - 1
                        if maskseq[pos] == 'X': continue    # Affects masked sequence: ignore.
                        mutations[prot].append(mutant)
                    mutations[prot]  = string.replace(string.join(mutations[prot],';'),'p.','')
                self.printLog('\r#MASK','%s Masked %s %s input sequences.' % (slimcore.maskText(),rje.iStr(stot),runid))
                self.printLog('#MUT','Generated mutation lists for %s proteins.' % rje.iLen(mutations))
                self.mutFreqDict(mutations)
                ## ~ [2c] Perform searches and populate odb and edb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# Cannot use slim.misMatch({1:1} as position needs to be inverted.


                sx = 0.0; stot = len(self.slims()) * len(mdb.indexKeys(self.getStr('ProtField')))
                for prot in mdb.indexKeys(self.getStr('ProtField')):
                    (name,sequence) = seqlist.getSeq(seqdict[prot])
                    sname = string.split(name)[0]
                    maskseq = slimcore.maskSequence(sequence,name,log=self.getBool('LogMask'),screen=False)

                    #!# Strip away the motif extra _1a etc. before combining.
                    #!# Do at the SLiMSearch stage! Will need to make a dictionary of motifs to keep analysis together?
                    #!# Make a list (or dict) of lists of SLiM objects. Add as rje_slimlist method.
                    # self.obj['SlimList'].slimDict()   {core:[slims]}

                    for slim in self.slims():
                        self.progLog('\r#SEARCH','Searching SLiMs: %.2f%%' % (sx/stot)); sx += 100.0
                        # Scan wildtype
                        wtocc = []
                        for occ in slim.searchSequence(sequence=maskseq): wtocc.append(occ['Pos'])    # [{Pos,Variant,Match,ID,MisMatch}]
                        wtocc = rje.sortUnique(wtocc)
                        mutocc = {}
                        # Scan mutants
                        mutations = []
                        for mutant in mdb.indexDataList(self.getStr('ProtField'),prot,self.getStr('MutField')):
                            (wt,pos,mut) = rje.matchExp('(\D)(\d+)(\D)',mutant)
                            pos = int(pos) - 1
                            if maskseq[pos] == 'X': # Affects masked sequence: ignore.
                                continue    # Exclude masked mutations from counts.
                                  #for occpos in wtocc[0:]:     # Will not affect wt
                                #    if occpos not in mutocc: mutocc[occpos] = []
                                #    if mutant not in mutocc[occpos]: mutocc[occpos].append(mutant)
                            elif maskseq[pos] != wt: raise ValueError('"OK" mutation has sequence mismatch %s %d: wt %s vs prot %s' % (prot,pos+1,wt,maskseq[pos]))
                            else:
                                for occ in slim.searchSequence(sequence=rje.strSub(maskseq,pos,pos,mut)):
                                    occpos = occ['Pos']
                                    if occpos not in mutocc: mutocc[occpos] = []
                                    if mutant not in mutocc[occpos]: mutocc[occpos].append(mutant)
                            mutations.append(string.split(mutant,'.')[-1])  # All unmasked mutations
                        if not mutations: continue  # No unmasked mutations. Do not waste processing.
                        # Now have list of wtpos and a mutpos dictionary
                        allpos = rje.sortUnique(wtocc+mutocc.keys())   # All unique positions to check for GOF/LOF
                        if not allpos: continue     # No occurrences. Do not waste processing.
                        entrycore = {'RunID':slimcore.getStr('RunID'),'Motif':slim.name(),'Pattern':slim.pattern(),
                                     'Desc':string.join(string.split(name)[1:]),'Uniprot':sname,
                                     'End_Pos':0,'UPC':0,   #!# Drop these!
                                     'Mutations':string.join(mutations,';')}
                        for occpos in allpos:
                            if occpos not in mutocc: mutocc[occpos] = []
                            entry = rje.combineDict({'Start_Pos':occpos,'Mut':len(mutocc[occpos]),'Mutants':string.replace(string.join(mutocc[occpos],';'),'p.',''),'WT':len(mutations)},entrycore)
                            if occpos not in wtocc: entry['WT'] = 0
                            if entry['WT'] != entry['Mut']: odb.addEntry(entry)
                self.printLog('\r#SEARCH','Searching for SLiMs complete.')
            odb.saveToFile('%s.MutOcc.csv' % basefile)
            return True
        except:
            self.errorLog('%s.slimSearch error' % self.prog())
            if self.dev(): raise
            return False
#########################################################################################################################
    ### <6> ### Analyse Methods                                                                                         #
#########################################################################################################################
    def analyse(self):  ### Analysis of SLiMProb results
        '''
        Analysis of SLiMProb results.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            inbase = rje.baseFile(self.getStr('MutFile'))
            wtbase = wtfile = '%s.wildtype' % self.basefile()
            if not rje.exists(wtfile): wtfile = '%s.wildtype' % inbase
            mutbase = mutfile = '%s.mutant' % self.basefile()
            if not rje.exists(mutfile): mutfile = '%s.mutant' % inbase
            if not self.getStrLC('ResFile'): resfile = '%s.slimmutant.occ.csv' % self.basefile()
            else:
                resfile = self.getStrLC('ResFile')
                rsplit = string.split(resfile,'.')
                if rsplit[-2] != 'occ': resfile = string.join(rsplit[:-1]+['occ']+rsplit[-1:],'.')
            if not rje.exists(resfile) and not self.dev():
                self.errorLog('Cannot find SLiMProb results file "%s"!' % resfile,printerror=False)
                raise IOError
            ## ~ [1a] Analysed Mutations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Replace with parsing from original MutFile? Make MutFreq dictionary at generate/SLiMSearch stage? #!#
            mutations = {}  # Dictionary of AccNum:Mutations
            self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s.fas' % wtfile])
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqin=%s.fas' % mutfile])
            for seq in seqlist.seqs():
                name = seqlist.shortName(seq)
                nsplit = string.split(name,'.')
                if nsplit[-2] == 'p': prot = string.join(nsplit[:-2],'.')
                else: prot = string.join(nsplit[:-1],'.')
                if prot not in mutations: mutations[prot] = []
                mutations[prot].append(nsplit[-1])
            for prot in mutations:
                mutations[prot].sort()
                mutations[prot] = string.join(mutations[prot],';')
            self.printLog('#MUT','Generated mutation lists for %s proteins.' % rje.iLen(mutations))

            ### ~ [2] Load Occ Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mutoccfile = '%s.MutOcc.csv' % db.basefile()
            odb = None
            if self.force() or not rje.exists(mutoccfile):
                odb = db.addTable(resfile,mainkeys=['Dataset','RunID','Motif','Seq','Start_Pos','End_Pos'],name='MutOcc')
                for entry in odb.entries():
                    if entry['Dataset'].startswith(mutbase): entry['Dataset'] = 'Mut'
                    elif entry['Dataset'].startswith(wtbase): entry['Dataset'] = 'WT'
                    else: self.debug(entry)
                ox = 0
                for entry in odb.entries():
                    if entry['Dataset'] == 'WT' and entry['Seq'] not in mutations: entry['Dataset'] = 'NULL'; ox += 1
                odb.remakeKeys()
                odb.dropEntriesDirect('Dataset',['Mut','WT'],inverse=True)
                if ox: self.printLog('#WTOCC','Removed %s WT Occ for proteins without mutations.' % rje.iStr(ox))
                self.printLog('#WTOCC','%s wildtype SLiMProb Occurrences.' % rje.iLen(odb.indexEntries('Dataset','WT')))
                self.printLog('#MUTOCC','%s mutant SLiMProb Occurrences.' % rje.iLen(odb.indexEntries('Dataset','Mut')))
                if self.getStrLC('RunID'): odb.dropEntriesDirect('RunID',[self.getStr('RunID')],inverse=True)
                if self.list['MotList']: odb.dropEntriesDirect('Motif',[self.list['MotList']],inverse=True)
                # Dataset,RunID,Masking,Motif,Seq,Start_Pos,End_Pos,Prot_Len,Pattern,Match,Variant,MisMatch,Desc
                odb.dropFields(['Match','Variant','MisMatch','Masking','Prot_Len'])
                odb.makeField('#Seq#','Uniprot'); odb.addField('Mutation',evalue='')
                for entry in odb.indexEntries('Dataset','Mut'):
                    msplit = string.split(entry['Seq'],'.')
                    if msplit[-2] == 'p': entry['Uniprot'] = string.join(msplit[:-2],'.')
                    else: entry['Uniprot'] = string.join(msplit[:-1],'.')
                    entry['Mutation'] = msplit[-1]

            ### ~ [3] Reshape Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                odb.addField('Occ',evalue=1)
                odb.compress(['Dataset','RunID','Motif','Uniprot','Start_Pos','End_Pos'],rules={'Mutation':'list','Seq':'str','Occ':'sum'},default='mean',best=[])
                odb.dropField('Seq')
                odb.reshapeWide('Dataset',reshape=['Occ','Mutation'],evalue=0)
                odb.renameField('Occ|Mut','Mut'); odb.renameField('Occ|WT','WT')
                odb.renameField('Mutation|WT','Mutations'); odb.renameField('Mutation|Mut','Mutants')
                #self.debug(odb.fields())
                if self.dev(): odb.saveToFile('%s.temp.csv' % self.basefile())
                prex = odb.entryNum()
                for entry in odb.entries()[0:]:
                    #self.debug('%s' % entry)
                    if not entry['Mutants']: entry['Mutants'] = ''
                    entry['Mutations'] = mutations[entry['Uniprot']]
                    if entry['WT'] == 1 and entry['Mut'] == len(string.split(entry['Mutations'],';')): odb.dropEntry(entry)
                    if entry['WT'] == 1: entry['WT'] = len(string.split(entry['Mutations'],';'))
                self.printLog('#DROP','Purged unaffected occurrences: %s -> %s entries.' % (rje.iStr(prex),rje.iStr(odb.entryNum())))
                odb.dropField('Unaffected')
                odb.saveToFile()

            mutvwtfile = '%s.Mut_vs_WT.csv' % db.basefile()
            if self.force() or not rje.exists(mutvwtfile):
                if not odb: odb = self.db('MutOcc',add=False)
                if not odb:
                    odb = db.addTable(mutoccfile,mainkeys=['RunID','Motif','Uniprot','Start_Pos','End_Pos'],name='MutOcc')
                    odb.dataFormat({'Mut':'int','WT':'int'})
                for entry in odb.entries():
                    if entry['WT'] > entry['Mut']: entry['WT'] -= entry['Mut']; entry['Mut'] -= entry['Mut']
                    else: entry['WT'] -= entry['WT']; entry['Mut'] -= entry['WT']

            ### ~ [4] Perform analyses ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                odb.dropFields(['Desc','Mutation'])
                odb.compress(['RunID','Motif'],default='sum')
                odb.keepFields(['RunID','Motif','Pattern','Mut','WT'])
                odb.setStr({'Name':'Mut_vs_WT'})
                #!# Add some calculations based on mutation-frequency-derived expectations #!#
                mutfreq = self.mutFreqDict(mutations)
                odb.addField('fMut')
                for entry in odb.entries(): entry['fMut'] = (float(entry['Mut']) / (entry['Mut'] + entry['WT']))
                odb.addField('eMut')
                emut = {}
                for pattern in odb.index('Pattern'):
                    pmut = 0.0; px = 0
                    for el in rje_slim.prestoFromCode(rje_slim.slimFromPattern(pattern,trimx=True)):
                        if el == 'X': continue
                        if el not in emut:
                            efrom = 0.0
                            eto = 0.0
                            for xfrom in mutfreq:
                                if xfrom == 'X': continue
                                for xto in mutfreq[xfrom]:
                                    if xto == 'X': continue
                                    elif xfrom in el and xto in el: continue
                                    elif xfrom in el: efrom += mutfreq[xfrom][xto]
                                    elif xto in el: eto += mutfreq[xfrom][xto]
                            if not efrom and not eto: continue
                            emut[el] = (eto / (efrom + eto))    # Proportion of mutations that go TO SLiM
                        pmut += emut[el]; px += 1
                    for entry in odb.indexEntries('Pattern',pattern): entry['eMut'] = pmut / px
                odb.addField('pMut'); odb.addField('pWT')
                #!# Add proglog
                ox = 0.0; otot = odb.entryNum()
                for entry in odb.entries():
                    self.progLog('\r#PROB','Calculating pMut and pWT enrichment: %.2f%%' % (ox/otot)); ox += 100.0
                    entry['pMut'] = rje.poisson(entry['Mut'],(entry['Mut']+entry['WT'])*entry['eMut'],exact=False,callobj=self,uselog=True)
                    entry['pWT'] = rje.poisson(entry['WT'],(entry['Mut']+entry['WT'])*(1-entry['eMut']),exact=False,callobj=self,uselog=True)
                self.printLog('\r#PROB','Calculated pMut and pWT enrichment (Poisson from MutFreq)')
                odb.saveToFile()
            else: mutfreq = self.mutFreqDict(mutations,save=False)

            ### ~ [5] SLiMPPI Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            odmifile = '%s.Mut_vs_WTDMI.csv' % db.basefile()
            #self.debug(odmifile)
            if rje.exists(odmifile) and not self.force(): return True
            maxlinks = 3
            if not self.getBool('SLiMPPI'): self.printLog('#PPI','No SLiMPPI analysis (slimppi=F).'); return True
            pingu = self.obj['PINGU']
            ppidb = pingu.db('PPISource')
            pfamdb = pingu.db('Pfam')
            dmidb = pingu.db('DMI')
            occdb = self.db('MutOcc',add=False)
            if not occdb: occdb = self.db().addTable(name='MutOcc',mainkeys=['RunID','Motif','Uniprot','Start_Pos','End_Pos'])   # RunID,Motif,Start_Pos,End_Pos,Pattern,Desc,UPC,Uniprot,Mut,Mutants,WT,Mutations
            occdb.setStr({'Name':'MutOccDMI'})
            if 'DomLink' not in occdb.fields():
                occdb.addField('DomLink',evalue=0)
                # Motif field should match Motif in other tables
                # Uniprot is full name and should be reduced to accnum
                for entry in occdb.entries(): entry['Uniprot'] = string.split(entry['Uniprot'],'__')[-1]
                ## ~ [5a] Merge mutation gain/loss table with ELM-Pfam-PPI data ~~~~~~~~~~~~~~~~~~~~~~~ ##
                motlist = occdb.indexKeys('Motif'); mx = 0.0; mtot = len(motlist); dx = mtot
                for motif in motlist[0:]:
                    self.progLog('#LINK','Making Motif-Domain-Protein links: %.1f%%' % (mx/mtot)); mx += 100.0
                    if motif not in dmidb.index('Motif'): dx -= 1; continue
                    domlist = dmidb.indexDataList('Motif',motif,'Pfam')
                    # Start with all proteins containing relevant Pfam domains
                    domlinks = {0:[]}
                    for pfam in domlist: domlinks[0] += pfamdb.indexDataList('Pfam',pfam,'Uniprot')
                    domlinks[0] = rje.sortUnique(domlinks[0])
                    # Expand one PPI degree at a time
                    # Keep expanding until target protein hit or no more proteins added. (Or stop at 3?)
                    for n in range(maxlinks):
                        domlinks[n+1] = []
                        for entry in ppidb.entries():
                            if entry['HubUni'] not in domlinks[n]: continue
                            if n and entry['SpokeUni'] in domlinks[n-1]: continue
                            if entry['SpokeUni'] in domlinks[n+1]: continue
                            domlinks[n+1].append(entry['SpokeUni'])
                        domlinks[n+1].sort()
                    for entry in occdb.indexEntries('Motif',motif):
                        entry['DomLink'] = -1
                        for n in range(maxlinks):
                            if entry['Uniprot'] in domlinks[n+1]: entry['DomLink'] = n+1; break
                self.printLog('#LINK','Made Motif-Domain-Protein links for %s of %s motifs.' % (rje.iStr(dx),rje.iStr(mtot)))
            ## ~ [5b] Repeat above analysis by proximity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            odb = occdb
            odb.dataFormat({'WT':'int','Mut':'int'})
            for entry in odb.entries():
                if entry['WT'] > entry['Mut']: entry['WT'] -= entry['Mut']; entry['Mut'] -= entry['Mut']
                else: entry['WT'] -= entry['WT']; entry['Mut'] -= entry['WT']
            odb.dropFields(['Desc','Mutation']);
            #self.debug(odb.fields())
            odb.compress(['RunID','Motif','DomLink'],default='sum')
            odb.keepFields(['RunID','Motif','DomLink','Pattern','Mut','WT'])
            odb.setStr({'Name':'Mut_vs_WTDMI'})
            #!# Add some calculations based on mutation-frequency-derived expectations #!#
            odb.addField('fMut')
            for entry in odb.entries(): entry['fMut'] = (float(entry['Mut']) / (entry['Mut'] + entry['WT']))
            odb.addField('eMut')
            emut = {}
            for pattern in odb.index('Pattern'):
                pmut = 0.0; px = 0
                #self.debug(rje_slim.prestoFromCode(rje_slim.slimFromPattern(pattern,trimx=True)))
                for el in rje_slim.prestoFromCode(rje_slim.slimFromPattern(pattern,trimx=True)):
                    if el == 'X': continue
                    if el not in emut:
                        efrom = 0.0
                        eto = 0.0
                        for xfrom in mutfreq:
                            if xfrom == 'X': continue
                            for xto in mutfreq[xfrom]:
                                if xto == 'X': continue
                                elif xfrom in el and xto in el: continue
                                elif xfrom in el: efrom += mutfreq[xfrom][xto]
                                elif xto in el: eto += mutfreq[xfrom][xto]
                        if not efrom and not eto: continue
                        emut[el] = (eto / (efrom + eto))    # Proportion of mutations that go TO SLiM
                    pmut += emut[el]; px += 1
                for entry in odb.indexEntries('Pattern',pattern): entry['eMut'] = pmut / px
            odb.addField('pMut'); odb.addField('pWT')
            ox = 0.0; otot = odb.entryNum()
            for entry in odb.entries():
                self.progLog('\r#PROB','Calculating pMut and pWT enrichment: %.2f%%' % (ox/otot)); ox += 100.0
                entry['pMut'] = rje.poisson(entry['Mut'],(entry['Mut']+entry['WT'])*entry['eMut'],exact=False,callobj=self,uselog=True)
                entry['pWT'] = rje.poisson(entry['WT'],(entry['Mut']+entry['WT'])*(1-entry['eMut']),exact=False,callobj=self,uselog=True)
            self.printLog('\r#PROB','Calculated pMut and pWT enrichment (Poisson from MutFreq)')
            odb.saveToFile()
            return True
        except: self.errorLog('%s.analyse error' % self.prog()); return False
#########################################################################################################################
    def mutFreqDict(self,mutdict={},runid='SLiMMutant',save=True):  ### Returns optionally weighted mutation frequency dictionary
        '''
        Returns optionally weighted mutation frequency dictionary.
        >> mutdict:dict of {prot:[XnnnY]} aa subsitution mutations.
        >> runid:str ['SLiMMutant'] = RunID that identifies specific masking options etc.
        >> save:bool [True] = whether to save to file.
        << mutfreq: dict of {from:{to:X}} as raw numbers
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mutfreq = {'X':{'X':0}}
            mdb = self.db('MutFreq')
            if not mdb: mdb = self.db().addEmptyTable('MutFreq',['RunID','WT']+rje_seq.alph_protx,['RunID','WT'])
            ### ~ [1] Generate dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Generate from loaded data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not mutdict:
                if runid not in mdb.index('RunID'): raise ValueError('RunID "%s" not found in MutFreq table.' % runid)
                for entry in mdb.indexEntries('RunID',runid):
                    mutfreq[entry['WT']] = {}
                    for aa in rje_seq.alph_protx: mutfreq[entry['WT']][aa] = entry[aa]
                return mutfreq
            ## ~ [1b] Generate from given data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for prot in mutdict:
                if prot == '*':
                    entry = rje.combineDict({'RunID':runid,'WT':'*'},mutdict[prot])
                    mutfreq['*'] = mutdict[prot]    # This should contain total unmasked aa counts
                    mdb.addEntry(entry)
                    continue
                for mut in string.split(mutdict[prot],';'):
                    (xfrom, xto) = rje.matchExp('^(\D)\d+(\D)$',mut)
                    if xfrom not in mutfreq: mutfreq[xfrom] = {'X':0}
                    if xto not in mutfreq[xfrom]: mutfreq[xfrom][xto] = 0
                    if xto not in mutfreq['X']: mutfreq['X'][xto] = 0
                    for waa in ['X',xfrom]:
                        for maa in ['X',xto]: mutfreq[waa][maa] += 1
            ### ~ [2] Add to database table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for waa in rje_seq.alph_protx:
                mentry = {'WT':waa,'RunID':runid}
                for maa in rje_seq.alph_protx:
                    if waa in mutfreq and maa in mutfreq[waa]: mentry[maa] = mutfreq[waa][maa]
                    else: mentry[maa] = 0
                mdb.addEntry(mentry)
            if save: mdb.saveToFile()
            self.debug(mutfreq)
            return mutfreq
        except: self.errorLog('%s.mutFreqDict error' % self.prog()); return {}
#########################################################################################################################
### End of SECTION II: SLiMMutant Class                                                                                 #
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
    try: SLiMMutant(mainlog,cmd_list).run()

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
