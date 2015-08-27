#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       SLiMBench
Description:  Short Linear Motif prediction Benchmarking
Version:      1.9
Last Edit:    06/08/13
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMBench has two primary functions:

    1. Generating SLiM prediction benchmarking datasets from ELM (or other data in a similar format). This includes
    options for generating random and/or simulated datasets for ROC analysis etc.

    2. Assessing the results of SLiM predictions against a Benchmark. This program is designed to work with SLiMFinder
    and QSLiMFinder output, so some prior results parsing may be needed for other methods.

    Documentation for SLiMBench is currently under development. Please contact the author for more details.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sourcepath=PATH/    : Will look in this directory for input files if not found ['SourceData/']
    elmclass=FILE       : Download from ELM website of ELM classes ['elm_classes.tsv']
    elminstance=FILE    : Download from ELM website of ELM instances ['elm_instances.tsv']
    elmpfam=FILE        : Download from ELM website of ELM Pfam domain interactors ['elm_interaction_domains.tsv']
    uniprot=FILE        : File of downloaded UniProt entries (See rje_uniprot for more details) ['ELM.dat']

    ### ~ ELM BENCHMARK GENERATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genpath=PATH        : Output path for datasets generated with SLiMBench file generator [./SLiMBenchDatasets/]
    integrity=T/F       : Whether to quit by default if input integrity is breached [True]
    generate=T/F        : Whether to generate SLiMBench datasets from ELM input [False]
    slimmaker=T/F       : Whether to use SLiMMaker to "reduce" ELMs to more findable SLiMs [True]
    minupc=X            : Minimum number of UPC for ELM dataset [True]
    minic=X             : Min information content for a motif (1 fixed position = 1.0) [2.0]
    queries=T/F         : Whether to generate datasets with specific Query proteins [True]
    flankmask=LIST      : List of flanking mask options [none,win300,win100,flank5,site]
    searchini=LIST      : List of INI files containing search options (should have runid setting) []

    ### ~ ELM PPI/3DID BENCHMARK GENERATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    elmpfam=FILE        : Download from ELM website of ELM Pfam domain interactors ['elm_interaction_domains.tsv']
    pfamdata=FILE       : File mapping PFam domains onto genes/proteins (BioMart or HMM search) []
    xrefdata=FILE       : File of gene identifier cross-reference data from rje_genemap []
    3didsql=PATH        : Path to 3DID sql data. Use rje_mysql sqldump to extract 3DID DMI data. []
    dmidata=FILE        : File of 3DID DMI data ['3did.DMI.csv']
    pdbdata=FILE        : File mapping PDB identifiers onto genes/proteins []

    ### ~ RANDOM/SIMULATION BENCHMARK GENERATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    simulate=T/F        : Whether to generate simulated datasets using reduced ELMs (if found) [False]
    randomise=T/F       : Whether to generate randomised datasets (part of simulation if simulate=T) [False]
    randreps=X          : Number of replicates for each random (or simulated) datasets [10]
    simratios=LIST      : List of simulated ELM:Random rations [1,4,9,19]
    simcount=LIST       : Number of "TPs" to have in dataset [5,10]
    randir=PATH         : Output path for creation of randomised datasets [./SLiMBenchDatasets/Random/]
    randbase=X          : Base for random dataset name if simulate=F [ran]
    randsource=FILE     : Source for new sequences for random datasets [None]
    masking=T/F         : Whether to use SLiMCore masking for query selection [True]

    ### ~ BENCHMARK ASSESSMENT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    benchmark=T/F       : Whether to perfrom SLiMBench benchmarking assessment against motif file [False]
    datatype=X          : Type of data to be generated and/or benchmarked (elm/sim/simonly) [elm]
    resfiles=LIST       : List of (Q)SLiMFinder results files to use for benchmarking [*.csv]
    compdb=FILE         : Motif file to be used for benchmarking (default = reduced elmclass file) []
    benchbase=X         : Basefile for SLiMBench benchmarking output [slimbench]
    runid=LIST          : List of factors to split RunID column into (on '.') ['Program','Analysis']
    bycloud=X           : Whether to compress results into clouds prior to assessment (True/False/Both) [Both]
    sigcut=LIST         : Significance thresholds to use for assessment [0.05,0.01,0.001,0.0001]
    iccut=LIST          : Minimum IC for (Q)SLiMFinder results for benchmark assessment [2.0,2.1,3.0]
    slimlencut=LIST     : List of individual SLiM lengths to return results for (0=All) [0,3,4,5]
    noamb=T/F           : Filter out ambiguous patterns [False]
    # Add CompariMotif settings here for OT/TP etc.

    ### ~ GENERAL OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    force=T/F           : Whether to force regeneration of outputs (True) or assume existing outputs are right [False]
    backups=T/F         : Whether to (prompt if interactive and) generate backups before overwriting files [True]

See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seq, rje_seqlist, rje_slim, rje_slimcore, rje_slimlist, rje_uniprot, rje_zen
import comparimotif_V3 as comparimotif
import slimmaker, slimprob  #, slimsearch
# import rje_xref - need this for realistic test data?
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Functional version with benchmarking dataset generation.
    # 1.0 - Consolidation of "working" version with additional basic benchmarking analysis.
    # 1.1 - Added simulated dataset construction and benchmarking.
    # 1.2 - Added MinIC filtering to benchmark assessment. Sorted beginning/end of line for reduced ELMs.
    # 1.3 - Made SimCount a list rather than Integer. Sorted CompariMotif assessment issue.
    # 1.4 - Added ICCut and SLiMLenCut as lists and output columns.
    # 1.5 - Added Summary Results output table. Removed PropRes.
    # 1.6 - Added "simonly" to datatype - calculates both SN and FPR from "sim" data (ignores "ran") to check query bias.
    # 1.7 - Added Benchmarking of ELM datasets without queries.
    # 1.8 - Partially added Benchmarking dataset generation from PPI data and 3DID.
    # 1.9 - Added memsaver option. Replaced SLiMSearch with SLiMProb. Altered default IO paths.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [N] : Add running of actual (Q)SLiMFinder searches between generation and benchmarking.
    # [?] : Add generation of random datasets. (datatype=ran)
    # [Y] : Add generation of simulated datasets = random datasets with some non-random members. (datatype=sim)
    # [Y] : Sort out beginning/end of sequence problem.
    # [Y] : Check that ResMinIC filtering is working.
    # [Y] : Add some summary results table output too. (Reshaped etc.)
    # [Y] : Remove PropRes=T/F if no actual use for it.
    # [ ] : Add ELMdb.py script to download data directly from ELM.
    # [ ] : Consider how to handle ELMs that break rje_slimlist etc.
    # [ ] : Add 3DID benchmarking too.
    # [?] : Check that simulated data can return the query motif. (SLiMSearch)
    # [ ] : Consider splitting ELM instances by length and then trying to recombined with variable wildcards?
    # [ ] : Add output of graphs etc. for assessment results.
    # [ ] : Change summary to have maximum support values etc.
    # [ ] : Build an occurrence-based benchmarking for SLiMSearch-type analysis and also SLiMFinder coverage.
    # [ ] : Improve RunID and Dataset split error messages.
    # [X] : Add option to split Dataset into different things. No: add other datatypes as needed. (e.g. w/o query)
    # [?] : Try using rje_slimcore motifSeq() (in waves) as part of simulation process to check signal.
    # [ ] : Make sure that PPI and 3DID dataset generation is complete and test.
    # [ ] : Check for existence of motif splits file and suggest replacement. (Or use if force=F.)
    # [ ] : Test with memsaver=T.
    # [ ] : Test with SLiMSearch -> SLiMProb replacement.
    # [ ] : Output equivalence file from reduced runs to allow "fair" predictions.
    # [Y] : Add retrieval of UniProt entries directly from UniProt rest.
    # [ ] : Add possible option to save UniProt AccNum.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('SLiMBench', '1.9', 'August 2013', '2012')
    description = 'Short Linear Motif prediction Benchmarking'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show SLiMMaker commandline options?'): out.verbose(-1,4,text=slimmaker.__doc__)
            if rje.yesNo('Show SLiMProb commandline options?'): out.verbose(-1,4,text=slimprob.__doc__)
            if rje.yesNo('Show rje_slimlist commandline options?'): out.verbose(-1,4,text=rje_slimlist.__doc__)
            if rje.yesNo('Show rje_slimcore commandline options?'): out.verbose(-1,4,text=rje_slimcore.__doc__)
            if rje.yesNo('Show RJE_UniProt commandline options?'): out.verbose(-1,4,text=rje_uniprot.__doc__)
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
protected = ['Dataset','RunID','Masking','Build','Chance','RunTime','SeqNum','UPNum','AANum','MotNum','Rank','Sig',
             'Pattern','IC','Occ','Support','UP','ExpUP','Prob','Cloud','CloudSeq','CloudUP','IUP_mean','SA_mean']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMBench Class                                                                                         #
#########################################################################################################################
class SLiMBench(rje_obj.RJE_Object):     
    '''
    SLiMBench Class. Author: Rich Edwards (2012).

    Str:str
    - BenchBase = Basefile for SLiMBench benchmarking output [slimbench]
    - ByCloud = Whether to compress results into clouds prior to assessment (True/False/Both) [Both]
    - CompDB = Motif file to be used for benchmarking (default = reduced elmclass file) []
    - DataType = Type of data to be generated and/or benchmarked (elm/sim/simonly) [elm]
    - DMIData = File of 3DID DMI data ['3did.DMI.csv']
    - ELMClass = Download from ELM website of ELM classes ['elm_classes.tsv']
    - ELMInstance = Download from ELM website of ELM instances ['elm_instances.tsv']
    - ELMPfam = Download from ELM website of ELM Pfam domain interactors ['elm_interaction_domains.tsv']
    - GenPath = Output path for datasets generated with SLiMBench file generator [./SLiMBenchDatasets/]
    - PDBData = File mapping PDB identifiers onto genes/proteins []
    - PfamData = File mapping PFam domains onto genes/proteins (BioMart or HMM search) []
    - RanDir = Output path for creation of randomised datasets [SLiMBenchDatasets/Random/]
    - RandBase = Base for random dataset name [rand]
    - RandSource = Source for new sequences for random datasets [None]
    - SourcePath = Will look in this directory for input files if not found ['SourceData/']
    - Uniprot = File of downloaded UniProt entries (See rje_uniprot for more details) ['ELM.dat']
    - XRefData = File of gene identifier cross-reference data from rje_genemap []
    - 3DIDsql = Path to 3DID sql data. Use rje_mysql sqldump to extract 3DID DMI data. []

    Bool:boolean
    - Benchmark = Whether to perfrom SLiMBench benchmarking assessment [False]
    - Generate = Whether to generate SLiMBench datasets from ELM input [False]
    - Integrity = Whether to quit by default if input integrity is breached [True]
    - Masking = Whether to use SLiMCore masking for query selection [True]
    - NoAmb = Filter out ambiguous patterns [False]
    - Queries = Whether to generate datasets with specific Query proteins [True]
    - Randomise = Whether to generate randomised datasets (part of simulation if simulate=T) [False]
    - Simulate = Whether to generate simulated datasets using reduced ELMs (if found) [False]
    - SLiMMaker = Whether to use SLiMMaker to "reduce" ELMs to more findable SLiMs [True]

    Int:integer
    - MinUPC = Minimum number of UPC for ELM dataset [3]
    - RandReps = Number of replicates for each random (or simulated) datasets [10]

    Num:float
    - MinIC = Min information content for a motif (1 fixed position = 1.0) [2.0]
    - MinResIC = Minimum IC for (Q)SLiMFinder results for benchmark assessment [2.1] (Kept separate to keep 2 fixed pos ELMs for OT etc.)
    
    List:list
    - ByCloud = Whether to compress results into clouds prior to assessment (True/False) [True,False]
    - Dataset = List of headers to split dataset into. If blank, will use datatype defaults. []
    - FlankMask = List of flanking mask options [none,win300,win100,flank5,site]
    - ICCut = Minimum IC for (Q)SLiMFinder results for benchmark assessment [2.0,2.1,3.0]
    - ResFiles = List of (Q)SLiMFinder results files to use for benchmarking [*.csv]
    - RunID = List of factors to split RunID column into (on '.') ['Program','Settings']
    - SearchINI = List of INI files containing search options (should have runid setting) []
    - SimCount = Number of "TPs" to have in dataset [5,10]
    - SimRatios = List of simulated ELM:Random rations [1,4,9,19]
    - SigCut = Significance thresholds to use for assessment [0.05,0.01,0.001,0.0001]
    - SLiMLenCut = List of individual SLiM lengths to return results for (0=All) [0,3,4,5]

    Dict:dictionary
    - UniProt = dictionary of {AccNum:UniProtEntry} matching UniProt field added to ELM Instance table 

    Obj:RJE_Objects
    - CompDB = SLiMList object of Motif file to be used for benchmarking 
    - DB = rje_db.DBase object
    - SLiMMaker = slimmaker.SLiMMaker object
    - UniProt = rje_uniprot.UniProt object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BenchBase','ByCloud','CompDB','DataType','DMIData','ELMClass','ELMInstance','ELMPfam','GenPath',
                        'PDBData','PFanData','RanDir','RandBase','RandSource','SourcePath','Uniprot','3DIDsql']
        self.boollist = ['Benchmark','Generate','Integrity','Masking','NoAmb','Queries','Randomise','Simulate','SLiMMaker']
        self.intlist = ['MinUPC','RandReps']
        self.numlist = ['MinIC','MinResIC']
        self.listlist = ['ByCloud','FlankMask','ICCut','SimCount','ResFiles','RunID','SearchINI','SimRatios','SigCut','SLiMLenCut']
        self.dictlist = ['UniProt']
        self.objlist = ['DB','Uniprot']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'BenchBase':'slimbench','ByCloud':'Both','DataType':'elm',
                     'ELMClass':'elm_classes.tsv','ELMInstance':'elm_instances.tsv','Uniprot':'ELM.dat',
                     'ELMPfam':'elm_interaction_domains.tsv','DMIData':'3did.DMI.csv',
                     'GenPath':rje.makePath('./SLiMBenchDatasets/'),
                     'RanDir':rje.makePath('./SLiMBenchDatasets/Random/'),'RandBase':'rand',
                     'SourcePath':rje.makePath('./SourceData/')})
        self.setBool({'Generate':False,'Integrity':True,'Masking':True,'NoAmb':False,'Queries':True,'SLiMMaker':True})
        self.setInt({'MinUPC':3,'RandReps':10})
        self.setNum({'MinIC':2.0,'MinResIC':2.1})
        self.list['FlankMask'] = string.split('site,flank5,win50,win100,win300,none',',')
        self.list['RunID'] = ['Program','Analysis']
        self.list['SimRatios'] = [1,4,9,19]
        self.list['SimCount'] = [5,10]
        self.list['SigCut'] = [0.05,0.01,0.001,0.0001]
        self.list['SLiMLenCut'] = [0,3,4,5]
        self.list['ICCut'] = [2.0,2.1,3.0]
        self._cmdRead('resfiles=*.csv',type='glist',att='ResFiles')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['SLiMMaker'] = slimmaker.SLiMMaker(self.log,self.cmd_list)
        self.obj['UniProt'] = rje_uniprot.UniProt(self.log,['uniprot=ELM.dat']+self.cmd_list)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['ByCloud','DataType','RandBase'])   # Normal strings
                self._cmdReadList(cmd,'path',['GenPath','RanDir','SourcePath','3DIDsql'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['BenchBase','CompDB','DMIData','ELMClass','ELMInstance','ELMPfam','PDBData','PFanData','RandSource','Uniprot'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['Benchmark','Generate','Integrity','Masking','NoAmb','Queries','Randomise','Simulate','SLiMMaker'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinUPC','RandReps'])   # Integers
                self._cmdReadList(cmd,'float',['MinIC','MinResIC']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['RunID'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['SimRatios','SimCount','SLiMLenCut'])  
                self._cmdReadList(cmd,'nlist',['SigCut','ICCut'])  
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['ResFiles','SearchINI']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.str['DataType'] = self.str['DataType'].lower()  
        self.list['ByCloud'] = []
        if self.getStr('ByCloud').lower() in ['true','both']: self.list['ByCloud'].append(True)
        if self.getStr('ByCloud').lower() in ['false','both']: self.list['ByCloud'].append(False)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            generated = True            
            if self.getBool('Generate'): generated = self.generate()
            if generated and self.getBool('Simulate') or self.getBool('Randomise'): generated = self.simulator()
            if self.getBool('Benchmark'):
                if not generated and (self.i() < 0 or rje.yesNo('SLiMBench generation failed. Quit?')): return False
                if self.getBool('Generate'):
                    self.printLog('#RUN','SLiMBench not currently setup to run (Q)SLiMFinder.')
                    if self.i() < 0 or not rje.yesNo('Confirm (Y) when results copied to proceed, or N to Quit'): self.printLog('#RUN','Run benchmark=T on results later.'); return False
                return self.benchmark()
            if not (self.getBool('Generate') or self.getBool('Simulate') or self.getBool('Randomise') or self.getBool('Benchmark')):
                self.printLog('#RUN','No process chosen for SLiMBench run! Use generate=T, simulate=T, randomise=T and/or benchmark=T')
                return False
            return generated
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### SLiMBench Generator Methods                                                                             #
#########################################################################################################################
    def generate(self): ### Main SLiMBench Generator run method
        '''Main SLiMBench Generator run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupGenerator(): return False
            ### ~ [2] Generate basic ELM datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.fullELMDatasets()
            ### ~ [3] Perform SLiMMaker Reducation of ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.reducedELMs()
            #!# Output equivalencies for "fair" SLiMFinder run? #!#
            ### ~ [4] Run SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.slimProb()
            ### ~ [5] Select Subset of ELMs and generate Query Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Queries'):
                elmlist = []
                for entry in self.db('SLiMProb').entries():
                    if float(entry['IC']) >= self.getNum('MinIC') and float(entry['N_UPC']) >= self.getInt('MinUPC'): elmlist.append(entry['Motif'])
                elmlist.sort()
                self.printLog('#ELM','%d ELMs meet IC >= %.2f & UP >= %d cutoffs' % (len(elmlist),self.getNum('MinIC'),self.getInt('MinUPC')))
                self.queryDatasets(elmlist)
        except: self.errorLog('%s.generate error' % self)
#########################################################################################################################
    def sourceDataFile(self,str,ask=True,expect=True):   ### Returns source data file.                              #V1.9
        '''
        Returns source data file.
        >> str = Key for self.str dictionary corresponding to the source data file being sought.
        >> ask:bool [True] = Whether to ask for file if not found.
        >> expect:bool [True] = Whether source file is expected - raise error if not found.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for checkfile in [self.getStr(str), self.getStr('SourcePath') + self.getStr(str)]:
                if rje.exists(checkfile): break
            while checkfile and ask and not rje.exists(checkfile):
                btext = {True:'exit',False:'ignore'}
                checkfile = rje.choice('%s not found. Enter %s file path (blank to %s)' % (checkfile,str,btext))
            if rje.exists(checkfile): return checkfile
            elif expect and checkfile: raise IOError
            elif expect: os.exit(1)
            else: return ''
        except: self.errorLog('Problem locating %s source data file' % str,quitchoice=True); return ''
#########################################################################################################################
    def setupGenerator(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Load ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            rje.mkDir(self,self.getStr('GenPath'),True)
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            db.baseFile(elmbase)
            ## ~ [1a] Load ELM Classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmc = db.addTable(self.sourceDataFile('ELMClass'),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
            elmc.dataFormat({'#Instances':'int'})
            ## ~ [1b] Load ELM Instances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmi = db.addTable(self.sourceDataFile('ELMInstance'),mainkeys=['ELMIdentifier','ProteinName','Start','End'],datakeys='All',name='ELMInstance')
            elmi.dataFormat({'Start':'int','End':'int'})
            elmi.addField('Primary','ProteinName')
            for entry in elmi.entries():
                try: entry['Primary'] = string.split(entry['Accessions'])[0]
                except: self.printLog('#ERR','Warning! No Accession numbers for %s instance %s' % (entry['ELMIdentifier'],entry['ProteinName']))
            ## ~ [1c] Check and summaries loaded data integrity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            good_elm = []; errors_found = False
            for elm in elmi.indexKeys('ELMIdentifier'):
                icount = len(elmi.indexEntries('ELMIdentifier',elm))
                if elm not in elmc.data(): self.printLog('#ERR','ELM "%s" from Instance file not found in Class file' % elm); errors_found = True
                elif icount != elmc.data(elm)['#Instances']: self.printLog('#ERR','%d "%s" instance in Instance file vs. %d in Class file' % (icount,elm,elmc.data(elm)['#Instances'])); errors_found = True
                else: good_elm.append(elm)
            for elm in elmc.dataKeys():
                if elm in good_elm: continue
                elif not elmc.data(elm)['#Instances']: good_elm.append(elm)
                else: self.printLog('#ERR','No "%s" instance in Instance file vs. %d in Class file' % (elm,elmc.data(elm)['#Instances'])); errors_found = True
            if errors_found:
                if self.getBool('Integrity') and (self.i() < 0 or rje.yesNo('ELM Errors found. Quit?')):
                    self.printLog('#ERR','Setup aborted due to ELM file integrity errors')
                    return False
                elif not self.getBool('Integrity') and self.i() > 0 and rje.yesNo('ELM Errors found. Quit?',default='N'):
                    self.printLog('#ERR','Setup aborted due to ELM file integrity errors')
                    return False
            self.printLog('#ELM','ELM Class and Instance data matched for %d/%d ELMs.' % (len(good_elm),elmc.entryNum()))
                                                                           
            ### ~ [2] Load UniProt Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uni = self.obj['UniProt']
            ufile = self.sourceDataFile('UniProt',ask=False,expect=False)
            if not ufile: uni.setInfo({'DATOut':self.getStr('UniProt')})
            elmacc = elmi.indexKeys('ProteinName') + elmi.indexKeys('Primary') 
            while '' in elmacc: elmacc.remove('') 
            uni.readUniProt(acclist=elmacc,cleardata=False)
            elmaccdict = uni.accDictFromEntries(acc_list=elmacc)
            #self.deBug(elmaccdict)
            ## ~ [2a] Generate ELM UniProt mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmi.addField('UniProt'); missing = []
            unimap = self.dict['UniProt'] = {}
            unifound = 0; unimissed = 0
            for entry in elmi.entries():
                for acc in [entry['ProteinName'],entry['Primary']]:
                    if acc in elmaccdict:
                        u_entry = elmaccdict[acc]
                        u_name = u_entry.info['Name']
                        entry['UniProt'] = u_name
                        unimap[u_name] = u_entry
                if entry['UniProt']: unifound += 1
                else:
                    missing += [entry['ProteinName'],entry['Primary']]
                    self.printLog('#ERR','Could not find UniProt entry %s (%s)' % (entry['ProteinName'],entry['Primary']))
                    unimissed += 1
            self.printLog('#UNI','Mapped %d of %d instances onto UniProt entries' % (unifound,elmi.entryNum()))
            if unimissed: self.printLog('#MISS','UniProt sequences missing for %d entries!' % unimissed)
            if missing:
                missfile = '%s.missing.acc' % rje.baseFile(self.getStr('ELMInstance'))
                rje.backup(self,missfile)
                open(missfile,'w').write(string.join(missing,'\n'))
                self.printLog('#MISS','%d missing UniProt ID/AccNum output to %s' % (len(missing),missfile))
                if (self.getBool('Integrity') and (self.i() < 0 or rje.yesNo('ELM UniProt entries missing. Quit?'))) or (not self.getBool('Integrity') and self.i() > 0 and rje.yesNo('ELM UniProt entries missing. Quit?',default='N')):
                    self.printLog('#ERR','Setup aborted due to ELM file UniProt errors')
                    return False

            ### ~ [3] Generate ELM motif file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            motif_file = '%s.motifs.txt' % elmbase
            if rje.checkForFile(motif_file) and not self.force(): self.printLog('#MOTIF','%s file found (force=F).' % motif_file)
            else:
                rje.backup(self,motif_file)
                motif_out = []
                for elm in elmc.dataKeys():
                    entry = elmc.data(elm)
                    motif_out.append('%s  %s  # %s [%d ELM instances]' % (entry['ELMIdentifier'],entry['Regex'],entry['Description'],entry['#Instances']))
                open(motif_file,'w').write(string.join(motif_out,'\n'))
                self.printLog('#MOTIF','%s motif patterns output to %s' % (elmc.entryNum(),motif_file))

            return True     # Setup successful            
        except: self.errorLog('Problem during %s setupGenerator.' % self); return False  # Setup failed
#########################################################################################################################
    def fullELMDatasets(self):     ### Generates Full ELM Datasets
        '''Generates Full ELM Datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            datadir = rje.makePath('%sELM_Datasets/' % self.getStr('GenPath'))
            rje.mkDir(self,datadir,True)
            elmc = db.getTable('ELMClass')
            elmc.addField('FullSeqNum')
            elmi = db.getTable('ELMInstance')
            uni = self.obj['UniProt']; unimap = self.dict['UniProt']
            ### ~ [1] ~ Generate ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfilex = 0
            for elm in elmc.dataKeys():
                elmseqfile = '%s%s.fas' % (datadir,elm)
                if not self.force() and rje.checkForFile(elmseqfile):
                    elmc.data(elm)['FullSeqNum'] = rje_seq.SeqCount(self,elmseqfile)
                    seqfilex += 1
                else:
                    rje.backup(self,elmseqfile,appendable=False)
                    elmseqs = []
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        try: seq = unimap[entry['UniProt']].obj['Sequence']
                        except: self.printLog('#ERR','Failed to get sequence for %s:%s' % (elm,entry['ProteinName'])); continue
                        if seq not in elmseqs: elmseqs.append(seq)
                    elmc.data(elm)['FullSeqNum'] = len(elmseqs)
                    if elmseqs:
                        ELMSEQ = open(elmseqfile,'w')
                        for seq in elmseqs: ELMSEQ.write('>%s\n%s\n' % (seq.name(),seq.getSequence()))
                        ELMSEQ.close()
                        self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                        seqfilex += 1
                    else: self.printLog('#ELM','No sequences to output for %s' % elm)
            self.printLog('#DSETS','%d ELM sequence datasets' % seqfilex)
        except: self.errorLog('Problem during %s.fullELMDatasets().' % self); return False  
#########################################################################################################################
    def reducedELMs(self):  ### Performs ELM reduction using SLiMMaker
        '''Performs ELM reduction using SLiMMaker.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('SLiMMaker'): return
            db = self.obj['DB']
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            motif_file = '%s.reduced.motifs.txt' % elmbase
            datadir = rje.makePath('%sELM_Datasets_Reduced/' % self.getStr('GenPath'))
            rje.mkDir(self,datadir,True)
            instdir = rje.makePath('%sELM_Instances/' % self.getStr('GenPath'))
            rje.mkDir(self,instdir,True)
            elmc = db.getTable('ELMClass')
            elmc.addField('SLiMMaker')
            elmc.addField('ReducedInstance')
            elmc.addField('ReducedSeqNum')
            elmi = db.getTable('ELMInstance')
            elmi.addField('Instance')
            elmi.addField('Reduced')
            uni = self.obj['UniProt']; unimap = self.dict['UniProt']
            maker = self.obj['SLiMMaker']
            seqfilex = 0
            ### ~ [1] ~ Generate Reduced ELMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for elm in elmc.dataKeys():
                instances = []
                for entry in elmi.indexEntries('ELMIdentifier',elm):
                    try: seq = unimap[entry['UniProt']].obj['Sequence']
                    except:
                        self.printLog('#ERR','Failed to get sequence for %s:%s' % (elm,entry['ProteinName']))
                        entry['Instance'] = ''; entry['Reduced'] = 'N'; continue
                    entry['Instance'] = seq.getSequence()[entry['Start']-1:entry['End']]
                    if entry['Start'] == 1: entry['Instance'] = '^%s' % entry['Instance']
                    else: entry['Instance'] = 'X%s' % entry['Instance']
                    if entry['End'] == seq.seqLen(): entry['Instance'] = '%s$' % entry['Instance']
                    else: entry['Instance'] = '%sX' % entry['Instance']
                    instances.append(entry['Instance'])
                open('%s%s.instances.txt' % (instdir,elm),'w').write(string.join(instances,'\n'))
                #self.deBug('%s: %s' % (elmc.data(elm)['Regex'],instances))
                # Consider splitting by length and then recombine with wildcard spacers?! #
                maker.list['Peptides'] = instances
                slim = maker.run(iterate=True)[0]       #!# Make iterate an option? #!#
                if not slim: maker.list['Peptides'] = []
                elmc.data(elm)['SLiMMaker'] = slim
                elmc.data(elm)['ReducedInstance'] = len(maker.list['Peptides'])
                for entry in elmi.indexEntries('ELMIdentifier',elm):
                    if entry['Instance'] and entry['Instance'] in maker.list['Peptides']: entry['Reduced'] = 'Y'
                    else: entry['Reduced'] = 'N'
                #self.deBug('%s: %s' % (elmc.data(elm)['Regex'],slim))
                self.printLog('#MAKE','%s: %s -> "%s" = %d -> %d instances' % (elm,elmc.data(elm)['Regex'],slim,elmc.data(elm)['#Instances'],elmc.data(elm)['ReducedInstance']))
                ## ~ [1a] Generate Reduced ELM Dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elmseqfile = '%s%s.reduced.fas' % (datadir,elm)
                if not self.force() and rje.checkForFile(elmseqfile):
                    elmc.data(elm)['FullSeqNum'] = rje_seq.SeqCount(self,elmseqfile)
                    seqfilex += 1
                else:
                    rje.backup(self,elmseqfile,appendable=False)
                    elmseqs = []
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        if entry['Reduced'] == 'N': continue
                        try: seq = unimap[entry['UniProt']].obj['Sequence']
                        except: self.errorLog('Failed to get sequence for %s:%s' % (elm,entry['ProteinName']),quitchoice=False); continue
                        if seq not in elmseqs: elmseqs.append(seq)
                    elmc.data(elm)['ReducedSeqNum'] = len(elmseqs)
                    if elmseqs:
                        ELMSEQ = open(elmseqfile,'w')
                        for seq in elmseqs: ELMSEQ.write('>%s\n%s\n' % (seq.name(),seq.getSequence()))
                        ELMSEQ.close()
                        self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                        seqfilex += 1
                    else: self.printLog('#ELM','No sequences to output for %s' % elm)
            self.printLog('#DSETS','%d Reduced ELM sequence datasets' % seqfilex)
            ### ~ [2] ~ Create Reduced ELM Motif File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile(motif_file) and not self.force(): self.printLog('#MOTIF','%s file found (force=F).' % motif_file)
            else:
                rx = 0
                rje.backup(self,motif_file)
                motif_out = []
                for elm in elmc.dataKeys():
                    entry = elmc.data(elm)
                    if not entry['SLiMMaker']: continue
                    motif_out.append('%s  %s  # %s. ELM Definition: %s [%d -> %d ELM instances]' % (entry['ELMIdentifier'],entry['SLiMMaker'],entry['Description'],entry['Regex'],entry['#Instances'],entry['ReducedInstance']))
                    rx += 1
                open(motif_file,'w').write(string.join(motif_out,'\n'))
                self.printLog('#MOTIF','%d of %d motif patterns output to %s' % (rx,elmc.entryNum(),motif_file))
            return True
        except: self.errorLog('Problem during %s.reducedELMs().' % self); return False  
#########################################################################################################################
    def slimProb(self,datadir='ELM_Datasets_Reduced',slimkey='SLiMMaker'):     ### Performs SLiMProb of specific motifs against their datasets
        '''Performs SLiMProb of specific motifs against their datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            basefile = rje.makePath('%s%s' % (self.getStr('GenPath'),datadir),wholepath=True)
            ssfile = '%s.slimprob.occ.csv' % basefile
            if rje.checkForFile(ssfile):
                if not self.force():
                    self.printLog('#SEARCH','SLiMProb file %s found (force=F).' % ssfile)
                    db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Dataset','RunID','Motif'],datakeys='All',name='SLiMProb')
                    db.addTable(ssfile,mainkeys=['Dataset','RunID','Motif','Seq','Start_Pos'],datakeys='All',name='SLiMProbOcc')
                    return True
                else: rje.backup(self,ssfile,appendable=False)
            #self.deBug(rje.checkForFile(ssfile))
            datadir = rje.makePath(basefile)
            elmc = db.getTable('ELMClass')
            elmi = db.getTable('ELMInstance')
            ex = 0; rx = 0
            if not self.list['SearchINI']: self.list['SearchINI'] = ['None']
            #self.deBug('%s' % self.list['SearchINI'])
            ### ~ [1] ~ Perform SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for elmfas in glob.glob('%s*.fas' % datadir):
                elm = string.split(os.path.basename(elmfas),'.')[0]
                if elm not in elmc.dataKeys(): self.errorLog('Failed to find ELM "%s" in ELMClass table.' % elm,quitchoice=False); continue
                entry = elmc.data(elm)
                slim = entry[slimkey]
                if not slim: self.printLog('#SKIP','No "%s" SLiM to search for ELM "%s".' % (slimkey,elm)); continue
                ex += 1
                for ini in self.list['SearchINI']:
                    sscmd = ['minic=1.1','extras=0']+self.cmd_list+['resfile=%s' % ssfile,'append=T','motifs=','seqin=%s' % elmfas,'resdir=%sSLiMProb/' % datadir]
                    if ini.lower() in ['','none','n']: sscmd += ['runid=SLiMBench']
                    else: sscmd += rje.getCmdList(['ini=%s' % ini])
                    #self.deBug(sscmd)
                    ss = slimprob.SLiMProb(self.log,sscmd+['debug=F'])
                    #self.deBug(ss.cmd_list)
                    slimlist = ss.obj['SlimList']
                    slimlist.list['Motif'] = []
                    slimlist._addMotif(elm,slim,check=True)
                    if not slimlist.list['Motif']: self.printLog('#CUT','%s %s did not make the cut' % (elm,slim)); continue
                    ss.run()
                    rx += 1
            self.printLog('#SS','%d SLiMProb runs for %d ELMs (%s vs %s)' % (rx,ex,slimkey,datadir))
            ### ~ [2] ~ Load SLiMProb Results into Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Dataset','RunID','Motif'],datakeys='All',name='SLiMProb')
            db.addTable(ssfile,mainkeys=['Dataset','RunID','Motif','Seq','Start_Pos'],datakeys='All',name='SLiMProbOcc')
            return True
        except: self.errorLog('Problem during %s.slimProb().' % self); return False  
#########################################################################################################################
    def queryDatasets(self,elmlist,datadir='ELM_Datasets_Reduced'):  ### Produces ELM Query Datasets
        '''
        Produces ELM Query Datasets.
        >> datadir:path = Path to ELM dataset files
        >> elmlist:list = List of ELM datasets to make files for
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db('ELMInstance')
            indir = rje.makePath('%s%s' % (self.getStr('GenPath'),datadir))
            outdir = rje.makePath('%sELM_Datasets_WithQueries' % (self.getStr('GenPath')))
            rje.mkDir(self,outdir,True)
            self.list['FlankMask'] = string.split(string.join(self.list['FlankMask']).lower())
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','dna=F'])
            ### ~ [1] ~ Generate new ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qx = 0
            for elm in elmlist:
                #self.deBug(elm)
                for efile in glob.glob('%s%s.*fas' % (indir,elm)):
                    #self.deBug(efile)
                    self.printLog('#ELM','Process %s: %s' % (elm,efile))
                    seqlist.loadSeq(efile,nodup=True,clearseq=True,mode='list')     ### Loads sequences from file
                    #self.deBug(seqlist.list['Seq'])
                    eseqs = seqlist.list['Seq'][0:]
                    qseqs = seqlist.list['Seq'][0:]
                    ## ~ [1a] ~ Take each Query in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    while qseqs:
                        query = seqlist.getSeq(qseqs.pop(0))
                        #self.deBug(query)
                        qname = string.split(query[0])[0]
                        qlen = len(query[1])
                        #self.deBug(qname)
                        self.printLog('#QRY','%s - Query %s' % (elm,qname)); qx += 1
                        #self.deBug('%s: %s' % (elm,self.db('ELMClass').data(elm)['SLiMMaker']))
                        for flank in self.list['FlankMask']:
                            ## ~ [1b] ~ Establish Flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            #self.deBug(flank)
                            mask = []
                            for entry in db.indexEntries('ELMIdentifier',elm):
                                if entry['ProteinName'] in string.split(qname,'__') or entry['Primary'] in string.split(qname,'__'):
                                    occ = (int(entry['Start'])-1,int(entry['End']))
                                    if occ not in mask: mask.append(occ)
                            mask.sort()
                            #self.deBug(mask)
                            if not mask: self.printLog('#ERR','Cannot find any %s instances of %s in ELMInstances!' % (qname,elm)); continue
                            if flank == 'none': mask = [(0,qlen)]
                            elif flank == 'site': pass
                            elif flank[:5] == 'flank':
                                win = string.atoi(flank[5:])
                                newmask = []
                                for old in mask: newmask.append((max(0,old[0]-win),min(qlen,old[1]+win)))
                                mask = newmask
                            elif flank[:3] == 'win':
                                win = string.atoi(flank[3:])
                                #self.deBug(win)
                                masklen = 0
                                for old in mask: masklen += (old[1] - old[0])
                                while masklen < win and masklen < qlen:
                                    newmask = []
                                    for old in mask:
                                        start = max(0,old[0]-1)
                                        end = min(qlen,old[1]+1)
                                        if newmask and start < newmask[-1][1]: newmask[-1] = (newmask[-1][0],end)
                                        else: newmask.append((start,end))
                                    #self.deBug(newmask)
                                    mask = newmask
                                    masklen = 0
                                    for old in mask: masklen += (old[1] - old[0])
                                    #self.deBug(masklen)
                                    continue
                                    #!# Delete below #!#
                                    splits = len(mask) * 2
                                    #self.deBug(splits)
                                    splitlen = max(1,int(0.99+(win-masklen)/float(splits)))
                                    #self.deBug(splitlen)
                                    for old in mask:
                                        if old[1] == qlen: start = max(0,old[0]-2*splitlen)
                                        else: start = max(0,old[0]-splitlen)
                                        if old[0] == 0: end = min(qlen,old[1]+2*splitlen)
                                        else: end = min(qlen,old[1]+splitlen)
                                        if newmask and start < newmask[-1][1]: newmask[-1] = (newmask[-1][0],end)
                                        else: newmask.append((start,end))
                                        #self.deBug(newmask)
                                    mask = newmask
                                    masklen = 0
                                    for old in mask: masklen += (old[1] - old[0])
                                    #self.deBug(masklen)
                            #self.deBug(mask)
                            ## ~ [1c] ~ Mask Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            maskseq = query[1].lower()
                            for case in mask:
                                #self.deBug(case); self.deBug(maskseq.upper()[case[0]:case[1]])
                                maskseq = maskseq[:case[0]] + maskseq.upper()[case[0]:case[1]] + maskseq[case[1]:]
                            #self.deBug(maskseq)
                            ## ~ [1d] ~ Output seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            outfile = '%s%s.%s.%s.fas' % (outdir,elm,qname,flank)
                            OUT = open(outfile,'w')
                            OUT.write('>%s\n%s\n' % (query[0],maskseq))
                            for seq in eseqs[1:]: OUT.write('>%s\n%s\n' % seqlist.getSeq(seq))
                            OUT.close()
                            self.printLog('#FAS','%d sequences output to %s' % (len(eseqs),outfile))
                        eseqs.append(eseqs.pop(0))
            self.printLog('#QFAS','Fasta files made for %d Queries (%d ELMs)' % (qx,len(elmlist)))
            return True                                
        except: self.errorLog('Problem during %s.queryDatasets().' % self); return False  
#########################################################################################################################
    def randomiser(self,queries=[],seqnum=0,basefile='rand',start=0,end=0):   ### Generates random datasets using query masking
        '''
        Generates random datasets using query masking.
        >> queries:list of Sequence objects = First sequence will be the query for region masking.
        >> seqnum:int = Total number of sequences - add random sequences from self.obj['SeqList'] until complete
        >> basefile:str ['rand'] = Base for output files: basefile.region.fas
        >> start:int [0] = Start position of motif for region masking
        >> end:int [0] = Start position of motif for region masking
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            randsource = self.obj['SeqList']
            randseq = randsource.seqs()[0:]
            for qry in queries: randseq.remove(qry)
            while randseq and len(queries) < seqnum:
                r = random.randint(0,len(randseq)-1)
                queries.append(randseq.pop(r))
            if len(queries) < seqnum: raise ValueError
            ### ~ [1] Output a file for each region masking used ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for flank in self.list['FlankMask']:
                ## ~ [1a] ~ Establish occ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qry = queries[0]
                (qname,qseq) = randsource.getSeq(qry,'tuple')
                qlen = len(qseq)    #randsource.seqLen(qry)
                if not start:
                    start = random.randint(0,qlen-10)
                    end = start + 10
                mask = [(start,end)]
                ## ~ [1b] ~ Set up Region ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if flank == 'none': mask = [(0,qlen)]
                elif flank == 'site': pass
                elif flank[:5] == 'flank':
                    win = string.atoi(flank[5:])
                    newmask = []
                    for old in mask: newmask.append((max(0,old[0]-win),min(qlen,old[1]+win)))
                    mask = newmask
                elif flank[:3] == 'win':
                    win = string.atoi(flank[3:])
                    #self.deBug(win)
                    masklen = 0
                    for old in mask: masklen += (old[1] - old[0])
                    while masklen < win and masklen < qlen:
                        newmask = []
                        for old in mask:
                            start = max(0,old[0]-1)
                            end = min(qlen,old[1]+1)
                            if newmask and start < newmask[-1][1]: newmask[-1] = (newmask[-1][0],end)
                            else: newmask.append((start,end))
                        #self.deBug(newmask)
                        mask = newmask
                        masklen = 0
                        for old in mask: masklen += (old[1] - old[0])
                ## ~ [1c] ~ Mask Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                maskseq = qseq.lower()  #qry.info['Sequence'].lower()
                for case in mask:
                    #self.deBug(case); self.deBug(maskseq.upper()[case[0]:case[1]])
                    maskseq = maskseq[:case[0]] + maskseq.upper()[case[0]:case[1]] + maskseq[case[1]:]
                #self.deBug(maskseq)
                ## ~ [1d] ~ Output seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                outfile = '%s.%s.%s.fas' % (basefile,string.split(qname)[0],flank)
                OUT = open(outfile,'w')
                OUT.write('>%s\n%s\n' % (qname,maskseq))
                for seq in queries[1:]:
                    (seqname,seqseq) = randsource.getSeq(seq,'tuple')
                    OUT.write('>%s\n%s\n' % (seqname,seqseq.upper()))
                OUT.close()
                self.printLog('#FAS','%d sequences output to %s' % (len(queries),outfile))
            return True
        except: self.errorLog('%s.randomiser error' % self); return False
#########################################################################################################################
    def simulator(self):   ### Generates simulated datasets
        '''Generates simulated datasets.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~ RANDOMISE/SIMULATE SETUP ~~~~~~~~~~~~~~~~~~ #')
            core = rje_slimcore.SLiMCore(self.log,['logmask=T','maskpickle=T','resdir=%s' % self.getStr('GenPath')]+self.cmd_list+['efilter=F','slimbuild=F','walltime=0'])
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F','seqin=%s' % self.getStr('RandSource')]
            randsource = self.obj['SeqList'] = rje_seqlist.SeqList(self.log,seqcmd+['mode=file'])
            core.dict['MotifSeq'] = {}  #{Pattern:File}
            rje.mkDir(self,self.getStr('RanDir'))
            ## ~ [0a] Setup DB Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.baseFile(self.getStr('RandBase'))
            db = self.obj['DB']
            db.baseFile(self.getStr('RandBase'))
            for table in db.tables(): db.deleteTable(table)     # Clear db before proceeding
            ## ~ [0b] Load ELM Classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmc = db.addTable(self.getStr('ELMClass'),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
            elmc.dataFormat({'#Instances':'int'})
            ## ~ [0c] Load Reduced Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            reduced = '%s.reduced.motifs.txt' % elmbase
            genpath = rje.makePath('%sSourceMatch/' % self.getStr('GenPath'))
            rje.mkDir(self,genpath)
            sourcebase = rje.baseFile(self.getStr('RandSource'),True)
            if self.getBool('Simulate'):
                if not rje.checkForFile(reduced):
                    self.printLog('#MOTIF','Reduced Motif file "%s" not found: no simulation.' % reduced)
                    self.bool['Simulate'] = False
                else:
                    slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list)
                    slimlist.loadMotifs(reduced)    # Will automatically apply minic filter
                    for slim in slimlist.slims():
                        slimfas = '%s%s.%s.fas' % (genpath,slim.getStr('Name'),sourcebase)
                        if not self.force() and rje.checkForFile(slimfas):
                            self.printLog('#FAS','%s found. (Force=F)' % slimfas)
                            continue
                        core.dict['MotifSeq'][slim.pattern()] = slimfas  #{Pattern:File}                                                
            ## ~ [0d] Perform randomised data generation only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getBool('Simulate'):
                for rep in range(self.getInt('RandReps')):
                    r = random.randint(0,randsource.seqNum()-1)
                    qry = randsource.seqs()[r]
                    for ratio in self.list['SimRatios']:
                        for simcount in self.list['SimCount']:
                            seqnum = simcount*(1+ratio)
                            rbase = '%s%s.r%s.n%d' % (self.getStr('RanDir'),self.getStr('RandBase'),rje.preZero(rep+1,self.getInt('RandReps')),seqnum)
                            if not self.randomiser(queries=[qry],seqnum=seqnum,basefile=rbase): return False
                return True
            ### ~ [1] Run SLiMCore MotifSeq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if core.dict['MotifSeq']:
                #!# Add log output for progress #!#
                core.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
                core.setupBasefile()
                core.stat['StartTime'] = time.time()
                core.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
                if core.opt['Masked']: core.obj['SeqList'].saveFasta(seqfile='%s.%s.masked.fas' % (core.info['Basefile'],core.maskText()))
                core.motifSeq()
            ### ~ [2] Run SLiMProb on MotifSeq files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            searchdir = rje.makePath('%sSourceSearch/' % self.getStr('GenPath'))
            searchbase = '%s%s.%s.randsource' % (searchdir,self.baseFile(),sourcebase)
            ssfile = '%s.slimprob.occ.csv' % searchbase
            if rje.checkForFile(ssfile) and not self.force():
                self.printLog('#SEARCH','SLiMProb file %s found (force=F).' % ssfile)
            else:
                rje.backup(self,ssfile,appendable=False)
                rje.mkDir(self,'%sSLiMProb/' % searchdir)
                ## ~ [2a] Perform SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for slim in slimlist.slims():
                    slimfas = '%s%s.%s.fas' % (genpath,slim.getStr('Name'),sourcebase)
                    if open(slimfas,'r').read()[:1] != '>': continue
                    sscmd = ['extras=0']+self.cmd_list+['resfile=%s' % ssfile,'append=T','motifs=','seqin=%s' % slimfas,'resdir=%sSLiMProb/' % searchdir,'occupc=T']
                    ss = slimprob.SLiMProb(self.log,sscmd+['debug=F'])
                    ss.info['Basefile'] = ''
                    ss.obj['SlimList'].list['Motif'] = [slim]
                    ss.run()
            ## ~ [2b] ~ Load SLiMProb Results into Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sdb = db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Motif'],datakeys='All',name='SourceSearch')
            odb = db.addTable(ssfile,mainkeys=['Motif','Seq','Start_Pos'],datakeys='All',name='SourceOcc')
            ### ~ [3] Work through each OK ELM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['SimCount'].sort()
            seqdict = randsource.seqNameDic()
            dx = 0
            for elm in sdb.dataKeys():
                sentry = sdb.data(elm)
                if float(sentry['IC']) < self.getNum('MinIC'):
                    self.printLog('#IC','%s IC < %d: rejected!' % (elm,self.getNum('MinIC')))
                    continue
                for rep in range(self.getInt('RandReps')):
                    upc = []
                    queries = []
                    (start,end) = (0,0)
                    for simcount in self.list['SimCount']:
                        ## Check UPC
                        if int(sentry['N_UPC']) < simcount:
                            self.printLog('#UPC','%s has too few UPC for %d sim: rejected!' % (elm,simcount))
                            continue
                        for ratio in self.list['SimRatios']:
                            if queries: queries = queries[:1]; upc = upc[:1]
                            occ = odb.indexEntries('Motif',elm)
                            while occ and len(queries) < simcount:
                                oentry = occ.pop(random.randint(0,len(occ)-1))
                                if oentry['UPC'] in upc: continue
                                if not upc: (start,end) = (int(oentry['Start_Pos'])-1,int(oentry['End_Pos']))
                                upc.append(oentry['UPC'])
                                queries.append(seqdict[oentry['Seq']])
                            if len(queries) < simcount: raise ValueError
                            #r = random.randint(0,randsource.seqNum()-1)
                            #qry = randsource.seq[r]
                            qry = queries[0]
                            seqnum = simcount*(1+ratio)
                            ## Simulated Dataset output ##
                            rbase = '%s%s.sim.r%s.p%d.n%d' % (self.getStr('RanDir'),elm,rje.preZero(rep+1,self.getInt('RandReps')),simcount,seqnum)
                            if not self.randomiser(queries,seqnum,rbase,start,end): return False
                            ## Matching Random Dataset output ##
                            rbase = '%s%s.ran.r%s.p%d.n%d' % (self.getStr('RanDir'),elm,rje.preZero(rep+1,self.getInt('RandReps')),simcount,seqnum)
                            if not self.randomiser([qry],seqnum,rbase,start,end): return False
                            dx += 1
            self.printLog('#OUT','%s sets of simulated and randomised data output to %s' % (rje.iStr(dx),self.getStr('RanDir')))
            return True
        except: self.errorLog('%s.simulator error' % self); return False
#########################################################################################################################
    def dmi(self):   ### Make 3DID DMI benchmarks   #!# Currently taken purely from rje_misc
        '''Make 3DID DMI benchmarks.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = rje_db.Database(self.log,self.cmd_list)
            base = db.info['Basefile'] = '3did.DMI'
            dmi = db.addTable('3did.DMI.csv',name='DMI',mainkeys=['ID'])
            for entry in dmi.entries(): entry['PDB'] = entry['PDB'].upper()
            dmi.makeField('#PDB#:#dch#','DomPDB'); dmi.makeField('#PDB#:#pch#','MotPDB')
            pdb = db.addTable('3did.DMI.pdb.ens_HUMAN.mapping.tdt',name='PDB',mainkeys=['Query'])
            #pdb.dropEntriesDirect('Hit',['None'])
            for entry in pdb.entries(): entry['Query'] = string.split(entry['Query'],'|')[0]
            xref = db.addTable('../../HGNC/genemap.111005.data.tdt',name='HGNC',mainkeys=['Gene'])
            ### ~ [1] ~ Join Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            temp = db.joinTables(name='Temp',join=[(dmi,'DomPDB'),(pdb,'Query',['Hit'])],newkey=['ID'],empties=False)
            temp.renameField('Hit','DomSeq')
            temp2 = db.joinTables(name='Temp2',join=[(temp,'DomSeq'),(xref,'EnsLoci',['Gene'])],newkey=['ID','Gene'],empties=False)
            map = db.joinTables(name='PDBMap',join=[(temp2,'MotPDB'),(pdb,'Query',['Hit'])],newkey=['ID','Gene'],empties=False)
            map.renameField('Hit','MotSeq'); map.renameField('Gene','Hub')
            map.saveToFile()
            ### ~ [2] ~ Make sequence query files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map.dropEntriesDirect('DomSeq',['None'])
            map.dropEntriesDirect('MotSeq',['None'])
            map.makeField('#MotSeq#','qry')
            for entry in map.entries(): entry['qry'] = string.split(entry['qry'],'_')[-1]
            map.makeField('#peptide#.#Hub#.#qry#','dataset')
            map.dataFormat({'pstart':'int','pend':'int'})
            ensloci = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=../../EnsEMBL/ens_HUMAN.loci.fas','autoload=T','seqmode=file'])
            ensdict = ensloci.seqNameDic()
            SEQOUT = open('%s.pseq.fas' % base,'w')
            flank = 5
            ppidir = '../../Pingu/PPIFas/'
            rje.mkDir(self,'3DID_PPIFas/',log=True)
            for dataset in map.index('dataset'):
                pseq = ''; qry = ''; hub = ''; ex = 0
                for entry in map.indexEntries('dataset',dataset):
                    #while len(seq) < entry['pend']: seq += 'x'
                    qry = entry['MotSeq']; hub = entry['Hub']; ex += 1
                pseq = string.join(map.indexDataList('dataset',dataset,'pseq'))
                self.printLog('#DSET','%s: %d entries' % (dataset,ex))
                SEQOUT.write('>%s %s\n%s\n' % (dataset,qry,pseq))
                ppifas = ppidir + hub + '.ppi.fas'
                outfas = '3DID_PPIFas/%s.ppi.fas' % dataset
                if not os.path.exists(ppifas):
                    self.printLog('#NOPPI','No PPI file for hub %s' % hub)
                    continue
                ppiseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s' % ppifas,'autoload=T','seqnr=F','accnr=F','usecase=T'])
                ppidict = ppiseq.seqNameDic()
                if qry in ppidict: qseq = ppidict[qry]
                elif qry in ensdict:
                    (name,sequence) = ensloci.getSeq(ensdict[qry],'tuple')
                    self.printLog('#QRY','Adding Query %s to %s' % (qry,dataset))
                    qseq = ppiseq._addSeq(name,sequence)                    
                else:
                    self.printLog('#NOQRY','Cannot find query sequence %s' % qry)
                    continue
                if qseq in ppiseq.seq: ppiseq.seq.remove(qseq)
                else: self.printLog('#QRY','Adding Query %s to %s' % (qry,dataset))
                newseq = qseq.info['Sequence'].lower()
                uppseq = qseq.info['Sequence'].upper()
                for frag in map.indexDataList('dataset',dataset,'pseq'):
                    x = 0
                    self.deBug(newseq)
                    while uppseq.find(frag,x) >= 0:
                        i = uppseq.find(frag,x)
                        newseq = newseq[:max(i-flank,0)] + uppseq[max(i-flank,0):i+len(frag)+flank] + newseq[i+len(frag)+flank:]
                        self.deBug(newseq)
                        x = i + 1
                PPI = open(outfas,'w')
                PPI.write('>%s\n%s\n' % (qseq.info['Name'],newseq))
                for seq in ppiseq.seq: PPI.write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence'].upper()))
                PPI.close()
                self.printLog('#FAS','%s: %d seq' % (outfas,ppiseq.seqNum()+1))
                    
            SEQOUT.close()

            
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <4> ### SLiMBench Benchmarking Methods                                                                          #
#########################################################################################################################
    def benchmark(self): ### Main SLiMBench Benchmarking run method                                                 #V1.0
        '''Main SLiMBench Benchmarking run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            if not self.setupBenchmarking(): return False
            ### ~ [2] Load and Process (Q)SLiMFinder search data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.loadSLiMPredictions(): return False
            if self.getBool('NoAmb'): self.filterAmb(save=self.getBool('Test'))
            ### ~ [3] Perform CompariMotif search of SLiM Predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.compariMotif(): return False
            ### ~ [4] Rate SLiM Predictions as TP/FP/OT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.ratings(): return False
            ### ~ [5] Summary Performance Table (ELM only?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.summaryTables(): return False  #!# Add some summary results table output too. (Reshaped etc.) #!#
            ### ~ [6] Generate Sn/Sp etc. stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            afile = '%s.assessment.tdt' % self.getStr('BenchBase')
            if not self.force() and os.path.exists(afile):
                self.printLog('#BENCH','%s found! (force=False)' % afile)
                keep_unique = self.list['RunID'] + ['Region','ICCut','LenCut','SigCut','ByCloud']
                if not self.getBool('Queries'): keep_unique.remove('Region')
                if self.getStr('DataType')[:3] == 'sim': keep_unique += ['PosNum','SeqNum']
                db.addTable(afile,mainkeys=keep_unique,datakeys='All',name='results',replace=True)
            elif self.getBool('MemSaver'):
                fullresdb = db.copyTable('results','full_results')
                self.db('results').makeField('#%s#' % string.join(self.list['RunID'],'#|#'),'MemSaver')
                splitres = db.splitTable('results',field,asdict=True,keepfield=False)
                self.db().deleteTable('results')
                for ic_cut in self.list['ICCut']:
                    self.num['MinResIC'] = ic_cut
                    for len_cut in self.list['SLiMLenCut']:
                        self.int['ResLen'] = len_cut
                        for table in splitres.values():
                            self.db().copyTable(table,'results')
                            if not self.filterMinICLen(ic_cut,len_cut,save=self.getBool('Test') or self.getBool('DeBug') or self.v()>1): return False
                            for bycloud in self.list['ByCloud']:
                                for sig_cut in self.list['SigCut']:
                                    if not self.assessment(sig_cut,bycloud,test=self.getBool('Test')): return False
                            self.db().deleteTable('results')
                db.addTable(afile,mainkeys=keep_unique,datakeys='All',name='results',replace=True)
            else:
                fullresdb = self.db().copyTable('results','full_results')
                for ic_cut in self.list['ICCut']:
                    self.num['MinResIC'] = ic_cut
                    for len_cut in self.list['SLiMLenCut']:
                        self.int['ResLen'] = len_cut
                        self.db().deleteTable('results')
                        self.db().copyTable('full_results','results')
                        if not self.filterMinICLen(ic_cut,len_cut,save=self.getBool('Test') or self.getBool('DeBug') or self.v()>1): return False
                        for bycloud in self.list['ByCloud']:
                            for sig_cut in self.list['SigCut']:
                                if not self.assessment(sig_cut,bycloud,test=self.getBool('Test')): return False
                self.db().deleteTable('results')
                self.db().copyTable('full_results','results')
        except: self.errorLog('%s.benchmark error' % self)
#########################################################################################################################
    def setupBenchmarking(self):    ### Main class setup method.                                                    #V1.0
        '''Main class setup method.'''
        try:### ~ [1] Load ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.baseFile(self.getStr('BenchBase'))
            db = self.obj['DB']
            db.baseFile(self.getStr('BenchBase'))
            for table in db.tables(): db.deleteTable(table)     # Clear db before proceeding
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            ## ~ [1a] Load ELM Classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmc = db.addTable(self.getStr('ELMClass'),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
            elmc.dataFormat({'#Instances':'int'})
            #i# See setupGenerator() method if ELM instance and/or UniProt mapping data is needed.
            ### ~ [2] Load ELM motif file for comparison ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Establish file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            if not rje.checkForFile(self.getStr('CompDB')): self.str['CompDB'] = '%s.reduced.motifs.txt' % elmbase
            if not rje.checkForFile(self.getStr('CompDB')): self.str['CompDB'] = self.getStr('ELMClass')
            if not rje.checkForFile(self.getStr('CompDB')): self.printLog('#MOTIF','Motif input file (and ELM class file) not found'); raise IOError
            ## ~ [2b] Load Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['CompDB'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
            self.obj['CompDB'].loadMotifs(self.str['CompDB'])
            return True     # Setup successful            
        except: self.errorLog('Problem during %s setupBenchmarking.' % self); return False  # Setup failed
#########################################################################################################################
    def loadSLiMPredictions(self):  ### Load SLiM Predictions into Database object and establish runs etc.          #V1.0
        '''Load SLiM Predictions into Database object and establish runs etc.'''
        try:### ~ [0] Setup Header Splitting Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('DataType')[:3] == 'sim': dsethead = ['Motif','RType','Rep','PosNum','Query','Region']
            elif self.getBool('Queries'): dsethead = ['Motif','Query','Region']
            else: dsethead = ['Motif']
            for head in self.list['RunID']:
                if head in protected + dsethead:
                    self.errorLog('Cannot have protected/dataset field "%s" in RunID split.' % head, printerror=False)
                    raise ValueError
            ### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ LOAD PREDICTIONS ~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            db.setStr({'Delimit':self.getStr('Delimit')})
            rfile = '%s.results.%s' % (db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and os.path.exists(rfile):
                if self.getStr('DataType')[:3] == 'sim':
                    return db.addTable(rfile,mainkeys=['Motif','Query','Region','RType','Rep','PosNum','SeqNum'] + self.list['RunID'] + ['Pattern'],datakeys='All',name='results',replace=True)
                else: return db.addTable(rfile,mainkeys=dsethead + self.list['RunID'] + ['Pattern'],datakeys='All',name='results',replace=True)
            ## ~ [1a] Load and merge results files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            resdb = None
            for file in self.list['ResFiles']:
                if not resdb: resdb = db.addTable(file,mainkeys=['Dataset','RunID','Pattern'],datakeys='All',name='results')
                else: db.mergeTables(resdb,db.addTable(file,mainkeys=['Dataset','RunID','Pattern'],datakeys='All',name=file))
            resdb.setStr({'Delimit':self.getStr('Delimit')})
            ## ~ [1b] Split on RunID and report stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ex = 0.0; etot = resdb.entryNum()
            resdb.list['Fields'] = dsethead + self.list['RunID'] + resdb.list['Fields']
            newkeys = dsethead + self.list['RunID'] + ['Pattern']
            if self.getStr('DataType')[:3] == 'sim': newkeys += ['SeqNum']
            for entry in resdb.entries():
                self.progLog('\r#SPLIT','Splitting Dataset and RunID fields: %.2f%%' % (ex/etot)); ex += 100.0
                if self.getStr('DataType')[:3] == 'sim':
                    try:
                        [entry['Motif'],entry['RType'],entry['Rep'],entry['PosNum'],n,entry['Query'],entry['Region']] = string.split(entry['Dataset'],'.')
                        entry['PosNum'] = string.atoi(entry['PosNum'][1:])
                    except: self.errorLog('Simulated Dataset Name Error: %s' % entry['Dataset']); raise
                elif self.getBool('Queries'):
                    try: [entry['Motif'],entry['Query'],entry['Region']] = string.split(entry['Dataset'],'.')
                    except: self.errorLog('ELMBench (Queries) Dataset Name Error: %s' % entry['Dataset']); raise
                else:
                    try: entry['Motif'] = string.split(entry['Dataset'],'.')[0]
                    except: self.errorLog('ELMBench (No Queries) Dataset Name Error: %s' % entry['Dataset']); raise
                runid = string.split(entry['RunID'],'.')
                try:
                    for i in range(len(self.list['RunID'])): entry[self.list['RunID'][i]] = runid[i]
                except IndexError: self.errorLog('Check that RunID split list and results files match!'); return False
                except: raise
            self.printLog('\r#SPLIT','Split Dataset and RunID fields into %s.' % string.join(newkeys + self.list['RunID'],', '))
            for field in ['Dataset','RunID'] + newkeys: resdb.index(field)
            resdb.newKey(newkeys,True)
            resdb.saveToFile()
            return True
        except: self.errorLog('Problem during %s loadSLiMPredictions.' % self); return False  # Setup failed
#########################################################################################################################
    def compariMotif(self):  ### Perform CompariMotif search of Patterns against CompDB                             #V1.0
        '''Perform CompariMotif search of Patterns against CompDB.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ COMPARIMOTIF ~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            cfile = '%s.compare.tdt' % self.getStr('BenchBase')
            if not self.force() and os.path.exists(cfile): return db.addTable(cfile,mainkeys=['Name1','Name2'],name='comparimotif',replace=True)
            ### ~ [1] Save Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            patterns = self.db('results').indexKeys('Pattern')
            if '-' in patterns: patterns.remove('-')
            for badrun in '<>!':
                if badrun in patterns:
                    self.printLog('#WARN','Warning: %s results runs have "%s" as Pattern. Check run completion.' % (rje.iLen(self.db('results').index('Pattern')[badrun]),badrun))
                    patterns.remove(badrun)
                    for badkey in self.db('results').index('Pattern')[badrun]: self.db('results').dict['Data'].pop(badkey)
                    self.printLog('#DROP','%s "%s" results dropped. %s remain.' % (rje.iLen(self.db('results').index('Pattern')[badrun]),badrun,rje.iStr(self.db('results').entryNum())))
            pfile = '%s.patterns.txt' % self.getStr('BenchBase')
            rje.backup(self,pfile,appendable=False)
            open(pfile,'w').write(string.join(patterns,'\n'))
            self.printLog('#MOT','%s patterns written to %s for CompariMotif search' % (rje.iLen(patterns),pfile))
            ### ~ [2] Perform CompariMotif Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.str['CompDB']): self.errorLog('CompariMotif CompDB %s not found!' % self.str['CompDB'],printerror=False); raise IOError
            self.printLog('#COMPDB','Will use %s for CompariMotif SearchDB' % self.str['CompDB'])
            comp = comparimotif.CompariMotif(self.log,['xgmml=F','normcut=0.4']+self.cmd_list+['motifs=%s' % pfile,'searchdb=%s' % self.str['CompDB'],'i=0','basefile=%s' % self.getStr('BenchBase'),'motdesc=0','motific=T','outstyle=normal'])
            comp.run()
            ### ~ [3] Read in Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return db.addTable(cfile,mainkeys=['Name1','Name2'],name='comparimotif')
        except: self.errorLog('Problem during %s compariMotif.' % self); return False  # compariMotif failed
#########################################################################################################################
    def ratings(self):  ### Compares CompariMotif search results to Datasets and rates as TP,FP or OT               #V1.0
        '''Perform CompariMotif search of Patterns against CompDB.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ RATINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            if self.getBool('Queries'): newkeys = ['Motif','Query','Region'] + self.list['RunID'] + ['Pattern']
            else: newkeys = ['Motif'] + self.list['RunID'] + ['Pattern']
            if self.getStr('DataType')[:3] == 'sim': newkeys += ['RType','Rep','PosNum','SeqNum']
            if not self.force() and os.path.exists(rfile):
                return db.addTable(rfile,mainkeys=newkeys,datakeys='All',name='results',replace=True)
            ## ~ [0a] Make reduced ratings table from results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            temp = db.copyTable('results','temp')
            temp.dropEntriesDirect('Pattern',['>','<','-','!'])
            temp.keepFields(['Motif','Pattern'])
            temp.compress(['Motif','Pattern'])
            tdb = db.joinTables(name='tp',join=[(temp,'Pattern',['Motif','Pattern']),('comparimotif','Name1')],newkey=['Motif','Pattern','Name2'],cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
            tdb.dropFields(['Name1','Motif1'])
            tdb.dataFormat({'MatchPos':'int','MatchIC':'num','NormIC':'num','Score':'num'})
            db.deleteTable(temp)
            prex = tdb.entryNum()
            #!# Output full mapping for some investigation into best settings #!#
            tdb.saveToFile('%s.cm_full.tdt' % self.getStr('BenchBase'))
            tdb.dropEntries(['MatchPos<2'],inverse=False,log=True,logtxt='Removing single position hits')   # Must have 2+ positions match
            tdb.dropEntries(['MatchIC<1.5'],inverse=False,log=True,logtxt='Removing low MatchIC hits')      # Must match 1 fixed + [3] or 2x[2]
            tdb.dropEntries(['NormIC<0.5'],inverse=False,log=True,logtxt='Removing low NormIC hits')        # Must match half smallest motif
            tdb.addField('Rating',evalue='OT')
            for entry in tdb.entries():
                #if entry['MatchIC'] == 2 and entry['NormIC'] < 1: tdb.dropEntry(entry)
                match = {}
                for field in ['Motif','Name2']:
                    if rje.matchExp('^(.+)_\d$',entry[field]): match[field] = rje.matchExp('^(.+)_\d$',entry[field])[0]
                    elif rje.matchExp('^(.+)_\d_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+)_\d_[a-z][a-z0-9]*$',entry[field])[0]
                    elif rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field])[0]
                    else: match[field] = entry[field]
                if match['Motif'] == match['Name2']: entry['Rating'] = 'TP'
                else:
                    ## Second, stronger filter for the OT motifs? 
                    if entry['MatchIC'] < 2.5 and entry['NormIC'] < 1: tdb.dropEntry(entry)
            self.printLog('#REM','Removed low quality hits: %s ratings entries reduced to %s' % (rje.iStr(prex),rje.iStr(tdb.entryNum())))
            tdb.index('Pattern')
            try: self.printLog('#TP','%s TP CompariMotif combinations' % rje.iLen(tdb.index('Rating')['TP']))
            except: self.printLog('#TP','No TP CompariMotif combinations found!')
            tdb.saveToFile()
            ### ~ [1] Rate Motifs using CompariMotif hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb = self.db('results')
            rdb.addField('Rating')
            ex = 0.0; etot = rdb.entryNum()
            for entry in rdb.entries():
                self.progLog('\r#RATE','Rating SLiM Predictions: %.1f%%' % (ex/etot)); ex += 100.0
                tkey= '%s\t%s\t%s' % (entry['Motif'],entry['Pattern'],entry['Motif'])
                if self.getStr('DataType') == 'ran' or (self.getStr('DataType') == 'sim' and entry['RType'] == 'ran'):
                    if tkey in tdb.data(): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'TN'
                    else: entry['Rating'] = 'FP'
                else:
                    if tkey in tdb.data(): entry['Rating'] = tdb.data(tkey)['Rating']
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'FN'
                    else: entry['Rating'] = 'FP'
            self.printLog('\r#RATE','Rating SLiM Predictions complete.')
            for rating in ['TP','OT','FP','FN','TN']:
                if rating in rdb.index('Rating'): self.printLog('#%s' % rating,'%s %s results' % (rje.iLen(rdb.index('Rating')[rating]),rating))
                else: self.printLog('#%s' % rating,'No %s results!' % (rating))
            ### ~ [2] Tidy up table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for field in ['Dataset','RunID','Chance','RunTime']:
                if field in rdb.fields() and (self.i() < 1 or rje.yesNo('Removed field "%s" from Ratings table?' % field)): rdb.dropField(field)
            ex = 0.0; etot = rdb.entryNum()
            rdb.addField('Coverage',after='Support'); rdb.addField('Enrichment',after='Coverage')
            rdb.addField('CloudCoverage',after='Cloud')
            for entry in rdb.entries():
                self.progLog('\r#TIDY','Tidying Ratings Table: %.1f%%' % (ex/etot)); ex += 100.0
                try: entry['Enrichment'] = rje.expectString(float(entry['UP'])/float(entry['ExpUP']))
                except: entry['Enrichment'] = '-'
                if entry['Occ']:
                    entry['CloudCoverage'] = '%s/%s (%s/%s)' % (entry['CloudSeq'],entry['SeqNum'],entry['CloudUP'],entry['UPNum'])
                    entry['Coverage'] = '%s/%s (%s/%s)' % (entry['Support'],entry['SeqNum'],entry['UP'],entry['UPNum'])
                    entry['Support'] = '%s/%s/%s' % (entry['Occ'],entry['Support'],entry['UP'])
                else:
                    entry['CloudCoverage'] = '0/%s (0/%s)' % (entry['SeqNum'],entry['UPNum'])
                    entry['Coverage'] = '0/%s (0/%s)' % (entry['SeqNum'],entry['UPNum'])
                    entry['Support'] = '0/0/0'
            self.printLog('\r#TIDY','Tidyied Ratings Table: %.1f%%' % (ex/etot),log=False)
            rdb.dropFields(['UPNum','AANum','MotNum','Occ','UP','ExpUP','Prob','CloudSeq','CloudUP'])
            rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during %s ratings.' % self); return False  # rating failed
#########################################################################################################################
    def summaryTables(self):    ### Produces summary table of best TP results                                       #V1.5
        '''Reduces results to those with minimum IC or greater and relevant SLiMLen.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SUMMARY TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            sfile = '%s.summary.tdt' % (self.getStr('BenchBase'))
            sdb = db.copyTable('results','summ_temp')
            # Motif	Query	Region	Program	Masking	Pattern	Build	SeqNum	Rank	Sig	IC	Support	Coverage	Enrichment	Cloud	CloudCoverage	IUP_mean	SA_mean	Rating
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            reduced = '%s.reduced.motifs.txt' % elmbase
            slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list)
            slimlist.loadMotifs(reduced)    # Will automatically apply minic filter
            edb = self.db('ELMClass')
            edb.addField('Reduced')
            for slim in slimlist.slims():
                if slim.getStr('Name') in edb.index('ELMIdentifier'): edb.data(slim.getStr('Name'))['Reduced'] = slim.pattern()
            ### ~ [1] Remove FP and OT motif predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.dropFields(['Build','SeqNum','Cloud','IUP_mean','SA_mean'])
            for entry in sdb.entries():
                if entry['Enrichment'] == '-': entry['Enrichment'] = 0
                if not entry['IC']: entry['IC'] = 0
            sdb.dataFormat({'Rank':'int','IC':'num','Sig':'num','Enrichment':'num'})
            sdb.addField('Q',evalue=0)
            sdb.addField('N',evalue=1)
            for entry in sdb.indexEntries('Rating','FP') + sdb.indexEntries('Rating','OT'):
                entry['Sig'] = 0.1
                entry['Pattern'] = '-'
                entry['Coverage'] = '0/%s (0/%s)' % rje.matchExp('\d+/(\d+) \(\d+/(\d+)\)',entry['Coverage'])
                entry['CloudCoverage'] = '0/%s (0/%s)' % rje.matchExp('\d+/(\d+) \(\d+/(\d+)\)',entry['CloudCoverage'])
                entry['Support'] = '0/0/0'
            for entry in sdb.indexEntries('Rating','TP'): entry['Q'] = 1
            ### ~ [2] Compress Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Queries'):
                sdb.compress(['Motif','Query','Region','Program','Masking'],rules={'Sig':'min','Enrichment':'max','Q':'max','N':'max'},default='str',best=['Sig','IC','Enrichment','Q'])
                qdb = db.copyTable(sdb,'qtemp'); sdb.dropFields(['N'])
                qdb.dropFields(['Pattern','Rank','Sig','Coverage','CloudCoverage','Support'])
                sdb.compress(['Motif','Region','Program','Masking'],rules={'Sig':'min','Enrichment':'max'},default='str',best=['Sig','IC','Enrichment'])
                qdb.compress(['Motif','Region','Program','Masking'],rules={'Sig':'min','Enrichment':'max','Q':'sum','N':'sum'},default='mean')
                for entry in qdb.entries():
                    entry['Q'] = float(entry['Q']) / entry['N']
                    sdb.data(sdb.makeKey(entry))['Q'] = entry['Q']
                sdb.dropFields(['IC','Enrichment','Rating','Query'])
            else:
                sdb.compress(['Motif','Program','Masking'],rules={'Sig':'min','Enrichment':'max'},default='str',best=['Sig','IC','Enrichment'])
                sdb.dropFields(['IC','Enrichment','Rating'])
            ### ~ [3] Reshap Wide ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.reshapeWide('Program',['Pattern','Rank','Sig','Q','Coverage','CloudCoverage','Support'])
            ### ~ [4] Add ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Queries'): sdb_newkeys = ['Motif','Region','Masking']
            else: sdb_newkeys = ['Motif','Masking']
            sdb = db.joinTables(name='summary',join=[(sdb,'Motif'),(edb,'ELMIdentifier',['Regex','Reduced'])],newkey=sdb_newkeys,cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
            sdb.dropField('summ_temp_Masking')
            self.deBug(sdb.fields())
            sfields = sdb.fields()[0:]
            for field in sdb_newkeys + ['Reduced','Regex']: sfields.remove(field)
            sdb_newkeys.reverse()
            sdb.list['Fields'] = sdb_newkeys + ['Regex','Reduced'] + sfields
            sdb.saveToFile(sfile)
            return True
        except: self.errorLog('Problem during %s ratings.' % self); return False  # rating failed
#########################################################################################################################
    def filterAmb(self,save=False):  ### Reduces results to those without ambiguity                                 #V1.5
        '''Reduces results to those with minimum IC or greater and relevant SLiMLen.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ Ambiguity FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.noamb.tdt' % (self.getStr('BenchBase'))
            rdb = self.db('results')
            rdb.dataFormat({'Rank':'int'})
            rkeys = rdb.list['Keys'][0:]
            rkeys.remove('Pattern')
            rindex = '#%s#' % string.join(rkeys,'#|#')
            ### ~ [1] Screen low MinIC motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = rdb.entryNum(); ambx = 0
            rdb.makeField(rindex)
            for dset in rdb.index(rindex):
                rentries = rdb.indexEntries(rindex,dset)
                rank1 = True; maxrank = 0
                for entry in rentries[0:]:
                    if '[' in entry['Pattern'] or '{' in entry['Pattern']: 
                        if entry['Rank'] == 1: rank1 = False
                        rentries.remove(entry)
                        if rentries: rdb.dropEntry(entry)
                        else: entry['Pattern'] = '-'; entry['Cloud'] = ''; entry['Rank'] = 0; entry['Rating'] = 'FN'
                        ambx += 1
                    maxrank = max(maxrank,entry['Rank'])
                nextrank = 1
                while not rank1 and rentries and (len(rentries) > 1 or rentries[0]['Rank'] > 0):
                    nextrank += 1
                    if nextrank > maxrank:
                        try: raise ValueError
                        except:
                            print rentries
                            self.errorLog('Cannot find new Rank 1 motif for %s' % dset); raise
                    for entry in rentries[0:]:
                        if entry['Rank'] == nextrank: entry['Rank'] = 1; rank1 = True; break
            self.printLog('#AMB','%s ambiguous patterns filtered (noamb=T).' % rje.iStr(ambx))
            if ambx: rdb.remakeKeys()
            rdb.dropField(rindex)
            if save: rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during %s filterAmb.' % self); return False  # rating failed
#########################################################################################################################
    def filterMinICLen(self,minic,slimlen,save=False):  ### Reduces results to those with minimum IC or greater     #V1.4
        '''Reduces results to those with minimum IC or greater and relevant SLiMLen.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ MinIC/SLiMLen FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.l%d.i%2f.tdt' % (self.getStr('BenchBase'),slimlen,minic)
            rdb = self.db('results')
            rdb.dataFormat({'Rank':'int','IC':'num'})
            rkeys = rdb.list['Keys'][0:]
            rkeys.remove('Pattern')
            rindex = '#%s#' % string.join(rkeys,'#|#')
            ### ~ [1] Screen low MinIC motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prex = rdb.entryNum(); dumpx = 0; crapx = 0; ambx = 0
            rdb.makeField(rindex)
            for dset in rdb.index(rindex):
                rentries = rdb.indexEntries(rindex,dset)
                rank1 = True; maxrank = 0
                for entry in rentries[0:]:
                    if entry['Pattern'] in ['>','<','!']:
                        crapx += 1; elen_ok = False; entry['Pattern'] = '-'; entry['Rating'] = 'FN'
                    #elif self.getBool('NoAmb') and ('[' in entry['Pattern'] or '{' in entry['Pattern']):
                    #    ambx += 1; elen_ok = False
                    else:
                        elen = rje_slim.slimLen(entry['Pattern'])
                        elen_ok = slimlen in [0,elen]
                    maxrank = max(maxrank,entry['Rank'])
                    if entry['Rank'] > 0 and (entry['IC'] < minic or not elen_ok):
                        if entry['Rank'] == 1: rank1 = False
                        rentries.remove(entry)
                        if rentries: rdb.dropEntry(entry)
                        else: entry['Pattern'] = '-'; entry['Cloud'] = ''; entry['Rank'] = 0; entry['Rating'] = 'FN'
                        dumpx += 1
                nextrank = 1
                while not rank1 and rentries and (len(rentries) > 1 or rentries[0]['Rank'] > 0):
                    nextrank += 1
                    if nextrank > maxrank:
                        try: raise ValueError
                        except:
                            print rentries
                            self.errorLog('Cannot find new Rank 1 motif for %s' % dset); raise
                    for entry in rentries[0:]:
                        if entry['Rank'] == nextrank: entry['Rank'] = 1; rank1 = True; break
            if ambx: self.printLog('#AMB','%s ambiguous patterns filtered (noamb=T).' % rje.iStr(ambx))
            self.printLog('#MINIC','MinResIC filtering < %s: %s filtered' % (minic,rje.iStr(dumpx-ambx)))
            if crapx: self.printLog('#BAD','%s <, >  or ! patterns modified.' % rje.iStr(crapx))
            if dumpx: rdb.remakeKeys()
            rdb.dropField(rindex)
            if save: rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during %s filterMinICLen.' % self); return False  # rating failed
#########################################################################################################################
    def assessment(self,sig_cut=0.05,cloud=False,test=False):   ### Calculates statistics for different searches    #V1.0
        '''
        Calculates statistics for different searches.
        >> sig_cut:float [0.05] = Significance cut-off for analysis
        >> cloud:bool [False] = Whether to calculate stats on clouds (True) or individual motifs (False)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~ ASSESSMENT (IC=%.2f; Len=%d) ~~~~~~~~~~~~~~~~~~~~~~~ #' % (self.num['MinResIC'],self.int['ResLen']))
            self.printLog('#SIG','Significance cut-off = %s; ByCloud = %s; Queries = %s' % (sig_cut,cloud,self.getBool('Queries')))
            db = self.obj['DB']
            ## ~ [0a] Determine headers that will identify a unique benchmarking run ~~~~~~~~~~~~~~ ##
            keep_unique = self.list['RunID'][0:]
            if self.getBool('Queries'): keep_unique += ['Region']
            keep_unique += ['ICCut','LenCut','SigCut','ByCloud']
            if self.getStr('DataType')[:3] == 'sim': keep_unique += ['Rep','PosNum','SeqNum','RType']
            mdb = db.getTable('assessment')
            bdb = db.getTable('bymotif')
            if mdb: adb = db.copyTable('results','ass_temp')
            else: mdb = adb = db.copyTable('results','assessment')
            adb.addField('ICCut',evalue=self.num['MinResIC'])
            adb.addField('LenCut',evalue=self.int['ResLen'])
            adb.addField('SigCut',evalue=sig_cut)
            adb.addField('ByCloud',evalue='%s' % cloud)
            ## ~ [0b] Significance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            prex = adb.entryNum(); dumpx = 0
            for entry in adb.entries():
                if float(entry['Sig']) > sig_cut and int(entry['Rank']) > 0:
                    if int(entry['Rank']) > 1: adb.dropEntry(entry)
                    else: entry['Cloud'] = ''; entry['Rating'] = 'FN'
                    dumpx += 1
            self.printLog('#SIG','Significance filtering > %s: %s filtered' % (sig_cut,rje.iStr(dumpx + prex - adb.entryNum())))
            ## ~ [0c] Setup fields for compression ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for rating in ['TP','OT','FP','FN','TN','SN','FPX','FPXn']: adb.addField(rating,evalue=0)
            #i# TP, OT and FP are the weighted proportions of *motifs* (or clouds)
            #i# FN and TN are the weighted proportions of datasets returning no motifs  #!# These are useless?! #!#
            #i# SN is the Sensitivity - the proportion of total TPs returned (elm or sim)
            #i# FPX is the FPR per dataset rather than by motif (%datasets returning FP). FPXn has OT=FP. (elm or ran)
            for entry in adb.entries():
                try: entry[entry['Rating']] = 1
                except: print entry; raise
            adb.keepFields(keep_unique + ['Motif','Query','SeqNum','Rank','Pattern','Cloud','TP','OT','FP','FN','TN','N','SN','FPX','FPXn'])
            if test: adb.saveToFile('test.rated.tdt')
            ### ~ [1] Optional cloud reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cloud:
                if self.getBool('Queries'): adb.compress(keep_unique + ['Motif','Query','Cloud'],default='max',rules={'SigCut':'max','ICCut':'max','LenCut':'max'})
                else: adb.compress(keep_unique + ['Motif','Cloud'],default='max',rules={'SigCut':'max','ICCut':'max','LenCut':'max'})
                for entry in adb.entries():     # Keep the best rating for cloud and rescale N to clouds, not motifs
                    if entry['TP']: entry['OT'] = entry['FP'] = entry['FN'] = 0
                    elif entry['OT']: entry['FP'] = entry['FN'] = 0
                if test: adb.saveToFile('test.cloud.tdt')
            ### ~ [2] Compress by Query (by Dataset) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Setup Dataset stats and results counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##  
            for entry in adb.entries():
                entry['SN'] = entry['TP']
                entry['FPX'] = entry['FP']
                entry['FPXn'] = max(entry['OT'],entry['FP'])
                if self.getStr('DataType') == 'sim':
                    if entry['RType'] == 'sim':     # Do not count FPX
                        entry['FPX'] = entry['FPXn'] = 0
                        if not entry['TP']: entry['FN'] = 1
                    else:                           # Do not have any FN or count TP/OT/FP for PPV
                        entry['TP'] = entry['OT'] = entry['FP'] = entry['FN'] = 0
                elif self.getStr('DataType') == 'simonly':
                    if entry['RType'] == 'sim':     
                        if not entry['TP']: entry['FN'] = 1
                    else:                           # Do not have any FN or count TP/OT/FP for PPV
                        entry['FPX'] = entry['FPXn'] = 0    # Do not count FPX
                        entry['TP'] = entry['OT'] = entry['FP'] = entry['FN'] = 0
            adb.dropField('TN')
            adb.addField('N',evalue=1)  #i# N is counting the number of results contributing to each level
            if test: adb.saveToFile('test.%s.ratingcheck.tdt' % self.getStr('DataType'))
            ## ~ [2b] Compress motifs/clouds to a single entry per dataset ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            comp_rules = {'N':'sum','SN':'max','FPX':'max','FPXn':'max'}
            self.deBug(keep_unique)
            if self.getStr('DataType')[:3] == 'sim': adb.compress(['Motif'] + keep_unique,default='mean',rules=comp_rules)
            elif self.getBool('Queries'):
                adb.compress(['Motif','Query'] + keep_unique,default='mean',rules=comp_rules)
                if test: adb.saveToFile('test.queries.tdt')
            adb.dropFields(['Pattern','Cloud','Rank'])
            adb.addField('ResNum')  #i# Total number of results (motifs/clouds)
            for entry in adb.entries():
                entry['ResNum'] = entry['N']
                entry['N'] = 1  #i# N is now the number of datasets
                # Re-normalise to proportions of results
                if self.getStr('DataType')[:3] != 'sim' or entry['RType'] == 'sim':
                    normsum = float(sum([entry['TP'],entry['OT'],entry['FP'],entry['FN']]))
                    for rating in ['TP','OT','FP','FN']: entry[rating] = entry[rating] / normsum
            if test: adb.saveToFile('test.conv.tdt')
            ### ~ [3] Compress queries/replicates to a single entry per ELM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            comp_rules = {'N':'sum','ResNum':'sum','SeqNum':'mean'}
            if self.getStr('DataType')[:3] == 'sim': keep_unique.remove('Rep')
            self.deBug(keep_unique)
            adb.compress(['Motif'] + keep_unique,default='mean',rules=comp_rules)
            if self.getStr('DataType')[:3] == 'sim':    # Combine the 'ran' and 'sim' data
                keep_unique.remove('RType'); 
                adb.compress(['Motif'] + keep_unique,default='sum')
                adb.dropField('RType'); adb.dropField('Rep')
            if 'Query' in adb.fields(): adb.dropField('Query')
            adb.addField('DsetNum')  #i# Total number of datasets (motifs/clouds)
            for entry in adb.entries():
                if self.getBool('Queries') and self.getStr('DataType') == 'elm' and entry['SeqNum'] != entry['N']: self.printLog('#ERR','MISSING DATA! %s' % adb.makeKey(entry))
                # Re-normalise to proportions of results
                normsum = float(sum([entry['TP'],entry['OT'],entry['FP'],entry['FN']]))
                if normsum:
                    for rating in ['TP','OT','FP','FN']: entry[rating] = entry[rating] / normsum
                else: self.printLog('#ERR','Rating error for: %s' % entry)
                entry['DsetNum'] = entry['N']
                entry['N'] = 1
            if test: adb.saveToFile('test.motifs.tdt')
            self.deBug(adb.entries()[0])
            if bdb:
                for entry in adb.entries(): bdb.addEntry(entry)
            else:
                bdb = db.copyTable(adb,'bymotif')
                rje.backup(self,'%s.bymotif.tdt' % self.baseFile())
            bdb.saveToFile(backup=False)
            ### ~ [4] Compress to unique keys for assessment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            comp_rules = {'N':'sum','ResNum':'sum','DsetNum':'sum','SeqNum':'mean'}
            adb.compress(keep_unique,default='mean',rules=comp_rules)
            adb.dropField('Motif')
            benchscores = ['PPV','PPVp','PPVn']
            for score in benchscores: adb.addField(score)
            for entry in adb.entries():
                if 'PPV' in benchscores:
                    entry['PPV'] = rje.ratio(entry['TP'], (entry['TP'] + entry['FP']))    # Precision
                    entry['PPVp'] = rje.ratio((entry['TP'] + entry['OT']), (entry['TP'] + entry['OT'] + entry['FP']))    # Precision
                    entry['PPVn'] = rje.ratio(entry['TP'], (entry['TP'] + entry['OT'] + entry['FP']))    # Precision
            ### ~ [5] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb.newKey(keep_unique)
            if self.getBool('MemSaver'):
                if mdb != adb: mdb.saveToFile(backup=False,append=True)
                else: mdb.saveToFile(backup=True,append=False)                
            else:
                if mdb != adb: db.mergeTables(mdb,adb)
                else: rje.backup(self,'%s.assessment.tdt' % self.baseFile())
                mdb.saveToFile(backup=False)
            return True
        except: self.errorLog('Problem during %s assessment.' % self); return False  # assessment failed
#########################################################################################################################
### End of SECTION II: SLiMBench Class                                                                                  #
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
    try: SLiMBench(mainlog,cmd_list).run()

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
