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
Version:      2.14.0
Last Edit:    08/11/17
Citation:     Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]
Copyright (C) 2012  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMBench has two primary functions:

    1. Generating SLiM prediction benchmarking datasets from ELM (or other data in a similar format). This includes
    options for generating random and/or simulated datasets for ROC analysis etc.

    2. Assessing the results of SLiM predictions against a Benchmark. This program is designed to work with SLiMFinder
    and QSLiMFinder output, so some prior results parsing may be needed for other methods.

    If `generate=F benchmark=F`, SLiMBench will check and optionally download the input files but perform no additional
    processing or analysis.

    Please see the SLiMBench manual for more details.

Commandline:
    ### ~ SOURCE DATA OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    sourcepath=PATH/    : Will look in this directory for input files if not found [SourceData/]
    sourcedate=DATE     : Source file date (YYYY-MM-DD) to preferentially use [None]
    elmclass=FILE       : Download from ELM website of ELM classes [elm_classes.tsv]
    elminstance=FILE    : Download from ELM website of ELM instances [elm_instances.tsv]
    elminteractors=FILE : Download from ELM website of ELM interactors [elm_interactions.tsv]
    elmdomains=FILE     : Download from ELM website of ELM Pfam domain interactors [elm_interaction_domains.tsv]
    elmdat=FILE         : File of downloaded UniProt entries (See rje_uniprot for more details) ['ELM.dat']
    ppisource=X         : Source of PPI data. (See documentation for details.) (HINT/FILE) ['HINT']
    ppispec=LIST        : List of PPI files/species/databases to generate PPI datasets from [HUMAN,MOUSE,DROME,YEAST]
    ppid=X              : PPI source protein identifier type (gene/uni/none; will work out from headers if None) [None]
    randsource=FILE     : Source for random/simulated dataset sequences. If species, will extract from UniProt [HUMAN]
    randat=T/F          : Whether to use DAT file for random source [False]
    download=T/F        : Whether to download files directly from websites where possible if missing [True]
    integrity=T/F       : Whether to quit by default if source data integrity is breached [False]
    unipath=PATH        : Path to UniProt download. Will query website if "URL" [URL]

    ### ~ GENERAL/ELM BENCHMARK GENERATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genpath=PATH        : Output path for datasets generated with SLiMBench file generator [./SLiMBenchDatasets/]
    generate=T/F        : Whether to generate SLiMBench datasets from ELM input. [False]
    genspec=LIST        : Restrict ELM/OccBench datasets to listed species (restricts ELM instances) []
    slimmaker=T/F       : Whether to use SLiMMaker to "reduce" ELMs to more findable SLiMs [True]
    minupc=X            : Minimum number of UPC for benchmark dataset [3]
    maxseq=X            : Maximum number of sequences for benchmark datasets [0]
    minic=X             : Min information content for a motif (1 fixed position = 1.0) [2.0; 1.1 for OccBench]
    filterdir=X         : Directory suffix for filtered benchmarking datasets [_Filtered/]
    queries=T/F         : Whether to generate datasets with specific Query proteins [False]
    flankmask=LIST      : List of flanking mask options (used with queries and simbench) [none,win100,flank5,site]
    elmbench=T/F        : Whether to generate ELM datasets [True]
    ppibench=T/F        : Whether to generate ELM PPI datasets [True]
    domlink=T/F         : Link ELMs to PPI via Pfam domains (True) or (False) just use direct protein links [True]
    itype=X             : Interaction identifer for PPI datasets [first element of ppisource]
    dombench=T/F        : Whether to generate Pfam domain ELM PPI datasets [True]
    occbench=T/F        : Whether to generate ELM OccBench datasets [True]
    
    ### ~ RANDOM/SIMULATION BENCHMARK GENERATION OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    simbench=T/F        : Whether to generate simulated datasets using reduced ELMs (if found) [False]
    ranbench=T/F        : Whether to generate randomised datasets (part of simulation if simbench=T) [False]
    randreps=X          : Number of replicates for each random (or simulated) datasets [8]
    simcount=LIST       : Number of "TPs" to have in dataset [4,8,16]
    simratios=LIST      : List of simulated ELM:Random ratios [0,1,3,7,15,31]
    randir=PATH         : Output path for creation of randomised datasets [./SLiMBenchDatasets/Random/]
    randbase=X          : Base for random dataset name if simbench=F [ran]
    masking=T/F         : Whether to use SLiMCore masking for query selection [True]
    searchini=FILE      : INI file containing SLiMProb search options that restrict returned positives []
    maxseq=X            : Maximum number of randsource sequences for SLiM to hit (also maxaa and maxupc limits) [1000]

    ### ~ BENCHMARK ASSESSMENT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    benchmark=T/F       : Whether to perfrom SLiMBench benchmarking assessment against motif file [False]
    datatype=X          : Type of data to be generated and/or benchmarked (occ/elm/ppi/dom/sim/simonly) [elm]
    queries=T/F         : Whether to datasets have specific Query proteins [False]
    resfiles=LIST       : List of (Q)SLiMFinder results files to use for benchmarking [*.csv]
    balanced=T/F        : Whether to reduce benchmarking to datasets found for all RunIDs [True]
    compdb=FILE         : Motif file to be used for benchmarking [elmclass file] (reduced unless occ/ppi)
    occbenchpos=FILE    : File of all positive occurrences for OccBench [genpath/ELM_OccBench/ELM.full.ratings.csv]
    benchbase=X         : Basefile for SLiMBench benchmarking output [slimbench]
    runid=LIST          : List of factors to split RunID column into (on '.') [Program,Analysis]
    bycloud=X           : Whether to compress results into clouds prior to assessment (True/False/Both) [Both]
    sigcut=LIST         : Significance thresholds to use for assessment [0.1,0.05,0.01,0.001,0.0001]
    iccut=LIST          : Minimum IC for (Q)SLiMFinder results for elm/sim/ppi benchmark assessment [2.0,2.1,3.0]
    slimlencut=LIST     : List of individual SLiM lengths to return results for (0=All) [0,3,4,5]
    noamb=T/F           : Filter out ambiguous patterns [False]

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
import rje, rje_db, rje_obj, rje_ppi, rje_seq, rje_seqlist, rje_slim, rje_slimcore, rje_slimlist, rje_uniprot
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
    # 1.9 - Removed 3DID again: new ELM interaction_domains file has position-specific PPI details.
    # 2.0 - Major overhaul of input options to standardise/clarify. Implemented auto-downloads and PPI datasets.
    # 2.1 - Fixed memsaver=T unless in development mode (dev=T). Removed old Assessment. Tested with simbench analysis.
    # 2.2 - Replaced searchini=LIST with searchini=FILE and moved to SimBench commands.
    # 2.2 - Modified the FN/TN and ResNum calculations. No longer rate TP in random data as OT.
    # 2.3 - Changed the default to queries=F. SearchINI bug fix. Added occbench generation.
    # 2.4 - Improved error messages.
    # 2.5 - Basic OccBench assessment benchmarking. Added ELM Uniprot acclist output. (Download issues?)
    # 2.6 - Added ELM domain interactions table: http://www.elm.eu.org/infos/browse_elm_interactiondomains.tsv.
    # 2.6 - Fixed issues introduced with new SLiMCore V2.0 SLiMSuite code.
    # 2.7 - Reinstate filtering. (Not sure why disabled.) Add genspec=LIST to filter by species. Added domlink=T/F.
    # 2.8.0 - Implemented PPIBench benchmarking for datasets without Motifs in name.
    # 2.8.1 - Removed use of Protein name for ELM Uniprot entries due to problems mapping old IDs.
    # 2.9.0 - Added SLiMMaker ELM reduction table and output.
    # 2.9.1 - Enabled download only with generate=F benchmark=F.
    # 2.10.0 - Add generation of table mapping PPIBench dataset generation.
    # 2.10.1 - Updated ELM Source URLs.
    # 2.10.2 - Updated HINT Source URLs.
    # 2.11.0 - Fixed issue with ELM motifs file names (*.motifs, not *.motifs.txt). Updated some warning/error messages.
    # 2.11.1 - Switched mergesplits=F for SLiMProb run. (Not expecting it.)
    # 2.11.2 - Trying to complete implementation of PPIBench.
    # 2.11.3 - Fixed ppdb bug when making simbench without ppibench.
    # 2.12.0 - Added randat=T/F : Whether to use DAT file for random source [False]
    # 2.13.0 - Added balanced=T/F : Whether to reduce benchmarking to datasets found for all RunIDs [True]
    # 2.13.1 - Set tuplekeys=T for benchmark assessment runs.
    # 2.13.2 - Fixed tuplekeys=T bug for benchmark assessment runs.
    # 2.14.0 - Added wPPV = SN/(SN+FPR) for OccBench.
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
    # [X] : Add ELMdb.py script to download data directly from ELM.
    # [?] : Consider how to handle ELMs that break rje_slimlist etc.
    # [X] : Add 3DID benchmarking too.
    # [?] : Check that simulated data can return the query motif. (SLiMSearch)
    # [Y] : Change summary to have maximum support values etc.
    # [Y] : Improve RunID and Dataset split error messages.
    # [X] : Add option to split Dataset into different things. No: add other datatypes as needed. (e.g. w/o query)
    # [?] : Try using rje_slimcore motifSeq() (in waves) as part of simulation process to check signal.
    # [ ] : Make sure that PPI and 3DID dataset generation is complete and test.
    # [?] : Check for existence of motif splits file and suggest replacement. (Or use if force=F.)
    # [Y] : Test with memsaver=T.
    # [Y] : Test with SLiMSearch -> SLiMProb replacement.
    # [Y] : Output equivalence file from reduced runs to allow "fair" predictions.
    # [Y] : Add retrieval of UniProt entries directly from UniProt rest.
    # [Y] : Fix issues of failing PPI parsing causing issues later. Check and give option to switch off species.
    # [Y] : Add dombench option through PFam parsing.
    # [Y] : Solve problem of duplicate sequences in PPI (and presumably DomPPI) datasets
    # [ ] : Add an option to benchmark by Rank rather than Significance. (Not all methods have Sig.)
    # [ ] : Add generation of PPI datasets using both Reduced and Full ELM motifs.
    # [Y] : Need to remove the searchini=LIST option. (Replace with searchini=FILE.)
    # [ ] : Add option to generate DAT or AccList files for benchmarking rather than fasta files.
    # [ ] : Build an occurrence-based benchmarking for SLiMSearch-type analysis and also SLiMFinder coverage.
    # [Y] : Add occbench generation - generate fasta files for full ELM and ELM reduced.
    # [ ] : Add occbench assessment - will need to be given search files and motifs for background. Needs reference results for TN.
    # [ ] : SLiMMaker: Consider splitting ELM instances by length and then trying to recombined with variable wildcards?
    # [?] : Add output of graphs etc. for assessment results.
    # [ ] : Add summary table field options for splitting.
    # [ ] : Add an over-ride for dataset name splitting.
    # [Y] : Remove .ppi and .dpi from PPIBench dataset names. (Have code that accepts with & without.)
    # [Y] : Add new ELM domain interactions table: http://www.elm.eu.org/infos/browse_elm_interactiondomains.tsv = "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
    # [ ] : Consider adding Domain Name to DPI files (in addition to Pfam ID).
    # [ ] : Add CompariMotif settings for OT/TP etc. - How will they affect benchmarking?
    # [ ] : Update to handle run with mergesplits=T.
    # [ ] : Fixed ppi benchmark in non-Memsaver mode.
    # [ ] : Check domppi benchmark=T.
    # [ ] : Check SimBench generation with queries=F and region settings.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('SLiMBench', '2.14.0', 'November 2017', '2012')
    description = 'Short Linear Motif prediction Benchmarking'
    author = 'Dr Richard J. Edwards.'
    comments = ['Cite: Palopoli N, Lythgow KT & Edwards RJ. Bioinformatics 2015; doi: 10.1093/bioinformatics/btv155 [PMID: 25792551]',
                'Please report bugs to Richard.Edwards@UNSW.edu.au']
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
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
            if rje.yesNo('Show SLiMMaker commandline options?','N'): out.verbose(-1,4,text=slimmaker.__doc__)
            if rje.yesNo('Show SLiMProb commandline options?','N'): out.verbose(-1,4,text=slimprob.__doc__)
            if rje.yesNo('Show rje_slimlist commandline options?','N'): out.verbose(-1,4,text=rje_slimlist.__doc__)
            if rje.yesNo('Show rje_slimcore commandline options?','N'): out.verbose(-1,4,text=rje_slimcore.__doc__)
            if rje.yesNo('Show RJE_UniProt commandline options?','N'): out.verbose(-1,4,text=rje_uniprot.__doc__)
            if rje.yesNo('Show general commandline options?','N'): out.verbose(-1,4,text=rje.__doc__)
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
    - ELMDat = File of downloaded UniProt entries (See rje_uniprot for more details) ['ELM.dat']
    - ELMClass = Download from ELM website of ELM classes ['elm_classes.tsv']
    - ELMInstance = Download from ELM website of ELM instances ['elm_instances.tsv']
    - ELMInteractors = Download from ELM website of ELM interactors ['elm_interaction.tsv']
    - ELMDomains = Download from ELM website of ELM Pfam domain interactors ['elm_interaction_domains.tsv']
    - FilterDir = Directory suffix for filtered benchmarking datasets [_Filtered]
    - GenPath = Output path for datasets generated with SLiMBench file generator [./SLiMBenchDatasets/]
    - IType = Interaction identifer for PPI datasets [first element of ppisource]
    - OccBenchPos = File of all positive occurrences for OccBench [genpath/ELM_OccBench/ELM.full.occ.csv]
    - PPID = PPI source protein identifier type (gene/uni/none; will work out from headers if None) [None]
    - PPISource = Source of PPI data. (See documentation for details.) (HINT/FILE) ['HINT']
    - RanDir = Output path for creation of randomised datasets [SLiMBenchDatasets/Random/]
    - RandBase = Base for random dataset name [rand]
    - RandSource = Source for new sequences for random datasets [None]
    - SearchINI = INI file containing search options (should have runid setting) []
    - SourceDate = Source file date (YYYY-MM-DD) to preferentially use [None]
    - SourcePath = Will look in this directory for input files if not found ['SourceData/']

    Bool:boolean
    - Balanced=T/F        : Whether to reduce benchmarking to datasets found for all RunIDs [True]
    - Benchmark = Whether to perfrom SLiMBench benchmarking assessment [False]
    - DomBench = Whether to generate Pfam domain ELM PPI datasets [True]
    - DomLink = Link ELMs to PPI via Pfam domains (True) or (False) just use direct protein links [True]
    - Download = Whether to download files directly from websites where possible if missing [True]
    - ELMBench = Whether to generate ELM datasets [True]
    - Generate = Whether to generate SLiMBench datasets from ELM input [False]
    - Integrity = Whether to quit by default if input integrity is breached [False]
    - Masking = Whether to use SLiMCore masking for query selection [True]
    - NoAmb = Filter out ambiguous patterns [False]
    - OccBench = Whether to generate ELM OccBench datasets [True]
    - PPIBench = Whether to generate ELM PPI datasets [True]
    - Queries = Whether to generate datasets with specific Query proteins [True]
    - RanBench = Whether to generate randomised datasets (part of simulation if SimBench=T) [False]
    - RanDat=T/F : Whether to use DAT file for random source [False]
    - SimBench = Whether to generate simulated datasets using reduced ELMs (if found) [False]
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
    - FlankMask = List of flanking mask options [none,win100,flank5,site]
    - GenSpec = Restrict ELM/OccBench datasets to listed species (filters ELM instances) []
    - ICCut = Minimum IC for (Q)SLiMFinder results for benchmark assessment [2.0,2.1,3.0]
    - PPISpec = List of PPI files/species/databases to generate PPI datasets from ['HUMAN','YEAST']        
    - ResFiles = List of (Q)SLiMFinder results files to use for benchmarking [*.csv]
    - RunID = List of factors to split RunID column into (on '.') ['Program','Settings']
    - SimCount = Number of "TPs" to have in dataset [4,8,16]
    - SimRatios = List of simulated ELM:Random rations [0,1,3,7,15,31]
    - SigCut = Significance thresholds to use for assessment [0.1,0.05,0.01,0.001,0.0001]
    - SLiMLenCut = List of individual SLiM lengths to return results for (0=All) [0,3,4,5]

    Dict:dictionary
    - PPI = dictionary of {Species:PPI object}?
    - UniProt = dictionary of {AccNum:UniProtEntry} matching UniProt field added to ELM Instance table
    - UniSpec = dictionary of {Species:UniProt Object} if memsaver=F

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
        self.strlist = ['BenchBase','ByCloud','CompDB','DataType','DMIData','ELMClass','ELMDat','ELMInstance','IType',
                        'ELMInteractors','FilterDir','GenPath','OccBenchPos','PPISource','PPID','ELMDomains',
                        'RanDir','RandBase','RandSource','SourcePath','SourceDate','SearchINI']
        self.boollist = ['Balanced','Benchmark','DomBench','Download','ELMBench','OccBench','PPIBench','Generate','Integrity','Masking',
                         'NoAmb','Queries','RanBench','RanDat','SimBench','SLiMMaker','DomLink']
        self.intlist = ['MinUPC','RandReps']
        self.numlist = ['MinIC','MinResIC']
        self.listlist = ['GenSpec','ByCloud','FlankMask','ICCut','SimCount','ResFiles','RunID','SimRatios','SigCut','SLiMLenCut','PPISpec']
        self.dictlist = ['UniProt','UniSpec']
        self.objlist = ['DB','Uniprot']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'BenchBase':'slimbench','ByCloud':'Both','DataType':'elm',
                     'ELMClass':'elm_classes.tsv','ELMInstance':'elm_instances.tsv','ELMDat':'ELM.dat',
                     'ELMInteractors':'elm_interactions.tsv','DMIData':'3did.DMI.csv',
                     'ELMDomains':'elm_interaction_domains.tsv','FilterDir':'_Filtered',
                     'GenPath':rje.makePath('./SLiMBenchDatasets/'),
                     'PPISource':'HINT','SourceDate':'','RandSource':'HUMAN',
                     'RanDir':rje.makePath('./Random/'),'RandBase':'rand',
                     'SourcePath':rje.makePath('./SourceData/')})
        self.setBool({'Generate':False,'Integrity':False,'Masking':True,'NoAmb':False,'Queries':False,'SLiMMaker':True,
                      'ELMBench':True,'PPIBench':True,'Download':True,'DomBench':True,'OccBench':True,'DomLink':True,
                      'Balanced':True})
        self.setInt({'MinUPC':3,'RandReps':8})
        self.setNum({'MinIC':2.0,'MinResIC':2.1})
        self.cmd_list = ['minic=2.0','minupc=3'] + self.cmd_list    # Propagate defaults
        self.list['FlankMask'] = string.split('site,flank5,win100,none',',')
        self.list['PPISpec'] = ['HUMAN','MOUSE','DROME','YEAST']
        self.list['RunID'] = ['Program','Analysis']
        self.list['SimRatios'] = [0,1,3,7,15,31]
        self.list['SimCount'] = [4,8,16]
        self.list['SigCut'] = [0.1,0.05,0.01,0.001,0.0001]
        self.list['SLiMLenCut'] = [0,3,4,5]
        self.list['ICCut'] = [2.0,2.1,3.0]
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.cmd_list.insert(0,'basefile=slimbench')
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
        self.obj['SLiMMaker'] = slimmaker.SLiMMaker(self.log,self.cmd_list)
        self.obj['UniProt'] = rje_uniprot.UniProt(self.log,['unipath=url']+self.cmd_list)
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
                self._cmdReadList(cmd,'str',['ByCloud','DataType','FilterDir','IType','RandBase','PPID'])   # Normal strings
                self._cmdReadList(cmd,'date',['SourceDate'])   # Normal strings
                self._cmdReadList(cmd,'path',['GenPath','RanDir','SourcePath'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['BenchBase','CompDB','ELMClass','ELMInstance','ELMInteractors','ELMDomains','SearchINI','RandSource','ELMDat','OccBenchPos','PPISource'])  # String representing file path
                self._cmdReadList(cmd,'bool',['Benchmark','DomBench','DomLink','Download','ELMBench','OccBench','PPIBench','Generate','Integrity','Masking','NoAmb','Queries','RanBench','RanDat','SimBench','SLiMMaker'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinUPC','RandReps'])   # Integers
                self._cmdReadList(cmd,'float',['MinIC','MinResIC']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['GenSpec','RunID','PPISpec'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['SimRatios','SimCount','SLiMLenCut'])
                self._cmdReadList(cmd,'nlist',['SigCut','ICCut'])  
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['ResFiles']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.str['DataType'] = self.str['DataType'].lower()  
        self.list['ByCloud'] = []
        if self.getStr('RanDir') != os.path.abspath(self.getStr('RanDir')):
            if not os.path.abspath(self.getStr('RanDir')).startswith(os.path.abspath(self.getStr('GenPath'))):
                self.setStr({'RanDir':'%s%s' % (self.getStr('GenPath'),self.getStr('RanDir'))})
        if self.getStrLC('ByCloud') in ['t','true','both']: self.list['ByCloud'].append(True)
        elif self.getStrLC('ByCloud') in ['f','false','both']: self.list['ByCloud'].append(False)
        else: self.warnLog('ByCloud "%s" not recognised. Will set to "Both".' % self.getStr('ByCloud')); self.list['ByCloud'] = [True,False]
        if self.getStrLC('SourceDate') in ['none','today']: self.setStr({'SourceDate':rje.dateTime(dateonly=True)})
        if not self.getStrLC('OccBenchPos'): self.setStr({'OccBenchPos':'%sELM_OccBench/ELM.full.ratings.csv' % self.getStr('GenPath')})
        if not self.dev(): self.setBool({'MemSaver':True})
        if not self.getStrLC('IType'): self.setStr({'IType':string.split(rje.baseFile(self.getStrLC('PPISource'),strip_path=True),'.')[0]})
        if self.getBool('Benchmark') and not self.list['ResFiles']: self._cmdRead('resfiles=*.csv',type='glist',att='ResFiles')
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method                                                                             #V2.0
        '''Main run method.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Generate') and self.getBool('Benchmark'):
                self.printLog('#BENCH','Cannot benchmark if generate=T. Switching benchmark=F.')
            elif not self.getBool('Generate') and not self.getBool('Benchmark'):
                self.warnLog('Data download/check only: generate=F benchmark=F.')
                if rje.yesNo('Switch on SLiMBench dataset generation (generate=T)?','N'): self.setBool({'Generate':True})
                elif rje.yesNo('Switch on SLiMBench benchmarking (benchmark=T)?','N'): self.setBool({'Benchmark':True})
                else: return self.setupSourceData()
            ### ~ [1] ~ Generate Benchmark Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Generate'): return self.generate()    
            ### ~ [2] ~ Benchmark motif predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Benchmark'): return self.benchmark()
            ### ~ [3] ~ Pointless run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RUN','No SLiMBench run if generate=F benchmark=F.'); return False
        except KeyboardInterrupt: raise
        except: self.errorLog('SLiMBench has encountered a terminal error'); sys.exit(1)
#########################################################################################################################
    def sourceDataFile(self,str,ask=True,expect=True,check=True,force=True):   ### Returns source data file.        #V2.0
        '''
        Returns source data file.
        >> str = Key for self.str dictionary corresponding to the source data file being sought.
        >> ask:bool [True] = Whether to ask for file if not found.
        >> expect:bool [True] = Whether to download/ask for file name if missing
        >> check:bool [True] = Whether to check for file presence. If False will return the desired file download name.
        >> force:bool [True] = Whether to use self.force() to govern file regeneration.
        << returns filename if present or None if missing.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            force = force and self.force()
            ask = ask and self.i() >= 0
            if not self.db('Source'): self.db().addEmptyTable('Source',['Name','File','Status','Entries'],keys=['Name'],log=False)   # Store Source info
            sentry = self.db('Source').data(str) 
            if not sentry: sentry = self.db('Source').addEntry({'Name':str})
            ## ~ [0a] Filename lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            filename = os.path.split(self.getStr(str))[-1]
            #self.deBug(filename)
            sourcefile = self.getStr('SourcePath') + filename
            fileparts = os.path.splitext(filename)
            if fileparts[0].endswith(self.getStr('SourceDate')):
                datefile = sourcefile
            else:
                if rje.matchExp('^(.+)\.(\d\d\d\d-\d\d-\d\d)$',fileparts[0]):
                    fileparts = (rje.matchExp('^(.+)\.(\d\d\d\d-\d\d-\d\d)$',fileparts[0])[0],fileparts[1])
                datefile = '%s%s.%s%s' % (self.getStr('SourcePath'),fileparts[0],self.getStr('SourceDate'),fileparts[1])
            nowfile = '%s%s.%s%s' % (self.getStr('SourcePath'),fileparts[0],rje.dateTime(dateonly=True),fileparts[1])
            if not check: return nowfile
            lastfile = None
            for gfile in glob.glob('%s%s.*%s' % (self.getStr('SourcePath'),fileparts[0],fileparts[1])):
                if rje.matchExp('^%s%s\.(\d\d\d\d)-(\d\d)-(\d\d)%s$' % (self.getStr('SourcePath'),fileparts[0],fileparts[1]),gfile):
                    lastfile = gfile    # Possibly use the latest dated version
            if lastfile in [self.getStr(str), datefile, sourcefile, nowfile]: lastfile = None   # Only interested if different
            ## ~ [0b] Source URLs for downloads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sourceurl = {'ELMClass':'http://www.elm.eu.org/elms/elms_index.tsv',
                         'ELMInstance':'http://www.elm.eu.org/instances.tsv?q=*',
                         'ELMInteractors':'http://www.elm.eu.org/interactions/as_tsv',
                         'ELMDomains':'http://www.elm.eu.org/infos/browse_elm_interactiondomains.tsv',
                         'HINT.HUMAN':'http://hint.yulab.org/download/HomoSapiens/binary/hq/',
                         'HINT.YEAST':'http://hint.yulab.org/download/SaccharomycesCerevisiaeS288C/binary/hq/',
                         'HINT.MOUSE':'http://hint.yulab.org/download/MusMusculus/binary/hq/',
                         'HINT.DROME':'http://hint.yulab.org/download/DrosophilaMelanogaster/binary/hq/',
                         'HINT.CAEEL':'http://hint.yulab.org/download/CaenorhabditisElegans/binary/hq/',
                         'Uniprot.HUMAN':'http://www.uniprot.org/uniprot/?query=organism:9606+AND+reviewed:yes&format=txt',
                         'Uniprot.YEAST':'http://www.uniprot.org/uniprot/?query=organism:559292+AND+reviewed:yes&format=txt',
                         'Uniprot.MOUSE':'http://www.uniprot.org/uniprot/?query=organism:10090+AND+reviewed:yes&format=txt',
                         'Uniprot.DROME':'http://www.uniprot.org/uniprot/?query=organism:7227+AND+reviewed:yes&format=txt',
                         'Uniprot.CAEEL':'http://www.uniprot.org/uniprot/?query=organism:6239+AND+reviewed:yes&format=txt'
                         }

            ### ~ [1] Return file if it exists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            checked = [None]
            for checkfile in [self.getStr(str), datefile, sourcefile, nowfile, lastfile]:
                if checkfile in checked: continue
                checked.append(checkfile)
                if not rje.exists(checkfile): continue
                if checkfile != nowfile and force: self.printLog('#FORCE','Ignoring %s. (force=T)' % checkfile); continue
                if checkfile == lastfile and not rje.yesNo('%s not found. Use %s?' % (sourcefile,lastfile)): continue
                if checkfile not in [self.getStr(str),datefile,nowfile]:
                    if not self.getStr('SourceDate'):
                        newdate = rje.matchExp('^(\d\d\d\d-\d\d-\d\d)$',lastfile.split('.')[-2])[0]
                        if rje.yesNo('Set sourcedate=%s?' % newdate): self.setStr({'SourceDate':newdate})
                        else:
                            self.setStr({'SourceDate':rje.dateTime(dateonly=True)})
                            self.warnLog('Using %s rather than %s (sourcedate=%s)' % (checkfile,datefile,self.getStr('SourceDate')))
                    else: self.warnLog('Using %s rather than %s (sourcedate=%s)' % (checkfile,datefile,self.getStr('SourceDate')))
                if force and rje.yesNo('%s found but force=T. Regenerate?' % checkfile): self.printLog('#FORCE','Ignoring %s. (force=T)' % checkfile); continue
                sentry['Status'] = 'Found'
                sentry['File'] = checkfile
                if checkfile != self.getStr(str): self.printLog('#SOURCE','Set %s=%s.' % (str.lower(),checkfile))
                self.setStr({str:checkfile})
                return checkfile
            if not expect: return None

            ### ~ [2] Optional download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x# Special fudge for DROME no longer required. NOTE: HINT headers have now changed.
            #if self.getBool('Download') and str in ['HINT.DROME']:
            #    self.warnLog('Trying to correct for dodgy HINT.DROME download w/o headers')
            #    ## Special fudge
            #    open(nowfile,'w').write('%s\n' % string.join('Id_A Id_B Gene_A Gene_B Pubmedid,EvidenceCode,HT'.split(),'\t'))
            #    sourcefile = nowfile
            #    rje.urlToFile(sourceurl[str],nowfile,self,backupfile=False)
            #    sentry['Status'] = 'Download'
            if self.getBool('Download') and str in sourceurl:
                sourcefile = nowfile
                rje.urlToFile(sourceurl[str],nowfile,self)
                sentry['Status'] = 'Download'
            elif self.getBool('Download') and expect: self.warnLog('Downloads for %s not supported' % str)

            ### ~ [3] Ask for replacement file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while sourcefile and ask and not rje.exists(sourcefile):
                btext = {True:'exit',False:'ignore'}[expect and self.getBool('Integrity')]
                sourcefile = rje.choice('%s not found. Enter %s file path (blank to %s)' % (sourcefile,str,btext))
                sentry['Status'] = 'Manual'
            sentry['File'] = sourcefile
            if rje.exists(sourcefile):
                if sourcefile != self.getStr(str): self.printLog('#SET','Set %s=%s.' % (str.lower(),sourcefile))
                self.setStr({str:sourcefile})
                return sourcefile
            elif expect and sourcefile: raise IOError
            elif expect and self.getBool('Integrity'): raise KeyboardInterrupt
            self.warnLog('Problem locating %s source data file' % str,quitchoice=self.getBool('Integrity'))
            sentry['Status'] = 'Missing'
            return ''
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem locating %s source data file' % str,quitchoice=self.getBool('Integrity')); return ''
#########################################################################################################################
    def setupSourceData(self):    ### Main class setup method.                                                      #V2.0
        '''Main class setup method.'''
        try:### ~ [0] Setup Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            sdb = db.addEmptyTable('Source',['Name','File','Status','Entries'],keys=['Name'],log=False)   # Store Source info
            rje.mkDir(self,self.getStr('SourcePath'),True)
            db.baseFile('%sslimbench' % self.getStr('SourcePath'))
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP SOURCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            
            ### ~ [1] Load ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# At some point, add option to download and try again if problem with file #!#
            ## ~ [1a] Load ELM Classes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sentry = sdb.addEntry({'Name':'ELMClass'})
            elmc = db.addTable(self.sourceDataFile('ELMClass'),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
            elmc.dataFormat({'#Instances':'int'})
            sentry['Entries'] = elmc.entryNum()
            if not sentry['Entries']: self.warnLog('No %s entries loaded from %s!' % (sentry['Name'],sentry['File']),quitchoice=True)
            ## ~ [1b] Load ELM Instances ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sentry = sdb.addEntry({'Name':'ELMInstance'})
            elmi = db.addTable(self.sourceDataFile('ELMInstance'),mainkeys=['ELMIdentifier','Primary_Acc','Start','End'],datakeys='All',name='ELMInstance')
            elmi.dataFormat({'Start':'int','End':'int'})
            sentry['Entries'] = elmi.entryNum()     # Note that Primary_Acc can include splice variants
            if not sentry['Entries']: self.warnLog('No %s entries loaded from %s!' % (sentry['Name'],sentry['File']),quitchoice=True)
            ## ~ [1c] Load ELM interactors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sentry = sdb.addEntry({'Name':'ELMInteractors'})
            elmp = db.addTable(self.sourceDataFile('ELMInteractors'),mainkeys='#',datakeys='All',name='ELMInteractors')
            try:
                for entry in elmp.entries():
                    if entry['StartDomain'] == 'None': entry['StartDomain'] = '0'   # Not all domain entries have positions
                    if entry['StopDomain'] == 'None': entry['StopDomain'] = '0'
                    if entry['taxonomyElm'].find('\t'): [entry['taxonomyElm'],entry['taxonomyDomain']] = string.split(entry['taxonomyElm'],'\t')
            except: self.errorLog('Something went wrong trying to fix dodgy %s format' % sentry['File'])
            elmp.dataFormat({'StartElm':'int','StopElm':'int','StartDomain':'int','StopDomain':'int'})
            sentry['Entries'] = elmp.entryNum()
            if not sentry['Entries']: self.warnLog('No %s entries loaded from %s!' % (sentry['Name'],sentry['File']),quitchoice=True)
            sentry = sdb.addEntry({'Name':'ELMDomains'})
            elmd = db.addTable(self.sourceDataFile('ELMDomains'),mainkeys='#',datakeys='All',name='ELMDomains') # "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
            elmd.renameField('ELM identifier','Elm'); elmd.renameField('Interaction Domain Id','Domain')
            sentry['Entries'] = elmd.entryNum()
            if not sentry['Entries']: self.warnLog('No %s entries loaded from %s!' % (sentry['Name'],sentry['File']),quitchoice=True)
            ## ~ [1d] Check and summarise loaded data integrity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmdates = []
            for entry in sdb.entries(): 
                try: elmdates.append(rje.matchExp('\.(\d\d\d\d-\d\d-\d\d)\.',entry['File'])[0])
                except: elmdates.append('None')
            if not elmdates.count(elmdates[0]) == len(elmdates): self.warnLog('Date mismatch in ELM files: %s' % string.join(elmdates,'; '),quitchoice=True)
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
            sentry = sdb.addEntry({'Name':'ELMDat'})
            uni = self.obj['UniProt']
            ufile = self.sourceDataFile('ELMDat',expect=False,force=False)
            if ufile: uni.setInfo({'Name':self.getStr('ELMDat'),'UniPath':self.getStr('SourcePath')})
            else: uni.setInfo({'DATOut':self.sourceDataFile('ELMDat',check=False)})
            elmacc = elmi.indexKeys('Primary_Acc') #+ elmi.indexKeys('ProteinName')
            while '' in elmacc: elmacc.remove('')
            if self.force() and ufile and rje.yesNo('%s found but force=T. Regenerate?' % ufile):
                self.printLog('#FORCE','Ignoring %s. (force=T)' % ufile); ufile = None
            if not ufile:
                accfile = '%s.acc' % rje.baseFile(uni.getStr('DATOut'))
                open(accfile,'w').write('%s\n' % string.join(elmacc,'\n'))
                self.printLog('#ELMACC','%s ELM Uniprot accession numbers output to %s for manual download in case of failure.' % (rje.iLen(elmacc),accfile))
            uni.readUniProt(acclist=elmacc,cleardata=False)         # Extract entries for ELM proteins
            elmaccdict = uni.accDictFromEntries(acc_list=elmacc)    # Generates dictionary of {acc:UniProtEntry}
            self.sourceDataFile('ELMDat',expect=True,ask=False,force=False)     # This should perform integrity check
            sentry['Entries'] = len(elmaccdict)
            ## ~ [2a] Generate ELM UniProt mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmi.addField('UniProt'); missing = []
            unimap = self.dict['UniProt'] = {}
            unifound = 0; unimissed = 0
            for entry in elmi.entries():
                #for acc in [entry['ProteinName'],entry['Primary_Acc']][0:]:
                acc = string.split(entry['Primary_Acc'],'-')[0]
                if acc in elmaccdict:
                    u_entry = elmaccdict[acc]
                    u_name = u_entry.info['Name']
                    entry['UniProt'] = u_name
                    unimap[u_name] = u_entry
                if entry['UniProt']: unifound += 1
                else:
                    #missing += [entry['ProteinName'],entry['Primary_Acc']]
                    missing += [entry['Primary_Acc']]
                    self.printLog('#ERR','Could not find UniProt entry %s (%s)' % (entry['Primary_Acc'],entry['ProteinName']))
                    unimissed += 1
            self.printLog('#UNI','Mapped %d of %d instances onto UniProt entries' % (unifound,elmi.entryNum()))
            if unimissed: self.printLog('#MISS','UniProt sequences missing for %d entries!' % unimissed)
            if missing:
                if self.getBool('Download'):
                    missfile = '%s.missing.acc' % rje.baseFile(self.getStr('ELMInstance'))
                    #rje.backup(self,missfile)
                    open(missfile,'w').write(string.join(missing,'\n'))
                    self.printLog('#MISS','%d missing UniProt ID/AccNum output to %s' % (len(missing),missfile))
                if (self.getBool('Integrity') and (self.i() < 0 or rje.yesNo('ELM UniProt entries missing. Quit? (Switch integrity=F to avoid question in future.)'))) or (not self.getBool('Integrity') and self.i() > 1 and rje.yesNo('ELM UniProt entries missing. Quit?',default='N')):
                    self.printLog('#ERR','Setup aborted due to ELM file UniProt errors')
                    return False
                else: self.warnLog('ELM UniProt entries missing.')
            ## ~ [2b] Populated Instance field ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmi.addField('SeqAcc')
            elmi.addField('Instance')
            dx = 0
            for elm in elmc.dataKeys():
                for entry in elmi.indexEntries('ELMIdentifier',elm):
                    entry['Instance'] = ''
                    try: seq = unimap[entry['UniProt']].obj['Sequence']     
                    except:
                        self.printLog('#ERR','Failed to get sequence for %s:%s (%s)' % (elm,entry['Primary_Acc'],entry['ProteinName']))
                        entry['Reduced'] = 'N'; continue
                    if entry['Primary_Acc'].find('-') > 0:
                        isoform = '%s-%s' % (seq.getStr('AccNum'),entry['Primary_Acc'].split('-')[1])   # Might be secondary Acc
                        entry['Instance'] = seq.getSequence(ikey=isoform)[entry['Start']-1:entry['End']]
                        seqlen = len(seq.getSequence(ikey=isoform))
                        if not entry['Instance']: self.warnLog('Could not find %s isoform sequence!' % entry['Primary_Acc'])
                        else: entry['SeqAcc'] = isoform
                    if not entry['Instance']:
                        entry['Instance'] = seq.getSequence()[entry['Start']-1:entry['End']]
                        seqlen = seq.seqLen()
                        entry['SeqAcc'] = seq.getStr('AccNum')
                    if entry['Start'] == 1: entry['Instance'] = '^%s' % entry['Instance']
                    else: entry['Instance'] = 'X%s' % entry['Instance']
                    if entry['End'] >= seqlen: entry['Instance'] = '%s$' % entry['Instance']
                    else: entry['Instance'] = '%sX' % entry['Instance']
                    regex = string.replace(rje.strEscape(elmc.data(elm)['Regex'],'^$'),'[\^','[^')
                    if not rje.matchExp('(%s)' % regex,entry['Instance']):
                        self.warnLog('Cannot find %s %s in %s %s (%s-%s)' % (elm,elmc.data(elm)['Regex'],entry['Primary_Acc'],entry['Instance'],entry['Start'],entry['End']))
                        elmi.dropEntry(entry); dx += 1
            if dx: self.printLog('#ELM','%s ELM instances mapped to sequences; %s failed.' % (rje.iLen(elmi.entries()),rje.iStr(dx)))

            ### ~ [3] Generate ELM motif file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elmbase = '%s%s' % (self.getStr('SourcePath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            motif_file = '%s.motifs' % elmbase
            sentry = sdb.addEntry({'Name':'ELMMotifs','File':motif_file,'Status':'Found'})
            if rje.checkForFile(motif_file) and not self.force(): self.printLog('#MOTIF','%s file found (force=F).' % motif_file)
            else:
                rje.backup(self,motif_file)
                motif_out = []
                for elm in elmc.dataKeys():
                    entry = elmc.data(elm)
                    motif_out.append('%s  %s  # %s [%d ELM instances]' % (entry['ELMIdentifier'],entry['Regex'],entry['Description'],entry['#Instances']))
                open(motif_file,'w').write(string.join(motif_out,'\n'))
                self.printLog('#MOTIF','%s motif patterns output to %s' % (elmc.entryNum(),motif_file))
                sentry['Status'] = 'Made'

            ### ~ [4] Download and Process PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            supported_db = ['HINT']
            if self.getStr('PPISource') not in supported_db:
                ppifile = self.sourceDataFile('PPISource',expect=self.getBool('PPIBench'),ask=True)    # If this exists, use for everything!
                ppdb = self.db().addTable(ppifile,mainkeys=['Hub','Spoke'],name='PPISource')
                sdb.data('PPISource')['Entries'] = ppdb.entryNum()
            for spec in self.list['PPISpec']:
                ## ~ [4a] Setup new species-specific attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.setStr({'%s.%s' % (self.getStr('PPISource'),spec): '%s.%s.ppi.tdt' % (self.getStr('PPISource'),spec),
                             'PPI.%s' % spec: '%s.%s.pairwise.tdt' % (self.getStr('PPISource'),spec),
                             'DomPPI.%s' % spec: '%s.%s.domppi.tdt' % (self.getStr('PPISource'),spec),
                             'Uniprot.%s' % spec: 'uniprot.%s.dat' % spec})
                ## ~ [4b] Download/Check UniProt DAT and fasta files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ppidat = self.sourceDataFile('Uniprot.%s' % spec,expect=self.getBool('PPIBench'),ask=False)  # If this exists, will use
                if not ppidat and self.getBool('PPIBench'): self.warnLog('Something went wrong making/finding %s UniProt file.' % spec, quitchoice=self.getBool('Integrity'))
                domdb = None
                domfile = self.sourceDataFile('DomPPI.%s' % spec,expect=False,ask=False)  # If this exists, will use
                #self.debug(domfile)
                ppuni = None
                ## ~ [4c] Download/Check/Generate PPI Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getStr('PPISource') in supported_db:   # Download sources are all found in sourceDataFile() for SOURCE.SPEC
                    ppifile = self.sourceDataFile('PPI.%s' % spec,expect=False,ask=False)  # If this exists, will use
                    needsource = (self.getBool('PPIBench') and not ppifile) or (self.getBool('DomBench') and not domfile)  or (self.getBool('DomBench') and not ppifile)
                    ppisource = self.sourceDataFile('%s.%s' % (self.getStr('PPISource'),spec),expect=needsource,ask=False)  # If this exists, will use
                ## ~ [4d] Download/Check/Generate Pairwise PPI Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    ppdb = None
                    if ppifile:
                        if not ppisource: self.warnLog('%s exists but source %s not found.' % (self.getStr('PPI.%s' % spec),self.getStr('PPISource')),quitchoice=self.getBool('Integrity'))
                        ppdb = self.db().addTable(ppifile,mainkeys=['Hub','Spoke'],name='PPI.%s' % spec)
                    elif not ppifile and self.getBool('PPIBench') or self.getBool('DomBench'):  ## Generate Pairwise PPI file
                        ppifile = string.replace(ppisource,'ppi','pairwise')
                        ppi = rje_ppi.PPI(self.log,self.cmd_list)
                        ppi.loadPairwisePPI(ppisource)
                        ppdb = ppi.db('Edge')
                        ppdb.setStr({'Name':'PPI.%s' % spec})
                        if self.getStrLC('PPID') in ['gene'] or (ppi.getStrLC('HubField').startswith('gene') and not self.getStrLC('PPID').startswith('uni')):
                            for field in ['SpokeUni','HubUni']: ppdb.addField(field,after='Spoke',evalue='')
                            if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                            else:
                                ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=flybase,pfam'])   #!# Make sure to add more as needed #!#
                                ppuni.setStr({'Name':ppidat})
                                ppuni.baseFile(rje.baseFile(ppidat))
                                ppuni.readUniProt()
                                sdb.data('Uniprot.%s' % spec)['Entries'] = ppuni.entryNum()
                                if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                            for entry in ppuni.entries():
                                gene = entry.gene().upper()   #!# Need to add additional mapping for DROME etc. #!#
                                acc = entry.acc()
                                if 'FlyBase' in entry.dict['DB']:
                                    for fbg in entry.dict['DB']['FlyBase']:
                                        for field in ['Id_A','Gene_A']:
                                            if field in ppdb.fields():
                                                for pentry in ppdb.indexEntries(field,fbg): pentry['HubUni'] = acc
                                        for field in ['Id_B','Gene_B']:
                                            if field in ppdb.fields():
                                                for pentry in ppdb.indexEntries(field,fbg): pentry['SpokeUni'] = acc
                                #if gene in ppdb.index('Hub'):
                                for pentry in ppdb.indexEntries('Hub',gene): pentry['HubUni'] = acc
                                for pentry in ppdb.indexEntries('Spoke',gene): pentry['SpokeUni'] = acc
                        else:   # Assume Hub and Spoke ARE UniProt IDs
                            for field in ['SpokeUni','HubUni']: ppdb.makeField('#%s#' % field[:-3],field,after='Spoke')
                        #ppdb.dropEntriesDirect('Hub',[''])
                        #ppdb.dropEntriesDirect('Spoke',[''])
                        ppdb.saveToFile(ppifile)
                        self.db().list['Tables'].append(ppdb)
                    if ppdb: sdb.data('PPI.%s' % spec)['Entries'] = ppdb.entryNum()
                    if not domfile and self.getBool('PPIBench'):  ## Generate Pairwise PPI file
                        domfile = string.replace(ppisource,'ppi','domppi')  #!# Is this used for anything ever?
                ## ~ [4e] Generate Protein-Pfam Links ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if domfile: domdb = self.db().addTable(domfile,mainkeys=['Pfam','Spoke'],name='DomPPI.%s' % spec,expect=False) #!# Is this used for anything ever?
                else:
                    domdb = self.db().addTable(mainkeys=['Pfam','Spoke'],name='DomPPI.%s' % spec,expect=False) #!# Is this used for anything ever?
                    sdb.data('DomPPI.%s' % spec)['File'] = '%s.DomPPI.%s.tdt' % (db.basefile(),spec)
                    sdb.data('DomPPI.%s' % spec)['Status'] = 'Found'
                if self.getBool('DomLink') or (domfile and not domdb):
                    makedomdb = not domdb
                    if makedomdb:
                        domdb = self.db().addEmptyTable('DomPPI.%s' % spec,['Pfam','Spoke','SpokeUni'],keys=['Pfam','Spoke'])
                        sdb.data('DomPPI.%s' % spec)['File'] = '%s.DomPPI.%s.tdt' % (db.basefile(),spec)
                        sdb.data('DomPPI.%s' % spec)['Status'] = 'Generated from Uniprot Pfam links'
                    pfamdb = self.db().addTable(name='Pfam.%s' % spec,mainkeys=['Uniprot'],expect=False)
                    makepfamdb = not pfamdb
                    if makepfamdb: pfamdb = self.db().addEmptyTable('Pfam.%s' % spec,['Uniprot','Pfam'],keys=['Uniprot'])
                    elif not makedomdb: continue    # Both files already made and loaded!
                    if not ppuni:
                        if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                        else:
                            ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['dbparse=flybase,pfam'])   #!# Make sure to add more as needed #!#
                            ppuni.setStr({'Name':ppidat})
                            ppuni.readUniProt()
                            sdb.data('Uniprot.%s' % spec)['Entries'] = ppuni.entryNum()
                            if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                    pfx = 0; epx = 0; ux = 0; utot = ppuni.entryNum()
                    for entry in ppuni.entries():
                        self.progLog('#PFAM','Adding DomPPI for %s UniProt entries: %.2f%%' % (rje.iStr(utot),ux/utot)); ux += 100.0
                        acc = entry.acc()
                        if 'Pfam' not in entry.dict['DB']: continue
                        epx += 1
                        if makepfamdb: pfamdb.addEntry({'Uniprot':acc,'Pfam':string.join(rje.sortKeys(entry.dict['DB']['Pfam']),'|')})
                        if not makedomdb: pfx += len(entry.dict['DB']['Pfam']); continue
                        for pfam in entry.dict['DB']['Pfam']:
                            pfx += 1
                            if ppdb:
                                for pentry in ppdb.indexEntries('HubUni',acc): domdb.addEntry({'Pfam':pfam,'Spoke':pentry['Spoke'],'SpokeUni':pentry['SpokeUni']},warn=False)
                    if makedomdb:
                        domdb.saveToFile(domfile)
                        self.printLog('#PFAM','Added %s DomPPI for %s Pfam domains in %s of %s UniProt entries.' % (rje.iStr(domdb.entryNum()),rje.iStr(pfx),rje.iStr(epx),rje.iStr(utot)))
                    if makepfamdb: pfamdb.saveToFile('%sslimbench.Pfam.%s.tdt' % (self.getStr('SourcePath'),spec)) #?#
                sdb.data('DomPPI.%s' % spec)['Entries'] = domdb.entryNum()

            ### ~ [5] Randomiser Source Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SimBench') or self.getBool('RanBench'):
                randfas = self.sourceDataFile('RandSource',expect=False,ask=False)  # If this exists, will use
                #self.deBug(randfas)
                while not randfas and not self.getStrLC('RandSource'):
                    self.setStr({'RandSource':rje.choice('Need randsource=X, where X is a file or Uniprot species code:')})
                while not randfas:
                    datkey = 'Uniprot.%s' % self.getStr('RandSource').upper()   # Try treating as species
                    #self.deBug(datkey)
                    if datkey not in self.str: self.setStr({datkey:'%s.dat' % datkey})
                    #self.deBug(datkey)
                    randdat = self.sourceDataFile(datkey,expect=True,ask=False)
                    #self.deBug(randdat)
                    if randdat:
                        if self.getBool('RanDat'): self.setStr({'RandSource':randdat})
                        else:
                            self.setStr({'RandSource':rje.baseFile(randdat) + '.fas'})
                            rseq = rje_seq.SeqList(self.log,['seqin=%s' % randdat,'seqout=%s' % self.getStr('RandSource'),'autoload=T','gnspacc=T','memsaver=T','accnr=F','seqnr=F','autofilter=T'])
                            #!# autofilter=T fudge to get this to work (for some reason) #!#
                            #!# Replace with UniProt code at some point - option to output splice variants too? #!#
                            sdb.data('RandSource')['Status'] = 'Made'
                            sdb.data('RandSource')['File'] = self.getStr('RandSource')
                            #self.deBug(self.getStr('RandSource')); self.deBug(rje.exists(self.getStr('RandSource')))
                        randfas = self.sourceDataFile('RandSource',expect=True,ask=False)  # If this exists, will use
                    else:
                        self.str.pop(datkey)
                        self.setStr({'RandSource':rje.choice('Source for Random Datasets. (Blank to exit.)')})
                        randfas = self.sourceDataFile('RandSource',expect=True,ask=False)  # If this exists, will use
                if not rje.exists(self.getStr('RandSource')): self.warnLog('Something went wrong making %s fasta file.' % self.getStr('RandSource'), quitchoice=self.getBool('Integrity'))

            ### ~ [6] Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SOURCE SUMMARY ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            for skey in sdb.dataKeys():
                entry = sdb.data(skey)
                if entry['Entries'] == '': self.printLog('#SOURCE','%s (%s): %s' % (entry['Name'],entry['File'],entry['Status']))
                else: self.printLog('#SOURCE','%s %s (%s): %s' % (rje.iStr(entry['Entries']),entry['Name'],entry['File'],entry['Status']))
            #!# Add integrity/date check using sb entries #!#
            if 'Download' in self.db('Source').index('Status'): sdb.saveToFile()
            if not self.force() and 'Download' in self.db('Source').index('Status'): 
                self.setBool({'Force':rje.yesNo('New source download detected: switch force=T for SLiMBench?')})
            return True     # Setup successful            
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem during SLiMBench.setupSourceData().'); return False  # Setup failed
#########################################################################################################################
    ### <3> ### SLiMBench Generator Methods                                                                             #
#########################################################################################################################
    def generate(self): ### Main SLiMBench Generator run method                                                     #V2.0
        '''Main SLiMBench Generator run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupSourceData(): return False            
            rje.mkDir(self,self.getStr('GenPath'),True)
            ## ~ [1a] Check GenSpec setting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.list['GenSpec'] = rje.listUpper(self.list['GenSpec']); self.list['GenSpec'].sort()
            genspecfile = '%sgenspec.txt' % self.getStr('GenPath')
            if os.path.exists(genspecfile):
                genspec = rje.listUpper(rje.listFromCommand(genspecfile,checkfile=True))
                genspec.sort()
                if not self.list['GenSpec']:
                    self.warnLog('%s file found: will use for GenSpec filtering.' % genspecfile,quitchoice=self.i() >= 0)
                    self.list['GenSpec'] = genspec
                else:
                    if self.list['GenSpec'] != genspec:
                        if self.i() >= 0 & rje.yesNo('Replace %d genspec=LIST codes with %d codes from %s?' % (len(self.list['GenSpec']),len(genspec),genspecfile)):
                            self.list['GenSpec'] = genspec
                        else:
                            self.warnLog('%s and genspec=LIST mismatch!' % genspecfile,quitchoice=True)
                            genspec = self.list['GenSpec']
            else: genspec = self.list['GenSpec']
            if genspec:
                elmi = self.db('ELMInstance'); dx = 0
                for ekey in elmi.dataKeys():
                    entry = elmi.data(ekey)
                    spec = string.split(entry['ProteinName'],'_')[-1]
                    if spec not in genspec: elmi.dropEntry(entry); dx += 1
                if dx: self.printLog('#SPEC','%s instances dropped based on genspec=LIST: %s remain.' % (rje.iStr(dx), rje.iStr(elmi.entryNum())))
                open(genspecfile,'w').write('%s\n' % string.join(genspec,'\n'))
            ### ~ [2] Generate basic ELM datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.fullELMDatasets()      # Also makes OccBench full.fas
            if self.getBool('PPIBench'): self.ppiELMDatasets()
            if self.getBool('DomBench'): self.pfamELMDatasets()
            ### ~ [3] Perform SLiMMaker Reduction of ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.reducedELMs()          # Also makes OccBench reduced.fas
            ### ~ [4] Run SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] ELMBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('SLiMMaker'):
                if self.getBool('ELMBench'): self.slimProb()
            else:
                if self.getBool('ELMBench'): self.slimProb('ELM_Datasets','Regex')
            if self.getBool('OccBench'): self.slimProbOcc()
            ### ~ [5] Select Subset of ELMs and generate Query Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [5a] ELMBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('ELMBench'):
                if self.getBool('SLiMMaker'): datadir = 'ELM_Reduced'
                else: datadir = 'ELM_Datasets'
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ %s DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir))
                elmlist = []
                for entry in self.db('%s.SLiMProb' % datadir).entries():
                    if float(entry['IC']) >= self.getNum('MinIC') and float(entry['N_UPC']) >= self.getInt('MinUPC'): elmlist.append(entry['Dataset'])
                elmlist.sort()
                filtertxt = 'ELMs meet IC >= %.2f & UP >= %d cutoffs' % (self.getNum('MinIC'),self.getInt('MinUPC'))
                self.printLog('#ELM','%s: %d ELMs meet IC >= %.2f & UP >= %d cutoffs' % (datadir,len(elmlist),self.getNum('MinIC'),self.getInt('MinUPC')))
                self.filterDatasets(elmlist,datadir,filtertxt)
                if self.getBool('Queries'): self.queryDatasets(elmlist,datadir)
            ## ~ [5b] PPIBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('PPIBench') or self.getBool('DomBench'):
                slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list+['force=F'])
                slimlist.setNum({'MinIC':self.getNum('MinIC')})
                slimlist.loadMotifs(self.getStr('ELMClass'))
                minicelms = slimlist.nameList(remsplit=True)
            if self.getBool('PPIBench'):
                for spec in self.list['PPISpec']:
                    datadir = 'PPI.%s' % spec
                    upcdir = '%s%s/SLiMCore/' % (self.getStr('GenPath'),datadir)
                    self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ %s DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir))
                    slimcore = rje_slimcore.SLiMCore(self.log,self.cmd_list+['resdir=%s' % upcdir,'batch=%s%s/*.fas' % (self.getStr('GenPath'),datadir)])
                    slimcore.run()  # Generates UPC files in upcdir
                    elmlist = []   # NOTE: elmlist is the list of datasets (w/o .fas), not actual ELMs => Cannot apply MinIC setting to PPI Datasets at this stage
                    for dset in glob.glob('%s/*.fas' % datadir):
                        dsetbase = rje.baseFile(dset,strip_path=True)
                        ufile = '%s%s.upc' % (upcdir,dsetbase)
                        if not rje.exists(ufile):
                            if slimcore.getInt('MaxSeq') > 0: self.printLog('#UPC','%s not found. Assumed SeqNum > MaxSeq (%d)' % (ufile,slimcore.getInt('MaxSeq')))
                            else: self.warnLog('%s not found and no MaxSeq cutoff!' % ufile)
                            continue
                        udata = open(ufile,'r').readline()
                        self.printLog('\r#UPC',udata)
                        udata = rje.matchExp('#(\S.*)# (\d+) Seq; (\d+) UPC; (\S+) MST',udata)
                        if int(udata[1]) > slimcore.getInt('MaxSeq') > 0: continue
                        if int(udata[2]) < self.getInt('MinUPC'): continue
                        # Check for any ELMs with enough IC
                        if rje.listIntersect(minicelms,self.hub2elm(string.split(dsetbase,'.')[0])): elmlist.append(dsetbase)
                    #for entry in self.db('%s.SLiMProb' % datadir).entries():
                    #    if float(entry['IC']) >= self.getNum('MinIC') and float(entry['N_UPC']) >= self.getInt('MinUPC'): elmlist.append(entry['Dataset'])
                    #elmlist.sort()
                    filtertxt = 'Hubs interact with 1+ ELMs with IC >= %.2f; PPI UP >= %d' % (self.getNum('MinIC'),self.getInt('MinUPC'))
                    if slimcore.getInt('MaxSeq') > 0: filtertxt += '; PPI Seq <= %d' % slimcore.getInt('MaxSeq')
                    self.printLog('#PPI','%s: %s %s' % (datadir,rje.iLen(elmlist),filtertxt))
                    self.filterDatasets(elmlist,datadir,filtertxt)
                    if self.getBool('Queries'): self.queryDatasets(elmlist,datadir)
            ## ~ [5c] DomBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DomBench'):
                for spec in self.list['PPISpec']:
                    datadir = 'DomPPI.%s' % spec
                    upcdir = '%s%s/SLiMCore/' % (self.getStr('GenPath'),datadir)
                    self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ %s DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir))
                    slimcore = rje_slimcore.SLiMCore(self.log,self.cmd_list+['resdir=%s' % upcdir,'batch=%s%s/*.fas' % (self.getStr('GenPath'),datadir)])
                    slimcore.run()  # Generates UPC files in upcdir
                    elmlist = []   # NOTE: elmlist is the list of datasets (w/o .fas), not actual ELMs => Cannot apply MinIC setting to PPI Datasets at this stage
                    for dset in glob.glob('%s/*.fas' % datadir):
                        dsetbase = rje.baseFile(dset,strip_path=True)
                        ufile = '%s%s.upc' % (upcdir,dsetbase)
                        if not rje.exists(ufile):
                            if slimcore.getInt('MaxSeq') > 0: self.printLog('#UPC','%s not found. Assumed SeqNum > MaxSeq (%d)' % (ufile,slimcore.getInt('MaxSeq')))
                            else: self.warnLog('%s not found and no MaxSeq cutoff!' % ufile)
                            continue
                        udata = open(ufile,'r').readline()
                        self.printLog('\r#UPC',udata)
                        udata = rje.matchExp('#(\S.*)# (\d+) Seq; (\d+) UPC; (\S+) MST',udata)
                        if int(udata[1]) > slimcore.getInt('MaxSeq') > 0: continue
                        if int(udata[2]) < self.getInt('MinUPC'): continue
                        # Check for any ELMs with enough IC
                        if rje.listIntersect(minicelms,self.hub2elm(string.split(dsetbase,'.')[0])): elmlist.append(dsetbase)
                    #for entry in self.db('%s.SLiMProb' % datadir).entries():
                    #    if float(entry['IC']) >= self.getNum('MinIC') and float(entry['N_UPC']) >= self.getInt('MinUPC'): elmlist.append(entry['Dataset'])
                    #elmlist.sort()
                    filtertxt = 'Hubs interact with 1+ ELMs with IC >= %.2f & PPI UP >= %d' % (self.getNum('MinIC'),self.getInt('MinUPC'))
                    if slimcore.getInt('MaxSeq') > 0: filtertxt += '; PPI Seq <= %d' % slimcore.getInt('MaxSeq')
                    self.printLog('#PPI','%s: %s %s' % (datadir,rje.iLen(elmlist),filtertxt))
                    self.filterDatasets(elmlist,datadir,filtertxt)
                    if self.getBool('Queries'): self.queryDatasets(elmlist,datadir)
            ### ~ [6] Simulated and/or random datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SimBench') or self.getBool('RanBench'): return self.simulator()
            return True
        except KeyboardInterrupt: raise
        except: self.errorLog('%s.generate error' % self,quitchoice=self.getBool('Integrity'))
        return not self.getBool('Generate')
#########################################################################################################################
    def ppiELMDatasets(self):     ### Generates ELM PPI Datasets                                                    #V2.0
        '''Generates PPI ELM Datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            force = self.force()
            db = self.obj['DB']
            elmp = db.getTable('ELMInteractors')
            elmd = db.getTable('ELMDomains')
            itype = self.getStrLC('IType')
            pbdb = db.addEmptyTable('ppibench.%s' % itype,['ELM','Hub','Link'],['ELM','Hub'])
            #unimap = self.obj['UniProt'].accDictFromEntries()   #           self.dict['UniProt']
            ### ~ [1] ~ Generate ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['PPISpec']:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ %s PPI ELM DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % spec)
                datadir = rje.makePath('%sPPI.%s/' % (self.getStr('GenPath'),spec))
                if os.path.exists(datadir) and force: rje.deleteDir(self,datadir)
                rje.mkDir(self,datadir,True)
                seqfilex = []
                if self.db('PPISource'): pdb = self.db('PPISource')
                else: pdb = self.db('PPI.%s' % spec)
                pdb.index('HubUni')
                if self.getBool('DomLink'):
                    pfamdb = self.db('Pfam.%s' % spec)
                    pfamdb.index('Pfam',splitchar='|')
                ppidat = self.getStr('Uniprot.%s' % spec)
                if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                else:
                    ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=%s' % self.getStr('SourcePath')])
                    ppuni.setStr({'Name':ppidat})
                    ppuni.readUniProt()
                    if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                unimap = ppuni.accDictFromEntries()     # Map all accession numbers onto species UniProt entries
                ## ~ [1a] Take each hub in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Interested in pairs of Elm:InteractorDomain (NB. Might add Elm:Domain datasets at some point too)
                fullelmppi = []
                for elm in elmp.index('Elm'):
                    directppi = elmp.indexDataList('Elm',elm,'interactorDomain')   # interactorDomain = Uniprot entry
                    if self.getBool('DomLink'):
                        elmppi = []
                        for pfam in rje.sortUnique(elmd.indexDataList('Elm',elm,'Domain') + elmp.indexDataList('Elm',elm,'Domain')):
                            if pfam in pfamdb.index('Pfam'): elmppi += pfamdb.index('Pfam')[pfam]
                        elmppi = rje.sortUnique(elmppi)
                    else: elmppi = directppi[0:]
                    #self.deBug('%s ELMPPI: %s' % (elm,elmppi))
                    unihubx = 0     # Number of mapped UniProt hubs
                    seqfilex.append(0)
                    for acc in elmppi:  # This is now a hub protein: find in ppi
                        if acc not in unimap: continue    # Wrong species and/or failed mapping
                        if acc in directppi: pbdb.addEntry({'ELM':elm,'Hub':acc,'Link':'ELMInteractors'})
                        else: pbdb.addEntry({'ELM':elm,'Hub':acc,'Link':'DomLink'})
                        if acc in fullelmppi: continue      # This hub is already present!
                        else: fullelmppi.append(acc)        # Only output once!
                        #x#elmseqfile = '%s%s.%s.fas' % (datadir,elm,acc)   #i# Removed ELM from dataset name.
                        elmseqfile = '%s%s.%s.fas' % (datadir,acc,itype)     #!# Add option to map to Symbol?
                        unihubx += 1
                        if rje.exists(elmseqfile) and not self.force(): seqfilex[-1] += 1; continue
                        uni = unimap[acc]
                        hubacc = uni.accNum()
                        #self.deBug(hubacc)
                        if hubacc not in pdb.index('HubUni'):
                            self.printLog('#PPI','No %s PPI read for %s' % (spec,uni.shortName()))
                            continue
                        elmseqs = []
                        spokelist = pdb.indexDataList('HubUni',hubacc,'SpokeUni')
                        while '' in spokelist: spokelist.remove('')
                        self.printLog('#PPI','%s %s PPI parsed for %s' % (rje.iLen(spokelist),spec,uni.shortName()))
                        #self.deBug(spokelist)
                        for spoke in spokelist:
                            spacc = spoke.split('-')[0]     # PPI can have isoforms?
                            try: seq = (unimap[spacc].obj['Sequence'],spoke)   # SeqAcc can be isoform
                            except: self.printLog('#ERR','Failed to get sequence for %s.' % (spoke)); continue
                            if seq not in elmseqs: elmseqs.append(seq)
                        if elmseqs:
                            ELMSEQ = open(elmseqfile,'w')
                            for seq in elmseqs: ELMSEQ.write(seq[0].fasta(seq[1]))
                            ELMSEQ.close()
                            self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                            seqfilex[-1] += 1
                        else:
                            self.printLog('#PPI','No sequences to output for %s %s %s PPI' % (elm,spec,acc))
                            #self.deBug(elmp.indexEntries('Elm',elm))
                    self.printLog('#PPI','%s of %s %s interactors mapped to %s: %s PPI files made.' % (rje.iStr(unihubx),rje.iLen(elmppi),elm,spec,rje.iStr(seqfilex[-1])))
                self.printLog('#DSETS','%s ELM %s PPI datasets (force=%s)' % (rje.iStr(sum(seqfilex)),spec,self.force()))
                pbdb.saveToFile('%sppibench.%s.tdt' % (datadir,itype))
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem during SLiMBench.ppiELMDatasets().'); return False  
#########################################################################################################################
    def pfamELMDatasets(self):     ### Generates ELM-Pfam PPI Datasets                                              #V2.0
        '''Generates ELM-Pfam PPI Datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            force = self.force()
            db = self.obj['DB']
            elmd = db.getTable('ELMDomains')    # "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
            elmp = db.getTable('ELMInteractors')
            itype = 'dom%s' % self.getStrLC('IType')
            #unimap = self.obj['UniProt'].accDictFromEntries()   #           self.dict['UniProt']
            ### ~ [1] ~ Generate ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in self.list['PPISpec']:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ %s PFam PPI ELM DATASETS ~~~~~~~~~~~~~~~~~~~~~~~ #' % spec)
                datadir = rje.makePath('%sDomPPI.%s/' % (self.getStr('GenPath'),spec))
                if os.path.exists(datadir) and force: rje.deleteDir(self,datadir)
                rje.mkDir(self,datadir,True)
                seqfilex = []
                pdb = self.db('DomPPI.%s' % spec)
                pdb.index('Pfam')
                ppidat = self.getStr('Uniprot.%s' % spec)
                if spec in self.dict['UniSpec']: ppuni = self.dict['UniSpec'][spec]
                else:
                    ppuni = rje_uniprot.UniProt(self.log,self.cmd_list+['unipath=%s' % self.getStr('SourceData')])  #!# If not memsaver=T, store these #!#
                    ppuni.setStr({'Name':ppidat})
                    ppuni.readUniProt()
                    if not self.getBool('MemSaver'): self.dict['UniSpec'][spec] = ppuni
                unimap = ppuni.accDictFromEntries()     # Map all accession numbers onto species UniProt entries
                ## ~ [1a] Take each hub in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Interested in pairs of Elm:Domain
                fullelmppi = []
                for elm in rje.sortUnique(elmp.indexKeys('Elm') + elmd.indexKeys('Elm')):
                    elmppi = rje.sortUnique(elmd.indexDataList('Elm',elm,'Domain') + elmp.indexDataList('Elm',elm,'Domain'))
                    pfamx = 0     # Number of mapped Pfam hubs
                    seqfilex.append(0)
                    for acc in elmppi:  # This is now a Pfam domain: find in pdb
                        if acc in fullelmppi: continue      # This hub is already present!
                        else: fullelmppi.append(acc)        # Only output once!
                        #x#elmseqfile = '%s%s.%s.fas' % (datadir,elm,acc)
                        elmseqfile = '%s%s.%s.fas' % (datadir,acc,itype)
                        pfamx += 1
                        if rje.exists(elmseqfile) and not self.force(): seqfilex[-1] += 1; continue
                        if acc not in pdb.index('Pfam'):
                            self.printLog('#DPI','No %s DPI read for domain %s.' % (spec,acc))
                            continue
                        elmseqs = []
                        spokelist = pdb.indexDataList('Pfam',acc,'SpokeUni')
                        while '' in spokelist: spokelist.remove('')
                        self.printLog('#DPI','%s %s PPI parsed for domain %s.' % (rje.iLen(spokelist),spec,acc))
                        #self.deBug(spokelist)
                        for spoke in spokelist:
                            spacc = spoke.split('-')[0]     # PPI can have isoforms?
                            try: seq = (unimap[spacc].obj['Sequence'],spoke)   # SeqAcc can be isoform
                            except: self.printLog('#ERR','Failed to get sequence for %s.' % (spoke)); continue
                            if seq not in elmseqs: elmseqs.append(seq)
                        if elmseqs:
                            ELMSEQ = open(elmseqfile,'w')
                            for seq in elmseqs: ELMSEQ.write(seq[0].fasta(seq[1]))
                            ELMSEQ.close()
                            self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                            seqfilex[-1] += 1
                        else:
                            self.printLog('#DPI','No sequences to output for %s %s %s PPI' % (elm,spec,acc))
                            #self.deBug(elmp.indexEntries('Elm',elm))
                    self.printLog('#DPI','%s of %s %s interactors mapped to %s: %s DPI files made.' % (rje.iStr(pfamx),rje.iLen(elmppi),elm,spec,rje.iStr(seqfilex[-1])))
                self.printLog('#DSETS','%s ELM %s Domain PPI datasets (force=%s)' % (rje.iStr(sum(seqfilex)),spec,self.force()))
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem during SLiMBench.pfamELMDatasets().'); return False  
#########################################################################################################################
    def fullELMDatasets(self):     ### Generates Full ELM Datasets                                                  #V1.x
        '''Generates Full ELM Datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('ELMBench') and not self.getBool('OccBench'): return
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ FULL ELM DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            force = self.force()
            db = self.obj['DB']
            datadir = rje.makePath('%sELM_Datasets/' % self.getStr('GenPath'))
            if os.path.exists(datadir) and force: rje.deleteDir(self,datadir)
            rje.mkDir(self,datadir,True)
            elmc = db.getTable('ELMClass')
            elmc.addField('FullSeqNum')
            elmi = db.getTable('ELMInstance')
            unimap = self.dict['UniProt']
            ### ~ [1] ~ Generate ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfilex = 0
            for elm in elmc.dataKeys():
                elmseqfile = '%s%s.fas' % (datadir,elm)
                if not force and rje.checkForFile(elmseqfile):
                    elmc.data(elm)['FullSeqNum'] = rje_seq.SeqCount(self,elmseqfile)
                    seqfilex += 1
                else:
                    rje.backup(self,elmseqfile,appendable=False)
                    elmseqs = []
                    if not elmi.indexEntries('ELMIdentifier',elm): self.printLog('#ELM','No instance data for %s' % elm); continue
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        try: seq = (unimap[entry['UniProt']].obj['Sequence'],entry['SeqAcc'])
                        except: self.printLog('#ERR','Failed to get sequence for %s:%s' % (elm,entry['ProteinName'])); continue
                        if seq not in elmseqs: elmseqs.append(seq)
                    elmc.data(elm)['FullSeqNum'] = len(elmseqs)
                    if elmseqs:
                        ELMSEQ = open(elmseqfile,'w')
                        for seq in elmseqs: ELMSEQ.write(seq[0].fasta(seq[1]))
                        ELMSEQ.close()
                        self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                        seqfilex += 1
                    else:
                        self.printLog('#ELM','No sequences to output for %s' % elm)
                        #self.deBug(elmi.indexEntries('ELMIdentifier',elm))
            self.printLog('#DSETS','%d ELM sequence datasets (force=%s)' % (seqfilex,force))
            ### ~ [2] ~ Generate OccBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('OccBench'): return
            occdir = rje.makePath('%sELM_OccBench/' % self.getStr('GenPath'))
            if os.path.exists(occdir) and force: rje.deleteDir(self,occdir)
            rje.mkDir(self,occdir,True)
            occfile = occdir + 'ELM.full.fas'
            if rje.exists(occfile): self.printLog('#OCC','%s found (force=F).' % occfile); return
            occseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqout=%s' % occfile,'seqmode=list'])
            for elm in elmc.dataKeys():
                elmseqfile = '%s%s.fas' % (datadir,elm)
                if rje.exists(elmseqfile): occseq.loadSeq(elmseqfile,clearseq=False)
                else: self.warnLog('ELM file %s missing!' % elmseqfile)
            occseq.saveSeq()
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem during SLiMBench.fullELMDatasets().'); return False  
#########################################################################################################################
    def reducedELMs(self):  ### Performs ELM reduction using SLiMMaker                                              #V1.x
        '''Performs ELM reduction using SLiMMaker.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('SLiMMaker'): return
            if not (self.getBool('ELMBench') or self.getBool('OccBench') or self.getBool('SimBench') or self.getBool('RanBench')): return
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SLiMMaker ELM Reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            motif_file = '%s.reduced.motifs' % elmbase
            datadir = rje.makePath('%sELM_Reduced/' % self.getStr('GenPath'))
            if os.path.exists(datadir) and self.force(): rje.deleteDir(self,datadir)
            rje.mkDir(self,datadir,True)
            instdir = rje.makePath('%sELM_Instances/' % self.getStr('GenPath'))
            if os.path.exists(datadir) and self.force(): rje.deleteDir(self,instdir)
            rje.mkDir(self,instdir,True)
            unimap = self.dict['UniProt']
            maker = self.obj['SLiMMaker']
            ## ~ [0a] Add extra fields to ELM Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmc = db.getTable('ELMClass')
            elmc.addField('SLiMMaker')
            elmc.addField('ReducedInstance')
            elmc.addField('ReducedSeqNum')
            elmi = db.getTable('ELMInstance')
            elmi.addField('Reduced')

            ### ~ [1] ~ Generate Reduced ELMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdbfile = '%s.reduced.tdt' % elmbase
            if rje.checkForFile(motif_file) and not self.force():
                self.printLog('#MOTIF','%s file found (force=F).' % motif_file)
                ## ~ [1a] Populate from existing file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for rmotif in open(motif_file,'r').readlines():
                    rdata = rje.matchExp('(\S+)\s+(\S+)\s+# (\S.+). ELM Definition: (\S+) \[(\d+) -> (\d+) ELM instances\]',rmotif)
                    if rdata:
                        elm = rdata[0]
                        centry = elmc.data(elm)
                        centry['SLiMMaker'] = rdata[1]
                        centry['ReducedInstance'] = int(rdata[-1])
                        for entry in elmi.indexEntries('ELMIdentifier',elm):
                            regex = regex = string.replace(rje.strEscape(rdata[1],'^$'),'[\^','[^')
                            if rje.matchExp('(%s)' % regex,entry['Instance']): entry['Reduced'] = 'Y'
                            else: entry['Reduced'] = 'N'
                    else: self.deBug(rmotif)
            else:
                ## ~ [1b] ~ Generate Reduced ELMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                equiv = []
                for elm in elmc.dataKeys():
                    instances = []
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        if entry['Instance']: instances.append(entry['Instance'])
                    open('%s%s.instances.txt' % (instdir,elm),'w').write(string.join(instances,'\n'))
                    #self.deBug('%s: %s' % (elmc.data(elm)['Regex'],instances))
                    # Consider splitting by length and then recombine with wildcard spacers?! #
                    self.printLog('#SLIM','%s: %s (%s instances)' % (elm,elmc.data(elm)['Regex'],elmc.data(elm)['#Instances']))
                    maker.list['Peptides'] = instances; maker.list['Input'] = []
                    maker.setStr({'PeptAlign':elmc.data(elm)['Regex']})
                    maker.obj['PeptCluster'].setStr({'PeptAlign':elmc.data(elm)['Regex']})
                    slim = maker.run(iterate=True)[0]       #!# Make iterate an option? #!#
                    if slim and rje_slim.patternIC(slim) < self.getNum('MinIC'):
                        self.printLog('#MAKE','%s: %s -> "%s" -> IC filtered' % (elm,elmc.data(elm)['Regex'],slim))
                        slim = ''
                    if slim:
                        for el in rje_slim.slimFromPattern(slim).split('-'):
                            if len(el) > 1 and el not in equiv: equiv.append(el)
                    else: maker.list['Peptides'] = []
                    elmc.data(elm)['SLiMMaker'] = slim
                    elmc.data(elm)['ReducedInstance'] = len(maker.list['Peptides'])
                    maker.list['Peptides'] = string.split(string.replace(string.join(maker.list['Peptides']),'-',''))
                    self.debug(maker.list['Peptides'])
                    self.debug(elmi.indexDataList('ELMIdentifier',elm,'Instance'))
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        if entry['Instance'] and entry['Instance'] in maker.list['Peptides']: entry['Reduced'] = 'Y'
                        else: entry['Reduced'] = 'N'
                    self.debug(elmi.indexDataList('ELMIdentifier',elm,'Reduced'))
                    #self.deBug('%s: %s' % (elmc.data(elm)['Regex'],slim))
                    self.printLog('#MAKE','%s: %s -> "%s" = %d -> %d instances' % (elm,elmc.data(elm)['Regex'],slim,elmc.data(elm)['#Instances'],elmc.data(elm)['ReducedInstance']))
                ## ~ [1c] ~ Create Reduced ELM Motif File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
                ## ~ [1d] ~ Create equivalence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                equiv.sort()
                for eq in equiv[0:]:
                    if string.join(equiv).count(eq) > 1: equiv.remove(eq)
                equivfile = rje.baseFile(motif_file) + '.equiv.txt'
                rje.backup(self,equivfile,appendable=False)
                open(equivfile,'w').write(string.join(equiv,'\n'))
                #!# Remove redundant substrings! #!#
            if self.force() or not rje.checkForFile(rdbfile): elmc.saveToFile(rdbfile)

            ### ~ [2] ~ Generate Reduced ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('ELMBench') and not self.getBool('OccBench'): return
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ Reduced ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            seqfilex = 0; redelms = []  # File count and list of Reduced ELMs output
            for elm in elmc.dataKeys():
                elmseqfile = '%s%s.reduced.fas' % (datadir,elm)
                if not self.force() and rje.checkForFile(elmseqfile):
                    elmc.data(elm)['ReducedSeqNum'] = rje_seq.SeqCount(self,elmseqfile)
                    seqfilex += 1; redelms.append(elm)
                else:
                    rje.backup(self,elmseqfile,appendable=False)
                    elmseqs = []
                    if not elmi.indexEntries('ELMIdentifier',elm): self.printLog('#ELM','No instance data for %s' % elm); continue
                    for entry in elmi.indexEntries('ELMIdentifier',elm):
                        if entry['Reduced'] == 'N': continue
                        try: seq = (unimap[entry['UniProt']].obj['Sequence'],entry['SeqAcc'])
                        except: self.errorLog('Failed to get sequence for %s:%s' % (elm,entry['ProteinName']),quitchoice=False); continue
                        if seq not in elmseqs: elmseqs.append(seq)
                    elmc.data(elm)['ReducedSeqNum'] = len(elmseqs)
                    if elmseqs:
                        ELMSEQ = open(elmseqfile,'w')
                        for seq in elmseqs: ELMSEQ.write(seq[0].fasta(seq[1]))
                        ELMSEQ.close()
                        self.printLog('#SEQ','%s sequences output to %s' % (len(elmseqs),elmseqfile))
                        seqfilex += 1; redelms.append(elm)
                    else: self.printLog('#ELM','No sequences to output for %s' % elm)
            self.printLog('#DSETS','%d Reduced ELM sequence datasets (force=%s)' % (seqfilex,self.force()))
            elmc.saveToFile(rdbfile,backup=False)
            if len(redelms) != seqfilex: raise ValueError('Reduced ELM seqfile count (%d) and ELM list size (%d) mismatch!' % (seqfilex,len(redelms)))

            ### ~ [3] ~ Generate OccBench Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('OccBench'): return
            occdir = rje.makePath('%sELM_OccBench/' % self.getStr('GenPath'))
            occfile = occdir + 'ELM.reduced.fas'
            if rje.exists(occfile):
                if self.force(): rje.backup(self,occfile)
                else: self.printLog('#OCC','%s found (force=F).' % occfile); return
            occseq = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqout=%s' % occfile,'seqmode=list'])
            noseq = elmc.indexDataList('ReducedSeqNum',0,'ELMIdentifier')
            for elm in redelms:
                elmseqfile = '%s%s.reduced.fas' % (datadir,elm)
                if rje.exists(elmseqfile): occseq.loadSeq(elmseqfile,clearseq=False)
                elif elm not in noseq: self.warnLog('Reduced ELM file %s missing!' % elmseqfile)
            occseq.saveSeq()

            return True
        except: self.errorLog('Problem during %s.reducedELMs().' % self); raise; return False
#########################################################################################################################
    def slimProb(self,datadir='ELM_Reduced',slimkey='SLiMMaker'):     ### Performs SLiMProb of specific motifs against their datasets
        '''Performs SLiMProb of specific motifs against their datasets.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ %s %s SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir,slimkey))
            db = self.obj['DB']
            basefile = rje.makePath('%s%s' % (self.getStr('GenPath'),datadir),wholepath=True)
            ssfile = '%s.slimprob.occ.csv' % basefile
            if rje.checkForFile(ssfile):
                if not self.force():
                    self.printLog('#SEARCH','SLiMProb file %s found (force=F).' % ssfile)
                    db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Dataset','RunID','Motif'],datakeys='All',name='%s.SLiMProb' % datadir)
                    #db.addTable(ssfile,mainkeys=['Dataset','RunID','Motif','Seq','Start_Pos'],datakeys='All',name='%s.SLiMProbOcc' % datadir)   #!# Why? #!#
                    return True
                else: rje.backup(self,ssfile,appendable=False)
            #self.deBug(rje.checkForFile(ssfile))
            datadirpath = rje.makePath(basefile)
            elmc = db.getTable('ELMClass')
            ex = 0; rx = 0
            ### ~ [1] ~ Perform SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for elmfas in glob.glob('%s*.fas' % datadirpath):
                elm = string.split(os.path.basename(elmfas),'.')[0]
                if elm not in elmc.dataKeys(): self.errorLog('Failed to find ELM "%s" in ELMClass table.' % elm,quitchoice=False); continue
                entry = elmc.data(elm)
                slim = entry[slimkey]
                if not slim: self.printLog('#SKIP','No "%s" SLiM to search for ELM "%s".' % (slimkey,elm)); continue
                ex += 1
                sscmd = ['extras=0','maxsize=0'] + self.cmd_list + ['runid=%s' % slimkey]   # Note that self.cmdlist defaults to minic=2.0
                if rje.checkForFile(self.getStr('SearchINI')): sscmd += rje.getCmdList(['ini=%s' % self.getStr('SearchINI')])
                sscmd += ['resfile=%s' % ssfile,'masking=F','basefile=None','occupc=T','append=T','motifs=','seqin=%s' % elmfas,'resdir=%sSLiMProb/' % datadirpath]
                ss = slimprob.SLiMProb(self.log,sscmd+['mergesplits=F'])  #+['debug=F'])
                slimlist = ss.obj['SlimList']
                slimlist.list['Motif'] = []
                slimlist._addMotif(elm,slim,check=True)
                if not slimlist.list['Motif']: self.printLog('#CUT','%s %s did not make the cut' % (elm,slim)); continue
                ss.run()
                rx += 1
            self.printLog('#SLIM','%d SLiMProb runs for %d ELMs (%s vs %s)' % (rx,ex,slimkey,datadir))
            ### ~ [2] ~ Load SLiMProb Results into Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try:
                db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Dataset','RunID','Motif'],datakeys='All',name='%s.SLiMProb' % datadir)
                #db.addTable(ssfile,mainkeys=['Dataset','RunID','Motif','Seq','Start_Pos'],datakeys='All',name='%s.SLiMProbOcc' % datadir)
                return True
            except:
                self.deBug(ssfile)
                self.errorLog('%s vs %s SLiMProb failed' % (slimkey,datadir)); return False
        except: self.errorLog('Problem during %s.slimProb().' % self); return False  
#########################################################################################################################
    def slimProbOcc(self):  ### Performs SLiMProb of all motifs against OccBench datasets to generate positives file.
        '''Performs SLiMProb of all motifs against OccBench datasets to generate positives file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ OccBench SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            motif_file = {'reduced':'%s%s.reduced.motifs' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True)),
                          'full':'%s%s.motifs' % (self.getStr('SourcePath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))}
            elmi = db.getTable('ELMInstance')   # mainkeys=['ELMIdentifier','Primary_Acc','Start','End']
            ### ~ [1] ~ Perform SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for etype in ['full']:  # Cannot do ,'reduced']: as instances do not map!
                basefile = '%sELM_OccBench/ELM.%s' % (self.getStr('GenPath'),etype)
                ## ~ [1a] ~ Check for existing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elmfas = '%s.fas' % basefile
                if not rje.checkForFile(elmfas):
                    self.printLog('#OCCFAS','OccBench fasta file %s not found.' % elmfas); continue
                resfile = '%s.csv' % basefile
                run_slimprob = self.force() or not rje.checkForFile(resfile)
                if rje.checkForFile(resfile):
                    if not self.force(): self.printLog('#SEARCH','SLiMProb file %s found (force=F).' % resfile)
                    else: rje.backup(self,resfile,appendable=False)
                ## ~ [1b] ~ Run SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if run_slimprob:
                    # Note that self.cmdlist defaults to minic=2.0
                    sscmd = ['extras=0'] + self.cmd_list + ['maxsize=0','runid=%s' % etype,'masking=F','resfile=%s' % resfile,'basefile=None','occupc=F','efilter=F','append=F','motifs=%s' % motif_file[etype],'seqin=%s' % elmfas,'resdir=%s.SLiMProb/' % basefile]
                    ssi = sscmd.index('minic=2.0')
                    sscmd[ssi] = 'minic=1.1'  # Change default to 1.1
                    self.debug(sscmd)
                    ss = slimprob.SLiMProb(self.log,sscmd+['mergesplits=F'])  #+['debug=F'])
                    ss.run()
                    self.printLog('#SLIM','SLiMProb run for %s OccBench data.' % (etype))
            ### ~ [2] ~ Rate SLiMProb Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                ratfile =  '%s.ratings.csv' % basefile
                rkeys = ['Motif','Seq','Start_Pos','End_Pos']
                if rje.checkForFile(ratfile):
                    if not self.force(): self.printLog('#SEARCH','SLiMProb ratings file %s found (force=F).' % ratfile); continue
                    else: rje.backup(self,ratfile,appendable=False)
                occfile = '%s.occ.csv' % basefile
                odb = db.addTable(occfile,mainkeys=rkeys,datakeys=rkeys,name='%s.rating' % etype)
                elmconv = {}
                for entry in odb.entries():
                    if entry['Motif'] in elmi.index('ELMIdentifier'): continue
                    if entry['Motif'] in elmconv: entry['Motif'] = elmconv[entry['Motif']]; continue
                    elmcore = string.join(string.split(entry['Motif'],'_')[:-1],'_')
                    #self.debug('%s: %s' % (elmcore,elmcore in elmi.index('ELMIdentifier')))
                    if elmcore in elmi.index('ELMIdentifier'): elmconv[entry['Motif']] = elmcore; entry['Motif'] = elmcore
                odb.remakeKeys()
                for elm in elmconv: self.printLog('#SPLIT','Split ELM %s -> %s' % (elm,elmconv[elm]))
                #self.debug(odb.indexKeys('Motif'))
                #self.debug(elmi.indexKeys('ELMIdentifier'))

                ## ~ [2a] Make a 'Seq':'Primary_Acc' mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqmap = {}
                primaries = []
                elms = []
                for elm in elmi.index('ELMIdentifier'):
                    if elm not in odb.index('Motif'): self.printLog('#ELM','ELM %s not found in %s.' % (elm,occfile),screen=self.v()>0 or self.debugging())
                    else: elms.append(elm); primaries += elmi.indexDataList('ELMIdentifier',elm,'SeqAcc')
                primaries = rje.sortUnique(primaries)   # Unique list of AccNum and isoforms
                ix = 0; sx = len(primaries)
                for seqacc in primaries:
                    for seq in odb.indexKeys('Seq'):
                        if seq.endswith(seqacc): seqmap[seqacc] = seq; break    # ELM instance AccNum points to sequence name
                    if seqacc not in seqmap: self.warnLog('Cannot find %s in %s results!' % (seqacc,occfile),warntype='acc_missing',quitchoice=True,suppress=True)
                if self.dev() or self.debugging():
                    elmi.indexReport('ELMIdentifier','#ELMI')
                    odb.indexReport('Motif','#OCC')
                ## ~ [2b] Rate occ.csv results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                ## >> Direct matches through AccNum = TP
                ## >> Other matches to different splice variants = X -> remove! (Might be TP or FP. Who knows?!)
                odb.addField('Rating',evalue='FP')
                fx = 0
                for entry in elmi.entries():
                    elm = entry['ELMIdentifier']
                    if elm not in elms: continue
                    seq = seqmap[entry['SeqAcc']]                   # Sequence name from SLiMProb results file
                    seqacc = string.split(entry['SeqAcc'],'-')[0]   # Strip splice variant information to get pure AccNum
                    ix += 1
                    ## Correct variant first
                    ekey = odb.makeKey({'Motif':elm,'Seq':seq,'Start_Pos':entry['Start'],'End_Pos':entry['End']})
                    try: odb.data(ekey)['Rating'] = 'TP'
                    except:
                        self.bugPrint(ekey)
                        if seq in odb.index('Seq'): self.bugPrint(odb.index('Seq')[seq])
                        else: self.bugPrint(odb.index('Motif')[elm])
                        self.warnLog('Failed to find ELM instance: %s' % ekey,quitchoice=self.dev()); fx += 1
                        #self.debug(entry)
                    ## Other entries -> X
                    for xentry in odb.indexEntries('Motif',elm):    # This also removes other hits to the SAME protein
                        if seqacc in xentry['Seq'] and xentry['Rating'] != 'TP': xentry['Rating'] = 'X'
                self.printLog('#OCC','Mapped %s instances of %s ELMs in %s sequences for rating' % (rje.iStr(ix),rje.iLen(elms),rje.iStr(sx)))
                if fx:
                    self.printLog('#FAIL','Failed to find %d instances. Might be an artefact of motif splitting for Regex compatibility. Some variants might have been screened by minic=%.2f.' % (fx,self.getNum('MinIC')))
                    if self.getBool('Integrity'): self.errorLog('Failed to find %d of %d ELM instances in %s.' %(fx,ix,occfile),printerror=False,quitchoice=True)
                    else: self.warnLog('Failed to find %d of %d ELM instances in %s.' %(fx,ix,occfile),'integrity')
                prex = odb.entryNum()
                odb.dropEntriesDirect('Rating',['X'])
                if prex != odb.entryNum():
                    self.printLog('#AMBIG','Potentially ambiguous instances, i.e. in a known TP protein (or isoform thereof) identified.')
                    self.printLog('#AMBIG','Dropped %s potentially ambiguous ("X"-Rated) instances.' % rje.iStr(prex - odb.entryNum()))
                odb.indexReport('Rating')
                odb.saveToFile(ratfile)
                db.deleteTable(odb)
            return True
        except: self.errorLog('Problem during %s.slimProbOcc().' % self); return False
#########################################################################################################################
    def queryDatasets(self,elmlist,datadir='ELM_Reduced'):  ### Produces ELM Query Datasets
        '''
        Produces ELM Query Datasets.
        >> datadir:path = Path to ELM dataset files
        >> elmlist:list = List of ELM datasets to make files for
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ %s QUERY DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir))
            db = self.db('ELMInstance')
            indir = rje.makePath('%s%s' % (self.getStr('GenPath'),datadir))
            outdir = rje.makePath('%s%s_WithQueries' % (self.getStr('GenPath'),datadir))
            rje.mkDir(self,outdir,True)
            self.list['FlankMask'] = string.split(string.join(self.list['FlankMask']).lower())
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','dna=F'])
            ### ~ [1] ~ Generate new ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qx = 0; esx = 0
            for dset in elmlist[0:]:     # These are now datasets!
                elm = string.split(dset,'.')[0]
                #for efile in glob.glob('%s%s.*fas' % (indir,elm)):
                efile  = '%s%s.fas' % (indir,dset)
                if not rje.exists(efile): self.warnLog('%s missing!' % efile)
                else:
                    self.printLog('#ELM','Process %s: %s' % (elm,efile))
                    seqlist.loadSeq(efile,nodup=True,clearseq=True,mode='list')     ### Loads sequences from file
                    eseqx = seqlist.seqNum()
                    eseqs = seqlist.list['Seq'][0:]
                    qseqs = seqlist.list['Seq'][0:]
                    self.debug(qseqs)
                    qoutx = 0
                    ## ~ [1a] ~ Take each Query in turn ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    while qseqs:
                        query = seqlist.getSeq(qseqs.pop(0))
                        self.debug(query)
                        qname = string.split(query[0])[0]
                        qlen = len(query[1])
                        qocc = []
                        #self.debug(db.indexEntries('ELMIdentifier',elm))    # NB. Will not work for PPI data
                        if elm in db.index('ELMIdentifier'): qelms = [elm]    # NB. Will not work for PPI data
                        else: qelms = self.hub2elm(elm)
                        for qelm in qelms:
                            for entry in db.indexEntries('ELMIdentifier',qelm):
                                if entry['ProteinName'] in string.split(qname,'__') or entry['Primary_Acc'] in string.split(qname,'__'):
                                    occ = (int(entry['Start'])-1,int(entry['End']))
                                    if occ not in qocc: qocc.append(occ)
                        if not qocc:
                            if 'PPI' not in datadir: self.errorLog('Cannot find any %s instances of %s in ELMInstances!' % (qname,elm),printerror=False)
                            eseqs.append(eseqs.pop(0))
                            continue
                        qocc.sort()
                        self.printLog('#QRY','%s - Query %s (%d occ)' % (elm,qname,len(qocc))); qx += 1
                        outx = 0
                        for flank in self.list['FlankMask']:
                            ## ~ [1b] ~ Establish Flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            outfile = '%s%s.%s.%s.fas' % (outdir,rje.baseFile(efile,strip_path=True),qname,flank)
                            self.debug(outfile)
                            if rje.exists(outfile) and not self.force(): outx += 1; continue
                            mask = qocc[0:]
                            if flank == 'none': mask = [(0,qlen)]
                            elif flank == 'site': pass
                            elif flank[:5] == 'flank':
                                win = string.atoi(flank[5:])
                                newmask = []
                                for old in mask: newmask.append((max(0,old[0]-win),min(qlen,old[1]+win)))
                                mask = newmask
                            elif flank[:3] == 'win':
                                win = string.atoi(flank[3:])
                                masklen = 0
                                for old in mask: masklen += (old[1] - old[0])
                                while masklen < win and masklen < qlen:
                                    newmask = []
                                    for old in mask:
                                        start = max(0,old[0]-1)
                                        end = min(qlen,old[1]+1)
                                        if newmask and start < newmask[-1][1]: newmask[-1] = (newmask[-1][0],end)
                                        else: newmask.append((start,end))
                                    mask = newmask
                                    masklen = 0
                                    for old in mask: masklen += (old[1] - old[0])
                            ## ~ [1c] ~ Mask Query ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            maskseq = query[1].lower()
                            for case in mask:
                                maskseq = maskseq[:case[0]] + maskseq.upper()[case[0]:case[1]] + maskseq[case[1]:]
                            ## ~ [1d] ~ Output seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                            OUT = open(outfile,'w')
                            OUT.write('>%s\n%s\n' % (query[0],maskseq))
                            for seq in eseqs[1:]: OUT.write('>%s\n%s\n' % seqlist.getSeq(seq))
                            OUT.close()
                            outx += 1
                        if len(eseqs) != eseqx: self.warnLog('Only %d of %d sequences output in %s %s flank files!' % (len(eseqs),eseqx,outx,qname))
                        eseqs.append(eseqs.pop(0))
                        qoutx += outx
                    self.printLog('#QFAS','%s: %s sequences files from %d query x %d flank' % (elm,rje.iStr(qoutx),eseqx,len(self.list['FlankMask'])))
                    esx += qoutx
            self.printLog('#QFAS','%s Fasta files made for %d %s Queries (%d ELMs; %d QRegions)' % (rje.iStr(esx),qx,datadir,len(elmlist),len(self.list['FlankMask'])))
            return True                                
        except: self.errorLog('Problem during %s.queryDatasets().' % self); return False  
#########################################################################################################################
    def filterDatasets(self,elmlist,datadir='ELM_Reduced',filtertxt=''):  ### Produces filtered ELM/PPI Datasets
        '''
        Produces filtered ELM/PPI Datasets.
        >> datadir:path = Path to ELM dataset files
        >> elmlist:list = List of ELM datasets to make files for
        >> filtertxt:str = Optional text to put output to filter_info.txt outlining filtering settings.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if self.dev(): self.printLog('#DEV','Dataset filter for ELMBench and PPIBench is currently disabled.')
            #return False
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ %s FILTERED DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % (datadir))
            db = self.db('ELMInstance')
            indir = rje.makePath('%s%s' % (self.getStr('GenPath'),datadir))
            outdir = rje.makePath('%s%s%s' % (self.getStr('GenPath'),datadir,self.getStr('FilterDir')))
            rje.mkDir(self,outdir,True)
            if filtertxt: open('%sfilter_info.txt' % outdir,'w').write(filtertxt)
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','dna=F'])
            ### ~ [1] ~ Generate new ELM Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qx = 0; esx = 0
            for elm in elmlist:     # Note that ELMList is already filtered by UPC.
                #for efile in glob.glob('%s%s.*fas' % (indir,elm)):
                efile  = '%s%s.fas' % (indir,elm)
                if not rje.exists(efile): self.warnLog('%s missing!' % efile); continue
                outfile = '%s%s.fas' % (outdir,rje.baseFile(efile,strip_path=True))
                rje.fileTransfer(fromfile=efile,tofile=outfile,deletefrom=False,append=False); esx += 1
                #?# Modify this code in future to output flankmasked files (all sequences masked)?
            self.printLog('#FAS','%s filtered %s*.fas files (%d ELMs).' % (rje.iStr(esx),outdir,len(elmlist)))
            return True
        except: self.errorLog('Problem during %s.filterDatasets().' % self); return False
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
            if not self.list['FlankMask']: self.list['FlankMask'] = ['']
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
                if flank in ['','none']: mask = [(0,qlen)]
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
                if flank: outfile = '%s.%s.%s.fas' % (basefile,string.split(qname)[0],flank)
                else: outfile = '%s.%s.fas' % (basefile,string.split(qname)[0])
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
            if self.getBool('RanDat'):
                core = rje_slimcore.SLiMCore(self.log,['logmask=T','maskpickle=T','resdir=%s' % self.getStr('GenPath')]+self.cmd_list+['efilter=F','slimbuild=F','walltime=0','seqin=%s' % self.getStr('RandSource')])
                core.setupSeqIn()
                randfas = rje.baseFile(self.getStr('RandSource')) + '.fas'
                seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F','seqin=%s' % randfas]
            else:
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
            reduced = '%s.reduced.motifs' % elmbase
            genpath = rje.makePath('%sSourceMatch/' % self.getStr('GenPath'))
            rje.mkDir(self,genpath)
            sourcebase = rje.baseFile(self.getStr('RandSource'),True)
            if self.getBool('SimBench'):
                if not rje.checkForFile(reduced):
                    self.printLog('#MOTIF','Reduced Motif file "%s" not found: no simulation.' % reduced)
                    self.bool['SimBench'] = False
                else:
                    core.obj['SlimList'] = slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list)
                    slimlist.loadMotifs(reduced)    # Will automatically apply minic filter
                    for slim in slimlist.slims():
                        slimfas = '%s%s.%s.fas' % (genpath,slim.getStr('Name'),sourcebase)
                        if not self.force() and rje.checkForFile(slimfas):
                            self.printLog('#FAS','%s found. (Force=F)' % slimfas)
                            continue
                        core.dict['MotifSeq'][slim.pattern()] = slimfas  #{Pattern:File}                    
            ## ~ [0d] Perform randomised data generation only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getBool('SimBench'):
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
            #self.deBug(core.dict['MotifSeq'])
            if core.dict['MotifSeq']:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~ SLIMCORE MAPPING MOTIFS TO SEQUENCES ~~~~~~~~~~~~~~~~~~ #')
                self.printLog('#MSEQ','%d MotifSeq datasets to generate.' % len(core.dict['MotifSeq']))
                #!# Add log output for progress #!#
                if self.getBool('RanDat'): #i# Should already have DAT loaded
                    coreuni = core.obj['SeqList'].obj['UniProt']
                    #i# This should match the sequences in the seqcmd fasta file, so masking should work fine
                    #!# Might need to do some kind of check.
                else:
                    core.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
                core.setupBasefile()
                core.setNum({'StartTime':time.time()})
                core.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
                if core.getBool('Masked'): core.obj['SeqList'].saveFasta(seqfile='%s%s.%s.masked.fas' % (genpath,core.basefile(),core.maskText()))
                core.motifSeq()
            ### ~ [2] Run SLiMProb on MotifSeq files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~ SLIMPROB MAPPING MOTIFS TO SEQUENCES ~~~~~~~~~~~~~~~~~~ #')
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
                    if not rje.exists(slimfas) or open(slimfas,'r').read()[:1] != '>': continue
                    sscmd = ['extras=0','maxsize=0','maxupc=0','maxseq=1000']+self.cmd_list+['resfile=%s' % ssfile,'append=T','motifs=','seqin=%s' % slimfas,'resdir=%sSLiMProb/' % searchdir,'occupc=T']
                    ss = slimprob.SLiMProb(self.log,sscmd+['debug=F','mergesplits=F'])
                    ss.setStr({'Basefile':''})
                    ss.obj['SlimList'].list['Motif'] = [slim]
                    ss.run()
            ## ~ [2b] ~ Load SLiMProb Results into Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sdb = db.addTable(string.replace(ssfile,'occ.csv','csv'),mainkeys=['Motif'],datakeys='All',name='SourceSearch')
            odb = db.addTable(ssfile,mainkeys=['Motif','Seq','Start_Pos'],datakeys='All',name='SourceOcc')
            ### ~ [3] Work through each OK ELM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~ GENERATE RANDOM DATASETS ~~~~~~~~~~~~~~~~~~ #')
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
    ### <4> ### Pfam PPI Mapping Methods                                                                                #
#########################################################################################################################
    def pfam2elm(self): ### (Generates and) Returns the pfam2elm mapping dictionary
        '''(Generates and) Returns the pfam2elm mapping dictionary.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Pfam2ELM' in self.dict and self.dict['Pfam2ELM']: return self.dict['Pfam2ELM']
            self.dict['Pfam2ELM'] = pfam2elm = {}
            ### ~ [1] Generate dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for entry in self.db('ELMDomains').entries() + self.db('ELMInteractors').entries():
                if entry['Domain'] not in pfam2elm: pfam2elm[entry['Domain']] = []
                if entry['Elm'] not in pfam2elm[entry['Domain']]: pfam2elm[entry['Domain']].append(entry['Elm'])
            return pfam2elm
        except: self.errorLog('pfam2elm error'); raise
#########################################################################################################################
    def hub2elm(self,hub,elm=None): ### Returns the ELMs mapping via Pfam to a Hub, or True/False if ELM given
        '''Returns the ELMs mapping via Pfam to a Hub, or True/False if ELM given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pfam2elm = self.pfam2elm()
            ## [0a] Check for Pfam hub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if hub in pfam2elm:
                elmlist = pfam2elm[hub]
                if elm: return elm in elmlist
                else: return rje.sortUnique(elmlist)
            ### ~ [1] Generate ELM List/Map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elmlist = []
            for spec in self.list['PPISpec']:
                ptable = self.db('Pfam.%s' % spec)
                if hub not in ptable.data(): continue
                for domain in string.split(ptable.data(hub)['Pfam'],'|'):
                    if domain in pfam2elm: elmlist += pfam2elm[domain]
            if not elmlist: self.warnLog('Cannot find hub "%s" in Pfam-ELM links!' % hub)
            if elm: return elm in elmlist
            else: return rje.sortUnique(elmlist)
        except: self.errorLog('hub2elm error'); raise
#########################################################################################################################
    ### <5> ### SLiMBench Benchmarking Methods                                                                          #
#########################################################################################################################
    def benchmark(self): ### Main SLiMBench Benchmarking run method                                                 #V1.0
        '''Main SLiMBench Benchmarking run method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setupBenchmarking(): return False
            if self.getStr('DataType') == 'occ': return self.occBenchmark()
            ### ~ [2] Load and Process (Q)SLiMFinder search data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.loadSLiMPredictions(): return False
            if self.getBool('NoAmb'): self.filterAmb(save=True)
            ### ~ [3] Perform CompariMotif search of SLiM Predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.compariMotif(): return False
            ### ~ [4] Rate SLiM Predictions as TP/FP/OT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.ratings(): return False
            ### ~ [5] Summary Performance Table (ELM only?) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.summaryTables()    # Will not run without 2+ RunID fields
            ### ~ [6] Generate Sn/Sp etc. stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.assessment(): return False
            return True
        except: self.errorLog('%s.benchmark error' % self); return False
#########################################################################################################################
    def setupBenchmarking(self):    ### Main class setup method.                                                    #V1.0
        '''Main class setup method.'''
        try:### ~ [1] Load ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.baseFile(self.getStr('BenchBase'))
            db = self.obj['DB']
            db.baseFile(self.getStr('BenchBase'))
            for table in db.tables(): db.deleteTable(table)     # Clear db before proceeding
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            ## ~ [1a] Load ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('DataType') != 'occ':
                elmc = db.addTable(self.sourceDataFile('ELMClass',force=False),mainkeys=['ELMIdentifier'],datakeys='All',name='ELMClass')
                elmc.dataFormat({'#Instances':'int'})
            if self.getStr('DataType') == 'ppi':
                elmi = db.addTable(self.sourceDataFile('ELMInteractors',force=False),mainkeys='#',datakeys='All',name='ELMInteractors') # "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
                elmd = db.addTable(self.sourceDataFile('ELMDomains',force=False),mainkeys='#',datakeys='All',name='ELMDomains') # "ELM identifier"	"Interaction Domain Id"	"Interaction Domain Description"	"Interaction Domain Name"
                elmd.renameField('ELM identifier','Elm'); elmd.renameField('Interaction Domain Id','Domain')
                #pfamdb = None
                for spec in self.list['PPISpec']:
                    #if pfamdb: self.db().mergeTables(pfamdb,self.db().addTable('%sslimbench.Pfam.%s.tdt' % (self.getStr('GenPath'),spec),mainkeys=['Uniprot'],expect=True,name=spec))
                    #else: pfamdb = self.db().addTable('%sslimbench.Pfam.%s.tdt' % (self.getStr('GenPath'),spec),mainkeys=['Uniprot'],expect=True,name='Pfam')
                    self.db().addTable('%sslimbench.Pfam.%s.tdt' % (self.getStr('SourcePath'),spec),mainkeys=['Uniprot'],expect=True,name='Pfam.%s' % spec)
            #i# See setupSourceData() method if UniProt mapping data is needed.
            ### ~ [2] Load ELM motif file for comparison ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Establish file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.sourceDataFile('ELMClass',force=False),strip_path=True))
            if self.getStr('DataType') not in ['occ','ppi'] and not rje.checkForFile(self.getStr('CompDB')): self.str['CompDB'] = '%s.reduced.motifs' % elmbase
            if not rje.checkForFile(self.getStr('CompDB')): self.str['CompDB'] = self.getStr('ELMClass')
            if not rje.checkForFile(self.getStr('CompDB')): self.printLog('#MOTIF','CompDB motif input file (and ELM class file) not found'); raise IOError
            ## ~ [2b] Load Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['CompDB'] = rje_slimlist.SLiMList(self.log,self.cmd_list+['force=F'])
            if not self.obj['CompDB'].loadMotifs(self.str['CompDB']): return False
            self.printLog('#SETUP','Benchmarking setup successful.')
            return True     # Setup successful            
        except: self.errorLog('Problem during %s setupBenchmarking.' % self); return False  # Setup failed
#########################################################################################################################
    def loadSLiMPredictions(self):  ### Load SLiM Predictions into Database object and establish runs etc.          #V1.0
        '''Load SLiM Predictions into Database object and establish runs etc.'''
        try:### ~ [0] Setup Header Splitting Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('DataType')[:3] == 'sim': dsethead = ['Motif','RType','Rep','PosNum']#,'SeqNum']
            #X#elif self.getStr('DataType') == 'ppi': dsethead = ['Motif','PPI','IType']
            elif self.getStr('DataType') == 'ppi': dsethead = ['Hub','IType']   # IType will be 'PPI' if not found
            elif self.getStr('DataType') == 'elm': dsethead = ['Motif']
            else: self.errorLog('DataType "%s" not recognised.',printerror=False); raise ValueError
            if self.getBool('Queries'): dsethead += ['Query','Region']
            for head in self.list['RunID']:
                if head in protected + dsethead:
                    self.errorLog('Cannot have protected/dataset field "%s" in RunID split.' % head, printerror=False)
                    raise ValueError
            ### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ LOAD PREDICTIONS ~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            db.setBool({'TupleKeys':True})
            db.setStr({'Delimit':self.getStr('Delimit')})
            rfile = '%s.results.%s' % (db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and os.path.exists('%s.ratings.tdt' % db.baseFile()) and self.getBool('MemSaver'):
                self.printLog('#MEM','Ratings file found. Will skip reading results. (memsaver=T; force=F)')
                return True
            if not self.force() and os.path.exists(rfile):
                if self.getStr('DataType')[:3] == 'sim':
                    return db.addTable(rfile,mainkeys=['Motif','Query','Region','RType','Rep','PosNum','SeqNum'] + self.list['RunID'] + ['Pattern'],datakeys='All',name='results',replace=True)
                else: return db.addTable(rfile,mainkeys=dsethead + self.list['RunID'] + ['Pattern'],datakeys='All',name='results',replace=True)
            ## ~ [1a] Load and merge results files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            resdb = None
            for file in self.list['ResFiles']:
                if not resdb: resdb = db.addTable(file,mainkeys=['Dataset','RunID','Pattern'],datakeys='All',name='results')
                else: db.mergeTables(resdb,db.addTable(file,mainkeys=['Dataset','RunID','Pattern'],datakeys='All',name=file))
            resdb.dropEntriesDirect('Pattern',['<','>','!'])    # Consider adding > and ! : Why keep them?!
            resdb.setStr({'Delimit':self.getStr('Delimit')})
            ## ~ [1b] Check all Datasets present for all RunIDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('DataType')[:3] == 'sim':
                dcheck = 'DCheck'
                checkfield = 'RCheck'
                resdb.addFields([dcheck,checkfield])
                for entry in resdb.entries():
                    entry['DCheck'] = string.split(entry['Dataset'],'.')
                    entry['RCheck'] = '%s|%s' % (entry['RunID'],entry['DCheck'][1])
                    entry['DCheck'] = string.join(entry['DCheck'][:1]+entry['DCheck'][2:],'.')
            else:
                dcheck = 'Dataset'
                checkfield = 'RunID'
            runidlist = resdb.indexKeys(checkfield)
            dsetlist =  resdb.indexKeys(dcheck)
            incomplete = []
            missdx = 0
            for runid in runidlist:
                rundsets = resdb.indexDataList(checkfield,runid,dcheck)
                missing = rje.listDifference(dsetlist,rundsets) # Returns the elements of list1 that are not found in list 2.
                if missing:
                    missdx += 1
                    self.warnLog('RunID "%s" has %s of %s datasets missing.' % (runid,rje.iLen(missing),rje.iLen(dsetlist)),warntype='Missing Data',quitchoice=False,suppress=False,dev=False,screen=True)
                    incomplete = rje.listUnion(incomplete,missing)
            if incomplete:
                self.warnLog('%s of %s RunIDs have datasets missing: %s of %s datasets have missing data' % (rje.iStr(missdx),rje.iLen(runidlist),rje.iLen(incomplete),rje.iLen(dsetlist)),warntype='Missing Data',quitchoice=self.getBool('Balanced'))
                if self.getBool('Balanced'):
                    self.printLog('#DATA','Balanced=T: incomplete datasets will be removed.')
                    resdb.dropEntriesDirect(dcheck,incomplete)
                    self.printLog('#DATA','Results for %s of %s detected datasets loaded for all %s RunID' % (rje.iLen(resdb.indexKeys(dcheck)),rje.iLen(dsetlist),rje.iLen(runidlist)))
                else:
                    self.warnLog('Balanced=F: incomplete datasets may result in unfair methods comparisons.')
            else: self.printLog('#DATA','Results for all %s detected datasets loaded for all %s RunID' % (rje.iLen(dsetlist),rje.iLen(runidlist)))
            if self.getStr('DataType')[:3] == 'sim':
                resdb.dropField(dcheck)
                resdb.dropField(checkfield)
            self.debug(resdb.fields())
            ## ~ [1c] Add Motif-Dataset links for PPIBench ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('DataType') == 'ppi':
                resdb.addField('Motif',evalue='PPI')
                resdb.newKey(['Dataset','Motif','RunID','Pattern'])
                for ekey in resdb.dataKeys():
                    entry = resdb.data(ekey)
                    if '.' not in entry['Dataset']:
                        if entry['Dataset'].startswith('PF'): entry['Dataset'] = '%s.domppi' % entry['Dataset']
                        else: entry['Dataset'] = '%s.ppi' % entry['Dataset']
                    hub = string.split(entry['Dataset'],'.')[0]
                    motifs = self.hub2elm(hub)
                    if not motifs:
                        self.warnLog('No ELMs mapped to %s via Pfam domains!' % hub)
                        entry['Motif'] = '!na!'
                    else:
                        entry['Motif'] = motifs[0]
                        self.debug(entry)
                        for motif in motifs[1:]:
                            self.debug(motif)
                            resdb.addEntry(rje.combineDict({'Motif':motif},entry,overwrite=False))
            ## ~ [1d] Split on RunID and report stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ex = 0.0; etot = resdb.entryNum()
            resdb.list['Fields'] = dsethead + self.list['RunID'] + resdb.list['Fields']
            newkeys = dsethead + self.list['RunID'] + ['Pattern']
            if self.getStr('DataType')[:3] == 'sim': newkeys += ['SeqNum']; dsethead.insert(4,'SeqNum')
            if self.getStr('DataType') == 'ppi': newkeys.insert(2,'Motif')
            self.printLog('#SPLIT','Dataset => %s' % string.join(dsethead))
            self.printLog('#SPLIT','RunID => %s' % string.join(self.list['RunID']))
            for entry in resdb.entries():
                self.progLog('\r#SPLIT','Splitting Dataset and RunID fields: %.2f%%' % (ex/etot)); ex += 100.0
                dsetdata = string.split(entry['Dataset'],'.')
                if self.getStr('DataType') == 'elm' and 'reduced' in dsetdata: dsetdata.remove('reduced')
                try:
                    for dhead in dsethead:
                        if dhead in ['SeqNum']: dsetdata.pop(0) # SeqNum already a results field
                        else: entry[dhead] = dsetdata.pop(0)
                    if dsetdata:
                        self.errorLog('Too many fields in dataset name! Expect: %s (%d remain.)' % (string.join(dsethead,'.'),len(dsetdata)),printerror=False)
                        raise ValueError
                except: self.errorLog('%s Dataset Name Error: %s' % (self.getStr('DataType').upper(),entry['Dataset'])); raise
                runid = string.split(entry['RunID'],'.')
                try:
                    for i in range(len(self.list['RunID'])): entry[self.list['RunID'][i]] = runid[i]
                except IndexError:
                    self.errorLog('Cannot split %s into %s: check that RunID split list and results files match!' % (entry['RunID'],string.join(self.list['RunID'],'.')))
                    return False
                except:
                    self.errorLog('Failed to split %s into %s' % (entry['RunID'],string.join(self.list['RunID'],'.')))
                    raise
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
            if not self.force() and os.path.exists(cfile):
                if self.getBool('MemSaver'): self.printLog('#COMP','CompariMotif results found. (memsaver=T; force=F)'); return True
                return db.addTable(cfile,mainkeys=['Name1','Name2'],name='comparimotif',replace=True)
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
            if self.getBool('MemSaver'): return True
            return db.addTable(cfile,mainkeys=['Name1','Name2'],name='comparimotif')
        except: self.errorLog('Problem during %s compariMotif.' % self); return False  # compariMotif failed
#########################################################################################################################
    def ratings(self):  ### Compares CompariMotif search results to Datasets and rates as TP,FP or OT               #V1.0
        '''Perform CompariMotif search of Patterns against CompDB.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MemSaver'): return self.ratingsMemSaver()
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ RATINGS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            newkeys = self.list['RunID'][0:] + ['Motif']
            if self.getBool('Queries'): newkeys += ['Query','Region']
            if self.getStr('DataType')[:3] == 'sim': newkeys += ['RType','Rep','PosNum','SeqNum']
            if self.getStr('DataType') == 'ppi': newkeys += ['IType','Hub']
            newkeys += ['Pattern']
            # Re-order ratings keys to match assessment groupings 
            akeys = self.list['RunID'][0:]
            if self.getStr('DataType')[:3] == 'sim': akeys += ['PosNum','SeqNum','RType','Rep']
            elif self.getStr('DataType') == 'ppi': akeys += ['IType']
            #akeys += ['Region','ICCut','LenCut','SigCut','ByCloud']
            if self.getBool('Queries'): akeys.append('Region')
            for field in newkeys:
                if field not in akeys: akeys.append(field)
            newkeys = akeys
            if not self.force() and os.path.exists(rfile):
                return db.addTable(rfile,mainkeys=newkeys,datakeys='All',name='results',replace=True)
            ## ~ [0a] Setup ELM-Pfam-Hub mapping for PPIBench ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStr('DataType') == 'ppi':
                pfam2elm = self.pfam2elm()
                pfamdb = self.db('Pfam')
            ## ~ [0b] Make reduced ratings table from results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
                    if rje.matchExp('^(.+_\d)$',entry[field]): match[field] = rje.matchExp('^(.+_\d)$',entry[field])[0]
                    elif rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',entry[field])[0]
                    elif rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field])[0]
                    else: match[field] = entry[field]
                tp = False
                if self.getStr('DataType') == 'ppi': tp = self.hub2elm(match['Motif'],match['Name2'])
                    #if match['Motif'] in pfam2elm:  # Pfam domain hub
                    #    if match['Name2'] in pfam2elm[match['Motif']]: tp = True
                    #else:
                    #    pentry = pfamdb.data(match['Motif'])
                    #    if pentry:
                    #        for domain in string.split(pentry['Pfam'],'|'):
                    #            if domain in pfam2elm and match['Name2'] in pfam2elm[domain]: tp = True
                    #    else: self.warnLog('Cannot find hub "%s" in Pfam-ELM links!' % match['Motif'])
                elif match['Motif'] == match['Name2']: tp = True
                if tp: entry['Rating'] = 'TP'
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
                #tkey= '%s\t%s\t%s' % (entry['Motif'],entry['Pattern'],entry['Motif'])
                tkey= tdb.makeKey({'Motif':entry['Motif'],'Pattern':entry['Pattern'],'CMHit':entry['Motif']})
                if self.getStr('DataType') == 'ran' or (self.getStr('DataType') == 'sim' and entry['RType'] == 'ran'):
                    if tkey in tdb.data(): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'TN'
                    else: entry['Rating'] = 'FP'
                else:
                    if tkey in tdb.data(): entry['Rating'] = tdb.data(tkey)['Rating']
                    elif self.getStr('DataType') == 'ppi' and entry['Pattern'] in tdb.indexDataList('Motif',entry['Motif'],'Pattern'):
                        for tentry in tdb.indexEntries('Motif',entry['Motif']):
                            if tentry['Pattern'] == entry['Pattern']: entry['Rating'] = tentry['Rating']
                            else: continue
                            if entry['Rating'] == 'TP': break
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'FN'
                    else: entry['Rating'] = 'FP'
            self.printLog('\r#RATE','Rating SLiM Predictions complete.')
            for rating in ['TP','OT','FP','FN','TN']:
                if rating in rdb.index('Rating'): self.printLog('#%s' % rating,'%s %s results' % (rje.iLen(rdb.index('Rating')[rating]),rating))
                else: self.printLog('#%s' % rating,'No %s results' % (rating))
            ### ~ [2] Tidy up table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for field in ['Dataset','RunID','Chance','RunTime']:
                if field in rdb.fields() and (self.i() < 1 or rje.yesNo('Remove field "%s" from Ratings table?' % field)): rdb.dropField(field)
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
            ## Make compatible with memsaver 
            rdb.list['Keys'] = newkeys
            newfields = rdb.list['Keys'][0:]
            for field in rdb.fields():
                if field not in newfields: newfields.append(field)
            rdb.list['Fields'] = newfields
            rdb.remakeKeys()
            ## Output
            rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during %s ratings.' % self); return False  # rating failed
#########################################################################################################################
    def ratingsMemSaver(self):  ### Compares CompariMotif search results to Datasets and rates as TP,FP or OT       #V2.0
        '''Perform CompariMotif search of Patterns against CompDB.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ MEMSAVER RATINGS ~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            newkeys = self.list['RunID'][0:] + ['Motif']
            if self.getBool('Queries'): newkeys += ['Query','Region']
            if self.getStr('DataType')[:3] == 'sim': newkeys += ['RType','Rep','PosNum','SeqNum']
            if self.getStr('DataType') == 'ppi': newkeys.remove('Motif'); newkeys += ['IType','Hub']
            newkeys += ['Pattern']
            # Re-order ratings keys to match assessment groupings 
            akeys = self.list['RunID'][0:]
            if self.getStr('DataType')[:3] == 'sim': akeys += ['PosNum','SeqNum','RType','Rep']
            elif self.getStr('DataType') == 'ppi': akeys += ['IType']
            #akeys += ['Region','ICCut','LenCut','SigCut','ByCloud']
            if self.getBool('Queries'): akeys.append('Region')
            for field in newkeys:
                if field not in akeys: akeys.append(field)
            newkeys = akeys
            if not self.force() and os.path.exists(rfile):
                rdb = db.addTable(rfile,mainkeys=newkeys,datakeys='All',name='results',replace=True)
                #self.deBug(rdb.keys())
                return rdb
            ## ~ [0a] Setup ELM-Pfam-Hub mapping for PPIBench ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #if self.getStr('DataType') == 'ppi':
            #    pfam2elm = self.pfam2elm()
                #for entry in self.db('ELMDomains').entries() + self.db('ELMInteractors').entries():
                #    if entry['Domain'] not in pfam2elm: pfam2elm[entry['Domain']] = []
                #    if entry['Elm'] not in pfam2elm[entry['Domain']]: pfam2elm[entry['Domain']].append(entry['Elm'])
                #self.debug(pfam2elm['PF02747'])
            #    pfamdb = self.db('Pfam')
            ## ~ [0b] Make reduced ratings table from results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rdb = self.db('results')
            rdb.index('Pattern')
            ### ~ [1] Find True Positive Patterns from CompariMotif hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tfields = ['Motif','Pattern','CMHit','MatchPos','MatchIC','NormIC','Rating']
            tkeys = ['Motif','Pattern','CMHit']
            tfile = '%s.tp.tdt' % self.getStr('BenchBase')
            if not self.force() and rje.exists(tfile):
                tdb = db.addTable(tfile,mainkeys=tkeys,name='tp',replace=True)
            else:
                tdb = db.addEmptyTable('tp',['Motif','Pattern','CMHit','MatchPos','MatchIC','NormIC','Rating'],['Motif','Pattern','CMHit'])
                for field in []: tdb.addField(field,evalue=0)                                   
                cfile = '%s.compare.tdt' % self.getStr('BenchBase')
                CM = open(cfile,'r')
                CM.seek(0,2)
                cend = float(CM.tell())
                CM.seek(0)
                chead = cline = rje.readDelimit(CM.readline())
                pi = chead.index('Name1')   # Pattern index pos
                mi = chead.index('Name2')   # Motif index pos
                pos = chead.index('MatchPos'); posx = 0
                ic = chead.index('MatchIC'); icx = 0
                norm = chead.index('NormIC'); normx = 0
                qualx = 0
                tpx = 0; otx = 0
                ## ~ [1a] Process CompariMotif results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while CM:
                    self.progLog('\r#RATE','Mapping CompariMotif results: %.2f%%' % (100.0*CM.tell()/cend))
                    cline = rje.readDelimit(CM.readline())
                    if not cline: CM.close(); break
                    try: cdata = {'CMHit':cline[mi],'Pattern':cline[pi],'MatchPos':int(cline[pos]),'MatchIC':float(cline[ic]),'NormIC':float(cline[norm])}
                    except: self.deBug(cline); continue
                    #self.deBug(cdata)
                    if cdata['MatchPos'] < 2: posx += 1; continue   # Must have 2+ positions match
                    if cdata['MatchIC'] < 1.5: icx += 1; continue   # Must match 1 fixed + [3] or 2x[2]
                    if cdata['NormIC'] < 0.5: normx += 1; continue   # Must match half smallest motif
                    match = {}
                    field = 'CMHit'
                    if rje.matchExp('^(.+_\d)$',cdata[field]): match[field] = rje.matchExp('^(.+_\d)$',cdata[field])[0]
                    elif rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',cdata[field]): match[field] = rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',cdata[field])[0]
                    elif rje.matchExp('^(.+)_[a-z][a-z0-9]*$',cdata[field]): match[field] = rje.matchExp('^(.+)_[a-z][a-z0-9]*$',cdata[field])[0]
                    else: match[field] = cdata[field]
                    for entry in rdb.indexEntries('Pattern',cdata['Pattern']):
                        field = 'Motif'
                        if rje.matchExp('^(.+_\d)$',entry[field]): match[field] = rje.matchExp('^(.+_\d)$',entry[field])[0]
                        elif rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+_\d)_[a-z][a-z0-9]*$',entry[field])[0]
                        elif rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field]): match[field] = rje.matchExp('^(.+)_[a-z][a-z0-9]*$',entry[field])[0]
                        else: match[field] = entry[field]
                        #self.deBug(entry); self.deBug(match)
                        tp = False
                        #if self.getStr('DataType') == 'ppi': tp = self.hub2elm(match['Motif'],match['CMHit'])
                            #if match['Motif'] in pfam2elm:  # Pfam domain hub
                            #    if match['Name2'] in pfam2elm[match['Motif']]: tp = True
                            #else:
                            #    pentry = pfamdb.data(match['Motif'])
                            #    if pentry:
                            #        for domain in string.split(pentry['Pfam'],'|'):
                            #            if domain in pfam2elm and match['CMHit'] in pfam2elm[domain]: tp = True
                            #    else: self.warnLog('Cannot find hub "%s" in Pfam-ELM links!' % match['Motif'])
                        #el
                        if match['Motif'] == match['CMHit']: tp = True
                        if tp: tdb.addEntry(rje.combineDict({'Motif':entry['Motif'],'Rating':'TP'},cdata),warn=False); tpx += 1
                        else:
                            ## Second, stronger filter for the OT motifs? 
                            if cdata['MatchIC'] < 2.5 and cdata['NormIC'] < 1: qualx += 1; continue
                            tdb.addEntry(rje.combineDict({'Motif':entry['Motif'],'Rating':'OT'},cdata),overwrite=False); otx += 1
                self.printLog('\r#RATE','Mapped %s CompariMotif hits: %s TP; %s OT.' % (rje.iStr(tdb.entryNum()),rje.iStr(tpx),rje.iStr(otx)))
                self.printLog('#FILT','Filtered %s single position hits' % rje.iStr(posx))   # Must have 2+ positions match
                self.printLog('#FILT','Filtered %s low MatchIC (< 1.5) hits' % rje.iStr(icx))   # Must match 1 fixed + [3] or 2x[2]
                self.printLog('#FILT','Filtered %s low NormIC (< 0.5) hits' % rje.iStr(normx))   # Must match half smallest motif
                self.printLog('#FILT','Filtered %s low quality (MatchIC < 2.5 & NormIC < 1) OT hits' % rje.iStr(qualx))   # Must match half smallest motif
                tdb.index('Pattern')
                try: self.printLog('#TP','%s TP CompariMotif combinations' % rje.iLen(tdb.index('Rating')['TP']))
                except: self.printLog('#TP','No TP CompariMotif combinations found!')
                tdb.saveToFile()
            ### ~ [2] Rate Motifs using CompariMotif hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb.list['Fields'].append('Rating')
            ratex = {}
            for rating in ['TP','OT','FP','FN','TN']: ratex[rating] = 0
            ex = 0.0; etot = rdb.entryNum()
            for entry in rdb.entries():
                self.progLog('\r#RATE','Rating SLiM Predictions: %.1f%%' % (ex/etot)); ex += 100.0
                #tkey= '%s\t%s\t%s' % (entry['Motif'],entry['Pattern'],entry['Motif'])
                tkey= tdb.makeKey({'Motif':entry['Motif'],'Pattern':entry['Pattern'],'CMHit':entry['Motif']})
                if self.getStr('DataType') == 'ran' or (self.getStr('DataType') == 'sim' and entry['RType'] == 'ran'):
                    if tkey in tdb.data(): entry['Rating'] = tdb.data(tkey)['Rating']   #'OT'
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'TN'
                    else: entry['Rating'] = 'FP'
                else:
                    if tkey in tdb.data(): entry['Rating'] = tdb.data(tkey)['Rating']
                    elif self.getStr('DataType') == 'ppi' and entry['Pattern'] in tdb.indexDataList('Motif',entry['Motif'],'Pattern'):
                        for tentry in tdb.indexEntries('Motif',entry['Motif']):
                            if tentry['Pattern'] == entry['Pattern']: entry['Rating'] = tentry['Rating']
                            else: continue
                            if entry['Rating'] == 'TP': break
                    elif entry['Pattern'] in tdb.index('Pattern'): entry['Rating'] = 'OT'
                    elif entry['Pattern'] in '!<>-': entry['Rating'] = 'FN'
                    else: entry['Rating'] = 'FP'
                ratex[entry['Rating']] += 1
            self.printLog('\r#RATE','Rating SLiM Predictions complete.')
            for rating in ['TP','OT','FP','FN','TN']:
                if ratex[rating]: self.printLog('#%s' % rating,'%s %s results' % (rje.iStr(ratex[rating]),rating))
                else: self.printLog('#%s' % rating,'No %s results' % (rating))
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
            if len(rdb.list['Keys']) != len(newkeys): self.deBug(rdb.keys()); self.deBug(newkeys); raise ValueError
            rdb.list['Keys'] = newkeys
            newfields = rdb.list['Keys'][0:]
            for field in rdb.fields():
                if field not in newfields: newfields.append(field)
            rdb.list['Fields'] = newfields
            rdb.remakeKeys()
            #self.deBug(rdb.keys())
            rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during %s ratings.' % self); return False  # rating failed
#########################################################################################################################
    def summaryTables(self):    ### Produces summary table of best TP results                                       #V1.5
        '''Reduces results to those with minimum IC or greater and relevant SLiMLen.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(self.list['RunID']) < 2:
                self.printLog('#SUMM','Cannot produced summary table without 2+ elements to RunID')
                return False
            analfield = self.list['RunID'][1]   #'Analysis'  # Formerly 'Masking'.
            progfield = self.list['RunID'][0]       # May need to change if not used in RunID
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SUMMARY TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            sdb = self.db('results')
            for field in [analfield,progfield]:
                if field not in sdb.fields():
                    self.printLog('#SUMM','Cannot produced summary table without RunID field "%s" in compiled results.' % field)
                    return False
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            rkeys = sdb.keys()
            sfile = '%s.summary.tdt' % (self.getStr('BenchBase'))
            if not self.force() and rje.exists(sfile): return self.printLog('#SUMM','Summary table %s found (force=F)' % sfile)
            # Motif	Query	Region	Program	Masking	Pattern	Build	SeqNum	Rank	Sig	IC	Support	Coverage	Enrichment	Cloud	CloudCoverage	IUP_mean	SA_mean	Rating
            elmbase = '%s%s' % (self.getStr('GenPath'),rje.baseFile(self.getStr('ELMClass'),strip_path=True))
            reduced = '%s.reduced.motifs' % elmbase
            slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list+['force=F'])
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
            #!# Add the extra fields here for non-ELM data #!#
            if self.getBool('Queries'):
                sdb.compress(['Motif','Query','Region',progfield,analfield],rules={'Sig':'min','Enrichment':'max','Q':'max','N':'max'},default='str',best=['Sig','IC','Enrichment','Q'])
                qdb = db.copyTable(sdb,'qtemp'); sdb.dropFields(['N'])
                qdb.dropFields(['Pattern','Rank','Sig','Coverage','CloudCoverage','Support'])
                sdb.compress(['Motif','Region',progfield,analfield],rules={'Sig':'min','Enrichment':'max'},default='str',best=['Sig','IC','Enrichment'])
                qdb.compress(['Motif','Region',progfield,analfield],rules={'Sig':'min','Enrichment':'max','Q':'sum','N':'sum'},default='mean')
                for entry in qdb.entries():
                    entry['Q'] = float(entry['Q']) / entry['N']
                    sdb.data(sdb.makeKey(entry))['Q'] = entry['Q']
                sdb.dropFields(['IC','Enrichment','Rating','Query'])
            else:
                sdb.compress(['Motif',progfield,analfield],rules={'Sig':'min','Enrichment':'max'},default='str',best=['Sig','IC','Enrichment'])
                sdb.dropFields(['IC','Enrichment','Rating'])
            ### ~ [3] Reshap Wide ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.reshapeWide(progfield,['Pattern','Rank','Sig','Q','Coverage','CloudCoverage','Support'])
            ### ~ [4] Add ELM Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Queries'): sdb_newkeys = ['Motif','Region',analfield]
            else: sdb_newkeys = ['Motif',analfield]
            sdb = db.joinTables(name='summary',join=[(sdb,'Motif'),(edb,'ELMIdentifier',['Regex','Reduced'])],newkey=sdb_newkeys,cleanup=True,delimit='\t',empties=True,check=False,keeptable=True)
            sdb.dropField('summ_temp_Masking')
            #self.deBug(sdb.fields())
            sfields = sdb.fields()[0:]
            for field in sdb_newkeys + ['Reduced','Regex']: sfields.remove(field)
            sdb_newkeys.reverse()
            sdb.list['Fields'] = sdb_newkeys + ['Regex','Reduced'] + sfields
            sdb.saveToFile(sfile)
            db.addTable(rfile,mainkeys=rkeys,datakeys='All',name='results',replace=True)
            return True
        except: self.errorLog('Problem during %s summaryTables.' % self); return False  # failed
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
    def assessment(self):   ### Perform benchmarking assessment for results                                         #V2.0
        '''Perform benchmarking assessment for results.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ASSESSMENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.db()
            afile = '%s.assessment.tdt' % self.getStr('BenchBase')
            akeys = self.list['RunID'][0:]
            if self.getStr('DataType')[:3] == 'sim': akeys += ['PosNum','SeqNum']
            elif self.getStr('DataType') == 'ppi': akeys += ['IType']
            akeys += ['Region','ICCut','LenCut','ByCloud','SigCut']
            if not self.getBool('Queries'): akeys.remove('Region')
            setkeys = akeys[:-4]
            #if self.getStr('DataType')[:3] == 'sim': setkeys += ['RType','Rep']
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            rkeys = self.db('results').keys()
            db.deleteTable('results')
            self.list['SigCut'].sort()
            self.list['SigCut'].reverse()   # Need big -> small for repeated filtering
            #self.deBug(self.list['SigCut'])
            ## ~ [0a] Read in list of analysed searches or clear file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.force() and os.path.exists(afile):
                self.printLog('#BENCH','%s found! (force=False)' % afile)
                adb = db.addTable(afile,mainkeys=akeys,datakeys='All',name='assessment',replace=True)
                append = True
            else:
                rje.backup(self,afile,appendable=False); append = False
                ahead = akeys + ['TP','OT','FP','FN','TN','SN','FPX','FPXn', 'N','ResNum','DsetNum''PPV','PPVp','PPVn']
                adb = db.addEmptyTable('assessment',ahead,keys=akeys)
            ## ~ [0b] ByMotif file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bfile = '%s.bymotif.tdt' % self.getStr('BenchBase')
            dkeys = ['Motif'] + akeys[:-3]
            if self.force(): rje.backup(self,bfile,appendable=False)

            ### ~ [1] Read in Results, one search at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rend = rje.endPos(filename=rfile)
            rdb = db.openTable(rfile,mainkeys=rkeys,name='results')
            #self.deBug(rkeys)
            #self.deBug(rdb.fields())
            #self.deBug(setkeys)
            acheck = '#%s#' % string.join(setkeys+['ICCut','LenCut','ByCloud'],'#|#')
            adb.index(acheck,make=True)
            while rdb.readSet(setkeys):   ### Read in next set of data for assessment
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ASSESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
                try: rx = (rend - rdb.obj['File'].tell()) * 100.0
                except: rx = 0.0    # Should be end of file
                matchlist = string.join(rdb.list['MatchData'],'|')
                self.printLog('#RES','Results set %s read for analysis. %.2f%% remains.' % (matchlist,rx/rend))
                for ic_cut in self.list['ICCut']:
                    self.num['MinResIC'] = ic_cut
                    for len_cut in self.list['SLiMLenCut']:
                        self.int['ResLen'] = len_cut
                        for bycloud in self.list['ByCloud']:
                            matchlist = string.join(rdb.list['MatchData']+['%s' % ic_cut,'%s' % len_cut,'%s' % bycloud],'|')
                            if matchlist in adb.index(acheck):
                                self.printLog('#RES','Results set %s found in %s: skipping (force=F).' % (matchlist,afile))
                                continue
                            db.copyTable('results','assess')
                            if not self.filterMinICLen(ic_cut,len_cut,save=self.getBool('Test') or self.getBool('DeBug') or self.v()>1,table='assess'): return False
                            for sig_cut in self.list['SigCut']:
                                if not self.assessSearchMemSaver(sig_cut,bycloud,test=self.getBool('Test'),append=append): return False
                                append = True
            self.db().deleteTable('results')
        except: self.errorLog('Problem during SLiMBench assessment.'); return False  # assessment failed
#########################################################################################################################
    def filterMinICLen(self,minic,slimlen,save=False,table='results'):  ### Reduces results to those with minimum IC or greater     #V1.4
        '''Reduces results to those with minimum IC or greater and relevant SLiMLen.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ MinIC/SLiMLen FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            rfile = '%s.l%d.i%2f.tdt' % (self.getStr('BenchBase'),slimlen,minic)
            rdb = self.db(table)
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
                    if self.getStr('DataType') == 'ran' or (self.getStr('DataType') == 'sim' and entry['RType'] == 'ran'): negrate = 'TN'
                    else: negrate = 'FN'
                    if entry['Pattern'] in ['>','<','!']:
                        crapx += 1; elen_ok = False; entry['Pattern'] = '-'; entry['Rating'] = negrate
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
                        else: entry['Pattern'] = '-'; entry['Cloud'] = ''; entry['Rank'] = 0; entry['Rating'] = negrate
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
    def assessSearchMemSaver(self,sig_cut=0.05,cloud=False,append=True,test=False):   ### Calculates statistics for a search      #V2.0
        '''
        Calculates statistics for different searches.
        >> sig_cut:float [0.05] = Significance cut-off for analysis
        >> cloud:bool [False] = Whether to calculate stats on clouds (True) or individual motifs (False)
        >> append:bool [True] = Whether to append assessment file
        >> test:bool [False] = Special test option for downloading intermediate states.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~ ASSESS SEARCH (IC=%.2f; Len=%d) ~~~~~~~~~~~~~~~~~~~~~~~ #' % (self.num['MinResIC'],self.int['ResLen']))
            self.printLog('#SIG','Significance cut-off = %s; ByCloud = %s; Queries = %s' % (sig_cut,cloud,self.getBool('Queries')))
            db = self.obj['DB']
            adb = db.getTable('assess')   # Single data combination without sig filtering. Can use for each filter.
            if self.v() < 2: adb.setInt({'Verbose':-1})
            ## ~ [0a] Determine headers that will identify a unique benchmarking run ~~~~~~~~~~~~~~ ##
            keep_unique = self.list['RunID'][0:][0:]
            if self.getStr('DataType')[:3] == 'ppi': keep_unique += ['IType']
            if self.getBool('Queries'): keep_unique += ['Region']
            keep_unique += ['ICCut','LenCut','SigCut','ByCloud']
            if self.getStr('DataType')[:3] == 'sim': keep_unique += ['Rep','PosNum','SeqNum','RType']
            elif self.getStr('DataType')[:3] == 'ppi': keep_unique += ['Hub']
            ## ~ [0b] Significance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            prex = adb.entryNum(); dumpx = 0
            for entry in adb.entries():
                if float(entry['Sig']) > sig_cut and int(entry['Rank']) > 0:
                    if int(entry['Rank']) > 1: adb.dropEntry(entry)
                    else:
                        if self.getStr('DataType') == 'ran' or (self.getStr('DataType') == 'sim' and entry['RType'] == 'ran'): negrate = 'TN'
                        else: negrate = 'FN'
                        entry['Cloud'] = ''; entry['Rating'] = negrate
                    dumpx += 1
            self.printLog('#SIG','Significance filtering > %s: %s filtered' % (sig_cut,rje.iStr(dumpx + prex - adb.entryNum())))
            adb = db.copyTable('assess','assess_sig')
            ## ~ [0c] Setup fields for compression ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            adb.addField('ICCut',evalue=self.num['MinResIC'],log=self.v()>0)
            adb.addField('LenCut',evalue=self.int['ResLen'],log=self.v()>0)
            adb.addField('SigCut',evalue=sig_cut,log=self.v()>0)
            adb.addField('ByCloud',evalue='%s' % cloud,log=self.v()>0)
            for rating in ['TP','OT','FP','FN','TN','SN','FPX','FPXn']: adb.addField(rating,evalue=0,log=self.v()>0)
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
                    if entry['TP']: entry['OT'] = entry['FP'] = entry['FN'] = entry['TN'] = 0
                    elif entry['OT']: entry['FP'] = entry['FN'] = entry['TN'] = 0
                if test: adb.saveToFile('test.cloud.tdt')
            ### ~ [2] Compress by Query (by Dataset) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Setup Dataset stats and results counter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##  
            adb.addField('N',evalue=0)  #i# N is counting the number of results contributing to each level
            for entry in adb.entries():
                entry['SN'] = entry['TP']
                entry['FPX'] = entry['FP']
                entry['FPXn'] = max(entry['OT'],entry['FP'])
                entry['N'] = max(entry['TP'],entry['OT'],entry['FP'])
                if self.getStr('DataType') == 'sim':
                    if entry['RType'] == 'sim':     # Do not count FPX
                        entry['FPX'] = entry['FPXn'] = 0
                        #if not entry['TP']: entry['FN'] = 1
                    else:                           # Do not count TP/OT/FP for PPV
                        entry['SN'] = entry['TP'] = entry['OT'] = entry['FP'] = 0
                elif self.getStr('DataType') == 'simonly':
                    if entry['RType'] == 'sim':
                        pass
                        #if not entry['TP']: entry['FN'] = 1
                    else:                           # Do not count TP/OT/FP for PPV
                        entry['FPX'] = entry['FPXn'] = 0    # Do not count FPX
                        entry['SN'] = entry['TP'] = entry['OT'] = entry['FP'] = 0
            #adb.dropField('TN')
            #adb.dropField('FN')     # FN is not useful, so dropping in V2.2
            if test: adb.saveToFile('test.%s.ratingcheck.tdt' % self.getStr('DataType'))
            ## ~ [2b] Compress motifs/clouds to a single entry per dataset ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            comp_rules = {'N':'sum','SN':'max','FPX':'max','FPXn':'max','FN':'max','TN':'max'}
            #self.deBug(keep_unique)
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
                    normsum = float(sum([entry['TP'],entry['OT'],entry['FP']]))                 # Removed FN
                    if normsum:
                        for rating in ['TP','OT','FP']: entry[rating] = entry[rating] / normsum     # Removed FN
            if test: adb.saveToFile('test.conv.tdt')
            ### ~ [3] Compress queries/replicates to a single entry per ELM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            comp_rules = {'N':'sum','ResNum':'sum','SeqNum':'mean'}
            if self.getStr('DataType')[:3] == 'sim': keep_unique.remove('Rep')
            #self.deBug(keep_unique)
            adb.compress(['Motif'] + keep_unique,default='mean',rules=comp_rules)
            if self.getStr('DataType')[:3] == 'sim':    # Combine the 'ran' and 'sim' data
                keep_unique.remove('RType'); 
                adb.compress(['Motif'] + keep_unique,default='sum')
                adb.dropField('RType'); adb.dropField('Rep')
            elif self.getStr('DataType') == 'ppi':
                #!# Modify this at some point. Currently pointless as 'Motif' = 'Hub'!
                keep_unique.remove('Hub') # Combine different PPI datasets for an ELM
                adb.compress(['Motif'] + keep_unique,default='mean',rules=comp_rules)
                adb.dropField('Hub')
            if 'Query' in adb.fields(): adb.dropField('Query')
            adb.addField('DsetNum')  #i# Total number of datasets (motifs/clouds)
            for entry in adb.entries():
                if self.getBool('Queries') and self.getStr('DataType') == 'elm' and entry['SeqNum'] != entry['N']: self.printLog('#ERR','MISSING DATA! %s' % adb.makeKey(entry))
                # Re-normalise to proportions of results
                normsum = float(sum([entry['TP'],entry['OT'],entry['FP']]))                     # Removed FN
                if normsum:
                    for rating in ['TP','OT','FP']: entry[rating] = entry[rating] / normsum     # Removed FN
                elif not entry['FN']: self.printLog('#ERR','Rating error for: %s - TP+OT+FP+FN=0!' % entry); raise ValueError    # Something has gone v. wrong
                entry['DsetNum'] = entry['N']
                entry['N'] = 1
            if test: adb.saveToFile('test.motifs.tdt')
            #self.deBug(adb.entries()[0])
            adb.saveToFile('%s.bymotif.tdt' % self.baseFile(),backup=False,append=append)
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
            adb.saveToFile('%s.assessment.tdt' % self.baseFile(),backup=False,append=append)
            return True
        except: self.errorLog('Problem during %s assessSearch.' % self); return False  # assessment failed
#########################################################################################################################
    ### <6> ### OccBench Benchmarking Methods                                                                          #
#########################################################################################################################
    def occBenchmark(self): ### OccBench SLiMBench Benchmarking run method                                          #V2.4
        '''OccBench SLiMBench Benchmarking run method.'''
        try:### ~ [1] Load and Process SLiMProb search data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.loadSLiMOccResults(): return False
            ### ~ [2] Rate SLiM Predictions as TP/FP/OT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.occRatings(): return False
            ### ~ [3] Generate Sn/Sp etc. stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.occAssessment(): return False
            return True
        except: self.errorLog('%s.benchmark error' % self); return False
#########################################################################################################################
    def loadSLiMOccResults(self):  ### Load SLiM Predictions into Database object and establish runs etc.          #V1.0
        '''Load SLiM Predictions into Database object and establish runs etc.'''
        try:### ~ [0] Setup Header Splitting Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mainkeys = ['Motif','Seq','Start_Pos','End_Pos','RunID']
            for head in self.list['RunID']:
                if head in protected + mainkeys:
                    self.errorLog('Cannot have protected field "%s" in RunID split.' % head, printerror=False)
                    raise ValueError
            occfilter = []
            if self.getStrLC('OccFilter'): occfilter.append(self.getStr('OccFilter'))
            ### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ LOAD OCC DATA ~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            db.setStr({'Delimit':self.getStr('Delimit')})
            rfile = '%s.results.%s' % (db.baseFile(),rje.delimitExt(db.getStr('Delimit')))
            if not self.force() and os.path.exists('%s.ratings.tdt' % db.baseFile()) and self.getBool('MemSaver'):
                self.printLog('#MEM','Ratings file found. Will skip reading results. (memsaver=T; force=F)')
                return True
            if not self.force() and os.path.exists(rfile):
                return db.addTable(rfile,mainkeys=mainkeys[:-1] + self.list['RunID'],datakeys='All',name='results',replace=True)
            ## ~ [1a] Load and merge results files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            resdb = None    #i# NOTE: Currently only mainkeys fields read in.
            #!# Will want to add filter scores at some point?
            # occfilter=X : Maybe just one field to use in place of sigcut
            # invoccfilter=T/F  : Whether to inverse and use high score as good.
            for rfile in self.list['ResFiles']:
                if not resdb: resdb = db.addTable(rfile,mainkeys=mainkeys,datakeys=mainkeys+occfilter,name='results')
                else: db.mergeTables(resdb,db.addTable(rfile,mainkeys=mainkeys,datakeys=mainkeys+occfilter,name=rfile))
            resdb.setStr({'Delimit':self.getStr('Delimit')})
            ## ~ [1b] Split on RunID and report stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ex = 0.0; etot = resdb.entryNum()
            resdb.list['Fields'] += self.list['RunID']
            newkeys = mainkeys[:-1] + self.list['RunID']
            for entry in resdb.entries():
                self.progLog('\r#SPLIT','Splitting RunID field: %.2f%%' % (ex/etot)); ex += 100.0
                runid = string.split(entry['RunID'],'.')
                try:
                    for i in range(len(self.list['RunID'])): entry[self.list['RunID'][i]] = runid[i]
                except IndexError: self.errorLog('Check that RunID split list and results files match!'); return False
                except: raise
            self.printLog('\r#SPLIT','Split RunID fields into %s.' % string.join(newkeys,', '))
            for field in ['RunID'] + newkeys: resdb.index(field)
            resdb.newKey(newkeys,True)
            resdb.saveToFile()
            return True
        except: self.errorLog('Problem during %s loadSLiMPredictions.' % self); return False  # Setup failed
#########################################################################################################################
    def occRatings(self):   ### Combines full SLiMProb and loaded results and rates TP, FP, TN and FN.              #V2.4
        '''Combines full SLiMProb and loaded results and rates TP, FP, TN and FN.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~ OCC RATINGS ~~~~~~~~~~~~~~~~~~~ #')
            db = self.obj['DB']
            ## ~ [0a] Check for file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            rfields = self.list['RunID'][0:] + ['Motif','Seq','Start_Pos','End_Pos','Rating']
            if not self.force() and os.path.exists(rfile):
                rdb = db.addTable(rfile,mainkeys=rfields[:-1],datakeys='All',name='results',replace=True)
                return rdb
            rdb = db.getTable('results')
            rdb.addField('Rating',evalue='X')   # Will need to screen out ambiguous hits that were screened out from ELM
            ## ~ [0b] Load ELM Rating file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            odb = db.addTable(self.getStr('OccBenchPos'),mainkeys=['Motif','Seq','Start_Pos','End_Pos'],datakeys='All',name='occpos',replace=True)
            if not odb: raise IOError('OccBenchPos file %s not found!' % self.getStr('OccBenchPos'))
            ## ~ [0c] Convert results Motifs to account for split ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elmconv = {}
            for entry in rdb.entries():
                if entry['Motif'] in odb.index('Motif'): continue
                if entry['Motif'] in elmconv: entry['Motif'] = elmconv[entry['Motif']]; continue
                elmcore = string.join(string.split(entry['Motif'],'_')[:-1],'_')
                #self.debug('%s: %s' % (elmcore,elmcore in elmi.index('ELMIdentifier')))
                if elmcore in odb.index('Motif'): entry['Motif'] = elmconv[entry['Motif']] = elmcore
            prex = rdb.entryNum()
            rdb.remakeKeys(warnings=False)
            if rdb.entryNum() != prex:
                self.printLog('#SPLIT','Split motifs combined: %s redundant occurrences removed.' % (rje.iStr(prex-rdb.entryNum())))
            ## ~ [0d] Load/Filter ELMs and reduce tables accordingly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            odb.dropEntriesDirect('Motif',self.obj['CompDB'].nameList(remsplit=True),inverse=True)
            rdb.dropEntriesDirect('Motif',self.obj['CompDB'].nameList(remsplit=True),inverse=True)
            prex = rdb.entryNum()
            rdb.dropEntriesDirect('Motif',odb.indexKeys('Motif'),inverse=True)
            if rdb.entryNum() != prex:
                self.printLog('#BENCH','Removed %s occurrences: motif not in OccBenchPos file.' % (rje.iStr(prex-rdb.entryNum())))
            for motif in odb.indexKeys('Motif'):
                #!# Could consider filtering out motifs on the basis of maxocc=X setting? And minocc?
                if motif in rdb.indexKeys('Motif'):
                    self.printLog('#OCC','%s: %s OccBenchPos; %s results (combined).' % (motif, rje.iLen(odb.index('Motif')[motif]), rje.iLen(rdb.index('Motif')[motif])),screen=False)
                else:
                    self.printLog('#OCC','%s: %s OccBenchPos; No results (combined)!' % (motif, rje.iLen(odb.index('Motif')[motif])),screen=False)
                    self.warnLog('%s (%s OccBenchPos) has no results (combined)! Check it was used in searches.' % (motif, rje.iLen(odb.index('Motif')[motif])))
            missing = []
            for motif in rdb.indexKeys('Motif'):
                if motif not in odb.indexKeys('Motif'):
                    self.printLog('#MOTIF','Results motif "%s" not found in OccBenchPos: check input options/filtering. (Will be removed.)' % motif)
                    missing.append(motif)
            if missing: rdb.dropEntriesDirect('Motif',missing)
            ### ~ [1] Combine OccBench ratings with results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for runid in rdb.index('RunID'):
                allocc = odb.dataKeys()     # Full list of FP and TP occurrences
                px = 0; nx = 0; xx = 0; rentry = None; ptot = len(allocc); rx = 100.0
                for entry in rdb.indexEntries('RunID',runid):
                    self.progLog('\r#RATE','Rating %s occurrences: %.2f%%' % (runid,rx/ptot)); rx += 100.0
                    rentry = entry
                    ekey = odb.makeKey(entry)
                    if ekey in allocc: entry['Rating'] = odb.data(ekey)['Rating']; allocc.remove(ekey); px += 1
                    else:
                        #X#self.deBug('"Ambiguous" result: %s' % ekey)
                        xx += 1; ptot += 1
                for ekey in allocc:
                    self.progLog('\r#RATE','Rating %s occurrences: %.2f%%' % (runid,rx/ptot)); rx += 100.0
                    pentry = odb.data(ekey)                 # Positive occurrence missing from results (i.e. Negative)
                    nentry = rje.combineDict({},rentry)     # New negative entry
                    nentry['Rating'] = {'FP':'TN','TP':'FN'}[pentry['Rating']]
                    for field in ['Motif','Seq','Start_Pos','End_Pos']: nentry[field] = pentry[field]
                    rdb.addEntry(nentry); nx += 1
                self.printLog('#RATE','%s occurrences rated: %s positives; %s negatives; %s ambiguous (to remove).' % (runid,rje.iStr(px),rje.iStr(nx),rje.iStr(xx)))
            prex = rdb.entryNum()
            rdb.dropEntriesDirect('Rating',['X'])   #!# Add an option to reclassify as FP? #!#
            if prex != rdb.entryNum(): self.printLog('#AMBIG','Dropped %s potentially ambiguous ("X"-Rated) instances: occur in TP protein at non-TP location.' % rje.iStr(prex - rdb.entryNum()))
            rdb.dropField('RunID')
            rdb.saveToFile(rfile)
            return True
        except: self.errorLog('Problem during OccBench ratings.'); return False  # rating failed
#########################################################################################################################
    def occAssessment(self):   ### Perform benchmarking assessment for OccBench results                             #V2.4
        '''Perform benchmarking assessment for OccBench results. Currently simplistic without occfilter.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ASSESSMENT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.db()
            rdb = self.db('results')
            afile = '%s.assessment.tdt' % self.getStr('BenchBase')
            rje.backup(self,afile,appendable=False)
            akeys = self.list['RunID'][0:]
            bfile = '%s.bymotif.tdt' % self.getStr('BenchBase')
            rje.backup(self,bfile,appendable=False)
            bkeys = ['Motif'] + akeys
            rfile = '%s.ratings.tdt' % self.getStr('BenchBase')
            rkeys = self.db('results').keys()

            ### ~ [1] Convert Ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for field in ['TP','FP','FN','TN']: rdb.addField(field,evalue=0)
            rdb.addField('ResNum',evalue=1)
            for entry in rdb.entries(): entry[entry['Rating']] = 1
            rdb.dropField('Rating')

            ### ~ [2] Compress and rate by motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb.compress(bkeys+['Seq'],rules={},default='sum',best=[])  #!# Make default='mean' to normalise by sequence
            rdb.compress(bkeys,rules={},default='sum',best=[])
            rdb.dropFields(['Seq','Start_Pos','End_Pos'])
            for field in ['SN','FPR','PPV','wPPV']: rdb.addField(field)
            for entry in rdb.entries():
                try:
                    entry['SN'] = float(entry['TP']) / (float(entry['TP']) + float(entry['FN']))
                    entry['FPR'] = float(entry['FP']) / (float(entry['FP']) + float(entry['TN']))
                    entry['PPV'] = float(entry['TP']) / (float(entry['TP']) + float(entry['FP']))
                    entry['wPPV'] = float(entry['SN']) / (float(entry['SN']) + float(entry['FPR']))
                except: pass
                self.debug(entry)
            rdb.saveToFile(bfile)

            ### ~ [3] Compress across motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdb.addField('N',evalue=1)
            rdb.compress(akeys,rules={'SN':'mean','FPR':'mean','PPV':'mean','wPPV':'mean'},default='sum',best=[])
            rdb.deleteField('Motif')
            rdb.saveToFile(afile)
            return True
        except: self.errorLog('Problem during OccBench assessment.'); return False  # assessment failed
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
    except KeyboardInterrupt: mainlog.errorLog('User terminated.',quitchoice=False)
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
