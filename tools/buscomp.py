#!/usr/bin/python

# See below for name and description
# Copyright (C) 2019 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       BUSCOMP
Description:  BUSCO Compiler and Comparison tool
Version:      0.9.7
Last Edit:    03/12/20
Citation:     Edwards RJ (2019). F1000Research 8:995 (slides) (doi: 10.7490/f1000research.1116972.1)
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    BUSCOMP is designed to overcome some of the non-deterministic limitations of BUSCO to:

    1. compile a non-redundant maximal set of complete BUSCOs from a set of assemblies, and
    2. use this set to provide a "true" comparison of completeness between different assemblies of the same genome
    with predictable behaviour.

    For each BUSCO gene, BUSCOMP will extract the best "Single Complete" sequence from those available, using the
    `full_table_*.tsv` results table and `single_copy_busco_sequences/` directory of hit sequences. BUSCOMP ranks all
    the hits across all assemblies by Score and keeps the top-ranking hits. Ties are then resolved by Length, keeping
    the longest sequence. Ties for Score and Length will keep an arbitrary entry as the winner. Single complete hits
    are given preference over Duplicate hits, even if they have a lower score, because only Single hit have their
    sequences saved by BUSCO in the `single_copy_busco_sequences/` directory. This set of predicted gene sequences
    forms the "BUSCOMPSeq" gene set.

    BUSCOMP uses minimap2 to map BUSCOSeq predicted CDS sequences onto genome/transcriptome assemblies, including
    those not included in the original BUSCO compilation. This way, the compiled set of species-specific BUSCO
    sequences can also be used to generate a quick-and-dirty assessment of completeness for a new genome assembly.
    Hits are converted into percentage coverage stats, which are then used to reclassify the BUSCO gene on the basis
    of coverage and identity. BUSCOMP ratings are designed to mimic the original BUSCO ratings but have different
    definitions. In addition, two extra classes of low quality hit have been added: "Partial" and "Ghost".

    * **Complete**: 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
    * **Duplicated**: 95%+ Coverage in 2+ contigs/scaffolds.
    * **Fragmented**: 95%+ combined coverage but not in any single contig/scaffold.
    * **Partial**: 40-95% combined coverage.
    * **Ghost**: Hits meeting local cutoff but <40% combined coverage.
    * **Missing**: No hits meeting local cutoff.

    In addition to individual assembly stats, BUSCO and BUSCOMP ratings are compiled across user-defined groups of
    assemblies with various outputs to give insight into how different assemblies complement each other. Ratings are
    also combined with traditional genome assembly statistics (NG50 and LG50) based on a given `genomesize=X` to help
    identify the "best" assemblies. Details of settings, key results, tables and plots are output to an HTML report
    using Rmarkdown.

    NOTE: For HTML output, R must be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    NOTE: BUSCOMPSeq sequences can be provided with `buscofas=FILE` in place of compilation. This option has not been
    tested and might give some unexpected behaviours, as some of the quoted figures will still be based on the
    calculated BUSCOMPSeq data. Please report any unexpected behaviour.

    For full documentation of the BUSCOMP workflow, run with `dochtml=T` and read the `*.docs.html` file generated, or
    visit https://slimsuite.github.io/buscomp/.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    runs=DIRLIST    : List of BUSCO run directories (wildcards allowed) [run_*]
    fastadir=DIRLIST: List of directories containing genome fasta files (wildcards allowed) [./]
    fastaext=LIST   : List of accepted fasta file extensions that will be checked for in fastadir [fasta,fas,fsa,fna,fa]
    genomes=FILE    : File of Prefix and Genome fields for generate user-friendly output [*.genomes.tdt if found]
    restrict=T/F    : Restrict analysis to genomes with a loaded alias [False]
    runsort=X       : Output sorting order for genomes and groups (or "Genome","Prefix","Complete","Group") [Group]
    stripnum=T/F    : Whether to strip numbers ("XX_*") at start of Genome alias in output [True]
    groups=FILE     : File of Genome and Group fields to define Groups for compilation [*.groups.tdt]
    buscofas=FASFILE: Fasta file of BUSCO DNA sequences. Will combine and make (NR) if not given [None]
    buscomp=T/F     : Whether to run BUSCO compilation across full results tables [True]
    dupbest=T/F     : Whether to rate "Duplicated" above "Complete" when compiling "best" BUSCOs across Groups [False]
    buscompseq=T/F  : Whether to run full BUSCO comparison using buscofas and minimap2 [True]
    ratefas=FILELIST: Additional fasta files of assemblies to rate with BUSCOMPSeq (No BUSCO run) (wildcards allowed) []
    rmdreport=T/F   : Generate Rmarkdown report and knit into HTML [True]
    ggplot=T/F      : Whether to use ggplot code for plotting [True]
    fullreport=T/F  : Generate full Rmarkdown report including detailed tables of all ratings [True]
    missing=T/F     : Generate summary tables for sets of Missing genes for each assembly/group [True]
    dochtml=T/F     : Generate HTML BUSCOMP documentation (*.docs.html) instead of main run [False]
    summarise=T/F   : Include summaries of genomes in main `*.genomes.tdt` output [True]
    loadsummary=T/F : Use existing genome summaries including NG50 from `*.genomes.tdt`, if present [True]
    ### ~ Mapping/Classification options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    endextend=X     : Extend minimap2 hits to end of sequence if query region with X bp of end [0]
    minlocid=INT    : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [20]
    uniquehit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [True]
    mmsecnum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [3]
    mmpcut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]
    mapopt=CDICT    : Dictionary of additional minimap2 options to apply (caution: over-rides conflicting settings) []
    alnseq=T/F      : Whether to use alnseq-based processing (True) or (False) faster CS-Gstring processing [False]
    ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import hashlib, os, string, sys, time, glob
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
gitcodepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)))) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'../code/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_menu, rje_paf, rje_rmd, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Basic working version.
    # 0.2.0 - Functional version with basic RMarkdown HTML output.
    # 0.3.0 - Added ratefas=FILELIST: Additional fasta files of assemblies to rate with BUSCOMPSeq (No BUSCO run) [].
    # 0.4.0 - Implemented forking and tidied up output a little.
    # 0.5.0 - Updated genome stats and RMarkdown HTML output. Reorganised assembly loading and proeccessing. Added menus.
    # 0.5.1 - Reorganised code for clearer flow and documentation. Unique and missing BUSCO output added.
    # 0.5.2 - Dropped paircomp method and added Rmarkdown control methods. Updated Rmarkdown descriptions. Updated log output.
    # 0.5.3 - Tweaked log output and fixed a few minor bugs.
    # 0.5.4 - Deleted some excess code and tweaked BUSCO percentage plot outputs.
    # 0.5.5 - Fixed minlocid bug and cleared up minimap temp directories. Added LnnIDxx to BUSCOMP outputs.
    # 0.5.6 - Added uniquehit=T/F : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [False]
    # 0.6.0 - Added more minimap options, changed defaults and dev generation of a table changes in ratings from BUSCO to BUSCOMP.
    # 0.6.1 - Fixed bug that was including Duplicated sequences in the buscomp.fasta file. Added option to exclude from BUSCOMPSeq compilation.
    # 0.6.2 - Fixed bug introduced that had broken manual group review/editing.
    # 0.7.0 - Updated the defaults in the light of test analyses. Tweaked Rmd report.
    # 0.7.1 - Fixed unique group count bug when some genomes are not in a group. Fixed running with non-standard options.
    # 0.7.2 - Added loadsummary=T/F option to regenerate summaries and fixed bugs running without BUSCO results.
    # 0.7.3 - Fixed bugs calculating Complete BUSCO scores in a couple of places. Added text summaries to plots.
    # 0.7.4 - Added ggplot option. Added group plots to full reports.
    # 0.7.5 - Reinstated BUSCOMP contribution reports when re-running.
    # 0.7.6 - Added additional error-handling for CS parsing errors.
    # 0.7.7 - Fixed problems with buscompseq=F. Fixed stripnum and Rmd bugs. Added sequence name checking for duplicates.
    # 0.7.8 - Fixed a bug where BUSCOMP was not being compiled for assemblies without BUSCO data.
    # 0.7.9 - Added listing of numbers to BUSCOMP Missing charts.
    # 0.8.0 - Added alnseq=F as default PAF parsing mode for improved efficiency.
    # 0.8.1 - Set endextend=0 due to bug.
    # 0.8.2 - Fixed full RMD chart labelling bug. Fixed endextend bug and reinstated endextend=10 default.
    # 0.8.3 - Fixed Unique rating bug with no groups.
    # 0.8.4 - Set endextend=0 due to another bug.
    # 0.8.5 - Fixed BUSCO table loading bug introduced by Diploidocus. Added error catching for logbinomial bug.
    # 0.8.6 - Tweaked code to handle BUSCO v4 files, but not (yet) file organisation.
    # 0.8.7 - Fixing issues with prefix parsing from BUSCO directories and files.
    # 0.9.0 - Updated parsing of single_copy_busco_sequences/ to enable multiple directories with "$PREFIX" suffixes.
    # 0.9.1 - Updated parsing to enable BUSCO v4 results recognition. (run with -o $GENOME.busco)
    # 0.9.2 - Fixed some bugs when files missing.
    # 0.9.3 - Minor fixes to output and clearer error messages. Fixed formatting for Python 2.6 back compatibility for servers.
    # 0.9.4 - Added contig statistics and fixed group description loading bug.
    # 0.9.5 - Fixed Group BUSCOMP plot output bug.
    # 0.9.6 - Added CtgNum: Number of contigs (`SeqNum`+`GapCount`).
    # 0.9.7 - Fixed some Rmd bugs to fix output after summary table changes.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [?] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [?] : Add `groupfas=T` to output a `*.buscomp.$GROUP.fasta` file per group.
    # [Y] : Add an option to strip a numerical prefix added to the run directory for sorting: run_XX_$PREFIX
    # [Y] : Add minlocid to this and PAF to speed up the CDS mapping. Q: Should this be the same as the overal LocalID?
    # [Y] : Rscript -e 'library(rmarkdown); rmarkdown::render("/path/to/test.Rmd", "html_document")'
    # [Y] : Add some additional log markers of program stages
    # [Y] : Update docs to include details of min cutoffs
    # [ ] : Add variables to control ratings cutoffs.
    # [Y] : Add additional text to RMarkdown results file, explaining results and outputs.
    # [Y] : Add capacity to skip BUSCO compilation and still perform BUSCOMPSeq analysis.
    # [Y] : Update logging with increased use of bugLog in place of printLog.
    # [Y] : Add group compilation of BUSCOMP results as well as BUSCO results.
    # [Y] : Add "Unique" field to BUSCO(MP) output = Complete genes not found in any other assembly
    # [ ] : Modify use of Identical to be a special subcategory of Single Complete, which can then be compiled.
    # [ ] : Test summarise=F.
    # [Y] : Add restrict=T/F option to restrict output to genomes with loaded prefixes (for re-running on reduced data).
    # [Y] : Finish documentation in run() docstring.
    # [ ] : Make the left margin for BUSCO plots more responsive to Genome labels.
    # [ ] : Consider dropping the localnid setting in favour of a global minID threshold?
    # [Y] : minlocid threshold update inc. docs
    # [Y] : uniquehit threshold update inc. docs
    # [Y] : mmsecnum threshold update inc. docs
    # [Y] : Fix *.buscomp.tdt output bug if BUSCO sequences not given and buscofas used: runs=../example/fulltables/ buscofas=test1.buscomp.fasta
    # [Y] : Add option to include text ratings summaries in the bar plots.
    # [ ] : Add group plot output.
    # [ ] : Make the plot unit scaling more responsive to limits, e.g. Mb and straight counts.
    # [ ] : Add option for BLAST+ replacement of Minimap2
    # [Y] : Add checking of input sequence files: fasta format and unique names.
    # [Y] : Check/fix bug with buscompseq=F.
    # [Y] : Update GABLAM statistics to work directly from CS strings.
    # [ ] : Don't plot BUSCO/BUSCOMP stats if relevant analysis not performed.
    # [Y] : Fix issue of group descriptions not loading.
    # [ ] : Add contig stats to output.
    # [ ] : Add warning of recalculation if complete genome.tdt file found but alias file given. (Or fix re-use issue.)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('BUSCOMP', '0.9.7', 'December 2020', '2019')
    description = 'BUSCO Compiler and Comparison tool'
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
            if rje.yesNo('Show RJE_PAF (Minimap2) commandline options?',default='N'): out.verbose(-1,4,text=rje_paf.__doc__)
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
        if len(sys.argv) == 2 and sys.argv[1] in ['version','-version','--version']: rje.printf(info.version); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('{0} v{1}'.format(info.program,info.version)); sys.exit(0)
        if len(sys.argv) == 2 and sys.argv[1] in ['description','-description','--description']: rje.printf('%s: %s' % (info.program,info.description)); sys.exit(0)
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
dformats = {'Complete':'int','Duplicated':'int','Fragmented':'int','Missing':'int','N':'int',
            'Sequences':'bool',
            'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int',
            'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int',
            'Length':'int','QryLen':'int','HitLen':'int','Qry_AlnLen':'num','Qry_AlnID':'num',
            'Identical':'int','Single':'int','Partial':'int','Ghost':'int'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: BUSCOMP Class                                                                                           #
#########################################################################################################################
class BUSCOMP(rje_obj.RJE_Object):
    '''
    BUSCOMP Class. Author: Rich Edwards (2019).

    Str:str
    - BUSCOFas=FASFILE: Fasta file of BUSCO DNA sequences. Will combine and make (NR) if not given [None]
    - Genomes=FILE    : File of Prefix and Genome fields for generate user-friendly output [*.genomes.tdt if found]
    - Groups=FILE     : File of Genome and Group fields to define Groups for compilation [*.groups.tdt]

    Bool:boolean
    - BUSCOMP=T/F     : Whether to run BUSCO compilation across full results tables [True]
    - BUSCOMPSeq=T/F  : Whether to run full BUSCO comparison using buscofas and minimap2 [True]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - DupBest=T/F     : Whether to rate "Duplicated" above "Complete" when compiling "best" BUSCOs across Groups [False]
    - FullReport=T/F  : Generate full Rmarkdown report including detailed tables of all ratings [True]
    - GGPlot=T/F      : Whether to use ggplot code for plotting [True]
    - LoadSummary=T/F : Use existing genome summaries including NG50 from `*.genomes.tdt`, if present [True]
    - Missing=T/F     : Generate summary tables for sets of Missing genes for each assembly/group [True]
    - Restrict=T/F    : Restrict analysis to genomes with a loaded alias [False]
    - RmdReport=T/F   : Generate Rmarkdown report and knit into HTML [True]
    - StripNum=T/F    : Whether to strip numbers ("XX_*") at start of Genome alias in output [True]
    - Summarise=T/F   : Include summaries of genomes in main `*.genomes.tdt` output [True]
    - UniqueHit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [True]

    Int:integer
    - MinLocID=INT    : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [20]
    - MMSecNum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [3]

    Num:float
    - MMPCut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]

    File:file handles with matching str filenames
    
    List:list
    - FastaDir=DIRLIST: List of directories containing genome fasta files (wildcards allowed) [./]
    - FastaExt=LIST   : List of accepted fasta file extensions that will be checked for in fastadir [fasta,fas,fsa,fna,fa]
    - RateFas=FILELIST: Additional fasta files of assemblies to rate with BUSCOMPSeq (No BUSCO run) []
    - Ratings = List of rating hierarchy [Duplicated, Complete, Fragmented, Missing]
    - Runs=DIRLIST    : List of BUSCO run directories (wildcards allowed) [run_*]
    - RunSort = Will be a list of Genomes in the order that they are to be output. Loaded, or special generation.

    Dict:dictionary
    - Alias = dictionary of {prefix/genome: {genome table entry} } for easy mapping back on to core data
    - Genomes = dictionary of prefix:SeqList object containing genome
    - Groups = dictionary of {group: [list of genome table entries]}

    Obj:RJE_Objects
    - BUSCOMP:SeqList object containing compiled BUSCO hit sequences
    - DB:Database object containing main data
    - Rmd:Rmd object for RMarkdown output
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BUSCOFas','Genomes','Groups']
        self.boollist = ['BUSCOMP','BUSCOMPSeq','DocHTML','DupBest','GGPlot','FullReport','LoadSummary','Missing','Restrict','RmdReport','Summarise','UniqueHit']
        self.intlist = ['MinLocLen','MinLocID','MMSecNum']
        self.numlist = ['MMPCut']
        self.filelist = []
        self.listlist = ['FastaDir','FastaExt','RateFas','Ratings','Runs','RunSort']
        self.dictlist = ['Alias','Genomes','Groups']
        self.objlist = ['BUSCOMP','DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'BUSCOMP':True,'BUSCOMPSeq':True,'DocHTML':False,'FullReport':True,'GGPlot':True,'LoadSummary':True,
                    'Missing':True,'Restrict':False,'RmdReport':True,'StripNum':True,'Summarise':True,'UniqueHit':True})
        self.setInt({'MinLocLen':20,'MinLocID':0,'MMSecNum':3})
        self.setNum({'MMPCut':0.0})
        self.list['FastaDir'] = ['./']
        self.list['FastaExt'] = string.split('fasta,fas,fsa,fna,fa',sep=',')
        self.list['Ratings'] =  ['Duplicated', 'Complete', 'Fragmented', 'Missing']
        self.list['Runs'] = glob.glob('run_*')
        self.list['RunSort'] = ['Group']
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
        self.obj['Rmd'] = rje_rmd.Rmd(self.log,self.cmd_list)
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
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['BUSCOFas','Genomes','Groups'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BUSCOMP','BUSCOMPSeq','DocHTML','DupBest','FullReport','GGPlot','LoadSummary','Missing','Restrict','RmdReport','Summarise','UniqueHit'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinLocLen','MinLocID','MMSecNum'])   # Integers
                self._cmdReadList(cmd,'float',['MMPCut']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['FastaExt','RunSort'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['FastaDir','RateFas','Runs']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ## Tidy commands ##
        if not self.getStrLC('Basefile'): self.baseFile('buscomp'); self.cmd_list.append('buscomp')
#########################################################################################################################
    ### <2> ### Main Run Method                                                                                         #
#########################################################################################################################
    def run(self):  ### Main run method                                                                         # v0.5.2
        '''
        # BUSCOMP: BUSCO Compiler and Comparison tool

        BUSCOMP is designed to overcome some of the non-deterministic limitations of BUSCO to:

        1. compile a non-redundant maximal set of complete BUSCOs from a set of assemblies, and
        2. use this set to provide a "true" comparison of completeness between different assemblies of the same genome
        with predictable behaviour.

        For each BUSCO gene, BUSCOMP will extract the best "Single Complete" sequence from those available, using the
        `full_table_*.tsv` results table and `single_copy_busco_sequences/` directory of hit sequences. BUSCOMP ranks all
        the hits across all assemblies by Score and keeps the top-ranking hits. Ties are then resolved by Length, keeping
        the longest sequence. Ties for Score and Length will keep an arbitrary entry as the winner. Single complete hits
        are given preference over Duplicate hits, even if they have a lower score, because only Single hit have their
        sequences saved by BUSCO in the `single_copy_busco_sequences/` directory. This set of predicted gene sequences
        forms the "BUSCOMPSeq" gene set.

        BUSCOMP uses minimap2 to map BUSCOSeq predicted CDS sequences onto genome/transcriptome assemblies, including
        those not included in the original BUSCO compilation. This way, the compiled set of species-specific BUSCO
        sequences can also be used to generate a quick-and-dirty assessment of completeness for a new genome assembly.
        Hits are converted into percentage coverage stats, which are then used to reclassify the BUSCO gene on the basis
        of coverage and identity. BUSCOMP ratings are designed to mimic the original BUSCO ratings but have different
        definitions. In addition, two extra classes of low quality hit have been added: "Partial" and "Ghost".

        * **Complete**: 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
        * **Duplicated**: 95%+ Coverage in 2+ contigs/scaffolds.
        * **Fragmented**: 95%+ combined coverage but not in any single contig/scaffold.
        * **Partial**: 40-95% combined coverage.
        * **Ghost**: Hits meeting local cutoff but <40% combined coverage.
        * **Missing**: No hits meeting local cutoff.

        In addition to individual assembly stats, BUSCO and BUSCOMP ratings are compiled across user-defined groups of
        assemblies with various outputs to give insight into how different assemblies complement each other. Ratings are
        also combined with traditional genome assembly statistics (NG50 and LG50) based on a given `genomesize=X` to help
        identify the "best" assemblies. Details of settings, key results, tables and plots are output to an HTML report
        using Rmarkdown.

        **NOTE:** For HTML output, R must be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        **NOTE:** BUSCOMPSeq sequences can be provided with `buscofas=FILE` in place of compilation. This option has not been
        tested and might give some unexpected behaviours, as some of the quoted figures will still be based on the
        calculated BUSCOMPSeq data. Please report any unexpected behaviour.

        ---

        # Running BUSCOMP

        BUSCOMP is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/buscomp.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [BUSCOMP git repo](https://github.com/slimsuite/buscomp), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [BUSCOMP git repo](https://github.com/slimsuite/buscomp)
        for running on example data.

        For BUSCOMPSeq analysis, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
        environment `$PATH` or given to BUSCOMP with the `minimap2=PROG` setting.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        runs=DIRLIST    : List of BUSCO run directories (wildcards allowed) [run_*]
        fastadir=DIRLIST: List of directories containing genome fasta files (wildcards allowed) [./]
        fastaext=LIST   : List of accepted fasta file extensions that will be checked for in fastadir [fasta,fas,fsa,fna,fa]
        genomes=FILE    : File of Prefix and Genome fields for generate user-friendly output [*.genomes.tdt if found]
        restrict=T/F    : Restrict analysis to genomes with a loaded alias [False]
        runsort=X       : Output sorting order for genomes and groups (or "Genome","Prefix","Complete","Group") [Group]
        stripnum=T/F    : Whether to strip numbers ("XX_*") at start of Genome alias in output [True]
        groups=FILE     : File of Genome and Group fields to define Groups for compilation [*.groups.tdt]
        buscofas=FASFILE: Fasta file of BUSCO DNA sequences. Will combine and make (NR) if not given [None]
        buscomp=T/F     : Whether to run BUSCO compilation across full results tables [True]
        dupbest=T/F     : Whether to rate "Duplicated" above "Complete" when compiling "best" BUSCOs across Groups [False]
        buscompseq=T/F  : Whether to run full BUSCO comparison using buscofas and minimap2 [True]
        ratefas=FILELIST: Additional fasta files of assemblies to rate with BUSCOMPSeq (No BUSCO run) (wildcards allowed) []
        rmdreport=T/F   : Generate Rmarkdown report and knit into HTML [True]
        fullreport=T/F  : Generate full Rmarkdown report including detailed tables of all ratings [True]
        missing=T/F     : Generate summary tables for sets of Missing genes for each assembly/group [True]
        dochtml=T/F     : Generate HTML BUSCOMP documentation (*.docs.html) instead of main run [False]
        summarise=T/F   : Include summaries of genomes in main `*.genomes.tdt` output [True]
        ### ~ Mapping/Classification options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        minimap2=PROG   : Full path to run minimap2 [minimap2]
        endextend=X     : Extend minimap2 hits to end of sequence if query region with X bp of end [0]
        minlocid=INT    : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1]
        uniquehit=T/F   : Option to use *.hitunique.tdt table of unique coverage for GABLAM coverage stats [True]
        mmsecnum=INT    : Max. number of secondary alignments to keep (minimap2 -N) [3]
        mmpcut=NUM      : Minimap2 Minimal secondary-to-primary score ratio to output secondary mappings (minimap2 -p) [0]
        mapopt=CDICT    : Dictionary of additional minimap2 options to apply (caution: over-rides conflicting settings) []
        ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        ## Input Data: BUSCO runs and Assembly fasta files

        STEP 1 is to set up the input data. BUSCOMP will first identify and check the run directories (given with the
        `runs=DIRLIST` command). This command can either be a list of directories, or an expandable list with wildcard.
        A list of individual directories can be provided as a file if desired.

        **NOTE:** `runs=DIRLIST` wants the actual BUSCO run directories (`run_$PREFIX[.busco]/`), not the parent
        directory. These directories can optionally have a `.busco` suffix, which will be trimmed off. By default,
        BUSCOMP will look for `./run*`.

        From these, the run list will be extracted, consisting of a number of genome `$PREFIX` values. It is expected
        that the `$PREFIX` for each run will match the input genome file, found in a directory set by `fastadir=DIRLIST`.
        This file can also match the `Genome` alias for that genome (see below). If you have not used consistent naming
        across files, please check data has been loaded and mapped correctly during the Review phase (below).

        The presence of a `full_$PREFIX[.busco].tsv` table and `single_copy_busco_sequences/` directory within each run
        directory will also be checked. If `buscocompseq=F` then all of the `full_*.tsv` tables in the `runs`
        directories will be loaded for compilation, but no sequence searching will be performed. The presence of
        sequences available for compilation will be stored in the `Sequences` field of `*.genomes.tdt`.

        **BUSCO v4.x:** From `v0.9.1`, BUSCOMP should recognise BUSCO v4 output. Due to the reorganisation and altered
        naming strategy of version 4, the naming of results files is more strict. The main results directory (set by
        `-o` when you run BUSCO) should match the genome assembly prefix with an optional `.busco` suffix. (A leading
        `run_` is also permitted and will be ignored. For example, if `assembly.fasta` was analysed with BUSCO v4,
        BUSCOMP should be able to parse results generated using `-o assembly`, `-o assembly.busco`, `-o run_assembly`, or
        `-o run_assembly.busco`. If none of these settings were used, the results directory can be manually renamed. For
        additional sorting, a XX_ numerical prefix _can_ be used (see below), e.g. `run_01_assembly` will still look
        for `assembly.*` in the `fastadir` path.

        ### BUSCOMPSeq Analysis Only

        Additional assemblies can rated using the BUSCOMPSeq analysis (see below) without having BUSCO data analysed.
        These are given with the `ratefas=LIST` command (wildcards allowed), and will be added to the main Genomes table
        after BUSCO compilation has been performed. The default prefix for these files will be their basename (filename,
        minus path and extension). This can be linked to an alias and/or have a prefix number stripped, as with the
        run directories. When combined with `buscofas=FASFILE`, this enables the BUSCOMPSeq analysis to be performed in
        the absence of _any_ BUSCO results.

        ### Genome Prefixes

        Each genome analysed has a "Prefix" that can be used to help sorting (if `runsort=prefix`) and identify relevant
        BUSCO files where the genome (Fasta) name and BUSCO run name do not match. It is generally assumed that this
        Prefix will match the prefix of the original assembly file on which BUSCO was run. This Prefix is set once as the
        data is loaded and is never changed. (Sorting and visualisations can be made more informative altered using the
        Genome (a.k.a. Alias) and Description fields.)

        The assembly Prefix will be set as:

        1. If a `full_table_*.tsv` file is found, the core (indicated with `*`) will be used as the prefix.
        2. If no BUSCO results are available, the path- and extension-stripped fasta file name will be used for the prefix

        ### Genome Aliases

        Genome Aliases will be parsed initially from the loaded data:

        1. If a BUSCO run directory is given and it contains a single BUSCO run, the run directory name will be used for
        the `Genome` alias, minus the leading "`run_`".
        2. If the run directory contains multiple BUSCO runs, the `full_table_*.tsv` core used for the `Prefix` will also
        be used for the `Genome` alias.
        3. If no BUSCO results are available, the path- and extension-stripped fasta file name will be used for the
        `Genome` alias.

        If `Genome` contains a trailing `.busco`, this will be stripped. If `Genome` starts with preceding numbers
        (`XX_`), these will be kept for sorting but stripped for outputs and fasta matching (unless `stripnum=F`).

        **NOTE:** If the assembly fasta file does not match the name used for the BUSCO run file, then the `Genome` alias
        *must* be set to the `basefile` of the fasta file by appropriate naming of the `run_` directory, _i.e._ the
        `run_*` directory should contain a single BUSCO run and the directory name should match a `*.fasta` file (or
        other accepted `fastaext=LIST` extension) in one of the directories given by `fastadir=DIRLIST`.

        An optional alias table can be provided (`genomes=FILE`) that contains `Prefix` and `Genome` fields (others
        allowed) and can give alternative names for each genome to be used in other outputs. If running interactively,
        there will also be an option to check/update these aliases manually. Genome aliases are output as part of the
        main `*.genomes.tdt` output, which can also be used as subsequent `genomes=FILE` input (the default, if found).

        This table can also have a `Description` field, which will be used to update descriptions for output if found.
        (Empty descriptions will not be mapped.) Otherwise, descriptions can be parsed from a `description_$PREFIX.txt`
        file (exactly matching the `$PREFIX` parsed from a `full_table_*.tsv` file) if found. Alternatively, if the run
        directory has a `run_` prefix AND there is only one `full_table`, or it is a BUSCO v4 output directory, a
        description will be parsed from `description.txt`. In each case, only the first line will be used. If neither
        condition is met, the name of the genome fasta file will be used. Failing that, `Genome` will be used. (Note that
        for BUSCO v4, `description.txt` must be in the same directory as `full_table.tsv`.

        **NOTE:** The optional Alias table *can* be used to change a `Genome` name to something other than its Fasta file
        - this mapping is performed _after_ fasta files have been identified.

        **NOTE:** `Prefix` and `Genome` fields must be unique. A common `Prefix` and `Genome` value is permitted but
        these will be assumed to refer to the _same_ assembly. (In other words, do not run BUSCO on two assemblies
        and give each BUSCO run the name of the other assembly!)

        **NOTE:** Some `Genome` names may cause conflicts with running minimap2 and/or the Rmarkdown output.
        **Whitespace is not permitted** in `Genome` names (though can be in R labels) and will be stripped. It is
        recommended to keep `Genome` names simple and only containing standard characters (letters, numbers, dots and
        underscores).

        ---

        ## Grouping and BUSCO compilation

        BUSCOMP compilation will take a set of genomes and identify the best rating for each BUSCO across the set. By
        default, this is done for all input genomes and assigned to a "BUSCOMP" rating. Ratings are compiled according
        to the `ratings=LIST` hierarchy, which by default is: `Complete`, `Duplicated`, `Fragmented`, `Missing`. The
        "best" rating (first in the hierarchy) will be used for that group, even if it is found in only a single
        assembly in the group. If `dupbest=T`, the hierarchy is `Duplicated`, `Complete`, `Fragmented`, `Missing`.

        ### Loading Groups

        Additional subsets of genomes can also be grouped together for compilation, using the `groups=FILE` input. This
        must have `Genome` and `Group` fields. Once generated, a Group becomes a "Genome". Groups can therefore include
        other Groups in them, provided there is no circular membership. Genome mapping should also recognise Prefix
        values. The "BUSCOMP" group with *all* genomes will be added unless a conflicting "BUSCOMP" group is provided,
        or `buscomp=F`.

        **NOTE:** Genome groups are generated mapping onto the `Genome` (or `Prefix`) values _after_ any alias mapping
        but _before_ any manual editing of `Genome` aliases.

        ### Manual group review/editing

        When running interactively (`i>=0`), Group review and editing can be accessed from the main Genome menu. This
        will show the genomes currently in each group, and provides the following menu options:

        * `<A>dd` = Create a new Group. First, a group name will need to be chosen, then the `<E>dit` menu will be
        triggered.
        * `<D>elete` = Remove the Group. Note that this will not remove Group members from any "parent" groups that
        contained the deleted group when they were loaded.
        * `<R>ename` = Change the Group name.
        * `<E>dit` = Edit group membership. This will give options to add/remove each assembly in the analysis.
        * `E<x>it menu` = Return to main Genome Menu.
        * `<Q>uit` = Quit BUSCOMP.

        ### Saving and sorting

        Once group edits are complete, group data will be saved to the `*.groups.tdt` file. If `runsort=group` then
        genome sorting will also be performed following group edits.


        ## BUSCOMP Genome/Group menu

        Following initial setup of genomes and groups, a Genome/Group edit menu enables updating loaded information.
        (This menu will not be available if running in non-interactive mode, e.g. `i=-1`.) The primary use for this is
        to convert unwieldy genome and/or BUSCO run names into something more user-friendly for display. It can also be
        used for adding descriptions to assemblies, adding or modifying group membership, and sorting genomes using
        different criteria. Assemblies can also be dropped from analysis from this menu.

        The main menu offers five options:

        - `R` = Review genomes
        - `G` = Review/edit genome groups
        - `S` = Sort genomes
        - `Q` = Quit BUSCOMP
        - `X` = Exit Menu and Proceed [default]

        ### Genome review menu

        If running in interactive mode, users will have the option to review and manually update loaded data. Each
        loaded genome will be displayed in turn, in its current sorted order. (Sorting can be altered or updated via
        the `Sort genomes` option of the main genome/group menu.) Options will be presented to modify data, or quit
        the menu (or BUSCOMP program). Modifications permitted are:

        - `<A>lias` = Sets the `Genome` field for the entry. This can have a "silent" numeric prefix for sorting if
        `runsort=Genome` or `runsort=Group`. Most importantly, this will be the default label in tables and plots, and
        used for fields in compiled data.
        - `<D>escription` = This gives the option to set a manual description of the assembly, which is used in the
        final output.
        - `<R>emove` = Drop this assembly from the analysis. BUSCOMP will ask for confirmation. This cannot be undone.

        If genomes are edited and `runsort=Genome` or `runsort=Group`, they will be re-sorted upon exiting the menu.


        ### Output order

        An optional `runsort=LIST` parameter will determine the output order. Future releases of BUSCOMP will proceess
        a mixed list of Prefix, Genome and Group names. Any data missing from the list will be appended at the end in
        Genome order. (NOTE: Group names are a special kind of Genome in this context.) This sorting is preformed after
        all manual reviews and updates are complete. Currently, providing a list of genomes (or `runsort=Manual`) will
        trigger manual ordering via the sort menu.

        There are also four special sort modes that can be implemented:

        - Genome = sort alphabetically based on Genome alias and Group name. (NOTE: "silent" numerical prefixes will
        affect this sorting (see below).
        - Group = sort alphabetically using `Genome` sorting within each group with the Group summary immediately
        following the last Group member. [Default]
        - Prefix = sort alphabetically based on Prefix. (NOTE: "silent" numerical prefixes will
        affect this sorting (see below).
        - Complete = sort in order of most complete to least complete.
        - Missing = sort in order of fewest missing to most missing BUSCOs

        The `runsort=prefix` option is best combined with the addition of a number to the run directory name:

        * `run_01_OneGenome/` containing `full_OneGenome.tsv` and generated from `OneGenome.fasta`
        * `run_02_AnotherGenome.busco/` containing `full_AnotherGenome.busco.tsv` and generated from `AnotherGenome.fasta`

        This example would sort OneGenome before AnotherGenome if runsort=prefix, but AnotherGenome before OneGenome if
        `runsort=Genome`.

        ---

        ## Genome Summary Statistics

        The final setup step is to load any genomes for which Fasta files have been provided. Unless `summarise=F`,
        summary statistics will be calculated for each genome, enabling assessements of genome completeness and quality
        in combination with the BUSCO and BUSCOMP results. Summary statistics are calculated by `rje_seqlist`. The
        following statistics are calculated for each genome:

        * **SeqNum**: The total number of scaffolds/contigs in the assembly.
        * **TotLength**: The total combined length of scaffolds/contigs in the assembly.
        * **MinLength**: The length of the shortest scaffold/contig in the assembly.
        * **MaxLength**: The length of the longest scaffold/contig in the assembly.
        * **MeanLength**: The mean length of scaffolds/contigs in the assembly.
        * **MedLength**: The median length of scaffolds/contigs in the assembly.
        * **N50Length**: At least half of the assembly is contained on scaffolds/contigs of this length or greater.
        * **L50Count**: The smallest number scaffolds/contigs needed to cover half the the assembly.
        * **CtgNum**: Number of contigs (`SeqNum`+`GapCount`).
        * **N50Ctg**: At least half of the assembly is contained on contigs of this length or greater.
        * **L50Ctg**: The smallest number contigs needed to cover half the the assembly.
        * **NG50Length**: At least half of the genome is contained on scaffolds/contigs of this length or greater.
        This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
        * **LG50Count**: The smallest number scaffolds/contigs needed to cover half the the genome.
        This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
        * **GapLength**: The total number of undefined "gap" (`N`) nucleotides in the assembly.
        * **GapCount**: The total number of undefined "gap" (`N`) regions above mingap=X (default 10bp).
        * **GC**: The %GC content of the assembly.

        **NOTE:** `NG50Length` and `LG50Count` statistics use `genomesize=X` or the biggest assembly loaded.
        If BUSCOMP has been run more than once on the same data (_e.g._ to update descriptions or sorting),
        previously calculates statistics will be read from the `genomes=FILE` table, if loaded. This may cause some
        inconsistencies between reported NG50 and LG50 values and the given/calculated `genomesize=X`. If partial
        results are reloaded following a change in `genomesize=X`, this will give rise to inconsistencies between
        assemblies. When re-running, please make sure that a consistent genome size is used.
        If in doubt, run with `force=T` and force regeneration of statistics.

        ---

        ## BUSCO Compilation

        Once all the genomes are loaded, BUSCOMP compiles the BUSCO results. This is done for three purposes:

        1. Report comparative BUSCO statistics across assemblies and plot with other assembly statistics.
        2. Combine BUSCO results across groups of assemblies to summarise total BUSCO coverage.
        3. Identify and compile the best `Complete` BUSCO hits for BUSCOMPSeq analysis (below).

        Loaded assemblies with an identified full results `Table` will be processed and added to the `*.full.tdt` table.
        The best rating per BUSCO gene will then be used for summary BUSCO counts for each genome. Finally, each Group
        will have its combined BUSCO rating calculated. This is performed by choosing the "best" rating for each BUSCO
        across the Group's members, where "best" is determined by the `dupbest=T/F` setting:

        - By default (`dupbest=F`), the rating hierarchy is: 'Complete', 'Duplicated', 'Fragmented', 'Missing'.
        - If `dupbest=T`, the rating hierarchy is: 'Duplicated', 'Complete', 'Fragmented', 'Missing'.

        Where ratings are defined (quoting from the [BUSCO v3 User Guide](http://gitlab.com/ezlab/busco/raw/master/BUSCO_v3_userguide.pdf) as:

        * `Complete`: Single-copy hits where "BUSCO matches have scored within the expected range of scores and within the expected range of length alignments to the
        BUSCO profile."
        * `Duplicated`: As `Complete` but 2+ copies.
        * `Fragmented`: "BUSCO matches ... within the range of scores but not within the range of length alignments to the BUSCO profile."
        * `Missing`: "Either no significant matches at all, or the BUSCO matches scored below the range of scores for the BUSCO profile."

        Total BUSCO counts (`N`) and summed Ratings are added to the main `*.genomes.tdt` table for each assembly and
        group. A `*.busco.tdt` file is also generated that has the rating for each BUSCO gene (rows) in each
        assembly/group (columns).

        If output is being sorted using `runsort=Complete` or `runsort=Missing`, data will be sorted at this point.
        Raw BUSCO results and compiled numbers will then be output to the log file.

        The following tables will then be saved (see **Output files**, below, for details):

        * `*.genomes.tdt` = summary of loaded data, genome statistics and BUSCO ratings.
        * `*.full.tdt` = full BUSCO results across all assemblies.
        * `*.busco.tdt` = compiled BUSCO results showing the rating per genome/group for each BUSCO gene.

        ### BUSCOMPSeq sequence compilation

        The final step of the BUSCOMP BUSCO compilation is to extract the best Complete BUSCO sequences from those
        available. For all assemblies with BUSCO results and a `single_copy_busco_sequences/`, BUSCOMP ranks all the hits
        across all assemblies by `Score` and keeps the top-ranking hits. Ties are then resolved by `Length`, keeping the
        longest sequence. Ties for `Score` and `Length` will keep an arbitrary entry as the winner. `Single` complete
        hits are given preference over `Duplicate` hits, even if they have a lower score, because only `Single` hit
        have their sequences saved by BUSCO in the `single_copy_busco_sequences/` directory.

        If `buscofas=FASFILE` has been given, this will be used for BUSCOMPSeq searches. Otherwise, the best BUSCO
        seqences identified for each BUSCO gene will be saved as `*.buscomp.fasta` for BUSCOMPSeq analysis. The
        exception is if `buscofas=FASFILE` is pointing to the `*.buscomp.fasta` and `force=T`, or if the
        `buscofas=FILE` fasta file cannot be found.

        **NOTE:** The `buscofas=FASFILE` option has not been tested and might give some unexpected behaviours, as some
        of the quoted figures will still be based on the calculated BUSCOMPSeq data.

        ## BUSCOMPSeq Minimap2 searches

        Once the gene sequences for BUSCOMPSeq have been established, the next step is to search for them in the assembly
        fasta files using a more deterministic (and faster) approach. For this, minimap2 has been chosen. BUSCOMP will
        perform a minimap2 search on each assembly fasta file. (Use `minimap2=FILE` to tell BUSCOMP where to find
        minimap2 if the path is not an environment variable.)

        Minimap2 will run with `-x splice -uf -C5` and hits will be filtered by length, keeping those that meet the
        `minloclen=X` cutoff (default 20 bp). A percentage identity cutoff can also be applied with `minlocid=X`
        (default 0%). Hits are also reduced to unique coverage by starting with the local alignment with the largest
        number of identical positions,
        and then trimming any other local alignments that overlap with the same part(s) of the query (BUSCO) sequence.
        This is repeated until all local alignments are non-overlapping (and still meet the identity and length
        criteria). Local hits are then compiled into global alignment (GABLAM) statistics of coverage and identity.

        By default, Minimap2 searches will keep the **Top 3** secondary alignments for each sequence (the minimap2 `-N`
        commandline parameter). In limited tests, this appears to give a good trade-off between search speed and the
        identification of `Duplicated` BUSCO genes. This default can be adjusted with `mmsecnum=INT` to make runs faster
        or more sensitive, as required. This can be further modulated using `mmpcut=NUM`, which sets the Minimap2 minimal
        secondary-to-primary score ratio to output secondary mappings (the minimap2 `-p` commandline parameter).

        Minimap2 search results are saved in a `*.minimap/` directory, with the prefix `$BASEFILE.$CUTOFFS`, where
        `$BASEFILE` is the fasta file basename for assembly, and `$CUTOFFS` takes the form `NXLXXIDXX`, where `NX` is the
        max number of secondary hits (`mmpcut=NUM`), `LXX` is the
        length cutoff (`minloclen=X`) and `IDXX` is the identity cutoff (`minlocid=X`) - `L20ID0` by default. Unless
        `uniquehit=F` (see below), `$CUTOFFS` will have a `U` suffix, e.g. `L20ID0U`. The Minimap2 `*.paf` file will only
        have the `.NX` suffix. A `*.NX.paf.cmd` file with the actual Minimap2 command used will also be saved, so the
        `mmpcut=NUM` setting and any other Minimap2 options set using `mapopt=CDICT` can be checked if required.

        Existing files will be loaded and reused unless `force=T`.

        Minimap2 hits are first converted into BLAST-style "local" hits, which are then converted into
        [GABLAM](http://rest.slimsuite.unsw.edu.au/gablam)-style summary statistics. For the `*.hitsum.tdt` table,
        local hits are reduced to unique coverage of each Query, such that each part of the Query is covered by a single
        hit. Local hits are selected in order of the number of identical positions between Query and Hit, followed by
        Length in the case of a tie. Overlapping regions from other local hits are trimmed, and the process repeated with
        the next-best local hit, ranked on trimmed stats where necessary. For the `*.gablam.tdt` output, each Query-Hit
        pair is considered independently and local hits are reduced to (a) unique coverage of the Query for `Qry_*`
        output, and (b) unique coverage of the Hit for `Hit_*` output. (Note that for each Query-Hit pair, only the local
        hits between that pair are considered for "unique" local hit reduction - there may be overlapping alignments to
        different queries/hits.) Unless `uniquehit=F`, hits will first be reduced to be non-overlapping across _all_
        queries before being coverted to Query-Hit pair GABLAM coverage stats.

        **NOTE:** Re-using results does NOT robustly check whether the BUSCOMPSeq data has changed, so this directory
        should be deleted if re-running with different sequences. BUSCOMP will save the name of the BUSCO file along
        with the md5sum of its contents to `*.NX.input.md5`, which will be checked if present and re-using data. Some basic
        checks should also be performed during the following results compilation stage, but please note that BUSCO IDs
        are used for sequence identifiers and so changes in BUSCO hit sequences will not be identified if files have
        been inappropriatley copied/moved between runs etc. - please look out for unexpected behaviour outside a "clean"
        run.

        **NOTE:** Minimap2 only works with high sequence identity. BUSCOMP is designed to be used on multiple assemblies
        _from the same species_. If using cross-species BUSCO analysis, odd results might be obtained, biased by the
        evolutionary distance from the species with the most BUSCOMP sequences. Under these circumstances, it might be
        better to swap out minimap2 for BLAST+. This can be achieved by independently running GABLAM using the
        BUSCOMPSeq fasta file as queries and each assembly as the search database, then replacing the `*.gablam.tdt` and
        `*.hitsum.tdt` files before re-running BUSCOMP. Please contact the author if you want help with this.

        ### BUSCOMPSeq Minimap2 search compilation.

        Minimap2 searches of the compiled BUSCOMP sequences are then compiled in similar vein to the original BUSCO
        results. The primary difference is that the search results need to first be converted into BUSCO-style ratings.
        This rating is explicitly more "quick and dirty" than the original BUSCO ratings, and should be considered
        complementary rather than superior. It aims to provide a quick, consistent assessment, but does have fewer checks
        and balances as a result.

        BUSCOMP uses results from both the `*.hitsum.tdt` table (overall coverage in assembly) and the `*.gablam.tdt`
        table (coverage in a single contig/scaffold) to rate sequences as:

        * **Complete**: 95%+ Coverage in a single contig/scaffold. (Note: unlike BUSCO, accuracy/identity is not considered.
        This will make it more prone to misidentification of closely related paralogues.)
        * **Duplicated**: 95%+ Coverage in 2+ contigs/scaffolds.
        * **Fragmented**: 95%+ combined coverage but not in any single contig/scaffold. (Note: as with BUSCOMP `Complete`
        ratings, these might include part of closely related paralogues.)
        * **Partial**: 40-95% combined coverage. (Note: these might be marked as "Fragmented" by BUSCO, which does not
        discriminate between "split" full-length hits, and single partial hits.)
        * **Ghost**: Hits meeting local cutoff but <40% combined coverage.
        * **Missing**: No hits meeting local cutoff.

        When compiling the results for all BUSCO genes, Single/Duplicated Complete hits will also be rated as
        **Identical** if they have 100% coverage and 100% identity in at least one contig/scaffold.

        Once all the individual assemblies have been rated for the full set of assemblies, results are compiled across
        Groups as described for the original BUSCO results (above). Because no genes receive an individual "Identical"
        rating, Groups will _not_ have a summary count for Identical hits.

        Individual gene ratings for each genome and group are output to `*.LnnIDxx.buscomp.tdt`, where `LnnIDxx` records
        the `minloclen=nn` and `minlocid=xx` settings. Compiled ratings are output to `*.LnnIDxx.ratings.tdt`.

        ### BUSCOMP Summary

        The final step of the BUSCOMP compilation is to summarise the findings in the log afile, akin to the BUSCO
        summary file. This will first generate a one line summary of the percentages, along with the original number of
        complete BUSCOs and the number of BUSCOMP sequences contributed by that assembly (i.e. the number with the best
        score of all the BUSCO searches.) This is followed by a more detailed breakdown of each category. For example:

        ```
        #BUSCO	00:00:41	Results: C:89.2%[S:87.8%,D:1.4%,I:23.3%],F:5.8%,P:3.8%,G:1.6%,M:0.9%,n:3736 - canetoad_v2.2 (3194 Complete BUSCOs; 102 BUSCOMP Seqs)
        #INFO	00:00:41	870 Identical BUSCOs (I)  [100% complete and identical]
        #INFO	00:00:41	3333 Complete BUSCOs (C)  [95%+ coverage in a single contig/scaffold]
        #INFO	00:00:41	3281 Complete and single-copy BUSCOs (S)  [1 Complete hit]
        #INFO	00:00:41	52 Complete and duplicated BUSCOs (D)  [2+ Complete hits]
        #INFO	00:00:41	217 Fragmented BUSCOs (F)  [95%+ coverage spread over 2+ contigs/scaffolds]
        #INFO	00:00:41	143 Partial BUSCOs (P)  [40-95% coverage]
        #INFO	00:00:41	61 Ghost BUSCOs (G)  [<40% coverage]
        #INFO	00:00:41	34 Missing BUSCOs (M)  [No hits]
        #INFO	00:00:41	3736 Total BUSCO gene hits searched
        ```

        ## BUSCO versus BUSCOMP comparisons

        There is a risk that performing a low stringency search will identify homologues or pseudogenes of the desired BUSCO gene in error.
        If there is a second copy of a gene in the genome that is detectable by the search then we would expect the same
        genes that go from `Missing` to `Complete` in some genomes to go from `Single` to `Duplicated` in others.

        To test this, data is reduced for each pair of genomes to BUSCO-BUSCOMP rating pairs of:

        * `Single`-`Single`
        * `Single`-`Duplicated`
        * `Missing`-`Missing`
        * `Missing`-`Single`

        This is then converted in to `G`ain ratings (`Single`-`Duplicated` & `Missing`-`Single`) or `N`o Gain ratings
        (`Single`-`Single` & `Missing`-`Missing`). The `Single`-`Duplicated` shift in one genome is then used to set the expected `Missing`-`Single`
        shift in the other, and assess the probability of observing the `Missing`-`Single` shift using a cumulative binomial
        distribution, where:

        * `k` is the number of observed `GG` pairs (`Single`-`Duplicated` _and_ `Missing`-`Single`)
        * `n` is the number of `Missing`-`Single` `G`ains in the focal genome (`NG`+`GG`)
        * `p` is the proportion of `Single`-`Duplicated` `G`ains in the background genome (`GN`+`GG` / (`GN`+`GG`+`NN`+`NG`))
        * `pB` is the probability of observing `k+` `Missing`-`Single` gains, given `p` and `n`

        This is output to `*.gain.tdt`, where each row is a Genome and each field gives the probability of the row
        genome's `Missing`-`Single` gains, given the column genome's `Single`-`Duplicated` gains.

        ---

        # Output files

        Contents of the main output files are given below. In addition to the tab-delimited text output, a summary HTML
        (and source RMarkdown) file will be generated, assuming R is installed. If `buscompseq=T` then a fasta file of
        the "best" (top-scoring) Single Complete BUSCO genes will also be generate (`*.buscomp.fasta`), unless it is
        provided using `buscofas=FASFILE`.

        ## BUSCOMPSeq Fasta files

        Unless `buscompseq=F` or `buscofas=FASFILE` is provided, nucleotide and protein sequences for the "best" BUSCO
        gene hits will be output to `*.buscomp.fasta` (nucleotide) and `*.buscomp.faa` (protein). The `*.buscomp.faa`
        protein file is not used by BUSCOMP but is provided as a useful set of (near-)complete protein sequences. Each
        fasta file is made from concatenating individual fasta files from the `single_copy_busco_sequences/` directories.
        As such, they will be named as with the standard BUSCO output:

        ```
        >$BUSCOID $FASTAFILE:$CONTIG:$START-$END
        ```

        **NOTE:** The `$FASTAFILE` path given in this file will be the one from the original BUSCO run, not the path to
        the genome file used by BUSCOMP if it has subsequently been moved/renamed.


        ## Data Tables

        ### genomes.tdt

        The `*.genomes.tdt` table is the main summary table for the input genomes, their BUSCO rating summaries and genome statistics.

        - `#` = Sorting order for output
        - `Directory` = Directory for this BUSCO run
        - `Prefix` = Prefix identified from BUSCO directory (or full table), corresponding to BUSCO output files
        - `Genome` = Optional alternative name (or "alias") to be used in outputs
        - `Fasta` = path to genome file, if found
        - Genome statistics = genome summary statistics if genome found and `summarise=T`. (See above.)
        - `Sequences` (True/False) = Whether `single_copy_busco_sequences/` directory found
        - `Table` = Path to `full_*.tsv` table
        - BUSCO Ratings summary will be added (assuming `Table` is `True`) (See above.)

        ### groups.tdt

        The `*.groups.tdt` table has the mappings between `Genome` and `Group` identifiers to make the groups for combined ratings (see
        above).

        - `#` = Arbitrary unique key for table
        - `Genome` = Genome alias or Prefix
        - `Group` = Name for grouped rating

        ### full.tdt

        The `*.full.tdt` table is a compiled version of the all the individual BUSCO full results tables.

        - `Genome` = Genome alias
        - `#` = Arbitrary unique key element for table (adjusting for duplicates)
\       - `BuscoID` = BUSCO gene identifier
        - `Status` = BUSCO rating
        - `Contig` = Contig containing BUSCO hit
        - `Start` = Start position
        - `End` = End position
        - `Score` = BUSCO score
        - `Length` = Length of BUSCO hit

        ### busco.tdt

        The `*.busco.tdt` table has the compiled ratings for individual BUSCO genes in each genome/group

        - `BuscoID` = BUSCO EOG gene identifier
        - Genomes and Groups fields have the rating of that gene in that Genome or Group

        ### buscoseq.tdt

        The `*.buscoseq.tdt` table has the set of best BUSCO genes used for the compiled BUSCOMPSeq sequence set. In
        addition to the BUSCO stats for all the BUSCOMP sequences, BUSCOMP ratings for each Genome are output in this
        file.

        - `Genome` = Genome alias
        - `BuscoID` = BUSCO gene identifier
        - `Status` = BUSCO rating
        - `Contig` = Contig containing BUSCO hit
        - `Start` = Start position
        - `End` = End position
        - `Score` = BUSCO score
        - `Length` = Length of BUSCO hit
        - Genomes fields have the BUSCOMP rating of that gene in that Genome

        ### buscomp.tdt

        The `*.LnnIDxx.buscomp.tdt` table is the same as the `*.busco.tdt` but with revised ratings based on BUSCOMP analysis.

        - `BuscoID` = BUSCO EOG gene identifier
        - Genomes and Groups fields have the rating of that gene in that Genome or Group

        ### ratings.tdt

        The `*.LnnIDxx.ratings.tdt` table has the compiled BUCOMP summary ratings per genome/group.

        - `Genome` = Genome alias or group name
        - `N` = Number of BUSCOMP genes
        - `Identical` = 100% coverage and 100% identity in at least one contig/scaffold. These will also be rated as
        `Complete` or `Duplicated`. Groups do not have `Identical` ratings.
        - `Complete` = 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
        - `Duplicated` = 95%+ Coverage in 2+ contigs/scaffolds.
        - `Fragmented` = 95%+ combined coverage but not in any single contig/scaffold.
        - `Partial` = 40-95% combined coverage.
        - `Ghost` = Hits meeting local cutoff but <40% combined coverage.
        - `Missing` = No hits meeting local cutoff.

        ### changes.tdt

        The `*.LnnIDxx.changes.tdt` table has counts of change ratings (BUSCO to BUSCOMP) for each genome, in addition
        to a Total count across all genomes.

        - `BUSCO` = BUSCO rating.
        - `BUSCOMP` = BUSCOMP rating.
        - Genome fields have the count of the number of genes with this combination of ratings.
        - `Total` = Summed count over all genomes.

        ### changes.full.tdt

        The `*.LnnIDxx.changes.full.tdt` table is the complete set of ratings changes (BUSCO to BUSCOMP) for each gene
        and genome.

        - `BuscoID` = BUSCO EOG gene identifier
        - Genome fields have the ratings change for that gene, using the first letters of the ratings.

        <small>**C**omplete, **D**uplicated, **F**ragmented, **P**artial, **G**host, **M**issing</small>

        ### unique.tdt

        The `*.LnnIDxx.unique.tdt` table has `Genome` and `Group` identifiers for each `BuscoID` where that gene is only
        `Complete` (or `Duplicated`) in that genome/group.

        - `BuscoID` = BUSCO EOG gene identifier
        - `BUSCO` = Genome/Group that uniquely has a BUSCO `Complete` rating.
        - `BUSCOMP` = Genome/Group that uniquely has a BUSCOMP `Complete` rating.

        ### rdata.tdt

        The `*.LnnIDxx.rdata.tdt` table has the BUSCOMP summary ratings (converted into percentage values) and sequence
        statistics for each `Genome`/`Group`, in addition to:

        - `BUSCO` = BUSCO `Complete` (`Single` and `Duplicated`) percentage
        - `NoBUSCO` = BUSCO `Missing` percentage
        - `UniqBUSCO` = Number of unique BUSCO `Complete` genes
        - `UniqBUSCOMP` = Number of unique BUSCOMP `Complete` genes
        - `col` = colour field for R plotting.
        - `pch` = point type for R plotting.
        - `label` = Genome/Group label for R plotting.
        - `plot` = TRUE/FALSE, whether to plot the Genome/Group statistics.
        - `best` = any categories for which this Genome is rated "best".


        ## Minimap2 output

        Minimap2 will be run for each genome with a `Fasta` file, generating:

        - `*.paf`
        - `*.paf.cmd`
        - `*.input.md5`
        - `*.L20ID80.local.tdt`
        - `*.L20ID80.hitunique.tdt`
        - `*.L20ID80.qryunique.tdt`
        - `*.L20ID80.hitsum.tdt`
        - `*.L20ID80.gablam.tdt`

        _Output details will be added here._

        ---

        # BUSCOMP Report (RMarkdown HTML output)

        The final step in BUSCOMP analysis is to generate a summary report. This is produced as an RMarkdown document
        (`*.Rmd`) and (assuming the path to pandoc is set) is then "knitted" to make an HTML file (`*.html`). The
        RMarkdown file is retained to make it easy for the user to modify the content and convert the output into a
        report of their own. HTML output can be switched off with `rmdreport=F`. Unless `fullreport=F`, a larger report
        with full BUSCO results tables and comparative plots of Missing BUSCO genes will be generated (`*.full.html').

        ## Compilation of ratings and genomes tables for summary plots and tables.

        Prior to generation of the document, results are first compiled into an overview stats file with additional plot
        attributes, before the "best" assemblies are identified under various criteria. Finally, there is an option to
        modify some of the plotting attributes before the report is generated.

        The main BUSCOMP `*.ratings.tdt` output is combined with key genome statistics and BUSCO ratings from the
        `*.genomes.tdt` table. BUSCOMP ratings are converted into percentage values. BUSCO ratings are converted into
        percentage values and reduced to `BUSCO` (Single and Duplicated Complete BUSCOs) and `NoBUSCO` (Missing BUSCOs).
        (Full BUSCO results are plotted directly from the `*.genomes.tdt`.)

        After results are compiled, additional plotting fields are generated:

        * `col` = Plotting colour. (Default, "red". Genomes with 1+ "best" ratings, "blue")
        * `pch` = Point type. (Default, 21 [circle]. Genomes with "best" rating will be 22 [square] for best BUSCO(MP)
        ratings, 24 [triangle] for best contiguity ratings, or 23 [diamond] for best in both categories.
        * `label` = Additional text field to be used for labels. In interactive mode, the option will be given to leave
        this blank for unlabelled plots.
        * `plot` = Whether or not to include a genome in the plot. (Default, TRUE)

        In interactive mode, the option will be provided to edit plotting attributes for each assembly, following the
        calculation of "best" assemblies (below). Compiled data are then saved to `*.LnnIDxx.rdata.tdt` for summary plots
        and tables in the RMarkdown output.

        ## Identifying the "best" assemblies

        There are many ways of assessing genome assembly quality, but they can be broadly broken down into four aspects:

        1. **Coverage.** How much of the genome is included in the assembly.

        2. **Contiguity.** How many fragments is the genome in, and how big are they.

        3. **Accuracy.** How accurate is the assembly in terms of sequence errors.

        4. **Redundancy.** How much of the genome has been included multiple times.

        Standard reference-free genome statistics (e.g. number, length and gappiness of scaffolds), can only give limited
        insights. Whilst assemblies smaller than the genome size are clearly missing bits, assembly size could be
        artificially inflated by redundant regions, making it impossible to assess Coverage. Scaffold sizes can give an
        indicator of Contiguity, but are also prone to biases introduced by filtering of small scaffolds, for example.
        If an estimated Genome Size can be provided, the most useful measures in this respect are `NG50` and `LG50`.
        These statistics sort scaffolds by length (big to small) and identify the contig/scaffold size at which half the
        *genome* (not the *assembly*) is covered by contigs/scaffolds of that size or bigger. This is the `NG50` value
        (called `NG50Length` in BUSCOMP), with `LG50` (`LG50Count`) being the number of contigs/scaffolds needed to
        cover half the genome. (`N50` and `L50` are the same statistics only for the assembly.) These statistics can
        still be mislead by redundancy in the assembly. None of theses statistics speak to sequence accuracy.

        The power of BUSCO is that it speaks to all four aspects of quality. `Complete` and `Fragmented` genes give
        an indication of Coverage, Continuity and Accuracy; `Duplicated` genes give an indication of Redundancy;
        `Missing` genes give an indication of Coverage and Accuracy. However, the weakness is that these different
        aspects cannot easily be disentangled. This makes side-by-side comparisons of different assemblies challenging,
        as it is not always clear what a difference is indicating.

        BUSCOMP is designed on the principle that **Coverage** and **Contiguity** are the two most important aspects of
        assembly quality. *Accuracy* can, in principle, be improved by additional error-correction steps (including
        manual curation). Suspected *Redundancy* can likewise be identified a flagged. *Coverage* and *Contiguity*, in
        contrast, can only be improved by further assembly - either trying again from scratch, or employing a
        scaffolding or gap-filling alogrithm.

        BUSCOMP has therefore identified seven statistics that can be used to rank assemblies on inferred Completeness
        or Contiguity:

        * **Completeness**:
            1. `Complete` = Maximum number of Complete BUSCOMP ratings.
            2. `BUSCO` = Maximum number of Complete BUSCO ratings.
            3. `Missing` = Smallest number of Missing BUSCOMP ratings.
            4. `NoBUSCO` = Smallest number of Missing BUSCO ratings.
        * **Contiguity**:
            1. `MaxLength` = Maximum length of contigs/scaffolds.
            2. `NG50Length` = Half the genome lies on contigs/scaffolds at least this size. (See above.)
            3. `LG50Count` = Half the genome can be covered by this many contigs/scaffolds. (See above.)

        Individual assemblies are rated as "best" under all seven criteria, with ties allowed. Each assembly can
        be best under multiple criteria and each criterion can have several best assemblies. BUSCOMP will report all
        such combinations.

        ---

        ## BUSCOMP Report sections

        The BUSCOMP HTML reports (full and summary) consist of the following sections:

        ### 1. Run summary

        Details of the BUSCOMP run, including the version, directory and commands.

        - 1.1 BUSCOMP summary: A brief overview of the "best" assembly ratings.

        ### 2. Genome summary

        This section contains summary details of the genome assemblies analysed, including the loaded data, the summary
        statistics for the assemblies, and coverage/contiguity assessment plots:

        - 2.1 Genome statistics
        - 2.2 Genome coverage assessment plots
        - 2.3 Genome contiguity assessment plots

        ### 3. BUSCO Ratings

        The compilation of BUSCO ratings, including details of the BUSCOMP Groups is given in this section. In addition
        to the summary ratings for each genome/group, the full report will have a table of all the individual gene
        ratings:

        - 3.1 Genome Groups
        - 3.2 BUSCO Summary
        - 3.3 BUSCO Gene Details

        ### 4. BUSCOMP Ratings

        Next, BUSCOMP ratings using the compiled BUSCOSeq dataset are reported. In addition
        to the summary ratings for each genome/group, the full report will have a table of all the individual gene
        ratings:

        - 4.1 BUSCOMP re-rating of genomes
        - 4.2 BUSCOMP re-rating full results

        ### 5. BUSCO and BUSCOMP Comparisons

        The final report section features direct comparisons of the BUSCO and BUSCOMP ratings, reporting
        changes in ratings between BUSCO and BUSCOMP. _Unique_ Complete genes are also identified: those rated as
        `Complete` in only a single genome, or multiple genomes in a single Group.
        The full report also has a series of summary ratings plots for subsets of BUSCO genes that are rated as
        `Missing` in a genome or group.

        - 5.1 BUSCO to BUSCOMP ratings changes
        - 5.2 Unique BUSCO and BUSCOMP Complete genes
        - 5.3 Ratings for Missing BUSCO genes

        ### 6. Appendix: BUSCOMP run details

        Details of the BUSCOMP run in terms of when, where and how (e.g. commandline options) are found in the Appendix.
        Any errors or warnings generated during the BUSCO run are reported here. Check the `*.log` file generated for
        details.

        - 6.1 BUSCOMP errors
        - 6.2 BUSCOMP warnings

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return self.docHTML()
            ## ~ [1a] ~ Setup the BUSCO runs for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# This will create the basic `genomes` and `alias` tables, and mark assemblies for buscompseq.
            if not self.setup(): return False
            ## ~ [1b] ~ Add additional genomes for rating ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# This will update the `genomes` table with any additional fasta files to process
            if not self.addRateFas(): return False
            ## ~ [1c] ~ BUSCOMP Alias Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# This will compare loaded data to an alias table (if given) and update entries.
            if not self.updateGenomesFromAlias(): return False
            ## ~ [1d] ~ BUSCOMP Genome Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# This will load the groups table and perform initial default group assignments (non-interactive)
            if not self.genomeGroups(review=False): return False
            ## ~ [1e] ~ Manual review of data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# In interactive mode, edit options for assemblies and groups will be offered.
            if not self.genomeMenu(): return False   #x#self.review()
            ## ~ [1f] ~ Load and summarise genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.loadGenomes()

            ### ~ [2] ~ Load and compile BUSCOs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.compileBUSCO(): return False
            #i# If sorting on BUSCO results, sort following calculations
            self.sortGenomes(restrict=['Complete','Missing'],quicksort=True)
            if self.db('full'):
                self.rawSummary()  # Summarise stats from self.db('genomes')
            ## ~ [2a] Main table and outputs from raw BUSCO results ~~~~~~~~~~~~~~~~~~~~~~ ##
            self.db('genomes').saveToFile(sfdict={'GC':4})
            if self.db('full'):
                self.db('full').newKey(['BuscoID','Genome','#'])
                self.db('full').saveToFile()
                self.db('busco').saveToFile()
            ## ~ [2b] Compile BUSCO sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            faswarn = rje.exists(self.getStr('BUSCOFas'))
            if self.db('full'):
                if not self.compileBUSCOSeqs(): return False     # - save as *.buscomp.fasta
            else: faswarn = False

            ### ~ [3] ~ BUSCOMPSeq searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Perform PAF Minimap2 searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('BUSCOMPSeq'):
                # Perform PAF Minimap2 searches
                # - Create $BASEFILE.minimap/ and run rje_paf within
                self.minimapSearches()
            ## ~ [3b] Compile PAF Minimap2 searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Compile PAF searches
                # - Load all GABLAM and HitSum tables
                # - Reclassify
                # - Output to *.buscomp.tdt
                self.compileMinimap()
            ## ~ [3c] Output compiled tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.db('buscoseq'): self.db('buscoseq').saveToFile()
                self.baseFile('%s.%s' % (self.baseFile(),self.pafBase()))
                if self.db('buscomp'): self.db('buscomp').saveToFile()
                self.db('ratings').newKey(['#'])
                self.db('ratings').saveToFile()
                self.db('ratings').newKey(['Genome'])
                #if faswarn:
                #    self.warnLog('Reported "BUSCOMP Seqs" will refer to compiled best sequences, not necessarily those in given buscofas=FASFILE.')
                self.compSummary(full=True)  # Summarise stats from self.db('genomes')
                self.compSummary(full=False)  # Summarise stats from self.db('genomes')
            ## ~ [3d] Changes in BUSCO[MP] ratings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.buscompChanges()
            ## ~ [3e] Identification of Unique BUSCO[MP]s ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.uniqueBUSCOs()
            # Add unique field = Missing
            # cycle through each genome/group
            # - if not Missing and unique Missing -> Genome
            # - if not Missing and unique not Missing -> ''
            # Do separately for genomes and groups. Use genome if possible, else group
            # Summarise counts per genome/group
            # Q. Add genomes to group?

            ### ~ [4] ~ Final HTML output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Tidy up Genome field for Rmd output. NOTE: This might mess up repeat calling of earlier methods.
            for gentry in self.db('genomes').entries():
                gentry['Genome'] = self.genomeField(gentry)
            ## ~ [4a] Rmd output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupRData()   # Compile ratings and genomes tables for plots
            self.buscompRmd(fullreport=False)
            if self.getBool('FullReport'): self.buscompRmd()

            return True
        #except SystemExit: return False
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def docHTML(self):  ### Generate the BUSCOMP Rmd and HTML documents.                                        # v0.5.2
        '''Generate the BUSCOMP Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='BUSCOMP Documentation',author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2019 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown BUSCOMP documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### BUSCOMP Setup and Loading Methods
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                            # v0.5.1
        '''
        STEP 1 is to set up the input data. BUSCOMP will first identify and check the run directories (given with the
        `runs=DIRLIST` command). This command can either be a list of directories, or an expandable list with wildcard.
        A list of individual directories can be provided as a file if desired. NOTE: `runs` wants the actual BUSCO run
        directories (`run_$PREFIX[.busco]/`), not the parent directory. These directories can optionally have a `.busco`
        suffix, which will be trimmed off. By default, BUSCOMP will look for `./run*`.

        From these, the run list will be extracted, consisting of a number of genome $PREFIX values. It is expected that
        the $PREFIX for each run will match the input genome file, found in a directory set by `fastadir=DIRLIST`. This
        file can also match the `Genome` alias for that genome (see below).

        The presence of a `full_$PREFIX[.busco].tsv` table and `single_copy_busco_sequences/` directory within each run
        directory will also be checked. If `buscocompseq=F` then all of the `full_*.tsv` tables in the `runs`
        directories will be loaded for compilation, but no sequence searching will be performed. The presence of
        sequences available for compilation will be stored in the `Sequences` field of `*.genomes.tdt`.

        ### Genome Prefixes

        Each genome analysed has a "Prefix" that is used to help sorting (if `runsort=prefix`) and identify relevant
        BUSCO files where the genome (Fasta) name and BUSCO run name do not match. This Prefix is set once as the data
        is loaded and is never changed. (Sorting and visualisations can be made more informative altered using the Genome
        (a.k.a. Alias) and Description fields.)

        The assembly Prefix will be set as:

        1. If a `full_table_*.tsv` file is found, the core (indicated with `*`) will be used as the prefix.
        2. If no BUSCO results are available, the path- and extension-stripped fasta file name will be used for the prefix

        ### Genome Aliases

        Genome Aliases will be parsed initially from the loaded data:

        1. If a BUSCO run directory is given and it contains a single BUSCO run, the run directory name will be used for
        the `Genome` alias, minus the leading "`run_`".
        2. If the run directory contains multiple BUSCO runs, the `full_table_*.tsv` core used for the `Prefix` will also
        be used for the `Genome` alias.
        3. If no BUSCO results are available, the path- and extension-stripped fasta file name will be used for the
        `Genome` alias.

        If `Genome` contains a trailing `.busco`, this will be stripped. If `Genome` starts with preceding numbers
        (`XX_`), these will be kept for sorting but stripped for outputs and fasta matching (unless `stripnum=F`).

        **NOTE:** If the assembly fasta file does not match the name used for the BUSCO run file, then the `Genome` alias
        *must* be set to the `basefile` of the fasta file by appropriate naming of the `run_` directory, _i.e._ the
        `run_*` directory should contain a single BUSCO run and the directory name should match a `*.fasta` file (or
        other accepted `fastaext=LIST` extension) in one of the directories given by `fastadir=DIRLIST`.

        An optional alias table can be provided (`genomes=FILE`) that contains `Prefix` and `Genome` fields (others
        allowed) and can give alternative names for each genome to be used in other outputs. If running interactively,
        there will also be an option to check/update these aliases manually. Genome aliases are output as part of the
        main `*.genomes.tdt` output, which can also be used as subsequent `genomes=FILE` input (the default, if found).
        This table can also have a `Description` field, which will be used to update descriptions for output if found.
        (Empty descriptions will not be mapped.)

        **NOTE:** The optional Alias table *can* be used to change a `Genome` name to something other than its Fasta file
        - this mapping is performed _after_ fasta files have been identified.

        **NOTE:** `Prefix` and `Genome` fields must be unique. A common `Prefix` and `Genome` value is permitted but
        these will be assumed to refer to the _same_ assembly. (In other words, do not run BUSCO on two assemblies
        and give each BUSCO run the name of the other assembly!)

        ## Table: genomes.tdt
        - # = Sorting order for output
        - Directory = Directory for this BUSCO run
        - Prefix = Prefix identified from BUSCO directory (or full table), corresponding to BUSCO output files
        - Genome = Optional alternative name (or "alias") to be used in outputs
        - Fasta = path to genome file, if found
        - Genome statistics = genome summary statistics if genome found and `summarise=T`
        - Sequences (True/False) = Whether `single_copy_busco_sequences/` directory found
        - Table = Path to `full_*.tsv` table
        - Ratings summary will be added (assuming `Table` is `True`)

        ## Table: alias (not saved: loaded from *.genomes.tdt)
        - Prefix = BUSCO run prefix. Where BUSCO run not performed, generic group/run identifier
        - Genome = unique user-friendly label (alias) for assembly

        NOTE: Not entirely clear how this is used in the code.

        '''
        try:### ~ [1] Setup DB ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ SETUP BUSCOMP RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('SETUP BUSCOMP RUN',line='=')
            db = self.db()

            ### ~ [2] Setup Runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Setup table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            genfields = ['#','Directory','Prefix','Genome','Fasta','Sequences','Table','Description']
            #i# NOTE: if summarise=T, will add Genome Statistics after 'Fasta'
            gdb = db.addEmptyTable('genomes',genfields,['#'],log=self.debugging())
            ## ~ [2b] Read runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dirx = 0; nonx = 0;
            for runpath in self.list['Runs']:
                # Identify and check directory
                runtype = self.runType(runpath)
                self.printLog('#RUNDIR','%s run type: %s' % (runpath,runtype))
                if not runtype:
                    self.warnLog('No BUSCO runs found for run path "%s"' % runpath)
                    continue
                #if not os.path.isdir(runpath): continue
                runsplit = os.path.split(runpath)
                if not runsplit[-1]: runsplit = os.path.split(runsplit[0])
                rundir = runsplit[-1]
                dirx += 1

                # V0.9.0 updated and simplified input data parsing strategy: BUSCO v3
                # Input BUSCO runs are provided by the directories given by `runs=DIRLIST`. These must have a `full_table_*.tsv` files
                # for each BUSCO genome. The `*` component of the name is used as the genome `Prefix` - a unique identifier
                # used for each genome. By default, this is also used to set the unique `Genome` identifier. If the run directory
                # has a `run_` prefix AND there is only one `full_table`, the `Genome` identifier will instead be parsed from
                # the run directory name. This will have any `.busco` suffix removed. Fasta files must match either the `Prefix` or `Genome` identifier loaded.
                # Either or both of the `Prefix` and `Genome` identifiers can have leading numbers (`NN_*`) for sorting purposes.
                # These will be ignored when finding genome fasta files and/or printing to output.

                # If there is a matching `single_copy_busco_sequences_$PREFIX` directory (exactly matching the `$PREFIX` parsed
                # from a `full_table_*.tsv` file), single copy BUSCO sequences will be used for the BUSCOMP compilation.
                # Alternatively, if the run directory has a `run_` prefix AND there is only one `full_table`, sequences will
                # be extracted from the `single_copy_busco_sequences` directory if present.

                # Similarly, a description for the genome will be parsed from a `description_$PREFIX.txt` file (exactly matching the `$PREFIX` parsed
                # from a `full_table_*.tsv` file) if found. Alternatively, if the run directory has a `run_` prefix AND there is only one `full_table`,
                # a description will be parsed from `description.txt`. In each case, only the first line will be used.
                # If neither condition is met, the name of the genome fasta file will be used. Failing that, `Genome` will be used.

                # The final `Genome` and `Description` values can be updated by mapping the `Prefix` or fasta file onto data loaded from the `genomes=FILE` table.
                # = updateGenomesFromAlias()

                # Extract run file(s)
                self.progLog('#RUNDIR','Processing %s...' % rundir)
                prefix = ''
                runlist = glob.glob('%s/full_table_*.tsv' % runpath)

                #i# NOTE: This might be very confusing because prefix refers to the prefix of the assembly files, whereas
                #i# genome refers to the alias being used for the genome, NOT the prefix of the genome assembly. (Unless
                #i# the two are the same!

                # V4 conversion
                if runtype == 'V4':
                    runpath = glob.glob('%s/run_*' % runpath)[0]
                    runlist = ['%s/full_table.tsv' % runpath]

                # Process run file(s)
                for runfile in runlist:
                    #i# Single BUSCO sequence directory
                    seqdir = rje.makePath(runpath) + 'single_copy_busco_sequences'
                    if runtype == 'V4': seqdir = rje.makePath(runpath) + 'busco_sequences/single_copy_busco_sequences'
                    #?# Replace Sequences T/F with number of *.fna files
                    rentry = {'#':gdb.entryNum()+1,'Table':runfile,'Sequences':runtype in ['Single','V4'] and rje.exists(seqdir)}
                    if rentry['Sequences']: self.printLog('#SEQDIR','Found SC BUSCO sequence directory: %s' % seqdir)
                    (rentry['Directory'],tsv) = os.path.split(runfile)
                    if runtype == 'V4':
                        prefix = rundir
                        if rundir.startswith('run_'): prefix = rundir[4:]
                    else:
                        prefix = tsv[len('full_table_'):-4]
                        prefseqdir = rje.makePath(runpath) + 'single_copy_busco_sequences_' + prefix
                        if rje.exists(prefseqdir):
                            rentry['Sequences'] = True
                            self.printLog('#SEQDIR','Found SC BUSCO sequence directory: %s' % prefseqdir)
                    #self.debug(rentry)
                    rentry['Prefix'] = prefix
                    # Set genome value: chop off .busco at end!
                    genome = prefix
                    # If runpath starts with run and has one tsv then use as alias.
                    #!# Add stripaliasnum option for output.
                    if runtype == 'Single' and rundir.startswith('run_'):   #X# and not prefix in adb.index('Prefix',force=True,log=False):
                        genome = rundir[4:]
                    if genome.endswith('.busco'): genome = genome[:-6]
                    rentry['Genome'] = string.join(string.split(genome),'')
                    rentry['Fasta'] = self._findFasta(prefix,genome)
                    descfile = rje.makePath(runpath) + 'description.txt'
                    prefdescfile = rje.makePath(runpath) + 'description_' + prefix + '.txt'
                    if rje.exists(prefdescfile):
                        rentry['Description'] = rje.chomp(open(prefdescfile,'r').readline())
                    elif runtype in ['Single','V4'] and rje.exists(descfile):
                        rentry['Description'] = rje.chomp(open(descfile,'r').readline())
                    elif rentry['Fasta']: rentry['Description'] = rentry['Fasta']
                    else:
                        self.warnLog('No Fasta file found for %s' % self.genomeString(rentry))
                        rentry['Description'] = rentry['Genome'] + ' (no Fasta)'
                    if not rentry['Sequences'] and runtype in ['Single','V4']:
                        self.warnLog('No SC BUSCO sequence directory identified for %s' % genome)
                    self.debug('%s' % rentry)
                    gdb.addEntry(rentry)
                # Summarise entry
                self.printLog('#RUNDIR','Processed %s: %d runs (%d total)' % (rundir,len(runlist),gdb.entryNum()))

            if not self.list['Runs']:
                self.printLog('#RUNDIR','No run directories provided: no BUSCO results loaded.')
            elif not gdb.entryNum():
                self.warnLog('Processed %d run directories but no BUSCO runs loaded. Check input setting.' % (len(self.list['Runs'])))

            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def runType(self,runpath):  ### Returns the type of BUSCO run directory given by runpath
        '''
        Returns the type of BUSCO run directory given by runpath.
        - Single = Single V3 BUSCO run (default expectation)
        - Collated = Multiple BUSCO full_table*tsv files in a single directory. (Can be V3 or V4)
        - V4 = Single V4 BUSCO run directory.
        - None = No BUSCO results detected.
        :param runpath:
        :return: runtype [str]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not os.path.isdir(runpath): return None
            runlist = glob.glob('%s/full_table*tsv' % runpath)
            if not runlist:
                self.debug('%s/run_*/full_table.tsv' % runpath)
                self.debug(glob.glob('%s/run_*/full_table.tsv' % runpath))
            if len(runlist) == 1 and runlist[0] == '%s/full_table.tsv' % runpath:
                self.warnLog('Possible internal v4 results directory: %s' % runpath)
                return None
            elif len(runlist) == 1: return 'Single'
            elif not runlist and len(glob.glob('%s/run_*/full_table.tsv' % runpath)) == 1:
                return 'V4'
            elif not runlist and len(glob.glob('%s/run_*/full_table.tsv' % runpath)) > 1:
                self.warnLog('Possible v4 results directory but multiple run_*/full_table.tsv files: %s' % runpath)
                return None
            elif not runlist:
                self.warnLog('No BUSCO results tables found in directory: %s' % runpath)
                return None
            else: return 'Collated'
        except: self.errorLog('Problem during %s runType(%s).' % (self.prog(),runpath)); return None  # Setup failed
#########################################################################################################################
    def addRateFas(self):    ### Adds additional fasta files for BUSCOMP rating.                                # v0.5.1
        '''
        Adds additional fasta files for BUSCOMP rating. Any already loaded with BUSCO data will be skipped.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['RateFas']: self.headLog('Adding fasta files for rating')
            gdb = self.db('genomes')
            ## ~ [1b] Add genomes if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gaddx = 0
            for ratefas in self.list['RateFas']:
                prefix = rje.baseFile(ratefas,strip_path=True)
                #x#genome = self.alias(prefix)
                genome = string.join(string.split(prefix),'')
                gentry = {'Prefix':prefix,'Genome':genome,'Fasta':ratefas,'Sequences':False}
                if ratefas in gdb.dataList(gdb.entries(),'Fasta'):
                    self.printLog('#SKIP','RateFas file "%s" skipped - already in input' % ratefas)
                    continue
                if not rje.exists(ratefas):
                    self.warnLog('RateFas file "%s" skipped - file not found!' % ratefas)
                    continue
                gentry['Description'] = gentry['Fasta'] + ' (no BUSCO)'
                gdb.addEntry(gentry)
                gaddx += 1
            if gaddx:
                self.printLog('#RATE','Added %s additional genomes for BUSCOMPSeq rating' % rje.iStr(gaddx))
            if not gdb.entryNum():
                self.errorLog('No genomes loaded for BUSCO complitation or rating. Check input, including ini=INFILE settings.',printerror=False)
                return False
            return True     # Setup successful
        except: self.errorLog('Problem during %s addRateFas.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def genomeTableFormat(self,table):
        iformats = {}
        for ifield in string.split('SeqNum	TotLength	MinLength	MaxLength	N50Length	L50Count CtgNum N50Ctg	L50Ctg	NG50Length	LG50Count GapLength GapCount') + ['UniqBUSCO','UniqBUSCOMP']:
            iformats[ifield] = 'int'
        for nfield in string.split('MeanLength MedLength GC'):
            iformats[nfield] = 'num'
        table.dataFormat(iformats)
#########################################################################################################################
    def updateGenomesFromAlias(self):   ### Compare loaded data to an alias table (if given) and update entries.# v0.5.1
        '''
        Compare loaded data to an alias table and update entries.
        This will also generate the self.dict['Alias'] dictionary.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ MAPPING GENOME ALIASES ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('MAPPING GENOME ALIASES',line='=')
            db = self.db()
            gdb = self.db('genomes')
            ## ~ [1a] Setup aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('Genomes'): self.setStr({'Genomes':'%s.genomes.tdt' % self.baseFile()})
            adb = db.addTable(self.getStr('Genomes'),mainkeys=['Genome'],name='alias',expect=False,ignore=[])
            if not adb: adb = db.addEmptyTable(fields=['Prefix','Genome'],keys=['Genome'],name='alias')
            self.dict['Alias'] = {}
            prefixes = adb.index('Prefix').keys()
            if 'Fasta' in adb.fields(): adb.index('Fasta')
            if 'Directory' in adb.fields(): adb.index('Directory')

            ### ~ [2] Map and Set aliases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gentry in gdb.entries(sorted=True):
                # Generate full list of aliases to check
                prefcheck = [gentry['Prefix'],gentry['Genome']]
                if gentry['Fasta']: prefcheck.append(rje.stripPath(gentry['Fasta']))
                for fcore in [gentry['Prefix'],gentry['Genome']]:
                    if self.getBool('StripNum') and rje.matchExp('^(\d+)_(\S+)',fcore):
                        prefcheck.append(rje.matchExp('^(\d+)_(\S+)',fcore)[1])
                for prefix in prefcheck[0:]: prefcheck.append('%s.busco' % prefix)
                # Identify alias if present
                alias = None
                desc = None
                for prefix in prefcheck:
                    if prefix not in prefixes: continue
                    prentries = adb.indexEntries('Prefix',prefix)
                    if not prentries: continue
                    if len(prentries) > 1: raise ValueError('Alias table has multiple Prefix entries for "%s"' % prefix)
                    if alias and alias != prentries[0]['Genome']:
                        raise ValueError('Alias table has multiple different mappings for %s' % self.genomeString(gentry))
                    if not alias:
                        alias = prentries[0]['Genome']
                        self.printLog('#ALIAS','%s: "%s" -> "%s"' % (self.genomeString(gentry),prefix,alias))
                    if 'Description' in adb.fields() and not desc: desc = prentries[0]['Description']
                # Map via Fasta or Directory
                for mfield in ['Fasta','Directory']:
                    if not gentry[mfield]: continue
                    if not mfield in adb.fields(): continue
                    if gentry[mfield] not in adb.index(mfield): continue
                    mentries = adb.indexEntries(mfield,gentry[mfield])
                    if len(mentries) > 1: raise ValueError('Multiple %s mappings for %s' % (mfield,self.genomeString(gentry)))
                    mentry = mentries[0]
                    if alias and mentry['Genome'] != alias: raise ValueError('Alias/%s mapping mismatch for %s' % (mfield,self.genomeString(gentry)))
                    alias = mentry['Genome']
                    if 'Description' in adb.fields() and not desc: desc = mentry['Description']

                # Update genome entry
                if alias:
                    gentry['Genome'] = alias
                    self.printLog('#ALIAS','Genome updated: %s' % self.genomeString(gentry))
                elif self.getBool('Restrict'):
                    gdb.dropEntry(gentry)
                    self.printLog('#DROP','Genome without alias dropped from analysis (restrict=T): %s' % self.genomeString(gentry))
                if desc:
                    gentry['Description'] = desc
                    self.printLog('#DESC','Description updated: %s = "%s"' % (self.genomeString(gentry),desc))
                # Update alias dictionary - checks for redundancy
                #?# Can we get away without using self.dict['Alias']?
                for afield in ['Prefix','Genome']:
                    if afield == 'Genome' and gentry[afield] == gentry['Prefix']: continue
                    if gentry[afield] in self.dict['Alias']:
                        raise ValueError('Assembly "%s" value "%s" already in use for %s!' % (afield, gentry[afield], self.genomeString(self.dict['Alias'][gentry[afield]])))
                    self.dict['Alias'][gentry[afield]] = gentry

            return True     # No rules broken
        except: self.errorLog('Problem during %s updateGenomesFromAlias.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def genomeMenu(self):   ### Main genome and group review/setup method.                                      # v0.5.1
        '''
        ## BUSCOMP Genome/Group menu

        Following initial setup of genomes and groups, a Genome/Group edit menu enables updating loaded information.
        (This menu will not be available if running in non-interactive mode, e.g. `i=-1`.) The primary use for this is
        to convert unwieldy genome and/or BUSCO run names into something more user-friendly for display. It can also be
        used for adding descriptions to assemblies, adding or modifying group membership, and sorting genomes using
        different criteria. Assemblies can also be dropped from analysis from this menu.

        The main menu offers five options:

        - `R` = Review genomes
        - `G` = Review/edit genome groups
        - `S` = Sort genomes
        - `Q` = Quit BUSCOMP
        - `X` = Exit Menu and Proceed [default]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headtext = 'BUSCOMP Genome/Group menu'
            edited = False
            gdb = self.db('genomes')
            self.sortGenomes(quicksort=True)
            if self.i() < 0: return True
            default = 'R'
            nextdefault = {'R':'G','G':'S','S':'X','X':'X'}

            ### ~ [2] Review menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            menulist = [('R','Review genomes','return','R'),
                        ('G','Review/edit genome groups','return','G'),
                        ('S','Sort genomes','return','S'),
                        ('Q','Quit BUSCOMP','return','Q'),
                        ('X','Exit Menu and Proceed','return','X')]
            while gdb.entryNum():
                choice =  rje_menu.menu(self,headtext,menulist,choicetext='Please select:',changecase=True,default=default,jointxt=' = ',confirm=False)
                if choice == 'R': edited = self.review()
                if choice == 'G': self.genomeGroups()
                if choice == 'S': self.sortGenomes(quicksort=False)
                if choice == 'X' and rje.yesNo('Exit genome/group menu and proceed?'): break
                if choice == 'Q' and rje.yesNo('Quit BUSCOMP?'): os.sys.exit(0)
                if choice in 'RGS': default = nextdefault[default]

            if not gdb.entryNum():
                raise IOError('No genomes left to process!')

            return True
        except: self.errorLog('Problem during %s genomeMenu.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def review(self):    ### Options to review and edit loaded assemblies.                                      # v0.5.1
        '''
        ### Genome review menu

        If running in interactive mode, users will have the option to review and manually update loaded data. Each
        loaded genome will be displayed in turn, in its current sorted order. (Sorting can be altered or updated via
        the `Sort genomes` option of the main genome/group menu.) Options will be presented to modify data, or quit
        the menu (or BUSCOMP program). Modifications permitted are:

        - `<G>enome alias` = Sets the `Genome` field for the entry. This can have a "silent" numeric prefix for sorting if
        `runsort=Genome` or `runsort=Group`. Most importantly, this will be the default label in tables and plots, and
        used for fields in compiled data.
        - `<D>escription` = This gives the option to set a manual description of the assembly, which is used in the
        final output.
        - `<R>emove` = Drop this assembly from the analysis. BUSCOMP will ask for confirmation. This cannot be undone.

        If genomes are edited and `runsort=Genome` or `runsort=Group`, they will be re-sorted upon exiting the menu.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db('genomes')
            grpdb = self.db('groups')
            #x# Replaced by menu: #self.sortGenomes(quicksort=True)

            ### ~ [2] Review ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edited = False
            if self.i() >= 0 and rje.yesNo('Review BUSCO datasets?'):
                #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ MANUAL BUSCO DATASET REVIEW ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
                self.headLog('MANUAL BUSCO DATASET REVIEW',line='=')
                i = 0; sum = True
                while gdb.entryNum():
                    default = 'X'
                    gentry = gdb.entries(sorted=True)[i]
                    gtext = ''
                    if sum: gtext += '\n' + gdb.entrySummary(gentry)
                    sum = True
                    gtext += '\n'
                    if i > 0: gtext += '<P>revious | '
                    if i+1 < gdb.entryNum(): gtext += '<N>ext | '; default = 'N'
                    gtext += '<G>enome alias | <D>escription | <R>emove | E<x>it Review | <Q>uit'
                    if gentry['Sequences']: gtext = string.replace(gtext,'<R>','Ignore <S>equences | <R>')
                    ## Loop through choices ##
                    choice = rje.choice(gtext,default=default).upper()
                    if choice == 'P': i -= 1
                    elif choice == 'N': i += 1
                    elif choice == 'G':
                        oldgenome = gentry['Genome']
                        gentry['Genome'] = rje.choice('Genome Alias?:',default=gentry['Genome'],confirm=True,whitespace=False)
                        while gentry['Genome'] not in [oldgenome,gentry['Prefix']] and gentry['Genome'] in self.dict['Alias']:
                            gentry['Genome'] = rje.choice('Genome Alias already in use. New alias?:',default=oldgenome,confirm=True,whitespace=False)
                        edited = gentry['Genome'] != oldgenome
                        if gentry['Genome'] != oldgenome: self.dict['Alias'][gentry['Genome']] = self.dict['Alias'].pop(oldgenome)
                        if gentry['Genome'] != oldgenome and grpdb and oldgenome in grpdb.index('Genome'):
                            for grpentry in grpdb.indexEntries('Genome',oldgenome): grpentry['Genome'] = gentry['Genome']
                            grpdb.index('Genome',force=True,log=False)
                    elif choice == 'D':
                        gentry['Description'] = rje.choice('Enter description',default=gentry['Description'],confirm=True)
                    #!# Add option to re-add sequences and scan for folder #!#
                    elif gentry['Sequences'] and choice == 'S' and rje.yesNo('Drop %s sequences from BUSCOMPSeq compilation?' % gentry['Genome']) and rje.yesNo('Are you sure? This cannot be undone.'):
                        gentry['Sequences'] = False
                    elif choice == 'R' and rje.yesNo('Drop %s from BUSCOMP analysis?' % gentry['Genome']):
                        gdb.dropEntry(gentry)
                        if gentry['Prefix'] in self.dict['Alias']: self.dict['Alias'].pop(gentry['Prefix'])
                        if gentry['Genome'] in self.dict['Alias']: self.dict['Alias'].pop(gentry['Genome'])
                        for group in self.dict['Groups']:
                            if gentry in self.dict['Groups'][group]: self.dict['Groups'][group].remove(gentry)
                    elif choice == 'X' and rje.yesNo('End review of genomes?'): break
                    elif choice == 'Q' and rje.yesNo('Quit BUSCOMP?'): os.sys.exit(0)
                    else: sum = False
            if edited: self.sortGenomes(quicksort=True,restrict=['Genome','Group'])

            return edited     # Return whether or not genomes have been edited
        except: self.errorLog('Problem during %s review.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def sortGenomes(self,quicksort=False,restrict=[]):    ### Reorder genomes prior to output.                  # v0.5.1
        '''
        >> quicksort:bool [False] = Whether to only conduct Genome/Prefix/Complete/Missing sorting
        >> restrict:list [] = Restrict sorting to these methods only.

        ### Output order

        An optional `runsort=LIST` parameter will determine the output order. Future releases of BUSCOMP will proceess
        a mixed list of Prefix, Genome and Group names. Any data missing from the list will be appended at the end in
        Genome order. (NOTE: Group names are a special kind of Genome in this context.) This sorting is preformed after
        all manual reviews and updates are complete. Currently, providing a list of genomes (or `runsort=Manual`) will
        trigger manual ordering via the sort menu.

        There are also four special sort modes that can be implemented:

        - Genome = sort alphabetically based on Genome alias and Group name. (NOTE: "silent" numerical prefixes will
        affect this sorting (see below).
        - Group = sort alphabetically using `Genome` sorting within each group with the Group summary immediately
        following the last Group member. [Default]
        - Prefix = sort alphabetically based on Prefix. (NOTE: "silent" numerical prefixes will
        affect this sorting (see below).
        - Complete = sort in order of most complete to least complete.
        - Missing = sort in order of fewest missing to most missing BUSCOs

        The `runsort=prefix` option is best combined with the addition of a number to the run directory name:

        * `run_01_OneGenome/` containing `full_OneGenome.tsv` and generated from `OneGenome.fasta`
        * `run_02_AnotherGenome.busco/` containing `full_AnotherGenome.busco.tsv` and generated from `AnotherGenome.fasta`

        This example would sort OneGenome before AnotherGenome if runsort=prefix, but AnotherGenome before OneGenome if
        `runsort=Genome`.

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gdb = self.db('genomes')
            quicksortlist = ['Genome','Prefix','Complete','Missing','Group']
            runsort = self.list['RunSort']
            if len(runsort) == 1: runsort = rje.strSentence(runsort[0].lower())
            if runsort not in quicksortlist + ['Manual']:
                self.printLog('#DEV','Custom sorting not yet implemented; set runsort=Manual')
                runsort = 'Manual'
                self.list['RunSort'] = ['Manual']
            if runsort in ['Groups','Genomes']: runsort = runsort[:-1]
            #self.printLog('#SORT',string.join(self.list['RunSort']),', ')
            ## ~ [1a] Sort Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            headtext = 'Select sorting option'
            sortopt = {'Group':'G','Genome':'L','Prefix':'P','Complete':'C','Missing':'M','Manual':'E'}
            menulist = [('G','Group (Genome label then group names)','return','Group'),
                        ('L','Genome label','return','Genome'),
                        ('P','Run prefix','return','Prefix'),
                        ('C','BUSCO complete','return','Complete'),
                        ('M','BUSCO missing','return','Missing'),
                        ('E','Edit sort (manual)','return','Manual'),
                        ('Q','Quit BUSCOMP','return','Quit')]
            if self.i() >=0 and not quicksort:
                choice =  rje_menu.menu(self,headtext,menulist,choicetext='Please select:',changecase=True,default=sortopt[runsort],jointxt=' = ',confirm=True)
                if choice == 'Q' and rje.yesNo('Quit BUSCOMP?'): os.sys.exit(0)
                self.list['RunSort'] = [choice]
                runsort = choice

            ### ~ [2] QuickSort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if restrict: quicksortlist = rje.listIntersect(quicksortlist,restrict)
            if runsort in quicksortlist and runsort in gdb.fields():
                #gdb.rankField(runsort,'RankSort',rev=runsort in ['Complete'],unique=True)
                gdb.rankField(runsort,'RankSort',rev=runsort in ['Complete'],unique=True,warn=False)
                for gentry in gdb.entries(): gentry['#'] = gentry['RankSort']
                gdb.dropField('RankSort',log=self.debugging())
                gdb.remakeKeys()
                if quicksort: return True

            ### ~ [3] Group sorting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if runsort == 'Group':
                # Generate sorted list of genomes
                #!# Add checks for mutliple matching genome names at some point! (Could introduce with edits) #!#
                genomes = gdb.dataList(gdb.entries(),'Genome',sortunique=True)
                grpdb = self.db('groups')
                groups = rje.sortKeys(grpdb.index('Group'))
                for group in groups:
                    if group in genomes: genomes.remove(group)
                    groupi = 0
                    for genome in grpdb.indexDataList('Group',group,'Genome'):
                        if genome in genomes:
                            groupi = max(groupi, genomes.index(genome))
                        else:
                            self.warnLog('Group "%s" genome "%s" not found in genomes list' % (group,genome))
                    genomes.insert(groupi+1,group)
                # Regenerate sorted list
                newdata = {}
                for gentry in gdb.entries():
                    gentry['#'] = genomes.index(gentry['Genome']) + 1
                    newdata[gentry['#']] = gentry
                gdb.dict['Data'] = newdata
                self.printLog('#SORT','Genomes sorted by genome name and group')
                if quicksort: return True

            ### ~ [4] Manual sorting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if quicksort: return False
            sorted = gdb.orderedDataList('Genome')
            sortxt = '\n\nGenome Sort order (%s sorting):\n' % runsort
            for genome in sorted: sortxt += '%3s: %s\n' % (sorted.index(genome)+1,genome)
            sortxt += '\n'
            if runsort == 'Manual':
                self.printLog('#DEV','Manual run sorting not yet implemented.')
            elif not quicksort:
                self.verbose(0,1,sortxt)

            return False     # Sorting not changed
            #return True     # Sorting (maybe) changed
        except SystemExit: raise
        except: self.errorLog('Problem during %s sortGenomes.' % self.prog()); raise
#########################################################################################################################
    def loadGenomes(self):    ### Load genomes into SeqList objects and add summaries to genomes table.         # v0.5.1
        '''
        Load genomes into SeqList objects and add summaries to genomes table.

        ## Genome Summary Statistics

        The final setup step is to load any genomes for which Fasta files have been provided. Unless `summarise=F`,
        summary statistics will be calculated for each genome, enabling assessements of genome completeness and quality
        in combination with the BUSCO and BUSCOMP results. Summary statistics are calculated by `rje_seqlist`. The
        following statistics are calculated for each genome:

        * **SeqNum**: The total number of scaffolds/contigs in the assembly.
        * **TotLength**: The total combined length of scaffolds/contigs in the assembly.
        * **MinLength**: The length of the shortest scaffold/contig in the assembly.
        * **MaxLength**: The length of the longest scaffold/contig in the assembly.
        * **MeanLength**: The mean length of scaffolds/contigs in the assembly.
        * **MedLength**: The median length of scaffolds/contigs in the assembly.
        * **N50Length**: At least half of the assembly is contained on scaffolds/contigs of this length or greater.
        * **L50Count**: The smallest number scaffolds/contigs needed to cover half the the assembly.
        * **NG50Length**: At least half of the genome is contained on scaffolds/contigs of this length or greater.
        This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
        * **LG50Count**: The smallest number scaffolds/contigs needed to cover half the the genome.
        This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
        * **GapLength**: The total number of undefined "gap" (`N`) nucleotides in the assembly.
        * **GC**: The %GC content of the assembly.

         **NOTE:** `NG50Length` and `LG50Count` statistics use `genomesize=X` or the biggest assembly loaded.
         If BUSCOMP has been run more than once on the same data (_e.g._ to update descriptions or sorting),
         previously calculates statistics will be read from the `genomes=FILE` table, if loaded. This may cause some
         inconsistencies between reported NG50 and LG50 values and the given/calculated `genomesize=X`. If partial
         results are reloaded following a change in `genomesize=X`, this will give rise to inconsistencies between
         assemblies. When re-running, please make sure that a consistent genome size is used.
         If in doubt, run with `force=T` and force regeneration of statistics.

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            statfields = string.split('SeqNum TotLength MinLength MaxLength MeanLength MedLength N50Length L50Count CtgNum N50Ctg L50Ctg NG50Length LG50Count GapLength GapCount GC')
            gdb = self.db('genomes')
            seqobj = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=F','seqin=None'])
            if not 'SeqNum' in gdb.list['Fields']:
                gdb.list['Fields'] += statfields
            self.genomeTableFormat(gdb)
            ## ~ [1a] Fasta files to process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fasgenomes = [] # List of genomes to calculate stats for - will be empty if all found in alias table
            adb = self.db('alias')
            adata = []      # List of loaded Fasta files
            if adb and 'SeqNum' in adb.fields() and 'Fasta' in adb.fields():
                adb.dataFormat(dformats)
                adata = adb.orderedDataList('Fasta',empties=False)
                adb.index('Fasta')
            if adb: self.genomeTableFormat(adb)
            for gentry in gdb.entries(sorted=True):
                if gentry['Fasta']:
                    if gentry['Fasta'] in adata and not self.force():    # Data should already exist
                        aentry = adb.indexEntries('Fasta',gentry['Fasta'])[0]
                        for field in statfields:
                            if field in aentry:
                                gentry[field] = aentry[field]
                            else: gentry[field] = ''
                        adata.remove(gentry['Fasta'])
                    else: fasgenomes.append(gentry['Fasta'])

            ### ~ [2] Longest assembly as genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #seqlist = self.dict['BUSCOMP'] = {}   # SeqList object containing compiled BUSCO hit sequences
            self.setInt({'GenomeSize':seqobj.getNum('GenomeSize')})
            if not seqobj.getNum('GenomeSize') > 0:
                self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ ASSEMBLY SIZES ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
                maxlen = 0
                lenfas = []
                for gentry in self.db('genomes').entries():
                    if not gentry['Fasta']: continue  # Cannot summarise
                    if gentry['Fasta'] in fasgenomes: lenfas.append(gentry['Fasta'])
                    else: maxlen = max(maxlen,gentry['TotLength'])
                if maxlen:
                    self.warnLog('Using loaded stats data - might need to check for genome size consistency for NG50 and LG50.')
                for fasfile in lenfas:
                    genlen = 0
                    self.progLog('#FASTA','Reading sequences from %s' % fasfile)
                    IN = open(fasfile,'r')
                    iline = IN.readline()
                    while iline:
                        if iline[:1] != '>': genlen += len(rje.chomp(iline))
                        iline = IN.readline()
                    IN.close()
                    maxlen = max(maxlen,genlen)
                    self.printLog('#FASTA','Read sequences from %s: %s bp' % (fasfile,rje.iStr(genlen)))
                self.printLog('#GSIZE','Using biggest assembly for genome size: %s' % rje_seqlist.dnaLen(maxlen))
                self.cmd_list.append('genomesize=%d' % maxlen)
                self.setInt({'GenomeSize':maxlen})

            ### ~ [3] Load Genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ LOAD GENOMES ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('LOAD GENOMES',line='=')
            for gentry in self.db('genomes').entries():
                if not gentry['Fasta']: continue  # Cannot summarise
                if 'SeqNum' in gentry and gentry['SeqNum'] and not self.force() and self.getBool('LoadSummary'):
                    self.printLog('#FASTA','%s fasta summarised (force=F)' % self.genomeString(gentry))
                    continue # Already summarised
                self.printLog('#FASTA','Loading sequences for %s' % self.genomeString(gentry))
                seqcmd = ['seqin=%s' % gentry['Fasta'],'autoload','seqmode=file','summarise=F','dna','autofilter=F']
                try:
                    #seqobj = seqlist[gentry['Genome']] = rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd)
                    seqobj = rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd)
                    seqobj.checkNames()
                    if self.getBool('Summarise'):
                        sumdata = seqobj.summarise(save=False)
                        sumdata['GC'] = sumdata.pop('GCPC')
                        gentry.update(sumdata)
                except:
                    self.errorLog('Something went wrong processing %s' % gentry['Fasta'])

            gdb.fillBlanks(log=False)

            return True     # Setup successful
        except: self.errorLog('Problem during %s loadGenomes.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <4> ### BUSCOMP Data Mapping Methods
#########################################################################################################################
    def _findFasta(self,prefix,alias): ### Find genome fasta file from given prefix and alias                   # v0.5.1
        '''
        Find genome fasta file from given prefix and alias by looking in FastaDir for with possible file extensions.
        >> prefix:str = possible $PREFIX.$EXT fasta file
        >> alias:str = possible $ALIAS.$EXT fasta file
        << fasta:str = Fasta file. None if not found.
        '''
        try:### ~ [1] Look for file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fcores = [prefix, alias]
            for fcore in [prefix, alias]:
                if self.getBool('StripNum') and rje.matchExp('^(\d+)_(\S+)',fcore):
                    fcores.append(rje.matchExp('^(\d+)_(\S+)',fcore)[1])
            for fcore in fcores[0:]:
                if fcore.endswith('.busco'): fcores.append(fcore[:-6])
            for fcore in fcores:
                for fpath in self.list['FastaDir']:
                    fpath = rje.makePath(fpath)
                    fasfile = fpath + fcore
                    if rje.exists(fasfile): return fasfile
                    for fext in self.list['FastaExt']:
                        fasfile = fpath + fcore + '.' + fext
                        if rje.exists(fasfile): return fasfile
            return None
        except: self.errorLog('%s._findFasta error' % self.prog())
#########################################################################################################################
    def aliasEntry(self,alias,expect=True):  ### Returns genomes table entry from alias (or prefix).             # v0.5.1
        '''
        Returns genomes table entry from alias (or prefix). Will also try add .busco and stripping number prefixes.
        >> alias:str = alias (or prefix) to fish out of genomes table.
        >> expect:bool [True] = whether to expect an entry and raise error if missing: else return None
        << gentry:dict = Genomes Table entry (or None)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not alias:
                if expect: raise ValueError('No alias given to aliasEntry()')
                return None
            ### ~ [2] Find entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alias in self.dict['Alias']: return self.dict['Alias'][alias]
            buscadd = '%s.busco' % alias
            if buscadd in self.dict['Alias']: return self.dict['Alias'][buscadd]
            for fcore in [alias,buscadd]:
                if self.getBool('StripNum') and rje.matchExp('^(\d+)_(\S+)',fcore):
                    astrip = rje.matchExp('^(\d+)_(\S+)',fcore)[1]
                    if astrip in self.dict['Alias']: return self.dict['Alias'][astrip]
            ## Try prefix-stripping alias keys too
            if self.getBool('StripNum'):
                for akey in rje.sortKeys(self.dict['Alias']):
                    if rje.matchExp('^(\d+)_(\S+)',akey) and alias == rje.matchExp('^(\d+)_(\S+)',akey)[1]:
                        return self.dict['Alias'][akey]
            if expect:
                self.debug('%s vs %s' % (alias,rje.sortKeys(self.dict['Alias'])))
                raise ValueError('Genome/Prefix "%s" has no mapping!' % alias)
            return None
        except: self.errorLog('%s.aliasEntry() error' % self.prog()); raise
#########################################################################################################################
    def genomeString(self,gentry):  ## Returns genome and prefix as string                                      # v0.5.1
        '''Returns genome and prefix as string.'''
        genome = gentry['Genome']
        if gentry['Prefix'] != genome: genome += ' (%s)' % gentry['Prefix']
        return genome
#########################################################################################################################
    def genomeField(self,gentry,raw=False): ### Returns the string to be used for genome fields                 # v0.5.1
        '''
        Returns the string to be used for genome fields. Will strip leading numbers if appropriate.
        '''
        if raw: genome = gentry
        else: genome = gentry['Genome']
        if self.getBool('StripNum') and rje.matchExp('^(\d+)_(\S+)',genome):
            genome = rje.matchExp('^(\d+)_(\S+)',genome)[1]
        return genome
#########################################################################################################################
    def genomeData(self,genome,field=None): ### Returns the genomes table entry (or value) for genome           # v0.6.2
        '''
        Returns the genomes table entry (or value) for genome.
        :param genome:str = genome (or group) name
        :param field:str [None] = return the value for that field, or whole entry if None
        :return: value for that field, or whole entry if None
        '''
        gentry = self.aliasEntry(genome)
        if field: return gentry[field]
        return gentry
#########################################################################################################################
    ### <5> ### BUSCOMP Grouping Methods
#########################################################################################################################
    def genomeGroups(self,review=True): ### Set up and/or edit groups membership                                # v0.5.1
        '''
        Set up and/or edit groups membership.
        >> review:bool [True] = Whether to give option to review/edit groups.

        ## Grouping and BUSCO compilation

        BUSCOMP compilation will take a set of genomes and identify the best rating for each BUSCO across the set. By
        default, this is done for all input genomes and assigned to a "BUSCOMP" rating. Ratings are compiled according
        to the `ratings=LIST` hierarchy, which by default is: `Duplicated`, `Complete`, `Fragmented`, `Missing`. The
        "best" rating (first in the hierarchy) will be used for that group, even if it is found in only a single
        assembly in the group.

        ### Loading Groups

        Additional subsets of genomes can also be grouped together for compilation, using the `groups=FILE` input. This
        must have `Genome` and `Group` fields. Once generated, a Group becomes a "Genome". Groups can therefore include
        other Groups in them, provided there is no circular membership. Genome mapping should also recognise Prefix
        values. The "BUSCOMP" group with *all* genomes will be added unless a conflicting "BUSCOMP" group is provided,
        or `buscomp=F`.

        **NOTE:** Genome groups are generated mapping onto the `Genome` (or `Prefix`) values _after_ any alias mapping
        but _before_ any manual editing of `Genome` aliases.

        ### Manual group review/editing

        When running interactively (`i>=0`), Group review and editing can be accessed from the main Genome menu. This
        will show the genomes currently in each group, and provides the following menu options:

        * `<A>dd` = Create a new Group. First, a group name will need to be chosen, then the `<E>dit` menu will be
        triggered.
        * `<D>elete` = Remove the Group. Note that this will not remove Group members from any "parent" groups that
        contained the deleted group when they were loaded.
        * `<R>ename` = Change the Group name.
        * `<E>dit` = Edit group membership. This will give options to add/remove each assembly in the analysis.
        * `E<x>it menu` = Return to main Genome Menu.
        * `<Q>uit` = Quit BUSCOMP.

        ### Saving and sorting

        Once group edits are complete, group data will be saved to the `*.groups.tdt` file. If `runsort=group` then
        genome sorting will also be performed following group edits.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            gdb = self.db('genomes')
            adb = self.db('alias')
            if adb and 'Genome' not in adb.fields(): adb = None
            if adb and 'Description' not in adb.fields(): adb = None
            grpdb = self.db('groups')
            genomes = gdb.index('Genome',force=True)
            grpsave = False

            ### ~ [2] Load/Create Groups table  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Load groups file if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('Groups'): self.setStr({'Groups':'%s.groups.tdt' % self.baseFile()})
            if not grpdb:
                #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ SETUP GENOME GROUPS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
                self.headLog('SETUP GENOME GROUPS',line='=')
                grpdb = db.addTable(self.getStr('Groups'),mainkeys=['Auto'],name='groups',expect=False,ignore=[])
                self.dict['Groups'] = {}
            ## ~ [2b] Check/update/expand groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if grpdb:
                checking = True
                grpx = grpdb.entryNum()
                while checking:
                    groups = grpdb.index('Group',force=True,log=False)
                    checking = False  # Switch to True if changes made in case Groups lost
                    for gentry in grpdb.entries():
                        if gentry['Genome'] in genomes: continue
                        if gentry['Genome'] in groups: continue   #i# Groups handled below
                        aentry = self.aliasEntry(gentry['Genome'],expect=False)
                        if aentry: gentry['Genome'] = aentry['Genome']
                        else:
                            self.warnLog('Alias/Prefix/Group "%s" not found: dropping group entry' % gentry['Genome'])
                            grpdb.dropEntry(gentry)
                    ## Expand groups of groups
                    while rje.listIntersect(grpdb.index('Genome',force=True,log=False).keys(),grpdb.index('Group',force=True,log=False).keys()):
                        # List of groups that are the 'Genome' in another group
                        groupgroups = rje.listIntersect(grpdb.index('Genome').keys(),grpdb.index('Group').keys())
                        for expandgroup in groupgroups:
                            expandgenomes = grpdb.indexDataList('Group',expandgroup,'Genome')
                            parentgroups = grpdb.indexDataList('Genome',expandgroup,'Group')
                            for group in parentgroups:
                                parentgenomes = grpdb.indexDataList('Group',parentgroup,'Genome')
                                for genome in rje.listDifference(expandgenomes,parentgenomes):
                                    if genome == group: self.warnLog('Circular group membership detected: %s assigned as subgroup of self!' % group)
                                    else: grpdb.addEntry({'Genome':genome,'Group':group})
                            grpdb.dropEntriesDirect('Group',[expandgroup])
                            checking = True
                ## Identify whether table needs saving
                grpsave = grpdb.entryNum() != grpx
            ## ~ [2c] Generate group prefixes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['GroupPrefix'] = grprefix = {}
            groups = {}
            if grpdb:
                groups = grpdb.index('Group',force=True,log=False)
                for group in rje.sortKeys(groups):
                    self.printLog('#GRP','%s: %s' % (group,string.join(grpdb.indexDataList('Group',group,'Genome'),'|')))
                    grprefix[group] = string.join(grpdb.indexDataList('Group',group,'Genome'),'|')
                    if adb and group in adb.index('Genome'):
                        desc = adb.indexDataList('Genome',group,'Description')[0]
                        self.dict['Alias'][group] = gdb.addEntry({'Prefix': grprefix[group],'Genome':group,'Description': desc})
                groups = grpdb.index('Group',force=True,log=False)
            ## ~ [2d] Create empty group table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                grpdb = db.addEmptyTable('groups',['#','Genome','Group'],['#'])
            ## ~ [2e] Add Group BUSCOMP if not present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# NOTE: Adding a "BUSCOMP" group will prevent addition of a BUSCOMP group with all genomes in it
            if not 'BUSCOMP' in groups:
                gx = grpdb.entryNum()
                for genome in genomes:
                    grpdb.addEntry({'Genome':genome,'Group':'BUSCOMP'})
                grprefix['BUSCOMP'] = 'AllGenomes'
                desc = 'All assemblies'
                group = 'BUSCOMP'
                self.printLog('#GROUP','Group BUSCOMP added: %d -> %d Genome-Group pairs' % (gx,grpdb.entryNum()))
                self.dict['Alias'][group] = gdb.addEntry({'Prefix': grprefix[group],'Genome':group,'Description': desc})
                #grpsave = True
                #x#gdb.addEntry({'Prefix':'AllGenomes','Genome':'BUSCOMP'})
            ## ~ [2f] Generate Groups dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if grpdb and not self.dict['Groups']:
                for group in grpdb.index('Group'):
                    # Pull out genome dictionaries so that any changes will be mapped across
                    self.dict['Groups'][group] = gdb.indexEntries('Genome',grpdb.indexDataList('Group',group,'Genome'))

            ### ~ [3] Group Edit/Review ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3c] Manual group review ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add option to add/edit groups #!#
            #grpdb.addEntry({'Prefix': string.join(grpdb.indexDataList('Group',group,'Genome'),'|'),'Genome':group})
            if review and self.i() >= 0: # and rje.yesNo('Review/Edit BUSCOMP groups?'):
                while self.i() >= 0:
                    gtext = '\n%d BUSCOMP Groups:' % grpdb.entryNum()
                    grpnum = grpdb.entryNum()
                    groups = rje.sortKeys(grpdb.index('Group'))
                    default = 'X'
                    i = 1
                    for group in groups:
                        #gtext += '\n%s. %2s: %s (%s)' % (i,group,string.join(grpdb.indexDataList('Group',group,'Genome'),', '),gdb.data(group)['Description'])
                        gtext += '\n%s. %2s: %s (%s)' % (i,group,string.join(grpdb.indexDataList('Group',group,'Genome'),', '),self.genomeData(group,'Description'))
                        i += 1
                    gtext += '\n\n<A>dd | <D>elete | <R>ename | <E>dit | E<x>it menu | <Q>uit'
                    ## Loop through choices ##
                    choice = rje.choice(gtext,default=default).upper()
                    if choice == 'D':
                        while i >= grpnum:
                            i = rje.getInt('Group to Delete? (0 to cancel)',default=0)
                        if i > 0 and rje.yesNo('Delete group "%s"' % groups[i-1]):
                            grpdb.dropEntriesDirect('Group',[groups[i-1]])
                            self.dict['Groups'].pop(groups[i-1])
                    elif choice == 'R':
                        while i >= grpnum:
                            i = rje.getInt('Group to Rename? (0 to cancel)',default=0)
                        if i > 0:
                            newgroup = rje.choice('New Group Name?:',default=groups[i-1],confirm=True)
                            while newgroup != groups[i-1] and newgroup in groups:
                                newgroup = rje.choice('Group name in use! New Group Name?:',default=groups[i-1],confirm=True,whitespace=False)
                            if newgroup != groups[i-1]:
                                for entry in grpdb.indexEntries('Group',groups[i-1]):
                                    entry['Group'] = newgroup
                                grpdb.index('Group',force=True,log=False)
                                self.dict['Groups'][newgroup] = self.dict['Groups'].pop(groups[i-1])
                                self.dict['Alias'][newgroup] = self.dict['Alias'].pop(groups[i-1])
                    elif choice == 'E':
                        while i > len(groups):
                            i = rje.getInt('Group to Edit? (0 to cancel)',default=0)
                        if i > 0:
                            self.editGroup(groups[i-1])
                    elif choice == 'A':
                        newgroup = rje.choice('Group Name?:',default='',confirm=True)
                        while newgroup in groups or newgroup in self.dict['Alias']:
                            newgroup = rje.choice('Group name in use! New Group Name?:',default='',confirm=True,whitespace=False)
                        if newgroup:
                            self.dict['Groups'][newgroup] = []
                            self.dict['Alias'][newgroup] = gdb.addEntry({'Genome':newgroup,'Prefix':'','Description':''})
                            self.editGroup(newgroup)
                    elif choice == 'X' and rje.yesNo('End review of groups?'): break
                    elif choice == 'Q' and rje.yesNo('Quit BUSCOMP?'): os.sys.exit(0)
                grpsave = grpdb.entryNum()

            ### ~ [4] Finish group processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Check Groups for errors and missing Genomes etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            grpsave = self.checkGroupGenomes() or grpsave
            ## ~ [4b] Save groups to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if grpsave: grpdb.saveToFile(savefields=['Genome','Group'])
            ## ~ [4c] Sort genomes if using Group sorting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if grpsave: self.sortGenomes(quicksort=True,restrict=['Group'])

            return True     # Setup successful
        except: self.errorLog('Problem during %s genomeGroups.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def checkGroupGenomes(self): ### Checks group membership and presence in genomes table                      # v0.5.1
        '''Checks group membership and presence in genomes table. Returns whether any groups were removed.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            grpdb = self.db('groups')
            gdb = self.db('genomes')
            adb = self.db('alias')
            grprefix = self.dict['GroupPrefix']
            #!# This method needs to be added and consolidate some other methods elsewhere.
            #i# Groups should all be sorted out following the loading of genomes, e.g. at the start.
            genomes  = gdb.index('Genome',force=True,log=False)
            groups = rje.sortKeys(grpdb.index('Group',force=True,log=False))

            ### ~ [2] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Check Groups are in the Genome Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for group in groups:
                if group not in grprefix:
                    grprefix[group] = string.join(grpdb.indexDataList('Group',group,'Genome'),'|')
                    self.deBug('Had to add grpprefix for %s' % group)
                if group not in genomes:
                    if group in self.dict['Alias']: raise ValueError('Cannot use "%s" for group name: already in use for group or genome!' % group)
                    desc = ''
                    if adb and 'Description' in adb.fields() and group in adb.index('Genome'):
                        desc = adb.indexEntries('Genome',group)[0]['Description']
                        self.printLog('#DESC','Description updated: %s = "%s"' % (group,desc))
                    if not desc: desc = grprefix[group]
                    gentry = gdb.addEntry({'Prefix': grprefix[group],'Genome':group,'Description': desc})
                    self.dict['Alias'][group] = gentry
                    self.printLog('#GROUP','Added group %s to Genome table (%s)' % (gentry['Genome'],gentry['Prefix']))
            ## ~ [2b] Cleanup empty groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dx = 0
            for entry in grpdb.entries():
                if entry['Genome'] in genomes: continue
                grpdb.dropEntry(entry)
                dx += 1
            if dx:
                self.warnLog('%s Group entries dropped as genome not in Genomes table' % rje.iStr(dx))
                groups = rje.sortKeys(grpdb.index('Group',force=True,log=False))
            ## ~ [2c] Cleanup missing groups from Genome Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gentry in gdb.entries():
                if gentry['Directory']: continue
                if gentry['Fasta']: continue
                if gentry['Genome'] in groups: continue
                self.warnLog('Genome entry "%s" not found in directories, fasta list or groups: dropping' % gentry['Genome'])
                gdb.dropEntry(gentry)
            #i# Returns whether any groups were removed
            ## ~ [2d] Cleanup self.dict['Groups'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for group in rje.sortKeys(self.dict['Groups']):
                for gentry in self.dict['Groups'][group][0:]:
                    if gentry not in gdb.entries():
                        self.dict['Groups'][group].remove(gentry)
                        self.devLog('Dropped %s %s from dict' % (group,gentry['Genome']))
                if not self.dict['Groups'][group]:
                    self.dict['Groups'].pop(group)
                    self.devLog('Dropped empty %s from dict' % (group))

            return dx
        except: self.errorLog('Problem during %s checkGroupGenomes.' % self.prog()); raise  # Setup failed
#########################################################################################################################
    def editGroup(self,group): ### Edits group membership                                                       # v0.5.1
        '''Edits group membership.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            grpdb = self.db('groups')
            gdb = self.db('genomes')
            grpentry = self.aliasEntry(group)
            grprefix = self.dict['GroupPrefix']
            updated = False
            while group:
                genomes = grpdb.indexDataList('Group',group,'Genome')
                gtext = '\nGroup name: %s\n' % group
                i = 1
                for gentry in gdb.entries(sorted=True):
                    gtext += '|-- %2s' % i
                    if gentry['Genome'] in genomes:
                        gtext += ' [+] '
                    elif gentry['Genome'] in self.dict['Groups']: gtext += ' [X] '
                    else: gtext += ' [ ] '
                    gtext += '%s (%s)\n' % (gentry['Genome'],gentry['Prefix'])
                    i += 1
                gtext += '\nEnter genome number to toggle membership; 0 to quit'
                choice = rje.getInt(gtext,default=0)
                if choice == 0:
                    if group not in grprefix:
                        grprefix[group] = string.join(grpdb.indexDataList('Group',group,'Genome'),'|')
                    if not grpentry['Description']: grpentry['Description'] = grprefix[group]
                    grpentry['Description'] = rje.choice('Group description?', default=grpentry['Description'], confirm=True)
                    break
                if choice <= len(gdb.entries(sorted=True)):
                    updated = True
                    genome = gdb.entries(sorted=True)[choice-1]['Genome']
                    if genome in self.dict['Groups']: continue
                    gentry = self.aliasEntry(genome)
                    if gentry in self.dict['Groups'][group]: self.dict['Groups'][group].remove(gentry)
                    else: self.dict['Groups'][group].append(gentry)
                    if genome in genomes:
                        genomes.remove(genome)
                        for entry in grpdb.indexEntries('Group',group):
                            if entry['Genome'] == genome: grpdb.dropEntry(entry)
                    else:
                        grpdb.addEntry({'Group':group,'Genome':genome})
                    grpdb.index('Group',force=True,log=False)
            if updated:
                #grprefix[group] = rje.choice('Group description',default=string.join(grpdb.indexDataList('Group',group,'Genome'),'|'),confirm=True)
                grprefix[group] = string.join(grpdb.indexDataList('Group',group,'Genome'),'|')
            if group in gdb.index('Genome',force=True,log=False):
                #if updated: gdb.data(gdb.index('Genome')[group][0])['Prefix'] = grprefix[group]
                if updated: self.genomeData(group)['Prefix'] = grprefix[group]
            else: gdb.addEntry({'Genome':group,'Prefix':grprefix[group]})
        except: self.errorLog('Problem during %s editGroup.' % self.prog())
        try: return updated
        except: return False
#########################################################################################################################
    ### <6> ### BUSCOMP Compilation Methods                                                                             #
#########################################################################################################################
    def bugShush(self): return (self.debugging() and self.quiet()) or self.silence()
#########################################################################################################################
    def compileBUSCO(self): ### Generic Load and compiles full BUSCO tables                                     # v0.5.1
        '''
        ## BUSCO Compilation

        Once all the genomes are loaded, BUSCOMP compiles the BUSCO results. This is done for three purposes:

        1. Report comparative BUSCO statistics across assemblies and plot with other assembly statistics.
        2. Combine BUSCO results across groups of assemblies to summarise total BUSCO coverage.
        3. Identify and compile the best `Complete` BUSCO hits for BUSCOMPSeq analysis (below).

        Loaded assemblies with an identified full results `Table` will be processed and added to the `*.full.tdt` table.
        The best rating per BUSCO gene will then be used for summary BUSCO counts for each genome. Finally, each Group
        will have its combined BUSCO rating calculated. This is performed by choosing the "best" rating for each BUSCO
        across the Group's members, where "best" is determined by the `dupbest=T/F` setting:

        - By default (`dupbest=F`), the rating hierarchy is: 'Complete', 'Duplicated', 'Fragmented', 'Missing'.
        - If `dupbest=T`, the rating hierarchy is: 'Duplicated', 'Complete', 'Fragmented', 'Missing'.

        Total BUSCO counts (`N`) and summed Ratings are added to the main `*.genomes.tdt` table for each assembly and
        group. A `*.busco.tdt` file is also generated that has the rating for each BUSCO gene (rows) in each
        assembly/group (columns).

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ COMPILE BUSCO RESULTS ~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            self.headLog('COMPILE BUSCO RESULTS',line='=')
            db = self.db()
            gdb = self.db('genomes')
            for field in ['N']+self.list['Ratings']: gdb.addField(field,evalue=0)
            fullhead = ['BuscoID','Status','Contig','Start','End','Score','Length']
            #?# Add option to upload and append/update existing table? Can just re-run for now - it's fast!
            fulldb = None
            buscodb = None

            ### ~ [2] Cycle through tables and load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gentry in gdb.entries(sorted=True):
                if not gentry['Table']: continue
                logtext = '%s [%s]' % (gentry['Genome'],gentry['Table'])
                self.progLog('\r#BUSCO','Processing %s...' % logtext)

                # - Add to Table:full (+Genome field)
                self.bugShush()
                fdb = db.addTable(gentry['Table'],mainkeys='auto',headers=fullhead,expect=True,name=gentry['Genome'])
                logtext += ' (%s entries)' % fdb.entryNum()
                self.talk(); self.progLog('\r#BUSCO','Processing %s...' % logtext); self.bugShush()
                fdb.fillBlanks(blank=0,fields=['Score','Length'],fillempty=True,prog=True,log=True)
                fdb.dataFormat({'Score':'num','Length':'int'})
                fdb.addField('Genome',evalue=self.genomeField(gentry))
                fdb.newKey(['Genome','#'],startfields=True,strict=True)
                if not fulldb: fulldb = db.copyTable(fdb,'full')
                else:
                    for entry in fdb.entries(): fulldb.addEntry(entry)    #i# This will break dictionary connection

                # - Rename rating with Genome
                fdb.renameField('Status',self.genomeField(gentry))
                # - Compress to best rating per EOG
                fdb.rankFieldByIndex('BuscoID','Score',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False,warn=True)
                fdb.dropEntriesDirect('Rank',[1],inverse=True)
                fdb.dropField('Rank',log=self.debugging())
                fdb.rankFieldByIndex('BuscoID','Length',newfield='Rank',rev=True,absolute=True,lowest=True,unique=True,warn=True)
                fdb.dropEntriesDirect('Rank',[1],inverse=True)
                fdb.dropField('Rank',log=self.debugging())
                fdb.newKey('BuscoID')
                fdb.keepFields(['Genome','BuscoID',self.genomeField(gentry)],log=self.debugging())
                logtext = logtext[:-1] + '-> %s BUSCOs)' % fdb.entryNum()
                self.talk(); self.progLog('\r#BUSCO','Processing %s' % logtext); self.bugShush()

                # - Add to Table:busco
                if not buscodb:
                    buscodb = db.copyTable(fdb,'busco')
                    buscodb.dropField('Genome',log=self.debugging())
                else:
                    buscodb.addField(self.genomeField(gentry))
                    for eog in buscodb.dataKeys():
                        if fdb.data(eog):
                            buscodb.data(eog)[self.genomeField(gentry)] = fdb.data(eog)[self.genomeField(gentry)]
                        else:
                            self.warnLog('Cannot find %s results for %s [%s]' % (eog,self.genomeString(gentry),gentry['Table']),quitchoice=self.debugging())

                # - Add N & Ratings to Table:genomes ['N']+self.list['Ratings']
                fdb.addField('N',evalue=1)
                for field in self.list['Ratings']:
                    fdb.addField(field,evalue=0,log=self.debugging())
                for entry in fdb.entries(): entry[entry[self.genomeField(gentry)]] = 1
                fdb.compress(['Genome'],default='sum')
                gdict = fdb.data(self.genomeField(gentry))
                db.deleteTable(fdb)
                gdict.pop('Genome'); gdict.pop('BuscoID')
                gentry.update(gdict)
                self.talk(); self.printLog('\r#BUSCO','Processed %s.' % logtext)
            if not fulldb:
                self.printLog('#BUSCO','No BUSCO results to process.')
                return True

            ### ~ [3] Process Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ranklist = self.list['Ratings']
            if not self.getBool('DupBest'): ranklist =  ['Complete', 'Duplicated', 'Fragmented', 'Missing']
            groups = rje.sortKeys(self.dict['Groups'])
            for group in groups:
                # Identify genomes in group
                gentries = []
                for genentry in self.dict['Groups'][group]:
                    if not genentry['Directory']:
                        self.printLog('#SKIP','Skipping Genome %s: no BUSCO data' % self.genomeString(genentry))
                    else: gentries.append(genentry)
                # Identify Group in Genome table
                try:
                    grpentry = self.aliasEntry(group,expect=False)
                    if not grpentry: raise ValueError('Group %s missing from Genomes table!' % group)
                except:
                    self.errorLog('Compile BUSCO Groups error',quitchoice=True)
                    continue
                # Add group to BUSCO table
                buscodb.addField(group,evalue=self.list['Ratings'][-1])
                for eog in buscodb.dataKeys():
                    for genentry in gentries:
                        gfield = self.genomeField(genentry)
                        if not gfield: raise ValueError('No GenomeField entry for: %s' % genentry)
                        try:
                            ranklist.index(buscodb.data(eog)[gfield])
                            ranklist.index(buscodb.data(eog)[group])
                        except:
                            self.errorLog('Problem with %s data - "%s" or "%s" classification not recognised: %s' % (eog, gfield, group, buscodb.data(eog)))
                        if ranklist.index(buscodb.data(eog)[gfield]) < ranklist.index(buscodb.data(eog)[group]):
                            buscodb.data(eog)[group] = buscodb.data(eog)[gfield]
                    grpentry['N'] += 1
                    grpentry[buscodb.data(eog)[group]] += 1
                self.printLog('#GRP','Compiled BUSCO stats for group %s (%d of %d genomes)' % (group,len(gentries),len(self.dict['Groups'][group])))

            return True

            # ### ~ [3] Process Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # ranklist = self.list['Ratings']
            # groups = rje.sortKeys(grpdb.index('Group'))
            # for group in groups:
            #     # Identify genomes in group
            #     genomes = grpdb.indexDataList('Group',group,'Genome')
            #     for genome in genomes[0:]:
            #         genentry = self.aliasEntry(genome)
            #         if not genentry['Directory']:
            #             self.printLog('#SKIP','Skipping Genome %s: no BUSCO data' % genome)
            #             genomes.remove(genome)
            #     # Identify Group in Genome table
            #     if group in gdb.index('Genome'):
            #         gentry = gdb.indexEntries('Genome',group)[0]
            #     else:
            #         try: raise ValueError('Group %s missing from Genomes table!' % group)
            #         except: self.errorLog('Compile BUSCO Groups error',quitchoice=True)
            #         continue
            #     # Add group to BUSCO table
            #     buscodb.addField(group,evalue=self.list['Ratings'][-1])
            #
            #     #!# NOTE: No longer allowing groups to be part of groups here - should be resolved earlier.
            #     #dontgoback = []
            #     #while rje.listIntersect(genomes,groups):
            #     #    for genome in genomes[0:]:
            #     #        if genome in dontgoback: raise ValueError('Circular group membership: %s and %s' % (group,genome))
            #     #        if genome in groups:
            #     #            genomes.pop(genome)
            #     #            genomes += grpdb.indexDataList('Group',genome,'Genome')
            #     #            dontgoback.append(genome)
            #
            #     for eog in buscodb.dataKeys():
            #         for genome in genomes:
            #             genentry = self.aliasEntry(genome)
            #             gfield = self.genomeField( genentry )
            #             if ranklist.index(buscodb.data(eog)[gfield]) < ranklist.index(buscodb.data(eog)[group]):
            #                 buscodb.data(eog)[group] = buscodb.data(eog)[gfield]
            #         gentry['N'] += 1
            #         gentry[buscodb.data(eog)[group]] += 1
            #     self.printLog('#GRP','Compiled BUSCO stats for group %s (%d genomes)' % (group,len(genomes)))
            #
            # return True
        except:
            self.talk()
            self.errorLog('%s.compileBUSCO error' % self.prog()); raise
#########################################################################################################################
    def rawSummary(self,log=True,full=False,screen=False):   ### Summarise stats from self.db('genomes')        # v0.5.1
        '''
        Summarise stats from self.db('genomes').
        >> log:bool [True] = Whether to output summary to log
        >> full:bool [False] = Whether to output individual rating numbers, like BUSCO output
        << str = Summary text as a string that could be wrapped in R code etc.

        BUSCO example -
        INFO    Results:
        INFO    C:17.5%[S:17.5%,D:0.0%],F:7.3%,M:75.2%,n:3950
        INFO    692 Complete BUSCOs (C)
        INFO    691 Complete and single-copy BUSCOs (S)
        INFO    1 Complete and duplicated BUSCOs (D)
        INFO    290 Fragmented BUSCOs (F)
        INFO    2968 Missing BUSCOs (M)
        INFO    3950 Total BUSCO groups searched
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            gdb = self.db('genomes')
            #i# Check data formatting. To include genomes w/o BUSCO analysis (zero values), switch skipblank=True
            gdb.dataFormat({'Complete':'int','Duplicated':'int','Fragmented':'int','Missing':'int','N':'int'},skipblank=True)
            txt = ''
            ### ~ [2] Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log:
                #if full: self.printLog('#==#','# ================ BUSCO RESULTS ================ #')
                #else: self.printLog('#==#','# ================ BUSCO SUMMARY ================ #')
                if full: self.headLog('BUSCO RESULTS',line='=')
                else: self.headLog('BUSCO SUMMARY',line='=')
            for gentry in gdb.entries(sorted=True):
                genome = self.genomeField(gentry)
                if gentry['N'] == '': self.printLog('#BUSCO','Skipping %s (no BUSCO results)' % genome); continue
                txt += '%s BUSCO Results:\n' % (genome)
                n = gentry['N'] / 100.0
                if n:
                    busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%],F:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete']+gentry['Duplicated'])/n,gentry['Complete']/n,gentry['Duplicated']/n,gentry['Fragmented']/n,gentry['Missing']/n,gentry['N'])
                else: busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%],F:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete']+gentry['Duplicated']),gentry['Complete'],gentry['Duplicated'],gentry['Fragmented'],gentry['Missing'],gentry['N'])
                txt += '        %s\n' % busco
                if log: self.printLog('#BUSCO','Results: %s - %s' % (busco,genome))
                if not full: txt += '\n'; continue
                for bustxt in ['%d Complete BUSCOs (C)' % (gentry['Complete']+gentry['Duplicated']),
                    '%d Complete and single-copy BUSCOs (S)' % gentry['Complete'],
                    '%d Complete and duplicated BUSCOs (D)' % +gentry['Duplicated'],
                    '%d Fragmented BUSCOs (F)' % gentry['Fragmented'],
                    '%d Missing BUSCOs (M)' % gentry['Missing'],
                    '%d Total BUSCO groups searched' % gentry['N']]:
                    txt += '        %s\n' % bustxt
                    if log: self.printLog('#INFO',bustxt)
                txt += '\n'

            if screen and not log: self.verbose(txt,v=0)
            return txt
        except: self.errorLog('Problem during %s sortGenomes.' % self.prog()); raise
#########################################################################################################################
    def compileBUSCOSeqs(self): ### Compile BUSCO sequences into *.buscomp.fasta per EOG by Score then Length.  # v0.5.2
        '''
        ### BUSCOMPSeq sequence compilation

        The final step of the BUSCOMP BUSCO compilation is to extract the best Complete BUSCO sequences from those
        available. For all assemblies with BUSCO results and a `single_copy_busco_sequences/`, BUSCOMP ranks all the hits
        across all assemblies by `Score` and keeps the top-ranking hits. Ties are then resolved by `Length`, keeping the
        longest sequence. Ties for `Score` and `Length` will keep an arbitrary entry as the winner. `Single` complete
        hits are given preference over `Duplicate` hits, even if they have a lower score, because only `Single` hit
        have their sequences saved by BUSCO in the `single_copy_busco_sequences/` directory.

        If `buscofas=FASFILE` has been given, this will be used for BUSCOMPSeq searches. Otherwise, the best BUSCO
        seqences identified for each BUSCO gene will be saved as `*.buscomp.fasta` for BUSCOMPSeq analysis. The
        exception is if `buscofas=FASFILE` is pointing to the `*.buscomp.fasta` and `force=T`, or if the
        `buscofas=FILE` fasta file cannot be found.

        **NOTE:** The `buscofas=FASFILE` option has not been tested and might give some unexpected behaviours, as some
        of the quoted figures will still be based on the calculated BUSCOMPSeq data.

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#==#','# ================ COMPILE BUSCO SEQUENCES ================ #')
            self.headLog('COMPILE BUSCO SEQUENCES',line='=')
            gdb = self.db('genomes')
            fdb = self.db().copyTable(self.db('full'),'buscoseq')
            genomes = fdb.index('Genome').keys()
            badgenome = []
            for entry in gdb.entries():
                if entry['Genome'] in genomes and not entry['Sequences']: badgenome.append(entry['Genome'])
            if badgenome:
                self.printLog('#DROP','Dropping %d genomes without sequence data for BUSCOMPSeq compilation' % len(badgenome))
                fdb.dropEntriesDirect('Genome',badgenome,log=True)

            ### ~ [2] Rank and extract ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # - Add Rank to Table:Full
            # - rank EOG by score (lowest)
            #i# Need to add clever filtering of Duplicates to keep the best Single score #!#
            self.progLog('#SEQ','Compiling best BUSCOs ...')
            for fentry in fdb.indexEntries('Status',['Duplicated','Fragmented']):
                fentry['Score'] /= 1000.0
            fdb.rankFieldByIndex('BuscoID','Score',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False,warn=True)
            fdb.dropEntriesDirect('Rank',[1],inverse=True)
            for fentry in fdb.indexEntries('Status',['Duplicated','Fragmented']):
                fentry['Score'] *= 1000.0
            fdb.dropField('Rank')
            # - extract rank 1
            # - rank by length (unique)
            fdb.rankFieldByIndex('BuscoID','Length',newfield='Rank',rev=True,absolute=True,lowest=True,unique=True,warn=True)
            fdb.dropEntriesDirect('Rank',[1],inverse=True)
            fdb.dropField('Rank')

            ### ~ [3] Extract to fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            buscofas = '%s.buscomp.fasta' % self.baseFile()
            buscofaa = '%s.buscomp.faa' % self.baseFile()
            if rje.exists(self.getStr('BUSCOFas')) and (buscofas != self.getStr('BUSCOFas') or not self.force()):
                buscofas = self.getStr('BUSCOFas')
                self.printLog('#FASTA','BUSCO fasta file given: %s. No BUSCOMPSeq fasta generation' % buscofas)
                return True
            elif self.getStrLC('BUSCOFas') and not self.force():
                try: raise IOError('Cannot find %s: will create and use %s' % (self.getStr('BUSCOFas'),buscofas))
                except: self.errorLog('BUSCOFas error')
                self.setStr({'BUSCOFas':buscofas})
            else: self.setStr({'BUSCOFas':buscofas})
            self.progLog('#SEQ','Generating %s and %s ...' % (buscofas,buscofaa))
            #?# - Read sequences into obj['BUSCOMP']
            rje.backup(self,buscofas)
            rje.backup(self,buscofaa)
            FASTA = open(buscofas,'w')
            FAA = open(buscofaa,'w')
            fdb.newKey(['BuscoID'])
            fdb.dropField('#')
            ex = 0; faax = 0
            for entry in fdb.entries(sorted=True):
                if entry['Status'] == 'Missing': entry['Genome'] = 'Missing'; continue
                self.progLog('#SEQ','Generating %s and %s ...' % (buscofas,buscofaa))
                eog = entry['BuscoID']
                genome = entry['Genome']
                gentry = self.aliasEntry(genome)
                prefix = gentry['Prefix']
                seqdir = rje.makePath(gentry['Directory']) + 'single_copy_busco_sequences'
                prefseqdir = rje.makePath(gentry['Directory']) + 'single_copy_busco_sequences_' + prefix
                v4seqdir = rje.makePath(gentry['Directory']) + 'busco_sequences/single_copy_busco_sequences'
                if rje.exists(prefseqdir): seqdir = prefseqdir
                elif rje.exists(v4seqdir): seqdir = v4seqdir
                fna = seqdir + '/%s.fna' % eog
                if rje.exists(fna):
                    #!# Fix this! 'Qry': 'EOG0907007V:tiger.wtdbg2v1.racon.fasta:ctg_NOTSC__NSCUV1ONT0006:601178-628876'
                    fnalines = open(fna,'r').readlines()
                    if string.split(fnalines[0],':',maxsplit=1)[0][1:] == eog:  # BUSCO v3
                        fnalines[0] = string.join(string.split(fnalines[0],':',maxsplit=1))
                    else:
                        fnalines[0] = '>{0} {1}'.format(eog,fnalines[0]) # BUSCO v4
                    for i in range(1,len(fnalines)):
                        if fnalines[i].startswith('>'):
                            self.warnLog('"Single copy" BUSCO %s has 2+ sequences in %s! (Keeping first.)' % (eog,fna))
                            fnalines = fnalines[:i]
                            break
                    FASTA.writelines(fnalines)
                    ex += 1
                else:
                    if entry['Status'] == 'Single':
                        self.warnLog('%s not found (rating=Single)' % fna)
                faa = seqdir + '/%s.faa' % eog
                if rje.exists(faa):
                    faalines = open(faa,'r').readlines()
                    if string.split(faalines[0],':',maxsplit=1)[0][1:] == eog:  # BUSCO v3
                        faalines[0] = string.join(string.split(faalines[0],':',maxsplit=1))
                    else:
                        faalines[0] = '>{0} {1}'.format(eog,faalines[0]) # BUSCO v4
                    for i in range(1,len(faalines)):
                        if faalines[i].startswith('>'):
                            self.warnLog('"Single copy" BUSCO %s has 2+ sequences in %s! (Keeping first.)' % (eog,faa))
                            faalines = faalines[:i]
                            break
                    FAA.writelines(faalines)
                    faax += 1
            FASTA.close()
            FAA.close()
            self.printLog('#SEQ','DNA sequences for %s best (of %s total) BUSCO hits output to %s' % (rje.iStr(ex),rje.iStr(fdb.entryNum()),self.getStr('BUSCOFas')))
            self.printLog('#SEQ','Protein sequences for %s best (of %s total) BUSCO hits output to %s' % (rje.iStr(faax),rje.iStr(fdb.entryNum()),buscofaa))

            return True     # Setup successful
        except: self.errorLog('Problem during %s compileBUSCOSeqs.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <7> ### BUSCOMPSeq Rating Methods                                                                               #
#########################################################################################################################
    def pafBase(self,genome=None):  ### Returns the base for PAF TDT files, incorporating cutoffs
        '''Returns the base for PAF TDT files, incorporating cutoffs.'''
        outdir = rje.makePath('%s.minimap/' % self.baseFile())
        minlocid = self.getInt('MinLocID')
        minloclen = self.getInt('MinLocLen')
        pafbase = 'N%dL%dID%d' % (self.getInt('MMSecNum'),minloclen,minlocid)
        if self.getBool('UniqueHit'): pafbase += 'U'
        if genome:
            pafbase = '%s%s.%s' % (outdir,genome,pafbase)
        return pafbase
#########################################################################################################################
    def minimapSearches(self):    ### Perform PAF Minimap2 searches.                                            # v0.5.2
        '''
        ## BUSCOMPSeq Minimap2 searches

        Once the gene sequences for BUSCOMPSeq have been established, the next step is to search for them in the assembly
        fasta files using a more deterministic (and faster) approach. For this, minimap2 has been chosen. BUSCOMP will
        perform a minimap2 search on each assembly fasta file. (Use `minimap2=FILE` to tell BUSCOMP where to find
        minimap2 if the path is not an environment variable.)

        Minimap2 will run with `-x splice -uf -C5` and hits will be filtered by length and percentage identity, keeping
        on those that meet the `minlocid=X` (default 0%) and `minloclen=X` (default 20 bp) cutoffs. Hits are also
        reduced to unique coverage by starting with the local alignment with the largest number of identical positions,
        and then trimming any other local alignments that overlap with the same part(s) of the query (BUSCO) sequence.
        This is repeated until all local alignments are non-overlapping (and still meet the identity and length
        criteria). Local hits are then compiled into global alignment (GABLAM) statistics of coverage and identity.

        Minimap2 search results are saved in a `*.minimap/` directory, with the prefix `$BASEFILE.$CUTOFFS`, where
        `$BASEFILE` is the fasta file basename for assembly, and `$CUTOFFS` takes the form `LXXIDXX`, where `LXX` is the
        length cutoff (`minloclen=X`) and `IDXX` is the identity cutoff (`minlocid=X`) - `L20ID80` by default.
        Existing files will be loaded and reused unless `force=T`.

        **NOTE:** Re-using results does NOT robustly check whether the BUSCOMPSeq data has changed, so this directory
        should be deleted if re-running with different sequences. BUSCOMP will save the name of the BUSCO file along
        with the md5sum of its contents to `*.input.md5`, which will be checked if present and re-using data. Some basic
        checks should also be performed during the following results compilation stage, but please note that BUSCO IDs
        are used for sequence identifiers and so changes in BUSCO hit sequences will not be identified if files have
        been inappropriatley copied/moved between runs etc. - please look out for unexpected behaviour outside a "clean"
        run.

        **NOTE:** Minimap2 only works with high sequence identity. BUSCOMP is designed to be used on multiple assemblies
        _from the same species_. If using cross-species BUSCO analysis, odd results might be obtained, biased by the
        evolutionary distance from the species with the most BUSCOMP sequences. Under these circumstances, it might be
        better to swap out minimap2 for BLAST+. This can be achieved by independently running GABLAM using the
        BUSCOMPSeq fasta file as queries and each assembly as the search database, then replacing the `*.gablam.tdt` and
        `*.hitsum.tdt` files before re-running BUSCOMP. Please contact the author if you want help with this.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('MINIMAP2 SEARCHES',line='=')
            outdir = rje.makePath('%s.minimap/' % self.baseFile())
            rje.mkDir(self,outdir)
            gdb = self.db('genomes')
            pafdefault = ['uniqueout=T','localaln=F','mockblast=F','uniquehit=T','alnseq=F','endextend=0',
                          'minlocid=%d' % self.getInt('MinLocID'),'minloclen=%d' % self.getInt('MinLocLen')]
            pafopt = {'p':self.getNum('MMPCut'),'N':self.getInt('MMSecNum')}
            buscofas = '%s.buscomp.fasta' % self.baseFile()
            if rje.exists(self.getStr('BUSCOFas')):
                buscofas = self.getStr('BUSCOFas')
            pafcmd = ['seqin=%s' % buscofas,'mapsplice=T','pafin=minimap']
            ## ~ [1a] md5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fasmd5 = ''
            if rje.exists(buscofas):
                fasmd5 = rje.file2md5(buscofas)
            else: self.warnLog('BUSCOMPSeq fasta file "%s" found - cannot check md5!' % buscofas)

            ### ~ [2] Perform Minimap searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ranx = 0; failx = 0; skipx = 0
            for entry in gdb.entries(sorted=True):
                if not entry['Fasta']: continue
                self.headLog('%s Minimap' % self.genomeField(entry),line='-')
                #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~~~~ %s Minimap ~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % entry['Genome'])
                # Using Fasta file basename not Genome for repeat analysis - allows copying and re-use with & w/o BUSCO
                prefix = rje.baseFile(entry['Fasta'],strip_path=True)
                basefile = self.pafBase(prefix) #'%s%s' % (outdir,entry['Genome'])
                purebase = '%s%s.N%d' % (outdir,prefix,self.getInt('MMSecNum'))
                #if rje.exists('%s%s.paf' % (outdir,entry['Prefix'])): os.rename('%s%s.paf' % (outdir,entry['Prefix']),'%s.paf' % purebase)

                paf = '%s.paf' % purebase
                hitfile = '%s.hitsum.tdt' % basefile
                gabfile = '%s.gablam.tdt' % basefile
                #!# Add md5sum check of files #!#
                md5file = '%s.input.md5' % purebase
                md5fas = md5hash = ''
                md5warn = []
                if rje.exists(md5file):
                    try:
                        [md5fas,md5hash] = string.split(open(md5file,'r').readline())[:2]
                        if md5fas != rje.stripPath(buscofas): md5warn.append('%s fasta filename mismatch' % md5file)
                        if md5hash != fasmd5: md5warn.append('%s md5 mismatch' % md5file)
                        if self.i() >= 0 and md5warn:
                            for warning in md5warn: self.warnLog(warning)
                            if rje.yesNo('Abort BUSCOMPSeq analysis?'):
                                raise ValueError('%s mismatch detected' % md5file)
                    except:
                        md5warn.append('Problem with %s format: cannot check md5 for input' % md5file)
                else: md5warn.append('Cannot find %s: cannot check md5 for input' % md5file)
                if rje.exists(hitfile) and rje.exists(gabfile) and not self.force():
                    self.printLog('#PAF','Using existing parsed %s results (force=F)' % paf)
                    skipx += 1
                    continue
                if not rje.exists(entry['Fasta']):
                    if entry['Fasta']: self.warnLog('Genome fasta file for %s not found: %s' % (entry['Genome'],entry['Fasta']))
                    failx += 1
                    continue
                runcmd = ['basefile=%s' % basefile, 'reference=%s' % entry['Fasta'], 'seqmode=file']
                if rje.exists(paf) and not self.force():
                    runcmd.append('pafin=%s' % paf)
                    self.printLog('#PAF','Using existing %s (force=F)' % paf)
                elif not rje.exists(buscofas):
                    raise IOError('Neither existing minimap results nor BUSCOMPSeq fasta file "%s" found!' % buscofas)
                else:
                    open(md5file,'w').write('%s %s\n' % (rje.stripPath(buscofas),fasmd5))
                pafobj = rje_paf.PAF(self.log,pafdefault+self.cmd_list+pafcmd+runcmd)
                pafobj.dict['MapOpt'] = rje.combineDict(pafobj.dict['MapOpt'],pafopt,overwrite=True)
                if pafobj.run(): ranx += 1
                else: failx += 1
                if not rje.exists(paf) and rje.exists('%s.paf' % basefile):
                    os.rename('%s.paf' % basefile,paf)
                    self.printLog('#PAF','Renamed %s.paf -> %s' % (basefile,paf))
                if not rje.exists('%s.cmd' % paf) and rje.exists('%s.paf.cmd' % basefile):
                    os.rename('%s.paf.cmd' % basefile,'%s.cmd' % paf)
                    self.printLog('#PAF','Renamed %s.paf.cmd -> %s.cmd' % (basefile,paf))
            self.printLog('#PAF','RJE_PAF Minimap2 runs complete for %d genomes; %d failed; %d skipped.' % (ranx,failx,skipx))

            return True     # Setup successful
        except SystemExit: raise    # Child
        except: self.errorLog('Problem during %s minimapSearches.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def md5warnings(self,md5warn,paf):  ### Raises md5 warning if re-using results
        '''
        Raises md5 warning if re-using results.
        '''
        for warning in md5warn: self.warnLog(warning)
        if md5warn: self.warnLog('Existing %s results may not be consistent with current sequence set.' % paf)
#########################################################################################################################
    def compileMinimap(self):    ### Compile PAF Minimap2 searches.                                             # v0.5.1
        '''
        ### BUSCOMPSeq Minimap2 search compilation.

        Minimap2 searches of the compiled BUSCOMP sequences are then compiled in similar vein to the original BUSCO
        results. The primary difference is that the search results need to first be converted into BUSCO-style ratings.
        BUSCOMP uses results from both the `*.hitsum.tdt` table (overall coverage in assembly) and the `*.gablam.tdt`
        table (coverage in a single contig/scaffold) to rate sequences as:

        * Complete: 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
        * Duplicated: 95%+ Coverage in 2+ contigs/scaffolds.
        * Fragmented: 95%+ combined coverage but not in any single contig/scaffold.
        * Partial: 40-95% combined coverage.
        * Ghost: Hits meeting local cutoff but <40% combined coverage.
        * Missing: No hits meeting local cutoff.

        When compiling the results for all BUSCO genes, Single/Duplicated Complete hits will also be rated as
        **Identical** if they have 100% coverage and 100% identity in at least one contig/scaffold.

        Once all the individual assemblies have been rated for the full set of assemblies, results are compiled across
        Groups as described for the original BUSCO results (above). Because no genes receive an individual "Identical"
        rating, Groups will _not_ have a summary count for Identical hits.

        Individual gene ratings for each genome and group are output to `*.buscomp.tdt`. Compiled ratings are output to
        `*.ratings.tdt`.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#==#','# ================ COMPILE MINIMAP2 SEARCHES ================ #')
            self.headLog('COMPILE MINIMAP2 SEARCHES',line='=')
            dataformats = {'AlnNum':'int','BitScore':'num','Expect':'num','Identity':'int',
                           'QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int',
                           'Length':'int','QryLen':'int','HitLen':'int','Qry_AlnLen':'num','Qry_AlnID':'num'}
            ## ~ [1a] Database tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            gendb = self.db('genomes')
            seqdb = self.db('buscoseq')
            if seqdb and not rje.exists(self.getStr('BUSCOFas')):
                compdb = db.copyTable(seqdb,'buscomp')
                compdb.keepFields(['BuscoID','Genome'])
            else: compdb = db.addEmptyTable('buscomp',fields=['BuscoID','Genome'],keys=['BuscoID'])
            ratefields = ['#','Genome','N','Identical','Complete','Single','Duplicated','Fragmented','Partial','Ghost','Missing']
            ratedb = db.addEmptyTable(fields=ratefields,keys=['Genome'],name='ratings')
            ## ~ [1b] BUSCOMPSeq fasta file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Count number of sequences in BUSCOFas file to check against hitsum file.
            fasx = 0
            buscofas = '%s.buscomp.fasta' % self.baseFile()
            faseog = []
            if rje.exists(self.getStr('BUSCOFas')):
                seqdb = None
                buscofas = self.getStr('BUSCOFas')
            if rje.exists(buscofas):
                for line in open(buscofas,'r').readlines():
                    if line.startswith('>'):
                        fasx += 1
                        eog = string.split(line)[0][1:]
                        if not compdb.data(eog):
                            compdb.addEntry({'BuscoID':eog,'Genome':'BUSCOFas'})
                        faseog.append(eog)

            ## ~ [1b] Rating cutoffs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Complete: 95%+ Coverage (Single/Duplicated based on numbers)
            complete = 95
            # Fragmented: 95%+ Coverage in HitSum but not GABLAM
            fragment = 95
            # Partial: 40%+ Coverage
            partial = 40
            # Ghost: <40% Coverage
            # Missing: No hits meeting local cutoff
            # Identical: Parallel counting of 100% matches. (Could be Single or Duplicated)

            ### ~ [2] Load all GABLAM and HitSum tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # - Load all GABLAM and HitSum tables
            for entry in gendb.entries(sorted=True):
                if not entry['Fasta']: continue
                #self.printLog('#RATING','# ~~~~~~~~ %s Minimap2 ratings ~~~~~~~~ #' % self.genomeField(entry))
                self.headLog('%s Minimap2 ratings' % self.genomeField(entry))
                if not rje.exists(entry['Fasta']):
                    self.warnLog('Genome fasta file for %s not found: %s' % (entry['Genome'],entry['Fasta']))
                    continue
                prefix = rje.baseFile(entry['Fasta'],strip_path=True)
                pafbase = self.pafBase(prefix)
                # The hitsum table gives the maximum coverage of the gene
                hitfile = '%s.hitsum.tdt' % pafbase
                hdb = db.addTable(hitfile,mainkeys=['Qry'],name='%s.hitsum' % entry['Genome'])
                if hdb.entryNum() != fasx:
                    self.warnLog('%s has %s entries but %s has %s sequences!' % (hitfile,rje.iStr(hdb.entryNum()),buscofas,rje.iStr(fasx)))
                # The GABLAM table gives the number and completeness on a single contig
                gabfile = '%s.gablam.tdt' % pafbase
                gdb = db.addTable(gabfile,mainkeys=['Qry','Hit'],name='%s.gablam' % entry['Genome'])

                #?# Warn if PAF might not have same cutoffs? Or encode in the PAF names, as with PAGSAT << This #!#

                #i# LocID filtering has now been removed from this section, as it is applied at the PAF generation stage
                #i# Overall HitSum and GABLAM local %identity can be lower than MinLocID
                # Reformat and filter
                for table in [hdb,gdb]:
                    table.dropEntriesDirect('Qry_AlnID','')
                    table.dataFormat(dataformats)

                # - Reclassify
                for field in ['Identical','Complete','Partial','Ghost']: hdb.addField(field,evalue=0)
                hdb.addField('Rating',evalue='Missing')
                hdb.addField('SeqRating',evalue='Missing')
                # Count individual hits
                gx = 0; gtot = gdb.entryNum()
                for gentry in gdb.entries():
                    self.progLog('\r#RATING','Rating %s hits: %.1f%%' % (entry['Genome'],gx/gtot)); gx += 100.0
                    #self.bugPrint('\n\n%s' % gdb.entrySummary(gentry))
                    coverage = gentry['Qry_AlnLen']
                    eog = gentry['Qry']
                    if coverage == 100.0 and gentry['Qry_AlnID'] == 100.0: hdb.data(eog)['Identical'] += 1
                    if coverage >= complete: hdb.data(eog)['Complete'] += 1
                    elif coverage >= partial: hdb.data(eog)['Partial'] += 1
                    else: hdb.data(eog)['Ghost'] += 1
                    #self.deBug(hdb.entrySummary(hdb.data(eog)))
                self.printLog('\r#RATING','Rating %s %s hits complete.' % (rje.iStr(gtot),entry['Genome']))
                # Convert to rating
                if seqdb: seqdb.addField(self.genomeField(entry),evalue='NULL')
                compdb.addField(self.genomeField(entry),evalue='NULL')
                hx = 0; htot = hdb.entryNum()
                rentry = {'Genome':self.genomeField(entry),'Identical':0,'N':htot,'#':entry['#']}
                for hentry in hdb.entries():
                    self.progLog('\r#RATING','Rating %s BUSCOs: %.1f%%' % (entry['Genome'],hx/htot)); hx += 100.0
                    hentry['Rating'] = 'Missing'
                    if hentry['Complete'] > 1: hentry['Rating'] = 'Duplicated'
                    elif hentry['Complete'] > 0: hentry['Rating'] = 'Complete'
                    elif hentry['Partial'] > 0 and hentry['Qry_AlnLen'] >= fragment: hentry['Rating'] = 'Fragmented'
                    elif hentry['Partial'] > 0: hentry['Rating'] = 'Partial'
                    elif hentry['Ghost'] > 0: hentry['Rating'] = 'Ghost'
                    hentry['SeqRating'] = hentry['Rating']
                    if hentry['Identical'] > 0:
                        rentry['Identical'] += 1
                        hentry['SeqRating'] = 'Identical'
                    #self.bugPrint(hentry)
                    if seqdb:
                        try:
                            seqdb.data(hentry['Qry'])[self.genomeField(entry)] = hentry['SeqRating']
                        except:
                            self.debug(seqdb.data(hentry['Qry']))
                    try:
                        compdb.data(hentry['Qry'])[self.genomeField(entry)] = hentry['Rating']
                    except:
                        self.debug(compdb.data(hentry['Qry']))
                self.printLog('\r#RATING','Rating %s %s BUSCOs complete.' % (rje.iStr(htot),entry['Genome']))

                # Ratings summary
                self.printLog('#SEQID','%s of %s identical to BUSCOSeq' % (rje.iStr(rentry['Identical']),rje.iStr(htot)))
                hdb.indexReport('Rating')
                for field in ['Complete','Duplicated','Fragmented','Partial','Ghost','Missing']:
                    rentry[field] = len(hdb.index('Rating').get(field,[]))
                rentry['Single'] = rentry['Complete']
                rentry['Complete'] += rentry['Duplicated']
                ratedb.addEntry(rentry)
                db.deleteTable(hdb)
                db.deleteTable(gdb)
                #!# Add a one-line summary like the original BUSCO
                self.printLog('#COMP','%s BUSCOMP Rating complete' % (self.genomeString(entry)))
                #self.deBug('^')

            ### ~ [3] Process Groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('COMPILE BUSCOMP GROUP RATINGS',line='=')
            ranklist = ['Complete','Duplicated','Fragmented','Partial','Ghost','Missing','NULL']
            groups = rje.sortKeys(self.dict['Groups'])
            for group in groups:
                # Identify genomes in group
                gentries = []
                for genentry in self.dict['Groups'][group]:
                    if not rje.exists(genentry['Fasta']):
                        self.printLog('#SKIP','Skipping Genome %s: no FASTA file' % self.genomeString(genentry))
                    else: gentries.append(genentry)
                # Identify Group in Genome table
                grpentry = {'Genome':group,'Identical':0,'Single':0,'N':0,'#':self.aliasEntry(group)['#']}
                for field in ranklist: grpentry[field] = 0
                # Add group to BUSCO table
                compdb.addField(group,evalue=ranklist[-1])
                for eog in compdb.dataKeys():
                    #self.debug('%s -> %s -> %s' % (eog,eog in faseog,seqdb.data(eog)['Status']))
                    if rje.exists(self.getStr('BUSCOFas')):
                        if eog not in faseog: continue
                    elif seqdb.data(eog)['Status'] != 'Complete': continue
                    for genentry in gentries:
                        gfield = self.genomeField( genentry )
                        if ranklist.index(compdb.data(eog)[gfield]) < ranklist.index(compdb.data(eog)[group]):
                            compdb.data(eog)[group] = compdb.data(eog)[gfield]
                    grpentry['N'] += 1
                    grpentry[compdb.data(eog)[group]] += 1
                    if compdb.data(eog)[group] == 'Duplicated': grpentry['Complete'] += 1
                    if compdb.data(eog)[group] == 'Complete': grpentry['Single'] += 1
                self.printLog('#GRP','Compiled BUSCOMP stats for group %s (%d of %d genomes)' % (group,len(gentries),len(self.dict['Groups'][group])))
                #!# Adjust for loaded fasta: not tested yet
                #self.debug(grpentry)
                if fasx and grpentry['N'] != fasx:
                    extras =  grpentry['N'] - fasx
                    grpentry['Missing'] -= extras
                    grpentry['N'] -= extras
                #self.debug(grpentry)
                ratedb.addEntry(grpentry)


            return True     # Setup successful
        except: self.errorLog('Problem during %s compileMinimap.' % self.prog()); raise
#########################################################################################################################
    def compSummary(self,log=True,full=False,screen=False):   ### Summarise stats from self.db('ratings')       # v0.5.1
        '''
        Summarise stats from self.db('ratings').
        >> log:bool [True] = Whether to output summary to log
        >> full:bool [False] = Whether to output individual rating numbers, like BUSCO output
        << str = Summary text as a string that could be wrapped in R code etc.

        ### BUSCOMP Summary

        The final step of the BUSCOMP compilation is to summarise the findings in the log afile, akin to the BUSCO
        summary file. This will first generate a one line summary of the percentages, along with the original number of
        complete BUSCOs and the number of BUSCOMP sequences contributed by that assembly (i.e. the number with the best
        score of all the BUSCO searches.) This is followed by a more detailed breakdown of each category. For example:

        ```
        #BUSCO	00:00:41	Results: C:89.2%[S:87.8%,D:1.4%,I:23.3%],F:5.8%,P:3.8%,G:1.6%,M:0.9%,n:3736 - canetoad_v2.2 (3194 Complete BUSCOs; 102 BUSCOMP Seqs)
        #INFO	00:00:41	870 Identical BUSCOs (I)  [100% complete and identical]
        #INFO	00:00:41	3333 Complete BUSCOs (C)  [95%+ coverage in a single contig/scaffold]
        #INFO	00:00:41	3281 Complete and single-copy BUSCOs (S)  [1 Complete hit]
        #INFO	00:00:41	52 Complete and duplicated BUSCOs (D)  [2+ Complete hits]
        #INFO	00:00:41	217 Fragmented BUSCOs (F)  [95%+ coverage spread over 2+ contigs/scaffolds]
        #INFO	00:00:41	143 Partial BUSCOs (P)  [40-95% coverage]
        #INFO	00:00:41	61 Ghost BUSCOs (G)  [<40% coverage]
        #INFO	00:00:41	34 Missing BUSCOs (M)  [No hits]
        #INFO	00:00:41	3736 Total BUSCO gene hits searched
        ```
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            busdb = self.db('buscomp')
            gdb = self.db('ratings')
            rawdb = self.db('busco')
            for gentry in gdb.entries(sorted=True): busdb.index(gentry['Genome'])
            if self.db('buscoseq'):
                seqdb = db.copyTable(self.db('buscoseq'),'best',add=False)
                seqdb.dropEntriesDirect('Status','Complete',inverse=True)
                seqdb.index('Genome')
            else: seqdb = None
            txt = ''
            ### ~ [2] Summary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if log:
                #if full: self.printLog('#==#','# ================ BUSCOMP RESULTS ================ #')
                #else: self.printLog('#==#','# ================ BUSCOMP SUMMARY ================ #')
                if full: self.headLog('BUSCOMP RESULTS',line='=')
                else: self.headLog('BUSCOMP SUMMARY',line='=')
            for gentry in gdb.entries(sorted=True):
                genome = self.genomeField(gentry)
                if full: self.headLog('%s BUSCOMP' % genome)
                if rawdb and genome in rawdb.fields():
                    rnum = len(rawdb.indexEntries(genome,'Complete')) + len(rawdb.indexEntries(genome,'Duplicated'))
                    N = rawdb.entryNum()
                    rperc = 0.0
                    if N: rperc = 100.0 * rnum / N
                    if seqdb: bnum = len(seqdb.indexEntries('Genome',genome))
                    else: bnum = 0
                    bperc = 0.0
                    if bnum: bperc = 100.0 * bnum / gentry['N']
                    seqtxt = '[No BUSCO or BUSCOMP compilation]'
                    if N and gentry['N'] and seqdb: #not rje.exists(self.getStr('BUSCOFas')):
                        seqtxt = '[%d (%.2f%%) Complete BUSCOs; %d (%.2f%%) BUSCOMP Seqs]' % (rnum,rperc,bnum,bperc)
                    elif N:
                        seqtxt = '[%d (%.2f%%) Complete BUSCOs; No BUSCOMP compilation]' % (rnum,rperc)
                    elif gentry['N'] and not rje.exists(self.getStr('BUSCOFas')):
                        seqtxt = '[No BUSCO analysis; %d (%.2f%%) BUSCOMP Seqs]' % (bnum,bperc)
                    #seqtxt = '(%d Complete BUSCOs; %d BUSCOMP Seqs)' % (len(rawdb.indexEntries(genome,'Complete')),len(seqdb.indexEntries('Genome',genome)))
                else:
                    seqtxt = '[No BUSCO analysis]'
                txt += '%s BUSCOMP Results %s:\n' % (genome,seqtxt)
                n = gentry['N'] / 100.0
                if genome in self.dict['Groups']:
                    if n:
                        busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%],F:%.1f%%,P:%.1f%%,G:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete'])/n,gentry['Single']/n,gentry['Duplicated']/n,gentry['Fragmented']/n,gentry['Partial']/n,gentry['Ghost']/n,gentry['Missing']/n,gentry['N'])
                    else: busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%],F:%.1f%%,P:%.1f%%,G:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete']),gentry['Single'],gentry['Duplicated'],gentry['Fragmented'],gentry['Partial'],gentry['Ghost'],gentry['Missing'],gentry['N'])
                else:
                    if n:
                        busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%,I:%.1f%%],F:%.1f%%,P:%.1f%%,G:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete'])/n,gentry['Single']/n,gentry['Duplicated']/n,gentry['Identical']/n,gentry['Fragmented']/n,gentry['Partial']/n,gentry['Ghost']/n,gentry['Missing']/n,gentry['N'])
                    else: busco = 'C:%.1f%%[S:%.1f%%,D:%.1f%%,I%.1f%%],F:%.1f%%,P:%.1f%%,G:%.1f%%,M:%.1f%%,n:%d' % ((gentry['Complete']),gentry['Single'],gentry['Duplicated'],gentry['Identical'],gentry['Fragmented'],gentry['Partial'],gentry['Ghost'],gentry['Missing'],gentry['N'])
                txt += '        %s\n' % busco
                if log: self.printLog('#BUSCO','Results: %s - %s %s' % (busco,genome,seqtxt))
                if not full: txt += '\n'; continue
                businfo = ['%d Identical BUSCOs (I)  [100%% complete and identical]' % (gentry['Identical']),
                    '%d Complete BUSCOs (C)  [95%%+ coverage in a single contig/scaffold]' % (gentry['Complete']),
                    '%d Complete and single-copy BUSCOs (S)  [1 Complete hit]' % gentry['Single'],
                    '%d Complete and duplicated BUSCOs (D)  [2+ Complete hits]' % +gentry['Duplicated'],
                    '%d Fragmented BUSCOs (F)  [95%%+ coverage spread over 2+ contigs/scaffolds]' % gentry['Fragmented'],
                    '%d Partial BUSCOs (P)  [40-95%% coverage]' % gentry['Partial'],
                    '%d Ghost BUSCOs (G)  [<40%% coverage]' % gentry['Ghost'],
                    '%d Missing BUSCOs (M)  [No hits]' % gentry['Missing'],
                    '%d Total BUSCO gene hits searched' % gentry['N']]
                #i# Currently no Identical ratings for groups #i#
                if genome in self.dict['Groups']: businfo.pop(0)
                for bustxt in businfo:
                    txt += '        %s\n' % bustxt
                    if log: self.printLog('#INFO',bustxt)
                txt += '\n'

            if screen and not log: self.verbose(txt,v=0)
            return txt
        except: self.errorLog('Problem during %s sortGenomes.' % self.prog()); raise
#########################################################################################################################
    def groupDifference(self,group1,group2):    ### Returns entries in Group1 not found in Group2
        '''Returns entries in Group1 not found in Group2.'''
        diff = group1[0:]
        for el in group2:
            if el in diff: diff.remove(el)
        return diff
#########################################################################################################################
    def uniqueBUSCOs(self): ### Establish genomes and groups with unique Complete BUSCO genes                   # v0.5.1
        '''
        ### Unique BUSCO genes

        Once BUSCOMP compilation is complete, BUSCOMP will identify "unique" BUSCO hits. These are `Complete` or
        `Duplicated` only in single genome. If present in multiple genomes but only a single Group, that Group will be
        assigned Unique status. (Redundant groups, wholly contained within another group, will be dropped from this
        analysis unless manually retained (interactive mode only). The `BUSCOMP` group woll also be excluded from this
        analysis.) Genes `Complete` in multiple genomes and Groups (or in the absence of groups) will have no Unique
        rating. Genes `Complete` in no genomes are rated "Incomplete" in this analysis.

        The full table of unique ratings is output to `*.unique.tdt` and the number of unique ratings per genome and
        group are reported in the log.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.printLog('#==#','# ================ UNIQUE COMPLETE BUSCOs ================ #')
            self.headLog('UNIQUE COMPLETE BUSCOs',line='=')
            db = self.db()
            busdb = self.db('buscomp')
            gdb = self.db('genomes')
            rawdb = self.db('busco')
            uniqdb = db.addEmptyTable(name='unique',keys=['BuscoID'],fields=['BuscoID','BUSCO','BUSCOMP'])
            ## ~ [1a] Groups for Unique assessment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ugroups = rje.sortKeys(self.dict['Groups'])
            if 'BUSCOMP' in ugroups: ugroups.remove('BUSCOMP')
            redundant = []
            for group in ugroups[0:]:
                for group2 in ugroups[0:]:
                    if group == group2: continue
                    if self.groupDifference(self.dict['Groups'][group],self.dict['Groups'][group2]): continue  # group has entries not in group2
                    if self.dict['Groups'][group] == self.dict['Groups'][group2]:
                        if group2 in redundant: continue
                    redundant.append(group)
            if self.i() >= 0:
                for group in ugroups[0:]:
                    default = {True:'N',False:'Y'}
                    if not rje.yesNo('Include group "%s" in unique BUSCO assessment?' % group,default=default[group in redundant]):
                        ugroups.remove(group)
            else: ugroups = self.groupDifference(ugroups,redundant)
            #i# ugroups is a list of groups not contained in other groups
            #i# ugenomes is a list of entries corresponding to groups and genomes not in groups
            ugenomes = gdb.entries()
            self.debug('%s' % ugenomes)
            for ugroup in ugroups:
                self.bugPrint(ugroup)
                # Remove entries in ugroup from unique genome list
                ugenomes = self.groupDifference(ugenomes,self.dict['Groups'][ugroup])
                self.debug('%s' % ugenomes)

            ### ~ [2] Identify unique BUSCOs and BUSCOMPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rawdb: eoglist = rawdb.dataKeys()
            elif busdb: eoglist = busdb.dataKeys()
            else: self.printLog('#UNIQ','No BUSCO or BUSCOMP data for Unique Complete assessment'); return
            for eog in eoglist:
                buscomp = busco = 'Incomplete'
                ubusco = []
                ubuscomp = []

                # Genomes
                for gentry in gdb.entries(sorted=True):
                    genome = self.genomeField(gentry)
                    if rawdb and genome in rawdb.fields() and rawdb.data(eog)[genome] in ['Complete','Duplicated']:
                        if gentry in ugroups: ubusco.append(gentry)
                        if genome in self.dict['Groups']: pass
                        elif busco == 'Incomplete': busco = genome
                        else: busco = 'Multiple'
                    if busdb and genome in busdb.fields() and busdb.data(eog) and busdb.data(eog)[genome] in ['Complete','Duplicated']:
                        if gentry in ugroups: ubuscomp.append(gentry)
                        if genome in self.dict['Groups']: pass
                        elif buscomp == 'Incomplete': buscomp = genome
                        else: buscomp = 'Multiple'
                    self.deBug('%s:%s -> %s' % (eog,genome,buscomp))

                # Groups
                if busco == 'Multiple' and len(ubusco) == 1: busco = self.genomeField(ubusco[0])
                if buscomp == 'Multiple' and len(ubuscomp) == 1: buscomp = self.genomeField(ubusco[0])

                uniqdb.addEntry({'BuscoID':eog,'BUSCO':busco,'BUSCOMP':buscomp})
            uniqdb.saveToFile()

            ### ~ [3] Report numbers of unique BUSCOs and BUSCOMPs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniqdb.index('BUSCO')
            uniqdb.index('BUSCOMP')
            self.list['UniqText'] = []
            for gentry in gdb.entries(sorted=True):
                genome = self.genomeField(gentry)
                if genome in self.dict['Groups']: continue
                busx = cmpx = 0
                if genome in uniqdb.index('BUSCO'): busx = len(uniqdb.index('BUSCO')[genome])
                if genome in uniqdb.index('BUSCOMP'): cmpx = len(uniqdb.index('BUSCOMP')[genome])
                utext = '%s unique Complete genes: %s BUSCO; %s BUSCOMP' % (genome, rje.iStr(busx), rje.iStr(cmpx))
                self.printLog('#UNIQ',utext)
                self.list['UniqText'].append(utext)
                uniqdb.addEntry({'BuscoID':genome,'BUSCO':busx,'BUSCOMP':cmpx})
            for genome in ugroups:
                busx = cmpx = 0
                if genome in uniqdb.index('BUSCO'): busx = len(uniqdb.index('BUSCO')[genome])
                if genome in uniqdb.index('BUSCOMP'): cmpx = len(uniqdb.index('BUSCOMP')[genome])
                #self.list['UniqText'].append(self.printLog('#UNIQ','%s Group unique Complete genes: %s BUSCO; %s BUSCOMP' % (genome, rje.iStr(busx), rje.iStr(cmpx))))
                utext = '%s unique Complete genes: %s BUSCO; %s BUSCOMP' % (genome, rje.iStr(busx), rje.iStr(cmpx))
                self.printLog('#UNIQ',utext)
                self.list['UniqText'].append(utext)
                uniqdb.addEntry({'BuscoID':genome,'BUSCO':busx,'BUSCOMP':cmpx})

        except: self.errorLog('Problem during %s uniqueBUSCOs.' % self.prog()); raise
#########################################################################################################################
    def buscompChanges(self):   ### Generate a table changes in ratings from BUSCO to BUSCOMP                   # v0.6.0
        '''
        Generate a table changes in ratings from BUSCO to BUSCOMP.

        ## BUSCO versus BUSCOMP comparisons

        There is a risk that performing a low stringency search will identify homologues or pseudogenes of the desired BUSCO gene in error.
        If there is a second copy of a gene in the genome that is detectable by the search then we would expect the same
        genes that go from `Missing` to `Complete` in some genomes to go from `Single` to `Duplicated` in others.

        To test this, data is reduced for each pair of genomes to BUSCO-BUSCOMP rating pairs of:

        * `Single`-`Single`
        * `Single`-`Duplicated`
        * `Missing`-`Missing`
        * `Missing`-`Single`

        This is then converted in to `G`ain ratings (`Single`-`Duplicated` & `Missing`-`Single`) or `N`o Gain ratings
        (`Single`-`Single` & `Missing`-`Missing`). The `Single`-`Duplicated` shift in one genome is then used to set the expected `Missing`-`Single`
        shift in the other, and assess the probability of observing the `Missing`-`Single` shift using a cumulative binomial
        distribution, where:

        * `k` is the number of observed `GG` pairs (`Single`-`Duplicated` _and_ `Missing`-`Single`)
        * `n` is the number of `Missing`-`Single` `G`ains in the focal genome (`NG`+`GG`)
        * `p` is the proportion of `Single`-`Duplicated` `G`ains in the background genome (`GN`+`GG` / (`GN`+`GG`+`NN`+`NG`))
        * `pB` is the probability of observing `k+` `Missing`-`Single` gains, given `p` and `n`

        This is output to `*.gain.tdt`, where each row is a Genome and each field gives the probability of the row
        genome's `Missing`-`Single` gains, given the column genome's `Single`-`Duplicated` gains.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            busdb = self.db('busco')
            comdb = self.db('buscomp')
            if not busdb or not comdb: return False
            self.headLog('BUSCOMP RATING CHANGES',line='=')
            eogdb = self.db().copyTable(busdb,'changes.full')
            for genome in busdb.fields()[1:]:
                if genome in self.dict['Groups'] or genome not in comdb.fields()[1:]:
                    eogdb.dropField(genome)
            newdb = self.db().addEmptyTable('changes',['BUSCO','BUSCOMP'],['BUSCO','BUSCOMP'])
            for busco in ['Complete','Duplicated','Fragmented','Missing']:
                for buscomp in ['Complete','Duplicated','Fragmented','Partial','Ghost','Missing','NULL']:
                    newdb.addEntry({'BUSCO':busco,'BUSCOMP':buscomp})
            newdb.addFields(eogdb.fields()[1:]+['TOTAL'],evalue=0)
            ### ~ [2] Ratings changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for eog in busdb.dataKeys():
                bentry = busdb.data(eog)
                centry = comdb.data(eog)
                eentry = eogdb.data(eog)
                for genome in eogdb.fields()[1:]:
                    if centry:
                        eentry[genome] = '%s%s' % (bentry[genome][0],centry[genome][0])     # First letter of changes
                        nentry = newdb.data((bentry[genome],centry[genome]))
                    else:
                        eentry[genome] = '%sN' % (bentry[genome][0])     # First letter of changes
                        nentry = newdb.data((bentry[genome],'NULL'))
                    #newdb.addEntry({'BuscoID':eog,'Genome':genome,'BUSCO':bentry[genome],'BUSCOMP':centry[genome]})
                    try: nentry[genome] += 1
                    except: self.errorLog('%s %s' % (bentry[genome],centry[genome]))
                    if genome not in self.dict['Groups']:
                        #newdb.addEntry({'BuscoID':'%s.%s' % (eog,genome),'Genome':'TOTAL','BUSCO':bentry[genome],'BUSCOMP':centry[genome]})
                        nentry['TOTAL'] += 1
            #newdb.saveToFile('%s.full-changes.tdt' % self.baseFile())
            ## ~ [2a] Compress to count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #newdb.addField('N',evalue=1)
            #newdb.compress(['Genome','BUSCO','BUSCOMP'],default='sum')
            #newdb.dropField('BuscoID')
            newdb.dropEntriesDirect('TOTAL',[0])
            eogdb.saveToFile()
            newdb.saveToFile()

            ### ~ [3] Gain probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#GAIN','Calculating rating gain probabilities...')
            # There is a risk that performing a low stringency search will identify homologues or pseudogenes of the desired BUSCO gene in error. If there is a second copy of a gene in the genome that is detectable by the search then we would expect the same genes that go from Missing to Complete in some genomes to go from Single to Duplicated in others.
            #
            # To test this, data is reduced for each pair of genomes to:
            # Single-Single & Single-Duplicated vs Missing-Missing & Missing-Single
            # - NoGain & Gain vs NoGain & Gain
            #
            # We can then look for an association between Gain using a chi square test, or use a binomial test against the Gain-Gain ratings, where:
            # k is the number of observed G-G
            # n is the number of G in genome 1
            # p is the proportion of Gain in genome 2
            calcdb = self.db().addEmptyTable('gain.calc',['BUSCO','BUSCOMP'],['BUSCO','BUSCOMP'])
            gaindb = self.db().addEmptyTable('gain',['Genome']+eogdb.fields()[1:],['Genome'])
            for genome1 in eogdb.fields()[1:]:
                for genome2 in eogdb.fields()[1:]:
                    calcdb.addEntry({'BUSCO':genome1,'BUSCOMP':genome2})
            calcdb.addFields(['NN','NG','GN','GG','k','p','n','pB'],evalue=0)
            ## ~ [3a] Count changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in eogdb.entries():
                for genome1 in eogdb.fields()[1:]:
                    for genome2 in eogdb.fields()[1:]:
                        gentry = calcdb.data((genome1,genome2))
                        if entry[genome1] == 'CC' and entry[genome2] == 'MM': gentry['NN'] += 1
                        if entry[genome1] == 'CC' and entry[genome2] == 'MC': gentry['NG'] += 1
                        if entry[genome1] == 'CD' and entry[genome2] == 'MC': gentry['GG'] += 1
                        if entry[genome1] == 'CD' and entry[genome2] == 'MM': gentry['GN'] += 1
            ## ~ [3b] Binomial stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gentry in calcdb.entries():
                gentry['k'] = gentry['GG']
                gentry['n'] = gentry['GG'] + gentry['NG']
                xG = float(gentry['GN']+gentry['GG'])
                if xG:
                    gentry['p'] = xG / (xG+gentry['NN']+gentry['NG'])
                    try:
                        gentry['pB'] = rje.logBinomial(gentry['k'],gentry['n'],gentry['p'],exact=False,callobj=self)
                    except:
                        self.errorLog('Problem calculating probability of {0}+ GG from {1} x p(*G)={2}'.format(gentry['k'],gentry['n'],gentry['p']))
                        gentry['p'] = 0.0
                        gentry['pB'] = 1.0
                else:
                    gentry['p'] = 0.0
                    gentry['pB'] = 1.0
            if self.dev() or self.debugging(): calcdb.saveToFile()
            ## ~ [3c] Probability matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Genome is the BUSCOMP MC v MM changes and the other fields are the probabilities based on the BUSCO
            #i# CC and CD changes in those other genomes.
            for genome1 in eogdb.fields()[1:]:
                gentry = {'Genome':genome1}
                for genome2 in eogdb.fields()[1:]:
                    gentry[genome2] = calcdb.data((genome2,genome1))['pB']
                gaindb.addEntry(gentry)
            gaindb.saveToFile(sfdict={'*':3})

        except: self.errorLog('Problem during %s uniqueBUSCOs.' % self.prog()); raise
#########################################################################################################################
    ### <8> ### BUSCOMP RMarkdown output Methods                                                                        #
#########################################################################################################################
#i# Generate some cut-down tables for Rmd processing (tmp.BASEFILE.*.tdt) then delete at end
# Summarise run details (date, command etc.)
# Summarise genomes
# Plot main genome stats if summarise=T
# Summarise group construction
# BUSCO Summary and compilation
# BUSCOSeq details and ratings for each genome
# BUSCOSeq re-rating summaries and charts
# BUSCOSeq complete table
#########################################################################################################################
    def setupRData(self):   # Compile ratings and genomes tables for plots                                      # v0.5.1
        '''
        ## RMarkdown HTML output

        The final step in BUSCOMP analysis is to generate a summary report. This is produced as an RMarkdown document and
        (assuming the path to pandoc is set) is then "knitted" to make an HTML file. The RMarkdown file is retained to
        make it easy for the user to modify the content and convert the output into a report of their own. Prior to
        generation of the document, results are first compiled into an overview stats file with additional plot
        attributes, before the "best" assemblies are identified under various criteria. Finally, there is an option to
        modify some of the plotting attributes before the report is generated.

        ### Compilation of ratings and genomes tables for summary plots and tables.

        The main BUSCOMP `*.ratings.tdt` output is combined with key genome statistics and BUSCO ratings from the
        `*.genomes.tdt` table. BUSCOMP ratings are converted into percentage values. BUSCO ratings are converted into
        percentage values and reduced to `BUSCO` (Single and Duplicated Complete BUSCOs) and `NoBUSCO` (Missing BUSCOs).
        (Full BUSCO results are plotted directly from the `*.genomes.tdt`.)

        After results are compiled, additional plotting fields are generated:

        * `col` = Plotting colour. (Default, "red". Genomes with 1+ "best" ratings, "blue")
        * `pch` = Point type. (Default, 21 [circle]. Genomes with "best" rating will be 22 [square] for best BUSCO(MP)
        ratings, 24 [triangle] for best contiguity ratings, or 23 [diamond] for best in both categories.
        * `label` = Additional text field to be used for labels. In interactive mode, the option will be given to leave
        this blank for unlabelled plots.
        * `plot` = Whether or not to include a genome in the plot. (Default, TRUE)

        In interactive mode, the option will be provided to edit plotting attributes for each assembly, following the
        calculation of "best" assemblies (below). Compiled data are then saved to `*.rdata.tdt` for summary plots and
        tables in the RMarkdown output.

        ### Identifying the "best" assemblies

        There are many ways of assessing genome assembly quality, but they can be broadly broken down into four aspects:

        1. **Coverage.** How much of the genome is included in the assembly.

        2. **Contiguity.** How many fragments is the genome in, and how big are they.

        3. **Accuracy.** How accurate is the assembly in terms of sequence errors.

        4. **Redundancy.** How much of the genome has been included multiple times.

        Standard reference-free genome statistics (e.g. number, length and gappiness of scaffolds), can only give limited
        insights. Whilst assemblies smaller than the genome size are clearly missing bits, assembly size could be
        artificially inflated by redundant regions, making it impossible to assess Coverage. Scaffold sizes can give an
        indicator of Contiguity, but are also prone to biases introduced by filtering of small scaffolds, for example.
        If an estimated Genome Size can be provided, the most useful measures in this respect are `NG50` and `LG50`.
        These statistics sort scaffolds by length (big to small) and identify the contig/scaffold size at which half the
        *genome* (not the *assembly*) is covered by contigs/scaffolds of that size or bigger. This is the `NG50` value
        (called `NG50Length` in BUSCOMP), with `LG50` (`LG50Count`) being the number of contigs/scaffolds needed to
        cover half the genome. (`N50` and `L50` are the same statistics only for the assembly.) These statistics can
        still be mislead by redundancy in the assembly. None of theses statistics speak to sequence accuracy.

        The power of BUSCO is that it speaks to all four aspects of quality. `Complete` and `Fragmented` genes give
        an indication of Coverage, Continuity and Accuracy; `Duplicated` genes give an indication of Redundancy;
        `Missing` genes give an indication of Coverage and Accuracy. However, the weakness is that these different
        aspects cannot easily be disentangled. This makes side-by-side comparisons of different assemblies challenging,
        as it is not always clear what a difference is indicating.

        BUSCOMP is designed on the principle that **Coverage** and **Contiguity** are the two most important aspects of
        assembly quality. *Accuracy* can, in principle, be improved by additional error-correction steps (including
        manual curation). Suspected *Redundancy* can likewise be identified a flagged. *Coverage* and *Contiguity*, in
        contrast, can only be improved by further assembly - either trying again from scratch, or employing a
        scaffolding or gap-filling alogrithm.

        BUSCOMP has therefore identified seven statistics that can be used to rank assemblies on inferred Completeness
        or Contiguity:

        * **Completeness**:
            1. `Complete` = Maximum number of Complete BUSCOMP ratings.
            2. `BUSCO` = Maximum number of Complete BUSCO ratings.
            3. `Missing` = Smallest number of Missing BUSCOMP ratings.
            4. `NoBUSCO` = Smallest number of Missing BUSCO ratings.
        * ** Contiguity**:
            1. `MaxLength` = Maximum length of contigs/scaffolds.
            2. `NG50Length` = Half the genome lies on contigs/scaffolds at least this size. (See above.)
            3. `LG50Count = Half the genome can be covered by this many contigs/scaffolds. (See above.)

        Individual assemblies are rated as "best" under all seven criteria, with ties allowed. Each assembly can
        be best under multiple criteria and each criterion can have several best assemblies. BUSCOMP will report all
        such combinations.

        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SETUP RMARKDOWN DATA',line='=')
            statfields = string.split('SeqNum TotLength MinLength MaxLength MeanLength MedLength N50Length L50Count CtgNum N50Ctg L50Ctg NG50Length LG50Count GapLength GapCount GC')
            statfields = rje.listIntersect(statfields,self.db('genomes').fields())
            if self.db('ratings'):
                rdb = self.db().copyTable('ratings','rdata')
            else:
                rdb = self.db().copyTable('genomes','rdata')
                rdb.keepFields(['#','Genome'])
                rdb.newKey(['Genome'])
                for field in string.split('N Identical	Complete	Single	Duplicated	Fragmented	Partial	Ghost	Missing'):
                    rdb.addField(field,evalue=0)
            iformats = {}
            for ifield in string.split('SeqNum	TotLength	MinLength	MaxLength	N50Length L50Count CtgNum N50Ctg L50Ctg NG50Length LG50Count GapLength GapCount') + ['UniqBUSCO','UniqBUSCOMP']:
                iformats[ifield] = 'int'
            for nfield in string.split('MeanLength MedLength GC'):
                iformats[nfield] = 'num'
            #i# BUSCO is the percentage of complete BUSCOs
            #i# NoBUSCO is the percentage of missing BUSCOs
            rdb.addFields(statfields+['BUSCO','NoBUSCO','UniqBUSCO','UniqBUSCOMP'])
            uniqdb = self.db('unique')
            ## ~ [1a] Update data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for gentry in self.db('genomes').entries():
                genome = gentry['Genome']
                if genome not in rdb.dataKeys():
                    rentry = rje.combineDict({},gentry)
                    for field in string.split('N Identical	Complete	Single	Duplicated	Fragmented	Partial	Ghost	Missing'): rentry[field] = 0
                    rdb.addEntry(rentry)
                rentry = rdb.data(genome)
                for field in statfields: rentry[field] = gentry[field]
                if gentry['N']:
                    rentry['BUSCO'] = 100.0 * (gentry['Complete'] + gentry['Duplicated'])  / gentry['N']
                    rentry['NoBUSCO'] = 100.0 * gentry['Missing'] / gentry['N']
                else: rentry['BUSCO'] = rentry['NoBUSCO'] = ''
                if rentry['N']:
                    for field in string.split('Identical	Complete	Single	Duplicated	Fragmented	Partial	Ghost	Missing'):
                        rentry[field] = 100.0 * rentry[field] / rentry['N']
                genome = self.genomeField(gentry)
                uentry = uniqdb.data(genome)
                if uentry:
                    rentry['UniqBUSCO'] = uentry['BUSCO']
                    rentry['UniqBUSCOMP'] = uentry['BUSCOMP']
                else: rentry['UniqBUSCO'] = rentry['UniqBUSCOMP'] = -1
            rdb.dataFormat(iformats)

            ### ~ [2] Add plotting fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Point colour ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pcol = 'red'
            if self.i() >= 0: pcol = rje.choice('Default genome stat plot colour?: ',default=pcol,confirm=True)
            rdb.addField('col',evalue=pcol)
            ## ~ [2b] Point shape ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pch = '21'
            if self.i() >= 0: pch = rje.getInt('Default genome stat point type?: ',default=pch,confirm=True)
            rdb.addField('pch',evalue=pch)
            ## ~ [2c] Point label ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.i() < 0 or rje.yesNo('Label points?'):
                rdb.makeField('#Genome#','label')
            else: rdb.addField('label',evalue='')
            rdb.addField('plot',evalue='TRUE')

            ### ~ [3] Best Assemblies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Best'] = {}
            ## ~ [3a] Identify the the best assemblies under different criteria ~~~~~~~~~~~~~~~~~~~ ##
            # Complete	Missing	LG50Count	MaxLength	NG50Length	BUSCO	NoBUSCO
            #i# Indices will generate {field:{values:[genomes]}}
            #for field in rdb.fields(): rdb.index(field)
            rdb.addField('best')
            sumentries = []
            busentries = []
            for entry in rdb.entries():
                entry['best'] = []
                gentry = self.aliasEntry(entry['Genome'],expect=True)
                if gentry['Fasta']: sumentries.append(entry)
                if gentry['Directory']: busentries.append(entry)
            #i# Sequence summary data & key BUSCO counts
            maxsum = string.split('Complete MaxLength NG50Length')
            minsum = string.split('Missing LG50Count')# SeqNum')
            maxbus = string.split('BUSCO')
            minbus = string.split('NoBUSCO')
            for field in maxsum + minsum + maxbus + minbus:
                rdb.index(field)
                rval = []
                if field in maxsum + minsum: rval = rdb.dataList(sumentries,field,sortunique=True,empties=False)
                if field in maxbus + minbus: rval = rdb.dataList(busentries,field,sortunique=True,empties=False)
                while rval and field in ['LG50Count'] and rval[0] <= 0: rval = rval[1:]
                self.devLog('#RVAL','%s %s' % (field,rval))
                self.dict['Best'][field] = {}
                if not rval: continue
                if field in maxsum + maxbus: self.dict['Best'][field] = {rval[-1]:[]}
                if field in minsum + minbus: self.dict['Best'][field] = {rval[0]:[]}
            #i# Update with genomes
            for field in string.split('NG50Length LG50Count MaxLength  Complete Missing BUSCO NoBUSCO'):
                if not self.dict['Best'][field]: continue
                for rval in self.dict['Best'][field]:
                    if field in string.split('Complete Missing BUSCO NoBUSCO'):
                        self.printLog('#BEST','Best %s: %s%%' % (field,rje.dp(rval,2)))
                    else:
                        self.printLog('#BEST','Best %s: %s' % (field,rval))
                    for genome in rdb.index(field)[rval]:
                        if genome in self.dict['Groups']: continue
                        self.dict['Best'][field][rval].append(genome)
                        rdb.data(genome)['best'].append(field)

            # maxlist = string.split('Complete MaxLength NG50Length BUSCO')
            # for field in maxlist:
            #     maxval = max(rdb.index(field).keys())
            #     self.dict['Best'][field] = {maxval:[]}
            #     for genome in rdb.index(field)[maxval]:
            #         self.dict['Best'][field][maxval].append(genome)
            #         rdb.data(genome)['best'].append(field)
            # minlist = string.split('Missing LG50Count NoBUSCO')
            # for field in minlist:
            #     minval = -1
            #     while minval < 0 and rdb.index(field):
            #         if minval in rdb.index(field): rdb.index(field).pop(minval)
            #         minval = min(rdb.index(field).keys())
            #     self.dict['Best'][field] = {minval:[]}
            #     for genome in rdb.index(field)[minval]:
            #         self.dict['Best'][field][minval].append(genome)
            #         rdb.data(genome)['best'].append(field)

            complist = string.split('Complete BUSCO Missing NoBUSCO')
            contlist = string.split('MaxLength NG50Length LG50Count')
            for entry in rdb.entries():
                if rje.listIntersect(entry['best'],complist) and rje.listIntersect(entry['best'],contlist): entry['pch'] = '23'
                elif rje.listIntersect(entry['best'],complist): entry['pch'] = '22'
                elif rje.listIntersect(entry['best'],contlist): entry['pch'] = '24'
                if entry['best']: entry['col'] = 'blue'
                entry['best'] = string.join(entry['best'],'|')

            ### ~ [4] Review ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.yesNo('Edit plotting attributes?',default='N',log=False):
                i = 0; sum = True
                while rdb.entryNum():
                    default = 'X'
                    gentry = rdb.entries(sorted=True)[i]
                    gtext = ''
                    if sum: gtext += '\n' + rdb.entrySummary(gentry)
                    sum = True
                    gtext += '\n'
                    if i > 0: gtext += '<P>revious | '
                    if i+1 < rdb.entryNum(): gtext += '<N>ext | '; default = 'N'
                    gtext += '<C>olour | Point <T>ype | <L>abel | '
                    if gentry['plot'] == 'TRUE': gtext += '<D>on\'t plot | '
                    else: gtext += '<D>o plot | '
                    gtext += 'E<x>it edit menu'
                    ## Loop through choices ##
                    choice = rje.choice(gtext,default=default).upper()
                    if choice == 'P': i -= 1
                    elif choice == 'N': i += 1
                    elif choice == 'C':
                        gentry['col'] = rje.choice('Colour?:',default=gentry['col'],confirm=True)
                    elif choice == 'T':
                        gentry['pch'] = rje.choice('Point type?:',default=gentry['pch'],confirm=True)
                    elif choice == 'L':
                        gentry['label'] = rje.choice('Point label (blank for none)?:',default='',confirm=True)
                    elif choice == 'D' and rje.yesNo('Drop %s from BUSCOMP plots?' % gentry['Genome']):
                        if gentry['plot'] == 'TRUE': gentry['plot'] = 'FALSE'
                        else: gentry['plot'] = 'TRUE'
                    elif choice == 'X' and rje.yesNo('End plotting attribute edit?'): break
                    else: sum = False

            ### ~ [5] Sort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            runsort = self.list['RunSort']
            if len(runsort) == 1: runsort = rje.strSentence(runsort[0].lower())
            if runsort in ['Complete','Missing']:
                rdb.rankField(runsort,'#',rev=runsort in ['Complete'],unique=True,warn=False)
                rdb.remakeKeys()

            ### ~ [6] Save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sfdict = {}
            for field in string.split('Identical Complete Single Duplicated Fragmented Partial Ghost Missing GC MeanLength BUSCO NoBUSCO'):
                sfdict[field] = 4
            rdb.saveToFile(sfdict=sfdict)


        except: self.errorLog('Problem during %s setupRData.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def rGenomeField(self,genome):   ### Convert genome field into R-compatible genome field                     # v0.5.1
        '''Convert genome field into R-compatible genome field.'''
        if rje.matchExp('^(\d+)',genome): genome = 'X%s' % genome
        for rep in '- +?<>': genome = genome.replace(rep,'.')
        return genome
#########################################################################################################################
    def buscompRmd(self,fullreport=True):   ### Generate Rmd and HTML output                                     # v0.5.2
        '''Generate Rmd and HTML output.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            summarise = self.getBool('Summarise')
            buscompseq = self.getBool('BUSCOMPSeq')
            gdb = self.db('genomes')
            rbase = self.baseFile()
            if fullreport:
                rbase = '%s.full' % self.baseFile()
                self.headLog('FULL RMARKDOWN OUTPUT',line='=')
            else:
                self.headLog('SUMMARY RMARKDOWN OUTPUT',line='=')
            #idbase = '%s.L%dID%d' % (self.baseFile(),self.getInt('MinLocLen'),self.getInt('MinLocID'))
            idbase = self.baseFile()
            self.baseFile(rje.baseFile(idbase))
            info = self.log.obj['Info']
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            if fullreport:
                rtxt = rmd.rmdHead(title='BUSCOMP Full Report',author='%s BUSCOMP Analysis' % self.baseFile(),setup=True)
            else:
                rtxt = rmd.rmdHead(title='BUSCOMP Summary Report',author='%s BUSCOMP Analysis' % self.baseFile(),setup=True)
            kable = {True:None, False:True}[fullreport]
            if self.getBool('GGPlot'):
                rtxt += '''
```{r libraries, echo=FALSE}
library(ggplot2)
library(ggrepel)
```

'''

            ## ~ [1a] Table of Contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #toc = '## Contents\n\n'
            toc = '---\n\n**Report contents:**\n\n'
            toc += '* <a href="#Top">Run summary</a>\n'
            toc += '* <a href="#Summary">BUSCOMP summary</a>\n'
            toc += '* <a href="#Genomes">Genome summary</a>\n'
            #toc += '* <a href="#Groups">BUSCOMP Groups</a>\n'
            toc += '* <a href="#BUSCO">BUSCO Ratings</a>\n'
            toc += '* <a href="#BUSCOFull">BUSCO full results compilation</a>\n'
            toc += '* <a href="#BUSCOSeq">BUSCOMP Sequence details and rating</a>\n'
            toc += '* <a href="#BUSCOMP">BUSCOMP re-rating of genomes</a>\n'
            toc += '* <a href="#BUSCOSeqFull">BUSCOMP re-rating full results</a>\n'
            toc += '* <a href="#BUSCOMPUnique">Unique BUSCO and BUSCOMP Complete genes</a>\n'
            toc += '* <a href="#BUSCOMissing">Ratings for Missing BUSCO genes</a>\n'
            toc += '* <a href="#Appendix">Appendix: BUSCOMP run details</a>\n'
            toc += '\n---\n\n'
            #rtxt += '<hr>\n%s\n<hr>\n\n' % toc

            ## ~ [1b] Summarise run details (date, command etc.) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rtxt += '<a name="Top" />\n\n# BUSCOMP Run Summary\n\n'
            rtxt += '    %s\n\n' % self.log.runDetails()
            rtxt += '\n\nSee the <a href="#Appendix">run details appendix</a> end of this document for details of '
            rtxt += 'the <a href="%s">log file</a>' % (self.log.info['LogFile'])
            if self.log.list['ErrorLog'] or self.log.list['WarnLog']: rtxt += ', '
            else: rtxt += ' and '
            rtxt += '<a href="#Appendix">commandline parameters</a>'
            if self.log.list['ErrorLog'] or self.log.list['WarnLog']:
                rtxt += ' and runtime <a href="#Errors">BUSCOMP errors and warnings</a>.\n\n'
            else:
                rtxt += '. No runtime <a href="#Errors">BUSCOMP errors or warnings</a> reported.\n\n'
            rtxt += '**NOTE:** To edit this document, open `%s.Rmd` in RStudio, edit and re-knit to HTML.\n\n' % rbase

            ## ~ [1c] Summarise BUSCOMP result ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rdb = self.db('rdata')
            gx = 0; fx = 0; bx = 0; bestx = 0
            bestxt = []
            for entry in rdb.entries():
                if entry['Genome'] in self.dict['Groups']: continue
                if entry['BUSCO'] > 0 or entry['NoBUSCO']:
                    if entry['TotLength'] > 0: bx += 1
                    else: gx += 1
                elif entry['TotLength'] > 0: fx += 1
                if entry['best']:
                    bestx += 1
                    bestxt.append('* `%s`: %s.' % (entry['Genome'],string.replace(entry['best'],'|',', ')))
            #gx -= bx; fx -= bx
            rtxt += '<a name="Summary" />\n\n## BUSCOMP Results Summary\n\n'
            rtxt += 'Assemblies can be assessed on a number of criteria, but the main ones (in the absence of a reference "truth" genome) are either to judge contiguity or completeness.\n'
            rtxt += 'NG50 and LG50 values are based on a genome size of %s.\n' % (rje_seqlist.dnaLen(self.getNum('GenomeSize'),dp=0,sf=3))
            rtxt += 'If the `genomesize=X` parameter was not set (see command list in <a href="#Appendix">appendix</a>), '
            rtxt += 'this will be based on the longest assembly (see sequence stats, below).\n\n'

            rtxt += 'Of the %d assemblies analysed (%s BUSCO; %s fasta; %d both),\n' % (gx + fx + bx,bx+gx,bx+fx,bx)
            rtxt += '%d genomes were rated as the "best" by at least one criterion:\n\n%s\n' % (bestx, string.join(bestxt,'\n'))

            #
            if summarise:
                rtxt += '\nBest assemblies by assembly contiguity critera:\n\n'
                rtxt += '* **NG50Length.** Longest NG50 contig/scaffold length (%s bp): `%s`\n' % (rje.iStr(self.dict['Best']['NG50Length'].keys()[0]),string.join(self.dict['Best']['NG50Length'].values()[0],'`, `'))
                rtxt += '* **LG50Count.** Smallest LG50 contig/scaffold count (%s): `%s`\n' % (rje.iStr(self.dict['Best']['LG50Count'].keys()[0]),string.join(self.dict['Best']['LG50Count'].values()[0],'`, `'))
                rtxt += '* **MaxLength.** Maximum contig/scaffold length (%s bp): `%s`\n' % (rje.iStr(self.dict['Best']['MaxLength'].keys()[0]),string.join(self.dict['Best']['MaxLength'].values()[0],'`, `'))
                #rtxt += '* **SeqNum.** Minimum contig/scaffold count (%s): `%s`\n' % (rje.iStr(self.dict['Best']['SeqNum'].keys()[0]),string.join(self.dict['Best']['SeqNum'].values()[0],'`, `'))
            else: rtxt += '\n**NOTE:** `summarise=F` - no contiguity assessments.\n\n'

            rtxt += '\nBest assemblies by completeness critera:\n\n'
            if buscompseq:
                rtxt += '* **Complete.** Most Complete (Single & Duplicated) BUSCOMP sequences (%.1f %%): `%s`\n' % (self.dict['Best']['Complete'].keys()[0],string.join(self.dict['Best']['Complete'].values()[0],'`, `'))
                rtxt += '* **Missing.** Fewest Missing BUSCOMP sequences (%.1f %%): `%s`\n' % (self.dict['Best']['Missing'].keys()[0],string.join(self.dict['Best']['Missing'].values()[0],'`, `'))
            else: rtxt += '\n**NOTE:** `buscompseq=F` - no BUSCOMP compilation assessments.\n\n'
            if self.db('busco'):
                rtxt += '* **BUSCO.** Most Complete (Single & Duplicated) BUSCO sequences (%.1f %%): `%s`\n' % (self.dict['Best']['BUSCO'].keys()[0],string.join(self.dict['Best']['BUSCO'].values()[0],'`, `'))
                rtxt += '* **NoBUSCO.** Fewest Missing BUSCO sequences (%.1f %%): `%s`\n' % (self.dict['Best']['NoBUSCO'].keys()[0],string.join(self.dict['Best']['NoBUSCO'].values()[0],'`, `'))
            else: rtxt += '\n**NOTE:** No BUSCO results loaded.\n\n'


            ### ~ [2] Summarise genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Load the genomes table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            genfile = '%s.genomes.tdt' % self.baseFile()
            codename = 'genomesetup'
            rtxt += '<a name="Genomes" />\n\n# Genome Summary\n\n'
            rcode = rmd.rmdTable('%s.rdata.tdt' % idbase,name='rdata',codesection=False,loadtable=True,showtable=False,delim='tab')
            rcode += 'rdata = rdata[order(rdata[,1]),]\n\n'
            rcode += rmd.rmdTable(genfile,name='gentable',codesection=False,loadtable=True,showtable=False,delim='tab',kable=None,rows=10,cols=10)
            rcode += 'gentable$Genome = rdata$Genome\n'
            rcode += 'rownames(gentable) = gentable[,1]\n'
            rcode += 'gensum = gentable[! is.na(gentable$N),2:6]\n'
            rcode += 'gensum = gensum[! (gensum$Directory == "" & gensum$Fasta == ""),]\n'
            #i# v0.9.4 changes
            # old statfields = string.split('SeqNum TotLength MinLength MaxLength MeanLength MedLength N50Length L50Count NG50Length LG50Count GapLength GC')
            # -> statfields = string.split('SeqNum TotLength MinLength MaxLength MeanLength MedLength N50Length L50Count CtgNum N50Ctg L50Ctg NG50Length LG50Count GapLength GapCount GC')
            #i# Scan number of fields in gentable and pick appropriate field - update to use tidyverse at some point
            rcode += 'nfield = which(colnames(gentable)=="N")\n'
            rcode += 'genstat = gentable[! is.na(gentable$SeqNum),c(4,9:nfield-1)]\n'
            rcode += 'busco = gentable[,c(4,nfield:ncol(gentable))]\n'
            rcode += 'busco$Single = busco$Complete\n'
            rcode += 'busco$Complete = busco$Single + busco$Duplicated\n'
            rcode += 'busco = busco[,c(1,2,4,7,3,5:6)]\n'
            rmd.list['CodeChunks'].append(codename)
            rtxt += '```{r %s, echo=FALSE}\n%s\n```\n\n' % (codename,rcode)
            ## ~ [2b] Summarise genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rtxt += 'The following genomes and BUSCO results were analysed by BUSCOMP:\n\n'
            for gentry in gdb.entries(sorted=True):
                gtype = []
                if gentry['Directory']: gtype.append('<code>BUSCO</code>')
                if gentry['Fasta']: gtype.append('<code>Fasta</code>')
                if not gtype: continue #gtype = ['<code>Group</code>']
                rtxt += '* **%s**. [%s] %s\n' % (gentry['Genome'],string.join(gtype,'|'),gentry['Description'])
            rtxt += '\nDetails of the directories and files are below:\n\n'
            rcode = rmd.rmdTable(genfile,name='gensum',codesection=True,loadtable=False,showtable=True,delim='tab',kable=kable,rows=10,cols=10)
            rtxt += rcode
            rtxt += '\nGenomes with a **Directory** listed had BUSCO results available.\n'
            rtxt += 'If **Sequences** is `True`, these would be have been compiled to generate the BUSCOMP sequence set (unless `buscompseq=F`, or alternative sequences were provided with `buscofas=FASFILE`).\n'
            rtxt += 'Genomes with a **Fasta** listed had sequence data available for BUSCOMP searches.\n\n'
            ## ~ [2b] Genome stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Summarise'):
                rtxt += genomeStats  #'`RJE_SEQ` summary statistics were also calculated for these genomes:\n\n'
                rtxt += '\n\n'
                rcode = rmd.rmdTable(genfile,name='genstat',codesection=True,loadtable=False,showtable=True,delim='tab',kable=True,rows=10,cols=10)
                rtxt += rcode
                rtxt += '\n**NOTE:** `NG50Length` and `LG50Count` statistics use `genomesize=X` or the biggest assembly loaded (%s). ' % rje_seqlist.dnaLen(self.getNum('GenomeSize'))
                rtxt += 'If BUSCOMP has been run more than once on the same data (_e.g._ to update descriptions or sorting), '
                rtxt += 'please make sure that a consistent genome size is used, or these values may be wrong. '
                rtxt += 'If in doubt, run with `force=T` and force regeneration of statistics.\n\n'

                # Plot main genome stats if summarise=T
                rcode = 'pdata = rdata[rdata$plot == TRUE & ! is.na(rdata$MeanLength),]\n\n'

                #!# Future: make this conversion reactive to max value
                rcode += 'pdata$TotLength = pdata$TotLength / 1e9\n'
                rcode += 'pdata$MaxLength = pdata$MaxLength / 1e6\n'
                rcode += 'pdata$N50Length = pdata$N50Length / 1e6\n'
                #!# Make this responsive to max SeqNum: rcode += 'pdata$SeqNum = pdata$SeqNum / 1e3\n'
                if 'NG50Length' in self.db('genomes').fields():
                    rcode += 'pdata$NG50Length = pdata$NG50Length / 1e6\n'
                if 'NG50Ctg' in self.db('genomes').fields():
                    rcode += 'pdata$NG50Ctg = pdata$NG50Ctg / 1e6\n'
                # Output setup code - want to make individual plots
                codename = 'plotsetup'
                rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (codename,rcode)
                rcode = ''

                # Replace plots with separate code chunks for...
                rtxt += '\n## Genome coverage assessment plots\n\n'

                # PLOT 1. General assembly stats = No. scaffolds (1000s) versus Assembly size
                rtxt += 'In general, a good assembly will be approx. the same size as the genome and in as few pieces '
                rtxt += 'as possible. Any assembly smaller than the predicted genome size is clearly missing coverage. '
                rtxt += 'Assemblies bigger than the genome size might still be missing chunks of the genome if redundancy/duplication is a problem. '
                rtxt += 'In the following plot, the grey line marks the given genome size of %s.\n\n' %  (rje_seqlist.dnaLen(self.getNum('GenomeSize'),0))
                genGB = self.getNum('GenomeSize') / 1e9
                #!# Make this responsive to max SeqNum: yfield = 'SeqNum'; ylab = 'No. contigs/scaffolds (1000s)'
                yfield = 'SeqNum'; ylab = 'No. contigs/scaffolds'
                if max(gdb.dataList(gdb.entries(),'SeqNum')) > 10000:
                    yfield = 'SeqNum/1000'; ylab = 'No. contigs/scaffolds (1000s)'
                xfield = 'TotLength'; xlab = 'Assembly Size (Gb)'
                title = 'Contig/scaffold count versus total assembly size'
                plims = 'xlim=c(0,max(pdata$%s,%s)),ylim=c(0,max(pdata$%s))' % (xfield,genGB,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab,xmax=genGB,vline=genGB,xlog=False,ylog=False)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",%s)\n' % (xfield,yfield,xlab,ylab,title,plims)   # Removed ,log="xy"
                    #!# Make this dashed!
                    rcode += 'abline(v=%s, col="grey")\n' % genGB
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                rtxt += '```{r numVsize, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)

                # PLOT 2a. Completeness = NoBUSCO vs Assembly Size
                rtxt += 'A better indicator of the overall coverage of the genome is the number of `Missing` BUSCO '
                rtxt += 'genes. As BUSCO is highly dependent on the accuracy of the sequence and the gene models it makes, '
                rtxt += 'the `Missing` BUSCOMP ratings arguably give a more consistent proxy for genome completeness. '
                rtxt += 'NOTE: this says nothing about the fragmentation or completeness of the genes themselves.\n\n'
                yfield = 'NoBUSCO'; ylab = 'BUSCO Missing (%)'
                title = 'Missing BUSCO versus total assembly size'
                plims = 'xlim=c(0,max(pdata$%s,%s)),ylim=c(0,max(pdata$%s))' % (xfield,genGB,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab,xmax=genGB,vline=genGB)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    #!# Make this dashed!
                    rcode += 'abline(v=%s, col="grey")\n' % genGB
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                if self.db('busco'):
                    rtxt += '```{r nobuscoVsize, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)
                else: rtxt += '\n**NOTE:** No BUSCO results loaded.\n\n'

                # PLOT 2b. Completeness = Missing vs Assembly Size
                yfield = 'Missing'; ylab = 'BUSCOMP Missing (%)'
                title = 'Missing BUSCOMP versus total assembly size'
                plims = 'xlim=c(0,max(pdata$%s,%s)),ylim=c(0,max(pdata$%s))' % (xfield,genGB,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab,xmax=genGB,vline=genGB)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    #!# Make this dashed!
                    rcode += 'abline(v=%s, col="grey")\n' % genGB
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                rtxt += '```{r missingVsize, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)

                # PLOT 3. BUSCO versus BUSCOMP 1 = Missing vs NoBUSCO
                #rtxt += 'A better indicator of the overall coverage of the genome is the number of `Missing` BUSCO'
                #rtxt += 'genes. As BUSCO is highly dependent on the accuracy of the sequence and the gene models it makes,'
                #rtxt += 'the `Missing` BUSCOMP ratings arguably give a more consistent proxy for genome completeness.'
                #rtxt += 'NOTE: this says nothing about the fragmentation or completeness of the genes themselves.\n\n'
                yfield = 'Missing'; ylab = 'BUSCOMP Missing (%)'
                xfield = 'NoBUSCO'; xlab = 'BUSCO Missing (%)'
                title = 'Missing BUSCOMP versus Missing BUSCO'
                plims = 'xlim=c(0,max(pdata$%s)),ylim=c(0,max(pdata$%s))' % (xfield,yfield)
                #!# Add calculation of size to make square boxes.
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                if self.db('busco'):
                    rtxt += '```{r nobuscompVnobusco, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)
                else: rtxt += '\n**NOTE:** No BUSCO results loaded.\n\n'

                rtxt += '\n## Genome contiguity assessment plots\n\n'

                # PLOT 4. General assembly stats = LG50 vs NG50
                rtxt += 'In general, a good assembly will be in fewer, bigger pieces. This is approximated using NG50 and LG50, '
                rtxt += 'which are the min. length and number of contigs/scaffolds required to cover at least half the genome. '
                rtxt += 'These stats use the given genome size of %s.\n\n'  %  (rje_seqlist.dnaLen(self.getNum('GenomeSize'),0))
                yfield = 'LG50Count'; ylab = 'LG50 contig/scaffold count'
                xfield = 'NG50Length'; xlab = 'NG50 (Mb)'
                title = 'LG50 count versus NG50 length'
                plims = 'xlim=c(0,max(pdata$%s)),ylim=c(0,max(pdata$%s))' % (xfield,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                rtxt += '```{r lg50Vng50, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)

                # PLOT 5b. Contiguity = Complete vs NG50
                yfield = 'Complete'; ylab = 'BUSCOMP Complete (%)'
                title = 'Complete BUSCOMP versus NG50'
                plims = 'xlim=c(0,max(pdata$%s)),ylim=c(0,max(pdata$%s))' % (xfield,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                rtxt += '```{r completeVng50, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)

                # PLOT 5c. Contiguity = BUSCO vs NG50
                yfield = 'BUSCO'; ylab = 'BUSCO Complete (%)'
                title = 'Complete BUSCO versus NG50'
                plims = 'xlim=c(0,max(pdata$%s)),ylim=c(0,max(pdata$%s))' % (xfield,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                if self.db('busco'):
                    rtxt += '```{r buscoVng50, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)
                else: rtxt += '\n**NOTE:** No BUSCO results loaded.\n\n'

                # PLOT 3. BUSCO versus BUSCOMP 2 = Complete vs BUSCO
                yfield = 'Complete'; ylab = 'BUSCOMP Complete (%)'
                xfield = 'BUSCO'; xlab = 'BUSCO Complete (%)'
                title = 'Complete BUSCOMP versus Complete BUSCO'
                plims = 'xlim=c(0,max(pdata$%s)),ylim=c(0,max(pdata$%s))' % (xfield,yfield)
                if self.getBool('GGPlot'):
                    rcode = ggPlotCode(xfield,yfield,title,xlab,ylab)
                else:
                    rcode = 'plot(pdata$%s,pdata$%s,xlab="%s",ylab="%s",col=pdata$col,pch=pdata$pch,main="%s",xlog=TRUE,ylog=TRUE,%s)\n' % (xfield,yfield,xlab,ylab,title,plims)
                    rcode += 'text(pdata$%s,pdata$%s,pdata$label,pos=3,col=pdata$col)\n\n' % (xfield,yfield)
                if self.db('busco'):
                    rtxt += '```{r buscompVbusco, echo=FALSE, fig.width=12, fig.height=8}\n%s\n```\n\n' % (rcode)
                else: rtxt += '\n**NOTE:** No BUSCO results loaded.\n\n'


            rtxt += '**NOTE:** To modify these plots and tables, edit the `*.genomes.tdt` and `*.NxLxxIDxx.rdata.tdt` files and re-knit the `*.NxLxxIDxx.Rmd` file.\n\n'




            ### ~ [4] BUSCO Summary and compilation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            grpdb = self.db('groups')
            #groups = rje.sortKeys(grpdb.index('Group'))
            groups = []
            for gentry in gdb.entries(sorted=True):
                if gentry['Genome'] in self.dict['Groups']: groups.append(gentry['Genome'])
            ### Group rows data ###
            grouprows = {}
            if groups:
                for group in groups:
                    grouprows[group] = []
                    genomes = grpdb.indexDataList('Group',group,'Genome')
                    for genome in genomes:
                        gentry = self.aliasEntry(genome)
                        grouprows[group].append('%s' % gentry['#'])
                    grouprows[group].append('%s' % self.aliasEntry(group)['#'])
                    grouprows[group] = 'c(%s)' % string.join(grouprows[group],',')

            rtxt += '<a name="BUSCO" />\n\n# BUSCO Ratings\n\n'
            if self.db('busco'):
                rtxt += 'Compiled BUSCO results for %d assemblies and %d groups have been saved in `%s`. %s\n\n' % (gdb.entryNum()-len(groups),len(groups),genfile,buscoText)
                #rcode = 'busco$Genome = rdata$Genome\n'
                rcode = rmd.rmdTable(genfile,name='busco',codesection=True,loadtable=False,showtable=True,delim='tab',kable=True,rows=10,cols=10)
                rtxt += rcode
                #!# Add chart
                codename = 'buscochart'
                rmd.list['CodeChunks'].append(codename)
                rcode = 'sumdata = busco[! is.na(busco$N) & busco$N > 0,c(1,4:7)]\n'
                rcode += 'colnames(sumdata) = c("Dataset", "Complete", "Duplicated", "Fragmented", "Missing")\n'
                rcode += buscoPercPlot
                rcode += 'buscoPercPlot(sumdata,"BUSCO Rating Summary",maketext=%s)\n' % str(fullreport).upper()
                rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(self.db('genomes').entryNum()/2),rcode)
                ## ~ [2c] Summarise group construction  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rtxt += '<a name="Groups" />\n\n## Genome Groups\n\n'
                if groups:
                    rtxt += '`BUSCOMP` compiled the following groups of genomes (where BUSCO data was loaded), keeping the "best" rating for each BUSCO gene\n'
                    rtxt += ' across the group:\n\n'
                    for group in groups:
                        genomes = grpdb.indexDataList('Group',group,'Genome')
                        rtxt += '* **%s**: ' % group
                        for genome in genomes:
                            gentry = self.aliasEntry(genome)
                            if gentry['N']: rtxt += '`%s` ' % self.genomeField(genome,raw=True)
                            else: rtxt += '(`%s`) ' % self.genomeField(genome,raw=True)
                        rtxt = rtxt[:-1] + '\n'
                        #rtxt += '* **%s**: `%s`\n' % (group,string.join(genomes,'`,`'))
                    rtxt += '\n'
                else:
                    rtxt += 'No groups for `BUSCOMP` compilation of BUSCO ratings.\n\n'
                # BUSCO Summary
                rtxt += '<a name="BUSCOSummary" />\n\n## BUSCO Summary\n\n'
                rtxt += '```\n'
                rtxt += self.rawSummary(log=False)
                rtxt += '```\n'

                # BUSCOSeq details and ratings for each genome
                rtxt += '<a name="BUSCOFull" />\n\n## BUSCO Gene Details\n\n'
                dbfile = '%s.busco.tdt' % self.baseFile()
                rtxt += 'Full BUSCO results with ratings for each gene have been compiled in `%s`' % dbfile
                if fullreport:
                    rtxt += ':\n\n'
                    rcode = rmd.rmdTable(dbfile,name='buscofull',codesection=True,loadtable=True,showtable=True,delim='tab',kable=None,rows=10,cols=10)
                    rtxt += rcode
                else: rtxt += '.\n\n'


                ## ~ [2c+] Group charts  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if groups and fullreport:
                    rtxt += '<a name="GroupBUSCO" />\n\n## Genome Group BUSCO charts\n\n'
                    for group in groups:
                        gnum = len(grpdb.indexDataList('Group',group,'Genome')) + 1
                        codename = '%sbuscochart' % string.replace(group,' ','')
                        rmd.list['CodeChunks'].append(codename)
                        rcode = ''
                        rcode += '# %s %s\n' % (group,grouprows[group])
                        rcode += 'sumdata = busco[%s,]\n' % grouprows[group]
                        rcode += 'sumdata = sumdata[! is.na(sumdata$N) & sumdata$N > 0,c(1,4:7)]\n'
                        rcode += 'colnames(sumdata) = c("Dataset", "Complete", "Duplicated", "Fragmented", "Missing")\n'
                        rcode += 'buscoPercPlot(sumdata,"%s BUSCO Rating Summary",maketext=%s)\n' % (group,str(fullreport).upper())
                        rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(gnum/2),rcode)

            else:
                rtxt += 'No BUSCO results loaded.\n\n'

            # BUSCOSeq details and ratings for each genome
            rtxt += '<a name="BUSCOSeq" />\n\n# BUSCOMP Ratings\n\n'
            dbfile = '%s.buscoseq.tdt' % self.baseFile()
            if os.path.exists(dbfile):
                rtxt += 'The best complete BUSCO hit results (based on `Score` and `Length`) have been compiled in `%s`.\n' % dbfile
                rtxt += 'The `Genome` field indicates the assembly with the best hit, which is followed by details of that hit '
                rtxt += '(`Contig`, `Start`, `End`, `Score`, `Length`).\n'
                rtxt += '''BUSCOMP ratings for each assembly are then given in subsequent fields:
    
    * `Identical`: 100% coverage and 100% identity in at least one contig/scaffold.
    * `Complete`: 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
    * `Duplicated`: 95%+ Coverage in 2+ contigs/scaffolds.
    * `Fragmented`: 95%+ combined coverage but not in any single contig/scaffold.
    * `Partial`: 40-95% combined coverage.
    * `Ghost`: Hits meeting local cutoff but <40% combined coverage.
    * `Missing`: No hits meeting local cutoff.
    
    '''
                if fullreport:
                    rtxt += '\n'
                    rcode = rmd.rmdTable(dbfile,name='buscoseq',codesection=True,loadtable=True,showtable=True,delim='tab',kable=None,rows=10,cols=10)
                    rtxt += rcode
            elif buscompseq:
                rtxt += 'BUSCOMP ratings for sequences in `%s` are calculated for each assembly:\n' % self.getStr('BUSCOFas')
                rtxt += '''
    * `Identical`: 100% coverage and 100% identity in at least one contig/scaffold.
    * `Complete`: 95%+ Coverage in a single contig/scaffold. (Note: accuracy/identity is not considered.)
    * `Duplicated`: 95%+ Coverage in 2+ contigs/scaffolds.
    * `Fragmented`: 95%+ combined coverage but not in any single contig/scaffold.
    * `Partial`: 40-95% combined coverage.
    * `Ghost`: Hits meeting local cutoff but <40% combined coverage.
    * `Missing`: No hits meeting local cutoff.
    
    '''


            # BUSCOSeq re-rating summaries and charts
            rtxt += '\n<a name="BUSCOMP" />\n\n## BUSCOSeq Rating Summary\n\n'
            if buscompseq:
                dbfile = '%s.ratings.tdt' % idbase
                rtxt += 'BUSCOMP ratings (see above) are compiled to summary statistics in `%s`. ' % dbfile
                rtxt += 'Note that `Identical` ratings in this table will also be rated as `Complete`, which in turn are `Single` or `Duplicated`.\n'
                rtxt += 'Percentage summaries are plotted below, along with a BUSCO-style one-line summary per assembly/group.\n\n'
                rtxt += '**NOTE:** Group summaries do not include `Identical` ratings.\n\n'
                rcode = rmd.rmdTable(dbfile,name='ratings',codesection=True,loadtable=True,showtable=True,delim='tab',kable=True,rows=10,cols=10)
                rtxt += rcode
                #!# Add chart
                codename = 'buscoseqchart'
                rmd.list['CodeChunks'].append(codename)
                rcode = 'sumdata = ratings[,c(1,5:10)+1]\n'
                rcode += buscompSeqPercPlot
                rcode += 'buscompSeqPercPlot(sumdata,"BUSCOSeq Rating Summary",maketext=%s)\n' % str(fullreport).upper()
                rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(self.db('ratings').entryNum()/2),rcode)
                # BUSCOSeq Summary
                rtxt += '```\n'
                rtxt += self.compSummary(log=False)
                rtxt += '```\n'

                # BUSCOSeq complete table
                rtxt += '<a name="BUSCOSeqFull" />\n\n## BUSCOSeq Full Results Table\n\n'
                dbfile = '%s.buscomp.tdt' % idbase
                rtxt += 'Full BUSCOMP results with ratings for each gene in every assembly and group have been compiled in `%s`' % dbfile
                if fullreport:
                    rtxt += ':\n\n'
                    rcode = rmd.rmdTable(dbfile,name='buscomp',codesection=True,loadtable=True,showtable=True,delim='tab',kable=None,rows=10,cols=10)
                    rtxt += rcode
                else:
                    rtxt += '.\n\n'


                ## ~ [2c+] Group charts  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if groups and fullreport and self.db('buscomp'):
                    rtxt += '<a name="GroupBUSCOMP" />\n\n## Genome Group BUSCOMP charts\n\n'
                    for group in groups:
                        gnum = len(grpdb.indexDataList('Group',group,'Genome')) + 1
                        codename = '%sbuscompchart' % string.replace(group,' ','')
                        rmd.list['CodeChunks'].append(codename)
                        rcode = ''
                        rcode += '# %s %s\n' % (group,grouprows[group])
                        rcode += 'sumdata = ratings[%s,c(1,5:10)+1]\n' % grouprows[group]
                        rcode += 'buscompSeqPercPlot(sumdata,"%s BUSCOSeq Rating Summary",maketext=%s)\n' % (group,str(fullreport).upper())
                        rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(gnum/2),rcode)

            else: rtxt += 'No BUSCOSeq Rating analysis (`buscompseq=F`).'

            rtxt += '<a name="BUSCOMParisons" />\n\n# BUSCO and BUSCOMP Comparisons\n\n'

            if self.db('busco') and self.db('buscomp') and buscompseq:

                # Changes of Ratings
                rtxt += '<a name="BUSCOMPChanges" />\n\n## BUSCO to BUSCOMP Rating Changes\n\n'
                rtxt += 'Ratings changes from BUSCO to BUSCOMP (where `NULL` ratings indicate no BUSCOMP sequence):\n\n'
                dbfile = '%s.changes.tdt' % idbase
                rtxt += rmd.rmdTable(dbfile,name='changes',codesection=True,loadtable=True,showtable=True,delim='tab',kable=True,rows=10,cols=10)
                rtxt += '\n\n'

                if fullreport:
                    rtxt += 'Full table of Ratings changes from by gene:\n\n'
                    dbfile = '%s.changes.full.tdt' % idbase
                    rtxt += rmd.rmdTable(dbfile,name='fullchanges',codesection=True,loadtable=True,showtable=True,delim='tab',kable=None,rows=10,cols=10)
                    rtxt += '\n\n'
                    rtxt += '<small>**C**omplete, **D**uplicated, **F**ragmented, **P**artial, **G**host, **M**issing, **N**ULL (no BUSCOMP sequence)</small>\n\n'

                rtxt += '### BUSCOMP Gain test\n\n'
                rtxt += '''There is a risk that performing a low stringency search will identify homologues or pseudogenes of the desired BUSCO gene in error.
If there is a second copy of a gene in the genome that is detectable by the search then we would expect the same
genes that go from `Missing` to `Complete` in some genomes to go from `Single` to `Duplicated` in others.

To test this, data is reduced for each pair of genomes to BUSCO-BUSCOMP rating pairs of:

* `Single`-`Single`
* `Single`-`Duplicated`
* `Missing`-`Missing`
* `Missing`-`Single`

This is then converted in to `G`ain ratings (`Single`-`Duplicated` & `Missing`-`Single`) or `N`o Gain ratings
(`Single`-`Single` & `Missing`-`Missing`). The `Single`-`Duplicated` shift in one genome is then used to set the expected `Missing`-`Single`
shift in the other, and assess the probability of observing the `Missing`-`Single` shift using a cumulative binomial
distribution, where:

* `k` is the number of observed `GG` pairs (`Single`-`Duplicated` _and_ `Missing`-`Single`)
* `n` is the number of `Missing`-`Single` `G`ains in the focal genome (`NG`+`GG`)
* `p` is the proportion of `Single`-`Duplicated` `G`ains in the background genome (`GN`+`GG` / (`GN`+`GG`+`NN`+`NG`))
* `pB` is the probability of observing `k+` `Missing`-`Single` gains, given `p` and `n`

This is output to `*.gain.tdt`, where each row is a Genome and each field gives the probability of the row
genome's `Missing`-`Single` gains, given the column genome's `Single`-`Duplicated` gains:

'''
                dbfile = '%s.gain.tdt' % idbase
                rtxt += rmd.rmdTable(dbfile,name='gain',codesection=True,loadtable=True,showtable=True,delim='tab',kable=kable,rows=15,cols=15)
                rtxt += '\n\n'
                #!# Add report of the lowest probabilities?
                rtxt += 'Low probabilities indicate that BUSCOMP might be rating paralogues or pseudogenes and not functional orthologues of the BUSCO gene.\n'
                rtxt += 'Note that there is *no* correction for multiple testing, nor any adjustment for lack of independence between samples.\n\n'
            else: rtxt += 'No comparison without both BUSCO and BUSCOMP rating.\n\n'

            # BUSCOMPUnique unique ratings
            rtxt += '<a name="BUSCOMPUnique" />\n\n## Unique BUSCO and BUSCOMP Complete Genes\n\n'
            rtxt += 'BUSCO and BUSCOMP `Complete` ratings were compared for each BUSCO gene to identify those genes unique '
            rtxt += 'to either a single assembly or a group of assemblies. The `BUSCOMP` group is excluded from this analysis, '
            rtxt += 'as (typically) are other redundant groups wholly contained within another group. (Inclusion of such groups '
            rtxt += 'is guaranteed to result in 2+ groups containing any `Complete` BUSCOs they have.)\n'
            rtxt += '\n```\n%s\n```\n\n' % string.join(self.list['UniqText'],'\n')
            rcode = 'rownames(rdata) = rdata$Genome\n'
            rcode += 'pdata = rdata[rev(order(rdata[,1])),]\n'
            rcode += 'pdata = pdata[pdata$UniqBUSCO >= 0,]\n'
            rcode += 'par(mar=c(5,15,4,1)+0.1)\n'
            rcode += 'barplot(t(as.matrix(pdata[,c("UniqBUSCOMP","UniqBUSCO")])),horiz=TRUE,axes=TRUE,col=c("red","blue"),main="Unique Complete BUSCO Genes",xlab="Number of Unique Genes",names.arg=pdata$Genome,las=1,beside = TRUE,legend=TRUE)\n'
            rcode += 'par(mar=c(5,4,4,2)+0.1)\n'
            rtxt += '```{r uniqueplots, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (2+int(self.db('genomes').entryNum()/2),rcode)

            # Ratings of Missing
            rtxt += '<a name="BUSCOMissing" />\n\n## Ratings for Missing BUSCO genes\n\n'
            rtxt += 'In addition to the unique ratings (above), it can be useful to know how genes `Missing` from one '
            rtxt += 'assembly/group are rated in the others. These plots are generated for each assembly/group in turn. '
            rtxt += 'The full BUSCO (`*.busco.tdt`) and BUSCOMP (`*.LnnIDxx.buscomp.tdt`) tables are reduced to the subset of genes '
            rtxt += 'that are missing in the assembly/group of interest, and then the summary ratings recalculated for that subset.\n\n'
            rtxt += '''In each case, three plots are made (assuming both BUSCO and BUSCOMP data is available):
            
1. BUSCO ratings for missing BUSCO genes.
2. BUSCOMP ratings for missing BUSCO genes. As well as being more relaxed than pure BUSCO results, this will indicate
when BUSCOMP has found a gene in the focal assembly/group where BUSCO did not.
3. BUSCOMP ratings for missing BUSCOMP genes. It is expected that assemblies will be much more similar in terms of BUSCOMP
coverage.

'''
            if fullreport:
                rcode = missingPlot + missingSeqPlot
                if self.db('busco'):
                    rcode += 'buscofull <- read.delim("%s.busco.tdt", header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)\n' % self.baseFile()
                if self.db('buscomp'):
                    rcode += 'buscomp <- read.delim("%s.buscomp.tdt", header = TRUE, stringsAsFactors = FALSE, comment.char = "", fill = TRUE)\n' % idbase
                rtxt += '```{r missingfunctions, echo=FALSE}\n%s\n```\n\n' % (rcode)
                for gentry in gdb.entries(sorted=True):
                    genome = gentry['Genome']
                    if self.db('busco'):
                        rtxt += '## Missing %s BUSCO genes\n\n' % genome
                        if gentry['Directory'] or genome in self.dict['Groups']:
                            rtxt += 'BUSCO ratings for `Missing` %s BUSCO genes:\n\n' % genome
                            genome = self.rGenomeField(genome)
                            rcode = 'missingPlot(buscofull,"%s")' % genome
                            codename = 'missing.%s' % genome
                            rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(self.db('genomes').entryNum()/2),rcode)
                        if self.db('buscomp'):
                            if gentry['Directory'] or genome in self.dict['Groups']:
                                rtxt += 'BUSCOMP ratings for `Missing` %s BUSCO genes:\n\n' % genome
                                genome = self.rGenomeField(genome)
                                rcode = 'missingSeqPlot(buscomp,buscofull,"%s")' % genome
                                codename = 'nobuscoseq.%s' % genome
                                rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(self.db('genomes').entryNum()/2),rcode)
                    if self.db('buscomp'):
                        if gentry['Fasta'] or genome in self.dict['Groups']:
                            rtxt += 'BUSCOMP ratings for `Missing` %s BUSCOMP genes:\n\n' % genome
                            genome = self.rGenomeField(genome)
                            rcode = 'missingSeqPlot(buscomp,buscomp,"%s")' % genome
                            codename = 'missingseq.%s' % genome
                            rtxt += '```{r %s, echo=FALSE, fig.width=12, fig.height=%d}\n%s\n```\n\n' % (codename,2+int(self.db('genomes').entryNum()/2),rcode)
            else:
                rtxt += 'Plots for `Missing` BUSCO genes are only produced for `fullreport=T` output.\n\n'

            ### [X] Errors and warnings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rtxt += '<a name="Appendix" />\n\n# Appendix: BUSCOMP run details\n\n'
            rtxt += '    %s\n\n' % self.log.runDetails()
            rtxt += 'This analysis was run in:\n\n'
            rtxt += '    %s\n\n' %  os.path.abspath(os.curdir)
            rtxt += '* Log file: <a href="%s">`%s`</a>\n' % (self.log.info['LogFile'],self.log.info['LogFile'])
            #rtxt += '## BUSCOMP run details\n\n'
            argcmd = rje.longCmd(sys.argv[1:])
            rtxt += '* Commandline arguments: `%s`\n' % string.join(argcmd,'` `')
            rtxt += '* Full Command List: `%s`\n\n' % string.join(string.split(rje.argString(rje.tidyArgs(self.cmd_list))),'` `')
            rtxt += '<a name="Errors" />\n\n'
            rtxt += '## BUSCOMP errors\n\n'
            errors = self.log.list['ErrorLog'][0:]     # List of log error messages.
            if errors:
                #rtxt += 'See run log for further details:\n' + string.join(['']+errors,'\n* ') + '\n\n'
                rtxt += 'See run log for further details:\n\n```\n' + string.join(errors,'\n') + '\n```\n\n'
            else: rtxt += 'BUSCOMP returned no runtime errors.\n\n'

            rtxt += '## BUSCOMP warnings\n\n'
            warnings = self.log.list['WarnLog'][0:]    # List of log warning messages. Stored when silent=True or memsaver=False.
            if warnings:
                #!# Replace with '* `%s`\n' % warning
                #rtxt += 'See run log for further details:\n' + string.join(['']+warnings,'\n* \\') + '\n\n'
                rtxt += 'See run log for further details:\n\n```\n' + string.join(warnings,'\n') + '\n```\n\n'
            else: rtxt += 'BUSCOMP returned no runtime warnings.\n\n'

            ### [X] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rtxt += toc
            rtxt += '\n<small>Output generated by BUSCOMP v%s &copy; 2019 Richard Edwards | richard.edwards@unsw.edu.au</small>\n' % info.version
            self.deBug(rtxt)
            ## [Xa] Output and convert rtxt to html ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.baseFile(idbase)
            rmdfile = '%s.Rmd' % rbase
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown BUSCOMP summary output to %s' % rmdfile)
            if self.getBool('RmdReport'):
                rmd.rmdKnit(rmdfile,stdout=self.debugging())
            else: self.printLog('#RMD','Set rmdreport=T to knit %s to HTML' % rmdfile)

        except: self.errorLog('Problem during %s buscompRmd.' % self.prog()); return False  # Setup failed
#########################################################################################################################
### End of SECTION II: BUSCOMP Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: BUSCOMP R CODE                                                                                         #
#########################################################################################################################
def ggPlotCode(x,y,title,xlab,ylab,xmax=0,ymax=0,vline=None,xlog=False,ylog=False):

    plotcode = 'ggdata = pdata\n'
    #plotcode += 'numVsizeGG <- ggplot(ggdata, aes(x=%s, y=%s, colour=col, label=label, size=1)) +' % (x,y)
    plotcode += 'ggplot(ggdata, aes(x=%s, y=%s, colour=col, label=label, size=1)) +' % (x,y)
    plotcode += '    geom_point(size = 3,shape=ggdata$pch, colour=ggdata$col, na.rm=TRUE) +'
    plotcode += '    theme_bw() +'
    plotcode += '    theme(text = element_text(size=15)) +'
    plotcode += '    theme(plot.title = element_text(hjust = 0.5)) +'
    plotcode += '    theme(panel.border = element_rect(colour="black", size=0.5)) +'
    plotcode += '    ggtitle("%s") +' % title
    plotcode += '    guides(color = "none") +'
    plotcode += '    scale_shape(solid = FALSE) +'
    plotcode += '    xlab("%s") +' % xlab
    plotcode += '    ylab("%s") +' % ylab
    plotcode += '    xlim(0,max(ggdata$%s,%s)) +' % (x,xmax)
    if xlog: plotcode += '    coord_trans(x="log10") + '
    if ylog: plotcode += '    coord_trans(y="log10") + '
    plotcode += '    ylim(0,max(ggdata$%s,%s)) +' % (y,ymax)
    plotcode += '    geom_label_repel(aes(label=label), size = 5, colour=ggdata$col, na.rm=TRUE) +'
    if vline:
        plotcode += '    geom_vline(xintercept = %s, color = "grey", linetype="dashed", size=0.5) +' % vline

    plotcode = plotcode[:-2]
    plotcode = plotcode.replace(' +',' +\n')

    return plotcode


#########################################################################################################################
buscoPlot = '''# BUSCO Plot of summary data table (Dataset, Complete, Duplicated, Fragmented, Missing)
buscoPlot = function(sumdata,title=""){
  sumdata = sumdata[nrow(sumdata):1,]
  rownames(sumdata) = sumdata$Dataset
  N = sum(sumdata[1,2:5])
  my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
  par(mar=c(5,15,4,0)+0.1)
  barplot(t(as.matrix(sumdata[2:5])),horiz=TRUE,legend=TRUE,axes=TRUE,las=1,args.legend=c(x="topleft"),col=my_colors,main=title,xlab=paste0("BUSCO Count (n=",N,")"),xlim=c(0,N))
  par(mar=c(5,4,4,2)+0.1)
}
'''

buscoPercPlot = '''# BUSCO Plot of summary data table as % (Dataset, Complete, Duplicated, Fragmented, Missing)
buscoPercPlot = function(sumdata,title="",maketext=TRUE){
  sumdata = sumdata[nrow(sumdata):1,]
  if(maketext){
    sumdata$text = ""
    for(i in 1:nrow(sumdata)){
      sumdata$text[i] = paste0("C:",sumdata$Complete[i]+sumdata$Duplicated[i]," [S:",sumdata$Complete[i],", D:",sumdata$Duplicated[i],"], F:",sumdata$Fragmented[i],", M:",sumdata$Missing[i],", n:",sum(sumdata[i,2:5]))
    }
  }
  rownames(sumdata) = sumdata$Dataset
  N = sum(sumdata[1,2:5])
  my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F04442")
  par(mar=c(5,12,4,1)+0.1)
  barplot(t(as.matrix(sumdata[2:5]))*100/N,horiz=TRUE,legend=TRUE,axes=TRUE,las=1,args.legend=c(x="topright"),col=my_colors,main=title,xlab=paste0("BUSCO Percentage (n=",N,")"),xlim=c(0,119))
  text(0,0.6+0:(dim(sumdata)[1]-1)*1.2,sumdata$text,pos=4)
  par(mar=c(5,4,4,2)+0.1)
}
'''

#!# Need to pick better colours for Partial and Ghost
buscompSeqPercPlot = '''# BUSCO Plot of summary data table as % (Genome, Single, Duplicated, Fragmented, Partial, Ghost, Missing)
buscompSeqPercPlot = function(sumdata,title="",maketext=TRUE){
  sumdata = sumdata[nrow(sumdata):1,]
  sumdata = sumdata[! is.na(sumdata$Genome),]
  rownames(sumdata) = sumdata$Genome
  if(maketext){
    sumdata$text = ""
    for(i in 1:nrow(sumdata)){
      sumdata$text[i] = paste0("C:",sumdata$Single[i]+sumdata$Duplicated[i]," [S:",sumdata$Single[i],", D:",sumdata$Duplicated[i],"], F+P:",sumdata$Fragmented[i]+sumdata$Partial[i],", G+M:",sumdata$Ghost[i]+sumdata$Missing[i],", n:",sum(sumdata[i,2:7]))
    }
  }
  N = sum(sumdata[1,2:7])
  my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F0E442", "#F04442", "#F04442")
  par(mar=c(5,12,4,1)+0.1)
  barplot(t(as.matrix(sumdata[2:7]))*100/N,horiz=TRUE,legend=TRUE,axes=TRUE,las=1,args.legend=c(x="topright"),col=my_colors,main=title,xlab=paste0("BUSCOMP Percentage (n=",N,")"),xlim=c(0,119))
  text(0,0.6+0:(dim(sumdata)[1]-1)*1.2,sumdata$text,pos=4)
  par(mar=c(5,4,4,2)+0.1)
}

buscompSeqPercPlotNA = function(sumdata,title="",maketext=TRUE){
  sumdata = sumdata[nrow(sumdata):1,]
  rownames(sumdata) = sumdata$Genome
  if(maketext){
    sumdata$text = ""
    for(i in 1:nrow(sumdata)){
      sumdata$text[i] = paste0("C:",sumdata$Complete[i]+sumdata$Duplicated[i]," [S:",sumdata$Complete[i],", D:",sumdata$Duplicated[i],"], F+P:",sumdata$Fragmented[i]+sumdata$Partial[i],", G+M:",sumdata$Ghost[i]+sumdata$Missing[i],", n:",sum(sumdata[i,2:7]))
    }
  }
  N = sum(sumdata[1,2:8])
  my_colors <- c("#56B4E9", "#3492C7", "#F0E442", "#F0E442", "#F04442", "#F04442", "white")
  par(mar=c(5,12,4,1)+0.1)
  barplot(t(as.matrix(sumdata[2:8]))*100/N,horiz=TRUE,legend=TRUE,axes=TRUE,las=1,args.legend=c(x="topright"),col=my_colors,main=title,xlab=paste0("BUSCOMP Percentage (n=",N,")"),xlim=c(0,119))
  text(0,0.6+0:(dim(sumdata)[1]-1)*1.2,sumdata$text,pos=4)
  par(mar=c(5,4,4,2)+0.1)
}
'''

missingPlot = '''# Function to plot status of missing BUSCOs
missingPlot = function(plotdf,genome){
    missdf = plotdf[plotdf[[genome]]=="Missing",]
    df = data.frame(Dataset=colnames(missdf)[-1], Complete = 0, Duplicated = 0, Fragmented = 0, Missing = 0)
    rownames(df) = colnames(missdf)[-1] 
    for(gen in colnames(missdf)[-1]){
        missdf[,gen] = factor(missdf[,gen],levels=colnames(df)[-1])
        #levels(missdf[,gen]) = colnames(df)[-1]
        gtab = table(missdf[,gen])
        for(field in colnames(df)[-1]){
            df[gen,field] = gtab[[field]]
        }
    }
    buscoPercPlot(df,title=paste("Missing",genome,"BUSCOs"))
}
'''

missingSeqPlot = '''# Function to plot status of missing BUSCOs
missingSeqPlot = function(plotdf,ratedf,genome){
    misseog = ratedf[ratedf[[genome]]=="Missing" | ratedf[[genome]]=="NULL" | is.na(ratedf[[genome]]),]$BuscoID
    missdf = plotdf[plotdf$BuscoID %in% misseog,]
    df = data.frame(Genome=colnames(missdf)[-c(1:2)], Complete = 0, Duplicated = 0, Fragmented = 0, Partial = 0, Ghost = 0, Missing = 0, NULL = 0)
    rownames(df) = colnames(missdf)[-c(1:2)] 
    for(gen in colnames(missdf)[-c(1:2)]){
        missdf[,gen] = factor(missdf[,gen],levels=colnames(df)[-1])
        #levels(missdf[,gen]) = colnames(df)[-1]
        gtab = table(missdf[,gen])
        for(field in colnames(df)[-1]){
            df[gen,field] = gtab[[field]]
        }
    }
    if(identical(plotdf,ratedf)){ 
        ptitle=paste("Missing",genome,"BUSCOMPs") 
    }else{ 
        ptitle=paste("Missing",genome,"BUSCOs: BUSCOMP ratings")
    } 
    buscompSeqPercPlotNA(df,title=ptitle)
}
'''


genomeStats = '''## Genome statistics

The following genome statistics were also calculated by `RJE_SeqList` for each genome (table, below):

* **SeqNum**: The total number of scaffolds/contigs in the assembly.
* **TotLength**: The total combined length of scaffolds/contigs in the assembly.
* **MinLength**: The length of the shortest scaffold/contig in the assembly.
* **MaxLength**: The length of the longest scaffold/contig in the assembly.
* **MeanLength**: The mean length of scaffolds/contigs in the assembly.
* **MedLength**: The median length of scaffolds/contigs in the assembly.
* **N50Length**: At least half of the assembly is contained on scaffolds/contigs of this length or greater.
* **L50Count**: The smallest number scaffolds/contigs needed to cover half the the assembly.
* **CtgNum**: Number of contigs (`SeqNum`+`GapCount`).
* **N50Ctg**: At least half of the assembly is contained on contigs of this length or greater.
* **L50Ctg**: The smallest number contigs needed to cover half the the assembly.
* **NG50Length**: At least half of the genome is contained on scaffolds/contigs of this length or greater.
This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
* **LG50Count**: The smallest number scaffolds/contigs needed to cover half the the genome.
This is based on `genomesize=X`. If no genome size is given, it will be relative to the biggest assembly.
* **GapLength**: The total number of undefined "gap" (`N`) nucleotides in the assembly.
* **GapCount**: The total number of undefined "gap" (`N`) regions in the assembly.
* **GC**: The %GC content of the assembly.

'''

bestText = '''There are many ways of assessing genome assembly quality, but they can be broadly broken down into four aspects:

1. **Coverage.** How much of the genome is included in the assembly.

2. **Contiguity.** How many fragments is the genome in, and how big are they.

3. **Accuracy.** How accurate is the assembly in terms of sequence errors.

4. **Redundancy.** How much of the genome has been included multiple times.

Standard reference-free genome statistics (e.g. number, length and gappiness of scaffolds), can only give limited
insights. Whilst assemblies smaller than the genome size are clearly missing bits, assembly size could be
artificially inflated by redundant regions, making it impossible to assess Coverage. Scaffold sizes can give an
indicator of Contiguity, but are also prone to biases introduced by filtering of small scaffolds, for example.
If an estimated Genome Size can be provided, the most useful measures in this respect are `NG50` and `LG50`.
These statistics sort scaffolds by length (big to small) and identify the contig/scaffold size at which half the
*genome* (not the *assembly*) is covered by contigs/scaffolds of that size or bigger. This is the `NG50` value
(called `NG50Length` in BUSCOMP), with `LG50` (`LG50Count`) being the number of contigs/scaffolds needed to
cover half the genome. (`N50` and `L50` are the same statistics only for the assembly.) These statistics can
still be mislead by redundancy in the assembly. None of theses statistics speak to sequence accuracy.

The power of BUSCO is that it speaks to all four aspects of quality. `Complete` and `Fragmented` genes give
an indication of Coverage, Continuity and Accuracy; `Duplicated` genes give an indication of Redundancy;
`Missing` genes give an indication of Coverage and Accuracy. However, the weakness is that these different
aspects cannot easily be disentangled. This makes side-by-side comparisons of different assemblies challenging,
as it is not always clear what a difference is indicating.

BUSCOMP is designed on the principle that **Coverage** and **Contiguity** are the two most important aspects of
assembly quality. *Accuracy* can, in principle, be improved by additional error-correction steps (including
manual curation). Suspected *Redundancy* can likewise be identified a flagged. *Coverage* and *Contiguity*, in
contrast, can only be improved by further assembly - either trying again from scratch, or employing a
scaffolding or gap-filling alogrithm.

BUSCOMP has therefore identified seven statistics that can be used to rank assemblies on inferred Completeness
or Contiguity:

* **Completeness**:
    1. `Complete` = Maximum number of Complete BUSCOMP ratings.
    2. `BUSCO` = Maximum number of Complete BUSCO ratings.
    3. `Missing` = Smallest number of Missing BUSCOMP ratings.
    4. `NoBUSCO` = Smallest number of Missing BUSCO ratings.
* ** Contiguity**:
    1. `MaxLength` = Maximum length of contigs/scaffolds.
    2. `NG50Length` = Half the genome lies on contigs/scaffolds at least this size. (See above.)
    3. `LG50Count = Half the genome can be covered by this many contigs/scaffolds. (See above.)

Individual assemblies are rated as "best" under all seven criteria, with ties allowed. Each assembly can
be best under multiple criteria and each criterion can have several best assemblies. BUSCOMP will report all
such combinations.
'''

buscoText = '''BUSCO ratings are defined (quoting from the [BUSCO v3 User Guide](http://gitlab.com/ezlab/busco/raw/master/BUSCO_v3_userguide.pdf) as:

* `Complete`: Single-copy hits where "BUSCO matches have scored within the expected range of scores and within the expected range of length alignments to the
BUSCO profile."
* `Duplicated`: As `Complete` but 2+ copies.
* `Fragmented`: "BUSCO matches ... within the range of scores but not within the range of length alignments to the BUSCO profile."
* `Missing`: "Either no significant matches at all, or the BUSCO matches scored below the range of scores for the BUSCO profile."
'''
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
    try: BUSCOMP(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
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
