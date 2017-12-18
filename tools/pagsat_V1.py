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
# Author contact: <seqsuite@gmail.com> / School of Biotechnology and Biomolecular Sciences, UNSW, Sydney, Australia.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       PAGSAT
Description:  Pairwise Assembled Genome Sequence Analysis Tool
Version:      1.12.0
Last Edit:    26/09/16
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for the assessment of an assembled genome versus a suitable reference. For optimal results, the
    reference genome will be close to identical to that which should be assembled. However, comparative analyses should
    still be useful when different assemblies are run against a related genome - although there will not be the same
    expectation for 100% coverage and accuracy, inaccuracies would still be expected to make an assembly less similar
    to the reference.

    Main input for PAGSAT is an assembled genome in fasta format (`assembly=FILE`) and a reference genome in fasta format
    (`refgenome=FILE` or `reference=FILE`) with corresponding `*.gb` or `*.gbk` genbank download for feature extraction.

Output:
    Main output is a number of delimited text files and PNG graphics made with R. Details to follow.

Commandline:
    ### ~ Input/Setup Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    assembly=FILE   : Fasta file of assembled contigs to assess [None]
    refgenome=FILE  : Fasta file of reference genome for assessment (also *.gb for full functionality) [None]
    spcode=X        : Species code for reference genome (if not already processed by rje_genbank) [None]
    minqv=X         : Minimum mean QV score for assembly contigs (read from *.qv.csv) [20]
    mincontiglen=X  : Minimum contig length to retain in assembly (QV filtering only) [1000]
    casefilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    ### ~ Reference vs Assembly Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minlocid=X      : Minimum percentage identity for local hits mapping to chromosome coverage [95.0]
    minloclen=X     : Mininum length for local hits mapping to chromosome coverage [250]
    genesummary=T/F : Whether to include reference gene searches in summary data [True]
    protsummary=T/F : Whether to include reference protein searches in summary data [True]
    tophitbuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]
    diploid=T/F     : Whether to treat assembly as a diploid [False]
    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X      : Basename for output files and directories. [assembly+ref]
    chromalign=T/F  : Whether to perform crude chromosome-contig alignment [True]
    rgraphics=T/F   : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
    dotplots=T/F    : Whether to use gablam.r to output dotplots for all ref vs assembly. [False]
    report=T/F      : Whether to generate HTML report [True]
    genetar=T/F     : Whether to tar and zip the GeneHits/ and ProtHits/ folders (if generated & Mac/Linux) [True]
    ### ~ Comparison Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    compare=FILES   : Compare assemblies selected using a list of *.Summary.tdt files (wildcards allowed). []
    fragcov=LIST    : List of coverage thresholds to count min. local BLAST hits (checks integrity) [50,90,95,99]
    chromcov=LIST   : Report no. of chromosomes covered by a single contig at different %globID (GABLAM table) [95,98,99]
    ### ~ Assembly Tidy/Edit Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    tidy=T/F        : Execute semi-automated assembly tidy/edit mode to complete draft assembly [False]
    newacc=X        : New base for edited contig accession numbers (None will keep old accnum) [None]
    newchr=X        : Code to replace "chr" in new sequence names for additional PAGSAT compatibility [ctg]
    orphans=T/F     : Whether to include and process orphan contigs [True]
    chrmap=X        : Contig:Chromosome mapping mode for assembly tidy (unique/align) [unique]
    joinsort=X      : Whether to sort potential chromosome joins by `Length` or `Identity` [Identity]
    joinmerge=X     : Merging mode for joining chromosomes (consensus/end) [end]
    joinmargin=X    : Number of extra bases allowed to still be considered an end local BLAST hit [10]
    snapper=T/F     : Run Snapper on ctidX/haploid output following PAGSAT Tidy. (Re-Quiver recommended first.) [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_genbank, rje_html, rje_obj, rje_seqlist, rje_sequence, rje_synteny, rje_tree, rje_tree_group, rje_xref
import rje_blast_V2 as rje_blast
import rje_dismatrix_V3 as rje_dismatrix
import gablam, snapper
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0.0 - Initial working version for based on rje_pacbio assessment=T.
    # 1.1.0 - Fixed bug with gene and protein summary data. Removed gene/protein reciprocal searches. Added compare mode.
    # 1.1.1 - Added PAGSAT output directory for tidiness!
    # 1.1.2 - Renamed the PacBio class PAGSAT.
    # 1.2.0 - Tidied up output directories. Added QV filter and Top Gene/Protein hits output.
    # 1.2.1 - Added casefilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    # 1.3.0 - Added tophitbuffer=X and initial synteny analysis for keeping best reference hits.
    # 1.4.0 - Added chrom-v-contig alignment files along with *.ordered.fas.
    # 1.4.1 - Made default chromalign=T.
    # 1.4.2 - Fixed casefilter=F.
    # 1.5.0 - diploid=T/F     : Whether to treat assembly as a diploid [False]
    # 1.6.0 - mincontiglen=X  : Minimum contig length to retain in assembly [1000]
    # 1.6.1 - Added diploid=T/F to R PNG call.
    # 1.7.0 - Added tidy=T/F option. (Development)
    # 1.7.1 - Updated tidy=T/F to include initial assembly.
    # 1.7.2 - Fixed some bugs introduced by changing gablam fragment output.
    # 1.7.3 - Added circularise sequence generation.
    # 1.8.0 - Added orphan processing and non-chr naming of Reference.
    # 1.9.0 - Modified the join sorting and merging. Added better tracking of positions when trimming.
    # 1.9.1 - Added joinmargin=X    : Number of extra bases allowed to still be considered an end local BLAST hit [10]
    # 1.10.0 - Added weighted tree output and removed report warning.
    # 1.10.1 - Fixed issue related to having Description in GABLAM HitSum tables.
    # 1.10.2 - Tweaked haploid core output.
    # 1.10.3 - Fixed tidy bug for RevComp contigs and switched joinsort default to Identity. (Needs testing.)
    # 1.10.4 - Added genetar option to tidy out genesummary and protsummary output. Incorporated rje_synteny.
    # 1.10.5 - Set gablamfrag=1 for gene/protein hits.
    # 1.11.0 - Consolidated automated tidy mode and cleaned up some excess code.
    # 1.11.1 - Added option for running self-PAGSAT of ctidX contigs versus haploid set. Replaced ctid "X" with "N".
    # 1.11.2 - Fixed Snapper run choice bug.
    # 1.11.3 - Added reference=FILE as alias for refgenome=FILE. Fixed orphan delete bug.
    # 1.12.0 - Tidying up and documenting outputs. Changed default minloclen=250 and minlocid=95. (LTR identification.)
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add reduced functionality if genbank file not given (i.e. no features).
    # [Y] : Add HTML report output.
    # [ ] : Add indels to QAssemble
    # [?] : Have a min minloclen (=localcut) but then try increasing by 1kb chunks without increasing "N" periods?
    # [Y] : Calculate the difference from Reference as a number for comparisons.
    # [?] : Pull out "N" regions of the reference - QAssemble back against the pre-assembly and subreads.
    # [ ] : Add (interactive?) reformatting of *.full.fas for refgenome input.
    # [Y] : Include minloclen in relevant file names - GABLAM and default basefile.
    # [Y] : Consider using minloclen for localcut=X GABLAM Cut-off length for local alignments contributing to global stats.
    # [Y] : Option to switch off Gene and Protein searches for increased speed. (Need to edit R summGraph() too)
    # [ ] : Improve gene search to extend hits to full length of genes/proteins.
    # [ ] : Separate minloclen and localcut?
    # [ ] : Add thumbnails=T/F for report and R graphics
    # [ ] : Add maxcontig=X for a single summary page. (Ask to proceed if more contigs? Or split?)
    # [ ] : Add summary stats for assembly and reference? (Could load them into SeqList objects)
    # [Y] : Add reading of unitig coverage and quality scores from quiver (*.qv.csv files)
    # [Y] : Need to tidy up outputs. Generate a *.PAGSAT/ directory for most outputs: update R script accordingly.
    # [N] : blastdir=PATH   : Path for blast results file (unless keepblast=F) [./assembly.BLAST/]
    # [Y] : Add separate GABLAM directory for feature searches? (./assembly.GABLAM/) (These ignore cutoffs, right?)
    # [Y] : Add (and distinguish) minlocid to output file names as well as minloclen. (Make integer? *.LXXX.IDXX.*)
    # [X] : Contemplate setting softmask=F for BLAST searches. (Why are some genes missing?!)
    # [ ] : Consider replacing GABLAM fragfas with own gene extraction algorithm that extends ORFs & combines exons.
    # [Y] : Move R graphics such that it can be re-run in isolation (if summary.png missing).
    # [ ] : Need to add Ty element searching and mapping. Probably need to do this prior to assembly?
    # [ ] : Add Ty element map to the contig output somehow? Or just a separate map?
    # [Y] : Sort out Full, Tandem, Partial etc. TopHits classification.
    # [ ] : Add an "ends" mode for tidying where only unique-mapped ends are considered for joining.
    # [ ] : Will need to modify contig chr naming for "ends" mode - some ambiguous, some not.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('PAGSAT', '1.12.0', 'September 2016', '2015')
    description = 'Pairwise Assembled Genome Sequence Analysis Tool'
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
def setupProgram(extracmd=[]): ### Basic Setup of Program when called from commandline.
    '''
    Basic Setup of Program when called from commandline:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:### ~ [1] ~ Initial Command Setup & Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        info = makeInfo()                                   # Sets up Info object with program details
        cmd_list = rje.getCmdList(sys.argv[1:]+extracmd,info=info)   # Reads arguments and load defaults from program.ini
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
### SECTION II: PAGSAT Class                                                                                            #
#########################################################################################################################
class PAGSAT(rje_obj.RJE_Object):
    '''
    PAGSAT Class. Author: Rich Edwards (2015).

    Str:str
    - Assembly=FILE   : Fasta file of assembled contigs to assess [None]
    - BaseBase        : Path-trimmed Basefile for outputs that do not involve data filtering etc.
    - ChrMap=X        : Contig:Chromosome mapping mode for assembly tidy (unique/align) [unique]
    - CutBase         : Path-trimmed Basefile for outputs that include data filtering etc.
    - GABLAMDir       : Parent directory for all BLAST and GABLAM searches.
    - JoinMerge=X     : Merging mode for joining chromosomes (consensus/long/end) [end]
    - JoinSort=X      : Whether to sort potential chromosome joins by `Length` or `Identity` [Length]
    - NewAcc=X        : New base for edited contig accession numbers (None will keep old accnum) [None]
    - NewChr=X        : Code to replace "chr" in new sequence names for additional PAGSAT compatibility [ctg]
    - RefBase=X       : Basefile for reference genome for assessment (*.gb) [None]
    - RefGenome=FILE  : Fasta file of reference genome for assessment (also *.gb for full functionality) [None]
    - ResDir          : Results directory = BASEFILE.PAGSAT/

    Bool:boolean
    - CaseFilter=T/F  : Whether to filter leading/trailing lower case (low QV) sequences [True]
    - ChromAlign=T/F  : Whether to align chromosomes with contigs (slow!) [False]
    - Diploid=T/F     : Whether to treat assembly as a diploid [False]
    - DotPlots=T/F    : Whether to use gablam.r to output dotplots for all ref vs assembly. [False]
    - Features=T/F    : Whether to expect a Features table (from Genbank processing) [True]
    - GeneSummary=T/F : Whether to include reference gene searches in summary data [True]
    - GeneTar=T/F     : Whether to tar and zip the GeneHits/ and ProtHits/ folders (if generated & Mac/Linux) [True]
    - Orphans=T/F     : Whether to include and process orphan contigs [True]
    - ProtSummary=T/F : Whether to include reference protein searches in summary data [True]
    - RGraphics=T/F   : Whether to generate PNG graphics using R [True]
    - Report=T/F      : Whether to generate HTML report [True]
    - Snapper=T/F     : Run Snapper on ctidX output following PAGSAT Tidy [False]
    - Tidy=T/F        : Execute semi-automated assembly tidy/edit mode to complete draft assembly [False]

    Int:integer
    - JoinMargin=X    : Number of extra bases allowed to still be considered an end local BLAST hit [10]
    - MinContigLen=X  : Minimum contig length to retain in assembly [1000]
    - MinLocLen=X     : Mininum length for local hits mapping to chromosome coverage [100]
    - MinQV=X         : Minimum mean QV score for assembly contigs (read from *.qv.csv) [20]

    Num:float
    - MinLocID=X      : Minimum percentage identity for local hits mapping to chromosome coverage [0.99]
    - TopHitBuffer=X  : Percentage identity difference to keep best hits for reference genes/proteins. [1.0]

    File:file handles with matching str filenames
    
    List:list
    - ChromCov=LIST   : Report no. of chromosomes covered by a single contig at different %globID (GABLAM table) [95,98,99]
    - Compare=FILES   : Special mode to compare a list of *.Summary.tdt files (wildcards allowed). []
    - FragCov=LIST    : List of coverage thresholds to count min. local BLAST hits (checks integrity) [50,90,95,99]

    Dict:dictionary
    - QV              : Dictionary of contig:QV score read from *.qv.csv

    Obj:RJE_Objects
    - Features  : Reference Features Database table (reused for Compare=T).
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Assembly','BaseBase','CaseFilter','ChrMap','CutBase','GABLAMDir','JoinMerge','JoinSort','NewAcc','NewChr','RefBase','RefGenome','ResDir']
        self.boollist = ['ChromAlign','Diploid','DotPlots','Features','GeneSummary','GeneTar','Orphans','ProtSummary','RGraphics','Report','Snapper','Tidy']
        self.intlist = ['JoinMargin','MinContigLen','MinLocLen','MinQV']
        self.numlist = ['MinLocID','TopHitBuffer']
        self.filelist = []
        self.listlist = ['ChromCov','Compare','FragCov']
        self.dictlist = ['QV']
        self.objlist = ['DB','Features']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'ChrMap':'unique','GABLAMDir':rje.makePath('GABLAM/'),'JoinMerge':'end','JoinSort':'Identity','NewChr':'ctg','ResDir':rje.makePath('PAGSAT/')})
        self.setBool({'CaseFilter':True,'ChromAlign':True,'Diploid':False,'DotPlots':False,'Features':True,'Orphans':True,
                      'GeneSummary':True,'GeneTar':True,'ProtSummary':True,'RGraphics':True,'Report':True,'Snapper':False,'Tidy':False})
        self.setInt({'JoinMargin':10,'MinLocLen':250,'MinQV':20,'MinContigLen':1000})
        self.setNum({'MinLocID':95.0,'TopHitBuffer':1.0})
        self.list['ChromCov'] = [95,98,99]
        self.list['FragCov'] = [50,90,95,99]
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'str',['ChrMap','JoinMerge','JoinSort','NewAcc','NewChr'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Assembly','RefGenome'])  # String representing file path
                self._cmdRead(cmd,type='file',att='RefGenome',arg='reference')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['CaseFilter','ChromAlign','Diploid','DotPlots','GeneSummary','GeneTar',
                                              'Orphans','ProtSummary','RGraphics','Report','Snapper','Tidy'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['JoinMargin','MinContigLen','MinLocLen','MinQV'])   # Integers
                self._cmdReadList(cmd,'float',['TopHitBuffer']) # Floats
                self._cmdReadList(cmd,'perc',['MinLocID'])
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',[])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['ChromCov','FragCov'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Compare']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.win32() and self.getBool('GeneTar'):
            self.printLog('#WIN32','Cannot use targz on Windows: GeneTar=Fals.')
            self.setBool({'GeneTar':False})
        if self.getBool('Tidy') and (self.getBool('GeneSummary') or self.getBool('ProtSummary')):
            self.printLog('#TIDY','PAGAT Tidy mode selected with: GeneSummary:%s; ProtSummary:%s' % (self.getBool('GeneSummary'),self.getBool('ProtSummary')))
            self.printLog('#TIDY','Switching genesummary=F protsummary=F will make Tidy Quicker.')
            if self.yesNo('Switch off GeneSummary and ProtSummary for faster PASGAT Tidy run? (Will Reactivate for follow-up PAGSAT on tidied assembly.)'):
                self.setBool({'GeneSummary':False,'ProtSummary':False})
                self.printLog('#TIDY','Set: GeneSummary=%s; ProtSummary=%s' % (self.getBool('GeneSummary'),self.getBool('ProtSummary')))
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','Run')
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['Compare']: return self.compare()
            if self.getBool('Tidy'): return self.tidy()
            elif self.getBool('Report'): return self.report()
            else: return self.assessment()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''
        Main class setup method. This sets up the main Database object.

        In Compare mode, nothing else is done beyond checking/setting basefile.

        Other modes go through a series of setup stages:
        1. Set up Reference Genome.
        2. Set up Assembly with QV filtering if required.
        3. Report on key settings and setup output directories etc.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','Setup')
            ## ~ [0a] Database Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# First, the Database object is setup. This is used to control all output tables and integrate data.
            self.printLog('#~~#','## ~~~~~ PAGSTAT Setup ~~~~~ ##')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            ## ~ [0b] Check for Compare Mode and cancel rest of setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Compare']:
                if not self.basefile(return_none=''): self.basefile('pagsat')
                return True

            ### ~ [1] Set up Reference Genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            checkfiles = []     # List of files that need to exist to proceed
            ## ~ [1a] Set up Reference Genome basefile (RefBase) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# RefBase will include the path to the reference genomes
            if rje.exists('%s.gb' % self.getStr('RefGenome')): self.setStr({'RefBase':self.getStr('RefGenome')})
            elif rje.exists('%s.gbk' % self.getStr('RefGenome')): self.setStr({'RefBase':self.getStr('RefGenome')})
            else: self.setStr({'RefBase':rje.baseFile(self.getStr('RefGenome'))})
            if rje.exists('%s.gbk' % self.getStr('RefBase')): gbfile = '%s.gbk' % self.getStr('RefBase')
            else: gbfile = '%s.gb' % self.getStr('RefBase')
            ## ~ [1b] Establish whether Genbank processing of RefGenome needs to be performed ~~~~~ ##
            rungb = False   # Whether to run rje_genbank on RefGenome
            for rfile in ['full.fas','gene.fas','prot.fas','Feature.tdt']:
                gfile = '%s.%s' % (self.getStr('RefBase'),rfile)
                self.printLog('#CHECK','%s: %s' % (gfile,{True:'Found.',False:'Missing!'}[os.path.exists(gfile)]))
                rungb = rungb or not os.path.exists(gfile)
                checkfiles.append(gfile)
            self.debug('Run Genbank: %s' % rungb)
            #i# NOTE: If *.gb (and *.gbk) file missing, Genbank will not be run but can still use manual features table.
            if not rje.exists(gbfile): rungb = False
            ## ~ [1b] Run Genbank to process sequence data if required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rungb:
                gcmd = ['protacc=locus_tag','details=product,gene_synonym,note,db_xref']   # Defaults
                gcmd += self.cmd_list   # Can over-ride/add. This include spcode=X
                gcmd += ['seqin=%s' % gbfile,'taxdir=','tabout=T','fasout=full,gene,cds,prot']
                rje_genbank.GenBank(self.log,gcmd).run()
                for cfile in checkfiles:
                    if not rje.exists(cfile): raise IOError('Cannot find %s! Genbank processing failure?' % cfile)
            ## ~ [1c] Check for existence of Features file. Could be manually created ~~~~~~~~~~~~~ ##
            ftfile = '%s.Feature.tdt' % self.getStr('RefBase')
            self.setBool({'Features':os.path.exists(ftfile)})       # Whether features table generated/found
            self.debug('Features: %s' % self.getBool('Features'))
            ## ~ [1d] Check sequence naming and offer reformatting options ~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Need to add reformatting here #!#
            if string.split(self.getStr('RefGenome'),'.')[-1] in ['gb','gbk']:  # Look for *.fas
                self.setStr({'RefGenome':'%s.fas' % self.getStr('RefBase')})
            if not self.getStr('RefGenome').endswith('.fas') or not rje.exists(self.getStr('RefGenome')):
                self.printLog('#NAMES','%s.full.fas sequence names will not be suitable.' % self.getStr('RefBase'))
                self.printLog('#NAMES','Please modify gene names in %s.full.fas (e.g. ChrX) and save as %s.fas.' % (self.getStr('RefBase'),self.getStr('RefBase')))
                self.printLog('#NAMES','Then re-run with refgenome=%s.fas.' % self.getStr('RefBase'))
                return False
            self.printLog('#REF','Reference genome fasta: %s' % self.getStr('RefGenome'))
            if not rje.exists(self.getStr('RefGenome')): raise IOError('Cannot find RefGenome: %s!' % self.getStr('RefGenome'))

            ### ~ [2] Assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('Assembly')): raise IOError('Cannot find Assembly: %s!' % self.getStr('Assembly'))
            assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=file'])
            ## ~ [2a] Check assembly gene names. Should be unique. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Should not match Reference either?
            #!# Add reformatting of assembly genes if not acceptable.
            #!# Add checking versus reference genes if required.
            assgenes = []   # These should be unique for later processing
            while assembly.nextSeq():
                agene = assembly.seqGene()
                if agene in assgenes:
                    self.printLog('#NAMES','%s sequence names will not be suitable.' % self.getStr('Assembly'))
                    self.printLog('#NAMES','Please modify gene names in %s to be unique (e.g. CtgX) and re-run.' % (self.getStr('Assembly')))
                    return False
                else: assgenes.append(agene)
            ## ~ [2b] Check for QV data and perform QV filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qvfile = '%s.qv.csv' % rje.baseFile(self.getStr('Assembly'))    # QV data file
            qvcut = self.getInt('MinQV')                                    # QV cutoff
            abase = rje.baseFile(self.getStr('Assembly'),strip_path=True)   # Assembly basefile
            qvbase = '%s.QV/%s.qv%d' % (abase,abase,qvcut)                  # QV-filtered assembly basefile
            qvfas = '%s.fas' % (qvbase)                                     # QV-filtered assembly fasta file
            if rje.exists(qvfile) and not self.force() and rje.exists(qvfas): self.setStr({'Assembly':qvfas})
            if rje.exists(qvfile) and (self.force() or not rje.exists(qvfas)):
                qdb = self.db().addTable(qvfile,['contig_id'],name='QV')
                qdb.dataFormat({'mean_coverage':'num','mean_qv':'num'})
                qdb.addField('Contig'); qdb.addField('Seq'); qdb.addField('Len')
                assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=file','usecase=T'])
                qx = 0
                while assembly.nextSeq():
                    seqname = assembly.seqName()
                    sname = assembly.shortName()
                    for qentry in qdb.entries():
                        if '%s|' % qentry['contig_id'] in seqname:
                            if qentry['Seq']: raise ValueError('Multiple QV sequence mapping error! (%s)' % seqname)
                            qentry['Seq'] = assembly.obj['Current']; qentry['Contig'] = sname
                            qentry['Len'] = assembly.seqLen()
                            if sname in self.dict['QV']: raise ValueError('Multiple QV sequence mapping error! (%s)' % seqname)
                            self.dict['QV'][sname] = qentry['mean_qv']; qx += 1
                if qx != qdb.entryNum(): raise ValueError('Only %d of %d QV values mapped to sequences!' % (qx,qdb.entryNum()))
                if qx != assembly.seqNum(): raise ValueError('%d QV values but %d sequences!' % (qx,assembly.seqNum()))
                # Filter qdb on XCov and/or QV and see if entries lost
                lencut = self.getInt('MinContigLen')
                qvcutlen = 0
                for qentry in qdb.sortedEntries('Contig'):
                    if qentry['mean_qv'] < qvcut:
                        self.printLog('#MINQV','%s failed to meet QV>=%d: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],qvcut,qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                        qvcutlen += qentry['Len']
                        qdb.dropEntry(qentry)
                    elif qentry['Len'] < lencut:
                        self.printLog('#MINLEN','%s failed to meet Len>=%d: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],lencut,qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                        qvcutlen += qentry['Len']
                        qdb.dropEntry(qentry)
                    else:
                        self.printLog('#QV','%s: %.3f kb; XCov=%s; QV=%s' % (qentry['Contig'],qentry['Len']/1000.0,rje.sf(qentry['mean_coverage'],3),rje.sf(qentry['mean_qv'],3)))
                if qx != qdb.entryNum() or self.getBool('CaseFilter'):
                    self.printLog('#QV','%d of %d contigs (%.2f kb) fail to meet QV>=%d' % (qx-qdb.entryNum(),qx,qvcutlen/1000.0,qvcut))
                    rje.mkDir(self,qvbase)
                    if rje.exists(qvfas) and not self.force():
                        self.printLog('#QV','QVFas file found (force=F): %s' % qvfas)
                        QVFAS = None
                    else:
                        # Make new Assembly file in new directory. BASEFILE.Assembly.QVX
                        rje.backup(self,qvfas)
                        QVFAS = open(qvfas,'w')
                    sx = 0; trimx = 0
                    for seq in qdb.indexKeys('Seq'):
                        (name,sequence) = assembly.getSeq(seq)
                        if self.getBool('CaseFilter'):
                            i = 0; slen = j = len(sequence)
                            while i < len(sequence) and sequence[i] == sequence[i].lower(): i += 1
                            while j and sequence[j-1] == sequence[j-1].lower(): j -= 1
                            name = name + ' QVTrimmed:%d-%d' % (i+1,j)
                            sequence = sequence[i:j]
                            trimx += slen - len(sequence)
                            if len(sequence) != slen: self.printLog('#TRIM','%s: %s bp QV trimmed.' % (assembly.shortName(seq),rje.iStr(slen-len(sequence))))
                            if len(sequence) < lencut: self.printLog('#TRIM','%s: %s bp QV trimmed => now too short!' % (assembly.shortName(seq),rje.iStr(slen-len(sequence)))); continue
                            if not sequence: self.warnLog('No quality (upper case) sequence for %s!' % name); continue
                        if QVFAS: QVFAS.write('>%s\n%s\n' % (name, sequence)); sx += 1
                    if self.getBool('CaseFilter'):
                        self.printLog('#TRIM','Total assembly QV trimmed: %.3f kb' % (trimx/1000.0))
                    if QVFAS:
                        QVFAS.close()
                        self.printLog('\r#OUT','%s of %s sequences output to %s.' % (sx,assembly.seqNum(),qvfas))
                        # Save qdb too!
                        qdb.saveToFile('%s.qv.csv' % qvbase)
                    self.setStr({'Assembly':qvfas})
                else: self.printLog('#QV','%d of %d contigs meet QV>=%d' % (qdb.entryNum(),qx,qvcut))
            else: self.warnLog('%s not found: no QV filtering/reporting/' % qvfile)

            ### ~ [3] Report on Key settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#INPUT','Assembly file: %s' % self.getStr('Assembly'))
            self.printLog('#INPUT','Reference genome file: %s' % self.getStr('RefGenome'))
            self.printLog('#INPUT','Reference genbank file: %s.gb' % self.getStr('RefBase'))
            minloclen = self.getInt('MinLocLen')
            self.printLog('#PARAM','Min. local BLAST alignment length: %sbp ("L%d")' % (rje.iStr(minloclen),minloclen))
            minlocid = self.getNum('MinLocID')
            if minlocid > 99: minlocid = minlocid * 10
            self.printLog('#PARAM','Min. local BLAST alignment %%identity: %s%% ("ID%d")' % (rje.sf(self.getNum('MinLocID'),3),int(minlocid)))
            if not self.baseFile(return_none=None):
                self.baseFile('%s.%s' % (rje.baseFile(self.getStr('Assembly'),strip_path=True),rje.baseFile(self.getStr('RefGenome'),strip_path=True)))
            basedir = rje.makePath('%s.PAGSAT/' % self.baseFile())
            gabdir = rje.makePath('%s.GABLAM/' % self.baseFile())
            rje.mkDir(self,basedir)
            rje.mkDir(self,gabdir)
            self.setStr({'GABLAMDir':gabdir,'ResDir':basedir,'BaseBase':os.path.basename(self.baseFile())})
            self.printLog('#GABDIR',self.getStr('GABLAMDir'))   # Directory for GABLAM and BLAST output
            self.printLog('#PAGDIR',self.getStr('ResDir'))      # Directory for PAGSAT output
            # Primary basefile is PAGSAT directory with additional cut-off information added
            self.baseFile('%s%s.L%dID%d' % (basedir,os.path.basename(self.baseFile()),minloclen,int(minlocid)))
            self.setStr({'CutBase':os.path.basename(self.baseFile())})
            self.printLog('#BASE',self.getStr('BaseBase'))      # Root basefile
            self.printLog('#PAGOUT','%s.*' % self.baseFile())   # Primary output basefile
            if self.baseFile() != self.fileBase(): raise ValueError('Dev problem! Please report. Avoid use of basefile=X.')
            self.db().basefile(self.basefile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def fileBase(self,resdir='Res',base='Cut',extra=None):   ### Returns appropriate file output basename
        '''
        Returns appropriate file output basename. Default should return same as self.baseFile()
        @param dir:str ['Res'] = Whether output directory is 'GABLAM' or 'Res'.
        @param base:str ['Cut'] = Whether file basename included cutoffs ('Cut') or not ('Base')
        @param extra:str [None] = Additional filename element to add to basefile.
        @return: path constructed from GABLAMDir/ResDir and CutBase/BaseBase.
        '''
        filebase = '%s%s' % (self.getStr('%sDir' % resdir),self.getStr('%sBase' % base))
        if extra: filebase = '%s.%s' % (filebase,extra)
        return filebase
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        coverage = main results table
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return ['coverage']
#########################################################################################################################
    ### <3> ### PAGSAT GABLAM Methods                                                                                   #
#########################################################################################################################
    def runGABLAM(self,gabcmd,gtype='Reference',blastgz=True): ### Runs GABLAM with given commands and basefile, managing *blast file.
        '''
        Runs GABLAM with given commands and basefile, managing *blast file. Adds minloclen for actual run.
        >> gabcmd:list = List of GABLAM run commands. Will extract basefile from this list.
        >> blastgz:bool [True] = whether to make/process a general BLAST file if possible. (e.g. rename and g(un)zip)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','runGABLAM')
            basefile = None
            for gcmd in gabcmd:
                if gcmd.startswith('basefile='): basefile = string.split(gcmd,'=',1)[1]
            self.printLog('#~~#','## ~~~~~ %s GABLAM ~~~~~ ##' % basefile)
            # Set blastbase = basefile name for blast file
            if gtype == 'Self': blastbase = string.join(string.split(basefile,'.')[:-1],'.')
            else: blastbase = self.fileBase('GABLAM','Base',gtype)
            #if basefile.endswith('.%d' % self.getInt('MinLocLen')): blastbase = string.join(string.split(basefile,'.')[:-1],'.')
            #else: blastbase = basefile
            ### ~ [1] Run GABLAM, processing BLAST file as required ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blastfile = '%s.blast' % blastbase
            # Unzip file if required
            if rje.exists('%s.gz' % blastfile) and blastgz:
                os.system('gunzip %s.gz' % blastfile)
                self.printLog('#GUNZIP','%s unzipped.' % blastfile)
            # Rename BLAST file if found
            if rje.exists(blastfile) and blastbase != basefile:
                os.rename(blastfile,'%s.blast' % basefile)
                self.printLog('#BLAST','%s -> %s.blast' % (blastfile,basefile))
            elif not rje.exists(blastfile): self.printLog('#BLAST','%s not found.' % (blastfile))
            # Run GABLAM
            gablam.GABLAM(self.log,gabcmd).run()
            # Rename BLAST file if kept
            if rje.exists('%s.blast' % basefile) and blastbase != basefile:
                os.rename('%s.blast' % basefile,blastfile)
                self.printLog('#BLAST','%s.blast -> %s' % (basefile,blastfile))
            # (Re)zip BLAST file if required
            if rje.exists(blastfile) and blastgz:
                os.system('gzip %s' % blastfile)
                self.printLog('#GZIP','%s (re)zipped.' % blastfile)
            return True
        except: self.errorLog('%s.runGABLAM() error' % self.prog()); return False
#########################################################################################################################
    def qAssembleGABLAM(self,gtype='Reference'):   ### Generate GABLAM analyses for assembly assessment.
        '''
        Generate GABLAM analyses for assembly assessment.
        >> gtype:str ['qassemble'] = type of GABLAM to run (Reference/Assembly/Genes/Proteins/Genes.Reciprocal/Proteins.Reciprocal)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','qAssembleGABLAM')
            ## Process reference gb file if required and find genome, genes and proteins
            refbase = self.getStr('RefBase')    # Used for searches of genes and proteins
            gdir = self.getStr('GABLAMDir')     # GABLAMDir. Also used as BLAST directory
            if gtype in ['Reference','Assembly']:  # These searches used cutoff data
                gbase = self.fileBase('GABLAM','Cut',gtype)
            else: gbase = self.fileBase('GABLAM','Base',gtype)     # Other searches do not need cutoffs.

            #basefile = self.baseFile()
            #gdir = '%s.GABLAM/' % basefile
            #if gtype not in ['Assembly','Reciprocal'] and basefile.endswith('.%d' % self.getInt('MinLocLen')):
            #    gdir = '%s.GABLAM/' % string.join(string.split(basefile,'.')[:-1],'.')
            #    gbase = '%s%s' % (gdir,string.join(string.split(self.baseFile(strip_path=True),'.')[:-1],'.'))
            #else: gbase = '%s%s' % (gdir,self.baseFile(strip_path=True))

            self.printLog('#GDIR',gdir)
            rje.mkDir(self,gdir,log=True)
            if not self.force():
                runfound = True
                if gtype == 'Reference' and not rje.exists('%s.hitsum.tdt' % gbase): runfound = False
                if gtype == 'Reference' and not rje.exists('%s.gablam.tdt' % gbase): runfound = False
                if gtype == 'Reference' and not rje.exists('%s.local.tdt' % gbase): runfound = False
                if gtype != 'Reference' and not rje.exists('%s.hitsum.tdt' % (gbase)): runfound = False
                self.printLog('#GABLAM','%s GABLAM %s.* found: %s' % (gtype,gbase,runfound))
                if runfound: return True
            ## Might want to set forks for speed up
            if self.getInt('Forks') < 2: self.printLog('#INFO','Consider setting forks=X to speed up multiple BLASTs.')
            ## GABLAM defaults that can be over-ridden by commandline
            gabdefault = ['keepblast=T','fullblast=T','blastdir=%s' % gdir,'dismat=F','distrees=F','disgraph=F','dotplots=F']
            ## GABLAM options that must be set
            gabcmd = gabdefault + self.cmd_list + ['qassemble=T','outstats=GABLAM','qryacc=F','percres=T','localcut=%d' % self.getInt('MinLocLen'),'basefile=%s' % gbase]
            fascmd = ['fasout=T','fragfas=T','combinedfas=T','localcut=0']

            ### ~ [1] Perform different GABLAM searches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # 1. QAssemble GABLAM of the reference genome against the assembly to get full genome coverage.
            if gtype == 'Reference':
                gcmd = ['seqin=%s' % self.getStr('RefGenome'),'searchdb=%s' % self.getStr('Assembly'),'dna=T','blastp=blastn','dotplots=%s' % self.getBool('DotPlots'),'dotlocalmin=%s' % self.getInt('MinLocLen')]
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 2. QAssemble GABLAM of the assembly versus reference genome to assess assembly accuracy in terms of excess sequence.
            if gtype == 'Assembly':
                gcmd = ['seqin=%s' % self.getStr('Assembly'),'searchdb=%s' % self.getStr('RefGenome'),'dna=T','blastp=blastn']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 3. QAssemble GABLAM of the reference genes (from Genbank annotation) to assess accuracy in terms of annotated features.
            if gtype == 'Genes':
                fasdir = '%s.GeneHits/' % gbase
                gcmd = fascmd + ['seqin=%s.gene.fas' % refbase,'searchdb=%s' % self.getStr('Assembly'),'dna=T','blastp=blastn','fasdir=%s.GeneHits/'% gbase]
                #gcmd.append('gablamfrag=1') #!# No: added FragMerge instead!
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
                if rje.exists(fasdir) and self.getBool('GeneTar'):
                    rje.targz(self,fasdir)
                    rje.deleteDir(self,fasdir,contentsonly=False,confirm=self.dev() or self.debugging() or self.i() > 1,report=True)
            # Perform Fragments search of hit genes
            if gtype == 'Genes.Fragments':
                gcmd = ['localcut=0','seqin=%s.gene.fas' % refbase,'searchdb=%s.fas' % gbase[:-10],'dna=T','blastp=blastn']
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Reciprocal search of hit genes
            if gtype == 'Genes.Reciprocal':
                gcmd = ['localcut=0','searchdb=%s.gene.fas' % refbase,'seqin=%s.fas' % gbase[:-11],'dna=T','blastp=blastn']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
            # 4. QAssemble GABLAM of the reference proteins (from Genbank annotation) to assess accuracy in terms of proteome coverage.
            if gtype == 'Proteins':
                fasdir = '%s.ProtHits/' % gbase
                gcmd = fascmd + ['seqin=%s.prot.fas' % refbase,'searchdb=%s' % self.getStr('Assembly'),'blastp=tblastn','fasdir=%s.ProtHits/'% gbase]
                #gcmd.append('gablamfrag=1') #!# No: added FragMerge instead!
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
                if rje.exists(fasdir) and self.getBool('GeneTar'):
                    rje.targz(self,fasdir)
                    rje.deleteDir(self,fasdir,contentsonly=False,confirm=self.dev() or self.debugging() or self.i() > 1,report=True)
            # Perform Fragments search of hit protein-coding sequences by proteins
            if gtype == 'Proteins.Fragments':
                gcmd = ['localcut=0','seqin=%s.prot.fas' % refbase,'searchdb=%s.fas' % gbase[:-10],'blastp=tblastn']
                self.runGABLAM(gabcmd + gcmd,gtype)
            # Perform Reciprocal search of hit protein-coding sequences versus proteins
            #Q# Should these be translated first?
            if gtype == 'Proteins.Reciprocal':
                gcmd = ['localcut=0','searchdb=%s.prot.fas' % refbase,'seqin=%s.fas' % gbase[:-11],'blastp=blastx']
                #gablam.GABLAM(self.log,gabcmd + gcmd).run()
                self.runGABLAM(gabcmd + gcmd,gtype)
        except: self.errorLog('%s.qAssembleGABLAM error' % self.prog())
#########################################################################################################################
    def selfQAssemble(self,genome=None):    ### Runs QAssemble GABLAM against self.
        '''
        Runs QAssemble GABLAM against self.
        >> genome:str = Genome to self QAssemble. Uses self.getStr('RefGenome') if None.
        << selfdb = Database Object with HitSum and Local tables loaded.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','selfQAssemble')
            if not genome: genome = self.getStr('RefGenome')
            genbase = '%s.L%dID%d' % (rje.baseFile(genome),self.getInt('MinLocLen'),self.getInt('MinLocID'))
            #i# This is using its own Database object
            selfdb = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % genbase])
            ## ~ [1] Check for existing run data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.force():
                hdb = selfdb.addTable(mainkeys=['Qry'],name='hitsum',expect=False)
                gdb = selfdb.addTable(mainkeys=['Qry','Hit'],name='gablam',expect=False)
                ldb = selfdb.addTable(mainkeys=['Qry','Hit','AlnNum'],name='local',expect=False)
                if hdb and ldb and gdb: return selfdb
            ## ~ [2] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Might want to set forks for speed up
            if self.getInt('Forks') < 2: self.printLog('#INFO','Consider setting forks=X to speed up multiple BLASTs.')
            ## Process reference gb file if required and find genome, genes and proteins
            blastdir = rje.makePath(os.path.dirname(genbase))
            gabdefault = ['keepblast=T','fullblast=T','blastdir=%s' % blastdir,'dismat=F','distrees=F','disgraph=F','dotplots=F']
            ## GABLAM options that must be set
            gabcmd = gabdefault + self.cmd_list + ['qassemble=T','outstats=GABLAM','qryacc=F','percres=T']
            gabcmd += ['seqin=%s' % genome,'dna=T','blastp=blastn','basefile=%s' % genbase,'selfhit=T','selfsum=T']
            #gablam.GABLAM(self.log,gabcmd).run()
            self.runGABLAM(gabcmd,'Self')
            hdb = selfdb.addTable(mainkeys=['Qry'],name='hitsum',expect=False)
            gdb = selfdb.addTable(mainkeys=['Qry','Hit'],name='gablam',expect=False)
            ldb = selfdb.addTable(mainkeys=['Qry','Hit','AlnNum'],name='local',expect=False)
            if hdb and ldb and gdb: return selfdb
            else: return None
        except: self.errorLog('%s.selfQAssemble error' % self.prog()); return None
#########################################################################################################################
    def disMatrixOut(self,gtables=[]): ### Output BLAST GABLAM distance matrix and Tree.
        '''
        Output BLAST GABLAM distance matrix and Tree.
        >> gtables:list [] = List of GABLAM summary tables to combine into dismatrix and tree.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','disMatrixOut')
            gablam = rje_dismatrix.DisMatrix(self.log,self.cmd_list+['basefile=%s' % self.baseFile()])
            gablam.setStr({'Name':'%s GABLAM' % self.baseFile()})
            diskey = 'Qry_AlnID'
            namedict = {}
            spcode = None

            ### ~ [1] ~ Generate and output DisMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for table in gtables:
                #self.debug(self.printLog('#DIS','Adding distances from %s' % table.name()))
                if not table: raise ValueError('Problem with PacBio dismatrix tables: %s' % gtables)
                table.dataFormat({diskey:'float','Hit_AlnID':'float'})
                for entry in table.entries():
                    if not spcode: spcode = string.split(entry['Qry'],'_')[1]
                    gablam.addDis(entry['Qry'],entry['Hit'],100.0-entry[diskey])
                    # Ref vs Assembly needs to be added both ways round
                    if not gablam.getDis(entry['Hit'],entry['Qry']): gablam.addDis(entry['Hit'],entry['Qry'],100.0-entry['Hit_AlnID'])
                    #self.debug('%s vs %s (%s & %s)' % (entry['Qry'],entry['Hit'],gablam.getDis(entry['Qry'],entry['Hit']),gablam.getDis(entry['Hit'],entry['Qry'])))
                    # Add to namedict
                    if entry['Qry'] not in namedict:
                        slen = float(entry['QryLen'])
                        if slen > 1e6: namedict[entry['Qry']] = '%s (%.2f Mb)' % (entry['Qry'],slen/1e6)
                        elif slen > 1e3: namedict[entry['Qry']] = '%s (%.2f kb)' % (entry['Qry'],slen/1e3)
                        else: namedict[entry['Qry']] = '%s (%.3f kb)' % (entry['Qry'],slen/1e3)
                        if entry['Qry'] in self.dict['QV']: namedict[entry['Qry']] = '%s; QV=%s)' % (namedict[entry['Qry']][:-1],rje.sf(self.dict['QV'][entry['Qry']],3))
                    if entry['Hit'] not in namedict:
                        slen = float(entry['HitLen'])
                        if slen > 1e6: namedict[entry['Hit']] = '%s (%.2f Mb)' % (entry['Hit'],slen/1e6)
                        elif slen > 1e3: namedict[entry['Hit']] = '%s (%.2f kb)' % (entry['Hit'],slen/1e3)
                        else: namedict[entry['Hit']] = '%s (%.3f kb)' % (entry['Hit'],slen/1e3)
                        if entry['Hit'] in self.dict['QV']: namedict[entry['Hit']] = '%s; QV=%s)' % (namedict[entry['Hit']][:-1],rje.sf(self.dict['QV'][entry['Hit']],3))
            gablam.saveMatrix(filename='%s.dismatrix.csv' % self.baseFile())

            ### ~ [2] ~ Make DisMatrix symmetrical on mean difference and generate tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gablam.forceSymmetry(method='mean',missing=100.0)
            self.setNum({'MST':gablam.MST()})
            ## ~ [2a] ~ Generate tree outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            upgma = rje_tree.Tree(self.log,['savetype=none','treeformats=nwk,text,png,r']+self.cmd_list+['autoload=F'])
            nsftree = gablam.upgma()
            upgma.buildTree(nsftree,type='nsf',postprocess=False)
            self.printLog('#TREE','Total TreeLen=%s; MST=%s' % (rje.sf(upgma.treeLen(),3),rje.sf(gablam.MST(),3)))
            for node in upgma.node:
                if node.info['Name'] in namedict: node.info['Name'] = namedict[node.info['Name']]
            upgma.basefile(self.basefile())
            upgma.info['GroupSpecies'] = spcode
            self.printLog('#SPEC','Looking for tree duplications based on %s' % spcode)
            upgma.opt['QueryGroup'] = True
            rje_tree_group._dupGroup(upgma,useseq=False) #.findDuplicationsNoSeq(spcode)
            upgma.saveTrees()
            #gablam.savePNG(upgma)
        except: self.errorLog('Major problem with %s.disMatrixOut' % self)
#########################################################################################################################
    ### <4> ### PAGSAT Assessment Methods                                                                               #
#########################################################################################################################
    def assessment(self,qassemble=True):   ### Generate GABLAM analyses and summarise assembly assessment.
        '''Generate GABLAM analyses and summarise assembly assessment.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','assessment')
            self.printLog('#~~#','## ~~~~~ PAGSTAT Assessment ~~~~~ ##')
            db = self.db()
            if not self.force():
                wanted = ['Summary.tdt',    # Summary delimited text file
                          'png']            # Summary tree of assembly vs reference chromosomes
                complete = True
                for wext in wanted:
                    wfile = '%s.%s' % (self.baseFile(),wext)
                    self.printLog('#CHECK','%s: %s.' % (wfile,os.path.exists(wfile)))
                    complete = complete and os.path.exists(wfile)
                if complete:
                    self.printLog('#SKIP','Assessment run found (force=F).')
                    if self.getBool('RGraphics'): return self.rGraphics()
                    return True
            minloclen = self.getInt('MinLocLen')    #!# Add to basefile #!#
            minid = self.getNum('MinLocID') / 100.0
            ## ~ [0a] Load Sequence files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Read in reference and assembly to seqlists
            if self.dev():
                self.obj['RefSeq'] = refseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('RefGenome'),'autoload=T','seqmode=filedb'])
                self.obj['Assembly'] = assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=filedb'])
            else:
                self.obj['RefSeq'] = refseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('RefGenome'),'autoload=T','seqmode=file'])
                self.obj['Assembly'] = assembly = rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % self.getStr('Assembly'),'autoload=T','seqmode=file'])
            # Generate summary tables: db('sequences'): ['name','desc','gene','spec','accnum','length'],['name']
            self.obj['RefSeq'].summarise(sumdb=True)
            self.obj['Assembly'].summarise(sumdb=True)
            ## ~ [0b] Reference All-by-all (including self-assessment benchmark) ~~~~~~~~~~~~~~~~~~ ##
            # NOTE: Cannot do a combined reference+assembly self-GABLAM as QAssemble stats would be messed up.
            # >> This is actually better in some ways as reference self-search can be re-used!
            selfdb = self.selfQAssemble()
            if not selfdb: raise ValueError('Reference self-GABLAM failure.')
            for table in ['hitsum','local','gablam']:
                tdb = selfdb.getTable(table)
                tdb.setStr({'Name':'Self.%s' % table})
                tdb.baseFile(self.baseFile())
                db.list['Tables'].append(tdb)
            ## ~ [0c] Assembly All-by-all for dismatrix and tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aselfdb = self.selfQAssemble(self.getStr('Assembly'))
            if not aselfdb: raise ValueError('Assembly self-GABLAM failure.')
            ## ~ [0d] Load/Generate GABLAM Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #gdir = '%s.GABLAM/' % self.baseFile()
            gdir = self.getStr('GABLAMDir')
            gbase = self.fileBase('GABLAM','Cut')
            db.baseFile(gbase)
            # Look to add QAssemble GABLAM tables: run qAssembleGABLAM() if force/missing.
            self.qAssembleGABLAM('Reference')
            qrefdb = self.db(table='Reference.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
            if not qrefdb: raise IOError('Failed to generate %s.Reference.hitsum.tdt' % self.db().baseFile())
            gdb = self.db(table='Reference.gablam',add=True,forcecheck=True,mainkeys=['Qry','Hit'])
            # Read in local genome alignment results
            locdb = self.db(table='Reference.local',add=True,forcecheck=True,mainkeys=['Qry','Hit','AlnNum'])
            #if qassemble and not locdb: self.qAssembleGABLAM('Assembly'); return self.assessment(False)
            #if qassemble and not (qrefdb and gdb): self.qAssembleGABLAM('Assembly'); return self.assessment(False)
            self.qAssembleGABLAM('Assembly')
            qassdb = self.db(table='Assembly.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
            #if qassemble and not qassdb: self.qAssembleGABLAM('Reciprocal'); return self.assessment(False)
            #qassdb.setStr({'Name':'Assembly.hitsum'})

            # Read in gene and protein results
            ftbase = self.fileBase('GABLAM','Base')
            db.baseFile(ftbase)
            #i# NOTE: Reciprocal Gene and Protein searches have been disabled for now as not that useful.
            genedb = protdb = grepdb = prepdb = None
            if self.getBool('GeneSummary'):
                self.qAssembleGABLAM('Genes')
                genedb = self.db(table='Genes.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not genedb: self.qAssembleGABLAM('Genes'); return self.assessment(False)
                self.qAssembleGABLAM('Genes.Fragments')    # Used only for TopHits Analysis
                #self.qAssembleGABLAM('Genes.Reciprocal')
                #grepdb = self.db(table='Genes.Reciprocal.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not grepdb: self.qAssembleGABLAM('Genes.Reciprocal'); return self.assessment(False)
            if self.getBool('ProtSummary'):
                self.qAssembleGABLAM('Proteins')            # Used only for TopHits Analysis
                protdb = self.db(table='Proteins.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not protdb: self.qAssembleGABLAM('Proteins'); return self.assessment(False)
                self.qAssembleGABLAM('Proteins.Fragments')  # Used only for TopHits Analysis
                #self.qAssembleGABLAM('Proteins.Reciprocal')
                #prepdb = self.db(table='Proteins.Reciprocal.hitsum',add=True,forcecheck=True,mainkeys=['Qry'])
                #if qassemble and not prepdb: self.qAssembleGABLAM('Proteins.Reciprocal'); return self.assessment(False)
            db.baseFile(self.baseFile())
            for table in db.tables(): table.baseFile(self.baseFile())

            # Plot both chr-contig and chr-chr in same table
            # Generate stats on duplication (contigs > self) and fragmentation (no local hits?) for each chromosome
            # ? Can we also map genes and proteins to chromosomes and make the summary per chromosome as well as total?
            # ? Should we also perform the gene/protein searches against both Reference and Assembly for benchmark?


            ### ~ [1] Order assembly contigs and reference by matches (based on biggest local alignments) ~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSAT Contig Ordering ~~~~~ ##')
            locdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            locdb.dropField('Positives')
            locdb.makeField('Identity/Length','Local')
            locdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor local hits')
            # Create lists of contigs and chromosomes from seqlists
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            #self.debug(chrdict)
            chrom = refseq.names()
            contigs = assembly.names()
            mapping = []
            # Work through local alignments in identity (count) order for primary pairings
            for lentry in locdb.sortedEntries('Identity',reverse=True):
                if not contigs: break
                if lentry['Hit'] not in contigs: continue
                # If contig OR chromosome still in list:
                # a. Add to assembly list: (chr,start,end,contig,start,end)
                mapping.append((lentry['Qry'],lentry['QryStart'],lentry['QryEnd'],lentry['Hit'],lentry['SbjStart'],lentry['SbjEnd']))
                # b. Remove chr and contig from lists
                if lentry['Qry'] in chrom: chrom.remove(lentry['Qry'])
                if lentry['Hit'] in contigs: contigs.remove(lentry['Hit'])
            if chrom: self.printLog('#CHROM','%d reference chromosomes without primary contig hits: %s' % (len(chrom),string.join(chrom,', ')))
            else: self.printLog('#CHROM','No reference chromosomes without primary contig hits.')
            if contigs: self.printLog('#CHROM','%d assembly contigs without primary reference hits: %s' % (len(contigs),string.join(contigs,', ')))
            else: self.printLog('#CHROM','No assembly contigs without primary reference hits.')
            # Sort assembly list
            mapping.sort()
            c2cmap = mapping[0:]
            # Cycle through chromosomes in order and output matching contigs to *.ordered.fas
            ordfile = '%s.ordered.fas' % self.baseFile()
            if self.force() or not rje.exists(ordfile):
                ordered = []
                for chrom in refseq.names():
                    for pair in mapping[0:]:
                        if pair[0] != chrom: continue
                        ordered.append(seqdict[pair[3]])
                        mapping.remove(pair)
                # Add (and log) any unmatched contigs to *.ordered.fas
                for seq in contigs: ordered.append(seqdict[seq])
                assembly.saveSeq(ordered,seqfile=ordfile)

            ### ~ [2] Chromosome-Contig Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            adb = self.db(table='ChromAlign',add=True,forcecheck=True,mainkeys=['Chrom'])
            #!# Need to add a diploid=T mode: duplicate queries to W&C copies, run, then recombine #!#
            if self.getBool('ChromAlign') and not adb and self.getBool('Diploid'): self.chromAlign(c2cmap,True)
            if self.getBool('ChromAlign') and not adb: adb = self.chromAlign(c2cmap)

            ### ~ [3] Generate Summary Tables, Charts and Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Summary Tables ~~~~~ ##')
            ## ~ [1a] Generate Distance Matrix and Tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Think about how best to make the distance matrix and tree. Add an assembly self-QAssemble and combine?
            self.disMatrixOut([gdb,self.db('Self.gablam'),aselfdb.getTable('gablam')])
            ## ~ [1b] Generate Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Add summary columns etc.
            sumtab = []
            for table in [self.db('Self.hitsum'),qrefdb,qassdb,genedb,protdb,grepdb,prepdb]:
                #!# Replace with tname list and self.db(tname)
                if not table: continue
                self.printLog('#TABLE','### ~~~~~~~~ %s ~~~~~~~~ ###' % table.name())
                table.setStr({'Name':string.replace(table.name(),'hitsum','Coverage')})
                table.dataFormat({'HitNum':'int','Positives':'int','EVal':'num','Length':'int','Coverage':'int','Identity':'int'})
                table.makeField('Length-Coverage','Missing',evalue=0)
                table.makeField('Coverage-Identity','Errors',evalue=0)
                table.dataFormat({'Missing':'int','Errors':'int'})
                table.addField('Perfect',evalue=0)
                for entry in table.entries():
                    if not entry['Missing'] and not entry['Errors']: entry['Perfect'] = 1
                table.saveToFile()
                sumtab.append(db.copyTable(table,string.replace(table.name(),'.Coverage','')))
            # Output summary statistics
            covdb = None
            if adb: fullsum = sumtab + db.splitTable(adb,'Type')
            else: fullsum = sumtab
            for table in fullsum:
                if table not in sumtab:
                    table.dropFields(['Insertions','Deletions']); table.renameField('Chrom','Qry')
                    table.setStr({'Name':string.split(table.name(),'_')[-1]+'Align'})
                else: table.dropFields(['MaxScore','EVal','Description'])
                table.addField('N',evalue=1)
                table.addField('Summary',evalue=table.name())
                table.compress(['Summary'],default='sum')
                if table not in sumtab: table.list['Fields'] = covdb.fields()
                else: table.list['Fields'] = ['Summary'] + table.fields()[:-1]
                table.dropField('Qry')
                if covdb: db.mergeTables(covdb,table)
                else: covdb = table; table.setStr({'Name':'Summary'})
            covdb.saveToFile()
            # Generate Excel file
            if self.dev() and False:
                try:
                    import XlsxWriter
                    #!# Add code here
                except: self.errorLog('XlsxWriter not on system: no Excel output')



            ### ~ [4] Generate coverage plotting data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ PAGSTAT Coverage Data ~~~~~ ##')
            ## ~ [4a] Coverage of Reference Chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chrdb = db.addEmptyTable('covplot.chrom',['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom'],['Chrom','Pos'])
            selfdb = self.db('Self.local')  # This is the chromosome versus chromosome plot
            selfdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            selfdb.dropField('Positives')
            selfdb.makeField('Identity/Length','Local')
            selfdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor (self) local hits')
            for chrom in locdb.index('Qry'):    # Should be the same for both.
                chromlen = refseq.seqLen(chrdict[chrom])
                posdict = {chromlen:{}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for pos in locdb.indexDataList('Qry',chrom,'QryStart') + selfdb.indexDataList('Qry',chrom,'QryStart'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos-1) not in posdict: posdict[pos-1] = {}
                for pos in locdb.indexDataList('Qry',chrom,'QryEnd') + selfdb.indexDataList('Qry',chrom,'QryEnd'):
                    if pos not in posdict: posdict[pos] = {}
                    if pos < chromlen and (pos+1) not in posdict: posdict[pos+1] = {}
                poslist = rje.sortKeys(posdict)
                for lentry in locdb.indexEntries('Qry',chrom) + selfdb.indexEntries('Qry',chrom):
                    hit = lentry['Hit']
                    for pos in poslist:
                        if pos < lentry['QryStart']: continue
                        if pos > lentry['QryEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['SbjStart'])
                for pos in poslist:
                    pentry = {'Chrom':chrom,'Pos':pos,
                              'ContigNum':0,'HitNum':0,'Class':'N','Contigs':[],
                              'ChromHit':0,'ChromNum':0,'RefChrom':[]}
                    for contig in posdict[pos]:
                        if contig in seqdict:   # Contig
                            pentry['ContigNum'] += 1
                            pentry['HitNum'] += len(posdict[pos][contig])
                            pentry['Contigs'].append(string.split(contig,'_')[0])
                        else:
                            pentry['ChromNum'] += 1
                            pentry['ChromHit'] += len(posdict[pos][contig])
                            pentry['RefChrom'].append(string.split(contig,'_')[0])
                    pentry['Contigs'].sort()
                    pentry['Contigs'] = string.join(pentry['Contigs'],';')
                    pentry['RefChrom'].sort()
                    pentry['RefChrom'] = string.join(pentry['RefChrom'],';')
                    # Class
                    if pentry['HitNum'] == 1: pentry['Class'] = 'U'
                    elif pentry['ContigNum'] == 1: pentry['Class'] = 'C'
                    elif pentry['ContigNum'] == pentry['HitNum'] == 2: pentry['Class'] = 'D'
                    elif pentry['ContigNum'] > 1: pentry['Class'] = 'M'
                    chrdb.addEntry(pentry)
            chrdb.saveToFile()
            ## ~ [4b] Coverage of Assembled Chromosomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Rationalise titles of these tables?
            # Make sure that SbjStart < SbjEnd. Store direction for 3c below.
            locdb.addField('Dirn',evalue='Fwd')
            for lentry in locdb.entries():
                if lentry['SbjStart'] > lentry['SbjEnd']:
                    (lentry['SbjStart'],lentry['SbjEnd']) = (lentry['SbjEnd'],lentry['SbjStart'])
                    lentry['Dirn'] = 'Rev'
            # This is now basically the same but swapping chromosomes and contigs
            chrdb = db.addEmptyTable('covplot.contig',['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom'],['Chrom','Pos'])
            chrdict = refseq.seqNameDic()   # Chromosomes to sequences
            ctgdict = seqdict               # Contigs to sequences
            selfdb = aselfdb.getTable('local')  # Assembly versus self
            selfdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            selfdb.dropField('Positives')
            selfdb.makeField('Identity/Length','Local')
            selfdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor (self) local hits')
            for chrom in locdb.index('Hit'):    # Should be the same for both.
                # First, make a dictionary of all the boundaries of local alignments
                # e.g. Every start and end position is a boundary as is every position BEFORE a start, or AFTER and end
                # This means that the change in numbers will occur between Starts and the position before, or ends and the position after
                posdict = {assembly.seqLen(ctgdict[chrom]):{}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for pos in locdb.indexDataList('Hit',chrom,'SbjStart') + selfdb.indexDataList('Qry',chrom,'QryStart'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos-1) not in posdict: posdict[pos-1] = {}
                for pos in locdb.indexDataList('Hit',chrom,'SbjEnd') + selfdb.indexDataList('Qry',chrom,'QryEnd'):
                    if pos not in posdict: posdict[pos] = {}
                    if (pos+1) not in posdict: posdict[pos+1] = {}
                # The boundary dictionary is a dictionary of dictionaries.
                # The value dictionaries contain contig keys (that hit this chrom) and a list of hit positions
                #?# Not sure it cycles through both hits and queries for Hit chrom. Does this not count them twice?! (Clearly NOT!)
                #?# Is one of the lists just empty?!
                #!# Tidy this up a bit
                poslist = rje.sortKeys(posdict)
                for lentry in locdb.indexEntries('Hit',chrom):
                    hit = lentry['Qry']
                    for pos in poslist:
                        if pos < lentry['SbjStart']: continue
                        if pos > lentry['SbjEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['QryStart'])
                for lentry in selfdb.indexEntries('Qry',chrom):
                    hit = lentry['Hit']
                    for pos in poslist:
                        if pos < lentry['QryStart']: continue
                        if pos > lentry['QryEnd']: break    # Next entry
                        if hit not in posdict[pos]: posdict[pos][hit] = []
                        posdict[pos][hit].append(lentry['SbjStart'])
                for pos in poslist:
                    pentry = {'Chrom':chrom,'Pos':pos,
                              'ContigNum':0,'HitNum':0,'Class':'N','Contigs':[],
                              'ChromHit':0,'ChromNum':0,'RefChrom':[]}
                    # Cycle through and update number of hit chromosomes and total BLAST hits.
                    for contig in posdict[pos]:
                        if contig in chrdict:   # Chromosome
                            pentry['ContigNum'] += 1
                            pentry['HitNum'] += len(posdict[pos][contig])
                            pentry['Contigs'].append(string.split(contig,'_')[0])
                        else:
                            pentry['ChromNum'] += 1
                            pentry['ChromHit'] += len(posdict[pos][contig])
                            pentry['RefChrom'].append(string.split(contig,'_')[0])
                    pentry['Contigs'].sort()
                    pentry['Contigs'] = string.join(pentry['Contigs'],';')
                    pentry['RefChrom'].sort()
                    pentry['RefChrom'] = string.join(pentry['RefChrom'],';')
                    # Class
                    if pentry['HitNum'] == 1: pentry['Class'] = 'U'
                    elif pentry['ContigNum'] == 1: pentry['Class'] = 'C'
                    elif pentry['ContigNum'] == pentry['HitNum'] == 2: pentry['Class'] = 'D'
                    elif pentry['ContigNum'] > 1: pentry['Class'] = 'M'
                    chrdb.addEntry(pentry)
            chrdb.saveToFile()

            ## ~ [4c] Directional coverage of Assembly Contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ctgdb = db.addEmptyTable('dirplot.contig',['Contig','Pos','FwdHit','FwdChrom','RevHit','RevChrom','Chroms','Class'],['Contig','Pos'])
            for chrom in locdb.index('Hit'):
                posdict = {'Fwd':{assembly.seqLen(seqdict[chrom]):{}},'Rev':{assembly.seqLen(seqdict[chrom]):{}}}    # Dictionary of all local alignment boundaries : {contig:[starts]}
                for lentry in locdb.indexEntries('Hit',chrom):
                    for dirn in ['Fwd','Rev']:
                        pos = lentry['SbjStart']
                        if pos not in posdict[dirn]: posdict[dirn][pos] = {}
                        if (pos-1) not in posdict[dirn]: posdict[dirn][pos-1] = {}
                        pos = lentry['SbjEnd']
                        if pos not in posdict[dirn]: posdict[dirn][pos] = {}
                        if (pos+1) not in posdict[dirn]: posdict[dirn][pos+1] = {}
                poslist = rje.sortUnique(rje.sortKeys(posdict['Fwd']) + rje.sortKeys(posdict['Rev']),num=True)
                for lentry in locdb.indexEntries('Hit',chrom):
                    dirn = lentry['Dirn']
                    hit = lentry['Qry']
                    for pos in poslist:
                        if pos < lentry['SbjStart']: continue
                        if pos > lentry['SbjEnd']: break    # Next entry
                        if hit not in posdict[dirn][pos]: posdict[dirn][pos][hit] = []
                        posdict[dirn][pos][hit].append(lentry['QryStart'])
                for pos in poslist:
                    pentry = {'Contig':chrom,'Pos':pos,'Class':'N','Chroms':[]}
                    for dirn in ['Fwd','Rev']:
                        pentry['%sHit' % dirn] = 0; pentry['%sChrom' % dirn] = 0
                        if pos in posdict[dirn]:
                            for contig in posdict[dirn][pos]:
                                pentry['%sHit' % dirn] += len(posdict[dirn][pos][contig])
                                pentry['%sChrom' % dirn] += 1
                                pchrom = string.split(contig,'_')[0]
                                if pchrom not in pentry['Chroms']: pentry['Chroms'].append(pchrom)
                    pentry['Chroms'].sort()
                    # Class
                    if (pentry['FwdHit']+pentry['RevHit']) == 1: pentry['Class'] = 'U'
                    elif len(pentry['Chroms']) == 1: pentry['Class'] = 'C'
                    elif len(pentry['Chroms']) > 1: pentry['Class'] = 'M'
                    pentry['Chroms'] = string.join(pentry['Chroms'],';')
                    ctgdb.addEntry(pentry)
            ctgdb.saveToFile()

            ## ~ [4c] Table of unique mapping based on sorted ctgdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Needs tuplekeys=T for sorting.
            #i# ctgdb = db.addEmptyTable('dirplot.contig',['Contig','Pos','FwdHit','FwdChrom','RevHit','RevChrom','Chroms','Class'],['Contig','Pos'])
            mapdb = db.addEmptyTable('mapping.contig',['Contig','Chrom','Fwd','Rev'],['Contig','Chrom'])
            mapchr = {}
            for chrname in refseq.names(): mapchr[string.split(chrname,'_')[0]] = chrname
            for contig in ctgdb.index('Contig'):
                ckeys = ctgdb.index('Contig')[contig][0:]
                while ckeys:    # These should be in pairs!
                    cstart = cend = ctgdb.data(ckeys.pop(0))
                    if cstart['Class'] != 'U': continue
                    while ckeys and ctgdb.data(ckeys[0])['Chroms'] == cstart['Chroms'] and ctgdb.data(ckeys[0])['Class'] == 'U':
                        cend = ctgdb.data(ckeys.pop(0))
                    self.debug(cstart)
                    if (len(string.split(cstart['Chroms'],';')) + cstart['FwdChrom'] + cstart['RevChrom']) != 2:
                        raise ValueError(cstart)
                    if (len(string.split(cend['Chroms'],';')) + cend['FwdChrom'] + cend['RevChrom']) != 2:
                        raise ValueError(cend)
                    # Now have a region mapping to a single chromosome
                    chrom = mapchr[cstart['Chroms']]
                    cpair = (contig,chrom)
                    mentry = mapdb.data(cpair)
                    if not mentry: mentry = mapdb.addEntry({'Contig':contig,'Chrom':chrom,'Fwd':0,'Rev':0})
                    mentry['Fwd'] += (cend['Pos'] - cstart['Pos'] + 1) * cstart['FwdChrom']
                    mentry['Rev'] += (cend['Pos'] - cstart['Pos'] + 1) * cstart['RevChrom']
            mapdb.saveToFile()

            ### ~ [5] Generate additional summary of Gene and Protein Top Hits using Ref features ~~~~~~~~~~~~~~~~~~~ ###
            self.topHits()

            ### ~ [6] Generate detailed chromosome to contig mappings from local alignments ~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('RGraphics'): return self.rGraphics()
        except: self.errorLog('%s.assessment error' % self.prog())
#########################################################################################################################
    def chromAlign(self,c2cmap,diploid=False): ### Generates PAGSAT Chrom Alignment tables
        '''Returns whether *Plots/*.summary.png has been made.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','chromAlign')
            self.printLog('#~~#','## ~~~~~ PAGSAT Chromosome Alignments (Diploid: %s) ~~~~~ ##' % diploid)
            db = self.db()
            alndir = rje.makePath('%s.ALN/' % self.baseFile())
            rje.mkDir(self,alndir)
            minloclen = self.getInt('MinLocLen')
            minid = self.getNum('MinLocID') / 100.0
            refseq = self.obj['RefSeq']
            assembly = self.obj['Assembly']
            ### ~ [1] Reference vs Assembly local BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Read BLAST into local database table including alignments
            blastfile = '%s.blast' % self.fileBase('GABLAM','Base','Reference')
            # Unzip file if required
            blastgz = rje.exists('%s.gz' % blastfile)
            if blastgz:
                os.system('gunzip %s.gz' % blastfile)
                self.printLog('#GUNZIP','%s unzipped.' % blastfile)
            if not rje.exists(blastfile): raise IOError('#BLAST','%s not found.' % (blastfile))
            # Better to re-read in the original BLAST
            blast = rje_blast.blastObj(self.log,self.cmd_list+['blastp=blastn'])
            blast.readBLAST(resfile=blastfile,clear=False,gablam=False,unlink=False,local=True,screen=True,log=False,keepaln=True)
            # (Re)zip BLAST file if required
            if blastgz:
                os.system('gzip %s' % blastfile)
                self.printLog('#GZIP','%s (re)zipped.' % blastfile)
            ## Use data from:
            bdb = blast.db('Local')
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            bdb.dataFormat({'Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int','Length':'int'})
            bdb.dropEntries(['Identity<%f' % minid,'Length<%d' % minloclen],inverse=False,log=True,logtxt='Removing poor local hits')
            # If diploid mode, duplicate queries here >> -A and -B
            if self.getBool('Diploid') and diploid:
                for entry in bdb.entries():
                    newentry = rje.combineDict({},entry)
                    newentry['Query'] = '%s-B' % entry['Query']
                    bdb.addEntry(newentry)
                    entry['Query'] = '%s-A' % entry['Query']
                bdb.remakeKeys()
            ### ~ [2] Chromosome-Contig Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bentries = bdb.sortedEntries('Identity',reverse=True)   # List of all entries (sorted) to process
            aentries = {}                                           # Dictionary of final entries for alignments
            alignpos = {}; ax = 0   # Dictionary of {chr/ctg:[(start,stop) list of positions included in local aln]}
            ## Process local alignments, keeping good ones and modifying the rest
            while bentries:
                entry = bentries.pop(0)     # This is best remaining hit
                ax += 1
                self.progLog('\r#LOCALN','Processing local alignments: %s -> %s' % (rje.iLen(bentries),rje.iLen(aentries)))
                entry['Qry'] = entry['Query']; entry['Sbj'] = entry['Hit']
                self.bugPrint('%s & %s' % (entry['Query'],entry['Sbj']))
                qh = 'Qry'
                if entry[qh] not in aentries: aentries[entry[qh]] = {}
                aentries[entry[qh]][entry['%sStart' % qh]] = entry
                # Update alignpos
                for qh in ['Qry','Sbj']:
                    if entry[qh] not in alignpos: alignpos[entry[qh]] = []
                    alignpos[entry[qh]].append((min(entry['%sStart' % qh],entry['%sEnd' % qh]),max(entry['%sStart' % qh],entry['%sEnd' % qh])))
                    alignpos[entry[qh]].sort()
                    i = 1
                    while i < len(alignpos[entry[qh]]):
                        if alignpos[entry[qh]][i][0] <= alignpos[entry[qh]][i-1][1]+1:
                            alignpos[entry[qh]][i-1] = (alignpos[entry[qh]][i-1][0],max(alignpos[entry[qh]][i-1][1],alignpos[entry[qh]][i][1]))
                            alignpos[entry[qh]].pop(i)
                        else: i += 1
                # Adjust/Filter remaining entries
                ex = 0
                while ex < len(bentries):
                    aentry = bentries[ex]
                    # Skip or delete if Query/Hit pair do not match
                    if entry['Query'] != aentry['Query'] and entry['Hit'] != aentry['Hit']: ex += 1; continue
                    if entry['Query'] != aentry['Query'] and entry['Hit'] == aentry['Hit']: bentries.pop(ex); continue
                    # Check for overlapping regions in Query
                    #??# What if in middle! :-(
                    #!# Need to fix this! #!#
                    dump = False
                    for (qstart,qend) in alignpos[entry['Query']]:
                        if qend < aentry['QryStart'] or qstart > aentry['QryEnd']: continue
                        if qstart <= aentry['QryStart'] and qend >= aentry['QryEnd']: dump = True; break
                        # Trim according to query positions
                        apos = 0                        # Position in alignment
                        qpos = aentry['QryStart'] - 1   # Position in query (prior to apos)
                        sdir = {True:1,False:-1}[aentry['SbjStart']<aentry['SbjEnd']]
                        spos = aentry['SbjStart'] - sdir  # Position in subject (prior to apos)
                        trimstart = aentry['QryEnd'] >= qend >= aentry['QryStart']  # and aentry['QryEnd'] >= qstart
                        trimend = aentry['QryStart'] <= qstart <= aentry['QryEnd'] #and aentry['QryStart'] <= qend
                        if (trimstart and trimend): self.bugPrint('Q: (%s,%s) vs (%s,%s)!' % (qstart,qend,aentry['QryStart'],aentry['QryEnd']))
                        while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                            if aentry['QrySeq'][apos] != '-': qpos += 1
                            if aentry['SbjSeq'][apos] != '-': spos += sdir
                            if trimstart and qpos > qend:   # Trim!
                                for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                aentry['QryStart'] = qpos
                                aentry['SbjStart'] = spos
                                apos = 1; trimstart = False
                            elif trimend and qpos == qstart - 1:
                                for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                aentry['QryEnd'] = qpos
                                aentry['SbjEnd'] = spos
                                break
                            else: apos += 1
                    if dump: bentries.pop(ex); continue
                    if entry['Query'] == aentry['Query'] and entry['Hit'] != aentry['Hit']: ex += 1; continue   # Can have 2+ contigs
                    # Check for overlapping regions in Hit
                    for (hstart,hend) in alignpos[entry['Hit']]:
                        if hend < min(aentry['SbjStart'],aentry['SbjEnd']) or hstart > max(aentry['SbjStart'],aentry['SbjEnd']): continue
                        if hstart <= min(aentry['SbjStart'],aentry['SbjEnd']) and hend >= max(aentry['SbjStart'],aentry['SbjEnd']): dump = True; continue
                        # Trim according to subject positions
                        apos = 0                        # Position in alignment
                        qpos = aentry['QryStart'] - 1   # Position in query (prior to apos)
                        sdir = {True:1,False:-1}[aentry['SbjStart']<aentry['SbjEnd']]
                        spos = aentry['SbjStart'] - sdir  # Position in subject (prior to apos)
                        if sdir > 0:    # Fwd hit
                            trimstart = aentry['SbjEnd'] >= hend >= aentry['SbjStart'] #and aentry['SbjEnd'] >= hstart
                            trimend = aentry['SbjStart'] <= hstart <= aentry['SbjEnd'] #and aentry['SbjStart'] <= hend
                            #if (trimstart and trimend): self.warnLog('(%s,%s) vs %s!' % (hstart,hend,aentry))
                            if (trimstart and trimend): self.bugPrint('Fwd: (%s,%s) vs (%s,%s)!' % (hstart,hend,aentry['SbjStart'],aentry['SbjEnd']))
                            while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                                if aentry['QrySeq'][apos] != '-': qpos += 1
                                if aentry['SbjSeq'][apos] != '-': spos += sdir
                                if trimstart and spos > hend:   # Trim!
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                    aentry['QryStart'] = qpos
                                    aentry['SbjStart'] = spos
                                    apos = 1; trimstart = False
                                elif trimend and spos == hstart - 1:
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                    aentry['QryEnd'] = qpos
                                    aentry['SbjEnd'] = spos
                                    break
                                else: apos += 1
                        else:
                            trimstart = aentry['SbjEnd'] <= hstart <= aentry['SbjStart'] #and aentry['SbjEnd'] <= hend
                            trimend = aentry['SbjStart'] >= hend >= aentry['SbjEnd'] #and aentry['SbjStart'] >= hstart
                            #if (trimstart and trimend): self.warnLog('(%s,%s) vs %s!' % (hstart,hend,aentry))
                            if (trimstart and trimend): self.bugPrint('Bwd: (%s,%s) vs (%s,%s)!' % (hstart,hend,aentry['SbjStart'],aentry['SbjEnd']))
                            while apos < len(aentry['AlnSeq']) and (trimstart or trimend):
                                if aentry['QrySeq'][apos] != '-': qpos += 1
                                if aentry['SbjSeq'][apos] != '-': spos += sdir
                                if trimstart and spos == hstart - 1:   # Trim!
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][apos:]
                                    aentry['QryStart'] = qpos
                                    aentry['SbjStart'] = spos
                                    apos = 1; trimstart = False
                                elif trimend and spos == hend + 1:
                                    for qh in ['Qry','Sbj','Aln']: aentry['%sSeq' % qh] = aentry['%sSeq' % qh][:apos+1]
                                    aentry['QryEnd'] = qpos
                                    aentry['SbjEnd'] = spos
                                    break
                                else: apos += 1
                    if dump: bentries.pop(ex); continue
                    ex += 1
            #self.debug(alignpos)
            self.printLog('\r#LOCALN','Processed local alignments -> %s pairwise aln.' % (rje.iStr(ax)))
            #bdb.dict['Data'] = {}
            #for entry in aentries: bdb.addEntry(entry)
            #bdb.newKey(['Query','QryStart'])
            #self.debug(bdb.dataKeys())
            ## Generate localn files:
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            adb = self.db('ChromAlign')
            if not adb: adb = db.addEmptyTable('ChromAlign',['Chrom','Ctid','Type','Length','Identity','Coverage','Insertions','Deletions'],['Chrom','Ctid','Type'])
            ddb = db.addEmptyTable('ChromAlignLoc', ['Query','Hit','Ctid','AlnID','Length','Identity','QryStart','QryEnd','SbjStart','SbjEnd'],['Query','Hit','Ctid','AlnID'])
            for chrom in rje.sortKeys(aentries):
                self.debug(chrom)
                if self.getBool('Diploid') and diploid: dentry = {'Query':chrom[:-2],'Ctid':chrom[-1]}
                else: dentry = {'Query':chrom,'Ctid':'H'}
                for entry in aentries[chrom].values():
                    #self.debug(entry)
                    #if self.getBool('Diploid') and diploid: dentry = {'Query':chrom[:-2],'Ctid':entry['Query'][-1]}
                    #else: dentry = {'Query':chrom,'Ctid':'H'}
                    ddentry = rje.combineDict({},dentry,overwrite=False)
                    ddentry = rje.combineDict(ddentry,entry,overwrite=False)
                    ddb.addEntry(ddentry)
                    #self.debug(ddentry)
                if self.getBool('Diploid') and diploid: chrseq = refseq.getSeq(chrdict[chrom[:-2]])
                else: chrseq = refseq.getSeq(chrdict[chrom])
                contigs = []
                poslist = rje.sortKeys(aentries[chrom])
                entry = aentries[chrom][poslist[0]]
                qname = string.split(entry['Query'],'__')[0]
                if self.getBool('Diploid') and diploid: qname += chrom[-2:]
                #self.debug(qname)
                aname = string.split(qname,'_')[0] + '_ALIGN'
                hname = string.split(qname,'_')[0] + '_' + string.split(entry['Hit'],'_')[1]
                qentry = {'Chrom':qname,'Type':'Reference','Ctid':dentry['Ctid']}; hentry = {'Chrom':hname,'Type':'Assembly','Ctid':dentry['Ctid']}
                qpos = 0    # Last position covered
                qseq = aseq = hseq = ''
                for spos in poslist:
                    entry = aentries[chrom][spos]
                    if entry['Hit'] not in contigs: contigs.append(entry['Hit'])
                    qname += ' %s(%s..%s)' % (string.split(entry['Qry'],'_')[-1],entry['QryStart'],entry['QryEnd'])
                    hname += ' %s(%s..%s)' % (string.split(entry['Hit'],'_')[-1],entry['SbjStart'],entry['SbjEnd'])
                    addseq = chrseq[1][qpos:entry['QryStart']]
                    if addseq:
                        qseq += addseq
                        aseq += '-' * len(addseq)
                        hseq += '.' * len(addseq)
                    if qseq: qseq += '.....'; aseq += 'XXXXX'; hseq += '.....'  # Small spacer between alignments
                    qseq += entry['QrySeq']
                    aseq += entry['AlnSeq']
                    hseq += entry['SbjSeq']
                    qpos = entry['QryEnd']
                addseq = chrseq[1][qpos:]
                if addseq:
                    qseq += addseq
                    aseq += '-' * len(addseq)
                    hseq += '.' * len(addseq)
                aseq = string.replace(aseq,' ','x')
                qentry['Length'] = len(chrseq[1])
                hentry['Length'] = 0
                for ctg in contigs: hentry['Length'] += assembly.seqLen(seqdict[ctg])
                qentry['Identity'] = hentry['Identity'] = string.count(aseq,'|')
                qentry['Insertions'] = hentry['Deletions'] = string.count(hseq,'-')
                qentry['Deletions'] = hentry['Insertions'] = string.count(qseq,'-')
                qentry['Coverage'] = hentry['Coverage'] = string.count(aseq,'x') + qentry['Identity'] - qentry['Insertions'] - qentry['Deletions']
                if not diploid:
                    ALNFAS = open('%s%s.aln.fas' % (alndir,chrom),'w')
                    ALNFAS.write('>%s\n%s\n' % (qname,qseq))
                    ALNFAS.write('>%s\n%s\n' % (aname,aseq))
                    ALNFAS.write('>%s\n%s\n' % (hname,hseq))
                    ALNFAS.close()
                    adb.addEntry(qentry)
                if diploid == self.getBool('Diploid'):
                    adb.addEntry(hentry)
            if diploid == self.getBool('Diploid'):
                ddb.saveToFile()
            if not diploid:
                adb.makeField('Length-Coverage','Missing',evalue=0)
                adb.makeField('Coverage-Identity','Errors',evalue=0)
                adb.addField('Perfect',evalue=0)
                for entry in adb.entries():
                    if not entry['Missing'] and not entry['Errors']: entry['Perfect'] = 1
                adb.compress(['Chrom'],default='sum')
                adb.dropField('Ctid')
                adb.saveToFile()
            #else: return adb

            # Build up and combine local alignments
            # >> Build up best -> Worst (BitScore), keeping if new sequences involved and trimming ends where appropriate
            # >> Fill in the gaps by splitting the shortest sequence in the middle and inserting gaps
            # Q. How to deal with inversions? RevComp and lower case?

            # Create lists of contigs and chromosomes from seqlists
            chrdict = refseq.seqNameDic()
            seqdict = assembly.seqNameDic()
            # Cycle through chromosomes in order and output matching contigs to *.ordered.fas
            for chrom in rje.sortKeys(chrdict): # refseq.names():
                #self.debug('%s -> %s' % (chrom,chrdict[chrom]))
                ALNFAS = open('%s%s.fas' % (alndir,chrom),'w')
                ALNFAS.write(refseq.fasta(chrdict[chrom]))
                for pair in c2cmap[0:]:
                    if pair[0] != chrom: continue
                    c2cmap.remove(pair)
                    if pair[-2] > pair[-1]:     # reverse!
                        (name,sequence) = assembly.getSeq(seqdict[pair[3]])
                        name += 'Reverse Complement'
                        #self.bugPrint(sequence[:100])
                        sequence = rje_sequence.reverseComplement(sequence)
                        #self.deBug(sequence[-100:])
                        ALNFAS.write('>%s\n%s\n' % (name,sequence))
                    else: ALNFAS.write(assembly.fasta(seqdict[pair[3]]))
                ALNFAS.close()
                #!# The following takes a long time and should be replaced with a better genome alignment program! (Exonerate?)
                #rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s%s.fas' % (alndir,chrom),'autoload=T','autofilter=F','dna=T','align=F']).align(outfile='%s%s.aln.fas' % (alndir,chrom))
            return adb
        except: self.errorLog('%s.chromAlign error' % self.prog()); raise
#########################################################################################################################
    def rGraphics(self):    ### Generates PAGSAT R Graphics
        '''Returns whether *Plots/*.summary.png has been made.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','rGraphics')
            self.printLog('#~~#','## ~~~~~ PAGSTAT R Graphics ~~~~~ ##')
            basename = self.baseFile(strip_path=True)
            rsfile ='%s.Plots/%s.summary.png' % (self.baseFile(),basename) # Summary reference figure.
            if not self.force() and os.path.exists(rsfile):
                self.printLog('#CHECK','%s: Found. (Force=F)' % rsfile)
                return True
            rcmd = '%s --no-restore --no-save --args "pagsat_V1" "%s"' % (self.getStr('RPath'),self.baseFile())
            rcmd += ' "refbase=%s"' % rje.baseFile(self.getStr('RefGenome'))
            rcmd += ' "minloclen=%d"' % self.getInt('MinLocLen')
            for boolvar in ['Diploid','GeneSummary','ProtSummary','ChromAlign','Features']:
                if not self.getBool(boolvar): rcmd += ' "%s=F"' % boolvar.lower()
            rdir = '%slibraries/r/' % slimsuitepath
            rcmd += ' "rdir=%s" < "%srje.r" > "%s.r.tmp.txt"' % (rdir,rdir,self.baseFile())
            self.printLog('#RPNG',rcmd)
            problems = os.popen(rcmd).read()
            if problems:
                for ptxt in problems: self.warnLog(ptxt)
            # Optional cleanup of *.r.tmp.txt ?
            return os.path.exists(rsfile)
        except: self.errorLog('%s.rGraphics error' % self.prog()); return False
#########################################################################################################################
    def topHits(self):    ### Generates output for Gene/Protein TopHits analysis
        '''Returns the Reference Features table, if given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','topHits')
            if not self.getBool('Features'):
                self.printLog('#TOPHIT','No Gene/Protein TopHit analysis without Genbank features table.')
                return False
            self.printLog('#~~#','## ~~~~~ PAGSTAT TopHits/Synteny Analysis ~~~~~ ##')
            db = self.db()
            if not self.force() and rje.exists('%s.Genes.TopHits.tdt' % db.baseFile()) and rje.exists('%s.Proteins.TopHits.tdt' % db.baseFile()):
                self.printLog('#HITS','Genes.TopHits and Proteins.TopHits tables found (force=F).'); return True
            tophitbuffer = self.getNum('TopHitBuffer')
            if tophitbuffer < 0: self.warnLog('Cannot have TopHitBuffer < 0.0!'); tophitbuffer = 0.0
            refseq = self.obj['RefSeq']
            ftdict = db.splitTable(self.ftdb(),'feature',asdict=True,keepfield=False,splitchar=None,values=['gene','CDS'])    # Reference features tables
            ftdict['Genes'] = ftdict.pop('gene')
            #self.debug(ftdict['Genes'].entries()[:10])
            ftdict['Proteins'] = ftdict.pop('CDS')
            #self.debug(ftdict['Proteins'].index('note').keys()[:10])
            acc2chr = {}
            for seq in refseq.seqs(): acc2chr[string.split(refseq.seqAcc(seq),'.')[0]] = refseq.seqGene(seq)
            ### ~ [1] Process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for gtype in ['Genes','Proteins']:
                ## ~ [1a] Load GABLAM data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.getBool('%sSummary' % gtype[:4]): continue
                gdb = db.addTable(filename='%s.%s.Fragments.gablam.tdt' % (self.fileBase('GABLAM','Base'),gtype),mainkeys=['Qry','Hit'],name='%s.gablam' % gtype,expect=True)
                if not gdb: self.warnLog('%s.%s.Fragments.gablam.tdt missing!' % (self.fileBase('GABLAM','Base'),gtype)); continue
                ## ~ [1b] Generate TopHits and Synteny predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                synteny = rje_synteny.Synteny(self.log,self.cmd_list)   #i# Uses own tophitbuffer=X [1.0]
                synteny.baseFile(db.baseFile())
                synteny.obj['DB'] = db
                synteny.obj['RefSeq'] = self.obj['RefSeq']
                synteny.topHitSynteny(ftdict[gtype],gdb,proteins=gtype=='Proteins',tabname='%s.TopHits' % gtype)
        except: self.errorLog('%s.topHits error' % self.prog())
#########################################################################################################################
    ### <5> ### PAGSAT Report Methods                                                                                   #
#########################################################################################################################
    def ftdb(self): ### Returns the Reference Features table, if given
        '''Returns the Reference Features table, if given.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','ftdb')
            if self.obj['Features']: return self.obj['Features']
            if not self.getBool('Features'): return False
            db = self.db()
            ### ~ [1] Load Reference Feature Table (if refgenome given) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('RefGenome'):
                ftdb = db.addTable('%s.Feature.tdt' % rje.baseFile(self.getStr('RefGenome')),mainkeys=['locus','feature','position'],name='features')
                ftdb.dataFormat({'start':'int','end':'int'})
            else: self.printLog('#REFFT','Cannot load/report/assess features without refgenome=FILE.'); ftdb = None
            self.obj['Features'] = ftdb
            return self.obj['Features']
        except: self.errorLog('%s.ftdb error' % self.prog()); return None
#########################################################################################################################
    def report(self):   ### Generates HTML reports of PAGSAT assessment
        '''Generates HTML reports of PAGSAT assessment.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','report')
            basename = self.baseFile(strip_path=True)
            wanted = ['Summary.tdt',    # Summary delimited text file
                      'Reference.Coverage.tdt', # Summary of the combined coverage of reference chromosomes
                      'Assembly.Coverage.tdt',  # Summary of the combined coverage of assembly contigs
                      'png',            # Summary tree of assembly vs reference chromosomes
                      'Plots/%s.summary.png' % basename] # Summary reference figure.
            for wext in wanted:
                wfile = '%s.%s' % (self.baseFile(),wext)
                if not os.path.exists(wfile):
                    self.printLog('#CHECK','%s: Missing! Will generate.' % wfile)
                    self.assessment(); break
                self.printLog('#CHECK','%s: Found.' % wfile)
            for wext in wanted:
                wfile = '%s.%s' % (self.baseFile(),wext)
                if not os.path.exists(wfile): raise IOError('Cannot find %s!' % wfile)
            ## ~ [0a] Setup HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = rje_html.HTML(self.log,self.cmd_list)
            hfile = '%s.report.html' % self.baseFile()
            ## ~ [0b] Dictionary of links and descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            links = {'Summary Table':'summarytable','Reference Coverage':'refcov','Assembly Coverage':'asscov'}
            desc = {'tree':'Summary tree of chromosomes vs contigs (% global identity)',
                    'summary':'Summary plot of assembly against reference chromosomes',
                    'summarytable':'Summary table of assembly against reference chromosomes',
                    'contents':'Report contents and quick links',
                    'refcov':'Summary table of assembled reference coverage',
                    'asscov':'Summary table of assembled assembly coverage',
                    'chromalign':'Summary plot of aligned contigs against reference chromosomes',
                    'assembly':'Summary plot of reference chromosomes against assembly contigs'}
            ## ~ [0c] Individual PNG files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            chrdb = self.db('Reference.Coverage',add=True)
            chrdb.index('Qry')
            pngfiles = rje.getFileList(self,'%s.Plots/' % self.baseFile(),['%s.covplot.*.png' % basename])
            trimx = len('%s.covplot.' % basename)
            chrpng = {}
            ctgpng = {}
            for pfile in pngfiles:
                ctg = rje.stripPath(pfile)[trimx:-4]
                #self.debug(ctg)
                #!# NB. This does not work if genesummary=F! Can we read it a better way? (Read a TDT?)
                #if rje.exists(string.replace(pfile,'covplot','genehits')): chrpng[ctg] = pfile
                #else: ctgpng[ctg] = pfile
                if ctg in chrdb.index('Qry'): chrpng[ctg] = pfile
                else: ctgpng[ctg] = pfile

            ### ~ [1] Initial Summary and quick links to sections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hbody = ['<a name="head"><h1>%s Report</h1></a>' % basename,'']
            hbody += ['<a name="contents"><h2>Contents</h2></a>','']
            hlink = ['<p>','Quick Links:']
            for linkid in ['Contents','Summary Table','Reference Coverage','Assembly Coverage','Tree','Summary PNG','ChromAlign PNG','Assembly PNG']:
                if linkid in links: link = links[linkid]
                else: link = string.split(linkid)[0].lower()
                hlink.append('~ <a href="#%s" title="%s">%s</a>' % (link,desc[link],linkid))
            hlink += ['</p>','']
            hbody += hlink

            hbody += ['<p>','Reference Chromosomes:']
            for chrom in rje.sortKeys(chrpng):
                hbody.append('~ <a href="#%s" title="%s">%s</a>' % (chrom,chrom,chrom))
            hbody += ['</p>','']
            hbody += ['<p>','Assembly Contigs:']
            for chrom in rje.sortKeys(ctgpng):
                hbody.append('~ <a href="#%s" title="%s">%s</a>' % (chrom,chrom,chrom))
            hbody += ['</p>','']


            hbody += ['<a name="summarytable"><h2>%s Summary Table</h2></a>' % basename,'']
            sumtable = open('%s.Summary.tdt' % self.basefile()).read()
            hbody.append(rje_html.tableToHTML(sumtable,'\t',tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''))

            hbody += ['<a name="refcov"><h2>%s Reference Coverage</h2></a>' % basename,'']
            sumtable = open('%s.Reference.Coverage.tdt' % self.basefile()).read()
            hbody.append(rje_html.tableToHTML(sumtable,'\t',tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''))

            hbody += ['<a name="asscov"><h2>%s Assembly Coverage</h2></a>' % basename,'']
            sumtable = open('%s.Assembly.Coverage.tdt' % self.basefile()).read()
            hbody.append(rje_html.tableToHTML(sumtable,'\t',tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''))

            hbody += ['<a name="tree"><h2>%s Summary Tree</h2><a>' % basename,'']
            hbody.append('<a href="./%s.png"><img src="./%s.png" width="100%%" title="%s"></a>' % (basename,basename,desc['tree']))

            for png in ['summary','chromalign','assembly']:
                pngfile = '%s.Plots/%s.%s.png' % (self.baseFile(),basename,png)
                pnglink = './%s.Plots/%s.%s.png' % (basename,basename,png)
                hbody += ['<a name="%s"><h2>%s</h2><a>' % (png,desc[png]),'']
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,desc[png]))

            # SECTIONS:
            # Reference-based Summary

            #self.warnLog('PAGSAT.Report() only partially implemented. Sorry!')

            # Include features overlapping:
            # (a) the missing regions
            # (b) repeated regions (where assembly > reference) in *.covplot.chrom.tdt

            '''
            head MBG8150.SP16481.hcq.sgd.srt.1000.covplot.chrom.tdt
            Chrom	Pos	HitNum	ContigNum	Contigs	Class	ChromHit	ChromNum	RefChrom
            chrIII_S288C__BK006937.2	0	0	0		N	0	0
            chrIII_S288C__BK006937.2	1	2	2	hcq1;hcq11	D	1	1	chrIII
            chrIII_S288C__BK006937.2	10	3	3	hcq1;hcq11;hcq6	M	1	1	chrIII
            chrIII_S288C__BK006937.2	11	4	4	hcq0;hcq1;hcq11;hcq6	M	1	1	chrIII
            chrIII_S288C__BK006937.2	11225	2	2	hcq16;hcq7	D	2	2	chrIII;chrXI
            chrIII_S288C__BK006937.2	11226	1	1	hcq16	U	1	1	chrIII
            chrIII_S288C__BK006937.2	114	6	6	hcq0;hcq1;hcq11;hcq4;hcq5;hcq6	M	2	2	chrIII;chrXIV
            chrIII_S288C__BK006937.2	1149	9	9	hcq0;hcq1;hcq11;hcq16;hcq2;hcq4;hcq5;hcq6;hcq7	M	2	2	chrIII;chrXIV
            chrIII_S288C__BK006937.2	115	7	7	hcq0;hcq1;hcq11;hcq4;hcq5;hcq6;hcq7	M	2	2	chrIII;chrXIV
            '''


            ### ~ [2] Individual Chromosomes and Contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hbody += ['<h2>Reference Chromosomes</h2>']
            for chrom in rje.sortKeys(chrpng):
                hbody.append('<a name="%s"><h3>%s</h3></a>' % (chrom,chrom))
                pngfile = chrpng[chrom]
                pnglink = './%s.Plots/%s' % (basename,rje.stripPath(pngfile))
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,'%s coverage plot' % chrom))
                pngfile = string.replace(pngfile,'covplot','genehits')
                pnglink = './%s.Plots/%s' % (basename,rje.stripPath(pngfile))
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,'%s gene hits plot' % chrom))
            hbody += ['<h2>Assembly Chromosomes</h2>']
            for chrom in rje.sortKeys(ctgpng):
                hbody.append('<a name="%s"><h3>%s</h3></a>' % (chrom,chrom))
                pngfile = ctgpng[chrom]
                pnglink = './%s.Plots/%s' % (basename,rje.stripPath(pngfile))
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,'%s coverage plot' % chrom))


            ### ~ [X] Output HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            HTML = open(hfile,'w')
            HTML.write(html.htmlHead(title=basename,tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
            HTML.write(string.join(hbody,'\n'))
            HTML.write(html.htmlTail(tabber=False))
            HTML.close()
            self.printLog('#HTML','HTML report output: %s' % hfile)
        except: self.errorLog('%s.report error' % self.prog())
#########################################################################################################################
    ### <8> ### PAGSAT Assembly Tidying/Editing                                                                         #
#########################################################################################################################
    def tidy(self):   ### Semi-automated tidying and editing of assembly to generate final draft genome.
        '''Semi-automated tidying and editing of assembly to generate final draft genome.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','tidy')
            if self.getBool('Report'): self.report()
            self.printLog('#~~#','## ~~~~~ PAGSTAT Assembly Tidying/Editing ~~~~~ ##')
            db = self.db()
            basename = self.baseFile(strip_path=True)   # Text of basename for output to screen (not for file management)
            #if self.dev(): chrmap = 'unique'
            #else: chrmap = 'align'
            chrmap = self.getStrLC('ChrMap')
            if chrmap == 'unique': maptable = 'mapping.contig'
            elif chrmap == 'align': maptable = 'ChromAlignLoc'
            else:
                self.warnLog('Chromosome Mapping method "%s" not recognised: using "unique" hits.')
                maptable = 'mapping.contig'
            ## ~ [0a] Check/create PAGSAT files in *.PAGSAT/ directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Establish files needed for tidy=T
            wanted = ['%s.tdt' % maptable,    # Summary delimited text file
                      'Plots/%s.chromalign.png' % basename] # Summary reference figure.
            # Check for files
            if not rje.checkForFiles(wanted,'%s.' % self.baseFile(),log=self.log,cutshort=True,missingtext=' Will generate.'):
                self.assessment()
                # Check for files again. (Should have been made if missing during first check.)
                rje.checkForFiles(wanted,'%s.' % self.baseFile(),log=self.log,cutshort=True,ioerror=True)
            ## ~ [0b] Setup new directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Use self.fileBase() to return basefile for all input files (e.g. in *.PAGSAT/).
            #i# Use self.fileBase(resdir='Assemble') to return basefile for all output files.
            pagdir = self.getStr('ResDir')
            #Replace PAGSAT/ with TIDY/
            assdir = rje.makePath('%s.ASSEMBLE/' % string.join(string.split(pagdir,'.')[:-1],'.'))
            rje.mkDir(self,assdir)
            plotdir = rje.makePath('%s.PLOT/' % self.fileBase())    # Where to find PAGSAT plots.
            self.setStr({'AssembleDir':assdir,'PlotDir':plotdir})
            self.printLog('#ASSDIR',self.getStr('AssembleDir'))      # Directory for PAGSAT output
            ## ~ [0c] Setup HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            html = rje_html.HTML(self.log,self.cmd_list)
            hfile = '%s.assembly.html' % self.fileBase('Assemble')
            ## ~ [0d] Dictionary of links and descriptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            desc = {'chromosome':'Summary of contig to reference chromosome mapping',
                    'chromalign':'Summary plot of aligned contigs against reference chromosomes',
                    'reference':'Reference chromosome PAGSAT plots',
                    'assembly':'Assembly contig PAGSAT plots'}
            ## ~ [0e] Initial Summary and quick links to sections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hbody = ['<h1>%s Assembly Tidy</h1>' % basename,'']
            hbody += ['<p>','Quick Links:']
            for hid in ['Chromosome Map','ChromAlign PNG','Reference Plots','Assembly Plots']:
                link = string.split(hid)[0].lower()
                hbody.append('~ <a href="#%s" title="%s">%s</a>' % (link,desc[link],hid))
            hbody += ['</p>','']

            ### ~ [1] Tidy and assemble contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #?# What to do about tel/rRNA contigs etc. that do not really assemble - manual ID/renaming?
            db.baseFile('%s%s' % (assdir,basename))
            ## ~ [1a] Compress ChromAlignLoc tdt output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if maptable == 'mapping.contig':  # Generate mapping based on ('mapping.contig',['Contig','Chrom','Fwd','Rev'],['Contig','Chrom'])
                self.printLog('#CTGMAP','Contig:Chrom mapping from mapping.contig table.')
                # Query	Hit	Ctid	AlnID	Length	Identity	QryStart	QryEnd	SbjStart	SbjEnd
                # Ctid = H if haploid, else A/B. Think about adding C or N for initial mapping?
                cdb = db.addTable('%s.mapping.contig.tdt' % self.fileBase(),['Contig','Chrom'],name='chrmap',expect=True)
                cdb.dataFormat({'Fwd':'int','Rev':'int'})
                cdb.makeField(formula='Fwd+Rev',fieldname='Length')
                cdb.rankFieldByIndex('Contig','Length',newfield='Rank',rev=True,absolute=True,lowest=True)
                cdb.dropEntries(['Rank>1'],inverse=False)
                cdb.dropField('Rank')
                cdb.renameField('Chrom','Query')    #!# Should change the other Query/Hit to Contig and Chrom!
                cdb.renameField('Contig','Hit')
                # Reduce to best
                cdb.compress(['Hit'],default='max')
                cdb.makeField(formula='Fwd-Rev',fieldname='Dirn')
                cdb.addField('Ctid',evalue='H')
                cdb.saveToFile()
            else:
                self.printLog('#CTGMAP','Contig:Chrom mapping from ChromAlignLoc table.')
                # Query	Hit	Ctid	AlnID	Length	Identity	QryStart	QryEnd	SbjStart	SbjEnd
                # Ctid = H if haploid, else A/B.
                cdb = db.addTable('%s.ChromAlignLoc.tdt' % self.fileBase(),['Query','Hit','AlnID'],name='chrmap',expect=True)
                cdb.dataFormat({'AlnID':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                cdb.makeField(formula='SbjEnd-SbjStart',fieldname='Dirn')
                cdb.compress(['Query','Hit'],rules={'AlnID':'max','QryStart':'min','QryEnd':'max','SbjStart':'min','SbjEnd':'max'},default='sum')
                cdb.saveToFile()
            ## ~ [1b] Convert to ordered list of chromosomes per contig ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Add orphan circularisation to tidy
            chrmap = {}; ctg2name = {}; ctg2chrom = {}; newname = {}
            for entry in cdb.entries():
                chrom = string.split(entry['Query'],'_')[0] + entry['Ctid']
                ctg2name[chrom] = entry['Query']
                if chrom not in chrmap: chrmap[chrom] = []
                [ctg,spec,n,acc] =  string.split(entry['Hit'],'_')
                if self.getStrLC('NewAcc'): newacc = '%s.%s' % (self.getStr('NewAcc'),string.split(acc,'.')[1])
                else: newacc = acc
                if chrom.startswith('chr'): newname[entry['Hit']] = string.join([string.replace(chrom,'chr',self.getStr('NewChr')),spec,n,newacc],'_')
                elif chrom.startswith('mt'): newname[entry['Hit']] = string.join(['%sMT' % string.replace(chrom,'mt',self.getStr('NewChr')),spec,n,newacc],'_')
                else:
                    newname[entry['Hit']] = string.join([string.replace(chrom,'chr',self.getStr('NewChr')),spec,n,newacc],'_')
                    self.warnLog('Non chrN/mt sequence name: %s' % newname[entry['Hit']])
                if entry['Dirn'] < 0: newname[entry['Hit']] += ' RevComp'
                ctg2name[ctg] = entry['Hit']
                ctg2chrom[ctg] = chrom
                if maptable == 'ChromAlignLoc': chrmap[chrom].append((entry['QryStart'],ctg,entry['Dirn']))
                else: chrmap[chrom].append((-entry['Length'],ctg,entry['Dirn']))
            hbody += ['<a name="chromosome"><h2>%s Chromosome Map</h2>' % basename,'']
            chrmapurl = {}  # Dictionary of chrom/ctg name linked to the text of its mapped contigs/chrom
            for chrom in rje.sortKeys(chrmap):
                chrmap[chrom].sort()
                chrtxt = []
                hmap = []
                chrhtml = '<a href="#%s">%s</a>' % (ctg2name[chrom],chrom)
                for (x,ctg,dirn) in chrmap[chrom]:
                    chrmapurl[ctg2name[ctg]] = chrhtml
                    chrtxt.append(ctg)
                    hmap.append('<a href="#%s">%s</a>' % (ctg2name[ctg],ctg))
                    if dirn < 0: chrtxt[-1] += 'rev'; hmap[-1] += '(Rev)'
                self.printLog('#CHRMAP','%s: %s' % (chrom,string.join(chrtxt,'|')))
                chrmapurl[ctg2name[chrom]] = string.join(hmap,' | ')
                hbody += ['<a href="#%s"><b>%s:</b></a> %s<br>' % (ctg2name[chrom],chrom,string.join(hmap,' | '))]
            #!# Add orphan sequences here
            pagbase = string.join(string.split(basename,'.')[:-1],'.')
            hbody += ['<p>Please see <a href="../%s.PAGSAT/%s.report.html">PAGSAT Report</a> for contig plots of orphan contigs not found above.</p>' % (pagbase,basename)]
            #sumtable = open('%s.Summary.tdt' % self.basefile()).read()
            #hbody.append(rje_html.tableToHTML(sumtable,'\t',tabwidth='100%',tdwidths=[],tdalign=[],valign='center',thead=True,border=1,tabid=''))
            ## ~ [1c] Generate *.assemble.html output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for png in ['chromalign']:
                #pngfile = '%s.%s.png' % (self.fileBase(resdir='Plot'),png)
                pnglink = '../%s.PAGSAT/%s.Plots/%s.%s.png' % (pagbase,basename,basename,png)
                hbody += ['<a name="%s"><h2>%s</h2></a>' % (png,desc[png]),'']
                hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,desc[png]))
            hplots = ['<a name="reference"><h2>Reference PAGSAT Plots</h2><a>','']
            for png in cdb.indexKeys('Query'):
                desc[png] = 'Coverage plot of reference chromosome %s' % png
                pnglink = '../%s.PAGSAT/%s.Plots/%s.covplot.%s.png' % (pagbase,basename,basename,png)
                hplots += ['<a name="%s"><h2>%s</h2></a>' % (png,desc[png]),'']
                if png in chrmapurl: hplots += ['<p>Mapped to: %s</p>' % chrmapurl[png]]
                hplots.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,desc[png]))
            hplots += ['<a name="assembly"><h2>Reference PAGSAT Plots</h2><a>','']
            for png in cdb.indexKeys('Hit'):
                desc[png] = 'Coverage plot of assembly contig %s' % png
                pnglink = '../%s.PAGSAT/%s.Plots/%s.covplot.%s.png' % (pagbase,basename,basename,png)
                hplots += ['<a name="%s"><h2>%s</h2></a>' % (png,desc[png]),'']
                if png in chrmapurl: hplots += ['<p>Mapped to: %s</p>' % chrmapurl[png]]
                hplots.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pnglink,pnglink,desc[png]))
            hbody += hplots
            HTML = open(hfile,'w')
            HTML.write(html.htmlHead(title=basename,tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
            HTML.write(string.join(hbody,'\n'))
            HTML.write(html.htmlTail(tabber=False))
            HTML.close()
            self.printLog('#HTML','HTML summary output: %s' % hfile)

            ## ~ [1d] Load, rename and sort sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqfile = '%s%s.pagsat.fas' % (assdir,self.fileBase('Assembly'))
            revseqfile = '%s%s.reviewed.fas' % (assdir,self.fileBase('Assembly'))
            if rje.exists(seqfile) and not self.force() and (self.yesNo('%s found. Use existing file sequence assembly/names?' % seqfile)):
                seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'dna=T','autoload=F','seqmode=list']
                seqlist = rje_seqlist.SeqList(self.log,seqcmd)
                seqlist.loadSeq()
                self.printLog('#PAGMAP','PAGSAT-mapped contigs read from %s. Please check %s for incorrect mapping.' % (seqfile,hfile))
            elif rje.exists(revseqfile) and not self.force() and (self.yesNo('%s found. Use existing file sequence assembly/names?' % revseqfile)):
                seqcmd = self.cmd_list + ['seqin=%s' % revseqfile,'dna=T','autoload=F','seqmode=list']
                seqlist = rje_seqlist.SeqList(self.log,seqcmd)
                seqlist.loadSeq()
                self.printLog('#PAGMAP','PAGSAT-mapped contigs read from %s. Please check %s for incorrect mapping.' % (revseqfile,hfile))
            else:
                seqcmd = self.cmd_list + ['seqin=%s' % self.getStr('Assembly'),'dna=T','autoload=F','seqmode=list']
                seqlist = rje_seqlist.SeqList(self.log,seqcmd)
                seqlist.loadSeq()
                newseq = []
                while seqlist.list['Seq']:
                    (sname,sequence) = seqlist.list['Seq'].pop(0)
                    namedata = string.split(sname,maxsplit=1)
                    if len(namedata) == 1: short = namedata[0]; desc = ''
                    else: [short,desc] = namedata
                    if short not in newname:
                        self.warnLog('Contig %s not found in %s table.' % (short,maptable))
                        if not self.getBool('Orphans'):
                            self.printLog('#DEL','Deleted contig %s: not found in %s table.' % (short,maptable))
                            continue
                        newname[short] = string.join(['Orphan']+string.split(short,'_')[1:],'_')
                    if newname[short].endswith('RevComp'): sequence = rje_sequence.reverseComplement(sequence)
                    newseq.append(('%s %s' % (newname[short],desc),sequence))
                    self.printLog('#EDIT','%s -> %s' % (short,newname[short]))
                newseq.sort()
                seqlist.list['Seq'] = newseq
                seqlist.saveSeq(seqfile=seqfile)
                self.printLog('#PAGMAP','PAGSAT-mapped contigs output to %s. Please check %s for incorrect mapping.' % (seqfile,hfile))
            seqlist.setStr({'SeqIn':seqfile})

            ## ~ [1e] Option to review/accept/quit ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tdb = self.db().addEmptyTable('tidy',['Contig','Start','End','Desc','Seq'],['Contig','Start'])   # Table of final comparison data
            stepi = 0; steps = ['R','A','O','F']
            if self.i() <= 0: stepi = 1
            while steps:
                choice = self.choice('\n\n<R>eview/edit Contigs; <A>ssemble; <O>rphans; <F>inish Tidy; <Q>uit PAGSAT?',default=steps[stepi],confirm=True).upper()
                if choice in ['Q','QUIT']: break
                elif choice == 'A':
                    self.assemble(seqlist)
                    aseqfile = '%s%s.assemble.fas' % (assdir,self.fileBase('Assembly'))
                    seqlist.saveSeq(seqfile=aseqfile)
                    self.printLog('#PAGASS','PAGSAT Assembly cycle complete')
                    stepi = 2
                elif choice == 'R':
                    seqlist.edit()
                    seqlist.list['Seq'].sort()
                    rseqfile = '%s%s.reviewed.fas' % (assdir,self.fileBase('Assembly'))
                    seqlist.saveSeq(seqfile=rseqfile)
                    stepi = 1
                    #?# Remake HTML #?#
                    #!# Add improved tidy table checking using Seq field
                    if tdb.entries(): self.warnLog('Manual sequence edits may cause *.tidy.tdt mismatch.')
                elif choice == 'O': # <O>rphans
                    expand_orphans = self.yesNo('Treat any non-"%s" sequences as orphans?' % self.getStr('NewChr'))
                    olist = []; otxt = []
                    for seq in seqlist.list['Seq']:
                        (seqname,sequence) = seq
                        if seqname.lower().startswith('orphan'): olist.append(seq); otxt.append(seqname)
                        elif expand_orphans and not seqname.startswith(self.getStr('NewChr')): olist.append(seq); otxt.append(seqname)
                    if olist:
                        print '%s\n\n' % string.join(['\n%d Orphan sequences:' % len(olist)]+otxt,'\n - ')
                        if self.yesNo('Delete %d orphan sequences? (Review/Edit for individual changes.)' % len(olist)):
                            for seq in olist:
                                seqlist.list['Seq'].remove(seq)
                                self.printLog('#DEL','Deleted contig %s: Orphan contig.' % (string.split(seq[0])[0]))
                    elif self.getBool('Orphans'): print '\n\nNo Orphan contigs.\n\n'
                    else: self.printLog('#INFO','No Orphans permitted. Check earlier #DEL entries in log file.')
                    stepi = 3
                elif choice == 'F': # <F>inish
                    break
            if tdb.entries():
                tdb.dropField('Seq')
                tdb.saveToFile()    # Table of tidy joins to checking with read coverage

            ## ~ [1f] Save sequence files for additional processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,seqfile)
            seqlist.saveSeq(seqfile=seqfile)
            aseqfile = None     # This will get set by Haploid core output
            nonhapx = -1        # Number of non-haploid sequences
            if self.yesNo('Output CtidA(/H) sequences for "haploid core" PAGSAT run?'):
                ctidAseq = []
                nonhapx = 0
                for seq in seqlist.seqs():
                    if seqlist.seqGene(seq).endswith('A') or seqlist.seqGene(seq).endswith('H'):
                        #X#ctidAseq.append(seq)
                        (seqname,sequence) = seqlist.getSeq(seq)
                        seqname = '%s_%s__%s %s' % (seqlist.seqGene(seq)[:-1],seqlist.seqSpec(seq),seqlist.seqAcc(seq),seqlist.seqDesc(seq))
                        self.printLog('#CHR','%s -> %s' % (seqlist.shortName(seq),string.split(seqname)[0]))
                        ctidAseq.append((seqname,sequence))
                    else: nonhapx += 1
                aseqfile = '%s%s.haploid.fas' % (assdir,self.fileBase('Assembly'))
                seqlist.saveSeq(seqfile=aseqfile,seqs=ctidAseq)
            pseqfile = None     # This will get set by full contig output
            fulldefault = {True:'Y',False:'N'}[nonhapx != 0]
            if self.yesNo('Output renamed sequences with unitig numbers for full PAGSAT run? (Assumes X.Y accnum)',fulldefault):
                seqs = []
                for seq in seqlist.seqs():
                    (seqname,sequence) = seqlist.getSeq(seq)
                    acc = seqlist.seqAcc(seq)
                    utig = string.split(acc,'.')[-1]
                    seqname = '%s%s_%s__%s %s' % (seqlist.seqGene(seq),utig,seqlist.seqSpec(seq),acc,seqlist.seqDesc(seq))
                    self.printLog('#UTIG','%s -> %s' % (seqlist.shortName(seq),string.split(seqname)[0]))
                    seqs.append((seqname,sequence))
                pseqfile = '%s%s.ctidX.fas' % (assdir,self.fileBase('Assembly'))
                seqlist.saveSeq(seqfile=pseqfile,seqs=seqs,seqtuples=True)

            ### ~ [2] Option for running Snapper for CNV analysis etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Actually recommended to run Quiver first. This was added for a Year 2 intern analysis
            snapdefault = {True:'Y',False:'N'}[self.getBool('Snapper')]
            if pseqfile and self.yesNo('Run Snapper on "ctidX" full contig set?',snapdefault):
                self.printLog('#SNPMAP','Running SNAPPER: see %s.snapper.log' % pagbase)
                snapcmd = ['seqin=%s' % pseqfile,'reference=%s' % self.getStr('RefGenome'),'basefile=%s.Snapper/%s' % (pagbase,basename),'log=%s.snapper.log' % pagbase]
                (info,out,mainlog,cmd_list) = snapper.setupProgram(snapcmd)
                snapper.Snapper(mainlog,cmd_list).run()
                mainlog.endLog(info)
            elif aseqfile and self.yesNo('Run Snapper on "haploid core"?',snapdefault):
                self.printLog('#SNPMAP','Running SNAPPER: see %s.snapper.log' % pagbase)
                snapcmd = ['seqin=%s' % aseqfile,'reference=%s' % self.getStr('RefGenome'),'basefile=%s.Snapper/%s' % (pagbase,basename),'log=%s.snapper.log' % pagbase]
                (info,out,mainlog,cmd_list) = snapper.setupProgram(snapcmd)
                snapper.Snapper(mainlog,cmd_list).run()
                mainlog.endLog(info)

            ### ~ [3] Option for regenerating new plots etc. using PAGSAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if aseqfile and pseqfile and self.yesNo('Run self-PAGSAT of full contig set vs "haploid core" (tidy=F)?'):
                self.printLog('#PAGSAT','Running PAGSAT: see %s.diploid.log' % pagbase)
                pagcmd = ['assembly=%s' % pseqfile,'refgenome=%s' % aseqfile,'basefile=%s.diploid' % pagbase,'tidy=F',
                          'genesummary=F','protsummary=F']
                (info,out,mainlog,cmd_list) = setupProgram(pagcmd)
                PAGSAT(mainlog,cmd_list).run()
                mainlog.endLog(info)
            if aseqfile and self.yesNo('Run PAGSAT on "haploid core" vs Reference (tidy=F)?'):
                self.printLog('#PAGSAT','Running PAGSAT: see %s.haploid.log' % pagbase)
                pagcmd = ['assembly=%s' % aseqfile,'basefile=%s.haploid' % pagbase,'tidy=F']
                (info,out,mainlog,cmd_list) = setupProgram(pagcmd)
                PAGSAT(mainlog,cmd_list).run()
                mainlog.endLog(info)
            if pseqfile and self.yesNo('Run PAGSAT on "ctidX" full contig set vs Reference (tidy=F)?',fulldefault):
                self.printLog('#PAGSAT','Running PAGSAT: see %s.ctidX.log' % pagbase)
                pagcmd = ['assembly=%s' % pseqfile,'basefile=%s.ctidX' % pagbase,'tidy=F']
                (info,out,mainlog,cmd_list) = setupProgram(pagcmd)
                PAGSAT(mainlog,cmd_list).run()
                mainlog.endLog(info)

        except: self.errorLog('%s.tidy() error' % self.prog())
#########################################################################################################################
    def assembleChrom(self,seqlist,chrom,chromseq):     ### Run the assembly process for just one chromosome
        '''Run the assembly process for just one chromosome.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','assembleChrom')
            if self.i() > 0 and not rje.yesNo('Assess and assemble %s?' % chrom): return False
            tdb = self.db('tidy') # ['Contig','Start','End','Desc','Seq'],['Contig','Start'])
            seqdict = seqlist.makeSeqNameDic()  # Dictionary of shortname to sequence, updated for any edits
            sfile = '%s.%s.fas' % (self.fileBase('CtgGABLAM'),chrom)
            sbase = '%s.%s' % (self.fileBase('CtgGABLAM'),chrom)
            circularise = False
            if len(chromseq) == 1:
                try:
                    (seqname,sequence) = chromseq[0]
                    seqdesc = string.split(seqname,maxsplit=1)[1]
                except: seqdesc = ''
                self.debug(seqdesc)
                choicedef = {True:'N',False:'Y'}[seqdesc.startswith('Circle')]
                if self.yesNo('Check %s for circularity?' % chrom,default=choicedef): circularise = True
                else: return True   #!# Should check/add -A or -H suffix.
            ### ~ [1] Perform all-by-all BLAST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Save and perform BLAST
            seqlist.saveSeq(seqfile=sfile,seqs=chromseq,reformat='short',backup=False)
            #gcmd = self.cmd_list + ['seqin=%s' % sfile,'searchdb=%s' % sfile,'minloclen=%d' % self.getInt('MinLocLen'),'qryacc=F','dna=T','blastp=blastn','fullblast=T','basefile=%s' % gbase]
            bcmd = self.cmd_list + ['blasti=%s' % sfile,'blastd=%s' % sfile,'blastp=blastn','blastf=F','basefile=%s' % sbase,'blasto=%s.blast' % sbase]
            blast = rje_blast.blastObj(self.log,bcmd+['backups=F','gablamfrag=0'])
            blast.formatDB(fasfile=sfile,protein=False,force=True,log=True,checkage=None,details=False)
            blast.blast(wait=True,cleandb=True,use_existing=False,log=True)
            blast.readBLAST(clear=True,gablam=False,unlink=False,local=True,screen=True,log=True,keepaln=True)
            ## ~ [1a] Tidy and filter results tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # 'Run',['Run','Type','E-Value','DBase','InFile','BLASTCmd','Complexity Filter','Composition Statistics','SoftMask','GappedBLAST','OneLine','HitAln','DBLen','DBNum'],['Run']
            # 'Search',['Query','Length','Hits','MaxScore','TopE'],['Query']
            # 'Hit',['Query','Rank','Hit','Description','BitScore','E-Value','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Query','Hit']
            # 'Local',['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],['Query','Hit','AlnID']
            if not circularise:
                for table in ['Hit','Local']:    #,'GABLAM']:
                    blast.db(table).dropEntries(['Query==Hit'])   # No self-hits #!# This does not work. (Why?!)
            # blast.db('Local').dropEntries('Length<%d' % self.getInt('MinLocLen')) ?
            # Read Local table
            clocdb = blast.db('Local')
            clocdb.dropField('AlnSeq')
            clocdb.renameField('Query','Qry')
            # Pure check:
            sex = 0
            for entry in clocdb.entries():
                if entry['Qry'] == entry['Hit'] and not circularise:    #!# Something separate for circularisation #!#
                    #self.warnLog('%s Self-Hit detected and removed!' % entry['Qry'])
                    clocdb.dropEntry(entry); sex += 1
            self.printLog('#SELF','%s self-hits removed' % rje.iStr(sex))
            # Reformat data and add new fields
            clocdb.dataFormat({'AlnNum':'int','BitScore':'num','Expect':'num','Length':'int','Identity':'int','Positives':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
            clocdb.renameField('SbjStart','HitStart')
            clocdb.renameField('SbjEnd','HitEnd')
            clocdb.renameField('SbjSeq','HitSeq')
            clocdb.addFields(['QryLen','HitLen'])
            cqrydb = blast.db('Search')
            for entry in clocdb.entries():
                entry['QryLen'] = cqrydb.data(entry['Qry'])['Length']
                entry['HitLen'] = cqrydb.data(entry['Hit'])['Length']
            # Rate overlaps. (May be overkill but useful for clarity/checking.)
            clocdb.addFields(['QryType','HitType'])
            for entry in clocdb.entries():
                for qh in ['Qry','Hit']:
                    if entry['%sStart' % qh] <= self.getInt('JoinMargin') < entry['%sEnd' % qh]:
                        if entry['%sEnd' % qh] >= (entry['%sLen' % qh]-self.getInt('JoinMargin')): entry['%sType' % qh] = 'Full'
                        else: entry['%sType' % qh] = 'Start'
                    elif entry['%sEnd' % qh] >= (entry['%sLen' % qh]-self.getInt('JoinMargin')) >  entry['%sStart' % qh]: entry['%sType' % qh] = 'End'
                    elif entry['%sEnd' % qh] <= self.getInt('JoinMargin'):
                        if entry['%sStart' % qh] >= (entry['%sLen' % qh]-self.getInt('JoinMargin')): entry['%sType' % qh] = 'InvFull'
                        else: entry['%sType' % qh] = 'InvStart'
                    elif entry['%sStart' % qh] >= (entry['%sLen' % qh]-self.getInt('JoinMargin')): entry['%sType' % qh] = 'InvEnd'
                    elif entry['%sStart' % qh] > entry['%sEnd' % qh]: entry['%sType' % qh] = 'Inverted'
                    else: entry['%sType' % qh] = 'Internal'
            clocdb.indexReport('QryType',logstr='#QTYPE')
            clocdb.indexReport('HitType',logstr='#HTYPE')
            #if self.dev():
            savefields = clocdb.fields()
            savefields.remove('QrySeq'); savefields.remove('HitSeq')
            clocdb.saveToFile(savefields=savefields)
            cassdb = blast.db().copyTable(clocdb,newname='Assemble')
            # Check for bad RevComp signs = non-self 1-x / 1-y or x-L / y-L overlaps
            for qtype in clocdb.index('QryType'):
                if qtype.startswith('Inv'): raise ValueError('QryType should not be inverted!')
            # Reduce to terminal overlaps
            invseq = []
            clocdb.dropEntriesDirect('QryType',['Start','End'],inverse=True)
            for htype in clocdb.index('HitType'):
                if htype in ['InvStart','InvFull','InvEnd']:
                    self.warnLog('%s inverted terminal hits detected: possible local inversion or RevComp errors' % rje.iLen(clocdb.index('HitType')[htype]))
                    for centry in clocdb.indexEntries('HitType',htype):
                        self.printLog('#%s' % htype.upper()[:4],'%s: %s vs %s' % (htype,centry['Qry'],centry['Hit']))
                        invseq += [centry['Qry'],centry['Hit']]
            clocdb.dropEntriesDirect('HitType',['Start','End'],inverse=True)
            clocdb.dropEntries(['QryType==HitType'])    #!# Not working!
            eqx = 0
            for entry in clocdb.entries():
                if entry['QryType'] == entry['HitType']: clocdb.dropEntry(entry); eqx += 1
            if eqx: self.printLog('#TYPE','%s Start:Start or End:End hits removed.' % rje.iStr(eqx))
            while not clocdb.entryNum() and invseq:
                invseq = rje.sortUnique(invseq)
                self.printLog('#CJOIN','No possible contig joins: %d possible inverted sequences.' % len(invseq))
                ctext = '\n'
                for ci in range(len(invseq)): ctext += '<%d> : Reverse complement %s\n' % (ci+1,invseq[ci])
                ctext += '\n<0> : Continue without additional edits.'
                if self.i() >= 0: ji = rje.getInt(ctext,default=0,confirm=True)   #!# Add acceptable limits!
                else: ji = 0
                if not ji: break
                if ji < 0 or ji > len(invseq): continue
                sname = invseq[ji-1]
                iseq = seqdict[sname]
                si = seqlist.list['Seq'].index(iseq)
                (seqname,sequence) = seqlist.getSeq(iseq)
                chromseq.remove((seqname,sequence))
                seqname = string.split(seqname)
                if len(seqname) > 1 and seqname[1] == 'RevComp': seqname.pop(1)
                else: seqname.insert(1,'RevComp')
                seqname = string.join(seqname)
                sequence = rje_sequence.reverseComplement(sequence)
                self.printLog('#EDIT','%s -> %s' % (sname,seqname))
                seqdict[sname] = seqlist.list['Seq'][si] = (seqname,sequence)
                chromseq.insert(0,(seqname,sequence))
                return self.assembleChrom(seqlist,chrom,chromseq)
            if not clocdb.entryNum():
                self.printLog('#CJOIN','No possible contig joins: assigning chromatids.')
                qdb = blast.db('Search') # 'Search',['Query','Length','Hits','MaxScore','TopE'],['Query']
                hdb = blast.db('Hit')
                cassdb.dropEntriesDirect('HitType',['Internal','Inverted'])
                qsort = []
                for entry in qdb.entries(): qsort.append((entry['Length'],entry['Query']))
                qsort.sort(reverse=True)
                # The longest contig becomes chromatid A
                ctidA = qsort.pop(0)[1]
                #self.debug(ctidA)
                cassdb.dropEntriesDirect('Hit',[ctidA])
                cassdb.dropEntriesDirect('Qry',[ctidA],inverse=True)
                cassdb.index('Hit')
                #cassdb.makeField('#Qry#|#Hit#','QH')
                #cassdb.index('QH')
                ctid = {'A':[ctidA],
                        'B':[],  # List of contigs "neatly" contained in ctidA
                        'N':[]}  # List of "messy" contigs
                for (hlen,hit) in qsort:
                    overlaps = cassdb.indexDataList('Hit',hit,'HitType')
                    if 'Full' in overlaps or ('Start' in overlaps and 'End' in overlaps):   #?# Check End > Start?
                        ctid['B'].append(hit)
                    else: ctid['N'].append(hit)
                #self.debug(ctid)
                # Rename sequences
                for c in 'ABN':
                    self.printLog('#CTID%s' % c,string.join(ctid[c],'; '))
                    #self.debug(rje.sortKeys(seqdict))
                    for contig in ctid[c]:
                        cseq = seqdict[contig]
                        if seqlist.seqGene(cseq).endswith(c): continue
                        if seqlist.seqGene(cseq)[-1] in 'ABHN':
                            newgene = seqlist.seqGene(cseq)[:-1] + c
                        else: newgene = seqlist.seqGene(cseq) + c
                        newname = '%s_%s__%s %s' % (newgene,seqlist.seqSpec(cseq),seqlist.seqAcc(cseq),seqlist.seqDesc(cseq))
                        ci = seqlist.list['Seq'].index(cseq)
                        newseq = (newname,cseq[1])
                        seqlist.list['Seq'][ci] = newseq
                        seqdict.pop(contig)
                        seqdict[seqlist.shortName(newseq)] = newseq
                return True
            self.printLog('#CJOIN','%d(/2) possible %s contig joins' % (clocdb.entryNum(),chrom))
            ### ~ [X] Special Circularise assessment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #X# Circularise should not need special treatment here: sort out in join step below.
            #if circularise:
            #    self.printLog('#DEV','Circularisation not yet implemented!')
            #    #!# Add specific analysis looking for self start/end overlap
            #    return False
            ### ~ [3] Work through possible contig joins and repeat process until done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Changed this to be sorted by smallest overlap first: over-ruled several times on this basis.
            #?# Could add as an option: Length/Identity
            devopt_joinsort = self.getStrLC('JoinSort')   #'Length'
            cjoins = []     # List of (joinID,qryLen,entry) tuples
            for entry in clocdb.entries():
                if devopt_joinsort == 'identity':   # Add -ve identity to sort from big to small
                    cjoins.append((-float(entry['Identity'])/entry['Length'],entry['QryLen'],entry))
                elif devopt_joinsort == 'length':
                    cjoins.append((entry['Length'],entry['QryLen'],entry))
                else: raise ValueError('joinsort="%s" not recognised. (Length/Identity)' % self.getStr('JoinSort'))
                #if circularise: break   # Should only be one!
            cjoins.sort(reverse=False)  #devopt_joinsort in ['identity'])
            ctext = '\n'
            for ci in range(len(cjoins)):
                entry = cjoins[ci][-1]
                ctext += '<%d> : Join %s %s-%s with %s %s-%s (%s nt = %.2f%% identity)\n' % (ci+1,entry['Qry'],rje.iStr(entry['QryStart']),rje.iStr(entry['QryEnd']),entry['Hit'],rje.iStr(entry['HitStart']),rje.iStr(entry['HitEnd']),rje.iStr(entry['Length']),100.0*float(entry['Identity'])/entry['Length'])
            ctext += '<0> : Abort join for %s\n\nJoin choice?' % chrom
            if self.i() >= 0: ji = rje.getInt(ctext,default=1,confirm=True)
            else: ji = 1
            if not ji: self.warnLog('%s assembly aborted. Will need manual ctid assignment.' % chrom); return True #!# Might mess things up?
            ## ~ [3a] Perform join ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            entry = cjoins[ji-1][-1]
            # Make consensus
            if self.getStrLC('JoinMerge') == 'consensus':
                jname = 'Consensus:%s:%s-%s/%s:%s-%s (%snt = %.2f%% identity)' % (entry['Qry'],rje.iStr(entry['QryStart']),rje.iStr(entry['QryEnd']),entry['Hit'],rje.iStr(entry['HitStart']),rje.iStr(entry['HitEnd']),rje.iStr(entry['Length']),100.0*float(entry['Identity'])/entry['Length'])
                jsequence = seqlist.makeConsensus(entry['QrySeq'],[entry['HitSeq']])
                gapx = jsequence.count('-')
                jsequence = jsequence.replace('-','')
                #seqlist._addSeq(jname,jsequence)
                jseq = (jname,jsequence)
                self.printLog('#SEQ','Sequence added: %s = %s nt; %s gaps removed.' % (jname,rje.iLen(jsequence),rje.iStr(gapx)))
                # Trim Query and Hit
                if circularise:     # Trim off both ends!
                    if entry['Qry'] != entry['Hit']: raise ValueError
                    qh = 'Qry'
                    chromseq.remove(seqdict[entry[qh]])
                    (seqname,sequence) = seqdict[entry[qh]]
                    si = seqlist.list['Seq'].index(seqdict[entry[qh]])
                    if entry['%sType' % qh] == 'Start':
                        x = entry['%sEnd' % qh] + 1
                        y = entry['HitStart'] - 1
                    elif entry['%sType' % qh] == 'End':
                        x = entry['HitEnd'] + 1
                        y = entry['%sStart' % qh] - 1
                    else: raise ValueError(entry['%sType' % qh])
                    sequence = sequence[x-1:y]
                    self.printLog('#EDIT','%s -> Region %d to %d.' % (seqname,x,y))
                    seqname = '%s (Region %d to %d)' % (seqname,x,y)
                    seqdict[entry[qh]] = seqlist.list['Seq'][si] = (seqname,sequence)
                    # Perform join
                    qseq = seqdict[entry['Qry']]
                    joinseq = (jseq,qseq)
                    seqdesc = 'Circle[%s & %s]' % (joinseq[0][0],joinseq[1][0])
                    seqname = '%s %s' % (string.split(joinseq[1][0])[0],seqdesc)
                    sequence = joinseq[0][1] + joinseq[1][1]
                    qi = seqlist.list['Seq'].index(qseq)
                    seqdict[entry['Qry']] = seqlist.list['Seq'][qi] = (seqname,sequence)
                    self.printLog('#EDIT',seqdesc)
                    chromseq.insert(0,(seqname,sequence))
                else:
                    for qh in ['Qry','Hit']:
                        try: chromseq.remove(seqdict[entry[qh]])
                        except: self.warnLog('%s %s not found in %s chromseq dictionary?!' % (qh,entry[qh],chrom))
                        #!# This happens if manual editing of sequence between Assembly runs: make sure chromseq is regenerated at right point.
                        (seqname,sequence) = seqdict[entry[qh]]
                        si = seqlist.list['Seq'].index(seqdict[entry[qh]])
                        if entry['%sType' % qh] == 'Start':
                            x = entry['%sEnd' % qh] + 1
                            y = len(sequence)
                        elif entry['%sType' % qh] == 'End':
                            x = 1
                            y = entry['%sStart' % qh] - 1
                        else: raise ValueError(entry['%sType' % qh])
                        sequence = sequence[x-1:y]
                        self.printLog('#EDIT','%s -> Region %d to %d.' % (seqname,x,y))
                        seqname = '%s (Region %d to %d)' % (seqname,x,y)
                        seqdict[entry[qh]] = seqlist.list['Seq'][si] = (seqname,sequence)
                    # Perform join
                    qseq = seqdict[entry['Qry']]
                    hseq = seqdict[entry['Hit']]
                    if entry['QryType'] == 'End':
                        joinseq = (seqdict[entry['Qry']],jseq,seqdict[entry['Hit']])
                        seqdesc = 'Join[%s & %s & %s]' % (joinseq[0][0],joinseq[1][0],joinseq[2][0])
                        seqname = '%s %s' % (string.split(joinseq[0][0])[0],seqdesc)
                    elif entry['QryType'] == 'Start':
                        joinseq = (seqdict[entry['Hit']],jseq,seqdict[entry['Qry']])
                        seqdesc = 'Join[%s & %s & %s]' % (joinseq[0][0],joinseq[1][0],joinseq[2][0])
                        seqname = '%s %s' % (string.split(joinseq[2][0])[0],seqdesc)
                    else: raise ValueError(entry['QryType'])
                    sequence = joinseq[0][1] + joinseq[1][1] + joinseq[2][1]
                    qi = seqlist.list['Seq'].index(qseq)
                    hi = seqlist.list['Seq'].index(hseq)
                    seqdict[entry['Qry']] = seqlist.list['Seq'][qi] = (seqname,sequence)
                    self.printLog('#EDIT',seqdesc)
                    seqlist.list['Seq'].pop(hi)
                    seqdict.pop(entry['Hit'])
                    chromseq.insert(0,(seqname,sequence))
            else:   # Simple cut and stick join
                if self.getStrLC('JoinMerge') != 'end':
                    self.warnLog('JoinMerge=%s not recognised: will use "end" mode.' % self.getStr('JoinMerge'))
                # Trim Query and Hit
                if circularise:
                    if entry['Qry'] != entry['Hit']: raise ValueError
                    # Want to keep the end
                    qh = 'Qry'
                    chromseq.remove(seqdict[entry[qh]])
                    (seqname,sequence) = seqdict[entry[qh]]
                    si = seqlist.list['Seq'].index(seqdict[entry[qh]])
                    try: (sname,seqdesc) = string.split(seqname,maxsplit=1)
                    except: sname = seqname; seqdesc = ''

                    #?# Why did this previously care about RevComp#?#
                    #if seqdesc.startswith('RevComp'):
                    #    if entry['QryType'] == 'Start': y = entry['HitStart'] - 1
                    #    elif entry['QryType'] == 'End': y = entry['QryStart'] - 1
                    #    else: raise ValueError(entry['QryType'])
                    #    x = 1
                    #else:
                    #    if entry['QryType'] == 'Start': x = entry['QryEnd'] + 1
                    #    elif entry['QryType'] == 'End': x = entry['HitEnd'] + 1
                    #    else: raise ValueError(entry['QryType'])
                    #    y = len(sequence)

                    if entry['QryType'] == 'Start':     # Want to chop of the overlapping region from the query
                        x = entry['QryEnd'] + 1
                    elif entry['QryType'] == 'End':     # Want to chop of the overlapping region from the hit
                        x = entry['HitEnd'] + 1
                    else: raise ValueError(entry['QryType'])
                    y = len(sequence)                   # Always want to keep the end of the sequence

                    sequence = sequence[x-1:y]
                    self.printLog('#EDIT','%s -> Circle[Region %d to %d]' % (seqname,x,y))
                    seqname = '%s Circle[Region %d to %d|%s]' % (sname,x,y,seqdesc)
                    seqdict[entry[qh]] = seqlist.list['Seq'][si] = (seqname,sequence)
                    # Perform join
                    qseq = seqdict[entry['Qry']]
                    qi = seqlist.list['Seq'].index(qseq)
                    seqdict[entry['Qry']] = seqlist.list['Seq'][qi] = (seqname,sequence)
                    chromseq.insert(0,(seqname,sequence))
                    #!# Should check/add -A or -H suffix.
                else:
                    for qh in ['Qry','Hit']:
                        try: chromseq.remove(seqdict[entry[qh]])
                        except: self.warnLog('%s %s not found in %s chromseq dictionary?!' % (qh,entry[qh],chrom))
                        #!# This happens if manual editing of sequence between Assembly runs: make sure chromseq is regenerated at right point.
                        (seqname,sequence) = seqdict[entry[qh]]
                        si = seqlist.list['Seq'].index(seqdict[entry[qh]])

                        #?# Why did this previously care about RevComp#?#
                        #try: (sname,seqdesc) = string.split(seqname,maxsplit=1)
                        #except: sname = seqname; seqdesc = ''

                        #if seqdesc.startswith('RevComp'):
                        #    if entry['%sType' % qh] == 'Start': continue
                        #    elif entry['%sType' % qh] == 'End':
                        #        x = 1
                        #        y = entry['%sStart' % qh] - 1
                        #    else: raise ValueError(entry['%sType' % qh])
                        #else:
                        #    if entry['%sType' % qh] == 'Start':
                        #        x = entry['%sEnd' % qh] + 1
                        #        y = len(sequence)
                        #    elif entry['%sType' % qh] == 'End': continue
                        #    else: raise ValueError(entry['%sType' % qh])

                        # For the simple "cut" job, we trim the start of the second sequence and stick it on the first
                        if entry['%sType' % qh] == 'Start':             # This needs trimming
                            x = entry['%sEnd' % qh] + 1
                            y = len(sequence)
                        elif entry['%sType' % qh] == 'End': continue    # This stays full length
                        else: raise ValueError(entry['%sType' % qh])
                        # Edit only performed for "Start" sequence
                        sequence = sequence[x-1:y]
                        self.printLog('#EDIT','%s -> Region %d to %d.' % (seqname,x,y))
                        seqname = '%s (Region %d to %d)' % (seqname,x,y)
                        seqdict[entry[qh]] = seqlist.list['Seq'][si] = (seqname,sequence)

                    # Perform join
                    qseq = seqdict[entry['Qry']]
                    qacc = seqlist.seqAcc(qseq); qplus = 0
                    hseq = seqdict[entry['Hit']]
                    hacc = seqlist.seqAcc(hseq); hplus = 0
                    qend = True   # Join happens at end of query
                    if entry['QryType'] == 'End':
                        chopx = entry['HitEnd']  # Front chopped from Hit
                        joinseq = (seqdict[entry['Qry']],seqdict[entry['Hit']])
                        seqdesc = 'Join[%s & %s]' % (joinseq[0][0],joinseq[1][0])
                        seqname = '%s %s' % (string.split(joinseq[0][0])[0],seqdesc)
                        hplus = len(joinseq[0][1])  # Length that will need to be added to hacc tdb entries
                    elif entry['QryType'] == 'Start':
                        qend = False
                        chopx = entry['QryEnd']   # Front chopped from Qry
                        joinseq = (seqdict[entry['Hit']],seqdict[entry['Qry']])
                        seqdesc = 'Join[%s & %s]' % (joinseq[0][0],joinseq[1][0])
                        seqname = '%s %s' % (string.split(joinseq[1][0])[0],seqdesc)
                        qplus = len(joinseq[0][1])  # Length that will need to be added to qacc tdb entries
                    else: raise ValueError(entry['QryType'])
                    sequence = joinseq[0][1] + joinseq[1][1]
                    qi = seqlist.list['Seq'].index(qseq)
                    hi = seqlist.list['Seq'].index(hseq)
                    seqdict[entry['Qry']] = seqlist.list['Seq'][qi] = (seqname,sequence)
                    self.printLog('#EDIT',seqdesc)
                    seqlist.list['Seq'].pop(hi)
                    seqdict.pop(entry['Hit'])
                    # New sequence name is always the Qry (qryacc)
                    chromseq.insert(0,(seqname,sequence))
                    # Update tidy table
                    tdb.index('Contig',force=True)
                    for tentry in tdb.indexEntries('Contig',hacc):
                        if qend:
                            tentry['Start'] = max(1,tentry['Start']-chopx)
                            tentry['End'] = max(0,tentry['End']-chopx)
                            if not tentry['End']: tdb.dropEntry(tentry); continue   # Chopped off region!
                        tentry['Contig'] = qacc
                        tentry['Start'] += hplus
                        tentry['End'] += hplus
                        tentry['Desc'] = seqdesc
                        tentry['Seq'] = sequence
                    for tentry in tdb.indexEntries('Contig',qacc):
                        if not qend:
                            tentry['Start'] = max(1,tentry['Start']-chopx)
                            tentry['End'] = max(0,tentry['End']-chopx)
                            if not tentry['End']: tdb.dropEntry(tentry); continue   # Chopped off region!
                        tentry['Start'] += qplus
                        tentry['End'] += qplus
                        tentry['Desc'] = seqdesc
                        tentry['Seq'] = sequence
                    jentry = {'Contig':qacc,'Start':len(joinseq[0][1])+1}
                    jentry['End'] = jentry['Start'] + chopx - 1
                    jentry['Desc'] = seqdesc
                    jentry['Seq'] = sequence
                    tdb.addEntry(jentry)

            return self.assembleChrom(seqlist,chrom,chromseq)

        except: self.errorLog('%s.assembleChrom() error' % self.prog())
#########################################################################################################################
    def assemble(self,seqlist=None,seqfile=None):  ### Generates summary of statistics across multiple PAGSAT runs.
        '''Generates summary of statistics across multiple PAGSAT runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','assemble')
            if self.i() < 0:
                self.warnLog('Cannot perform PAGSTAT Manual Assembly if interactivity i<0: will use auto settings.')
            db = self.db()
            self.printLog('#~~#','## ~~~~~ PAGSTAT Manual Assembly ~~~~~ ##')
            if not seqlist:
                if not seqfile: raise ValueError('PAGSAT.assemble() needs seqlist or seqfile!')
                seqcmd = self.cmd_list + ['seqin=%s' % seqfile,'dna=T','autoload=F','seqmode=list']
                seqlist = rje_seqlist.SeqList(self.log,seqcmd)
                seqlist.loadSeq()
            elif not seqfile: seqfile = seqlist.getStr('SeqIn')
            #self.debug(seqfile)
            #x#assdir = rje.basePath(seqfile)  # Directory in which assembly lives
            assdir = self.getStr('AssembleDir')
            #self.debug(assdir)
            cgablamdir = rje.makePath('%sChromGABLAM/' % assdir)    # Where to generate GABLAM data.
            #self.debug(cgablamdir)
            rje.mkDir(self,cgablamdir)
            self.setStr({'CtgGABLAMDir':cgablamdir})
            ## ~ [0a] Initial split of sequences per chromosome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Replace this with custom assignment per chromosome, not per chromatid. (Or remove A/B from gene?)
            #splits = seqlist.splitSeq(basename=self.fileBase('CtgGABLAM'),splitseq='gene')
            splitseq = {}       # Dictionary of sequence tuples per chromosome
            seqdict = seqlist.seqNameDic()  # Dictionary of shortname to sequence
            for seq in seqlist.seqs():
                chrom = string.split(seq[0],'_')[0]
                if chrom[-1:] in ['H','A','B','N']: chrom = chrom[:-1]  # This should always be true but checking for future compatibility
                if chrom not in splitseq: splitseq[chrom] = []
                splitseq[chrom].append(seq)
            #self.debug(splitseq.keys())

            ### ~ [1] Cycle through each chromosome and try to iteratively assemble it ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            chrom2assemble = rje.sortKeys(splitseq)
            if 'Orphan' in chrom2assemble: chrom2assemble.remove('Orphan'); chrom2assemble.append('Orphan')
            #!# Make this a method for each chromosome?
            # Include localnfas output without alnseq
            for chrom in chrom2assemble:
                if chrom == 'Orphan' and self.yesNo('Process each Orphan separately?',default='Y'):
                    for seq in splitseq[chrom]:
                        sname = seqlist.shortName(seq)
                        self.assembleChrom(seqlist,sname,[seq])
                    continue
                self.assembleChrom(seqlist,chrom,splitseq[chrom])
                #?# Add chromatid assignment and merging/renaming (see below for ctid assignment) ??? Done elsewhere now?
                continue    # Trying to rationalise below in specific methods


            #!# Add option to look for joins in all contigs?



                    # If diploid=F, give option to (a) reject, or (b) combine chromatid B sequences with ctid A

            ## ~ [1x] Option to repeat assembly process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Regenerate seqlist.list['Seq'] from splitseq.values()
            # Resave with new filename (ask and check first)
            if self.i() > 0 and self.yesNo('Repeat manual assembly cycle?',default='N'): return self.assemble(seqlist)


            ## ~ [1h] Handle generation of diploid sequence: assemble track <C>ore then fill in track <D>iploid  ##


            ### ~ [2] Additional analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Option to run PAGSAT on new assembly
            # Option to run GABLAM -> SNPTable -> Snapper on new assembly


        except: self.errorLog('%s.assemble() error' % self.prog())
#########################################################################################################################
    ### <7> ### PAGSAT Comparison Methods                                                                               #
#########################################################################################################################
    def compare(self):  ### Generates summary of statistics across multiple PAGSAT runs.
        '''Generates summary of statistics across multiple PAGSAT runs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','compare')
            self.printLog('#~~#','## ~~~~~ PAGSTAT Compare mode (%d files) ~~~~~ ##' % len(self.list['Compare']))
            db = self.db()
            # This method essentially wants to read and combine the summary data from several runs, and then extract the
            # most useful information for assessing assessment quality, including:
            # - %coverage and %accuracy for assembly and reference.
            # - no. contigs
            # - optional gene/protein data if present
            ## ~ [0a] Load Reference Feature Table (if refgenome given) ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            repeatft = ['rRNA','mobile','LTR','centromere','telomere']   # self.list['RepeatFT']
            repeatft.sort()
            if repeatft: ftdb = self.ftdb()
            else: ftdb = None
            if ftdb:
                ftdb.dropEntriesDirect('feature',repeatft,inverse=True)
                self.printLog('#RPTFT','%s repeat features to exclude in "Uniq" outputs.' % rje.iStr(ftdb.entryNum()))
                if not ftdb.entryNum(): ftdb = None
            else:
                if not self.getStrLC('RefGenome'): self.printLog('#RPTFT','Cannot filter repeat features without refgenome=FILE.')
                if not repeatft: self.printLog('#RPTFT','Cannot filter repeat features without repeatft=LIST.')
            loc2chr = {}                # Will be a locus -> chromosome dictionary
            fragcov = self.list['FragCov']   # = [50,90,95,99] List of coverage thresholds to count (local table)
            chromcov = self.list['ChromCov'] # = [95,98,99] No. of chromosomes covered by a single contig (GABLAM table)
            fragcov.sort(); chromcov.sort()

            ### ~ [1] Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            compfields = ['Assembly','N','%AssCov','%AssAcc','%AssAlnCov','%AssAlnAcc','Multiplicity','Parsimony','%RefCov','%RefAcc','%GlobRefAcc','%RefAlnCov','%RefAlnAcc','%GlobAlnAcc','Missing','Errors','Extra','Duplicate','TreeLen','WtTreeLen']
            if ftdb: compfields += ['UniqCov','UniqCtg','UniqDup','RepeatFT']
            for chromx in chromcov: compfields.append('Chrom%d' % chromx)
            for fragx in fragcov:
                compfields.append('Frag%d' % fragx)
                if ftdb: compfields.append('UniqFrag%d' % fragx)
            compdb = db.addEmptyTable('compare',compfields,['Assembly'])   # Table of final comparison data
            if self.getBool('GeneSummary'): compdb.addFields(['%GeneCov','%GeneAcc'])#,'%GeneIntegrity'])
            if self.getBool('ProtSummary'): compdb.addFields(['%ProtCov','%ProtAcc'])#,'%ProtIntegrity'])
            for pfile in self.list['Compare']:
                #!# Can/should this whole process be moved into a function that can be run on a single dataset? #!#
                #!# Can then simply compile the datasets if found with the right headers.
                ## ~ [1a] Load and Process Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pbase = rje.baseFile(pfile,strip_path=True)
                if pbase.endswith('.Summary'): pbase = string.join(string.split(pbase,'.')[:-1],'.')     # Strip Summary
                basedir = rje.makePath('%s.PAGSAT/' % string.join(string.split(pbase,'.')[:-1],'.'))
                gabdir = rje.makePath('%s.GABLAM/' % string.join(string.split(pbase,'.')[:-1],'.'))
                self.setStr({'GABLAMDir':gabdir,'ResDir':basedir,'BaseBase':string.join(string.split(pbase,'.')[:-1],'.'),'CutBase':pbase})
                try: pdb = db.addTable(pfile,['Summary'],name=pbase,expect=True)
                except: self.errorLog('Cannot load PAGSAT Summary table "%s": check format' % pfile); continue
                pdb.dataFormat({'Length':'int','Coverage':'int','Identity':'int','Missing':'int','Errors':'int'})
                centry = {'Assembly':pbase}
                for entry in pdb.entries():
                    if entry['Summary'] == 'Reference':
                        centry['%RefCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Coverage
                        centry['%RefAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['%GlobRefAcc'] = 100.0 * entry['Identity'] / entry['Length']
                        centry['Missing'] = entry['Missing']
                        centry['Errors'] = entry['Errors']
                        centry['Duplicate'] = 0
                    if entry['Summary'] == 'ReferenceAlign':
                        centry['%RefAlnCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Coverage
                        centry['%RefAlnAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['%GlobAlnAcc'] = 100.0 * entry['Identity'] / entry['Length']
                        centry['Duplicate'] = 0
                    if entry['Summary'] == 'Assembly':
                        centry['%AssCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Validity
                        centry['%AssAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                        centry['Multiplicity'] = float(entry['Coverage']) / pdb.data('Reference')['Coverage']
                        centry['Parsimony'] = float(entry['Length']) / pdb.data('Reference')['Coverage']
                        centry['Extra'] = entry['Missing']
                        centry['N'] = entry['N']
                    if entry['Summary'] == 'AssemblyAlign':
                        centry['%AssAlnCov'] = 100.0 * entry['Coverage'] / entry['Length']     # Validity
                        centry['%AssAlnAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    if entry['Summary'] == 'Genes':
                        centry['%GeneCov'] = 100.0 * entry['Coverage'] / entry['Length']
                        centry['%GeneAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    if entry['Summary'] == 'Proteins':
                        centry['%ProtCov'] = 100.0 * entry['Coverage'] / entry['Length']
                        centry['%ProtAcc'] = 100.0 * entry['Identity'] / entry['Coverage']
                    #i# The Reciprocal Gene searches do not seem to very useful as summary data.
                    #if entry['Summary'] == 'Genes.Reciprocal':
                    #    centry['%GeneIntegrity'] = 100.0 * entry['Identity'] / entry['Length']
                    #if entry['Summary'] == 'Proteins.Reciprocal':
                    #    centry['%ProtIntegrity'] = 100.0 * entry['Identity'] / entry['Length']
                centry = compdb.addEntry(centry)
                #self.debug(centry)
                #self.debug('%s' % compdb.data(pbase))
                ## ~ [1b] Load and process CovPlot Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #MBG8150.SP16481.hcq.sgd.srt.1000.covplot.chrom.tdt
                cfile = string.replace(pfile,'Summary','covplot.chrom')
                if not rje.exists(cfile): self.warnLog('Could not locate %s' % cfile); continue
                try: cdb = db.addTable(cfile,['Chrom','Pos'],name='covplot',expect=True)
                except: self.errorLog('Cannot load PAGSAT chromosome coverage table "%s": check format' % cfile); continue
                cdb.dataFormat({'Pos':'int','ContigNum':'int','ChromNum':'int'})
                # ['Chrom','Pos','HitNum','ContigNum','Contigs','Class','ChromHit','ChromNum','RefChrom']
                covdata = {}    # Dict of {chrom:{pos:excess}}
                for entry in cdb.entries():
                    if entry['Chrom'] not in covdata: covdata[entry['Chrom']] = {}
                    covdata[entry['Chrom']][entry['Pos']] = entry['ContigNum'] - entry['ChromNum']
                centry['Duplicate'] = 0
                for chrom in covdata:
                    cpos = rje.sortKeys(covdata[chrom])
                    (x,i) = (0,0)
                    while i < len(cpos):
                        cx = covdata[chrom][cpos[i]]    # Contig hits - chrom hits
                        #i# I think this means that "Duplicate" only counts CNV increases on different chromosomes?
                        if cx > 0: centry['Duplicate'] +=  cx * (cpos[i] - x)
                        x = cpos[i]
                        i += 1
                ## ~ [1c] Unique coverage and duplication ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                bad_loc2chr = False
                if ftdb:        #!# NOTE: Unique stuff is not working!
                    centry['RepeatFT'] = string.join(repeatft,';')
                    # Note: locus=BK006937; Chrom=chrIII_S288C__BK006937.2
                    if not loc2chr:     # Only need to make once
                        for locus in ftdb.index('locus'):
                            loc2chr[locus] = None
                            for chrom in cdb.index('Chrom'):
                                if string.split(chrom,'__')[1].startswith(locus): loc2chr[locus] = chrom
                            if not loc2chr[locus]: self.warnLog('CovPlot missing %s?' % locus); bad_loc2chr = True
                        #self.debug(loc2chr)
                    # First, modify cdb to include R entries where Contig hits = Chrom hits = 0
                    for ft in ftdb.entries():
                        chrom = loc2chr[ft['locus']]
                        if not chrom or not cdb.indexEntries('Chrom',chrom): continue
                        fmin = ()   # (pos,centry) closest to 5' end of feature (to be at Pos=start-1)
                        fmax = ()   # (pos,centry) closest to 3' end of feature (to be at Pos=end+1)
                        for chrentry in cdb.indexEntries('Chrom',chrom)[0:]:
                            if chrentry['Pos'] < (ft['start']-1):
                                if not fmin or chrentry['Pos'] > fmin[0]: fmin = (chrentry['Pos'],chrentry)
                                continue    # No overlap, so keep
                            if chrentry['Pos'] > (ft['end']+1):
                                if not fmax or chrentry['Pos'] < fmax[0]: fmax = (chrentry['Pos'],chrentry)
                                continue    # No overlap, so keep
                            if not fmin or chrentry['Pos'] < fmin[0]: fmin = (chrentry['Pos'],chrentry)
                            if not fmax or chrentry['Pos'] > fmax[0]: fmax = (chrentry['Pos'],chrentry)
                            cdb.dropEntry(chrentry)
                        chrentry = fmin[1]
                        cdb.addEntry(rje.combineDict({'Chrom':chrom,'Pos':ft['start']-1},chrentry,overwrite=False))
                        cdb.addEntry({'Chrom':chrom,'Pos':ft['start'],'ContigNum':0,'Class':'R','ChromNum':0})
                        cdb.addEntry({'Chrom':chrom,'Pos':ft['end'],'ContigNum':0,'Class':'R','ChromNum':0})
                        chrentry = fmax[1]
                        cdb.addEntry(rje.combineDict({'Chrom':chrom,'Pos':ft['end']+1},chrentry,overwrite=False))
                    # Then calculate UniqDup as before
                    covdata = {}    # Dict of {chrom:{pos:excess}}
                    uniqdata = {}   # Dict of {chrom:{pos:class}}
                    uniqlen = totlen = 0; uniqcov = uniqctg = 0
                    for entry in cdb.entries():
                        if entry['Chrom'] not in covdata: covdata[entry['Chrom']] = {}; uniqdata[entry['Chrom']] = {}
                        try: covdata[entry['Chrom']][entry['Pos']] = entry['ContigNum'] - entry['ChromNum']
                        except: self.debug(entry)
                        uniqdata[entry['Chrom']][entry['Pos']] = entry['Class']
                    centry['UniqDup'] = 0
                    for chrom in covdata:
                        cpos = rje.sortKeys(covdata[chrom])
                        totlen += cpos[-1] - 1
                        uniqlen += cpos[-1] - 1
                        (x,i) = (0,0)
                        while i < len(cpos):
                            if uniqdata[chrom][cpos[i]] == 'R': uniqlen -= cpos[i] - x
                            elif uniqdata[chrom][cpos[i]] != 'N': uniqcov += cpos[i] - x    # Was C/U but I think
                            if uniqdata[chrom][cpos[i]] in ['C','U']: uniqctg += cpos[i] - x    # Need to annotate class ratings
                            cx = covdata[chrom][cpos[i]]    # Contig hits - chrom hits
                            if cx > 0: centry['UniqDup'] +=  cx * (cpos[i] - x)
                            x = cpos[i]
                            i += 1
                    centry['UniqCov'] = 100.0 * uniqcov / uniqlen
                    centry['UniqCtg'] = 100.0 * uniqctg / uniqlen
                db.deleteTable(cdb)
                ## ~ [1d] Load an process Tree file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                tfile = string.replace(pfile,'Summary.tdt','nwk')
                tree = rje_tree.Tree(self.log,self.cmd_list+['autoload=F'])
                tree.loadTree(tfile,postprocess=False)
                centry['TreeLen'] = rje.dp(tree.treeLen(),1)
                self.printLog('#TREE','Tree length = %s' % centry['TreeLen'])
                # Weighted tree length
                gfile = '%s.hitsum.tdt' % self.fileBase('GABLAM','Cut','Reference')
                if not rje.exists(gfile): self.warnLog('Could not locate %s' % gfile); continue
                try: gdb = db.addTable(gfile,['Qry'],name='rhitsum',expect=True)
                except: self.errorLog('Cannot load GABLAM table "%s": check format' % gfile); continue
                gdb.dataFormat({'Length':'int'})
                afile = '%s.hitsum.tdt' % self.fileBase('GABLAM','Cut','Assembly')
                if not rje.exists(afile): self.warnLog('Could not locate %s' % afile); continue
                try: adb = db.addTable(afile,['Qry'],name='ahitsum',expect=True)
                except: self.errorLog('Cannot load GABLAM table "%s": check format' % afile); continue
                adb.dataFormat({'Length':'int'})
                # Qry	Hit	Rank	Score	EVal	QryLen	HitLen
                # Qry = Reference; Hit = Assembly
                seqlen = {}
                for entry in gdb.entries() + adb.entries(): seqlen[entry['Qry']] = entry['Length']
                db.deleteTable(gdb); db.deleteTable(adb)
                centry['WtTreeLen'] = 0.0
                try:
                    treelen = 0.0
                    for branch in tree.branch:
                        blen = tree.pathLen([branch])
                        lens = []
                        for node in tree.branchClades(branch)[1]: lens.append(seqlen[node.shortName()])
                        treelen += blen * rje.meanse(lens)[0]
                    centry['WtTreeLen'] = rje.sf(treelen/1e6,4)
                    self.printLog('#WTLEN','Weighted tree length = %s' % centry['WtTreeLen'])
                except: self.errorLog('Problem generating weighted tree length!')
                ## ~ [1e] FragX and ChrX ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #!# NOTE: This needs to be improved with respect to checking options etc. #!#
                #lfile = '%s/%s.local.tdt' % (string.replace(pfile,'Summary.tdt','GABLAM'),pbase)
                # Not: MBG479.SP16499.hcq.qv20.sgd.srt.PAGSAT/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.GABLAM/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.local.tdt
                # -> MBG479.SP16499.hcq.qv20.sgd.srt.GABLAM/MBG479.SP16499.hcq.qv20.sgd.srt.L500ID800.Reference.local.tdt
                lfile = '%s.local.tdt' % self.fileBase('GABLAM','Cut','Reference')
                if not rje.exists(lfile): self.warnLog('Could not locate %s' % lfile); continue
                try: locdb = db.addTable(lfile,['Qry','Hit','AlnNum'],name='local',expect=True)
                except: self.errorLog('Cannot load local hit table "%s": check format' % lfile); continue
                locdb.dataFormat({'Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                for lentry in locdb.entries():
                    if lentry['SbjStart'] > lentry['SbjEnd']:
                        (lentry['SbjStart'],lentry['SbjEnd']) = (lentry['SbjEnd'],lentry['SbjStart'])
                # Use uniqlen and totlen calculated above to count number of local BLAST hits needed to exceed length thresholds
                #fragcov = [50,90,95,99]     # List of coverage thresholds to count (local table)
                for fragx in fragcov: centry['Frag%d' % fragx] = centry['UniqFrag%d' % fragx] = 0
                # Make lists of coverage (start, end), merging as required, sum up and compare to totlen
                # For uniqlen, start with a list of features before adding local hits
                covdict = {'Frag':{},'UniqFrag':{}}    # Dictionary of {chromosome:[(start,end)]
                covtot = {'Frag':{},'UniqFrag':{}}  # Dictionary of {Chromosome:total coverage}
                covlen = {'Frag':totlen,'UniqFrag':uniqlen}  # Dictionary of {Chromosome:total coverage}
                for qry in locdb.index('Qry'):
                    covdict['Frag'][qry] = []
                    covdict['UniqFrag'][qry] = []
                    covtot['Frag'][qry] = 0
                ucovdict = covdict['UniqFrag']
                for ft in ftdb.entries():
                    qry = loc2chr[ft['locus']]
                    if not qry: continue
                    ucovdict[qry].append((ft['start'],ft['end']))
                for qry in locdb.index('Qry'):
                    ucovdict[qry].sort()
                    x = 1
                    while x < len(ucovdict[qry]):
                        if ucovdict[qry][x][0] <= (ucovdict[qry][x-1][1] + 1):    # Merge
                            ucovdict[qry][x-1] = (ucovdict[qry][x-1][0],max(ucovdict[qry][x-1][1],ucovdict[qry][x][1]))
                            ucovdict[qry].pop(x)
                        else: x += 1
                    covtot['UniqFrag'][qry] = 0
                    for (i,j) in ucovdict[qry]: covtot['UniqFrag'][qry] += (j - i + 1)
                    #self.debug(ucovdict[qry])
                    #self.debug(covtot['UniqFrag'][qry])
                uniqlen = sum(covtot['UniqFrag'].values())
                covlen['UniqFrag'] = totlen - uniqlen
                # Add local hits in size order.
                hitx = 0    # Hit counter
                for lentry in locdb.sortedEntries('Identity',reverse=True):
                    hitx += 1
                    qry = lentry['Qry']
                    for c in ['Frag','UniqFrag']:
                        qfrag = covdict[c][qry]
                        qfrag.append((lentry['QryStart'],lentry['QryEnd']))
                        qfrag.sort()
                        x = 1
                        while x < len(qfrag):
                            if qfrag[x][0] <= (qfrag[x-1][1] + 1):    # Merge
                                qfrag[x-1] = (qfrag[x-1][0],max(qfrag[x-1][1],qfrag[x][1]))
                                qfrag.pop(x)
                            else: x += 1
                        covtot[c][qry] = 0
                        for (i,j) in qfrag: covtot[c][qry] += (j - i + 1)
                        # Assess coverage:
                        for fragx in fragcov:
                            if centry['%s%d' % (c,fragx)]: continue
                            if c == 'UniqFrag':
                                if 100.0 * (sum(covtot[c].values()) - uniqlen) / covlen[c] >= fragx: centry['%s%d' % (c,fragx)] = hitx
                            elif 100.0 * sum(covtot[c].values()) / covlen[c] >= fragx: centry['%s%d' % (c,fragx)] = hitx
                    if centry['Frag%d' % fragcov[-1]] and centry['UniqFrag%d' % fragcov[-1]]: break
                db.deleteTable(locdb)
                if bad_loc2chr: loc2chr = {}

                #chromcov = [50,95,98,99]    # No. of chromosomes covered by a single contig (GABLAM table)
                #gfile = string.replace(lfile,'local','gablam')
                gfile = '%s.gablam.tdt' % self.fileBase('GABLAM','Cut','Reference')
                if not rje.exists(gfile): self.warnLog('Could not locate %s' % gfile); continue
                try: gdb = db.addTable(gfile,['Qry','Hit'],name='gablam',expect=True)
                except: self.errorLog('Cannot load GABLAM table "%s": check format' % gfile); continue
                gdb.dropEntriesDirect('Rank',['1'],inverse=True)
                gxfield = 'Qry_AlnID'   # Could also use Qry_AlnLen
                gdb.dataFormat({gxfield:'float'})
                for chromx in chromcov:
                    centry['Chrom%d' % chromx] = 0
                    for gentry in gdb.entries():
                        if gentry[gxfield] >= chromx: centry['Chrom%d' % chromx] += 1

                db.deleteTable(gdb)

                self.debug(centry)
            compdb.saveToFile()



            pagfiles = []
            pagbase = []    # List of basefiles for PAGSAT results

            # Load in summary table, add assembly name and then combine with others
            # Reshape wide and then reshape long again!

            #Summary	HitNum	Length	Coverage	Identity	Positives	Missing	Errors	Perfect	N
            #Assembly	1173	13235834	13233590	13232765	13232765	2244	825	28	120
            #Reference	1190	12157104	12124470	12123837	12123837	32634	633	0	17
            #Self	268	12157104	12157104	12157104	12157104	0	0	17	17


            datatypes = ['Genes','Genes.Reciprocal','Proteins','Proteins.Reciprocal','Reference','Self']


        except: self.errorLog('%s.compare error' % self.prog())
#########################################################################################################################
### End of SECTION II: PAGSAT Class                                                                                     #
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
    try: PAGSAT(mainlog,cmd_list).run()

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
