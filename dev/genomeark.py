#!/usr/bin/python

# See below for name and description
# Copyright (C) 2016 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       GenomeARK
Description:  Genome Assembly, Read and Kmer analysis tool
Version:      0.2.0
Last Edit:    11/02/18
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The initial function of GenomeARK is to use information from kmer counts, read depth and sequence homology (using
    blastn+ for convenience) to assess assembled contigs for redundancy and duplication.

    General stats - have these fed as options or empirically determined - have a single pass and double pass method?
    - Can analyse just a subset of the assembly and/or over-ride empirical calculations.
    - RawDepth=X gives the haploid read depth generated from the raw data kmers
    - BAMDepth=LIST gives the haploid read depth for each BAM file. If -1 will use RawXDepth. If 0 will use median.


    Fork out each sequence
    - BLAST Hit Depth at different %identity cutoffs: need a localalnlen [250] and localalnid [40,90,99] settings
        -- Store best hits for single (non-self) hits?
        -- Filter out self hits with same start/end
    - Read Depths
        -- Can have multiple SAM/BAM files [Need a way to extract from BAM and convert to pileup?]
    - Assembly Kmer count
        -- jellyfish query ASSEMBLY.jf KMER -> KMER COUNT (e.g. AACACCCCAATGCTCGTAGCCCAAT 1)
    - Raw Kmer count
        -- jellyfish query ASSEMBLY.jf KMER

    Save and update results tables (unless force=T)
    - Sequences
        -- Calculate probable percentage: error, haploid, diploid, hetdup, homdup, repeat | av repeat CN
    - Regions
    - Windows for plots - Need to have sampling frequency and window size settings

    Option to mask sequences/sort based on ratings?
    - Lower case for predicted errors


    Initially just use one high %id cutoff for estimating assembly copy number based on BLAST hits - can see where this
    agrees or disagrees with assembly copy number based on assembly kmer counts.


    Read depths cannot say anything overall about the region but do say something about that specific contig w.r.t to how
    consistent it is with different genomic and assembly copy numbers.



    Want a separate repeat analysis pipeline using GABLAM results with lower %identity cutoff?



Internal Organisation:
    - rje_seq.SeqList stores the input contigs
    - rje_samtools.Samtools objects store the read depths (one per input file); popen BAM to SAM command (pull out contig)
        -- popen('samtools view input.bam contig') -> read into parse SAM code
    - rje_db.Database

Overview:
    Load contigs file.
    Setup output and identify processed contigs.
    Cycle through each contig and skip/fork.
    - skiplist=LIST will additionally skip certain contigs from analysis
    - contigs=LIST will give GenomeARK a subset of contigs to process.
    Compile results.

Setup:
    1. Load Genome Fasta file into SeqList object (seqmode=file).
    2. Generate list of contigs. = seqlist.names()
        - a. Compare to given contig list (if any), report errors and update.
        - b. Compare to skiplist (if any) and report errors.
    3. Backup and/or load main output.
        - a. Update skiplist if force=F

Output:
    ### ~ Main contig output table ~ ###
    * Contig = contig name.
    * Gap = Proportion predicted to be gapped (Ns in kmer)
    * Error = Proportion predicted to be sequencing/assembly error
    * Redundant = Proportion predicted to be redundant, i.e. more than 2x copies in assembly versus genome.
    * Hetero = Proportion predicted to be haploid coverage (two copies in assembly, one in genome)
    * Homo = Proportion predicted to be diploid coverage (one copy in assembly, one in genome)
    * Collapsed = one copy in assembly, 2+ copies in genome
    * CollapsedRep = 2+ copies in assembly but fewer than one copy in assembly per copy in genome
    * Repeat = Multiple copies in both genome and assembly (i.e. complex)
    * Unknown = Regions not mapping to any class with high probability
    * SCDep = Median BAM read depths for SC regions (AssemblyK=1, No local BLAST hits)
    * DCDep = Median BAM read depths for DC regions (AssemblyK=1 or 2, 1 local BLAST hits)
    * IDSC = Proportion single copy based on GABLAM (| sep list for each %id cutoff)
    * IDDC = Proportion double copy based on GABLAM (| sep list for each %id cutoff, n)
    * IDMC = Proportion multi copy based on GABLAM (| sep list for each %id cutoff, n)

    #!# Need to add classification of entire contig based on these criteria


    ### ~ Detailed contig position output table ~ ###
    * Contig = Contig name
    * Position = Position in contig
    * RawK = Kmer count in raw data
    * AssemblyK = Kmer count in assembly
    * NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
    * NContig = Number of non-self contigs with 1+ local BLAST hits (assembly copy number) (| separated list)
    * Dep = Normalised read depth for each BAM file (| separated list)
    * pErrK = Probability of Error (1X) based on raw kmer
    * pRedK = Probability of Redundant (0.5NX) based on raw kmer
    * pHetK = Probability of Heterozygous (1NX) based on raw kmer
    * pHomK = Probability of Homozygous (2NX) based on raw kmer
    * pRepK = Probability of Repeat (3NX+) based on raw kmer
    * RepCNK = Predicted repeat copy number based on raw kmer
    * pErrD = Probability of Error (1X) based on read depth
    * pRedD = Probability of Redundant (0.5NX) based on read depth
    * pHetD = Probability of Heterozygous (1NX) based on read depth
    * pHomD = Probability of Homozygous (2NX) based on read depth
    * pRepD = Probability of Repeat (3NX+) based on read depth


    #!# #?# Can we do something more sophisticated using BLAST local alignments and the kmers hit by a given kmer? #!#

    #!# Might want to only use one BAM file (PacBio if possible) - possible sequencing biases
    NOTE: Kmers based on raw sequencing data may have Illumina biases messing things up. Absence in raw sequencing data
    does not necessarily imply sequencing error: "Error" category may be an error in the *sequencing data* itself, not
    an error in the *assembly sequence*. ("Inconsistency" rather than "error".)

Fork plan:
    1. Save contig to query file.
    2. GABLAM against assembly (as fork then read in results?).
    3. Step through sequence, starting at (k+1)/2 then every sampling=X nucleotides.
        - a. Extract (k-1)/2 each side to make kmer. If kmer contains Ns then assign to gap.
        - b. Pull out kmer count for assembly and raw data.
        - c. Calculate BLAST depth (local and contig) at position for each %id cutoff.
        - d. Pull contig out of BAM files and calculate read depth for each position.
    4. Step through sequence again and:
        - a. Calculate probability of gap, error (1X), redundant (0.5NX), haploid (1NX), diploid (2NX), repeat (3NX+)
            - think about expanded: hetdup (3NX), homdup (4NX), repeat (5NX+) & across window (see below)
            - Based on raw kmers = copies in genome
        - b. Calculate copies in genome based on assembly kmers
        - c. Calculate genome:assembly ratio based on read depths
            - Probablilty of error (1X), redundant (0.5NX), haploid (1NX), diploid (2NX), multi (3NX)
        - d. Classify each based on high probability (plus "unknown") - need a probability threshold
        - e. Calculate copies in genome based on GABLAM local hits at each %identity
        - f. Calculate genome:assembly ratio from each BAM file using raw k and assembly k
        - g. Calculate genome copies using raw k and read depth (average over different BAM files)
    5. Combine statistics into per-contig summaries.
        - a. Classify contigs. #!# Need classification system #!#
    6. Optional masked/reformatted output of contigs.


    #?# What are we using the windows for?
    Windows are GOOD because they smooth out fluctuations.
    Windows are BAD because they cause problems at boundaries.
    Rather than windows, could use a simple flanking assignment system for the "Unknowns" if both neighbours agree?
    Or just not worry about it in the first instance?

    #i# For now, will sort the localidmin and use the biggest value for BLAST-based assembly copy number
    #localidmin=PERC  : Minimum local %identity thresholds for determing assembly copy numbers [99.0]


Probability assignments:


Commandline:
    ### ~ Input Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE    : Full genome assembly in fasta format. [assembly.fasta]
    maskfile=FASFILE : Repeat-Masked full genome assembly in fasta format. (Will use seqin if None) [None]
    bamlist=FILELIST : List of BAM files for read depth analysis (wildcards allowed). [ASSEMBLYBASE.bam]
    bamtypes=LIST    : Type of read data for BAM file (if known) (long/short) [None]
    bamdepth=NUMLIST : List of 1N X read depths for each BAM file. If -1 will use RawXDepth. If 0 will use median. [-1]
    assemblyjf=JFILE : Jellyfish file of kmers in the assembly file [ASSEMBLYBASE.jf]
    rawjf=JFILE      : Jellyfish file of kmers in the raw sequencing data [rawdata.jf]
    rawdepth=NUM     : Haploid read depth generated from the raw data kmers *Required* []
    contigs=LIST     : Optional subset of contigs to process. []
    skiplist=LIST    : Optional set of contigs to skip procesing []

    ### ~ Processing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    k=INT               : k used for the jellyfish runs (*must be odd*) [25]
    gcwin=INT           : Window size for GC content calculation (*must be odd*) [101]
    localmin=INT        : Minimum length of local alignment to output to local stats table [250]
    localidmin=NUMLIST  : List of minimum local %identity thresholds for local alignments [40,90,99]
    sampling=INT        : Sampling frequency for kmer and depth analysis [k]
    winsize=INT         : Window size up/downstream for sliding window analysis, e.g. 500 becomes 1001 bp window [500]
    blaste=X            : E-Value cut-off for BLAST searches (BLAST -e X) [1e-10]
    tophits=X           : Sets max number of BLAST hits returned (blastb and blastv) [1000]
    depcopy=X           : Type of BAM file (long/short/any) to use for depth-based CN estimation [any]
    recalculate=T/F     : Run special "recalculate" mode on existing output [False]
    depthcalc=T/F       : Run special "depth calculation" mode on first bamfile [False]

    ### ~ Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    basefile=X          : Sets base for output files [genomeark]
    contigfas=T/F       : Save processed contigs to fasta file. (Needed as reference for partial SAM reference) [False]
    localsam=T/F        : Save local hits data as SAM files in addition to TDT [False]
    pickup=T/F          : Whether to pick up and complete run (True) or assume fully (in)complete run (False) [False]

    ### ~ Forking Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    noforks=T/F     : Whether to avoid forks [False]
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
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
import rje, rje_db, rje_forker, rje_obj, rje_samtools, rje_seqlist
import rje_blast_V2 as rje_blast
import gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation with partial statstics output.
    # 0.1.0 - Initial "working" version with minimal unverified statistics and ratings.
    # 0.2.0 - Reorganised forking and recalculate mode.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [?] : Better name with Homology in there too
    # [ ] : Add some compilation statistics - median BLAST
    # [ ] : Add a randomcontig=X setting, that pulls out a random X contigs
    # [ ] : Add a separate kmer.tdt output for improved re-running.
    # [ ] : Add table headers to main sampling and contigark output.
    # [ ] : Add Recalculate option with samplingToContigARK() method.
    # [ ] : Replace bamlist=FILES with separate PacBio and Illumina BAM files (optional)? or...
    # [ ] : Add bamtypes=LIST (long/short) and depcopy=long/short/any to restrict depth-based CN estimation
    # [ ] : Implement a Forking/ subdirectory to store all the forked intermediates.
    # [Y] : Add pickup or checkpos option to check for complete files are just assume complete and skip
    # [Y] : Add skipping and appending to deplthcalc.
    # [ ] : Add forking and/or subsampling to depthcalc. (Sampling may cause issues where there are regions w/o reads)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('GenomeArk', '0.2.0', 'February 2018', '2018')
    description = 'Genome Assembly, Read and Kmer analysis tool'
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
### SECTION II: GenomeARK Class                                                                                         #
#########################################################################################################################
class GenomeARK(rje_forker.Forker):
    '''
    GenomeARK Class. Author: Rich Edwards (2015).

    Str:str
    - AssemblyJF=JFILE : Jellyfish file of kmers in the assembly file [ASSEMBLYBASE.jf]
    - DepCopy=X        : Type of BAM file (long/short/any) to use for depth-based CN estimation [any]
    - MaskFile=FASFILE : Repeat-Masked full genome assembly in fasta format. (Will use seqin if None) [None]
    - RawJF=JFILE      : Jellyfish file of kmers in the raw sequencing data [rawdata.jf]
    - SeqIn=FASFILE    : Full genome assembly in fasta format. [assembly.fasta]

    Bool:boolean
    - ContigFas=T/F       : Save processed contigs to fasta file. (Needed as reference for partial SAM reference) [False]
    - DepthCalc=T/F       : Run special "depth calculation" mode on first bamfile [False]
    - LocalSAM=T/F        : Save local hits data as SAM files in addition to TDT [False]
    - Pickup=T/F          : Whether to pick up and complete run (True) or assume fully (in)complete run (False) [False]
    - Recalculate=T/F     : Run special "recalculate" mode on existing output [False]

    Int:integer
    - GCWin=INT           : Window size for GC content calculation (*must be odd*) [101]
    - k =INT                 : k used for the jellyfish runs (*must be odd*) [25]
    - LocalMin=INT        : Minimum length of local alignment to output to local stats table [250]
    - Sampling=INT        : Sampling frequency for kmer and depth analysis [k]
    - TopHits=X           : Sets max number of BLAST hits returned (blastb and blastv) [1000]
    - WinSize=INT         : Window size up/downstream for sliding window analysis, e.g. 500 becomes 1001 bp window [500]

    Num:float
    - RawDepth=NUM     : Haploid read depth generated from the raw data kmers *Required* []

    File:file handles with matching str filenames
    - Out.Contigs = Main contigs table to append
    - Out.Sampling = Main position table to append

    List:list
    - BAMDepth=NUMLIST : List of 1N X read depths for each BAM file. If -1 will use RawXDepth. If 0 will use median. [-1]
    - BAMList=FILELIST     : List of BAM files for read depth analysis. [ASSEMBLYBASE.bam]
    - BAMTypes=LIST     : Type of read data for BAM file (if known) (long/short) []
    - Contigs=LIST     : Optional subset of contigs to process. []
    - LocalIDMin=NUMLIST  : List of minimum local %identity thresholds for local alignments [40,90,99]
    - SkipList=LIST    : Optional set of contigs to skip procesing []

    Dict:dictionary

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['AssemblyJF','DepCopy','MaskFile','RawJF','SeqIn','Out.Contigs','Out.Sampling']
        self.boollist = ['ContigFas','DepthCalc','LocalSAM','Pickup','Recalculate']
        self.intlist = ['GCWin','k','LocalMin','Sampling','TopHits','WinSize']
        self.numlist = ['RawDepth']
        self.filelist = ['Out.Contigs','Out.Sampling']
        self.listlist = ['BAMDepth','BAMList','BAMTypes','Contigs','LocalIDMin','SkipList']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'DepCopy':'any','RawJF':'rawdata.jf','SeqIn':'assembly.fasta'})
        self.setBool({'ContigFas':False,'LocalSAM':False,'Pickup':False,'Recalculate':False})
        self.setInt({'GCWin':101,'k':25,'LocalMin':250,'Sampling':0,'TopHits':1000,'WinSize':500})
        self.setNum({'RawDepth':0})
        self.list['LocalIDMin'] = [40,90,99]
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.baseFile('genomeark')
        self.cmd_list.insert(0,'tuplekeys=T')
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] ~ Read commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ###
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['DepCopy'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path
                self._cmdReadList(cmd,'file',['AssemblyJF','MaskFile','RawJF','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['ContigFas','DepthCalc','LocalSAM','Pickup','Recalculate'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['GCWin','k','LocalMin','Sampling','TopHits','WinSize'])   # Integers
                self._cmdReadList(cmd,'float',['RawDepth']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['BAMTypes','Contigs','SkipList'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'nlist',['BAMDepth'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['LocalIDMin'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['BAMList']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [2] ~ Tidy/update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.getStrLC('AssemblyJF'): self.setStr({'AssemblyJF':'%s.jf' % rje.baseFile(self.getStr('SeqIn'))})
        if not self.getStrLC('MaskFile'): self.setStr({'MaskFile':self.getStr('SeqIn')})
        if not self.list['BAMList']: self.list['BAMList'] = ['%s.bam' % rje.baseFile(self.getStr('SeqIn'))]
        if not rje.isOdd(self.getInt('GCWin')):
            self.setInt({'GCWin': self.getInt('GCWin')+1})
            self.warnLog('GCWin=INT must be odd: increased to %d' % self.getInt('GCWin'))
        if self.getInt('Sampling') <= 0: self.setInt({'Sampling':self.getInt('k')})
        if self.getBool('LocalSAM') and not self.getBool('ContigFas'):
            self.setBool({'ContigFas':True})
            self.printLog('#ARG','localsam=T: set contigfas=T for IGV reference.')
        self.cmd_list.append('debug=F')
        if not self.list['BAMTypes']: self.list['BAMTypes'] = ['None'] * len(self.list['BAMList'])
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        The initial function of GenomeARK is to use information from kmer counts, read depth and sequence homology (using
        blastn+ for convenience) to assess assembled contigs for redundancy and duplication.

        Three boolean options control the reuse of existing data:
        - If force=T, all data will be regenerated.
        - If pickup=F, it will be assumed that intermediate files are complete if discovered.
        - If recalculate=T, existing posark and contigark outputs will be regenerated, even if existing.

        Overview:
            Load contigs file.
            Setup output and identify processed contigs.
            Cycle through each contig and skip/fork.
            - skiplist=LIST will additionally skip certain contigs from analysis
            - contigs=LIST will give GenomeARK a subset of contigs to process.
            Compile results.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Set up objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.setup(): return False
            ## ~ [1b] ~ Special run modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DepthCalc'): return self.depthCalc()

            ### ~ [2] ~ Main forking of contigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Initial forking of contigs will do the basic kmer, depth and homology measurements
            #i# A second round of forking will perform the actual calculations and predictions
            #i# Recalculate mode will skip forking of contigs and re-process existing data
            stypes = ['rawk','assk','copy']
            for bami in range(len(self.list['BAMList'])):
                stypes.append('dep%d' % bami)
            if not self.getBool('Recalculate') and not self.forkContigs(stypes):
                raise ValueError('Forking did not end gracefully.')

            #!# Ideally, update the combine forking to be responsive to complete sets of data

            self.printLog('#COMB','Combining contig sampling data')
            if not self.getBool('Recalculate') and not self.forkContigs(['combine']):
                raise ValueError('Forking did not end gracefully.')
            ### ~ [3] ~ Assess or calculate single copy raw kmer frequencies and/or BAM depths ~~~~~~~~~~~~~~~~~~~~~~ ###
            #?# Add options for calculating SC numbers and also for doing the final calculations.
            #i# The '*.sampling.tdt' is loaded and used to estimate the kmer and read depths and
            self.singleCopyGenome()

            #!# Make the following into a method and use during setup to generate initial skiplist

            tfile = '%s.sampling.tdt' % self.baseFile()
            self.dict['ContigPos'] = {}     # Dictionary of filepos and Contig starts
            self.db().list['Tables'] = []
            self.printLog('#COMB','Paring contig sets from %s' % tfile)
            sdb = self.db().openTable(tfile,mainkeys=['Contig','Pos'],name='sampling')
            while sdb.obj['File']:
                fpos = sdb.obj['File'].tell()
                fcontig = sdb.readSet(['Contig'])
                if self.dev(): self.printLog('#SET','%s: %s = %d' % (tfile,fcontig[0],fpos))
                self.dict['ContigPos'][fcontig[0]] = fpos
            if sdb.obj['File']: sdb.obj['File'].close()
            self.db().deleteTable(sdb)

            ### ~ [4] ~ Collate/summarise results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Individual contig results will be collected during forking
            #i# Everything should be in a single table for all contigs but the region data will not be stored.
            #i# Each calculate fork should load the relevant sampling.tdt data using self.dict['ContigPos']
            if not self.forkContigs(['calculate']):
                raise ValueError('Forking did not end gracefully.')


        #?# Option to mask sequences/sort based on ratings?
        #?# - Lower case for predicted errors


        #?# Initially just use one high %id cutoff for estimating assembly copy number based on BLAST hits - can see where this
        #?# agrees or disagrees with assembly copy number based on assembly kmer counts.


        #?# Read depths cannot say anything overall about the region but do say something about that specific contig w.r.t to how
        #?# consistent it is with different genomic and assembly copy numbers.



        #?# Want a separate repeat analysis pipeline using GABLAM results with lower %identity cutoff?

        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''
        Main class setup method.
        Internal Organisation:
            - rje_seq.SeqList stores the input contigs
            - rje_samtools.Samtools objects store the read depths (one per input file); popen BAM to SAM command (pull out contig)
                -- popen('samtools view input.bam contig') -> read into parse SAM code
            - rje_db.Database
        General stats - have these fed as options or empirically determined - have a single pass and double pass method?
        - Optionally have different assembly and contigset
        - RawDepth=X gives the haploid read depth generated from the raw data kmers
        - BAMDepth=LIST gives the haploid read depth for each BAM file. If -1 will use RawXDepth. If 0 will use median.

        Three boolean options control the reuse of existing data:
        - If force=T, all data will be regenerated.
        - If pickup=F, it will be assumed that intermediate files are complete if discovered.
        - If recalculate=T, existing posark and contigark outputs will be regenerated, even if existing.

        NOTE: Need to update the following output descriptions.
        Output:
            ### ~ Main contig output table (*.contigark.tdt) ~ ###
            * Contig = contig name.
            * Len = contig length
            * Gap = Proportion predicted to be gapped (Ns in kmer)
            * Error = Proportion predicted to be sequencing/assembly error
            * Redundant = Proportion predicted to be redundant, i.e. more than 2x copies in assembly versus genome.
            * Hetero = Proportion predicted to be haploid coverage (two copies in assembly, one in genome)
            * Homo = Proportion predicted to be diploid coverage (one copy in assembly, one in genome)
            * Collapsed = one copy in assembly, 2+ copies in genome
            * CollapsedRep = 2+ copies in assembly but fewer than one copy in assembly per copy in genome
            * Repeat = Multiple copies in both genome and assembly (i.e. complex)
            * Unknown = Regions not mapping to any class with high probability
            * IDSC = Proportion single copy based on GABLAM (| sep list for each %id cutoff)
            * IDDC = Proportion double copy based on GABLAM (| sep list for each %id cutoff, n)
            * IDMC = Proportion multi copy based on GABLAM (| sep list for each %id cutoff, n)

            ### ~ Detailed contig position output table (*.sampling.tdt) ~ ###
            * Contig = Contig name
            * Position = Position in contig
            * RawK = Kmer count in raw data
            * AssemblyK = Kmer count in assembly
            * NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
            * NContig = Number of non-self contigs with 1+ local BLAST hits (assembly copy number) (| separated list)
            * Dep = Normalised read depth for each BAM file (| separated list)
            * pGapK = Probability of Gap (1+ N) based on raw kmer (either 0 or 1)
            * pErrK = Probability of Error (1X) based on raw kmer
            * pRedK = Probability of Redundant (0.5NX) based on raw kmer
            * pHetK = Probability of Heterozygous (1NX) based on raw kmer
            * pHomK = Probability of Homozygous (2NX) based on raw kmer
            * pRepK = Probability of Repeat (3NX+) based on raw kmer
            * RepCNK = Predicted repeat copy number based on raw kmer
            * pErrD = Probability of Error (1X) based on read depth
            * pRedD = Probability of Redundant (0.5NX) based on read depth
            * pHetD = Probability of Heterozygous (1NX) based on read depth
            * pHomD = Probability of Homozygous (2NX) based on read depth
            * pRepD = Probability of Repeat (3NX+) based on read depth
        '''
        try:### ~ [0] Setup General Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['LocalIDMin'].sort()
            if not self.list['LocalIDMin']: raise ValueError('Need at least one LocalIDMin values!')
            ## ~ [0a] Raw XDepth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while self.getNum('RawDepth') <= 0:
                if self.i() >= 0: self.setNum({'RawDepth':rje.getFloat('Enter haploid (1N) Xdepth for raw kmers',confirm=True)})
                else: raise ValueError('Need to set rawdepth=NUM.')
            ## ~ [0b] BAM XDepth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while len(self.list['BAMDepth']) < len(self.list['BAMList']):
                if self.list['BAMDepth'] and self.list['BAMDepth'][-1] == -1: self.list['BAMDepth'].append(-1)
                self.list['BAMDepth'].append(0)
            for i in range(len(self.list['BAMList'])):
                if self.list['BAMDepth'][i] == -1:
                    self.list['BAMDepth'][i] = self.getNum('RawDepth')
                    self.warnLog('Due to finite read lengths, XCoverage will not be the same for kmers and read mapping.')
                if self.list['BAMDepth'][i] == 0:
                    #!# Add a method for pulling out median read depth
                    self.printLog('#DEV','Development note: method for extracting median read depth not yet implemented.')
                    self.list['BAMDepth'][i] = rje.getFloat('Enter haploid (1N) Xdepth for %s' % self.list['BAMList'][i],confirm=True)

            ### ~ [1] Setup Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Load Genome Fasta file into SeqList object (seqmode=file) ~~~~~~~~~~~~~~~~~~~ ##
            seqcmd = ['seqmode=file','seqin=%s' % self.getStr('SeqIn'),'edit=F','autoload=T']
            seqlist = self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd)
            seqlist.seqNameDic()
            self.debug(self.getStr('MaskFile'))
            maskseq = self.obj['MaskSeq'] = rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd+['seqin=%s' % self.getStr('MaskFile')])
            maskseq.seqNameDic()
            contigs = seqlist.names()
            blast = rje_blast.blastObj(self.log,['blastf=F','gablamfrag=1']+self.cmd_list,type='dev')
            blast.formatDB(fasfile=self.getStr('SeqIn'),protein=False,force=False)
            ## ~ [1b] Load Genome Fasta file into SeqList object (seqmode=file) ~~~~~~~~~~~~~~~~~~~ ##
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db.setBasefile(self.baseFile())

            ### ~ [2] Check and/or generate contig list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Compare to given contig list (if any), report errors and update ~~~~~~~~~~~~~ ##
            if self.list['Contigs']:
                badcontigs = rje.listDifference(self.list['Contigs'],contigs)   #i# Returns the elements of list1 that are not found in list 2.
                if badcontigs: self.warnLog('%s of %s contigs=LIST contigs not found in assembly (seqin=%s)' % (rje.iLen(badcontigs),rje.iLen(self.list['Contigs']),self.getStr('SeqIn')),quitchoice=True)
                else: self.printLog('#CONTIG','%s contigs=LIST contigs found in assembly (seqin=%s)' % (rje.iLen(self.list['Contigs']),self.getStr('SeqIn')))
            else: self.list['Contigs'] = contigs[0:]
            ## ~ [2b] Compare to skiplist (if any) and report errors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['SkipList']:
                badcontigs = rje.listDifference(self.list['SkipList'],contigs)   #i# Returns the elements of list1 that are not found in list 2.
                if badcontigs: self.warnLog('%s of %s skiplist=LIST contigs not found in assembly (seqin=%s)' % (rje.iLen(badcontigs),rje.iLen(self.list['SkipList']),self.getStr('SeqIn')),quitchoice=True)
                else: self.printLog('#SKIP','%s skiplist=LIST contigs found in assembly (seqin=%s)' % (rje.iLen(self.list['SkipList']),self.getStr('SeqIn')))

            ### SPECIAL MODE QUICK EXIT ###
            if self.getBool('DepthCalc'): return True


            ### ~ [3] Backup and/or load main output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Three boolean options control the reuse of existing data:
            #i# - If force=T, all data will be regenerated.
            #i# - If pickup=F, it will be assumed that intermediate files are complete if discovered.
            #i# - If recalculate=T, existing posark and contigark outputs will be regenerated, even if existing.

            #>>> MODIFY >>>#

            ## ~ [3a] Main contig output table (*.contigark.tdt) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tfile = '%s.contigark.tdt' % self.baseFile()
            tfields = self.tableFields('contigark') # ['Contig','Len','Gap','Error','Redundant','Hetero','Homo','Collapsed','CollapsedRep','Repeat','Unknown','IDSC','IDDC','IDMC']
            tkeys = ['Contig']
            self.str['Out.Contigs'] = tfile
            if self.getBool('Recalculate') and rje.exists(tfile): rje.backup(self,tfile)
            if self.force() or not rje.exists(tfile):    #i# Make new file and table to append
                cdb = db.addEmptyTable('contigark',tfields,tkeys)
                cdb.saveToFile(tfile)

            #x# No longer update skiplist here
            elif not True:
                cdb = db.addTable(tfile,tkeys,name='contigark',expect=True)
                #i# Update skiplist if force=F
                skiplist = cdb.dataKeys()
                skiptxt = '%s loaded contigs added to %s skiplist=LIST contigs' % (rje.iLen(skiplist),rje.iLen(self.list['SkipList']))
                self.list['SkipList'] = rje.listUnion(self.list['SkipList'],skiplist)
                self.printLog('#SKIP','%s -> %s combined contigs to skip' % (skiptxt,rje.iLen(self.list['SkipList'])))

            #!# Need to sort out intermediate handling of sampling and posark files #!#
            #!# New forking process is not as simple as before #!#

            ## ~ [3b] Main contig sampling output table (*.sampling.tdt) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tfile = '%s.sampling.tdt' % self.baseFile()
            tfields = self.tableFields('sampling') #['Contig','Position','RawK','AssemblyK','NLocal','NContig','Dep','pErrK','pRedK','pHetK','pHomK','pRepK','RepCNK','pErrD','pRedD','pHetD','pHomD','pRepD']
            tkeys = ['Contig','Pos']
            self.str['Out.Sampling'] = tfile
            if self.force() or not rje.exists(tfile):    #i# Make new file and table to append
                sdb = db.addEmptyTable('sampling',tfields,tkeys)
                sdb.saveToFile(tfile)
            #else:
            #    sdb = db.addTable(tfile,tkeys,name='sampling',expect=True)
            ## ~ [3c] Main contig sampling output table (*.sampling.tdt) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tfile = '%s.posark.tdt' % self.baseFile()
            tfields = self.tableFields('posark') #['Contig','Position','RawK','AssemblyK','NLocal','NContig','Dep','pErrK','pRedK','pHetK','pHomK','pRepK','RepCNK','pErrD','pRedD','pHetD','pHomD','pRepD']
            tkeys = ['Contig','Pos']
            self.str['Out.PosARK'] = tfile
            if self.getBool('Recalculate') and rje.exists(tfile): rje.backup(self,tfile)
            if self.force() or not rje.exists(tfile):    #i# Make new file and table to append
                sdb = db.addEmptyTable('posark',tfields,tkeys)
                sdb.saveToFile(tfile)
            #else:
            #    sdb = db.addTable(tfile,tkeys,name='posark',expect=True)

            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def tableFields(self,output):   ### Returns fields for given output
        '''Returns fields for given output.'''
        if output == 'contigark':
            fields = ['Contig','Len','Gap','Error','Redundant','Hetero','Homo','Collapsed','CollapsedRep','Repeat','Unknown',
                      'SCDep','DCDep']
            for i in range(len(self.list['BAMList'])):
                if i:
                    fields.append('SCDep%d' % (i+1))
                    fields.append('DCDep%d' % (i+1))
            for locid in self.list['LocalIDMin']:
                fields.append('IDSC|%d' % locid)
                fields.append('IDDC|%d' % locid)
                fields.append('IDMC|%d' % locid)
        if output in ['sampling','posark']:
            fields = ['Contig','Pos','RawK','AssK']
            for locid in self.list['LocalIDMin']: fields.append('NLocal|%d' % locid)
            fields += ['NContig','Mask','Gap','GC','Dep']
            for i in range(len(self.list['BAMList'])):
                if i: fields.append('Dep%d' % (i+1))
            if output == 'posark':
                fields += ['pGapK','pErrK','pRedK','pHetK','pHomK','pRepK','RepCNK',
                           'pErrD','pRedD','pHetD','pHomD','pRepD','Rating','Confidence']
        try: return fields
        except: raise ValueError('Output %s not recognised!' % output)
#########################################################################################################################
    def forkContigs(self,stypes=[]):  ### Main forking method. Based of Forker run method.
        '''
        Main forking method. Cycle through each contig and skip/fork.
        - skiplist=LIST will additionally skip certain contigs from analysis
        - contigs=LIST will give GenomeARK a subset of contigs to process.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contigs = self.list['Contigs']
            skiplist = self.list['SkipList']
            if skiplist: skiplist = rje.listIntersect(contigs,self.list['SkipList'])
            self.list['ToFork'] = []
            cx = 0
            for contig in rje.listDifference(contigs,skiplist):
                cx += 1
                for ftype in stypes:
                    self.list['ToFork'].append((contig,ftype))
            self.printLog('#FORK','Fork list set up: %s of %s contigs (%s skipped)' % (rje.iStr(cx),rje.iLen(contigs),rje.iLen(skiplist)))
            forkx = len(self.list['ToFork'])    #!# Replace this with list of contigs.
            self.setBool({'LogFork':True})
            ### ~ [2] ~ Main forking cycle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Start initial set of contigs to forking
            self.setBool({'Child':False})   # Whether this is a child process
            self.list['Forked'] = []
            while len(self.list['Forked']) < self.getNum('Forks') and self.list['ToFork']: self.nextFork()
            #i# Cycle through forking cycle
            self.forking()  # Modify this to do the proper forking.
            self.printLog('#FORK','Forking of %s jobs completed.' % (rje.iStr(forkx)))
        except:  self.errorLog('%s.forkContigs() Error' % self.prog())
        if self.list['Forked']:
            self.warnLog('%s fork jobs remain unforked.' % rje.iLen(self.list['Forked']))
            return False
        return True
#########################################################################################################################
    def singleCopyGenome(self):     ### Method for assessing read depth profiles in GC bins for SC regions
        '''Method for assessing read depth profiles in GC bins for SC regions.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SCDEP','Single copy depth calculations not yet implemented')
        except:  self.errorLog('%s.singleCopyGenome() Error' % (self.prog()))
#########################################################################################################################
    ### <3> ### Forking Control Methods                                                                                 #
#########################################################################################################################
    #i# Under the GenomeARK redesign, each contig analysis element will be forked separately:
    #i# - contigForkKmer(fdict) will fork out kmer counts for a given jellyfish file.
    #i# - contigForkGABLAM(fdict) will fork a GABLAM of the contig versus a given genome file.
    #i# - contigForkDepth(fdict) will fork out read depth analysis for a given BAM file.
    #i# - contigForkCombine(fdict) will fork out combining the different outputs for a given contig and updating main tables.
    #i# In each case, fdict will store 'type':str (rawk/assk/copy/depN/combine) and (except combine) 'file':str.
#########################################################################################################################
    def needFork(self,fdict):  ### Checks for completeness of results file.
        '''Checks for completeness of results file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']
            forkfile = '%s.%s.tdt' % (fdict['FID'],fdict['type'])
            if self.force() or not rje.exists(forkfile): return True
            if rje.exists(forkfile) and not self.getBool('Pickup'): return False
            fdb = self.db().addTable(forkfile,['Contig','Pos'],name='fork')
            fdb.dataFormat({'Pos':'int'})
            k = self.getInt('k')
            w = (k-1)/2
            i = w  #i# Adjusted for 0<L, whereas Pos will be 1-L
            while (i+w) < len(sequence):
                pos = i + 1
                if (contig,pos) not in fdb.dataKeys():
                    self.db().deleteTable(fdb)
                    return True
                i += self.getInt('Sampling')
            self.db().deleteTable(fdb)
            return False
        except:  self.errorLog('%s.needFork(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigForkKmer(self,fdict):    ### Fork out kmer count extraction for a given jellyfish file.
        '''
        Fork out kmer count extraction for a given jellyfish file.
            -- jellyfish query JELLYFISHFILE.jf KMER -> KMER COUNT (e.g. AACACCCCAATGCTCGTAGCCCAAT 1)

        Kmer table:
            * Contig = Contig name
            * Pos = Position in contig
            * K = Kmer count in raw data

        fdict:
            * contig = contig
            * sequence = sequence
            * FID = basefile for output *.tdt & *.log (*.assemblyk or *.rawk)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']
            if fdict['type'] == 'rawk': jffile = self.getStr('RawJF'); ktype = 'raw'
            if fdict['type'] == 'assk': jffile = self.getStr('AssemblyJF'); ktype = 'assembly'
            ## ~ [1a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            sdb = db.addEmptyTable('kmer',['Contig','Pos','K'],['Contig','Pos'])

            ### ~ [4] Step through sequence, starting at (k+1)/2 then every sampling=X nucleotides ~~~~~~~~~~~~~~~~~~ ###
            k = self.getInt('k')
            w = (k-1)/2
            i = w  #i# Adjusted for 0<L, whereas Pos will be 1-L
            while (i+w) < len(sequence):
                #if self.dev():
                self.progLog('\r#KMERS','%s: |--%s--| (%s bp)' % (contig,rje.iStr(i),rje.iLen(sequence)),rand=0.1)
                pos = i + 1
                ## ~ [4a] Extract (k-1)/2 each side to make kmer. If kmer contains Ns then assign to gap. ~ ##
                x = i - w   # Start of kmer
                y = i + w   # End of kmer
                i += self.getInt('Sampling')
                kmer = sequence[x:y+1].upper()
                if len(kmer) != k: raise ValueError('kmer buggered!')
                nx = kmer.count('N')
                ## ~ [4b] Pull out kmer count for assembly and raw data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if nx:
                    #i# Give K a value of -1 if gapped region detected.
                    sentry = {'Contig':contig,'Pos':pos,'K':-1}
                else:
                    ajf = os.popen('jellyfish query %s %s' % (jffile,kmer)).read()
                    ajf = string.split(ajf)     #i# e.g. AACACCCCAATGCTCGTAGCCCAAT 1
                    if len(ajf[0]) != k: raise ValueError('jellyfish query problem for %s %s' % (ktype,kmer))
                    sentry = {'Contig':contig,'Pos':pos,'K':int(ajf[1])}
                ## ~ [4x] Add entry to table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sdb.addEntry(sentry)
            self.printLog('\r#KMERS','Processed %s %s kmers (%s bp)' % (rje.iStr(sdb.entryNum()),contig,rje.iLen(sequence)))

            #i# Save to file - name based on fork type - don't use ResFile or will be transferred at fork end
            if fdict['type'] == 'rawk': sdb.renameField('K','RawK')
            if fdict['type'] == 'assk': sdb.renameField('K','AssK')
            sdb.saveToFile('%s.%s.tdt' % (fdict['FID'],fdict['type']))

        except:  self.errorLog('%s.contigForkKmer(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigForkDepth(self,fdict):    ### Fork out BAM file read depth extraction for sampled kmers.
        '''
        Fork out BAM file read depth extraction for sampled kmers.

        Depth table:
            * Contig = Contig name
            * Pos = Position in contig
            * Dep = Normalised read depth for each BAM file (| separated list)

        fdict:
            * contig = contig
            * sequence = sequence
            * FID = basefile for output *.tdt & *.log (*.assemblyk or *.rawk)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']
            bami = string.atoi(fdict['type'][3:])   # Type is depN where N is the index of the BAMFile
            bamfile = self.list['BAMList'][bami]
            ## ~ [1a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            sdb = db.addEmptyTable('dep',['Contig','Pos','Dep'],['Contig','Pos'])

            ### ~ [1] Pull out BAM depth lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# NOTE: Calculations of state based on read depth moved to separate calculate method
            sam = rje_samtools.SAMtools(self.log,self.cmd_list)
            qc = [0]        # List of counts at different quality scores
            ri = 0          # RID Counter
            ridlist = []    # List of read IDs
            rdel = {}       # Dictionary of {rid:current deletion}
            dep = []        # List of depths across contig positions
            pcmd = 'samtools view %s %s -b | samtools mpileup -BQ0 -d10000000 -f %s -' % (bamfile,contig,self.getStr('SeqIn'))
            self.printLog('#BAM',pcmd)
            PILEUP = os.popen(pcmd)
            while PILEUP:
                #if self.dev():
                self.progLog('\r#BAMDEP','%s: %s %s depths' % (bamfile,rje.iLen(dep),contig),rand=0.001)
                line = PILEUP.readline()
                if not line: break
                bentry = sam.parsePileupLine(line,ri,ridlist,rdel,qc)
                while len(dep) < (bentry['Pos'] - 1): dep.append(0)
                if len(dep) == (bentry['Pos'] - 1): dep.append(bentry['N'])
                if ridlist: ri = max(ri,max(ridlist))
            PILEUP.close()
            self.printLog('\r#BAMDEP','%s: %s %s depths' % (bamfile,rje.iLen(dep),contig))

            ### ~ [4] Step through sequence, starting at (k+1)/2 then every sampling=X nucleotides ~~~~~~~~~~~~~~~~~~ ###
            k = self.getInt('k')
            w = (k-1)/2
            i = w  #i# Adjusted for 0<L, whereas Pos will be 1-L
            self.progLog('\r#DEPTH','%s %s depth for kmer sampling...' % (contig,bamfile))
            while (i+w) < len(sequence):
                #if self.dev():
                #self.progLog('\r#NLOC','%s: |--%s--| (%s bp)' % (contig,rje.iStr(i),rje.iLen(sequence)),rand=0.1)
                pos = i + 1
                ## ~ [4a] Extract (k-1)/2 each side to make kmer. If kmer contains Ns then assign to gap. ~ ##
                x = i - w   # Start of kmer
                y = i + w   # End of kmer
                i += self.getInt('Sampling')
                sentry = {'Contig':contig,'Pos':pos,'Dep':int(rje.mean(dep[x:y+1])+0.5)}
                sdb.addEntry(sentry)
            self.printLog('\r#DEPTH','%s %s depth for kmer sampling calculated.' % (contig,bamfile))


            #i# Save to file - name based on fork type - don't use ResFile or will be transferred at fork end
            sdb.renameField('Dep','Dep%d' % (bami+1))
            sdb.saveToFile('%s.%s.tdt' % (fdict['FID'],fdict['type']))

        except:  self.errorLog('%s.contigForkDepth(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigForkGABLAM(self,fdict):    ### Fork out kmer count extraction for a given jellyfish file.
        '''
        Fork out kmer count extraction for a given jellyfish file.
            -- jellyfish query JELLYFISHFILE.jf KMER -> KMER COUNT (e.g. AACACCCCAATGCTCGTAGCCCAAT 1)

        GABLAM Copy table:
            * Contig = Contig name
            * Position = Position in contig
            * NLocalX = Number of local hits at %identity cutoff
            * NContig = Median number of different contigs with BLAST homology across kmer
            * Mask = Proportion of kmer that is repeat masked
            * GC = GC frequency across GCWin
            * Gap = Proportion gaps (N) in unmasked assembly

        fdict:
            * contig = contig
            * sequence = sequence
            * maskseq = contig sequence for repeat masked assembly
            * FID = basefile for output *.tdt & *.log (*.assemblyk or *.rawk)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']    # Sequence of contig from SeqIn
            maskseq = fdict['maskseq']      # Sequence of contig from MaskFile
            basefile = fdict['FID']
            gcwin = (self.getInt('GCWin') - 1) / 2
            ## ~ [1a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            locfield = []
            for locidmin in self.list['LocalIDMin']:
                locfield.append('NLocal|%d' % locidmin)
            sdb = db.addEmptyTable('gab',['Contig','Pos']+locfield+['NContig','Mask','Gap','GC'],['Contig','Pos'])
            #i# Mask = 0/1 whether kmer has any N masking in MaskFile versus SeqIn
            #i# GC = GC content across GCWin

            ### ~ [1] Save contig to query file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qfile = '%s.contig.fas' % basefile
            open(qfile,'w').write('>%s\n%s\n' % (contig,sequence))

            ### ~ [2] GABLAM against assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# In future, look to fork this out and then read in results later, processing jf in meantime #!#
            gcmd = ['seqin=%s' % qfile,'searchdb=%s' % self.getStr('MaskFile'),'dna=T','blastp=blastn','basefile=%s' % basefile,
                    'fullres=T','hitsum=T','local=T','noforks=T','selfhit=F','selfsum=F','nrseq=F','qryacc=F','fullblast=T','blasta=4',
                    'localidmin=%f' % self.list['LocalIDMin'][0],'localmin=%d' % self.getInt('LocalMin'),'backups=F']
            gdefaults = ['blaste=1e-10','tophits=1000']
            if self.dev(): gdefaults.append('keepblast=T')
            else: gdefaults.append('keepblast=F')
            gablam.GABLAM(self.log,gdefaults+self.cmd_list+gcmd).gablam()
            ## ~ [2a] Generate GABLAM depth lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ldb = db.addTable('%s.local.tdt' % basefile,['Qry','Hit','AlnNum'],name='local',expect=True)
            ldb.dataFormat({'QryStart':'int','QryEnd':'int'})
            #!# ['Qry','Hit','AlnNum','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq']
            ldb.dropEntriesDirect('Hit',[contig])
            ldb.dropEntries(['Length<%d' % self.getInt('LocalMin')],logtxt='LocalMin',log=True)
            locdep = {}      # Number of non-self local BLAST hits (assembly repeat copy number?)
            gabdep = []      # Number of non-self contigs with 1+ local BLAST hits (assembly copy number)
            #i# Cycle through sorted percentage cutoffs
            hitlimx = 0     # Number of positions reaching the 'TopHits' limit
            ncontig = []      # Number of non-self contigs with 1+ local BLAST hits (assembly copy number)
            for i in range(len(sequence)): ncontig.append([])
            for locidmin in self.list['LocalIDMin']:
                #i# Reduce local data
                if locidmin > 0.0:
                    badidx = 0
                    for lentry in ldb.entries()[0:]:
                        if 100.0 * float(lentry['Identity']) / int(lentry['Length']) < locidmin:
                            badidx += 1
                            ldb.dropEntry(lentry)
                    self.printLog('#MINID','Dropped %s entries < LocalIDMin=%s%% -> %s entries' % (rje.iStr(badidx),rje.sf(locidmin),rje.iStr(ldb.entryNum())))
                #i# Set up depth lists
                locdep[locidmin] = nlocal = [0] * len(sequence)     # Number of non-self local BLAST hits (assembly repeat copy number?)
                #i# Calculate (can be made more efficient)
                lx = 0.0; ltot = ldb.entryNum()
                for lentry in ldb.entries():
                    self.progLog('\r#COPY','Calculating GABLAM-based copy number: %.1f%%' % (lx/ltot),rand=0.01); lx += 100.0
                    for i in range(lentry['QryStart']-1,lentry['QryEnd']):
                        nlocal[i] += 1

                    if locidmin == self.list['LocalIDMin'][-1]:
                        self.debug(lentry)
                        self.debug(range(lentry['QryStart']-1,lentry['QryEnd']))
                        for i in range(lentry['QryStart']-1,lentry['QryEnd']):
                            #self.bugPrint(i)
                            if lentry['Hit'] not in ncontig[i]:
                                ncontig[i].append(lentry['Hit'])
                                #if self.dev(): self.printLog('#GLOB','%d - %s - %d' % (i,string.join(ncontig[i]),len(ncontig[i])))
                        self.debug(ncontig[37800:37900])
                self.printLog('\r#COPY','Calculated GABLAM-based copy number at %s%%ID' % rje.sf(locidmin))
            for i in range(len(ncontig)):
                ncontig[i] = len(ncontig[i])
                if ncontig[i] >= (self.getInt('TopHits')-1): hitlimx += 1
            if hitlimx: self.warnLog('%s of %s positions hit the tophits=X BLAST limit!' % (rje.iStr(hitlimx),rje.iLen(ncontig)))
            #if self.dev(): self.printLog('#GLOBDEP',ncontig)

            ### ~ [4] Step through sequence, starting at (k+1)/2 then every sampling=X nucleotides ~~~~~~~~~~~~~~~~~~ ###
            k = self.getInt('k')
            w = (k-1)/2
            i = w  #i# Adjusted for 0<L, whereas Pos will be 1-L
            self.progLog('\r#COPY','%s GABLAM-based copy number for kmer sampling...' % contig)
            while (i+w) < len(sequence):
                #if self.dev():
                #self.progLog('\r#NLOC','%s: |--%s--| (%s bp)' % (contig,rje.iStr(i),rje.iLen(sequence)),rand=0.1)
                pos = i + 1
                sentry = {'Contig':contig,'Pos':pos,'Rating':'Unknown'}
                ## ~ [4a] Extract (k-1)/2 each side to make kmer. If kmer contains Ns then assign to gap. ~ ##
                x = i - w   # Start of kmer
                y = i + w   # End of kmer
                i += self.getInt('Sampling')
                kmer = sequence[x:y+1].upper()
                mask = maskseq[x:y+1].upper()
                sentry['Mask'] = mask.count('N') - kmer.count('N')
                if sentry['Mask']: sentry['Mask'] = sentry['Mask'] / float(k - kmer.count('N'))
                sentry['Gap'] = float(kmer.count('N')) / k
                gcseq = sequence[max(0,i-gcwin):i+gcwin]
                gcnonn = len(gcseq) - gcseq.count('N')
                if gcnonn: sentry['GC'] = float(gcseq.count('G')+gcseq.count('C')) / gcnonn
                else: sentry['GC'] = 0.0
                ## ~ [4c] Calculate BLAST depth (local and contig) at position for each %id cutoff ~ ##
                #* NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
                #* NContig = Number of non-self contigs with 1+ local BLAST hits (assembly copy number) (| separated list)
                # Cycle through sampling dictionaries and take mean
                for locidmin in self.list['LocalIDMin']:
                    sentry['NLocal|%d' % locidmin] = rje.median(locdep[locidmin][x:y+1],avtie=False)
                sentry['NContig'] = rje.median(ncontig[x:y+1],avtie=False)
                sdb.addEntry(sentry)
            self.printLog('\r#COPY','%s GABLAM-based copy number for kmer sampling calculated.' % contig)

            #i# Save to file - name based on fork type
            #X# sdb.saveToFile(fdict['ResFile'][fdict['type']]) - don't use ResFile or will be transferred at fork end
            sdb.saveToFile('%s.%s.tdt' % (fdict['FID'],fdict['type']))

        except:  self.errorLog('%s.contigForkGABLAM(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigForkCombine(self,fdict):    ### Combine data from different individual generator forks
        '''
        Combine data from different individual generator forks:
        - rawk
        - assk
        - copy
        - depN

        New sampling table:
            * Contig = Contig name
            * Pos = Position in contig
            * RawK = Kmer count in raw data
            * AssemblyK = Kmer count in assembly
            * NLocalX = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
            * NContig = Number of non-self contigs with 1+ local BLAST hits at highest %identity
            * Mask = Proportion of masked genes
            * Gap = Proportion of kmer that is gaps (N)
            * GC = GC content across GCWin
            * DepN = Normalised read depth for each BAM file
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            stypes = ['rawk','assk','copy']
            for bami in range(len(self.list['BAMList'])):
                stypes.append('dep%d' % bami)
            contig = fdict['contig']
            basefile = fdict['FID']
            ## ~ [1a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            sdb = None  #db.addEmptyTable('sampling',self.tableFields('sampling'),['Contig','Pos'])
            ## ~ [1b] Load input tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for stype in stypes:
                xdb = db.addTable('%s.%s.tdt' % (basefile,stype),mainkeys=['Contig','Pos'],name=stype,expect=True)
                xdb.dataFormat({'Pos':'int'})
                if not self.dev(): fdict['Cleanup'].append('%s.tdt' % stype)
                if not sdb: sdb = xdb; sdb.setStr({'Name':'sampling'})
                for field in xdb.fields():
                    if field not in sdb.fields(): sdb.addField(field)
            ## ~ [1c] Combine table data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for stype in stypes[1:]:
                self.progLog('#JOIN','Joining %s data on (Contig,Pos)...' % stype)
                xdb = self.db(stype)
                for skey in sdb.dataKeys():
                    sentry = sdb.data(skey)
                    for field in xdb.fields()[2:]: sentry[field] = xdb.data(skey)[field]
            self.printLog('#JOIN','Joined %s data on (Contig,Pos)' % string.join(stypes,'+'))

            ### ~ [5] Step through entries and calculate probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.renameField('Dep1','Dep')
            sdb.saveToFile(fdict['ResFile']['sampling'])

        except:  self.errorLog('%s.contigForkCombine(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigForkCalculate(self,fdict):    ### Perform actual rating calculations for a contig based on sampling data.
        '''
        Perform actual rating calculations for a contig based on sampling data:
        - BLAST Hit Depth at different %identity cutoffs: need a localalnlen [250] and localalnid [40,90,99] settings
            -- Store best hits for single (non-self) hits?
            -- Filter out self hits with same start/end
        - Read Depths
            -- Can have multiple SAM/BAM files [Need a way to extract from BAM and convert to pileup?]
        - Assembly Kmer count
            -- jellyfish query ASSEMBLY.jf KMER -> KMER COUNT (e.g. AACACCCCAATGCTCGTAGCCCAAT 1)
        - Raw Kmer count
            -- jellyfish query ASSEMBLY.jf KMER

        Save and update results tables (unless force=T)
        - Sequences
            -- Calculate probable percentage: error, haploid, diploid, hetdup, homdup, repeat | av repeat CN
        - Regions
        - Windows for plots - Need to have sampling frequency and window size settings

        Sampling table:
            * Contig = Contig name
            * Position = Position in contig
            * RawK = Kmer count in raw data
            * AssK = Kmer count in assembly
            * NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
            * NContig = Number of non-self contigs with 1+ local BLAST hits at highest %identity
            * Gap = Proportion of kmer that is gap
            * Dep = Normalised read depth for BAM file 1
            * DepN = Normalised read depth for BAM file 2+ [Optional]
            * pGapK = Probability of Gap (1+ N) based on raw kmer (either 0 or 1)
            * pErrK = Probability of Error (1X) based on raw kmer
            * pRedK = Probability of Redundant (0.5NX) based on raw kmer
            * pHetK = Probability of Heterozygous (1NX) based on raw kmer
            * pHomK = Probability of Homozygous (2NX) based on raw kmer
            * pRepK = Probability of Repeat (3NX+) based on raw kmer
            * RepCNK = Predicted repeat copy number based on raw kmer
            * pErrD = Probability of Error (1X) based on read depth
            * pRedD = Probability of Redundant (0.5NX) based on read depth
            * pHetD = Probability of Heterozygous (1NX) based on read depth
            * pHomD = Probability of Homozygous (2NX) based on read depth
            * pRepD = Probability of Repeat (3NX+) based on read depth

        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']
            basefile = fdict['FID']
            ## ~ [0a] Setup input/output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            #i# The sampling table contains the compiled raw counts etc. If running recalculate mode, this will need to
            #i# to be saved for a contig prior to forking.
            #!# Make sure required fields are renamed during recalculate!
            #X#sdb = db.addTable('%s.sampling.tdt' % basefile,['Contig','Pos'],name='sampling')
            tfile = '%s.sampling.tdt' % fdict['arkbase']
            sdb = self.db().openTable(tfile,mainkeys=['Contig','Pos'],name='sampling')
            sdb.obj['File'].seek(self.dict['ContigPos'][contig])
            sdb.readSet(['Contig'])
            if sdb.obj['File']: sdb.obj['File'].close()
            sdb.dataFormat({'Pos':'int','RawK':'int','AssK':'int','NContig':'int','Gap':'num','Mask':'num'})
            pdb = db.addEmptyTable('posark',self.tableFields('posark'),['Contig','Pos'])
            cdb = db.addEmptyTable('contigark',self.tableFields('contigark'),['Contig'])
            ## ~ [0b] Setup expected depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            kdepth = self.getNum('RawDepth')
            expected = {'K':{'Err':1,'Red':0.5*kdepth,'Het':kdepth,'Hom':2*kdepth,'Rep':3*kdepth}}
            deplist = ['Dep']      # BAM depth fields
            for i in range(len(self.list['BAMList'])):
                if i: dfield = 'Dep%d' % (i + 1); deplist.append(dfield)
                else: dfield = 'Dep'
                sdb.dataFormat({dfield:'num'})
                xdepth = self.list['BAMDepth'][i]   #i# This needs to be provided or calculated from full sampling data
                expected[dfield] = {'Err':1,'Red':0.5*xdepth,'Het':xdepth,'Hom':2*xdepth,'Rep':3*xdepth}

            ### ~ [1] Step through each sampled position of sequence and calculate stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = sdb.entryNum()
            for sentry in sdb.entries():
                self.progLog('#CALC','Calculating position stats for %s: %.1f%%' % (contig,sx/stot)); sx += 100.0
                pos = sentry['Pos']
                sentry['Rating'] = 'Unknown'
                for field in self.tableFields('sampling'):
                    if field not in sentry: sentry[field] = 0
                for p in ['Err','Red','Het','Hom','Rep']: #i# BAMFiles
                    sentry['p%sD' % p] = 1.0
                    sentry['p%sK' % p] = 0.0
                ## ~ [1a] Pull out kmer count for assembly and raw data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #* pGapK = Probability of Gap (1+ N) based on raw kmer (either 0 or 1)
                #* pErrK = Probability of Error (1X) based on raw kmer
                #* pRedK = Probability of Redundant (0.5NX) based on raw kmer
                #* pHetK = Probability of Heterozygous (1NX) based on raw kmer
                #* pHomK = Probability of Homozygous (2NX) based on raw kmer
                #* pRepK = Probability of Repeat (3NX+) based on raw kmer
                #* RepCNK = Predicted repeat copy number based on raw kmer
                sentry['pGapK'] = 0.0
                if sentry['Gap'] > 0:
                    sentry['pGapK'] = 1.0
                    sentry['RepCNK'] = 0.0
                else:
                    sentry['RepCNK'] = sentry['RawK'] / (2.0 * self.getNum('RawDepth'))
                    #i# Raw kmer copy probs
                    kprob = 0.0
                    for p in ['Err','Red','Het','Hom','Rep']:
                        sentry['p%sK' % p] = rje.logPoisson(sentry['RawK'],expected['K'][p],exact=True,callobj=self)
                        kprob += sentry['p%sK' % p]
                    #i# Convert to relative likelihood
                    if kprob:
                        for p in ['Err','Red','Het','Hom','Rep']:
                            sentry['p%sK' % p] /= kprob

                ## ~ [1b] Pull contig out of BAM files and calculate read depth for each position ~~ ##
                #* Dep = Normalised read depth for each BAM file (| separated list)
                #* pErrD = Probability of Error (1X) based on read depth
                #* pRedD = Probability of Redundant (0.5NX) based on read depth
                #* pHetD = Probability of Heterozygous (1NX) based on read depth
                #* pHomD = Probability of Homozygous (2NX) based on read depth
                #* pRepD = Probability of Repeat (3NX+) based on read depth
                # Cycle through sampling dictionaries and take mean

                #!# Redefine deplist using BAMTypes (short/long/None) and DepCopy (short/long/any)

                for dfield in deplist:
                    # Make a combined probability from all files
                    for p in ['Err','Red','Het','Hom','Rep']:
                        sentry['p%sD' % p] *= rje.logPoisson(int(sentry[dfield]),expected[dfield][p],exact=True,callobj=self)
                #i# Convert to relative likelihood
                dprob = 0.0
                for p in ['Err','Red','Het','Hom','Rep']: dprob += sentry['p%sD' % p]
                if dprob:
                    for p in ['Err','Red','Het','Hom','Rep']: sentry['p%sD' % p] /= dprob

                ## ~ [1c] Add Rating of section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Classify each based on high probability (plus "unknown") - need a probability threshold
                #* Unknown = Regions not mapping to any class with high probability = DEFAULT
                #* Gap = Proportion predicted to be gapped (Ns in kmer)
                #i# kmer has Ns and no contig matches -> Gap
                ncopy = sentry['NContig']
                if sentry['pGapK'] and ncopy == 0:
                    sentry['Rating'] = 'Gap'; sentry['Confidence'] = 'High'
                #* Error = Proportion predicted to be sequencing/assembly error
                if sentry['pErrK'] > 0.5 and ncopy == 0:
                    sentry['Rating'] = 'Error'; sentry['Confidence'] = 'High'
                #* Hetero = Proportion predicted to be haploid coverage (two copies in assembly, one in genome)
                #i# High kmer hetero probability and one contig matches
                elif sentry['pHetK'] > 0.5 and ncopy == 1:
                    sentry['Rating'] = 'Hetero'; sentry['Confidence'] = 'High'
                #* Homo = Proportion predicted to be diploid coverage (one copy in assembly, one in genome)
                #i# High kmer homo probability and no contig matches
                elif sentry['pHomK'] > 0.5 and ncopy == 0:
                    sentry['Rating'] = 'Homo'; sentry['Confidence'] = 'High'
                #* Collapsed = one copy in assembly, 2+ copies in genome
                elif sentry['pRepK'] > 0.5 and ncopy == 1:
                    sentry['Rating'] = 'Collapsed'; sentry['Confidence'] = 'High'
                #* CollapsedRep = 2+ copies in assembly but fewer than one copy in assembly per copy in genome
                elif sentry['pRepK'] > 0.5 and sentry['RepCNK'] > max(sentry['AssK'],ncopy+1):
                    sentry['Rating'] = 'CollapsedRep'; sentry['Confidence'] = 'High'
                #* Repeat = Multiple copies in both genome and assembly (i.e. complex)
                elif sentry['pRepK'] > 0.5 and ncopy > 0:
                    sentry['Rating'] = 'Repeat'; sentry['Confidence'] = 'High'
                #* Redundant = Proportion predicted to be redundant, i.e. more than 2x copies in assembly versus genome.
                elif sentry['RepCNK'] < 0.5 * min(sentry['AssK'],ncopy+1):
                    sentry['Rating'] = 'Redundant'; sentry['Confidence'] = 'High'

                #!# Other ideas/considerations
                ## ~ [- e. Calculate copies in genome based on GABLAM local hits at each %identity
                ## ~ [- f. Calculate genome:assembly ratio from each BAM file using raw k and assembly k
                ## ~ [- g. Calculate genome copies using raw k and read depth (average over different BAM files)

                ## ~ [1d] Add entry to table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pdb.addEntry(sentry)
            self.printLog('#CALC','Calculation of position stats for %s complete: %s sampled positions' % (contig,rje.iStr(stot)))

            ### ~ [2] Step through entries and calculate probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb.saveToFile(fdict['ResFile']['posark'])


            ### ~ [3] Combine statistics into per-contig summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #* Contig = contig name.
            #* Len = contig length
            #* SCDep = Median BAM read depths for SC regions (AssK=1, No local BLAST hits)
            #* DCDep = Median BAM read depths for DC regions (AssK=1 or 2, 1 local BLAST hits)
            #* IDSC = Proportion single copy based on GABLAM (| sep list for each %id cutoff)
            #* IDDC = Proportion double copy based on GABLAM (| sep list for each %id cutoff, n)
            #* IDMC = Proportion multi copy based on GABLAM (| sep list for each %id cutoff, n)
            ## ~ [5a] Convert sampling data into contig stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Split the contig BLAST depth counts and convert into copy numbers per ID cutoff.
            #i# Rejoin and tidy the ID* fields after compression and taking means.
            for field in ['Error','Redundant','Hetero','Homo','Collapsed','CollapsedRep','Repeat','Unknown']:
                sdb.addField(field,evalue=0)
            for i in range(len(self.list['LocalIDMin'])):
                locidmin = self.list['LocalIDMin'][i]
                for field in ['IDSC','IDDC','IDMC']:
                    sdb.addField('%s|%s' % (field,locidmin),evalue=0)
                for sentry in sdb.entries():
                    ncopy = sentry['NLocal|%d' % locidmin]
                    if ncopy < 1: field = 'IDSC'
                    elif ncopy == 1: field = 'IDDC'
                    else: field = 'IDMC'
                    sentry['%s|%s' % (field,locidmin)] = 1.0
                    sentry[sentry['Rating']] = 1.0
            #i# Keep most robust analysis as the BLAST-based assembly copy number
            #i# Compare to the Assembly Kmer copy number? Or Raw:Assembly Kmer ratio

            depfields = ['Pos','Copy']
            for i in range(len(self.list['BAMList'])): depfields.append('BAM%d' % i)
            depdb = self.db().addEmptyTable('%s.dep' % contig,depfields,['Pos'],log=True)
            locidmin = self.list['LocalIDMin'][0]
            for sentry in sdb.entries():
                ncopy = sentry['NLocal|%d' % locidmin]
                dentry = {'Pos':sentry['Pos']}

                for i in range(len(self.list['BAMList'])):
                    dfield = 'Dep'
                    if i: dfield = 'Dep%d' % (i+1)
                    dentry[dfield] = int(sentry[dfield])
                if sentry['AssK'] == 1 and ncopy == 0:         # Solid single copy
                    dentry['Copy'] = 'SC'
                    depdb.addEntry(dentry)
                elif 0 < sentry['AssK'] < 3 and ncopy < 2:     # Probable double copy
                    dentry['Copy'] = 'DC'
                    depdb.addEntry(dentry)
            depdb.compress(['Copy'],default='median')
            #bamdep = {}
            #for dentry in depdb.entries():
            #    dep = []
            #    for i in range(len(self.list['BAMList'])):
            #        dep.append('%.1f' % dentry['BAM%d' % i])
            #    bamdep[dentry['Copy']] = string.join(dep,'|')

            #sdb.dropFields(['NLocal','Dep'])

            ## ~ [5b] Compress the sampling entries, taking mean values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.compress(['Contig'])
            sentry = sdb.entries()[0] #i# should only be one!
            #for field in ['IDSC','IDDC','IDMC']:
                #sentry[field] = []
                #for locidmin in self.list['LocalIDMin']:
                #    sentry[field].append('%.3f' % sentry['%s|%s' % (field,locidmin)])
                #sentry[field] = string.join(sentry[field],'|')
            sentry['Len'] = len(sequence)
            for dentry in depdb.entries():
                for i in range(len(self.list['BAMList'])):
                    if i: sentry['%sDep' % dentry['Copy']] = dentry['BAM%d' % i]
                    else: sentry['%sDep%d' % (dentry['Copy'],i+1)] = dentry['BAM%d' % i]
                #bamdep[dentry['Copy']] = string.join(dep,'|')

            #sentry['SCDep'] = bamdep['SC']
            #sentry['DCDep'] = bamdep['DC']

            ## ~ [- a. Classify contigs. #!# Need classification system #!#
            ### ~ [6. Optional masked/reformatted output of contigs.

            cdb.addEntry(sentry)
            cdb.saveToFile(fdict['ResFile']['contigark'])
        except:  self.errorLog('%s.contigForkCalculate(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def contigFork(self,fdict):    ### Perform actual forked analysis for a contig
        '''
        Perform actual forked analysis for a contig.
        Fork out each sequence
        - BLAST Hit Depth at different %identity cutoffs: need a localalnlen [250] and localalnid [40,90,99] settings
            -- Store best hits for single (non-self) hits?
            -- Filter out self hits with same start/end
        - Read Depths
            -- Can have multiple SAM/BAM files [Need a way to extract from BAM and convert to pileup?]
        - Assembly Kmer count
            -- jellyfish query ASSEMBLY.jf KMER -> KMER COUNT (e.g. AACACCCCAATGCTCGTAGCCCAAT 1)
        - Raw Kmer count
            -- jellyfish query ASSEMBLY.jf KMER

        Save and update results tables (unless force=T)
        - Sequences
            -- Calculate probable percentage: error, haploid, diploid, hetdup, homdup, repeat | av repeat CN
        - Regions
        - Windows for plots - Need to have sampling frequency and window size settings

        Sampling table:
            * Contig = Contig name
            * Position = Position in contig
            * RawK = Kmer count in raw data
            * AssemblyK = Kmer count in assembly
            * NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
            * NContig = Number of non-self contigs with 1+ local BLAST hits at highest %identity
            * Dep = Normalised read depth for each BAM file (| separated list)
            * pGapK = Probability of Gap (1+ N) based on raw kmer (either 0 or 1)
            * pErrK = Probability of Error (1X) based on raw kmer
            * pRedK = Probability of Redundant (0.5NX) based on raw kmer
            * pHetK = Probability of Heterozygous (1NX) based on raw kmer
            * pHomK = Probability of Homozygous (2NX) based on raw kmer
            * pRepK = Probability of Repeat (3NX+) based on raw kmer
            * RepCNK = Predicted repeat copy number based on raw kmer
            * pErrD = Probability of Error (1X) based on read depth
            * pRedD = Probability of Redundant (0.5NX) based on read depth
            * pHetD = Probability of Heterozygous (1NX) based on read depth
            * pHomD = Probability of Homozygous (2NX) based on read depth
            * pRepD = Probability of Repeat (3NX+) based on read depth

        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            contig = fdict['contig']
            sequence = fdict['sequence']
            basefile = fdict['FID']
            ## ~ [1a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            cdb = db.addEmptyTable('contigark',self.tableFields('contigark'),['Contig'])
            sdb = db.addEmptyTable('sampling',self.tableFields('sampling'),['Contig','Pos'])

            ### ~ [1] Save contig to query file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qfile = '%s.contig.fas' % basefile
            open(qfile,'w').write('>%s\n%s\n' % (contig,sequence))

            ### ~ [2] GABLAM against assembly ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# In future, look to fork this out and then read in results later, processing jf in meantime #!#
            gcmd = ['seqin=%s' % qfile,'searchdb=%s' % self.getStr('SeqIn'),'dna=T','blastp=blastn','basefile=%s' % basefile,
                    'fullres=T','hitsum=T','local=T','noforks=T','selfhit=F','selfsum=F','nrseq=F','qryacc=F','fullblast=T','blasta=4',
                    'localidmin=%f' % self.list['LocalIDMin'][0],'localmin=%d' % self.getInt('LocalMin'),'backups=F']
            gdefaults = ['blaste=1e-10','tophits=1000']
            if self.dev(): gdefaults.append('keepblast=T')
            else: gdefaults.append('keepblast=F')
            gablam.GABLAM(self.log,gdefaults+self.cmd_list+gcmd).gablam()
            ## ~ [2a] Generate GABLAM depth lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ldb = db.addTable('%s.local.tdt' % basefile,['Qry','Hit','AlnNum'],name='local',expect=True)
            ldb.dataFormat({'QryStart':'int','QryEnd':'int'})
            #!# ['Qry','Hit','AlnNum','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq']
            ldb.dropEntriesDirect('Hit',[contig])
            ldb.dropEntries(['Length<%d' % self.getInt('LocalMin')],logtxt='LocalMin',log=True)
            locdep = {}      # Number of non-self local BLAST hits (assembly repeat copy number?)
            gabdep = []      # Number of non-self contigs with 1+ local BLAST hits (assembly copy number)
            #i# Cycle through sorted percentage cutoffs
            hitlimx = 0     # Number of positions reaching the 'TopHits' limit
            ncontig = []      # Number of non-self contigs with 1+ local BLAST hits (assembly copy number)
            for i in range(len(sequence)): ncontig.append([])
            for locidmin in self.list['LocalIDMin']:
                #i# Reduce local data
                if locidmin > 0.0:
                    badidx = 0
                    for lentry in ldb.entries()[0:]:
                        if 100.0 * float(lentry['Identity']) / int(lentry['Length']) < locidmin:
                            badidx += 1
                            ldb.dropEntry(lentry)
                    self.printLog('#MINID','Dropped %s entries < LocalIDMin=%s%% -> %s entries' % (rje.iStr(badidx),rje.sf(locidmin),rje.iStr(ldb.entryNum())))
                #i# Set up depth lists
                locdep[locidmin] = nlocal = [0] * len(sequence)     # Number of non-self local BLAST hits (assembly repeat copy number?)
                #i# Calculate (can be made more efficient)
                lx = 0.0; ltot = ldb.entryNum()
                for lentry in ldb.entries():
                    self.progLog('\r#COPY','Calculating GABLAM-based copy number: %.1f%%' % (lx/ltot),rand=0.01); lx += 100.0
                    for i in range(lentry['QryStart']-1,lentry['QryEnd']):
                        nlocal[i] += 1

                    if locidmin == self.list['LocalIDMin'][-1]:
                        self.debug(lentry)
                        self.debug(range(lentry['QryStart']-1,lentry['QryEnd']))
                        for i in range(lentry['QryStart']-1,lentry['QryEnd']):
                            #self.bugPrint(i)
                            if lentry['Hit'] not in ncontig[i]:
                                ncontig[i].append(lentry['Hit'])
                                #if self.dev(): self.printLog('#GLOB','%d - %s - %d' % (i,string.join(ncontig[i]),len(ncontig[i])))
                        self.debug(ncontig[37800:37900])
                self.printLog('\r#COPY','Calculated GABLAM-based copy number at %s%%ID' % rje.sf(locidmin))
            for i in range(len(ncontig)):
                ncontig[i] = len(ncontig[i])
                if ncontig[i] >= (self.getInt('TopHits')-1): hitlimx += 1
            if hitlimx: self.warnLog('%s of %s positions hit the tophits=X BLAST limit!' % (rje.iStr(hitlimx),rje.iLen(ncontig)))
            #if self.dev(): self.printLog('#GLOBDEP',ncontig)

            ### ~ [3] Pull out BAM depth lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sam = rje_samtools.SAMtools(self.log,self.cmd_list)
            bamdep = {}     # Dictionary of bamfile:list of read depths
            kdepth = self.getNum('RawDepth')
            expected = {'K':{'Err':1,'Red':0.5*kdepth,'Het':kdepth,'Hom':2*kdepth,'Rep':3*kdepth}}
            for i in range(len(self.list['BAMList'])):
                bamfile = self.list['BAMList'][i]
                xdepth = self.list['BAMDepth'][i]
                expected[bamfile] = {'Err':1,'Red':0.5*xdepth,'Het':xdepth,'Hom':2*xdepth,'Rep':3*xdepth}
                qc = [0]        # List of counts at different quality scores
                ri = 0          # RID Counter
                ridlist = []    # List of read IDs
                rdel = {}       # Dictionary of {rid:current deletion}
                dep = []        # List of depths across contig positions
                pcmd = 'samtools view %s %s -b | samtools mpileup -BQ0 -d10000000 -f %s -' % (bamfile,contig,self.getStr('SeqIn'))
                self.printLog('#BAM',pcmd)
                PILEUP = os.popen(pcmd)
                while PILEUP:
                    #if self.dev():
                    self.progLog('\r#BAMDEP','%s: %s %s depths' % (bamfile,rje.iLen(dep),contig),rand=0.001)
                    line = PILEUP.readline()
                    if not line: break
                    bentry = sam.parsePileupLine(line,ri,ridlist,rdel,qc)
                    while len(dep) < (bentry['Pos'] - 1): dep.append(0)
                    if len(dep) == (bentry['Pos'] - 1): dep.append(bentry['N'])
                    if ridlist: ri = max(ri,max(ridlist))
                PILEUP.close()
                bamdep[bamfile] = dep[0:]
                self.printLog('\r#BAMDEP','%s: %s %s depths' % (bamfile,rje.iLen(dep),contig))

            ### ~ [4] Step through sequence, starting at (k+1)/2 then every sampling=X nucleotides ~~~~~~~~~~~~~~~~~~ ###
            k = self.getInt('k')
            w = (k-1)/2
            i = w  #i# Adjusted for 0<L, whereas Pos will be 1-L
            px = 0
            while (i+w) < len(sequence):
                #if self.dev():
                self.progLog('\r#KMERS','%s: |--%s--| (%s bp)' % (contig,rje.iStr(i),rje.iLen(sequence)),rand=0.1)
                pos = i + 1; px +=1
                sentry = {'Contig':contig,'Pos':pos,'Rating':'Unknown'}
                for field in self.tableFields('sampling'):
                    if field not in sentry: sentry[field] = 0
                for p in ['Err','Red','Het','Hom','Rep']:
                    #i# BAMFiles
                    sentry['p%sD' % p] = 1.0
                ## ~ [4a] Extract (k-1)/2 each side to make kmer. If kmer contains Ns then assign to gap. ~ ##
                x = i - w   # Start of kmer
                y = i + w   # End of kmer
                i += self.getInt('Sampling')
                kmer = sequence[x:y+1].upper()
                if len(kmer) != k: raise ValueError('kmer buggered!')
                nx = kmer.count('N')
                ## ~ [4b] Pull out kmer count for assembly and raw data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #* pGapK = Probability of Gap (1+ N) based on raw kmer (either 0 or 1)
                #* pErrK = Probability of Error (1X) based on raw kmer
                #* pRedK = Probability of Redundant (0.5NX) based on raw kmer
                #* pHetK = Probability of Heterozygous (1NX) based on raw kmer
                #* pHomK = Probability of Homozygous (2NX) based on raw kmer
                #* pRepK = Probability of Repeat (3NX+) based on raw kmer
                #* RepCNK = Predicted repeat copy number based on raw kmer
                if nx: sentry['pGapK'] = 1.0
                else:
                    ajf = os.popen('jellyfish query %s %s' % (self.getStr('AssemblyJF'),kmer)).read()
                    ajf = string.split(ajf)     #i# e.g. AACACCCCAATGCTCGTAGCCCAAT 1
                    if len(ajf[0]) != k: raise ValueError('jellyfish query problem for assembly %s' % kmer)
                    sentry['AssemblyK'] = int(ajf[1])
                    ajf = os.popen('jellyfish query %s %s' % (self.getStr('RawJF'),kmer)).read()
                    ajf = string.split(ajf)     #i# e.g. AACACCCCAATGCTCGTAGCCCAAT 1
                    if len(ajf[0]) != k: raise ValueError('jellyfish query problem for raw %s' % kmer)
                    sentry['RawK'] = int(ajf[1])
                    sentry['RepCNK'] = sentry['RawK'] / (2.0 * self.getNum('RawDepth'))
                    #i# Raw kmer copy probs
                    kprob = 0.0
                    for p in ['Err','Red','Het','Hom','Rep']:
                        sentry['p%sK' % p] = rje.logPoisson(sentry['RawK'],expected['K'][p],exact=True,callobj=self)
                        kprob += sentry['p%sK' % p]
                    #i# Convert to relative likelihood
                    if kprob:
                        for p in ['Err','Red','Het','Hom','Rep']:
                            sentry['p%sK' % p] /= kprob
                ## ~ [4c] Calculate BLAST depth (local and contig) at position for each %id cutoff ~ ##
                #* NLocal = Number of non-self local BLAST hits (assembly repeat copy number?) (| separated list)
                #* NContig = Number of non-self contigs with 1+ local BLAST hits (assembly copy number) (| separated list)
                # Cycle through sampling dictionaries and take mean
                for locidmin in self.list['LocalIDMin']:
                    sentry['NLocal|%d' % locidmin] = rje.median(locdep[locidmin][x:y+1],avtie=False)
                #i# Using the median value across the kmer for now
                #x#ncopy = rje.mean(ncontig[x:y+1])
                #x#sentry['NContig'] = '%.1f' % ncopy
                sentry['NContig'] = ncopy = rje.median(ncontig[x:y+1],avtie=False)
                #x#sentry['NLocal'] = string.join(sentry['NLocal'],'|')
                #x# Use single value for NContig. sentry['NContig'] = string.join(sentry['NContig'],'|')

                ## ~ [4d] Pull contig out of BAM files and calculate read depth for each position ~~ ##
                #* Dep = Normalised read depth for each BAM file (| separated list)
                #* pErrD = Probability of Error (1X) based on read depth
                #* pRedD = Probability of Redundant (0.5NX) based on read depth
                #* pHetD = Probability of Heterozygous (1NX) based on read depth
                #* pHomD = Probability of Homozygous (2NX) based on read depth
                #* pRepD = Probability of Repeat (3NX+) based on read depth
                # Cycle through sampling dictionaries and take mean
                sentry['Dep'] = []
                for bamfile in self.list['BAMList']:
                    sentry['Dep'].append('%d' % int(rje.mean(bamdep[bamfile][x:y+1])+0.5))
                    # Make a combined probability from all files
                    for p in ['Err','Red','Het','Hom','Rep']:
                        sentry['p%sD' % p] *= rje.logPoisson(int(sentry['Dep'][-1]),expected[bamfile][p],exact=True,callobj=self)
                #i# Convert to relative likelihood
                dprob = 0.0
                for p in ['Err','Red','Het','Hom','Rep']: dprob += sentry['p%sD' % p]
                if dprob:
                    for p in ['Err','Red','Het','Hom','Rep']: sentry['p%sD' % p] /= dprob
                sentry['Dep'] = string.join(sentry['Dep'],'|')

                ## ~ [4e] Add Rating of section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #i# Classify each based on high probability (plus "unknown") - need a probability threshold
                #* Unknown = Regions not mapping to any class with high probability = DEFAULT
                #* Gap = Proportion predicted to be gapped (Ns in kmer)
                #i# kmer has Ns and no contig matches -> Gap
                if sentry['pGapK'] and ncopy == 0:
                    sentry['Rating'] = 'Gap'; sentry['Confidence'] = 'High'
                #* Error = Proportion predicted to be sequencing/assembly error
                if sentry['pErrK'] > 0.5 and ncopy == 0:
                    sentry['Rating'] = 'Error'; sentry['Confidence'] = 'High'
                #* Hetero = Proportion predicted to be haploid coverage (two copies in assembly, one in genome)
                #i# High kmer hetero probability and one contig matches
                elif sentry['pHetK'] > 0.5 and ncopy == 1:
                    sentry['Rating'] = 'Hetero'; sentry['Confidence'] = 'High'
                #* Homo = Proportion predicted to be diploid coverage (one copy in assembly, one in genome)
                #i# High kmer homo probability and no contig matches
                elif sentry['pHomK'] > 0.5 and ncopy == 0:
                    sentry['Rating'] = 'Homo'; sentry['Confidence'] = 'High'
                #* Collapsed = one copy in assembly, 2+ copies in genome
                elif sentry['pRepK'] > 0.5 and ncopy == 1:
                    sentry['Rating'] = 'Collapsed'; sentry['Confidence'] = 'High'
                #* CollapsedRep = 2+ copies in assembly but fewer than one copy in assembly per copy in genome
                elif sentry['pRepK'] > 0.5 and sentry['RepCNK'] > max(sentry['AssemblyK'],ncopy+1):
                    sentry['Rating'] = 'CollapsedRep'; sentry['Confidence'] = 'High'
                #* Repeat = Multiple copies in both genome and assembly (i.e. complex)
                elif sentry['pRepK'] > 0.5 and ncopy > 0:
                    sentry['Rating'] = 'Repeat'; sentry['Confidence'] = 'High'
                #* Redundant = Proportion predicted to be redundant, i.e. more than 2x copies in assembly versus genome.
                elif sentry['RepCNK'] < 0.5 * min(sentry['AssemblyK'],ncopy+1):
                    sentry['Rating'] = 'Redundant'; sentry['Confidence'] = 'High'

                #!# Other ideas/considerations
                ## ~ [- e. Calculate copies in genome based on GABLAM local hits at each %identity
                ## ~ [- f. Calculate genome:assembly ratio from each BAM file using raw k and assembly k
                ## ~ [- g. Calculate genome copies using raw k and read depth (average over different BAM files)

                ## ~ [4x] Add entry to table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sdb.addEntry(sentry)
            self.printLog('\r#KMERS','Processed %s %s kmers: %s sampling entries (%s bp)' % (rje.iStr(px),contig,rje.iStr(sdb.entryNum()),rje.iLen(sequence)))

            ### ~ [5] Step through entries and calculate probabilities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.saveToFile(fdict['ResFile']['sampling'])


            ### ~ [5] Combine statistics into per-contig summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #* Contig = contig name.
            #* Len = contig length
            #* SCDep = Median BAM read depths for SC regions (AssemblyK=1, No local BLAST hits)
            #* DCDep = Median BAM read depths for DC regions (AssemblyK=1 or 2, 1 local BLAST hits)
            #* IDSC = Proportion single copy based on GABLAM (| sep list for each %id cutoff)
            #* IDDC = Proportion double copy based on GABLAM (| sep list for each %id cutoff, n)
            #* IDMC = Proportion multi copy based on GABLAM (| sep list for each %id cutoff, n)
            ## ~ [5a] Convert sampling data into contig stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Split the contig BLAST depth counts and convert into copy numbers per ID cutoff.
            #i# Rejoin and tidy the ID* fields after compression and taking means.
            for field in ['Gap','Error','Redundant','Hetero','Homo','Collapsed','CollapsedRep','Repeat','Unknown']:
                sdb.addField(field,evalue=0)
            for i in range(len(self.list['LocalIDMin'])):
                locidmin = self.list['LocalIDMin'][i]
                for field in ['IDSC','IDDC','IDMC']:
                    sdb.addField('%s|%s' % (field,locidmin),evalue=0)
                for sentry in sdb.entries():
                    ncopy = sentry['NLocal|%d' % locidmin]
                    if ncopy < 1: field = 'IDSC'
                    elif ncopy == 1: field = 'IDDC'
                    else: field = 'IDMC'
                    sentry['%s|%s' % (field,locidmin)] = 1.0
                    sentry[sentry['Rating']] = 1.0
            #i# Keep most robust analysis as the BLAST-based assembly copy number
            #i# Compare to the Assembly Kmer copy number? Or Raw:Assembly Kmer ratio

            depfields = ['Pos','Copy']
            for i in range(len(self.list['BAMList'])): depfields.append('BAM%d' % i)
            depdb = self.db().addEmptyTable('%s.dep' % contig,depfields,['Pos'],log=True)
            locidmin = self.list['LocalIDMin'][0]
            for sentry in sdb.entries():
                ncopy = sentry['NLocal|%d' % locidmin]
                dentry = {'Pos':sentry['Pos']}
                dbam = string.split(sentry['Dep'],'|')
                for i in range(len(self.list['BAMList'])):
                    dentry['BAM%d' % i] = int(dbam[i])
                if sentry['AssemblyK'] == 1 and ncopy == 0:         # Solid single copy
                    dentry['Copy'] = 'SC'
                    depdb.addEntry(dentry)
                elif 0 < sentry['AssemblyK'] < 3 and ncopy < 2:     # Probable double copy
                    dentry['Copy'] = 'DC'
                    depdb.addEntry(dentry)
            depdb.compress(['Copy'],default='median')
            #bamdep = {}
            #for dentry in depdb.entries():
            #    dep = []
            #    for i in range(len(self.list['BAMList'])):
            #        dep.append('%.1f' % dentry['BAM%d' % i])
            #    bamdep[dentry['Copy']] = string.join(dep,'|')

            #sdb.dropFields(['NLocal','Dep'])

            ## ~ [5b] Compress the sampling entries, taking mean values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb.compress(['Contig'])
            sentry = sdb.entries()[0] #i# should only be one!
            #for field in ['IDSC','IDDC','IDMC']:
                #sentry[field] = []
                #for locidmin in self.list['LocalIDMin']:
                #    sentry[field].append('%.3f' % sentry['%s|%s' % (field,locidmin)])
                #sentry[field] = string.join(sentry[field],'|')
            sentry['Len'] = len(sequence)
            for dentry in depdb.entries():
                for i in range(len(self.list['BAMList'])):
                    if i: sentry['%sDep' % dentry['Copy']] = dentry['BAM%d' % i]
                    else: sentry['%sDep%d' % (dentry['Copy'],i+1)] = dentry['BAM%d' % i]
                bamdep[dentry['Copy']] = string.join(dep,'|')

            #sentry['SCDep'] = bamdep['SC']
            #sentry['DCDep'] = bamdep['DC']

            ## ~ [- a. Classify contigs. #!# Need classification system #!#
            ### ~ [6. Optional masked/reformatted output of contigs.

            cdb.addEntry(sentry)
            cdb.saveToFile(fdict['ResFile']['contigark'])
        except:  self.errorLog('%s.contigFork(%s) Error' % (self.prog(),contig))
#########################################################################################################################
    def forking(self):  ### Keeps forking out and processing jobs until no more jobs in self.list['Forked'].
        '''
        Keeps forking out and processing jobs until no more jobs in self.list['Forked'].
        '''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('PIDCheck') or self.dev(): pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])    # Set *.pid object to match log
        else: pidcheck = None
        #self.deBug(pidcheck)
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.list['Forked']:
            if pidcheck: PIDCHECK = open(pidcheck,'w')
            for fdict in self.list['Forked'][0:]:
                try:
                    pid = fdict['PID']
                    if pidcheck: PIDCHECK.write('%s: %s\n' % (self.list['Forked'].index(fdict),pid))
                    if string.split('%s' % pid)[0] == 'WAIT': status = 1
                    else: (status,exit_stat) = os.waitpid(pid,os.WNOHANG)
                except:
                    self.errorLog('!')
                    status = 1
                if status > 0:
                    self.list['Forked'].remove(fdict)
                    self.endFork(fdict)   # Fork has finished: can replace with processing
            if pidcheck:
                PIDCHECK.close()
                #self.deBug(open(pidcheck,'r').read())
            ## ~ [2a] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if time.time() - self.getNum('KillTime') > self.getNum('KillForks'):
                self.verbose(0,1,'\n%d seconds of main thread inactivity. %d forks still active!' % (self.getNum('KillForks'),len(self.list['Forked'])),1)
                for fdict in self.list['Forked']:
                    self.verbose(0,2,' => Fork %s, PID %d still Active!' % (fdict['ID'],fdict['PID']),1)
                if self.i() < 0 or rje.yesNo('Kill Main Thread?'):
                    raise ValueError('%d seconds of main thread inactivity. %d forks still active!' % (self.getNum('KillForks'),len(self.list['Forked'])))
                elif rje.yesNo('Kill hanging forks?'):
                    for fdict in self.list['Forked']:
                        self.printLog('#KILL','Killing Fork %s, PID %d.' % (fdict['ID'],fdict['PID']))
                        os.system('kill %d' % fdict['PID'])
                else: self.setNum({'KillTime':time.time()})
            ## ~ [2b] Sleep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            time.sleep(self.getNum('ForkSleep'))
        if pidcheck and not self.list['Forked'] and not self.list['ToFork'] and rje.exists(pidcheck): os.unlink(pidcheck)
#########################################################################################################################
    def nextFork(self):  ### Sets an new contig fork running.
        '''Sets an new contig fork running.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.list['ToFork']: self.setNum({'KillTime':time.time()}); return      # No more runs to fork
            fdict = {}      # Setup empty dictionary to fill, if jobs available
            self.list['Forked'].append(fdict)
            ## ~ [0a] ~ Check memory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            freemem = self.freeMem()
            fdict['Mem'] = freemem * 100.0
            if self.getNum('MemFree') > 0.0 and freemem < self.getNum('MemFree'):
                fdict['PID'] = 'WAIT - %.1f%% free memory.' % (fdict['Mem'])
                return
            self.startFork(fdict)
        except: self.errorLog('Forker.nextFork error')
#########################################################################################################################
    def startFork(self,fdict):  ### Sets a new fork going using the data in fdict.
        '''Sets a new fork going using the data in fdict.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (fdict['contig'],fdict['type']) = self.list['ToFork'].pop(0)
            fdict['arkbase'] = self.baseFile()
            fdict['sequence'] = self.obj['SeqList'].getSeq(self.obj['SeqList'].seqNameDic()[fdict['contig']])[1]
            fdict['maskseq'] = self.obj['MaskSeq'].getSeq(self.obj['MaskSeq'].seqNameDic()[fdict['contig']])[1]
            fdict['ID'] = 'Fork %d' % self.list['Forked'].index(fdict)
            # Fork ID, used as base for outputs
            fdict['FID'] = '%s.fork' % fdict['contig']
            fdict['Log'] = '%s%s.fork.%s.log' % (self.getStr('RunPath'),fdict['contig'],fdict['type'])
            fdict['ResFile'] = {}
            if 'type' in fdict:
                restype = fdict['type']
                if restype == 'copy':
                    if self.getBool('ContigFas'): fdict['ResFile']['contig.fas'] = '%s.contig.fas' % (fdict['FID'])
                    if self.getBool('LocalSAM'): fdict['ResFile']['local.sam'] = '%s.local.sam' % (fdict['FID'])
                if restype == 'combine':
                    for otype in ['hitsum','gablam','local','sampling']:
                        fdict['ResFile'][otype] = '%s.%s.tdt' % (fdict['FID'],otype)
                if restype == 'calculate':
                    for otype in ['posark','contigark']:
                        fdict['ResFile'][otype] = '%s.%s.tdt' % (fdict['FID'],otype)
                #fdict['ResFile'][restype] = '%s.%s.tdt' % (fdict['FID'],restype)
                #?# Don't want these in ResFile until the combine phase #?#
                #?# One future option is to combine right at the end after compiling all different table types?
            else:   #i# No type = old run method
                for restype in ['hitsum','gablam','local','contigark','sampling']:
                    fdict['ResFile'][restype] = '%s.%s.tdt' % (fdict['FID'],restype)
                if self.getBool('ContigFas'): fdict['ResFile']['contig.fas'] = '%s.contig.fas' % (fdict['FID'])
                if self.getBool('LocalSAM'): fdict['ResFile']['local.sam'] = '%s.local.sam' % (fdict['FID'])
                fdict['Cleanup'] = ['contig.fas','contig.fas.index']
            try: open(fdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting fork.'); return self.endJob(fdict)

            ### ~ [1] Special Debugging shortcut to single contig run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('NoForks') or self.getNum('Forks') < 2 or (self.dev() and self.debugging()):
                fdict['PID'] = '%s.%s' % (fdict['contig'],fdict['type'])
                #self.baseFile(fdict['FID'])
                self.runFork(fdict)
                self.list['Forked'].remove(fdict)
                self.endFork(fdict)   # Process finished "Fork"
                return False

            ### ~ [2] ~ Add Fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setNum({'KillTime':time.time()})
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                fdict['PID'] = cpid
                self.printLog('#FORK','Forking %s.%s as %s: %d remain; %.1f%% mem free' % (fdict['contig'],fdict['type'],cpid,len(self.list['ToFork']),fdict['Mem']))
                self.printLog('#FORK','%s: %s.%s' % (cpid,fdict['FID'],fdict['type']))
            else:                   # child process
                self.setBool({'Child':True})   # Whether this is a child process
                self.setInt({'Interactive':-1})
                self.log.stat['Interactive'] = -1
                self.log.info['LogFile'] = fdict['Log']
                self.baseFile(fdict['FID'])
                if self.needFork(fdict): self.runFork(fdict)
                else: self.printLog('#DONE','Skipping %s %s: complete' % (fdict['contig'],fdict['type']))
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('Forker.startFork error')
#########################################################################################################################
    def runFork(self,fdict):
        if fdict['type'] == 'copy': self.contigForkGABLAM(fdict)
        elif fdict['type'][:3] == 'dep': self.contigForkDepth(fdict)
        elif fdict['type'] in ['rawk','assk']: self.contigForkKmer(fdict)
        elif fdict['type'] == 'combine': self.contigForkCombine(fdict)
        elif fdict['type'] == 'calculate': self.contigForkCalculate(fdict)
        else: self.contigFork(fdict)
#########################################################################################################################
    def endFork(self,fdict):   ### Ends fork, tidies and sets new one running
        '''Ends fork, tidies and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'ResFile' in fdict:
                for resfile in fdict['ResFile']:
                    fromfile = fdict['ResFile'][resfile]
                    if not rje.exists(fromfile): self.warnLog('Results file %s missing!' % fromfile); continue
                    if fromfile.endswith('.tdt'):
                        tofile = '%s.%s.tdt' % (self.baseFile(),resfile)
                        if rje.exists(tofile):
                            open(tofile,'a').writelines(open(fromfile,'r').readlines()[1:])
                            if not self.dev(): os.unlink(fromfile)
                        else: rje.fileTransfer(fromfile,tofile)
                    else:
                        tofile =  '%s.%s' % (self.baseFile(),resfile)
                        rje.fileTransfer(fromfile,tofile)
            if 'Cleanup' in fdict:
                for clean in fdict['Cleanup']:
                    cleanfile = '%s.%s' % (fdict['FID'],clean)
                    if rje.exists(cleanfile): os.unlink(cleanfile)
            if 'Log' in fdict:
                rje.fileTransfer(fdict['Log'],self.log.info['LogFile'])
                if self.getBool('LogFork'):
                    self.printLog('#END','Fork %s ended: %s.%s log content and results transferred' % (fdict['PID'],fdict['contig'],fdict['type']))
                    self.printLog('#~~#','#~~#',timeout=False)
                #if self.dev(): self.deBug(fdict['Log'])
                #if self.dev(): self.deBug(rje.exists(fdict['Log']))
            elif 'PID' in fdict and string.split('%s' % fdict['PID'])[0] == 'WAIT': pass
            else: self.printLog('#END','Fork %s ended.' % fdict['PID'])
        except IOError:
            if self.getInt('IOError') == 1: self.errorLog('Forker.endFork IOError limit reached'); raise
            else: self.int['IOError'] -= 1; self.errorLog('Forker.endFork')
        except: self.errorLog('Forker.endFork error')
        self.nextFork()   # Carry on regardless
#########################################################################################################################
    ### <4> ### Special Class Methods                                                                                #
#########################################################################################################################
    def depthCalc(self):    ### Work through contigs in order and update depth calculations from first bamfile.
        '''Work through contigs in order and update depth calculations from first bamfile.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bamfile = self.list['BAMList'][0]
            self.printLog('#CALC','DepthCalc mode for: %s' % bamfile)
            contigs = self.list['Contigs']
            skiplist = self.list['SkipList']
            if skiplist: skiplist = rje.listIntersect(contigs,self.list['SkipList'])

            ## ~ [0a] Setup output tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db = self.db()
            db.list['Tables'] = []
            dfile = '%s.dep.tdt' % (self.baseFile())
            hfile = '%s.hist.tdt' % (self.baseFile())
            if rje.exists(dfile) and rje.exists(hfile) and self.getBool('Pickup') and not self.force():
                ddb = db.addTable(dfile,mainkeys=['Contig'],name='dep')
                ddb.dataFormat({'Reads':'int'})
                readx = sum(ddb.dataList(ddb.entries(),'Reads',sortunique=False))
                hdb = db.addTable(hfile,name='hist',mainkeys=['Dep'])
                hdb.dataFormat({'Dep':'int','Count':'int','Mode':'int'})
            else:
                ddb = db.addEmptyTable('dep',['Contig','Len','Null','Min','Mean','Median','Mode','Max','Sum','Raw','Reads'],['Contig'])
                ddb.saveToFile()
                hdb = db.addEmptyTable('hist',['Dep','Count','Mode'],['Dep'])
                hdb.saveToFile()
                readx = 0

            ### ~ [1] Pull out BAM depth lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cx = 0
            for contig in rje.listDifference(contigs,skiplist):
                cx += 1
                self.printLog('#CONTIG','Contig %d: %s' % (cx,contig))
                sequence = self.obj['SeqList'].getSeq(self.obj['SeqList'].seqNameDic()[contig])[1]
                maskseq = self.obj['MaskSeq'].getSeq(self.obj['MaskSeq'].seqNameDic()[contig])[1]

                ## ~ [1a] Read summary for contig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sam = rje_samtools.SAMtools(self.log,self.cmd_list)
                pcmd = 'samtools view %s %s' % (bamfile,contig)
                SAM = os.popen(pcmd)
                rid = 0          # Read counter (ID counter)
                ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                dentry = {'Contig':contig,'Len':len(sequence),'Null':0,'Raw':0,'Reads':0}
                for line in SAM:

                    ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if line.startswith('@'): continue
                    samdata = string.split(line)
                    if len(samdata) < 11: continue
                    self.progLog('\r#PARSE','Parsing %s: %s reads...' % (bamfile,rje.iStr(rid)),rand=0.1)
                    rid += 1
                    rname = samdata[0]
                    zdata = string.split(string.replace(rname,'/',' '))
                    try:
                        [smrt,zmw,pos] = zdata[:3]
                        [beg,end] = string.split(pos,'_')
                        rlen = int(end) - int(beg) + 1
                    except: rlen = len(samdata[9])
                    cigstr = samdata[5]
                    cigdata = rje_samtools.parseCigar(cigstr)
                    if 'S' in cigdata: mlen = rlen - cigdata['S']
                    else: mlen = rlen
                    dentry['Raw'] += mlen
                    dentry['Reads'] += 1
                readx += dentry['Reads']
                SAM.close()
                self.printLog('\r#PARSE','Parsed %s: %s reads; %s bp' % (bamfile,rje.iStr(rid),rje.iStr(dentry['Raw'])))

                ## ~ [1b] Pileup depths for contig ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                counts = [0]
                qc = [0]        # List of counts at different quality scores
                ri = 0          # RID Counter
                ridlist = []    # List of read IDs
                rdel = {}       # Dictionary of {rid:current deletion}
                pcmd = 'samtools view %s %s -b | samtools mpileup -BQ0 -d10000000 -f %s -' % (bamfile,contig,self.getStr('SeqIn'))
                self.printLog('#BAM',pcmd)
                PILEUP = os.popen(pcmd); bx = 0
                while PILEUP:
                    #if self.dev():
                    self.progLog('\r#BAMDEP','%s: %s %s bases' % (bamfile,rje.iStr(bx),contig),rand=0.001)
                    line = PILEUP.readline()
                    if not line: break
                    bentry = sam.parsePileupLine(line,ri,ridlist,rdel,qc); bx += 1
                    dep = bentry['N']
                    while dep >= len(counts): counts.append(0)
                    counts[dep] += 1

                    if ridlist: ri = max(ri,max(ridlist))
                PILEUP.close()
                dentry['Null'] = dentry['Len'] - bx
                self.printLog('\r#BAMDEP','%s: %s -> %sX depth' % (bamfile,contig,rje.iLen(counts[1:])))

                ## ~ [1c] Update statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                dentry['Sum'] = 0
                for dx in range(len(counts)):
                    dentry['Sum'] += counts[dx] * dx
                if sum(counts):
                    dentry['Mean'] = dentry['Sum'] / float(sum(counts))
                else: dentry['Mean'] = 0.0
                half = int((sum(counts)+0.5)/2.0)
                dentry['Null'] += counts[0]
                counts[0] = dentry['Null']
                dentry['Min'] = -1
                dentry['Mode'] = 0
                dentry['Median'] = 0
                sumx = 0
                for dx in range(len(counts)):
                    sumx += counts[dx]
                    if sumx >= half and not dentry['Median']: dentry['Median'] = dx
                    if counts[dx] >= counts[dentry['Mode']]: dentry['Mode'] = dx
                    if counts[dx] > 0:
                        dentry['Max'] = dx
                        if dentry['Min'] < 0: dentry['Min'] = dx
                ddb.addEntry(dentry)
                self.printLog('#DEPTH','%s: %s' % (contig,dentry))
                print(ddb.entrySummary(dentry))

                modex = 0
                for dx in range(len(counts)):
                    self.debug(dx in hdb.data())
                    if dx in hdb.data():
                        self.debug(hdb.data()[dx])
                        hdb.data()[dx]['Count'] += counts[dx]
                    elif counts[dx]: #i# Only output entries with counts
                        hdb.addEntry({'Dep':dx,'Count':counts[dx],'Mode':0})
                    if counts[dx]:
                        if not modex: modex = dx
                        elif dx in hdb.data() and hdb.data()[dx]['Count'] >= hdb.data()[modex]['Count']: modex = dx
                hdb.data()[dentry['Mode']]['Mode'] += 1
                self.printLog('#MODE','%s contigs (%s reads) => total mode = %dX' % (rje.iStr(ddb.entryNum()),rje.iStr(readx),modex))

                ddb.saveToFile(append=True,savekeys=[contig])
                hdb.saveToFile(append=False,backup=False)

        except:  self.errorLog('%s.depthCalc(%s) Error' % (self.prog(),bamfile))
#########################################################################################################################
### End of SECTION II: GenomeARK Class                                                                                  #
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
    try: GenomeARK(mainlog,cmd_list).run()

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
