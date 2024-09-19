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
Module:       depthsizer
Description:  Read-depth based genome size prediction
Version:      1.9.3
Last Edit:    01/08/24
Citation:     Chen SH et al. & Edwards RJ (2022): Mol. Ecol. Res. (doi: 10.1111/1755-0998.13574)
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:

    DepthSizer is a modified wrapper for the genome size estimate methods of Diploidocus. DepthSizer needs a genome assembly
    (fasta format, `seqin=FILE`), a set of long read (ONT, PacBio or HiFi) data for the assembly (`reads=FILELIST` and
    `readtype=LIST`) (or `readbp=INT`), and a BUSCO full table of results (`busco=TSVFILE`).

    DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid
    read depth) regions along with some poor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
    are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
    region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
    read depth.

    DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth.
    This is converted into an estimated single copy read depth using a density plot of BUSCO single copy genes.
    Legacy mode will calculate modal read depth for each BUSCO gene along with the overall modal read depth for all gene
    regions. Genome size is then estimated based on a crude calculation using the total combined sequencing length.
    This will be caculated from `reads=FILELIST` unless provided with `readbp=INT`.

    BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
    `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
    ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
    and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
    Output from BUSCOMP is also compatible with DepthSizer.

    **NOTE:** The current unadjusted genome size prediction appears to be an over-estimate. Please see documentation for
    more details of attempts to correct for contamination, read mapping and/or imbalanced insertion:deletion ratios
    etc.

    ---

Dependencies:

    Unless `bam=FILE` is given, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
    environment `$PATH` or given to DepthSizer with the `minimap2=PROG` setting, and [samtools](http://www.htslib.org/)
    needs to be installed. Unless `legacy=T depdensity=F`, R will also need be installed.

    To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For DepthSizer documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Main DepthSizer run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [gapspanner or $SEQIN basefile]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    bamcsi=T/F      : Use CSI indexing for BAM files, not BAI (needed for v long scaffolds) [False]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    dochtml=T/F     : Generate HTML DepthSizer documentation (*.docs.html) instead of main run [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
    ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
    adjustmode=X    : Map adjustment method to apply (None/CovBases/IndelRatio/MapBases/MapAdjust/MapRatio/OldAdjust/OldCovBases) [IndelRatio]
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    depchunk=INT    : Chunk input into minimum of INT bp chunks for temp depth calculation [1e6]
    deponly=T/F     : Cease execution following checking/creating BAM and fastdep/fastmp files [False]
    depfile=FILE    : Precomputed depth file (*.fastdep or *.fastmp) to use [None]
    covbases=T/F    : Whether to calculate predicted minimum genome size based on mapped reads only [True]
    mapadjust=T/F   : Whether to calculate mapadjust predicted genome size based on read length:mapping ratio [False]
    benchmark=T/F   : Activate benchmarking mode and also output the assembly size and mean depth estimate [False]
    legacy=T/F      : Whether to perform Legacy v1.0.0 (Diploidocus) calculations [False]
    depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth in legacy mode [True]
    depadjust=INT   : Advanced R density bandwidth adjustment parameter [12]
    seqstats=T/F    : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
    reduced=T/F     : Only generate/use fastmp for BUSCO-containing sequences (*.busco.fastmp) [True]
    fragmented=T/F  : Whether to use Fragmented as well as Complete BUSCO genes for SC Depth estimates [False]
    ### ~ Forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
    logfork=T/F     : Whether to log forking in main log [False]
    memsaver=T/F    : Whether to disable threading for the R script [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_rmd, rje_readcore, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 1.0.0 - Initial working version of DepthSizer based on Diploidocus v0.16.2.
    # 1.1.0 - Updated to use updated ReadCore and DepthCopy code, plus new Indel ratio calculation (dev=T).
    # 1.2.0 - Added CovBases lower depthsizer estimate output based solely on mapped reads. Fixed unmapped read bug.
    # 1.2.1 - Fixed major flaw in indelratio calculation.
    # 1.3.0 - Added adjust string and benchmark=T option for additional calculations. IndelRatio no longer dev only.
    # 1.3.1 - Tweaked some input checks and log output. Replaced indelratio sort -u with uniq for speed and memory.
    # 1.4.0 - Added seqstats=T/F : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
    # 1.4.1 - Added citation and fixed minor output typo.
    # 1.4.2 - Fixed bug that causes clashes with v5 full_table.bed files.
    # 1.5.0 - Add additional map adjustment variants:
    #       - MapAdjust2 = allbases, not covbases
    #       - MapBases = Use map bases, not covbases for min read volumne
    #       - MapRatio = Use mapbases adjusted by indelratio
    # 1.6.0 - Disable legacy mode using Diploidocus.
    # 1.6.1 - Bug fixes to underlying R script and related core codebase.
    # 1.6.2 - Updated citation to Mol Ecol Res paper.
    # 1.6.3 - Fixed R code bug. Added bamcsi=T/F to use CSI indexing.
    # 1.7.0 - Fixed a problem with lack of Duplicated BUSCOs. Added fragmented=T option.
    # 1.7.1 - Fixed inconsistency with output s.f.
    # 1.8.0 - Added reduced=T/F : Only generate/use fastmp for BUSCO-containing sequences (*.busco.fastmp) [True]
    # 1.9.0 - Added multi-threading to the R script and chunking of input sequences for depth calculation.
    # 1.9.1 - Fixed indelratio to cope with -ve strand BUSCO formatting.
    # 1.9.2 - Fixed fragmented mode to use different scdepth file.
    # 1.9.3 - Fixed bug with chunking of input sequences when reduced=T that missed some BUSCO genes.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [N] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [Y] : Update docstring.
    # [Y] : Work out why the indel ratio is not working properly.
    # [Y] : Rationalise the log outputs etc.
    # [Y] : Add Lower output based solely on mapped reads (+/- Adjust)
    # [ ] : Add the sequence, window and chromcheck settings to DepthSizer. (Also available in DepthKopy.)
    # [Y] : Check why genomesize is an option and possibly remove? (Legacy from Diploidocus origins?)
    # [ ] : Add option to use complete genome rather than BUSCOs for the SC depth.
    # [Y] : Add a reduced=T/F option that only generates fastmp for BUSCO-containing sequences (*.busco.fastmp)
    # [ ] : Add a forkdep=T/F option that only generates the BAM & fastmps file then stops, for efficient Nextflow etc.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('DepthSizer', '1.9.3', 'August 2024', '2021')
    description = 'Read-depth based genome size prediction'
    author = 'Dr Richard J. Edwards.'
    comments = ['Citation: Chen SH et al. & Edwards RJ (2022): Mol. Ecol. Res. (doi: 10.1111/1755-0998.13574)',
                'Please raise bugs or questions at https://github.com/slimsuite/depthsizer.',rje_obj.zen()]
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
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: rje.printf('Major Problem with cmdHelp()')
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
        if len(sys.argv) == 2 and sys.argv[1] in ['details','-details','--details']: rje.printf('%s v%s' % (info.program,info.version)); sys.exit(0)
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
    except: rje.printf('Problem during initial setup.'); raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: DepthSizer Class                                                                                       #
#########################################################################################################################
class DepthSizer(rje_readcore.ReadCore):
    '''
    DepthSizer Class. Author: Rich Edwards (2021).

    Str:str
    - AdjustMode=X        : Map adjustment method to apply (None/CovBases/IndelRatio/MapAdjust) [IndelRatio]

    Bool:boolean
    - Benchmark=T/F   : Activate benchmarking mode and also output the assembly size and mean depth estimate [False]
    - CovBases=T/F    : Whether to calculate predicted minimum genome size based on mapped reads only [True]
    - DepOnly=T/F     : Cease execution following checking/creating BAM and fastdep/fastmp files [False]
    - Fragmented=T/F  : Whether to use Fragmented as well as Complete BUSCO genes for SC Depth estimates [False]
    - Legacy=T/F      : Whether to perform Legacy v1.0.0 (Diploidocus) calculations [False]
    - Reduced=T/F     : Only generate/use fastmp for BUSCO-containing sequences (*.busco.fastmp) [True]

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['AdjustMode']
        self.boollist = ['Benchmark','CovBases','DepOnly','DocHTML','Fragmented','Legacy','MapAdjust','Reduced']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self._setReadCoreAttributes()   # See rje_readcore
        self.setStr({'AdjustMode':'IndelRatio'})
        self.setBool({'CovBases':True,'Legacy':False,'Reduced':True})
        self.setInt({})
        self.setNum({})
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
                self._readCoreCmd(cmd)  # Will set all the core commands recognised.
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['AdjustMode'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['Benchmark','CovBases','DepOnly','DocHTML','Fragmented','Legacy','MapAdjust','Reduced'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # DepthSizer: Read-depth based genome size prediction

        DepthSizer is an updated version of the genome size estimate methods of Diploidocus. DepthSizer needs a genome assembly
        (fasta format, `seqin=FILE`), a set of long read (ONT, PacBio or HiFi) data for the assembly (`reads=FILELIST` and
        `readtype=LIST`) (or `readbp=INT`), and a BUSCO/BUSCOMP full table of results (`busco=TSVFILE`).

        DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid
        read depth) regions along with some poor quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
        are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
        region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
        read depth.

        DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth.
        This is converted into an estimated single copy read depth using a smoothed density plot of BUSCO single copy genes.
        Genome size is then estimated based on a crude calculation using the total combined sequencing length.
        This will be calculated from `reads=FILELIST` unless provided with `readbp=INT`.

        BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table matching the BUSCO fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
        and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
        Output from BUSCOMP is also compatible with DepthSizer. DepthSizer has been tested with outputs from BUSCO v3 and v5.

        **NOTE:** The basic DepthSizer approach assumes that the raw long read data has a 1:1 correspondence to the
        genomic DNA being sequenced, i.e. there is no contamination (including plastids) and no bias towards insertion
        or deletion read errors. As a consequence, the default genome size prediction is expected to be an over-estimate.
        DepthSizer will also calculate an estimated lower bound, based on only those reads that map to the assembly (unless
        `covbases=F`) . An adjustment for read error profiles is made by calculating the ratio of read:genomic data for
        mapped read from the BAM CIGAR strings ((insertions+matches)/(deletions+matches)) and reported as the `IndelRatio`
        adjustment. The older `MapAjust` method, which uses mapped reads and mapped bases calculated from `samtools coverage`
        and `samtools fasta`) to try to correct for read mapping and imbalanced insertion:deletion ratios, can also be
        switched on with `mapadjust=T` (or `benchmark=T`). Benchmarking of the different adjustments is ongoing. Read
        volumes can also be manually adjusted with `readbp=INT`. All calculated sizes will be reported in the
        `*.gensize.tdt` output, but the adjustment method selected by `adjustmode=X` (None/CovBases/IndelRatio/MapAdjust,
        default `IndelRatio`) will be used for "the" genome size prediction.

        **Version 1.1.** The core depth calculation shifted in Version 1.1. `Legacy` mode will use the old code to
        calculate the modal read depth for each BUSCO gene along with the overall modal read depth for all gene
        regions. These are not recommended.

        **Version 1.8.** Version 1.8 introduced a new `reduced=T/F` mode, which only processes sequences that have BUSCO
        predictions. (Complete, Duplicated or Fragmented.) This is _on_ (`True`) by default, and substantially reduces
        the disk footprint and processing time for highly fragmented genomes. If the BUSCO Completeness is low, using the
        `fragmented=T` option (introduced in version 1.7, default `False`) will use `Fragmented` BUSCO genes as well as
        `Complete` genes to establish the single-copy read depth.

        ## Citation

        DepthSizer has been published as part of the Waratah genome paper:

        > Chen SH, Rossetto M, van der Merwe M, Lu-Irving P, Yap JS, Sauquet H, Bourke G, Amos TG, Bragg JG & Edwards RJ (2022).
        Chromosome-level de novo genome assembly of Telopea speciosissima (New South Wales waratah) using long-reads,
        linked-reads and Hi-C. Molecular Ecology Resources doi: [10.1111/1755-0998.13574](https://doi.org/10.1111/1755-0998.13574)

        Please contact the author if you have trouble getting the full text version, or read the bioRxiv preprint version:

        > Chromosome-level de novo genome assembly of Telopea speciosissima (New South Wales waratah) using long-reads,
        linked-reads and Hi-C. [bioRxiv 2021.06.02.444084](https://www.biorxiv.org/content/10.1101/2021.06.02.444084v2.full);
        doi: 10.1101/2021.06.02.444084.

        ---

        # Running DepthSizer

        DepthSizer is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/depthsizer.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [DepthSizer git repo](https://github.com/slimsuite/depthsizer), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [DepthSizer git repo](https://github.com/slimsuite/depthsizer)
        for running on example data.

        ## Dependencies

        Unless `bam=FILE` is given, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the
        environment `$PATH` or given to DepthSizer with the `minimap2=PROG` setting, and [samtools](http://www.htslib.org/)
        needs to be installed. Unless `legacy=T depdensity=F`, R will also need be installed.

        To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For DepthSizer documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main DepthSizer run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        basefile=FILE   : Root of output file names [gapspanner or $SEQIN basefile]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
        bamcsi=T/F      : Use CSI indexing for BAM files, not BAI (needed for v long scaffolds) [False]
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        dochtml=T/F     : Generate HTML DepthSizer documentation (*.docs.html) instead of main run [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~ Genome size prediction options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
        readbp=INT      : Total combined read length for depth calculations (over-rides reads=FILELIST) []
        adjustmode=X    : Map adjustment method to apply (None/CovBases/IndelRatio/MapBases/MapAdjust/MapRatio/OldAdjust/OldCovBases) [IndelRatio]
        quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
        depchunk=INT    : Chunk input into minimum of INT bp chunks for temp depth calculation [1e6]
        deponly=T/F     : Cease execution following checking/creating BAM and fastdep/fastmp files [False]
        depfile=FILE    : Precomputed depth file (*.fastdep or *.fastmp) to use [None]
        covbases=T/F    : Whether to calculate predicted minimum genome size based on mapped reads only [True]
        mapadjust=T/F   : Whether to calculate mapadjust predicted genome size based on read length:mapping ratio [False]
        benchmark=T/F   : Activate benchmarking mode and also output the assembly size and mean depth estimate [False]
        legacy=T/F      : Whether to perform Legacy v1.0.0 (Diploidocus) calculations [False]
        depdensity=T/F  : Whether to use the BUSCO depth density profile in place of modal depth in legacy mode [True]
        depadjust=INT   : Advanced R density bandwidth adjustment parameter [12]
        seqstats=T/F    : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
        reduced=T/F     : Only generate/use fastmp for BUSCO-containing sequences (*.busco.fastmp) [True]
        fragmented=T/F  : Whether to use Fragmented as well as Complete BUSCO genes for SC Depth estimates [False]
        ### ~ Forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
        logfork=T/F     : Whether to log forking in main log [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        # DepthSizer workflow and options

        The main inputs for DepthSizer genome size prediction are:

        * `seqin=FILE` : Input sequence assembly to tidy [Required].
        * `reads=FILELIST`: List of fasta/fastq files containing long reads. Wildcard allowed. Can be gzipped.
        * `readtype=LIST` : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        * `busco=TSVFILE` : BUSCO full table [`full_table_$BASEFILE.busco.tsv`] used for calculating single copy ("diploid") read depth.

        ## Step 1: BAM file (read mapping)

        The first step is to generate a BAM file by mapping `reads` on to `seqin` using [minimap2](https://github.com/lh3/minimap2).
        A pre-generated BAM file can be given instead using `bam=FILE`. There should be no secondary mapping of reads, as
        these will inflate read depths, so filter these out if they were allowed during mapping. Similarly, the BAM file
        should not contain unmapped reads. (These should be filtered during processing if present.) If no BAM file setting
        is given, the BAM file will be named `$BASEFILE.bam`, where `$BASEFILE` is set by `basefile=X`.

        ## Step 2: BUSCO(MP) results

        DepthSizer works on the principle that `Complete` BUSCO genes should represent predominantly single copy (diploid read
        depth) regions along with some poor-quality and/or repeat regions. Assembly artefacts and collapsed repeats etc.
        are predicted to deviate from diploid read depth in an inconsistent manner. Therefore, even if less than half the
        region is actually diploid coverage, the **modal** read depth is expected to represent the actual single copy
        read depth. This is estimated using a smoothed density distribution calculated using R `density()`.

        BUSCO single-copy genes are parsed from a BUSCO full results table, given by `busco=TSVFILE` (default
        `full_table_$BASEFILE.busco.tsv`). This can be replaced with any table with the fields:
        ['BuscoID','Status','Contig','Start','End','Score','Length']. Entries are reduced to those with `Status` = `Complete`
        and the `Contig`, `Start` and `End` fields are used to define the regions that should be predominantly single copy.
        [BUSCOMP](https://github.com/slimsuite/buscomp) v0.10.0 and above will generate a `*.complete.tsv` file that can
        be used in place of BUSCO results. This can enable rapid re-annotation of BUSCO genes following, for example,
        vector trimming with [Diploidocus](https://github.com/slimsuite/diploidocus). If `fragmented=T` then entries with
        `Status` = `Fragmented` are also used. This is useful when Completeness is low.

        ## Step 3: Single-copy read depth

        DepthSizer uses `samtools mpileup` (or `samtools depth` if `quickdepth=T`) to calculate the per-base read depth
        and extracts the smoothed modal read depth for all single-copy (`Complete` BUSCO genes) using the `density()`
        function of R. To avoid a minority of extremely deep-coverage bases disrupting the density profile, the depth
        range is first limited to the range from zero to 1000, or four time the pure modal read depth if over 1000. If
        the pure mode is zero coverage, zero is returned. The number of bins for the density function is set to be
        greater than 5 times the max depth for the calculation.

        By default, the density bandwidth smoothing parameter is set to `adjust=12`. This can be modified with
        `depadjust=INT`. The raw and smoothed profiles are output to `*.plots/*.raw.png` `*.plots/*.scdepth.png`
        to check smoothing if required. Additional checking plots are also output (see Outputs below).

        The full output of depths per position is output to `$BAM.fastmp` (or `$BAM.fastdep` if `quickdepth=T`). The
        single-copy is also output to `$BAM.fastmp.scdepth`. If `reduced=T` (the default) then the `fastmp` or `fastdep`
        file will have a `$BAM.busco.*` prefix and only include the sequences in the BUSCO table. By default, generation
        of the fastdep/fastmp data is performed by chunking up the assembly and creating temporary files in parallel
        (`tmpdir=PATH`). Sequences are batched in order such that each batch meets the minimum size criterion set by
        `depchunk=INT` (default 1Mbp). If `depchunk=0` then each sequence will be processed individually. This is not
        recommended for large, highly fragmented genomes. Unless `dev=T` or `debug=T`, the temporary files will be
        deleted once the final file is made. If DepthSizer crashed during the generation of the file, it should be
        possible to re-run and it will re-use existing temporary files.

        **NOTE:** To generate output that is compatible with [DepthKopy](https://github.com/slimsuite/depthkopy), run
        with `reduced=F`.

        ## Step 4: Read mapping adjustments

        The basic DepthSizer approach (`adjustmode=None`) assumes that the raw long read data has a 1:1 correspondence to
        the genomic DNA being sequenced, i.e. there is no contamination (including plastids) and no bias towards insertion
        or deletion read errors. As a consequence, the default genome size prediction is expected to be an over-estimate.

        DepthSizer will also calculate several adjusted estimation values that aim to provide the range in which the true
        genome size is expected to lie. Note that benchmarking and refinement of these adjustments is ongoing and will be
        expanded in future releases. For most use cases, it is anticipated that `CovBases` and `None` will provide the
        lower and upper bounds, whilst `IndelRatio` and `MapAdjust` should fall in between and be more accurate. (See
        notes below for each method.

        Currently, there are four `adjustmode` settings that can be output, in addition to two benchmarking calculations:

        * `None` : The purest DepthSizer mode makes no adjustment to total sequencing depth. (See Step 5.)
        * `CovBases` : This uses `samtools coverage` to calculate the total number of mapped bases as covered bases
        multiplied by the mean depth: `samtools coverage $BAM | grep -v coverage | awk  '{{sum += ($7 * $5)}} END {{print sum}}'`.
        Very big differences between `CovBases` and `None` may indicate a very incomplete assembly and/or an excess of
        contamination in the raw sequencing data. If BUSCO scores etc. indicate good completeness, it is advisable to
        carefully check the `read=FILELIST` data provided to DepthSizer.
        * `IndeRatio` : This mode extracts the CIGAR strings from the BAM file and sums up the insertions (`nI`),
        deletions (`nD`) and mapped bases (`nX`+`nM`+`n=`). The insertion:deletion ratio is then calculated as:
        (`I+X+M+=`)/(`D+X+M+=`). The goal here is to estimate whether the raw sequencing data is biased towards insertion
        or deletion errors. Insertion bias will inflate the apparent volume of sequencing. The total read volume is
        therefore adjusted by dividing by the indelratio. An insertion bias will decrease the estimated genome size,
        whereas a deletion bias will increase the prediction. To reduce issues caused by poor-quality regions of the
        assembly, mapped regions are reduced to `Complete` BUSCO genes (in a BED file, `$BED`):
        `samtools view -h -F 4 $BAM -L $BED | grep -v '^@' | awk '{print $6;}' | uniq`. The combined CIGAR counts are
        saved in `$BAM.indelratio.txt` to accelerate re-calculation.
        * `MapAdjust` : An earlier attempt to model insertion:deletion ratios, this combines the `CovBases` calculation
        of total assembly base coverage with `samtools fasta` to calculate the total number of bases in the mapped reads:
        `samtools view -hb -F 4 {0} | samtools fasta - | grep -v '^>' | wc | awk '{{ $4 = $3 - $2 }} 1' | awk '{{print $4}}'`
        MapAdjust calculates the general loss of raw sequencing during mapping as `CovBases/MapBases`. This aims to
        estimate the proportion of the raw sequencing data that contributed to the single copy read depth and is used as
        a multiplier for the total read volume. Extreme mapadjust ratios should be treated with caution and
        may indicate problems with the assembly and/or source data.
        * `Assembly` : In `benchmark=T` mode, the observed assembly size is output.
        * `MeanX` : In `benchmark=T` mode, the mean coverage is calculated as `CovBases`/`AssemblySize` and used in place
        of `scdepth` for the genome size estimation using the full sequencing volume.

        By default, DepthSizer will estimate genome sizes using `IndelRatio`, `CovBases` in addition to `None`. For
        speed, `CovBases` can be switched off with `covbases=F` and `IndelRatio` by setting `adjustmode=None`. The old
        `MapAdjust` calculation is not made by default, but can be switched on with `mapadjust=T`, `adjustmode=MapAdjust`,
        or `benchmark=T`. Setting `benchmark=T` will output all six estimates.

        **NOTE:** v1.5.0 expands the options to None/CovBases/IndelRatio/MapBases/MapAdjust/MapRatio/OldAdjust/OldCovBases.
        Details to follow.

        ## Step 5: Total read volume

        Genome size is estimated using the total combined sequencing length. This will be calculated from `reads=FILELIST`
        unless provided with `readbp=INT`. The number of bases for each input `$READFILE` is saved as
        `$READFILE.basecount.txt` and will be reloaded for future runs unless `force=T`.

        ## Step 6: Genome size prediction

        The final genome size is predicted based on the total (adjusted) combined sequencing length and the single-copy
        read depth, as: `EstGenomeSize`=`ReadBP`/`SCDepth`. Size predictions will be output for to `*.gensize.tdt` and
        as `#GSIZE` entries in the log file.

        # Outputs

        The main DepthSizer outputs are:

        * `*.gensize.tdt` = Main genome size prediction table.
        * `*.log` = DepthSizer log file with key steps and details of any errors or warnings generated.
        * `*.plots/` = Directory of PNG plots (see below)

        The primary output is the `*.gensize.tdt` table, which has the following fields:

        * `SeqFile` = Assembly file used for genome size prediction.
        * `DepMethod` = Depth estimation method used (`mpileup` or `depth`)
        * `Adjust` = Read mapping adjustment (see Step 4, above)
        * `ReadBP` = Total read volume (see Step 5, above)
        * `MapAdjust` = The relevant adjustment ratio.
        * `SCDepth` = The single copy read depth used in the prediction.
        * `EstGenomeSize` = Genome size prediction (bp).

        ## DepthSizer plots and additional tables

        The first time `SCDepth` is calculated (if not provided with `scdepth=NUM`), DepthSizer will also generate a
        number of plots for additional QC of results, which are output in `$BASE.plots/` as PNG files.

        First, the raw and smoothed read depth profiles will be output to:

        * `$BASE.plots/$BASE.raw.png` = raw depth profile
        * `$BASE.plots/$BASE.scdepth.png` = smoothed depth profile with SC depth marked

        In addition, violin plots will be generated for BUSCO `Complete` and `Duplicated` genes for the following `$STAT`
        values:

        * `$BASE.plots/$BASE.MeanX.png` = Mean depth of coverage
        * `$BASE.plots/$BASE.MedX.png` = Median depth of coverage
        * `$BASE.plots/$BASE.ModeX.png` = Pure Modal depth of coverage
        * `$BASE.plots/$BASE.DensX.png` = Smoothed density modal depth of coverage
        * `$BASE.plots/$BASE.CN.png` = Estimated copy number. (See [DepthKopy](https://github.com/slimsuite/depthkopy)
        for more details.)

        If `seqstats=T` then each assembly sequence will also be output in a `Sequences` violin plot for comparison.
        Each point in the plots is a separate gene or sequence. Values for individual genes/sequences are also output as
        a density scatter plot named `$BASE.plots/$BASE.$REGIONS.$STAT.png`, where

        * `$REGIONS` is the type of region plotted (`BUSCO` complete, `Duplicated`, or assembly `Sequences`).
        * `$STAT` is the output statistic: `MeanX`, `MedX`, `ModeX` or `DensX`.

        The BUSCO input file will also get a `*.regcnv.tsv` and `*.dupcnv.tsv` copy number prediction tables generated.
        See [DepthKopy](https://github.com/slimsuite/depthkopy) for more details.

        ## BAM file outputs

        In addition to the BAM file of mapped reads itself, DepthSizer may generate (or reload) the following files
        associated with the BAM generated/provided:

        * `$BAM.bai` = Samtools index
        * `$BAM.fastmp` = `samtools mpileup` read depths per position. This format apes KAT `*.cvg` format, with a
        fasta-style sequence header, followed on the next line by a depth value per position.
        * `$BAM.fastdep` = `samtools depth` read depths per position. This format apes KAT `*.cvg` format, with a
        fasta-style sequence header, followed on the next line by a depth value per position.
        * `$BAM.fast*.scdepth` = Single-copy read depth calculated for the above files.
        * `$BAM.indelratio.txt` = Compiled CIGAR strings for IndelRatio calculation.
        * `$BAM.mapratio.txt` = CovBases and MapBases numbers for MapAdjust calculation.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            self.setup()
            if self.getBool('Legacy'):
                return self.legacyDepthSizer()
            if self.getBool('DepOnly'):
                self.getFastDep()
                self.printLog('#EXIT','Terminating after fastdep creation (deponly=T).')
                return False
            return self.depthSizer()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def depthSizer(self):    ### Main Depth Copy method.
        '''
        Main Depth Copy method. This will generate the BAM if needed, generate a depth file, and then call the
        depthcopy.R script to perform the actual calculations.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DEPTHSIZER GENOME SIZE ESTIMATION', line='=')
            self.printLog('#MODE','Read mapping adjustment mode: {0}'.format(self.getStr('AdjustMode')))
            if self.getBool('Benchmark'):
                self.printLog('#MODE', 'Benchmarking mode: True')
            if not self.checkInput(busco=self.getNum('SCDepth')<=0,reads=self.getNum('ReadBP')<=0): return False
            self.seqinObj(summarise=True)
            ### ~ [2] Calculate key stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# v1.5.0 => None/CovBases/IndelRatio/MapBases/MapAdjust/MapRatio/OldAdjust/OldCovBases
            self.headLog('Calculate Single-Copy Read Depth', line='-')
            if not self.getSCDepth():     # This will generate the BAM file and depfile if needed
                self.printLog('#DEPTH','Cannot calculate genome size without single copy read depth')
                return False
            ## ~ [2a] IndelRatio ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Benchmark') or self.getStrLC('AdjustMode') in ['indelratio','mapratio']:
                self.indelRatio()  # Calculate with and without indel ratio correction
            ## ~ [2b] CovBases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Benchmark') or self.getBool('CovBases') or self.getStrLC('AdjustMode') == 'covbases':
                self.allBases()
            if self.getBool('Benchmark') or self.getStrLC('AdjustMode') == 'oldcovbases':
                self.covBases()
            if self.getBool('Benchmark') or self.getStrLC('AdjustMode') in ['mapbases','mapratio']:
                self.mapBases()
            ## ~ [2c] MapAdjust ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Benchmark') or self.getBool('MapAdjust') or self.getStrLC('AdjustMode') == 'mapadjust':
                self.mapAdjust()
            if self.getBool('Benchmark') or self.getStrLC('AdjustMode') == 'oldmapadjust':
                self.mapAdjust(allbases=False)
            ### ~ [3] Genome size prediction and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gensize = self.calculateGenomeSize(self.getStr('AdjustMode'),benchmark=self.getBool('Benchmark'))
            gdb = self.db('gensize')
            gfile = gdb.saveToFileName()
            self.infoLog('Genome size estimates output to: {0}'.format(gfile))
            self.infoLog('Actual genome size likely to fall within range of output predictions.')
            if self.getBool('Benchmark'):
                self.infoLog('Ignore "Assembly" and "MeanX" predictions: benchmarking only.')
            self.printLog('#FINAL','Estimated genome size (adjustmode={1}) = {0}'.format(rje_seqlist.dnaLen(gensize,dp=0,sf=4),self.getStr('AdjustMode')))
            if self.getStrLC('AdjustMode') == 'covbases':
                self.warnLog('Estimated genome size is likely to be an underestimate. See: {0}'.format(gfile))
            else:
                self.warnLog('Estimated genome size is likely to be an overestimate. See: {0}'.format(gfile))
            if self.getBool('Benchmark') or self.getBool('MapAdjust') or self.getStrLC('AdjustMode') == 'mapadjust':
                self.warnLog('If "MapAdjust" adjustment is very large, treat with caution.')
            return gensize
        except: self.errorLog('%s.depthSizer error' % self.prog())
#########################################################################################################################
    def legacyDepthSizer(self):    ### Legacy DepthSize method, using Diploidocus.
        '''
        Legacy DepthSize method, using Diploidocus.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.infoLog('Legacy mode disabled: run Diploidocus gensize mode with legacy=T')
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: DepthSizer Class                                                                                 #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: DepthSizer(mainlog,['basefile=depthsizer','reduced=T']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
