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
Module:       GapSpanner
Description:  Genome assembly gap long read support and reassembly tool
Version:      0.1.1
Last Edit:    09/09/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    GapSpanner is a wrapper for the gap-spanning and reassembly methods introduced by Diploidocus v0.12.0 and v0.13.0.
    GapSpanner needs a genome assembly (fasta format, `seqin=FILE`) and a set of long read (ONT, PacBio or HiFi) data
    for the assembly (`reads=FILELIST` and `readtype=LIST`).

    First, all gaps in the assembly at least 10bp (`mingap=INT`) are identified and saved in a table
    (`$SEQBASE.gaps.tdt`). Long read data is mapped onto the assembly using [Minimap2](https://github.com/lh3/minimap2)
    and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF
    file. In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST`
    (default 100bp, 1kb and 5kb) will also be calculated. (If the end of a sequence is reached before the end of the read,
    this will also count as flanking. Such regions can be identified using the `MaxFlank5` and `MaxFlank3` fields, which
    identify the maximum distance 5' and 3' that a read can span due to sequence length constraints.) Gap spanning output
    will be saved to `*.checkpos.tdt`. IDs for reads spanning each gap will also be saved to `$BASEFILE_spanid/`, with
    each gap will be named: `seqname.start-end`.

    The `gapass` and `gapfill` run modes will then attempt to re-assemble any gaps spanned by 2+ (using `mingapspan=INT`)
    the [Flye](https://github.com/fenderglass/Flye) assembler. First, long reads are re-mapped to generate BAM output and
    then [samtools](http://www.htslib.org/) (`samtools view` and `samtools fasta`) used to extract the reads that span
    the gap into a file for re-assembly. Gap assemblies will then be farmed out, running `forks=X`
    Flye processes at a time, with each Flye assembly using `subforks=INT` threads. This can take some time, and
    GapSpanner will terminate all running Flye assemblies if none finish within a 10 hour window (`killforks=INT`). This
    can be switched off with `killforks=0`. Assemblies are performed in `$BASEFILE__gapassemble/`.

    NOTE: Not all gaps will be successfully re-assembled at this point. See the log file for details.

    With `gapfill` run mode, re-assembled gap regions are compiled into `*.assembledgaps.fasta` single file and then
    mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
    must be at least 500bp. This stringency can be controlled with `minlocid=PERC` and `minloclen=INT`. Local hits
    are reduced to unique coverage of the assembly sequences and saved to `*.gapfiller.tdt`. Gaps are filled if one of
    the two conditions are met:

    1. A single local alignment spans an entire gap.
    2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

    In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
    re-assembly contig region. For a pair of hits, the region between the two hits is replaced. The updated sequences are
    output to `*.fillcheck.fasta` with modified sequences given a `-Xfix` suffix.

    Long reads are then mapped on to the updated assembly and the [Diploidocus](https://github.com/slimsuite/diploidocus)
    `regcheck` mode is run a second time to check for mapped reads that span the filled gaps. If the single copy read
    depth is given with `scdepth=X` then the estimated copy number of the replaced region will also be calculated.
    If no reads span the replaced region, the gap-filling will be reversed. (NOTE: this may cause issues when the full
    assembly is inserted, and future releases may instead use [DepthCharge](https://github.com/slimsuite/depthcharge) to
    identify bad gap-filling.) Such sequences will have an additional `-Xrev` suffix.

    The final gap-filled assembly is output to `*.gapfill.fasta`.

    Output will be saved to files with a prefix set by `basefile=X` (default named after the `seqin=FILE` prefix).
    For more details, see the Diploidocus documentation.

GapSpanner run modes:

    ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

    This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
    existing) and then runs the Diploidocus read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
    with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
    This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
    In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
    sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
    using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
    due to sequence length constraints.

    Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
    `seqname.start-end`.

    ---

    ### ~ Assembly gap re-assembly [runmode=gapass] ~ ###

    In addition to the `gapspan` analysis, reads identified as spanning each gap are extracted and assembled using `flye`
    in a `$BASEFILE__gapassemble/` output directory.

    ---

    ### ~ Re-assembled gap-filling [runmode=gapfill] ~ ###

    In addition to the `gapspan` and `gapass` outputs, re-assembled gap regions are compiled into a single file and then
    mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
    are reduced to unique coverage of the assembly sequences. Gaps are filled if one of the two conditions are met:

    1. A single local alignment spans an entire gap.
    2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

    In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
    re-assembly contig region. For a pair of hits, the region between the two hits is replaced.

    ---

    ## Dependencies

    [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the environment `$PATH` or given to
    GapSpanner with the `minimap2=PROG` setting. For `gapass` or `gapfill` run modes, [samtools](http://www.htslib.org/)
    and the [Flye](https://github.com/fenderglass/Flye) assembler need to be installed.

    To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For GapSpanner documentation, run with `dochtml=T` and read the `*.docs.html` file generated.


Commandline:
    ### ~ Main GapSpanner run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    runmode=X       : GapSpanner run mode (gapspan/gapass/gapfill) [gapspan]
    basefile=FILE   : Root of output file names [gapspanner or $SEQIN basefile]
    summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
    genomesize=INT  : Haploid genome size (bp) [0]
    paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    dochtml=T/F     : Generate HTML GapSpanner documentation (*.docs.html) instead of main run [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
    ### ~ Gaps spanning and reassembly options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
    subforks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
    mingapspan=INT  : Minimum number of reads spanning a gap in order to re-assemble [2]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [500]
    scdepth=NUM     : Single copy ("diploid") read depth. If zero, will not calculate CNV for filled gaps [0]
    ### ~ Re-assembly forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
    logfork=T/F     : Whether to log forking in main log [False]
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
import diploidocus
import rje, rje_obj, rje_rmd
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation. Wrapper for Diploidocus v0.14.0 gapspan modes.
    # 0.1.0 - Modified sequence summary output.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [X] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Fixed missing basefile bug.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('GapSpanner', '0.1.1', 'September 2021', '2021')
    description = 'Genome assembly gap long read support and reassembly tool'
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
### SECTION II: GapSpanner Class                                                                                        #
#########################################################################################################################
class GapSpanner(rje_obj.RJE_Object):
    '''
    GapSpanner Class. Author: Rich Edwards (2021).

    Str:str
    - SeqIn=FILE      : Input sequence assembly [None]

    Bool:boolean
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]

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
        self.strlist = ['SeqIn']
        self.boollist = ['DocHTML']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
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
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['SeqIn'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML'])  # True/False Booleans
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
        # GapSpanner: Genome assembly gap long read support and reassembly tool

        GapSpanner is a wrapper for the gap-spanning and reassembly methods introduced by Diploidocus v0.12.0 and v0.13.0.
        GapSpanner needs a genome assembly (fasta format, `seqin=FILE`) and a set of long read (ONT, PacBio or HiFi) data
        for the assembly (`reads=FILELIST` and `readtype=LIST`).

        First, all gaps in the assembly at least 10bp (`mingap=INT`) are identified and saved in a table
        (`$SEQBASE.gaps.tdt`). Long read data is mapped onto the assembly using [Minimap2](https://github.com/lh3/minimap2)
        and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF
        file. In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST`
        (default 100bp, 1kb and 5kb) will also be calculated. (If the end of a sequence is reached before the end of the read,
        this will also count as flanking. Such regions can be identified using the `MaxFlank5` and `MaxFlank3` fields, which
        identify the maximum distance 5' and 3' that a read can span due to sequence length constraints.) Gap spanning output
        will be saved to `*.checkpos.tdt`. IDs for reads spanning each gap will also be saved to `$BASEFILE_spanid/`, with
        each gap will be named: `seqname.start-end`.

        The `gapass` and `gapfill` run modes will then attempt to re-assemble any gaps spanned by 2+ (using `mingapspan=INT`)
        the [Flye](https://github.com/fenderglass/Flye) assembler. First, long reads are re-mapped to generate BAM output and
        then [samtools](http://www.htslib.org/) (`samtools view` and `samtools fasta`) used to extract the reads that span
        the gap into a file for re-assembly. Gap assemblies will then be farmed out, running `forks=X`
        Flye processes at a time, with each Flye assembly using `subforks=INT` threads. This can take some time, and
        GapSpanner will terminate all running Flye assemblies if none finish within a 10 hour window (`killforks=INT`). This
        can be switched off with `killforks=0`. Assemblies are performed in `$BASEFILE__gapassemble/`.

        NOTE: Not all gaps will be successfully re-assembled at this point. See the log file for details.

        With `gapfill` run mode, re-assembled gap regions are compiled into `*.assembledgaps.fasta` single file and then
        mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
        must be at least 500bp. This stringency can be controlled with `minlocid=PERC` and `minloclen=INT`. Local hits
        are reduced to unique coverage of the assembly sequences and saved to `*.gapfiller.tdt`. Gaps are filled if one of
        the two conditions are met:

        1. A single local alignment spans an entire gap.
        2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

        In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
        re-assembly contig region. For a pair of hits, the region between the two hits is replaced. The updated sequences are
        output to `*.fillcheck.fasta` with modified sequences given a `-Xfix` suffix.

        Long reads are then mapped on to the updated assembly and the [Diploidocus](https://github.com/slimsuite/diploidocus)
        `regcheck` mode is run a second time to check for mapped reads that span the filled gaps. If the single copy read
        depth is given with `scdepth=X` then the estimated copy number of the replaced region will also be calculated.
        If no reads span the replaced region, the gap-filling will be reversed. (NOTE: this may cause issues when the full
        assembly is inserted, and future releases may instead use [DepthCharge](https://github.com/slimsuite/depthcharge) to
        identify bad gap-filling.) Such sequences will have an additional `-Xrev` suffix.

        The final gap-filled assembly is output to `*.gapfill.fasta`.

        Output will be saved to files with a prefix set by `basefile=X` (default named after the `seqin=FILE` prefix).
        For more details, see the Diploidocus documentation.

        ---

        # Running GapSpanner

        GapSpanner is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/gapspanner.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [GapSpanner git repo](https://github.com/slimsuite/gapspanner), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [GapSpanner git repo](https://github.com/slimsuite/gapspanner)
        for running on example data.

        ## Dependencies

        [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the environment `$PATH` or given to
        GapSpanner with the `minimap2=PROG` setting. For `gapass` or `gapfill` run modes, [samtools](http://www.htslib.org/)
        and the [Flye](https://github.com/fenderglass/Flye) assembler need to be installed.

        To generate documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For GapSpanner documentation, run with `dochtml=T` and read the `*.docs.html` file generated.

        ## Commandline options

        A list of commandline options can be generated at run-time using the `-h` or `help` flags. Please see the general
        [SLiMSuite documentation](http://slimsuite.blogspot.com/2013/08/command-line-options.html) for details of how to
        use commandline options, including setting default values with **INI files**.

        ```
        ### ~ Main GapSpanner run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        runmode=X       : GapSpanner run mode (gapspan/gapass/gapfill) [gapspan]
        basefile=FILE   : Root of output file names [gapspanner or $SEQIN basefile]
        summarise=T/F   : Whether to generate and output summary statistics sequence data before and after processing [True]
        genomesize=INT  : Haploid genome size (bp) [0]
        paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
        bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        dochtml=T/F     : Generate HTML GapSpanner documentation (*.docs.html) instead of main run [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~ Gaps spanning and reassembly options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkflanks=LIST: List of lengths flanking check regions that must also be spanned by reads [0,100,1000,5000]
        subforks=INT    : Number of forks for assembly subproccesses during gapfill and gapass modes [1]
        mingapspan=INT  : Minimum number of reads spanning a gap in order to re-assemble [2]
        minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [0]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [500]
        scdepth=NUM     : Single copy ("diploid") read depth. If zero, will not calculate CNV for filled gaps [0]
        ### ~ Re-assembly forking options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        forks=X         : Number of parallel sequences to process at once [0]
        killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
        killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [False]
        logfork=T/F     : Whether to log forking in main log [False]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        ---

        # GapSpanner run modes

        ### ~ Assembly gap read-spanning analysis [runmode=gapspan] ~ ###

        This mode first identifies all the gaps in an assembly (`seqin=FILE`) (using SeqList `gapstats` or `$SEQBASE.gaps.tdt` if pre-
        existing) and then runs the Diploidocus read spanning analysis (`runmode=regcheck`) with `regcnv=F`. Long read data, given
        with the `reads=FILELIST` and `readtype=LIST` options, are mapped onto the assembly using minimap2 to generate a PAF file.
        This is then parsed and reads spanning each gap are identified based on their positions and the target start and end positions in the PAF file.
        In addition to absolute spanning of regions, reads spanning a region +/- distances set by `checkflanks=LIST` will also be calculated. If the end of a
        sequence is reached before the end of the read, this will also count as flanking. Such regions can be identified
        using the `MaxFlank5` and `MaxFlank3` fields, which identify the maximum distance 5' and 3' that a read can span
        due to sequence length constraints.

        Spanning `spanid` output is also generated for each gap and saved in `$BASEFILE_spanid`. Each gap will be named:
        `seqname.start-end`.

        ---

        ### ~ Assembly gap re-assembly [runmode=gapass] ~ ###

        In addition to the `gapspan` analysis, reads identified as spanning each gap are extracted and assembled using `flye`
        in a `$BASEFILE__gapassemble/` output directory.

        ---

        ### ~ Re-assembled gap-filling [runmode=gapfill] ~ ###

        In addition to the `gapspan` and `gapass` outputs, re-assembled gap regions are compiled into a single file and then
        mapped back on the original assembly using Minimap2, with tabulated hit output into `$BASEFILE__gapfill/`. Local hits
        are reduced to unique coverage of the assembly sequences. Gaps are filled if one of the two conditions are met:

        1. A single local alignment spans an entire gap.
        2. A pair of concordant local alignments from the same re-assembly contig (same orientation) flank an entire gap.

        In the case of a single spanning local alignment, the entire assembly region is replaced by the corresponding
        re-assembly contig region. For a pair of hits, the region between the two hits is replaced.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dipobj = diploidocus.Diploidocus(self.log,['dna=T','runmode=gapspan']+self.cmd_list)
            dipobj.setup()
            if not dipobj.getStrLC('RunMode') in ['gapspan','gapass','gapfill']:
                raise ValueError('RunMode "{0}" not recognised'.format(dipobj.getStr('RunMode')))
            self.printLog('#MODE',dipobj.getStrLC('RunMode'))
            dipobj.seqinObj(summarise=False)
            dipobj.gapSpan()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
                else: raise ValueError('seqin=FILE must be set')
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
### End of SECTION II: GapSpanner Class                                                                                 #
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
    try: GapSpanner(mainlog,cmd_list).run()

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
