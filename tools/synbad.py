#!/usr/bin/python

# See below for name and description
# Copyright (C) 2020 Richard J. Edwards <dr.r.edwards@icloud.com>
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
Module:       synbad
Description:  Synteny-based scaffolding adjustment
Version:      0.5.0
Last Edit:    20/01/21
GitHub:       https://github.com/slimsuite/synbad
Copyright (C) 2020  Richard J. Edwards - See source code for GNU License Notice

Function:
    SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
    between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

    Synbad will use or create:

    1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
    spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
    correspond to misassemblies.

    2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

    Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
    by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
    should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

    Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
    sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
    new fields with 5', 3' or Both.

    The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
    corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
    If the local hits are against two different sequences from the other genome, the two sequence names will be
    entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
    small gaps), `SynSpan` or `GapSpan` will have a value of zero.

    Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:

    * `Aligned` = Gap is found in the middle of a local alignment to the Hit
    * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
    * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
    * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
    * `Duplication` = Overlapping flanking hits on the same strand.
    * `Inversion` = Flanking hits are on alternative strands.
    * `Translocation` = `SynSpan` indicates matches are on different contigs.
    * `Terminal` = Gap is between a local alignment and the end of the query sequence.
    * `Null` = No mapping between genomes for that gap.

    If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
    `minreadspan=INT` reads span the gap.

    A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
    where possible.

Dependencies:
    SynBad needs Minimap2 installed. For `gapass` gap mode, Flye also needs to be installed. To generate
    documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

        export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

    For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.

Commandline:
    ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    basefile=X      : Prefix for output files [synbad]
    gablam=X        : Optional prefix for GABLAM search [defaults to basefile=X]
    gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
    maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [1]
    maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [10000]
    chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    ### ~ Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
    fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Breakpoint,Inversion,Translocation]
    minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
    ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
    reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
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
import rje, rje_obj, rje_db, rje_rmd, rje_seqlist
import diploidocus, gablam
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial working version without fragment=T implementation.
    # 0.1.1 - Added minlocid=PERC and minloclen=INT.
    # 0.1.2 - Added additional translocation skipping for SynTrans rating.
    # 0.1.3 - Modified the SynBad classification text.
    # 0.1.4 - Modified code to be able to run without long read mapping. Added dochtml output.
    # 0.2.0 - Added fragment=T output.
    # 0.3.0 - Added chromosome scaffold Translocation restriction.
    # 0.4.0 - Added an Duplication rating in place of Breakpoint for overlapping flanking hits; added top sequence pairs.
    # 0.5.0 - Added HiFi read type. Changed default to gapmode=gapass.
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
    # [Y] : Add chr1 and chr2 chromosome prefix identifers (PAFScaff prefixes) to restrict Translocations to chromosomes.
    # [ ] : Add GapGFF function from rje_genomics.
    # [ ] : Add fullcollapse=T/F option to collapse entire (pairwise) matching regions between gaps.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SynBad', '0.5.0', 'January 2021', '2020')
    description = 'Synteny-based scaffolding adjustment'
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
### SECTION II: New Class                                                                                               #
#########################################################################################################################
class SynBad(rje_obj.RJE_Object):
    '''
    Class. Author: Rich Edwards (2015).

    Str:str
    - BAM1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
    - BAM2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
    - Chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
    - Chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
    - GABLAM=X        : Optional prefix for GABLAM search [defaults to basefile=X]
    - GapMode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
    - Genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
    - Genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
    - PAF1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
    - PAF2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]

    Bool:boolean
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]
    - Fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]

    Int:integer
    - MaxSynSkip=INT  : Maximum number of local alignments to skip for SynTrans classification [1]
    - MaxSynSpan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [10000]
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
    - MinReadSpan=INT     : Min number of Span0 reads in gaps table to prevent fragmentation [1]

    Num:float
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]

    File:file handles with matching str filenames
    
    List:list
    - Reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - Reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    - ReadType2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database object
    - Genome1 = Genome 1 SeqList Object
    - Genome2 = Genome 2 SeqList Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM1','BAM2','Chr1','Chr2','GABLAM','GapMode','Genome1','Genome2','PAF1','PAF2']
        self.boollist = ['DocHTML','Fragment']
        self.intlist = ['MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan']
        self.numlist = ['MinLocID']
        self.filelist = []
        self.listlist = ['Reads1','Reads2','ReadType1','ReadType2']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'GapMode':'gapspan'})
        self.setBool({'Fragment':False})
        self.setInt({'MaxSynSkip':1,'MaxSynSpan':10000,'MinLocLen':1000,'MinReadSpan':1})
        self.setNum({'MinLocID':50.0})
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
                self._cmdReadList(cmd,'str',['Chr1','Chr2','GABLAM','GapMode'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['BAM1','BAM2','Genome1','Genome2','PAF1','PAF2'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['DocHTML','Fragment'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxSynSkip','MaxSynSpan','MinLocLen','MinReadSpan'])   # Integers
                self._cmdReadList(cmd,'perc',['MinLocID']) # Percentage
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ReadType1','ReadType2'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Reads1','Reads2']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # SynBad:  Synteny-based scaffolding adjustment

        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.

        SynBad will use or create:

        1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
        spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
        correspond to misassemblies.

        2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.

        Pairwise hits between the genomes are filtered according to the `minlocid=PERC` and `minloclen=INT` criteria, which
        by default limits hits to be at least 1kb and 50% identity. Note that this is applied after GABLAM has run, so it
        should be possible to re-run with more relaxed constraints and re-use the GABLAM tables.

        Next, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
        sorted and `QryGap` and `SbjGap` fields added. Any local alignments with flanking hits are then flagged in these
        new fields with the number of flanking 5', 3' gaps.

        The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
        corresponding local hits on the Qry and Sbj genomes. If there is also an inversion, `SynSpan` will be negative.
        If the local hits are against two different sequences from the other genome, the two sequence names will be
        entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
        small gaps), `SynSpan` or `GapSpan` will have a value of zero.

        Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:

        * `Aligned` = Gap is found in the middle of a local alignment to the Hit
        * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
        * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
        * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
        * `Duplication` = Overlapping flanking hits on the same strand.
        * `Inversion` = Flanking hits are on alternative strands.
        * `Translocation` = `SynSpan` indicates matches are on different contigs.
        * `Fragmentation` = `SynSpan` indicates matches are on different contigs, 1+ of which is not a chromosome scaffold.
        * `Terminal` = Gap is between a local alignment and the end of the query sequence.
        * `Null` = No mapping between genomes for that gap.

        If `chr1=X` and/or `chr2=X` chromosome scaffold prefixes are provided then `Translocation` will be restricted to
        matches between two different chromosome scaffolds. Matches including one or more non-chromosome scaffolds will
        be classed as `Fragmentation`.

        If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
        `minreadspan=INT` reads span the gap.

        A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
        where possible.

        ## Dependencies

        SynBad needs Minimap2 installed. For `gapass` gap mode, Flye also needs to be installed. To generate
        documentation with `dochtml`, R will need to be installed and a pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        For full documentation of the SynBad workflow, run with `dochtml=T` and read the `*.docs.html` file generated.


        ## Commandline options

        ```
        ### ~ Main SynBad run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        genome1=FILE    : Genome assembly used as the query in the GABLAM searches []
        genome2=FILE    : Genome assembly used as the searchdb in the GABLAM searches []
        basefile=X      : Prefix for output files [synbad]
        gablam=X        : Optional prefix for GABLAM search [defaults to basefile=X]
        gapmode=X       : Diploidocus gap run mode (gapspan/gapass) [gapspan]
        minloclen=INT   : Minimum length for aligned chunk to be kept (local hit length in bp) [1000]
        minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
        maxsynskip=INT  : Maximum number of local alignments to skip for SynTrans classification [1]
        maxsynspan=INT  : Maximum distance (bp) between syntenic local alignments to count as syntenic [10000]
        chr1=X          : PAFScaff-style chromosome prefix for Genome 1 to distinguish Translocation from Fragmentation []
        chr2=X          : PAFScaff-style chromosome prefix for Genome 2 to distinguish Translocation from Fragmentation []
        ### ~ Fragmentation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        fragment=T/F    : Whether to fragment the assembly at gaps marked as non-syntenic [False]
        fragtypes=LIST  : List of SynBad ratings to trigger fragmentation [Breakpoint,Inversion,Translocation]
        minreadspan=INT : Min number of Span0 reads in gaps table to prevent fragmentation [1]
        ### ~ Additional input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        bam1=FILE       : Optional BAM file of long reads mapped onto assembly 1 [$BASEFILE1.bam]
        paf1=FILE       : Optional PAF file of long reads mapped onto assembly 1 [$BASEFILE1.paf]
        reads1=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype1=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        bam2=FILE       : Optional BAM file of long reads mapped onto assembly 2 [$BASEFILE2.bam]
        paf2=FILE       : Optional PAF file of long reads mapped onto assembly 2 [$BASEFILE2.paf]
        reads2=FILELIST : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype2=LIST  : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        ### ~ Additional output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        ```

        ---

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return self.docHTML()
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.synBad()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if self.getStrLC('GapMode') not in ['gapspan','gapass','gapfill']:
                self.printLog('#MODE','GapMode "{0}" not recognised. Setting gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                self.setStr({'GapMode':'gapspan'})
            if self.getStrLC('GapMode') not in ['gapspan','gapass']:
                self.printLog('#MODE','GapMode "{0}" no longer supported. Please run Diploidocus separately or use gapmode=gapspan.'.format(self.getStrLC('GapMode')))
                if rje.yesNo('Set gapmode=gapspan?'):
                    self.setStr({'GapMode':'gapspan'})
                    self.printLog('#MODE','Setting gapmode=gapspan.')
                else: return False
            ### ~ [2] Check for files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.exists(self.getStr('Genome1')): raise IOError('Genome1 file "{0}" not found!'.format(self.getStr('Genome1')))
            if not rje.exists(self.getStr('Genome2')): raise IOError('Genome2 file "{0}" not found!'.format(self.getStr('Genome2')))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def docHTML(self):  ### Generate the Diploidocus Rmd and HTML documents.                                        # v0.1.0
        '''Generate the Diploidocus Rmd and HTML documents.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            info = self.log.obj['Info']
            prog = '%s V%s' % (info.program,info.version)
            rmd = rje_rmd.Rmd(self.log,self.cmd_list)
            rtxt = rmd.rmdHead(title='%s Documentation' % prog,author='Richard J. Edwards',setup=True)
            #!# Replace this with documentation text?
            rtxt += string.replace(self.run.__doc__,'\n        ','\n')
            rtxt += '\n\n<br>\n<small>&copy; 2020 Richard Edwards | richard.edwards@unsw.edu.au</small>\n'
            rmdfile = '%s.docs.Rmd' % self.baseFile()
            open(rmdfile,'w').write(rtxt)
            self.printLog('#RMD','RMarkdown SynBad documentation output to %s' % rmdfile)
            rmd.rmdKnit(rmdfile)
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <3> ### Main SynBad Methods                                                                                     #
#########################################################################################################################
    def synBad(self):   ### Main SynBad run method
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            basefile = self.baseFile()
            dcmd = ['bam=','paf=','reads=','readtype=ont']
            chr1 = self.getStrLC('Chr1')
            if chr1: chr1 = self.getStr('Chr1')
            chr2 = self.getStrLC('Chr2')
            if chr2: chr2 = self.getStr('Chr2')
            #!# Add checking of seqnames read in for gaps and local hits and warn if none match #!#

            ### ~ [1] Run Diploidocus on each genome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('DIPLOIDOCUS GAP ANALYSIS',line='=')
            # SynBad will use or create:
            #
            # 1. A table of gap positions for each assembly (seqname, start, end). This can optionally have long reads mapped and
            # spanning coverage calculated for each gap using Diploidocus. Gaps without spanning long reads are more likely to
            # correspond to misassemblies.
            #
            # ==> tigersnake.v2.7.pafscaff.checkpos.tdt <==
            # #       seqname start   end     seqlen  gaplen  gap     MaxFlank5       MaxFlank3       Span0   Span100 Span1000        Span5000
            # 209     NSCUCHR1.01     57638   57737   341225390       100     NSCUCHR1.01.57638-57737 57637   341167653       1       1       1       0
            # 2       NSCUCHR1.01     103648  104147  341225390       500     NSCUCHR1.01.103648-104147       103647  341121243       14      14      12      5

            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base1 = rje.baseFile(self.getStr('Genome1'),strip_path=True)
            cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base1)]
            rundip = False
            if self.getStrLC('BAM1'): cmd1.append('bam={0}'.format(self.getStr('BAM1')))
            if self.getStrLC('PAF1'): cmd1.append('paf={0}'.format(self.getStr('PAF1'))); rundip = True
            if self.list['Reads1']: cmd1.append('reads={0}'.format(','.join(self.list['Reads1']))); rundip = True
            if self.list['ReadType1']: cmd1.append('readtype={0}'.format(','.join(self.list['ReadType1'])))
            dip1 = diploidocus.Diploidocus(self.log,self.cmd_list+dcmd+cmd1)
            cdb1 = None
            if rundip:
                dip1.run()
                cdb1 = dip1.db('checkpos')
                cdb1.baseFile(basefile)
                cdb1.setStr({'Name':'qrygap'})
                self.db().list['Tables'].append(cdb1)
            else:
                self.printLog('#READS','No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base1))
                if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=base1,log=self.log):
                    seqcmd = self.cmd_list + cmd1 + ['summarise=T','gapstats=T','raw=F','dna=T']
                    self.obj['SeqList1'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
                cdb1 = db.addTable('%s.gaps.tdt' % base1,mainkeys=['seqname','start','end'],name='qrygap',ignore=[],expect=True)
                cdb1.addField('Span0',evalue=0)
            #checkpos1 = '{0}.checkpos.tdt'.format(base1)
            #cdb1 = db.addTable(checkpos1,mainkeys=['seqname','start','end'],name='qrygap',ignore=[],expect=True)
            cdb1.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            cdb1.addField('GapSpan',evalue='.')
            cdb1.addField('SynSpan',evalue='.')
            cdb1.addField('SynBad',evalue='Null')
            cdb1.index('seqname')

            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            base2 = rje.baseFile(self.getStr('Genome2'),strip_path=True)
            cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                    'runmode={0}'.format(self.getStrLC('GapMode')),
                    'basefile={0}'.format(base2)]
            rundip = False
            if self.getStrLC('BAM2'): cmd2.append('bam={0}'.format(self.getStr('BAM2')))
            if self.getStrLC('PAF2'): cmd2.append('paf={0}'.format(self.getStr('PAF2'))); rundip = True
            if self.list['Reads2']: cmd2.append('reads={0}'.format(','.join(self.list['Reads2']))); rundip = True
            if self.list['ReadType2']: cmd2.append('readtype={0}'.format(','.join(self.list['ReadType2'])))
            dip2 = diploidocus.Diploidocus(self.log,self.cmd_list+dcmd+cmd2)
            cdb2 = None
            if rundip:
                dip2.run()
                cdb2 = dip2.db('checkpos')
                cdb2.baseFile(basefile)
                cdb2.setStr({'Name':'hitgap'})
                self.db().list['Tables'].append(cdb2)
            else:
                self.printLog('#READS','No reads or PAF mapping provided for {0}: no read spanning analysis'.format(base2))
                if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=base2,log=self.log):
                    seqcmd = self.cmd_list + cmd2 + ['summarise=T','gapstats=T','raw=F','dna=T']
                    self.obj['SeqList2'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
                    # dip2.cmd_list.append('gapstats')    # Try to run automatically if possible
                    # seqin = dip2.seqinObj()
                    # if not rje.checkForFiles(filelist=['.gaps.tdt'],basename=base2,log=self.log):
                    #     seqin.setBool({'Raw':False,'GapStats':True,'DNA':True})
                    #     seqin.str['SeqType'] = 'dna'
                    #     seqin.summarise()
                cdb2 = db.addTable('%s.gaps.tdt' % base2,mainkeys=['seqname','start','end'],name='hitgap',ignore=[],expect=True)
                cdb2.addField('Span0',evalue=0)
            #checkpos2 = '{0}.checkpos.tdt'.format(base2)
            #cdb2 = db.addTable(checkpos2,mainkeys=['seqname','start','end'],name='hitgap',ignore=[],expect=True)
            cdb2.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            cdb2.addField('GapSpan',evalue='.')
            cdb2.addField('SynSpan',evalue='.')
            cdb2.addField('SynBad',evalue='Null')
            cdb2.index('seqname')

            ### ~ [2] GABLAM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GABLAM SEARCH',line='=')
            gabbase = basefile
            if self.getStrLC('GABLAM'): gabbase = self.getStr('GABLAM')
            self.printLog('#GABLAM','GABLAM output basefile: {0}'.format(gabbase))
            # 2. The qryunique and hitunique local hits tables from a GABLAM run using Minimap2.
            # ==> tiger.v.najna.XXXunique.tdt <==
            # Qry     Hit     AlnNum  BitScore        Expect  Length  Identity        Positives       QryStart        QryEnd  SbjStart        SbjEnd
            # In each case, Qry is Genome1 and Hit is Genome2
            quniq = '{0}.qryunique.tdt'.format(gabbase)
            huniq = '{0}.hitunique.tdt'.format(gabbase)
            ## ~ [2a] Run GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.force() or not rje.exists(quniq) or not rje.exists(huniq):
                gabcmd = ['seqin={0}'.format(self.getStr('Genome1')),'searchdb={0}'.format(self.getStr('Genome2')),'mapper=minimap','minlocid=0','minloclen=0','basefile={0}'.format(gabbase)]
                gabobj = gablam.GABLAM(self.log,self.cmd_list+gabcmd)
                gabobj.gablam()
            ## ~ [2b] Load tables and reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                ufile = '{0}.{1}unique.tdt'.format(gabbase,qh.lower())
                udb = db.addTable(ufile,mainkeys=['Qry','Hit','AlnNum'],name=qh.lower(),ignore=[],expect=True)
                #udb = db.addTable(ufile,mainkeys=[qh,'{0}Start'.format(qh),'{0}End'.format(qh)],name=qh.lower(),ignore=[],expect=True)
                udb.dataFormat({'AlnNum':'int','Length':'int','Identity':'int','QryStart':'int','QryEnd':'int','SbjStart':'int','SbjEnd':'int'})
                lenx = 0; idx = 0
                for entry in udb.entries():
                    if entry['Length'] < self.getInt('MinLocLen'): udb.dropEntry(entry); lenx += 1
                    elif (100.0 * entry['Identity'] / entry['Length']) < self.getNum('MinLocID'): udb.dropEntry(entry); idx += 1
                self.printLog('#MINCUT','Dropped %s entries < %s bp and %s < %.1f%% identity' % (rje.iStr(lenx),rje.iStr(self.getInt('MinLocLen')),rje.iStr(idx),self.getNum('MinLocID')))
                udb.addField('HitStart',evalue=0)
                udb.addField('HitEnd',evalue=0)
                udb.addField('Strand',evalue='+')
                for entry in udb.entries():
                    if entry['SbjStart'] > entry['SbjEnd']:
                        entry['HitStart'] = entry['SbjEnd']
                        entry['HitEnd'] = entry['SbjStart']
                        entry['Strand'] = '-'
                    else:
                        entry['HitStart'] = entry['SbjStart']
                        entry['HitEnd'] = entry['SbjEnd']
                udb.newKey([qh,'{0}Start'.format(qh),'{0}End'.format(qh)])
                udb.setFields(['Qry','QryStart','QryEnd','Hit','HitStart','HitEnd','Strand','AlnNum','Length','Identity'])
                udb.addField('Unaligned',evalue=0)
                udb.addField('QryGap',evalue='')
                udb.addField('HitGap',evalue='')
                udb.addField('Qry5',evalue='')
                udb.addField('Qry3',evalue='')
                udb.addField('Hit5',evalue='')
                udb.addField('Hit3',evalue='')

            ### ~ [3] Process tables, adding Gap and local hit distance information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('SYNBAD GAP MAPPING',line='=')
            self.printLog('#NOTE','Genome1 ({0}) is always Qry and Genome2 ({1}) is always Hit'.format(base1,base2))
            # First, all gap positions are read in along with the local hits tables. For each genome, the local hit tables are
            # sorted and `QryGap` and `HitGap` fields added. Any local alignments with flanking hits are then flagged in these
            # new fields with the number of flanking 5', 3' gaps indicated by < and > characters.
            #
            # The gap tables will also be updated with `GapSpan` and `SynSpan` fields that have the distance between the
            # corresponding local hits on the Qry and Hit genomes. If there is also an inversion, `SynSpan` will be negative.
            # If the local hits are against two different sequences from the other genome, the two sequence names will be
            # entered in the `SynSpan` field instead. If the gap is in the middle of local hit (likely to be true only for
            # small gaps), `SynSpan` or `GapSpan` will have a value of zero.
            ## ~ [3a] Work through query and hit tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #self.debug(self.db().tableNames())
            for qh in ('Qry','Hit'):
                qry = qh
                hit = {'Qry':'Hit','Hit':'Qry'}[qry]
                hitchr = {'Qry':chr1,'Hit':chr2}[hit]
                gapdb = self.db('{0}gap'.format(qh.lower()))
                #self.debug(gapdb)
                altdb = {cdb1:cdb2,cdb2:cdb1}[gapdb]
                locdb = self.db(qh.lower())
                locdb.index(qry)
                locdb.index(hit)
                # First, deal with QryGaps
                for seqname in gapdb.indexKeys('seqname'):
                    if seqname not in locdb.index(qry):
                        self.printLog('#ORPHAN','No alignments found for {0}'.format(seqname))
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = gapdb.index('seqname')[seqname][0:]
                    loc = [(seqname,0,0)] + locdb.index(qry)[seqname] + [(seqname,-1,-1)]
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        entry5 = None
                        entry3 = None
                        if loc[0][1] > 0:
                            entry5 = locdb.data(loc[0])
                            entry5['{0}Gap'.format(qry)] += '>'
                            if not entry5['{0}3'.format(qry)]: entry5['{0}3'.format(qry)] = []
                            entry5['{0}3'.format(qry)].append(gap[0])
                        if loc[1][1] > 0:
                            entry3 = locdb.data(loc[1])
                            entry3['{0}Gap'.format(qry)] += '<'
                            if not entry3['{0}5'.format(qry)]: entry3['{0}5'.format(qry)] = []
                            entry3['{0}5'.format(qry)].append(gap[0])
                        thisgap = gap.pop(0)
                        gentry = gapdb.data(thisgap)
                        if entry3 and thisgap[2] >= loc[1][1]: # Gap overlapping local alignment
                            gentry['GapSpan'] = 0
                            gentry['SynSpan'] = '{0}:0'.format(entry3[hit])
                            gentry['SynBad'] = 'Aligned'
                            continue
                        # At this point, loc[0][2] < thisgap[1] & thisgap[2] < loc[1][1]
                        if loc[1][1] == -1: gentry['GapSpan'] = gentry['seqlen'] - loc[0][2]
                        else: gentry['GapSpan'] = loc[1][1] - loc[0][2]
                        if not entry5 or not entry3:
                            gentry['SynBad'] = 'Terminal'
                            continue
                        # Assess synteny and classify
                        hitstart = '{0}Start'.format(hit)
                        hitend = '{0}End'.format(hit)
                        gap5 = {'+':entry5[hitend],'-':entry5[hitstart]}[entry5['Strand']]
                        gap3 = {'-':entry3[hitend],'+':entry3[hitstart]}[entry3['Strand']]
                        synspan = max(gap5,gap3) - min(gap5,gap3)
                        #i# Look for picking up synteny by skipping fragments from other sequences
                        #i# Where a local alignment switches the hit sequence, SynBad will look downstream for another
                        #i# aligment to the same sequence. A maximum of `maxsynskip=INT` alignments will be skipped, up to
                        #i# a maximum distance of `maxsynspan=INT` bp.
                        maxskip = self.getInt('MaxSynSkip')     # Maximum number of local alignments to skip for SynTrans classification
                        entryskip = None
                        egap = -1
                        if entry3[hit] != entry5[hit]:
                            for i in range(maxskip):
                                if len(loc) > (2+i) and loc[2+i][1] > 0:
                                    entryskip = locdb.data(loc[2+i])
                                    if not entryskip:
                                        self.warnLog('Problem finding local hit entry for {0}'.format(str(loc[i+2])))
                                        continue
                                    goodskip = True
                                    for field in [hit,'Strand']:
                                        if entryskip[field] != entry5[field]: goodskip = False
                                    if entry5['Strand'] == '+' and entryskip[hitstart] < entry5[hitend]: goodskip = False
                                    elif entry5['Strand'] == '-' and entryskip[hitend] > entry5[hitstart]: goodskip = False
                                    elif entry5['Strand'] == '+': egap = entryskip[hitstart] - entry5[hitend]
                                    elif entry5['Strand'] == '-': egap = entryskip[hitend] - entry5[hitstart]
                                    if egap > self.getInt('MaxSynSpan'): goodskip = False
                                    if loc[2+i][1] - loc[0][2] > self.getInt('MaxSynSpan'): goodskip = False
                                    if goodskip: break
                                    entryskip = None
                        # Gaps will then be classified according to the associated `GapSpan` and `SynSpan` values:
                        #
                        # * `Aligned` = Gap is found in the middle of a local alignment to the Hit
                        # * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
                        # * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
                        # * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
                        # * `Duplication` = Overlapping flanking hits on the same strand.
                        # * `Inversion` = Flanking hits are on alternative strands.
                        # * `Translocation` = `SynSpan` indicates matches are on different contigs.
                        # * `Terminal` = Gap is between a local alignment and the end of the query sequence.
                        # * `Null` = No mapping between genomes for that gap.
                        if entryskip:
                            gentry['SynSpan'] = '{0}:{1}:{2}'.format(entry5[hit],egap,entryskip[hit])
                            gentry['SynBad'] = 'Insertion'
                        elif entry3[hit] != entry5[hit]:
                            gentry['SynSpan'] = '{0}::{1}'.format(entry5[hit],entry3[hit])
                            if hitchr and (not entry5[hit].startswith(hitchr) or not entry3[hit].startswith(hitchr)):
                                gentry['SynBad'] = 'Fragmentation'
                            else:
                                gentry['SynBad'] = 'Translocation'
                        elif entry3['Strand'] != entry5['Strand']:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Inversion'
                        elif (entry3['Strand'] == '+' and gap3 < gap5) or (entry3['Strand'] == '-' and gap3 > gap5):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],-synspan)
                            gentry['SynBad'] = 'Duplication'
                        elif (synspan - gentry['GapSpan']) > self.getInt('MaxSynSpan'):
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Breakpoint'
                        else:
                            gentry['SynSpan'] = '{0}:{1}'.format(entry5[hit],synspan)
                            gentry['SynBad'] = 'Syntenic'

                # Then, deal with Hit Gaps
                for seqname in altdb.indexKeys('seqname'):
                    if seqname not in locdb.index(hit):
                        continue
                    self.progLog('\r#SYNBAD','SynBad mapping: {0}'.format(seqname))
                    gap = altdb.index('seqname')[seqname][0:]
                    loc = [(seqname,0,0)]
                    for entry in locdb.indexEntries(hit,seqname):
                        loc.append((entry[hit],entry['{0}Start'.format(hit)],entry['{0}End'.format(hit)],entry))
                    loc.sort()
                    loc.append((seqname,-1,-1))
                    while gap and len(loc) > 1 :
                        # Cycle to get the first gap in the list between two local alignment entries
                        while gap and gap[0][1] <= loc[0][2]: gap.pop(0)
                        while len(loc) > 1 and loc[1][2] < gap[0][1]: loc.pop(0)
                        if not gap or not len(loc) > 1: continue
                        # Should now have loc[0][2] < gap[0][1] & loc[1][2] > gap[0][1]
                        # Also want to have gap[0][2] < loc[1][1] for a gap properly flanked by two local alignments
                        thisgap = gap.pop(0)
                        if loc[0][1] > 0:
                            loc[0][3]['{0}Gap'.format(hit)] += '>'
                            if not loc[0][3]['{0}3'.format(hit)]: loc[0][3]['{0}3'.format(hit)] = []
                            loc[0][3]['{0}3'.format(hit)].append(thisgap)
                        if loc[1][1] > 0:
                            loc[1][3]['{0}Gap'.format(hit)] += '<'
                            if not loc[1][3]['{0}5'.format(hit)]: loc[1][3]['{0}5'.format(hit)] = []
                            loc[1][3]['{0}5'.format(hit)].append(thisgap)
                self.progLog('\r#SYNBAD','{0} SynBad mapping complete.               '.format(qh))
                self.printLog('\r#SYNBAD','{0} SynBad mapping complete.'.format(qh))
            ## ~ [3b] Update tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qrygapdb = self.db('qrygap')
            hitgapdb = self.db('hitgap')
            qrydb = self.db('qry')
            hitdb = self.db('hit')
            ex = 0.0; etot = qrydb.entryNum() + hitdb.entryNum()
            for entry in qrydb.entries() + hitdb.entries():
                self.progLog('\r#UPDATE','Updating flanking gap SynBad ratings: %.2f%%' % (ex/etot)); ex += 100
                for field in ('Qry5','Qry3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(qrygapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                for field in ('Hit5','Hit3'):
                    if entry[field]:
                        gaps = []
                        for gap in entry[field]: gaps.append(hitgapdb.data(gap)['SynBad'])
                        entry[field] = ';'.join(rje.sortUnique(gaps))
                    else: entry[field] = '-'
                if not entry['QryGap']: entry['QryGap'] = '.'
                if not entry['HitGap']: entry['HitGap'] = '.'
            self.printLog('\r#UPDATE','Updating flanking gap SynBad ratings complete.')

            ## ~ [3c] Compress tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                locdb = self.db(qh.lower())
                prev = None
                ex = 0.0; etot = locdb.entryNum()
                for ekey in locdb.datakeys()[0:]:
                    self.progLog('\r#CMPRSS','Compress %s table: %.2f%%' % (qh,ex/etot)); ex += 100
                    #self.bugPrint(prev); self.bugPrint(ekey)
                    entry = locdb.data(ekey)
                    #i# Will compress adjacent entries with the same sequence pairs and no flanking gaps in Qry or Hit
                    if entry['QryGap'] != '.' or entry['HitGap'] != '.': prev = None; continue
                    if not prev: prev = ekey; continue
                    pentry = locdb.data(prev)
                    #self.bugPrint(pentry); self.bugPrint(entry)
                    matched = True
                    for field in ['Qry','Hit','Strand']:
                        matched = matched and entry[field] == pentry[field]
                    if not matched: prev = ekey; continue
                    # except:
                    #     self.deBug(prev)
                    #     self.deBug('?: %s' % pentry)
                    #i# Check Hit is compatible
                    if entry['Strand'] == '+' and entry['HitStart'] < pentry['HitEnd']: prev = ekey; continue
                    if entry['Strand'] == '-' and entry['HitEnd'] > pentry['HitStart']: prev = ekey; continue
                    #i# Combine and compress entries
                    for field in ['Length','Identity']:
                        pentry[field] += entry[field]
                    if qh == 'Qry':
                        pentry['Unaligned'] += (entry['QryStart'] - pentry['QryEnd'])
                    elif pentry['Strand'] == '+':
                        pentry['Unaligned'] += (entry['HitStart'] - pentry['HitEnd'])
                    else:
                        pentry['Unaligned'] += (pentry['HitStart'] - entry['HitEnd'])
                    pentry['QryEnd'] = entry['QryEnd']
                    if pentry['Strand'] == '+': pentry['HitEnd'] = entry['HitEnd']
                    else: pentry['HitStart'] = entry['HitStart']
                    # self.deBug(pentry)
                    # self.bugPrint('Popping...')
                    # self.bugPrint(ekey)
                    # self.bugPrint(prev)
                    locdb.dict['Data'].pop(ekey)
                    locdb.dict['Data'].pop(prev)
                    prev = locdb.makeKey(pentry)
                    locdb.dict['Data'][prev] = pentry
                self.printLog('\r#CMPRSS','Compressed %s table: %s alignments -> %s' % (qh,rje.iStr(etot),rje.iStr(locdb.entryNum())))

            ## ~ [3d] Save tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for qh in ('Qry','Hit'):
                self.db('{0}gap'.format(qh.lower())).saveToFile()
                self.db(qh.lower()).dropField('AlnNum')
                self.db(qh.lower()).saveToFile()

            ## ~ [3e] Top-matched pairs only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pairs = []
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'top{0}'.format(qh.lower()),replace=True,add=True)
                topdb.keepFields(['Qry','Hit','Length']+topdb.keys())
                topdb.compress(['Qry','Hit'],default='sum')
                topdb.keepFields(['Qry','Hit','Length'])
                topdb.rankFieldByIndex(qh,'Length',newfield='Rank',rev=True,absolute=True,lowest=True,unique=False,warn=True,highest=False)
                topdb.dropEntriesDirect('Rank',[1],inverse=True,log=True,force=False)
                pairs += topdb.dataKeys()
            for qh in ('Qry','Hit'):
                topdb = db.copyTable(self.db(qh.lower()),'{0}pairs'.format(qh.lower()),replace=True,add=True)
                ex = 0.0; etot = topdb.entryNum()
                for ekey in topdb.datakeys()[0:]:
                    self.progLog('\r#PAIRS','Reducing %s table to top-aligned pairs: %.2f%%' % (qh,ex/etot)); ex += 100
                    entry = topdb.data(ekey)
                    if (entry['Qry'],entry['Hit']) not in pairs: topdb.dict['Data'].pop(ekey)
                self.printLog('\r#PAIRS','Reduced %s table to %s alignments between top sequence pairs.' % (qh,rje.iStr(topdb.entryNum())))
                topdb.saveToFile()

            ## ~ [3f] Summary reports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.headLog(base1,line='-')
            self.db('qrygap').indexReport('SynBad')
            self.headLog(base2,line='-')
            self.db('hitgap').indexReport('SynBad')

            ### ~ [4] Fragment Gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('GAP FRAGMENTATION',line='=')
            # If `fragment=T`, the assemblies will then be fragmented on gaps that are not Syntenic, unless more than
            # `minreadspan=INT` reads span the gap.
            if self.getBool('Fragment'):
                self.fragment()
            else:
                self.printLog('#FRAG','No fragmentation (fragment=F).')

            # A future release of Synbad will optionally re-arrange the two assemblies, incorporating gapass assemblies
            # where possible.

            return
        except: self.errorLog('%s.synBad error' % self.prog())
#########################################################################################################################
    def seqObjSetup(self):   ### Loads the two genomes into sequence list objects
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Genome1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList1' not in self.obj or not self.obj['SeqList1']:
                base1 = rje.baseFile(self.getStr('Genome1'),strip_path=True)
                cmd1 = ['seqin={0}'.format(self.getStr('Genome1')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base1)]
                seqcmd = self.cmd_list + cmd1 + ['summarise=F','gapstats=F','raw=F','dna=T']
                self.obj['SeqList1'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
            ## ~ [1b] Genome2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'SeqList2' not in self.obj or not self.obj['SeqList2']:
                base2 = rje.baseFile(self.getStr('Genome2'),strip_path=True)
                cmd2 = ['seqin={0}'.format(self.getStr('Genome2')),
                        'runmode={0}'.format(self.getStrLC('GapMode')),
                        'basefile={0}'.format(base2)]
                seqcmd = self.cmd_list + cmd2 + ['summarise=F','gapstats=F','raw=F','dna=T']
                self.obj['SeqList2'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
            return
        except: self.errorLog('%s.seqObjSetup error' % self.prog())
#########################################################################################################################
    def fragment(self):   ### SynBad fragmentation
        '''
        SynBad is a tool for comparing two related genome assemblies and identify putative translocations and inversions
        between the two that correspond to gap positions. These positions could indicate misplaced scaffolding.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seqObjSetup()
            seqobj = {}
            seqobj['qryfrag'] = self.obj['SeqList1']
            seqobj['hitfrag'] = self.obj['SeqList2']
            base1 = rje.baseFile(self.getStr('Genome1'),strip_path=True)
            base2 = rje.baseFile(self.getStr('Genome2'),strip_path=True)
            db = self.db()
            basefile = self.baseFile()
            qrygap = self.db('qrygap')
            if not qrygap:
                qrygap = db.addTable('%s.gaps.tdt' % base1,mainkeys=['seqname','start','end'],name='qrygap',ignore=[],expect=True)
            qrygap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            qryfrag = db.copyTable(qrygap,'qryfrag',replace=True,add=True)
            hitgap = self.db('hitgap')
            if not hitgap:
                hitgap = db.addTable('%s.gaps.tdt' % base2,mainkeys=['seqname','start','end'],name='hitgap',ignore=[],expect=True)
            hitgap.dataFormat({'seqlen':'int','start':'int','end':'int','Span0':'int'})
            hitfrag = db.copyTable(hitgap,'hitfrag',replace=True,add=True)
            ### ~ [1] Filter gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # * `Aligned` = Gap is found in the middle of a local alignment to the Hit
            # * `Syntenic` = Difference between positive `SynSpan` and `GapSpan` is `maxsynspan=INT` or less (default 10kb).
            # * `Insertion` = Achieved `Syntenic` rating by skipping upto `maxsynskip=INT` local alignments and max `maxsynspan=INT` bp in both Qry and Hit.
            # * `Breakpoint` = Difference between positive `SynSpan` and `GapSpan` is bigger than the `maxsynspan=INT` distance.
            # * `Inversion` = Negative `SynSpan` value.
            # * `Translocation` = `SynSpan` indicates matches are on different contigs.
            # * `Terminal` = Gap is between a local alignment and the end of the query sequence.
            # * `Null` = No mapping between genomes for that gap.
            self.headLog('SYNBAD GAP FILTER',line='=')
            for table in (qryfrag,hitfrag):
                table.addField('fragid',evalue='')
                for entry in table.entries():
                    if entry['Span0'] < self.getInt('MinReadSpan') and entry['SynBad'] in ['Breakpoint','Inversion','Translocation']:
                        entry['fragid'] = 'Filter'
                table.dropEntriesDirect('fragid',['Filter'],inverse=True)
                table.newKey(['seqname','start','end','fragid'])
                table.addField('syn5',evalue='')
                table.addField('syn3',evalue='')
            ### ~ [2] Fragment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('FRAGMENT ON GAPS',line='=')
            for frag in ('qryfrag','hitfrag'):
                table = self.db(frag)
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                seqlist = seqobj[frag]
                prev = {'seqname':None}
                fragname = ''; fx = 1
                for ekey in table.dataKeys():
                    entry = table.data(ekey)
                    #i# Add extra 5' fragment for each sequence
                    if entry['seqname'] != prev['seqname']:
                        prev = rje.combineDict({},entry)
                        prev['start'] = 1
                        prev['end'] = entry['start'] - 1
                        fragname = entry['seqname']; fx = 1
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        prev['fragid'] = '{0}.{1}'.format(fragname,fx); fx += 1
                        prev['syn5'] = 'Terminal'
                        prev['syn3'] = 'Terminal'    #i# May be over-written
                        prev = table.addEntry(prev,remake=False)
                    #i# Update details of this entry
                    prev['syn3'] = entry['SynBad']
                    prev['end'] = entry['start'] - 1
                    entry['start'] = entry['end'] + 1
                    entry['fragid'] = '{0}.{1}'.format(fragname,fx); fx += 1
                    entry['syn5'] = entry['SynBad']
                    entry['syn3'] = 'Terminal'                          #i# May be over-written
                    entry['end'] = entry['seqlen']                      #i# May be overwritten
                    prev = entry
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    if entry['seqname'] == 'NSCUSCAFF155': self.debug('?')
                #i# Add sequences without any gaps
                table.remakeKeys()
                seqnames = rje.sortKeys(table.index('seqname'))
                self.printLog('#FRAG','{0} fragments from {1} sequences with gaps'.format(rje.iStr(table.entryNum()),rje.iLen(seqnames)))
                sx = 0
                for seq in seqlist.seqs():
                    seqname = seqlist.shortName(seq)
                    if seqname not in seqnames:
                        fragname = seqname
                        if rje.matchExp('\S+_\S+__(\S+)',fragname): fragname = rje.matchExp('\S+_\S+__(\S+)',fragname)[0]
                        fragid = '{0}.1'.format(fragname)
                        table.addEntry({'fragid':fragid,'seqname':seqname,'start':1,'seqlen':seqlist.seqLen(seq),'end':seqlist.seqLen(seq),'syn5':'Terminal','syn3':'Terminal'})
                        sx += 1
                self.printLog('#SEQ','Added {0} full-length sequences without gaps'.format(rje.iStr(sx)))
                #i# Tidy up fragment table
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
                #table.newKey(['fragid'])
                table.setFields(['fragid','seqname','start','end','seqlen','syn5','syn3'])
                table.saveToFile()
                #self.debug('{0} {1}: {2}'.format(table,table.name(),table.fields()))
            ### ~ [3] Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                seqdict = seqlist.seqNameDic()
                fragfas = '{0}.{1}.fasta'.format(basefile,frag)
                rje.backup(self,fragfas)
                FRAGFAS = open(fragfas,'w'); ex = 0
                for entry in table.entrySort():
                    self.progLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(ex),fragfas)); ex += 1
                    #self.bugPrint(table.entrySummary(entry,collapse=True))
                    seq = seqdict[entry['seqname']]
                    (seqname,sequence) = seqlist.getSeq(seq)
                    FRAGFAS.write('>{0}\n'.format(entry['fragid']))
                    FRAGFAS.write(sequence[entry['start']-1:entry['end']]+'\n')
                FRAGFAS.close()
                self.printLog('\r#FASOUT','{0} fragments output to {1}'.format(rje.iStr(table.entryNum()),fragfas))

        except:
            self.errorLog('%s.fragment error' % self.prog())
            raise   # Delete this if method error not terrible
#########################################################################################################################
### End of SECTION II: SynBad Class                                                                                     #
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
    try: SynBad(mainlog,['basefile=synbad']+cmd_list).run()

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
