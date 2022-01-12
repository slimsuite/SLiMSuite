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
Module:       rje_readcore
Description:  Read mapping and analysis core module
Version:      0.7.1
Last Edit:    10/01/22
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module has very simple standalone functionality to check for the existence of a BAM file and generate it from
    one or more long read input files if not found. It is primarily for inheritance by other RJE tools that need to
    make use of long read mapping and/or simple depth/coverage stat wrapping. The core ReadCore object will have
    methods for populating and checking key input files and settings for a number of other SeqSuite tools.

    ## Dependencies

    For read mapping, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to the environment
    `$PATH` or given with the `minimap2=PROG` setting. For depth summaries, [samtools](http://www.htslib.org/) needs to
    be installed. Presence of dependencies will be checked when needed and an error raised if not found.

Commandline:
    ### ~ Core input options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [$SEQIN basefile]
    paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
    bam=FILE        : BAM file of long reads mapped onto assembly [$BASEFILE.bam]
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    ### ~ Depth and Copy Number options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    scdepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]
    busco=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    quickdepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    depfile=FILE    : Precomputed depth file (*.fastdep or *.fastmp) to use [None]
    regfile=FILE    : File of SeqName, Start, End positions (or GFF) for read coverage checking [None]
    checkfields=LIST: Fields in checkpos file to give Locus, Start and End for checking [SeqName,Start,End]
    gfftype=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    depadjust=INT   : Advanced R density bandwidth adjustment parameter [12]
    seqstats=T/F    : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
    cnmax=INT       : Max. y-axis value for CN plot (and mode multiplier for related depth plots) [4]
    ### ~ System options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    tmpdir=PATH     : Path for temporary output files during forking [./tmpdir/]
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    rscript=PROG    : Full path to run minimap2 [Rscript]
    samtools=PROG   : Full path to run minimap2 [samtools]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys, time
mypath = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + os.path.sep
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_forker, rje_obj, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Adding forking for fastdepth file generation.
    # 0.2.0 - Added CovBases lower depthsizer estimate output based solely on mapped reads.
    # 0.2.1 - Added unique sorting of CIGAR strings for indel ratio. Fixed end padding of zero-coverage depths.
    # 0.2.2 - Fixed major flaw in indelratio calculation.
    # 0.3.0 - Add benchmark=T/F option to the genome size prediction. Tidied CovBase and MapAdjust.
    # 0.3.1 - Tweaked some input checks and log output. Replaced indelratio sort -u with uniq for speed and memory.
    # 0.4.0 - Added seqstats=T/F : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]
    # 0.4.1 - Fixed bug that causes clashes with v5 full_table.bed files.
    # 0.5.0 - Add additional map adjustment variants:
    #       - MapAdjust2 = allbases, not covbases
    #       - MapBases = Use map bases, not covbases for min read volumne
    #       - MapRatio = Use mapbases adjusted by indelratio
    # 0.6.0 - Added support for multiple regfiles and setting max limit for CN graphics.
    # 0.7.0 - Added passing on of gfftype=LIST option to Rscript.
    # 0.7.1 - Fixed readtype recycle bug.
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
    # [ ] : Add module load if cannot find program.
    # [ ] : Add NGMLR to mappers.
    # [Y] : Try using total sequence length not covbases (samtools coverage $3 not $5) for CovBases calculation.
    # [ ] : Add bamcsi=T/F : Use CSI indexing for BAM files, not BAI (needed for v long scaffolds) [False]
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('ReadMap', '0.7.1', 'January 2022', '2021')
    description = 'Read mapping analysis module'
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
paf_defaults = {'N':'250','p':'0.0001','x':'asm20'}
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: ReadCore Class                                                                                           #
#########################################################################################################################
class ReadCore(rje_obj.RJE_Object):
    '''
    ReadCore Class. Author: Rich Edwards (2021).

    Str:str
    - BAM=FILE        : BAM file of reads mapped onto assembly [$BASEFILE.bam]
    - BUSCO=TSVFILE   : BUSCO full table [full_table_$BASEFILE.busco.tsv]
    - DepFile=FILE    : Precomputed depth file (*.fastdep or *.fastmp) to use [None]
    - PAF=FILE        : PAF file of reads mapped onto assembly [$BASEFILE.paf]
    - RegFile=FILE    : File of SeqName, Start, End positions (or GFF) for read coverage checking [None]
    - SeqIn=FILE      : Input sequence assembly (sortnr/diphap modes) []
    - TmpDir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]

    Bool:boolean
    - Minimap2        : Whether Minimap2 found on system
    - QuickDepth=T/F  : Whether to use samtools depth in place of mpileup (quicker but underestimates?) [False]
    - Rscript         : Whether Rscript found on system
    - Samtools        : Whether Samtools found on system
    - SeqStats=T/F    : Whether to output CN and depth data for full sequences as well as BUSCO genes [False]

    Int:integer
    - Adjust=INT   : Advanced R density bandwidth adjustment parameter [12]
    - CNMax=INT       : Max. y-axis value for CN plot (and mode multiplier for related depth plots) [4]

    Num:float
    - SCDepth=NUM     : Single copy ("diploid") read depth. If zero, will use SC BUSCO mode [0]

    File:file handles with matching str filenames
    
    List:list
    - CheckFields=LIST: Fields in checkpos file to give Locus, Start and End for checking [SeqName,Start,End]
    - GFFType=LIST    : Optional feature types to use if performing regcheck on GFF file (e.g. gene) ['gene']
    - Reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    - ReadType=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database object
    - Forker = Forking controller
    - SeqIn = rje_seqlist.SeqList object (genome assembly)
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BAM','BUSCO','DepFile','PAF','RegFile','SeqIn','TmpDir']
        self.boollist = ['QuickDepth','SeqStats']
        self.intlist = ['Adjust','CNMax']
        self.numlist = ['SCDepth']
        self.filelist = []
        self.listlist = ['CheckFields','GFFType','Reads','ReadType']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self._setReadCoreAttributes()
        self.setStr({})
        self.setBool({})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _setReadCoreAttributes(self):   ### Sets defaults for all ReadCore attributes
        '''
        Sets defaults for all ReadCore attributes.
        '''
        self.setStr({'BAM':'None','BUSCO':'None','DepFile':'None','PAF':'None','RegFile':'None','SeqIn':'None','TmpDir':'./tmpdir/',
                     'Minimap2':'minimap2','Samtools':'samtools'})
        self.setBool({'QuickDepth':False,'SeqStats':False})
        self.setInt({'Adjust':12,'CNMax':4})
        self.setNum({'AllBases':0,'CovBases':0,'MapAjust':0,'MapBases':0,'OldAdjust':0,'SCDepth':0})
        self.list['CheckFields'] = ['SeqName','Start','End']
        self.list['GFFType'] = ['gene']
        self.list['Reads'] = []
        self.list['ReadType'] = ['ont']
        #i# Check for core programs
        self.setBool({'Minimap2':True,'Samtools':True,'Rscript':True})
        #i# Set up Database object
        self.obj['DB'] = rje_db.Database(self.log, self.cmd_list + ['tuplekeys=T'])
        self.obj['SeqIn'] = None
        self.obj['Forker']  = rje_forker.Forker(self.log,['logfork=F','killmain=F']+self.cmd_list)
#########################################################################################################################
    def _readCoreCmd(self,cmd):     ### Sets Core Attributes from commandline
        '''
        Sets Core Attributes from commandline.
        '''
        ### Class Options (No need for arg if arg = att.lower()) ###
        self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
        self._cmdReadList(cmd,'str',['RegFile'])  # String representing directory path
        self._cmdReadList(cmd,'file',['BAM','BUSCO','DepFile','PAF','SeqIn'])  # String representing file path
        self._cmdReadList(cmd,'bool',['QuickDepth','SeqStats','Minimap2','Samtools','Rscript'])
        self._cmdReadList(cmd,'int',['Adjust','CNMax'])
        self._cmdReadList(cmd,'num',['SCDepth'])
        self._cmdReadList(cmd,'glist',['Reads'])
        self._cmdReadList(cmd,'list',['CheckFields','GFFType','ReadType'])
        self._cmdRead(cmd,'int','Adjust','depadjust')   # Integers
        self._cmdRead(cmd,type='str',att='RegFile',arg='regcheck')  # No need for arg if arg = att.lower()
        self._cmdRead(cmd,type='list',att='CheckFields',arg='reghead')  # No need for arg if arg = att.lower()
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)      # Delete if no forking
                self._readCoreCmd(cmd)     # Will set all the core commands recognised.
                ### Class Options (No need for arg if arg = att.lower()) ### 
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['TmpDir'])  # String representing directory path
                #self._cmdReadList(cmd,'file',['BAM','BUSCO','PAF','RegFile','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['QuickDepth'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['DepAdjust'])   # Integers
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
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# If samtools is found, will make a BAM file. Else make a PAF file. Maybe add toggle at some point?
            if self.getBool('Samtools'):
                self.getBamFile()
            else:
                self.getPAFFile()
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def getBamFile(self): return self.getBAMFile()
    def getBAMFile(self,make=True,expect=True):  ### Checks/Creates indexed BAM file and returns filename as string
        '''
        Checks/Creates indexed BAM file and returns filename as string.
        >> make:bool [True] = Whether to make if missing.
        >> expect:bool [True] = Whether to raise error if missing.
        :return: bamfile [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Check for existing indexed BAM file and index if found but not indexed.
            if not self.getBool('Samtools'):
                raise RuntimeError('Cannot find samtools!')
            bamfile = self.setBAMFile()
            if self.checkBAMFile(bamfile): return bamfile
            if not make:
                if expect and rje.exists(bamfile): raise IOError('Problem with BAM File: {0}!'.format(bamfile))
                elif expect: raise IOError('{0} not found!'.format(bamfile))
                else: return None
            ### ~ [2] Generate BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Minimap2'):
                raise RuntimeError('Cannot find minimap2!')
            self.printLog('#BAM','Generating BAM file: {0}'.format(bamfile))
            return self.longreadMinimap()   #i# Included BAM file and index checks
        except:
            self.errorLog('{0}.getBamFile() error'.format(self.prog()))
            if expect: raise
        return None
#########################################################################################################################
    def getPAFFile(self,make=True,expect=True):  ### Checks/creates PAF file and returns filename as string
        '''
        Checks for PAF file and returns filename as string.
        >> make:bool [True] = Whether to make if missing.
        >> expect:bool [True] = Whether to raise error if missing.
        :return: paffile [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Check for existing PAF file
            seqin = self.getStr('SeqIn')
            paffile = self.setPAFFile()
            if not self.needToRemake(paffile,seqin): return paffile
            if not make:
                if expect: raise IOError('{0} not found or too old!'.format(paffile))
                else: return None
            ### ~ [2] Generate PAF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Minimap2'):
                raise RuntimeError('Cannot find minimap2!')
            self.printLog('#PAF','Generating PAF file: {0}'.format(paffile))
            return self.longreadMinimap(paf=True)
        except:
            self.errorLog('{0}.getPAFFile() error'.format(self.prog())); return None
#########################################################################################################################
    def getFastDep(self):   ### Checks for fastdep file and generates if missing
        '''
        Checks for fastdep file and generates if missing
        :return: fastdep file [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            depfile = self.getStr('DepFile')
            if not self.force() and rje.exists(depfile):
                self.printLog('#DEP', 'Existing {0} file found (force=F)'.format(depfile))
                return depfile
            if not self.getBool('Samtools'):
                raise RuntimeError('Cannot find samtools!')
            seqin = self.seqinObj()
            bamfile = self.getBAMFile()
            self.printLog('#BAM','BAM File for read depth extraction: {0}'.format(bamfile))
            fastdep = '{0}.fastmp'.format(bamfile)
            if self.getBool('QuickDepth'):
                fastdep = '{0}.fastdep'.format(bamfile)
            if not self.force() and not self.needToRemake(fastdep,bamfile):
                self.printLog('#DEP','Existing {0} file found (force=F)'.format(fastdep))
                self.setStr({'DepFile': fastdep})
                return fastdep
            ### ~ [2] ~ Generate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Generating depth file', line='-')
            #!# Add forking via temp files.
            rje.backup(self,fastdep)
            if self.threads() > 1: return self.forkDepth(bamfile,fastdep)
            self.printLog('#FORK', 'Note: program can be accelerated using forks=INT.')
            OUT = open(fastdep,'w')
            ox = 0.0; otot = seqin.seqNum()
            for seq in seqin.seqs():
                sname = seqin.shortName(seq)
                seqlen = seqin.seqLen(seq)
                self.progLog('\r#DEP','Generating depth file "{0}": {1:.1f}%'.format(fastdep,ox/otot)); ox += 100.0
                depcmd = 'samtools view -b -h -F 0x100 {0} {1} | '.format(bamfile,sname)
                if self.getBool('QuickDepth'):
                    depcmd += 'samtools depth -a - | '
                else:
                    depcmd += 'samtools mpileup -BQ0 - 2> /dev/null | '
                depcmd += 'awk \'BEGIN { prev_chr="";prev_pos=0;} { if($1==prev_chr && prev_pos+1!=int($2)) {for(i=prev_pos+1;i<int($2);++i) {printf("%s\\t%d\\t0\\n",$1,i);}} print; prev_chr=$1;prev_pos=int($2);}\' | '
                if self.getBool('QuickDepth'):
                    depcmd += "awk '{print $3;}'"
                else:
                    depcmd += "awk '{print $4;}'"
                # i# Pad the end if required
                deplist = os.popen(depcmd).read().split()
                if len(deplist) < seqlen:
                    deplist += ['0'] * (seqlen - len(deplist))
                OUT.write('>{0}\n'.format(sname) + ' '.join(deplist) + '\n')
            self.printLog('\r#DEP', 'Generated depth file "{0}": {1} sequences'.format(fastdep, otot))
            OUT.close()
            self.setStr({'DepFile':fastdep})
            return fastdep
        except:
            self.errorLog('{0}.getFastDep() error'.format(self.prog())); return None
#########################################################################################################################
    def forkDepth(self,bamfile,fastdep,secondary=False,setstr=True):    ### Generates the fast depth file using forks and temp directory.
        '''
        Generates the fast depth file using forks and temp directory.
        >> bamfile:str = Source BAM file for depth parsing
        >> fastdep:str = output depth file.
        >> secondary:bool [False] = whether to include secondary alignments
        >> setstr:bool [True] = Whether to set self.str['DepFile']
        :return:
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.force() and rje.exists(fastdep):
                self.printLog('#DEP', 'Existing {0} file found (force=F)'.format(fastdep))
                return fastdep
            if not self.getBool('Samtools'):
                raise RuntimeError('Cannot find samtools!')
            forker = self.obj['Forker']
            basefile = self.baseFile(strip_path=True)
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            seqin = self.seqinObj() #rje_seqlist.SeqList(self.log,['summarise=T']+self.cmd_list+['autoload=T','seqmode=file'])
            ## ~ [1a] Temp directory for forked depths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tmpdir = rje.makePath(self.getStr('TmpDir'),wholepath=False)
            if not rje.exists(tmpdir): rje.mkDir(self,tmpdir)

            ### ~ [2] Cycle through and fork out the depth calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Setup forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cleanup = 0; skipped = 0
            forker.list['ToFork'] = []
            for sname in seqin.names():
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,sname,depmethod)
                if rje.exists(tmpfile):
                    #?# Add checking for completeness #?#
                    if not self.force(): skipped += 1; continue
                    else: os.unlink(tmpfile); cleanup += 1
                #i# Depth command to run
                depcmd = 'samtools view -b -h -F 0x100 {0} {1} | '.format(bamfile,sname)
                if secondary: depcmd = 'samtools view -b -h {0} {1} | '.format(bamfile,sname)
                if self.getBool('QuickDepth'):
                    depcmd += 'samtools depth -a - | '
                else:
                    depcmd += 'samtools mpileup -BQ0 - 2> /dev/null | '
                depcmd += 'awk \'BEGIN { prev_chr="";prev_pos=0;} { if($1==prev_chr && prev_pos+1!=int($2)) {for(i=prev_pos+1;i<int($2);++i) {printf("%s\\t%d\\t0\\n",$1,i);}} print; prev_chr=$1;prev_pos=int($2);}\' | '
                if self.getBool('QuickDepth'):
                    depcmd += "awk '{print $3;}'"
                else:
                    depcmd += "awk '{print $4;}'"
                depcmd += ' > {0}'.format(tmpfile)
                #i# Forks
                forker.list['ToFork'].append(depcmd)
            self.printLog('#DEPTH','{} sequences queued for forking ({} existing files deleted); {} existing results skipped'.format(rje.iLen(forker.list['ToFork']),rje.iStr(cleanup),rje.iStr(skipped)))
            ## ~ [2b] Fork out depth analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if forker.list['ToFork']:
                if self.getNum('Forks') < 2:
                    #i# Warn lack of forking
                    self.printLog('#FORK','Note: program can be accelerated using forks=INT.')
                    for forkcmd in forker.list['ToFork']:
                        self.printLog('#SYS',forkcmd)
                        os.system(forkcmd)
                elif forker.run():
                    self.printLog('#FORK','Forking of depth parsing completed.')
                else:
                    try:
                        self.errorLog('Depth forking did not complete',printerror=False,quitchoice=True)
                    except:
                        raise RuntimeError('Depth forking did not complete')

            ### ~ [3] Convert to fastdepth format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            OUT = open(fastdep, 'w')
            ox = 0.0
            otot = seqin.seqNum()
            for seq in seqin.seqs():
                sname = seqin.shortName(seq)
                seqlen = seqin.seqLen(seq)
                self.progLog('\r#DEP', 'Generating depth file "{0}": {1:.1f}%'.format(fastdep, ox / otot))
                ox += 100.0
                tmpfile = '{}{}.{}.{}.tmp'.format(tmpdir,basefile,sname,depmethod)
                #i# Pad the end if required
                deplist = open(tmpfile,'r').read().split()
                if len(deplist) < seqlen:
                    deplist += ['0'] * (seqlen - len(deplist))
                OUT.write('>{0}\n'.format(sname) + ' '.join(deplist) + '\n')
            self.printLog('\r#DEP', 'Generated depth file "{0}": {1} sequences'.format(fastdep, otot))
            OUT.close()
            if setstr:
                self.setStr({'DepFile':fastdep})
            return fastdep

        except:
            self.errorLog('{0}.forkDepth() error'.format(self.prog()))
            raise
#########################################################################################################################
    ### <3> ### Setup Methods                                                                                           #
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setBool({'Minimap2': self.getBool('Minimap2') and self.checkForProgram('minimap2'),
                          'Samtools': self.getBool('Samtools') and self.checkForProgram('samtools'),
                          'Rscript': self.getBool('Rscript') and self.checkForProgram('Rscript')})
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def seqinObj(self,summarise=False): ### Returns the a SeqList object for the SeqIn file
        '''
        Returns the a SeqList object for the SeqIn file.
        :return: self.obj['SeqIn']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['SeqIn']:
                seqcmd = self.cmd_list
                if summarise: seqcmd = ['summarise=T']+self.cmd_list
                self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
                sx = 0.0; stot = self.obj['SeqIn'].seqNum()
                for seq in self.obj['SeqIn'].seqs():
                    self.progLog('\r#CHECK','Checking sequences names: %.1f%%' % (sx/stot)); sx += 100.0
                    if '|' in self.obj['SeqIn'].shortName(seq):
                        raise ValueError('Pipe "|" characters found in seqin=FILE names: will break program. Please rename and try again.')
                self.printLog('\r#CHECK','Checking sequences names complete.')
        except ValueError:
            self.printLog('\r#CHECK','Checking sequences names aborted.')
            self.errorLog('{0} input sequence error'.format(self.prog())); raise
        except:
            self.errorLog('{0} seqinObj() error'.format(self.prog()))
        return self.obj['SeqIn']
#########################################################################################################################
    def checkInput(self,busco=False,regfile=False,reads=False):   ### Checks whether possible input files exist and are in correct format
        '''
        Checks whether possible input files exist and are in correct format.
        >> busco:bool [False] = Whether BUSCO file needed
        >> regfile:bool [False] = Whether RegFile file needed
        >> reads:bool [False] = Whether reads are required
        << returns True or False
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            #i# Files to check are:
            # - BUSCO
            # - RegFile
            ### ~ [2] Check BUSCO File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('BUSCO'):
                rje.checkForFiles([self.getStr('BUSCO')],basename='',log=self.log,ioerror=True)
                #!# Add format checks at some point
            elif busco:
                raise IOError('BUSCO file not given (busco=FILE)')
            ### ~ [3] Check RegFile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('RegFile'):
                regfilename = string.split(string.split(self.getStr('RegFile'),',')[0],':')[-1]
                rje.checkForFiles([regfilename],basename='',log=self.log,ioerror=True)
                ## ~ [3a] Check Fields ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if regfilename.endswith('.gff') or regfilename.endswith('.gff3'):
                    self.printLog('#GFF','GFF file recognised for RegFile')
                else:
                    #!# Add feature to recognise and change first field if not given #!#
                    if not len(self.list['CheckFields']) == 3:
                        raise ValueError('checkfields=LIST must have exactly 3 elements: SeqName, Start, End. %d found!' % len(self.list['CheckFields']))
                    [locusfield,startfield,endfield] = self.list['CheckFields']
                    cdb = db.addTable(regfilename,mainkeys='auto',name='check',expect=True)
                    if not cdb: raise IOError('Cannot find checkpos file "%s"' % self.getStr('CheckPos'))
                    if locusfield not in cdb.fields():
                        self.warnLog('Field "%s" not found in checkpos file: will use "%s for sequence name' % (locusfield,cdb.fields()[0]))
                        locusfield = cdb.fields()[0]
                    if startfield not in cdb.fields():
                        newstart = cdb.fields()[1]
                        if 'Start' in cdb.fields(): newstart = 'Start'
                        self.warnLog('Field "%s" not found in checkpos file: will use "%s for start position' % (startfield,newstart))
                        startfield = newstart
                    if endfield not in cdb.fields():
                        newfield = cdb.fields()[2]
                        if 'End' in cdb.fields(): newfield = 'End'
                        self.warnLog('Field "%s" not found in checkpos file: will use "%s for end position' % (endfield,newfield))
                        endfield = newfield
                    self.list['CheckFields'] = [locusfield,startfield,endfield]
            elif regfile:
                raise IOError('Region file not given (regfile=FILE)')
            ### ~ [4] Check Reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.checkReadList(expect=reads)
            return True
        except: self.errorLog('%s.checkInput error' % self.prog()); return False
#########################################################################################################################
    def minimap2(self):
        try: return self.getStr('Minimap2')
        except: return 'minimap2'
#########################################################################################################################
    def samtools(self):
        try: return self.getStr('Samtools')
        except: return 'samtools'
#########################################################################################################################
    ### <4> ### Read Mapping Methods                                                                                    #
#########################################################################################################################
    def setBAMFile(self,baseoveride=True):   ### Sets BAM File name from settings
        '''
        Sets BAM File name from settings. If bam=FILE is not set, will name after seqin=FILE. If seqin=FILE is not set,
        will name after basefile=FILE.
        >> baseoveride:bool [True] = Whether to over-ride SeqIn basefile if $BASE.bam is found.
        '''
        if not self.getStrLC('BAM'):
            basebam = self.baseFile(strip_path=True) + '.bam'
            if baseoveride and rje.exists(basebam):
                return basebam
            if self.getStrLC('SeqIn'):
                self.setStr({'BAM': rje.baseFile(self.getStr('SeqIn'),strip_path=True) + '.bam'})
            else:
                self.setStr({'BAM':basebam})
            self.printLog('#BAM','Set BAM file: {0}'.format(self.getStr('BAM')))
        return self.getStr('BAM')
#########################################################################################################################
    def setPAFFile(self,baseoveride=True):   ### Sets PAF File name from settings
        '''
        Sets PAF File name from settings. If paf=FILE is not set, will name after seqin=FILE. If seqin=FILE is not set,
        will name after basefile=FILE.
        >> baseoveride:bool [True] = Whether to over-ride SeqIn basefile if $BASE.paf is found.
        '''
        if not self.getStrLC('PAF'):
            basepaf = self.baseFile(strip_path=True) + '.paf'
            if baseoveride and rje.exists(basepaf):
                return basepaf
            if self.getStrLC('SeqIn'):
                self.setStr({'PAF': rje.baseFile(self.getStr('SeqIn'),strip_path=True) + '.paf'})
            else:
                self.setStr({'PAF':basepaf})
            self.printLog('#PAF','Set PAF file: {0}'.format(self.getStr('PAF')))
        return self.getStr('PAF')
#########################################################################################################################
    def checkBAMFile(self,bamfile,makeindex=True,bai=False,csi=False,needed=False): ### Checks for indexed BAM file
        '''
        Checks for indexed BAM file.
        :param bamfile: str = BAM file to check
        :param makeindex: bool [True] = Whether to check for index and make if needed
        :param bai: bool [False] = Whether to make sure the BAI file is present
        :param csi: bool [False] = Whether to make sure the CSI file is present
        :param needed: bool [False] = Whether files are needed (raises IOErrors if missing)
        :return: True/False whether (indexed?) BAM is found
        '''
        ### ~ [1] BAM file check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not rje.exists(bamfile):
            if needed: raise IOError('Cannot find BAM file "{0}" (bam=FILE)'.format(bamfile))
            return False
        if os.system("samtools quickcheck {0}".format(bamfile)):
            raise RuntimeError("samtools quickcheck failed for {0}".format(bamfile))
        if not makeindex: return True
        seqin = self.getStr('SeqIn')
        if rje.exists(seqin) and self.needToRemake(bamfile, seqin):
            if needed: raise IOError('SeqIn file "{0}" younger than BAM file "{1}"'.format(seqin,bamfile))
            self.printLog('#BAM', 'SeqIn younger than BAM file "{}": regenerating'.format(bamfile))
            return False
        ### ~ [2] Index BAM file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        baifile = '{}.bai'.format(bamfile)
        csifile = '{}.csi'.format(bamfile)
        ## ~ [2a] BAI Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if not csi:
            if self.needToRemake(baifile, bamfile):
                makebai = 'samtools index -b {0} {1}.bai'.format(bamfile, bamfile)
                logline = self.loggedSysCall(makebai, append=True, threaded=False,
                                             nologline='No stdout from samtools index')
                self.printLog('#BAI','Indexed BAM file (-b).')
            if rje.exists(baifile): return True
            elif needed: raise IOError('Cannot find BAI index file "{0}"'.format(baifile))
            else: return False
        ## ~ [2b] CSI Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.needToRemake(csifile, bamfile):
            makebai = 'samtools index -c {0} {1}.csi'.format(bamfile, bamfile)
            logline = self.loggedSysCall(makebai, append=True, threaded=False,
                                         nologline='No stdout from samtools index')
            self.printLog('#CSI','Indexed CSI file (-c).')
        if rje.exists(csifile): return True
        elif needed: raise IOError('Cannot find CSI index file "{0}"'.format(csifile))
        else: return False
#########################################################################################################################
    def samToBam(self,prefix):  ### Converts SAM file to sorted BAM file, then checks
        '''
        Converts SAM file to sorted and indexed BAM file, then checks
        >> prefix:str = Prefix for SAM file and BAM file.
        :return: bamfile [str]
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tmpfile = '{0}.tmp.bam'.format(prefix)
            bamfile = '{0}.bam'.format(prefix)
            maplog = '{0}.log'.format(prefix)
            ### ~ [2] Converting SAM to BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Try multithreaded ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#BAM', 'SAM to BAM conversion')
            try:
                sam2bam = 'samtools view -bo {0}.tmp.bam -@ {1} -S {2}.sam'.format(prefix,self.threads()-1,prefix)
                logline = self.loggedSysCall(sam2bam, maplog, append=True, nologline='No stdout from sam2bam')
                if not self.checkBAMFile(tmpfile,makeindex=False,needed=True): raise IOError('Problem with temp BAM file')
            except:
                self.printLog('#BAM','Converting SAM to BAM. Using a single thread due to missing data.')
                sam2bam = 'samtools view -bo {0}.tmp.bam -S {1}.sam'.format(prefix,prefix)
                logline = self.loggedSysCall(sam2bam,maplog,append=True,nologline='No stdout from sam2bam')
                if not self.checkBAMFile(tmpfile, makeindex=False, needed=True): raise IOError('Problem with temp BAM file')
            ## ~ [2c] Sorting BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#BAM','Sorting BAM file.')
            bamsort = 'samtools sort -@ {0} -o {1}.bam -m 6G {2}.tmp.bam'.format(self.threads()-1,prefix,prefix)
            logline = self.loggedSysCall(bamsort,maplog,append=True)
            if not self.checkBAMFile(bamfile, makeindex=False, needed=True): raise IOError('Problem with sorted BAM file')
            os.unlink('{0}.tmp.bam'.format(prefix))
            if self.debugging: self.printLog('#SAM','Keeping "{0}.sam" (debug=T): delete later'.format(prefix))
            else:
                os.unlink('{0}.sam'.format(prefix))
            return bamfile
        except:
            self.errorLog('{0}.samToBam() error'.format(self.prog())); raise
#########################################################################################################################
    def checkReadList(self,expect=True):    ### Checks and fixes reads=LIST and readtype=LIST
        '''Checks and fixes reads=LIST and readtype=LIST.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Check reads
            if not self.list['Reads']:
                if expect: raise IOError('No reads=FILELIST files found')
                self.printLog('#READS','No reads=FILELIST files given.')
                return False
            if not rje.checkForFiles(filelist=self.list['Reads'],basename='',log=self.log,cutshort=False,ioerror=expect):
                return False
            #i# Check read types
            if not self.list['ReadType']:
                self.warnLog('Read Type not given (pb/ont): check readtype=LIST. Will use "ont".')
                self.list['ReadType'] = ['ont'] * len(self.list['Reads'])
            elif len(self.list['ReadType']) == 1 and len(self.list['Reads']) != 1:
                self.printLog('#READS','Using "%s" as read type for all long reads' % self.list['ReadType'][0])
                self.list['ReadType'] = [self.list['ReadType'][0]] * len(self.list['Reads'])
            elif len(self.list['ReadType']) > len(self.list['Reads']):
                self.warnLog('reads=FILELIST < readtype=LIST length mismatch: will truncate readtype=LIST.')
                self.list['ReadType'] = self.list['ReadType'][:len(self.list['Reads'])]
            elif len(self.list['ReadType']) != len(self.list['Reads']):
                self.warnLog('reads=FILELIST vs readtype=LIST length mismatch: check readtype=LIST. Will cycle if needed.')
                rx = 0
                while len(self.list['ReadType']) < len(self.list['Reads']):
                    self.list['ReadType'].append(self.list['ReadType'][rx])
                    rx += 1
            self.printLog('#READS','{0} read files ({1})'.format(len(self.list['ReadType']),'/'.join(self.list['ReadType'])))
            return True
        except:
            self.errorLog('{0}.checkReadList() error'.format(self.prog())); raise
#########################################################################################################################
    def iMax(self):
        imax = 4
        seqinsize = os.path.getsize(self.getStr('SeqIn')) / 1e9
        while seqinsize > imax: imax += 1
        return imax
#########################################################################################################################
    def longreadMinimap(self,paf=False):  ### Performs long read versus assembly minimap2 and to converts BAM or PAF
        '''
        Performs long read versus assembly minimap2 and converts to BAM file
        >> paf:bool [False] = Whether to output PAF file rather than BAM file.
        :return: bamfile/None
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Need to sort out consistent use of force
            #?# Check and update self.list['ReadType'] and self.list['Reads'] in a dedicated method
            if not self.getBool('Minimap2'):
                raise RuntimeError('Cannot find minimap2!')
            seqin = self.getStr('SeqIn')
            outfile = self.setBAMFile()
            if paf:
                outfile = self.setPAFFile()
            if not self.force() and not self.needToRemake(outfile, seqin): return outfile
            ### ~ [2] ~ Generate individual Map files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.checkReadList(): raise IOError('Problem with reads for mapping')
            self.headLog('Map long reads ({0} output)'.format({True:'PAF',False:'BAM'}[paf]), line='-')
            rje.backup(self,outfile)
            bamlist = []; rx = 0
            for readfile in self.list['Reads']:
                rtype = self.list['ReadType'][rx]; rx +=1
                if rtype in ['pacbio','pac']: rtype = 'pb'
                if rtype in ['hifi','ccs']: rtype = 'hifi'
                if rtype not in ['ont','pb','hifi']:
                    self.warnLog('Read Type "%s" not recognised (pb/ont): check readtype=LIST. Will use "ont".' % rtype)
                    rtype = 'ont'
                prefix = '{0}.{1}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(readfile,strip_path=True))
                maplog = '{0}.log'.format(prefix)
                ## ~ [2a] Map reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                imax = 4
                seqinsize = os.path.getsize(self.getStr('SeqIn')) / 1e9
                while seqinsize > imax: imax += 1
                mapcmd = '{0} -t {1} -I {2}G -L --secondary=no'.format(self.minimap2(),self.threads(),imax)
                ## ~ [2b] ~ Check for existing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if paf:
                    mapfile = '{0}.paf'.format(prefix)
                    mapcmd = mapcmd + ' -o {0}.paf -x map-{1} {2} {3}'.format(prefix,rtype,seqin,readfile)
                else:
                    bamfile = '{0}.bam'.format(prefix)
                    mapfile = '{0}.sam'.format(prefix)
                    mapcmd = mapcmd + ' -o {0}.sam -ax map-{1} {2} {3}'.format(prefix,rtype,seqin,readfile)
                runmap = self.force() or self.needToRemake(mapfile,seqin) or self.needToRemake(mapfile,readfile) or not os.path.getsize(mapfile)
                ## ~ [2c] ~ Generate map file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Extra BAM check
                if not paf and self.checkBAMFile(bamfile):
                    runmap = False
                    mapfile = bamfile
                if runmap:
                    logline = self.loggedSysCall(mapcmd,maplog,append=True)
                else:
                    self.printLog('#MAP','Will use existing file "{0}" (force=F).'.format(mapfile))
                if not os.path.getsize(mapfile): raise IOError('{0} not generated correctly'.format(mapfile))
                ## ~ [2b] Converting SAM to BAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not paf and mapfile != bamfile:
                    mapfile = self.samToBam(prefix)
                bamlist.append(mapfile)
            ### ~ [3] ~ Merge individual Map files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if paf:
                paffile = outfile
                if len(bamlist) > 1:
                    bammerge = 'cat {0} > {1}'.format(' '.join(bamlist), paffile)
                    logline = self.loggedSysCall(bammerge, append=True)
                    if not rje.exists(paffile): raise IOError('Merged PAF file "%s" not generated' % paffile)
                    for sortbam in paflist: os.unlink(sortbam)
                else:
                    os.rename(bamlist[0], paffile)
                self.printLog('#PAF','{0} generated.'.format(paffile))
                return paffile
            else:
                if len(bamlist) > 1:
                    # samtools merge - merges multiple sorted input files into a single output.
                    bammerge = 'samtools merge -@ {0} {1} {2}'.format(self.threads()-1,outfile,' '.join(bamlist))
                    logline = self.loggedSysCall(bammerge,append=True)
                else:
                    os.rename(bamlist[0], outfile)
                if not self.checkBAMFile(outfile,makeindex=True,needed=True):
                    raise IOError('Merged BAM file "%s" not generated' % outfile)
                self.printLog('#BAM','{0} generated.'.format(outfile))
                return outfile
        except:
            self.errorLog('{0}.longreadMinimap({1}) error'.format(self.prog(),{True:'PAF',False:'BAM'}[paf]))
            raise
#########################################################################################################################
    def cigarSum(self,cstr,cigre):    ### Returns a dictionary of CIGAR element and sum based on cigre regex
        '''
        Returns a dictionary of CIGAR element and sum based on cigre regex.
        :param cstr: CIGAR string to parse
        :param cigre: Dictionary of regex
        :return: Dictionary of sums
        '''
        cigsum = {}
        for (x,reX) in cigre.items():
            cigsum[x] = 0
            if x in cstr:
                cigsum[x] = sum(map(int, reX.findall(cstr)))
        return cigsum
#########################################################################################################################
    def cigIndelRatio(self,cigsum,logsum=True): ### Returns indel ratio from cigsum dictionary generated by cigarSum()
        '''
        Returns indel ratio from cigsum dictionary generated by cigarSum().
        :param cigsum:
        :return: num
        '''
        base = 0
        for x in 'XM=':
            if x in cigsum: base +=  cigsum[x]
        inscount = base + cigsum['I']
        delcount = base + cigsum['D']
        indelratio = 0.0
        if delcount: indelratio = inscount / float(delcount)
        if logsum:
            self.printLog('#MAPINS','{0} mapped + {1} insertions = {2}'.format(dnaLen(base),dnaLen(cigsum['I']),dnaLen(inscount)))
            self.printLog('#MAPDEL','{0} mapped + {1} deletions = {2}'.format(dnaLen(base),dnaLen(cigsum['D']),dnaLen(delcount)))
            self.printLog('#INSDEL','Insertion:Deletion ratio = {0}'.format(rje.dp(indelratio,3)))
        return indelratio
#########################################################################################################################
    def indelRatio(self,density=False):   ### Checks or calculates the indel ratio for a BAM file using samtools
        '''
        Checks or calculates the indel ratio for a BAM file using samtools.
        >> density:bool [False] = Whether to use density mode
        << indelratio:num = Insertion:Deletion ratio parsed from cigar string.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cigs = 'IDXM='  # CIGAR string elements to parse
            cigre = {}
            for x in cigs:
                cigre[x] = re.compile('(\d+){0}'.format(x))
            #reD = re.compile('(\d+)D')
            #reI = re.compile('(\d+)I')
            bamfile = self.getBAMFile()
            ratiofile = bamfile + '.indelratio.txt'
            indelratio = 0.0
            if density:
                ratiofile = bamfile + '.indelratiovec.txt'
            elif not self.force() and rje.exists(ratiofile):
                try:
                    cstr = open(ratiofile,'r').readline().split()[0]
                    indelratio = self.cigIndelRatio(self.cigarSum(cstr,cigre))
                    if not indelratio: raise ValueError('{0} has incorrect format'.format(ratiofile))
                except:
                    self.errorLog('Problem with indelratio file {0}: deleting'.format(ratiofile))
                    os.unlink(ratiofile)
            ### ~ [2] Parse BAM CIGAR strings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile,'r').readline().split()[0]:
                if not self.getBool('Samtools'): raise IOError('Cannot find samtools')
                if not self.checkInput(busco=True): raise IOError('No busco=FILE table')
                bedfile = rje.baseFile(self.getStr('BUSCO'),strip_path=True) + '.bed'
                if bedfile == 'full_table.bed':
                    bedfile = self.baseFile() + '.full_table.bed'
                if self.force() or not rje.checkForFiles([bedfile],log=self.log):
                    os.system('grep Complete %s | awk \'{print $3"\t"$4"\t"$5;}\' > %s' % (self.getStr('BUSCO'),bedfile))
                    if rje.exists(bedfile):
                        self.printLog('#BUSCO','Converted BUSCO Complete to {0}'.format(bedfile))
                    else:
                        raise IOError('BED file generation failed.')
                self.printLog('#CIGAR', "samtools view -h -F 4 %s -L %s | grep -v '^@' | awk '{print $6;}' | uniq" % (bamfile,bedfile))
                if density:
                    RAT = open(ratiofile,'w')
                    CIG = os.popen("samtools view -h -F 4 %s -L %s | grep -v '^@' | awk '{print $6;}' | uniq" % (bamfile,bedfile))
                    cstr = CIG.readline(); cx = 1
                    while cstr:
                        self.progLog('\r#INDEL', 'Calculating indel ratio: {0} reads'.format(cx))
                        #delcount = sum(map(int, reD.findall(cstr)))
                        #inscount = sum(map(int, reI.findall(cstr)))
                        insrat = self.cigIndelRatio(self.cigarSum(cstr,cigre),logsum=False)
                        RAT.write('{0}\n'.format(insrat))
                        cstr = CIG.readline()
                        cx += 1
                    RAT.close()
                    CIG.close()
                    self.printLog('\r#INDEL', 'Calculated indel ratio: {0} reads'.format(rje.iStr(cx)))
                else:
                    self.progLog('\r#INDEL', 'Calculating indel ratio...')
                    cstr = os.popen("samtools view -h -F 4 %s -L %s | grep -v '^@' | awk '{print $6;}' | uniq" % (bamfile,bedfile)).read()
                    #delcount = sum(map(int, reD.findall(cstr)))
                    #inscount = sum(map(int, reI.findall(cstr)))
                    #open(ratiofile, 'w').write('{0} {1}\n'.format(inscount,delcount))
                    cigsum = self.cigarSum(cstr, cigre)
                    indelratio = self.cigIndelRatio(cigsum)
                    cigsumstr = ''
                    for (x,reX) in cigre.items():
                        cigsumstr += '{0}{1}'.format(cigsum[x],x)
                    if indelratio:
                        open(ratiofile, 'w').write('{0}\n'.format(cigsumstr))
                        self.printLog('\r#INDEL', 'Saved indel ratio data to {0}'.format(ratiofile))
                    else:
                        self.warnLog('Calculated indelratio=0.0. Something probably went wrong. Indel ratio data not saved.')
            elif density:
                self.printLog('\r#INDEL', 'Use indel ratio vectors from {0} (force=F)'.format(ratiofile))
            else:
                self.printLog('\r#INDEL','Use indel ratio from {0} (force=F)'.format(ratiofile))
            ### ~ [3] Calculate indel ratio if needed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if density:
                rdir = self.rDir()
                indelratio = float(rje.chomp(os.popen('Rscript {0}depmodepure.R {1}'.format(rdir, ratiofile)).readlines()[0]))

            self.printLog('#INDEL','{0} indel ratio = {1:.3f}'.format(bamfile,indelratio))
            self.setNum({'IndelRatio':indelratio})
            return indelratio
        except: self.errorLog('%s.indelRatio error' % self.prog()); raise
#########################################################################################################################
    def baseCount(self):    ### Reads/calculates total sequencing bases
        '''
        Reads/calculates total sequencing bases.
        << readbp:num = Total read basecount.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            readbp = self.getNum('ReadBP')
            if readbp: return readbp
            readbp = 0
            if not self.list['Reads']: self.setInt({'ReadBP':readbp}); return readbp
            ### ~ [2] Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for rfile in self.list['Reads']:
                countfile = '{}.basecount.txt'.format(rfile)
                if rje.exists(countfile) and not self.force():
                    filebp = int(rje.chomp(open(countfile,'r').readline()))
                    self.printLog('#READBP','%s: %s bases in %s' % (countfile,rje.iStr(filebp),rfile))
                    readbp += filebp
                    continue
                self.progLog('\r#READBP','Counting bases in {}...'.format(rfile))
                gzip = rfile.endswith('.gz')
                fastq = rfile.endswith('q.gz') or rfile.endswith('q')
                if gzip and fastq:
                    filebp = int(os.popen("zcat %s | grep '^+$' -B1 | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                elif fastq:
                    filebp = int(os.popen("grep '^+$' -B1 %s | grep -v '^+$' | grep -v '^--' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                elif gzip:
                    filebp = int(os.popen("zcat %s | grep -v '^>' | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                else:
                    filebp = int(os.popen("grep -v '^>' %s | wc | awk '{ $4 = $3 - $2 } 1' | awk '{print $4}'" % rfile).readline())
                self.printLog('#READBP','%s bases counted in %s' % (rje.iStr(filebp),rfile))
                open(countfile,'w').write('%d\n' % filebp)
                readbp += filebp
            self.setInt({'ReadBP':readbp})
            self.printLog('#READBP', 'Total base count (unadjusted): %s' % (rje.iStr(readbp)))
            return readbp
        except: self.errorLog('%s.baseCount error' % self.prog())
#########################################################################################################################
    def covBases(self): ### Partial mapadjust calculation just calculating covbases.
        '''
        Partial mapadjust calcuation just calculating covbases.
        >> allbases:bool [False] = Whether to use all bases, not just covered bases, for calculation.
        << covbases:num = Covbases portion of mapajdust ratio.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('CovBases'): return self.getNum('CovBases')
            bamfile = self.getBAMFile()
            ratiofile = bamfile + '.mapratio.txt'
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile, 'r').readline().split()[0]:
                ratiofile = bamfile + '.covbases.txt'
            ### ~ [2] CovBases Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile,'r').readline().split()[0]:
                if not self.getBool('Samtools'): raise IOError
                self.printLog('\r#ADJUST', 'Calculating CovBases adjustment (samtools coverage)...')
                covbases = float(os.popen("samtools coverage {0} | grep -v coverage | awk  '{{sum += ($7 * $5)}} END {{print sum}}'".format(bamfile)).read().split()[0])
                open(ratiofile,'w').write('{0}\n'.format(covbases))
            else:
                self.printLog('\r#ADJUST','Use coverage data from from {0} (force=F)'.format(ratiofile))
                ratio = open(ratiofile, 'r').read().split()
                covbases = float(ratio[0])
            self.printLog('\r#ADJUST','Reference base read coverage: {0}'.format(rje.iStr(int(covbases))))
            self.setNum({'CovBases': covbases})
            return covbases
        except: self.errorLog('%s.covBases error' % self.prog()); raise
#########################################################################################################################
    def allBases(self): ### Partial mapadjust calculation just calculating improved covbases.
        '''
        Partial mapadjust calcuation just calculating covbases.
        << covbases:num = Covbases portion of mapajdust ratio.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('AllBases'): return self.getNum('AllBases')
            bamfile = self.getBAMFile()
            ### ~ [2] CovBases Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ratiofile = bamfile + '.allbases.txt'
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile,'r').readline().split()[0]:
                if not self.getBool('Samtools'): raise IOError
                self.printLog('\r#ADJUST', 'Calculating CovBases adjustment (samtools coverage)...')
                covbases = float(os.popen("samtools coverage {0} | grep -v coverage | awk  '{{sum += ($7 * ($3 - $2 + 1))}} END {{print sum}}'".format(bamfile)).read().split()[0])
                open(ratiofile,'w').write('{0}\n'.format(covbases))
            else:
                self.printLog('\r#ADJUST','Use coverage data from from {0} (force=F)'.format(ratiofile))
                ratio = open(ratiofile, 'r').read().split()
                covbases = float(ratio[0])
            self.printLog('\r#ADJUST','Reference base read coverage: {0}'.format(rje.iStr(int(covbases))))
            self.setNum({'AllBases': covbases})
            return covbases
        except: self.errorLog('%s.allBases error' % self.prog()); raise
#########################################################################################################################
    def mapBases(self): ### Partial mapadjust calculation just calculating mapbases.
        '''
        Partial mapadjust calcuation just calculating mapbases.
        << mapbases:num = Mapbases portion of mapajdust ratio.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('MapBases'): return self.getNum('MapBases')
            bamfile = self.getBAMFile()
            spliti = 1
            ratiofile = bamfile + '.mapratio.txt'
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile, 'r').readline().split()[0]:
                ratiofile = bamfile + '.mapbases.txt'
                spliti = 0
            ### ~ [2] MapBases Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.force() or not rje.exists(ratiofile) or not open(ratiofile,'r').readline().split()[0]:
                if not self.getBool('Samtools'): raise IOError
                self.printLog('\r#ADJUST', 'Calculating MapAdjust ratio (samtools fasta)...')
                mapbases = float(os.popen("samtools view -hb -F 4 {0} | samtools fasta - | grep -v '^>' | wc | awk '{{ $4 = $3 - $2 }} 1' | awk '{{print $4}}'".format(bamfile)).read().split()[0])
                open(ratiofile,'w').write('{0}\n'.format(mapbases))
            else:
                self.printLog('\r#ADJUST','Use coverage data from from {0} (force=F)'.format(ratiofile))
                ratio = open(ratiofile, 'r').read().split()
                mapbases = float(ratio[spliti])
            self.printLog('\r#ADJUST','Mapped read bases: {0}'.format(rje.iStr(int(mapbases))))
            self.setNum({'MapBases': mapbases})
            return mapbases
        except: self.errorLog('%s.mapBases error' % self.prog()); raise
#########################################################################################################################
    def mapAdjust(self,allbases=True):    ### Legacy read mapping adjustment method
        '''
        Legacy read mapping adjustment method.
        >> allbases:bool [True] = Whether to use
        << readbp:num = Legacy mapadjust ratio.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('MapAdjust'): return self.getNum('MapAdjust')
            bamfile = self.getBAMFile()
            if allbases:
                covbases = self.allBases()
            else:
                covbases = self.covBases()
            mapbases = self.mapBases()
            ### ~ [2] MapAdjust Total Sequencing Bases ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not covbases or not mapbases: raise ValueError('Cannot have zero values for mapadjust. Something has gone wrong!')
            mapadjust = covbases / mapbases
            if allbases:
                self.printLog('#ADJUST', 'MapAdjust ratio: {0:.3f}'.format(mapadjust))
                self.setNum({'MapAdjust':mapadjust})
            else:
                self.printLog('#ADJUST', 'OldAdjust ratio: {0:.3f}'.format(mapadjust))
                self.setNum({'OldAdjust':mapadjust})
            return mapadjust
        except: self.errorLog('%s.mapAdjust error' % self.prog()); raise
#########################################################################################################################
    ### <5> ### SC Depth Methods                                                                                        #
#########################################################################################################################
    def getSCDepth(self):  ### Loads or calculates SC depth from depthcopy.R Rscript.
        '''
        Loads or calculates SC depth from depthcopy.R Rscript.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            depfile = self.getFastDep()         # This will generate the BAM file if needed
            scdfile = '{0}.scdepth'.format(depfile)
            if self.needToRemake(scdfile,depfile):
                optionstr = self.makeOptionStr(['depfile','busco','adjust','basefile'])
                self.callRscript(optionstr)
            if not rje.checkForFile(scdfile): raise IOError('Cannot find "{0}"'.format(scdfile))
            scdepth = float(rje.chomp(open(scdfile,'r').readline()))
            self.setNum({'SCDepth':scdepth})
            self.printLog('#SCDEP','Single copy read depth = {0:.2f}X'.format(scdepth))
            return scdepth
        except: self.errorLog('%s.getSCDepth error' % self.prog()); raise
#########################################################################################################################
    def makeOptionStr(self,cmds=None):  ### Returns the option string for the depthcopy.R Rscript.
        '''
        Returns the option string for the depthcopy.R Rscript.
        depfile=FILE [busco=FILE] [scdepth=INT] [regfile=FILE] [reghead=LIST] [gfftype=LIST] [winsize=INT] [winstep=NUM]
        >> cmds:list [None] = list of command dictionary keys
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            options = ['pngdir={0}.plots'.format(self.baseFile())]
            if not cmds: cmds = ['depfile','busco','scdepth','regfile','gfftype','winsize','winstep','adjust','cnmax','basefile']
            for cmd in cmds:
                val =  self.getData(cmd)
                if val and val != 'None':
                    options.append('{0}={1}'.format(cmd,val))
            for lcmd in rje.sortKeys(self.list):
                if lcmd.lower() in cmds:
                    options.append('{0}={1}'.format(lcmd.lower(), ','.join(self.list[lcmd])))
            if self.debugging(): options.append('debug=TRUE')
            if self.getBool('SeqStats'): options.append('seqstats=TRUE')
            optionstr = ' '.join(options)
            return optionstr
        except: self.errorLog('%s.callRscript error' % self.prog())
#########################################################################################################################
    def rDir(self,rscript='depthcopy.R'):
        if rje.exists(mypath+rscript): return mypath
        else: return '%slibraries/r/' % slimsuitepath
#########################################################################################################################
    def callRscript(self,optionstr=''):  ### Calls the depthcopy.R Rscript and parses output.
        '''
        Calls the depthcopy.R Rscript and parses output.
        >> optionstr:str = String of options to give to depthcopy.R script.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rdir = self.rDir()
            if not optionstr: optionstr = self.makeOptionStr()
            ### ~ [2] Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            complete = False
            rcmd = 'Rscript {0}depthcopy.R {1}'.format(rdir, optionstr)
            self.printLog('#RCMD',rcmd)
            RCMD = os.popen(rcmd)
            rline = RCMD.readline()
            while rline:
                if '] #' in rline:
                    rline = ' '.join(rline.split('] ')[1:])
                    logstr = rje.chomp(rline).split()
                    self.printLog('{0}'.format(logstr[0].upper()), ' '.join(logstr[1:]))
                elif rline[:1] == '[':
                    self.verbose(v=0,text=rje.chomp(rline),newline=0)
                else:
                    self.verbose(v=1,text=rje.chomp(rline),newline=0)
                #!# Parse scdepth and other key points to printLog
                complete = complete or 'DepthCopy.R finished' in rline
                rline = RCMD.readline()
            RCMD.close()
            return complete
        except: self.errorLog('%s.callRscript error' % self.prog())
#########################################################################################################################
    def calculateGenomeSize(self,adjust=None,save=True,benchmark=False):  ### Calculates genome size from stored values and reports
        '''
        Calculates genome size from stored values and reports
        >> adjust:str ['indel'] = Adjustment to use for stored and reported genome size ('none','indel','mapadjust')
        >> save:bool [True] = whether to save the genome size estimate table
        '''
        try:### ~ [0] Set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Calculate Genome Size', line='-')
            readbp = self.baseCount()
            if adjust == 'None': adjust = None
            gdb = self.db('gensize')
            if not gdb: gdb = self.db().addEmptyTable('gensize',['SeqFile','DepMethod','Adjust','ReadBP','MapAdjust','SCDepth','EstGenomeSize'],['SeqFile','DepMethod','Adjust'])
            depmethod = 'mpileup'
            if self.getBool('QuickDepth'): depmethod = 'depth'
            adjustments = [None,'IndelRatio','MapAdjust','MapRatio','CovBases','OldAdjust','AllBases','MapBases','OldCovBases']
            if self.getNum('IndelRatio') and self.getNum('MapBases'):
                mapadjust = self.getNum('MapBases') / self.getNum('IndelRatio')
                self.setNum({'MapRatio':mapadjust})
                self.printLog('#ADJUST', 'MapRatio mapping and indel adjustment: {0} -> {1}'.format(dnaLen(readbp),dnaLen(mapadjust)))
            if self.getNum('CovBases') and self.getNum('MapBases'):
                mapadjust = self.getNum('CovBases') / self.getNum('MapBases')
                self.setNum({'OldAdjust':mapadjust})
                self.printLog('#ADJUST', 'OldAdjust MapRatio: {0:.3f}'.format(mapadjust))
            if benchmark:
                adjustments += ['Assembly','MeanX']
                seqin = self.seqinObj()
                seqlen = 0
                for seq in seqin.seqs():
                    seqlen += seqin.seqLen(seq)
                self.setNum({'Assembly':seqlen,'MeanX':self.getNum('AllBases')/float(seqlen)})
            if not save: adjustments = [adjust]
            scdepth = self.getSCDepth()
            ### ~ [1] Calculate genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for nkey in adjustments:
                if nkey and not self.getNum(nkey): continue
                mapadjust = 1.0
                desc = 'Unadjusted'
                if nkey:
                    mapadjust = self.getNum(nkey)
                    desc = nkey
                if nkey in ['CovBases','AllBases','MapBases','MapRatio']:
                    adjreadbp = mapadjust
                    mapadjust = adjreadbp / readbp
                elif nkey == 'IndelRatio':
                    adjreadbp = readbp / mapadjust
                elif nkey == 'Assembly':
                    adjreadbp = scdepth * mapadjust
                    mapadjust = adjreadbp / readbp
                elif nkey == 'MeanX':
                    scdepth = mapadjust
                    adjreadbp = readbp
                    mapadjust = 1.0
                else:
                    adjreadbp = readbp * mapadjust
                if nkey:
                    self.printLog('#READBP', 'Total base count ({0} adjusted): {1}'.format(nkey,rje.iStr(adjreadbp)))
                estgensize = int(0.5+(float(adjreadbp) / scdepth))
                scprint = rje.dp(scdepth, 2)
                akey = nkey
                if nkey == 'AllBases': akey = 'CovBases'
                if nkey == 'CovBases': akey = 'OldCovBases'
                self.printLog('#GSIZE','{0} ({1}) estimated genome size ({2} at {3}X): {4}'.format(desc,depmethod,rje_seqlist.dnaLen(adjreadbp,dp=0,sf=3),scprint,rje_seqlist.dnaLen(estgensize,dp=0,sf=4)))
                gdb.addEntry({'SeqFile':os.path.basename(self.getStr('SeqIn')),'DepMethod':depmethod,
                                  'Adjust':akey,'ReadBP':readbp,'MapAdjust':mapadjust,
                                  'SCDepth':scprint,'EstGenomeSize':estgensize})
                if not nkey and not adjust:
                    self.setInt({'EstGenomeSize':estgensize})
                elif nkey and nkey.lower() == adjust.lower():
                    self.setInt({'EstGenomeSize':estgensize})
            if save: gdb.saveToFile()
            return self.getInt('EstGenomeSize')
        except:
            self.errorLog('{0}.calculateGenomeSize() error'.format(self.prog()))
            raise
#########################################################################################################################
### End of SECTION II: ReadCore Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def dnaLen(seqlen,dp=2,sf=3): return rje_seqlist.dnaLen(seqlen,dp,sf)
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
    try: ReadCore(mainlog,cmd_list).run()

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
