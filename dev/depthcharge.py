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
Module:       DepthCharge
Description:  Genome assembly quality control and misassembly repair
Version:      0.3.0
Last Edit:    24/07/24
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    DepthCharge is an assembly quality control and misassembly repair program. It uses mapped long read depth of
    coverage to charge through a genome assembly and identify coverage "cliffs" that may indicate a misassembly.
    If appropriate, it will then blast the assembly into fragment at those misassemblies.

    DepthCharge uses a genome assembly and PAF file of mapped reads as input. If no file is provided, minimap2 will
    be used to generate one.

    For each sequence, DepthCharge starts at the beginning of the sequence and scans through the PAF file for
    coverage to drop below the `mindepth=INT` threshold (default = 1 read). These positions are marked as "bad" and
    compressed into regions of adjacent bad positions. Regions at the start or end of a sequnece are labelled "end".
    Regions overlapping gaps are labelled "gap". Otherwise, regions are labelled "bad". All regions are output to
    `*.depthcharge.tdt` along with the length of each sequence (region type "all").

    Future versions will either fragment the assembly at "bad" regions (and "gap" regions if 'breakgaps=T`. If
    `breakmode=gap` then DepthCharge will replace bad regions with a gap (`NNNN...`) of length `gapsize=INT`. If
    `breakmode=report` then no additional processing of the assembly will be performed. Otherwise, the processed
    assembly will be saved as `*.depthcharge.fasta`.

Commandline:
    ### ~ Main DepthCharge run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Input sequence assembly [None]
    basefile=FILE   : Root of output file names [$SEQIN basefile]
    paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
    breakmode=X     : How to treat misassemblies (report/gap/fragment) [fragment]
    breakgaps=T/F   : Whether to break at gaps where coverage drops if breakmode=fragment [False]
    gapsize=INT     : Size of gaps to insert when breakmode=gap [100]
    mindepth=INT    : Minimum depth to class as OK [1]
    minspan=INT     : Minimum spanning bp at end of reads (trims from PAF alignments) [0]
    ### ~ PAF file generation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
    readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
    minimap2=PROG   : Full path to run minimap2 [minimap2]
    mapopt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm5]
    ### ~ Additional options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    dochtml=T/F     : Generate HTML DepthCharge documentation (*.docs.html) instead of main run [False]
    logfork=T/F     : Whether to log forking in main log [False]
    tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
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
import rje, rje_forker, rje_obj, rje_db, rje_rmd, rje_paf, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Removed endbuffer and gapbuffer in favour of straight overlap assignment.
    # 0.2.0 - Added HiFi read type.
    # 0.3.0 - Added minspan=INT : Minimum spanning bp at end of reads (trims from PAF alignments). Fixed forcing. [0]
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [?] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [ ] : Add minlen=INT setting to cull short outputs (and merge gaps?)
    # [ ] : Add sequence summary to follow sequence output.
    # [ ] : Add PE filter charge mode using TLEN filter.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('DepthCharge', '0.3.0', 'July 2024', '2021')
    description = 'Genome assembly quality control and misassembly repair'
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
        log = rje.setLog(info,out,cmd_list,py3warn=False)   # Sets up Log object for controlling log file output
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
class DepthCharge(rje_forker.Forker):
    '''
    DepthCharge Class. Author: Rich Edwards (2021).

    Str:str
    - BreakMode=X     : How to treat misassemblies (report/gap/fragment) [fragment]
    - PAF=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
    - SeqIn=FILE      : Input sequence assembly [None]

    Bool:boolean
    - BreakGaps=T/F   : Whether to break at gaps where coverage drops if breakmode=fragment [False]
    - DocHTML=T/F     : Generate HTML BUSCOMP documentation (*.info.html) instead of main run [False]

    Int:integer
    - EndBuffer=INT   : Buffer size for sequence ends and gaps to avoid breakage [500]
    - GapBuffer=INT   : Buffer size for sequence ends and gaps to avoid breakage [100]
    - GapSize=INT     : Size of gaps to insert when breakmode=gap [100]
    - MinDepth=INT    : Minimum depth to class as OK [1]
    - MinSpan=INT     : Minimum spanning bp at end of reads (trims from PAF alignments) [0]

    Num:float

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - PAF = rje_paf.PAF object controlling minimap2 read mapping and PAF generation.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        self.cmd_list = ['logfork=F','tuplekeys=T'] + self.cmd_list
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['BreakMode','PAF','SeqIn']
        self.boollist = ['BreakGaps','DocHTML']
        self.intlist = ['EndBuffer','GapBuffer','GapSize','MinDepth','MinSpan']
        self.numlist = []
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        self._setForkerAttributes()
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'BreakMode':'fragment'})
        self.setBool({'BreakGaps':False,'DocHTML':False})
        self.setInt({'EndBuffer':500,'GapBuffer':100,'GapSize':100,'MinDepth':1})
        self.setNum({})
        #i# Forker Defaults
        self.setBool({'RjePy':False,'LogFork':True,'KillMain':True})
        self.setInt({'IOLimit':50})
        self.setNum({'MemFree':0.0,'ForkSleep':0.0,'KillForks':36000,'KillTime':time.time()})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.obj['SeqIn'] = None
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self._forkerCmd()
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['BreakMode'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['PAF','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BreakGaps','DocHTML'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['EndBuffer','GapBuffer','GapSize','MinDepth','MinSpan'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.setStr({'BreakMode':self.getStrLC('BreakMode')})
        if self.getStrLC('BreakMode') not in ['report','fragment','gap']:
            self.warnLog('BreakMode "{0}" not recognised: setting to breakmode=report'.format(self.getStrLC('BreakMode')))
            self.setStr({'BreakMode':'report'})
        self.setInt({'Forks':max(1,self.getInt('Forks'))})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        # DepthCharge: genome assembly quality control and misassembly repair.

        DepthCharge is an assembly quality control and misassembly repair program. It uses mapped long read depth of
        coverage to charge through a genome assembly and identify coverage "cliffs" that may indicate a misassembly.
        If appropriate, it will then blast the assembly into fragment at those misassemblies.

        DepthCharge uses a genome assembly and PAF file of mapped reads as input. If no file is provided, minimap2 will
        be used to generate one.

        For each sequence, DepthCharge starts at the beginning of the sequence and scans through the PAF file for
        coverage to drop below the `mindepth=INT` threshold (default = 1 read). These positions are marked as "bad" and
        compressed into regions of adjacent bad positions. Regions at the start or end of a sequnece are labelled "end".
        Regions overlapping gaps are labelled "gap". Otherwise, regions are labelled "bad". All regions are output to
        `*.depthcharge.tdt` along with the length of each sequence (region type "all").

        Future versions will either fragment the assembly at "bad" regions (and "gap" regions if `breakgaps=T`. If
        `breakmode=gap` then DepthCharge will replace bad regions with a gap (`NNNN...`) of length `gapsize=INT`. If
        `breakmode=report` then no additional processing of the assembly will be performed. Otherwise, the processed
        assembly will be saved as `*.depthcharge.fasta`.

        ---

        # Running DepthCharge

        DepthCharge is written in Python 2.x and can be run directly from the commandline:

            python $CODEPATH/depthcharge.py [OPTIONS]

        If running as part of [SLiMSuite](http://slimsuite.blogspot.com/), `$CODEPATH` will be the SLiMSuite `tools/`
        directory. If running from the standalone [DepthCharge git repo](https://github.com/slimsuite/depthcharge), `$CODEPATH`
        will be the path the to `code/` directory. Please see details in the [DepthCharge git repo](https://github.com/slimsuite/depthcharge)
        for running on example data.

        ## Dependencies

        DepthCharge uses `grep` and `awk`. To generate documentation with `dochtml`, R will need to be installed and a
        pandoc environment variable must be set, e.g.

            export RSTUDIO_PANDOC=/Applications/RStudio.app/Contents/MacOS/pandoc

        If a PAF file is not provided, [minimap2](https://github.com/lh3/minimap2) must be installed and either added to
        the environment `$PATH` or given with the `minimap2=PROG` setting.

        For full documentation of the DepthCharge workflow, run with `dochtml=T` and read the `*.docs.html` file generated.


        ## Commandline options


        ```
        ### ~ Main DepthCharge run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FILE      : Input sequence assembly [None]
        basefile=FILE   : Root of output file names [$SEQIN basefile]
        paf=FILE        : PAF file of long reads mapped onto assembly [$BASEFILE.paf]
        breakmode=X     : How to treat misassemblies (report/gap/fragment) [fragment]
        breakgaps=T/F   : Whether to break at gaps where coverage drops if breakmode=fragment [False]
        gapsize=INT     : Size of gaps to insert when breakmode=gap [100]
        mindepth=INT    : Minimum depth to class as OK [1]
        minspan=INT     : Minimum spanning bp at end of reads (trims from PAF alignments) [0]
        ### ~ PAF file generation options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        reads=FILELIST  : List of fasta/fastq files containing reads. Wildcard allowed. Can be gzipped. []
        readtype=LIST   : List of ont/pb/hifi file types matching reads for minimap2 mapping [ont]
        minimap2=PROG   : Full path to run minimap2 [minimap2]
        mapopt=CDICT    : Dictionary of minimap2 options [N:100,p:0.0001,x:asm5]
        ### ~ Additional options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        dochtml=T/F     : Generate HTML Diploidocus documentation (*.docs.html) instead of main run [False]
        logfork=T/F     : Whether to log forking in main log [False]
        tmpdir=PATH     : Path for temporary output files during forking (not all modes) [./tmpdir/]
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        ```

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DocHTML'): return rje_rmd.docHTML(self)
            if not self.setup(): return False
            ### ~ [2] ~ DepthCharge ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Fork out processing of the PAF file for each input sequence.
            resfile = '{0}.depthcharge.tdt'.format(self.baseFile())
            ddb = self.depthChargeForker()  # depthcharge table - ['seqname','start','end','type']
            if not ddb: raise IOError('Generation of DepthCharge table failed')
            ddb.indexReport('type')
            breakup = 'bad' in ddb.index('type') or (self.getBool('BreakGaps') and 'gap' in ddb.index('type'))
            if breakup:
                ddb.printLog('#RESULT','Regions of bad coverage output to {0}'.format(resfile))
            elif 'gap' in ddb.index('type'):
                ddb.printLog('#RESULT','Gaps of bad coverage output to {0}'.format(resfile))
            else:
                ddb.printLog('#RESULT','No regions of bad coverage to output!')
            if self.getStrLC('BreakMode') == 'report' or not breakup: return True
            ### ~ [3] ~ Fragment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Fragment, insert gaps or just report regions
            fasfile = '{0}.depthcharge.fasta'.format(self.basefile())
            rje.backup(self,fasfile)
            FRAGFAS = open(fasfile,'w')
            seqin = self.seqinObj()
            seqx = 0
            for seq in seqin.seqs():
                (seqname,sequence) = seqin.getSeq(seq)
                sname = rje.split(seqname)[0]
                seqlen = len(sequence)
                regions = []
                for entry in ddb.indexEntries('seqname',sname):
                    if entry['type'] == 'bad' or (self.getBool('BreakGaps') and entry['type'] == 'gap'):
                        regions.append((entry['start'],entry['end']))
                if not regions:
                    FRAGFAS.write('>{0}\n{1}\n'.format(seqname,sequence)); seqx += 1
                    continue
                regions += [(0,0),(seqlen,seqlen)]
                regions.sort()
                fragx = 0
                newseq = ''
                while len(regions) > 1:
                    fragx += 1
                    if self.getStrLC('BreakMode') == 'gap':
                        if newseq: newseq += 'N' * self.getInt('GapSize')
                        newseq += sequence[regions[0][1]:regions[1][0]-1]
                    elif self.getStrLC('BreakMode') == 'fragment':
                        newname = '{0}.{1} {2}'.format(sname,fragx,seqname)
                        newseq = sequence[regions[0][1]:regions[1][0]-1]
                        FRAGFAS.write('>{0}\n{1}\n'.format(newname,newseq)); seqx += 1
                    regions.pop(0)
                if self.getStrLC('BreakMode') == 'gap':
                    newname = '{0}+{1}gaps {2}'.format(sname,fragx-1,seqname)
                    FRAGFAS.write('>{0}\n{1}\n'.format(newname,newseq)); seqx += 1
                    self.printLog('#ADDGAP','{0} gaps added to {1}'.format(fragx-1,sname))
                else:
                    self.printLog('#FRAG','{0} fragments of {1} output to {2}'.format(fragx,sname,fasfile))
            self.printLog('#FASOUT','{0} sequences output to {1}'.format(seqx,fasfile))
            # self.warnLog('BreakMode "{0}" not yet implemented!'.format(self.getStrLC('BreakMode')))
            return False
        except:
            self.errorLog(self.zen())
            return True   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''
        Main class setup method. This will load sequences into a SeqList object, gaps into a 'gaps' database table, and
        check or generate a PAF file from the mapped long reads.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            if not self.getStrLC('SeqIn'):
                raise ValueError('seqin=FILE must be set')
            if not rje.exists(self.getStr('SeqIn')):
                raise IOError('Unable to read seqin=FILE: "{0}"'.format(self.getStr('SeqIn')))
            seqbase = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if not self.getStrLC('Basefile'): self.baseFile(seqbase)
            if rje.checkForFiles(filelist=['.gaps.tdt'],basename=seqbase,log=self.log) and not self.force():
                self.cmd_list.append('gapstats=F')
            else:
                self.cmd_list.append('gapstats=T')
            seqin = self.seqinObj()
            gapdb = self.db().addTable('%s.gaps.tdt' % seqbase,mainkeys=['seqname','start','end'],name='gaps',ignore=[],expect=True)
            gapdb.dataFormat({'start':'int','end':'int'})
            ### ~ [2] PAF File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('PAF'):
                self.setStr({'PAF':self.baseFile()+'.paf'})
            pfile = self.getStr('PAF')
            if self.fullForce() or not rje.exists(pfile):
                paf = rje_paf.PAF(self.log,self.cmd_list)
                paf.longreadMinimapPAF(pfile)
            if not rje.exists(self.getStr('PAF')):
                raise IOError('Unable to read or create PAF file (fullforce={1}): {0}'.format(pfile,self.fullForce()))
            return True
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def seqinObj(self,summarise=True,gapstats=True): ### Returns the a SeqList object for the SeqIn file
        '''
        Returns the a SeqList object for the SeqIn file.
        :return: self.obj['SeqIn']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqbase = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if not self.obj['SeqIn']:
                seqcmd = self.cmd_list
                if summarise: seqcmd += ['summarise=T','dna=T','raw=F']
                gapstats = gapstats and (self.force() or not rje.exists('%s.gaps.tdt' % seqbase))
                if gapstats: seqcmd += ['gapstats']
                self.obj['SeqIn'] = rje_seqlist.SeqList(self.log,seqcmd+['autoload=T','seqmode=file','autofilter=F'])
                # sx = 0.0; stot = self.obj['SeqIn'].seqNum()
                # for seq in self.obj['SeqIn'].seqs():
                #     self.progLog('\r#CHECK','Checking sequences names: %.1f%%' % (sx/stot)); sx += 100.0
                #     if '|' in self.obj['SeqIn'].shortName(seq):
                #         raise ValueError('Pipe "|" characters found in seqin=FILE names: will break program. Please rename and try again.')
                # self.printLog('\r#CHECK','Checking sequences names complete.')
        except ValueError:
            self.printLog('\r#CHECK','Checking sequences names aborted.')
            self.errorLog('DepthCharge input sequence error'); raise
        except:
            self.errorLog('DepthCharge.seqinObj() error')
        return self.obj['SeqIn']
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=docs for program documentation and options. A plain text version is accessed with &rest=help.
        &rest=OUTFMT can be used to retrieve individual parts of the output, matching the tabs in the default
        (&rest=format) output. Individual `OUTFMT` elements can also be parsed from the full (&rest=full) server output,
        which is formatted as follows:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ... contents for OUTFMT section ...
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        ### Available REST Outputs:
        There is currently no specific help available on REST output for this program.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def depthChargeForker(self):  ### Main DepthCharge forking method
        '''
        Work through each sequence and fork it out for DepthCharge analysis.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.seqinObj()
            self.list['ToFork'] = seqin.list['Seq'][0:]
            resfile = '{0}.depthcharge.tdt'.format(self.baseFile())
            self.printLog('#MINX','Minimum X depth: mindepth={0}'.format(self.getInt('MinDepth')))
            self.printLog('#SPAN','Minimum bp spanning at end of reads: minspan={0}'.format(self.getInt('MinSpan')))
            if self.force(): rje.backup(self,resfile,appendable=False)
            elif rje.exists(resfile):
                self.warnLog('Results exist and force=F. Check that mindepth=INT and minspan=INT settings are unchanged.')
                ddb = self.db().addTable(resfile,['seqname','start','end','type'])
                ddb.dataFormat({'start':'int','end':'int'})
                complete = ddb.indexDataList('type','all','seqname')
                if complete:
                    cx = 0
                    for seq in self.list['ToFork'][0:]:
                        if seqin.shortName(seq) in complete: self.list['ToFork'].remove(seq); cx += 1
                    if cx: self.printLog('#SKIP','Skipping {0} previously processed sequences (force=F)'.format(rje.iStr(cx)))
                if not self.list['ToFork']:
                    self.printLog('#CHARGE','All sequences previously processed (force=F)')
                    return ddb
            while len(self.list['Forked']) < self.getNum('Forks') and self.list['ToFork']: self.nextFork()
            ### ~ [2] ~ Work through each sequence and fork out ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.forking()
            self.printLog('#FORK','Forking of %s jobs completed.' % (rje.iStr(seqin.seqNum())),log=self.getBool('LogFork'))
            ddb = self.db().addTable(resfile,['seqname','start','end','type'],replace=True)
            ddb.dataFormat({'start':'int','end':'int'})
            return ddb
        except: self.errorLog('%s.depthChargeForker error' % self.prog())
#########################################################################################################################
    def startFork(self,fdict):  ### Sets a new fork going using the data in fdict.
        '''Sets a new fork going using the data in fdict.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqin = self.seqinObj()
            fdict['seq'] = self.list['ToFork'].pop(0)
            (seqname,sequence) = seqin.getSeq(fdict['seq'])
            fdict['ID'] = 'Fork -%d' % (len(self.list['ToFork'])+1)
            fdict['FID'] = 'f_%s' % rje.randomString(6)
            fdict['Log'] = '%s%s.log' % (self.getStr('RunPath'),fdict['FID'])
            fdict['ResFile'] = ['depthcharge.tdt']
            try: open(fdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting fork.'); return self.endJob(fdict)
            ### ~ [2] ~ Add Fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setNum({'KillTime':time.time()})
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                fdict['PID'] = cpid
                self.printLog('\r#FORK','Forking seq as %s: %d remain; %.1f%% mem free' % (cpid,len(self.list['ToFork']),fdict['Mem']),log=self.getBool('LogFork'),screen=self.getBool('LogFork') or self.v() > 1)
                self.printLog('#FORK','%s seq: %s' % (cpid,fdict['seq']),log=self.getBool('LogFork'),screen=self.getBool('LogFork') or self.v() > 1)
            else:                   # child process
                self.baseFile(fdict['FID'])
                self.setInt({'Interactive':-1})
                self.log.info['ErrorLog'] = ''
                self.log.info['LogFile'] = fdict['Log']
                self.log.opt['Quiet'] = True       # When True, will not write to screen or log apart from errors.
                self.list['ResFile'] = fdict['ResFile']
                self.depthCharge(seqname,sequence)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('Forker.startFork error')
#########################################################################################################################
    def depthCharge(self,seqname,sequence):  ### Main DepthCharge method
        '''
        Main DepthCharge method:
        awk seqname and read in (start,end) tuples. sort().
        badpos = [] # list of bad positions for sequence
        span = [] # current list of 3' extensions
        pos = 1 # current position for consideration
        Until start>pos, pop from tuples. If end <=pos, dump, else add end to span. Sort span. Using mindepth set new pos else add pos to badpos and increment pos and clear span then repeat until end.
        Remove near termini or gaps then collapse and report and/or fragment.

        NOTE: v0.1.0 removed:
        endbuffer=INT   : Buffer size for sequence ends to avoid breakage [500]
        gapbuffer=INT   : Buffer size for sequence gaps to avoid breakage [100]
        '''
        sname = sequence
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sname = rje.split(seqname)[0]
            seqlen = len(sequence)
            reads = []      # list of (start,end) positions for seqname
            badpos = []     # list of bad positions for sequence
            pos = 1         # current position for consideration
            mindepth = self.getInt('MinDepth')
            ## ~ [1a] GapBuffer ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gaps = []
            if sname in self.db('gaps').index('seqname'):
                gaps = self.db('gaps').index('seqname')[sname]   # [(seqname,start,end)]
                gaps.sort()
            #endbuffer = max(0,self.getInt('GapBuffer'))
            # ends = range(1,endbuffer+1) + range(seqlen-endbuffer+1,seqlen+1)
            gapbuffer = []
            for gap in gaps:
                # gapbuffer += range(gap[1]-endbuffer,gap[1])
                # gapbuffer += range(gap[2]+1,gap[2]+endbuffer+1)
                #gapbuffer += range(gap[1]-endbuffer,gap[2]+endbuffer+1)
                gapbuffer += range(gap[1],gap[2]+1)
            #endbuffer = max(0,self.getInt('EndBuffer'))
            ### ~ [2] Read in (start,end) list for reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# grep seqname and use awk to read in (start,end) tuples -> sort().
            for pafline in os.popen("grep {0} {1} | awk '{{print $8,$9;}}'".format(sname,self.getStr('PAF'))).readlines():
                try:
                    [x,y] = rje.split(pafline)[:2]
                    [x,y] = [int(x),int(y)]
                    if self.getInt('MinSpan') > 0:
                        x += self.getInt('MinSpan')
                        y -= self.getInt('MinSpan')
                    if y > x:
                        reads.append((x,y))
                except: pass
            reads.sort()
            self.printLog('#READS','{0} reads parsed for {1}.'.format(rje.iLen(reads),sname))
            ### ~ [3] Charge! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# 1. Until start>pos, pop from tuples.
            #i# 2. If end <=pos, dump, else add end to span.
            #i# 3. Sort span. Using mindepth set new pos else add pos to badpos and increment pos and clear span then repeat until end.
            #i# 4. Remove near termini or gaps.
            span = []       # current list of 3' extensions
            while pos <= seqlen:
                while reads and reads[0][0] <= pos:
                    read = reads.pop(0)
                    if read[1] <= pos: continue
                    span.append(read[1])
                span.sort()
                if len(span) >= mindepth:
                    span = span[-mindepth:] # Keep the rest for the next position
                    pos = span.pop(0)
                else:
                    # if pos not in ends:
                    badpos.append(pos)
                    # Jump to next read
                    if reads: # NB reads[0][0] > pos
                        badpos += list(range(pos+1,reads[0][0]))
                    else:
                        badpos += list(range(pos+1,seqlen+1))
                    pos = badpos[-1] + 1
                    while span and span[0] <= pos: span.pop(0)
            px = len(badpos); self.progLog('#BAD','{0} bad positions for {1}.'.format(rje.iLen(badpos),sname))
            #i# 5. Collapse to report and/or fragment.
            #gaptuples = []
            #gappos = rje.listIntersect(badpos,gapbuffer)
            #badpos = rje.listDifference(badpos,gappos)    ### Returns the elements of list1 that are not found in list 2.
            badpos.sort()
            #gappos.sort()
            badtuples = []
            while badpos:
                i = j = badpos.pop(0)
                while badpos and badpos[0] in [j,j+1]: j = badpos.pop(0)
                badtuples.append((i,j))
            self.progLog('\r#BAD','{0} bad positions for {1} => {2} bad regions.'.format(rje.iStr(px),sname,rje.iLen(badtuples)))
            # while gappos:
            #     i = j = gappos.pop(0)
            #     while gappos and gappos[0] in [j,j+1]: j = gappos.pop(0)
            #     gaptuples.append((i,j))
            #i# 6. Save to *.depthcharge.tdt : seqname, start, end, type
            ddb = self.db().addEmptyTable('depthcharge',['seqname','start','end','type'],['seqname','type','start','end'],log=False)
            btype = 'all'; i = 1; j = seqlen; bx =0; ex = 0; gx = 0
            ddb.dict['Data'][(sname,i,j,btype)] = {'seqname':sname,'start':i,'end':j,'type':btype}
            for (i,j) in badtuples:
                #if j <= endbuffer or i >= (seqlen-endbuffer+1):
                if i == 1 or j == seqlen:   # End of sequence
                    btype = 'end'; ex += 1
                elif rje.listIntersect(range(i,j+1),gapbuffer):     # Overlaps a gap
                    btype = 'gap'; gx += 1
                else:
                    btype = 'bad'; bx += 1
                ddb.dict['Data'][(sname,i,j,btype)] = {'seqname':sname,'start':i,'end':j,'type':btype}
            self.printLog('\r#BAD','{0} bad positions for {1} => {2} bad regions; {3} gaps; {4} ends.'.format(rje.iStr(px),sname,rje.iStr(bx),rje.iStr(gx),rje.iStr(ex)))
            # btype = 'gap'
            # for (i,j) in gaptuples:
            #     ddb.dict['Data'][(sname,i,j,btype)] = {'seqname':sname,'start':i,'end':j,'type':btype}
            ddb.saveToFile()    #'{0}.{1}'.format(fdict['FID'],self.list['ResFile'][0]))

            return True
        except:
            self.errorLog('%s.depthCharge(%s) error' % (self.prog(),sname))
            return False
#########################################################################################################################
### End of SECTION II: DepthCharge Class                                                                                #
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
    try: DepthCharge(mainlog,cmd_list).run()

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
