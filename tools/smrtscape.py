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
Module:       SMRTSCAPE
Description:  SMRT Subread Coverage & Assembly Parameter Estimator
Version:      1.10.1
Last Edit:    26/05/16
Copyright (C) 2015  Richard J. Edwards - See source code for GNU License Notice

Function:
    SMRTSCAPE (SMRT Subread Coverage & Assembly Parameter Estimator) is tool in development as part of our PacBio
    sequencing projects for predicting and/or assessing the quantity and quality of useable data required/produced for
    HGAP3 de novo whole genome assembly. The current documentation is below. Some tutorials will be developed in the
    future - in the meantime, please get in touch if you want to use it and anything isn't clear.

    The main functions of `SMRTSCAPE` are:

    1. Estimate Genome Coverage and required numbers of SMRT cells given predicted read outputs. NOTE: Default settings
    for SMRT cell output are not reliable and you should speak to your sequencing provider for their up-to-date figures.

    2. Summarise the amount of sequence data obtained from one or more SMRT cells, including unique coverage (one read
    per ZMW).

    3. Calculate predicted coverage from subread data for difference length and quality cutoffs.

    4. Predict HGAP3 length and quality settings to achieve a given coverage and accuracy.

	SMRTSCAPE `coverage=T` mode can be run from the EdwardsLab server at:
	<http://www.slimsuite.unsw.edu.au/servers/pacbio.php>

Commandline:
    ### ~ General Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    genomesize=X    : Genome size (bp) [0]
    ### ~ Genome Coverage Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    coverage=T/F    : Whether to generate coverage report [False]
    avread=X        : Average read length (bp) [20000]
    smrtreads=X     : Average assemble output of a SMRT cell [50000]
    smrtunits=X     : Units for smrtreads=X (reads/Gb/Mb) [reads]
    errperbase=X    : Error-rate per base [0.14]
    maxcov=X        : Maximum X coverage to calculate [100]
    bysmrt=T/F      : Whether to output estimated  coverage by SMRT cell rather than X coverage [False]
    xnlist=LIST     : Additional columns giving % sites with coverage >= Xn [1+`minanchorx`->`targetxcov`+`minanchorx`]
    ### ~ SubRead Summary Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    summarise=T/F   : Generate subread summary statistics including ZMW summary data [False]
    seqin=FILE      : Subread sequence file for analysis [None]
    batch=FILELIST  : Batch input of multiple subread fasta files (wildcards allowed) if seqin=None []
    targetcov=X     : Target percentage coverage for final genome [99.999]
    targeterr=X     : Target errors per base for preassembly [1/genome size]
    calculate=T/F   : Calculate X coverage and target X coverage for given seed, anchor + RQ combinations [False]
    minanchorx=X    : Minimum X coverage for anchor subreads [6]
    minreadlen=X    : Absolute minimum read length for calculations (use minlen=X to affect summary also) [500]
    rq=X,Y          : Minimum (X) and maximum (Y) values for read quality cutoffs [0.8,0.9]
    rqstep=X        : Size of RQ jumps for calculation (min 0.001) [0.01]
    ### ~ Preassembly Fragmentation analysis Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    preassembly=FILE: Preassembly fasta file to assess/correct over-fragmentation (use seqin=FILE for subreads) [None]
    ### ~ Assembly Parameter Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    parameters=T/F  : Whether to output predicted "best" set of parameters [False]
    targetxcov=X    : Target 100% X Coverage for pre-assembly [3]
    xmargin=X       : "Safety margin" inflation of X coverage [1]
    mapefficiency=X : [Adv.] Efficiency of mapping anchor subreads onto seed reads for correction [1.0]
    xsteplen=X      : [Adv.] Size (bp) of increasing coverage steps for calculating required depths of coverage [1e5]
    parseparam=FILES: Parse parameter settings from 1+ assembly runs []
    paramlist=LIST  : List of parameters to retain for parseparam output (file or comma separated, blank=all) []
    predict=T/F     : Whether to add XCoverage prediction and efficiency estimation from parameters and subreads [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import math, os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_tree
import rje_dismatrix_V3 as rje_dismatrix
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 1.0.0 - Initial working version for server.
    # 1.1.0 - Added xnlist=LIST : Additional columns giving % sites with coverage >= Xn [10,25,50,100].
    # 1.2.0 - Added assessment -> now PAGSAT.
    # 1.3.0 - Added seed and anchor read coverage generator (calculate=T).
    # 1.3.1 - Deleted assessment function. (Now handled by PAGSAT.)
    # 1.4.0 - Added new coverage=T function that incorporates seed and anchor subreads.
    # 1.5.0 - Added parseparam=FILES with paramlist=LIST to parse restricted sets of parameters.
    # 1.6.0 - New SMRTSCAPE program building on PacBio v1.5.0. Added predict=T/F option.
    # 1.6.1 - Updated parameters=T to incorporate that the seed read counts as X=1.
    # 1.7.0 - Added *.summary.tdt output from subread summary analysis. Added minreadlen.
    # 1.8.0 - preassembly=FILE: Preassembly fasta file to assess/correct over-fragmentation (use seqin=FILE for subreads)
    # 1.9.0 - Updated empirical preassembly mapefficiency calculation.
    # 1.10.0 - Added batch processing of subread files.
    # 1.10.1 - Fixed bug in batch processing.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Add REST outputs to restSetup() and restOutputOrder()
    # [Y] : Add to SLiMSuite or SeqSuite.
    # [ ] : Improved error estimation.
    # [Y] : Option to do a per-SMRT cell analysis.
    # [ ] : Add costing.
    # [Y] : Add subread summary and seed read length estimation.
    # [ ] : Reading of genome size from refgenome file if genomesize=0?
    # [Y] : Update parameters to use calculate.
    # [?] : Implement rqmean=T/F      : Whether to use mean RQ instead of min RQ for parameters=T calculations [False]
    # [Y] : Add a MapEfficiency Setting.
    # [Y] : Add predict=T setting to use specific settings and data to predict performance.
    # [Y] : Should `targetcov` be squared prior to analysis? (Needs seed AND anchor?)
    # [?] : Make sure that auto-generated seed lengths are recognised and calculated. Added but not convinced it's right.
    # [Y] : Add other key parameters to predict.tdt. (Check right udb/zdb being used for coverage and add docs.)
    # [ ] : Add option to read/generate ZMW, Unique and RQ data with separate basefile.
    # [ ] : Add an option to use restricted sets of SMRT cells.
    # [Y] : Consider reducing default xsteplen to 1e5?
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SMRTSCAPE', '1.10.1', 'May 2016', '2015')
    description = 'SMRT Subread Coverage & Assembly Parameter Estimator'
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
### SECTION II: SMRTSCAPE Class                                                                                         #
#########################################################################################################################
class SMRTSCAPE(rje_obj.RJE_Object):
    '''
    SMRTSCAPE Class. Author: Rich Edwards (2015).

    Str:str
    - Preassembly=FILE: Preassembly fasta file to assess/correct over-fragmentation (use seqin=FILE for subreads) [None]
    - SeqIn=FILE      : Subreads sequence input file
    - SMRTUnits=X     : Units for smrtreads=X (reads/Gb/Mb) [reads]

    Bool:boolean
    - BySMRT=T/F      : Whether to output estimated by SMRT cell rather than X coverage [False]
    - Calculate=T/F   : Calculate X coverage and target X coverage for given seed, anchor + RQ combinations [False]
    - Coverage=T/F    : Whether to generate coverage report [False]
    - FastSeedX=T/F   : Whether to use Fast SeedX estimator rather than more accurate "slow" method [True]
    - Parameters=T/F  : Whether to output
    - Predict=T/F     : Whether to add XCoverage prediction and efficiency estimation from parameters and subreads [False]
    - RQMean=T/F      : Whether to use mean RQ instead of min RQ for calculations [False]
    - Summarise=T/F   : Generate subread summary statistics including ZMW summary data [False]

    Int:integer
    - MaxCov=X        : Maximmum X coverage to calculate [100]
    - MinAnchorX=X    : Minimum X coverage for anchor subreads [10]
    - MinReadLen=X    : Absolute minimum read length for calculations (use minlen=X to affect summary also) [500]

    Num:float
    - AvRead=X        : Average read length (bp) [20000]
    - ErrPerBase=X    : Error-rate per base [0.14]
    - GenomeSize=X    : Genome size (bp) [4e9]
    - MapEfficiency=X : Efficiency of mapping anchor subreads onto seed reads for correction [1.0]
    - MaxRQ=X,Y          : Minimum (X) and maximum (Y) values for read quality cutoffs [0.8,0.9]
    - MinRQ=X,Y          : Minimum (X) and maximum (Y) values for read quality cutoffs [0.8,0.9]
    - RQStep=X        : Size of RQ jumps for calculation (min 0.001) [0.01]
    - SMRTReads=X     : Average number of reads of a SMRT cell [50000]
    - TargetCov=X     : Target coverage for final genome [99.999]
    - TargetErr=X     : Target errors per base for preassembly [1e-6]
    - TargetXCov=X    : Target 100% X Coverage for pre-assembly [3]
    - XMargin=X       : "Safety margin" inflation of X coverage [1]
    - XStepLen=X      : Size (bp) of increasing coverage steps for calculating required depths of coverage [1e6]

    File:file handles with matching str filenames
    
    List:list
    - Accuracy : list of %accuracy for each xcoverage level
    - Batch=FILELIST  : Batch input of multiple subread fasta files (wildcards allowed) []
    - ParamList=LIST  : List of parameters to retain for parseparam output (file or comma separated, blank=all) []
    - ParseParam=FILES: Parse parameter settings from 1+ assembly runs []
    - TargetXDepth= []    # List of target x (index) and required depth (value)
    - XCovLimits : List of total amount of sequence for each X coverage at TargetCov %coverage []
    - XnList=LIST     : Additional columns giving % sites with coverage >= Xn [10,25,50,100]

    Dict:dictionary
    - PercXDepth : {XCov:[%coverage at each Xdepth]}

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Preassembly','SeqIn','SMRTUnits']
        self.boollist = ['BySMRT','Calculate','Coverage','FastSeedX','Parameters','Predict','RQMean','Summarise']
        self.intlist = ['MaxCov','MinAnchorX','MinReadLen','TargetXCov','XMargin']
        self.numlist = ['AvRead','ErrPerBase','GenomeSize','MapEfficiency','MaxRQ','MinRQ','RQStep','SMRTReads','TargetCov','TargetErr','XStepLen']
        self.filelist = []
        self.listlist = ['Accuracy','Batch','ParamList','ParseParam','TargetXDepth','XCovLimits','XnList']
        self.dictlist = ['PercXDepth']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'SMRTUnits':'reads'})
        self.setBool({'BySMRT':False,'Calculate':False,'Coverage':False,'FastSeedX':True,'Parameters':False,'Summarise':False})
        self.setInt({'MaxCov':100,'MinAnchorX':6,'MinReadLen':500,'TargetXCov':3,'XMargin':1})
        self.setNum({'AvRead':20000,'ErrPerBase':0.14,'MapEfficiency':1.0,'MaxRQ':0.9,'MinRQ':0.8,'RQStep':0.01,
                     'GenomeSize':0,'SMRTReads':50000,'TargetCov':99.999,'TargetErr':-1,'XStepLen':1e5})
        self.list['XnList'] = []    #10,25,50,100]
        self.dict['PercXDepth'] = {0:[1.0]}
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
                self._cmdReadList(cmd,'str',['SMRTUnits'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Preassembly','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['BySMRT','Calculate','Coverage','FastSeedX','Parameters','Predict','RQMean','Summarise'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MaxCov','MinAnchorX','MinReadLen','TargetXCov','XMargin'])   # Integers
                self._cmdReadList(cmd,'float',['AvRead','ErrPerBase','GenomeSize','MapEfficiency','RQStep','SMRTReads','TargetErr','XStepLen']) # Floats
                self._cmdReadList(cmd,'perc',['TargetCov'])
                self._cmdRead(cmd,'fmax','MaxRQ','rq')   # Integer value part of min,max command
                self._cmdRead(cmd,'fmin','MinRQ','rq')   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ParamList','XnList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['Batch','ParseParam']) # List of files using wildcards and glob
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
            if self.getStrLC('Preassembly'): return self.preassemblyFrag()
            if self.getBool('Predict'): return self.predict()       # Includes parseParam and summarise if required.
            elif self.list['ParseParam']: return self.parseParam()  # Standalone function.
            if self.getBool('Parameters'): self.parameters()    # Includes calculate if required.
            elif self.getBool('Calculate'): self.calculate()    # Includes summarise if required.
            elif self.getBool('Summarise'): self.summarise()
            if self.getBool('Coverage'): self.coverage()        # Does not require summarise but will use if present.
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.baseFile(return_none=None) and self.getStrLC('SeqIn'):
                basename = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
                if basename.endswith('.subreads'): basename = os.path.splitext(basename)[0]
                self.baseFile(basename)
            elif not self.baseFile(return_none=None) and self.getStrLC('Preassembly'):
                basename = rje.baseFile(self.getStr('Preassembly'),strip_path=True)
                if basename.endswith('.preassembly'): basename = os.path.splitext(basename)[0]
                self.baseFile(basename)
            elif not self.baseFile(return_none=None): self.baseFile('smrtscape')
            if not self.list['ParseParam'] and not self.getBool('Predict'):
                while not self.getInt('GenomeSize'):
                    if self.i() >= 0: self.setNum({'GenomeSize':int(float(rje.choice('Enter Genome Size (bp)',confirm=True)))})
                    else: raise ValueError('Need to set genomesize=X')
                self.printLog('#GSIZE','Genome Size: %s bp' % rje.iStr(self.getInt('GenomeSize')))
                if self.getNum('TargetErr') < 0: self.setNum({'TargetErr':1.0/self.getInt('GenomeSize')})
                self.printLog('#PCERR','Target Error Rate = %s errors per base.' % rje.expectString(self.getNum('TargetErr')))
                if self.getStrLC('REST') and not self.basefile(return_none=''): self.basefile('smrtscape')
                self.list['Accuracy'] = [0,1.0 - self.getNum('ErrPerBase')]
            elif not self.baseFile(return_none=None): self.baseFile('smrtscape')
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            self.db().basefile(self.basefile())
            self.printLog('#BASE','Basename set for input/output: %s.*' % self.basefile())
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
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
    ### <3> ### Summarise Subreads Methods                                                                              #
#########################################################################################################################
    def batchSummarise(self):    ### Generate subread summary.
        '''Generate subread summary.'''
        try:### ~ [0] Setup SeqList and basic stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()  # Database object
            seqbatch = []   # List of SeqList objects
            self.printLog('#BATCH','%s sequence files to process.' % rje.iLen(self.list['Batch']))
            for seqfile in self.list['Batch']:
                seqcmd = self.cmd_list + ['seqmode=file','autoload=T','summarise=F','seqin=%s' % seqfile,'autofilter=F']
                seqbatch.append(rje_seqlist.SeqList(self.log,seqcmd))
            self.printLog('#BATCH','%s sequence files to summarise.' % rje.iLen(seqbatch))
            if not seqbatch: return IOError('No batch input fasta files found!')
            cells = []  # List of SMRT Cell identifiers. Add dictionary converter for names. Or ask?
            zdb = self.db().addEmptyTable('zmw',['SMRT','ZMW','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW','RN'])
            ### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db().addEmptyTable('summary',['SMRT','SeqNum','TotLength','MinLength','MaxLength','MeanLength','MedLength','N50Length','XCoverage'],['SMRT'])
            for seqlist in seqbatch:
                seqdata = seqlist.summarise()   # This will summarise read lengths etc.
                basename = rje.baseFile(seqlist.getStr('SeqIn'),strip_path=True)
                if basename.endswith('.subreads'): basename = os.path.splitext(basename)[0]
                seqdata['SMRT'] = basename
                seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
                sdb.addEntry(seqdata)
                # SMRT = SMRT cell identifier
                # ZMW = ZMW identifier
                # RN = subread number for that ZMW
                # Len = length of actual sequence read in
                # Pos = Position information from header
                # RQ = RQ Value
                # Partition out ZMWs and make histograms of lengths and numbers
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/9/0_3967 RQ=0.784
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/11/0_20195 RQ=0.868
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/12/0_18976 RQ=0.823
                #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/14/5776_27852 RQ=0.807
                #...
                #>m150625_043200_42272_c100792502550000001823157609091583_s1_p0/163464/0_1712 RQ=0.854
                prevzmw = None; rn = 0; sx = 0.0; stot = seqlist.seqNum()
                for seq in seqlist.seqs():
                    self.progLog('\r#SUB','Processing subreads: %.2f%%' % (sx/stot)); sx += 100.0
                    (name,sequence) = seqlist.getSeq(seq)
                    [smrt,zmw,pos,rq] = string.split(string.replace(name,'/',' '))
                    if smrt not in cells: cells.append(smrt)
                    smrt = cells.index(smrt)
                    if zmw != prevzmw: prevzmw = zmw; rn = 0
                    rn += 1
                    rq = rje.matchExp('RQ=(\S+)',rq)[0]
                    zdb.addEntry({'SMRT':smrt,'ZMW':zmw,'RN':rn,'Len':len(sequence),'Pos':pos,'RQ':rq,'Seq':seq})
                self.printLog('\r#SUB','Processed %s subreads (-> %d SMRT cells)' % (rje.iStr(stot),len(cells)))
            zdb.dataFormat({'RN':'int','Len':'int','RQ':'float'})
            zdb.saveToFile()
            zdb.index('RN')
            for smrt in cells:
                self.printLog('#~~#','# ~~~~~~~~~~~ SEQUENCE SUMMARY FOR %s ~~~~~~~~~~~ #' % smrt)
                self.printLog('#SMRT','SMRT %d = %s' % (cells.index(smrt),smrt))
                seqdata = self.summariseSeqLen(zdb.indexDataList('SMRT',cells.index(smrt),'Len',sortunique=False))
                seqdata['SMRT'] = smrt
                sdb.addEntry(seqdata)

            basename = self.baseFile()
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ COMBINED PACBIO SUBREAD SUMMARY FOR %s ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % basename)
            seqlen = zdb.dataList(zdb.entries(),'Len',sortunique=False)
            seqdata = self.summariseSeqLen(seqlen)
            sdb.addEntry(seqdata)
            ### >>> Rest as summarise() >>>

            # Make a ZMW unique table by collapsing on ZMW and keeping max length or quality where tied
            self.printLog('#~~#','# ~~~~~~~~~~~ PACBIO Unique ZMW subreads (%s) ~~~~~~~~~~~ #' % basename)
            udb = db.copyTable(zdb,'unique')
            #udb.dropField('Pos')
            udb.compress(['SMRT','ZMW'],default='max',best=['Len','RQ','RN'])
            udb.saveToFile()
            #udb.indexReport('RN','Best ZMW reads from pass (RN)')
            udb.index('RN')
            zrn = rje.sortKeys(zdb.index('RN')); rmax = max(zrn)
            for rn in zrn:
                zn = rje.iLen(zdb.index('RN')[rn])
                ptxt = 'Read %d longest in ' % rn
                if rn in udb.index('RN'):
                    un = rje.iLen(udb.index('RN')[rn])
                    if rn < rmax: px = len(zdb.index('RN')[rn+1])
                    else: px = 0
                    ptxt += '%s of %s ZMW; %s ZMW with %d+ passes.' % (un,zn,rje.iStr(px),rn+1)
                else: ptxt += '0 ZMW.'; continue
                self.printLog('#RN',ptxt)
                #if rn in udb.index('RN'): self.printLog('#RN','Read %d longest in %s of %s ZMW with %d+ passes.' % (rn,rje.iLen(udb.index('RN')[rn]),rje.iLen(zdb.index('RN')[rn]),rn))
                #else: self.printLog('#RN','Read %d longest in 0 of %s ZMW with %d+ passes.' % (rn,rje.iLen(zdb.index('RN')[rn]),rn))
            bestseq = udb.dataList(udb.entries(),'Seq',sortunique=False,empties=False)
            self.printLog('#~~#','# ~~~~~~~~~~~ SEQUENCE SUMMARY FOR %s UNIQUE ~~~~~~~~~~~ #' % basename)
            seqlen = udb.dataList(udb.entries(),'Len',sortunique=False)
            seqdata = self.summariseSeqLen(seqlen)
            seqdata['SMRT'] = '%s.unique' % basename
            sdb.addEntry(seqdata)
            for smrt in cells:
                self.printLog('#~~#','# ~~~~~~~~~~~ SEQUENCE SUMMARY FOR %s.unique ~~~~~~~~~~~ #' % smrt)
                seqdata = self.summariseSeqLen(udb.indexDataList('SMRT',cells.index(smrt),'Len',sortunique=False))
                seqdata['SMRT'] = '%s.unique' % smrt
                sdb.addEntry(seqdata)
            sdb.saveToFile()

            # Output number and percentage of subreads and longest reads at each rq. Fields: rq, subreads, bestreads
            self.printLog('#~~#','# ~~~~~~~~~~~ PACBIO ZMW subread RQ (%s) ~~~~~~~~~~~ #' % basename)
            if self.getNum('TargetErr') <= 0:
                self.setNum({'TargetErr':1.0/self.getNum('GenomeSize')})
                self.printLog('#RQERR','Set target RQ error per base to 1/genome size = %s' % rje.expectString(self.getNum('TargetErr')))
            # Report on read quality and optionally filter?
            rqz = {}; rqzlen = {}
            sumrq = 0   # Sum of RQ * subreads
            for rq in zdb.index('RQ'):
                rqz[rq] = len(zdb.index('RQ')[rq])
                rqzlen[rq] = sum(zdb.indexDataList('RQ',rq,'Len',sortunique=False))
                sumrq += (rq * rqzlen[rq])
            rqzlentot = sum(rqzlen.values())
            self.printLog('#TOTAL','Total length = %.2fMb' % (rqzlentot/1e6))
            maxx = float(rqzlentot) / self.getNum('GenomeSize')
            self.printLog('#MAXX','Max. XCoverage = %.1f' % maxx)
            rqzfreq = rje.dictFreq(rqz,total=False,newdict=True)

            rqu = {}; rqulen = {}
            for rq in udb.index('RQ'):
                rqu[rq] = len(udb.index('RQ')[rq])
                rqulen[rq] = sum(udb.indexDataList('RQ',rq,'Len',sortunique=False))
            rqulentot = sum(rqulen.values())
            rqufreq = rje.dictFreq(rqu,total=False,newdict=True)

            qdb = db.addEmptyTable('rq',['RQ','xerr','subreads','unique','f.subreads','f.unique','cum.subreads','cum.unique','x.subreads','x.unique','MeanRQ','Mean.XErr'],['RQ'])
            for rq in rje.sortKeys(rqz):
                meanrq = sumrq / rqzlentot
                self.progLog('\r#RQ','Processing RQ=%s (Mean=%.3f) ' % (rq,meanrq))
                if rq > 1: raise ValueError('RQ = %s' % rq)
                x = 1
                while ((1-rq) ** x) > self.getNum('TargetErr'): x += 1
                meanx = 1
                while ((1-meanrq) ** meanx) > self.getNum('TargetErr'): meanx += 1
                rentry = {'RQ':rq,'xerr':x,'subreads':rqz[rq],'unique':0,'f.subreads':rqzfreq[rq],'f.unique':0,
                          'cum.subreads':rqzlentot,'x.subreads':rqzlentot/self.getNum('GenomeSize'),
                          'MeanRQ':meanrq,'Mean.XErr':meanx}
                rqzlentot -= rqzlen[rq]
                sumrq -= (rq * rqzlen[rq])
                if rq in rqu:
                    rentry['unique'] = rqu[rq]; rentry['f.unique'] = rqufreq[rq]
                    rentry['cum.unique'] = rqulentot; rentry['x.unique'] = rqulentot/self.getNum('GenomeSize')
                    rqulentot -= rqulen[rq]
                #self.printLog('#RQ','RQ=%s: %s (%.2f%%) subreads; %s (%.2f%%) unique.' % (rq,rqz[rq],100.0*rqzfreq[rq],rentry['unique'],100.0*rentry['f.unique']))
                qdb.addEntry(rentry)
            self.printLog('\r#RQ','Processing of RQ<=%s complete.' % rq)
            qdb.saveToFile()

            #i# Old notes and contour data now in self.contours(). Replaced by calculate()?
        except: self.errorLog('%s.batchSummarise error' % self.prog())
#########################################################################################################################
    def summariseSeqLen(self,seqlen):   ### Generate summary data from list of sequence lengths
        '''
        Generate summary data from list of sequence lengths.
        @param seqlen: list of sequence lists
        @return: dictionary of summary data
        '''
        try:### ~ [0] Setup SeqList and basic stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlen.sort()
            seqnum = len(seqlen)
            sumlen = sum(seqlen)
            seqdata = {'SMRT':'Total','SeqNum':seqnum,'TotLength':sumlen,
                     'MinLength':min(seqlen),'MaxLength':max(seqlen)}
            seqdata['MeanLength'] = seqdata['TotLength']/float(seqdata['SeqNum'])
            self.printLog('#SUM','Total number of sequences: %s' % rje.iLen(seqlen))
            self.printLog('#SUM','Total length of sequences: %s' % rje.iStr(sumlen))
            self.printLog('#SUM','Min. length of sequences: %s' % rje.iStr(seqlen[0]))
            self.printLog('#SUM','Max. length of sequences: %s' % rje.iStr(seqlen[-1]))
            # Mean & Median sequence lengths
            meanlen = float(sumlen)/len(seqlen)
            meansplit = string.split('%.2f' % meanlen,'.')
            self.printLog('#SUM','Mean length of sequences: %s.%s' % (rje.iStr(meansplit[0]),meansplit[1]))
            if rje.isOdd(len(seqlen)): median = seqlen[len(seqlen)/2]
            else: median = sum(seqlen[len(seqlen)/2:][:2]) / 2.0
            self.printLog('#SUM','Median length of sequences: %s' % (rje.iStr(median)))
            seqdata['MedLength'] = median
            ## N50 calculation
            n50len = seqdata['TotLength'] / 2.0
            n50 = seqlen[0:]
            while n50len > 0 and n50: n50len -= n50.pop(-1)
            if n50:
                self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(n50[-1]))
                seqdata['N50Length'] = n50[-1]
            else:
                self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(seqlen[-1]))
                seqdata['N50Length'] = seqlen[-1]
            seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
            return seqdata
        except: self.errorLog('%s.summarise error' % self.prog())
#########################################################################################################################
    def summarise(self):    ### Generate subread summary.
        '''Generate subread summary.'''
        try:### ~ [0] Setup SeqList and basic stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['Batch'] and not self.getStrLC('SeqIn'): return self.batchSummarise()
            db = self.db()
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','summarise=F'])   # This will summarise read lengths etc.
            basename = rje.baseFile(seqlist.getStr('SeqIn'),strip_path=True)
            # >>> Now done in setup()
            #if self.getStrLC('Basefile'): basename = self.getStr('Basefile')
            #else:
            #    basename = rje.baseFile(seqlist.getStr('SeqIn'),strip_path=True)
            #    if basename.endswith('.subreads'): basename = os.path.splitext(basename)[0]
            #db.baseFile(basename)
            # <<<
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ PACBIO SUBREAD SUMMARY FOR %s ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % basename)
            ### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sdb = self.db().addEmptyTable('summary',['SMRT','SeqNum','TotLength','MinLength','MaxLength','MeanLength','MedLength','N50Length','XCoverage'],['SMRT'])
            #!# Add XCoverage #!#
            #!# For cutoffs.tdt, need optimisation when insufficient for full coverage - like calculate=T.
            seqdata = seqlist.summarise()
            seqdata['SMRT'] = basename
            seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
            sdb.addEntry(seqdata)
            # SMRT = SMRT cell identifier
            # ZMW = ZMW identifier
            # RN = subread number for that ZMW
            # Len = length of actual sequence read in
            # Pos = Position information from header
            # RQ = RQ Value
            cells = []  # List of SMRT Cell identifiers. Add dictionary converter for names. Or ask?
            zdb = self.db().addEmptyTable('zmw',['SMRT','ZMW','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW','RN'])
            # Partition out ZMWs and make histograms of lengths and numbers
            #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/9/0_3967 RQ=0.784
            #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/11/0_20195 RQ=0.868
            #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/12/0_18976 RQ=0.823
            #>m150625_001530_42272_c100792502550000001823157609091582_s1_p0/14/5776_27852 RQ=0.807
            #...
            #>m150625_043200_42272_c100792502550000001823157609091583_s1_p0/163464/0_1712 RQ=0.854
            prevzmw = None; rn = 0; sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#SUB','Processing subreads: %.2f%%' % (sx/stot)); sx += 100.0
                (name,sequence) = seqlist.getSeq(seq)
                [smrt,zmw,pos,rq] = string.split(string.replace(name,'/',' '))
                if smrt not in cells: cells.append(smrt)
                smrt = cells.index(smrt)
                if zmw != prevzmw: prevzmw = zmw; rn = 0
                rn += 1
                rq = rje.matchExp('RQ=(\S+)',rq)[0]
                zdb.addEntry({'SMRT':smrt,'ZMW':zmw,'RN':rn,'Len':len(sequence),'Pos':pos,'RQ':rq,'Seq':seq})
            self.printLog('\r#SUB','Processed %s subreads (%d SMRT cells)' % (rje.iStr(stot),len(cells)))
            for smrt in cells:
                self.printLog('#SMRT','SMRT %d = %s' % (cells.index(smrt),smrt))
                seqdata = seqlist.summarise(zdb.indexDataList('SMRT',cells.index(smrt),'Seq',sortunique=False),smrt)
                seqdata['SMRT'] = smrt
                seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
                sdb.addEntry(seqdata)
            zdb.saveToFile()
            zdb.dataFormat({'RN':'int','Len':'int','RQ':'float'})
            zdb.index('RN')

            # Make a ZMW unique table by collapsing on ZMW and keeping max length or quality where tied
            self.printLog('#~~#','# ~~~~~~~~~~~ PACBIO Unique ZMW subreads (%s) ~~~~~~~~~~~ #' % basename)
            udb = db.copyTable(zdb,'unique')
            #udb.dropField('Pos')
            udb.compress(['SMRT','ZMW'],default='max',best=['Len','RQ','RN'])
            udb.saveToFile()
            #udb.indexReport('RN','Best ZMW reads from pass (RN)')
            udb.index('RN')
            zrn = rje.sortKeys(zdb.index('RN')); rmax = max(zrn)
            for rn in zrn:
                zn = rje.iLen(zdb.index('RN')[rn])
                ptxt = 'Read %d longest in ' % rn
                if rn in udb.index('RN'):
                    un = rje.iLen(udb.index('RN')[rn])
                    if rn < rmax: px = len(zdb.index('RN')[rn+1])
                    else: px = 0
                    ptxt += '%s of %s ZMW; %s ZMW with %d+ passes.' % (un,zn,rje.iStr(px),rn+1)
                else: ptxt += '0 ZMW.'; continue
                self.printLog('#RN',ptxt)
                #if rn in udb.index('RN'): self.printLog('#RN','Read %d longest in %s of %s ZMW with %d+ passes.' % (rn,rje.iLen(udb.index('RN')[rn]),rje.iLen(zdb.index('RN')[rn]),rn))
                #else: self.printLog('#RN','Read %d longest in 0 of %s ZMW with %d+ passes.' % (rn,rje.iLen(zdb.index('RN')[rn]),rn))
            bestseq = udb.dataList(udb.entries(),'Seq',sortunique=False,empties=False)
            seqdata = seqlist.summarise(bestseq,basename='%s UNIQUE' % basename)
            seqdata['SMRT'] = '%s.unique' % basename
            seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
            sdb.addEntry(seqdata)
            for smrt in cells:
                seqdata = seqlist.summarise(udb.indexDataList('SMRT',cells.index(smrt),'Seq',sortunique=False),'%s.unique' % smrt)
                seqdata['SMRT'] = '%s.unique' % smrt
                seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
                sdb.addEntry(seqdata)
            sdb.saveToFile()

            # Output number and percentage of subreads and longest reads at each rq. Fields: rq, subreads, bestreads
            self.printLog('#~~#','# ~~~~~~~~~~~ PACBIO ZMW subread RQ (%s) ~~~~~~~~~~~ #' % basename)
            if self.getNum('TargetErr') <= 0:
                self.setNum({'TargetErr':1.0/self.getNum('GenomeSize')})
                self.printLog('#RQERR','Set target RQ error per base to 1/genome size = %s' % rje.expectString(self.getNum('TargetErr')))
            # Report on read quality and optionally filter?
            rqz = {}; rqzlen = {}
            sumrq = 0   # Sum of RQ * subreads
            for rq in zdb.index('RQ'):
                rqz[rq] = len(zdb.index('RQ')[rq])
                rqzlen[rq] = sum(zdb.indexDataList('RQ',rq,'Len',sortunique=False))
                sumrq += (rq * rqzlen[rq])
            rqzlentot = sum(rqzlen.values())
            self.printLog('#TOTAL','Total length = %.2fMb' % (rqzlentot/1e6))
            maxx = float(rqzlentot) / self.getNum('GenomeSize')
            self.printLog('#MAXX','Max. XCoverage = %.1f' % maxx)
            rqzfreq = rje.dictFreq(rqz,total=False,newdict=True)

            rqu = {}; rqulen = {}
            for rq in udb.index('RQ'):
                rqu[rq] = len(udb.index('RQ')[rq])
                rqulen[rq] = sum(udb.indexDataList('RQ',rq,'Len',sortunique=False))
            rqulentot = sum(rqulen.values())
            rqufreq = rje.dictFreq(rqu,total=False,newdict=True)

            qdb = db.addEmptyTable('rq',['RQ','xerr','subreads','unique','f.subreads','f.unique','cum.subreads','cum.unique','x.subreads','x.unique','MeanRQ','Mean.XErr'],['RQ'])
            for rq in rje.sortKeys(rqz):
                meanrq = sumrq / rqzlentot
                self.progLog('\r#RQ','Processing RQ=%s (Mean=%.3f) ' % (rq,meanrq))
                if rq > 1: raise ValueError('RQ = %s' % rq)
                x = 1
                while ((1-rq) ** x) > self.getNum('TargetErr'): x += 1
                meanx = 1
                while ((1-meanrq) ** meanx) > self.getNum('TargetErr'): meanx += 1
                rentry = {'RQ':rq,'xerr':x,'subreads':rqz[rq],'unique':0,'f.subreads':rqzfreq[rq],'f.unique':0,
                          'cum.subreads':rqzlentot,'x.subreads':rqzlentot/self.getNum('GenomeSize'),
                          'MeanRQ':meanrq,'Mean.XErr':meanx}
                rqzlentot -= rqzlen[rq]
                sumrq -= (rq * rqzlen[rq])
                if rq in rqu:
                    rentry['unique'] = rqu[rq]; rentry['f.unique'] = rqufreq[rq]
                    rentry['cum.unique'] = rqulentot; rentry['x.unique'] = rqulentot/self.getNum('GenomeSize')
                    rqulentot -= rqulen[rq]
                #self.printLog('#RQ','RQ=%s: %s (%.2f%%) subreads; %s (%.2f%%) unique.' % (rq,rqz[rq],100.0*rqzfreq[rq],rentry['unique'],100.0*rentry['f.unique']))
                qdb.addEntry(rentry)
            self.printLog('\r#RQ','Processing of RQ<=%s complete.' % rq)
            qdb.saveToFile()

            #i# Old notes and contour data now in self.contours(). Replaced by calculate()?
        except: self.errorLog('%s.summarise error' % self.prog())
#########################################################################################################################
    def zdb(self):  ### Returns ZMW Table
        '''
        Returns ZMW Table.
        ['SMRT','ZMW','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW','RN']
        '''
        zdb = self.db('zmw',add=False)
        if zdb: return zdb
        zdb = self.db('zmw',add=True,mainkeys=['SMRT','ZMW','RN'])
        if not zdb: # Try based on SeqIn
            self.db().baseFile(rje.baseFile(self.getStr('SeqIn')))
            zdb = self.db('zmw',add=True,mainkeys=['SMRT','ZMW','RN'])
            self.db().baseFile(self.baseFile())
            if zdb: self.printLog('#ZMW','Could not find %s.zmw.tdt: used %s.zmw.tdt.' % (self.basefile(),rje.baseFile(self.getStr('SeqIn'))))
        if not zdb: self.summarise(); zdb = self.db('zmw')
        if not zdb: raise IOError('Cannot find *.zmw.tdt table!')
        if not zdb.formatted(): zdb.dataFormat({'RN':'int','Len':'int','RQ':'float'})
        zdb.index('RN')
        return zdb
#########################################################################################################################
    def udb(self,minrq=0.0,force=False):  ### Returns Unique ZMW Table
        '''
        Returns Unique ZMW Table.
        >> minrq:float [0] = Min. RQ cutoff to apply prior unique read generation.
        >> force:bool [False] = Whether to force usage of ZMW table.
        ['SMRT','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW']
        '''
        udb = self.db('unique',add=False)
        if force and udb: self.db().deleteTable(udb); udb = None
        else:
            if not udb:
                udb = self.db('unique',add=True,mainkeys=['SMRT','ZMW'])
        if not udb: # Try based on SeqIn
            udb = self.db().copyTable(self.zdb(),'unique')
            udb.dataFormat({'Len':'int','RQ':'float'})
            if minrq and min(udb.indexKeys('RQ')) < minrq: udb.dropEntries(['RQ<%f' % minrq])
            udb.compress(['SMRT','ZMW'],default='max',best=['Len','RQ','RN'])
            udb.index('RN')
        else:
            if not udb.formatted(): udb.dataFormat({'Len':'int','RQ':'float'})
            if minrq and  min(udb.indexKeys('RQ')) < minrq: udb.dropEntries(['RQ<%f' % minrq])
        return udb
#########################################################################################################################
    def rqdb(self):  ### Returns RQ Table from summarise
        '''
        Returns RQ Table from summarise.
        ['RQ','xerr','subreads','unique','f.subreads','f.unique','cum.subreads','cum.unique','x.subreads','x.unique','MeanRQ','Mean.XErr'],['RQ']
        '''
        rqdb = self.db('rq',add=False)
        if rqdb: return rqdb
        rqdb = self.db('rq',add=True,mainkeys=['RQ'])
        if not rqdb: # Try based on SeqIn
            self.db().baseFile(rje.baseFile(self.getStr('SeqIn')))
            rqdb = self.db('rq',add=True,mainkeys=['RQ'])
            self.db().baseFile(self.baseFile())
            if rqdb: self.printLog('#RQ','Could not find %s.rq.tdt: used %s.rq.tdt.' % (self.basefile(),rje.baseFile(self.getStr('SeqIn'))))
        if not rqdb: self.summarise(); rqdb = self.db('rq')
        if not rqdb: raise IOError('Cannot find *.rq.tdt table!')
        if not rqdb.formatted(): rqdb.dataFormat({'RQ':'int'}) # add formatting as required
        return rqdb
#########################################################################################################################
    def xCovLimits(self): ### Generates/returns XCovLimits list of summed read lengths generating whole genome X Coverage
        '''Generates/returns XCovLimits list of summed read lengths generating whole genome X Coverage.'''
        # This is basically targetXDepth but without the conversion to Xdepth
        if self.list['XCovLimits']: return self.list['XCovLimits']
        zdb = self.zdb()
        xlen = steplen = self.getNum('XStepLen')
        plimit =  math.sqrt(self.getPerc('TargetCov'))      # Needs square root for anchor x seed?
        self.list['XCovLimits'] = xcovlimits = [0]          # Target summed sequence lengths for each index X coverage
        self.progLog('\r#COV','Prepping seq len...')
        maxlen = sum(zdb.dataList(zdb.entries(),'Len',sortunique=False,empties=True))
        # Setup targets for each XCoverage
        while xlen <= maxlen:
            self.progLog('\r#COV','Setting up XCovLimits: %dX (%s bp)' % (len(xcovlimits),rje.iStr(xlen)))
            if xlen/self.getNum('GenomeSize') < len(xcovlimits): xlen += steplen; continue
            while rje.poisson(len(xcovlimits),xlen/self.getNum('GenomeSize'),exact=False,callobj=self) >= plimit:
                xcovlimits.append(xlen)
            xlen += steplen
        self.printLog('\r#COV','Set up XCovLimits complete: %dX (%s bp)' % (len(xcovlimits),rje.iStr(maxlen)))
        return self.list['XCovLimits']
#########################################################################################################################
    def pcXDepth(self,xcov,xdepth=-1):  ### (Calculates and) returns the proportion of coverage at Xdepth for given Xcov.
        '''
        (Calculates and) returns the proportion of coverage at Xdepth for given Xcov.
        >> xcov:int = the overall X coverage of the genome.
        >> xdepth:int = the specific X depth for which to return proportional coverage. Returns list if < 0.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xcov < 0: raise ValueError('Cannot have XCoverage < 0!')
            xdepth = int(xdepth)
            xdict = self.dict['PercXDepth']
            ### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xcov not in xdict:
                xdict[xcov] = []
                bases = int(self.getNum('GenomeSize'))
                while bases > 1:
                    xdict[xcov].append(rje.logPoisson(len(xdict[xcov]),xcov,exact=True,callobj=self))
                    bases -= self.getNum('GenomeSize') * xdict[xcov][-1]
            ### ~ [2] Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if xdepth < 0: return xdict[xcov]
            elif xdepth > len(xdict[xcov]): return 0.0
            else: return xdict[xcov][xdepth]
        except: self.errorLog('%s.pcXDepth error' % self.prog()); raise
#########################################################################################################################
    ### <4> ### PacBio Coverage Cutoff Calculation Methods                                                              #
#########################################################################################################################
    def calculate(self):    ### Calculate cutoffs from subread summary data.
        '''Calculate cutoffs from subread summary data.'''
        try:### ~ [0] Setup Objects for calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            efficiency = self.getNum('MapEfficiency')
            db = self.db()
            zdb = self.zdb(); udb = self.udb()
            zdb.dropEntries(['Len<%d' % self.getInt('MinReadLen')],logtxt='Min. read length filter')
            udb.dropEntries(['Len<%d' % self.getInt('MinReadLen')],logtxt='Min. read length filter')
            xcovlimits = self.xCovLimits()
            self.list['RQCutoffs'] = rqcutoffs = []
            self.printLog('#MINRQ','Min. RQ=%.3f' % self.getNum('MinRQ'))
            if self.getNum('MaxRQ') < self.getNum('MinRQ'): self.setNum({'MaxRQ':self.getNum('MinRQ')})
            self.printLog('#MAXRQ','Max. RQ=%.3f' % self.getNum('MaxRQ'))
            rqstep = max(0.001,self.getNum('RQStep'))
            rq = self.getNum('MinRQ')
            while rq <= self.getNum('MaxRQ'):
                rqcut = rje.dp(rq,3)    # Reduce to 3dp
                if rqcut not in zdb.indexKeys('RQ'):
                    self.warnLog('Failed to find RQ=%s in ZMW file.' % rqcut)
                elif rqcut not in rqcutoffs: rqcutoffs.append(rqcut)
                rq += rqstep

            #?# Make this more like parameters:
            # Swap udb /zdb preference
            # Only retain SeedX values for SeedMinX = TargetXCov + XMargin
            # Only retain AnchorX values for AnchorMinX = TargetX + XMargin
            # Fill in SeedX and AnchorX from lengths and zentries
            ### ~ [1] Calculate unique X and coverage for given seed, anchor + RQ combination ~~~~~~~~~~~~~~~~~~~~~~~ ###
            pdb = db.addEmptyTable('cutoffs',['RQ','SeedMinX','SeedLen','SeedX','AnchorMinX','AnchorLen','AnchorX'],['RQ','SeedMinX','AnchorMinX'])
            # RQ = Read Quality Cutoff
            # SeedLen = Seed length Cutoff
            # SeedX = Total Seed read XCoverage (all reads)
            # SeedMinX = Min seed XCoverage for TargetCov % of genome (unique reads)
            # AnchorX = Total Anchor read XCoverage (all reads)
            # AnchorMinX = Min anchor XCoverage for TargetCov % of genome (unique reads)
            # Needs xcovlimits, rqcutoffs, zdb and udb from above
            zsorted = zdb.sortedEntries('Len',reverse=True)    # List of all entries, sorted by length
            usorted = udb.sortedEntries('Len',reverse=True)    # List of all entries, sorted by length
            for rq in rqcutoffs:
                # RQ and TargetError XCoverage
                if rq > 1: raise ValueError('RQ = %s' % rq)
                x = 1
                while ((1-rq) ** x) > self.getNum('TargetErr'): x += 1
                anchx = x
                targetx = min(x,self.getInt('MinAnchorX'))
                rqentry = False     # Whether anything has been output at this RQ
                #?# Should this be saved as a table somewhere?
                self.printLog('#XERR','Need %dX per base for %s per base error rate @RQ=%s' % (x,rje.expectString(self.getNum('TargetErr')),rq))
                #self.printLog('#XERR','Need %dX per base @RQ=%s for %s error.' % (targetx,rq,self.getNum('TargetErr')))
                if targetx > len(xcovlimits): continue
                # First, establish seed read lengths for each SeedX value
                seedx = [0]     # List of seed read lengths for different SeedX values
                rqzsum = 0      # Sum of all reads
                rqusum = 0      # Sum of all unique reads
                zentries = zsorted[0:]
                uentries = usorted[0:]
                urange = [self.getInt('TargetXCov'),self.getInt('TargetXCov')+self.getInt('XMargin')]
                while uentries: # First establish unique seed read coverage
                    self.progLog('\r#SEED','Establishing seed/anchor lengths/coverage %dX@RQ=%s...' % (len(seedx),rq))
                    while uentries and xcovlimits[len(seedx)] > rqusum: # Need more sequence
                        uentry = uentries.pop(0)
                        if uentry['RQ'] < rq: continue
                        rqusum += (uentry['Len'] * efficiency)
                        while zentries and zentries[0]['Len'] >= uentry['Len']:
                            zentry = zentries.pop(0)
                            if zentry['RQ'] < rq: continue
                            rqzsum += zentry['Len']
                    # Check for end where there is no more sequence available
                    if xcovlimits[len(seedx)] > rqusum and not uentries: break
                    # Update seedlen and calculate seedx
                    seedx.append(uentry['Len'])
                    if len(seedx) - 1 < urange[0]: continue #?# Should we include these?
                    if len(seedx) - 1 > urange[1]: break    # Enough output. Move along.
                    pcore = {'RQ':rq,'SeedLen':uentry['Len'],'SeedX':rqzsum/self.getNum('GenomeSize'),'SeedMinX':len(seedx)-1}
                    # Process anchor lengths
                    anchorx = [0]
                    aqzsum = 0; aqusum = 0
                    azentries = zentries[0:]     # Work through remaining subreads
                    auentries = uentries[0:]     # Work through remaining subreads
                    arange = [targetx,anchx+self.getInt('XMargin')]
                    while auentries:
                        while auentries and xcovlimits[len(anchorx)] > aqusum:
                            auentry = auentries.pop(0)
                            if auentry['RQ'] < rq: continue
                            aqusum += auentry['Len']
                            while azentries and azentries[0]['Len'] >= auentry['Len']:
                                azentry = azentries.pop(0)
                                if azentry['RQ'] < rq: continue
                                aqzsum += azentry['Len']
                        # Check for end
                        if not auentries and xcovlimits[len(anchorx)] > aqusum: break
                        anchorx.append(auentry['Len'])
                        #?# if len(anchorx) <= self.getInt('MinAnchorX'): continue
                        if len(anchorx) - 1 < arange[0]: continue #?# Should we include some of these?
                        if len(anchorx) - 1 > arange[1]: break    # Enough output. Move along.
                        # Generate entry
                        pentry = rje.combineDict({'AnchorLen':azentry['Len'],'AnchorX': aqzsum/self.getNum('GenomeSize'),'AnchorMinX':len(anchorx)-1},pcore)
                        rqentry = pdb.addEntry(pentry)
                if rqentry: continue   # Move along to next RQ

                # Now need to calculate the best we can do given the lack of data
                if rqusum < self.getNum('XStepLen'): continue   # So little data, not even worth it!
                best_seedx = self.seedXForBestCoverage(rqusum/self.getNum('GenomeSize'))
                best_seedcov = best_seedx * self.getNum('GenomeSize')
                rqzsum = 0      # Sum of all reads
                rqusum = 0      # Sum of all unique reads
                zentries = zsorted[0:]
                uentries = usorted[0:]
                seedx = [0]
                self.progLog('\r#SEED','Establishing seed/anchor lengths/coverage %.2fX@RQ=%s...' % (best_seedx,rq))
                while uentries and best_seedcov > rqusum: # Need more sequence
                    uentry = uentries.pop(0)
                    if uentry['RQ'] < rq: continue
                    rqusum += (uentry['Len'] * efficiency)
                    while zentries and zentries[0]['Len'] >= uentry['Len']:
                        zentry = zentries.pop(0)
                        if zentry['RQ'] < rq: continue
                        rqzsum += zentry['Len']
                    if xcovlimits[len(seedx)] <= rqusum: seedx.append(uentry['Len'])
                # Update seedlen and calculate seedx
                pcore = {'RQ':rq,'SeedLen':uentry['Len'],'SeedX':rqzsum/self.getNum('GenomeSize'),'SeedMinX':len(seedx)-1}
                if not pcore['SeedX']: continue     # Just not enough data!
                if not pcore['SeedMinX']:   # Calculate percentage coverage instead
                    pcore['SeedMinX'] = rje.logPoisson(1, rqusum/self.getNum('GenomeSize'),exact=False,callobj=self)
                # Process anchor lengths
                anchorx = [0]
                aqzsum = 0; aqusum = 0
                azentries = zentries[0:]     # Work through remaining subreads
                auentries = uentries[0:]     # Work through remaining subreads
                arange = [targetx,anchx+self.getInt('XMargin')]
                while auentries:
                    while auentries and xcovlimits[len(anchorx)] > aqusum:
                        auentry = auentries.pop(0)
                        if auentry['RQ'] < rq: continue
                        aqusum += auentry['Len']
                        while azentries and azentries[0]['Len'] >= auentry['Len']:
                            azentry = azentries.pop(0)
                            if azentry['RQ'] < rq: continue
                            aqzsum += azentry['Len']
                    # Check for end if lower X coverage have been output
                    if not auentries and xcovlimits[len(anchorx)] > aqusum and rqentry: break
                    if xcovlimits[len(anchorx)] <= aqusum: anchorx.append(auentry['Len'])
                    #?# if len(anchorx) <= self.getInt('MinAnchorX'): continue
                    if len(anchorx) < arange[0] and (rqentry or auentries): continue
                    if len(anchorx) > arange[1]: break    # Enough output. Move along.
                    # Generate entry
                    pentry = rje.combineDict({'AnchorLen':azentry['Len'],'AnchorX': aqzsum/self.getNum('GenomeSize'),'AnchorMinX':len(anchorx)-1},pcore)
                    if pentry['AnchorMinX'] < self.getInt('MinAnchorX'):   # Calculate percentage coverage instead
                        pentry['AnchorMinX'] = rje.logPoisson(self.getInt('MinAnchorX'), aqusum/self.getNum('GenomeSize'),exact=False,callobj=self)
                    rqentry = pdb.addEntry(pentry)

            self.printLog('\r#SEED','Established seed lengths and anchor lengths for different RQ.')
            pdb.saveToFile()
        except: self.errorLog('%s.calculate error' % self.prog())
#########################################################################################################################
    ### <5> ### PacBio Parameter Optimisation Methods                                                                   #
#########################################################################################################################
    def parameters(self):   ### Generate predicted optimum assembly settings from summary tables.
        '''Generate predicted optimum assembly settings from summary tables.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Load data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb = self.db(table='cutoffs',add=True,forcecheck=True,mainkeys=['RQ','SeedMinX','AnchorMinX'])
            if not cdb:
                self.calculate()
                cdb = self.db(table='cutoffs',add=True,forcecheck=False,mainkeys=['RQ','SeedMinX','AnchorMinX'])
                if not cdb: raise IOError('Unable to find/create "cutoffs" table!')
            if not cdb.formatted(): cdb.dataFormat({'RQ':'num','SeedMinX':'num','AnchorMinX':'num','SeedLen':'int','AnchorLen':'int'})
            if not cdb.entries(): self.warnLog('No *.cutoffs.tdt data! Too stringent for data quantity?'); return False
            #?# Add MeanRQ implementation. Is it really useful. Would need a meanRQ(rq,minlen,maxlen) function and then
            #?# work through all entries. Replace rqtargetx with entry['RQTargetX'].
            if self.getBool('RQMean'): self.warnLog('Sorry! RQMean=T not yet implemented!')
            ## ~ [0b] Establish targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Use targetxcov=X and xmargin=X to establish the required data combinations
            targetxcov = self.getInt('TargetXCov')      # Target 100% X Coverage for pre-assembly
            xmargin = self.getInt('XMargin')            # "Safety margin" inflation of X coverage
            if (targetxcov+xmargin) not in cdb.index('SeedMinX'):
                #?# Can we just pick the best one? Will it always be the lowest RQ? (All data) #?#
                self.warnLog('No *.cutoffs.tdt data at SeedMinX=%d! Too stringent for data quantity?' % (targetxcov+xmargin))
                self.printLog('#~~#','# ~~~~~~~~~~ Maximise genome coverage from insufficient data ~~~~~~~~~~~~~~ ##')
                best = (0.0,None)   # (Coverage, cdb entry)
                for centry in cdb.entries():
                    cov = centry['SeedMinX'] * centry['AnchorMinX']
                    #self.bugPrint('%s >> %.4f' % (centry,cov))
                    if cov > best[0] or not best[1]: best = (cov,centry)
                    elif cov < best[0]: continue
                    elif centry['RQ'] < best[1]['RQ']: best = (cov,centry)   # Go for more data!
                entry = best[1]
                rq = entry['RQ']
                self.printLog('#PARAM','Recommended RQ cutoff = %s' % rq)
                self.printLog('#PARAM','Recommended min. subread length = %d' % (entry['AnchorLen']))
                self.printLog('#PARAM','Recommended seed length = %d' % (entry['SeedLen']))
                return True
            cdb.dropEntriesDirect('SeedMinX',[targetxcov+xmargin],inverse=True)
            rqtargetx = {}
            for rq in cdb.indexKeys('RQ'):
                x = 1
                while ((1-rq) ** x) > self.getNum('TargetErr'): x += 1
                self.printLog('#XERR','Need %dX per base for %s per base error rate @RQ=%s' % (x,rje.expectString(self.getNum('TargetErr')),rq))
                targetx = x + xmargin - 1   # Note that the seed read counts as X=1.
                rqtargetx[rq] = x
                for entry in cdb.indexEntries('RQ',rq):
                    if entry['AnchorMinX'] != targetx: cdb.dropEntry(entry)
            ### ~ [1] Output selected parameter settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] MaxLen parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~ Maximise read lengths ~~~~~~~~~~~~~~ ##')
            rq = min(cdb.index('RQ',force=True))
            if len(cdb.indexEntries('RQ',rq)) > 1: self.warnLog('Something strange is afoot!')
            entry = cdb.indexEntries('RQ',rq)[0]
            self.printLog('#PARAM','Recommended RQ cutoff = %s' % rq)
            self.printLog('#PARAM','Recommended min. subread length (%dX + %dX - 1X) = %d' % (rqtargetx[rq],xmargin,entry['AnchorLen']))
            self.printLog('#PARAM','Recommended seed length (%dX + %dX) = %d' % (targetxcov,xmargin,entry['SeedLen']))
            ## ~ [1b] MaxRQ parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','# ~~~~~~~~~~ Maximise read quality ~~~~~~~~~~~~~~ ##')
            rq = max(cdb.index('RQ',force=True))
            if len(cdb.indexEntries('RQ',rq)) > 1: self.warnLog('Something strange is afoot!')
            entry = cdb.indexEntries('RQ',rq)[0]
            self.printLog('#PARAM','Recommended RQ cutoff = %s' % rq)
            self.printLog('#PARAM','Recommended min. subread length (%dX + %dX - 1X) = %d' % (rqtargetx[rq],xmargin,entry['AnchorLen']))
            self.printLog('#PARAM','Recommended seed length (%dX + %dX) = %d' % (targetxcov,xmargin,entry['SeedLen']))
            return True
        except: self.errorLog('%s.parameters error' % self.prog()); return False
#########################################################################################################################
    def stripWierd(self,text):
        maxord = max(ord(char) for char in text)
        if maxord >= 128:
            ascii = ''
            for char in text:
                if ord(char) < 128: ascii += char
            return ascii
        return text
#########################################################################################################################
    ### <6> ### PacBio Genome Coverage Methods                                                                          #
#########################################################################################################################
    def coverageSetup(self):    ### Sets up specific attributes for coverage function.
        '''Sets up specific attributes for coverage function.'''
        try:### ~ [0] ~ Setup Coverage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.baseFile(return_none=None): self.baseFile('pacbio')
            ## ~ [0a] Generate data from summarise *.unique.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            udb = self.db('unique',add=True,mainkeys=['SMRT','ZMW'])
            if udb:
                self.printLog('#SMRT','Calculating SMRT Output from unique data.')
                udb.dropFields(['Pos'])
                if not udb.formatted(): udb.dataFormat({'Len':'int','RQ':'num','Seq':'int'})
                for entry in udb.entries(): entry['Seq'] = 1
                udb.compress(['SMRT','RQ'],default='sum')
                udb.dropFields(['ZMW','RN'])
                udb.makeField('Len*RQ','RQSum')         # Calculate sum of quality bases for mean RQ calculation
                udb.compress(['SMRT'],default='sum')    # Compress to data per SMRT cell
                for entry in udb.entries(): entry['RQ'] = float(entry['RQSum']) / entry['Len']
                for entry in udb.entries(): entry['RQSum'] = 'All'
                udb.compress(['RQSum'],default='mean',rules={'Len':'mean','RQ':'mean','Seq':'mean'})     # Compress to data per SMRT cell
                udb.dropFields(['SMRT'])
                udata = udb.entries()[0]
                self.printLog('#SMRT','Mean Total Length = %.2fMb' % (udata['Len']/1e6))
                self.printLog('#SMRT','Average Read Length = %.2fkb' % (udata['Len']/udata['Seq']/1e3))
                self.printLog('#SMRT','Mean RQ = %.3f' % udata['RQ'])
                self.printLog('#SMRT','Mean SMRTReads = %s' % rje.iStr(udata['Seq']))
                self.setStr({'SMRTUnits':'reads'})
                self.setNum({'SMRTReads':udata['Seq'],'AvRead':udata['Len']/udata['Seq'],'ErrPerBase':udata['RQ']})
                self.list['Accuracy'] = [0,1.0 - self.getNum('ErrPerBase')]     # Update
            ## ~ [0b] SMRTReads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while self.getStrLC('SMRTUnits') not in ['reads','gb','mb']:
                txt = 'SMRTUnits "%s" not recognised'
                if self.getNum('SMRTReads') < 10: smrtunits = 'Gb'
                elif self.getNum('SMRTReads') > 10000: smrtunits = 'reads'
                else: smrtunits = 'Mb'
                if self.i() < 0 or rje.yesNo('%s: switch to (%s) %s?' % (txt,self.getNum('SMRTReads'),smrtunits)):
                    self.setStr({'SMRTUnits':smrtunits})
                elif self.i() >0: self.setStr({'SMRTUnits':rje.choice('SMRTUnits (reads/Gb/Mb)?')})
                self.printLog('#UNITS','%s => %s' % (txt,self.getStr('SMRTUnits')))
            if self.getStrLC('SMRTUnits') in ['gb','mb']:
                smrttotal = self.getNum('SMRTReads') * {'gb':1e9,'mb':1e6}[self.getStrLC('SMRTUnits')]
                txt =  '%s %s @ %.3f kb/read' % (self.getNum('SMRTReads'),self.getStr('SMRTUnits'),self.getNum('AvRead')/1000.0)
                self.setNum({'SMRTReads':smrttotal/self.getNum('AvRead')})
                txt += ' => %s reads' % rje.iStr(int(self.getNum('SMRTReads')))
                self.printLog('#READS',txt)
            ## ~ [0c] XnList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xnlist = range(self.getInt('MinAnchorX')+1,self.getInt('TargetXCov')+self.getInt('MinAnchorX')+1)
            for xn in self.list['XnList']:
                if xn == '': continue
                try:
                    ixn = int(xn)
                    if xn not in [ixn,'%d' % ixn]: self.printLog('#XN','"%s" -> %dX' % (xn,ixn))
                    if ixn == 0: self.printLog('#XN','No point in 0X output: use 1-%Coverage.')
                    elif ixn == 1: self.printLog('#XN','No point in 1X output: use %Coverage.')
                    else: xnlist.append(ixn)
                except: self.errorLog('Could not process %s as part of XnList. (Integers only.)' % xn)
            xnlist.sort()
            if xnlist: self.printLog('#XN','XnList: %sX.' % string.join(string.split('%s' % xnlist,','),'X, ')[1:-1])
            self.list['XnList'] = xnlist
            return True     # Setup successful
        except: self.errorLog('Problem during %s coverageSetup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def coverage(self): ### Calculates estimated % coverage and accuracy of genome sequencing.
        '''Calculates estimated % coverage and accuracy of genome sequencing.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] General Settings setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.coverageSetup()    # Sets up smrtread stats and xnlist
            self.targetXDepth()     # Sets up target X depth lists
            self.accuracy(self.getInt('MaxCov'))
            ## ~ [1b] Setup table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('BySMRT'): ckey = 'SMRT'
            else: ckey = 'XCoverage'
            cfields = ['XCoverage','SMRT','%Coverage','%Accuracy','%Complete','SeedX']
            # XCoverage = Total mean depth of coverage
            # SMRT = Total number of SMRT cells
            # %Coverage = Estimated %coverage of assembly
            # %Accuracy = Estimated %accuracy of assembled genome
            # %Complete = Product of %coverage and %accuracy
            # SeedX = Estimated optimal depth of coverage for seed reads
            for xn in self.list['XnList']: cfields.append('%%X%d' % xn)
            cdb = self.db().addEmptyTable('coverage',cfields,[ckey])
            ## ~ [1c] Setup stats for coverage calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cov_per_base_per_read = self.getNum('AvRead') / self.getNum('GenomeSize')
            # Calculate reads = number of reads per cycle of calculations
            if self.getBool('BySMRT'): reads = self.getInt('SMRTReads')                 # If going per SMRT cell
            else: reads = int(0.5 + self.getNum('GenomeSize') / self.getNum('AvRead'))  # if going per X coverage

            ### ~ [2] Calculate stats for each round of sequencing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xcoverage = 0      # Total X coverage of sequencing
            while xcoverage < self.getNum('MaxCov'):    # Cycle until limits reached
                self.progLog('\r#XCOV','Calculating coverage stats: %.1f%%' % (100.0 * xcoverage / self.getNum('MaxCov')))
                # Add sequence
                xcoverage += self.getNum('AvRead') * reads / self.getNum('GenomeSize')
                # Number of SMRT cells
                smrt = (self.getNum('GenomeSize') * xcoverage) / (self.getNum('AvRead') * self.getNum('SMRTReads'))
                # SeedX for determining coverage
                seedx = self.seedXForBestCoverage(xcoverage)
                anchx = xcoverage - seedx
                pccov = 100.0 * rje.logPoisson(1,seedx,exact=False,callobj=self) * rje.logPoisson(self.getInt('MinAnchorX'),xcoverage-seedx,exact=False,callobj=self)
                pcacc = 0.0
                # Calculate X coverage counts using binomial
                #self.dict['XCovPerBase'] = {}   # Dictionary of {TotalX:[List where index is X coverage and number is proportion of reads]}
                bases = int(self.getNum('GenomeSize'))
                xcov = []   # List where index is X coverage and number is proportion of reads
                while bases > 1:
                    xcov.append(rje.logPoisson(len(xcov),anchx,exact=True,callobj=self))
                    bases -= self.getNum('GenomeSize') * xcov[-1]
                    if len(xcov) > reads: raise ValueError('XCoverage cannot exceed read count!')
                covsum = sum(xcov[self.getInt('MinAnchorX'):])
                #self.debug(covsum)
                if covsum:
                    for x in range(self.getInt('MinAnchorX'),len(xcov)):
                        #self.debug(self.accuracy(x+1))   # Seed too, so x+1!
                        pcacc += xcov[x] * self.accuracy(x+1)   # Seed too, so x+1!
                    pcacc = 100.0 * pcacc / covsum
                elif pccov: raise ValueError('Somehow have anchor coverage without 0% genome coverage.')
                # XnList
                centry = {'XCoverage':rje.sf(xcoverage,3),'SMRT':rje.sf(smrt,3),'%Coverage':rje.sf(pccov,5),
                          '%Accuracy':rje.sf(pcacc,5),'%Complete':rje.sf(pccov*pcacc/100.0,5),'SeedX':rje.sf(seedx,3)}
                bases = int(self.getNum('GenomeSize'))
                xcov = []   # List where index is X coverage and number is proportion of reads
                while bases > 1:
                    xcov.append(rje.logPoisson(len(xcov),xcoverage,exact=True,callobj=self))
                    bases -= self.getNum('GenomeSize') * xcov[-1]
                    if len(xcov) > reads: raise ValueError('XCoverage cannot exceed read count!')
                #self.bugPrint(xcov)
                for xn in self.list['XnList']:
                    #self.bugPrint('%d >> %s = %s' % (xn,xcov[xn:],sum(xcov[xn:])))
                    if xn <= len(xcov): centry['%%X%d' % xn] = rje.sf(100.0*sum(xcov[xn:]),5)
                    else: centry['%%X%d' % xn] = 0.000
                cdb.addEntry(centry)
                #self.debug(centry)
            self.printLog('\r#XCOV','Calculated coverage stats upto %dX coverage.' % self.getInt('MaxCov'))

            ### ~ [4] Save results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for xkey in cdb.dataKeys():
                cdb.dict['Data'][float(xkey)] = cdb.dict['Data'].pop(xkey)
            cdb.saveToFile()

            return
        except: self.errorLog('%s.coverage error' % self.prog())
#########################################################################################################################
    def accuracy(self,xdepth):   ### Calculate accuracy (if required) at xdepth and returns
        '''Calculate accuracy (if required) at xdepth and returns.'''
        try:### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while len(self.list['Accuracy']) <= xdepth:
                if not int((1.0 - self.list['Accuracy'][-1]) * self.getNum('GenomeSize')): # Too few errors to worry!
                    self.list['Accuracy'].append(1.0)
                    continue
                xcov = len(self.list['Accuracy'])
                majority = int(0.33*xcov/1.33) + 1        # Number of correct reads needed for majority
                #self.debug(majority)
                try: self.list['Accuracy'].append(rje.logBinomial(majority,xcov,1.0 - self.getNum('ErrPerBase'),exact=False,callobj=self))
                except: self.list['Accuracy'].append(rje.logPoisson(majority,xcov*(1.0 - self.getNum('ErrPerBase')),exact=False,callobj=self))
                #self.debug(self.list['Accuracy'])
            return self.list['Accuracy'][xdepth]
        except: self.errorLog('%s.accuracy error' % self.prog())
#########################################################################################################################
    def seedXForBestCoverage(self,xcoverage):   ### Calculate optimal seedX for best %coverage given Xcoverage
        '''Calculate optimal seedX for best %coverage given Xcoverage.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seeddepth = self.getInt('TargetXCov')   # This is deemed required for optimal assembly
            anchdepth = self.getInt('MinAnchorX')   # This is absolutely required for seed read correction/retention
            efficiency = self.getNum('MapEfficiency')   # This is the proportion of anchor reads that map OK
            self.printLog('#BEST','Calculating best SeedX for min. %dX seed and %dX anchor coverage at %sX depth (Map Efficiency=%.3f)' % (seeddepth,anchdepth,rje.sf(xcoverage,3),efficiency))
            # Use: rje.poisson(xdepth,xcoverage,exact=False,callobj=self) >= targetcov
            ### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # First, calculate the total depth of coverage needed for target coverage upto anchdepth + seeddepth
            maxlen = max((anchdepth + seeddepth + 1),xcoverage) * self.getNum('GenomeSize')
            targetxdepth = self.targetXDepth(maxlen,anchdepth + seeddepth + 1)  # List of required depth for each Target X index
            # See if maximum requirements are met = easy!
            # Efficiency is purely a seed -> data-corrected correction.
            seedcov = targetxdepth[seeddepth] / efficiency              # XDepth for seed allowing for efficiency
            anchcov = (xcoverage - seedcov)                             # XDepth for anchors
            if targetxdepth[anchdepth] <= anchcov:
                self.printLog('#BEST','SeedX=%s; AnchorX=%s => min. %dX seed and %dX anchor coverage (Map Efficiency=%.3f)' % (rje.sf(seedcov,3),rje.sf(anchcov,3),seeddepth,anchdepth,efficiency))
                return seedcov
            # Otherwise, see if anchdepth can be met
            seedx = 0
            seedcov = targetxdepth[seedx] / efficiency              # XDepth for seed allowing for efficiency
            anchcov = (xcoverage - seedcov)                         # XDepth for anchors
            while targetxdepth[anchdepth] <= anchcov and seedx <= seeddepth: # Can meet anchor, so try increasing seed
                nextseed = targetxdepth[seedx+1] / efficiency       # XDepth for seed allowing for efficiency
                nextanch = (xcoverage - nextseed)                   # XDepth for anchors
                if targetxdepth[anchdepth] > nextanch: break        # Cannot meet anchor coverage at increased seedX
                # Icrement seedx and loop
                seedx += 1
                seedcov = nextseed
                anchcov = nextanch
            if seedx:
                self.printLog('#BEST','SeedX=%s; AnchorX=%s => min. %dX seed and %dX anchor coverage (Map Efficiency=%.3f)' % (rje.sf(seedcov,3),rje.sf(anchcov,3),seedx,anchdepth,efficiency))
                return seedcov
            # Finally, try to optimise poor xcoverage where:
            # rje.poisson(1,seedx,exact=False,callobj=self) = rje.poisson(anchdepth,xcoverage-seedx,exact=False,callobj=self)
            xstep = self.getNum('XStepLen') / self.getNum('GenomeSize')
            seedx = 0.0
            bestx = (0.0,1.0)   # (seedx, probability difference)
            while bestx[0] < xcoverage:     # Work through until optimal balance (or max xcoverage) reached
                seedx += xstep
                pdif = rje.modulus(rje.poisson(1,seedx*efficiency,exact=False,callobj=self) - rje.poisson(anchdepth,(xcoverage-seedx),exact=False,callobj=self))
                if pdif > bestx[1]: break
                bestx = (seedx,pdif)
            anchcov = (xcoverage - bestx[0])
            self.printLog('#BEST','SeedX=%s (%s 1X coverage); AnchorX=%s (%s %dX coverage). (Map Efficiency=%.3f)' % (rje.sf(bestx[0],3),rje.sf(rje.poisson(1,bestx[0]*efficiency,exact=False,callobj=self),3),rje.sf(anchcov,3),rje.sf(rje.poisson(anchdepth,(xcoverage-bestx[0]),exact=False,callobj=self),3),anchdepth,efficiency))
            return bestx[0]
        except: self.errorLog('%s.seedXForBestCoverage error' % self.prog())
#########################################################################################################################
    def targetXDepth(self,maxlen=0,maxdepth=0): ### Calculate list of target x (index) and required depth (value).
        '''Calculate list of target x (index) and required depth (value).'''
        try:### ~ [1] Calculate list of target x (index) and required depth (value) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['TargetXDepth']:
                if self.list['TargetXDepth'][-1] * self.getNum('GenomeSize') > maxlen: return self.list['TargetXDepth']
            else: self.list['TargetXDepth'] = [0]
            xcovlimits = self.list['TargetXDepth']    # List of target x (index) and required depth (value)
            xlen = steplen = self.getNum('XStepLen')
            if not maxlen: maxlen = max(self.getInt('TargetXCov')+self.getInt('MinAnchorX'),self.getInt('MaxCov')) * self.getNum('GenomeSize')
            targetcov = math.sqrt(self.getPerc('TargetCov'))  # Needs square root for anchor x seed
            # Setup TargetXDepth for each XCoverage
            while xlen <= maxlen or len(xcovlimits) <= maxdepth:
                self.progLog('\r#XDEPTH','Setting up TargetXDepth: %dX (%s bp)' % (len(xcovlimits),rje.iStr(xlen)))
                if xlen/self.getNum('GenomeSize') < len(xcovlimits): xlen += steplen; continue
                while rje.poisson(len(xcovlimits),xlen/self.getNum('GenomeSize'),exact=False,callobj=self) >= targetcov:
                    xcovlimits.append(xlen/self.getNum('GenomeSize'))
                xlen += steplen
            self.printLog('\r#XDEPTH','Set up TargetXDepth for up to %dX (%s bp).' % (len(xcovlimits)-1,rje.iStr(xlen)))
            return self.list['TargetXDepth']
        except: self.errorLog('%s.targetXDepth error' % self.prog())
#########################################################################################################################
    ### <7> ### PacBio Assembly Parameter Parsing Methods                                                               #
#########################################################################################################################
    def parseParam(self):   ### Parses assembly parameters, reorganises and saves
        '''Parses assembly parameters, reorganises and saves.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db(); pdb = None
            if not self.force(): pdb = self.db('settings',add=True,mainkeys=['Setting'])
            if not pdb: pdb = db.addEmptyTable('settings',['Setting','Prefix','Suffix','Variable'],['Setting'])
            if self.list['ParamList']: self.printLog('#PARAM','Restricting ParseParam to %d parameters.' % len(self.list['ParamList']))
            ### ~ [1] Parse parameter files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pfile in self.list['ParseParam']:
                pname = rje.baseFile(pfile,strip_path=True)
                self.printLog('#PARSE','Parsing %s parameters.' % pname)
                if pname in pdb.fields():
                    if not self.force():
                        self.printLog('#SKIP','Skipping previously parsed %s parameter settings.' % pname); continue
                    self.warnLog('Replacing contents of %s parameter settings.' % pname)
                    pdb.dropField(pname)
                pdb.addField(pname)
                for pline in open(pfile,'r').readlines():
                    pline = self.stripWierd(pline)
                    pdata = string.split(rje.chomp(pline),' = ',1)
                    if len(pdata) != 2: continue
                    if self.list['ParamList'] and pdata[0] not in self.list['ParamList']: continue
                    entry = pdb.data(pdata[0])
                    if entry:
                        entry[pname] = pdata[1]
                    else:
                        entry = {'Setting':pdata[0],pname:pdata[1],'Variable':'FALSE'}
                        try: [entry['Prefix'],entry['Suffix']] = string.split(pdata[0],'.')  # Will raise error if wrong
                        except:
                            if pdata[0] != pdata[0].upper(): self.warnLog('Cannot parse %s' % pdata[0])
                            continue
                        pdb.addEntry(entry)
            pfields = pdb.fields()[4:]
            for entry in pdb.entries():
                for p1 in range(len(pfields)-1):
                    for p2 in range(1,len(pfields)):
                        if entry[pfields[p1]] != entry[pfields[p2]]: entry['Variable'] = 'TRUE'; break
            pdb.saveToFile()
            return pdb
        except: self.errorLog('%s.parseParam error' % self.prog()); return False
#########################################################################################################################
    ### <8> ### PacBio Parameter Prediction Methods                                                                     #
#########################################################################################################################
    def predict(self):   ### Predicts coverage from parsed assembly parameters and compares to pre-assembly if possible.
        '''Predicts coverage from parsed assembly parameters and compares to pre-assembly if possible.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Load/create data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = self.parseParam()     # ['Setting','Prefix','Suffix','Variable']+[Assemblies],['Setting']
            zdb = self.zdb()            # ['SMRT','ZMW','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW','RN']
            self.printLog('#~~#','# ~~~~~~~~~~~ Predict PreAssembly Data (%d runs) ~~~~~~~~~~~ #' % len(self.list['ParseParam']))
            db = self.db(); xdb = None
            if not self.force(): xdb = self.db('predict',add=True,mainkeys=['Assembly'])
            xfields = ['Assembly','minSubReadLength','readScore','minLongReadLength','minCorCov','ovlErrorRate']
            xfields += ['SeedN','SeedX','AnchorX','SeedMinX','AnchorMinX','PreCov','CorPreCov']
            # SeedX','AnchorX','SeedMinX','AnchorMinX'
            # PreCov = predicted coverage of pre-assembly
            # CorPreCov = corrected predicted coverage of pre-assembly given mapefficiency
            xfields += ['PreN','PreX','PreMinX','PreMapEfficiency','AncMapEfficiency']
            if not xdb: xdb = db.addEmptyTable('predict',xfields,['Assembly'])
            else: xdb.list['Fields'] = xfields[0:]
            ## ~ [0b] Generate XCov limits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            anchdata = pdb.data('p_filter.minSubReadLength')
            rqdata = pdb.data('p_filter.readScore')
            seeddata = pdb.data('p_preassemblerdagcon.minLongReadLength')
            depthdata = pdb.data('p_preassemblerdagcon.minCorCov')
            eratedata = pdb.data('p_assembleunitig.ovlErrorRate')
            gsizedata = pdb.data('p_assembleunitig.genomeSize')
            compdata = pdb.data('p_preassemblerdagcon.computeLengthCutoff')
            xcovdata = pdb.data('p_assembleunitig.xCoverage')
            ## ~ [0c] Generate XCov limits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xlimits = self.xCovLimits()
            ### ~ [1] Calculate X stats from parameter settings and genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pfile in self.list['ParseParam']:
                pbase = rje.baseFile(pfile,strip_path=False)
                pname = rje.baseFile(pfile,strip_path=True)
                ## ~ [1a] Parameter data predictions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rq = float(rqdata[pname])
                udb = self.udb(rq,force=True)            # ['SMRT','RN','Len','Pos','RQ','Seq'],['SMRT','ZMW']
                anchlen = int(anchdata[pname])
                seedlen = int(seeddata[pname])
                if compdata[pname].upper() == 'TRUE':
                    seedtot = float(gsizedata[pname]) * float(xcovdata[pname])
                    seedsum = 0
                    zsorted = zdb.sortedEntries('Len',reverse=True)
                    while len(zsorted) > 1 and seedsum < seedtot: seedsum += zsorted.pop(0)['Len']
                    seedlen = zsorted.pop(0)['Len']
                    self.printLog('#SEED','Calculated ZMW seed length based on %sX of %s bp genome = %s.' % (xcovdata[pname],rje.sf(float(gsizedata[pname]),3),seedlen))
                    seedsum = 0
                    usorted = udb.sortedEntries('Len',reverse=True)
                    while len(usorted) > 1 and seedsum < seedtot: seedsum += usorted.pop(0)['Len']
                    seedlen = usorted.pop(0)['Len']
                    self.printLog('#SEED','Calculated unique seed length based on %sX of %s bp genome = %s.' % (xcovdata[pname],rje.sf(float(gsizedata[pname]),3),seedlen))
                xdepth = int(depthdata[pname])
                # Depth of coverage
                ax = sx = 0.0
                sn = 0
                for entry in zdb.entries():
                    if entry['RQ'] < rq: continue
                    if entry['Len'] >= seedlen: sx += entry['Len']; sn += 1
                    elif entry['Len'] >= anchlen: ax += entry['Len']
                # MinX Coverage
                aux = sux = 0.0
                for entry in udb.entries():
                    if entry['RQ'] < rq: continue
                    if entry['Len'] >= seedlen: sux += entry['Len']
                    elif entry['Len'] >= anchlen: aux += entry['Len']
                axcov = sxcov = 0
                while aux > xlimits[axcov]: axcov += 1
                while sux > xlimits[sxcov]: sxcov += 1
                prexdepth = self.pcXDepth(sx/self.getNum('GenomeSize'))
                precov = sum(prexdepth[xdepth:])
                xentry = {'Assembly':pname,'minSubReadLength':anchlen,'readScore':rq,'minLongReadLength':seedlen,'minCorCov':xdepth,
                          'SeedX':sx/self.getNum('GenomeSize'),'AnchorX':ax/self.getNum('GenomeSize'),
                          'ovlErrorRate':eratedata[pname],'SeedN':sn,
                          'SeedMinX':sxcov,'AnchorMinX':axcov,'PreCov':precov,'CorPreCov':precov*self.getNum('MapEfficiency')}
                ## ~ [1b] Load in pre-assembly data if possible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                prefile = '%s.preassembly.fasta' % pbase
                self.printLog('#PREX','Pre-assembly sequence file %s: %s' % (prefile,rje.exists(prefile)))
                #self.debug(rje.exists(prefile))
                if rje.exists(prefile):
                    seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','summarise=F','seqin=%s' % prefile])
                    seqdata = seqlist.summarise()   # This will summarise read lengths etc.
                    xentry['PreN'] = seqdata['SeqNum']
                    xentry['PreX'] = float(seqdata['TotLength']) / self.getNum('GenomeSize')
                    pxcov = 0
                    while seqdata['TotLength'] > xlimits[pxcov]: pxcov += 1
                    xentry['PreMinX'] = pxcov
                    #xentry['PreMapEfficiency'] = xentry['PreX'] / xentry['SeedX']
                    xentry['PreMapEfficiency'] = math.pow(xentry['PreX'] / xentry['SeedX'],1.0/xdepth)
                    # Adding AncMapEfficiency = empirical map efficiency based mean anchor depth of coverage
                    xentry['AncMapEfficiency'] = math.pow(xentry['PreX'] / xentry['SeedX'],1.0/ax)
                ## ~ [1c] Add entry to predict table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                xdb.addEntry(xentry)
            ### ~ [2] Calculate X stats from parameter settings and genome size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xdb.saveToFile()
        except: self.errorLog('%s.predict error' % self.prog()); return False
#########################################################################################################################
    ### <9> ### PacBio Parameter Prediction Methods                                                                     #
#########################################################################################################################
    def preassemblyFrag(self):  ### Assesses and tries to correct the preassembly fragmentation issue.
        '''Assesses and tries to correct the preassembly fragmentation issue.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ PACBIO PREASSEMBLY FRAGMENTATION ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #')
            db = self.db()
            if not rje.exists(self.getStr('Preassembly')): raise IOError('Preassembly file "%s" not found!' % self.getStr('Preassembly'))
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','summarise=F','seqin=%s' % self.getStr('Preassembly')])   # This will summarise read lengths etc.
            basename = rje.baseFile(self.getStr('Preassembly'),strip_path=True)
            fdb = self.db('fragment',add=True,mainkeys=['Assembly','Stage'])
            if not fdb: fdb = db.addEmptyTable('fragment',['Assembly','Stage','SeqNum','TotLength','MinLength','MaxLength','MeanLength','MedLength','N50Length','XCoverage','AdjN'],['Assembly','Stage'])
            ### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqdata = seqlist.summarise()
            seqdata['Assembly'] = basename; seqdata['Stage'] = 'preassembly'
            seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
            seqdata['AdjN'] = 0
            ##~ ~ [2] Adjacent fragment analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fixfas = '%s.smrtscape.fas' % basename
            rje.backup(self,fixfas)
            FIXFAS = open(fixfas,'w'); fx = 0
            prevseq = [None,-1,-1,'']
            while seqlist.nextSeq():
                (name,sequence) = seqlist.getSeq()
                (seed,i,j) = rje.matchExp('^(\S+)/(\d+)_(\d+)$',name)
                i = int(i); j = int(j)
                if seed == prevseq[0] and i == prevseq[2] + 1:  # 0bp fragmentation
                    prevseq[2] = j
                    prevseq[3] += sequence
                    seqdata['AdjN'] += 1
                else:
                    if prevseq[0]: FIXFAS.write('>%s/%d_%d\n%s\n' % (prevseq[0],prevseq[1],prevseq[2],prevseq[3])); fx += 1
                    prevseq = [seed,i,j,sequence]
            if prevseq[0]: FIXFAS.write('>%s/%d_%d\n%s\n' % (prevseq[0],prevseq[1],prevseq[2],prevseq[3])); fx += 1
            FIXFAS.close()
            self.printLog('#FAS','%s preassembly reads -> %s output to %s' % (rje.iStr(seqlist.seqNum()),rje.iStr(fx),fixfas))
            fdb.addEntry(seqdata)
            ### ~ [3] Calculate FixFas data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fixseq = rje_seqlist.SeqList(self.log,self.cmd_list+['seqmode=file','autoload=T','summarise=F','seqin=%s' % fixfas])   # This will summarise read lengths etc.
            seqdata = fixseq.summarise()
            seqdata['Assembly'] = basename; seqdata['Stage'] = 'smrtscape'
            seqdata['XCoverage'] = seqdata['TotLength'] / self.getNum('GenomeSize')
            seqdata['AdjN'] = 0
            fdb.addEntry(seqdata)
            ### ~ [4] Finish and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb.saveToFile()
            return fdb
        except: self.errorLog('%s.preassemblyFrag error' % self.prog()); return False
#########################################################################################################################
### End of SECTION II: SMRTSCAPE Class                                                                                  #
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
    try: SMRTSCAPE(mainlog,cmd_list).run()

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
