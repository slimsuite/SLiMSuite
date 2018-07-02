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
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       SLiMFarmer
Description:  SLiMSuite HPC job farming control program
Version:      1.9.0
Last Edit:    30/05/18
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to control and execute parallel processing jobs on an HPC cluster using PBS and QSUB. If
    qsub=T it will generate a job file and use qsub to place that job in the queue using the appropriate parameter
    settings. If `slimsuite=T` and `farm=X` gives a recognised program (below) or `hpcmode` is not `fork` then the qsub
    job will call SLiMFarmer with the same commandline options, plus `qsub=F i=-1 v=-1`. If `seqbyseq=T`, this will be
    run in a special way. (See SeqBySeq mode.) Otherwise `slimsuite=T` indicates that `farm=X` is a SLiMSuite program,
    for which the python call and `pypath` will be added. If this program uses forking then it should parallelise over a
    single multi-processor node. If `farm=X` contains a `/` path separator, this will be added to `pypath`, otherwise it
    will be assumed that `farm` is in `tools/`. If `slimsuite=F` then farm should be a program call to be queued in the
    PBS job file instead.

    Currently recognised SLiMSuite programs for farming: SLiMFinder, QSLiMFinder, SLiMProb, SLiMCore.

    Currently recognised SLiMSuite programs for rsh mode only qsub farming: GOPHER, SLiMSearch, UniFake.

    NOTE: Any commandline options that need bracketing quotes will need to be placed into an ini file. This can either
    be the ini file used by SLiMFarmer, or a `jobini=FILE` that will only be used by the farmed programs. Note that
    commands in `slimfarmer.ini` will not be passed on to other SLiMSuite programs unless `ini=slimfarmer.ini` is given
    as a commandline argument.

    The runid=X setting is important for SLiMSuite job farming as this is what separates different parameter setting
    combinations run on the same data and is also used for identifying which datasets have already been run. Running
    several jobs on the same data using the same SLiMSuite program but with different parameter settings will therefore
    cause problems. If runid is not set, it will default to the job=X setting.

    The hpcmode=X setting determines the method used for farming out jobs across the nodes. hpcmode=rsh uses rsh to spawn
    the additional processes out to other nodes, based on a script written for the IRIDIS HPC by Ivan Wolton.
    hpcmode=fork will restrict analysis to a single node and use Python forking to distribute jobs. This can be used even
    on a single multi-processor machine to fork out SLiMSuite jobs. basefile=X will set the log, RunID, ResFile, ResDir
    and Job: RunID and Job will have path stripped; ResFile will have .csv appended.

    Initially, it will call other programs but, in time, it is envisaged that other programs will make use of SLiMFarmer
    and have parallelisation built-in.

SeqBySeq Mode:
    In SeqBySeq mode, the program assumes that seqin=FILE and basefile=X are given and farm=X states the Python program
    to be run, which should be SLiMSuite program. (The SLiMSuite subdirectory will also need to be given unless
    slimsuite=F, in which case the whole path to the program should be given. pypath=PATH can set an alternative path.)

    Seqin will then be worked through in turn and each sequence farmed out to the farm program. Outputs given by OutList
    are then compiled, as is the Log, into the correct basefile=X given. In the case of *.csv and *.tdt files, the header
    row is copied for the first file and then excluded for all subsequent files. For all other files extensions, the
    whole output is copied.

Commandline:
    ### ~ Basic QSub Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    qsub=T/F        : Whether to execute QSub PDB job creation and queuing [False]
    jobini=FILE     : Ini file to pass to the farmed HPC jobs with SLiMFarmer options. Overrides commandline. [None]
    slimsuite=T/F   : Whether program is an RJE *.py script (adds log processing) [True]
    nodes=X         : Number of nodes to run on [1]
    ppn=X           : Processors per node [16]
    walltime=X      : Walltime for qsub job (hours) [12]
    vmem=X          : Virtual Memory limit for run (GB) [126]
    job=X           : Name of job file (.job added) [slimfarmer]

    ### ~ Advanced QSub Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    hpc=X           : Name of HPC system ['katana']
    pypath=PATH     : Path to python modules [slimsuite home directoy]
    qpath=PATH      : Path to change directory too [current path]
    pause=X         : Wait X seconds before attempting showstart [5]
    email=X         : Email address to email job stats to at end ['']
    mailstart=T/F   : Whether to email user at start of run [False]
    depend=LIST     : List of job ids to wait for before starting job (dependhpc=X added) []
    dependhpc=X     : Name of HPC system for depend ['blue30.iridis.soton.ac.uk']
    report=T/F      : Pull out running job IDs and run showstart [False]
    modules=LIST    : List of modules to add in job file e.g. blast+/2.2.31,clustalw []
    modpurge=T/F    : Whether to purge loaded modules in qsub job file prior to loading [True]
    precall=LIST    : List of additional commands to run between module loading and program call []
    daisychain=X    : Chain together a set of qsub runs of the same call that depend on the previous job.

    ### ~ Main SLiMFarmer Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    farm=X          : Execute a special SLiMFarm analysis on HPC [batch]
                        - batch will farm out a batch list of commands read in from subjobs=LIST
                        - gopher/slimfinder/qslimfinder/slimprob/slimcore/slimsearch/unifake = special SLiMSuite HPC.
                        - if seqbyseq=T, farm=X will specify the program to be run (see docs)
                        - otherwise, farm=X will be executed as a system call in place of SLiMFarmer
    hpcmode=X       : Mode to be used for farming jobs between nodes (rsh/fork) [fork]
    forks=X         : Number of forks to be used when hpcmode=fork and qsub=F. [1]
    jobini=FILE     : Ini file to pass to the farmed SLiMSuite run. (Also used for SLiMFarmer options if qsub=T.) [None]
    jobforks=X      : Number of forks to pass to farmed out run if >0 [0]

    ### ~ Standard HPC Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    subsleep=X      : Sleep time (seconds) between cycles of subbing out jobs to hosts [1]
    subjobs=LIST    : List of subjobs to farm out to HPC cluster []
    iolimit=X       : Limit of number of IOErrors before termination [50]
    memfree=X       : Min. proportion of node memory to be free before spawning job [0.1]
    test=T/F        : Whether to produce extra output in "test" mode [False]
    keepfree=X      : Number of processors to keep free on head node [1]

    ### ~ SeqBySeq Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqbyseq=T/F    : Activate seqbyseq mode - assumes basefile=X option used for output [False]
    seqin=FILE      : Input sequence file to farm out [None]
    basefile=X      : Base for output files - compiled from individual run results [None]
    outlist=LIST    : List of extensions of outputs to add to basefile for output (basefile.*) []
    pickhead=X      : Header to extract from OutList file and used to populate AccNum to skip []

    ### ~ SLiMSuite Farming Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    runid=X         : Text identifier for SLiMSuite job farming [`job`]
    resfile=FILE    : Main output file for SLiMSuite run [`farm`.csv]
    pickup=T/F      : Whether to pickup previous run based on existing results and RunID [True]
    sortrun=T/F     : Whether to sort input files by size and run big -> small to avoid hang at end [True]
    loadbalance=T/F : Whether to split SortRun jobs equally between large & small to avoid memory issues [True]
    basefile=X      : Set the log, RunID, ResFile, ResDir and Job to X [None].

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_hpc, rje_obj, rje_iridis, rje_qsub
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Functional version using rje_qsub and rje_iridis to fork out SLiMSuite runs.
    # 1.1 - Updated to use rje_hpc.JobFarmer and incorporate main SLiMSuite farming within SLiMFarmer class.
    # 1.2 - Implemented the slimsuite=T/F option and got SLiMFarmer qsub to work with GOPHER forking.
    # 1.3 - Modified default vmem request to 127GB from 64GB.
    # 1.4 - Added modules=LIST : List of modules to add in job file [clustalo,mafft]
    # 1.4.1 - Fixed farm=batch mode for qsub=T.
    # 1.4.2 - Fixed log transfer issues due to new #VIO line. Better handling of crashed runs.
    # 1.4.3 - Added recognition of missing slimsuite programs and switching to slimsuite=F.
    # 1.4.4 - Modified default vmem request to 126GB from 127GB.
    # 1.4.5 - Updated BLAST loading default to 2.2.31
    # 1.5.0 - mailstart=T/F : Whether to email user at start of run [False]
    # 1.6.0 - modpurge=T/F : Whether to purge loaded modules in qsub job file prior to loading [True]
    # 1.7.0 - precall=LIST : List of additional commands to run between module loading and program call []
    # 1.8.0 - jobforks=X : Number of forks to pass to farmed out run if >0 [0]
    # 1.9.0 - daisychain=X : Chain together a set of qsub runs of the same call that depend on the previous job.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Incorporate rje_iridis (rename rje_hpc) and rje_qsub.
    # [ ] : Test simple farming out of batch job by using farm=batch subjobs=LIST.
    # [Y] : Replace core rje_iridis with rje_hpc. Inherit.
    # [Y] : Move rje_iridis iRun code and functions to SLiMFarmer.
    # [Y] : Fix RunID issue. (Not working.)
    # [ ] : Test seqbyseq=T mode.
    # [ ] : Add farm=['gopher','slimsearch','unifake']
    # [ ] : Check/sort out job naming oddity.
    # [Y] : Add recognition of missing slimsuite programs and switching to slimsuite=F.
    # [ ] : Add running of several commands in serial.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SLiMFarmer', '1.9.0', 'May 2018', '2014')
    description = 'SLiMSuite HPC job farming control program'
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
### SECTION II: SLiMFarmer Class                                                                                        #
#########################################################################################################################
class SLiMFarmer(rje_hpc.JobFarmer):
    '''
    SLiMFarmer Class. Author: Rich Edwards (2014).

    Str:str
    - Farm = Execute a special SLiMFarm analysis on HPC [batch]
    - HPC = Name of HPC system ['IRIDIS4']
    - HPCMode = Mode to be used for farming jobs between nodes (rsh/fork) [fork]
    - Job = Name of job file (.job added) [slimfarmer]
    - JobINI = INI file to pass to the HPC job with SLiMFarmer options [None]
    - PyPath = Path to python modules [slimsuite home directoy]
    - ResFile = Main output file for SLiMSuite run [`farm`.csv]
    - RunID = Text identifier for SLiMSuite job farming [`job`]

    Bool:boolean
    - LoadBalance = Whether to split SortRun jobs equally between large & small to avoid memory issues [True]
    - Pickup = Whether to pickup previous run based on existing results and RunID [True]
    - QSub = Whether to execute QSub PDB job creation and queuing [False]
    - SeqBySeq = Activate seqbyseq mode - assumes basefile=X option used for output [False]
    - SLiMSuite = Whether program is an RJE *.py script (adds log processing) [True]
    - SortRun = Whether to sort input files by size and run big -> small to avoid hang at end [True]

    Int:integer
    - DaisyChain=X    : Chain together a set of qsub runs of the same call that depend on the previous job.
    - JobForks = Number of forks to pass to farmed out run if >0 [0]
    - KeepFree = Number of processors to keep free on head node [1]
    - Modules = List of modules to add in job file []

    Num:float
    
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
        self.strlist = ['Farm','HPC','HPCMode','Job','JobINI','PickHead','PyPath','StartFrom','ResFile','RunID','ResDir']
        self.boollist = ['QSub','SLiMSuite','SeqBySeq','Pickup','RjePy','SortRun','LoadBalance']
        self.intlist = ['DaisyChain','IOLimit','JobForks','KeepFree','SubSleep']
        self.numlist = ['MemFree']
        self.listlist = ['Hosts','SubJobs','OutList','Modules']
        self.dictlist = ['Running']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'Farm':'batch', 'HPC':'katana', 'HPCMode':'fork', 'Job':'', 'JobINI':'',
                     'PyPath':slimsuitepath})
        self.setBool({'QSub':False,'SLiMSuite':True,'SeqBySeq':False,'Pickup':True,'SortRun':True,'LoadBalance':True,
                      'RjePy':True})
        self.setInt({'IOLimit':50,'JobForks':0,'KeepFree':1,'SubSleep':1})
        self.setNum({})
        self.list['Modules'] = []
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self.setInt({'Forks':1})
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
                self._cmdReadList(cmd,'str',['Farm','HPC','HPCMode','Job','RunID'])   # Normal strings
                self._cmdReadList(cmd,'path',['PyPath'])  # String representing directory path
                self._cmdReadList(cmd,'file',['JobINI','ResFile'])  # String representing file path
                self._cmdReadList(cmd,'basefile',['ResFile'])  # String representing file path
                self._cmdReadList(cmd,'basepath',['ResDir'])  # String representing file path
                self._cmdReadList(cmd,'bool',['QSub','SLiMSuite','SeqBySeq','Pickup','SortRun','LoadBalance'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['DaisyChain','JobForks','KeepFree'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['Modules'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getStrLC('ResFile') and not self.getStr('ResFile').endswith('.csv'): self.str['ResFile'] += '.csv'
        self._hpcCmdList()
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [0] ~ Setup Basefile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setupBasefile()
            ### ~ [1] ~ QSub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('QSub'): return self.qSub()
            ### ~ [2] ~ SLiMFarming jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.slimFarm()
        except SystemExit: raise    # Child
        except:
            self.errorLog('Major Error during SLiMFarmer run.')
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupBasefile(self):    ### Uses Basefile and other parameters to set defaults
        '''Uses Basefile and other parameters to set defaults.'''
        try:### ~ [0] ~ Setup Basefile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Basefile'):
                if not self.getStrLC('Job'): self.baseFile('slimfarmer'); self.setStr({'Job':self.getStr('Basefile')})
                else: self.baseFile(self.getStr('Job'))
            if not self.getStrLC('Job'): self.setStr({'Job':self.basefile()})
            if not self.getStrLC('RunID'): self.setStr({'RunID':self.getStr('Job')})
            if not self.getStrLC('ResFile'): self.setStr({'ResFile':'%s.csv' % self.getStrLC('Farm')})
        except: self.errorLog('Problem during SLiMFarmer.setupBasefile().'); raise
#########################################################################################################################
    ### <3> ### QSub Class Methods                                                                                      #
#########################################################################################################################
    def qSub(self): ### Generates PBS job and places in queue with QSUB using rje_qsub.py
        '''
        Generates PBS job and places in queue with QSUB using rje_qsub.QSub. The following parameters are passed on:
        - job=X         : Name of job file (.job added) [qsub]
        - qpath=PATH    : Path to change directory too [current path]
        - nodes=X       : Number of nodes to run on [4]
        - ppn=X         : Processors per node [12]
        - walltime=X    : Walltime for qsub job (hours) [60]
        - depend=LIST   : List of job ids to wait for before starting job (dependhpc=X added) []
        - dependhpc=X   : Name of HPC system for depend ['blue30.iridis.soton.ac.uk']
        - pause=X       : Wait X seconds before attempting showstart [5]
        - report=T/F    : Pull out running job IDs and run showstart [False]
        - email=X       : Email address to email job stats to at end ['']
        - hpc=X         : Name of HPC system for depend ['IRIDIS4']

        SLiMFarmer constructs the following:
        - program=X     : Program call for Qsub (and options). Made from farm=X and slimsuite=T/F.
        - rjepy=T/F     : Whether program is an RJE *.py script (adds python PyPath/) [True]
        - pypath=PATH   : Path for RJE Python scripts [/rhome/re1u06/Serpentry/]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#JOB',self.getStr('Job'))
            qcmd = ['nodes=1','ppn=16','walltime=12','vmem=126','job=slimfarmer'] + self.cmd_list + ['rjepy=F','job=%s' % self.getStr('Job')]
            qsub = rje_qsub.QSub(self.log,qcmd)
            farm = string.split(self.getStr('Farm'))[0]
            if self.getBool('SLiMSuite') or farm == 'batch':
                if farm == 'batch':
                    program = 'python %stools/slimfarmer.py %s' % (self.getStr('PyPath'),string.join(sys.argv[1:]))
                    self.setBool({'SLiMSuite':True})
                elif self.getStrLC('HPCMode').startswith('fork') and farm not in ['slimfinder','qslimfinder','slimprob','slimcore']:
                    if farm.endswith('.py'): farm = farm[:-3]
                    if '/' in farm:
                        program = 'python %s%s.py %s' % (self.getStr('PyPath'),farm,string.join(sys.argv[1:]))
                        progpath = '%s%s.py' % (self.getStr('PyPath'),farm)
                        if not os.path.exists(progpath):
                            if self.i() < 0 or rje.yesNo('Cannot find "%s". Switching SLiMSuite=False?' % progpath):
                                self.warnLog('Cannot find "%s". Switching SLiMSuite=False.' % progpath)
                                self.setBool({'SLiMSuite':False})
                                program = self.getStr('Farm')
                    else:
                        program = 'python %stools/%s.py %s' % (self.getStr('PyPath'),farm,string.join(sys.argv[1:]))
                        progpath = '%stools/%s.py' % (self.getStr('PyPath'),farm)
                        if not os.path.exists(progpath):
                            if self.i() < 0 or rje.yesNo('Cannot find "%s". Switching SLiMSuite=False?' % progpath):
                                self.warnLog('Cannot find "%s". Switching SLiMSuite=False.' % progpath)
                                self.setBool({'SLiMSuite':False})
                                program = self.getStr('Farm')
                else: program = 'python %stools/slimfarmer.py %s' % (self.getStr('PyPath'),string.join(sys.argv[1:]))
                if self.getBool('SLiMSuite'):
                    if self.getInt('JobForks'): program += ' forks=%d qsub=F i=-1 v=-1 newlog=F' % self.getInt('JobForks')
                    else: program += ' forks=%d qsub=F i=-1 v=-1 newlog=F' % qsub.getInt('PPN')
            else: program = self.getStr('Farm')
            #x#if self.getStrLC('JobINI'): program += ' ini=%s' % self.getStr('JobINI')
            self.printLog('#QSUB',program)
            qsub.setStr({'Program':program,'PyPath':self.getStr('PyPath')})
            if self.getStrLC('HPCMode').startswith('fork'): qsub.setInt({'Nodes':1})
            while qsub.getNum('Walltime') <= 0:
                if self.i() < 0: raise ValueError('Cannot run QSub with walltime <= 0!')
                self.errorLog('Cannot run QSub with walltime <= 0!',printerror=False)
                qsub.setNum({'Walltime':rje.getFloat('New walltime (hours)?',default='12',confirm=True)})
            ### ~ [2] Run QSub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dev(): return os.system(program)
            ## ~ [2a] Special daisychain mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Chain together a set of qsub runs of the same call that depend on the previous job. [0]
            if self.getInt('DaisyChain') > 0:
                qid = ''
                for q in range(self.getInt('DaisyChain')):
                    job = '%s.%sof%d' % (self.getStr('Job'),rje.preZero(q+1,self.getInt('DaisyChain')),self.getInt('DaisyChain'))
                    qsub.setStr({'Job':job})
                    qid = qsub.qsub()
                    qsub.list['Depend'] = [qid]
                return qid
            ## ~ [2b] Normal qsub ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            return qsub.qsub()
        except: self.errorLog('SLiMFarmer.qSub() error')
#########################################################################################################################
    ### <4> ### SLiMFarmer Methods                                                                                      #
#########################################################################################################################
    def slimFarm(self): ### Farms out a job using rje_iridis.py
        '''
        Farms out a job using rje_iridis.IRIDIS. The following parameters are passed on:
        subsleep=X      : Sleep time (seconds) between cycles of subbing out jobs to hosts [1]
        subjobs=LIST    : List of subjobs to farm out to HPC cluster []
        iolimit=X       : Limit of number of IOErrors before termination [50]
        memfree=X       : Min. proportion of node memory to be free before spawning job [0.1]
        test=T/F        : Whether to produce extra output in "test" mode [False]
        keepfree=X      : Number of processors to keep free on head node [1]
        seqbyseq=T/F    : Activate seqbyseq mode - assumes basefile=X option used for output [False]
        seqin=FILE      : Input sequence file to farm out [None]
        basefile=X      : Base for output files - compiled from individual run results [None]
        outlist=LIST    : List of extensions of outputs to add to basefile for output (basefile.*) []
        pickup=X        : Header to extract from OutList file and used to populate AccNum to skip []
        runid=X         : Text identifier for iX run [SLiMFarmer]
        resfile=FILE    : Main output file for iX run [islimfinder.csv]
        sortrun=T/F     : Whether to sort input files by size and run big -> small to avoid hang at end [True]
        loadbalance=T/F : Whether to split SortRun jobs equally between large & small to avoid memory issues [True]

        SLiMFarmer constructs the following:
        irun=X          : Exectute a special iRun analysis on Iridis (gopher/slimfinder/qslimfinder/slimsearch/unifake) []
        iini=FILE       : Ini file to pass to the called program [None]
        pypath=PATH     : Path to python modules ['/home/re1u06/Serpentry/']
        rjepy=T/F       : Whether program is an RJE *.py script (adds log processing) [True]
        rsh=T/F         : Whether to use rsh to run jobs on other nodes [True]

        Key SLiMFarmer parameters:
        - farm=X          : Execute a special SLiMFarm analysis on HPC [batch]
                            - batch will farm out a batch list of commands read in from subjobs=LIST
                            - gopher/slimfinder/qslimfinder/slimsearch/unifake = special SLiMSuite HPC runs
                            - if seqbyseq=T, farm=X will specify the program to be run (see docs)
                            - otherwise, farm=X will be executed as a system call in place of SLiMFarmer unless
        - hpcmode=X       : Mode to be used for farming jobs between nodes (rsh/forks) [forks]
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('HPCMode').startswith('fork'):
                self.printLog('#FORK','HPC fork Mode: %d forks.' % self.getInt('Forks'))
                if not self.getInt('Forks'): raise ValueError('Cannot run HPC fork mode with 0 forks!')
            self.setup()
            if self.getBool('SeqBySeq'): return self.seqBySeq()
            farm = self.getStr('Farm')
            if farm in ['slimfinder','qslimfinder','slimprob','slimcore']: return self.farmSLiMJobs()
            ## ~ [1a] Setup rje_iridis.IRIDIS commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hpcmd = ['rjepy=T','rsh=F']     # Set up an extended cmdlist which will be propagated to special runs
            if farm in ['gopher','slimsearch','unifake']:   #!# Add these to SLiMFarmer
                hpcmd += ['irun=%s' % farm,'pypath=%stools/' % self.getStr('PyPath'),'iini=%s' % self.getStr('JobINI')]
                if self.getInt('JobForks'): hpcmd += ['forks=%d' % (self.getInt('JobForks'))]
            elif farm == 'batch': hpcmd[0] = 'rjepy=F'
            else:
                self.printLog('#SYS',farm)
                return os.system(farm)
            if self.getStrLC('HPCMode') == 'rsh': hpcmd[1] = 'rsh=T'
            if self.getStrLC('HPCMode').startswith('fork'): hpcmd.append('forks=%d' % self.getInt('Forks'))
            hpcmd += ['runid=%s' % self.getStr('RunID'), 'resfile=%s' % self.getStr('ResFile')]
            self.printLog('#RUNID',self.getStr('RunID'))
            self.printLog('#RES',self.getStr('ResFile'))
            self.deBug(hpcmd)
            ### ~ [2] Run HPC Controller ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hpc = rje_iridis.IRIDIS(self.log,self.cmd_list+hpcmd)
            hpc.iridisRun()
        except: self.errorLog('SLiMFarmer.slimFarm() error')
#########################################################################################################################
    def farmSLiMJobs(self): ### Farms out SLiMSuite jobs using modified rje_hpc methods
        '''Farms out SLiMSuite jobs using modified rje_hpc methods.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Extra Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            farm = self.getStr('Farm')
            self.setBool({'SLiMCore':farm == 'slimcore'})
            self.list['JRan'] = []      # List of random job numbers for checking and cleanup
            self.list['Batch'] = []     # List of files for batch running.
            ## ~ [0b] ~ Extra Commandline Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for cmd in self.cmd_list:
                try: self._cmdReadList(cmd,'list',['Batch'])
                except: self.errorLog('Problem with cmd:%s' % cmd)
            if not self.list['Batch']: return self.printLog('#EXIT','No batch=LIST for %s farming.' % farm)
            ### ~ [1] ~ Pickup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#RUNID','RunID: %s' % self.getStr('RunID'))
            self.printLog('#RES','ResFile: %s' % self.getStr('ResFile'))
            if self.getStrLC('ResDir'): self.printLog('#RES','ResDir: %s' % self.getStr('ResDir'))
            pickup = []
            try:
                if self.getBool('Pickup') and not self.getBool('SLiMCore') and rje.checkForFile(self.getStr('ResFile')):
                    try: pickup = rje.dataDict(self,self.getStr('ResFile'),['RunID'],['Dataset'],lists=True)
                    except: pickdat = {}
                    if self.getStr('RunID') in pickup:
                        pickup = pickup[self.getStr('RunID')]['Dataset']
                    else:
                        self.printLog('#PICKUP','RunID %s not found: %s' % (self.getStr('RunID'),rje.sortKeys(pickup)))
                        pickup = []
                    self.printLog('#PICKUP','%d %s datasets to skip' % (len(pickup),self.getStr('RunID')))
            except: self.errorLog('Pickup Failed'); pickup = []
            ### ~ [1] ~ Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
            self.printLog('\r#BATCH','%s batch files to run.' % rje.iLen(batchfiles))
            sizedict = {}
            newbatch = []       # New list of batch files in size order, starting with the biggest
            for file in batchfiles:
                if os.path.split(rje.baseFile(file))[1] in pickup:
                    self.printLog('#PICKUP','Skipping %s batch file %s' % (self.getStr('RunID'),file))
                    continue
                if self.getBool('SortRun'):
                    fsize = os.path.getsize(file)
                    if fsize not in sizedict: sizedict[fsize] = []
                    sizedict[fsize].append(file)
                else: newbatch.append(file)
            if self.getBool('SortRun'):
                sizes = rje.sortKeys(sizedict,revsort=True)
                for fsize in sizes: newbatch += sizedict.pop(fsize)
            self.list['NewBatch'] = newbatch
            self.int['BatchNum'] = len(self.list['NewBatch'])
            self.printLog('#JOBS','%s Batch jobs to run.' % rje.iStr(self.int['BatchNum']))
            ### ~ [2] ~ Run Jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.runSLiMJobs()
            self.cleanUpSLiMJobs()
            return True
        except SystemExit: raise    # Child
        except: self.errorLog('Problem with SLiMFarmer.farmSLiMJobs()'); self.cleanUpSLiMJobs()
        return False
#########################################################################################################################
    def runSLiMJobs(self):  ### Runs all the jobs in self.list['SubJobs']
        '''Runs all the jobs in self.list['SubJobs'].'''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkn = self.nprocs() - self.int['KeepFree']
        funny5min = False; big_problem = False
        nodes = rje.sortUnique(self.list['Hosts'])
        self.int['PPN'] = ppn = self.nprocs() / len(nodes)
        jlist = range(self.int['KeepFree'],self.nprocs())
        self.printLog('#PPN','%d ppn on %d host nodes (%d jobs + %d head)' % (ppn,len(nodes),len(jlist),self.int['KeepFree']))
        #for j in range(1,self.nprocs()): self.nextSLiMJob(j)    # Skip first node
        for p in range(ppn):
            for n in range(len(nodes)):
                j = n * ppn + p
                if j in jlist: self.nextSLiMJob(j)    # Skip first node
            time.sleep(self.int['SubSleep'])
        pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.dict['Running'] or self.list['NewBatch']:
            pidchecklist = []
            for j in rje.sortKeys(self.dict['Running']):
                if not self.dict['Running'][j]: self.dict['Running'].pop(j); continue   # No more jobs
                pid = '????'
                try:
                    pid = self.dict['Running'][j]['PID']
                    if string.split('%s' % pid)[0] == 'WAIT':
                        pidchecklist.append('%s: %s' % (j,pid))
                        status = 1
                    else:
                        pidchecklist.append('%s: %s (%s) = %s' % (j,pid,self.dict['Running'][j]['JRan'],rje.hhmmss(time.time() - self.dict['Running'][j]['Start'])))
                        (status,exit_stat) = os.waitpid(pid,os.WNOHANG)
                except:
                    status = 1
                    self.errorLog('Job %s problem.' % j)
                    try: self.printLog('#ERR','Job %s (pid:%s) problem after %s. (%s.*)' % (j,pid,rje.hhmmss(time.time() - self.dict['Running'][j]['Start']),self.dict['Running'][j]['JRan']))
                    except: pass
                if status > 0: self.endSLiMJob(j)       # subjob on processor j has finished: can replace with processing
            PIDCHECK = open(pidcheck,'w')
            PIDCHECK.write(string.join(pidchecklist+[''],'\n'))
            PIDCHECK.close()
            if self.list['NewBatch'] and len(self.dict['Running']) < checkn:    ### Why are jobs missing?
                if big_problem:
                    self.printLog('#ERR','Still only %d of %d jobs running but %d remain. Self-destructing.' % (len(self.dict['Running']),checkn,len(self.list['NewBatch'])))
                    self.cleanUpSLiMJobs()
                    raise ValueError
                elif funny5min:
                    self.printLog('#JOBS','Still only %d of %d jobs running but %d remain. Will try starting some more.' % (len(self.dict['Running']),checkn,len(self.list['NewBatch'])))
                    for j in range(1,self.nprocs()):
                        if j not in self.dict['Running']: self.nextSLiMJob(j)    # Skip first node
                    big_problem = True
                else:
                    self.printLog('#ERR','Only %d of %d jobs running but %d remain. Funny 5 minutes?' % (len(self.dict['Running']),checkn,len(self.list['NewBatch'])))
                    funny5min = True
                    time.sleep(300)
                    continue
            else: funny5min = False; big_problem = False
            time.sleep(self.int['SubSleep'])
#########################################################################################################################
    def cleanUpSLiMJobs(self):  ### Tries to clean up random run files that are still hanging about
        '''Tries to clean up random run files that are still hanging about.'''
        try:### ~ [1] Remove jobs without files from self.list['JRan'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            precleanx = len(self.list['JRan'])
            for jran in self.list['JRan'][0:]:
                jlog = '%s%s.log' % (self.str['RunPath'],jran)
                jcsv = '%s%s.csv' % (self.str['RunPath'],jran)
                if not os.path.exists(jlog) and not os.path.exists(jcsv): self.list['JRan'].remove(jran)
            self.printLog('#CLEAN','%s of %s job runs have files left' % (rje.iLen(self.list['JRan']),rje.iStr(precleanx)))
            if not self.list['JRan']: return
            ### ~ [2] Try cleaning up the files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tranx = 0; pretranx = len(self.list['JRan'])
            for jran in self.list['JRan'][0:]:
                jlog = '%s%s.log' % (self.getStr('RunPath'),jran)
                jcsv = '%s%s.csv' % (self.str['RunPath'],jran)
                jocc = '%s%s.occ.csv' % (self.str['RunPath'],jran)
                socc = '%s.occ.csv' % self.getStr('ResFile')[:-4]
                if os.path.exists(jlog):
                    jloglines = open(jlog,'r').readlines()
                    if jloglines:
                        if string.split(jloglines[-1])[0] == '#LOG' and 'End:' in string.split(jloglines[-1]):
                            self.printLog('#END','Random job %s looks finished' % jran)
                            if os.path.exists(jcsv):
                                if os.path.exists(self.getStr('ResFile')): open(self.getStr('ResFile'),'a').writelines(open(jcsv,'r').readlines()[1:])
                                else: open(self.getStr('ResFile'),'w').writelines(open(jcsv,'r').readlines())     # Headers too
                            if os.path.exists(jocc):
                                if os.path.exists(socc): open(socc,'a').writelines(open(jocc,'r').readlines()[1:])
                                else: open(socc,'w').writelines(open(jocc,'r').readlines())     # Headers too
                            loglines = open(jlog,'r').readlines()[4:]
                            if len(loglines) > 1:
                                if loglines[1].startswith('#VIO'): loglines = loglines[:1] + loglines[2:]
                                if loglines[-1].startswith('#LOG'): loglines = loglines[:-1]
                                open(self.log.info['LogFile'],'a').writelines(loglines)
                            else: self.printLog('#ERR','%s log content error' % jran)
                            self.printLog('#END','%s results & log content transferred' % jran); tranx += 1
                        else: self.printLog('#ERR','%s only reached: %s' % (jran,jloglines[-1]))
                    os.unlink(jlog)
                    self.printLog('#DEL','%s deleted in cleanup' % jlog)
                if os.path.exists(jcsv):
                    os.unlink(jcsv)
                    self.printLog('#DEL','%s deleted in cleanup' % jcsv)
            self.printLog('#CLEAN','%s of %s job runs had results & logs transferred' % (rje.iStr(tranx),rje.iLen(self.list['JRan'])))
        except: self.errorLog('SLiMFarmer.cleanUpSLiMJobs error. %s runs not checked.' % rje.iLen(self.list['JRan']))
#########################################################################################################################
    def endSLiMJob(self,host_id):   ### Replaced endJob to set a new SLiMFinder going
        '''Replaced endJob to set a new SLiMSuite job going.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            waiting = False
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                if 'Log' in jdict:
                    if os.path.exists(jdict['Res']):
                        if os.path.exists(self.getStr('ResFile')): open(self.getStr('ResFile'),'a').writelines(open(jdict['Res'],'r').readlines()[1:])
                        else: open(self.getStr('ResFile'),'w').writelines(open(jdict['Res'],'r').readlines())     # Headers too
                        os.unlink(jdict['Res'])
                    jfile = '%s.occ.csv' % jdict['Res'][:-4]
                    sfile = '%s.occ.csv' % self.getStr('ResFile')[:-4]
                    if os.path.exists(jfile):
                        if os.path.exists(sfile): open(sfile,'a').writelines(open(jfile,'r').readlines()[1:])
                        else: open(sfile,'w').writelines(open(jfile,'r').readlines())     # Headers too
                        os.unlink(jfile)
                    etxt = 'Job on processor %d ended' % host_id
                    dset = os.path.basename(jdict['Next'])
                    if os.path.exists(jdict['Log']):
                        open(self.log.info['LogFile'],'a').writelines(open(jdict['Log'],'r').readlines()[5:-1])
                        etxt = '%s: %s results & log content transferred' % (etxt,dset)
                        os.unlink(jdict['Log'])
                    else:
                        etxt = '%s: %s (%s) failure! Will rerun %s.' % (etxt,dset,jdict['JRan'],jdict['Next'])
                        self.list['NewBatch'].insert(0,jdict['Next'])
                        time.sleep(self.int['SubSleep'])
                    self.printLog('#END',etxt)
                    self.printLog('#~~#','#~~#',timeout=False)
                #x#if 'DAT' in jdict: os.system('rm %s*' % jdict['DAT'])
                elif 'PID' in jdict and string.split('%s' % jdict['PID'])[0] == 'WAIT': waiting = True
                else: self.printLog('#END','Job on processor %d ended.' % host_id)
        except IOError:
            if self.int['IOError'] == 1: self.errorLog('iSLiMFinder.endSLiMJob IOError limit reached'); raise
            else: self.int['IOError'] -= 1; self.errorLog('iSLiMFinder.endSLiMJob')
        except: self.errorLog('iSLiMFinder.endSLiMJob error')
        self.nextSLiMJob(host_id)   # Carry on regardless
        if waiting:
            try:
                if string.split('%s' % self.dict['Running'][host_id]['PID'])[0] == 'WAIT':
                    if self.test(): self.printLog('#WAIT','Resumed after waiting for memory. %d sec sleep!' % self.int['SubSleep'])
                    time.sleep(self.int['SubSleep'])
            except:
                if self.opt['Test']: self.printLog('#WAIT','Resumed after waiting for memory but no new PID!')
#########################################################################################################################
    def nextSLiMJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            node = self.list['Hosts'][host_id]
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            freemem = self.freeMem(node)
            if self.num['MemFree'] > freemem:
                jdict['PID'] = 'WAIT - %.1f%% %s mem (host %d)' % (freemem*100.0,node,host_id)
                if self.test(): self.printLog('#JOB',jdict['PID'])
                return
            if self.list['NewBatch']:
                if self.nprocs() != self.int['PPN'] and host_id < self.int['PPN']: next = self.list['NewBatch'].pop(-1)     # Keep head node running small jobs
                elif self.bool['SortRun'] and self.bool['LoadBalance'] and not len(self.list['NewBatch']) % 2: next = self.list['NewBatch'].pop(-1)
                else: next = self.list['NewBatch'].pop(0)   # SLiMFinder this dataset
            else: return                                                    # Out of datasets: stop
            jran = 'i_%s_%s' % (host_id,rje.randomString(6))
            while jran in self.list['JRan'] or os.path.exists('%s%s.log' % (self.str['RunPath'],jran)):
                jran = 'i_%s_%s' % (host_id,rje.randomString(6))
            self.list['JRan'].append(jran)
            jdict['JRan'] = jran; jdict['Next'] = next
            jdict['Log'] = '%s%s.log' % (self.str['RunPath'],jran)
            jdict['Res'] = '%s%s.csv' % (self.str['RunPath'],jran)
            jdict['Start'] = time.time()
            try: open(jdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endSLiMJob(host_id)
            #x#initial_cmds = 'cd ' + self.str['RunPath'] + ' ; echo %s as %s on `hostname` ; setenv IUPred_PATH /home/re1u06/Bioware/iupred/ ;' % (next,jran)
            initial_cmds = 'cd ' + self.str['RunPath'] + ' ; echo %s as %s on `hostname` ;' % (next,jran)  # bash: setenv: command not found
            if next[-3:] == 'acc':
                jdict['DAT'] = '%sdat' % next[:-3]
                while os.path.exists(jdict['DAT']):   # No need to run - skip onto the next one
                    if self.list['NewBatch']: next = self.list['NewBatch'].pop(0)
                    else: return
                    jdict['DAT'] = '%sdat' % next[:-3]
                    initial_cmds = 'cd ' + self.str['RunPath'] + ' ; echo %s as %s on `hostname` ; ' % (next,jran)
                job = 'python %srje_uniprot.py extract=%s datout=%s i=-1 v=-1 log=%s' % (self.str['PyPath'],next,jdict['DAT'],jdict['Log'])
                #initial_cmds = '%s python %srje_uniprot.py extract=%s datout=%s i=-1 v=-1 log=%s ; ' % (initial_cmds,self.info['PyPath'],next,jdict['DAT'],jdict['Log'])
                #next = jdict['DAT']
            else:
                if self.getBool('SLiMCore'): job = 'python %slibraries/rje_slimcore.py %s' % (self.str['PyPath'],string.join(sys.argv[1:]))
                else: job = 'python %stools/%s.py %s' % (self.str['PyPath'],self.str['Farm'],string.join(sys.argv[1:]))
                if self.str['JobINI']: job = '%s ini=%s' % (job,self.str['JobINI'])
                job = '%s pickup=F batch= basefile= seqin=%s i=-1 v=-1 runid=%s' % (job,next,self.getStr('RunID'))
                if self.getStrLC('ResDir'): job = '%s resdir=%s' % (job,self.getStr('ResDir'))
                job = '%s resfile=%s log=%s' % (job,jdict['Res'],jdict['Log'])
                if self.getInt('JobForks'): job = '%s forks=%d' % (job,self.getInt('JobForks'))
                if self.str['Farm'].startswith('slim'): job = '%s noforks=T' % job
            rsh = "rsh %s '%s %s'" % (self.list['Hosts'][host_id],initial_cmds,job)
            if self.getBool('RSH'): self.printLog('#RSH',rsh)
            else: self.printLog('#FORK',job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                remaining = len(self.list['NewBatch']); runnum = self.int['BatchNum'] - remaining
                if self.getNum('MemFree') >= 0.0: freemem = '; %.1f%% mem free.' % (freemem*100.0)
                else: freemem = '.'
                self.printLog('#BATCH','Running %s as %s (pid:%s) [%d::%s]: %d run; %d remain%s' % (os.path.basename(next),jran,cpid,host_id,self.list['Hosts'][host_id],runnum,remaining,freemem))
            else:                   # child process
                if self.getBool('RSH'): os.system(rsh)
                else: os.system(job)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('SLiMFarmer.nextSLiMJob error')
#########################################################################################################################
### End of SECTION II: SLiMFarmer Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def iRun(irun,itype): return rje_iridis.iRun(irun,itype)
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
    try: SLiMFarmer(mainlog,cmd_list).run()

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
