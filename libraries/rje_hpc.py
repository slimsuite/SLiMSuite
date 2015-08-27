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
Module:       rje_hpc
Description:  High Performance Computing job farming
Version:      1.1
Last Edit:    17/06/14
Copyright (C) 2014  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains the generic code for farming jobs out over multiple processors and/or nodes. It is based on
    rje_iridis Version 1.10, updated to the new RJE_Object structure and made more generic. The recommended tool for
    actual job farming is SLiMFarmer, which inherits the rje_hpc.JobFarmer object and adds extra bells and whistles.

Commandline:
    ### ~ Generic HPC Job Farming pptions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    subsleep=X      : Sleep time (seconds) between cycles of subbing out jobs to hosts [1]
    subjobs=LIST    : List of subjobs to farm out to HPC cluster []
    hpcmode=X       : Mode to be used for farming jobs between nodes (rsh/fork) [fork]
    iolimit=X       : Limit of number of IOErrors before termination [50]
    memfree=X       : Min. proportion of node memory to be free before spawning job [0.0]
    keepfree=X      : Number of processors to keep free on head node [1]

    ### ~ SLiMSuite and SeqSuite SeqBySeq program options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqbyseq=T/F    : Activate seqbyseq mode - assumes basefile=X option used for output [False]
    farm=X          : Program to farm out using seqbyseq mode. []
    seqin=FILE      : Input sequence file to farm out [None]
    basefile=X      : Base for output files - compiled from individual run results [None]
    outlist=LIST    : List of extensions of outputs to add to basefile for output (basefile.*) []
    pickhead=X      : Header to extract from OutList file and used to populate AccNum to skip []
    startfrom=X     : Sequence ID at which to begin the SeqBySeq farming [None]
    rjepy=T/F       : Whether program is an RJE *.py script (adds log processing) [True]
    jobini=FILE     : Ini file to pass to the called program [None]
    pypath=PATH     : Path to python modules ['/home/re1u06/Serpentry/']

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
import rje, rje_obj, rje_seq
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - Initial Compilation based on rje_iridis V1.10.
    # 1.1 - Disabled memory checking in Windows and OSX.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [Y] : Update check memory functions to check local memory when forking.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_HPC', '1.0', 'February 2014', '2014')
    description = 'High Performance Computing job farming'
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
### SECTION II: JobFarmer Class                                                                                         #
#########################################################################################################################
class JobFarmer(rje_obj.RJE_Object):
    '''
    JobFarmer Class. Author: Rich Edwards (2014).

    Str:str
    - Farm = Program to farm out using seqbyseq mode. [None]
    - HPCMode = Mode to be used for farming jobs between nodes (rsh/fork) [fork]
    - JobINI = Ini file to pass to the called program [None]
    - PickHead = Header to extract from OutList file and used to populate AccNum to skip []
    - PyPath = path to python modules [slimsuite home directory]
    - StartFrom = sequence ID at which to begin the SeqBySeq farming [None]

    Bool:boolean
    - RjePy = Whether program is an RJE *.py script (adds log processing) [True]
    - SeqBySeq = Activate seqbyseq mode - assumes basefile=X option used for output [False]

    Int:integer
    - IOLimit = Limit of number of IOErrors before termination [50]
    - KeepFree = Number of processors to keep free on head node [1]
    - SubSleep = Sleep time (seconds) between cycles of subbing out jobs to hosts [1]

    Num:float
    - MemFree = Min. proportion of node memory to be free before spawning job [0.0]

    List:list
    - Hosts = Host processors
    - OutList = List of extensions of outputs to add to basefile for output []
    - SubJobs = List of jobs to farms out to nodes

    Dict:dictionary
    - Running = List of running subjobs (Dictionaries of {HostIndex:{'PID'/'Job'/'Log'}})

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Farm','HPCMode','JobINI','PickHead','PyPath','StartFrom']
        self.boollist = ['RjePy','SeqBySeq']
        self.intlist = ['IOLimit','KeepFree','SubSleep']
        self.numlist = ['MemFree']
        self.listlist = ['Hosts','SubJobs','OutList']
        self.dictlist = ['Running']
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'HPCMode':'fork','PyPath':slimsuitepath})
        self.setBool({'RjePy':True})
        self.setInt({'IOLimit':50,'KeepFree':1,'SubSleep':1})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self): self._hpcCmdList()  ### Sets Attributes from commandline
#########################################################################################################################
    def _hpcCmdList(self):     ### Sets Attributes from commandline. Kept separate to permit inheritance.
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
                self._cmdReadList(cmd,'str',['Farm','HPCMode','JobINI','PickHead','StartFrom'])   # Normal strings
                self._cmdReadList(cmd,'path',['PyPath'])  # String representing directory path
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path
                self._cmdReadList(cmd,'bool',['RjePy','SeqBySeq'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['IOLimit','KeepFree','SubSleep'])   # Integers
                self._cmdReadList(cmd,'float',['MemFree']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if 'JobINI' not in self.str or not self.getStrLC('JobINI'): self.setStr({'JobINI':''})
        if self.getNum('MemFree') > 1: self.num['MemFree'] /= 100.0
        self.setStr({'HPCMode':self.getStrLC('HPCMode')[:4]})
#########################################################################################################################
    def nprocs(self): return len(self.list['Hosts'])                                                                #V1.0
    def rsh(self): return self.getStr('HPCMode') == 'rsh'                                                           #V1.0
    def forking(self): return self.getStr('HPCMode') == 'fork' and self.getInt('Forks')                             #V1.0
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method                                                                             #V1.0
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ SeqBySeq Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SeqBySeq'): return self.seqBySeq()
            ### ~ [2] ~ Run Jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.test(): self.set_subjobs()
            self.printLog('#JOBS','Running %d subjobs on %d hosts' % (len(self.list['SubJobs']),self.nprocs()))
            self.runJobs()
        except SystemExit: raise    # Child exit
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                                #V1.0
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Set job directory to RunPath if given, else directory from which job was submitted ~ ##
            try: jobdir = rje.makePath(os.environ['PBS_O_WORKDIR'])
            except: jobdir = None
            if self.getStr('RunPath') == rje.makePath(os.path.abspath(os.curdir)) and jobdir: self.setStr({'RunPath':jobdir})
            os.chdir(self.getStr('RunPath'))
            ## ~ [1b] ~ Read list of node names in file $PBS_NODEFILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setHosts()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def setHosts(self): ### Sets up self.list['Hosts']                                                              #V1.0
        '''Sets up self.list['Hosts'].'''
        try: self.list['Hosts'] = self.loadFromFile(os.environ['PBS_NODEFILE'],chomplines=True)
        except: self.list['Hosts'] = []; self.printLog('#HPC','HPC PBS not detected.')
        if not self.list['Hosts']:
            if self.forking():
                self.list['Hosts'] = ['fork'] * self.getInt('Forks')
            else:
                self.errorLog('No Hosts identified for HPC run',printerror=False)
                sys.exit()
        self.printLog('#HOSTS','%d hosts: %s.' % (self.nprocs(),string.join(self.list['Hosts'],'; ')))
#########################################################################################################################
    def checkMemory(self,node): ### Checks free memory on node and returns True/False if OK to run                  #V1.0
        '''Checks free memory on node and returns True/False if OK to run.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('MemFree') <= 0.0: return True
            ### ~ [2] Check details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.freeMem(node) < self.getNum('MemFree'): return False
            return True
        except: self.errorLog('Problem with checkMemory(%s)' % node); return False
#########################################################################################################################
    def freeMem(self,node): ### Uses RSH to read free memory from node.                                             #V1.0
        '''Uses RSH to read free memory from node.'''
        if self.getNum('MemFree') < 0.0: return 0.0
        if self.getBool('Win32') or self.getBool('OSX'): return 0.0
        if self.rsh():
            try: memdata = string.split(os.popen('rsh %s free' % node).readlines()[1])
            except: self.errorLog('Unable to read %s free memory' % node,printerror=False); return 0.0
        elif self.getStr('HPCMode') == 'fork':
            try: memdata = string.split(os.popen('free').readlines()[1])
            except: self.errorLog('Unable to read free memory',printerror=False); return 0.0
        else: raise ValueError('HPC Mode "%s" not recognised!' % self.getStr('HPCMode'))
        return string.atof(memdata[3]) / string.atof(memdata[1])
#########################################################################################################################
    ### <3> ### Job Farming Methods                                                                                     #
#########################################################################################################################
    def runJobs(self):  ### Runs all the jobs in self.list['SubJobs']                                               #V1.0
        '''Runs all the jobs in self.list['SubJobs'].'''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for j in range(self.getInt('KeepFree'),self.nprocs()): self.nextJob(j)    # Skip first node(s)
        pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.dict['Running']:
            PIDCHECK = open(pidcheck,'w')
            for j in rje.sortKeys(self.dict['Running']):
                if not self.dict['Running'][j]: self.dict['Running'].pop(j); continue   # No more jobs
                try:
                    pid = self.dict['Running'][j]['PID']
                    PIDCHECK.write('%s: %s\n' % (j,pid))
                    if string.split('%s' % pid)[0] == 'WAIT': status = 1
                    else: (status,exit_stat) = os.waitpid(pid,os.WNOHANG)
                except: status = 1
                if status > 0: self.endJob(j)       # subjob on processor j has finished: can replace with processing
            PIDCHECK.close()
            time.sleep(self.getInt('SubSleep'))
#########################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index                                #V1.0
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            node = self.list['Hosts'][host_id]
            freemem = self.freeMem(node)
            if self.getBool('SeqBySeq'): return self.nextSeqJob(host_id)
            jdict = self.dict['Running'][host_id] = {}      # Setup empty dictionary to fill, if jobs available
            if self.getNum('MemFree') > freemem:
                jdict['PID'] = 'WAIT - %.1f%% %s mem' % (freemem*100.0,node)
                return
            if self.list['SubJobs']: job = self.list['SubJobs'].pop(0)
            else: return
            if self.getBool('RjePy'): jdict['Log'] = '%si_%s.log' % (self.getStr('RunPath'),rje.randomString(6))
            if self.getStr('JobINI'): job = '%s ini=%s' % (job,self.getStr('JobINI'))
            if 'Log' in jdict:
                job = '%s log=%s' % (job,jdict['Log'])
                try: open(jdict['Log'],'w')
                except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#JOB','Running job as %s: %d remain; %.1f%% mem free' % (cpid,len(self.list['SubJobs']),freemem*100.0))
            else:                   # child process
                if self.rsh(): os.system(rsh)
                else: os.system(job)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('JobFarmer.nextJob error')
#########################################################################################################################
    def endJob(self,host_id):   ### Ends job on host and sets new one running
        '''Ends job on host and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('SeqBySeq'): return self.endSeqJob(host_id)
            jdict = self.dict['Running'].pop(host_id)
            if 'Log' in jdict:
                open(self.log.info['LogFile'],'a').writelines(open(jdict['Log'],'r').readlines()[5:-1])
                self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                self.printLog('#~~#','#~~#',timeout=False)
                os.unlink(jdict['Log'])
            elif 'PID' in jdict and string.split('%s' % jdict['PID'])[0] == 'WAIT': pass
            else: self.printLog('#END','Job on processor %d ended.' % host_id)
        except IOError:
            if self.getInt('IOError') == 1: self.errorLog('JobFarmer.endJob IOError limit reached'); raise
            else: self.int['IOError'] -= 1; self.errorLog('JobFarmer.endJob')
        except: self.errorLog('JobFarmer.endJob error')
        self.nextJob(host_id)   # Carry on regardless
#########################################################################################################################
    ### <4> ### SeqBySeq Run Methods                                                                                    #
#########################################################################################################################
    def seqBySeq(self):     ### Runs in SeqBySeq Mode                                                               #V1.0
        '''
        In SeqBySeq mode, the program assumes that seqin=FILE and basefile=X are given and farm states the program to be run.
        Seqin will then be worked through in turn and each sequence farmed out to the farm program. Outputs given by OutList
        are then compiled, as is the Log, into the correct basefile=X given. In the case of *.csv and *.tdt files, the header
        row is copied for the first file and then excluded for all subsequent files. For all other files extensions, the
        whole output is copied.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('Farm')[-3:] == '.py': self.str['Farm'] = self.str['Farm'][:-3]
            self.list['Seq'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','accnr=F','seqnr=F']).seq[0:]
            while self.getStrLC('StartFrom') and self.list['Seq']:
                if self.list['Seq'][0].shortName() != self.getStr('StartFrom'): self.list['Seq'] = self.list['Seq'][1:]
                else: self.str['StartFrom'] = ''
            self.printLog('#SEQ','%s query sequences to farm out' % rje.integerString(len(self.list['Seq'])))
            self.list['Pickup'] = self.pickupList()
            ### ~ [2] ~ Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.runJobs()
            return True
        except SystemExit: raise    # Child
        except: self.errorLog('JobFarmer.seqBySeq error')
        return False
#########################################################################################################################
    def nextSeqJob(self,host_id):  ### Sets an new job running on host with given index                             #V1.0
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.list['Seq']: seq = self.list['Seq'].pop(0)  # UniFake this sequence
            else: return                                        # Out of sequences: stop
            if seq.info['AccNum'] in self.list['Pickup']: return self.nextSeqJob(host_id)   # Skip this sequence
            jran = 'i_%s' % rje.randomString(6)
            jdict['Log'] = '%s%s.log' % (self.getStr('RunPath'),jran)
            jdict['Qry'] = '%s%s.qry' % (self.getStr('RunPath'),jran)
            for out in self.list['OutList']: jdict[out] = '%s%s.%s' % (self.getStr('RunPath'),jran,out)
            open(jdict['Qry'],'w').write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
            job = 'python %s%s.py' % (self.getStr('PyPath'),self.getStr('Farm'))
            if self.getStr('JobINI'): job = '%s ini=%s' % (job,self.getStr('JobINI'))
            job = '%s seqin=%s i=-1 v=-1 basefile=%s' % (job,jdict['Qry'],jran)
            initial_cmds = 'cd ' + self.getStr('RunPath') + ' ; echo %s on `hostname` as %s ; ' % (seq.shortName(),jran)
            if self.rsh():
                job = '%s %s ; echo Finishing on `hostname`' % (initial_cmds,job)
                job = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: cpid = os.fork()        # Fork child process
            except OSError:
                self.errorLog('JobFarmer.nextJob error. Will sleep for 10 mins then continue.')
                time.sleep(600)
                self.printLog('#YAWN','Re-awakening JobFarmer nextJob. Fingers crossed.')
                jdict['Error'] = 'JobFarmer.nextJob OSError.'     #!# Make sure this node is retried.
                return
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#SEQ','Running %s as %s [%d::%s]: %d remain' % (seq.shortName(),cpid,host_id,self.list['Hosts'][host_id],len(self.list['Seq'])))
            else:                   # child process
                os.system(job)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('JobFarmer.nextSeqJob error')
#########################################################################################################################
    def endSeqJob(self,host_id):   ### Ends job on host and sets new one running                                    #V1.0
        '''Ends job on host and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                for out in self.list['OutList']:
                    resfile = '%s.%s' % (self.baseFile(),out)
                    if os.path.exists(jdict[out]):
                        if os.path.exists(resfile) and string.split(resfile,'.')[-1] in ['tdt','csv']:
                            open(resfile,'a').writelines(open(jdict[out],'r').readlines()[1:])
                        else:
                            open(resfile,'a').writelines(open(jdict[out],'r').readlines())
                        os.unlink(jdict[out])
                    else: self.errorLog('Output %s file "%s" missing' % (out,jdict[out]),printerror=False)
                if os.path.exists(jdict['Log']):
                    loglines = open(jdict['Log'],'r').readlines()[5:-1]
                    if loglines: open(self.log.info['LogFile'],'a').writelines(loglines)
                    self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                    self.printLog('#~~#','#~~#',timeout=False)
                    os.unlink(jdict['Log'])
                elif 'PID' in jdict and string.split('%s' % jdict['PID'])[0] == 'WAIT': pass
                else: self.printLog('#END','Job on processor %d ended: already complete' % host_id)
                os.system('rm %s*' % jdict['Qry'])
        except: self.errorLog('JobFarmer.endJob error')
        self.nextJob(host_id)   # Carry on regardless
#########################################################################################################################
    def pickupList(self):   ### Generates Pickup List from file(s)                                                  #V1.0
        '''Generates Pickup List from file(s).'''
        if not self.getStrLC('PickHead').lower(): return []
        picklist = []
        for out in self.list['OutList']:
            resfile = '%s.%s' % (self.baseFile(),out)
            try: pickdat = rje.dataDict(self,resfile,[self.getStr('PickHead')])
            except: pickdat = {}
            picklist = picklist + pickdat.keys()
        return picklist
#########################################################################################################################
    ### <5> ### Testing Methods                                                                                         #
#########################################################################################################################
    def set_subjobs(self):  ### Sets temporary test list of subjobs to run                                          #V1.0
        '''Sets up list of subjob commands to be executed in array subjobs.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['SubJobs']: return     # Already given jobs, so just carry on with these
            ### ~ [2] ~ Make a test set of Zen jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jobnum = 2 * self.nprocs() + 3  # Default number of jobs
            initial_cmds = 'cd ' + self.getStr('RunPath') + ' ; echo `pwd` on `hostname` ; '
            for nj in range(jobnum):
                self.list['SubJobs'].append('%s python %slibraries/rje_zen.py wisdoms=%d' % (initial_cmds,slimsuitepath,random.randint(1,20)))
        except: self.errorLog('JobFarmer.set_subjobs error')
#########################################################################################################################
### End of SECTION II: JobFarmer Class                                                                                  #
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
    try: JobFarmer(mainlog,cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
