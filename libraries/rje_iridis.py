#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2008 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       RJE_IRIDIS
Description:  Parallel processing on IRIDIS
Version:      1.10.2
Last Edit:    29/11/15
Copyright (C) 2008  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to control and execute parallel processing jobs on the IRIDIS cluster based on the script
    written by Ivan Wolton. Initially, it will call other programs but, in time, it is envisaged that other programs will
    make use of this module and have parallelisation built-in.

    In SeqBySeq mode, the program assumes that seqin=FILE and basefile=X are given and irun states the program to be run.
    Seqin will then be worked through in turn and each sequence farmed out to the irun program. Outputs given by OutList
    are then compiled, as is the Log, into the correct basefile=X given. In the case of *.csv and *.tdt files, the header
    row is copied for the first file and then excluded for all subsequent files. For all other files extensions, the
    whole output is copied.

Commandline:
    ### ~ STANDARD RUN OPTIONS ~ ###
    irun=X          : Exectute a special iRun analysis on Iridis (gopher/slimfinder/qslimfinder/slimsearch/unifake) []
    iini=FILE       : Ini file to pass to the called program [None]
    pypath=PATH     : Path to python modules ['/home/re1u06/Serpentry/']
    rjepy=T/F       : Whether program is an RJE *.py script (adds log processing) [True]
    subsleep=X      : Sleep time (seconds) between cycles of subbing out jobs to hosts [1]
    subjobs=LIST    : List of subjobs to farm out to IRIDIS cluster []
    iolimit=X       : Limit of number of IOErrors before termination [50]
    memfree=X       : Min. proportion of node memory to be free before spawning job [0.0]
    test=T/F        : Whether to produce extra output in "test" mode [False]
    keepfree=X      : Number of processors to keep free on head node [1]
    rsh=T/F         : Whether to use rsh to run jobs on other nodes [True]

    ### ~ SEQBYSEQ OPTIONS ~ ###
    seqbyseq=T/F    : Activate seqbyseq mode - assumes basefile=X option used for output [False]
    seqin=FILE      : Input sequence file to farm out [None]
    basefile=X      : Base for output files - compiled from individual run results [None]
    outlist=LIST    : List of extensions of outputs to add to basefile for output (basefile.*) []
    pickup=X        : Header to extract from OutList file and used to populate AccNum to skip []

    ### ~ SPECIAL iRUN OPTIONS ~ ###
    runid=X         : Text identifier for iX run [None]
    resfile=FILE    : Main output file for iX run [islimfinder.csv]
    sortrun=T/F     : Whether to sort input files by size and run big -> small to avoid hang at end [True]
    loadbalance=T/F : Whether to split SortRun jobs equally between large & small to avoid memory issues [True]
    
Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../legacy/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_zen
try: import rje_seq, rje_slimlist, rje_uniprot
except: pass
try:
    import gopher_V2 as gopher
except: gopher = None
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Added additional functions to call other programs
    # 1.1 - Added UniFake.
    # 1.2 - Added generic seqbyseq option
    # 1.3 - Modified for IRIDIS3.
    # 1.4 - Added catching of IOErrors.
    # 1.5 - Added QSLiMFinder iRun
    # 1.6 - Modified iSLiMFinder job processing to try to catch errors better. (Not sure what is happening.)
    # 1.7 - Added memory checking before a run is spawned.
    # 1.8 - Added load balance option for SortRun: splits jobs equally between large and small input (& ends in middle).
    # 1.9 - Added scanning of legacy folder - moving GOPHER_V2!
    # 1.10- Modified freemem setting to run on Katana. Made rsh optional. Removed defunct IRIDIS3 option.
    # 1.10.1 - Attempted to fix SLiMFarmer batch run problem. (Should not be setting irun=batch!)
    # 1.10.2 - Trying to clean up unknown 30s pause. Might be freemem issue?
    # 1.10.3 - Fix issues with batch farming of subjobs splitting on commas.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Catch when something goes wrong with a BATCH job. Check for output somehow? Report transfer when nothing.
    # [ ] : Add memory check to iGopher, iSLiMSearch and iUniFake once check with iSLiMFinder.
    # [Y] : Add some kind of load balancing for SortRun to prevent memory problems.
    # [ ] : Replace GOPHER_V2 with gopher (V3) and test.
    # [Y] : Add forking without rsh to allow more streamlined runs on a single multi-core node.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_IRIDIS', '1.10.3', 'July 2019', '2008')
    description = 'Parallel processing on IRIDIS'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            rje.printf('\n\nHelp for {0} {1}: {2}\n'.format(info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
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
### SECTION II: IRIDIS Class                                                                                            #
#########################################################################################################################
class IRIDIS(rje.RJE_Object):     
    '''
    IRIDIS Class. Author: Rich Edwards (2008).

    Info:str
    - iIni = Ini file to pass to the called program [None]
    - iRun = Exectute a special iRun analysis on Iridis
    - PickUp = Header to extract from OutList file and used to populate AccNum to skip []
    - PyPath = path to python modules ['/home/re1u06/Serpentry/']
    
    Opt:boolean
    - RjePy = Whether program is an RJE *.py script (adds log processing) [True]
    - RSH = Whether to use rsh to run jobs on other nodes [True]
    - SeqBySeq = Activate seqbyseq mode - assumes basefile=X option used for output [False]
    - Test = Whether to produce extra output in "test" mode [False]
    
    Stat:numeric
    - IOLimit = Limit of number of IOErrors before termination [50]
    - KeepFree = Number of processors to keep free on head node [1]
    - MemFree = Min. proportion of node memory to be free before spawning job [0.0]
    - SubSleep = Sleep time (seconds) between cycles of subbing out jobs to hosts [1]

    List:list
    - Hosts = Host processors
    - OutList = List of extensions of outputs to add to basefile for output []
    - SubJobs = List of jobs to farms out to nodes

    Dict:dictionary    
    - Running = List of running subjobs (Dictionaries of {HostIndex:{'PID'/'Job'/'Log'}})
    
    Obj:RJE_Objects
    '''
    def nprocs(self): return len(self.list['Hosts'])
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['iRun','PyPath','StartFrom','iIni','Pickup']
        self.optlist = ['RjePy','SeqBySeq','Test','RSH']
        self.statlist = ['IOLimit','KeepFree','MemFree','SubSleep']
        self.listlist = ['Hosts','SubJobs','OutList']
        self.dictlist = ['Running']
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'PyPath':rje.makePath('/home/re1u06/Serpentry/')})
        self.setStat({'SubSleep':1.0,'IOLimit':50,'MemFree':0.0,'KeepFree':1})
        self.setOpt({'RjePy':True,'RSH':True})
        ### Other Attributes ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)      # General Options ###
                self._cmdReadList(cmd,'info',['iRun'])
                self._cmdRead(cmd,'info','iRun','farm')
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self._iridisCmdList()               # Allows easy incorporation when inheriting class
#########################################################################################################################
    def _iridisCmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._forkCmd(cmd)
                self._cmdReadList(cmd,'path',['PyPath'])
                self._cmdReadList(cmd,'file',['iIni'])
                self._cmdReadList(cmd,'info',['StartFrom','Pickup'])
                self._cmdRead(cmd,'info','Pickup','pickhead')
                self._cmdReadList(cmd,'int',['IOError','KeepFree'])
                self._cmdReadList(cmd,'stat',['SubSleep','MemFree'])  
                self._cmdReadList(cmd,'list',['OutList'])
                self._cmdReadList(cmd,'flist',['SubJobs'])
                self._cmdReadList(cmd,'opt',['RjePy','SeqBySeq','Test','RSH'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        if 'iIni' not in self.info or self.info['iIni'].lower() in ['none']: self.info['iIni'] = ''
        if self.stat['MemFree'] > 1: self.stat['MemFree'] /= 100.0
#########################################################################################################################
    ### <2> ### Main Run Methods                                                                                        #
#########################################################################################################################
    def iridisRun(self):  ### Main run method
        '''Main Run Method'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Set job directory to RunPath if given, else directory from which job was submitted ~ ##
            try: jobdir = rje.makePath(os.environ['PBS_O_WORKDIR'])
            except: jobdir = None
            if self.info['RunPath'] == rje.makePath(os.path.abspath(os.curdir)) and jobdir: self.info['RunPath'] = jobdir
            os.chdir(self.info['RunPath'])
            ## ~ [1b] ~ Read list of node names in file $PBS_NODEFILE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setHosts()
            if self.opt['SeqBySeq']: return self.seqBySeq()
            if self.info['iRun'] not in ['','None','batch']:
                self.printLog('#FARM','Performing %s run.' % self.getStr('iRun'))
                return iRun(self,self.info['iRun'])
            ## ~ [1c] ~ Set subjobs if not already given ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.set_subjobs()
            self.printLog('#JOBS','Running %d subjobs on %d hosts' % (len(self.list['SubJobs']),self.nprocs()))
            self.printLog('#TIME','%s second sleep between subjob cycles' % (self.getStat('SubSleep')))
            ### ~ [2] ~ Run Jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.runJobs()
        except SystemExit: raise    # Child exit
        except:
            self.log.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setHosts(self): ### Sets up self.list['Hosts']
        '''Sets up self.list['Hosts'].'''
        try: self.list['Hosts'] = self.loadFromFile(os.environ['PBS_NODEFILE'],chomplines=True)
        except: self.list['Hosts'] = []; self.printLog('#HPC','HPC PBS not detected.')
        if not self.list['Hosts']:
            if not self.getBool('RSH') and self.getInt('Forks'):
                self.list['Hosts'] = ['fork'] * self.getInt('Forks')
            else:
                self.errorLog('No Hosts identified for IRIDIS run',printerror=False)
                sys.exit()
        self.printLog('#HOSTS','%d hosts: %s.' % (self.nprocs(),rje.join(self.list['Hosts'],'; ')))
#########################################################################################################################
    def checkMemory(self,node): ### Checks free memory on node and returns True/False if OK to run
        '''Checks free memory on node and returns True/False if OK to run.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MemFree'] <= 0.0: return True
            ### ~ [2] Check details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.freeMem(node) < self.stat['MemFree']: return False
            return True
        except: self.errorLog('Problem with checkMemory(%s)' % node); return False
#########################################################################################################################
    def freeMem(self,node):
        try: memdata = rje.split(os.popen('rsh %s free' % node).readlines()[1])
        except: self.errorLog('Unable to read %s free memory' % node,printerror=False); return 0.0
        return rje.atof(memdata[3]) / rje.atof(memdata[1])
#########################################################################################################################
    def runJobs(self):  ### Runs all the jobs in self.list['SubJobs']
        '''Runs all the jobs in self.list['SubJobs'].'''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for j in range(self.stat['KeepFree'],self.nprocs()): self.nextJob(j)    # Skip first node
        pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.dict['Running']:
            PIDCHECK = open(pidcheck,'w')
            for j in rje.sortKeys(self.dict['Running']):
                if not self.dict['Running'][j]: self.dict['Running'].pop(j); continue   # No more jobs
                try:
                    pid = self.dict['Running'][j]['PID']
                    PIDCHECK.write('%s: %s\n' % (j,pid))
                    if rje.split('%s' % pid)[0] == 'WAIT': status = 1
                    else: (status,exit_stat) = os.waitpid(pid,os.WNOHANG)
                except: status = 1
                if status > 0: self.endJob(j)       # subjob on processor j has finished: can replace with processing
            PIDCHECK.close()
            time.sleep(self.stat['SubSleep'])
#########################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            node = self.list['Hosts'][host_id]
            if self.opt['SeqBySeq']: return self.nextSeqJob(host_id)
            jdict = self.dict['Running'][host_id] = {}      # Setup empty dictionary to fill, if jobs available
            if self.stat['MemFree'] > 0.0:
                freemem = self.freeMem(node)
                if freemem < self.stat['MemFree']:
                    jdict['PID'] = 'WAIT - %.1f%% %s mem' % (freemem*100.0,node)
                    return
            if self.list['SubJobs']: job = self.list['SubJobs'].pop(0)
            else: return
            if self.opt['RjePy']: jdict['Log'] = '%si_%s.log' % (self.info['RunPath'],rje.randomString(6))
            if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
            if 'Log' in jdict:
                job = '%s log=%s' % (job,jdict['Log'])
                try: open(jdict['Log'],'w')
                except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                if self.stat['MemFree'] > 0.0:
                    self.printLog('#JOB','Running job as %s: %d remain; %.1f%% mem free' % (cpid,len(self.list['SubJobs']),freemem*100.0))
                else: self.printLog('#JOB','Running job as %s: %d remain' % (cpid,len(self.list['SubJobs'])))
            else:                   # child process
                if self.getBool('RSH'): os.system(rsh)
                else: os.system(job)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('IRIDIS.nextJob error')
#########################################################################################################################
    def endJob(self,host_id):   ### Ends job on host and sets new one running
        '''Ends job on host and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SeqBySeq']: return self.endSeqJob(host_id)
            jdict = self.dict['Running'].pop(host_id)
            if 'Log' in jdict:
                open(self.log.info['LogFile'],'a').writelines(open(jdict['Log'],'r').readlines()[5:-1])
                self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                self.printLog('#~~#','#~~#',timeout=False)
                os.unlink(jdict['Log'])
            elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': pass
            else: self.printLog('#END','Job on processor %d ended.' % host_id)
        except IOError:
            if self.stat['IOError'] == 1: self.errorLog('IRIDIS.endJob IOError limit reached'); raise
            else: self.stat['IOError'] -= 1; self.errorLog('IRIDIS.endJob')
        except: self.errorLog('IRIDIS.endJob error')
        self.nextJob(host_id)   # Carry on regardless
#########################################################################################################################
    ### <3> ### SeqBySeq Run Methods                                                                                    #
#########################################################################################################################
    def seqBySeq(self):     ### Runs in SeqBySeq Mode
        '''
        In SeqBySeq mode, the program assumes that seqin=FILE and basefile=X are given and irun states the program to be run.
        Seqin will then be worked through in turn and each sequence farmed out to the irun program. Outputs given by OutList
        are then compiled, as is the Log, into the correct basefile=X given. In the case of *.csv and *.tdt files, the header
        row is copied for the first file and then excluded for all subsequent files. For all other files extensions, the
        whole output is copied.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['iRun'][-3:] == '.py': self.info['iRun'] = self.info['iRun'][:-3]
            if self.info['StartFrom'].lower() == 'none': self.info['StartFrom'] = ''
            self.list['Seq'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','accnr=F','seqnr=F']).seq[0:]
            while self.info['StartFrom'] and self.list['Seq']:
                if self.list['Seq'][0].shortName() != self.info['StartFrom']: self.list['Seq'] = self.list['Seq'][1:]
                else: self.info['StartFrom'] = ''
            self.printLog('#SEQ','%s query sequences to farm out' % rje.integerString(len(self.list['Seq'])))
            self.list['Pickup'] = self.pickupList()
            ### ~ [2] ~ Run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.runJobs()
            return True        
        except SystemExit: raise    # Child
        except: self.errorLog('IRIDIS.seqBySeq error')
        return False            
#########################################################################################################################
    def nextSeqJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.list['Seq']: seq = self.list['Seq'].pop(0)  # UniFake this sequence
            else: return                                        # Out of sequences: stop
            if seq.info['AccNum'] in self.list['Pickup']: return self.nextSeqJob(host_id)   # Skip this sequence
            jran = 'i_%s' % rje.randomString(6)
            jdict['Log'] = '%s%s.log' % (self.info['RunPath'],jran)
            jdict['Qry'] = '%s%s.qry' % (self.info['RunPath'],jran)
            for out in self.list['OutList']: jdict[out] = '%s%s.%s' % (self.info['RunPath'],jran,out)
            open(jdict['Qry'],'w').write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s on `hostname` as %s ; ' % (seq.shortName(),jran)
            if self.getBool('RSH'): job = '%s python %s%s.py' % (initial_cmds,self.info['PyPath'],self.info['iRun'])
            else: job = 'python %s%s.py' % (self.info['PyPath'],self.info['iRun'])
            if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
            job = '%s seqin=%s i=-1 v=-1 basefile=%s' % (job,jdict['Qry'],jran)
            if self.getBool('RSH'): job = '%s ; echo Finishing on `hostname`' % job
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: cpid = os.fork()        # Fork child process
            except OSError:
                self.errorLog('IRIDIS.nextJob error. Will sleep for 10 mins then continue.')
                time.sleep(600)
                self.printLog('#YAWN','Re-awakening IRIDIS nextJob. Fingers crossed.')
                jdict = {'Error':'IRIDIS.nextJob OSError.'}     # Make sure this node is retried.
                return
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#SEQ','Running %s as %s [%d::%s]: %d remain' % (seq.shortName(),cpid,host_id,self.list['Hosts'][host_id],len(self.list['Seq'])))
            else:                   # child process
                if self.getBool('RSH'): os.system(rsh)
                else: os.system(job)
                os._exit(0) 
        except SystemExit: raise    # Child
        except: self.errorLog('IRIDIS.nextSeqJob error')
#########################################################################################################################
    def endSeqJob(self,host_id):   ### Ends job on host and sets new one running
        '''Ends job on host and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                for out in self.list['OutList']:
                    resfile = '%s.%s' % (self.info['Basefile'],out)
                    if os.path.exists(jdict[out]):
                        if os.path.exists(resfile) and rje.split(resfile,'.')[-1] in ['tdt','csv']: open(resfile,'a').writelines(open(jdict[out],'r').readlines()[1:])
                        else: open(resfile,'a').writelines(open(jdict[out],'r').readlines())
                        os.unlink(jdict[out])
                    else: self.errorLog('Output %s file "%s" missing' % (out,jdict[out]),printerror=False)
                if os.path.exists(jdict['Log']):
                    loglines = open(jdict['Log'],'r').readlines()[5:-1]
                    if loglines: open(self.log.info['LogFile'],'a').writelines(loglines)
                    self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                    self.printLog('#~~#','#~~#',timeout=False)
                    os.unlink(jdict['Log'])
                elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': pass
                else: self.printLog('#END','Job on processor %d ended: already complete' % host_id)
                os.system('rm %s*' % jdict['Qry'])
        except: self.errorLog('IRIDIS.endJob error')
        self.nextJob(host_id)   # Carry on regardless        
#########################################################################################################################
    def pickupList(self):   ### Generates Pickup List from file(s)
        '''Generates Pickup List from file(s).'''
        if self.info['Pickup'].lower() in ['','none']: return []
        picklist = []
        for out in self.list['OutList']:
            resfile = '%s.%s' % (self.info['Basefile'],out)
            try: pickdat = rje.dataDict(self,resfile,[self.info['Pickup']])
            except: pickdat = {}
            picklist = picklist + list(pickdat.keys())
        return picklist
#########################################################################################################################
    ### <4> ### Testing Methods                                                                                         #
#########################################################################################################################
    def set_subjobs(self):  ### Sets temporary test list of subjobs to run
        '''Sets up list of subjob commands to be executed in array subjobs.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['SubJobs']: return     # Already given jobs, so just carry on with these
            ### ~ [2] ~ Make a test set of Zen jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jobnum = 2 * self.nprocs() + 3  # Default number of jobs
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo `pwd` on `hostname` ; '
            for nj in range(jobnum):
                self.list['SubJobs'].append('%s python %srje_zen.py wisdoms=%d' % (initial_cmds,self.info['PyPath'],random.randint(1,20)))
        except: self.errorLog('IRIDIS.set_subjobs error')
#########################################################################################################################
### End of SECTION II: IRIDIS Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def iRun(irun,itype):    ### Executes a special IRIDIS running using IRIDIS object irun
    '''
    Executes a special IRIDIS running using IRIDIS object irun. Generally will use ini file for options.
    >> irun:IRIDIS Object controlling parallel runs
    >> type:str = identifier for special IRIDIS run
    '''
    try:### ~ [1] ~ Select iRun method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if itype.lower() in ['gopher','igopher']: i = iGopher(irun.log,irun.cmd_list)
        elif itype.lower() == 'slimfinder': i = iSLiMFinder(irun.log,irun.cmd_list)
        elif itype.lower() == 'qslimfinder': i = iSLiMFinder(irun.log,irun.cmd_list+['qslimfinder=T'])
        elif itype.lower() == 'slimsearch': i = iSLiMSearch(irun.log,irun.cmd_list)
        elif itype.lower() == 'unifake': i = iUniFake(irun.log,irun.cmd_list)
        else: raise ValueError('%s is not a valid rje_iridis.IRun() type!' % itype)
        i.info['RunPath'] == rje.makePath(os.path.abspath(os.curdir))
        i.list['Hosts'] = irun.list['Hosts']
        return i.run()
    except SystemExit: raise
    except: irun.errorLog('Problem with iRun(%s)' % itype)
    return False
#########################################################################################################################
class iSLiMFinder(IRIDIS):  ### Farms a SLiMFinder run out over hosts
    '''Farms a SLiMFinder run out over hosts.'''
    #####################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self.opt['SLiMCore'] = False
        self.opt['SLiMSearch'] = False
        self.opt['QSLiMFinder'] = False
        self.opt['MakeDAT'] = False
        self.opt['Pickup'] = True
        self.info['RunID'] = 'None'
        self.info['ResFile'] = 'islimfinder.csv'	# Make sure resfile=X is set on commandline
        self.opt['SortRun'] = True
        self.opt['LoadBalance'] = True
        self.list['JRan'] = []      # List of random job numbers for checking and cleanup
        self.list['Batch'] = []     # List of files for batch running.
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)      # General Options ###
                self._cmdReadList(cmd,'list',['Batch'])
                self._cmdReadList(cmd,'file',['ResFile'])
                self._cmdReadList(cmd,'info',['RunID'])
                self._cmdReadList(cmd,'opt',['SLiMCore','QSLiMFinder','SLiMSearch','MakeDAT','Pickup','SortRun','LoadBalance'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self._iridisCmdList()               # Allows easy incorporation when inheriting class
    #####################################################################################################################
    def run(self):
        try:### ~ [0] ~ Pickup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pickup = []
            try:
                if self.opt['Pickup'] and not self.opt['SLiMCore'] and rje.checkForFile(self.info['ResFile']):
                    try: pickup = rje.dataDict(self,self.info['ResFile'],['RunID'],['Dataset'],lists=True)
                    except: pickdat = {}
                    if self.info['RunID'] in pickup:
                        pickup = pickup[self.info['RunID']]['Dataset']
                    else:
                        self.printLog('#PICKUP','RunID %s not found: %s' % (self.info['RunID'],rje.sortKeys(pickup)))
                        pickup = []
                    self.printLog('#PICKUP','%d %s datasets to skip' % (len(pickup),self.info['RunID']))
                elif self.opt['Pickup'] and self.opt['SLiMCore'] and rje.checkForFile(self.info['Pickup']):
                    pickup = self.loadFromFile(self.info['Pickup'],chomplines=True)
                    self.printLog('#PICKUP','%d %s datasets to skip' % (len(pickup),self.info['Pickup']))
            except: self.errorLog('Pickup Failed'); pickup = []
            ### ~ [1] ~ Batch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
            self.printLog('\r#BATCH','%s batch files to run.' % rje.iLen(batchfiles))
            sizedict = {}
            newbatch = []       # New list of batch files in size order, starting with the biggest
            for file in batchfiles:
                if os.path.split(rje.baseFile(file))[1] in pickup:
                    self.printLog('#PICKUP','Skipping %s batch file %s' % (self.info['RunID'],file))
                    continue
                if self.opt['SortRun']:
                    fsize = os.path.getsize(file)
                    if fsize not in sizedict: sizedict[fsize] = []
                    sizedict[fsize].append(file)
                else: newbatch.append(file)
            if self.opt['SortRun']:
                sizes = rje.sortKeys(sizedict,revsort=True)
                for fsize in sizes: newbatch += sizedict.pop(fsize)
            self.list['NewBatch'] = newbatch
            self.stat['BatchNum'] = len(self.list['NewBatch'])
            self.printLog('#JOBS','%s Batch jobs to run.' % rje.integerString(self.stat['BatchNum']))
            ### ~ [2] ~ Run Jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.runJobs()
            self.cleanUp()
            return True        
        except SystemExit: raise    # Child
        except: self.errorLog('Problem with iSLiMFinder.run()'); self.cleanUp()
        return False
    #####################################################################################################################
    def runJobs(self):  ### Runs all the jobs in self.list['SubJobs']
        '''Runs all the jobs in self.list['SubJobs'].'''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        checkn = self.nprocs() - self.stat['KeepFree']
        funny5min = False; big_problem = False
        nodes = rje.sortUnique(self.list['Hosts'])
        self.stat['PPN'] = ppn = self.nprocs() / len(nodes)
        jlist = range(self.stat['KeepFree'],self.nprocs())
        self.printLog('#PPN','%d ppn on %d host nodes (%d jobs + %d head)' % (ppn,len(nodes),len(jlist),self.stat['KeepFree']))
        #for j in range(1,self.nprocs()): self.nextJob(j)    # Skip first node
        for p in range(ppn):
            for n in range(len(nodes)):
                j = n * ppn + p
                if j in jlist: self.nextJob(j)    # Skip first node
            time.sleep(self.stat['SubSleep'])
        pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.dict['Running'] or self.list['NewBatch']:
            pidchecklist = []
            for j in rje.sortKeys(self.dict['Running']):
                if not self.dict['Running'][j]: self.dict['Running'].pop(j); continue   # No more jobs
                pid = '????'
                try:
                    pid = self.dict['Running'][j]['PID']
                    if rje.split('%s' % pid)[0] == 'WAIT':
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
                if status > 0: self.endJob(j)       # subjob on processor j has finished: can replace with processing
            PIDCHECK = open(pidcheck,'w')
            PIDCHECK.write(rje.join(pidchecklist+[''],'\n'))
            PIDCHECK.close()
            if self.list['NewBatch'] and len(self.dict['Running']) < checkn:    ### Why are jobs missing?
                if big_problem:
                    self.printLog('#ERR','Still only %d of %d jobs running but %d remain. Self-destructing.' % (len(self.dict['Running']),checkn))
                    self.cleanUp()
                    raise ValueError
                elif funny5min:
                    self.printLog('#JOBS','Still only %d of %d jobs running but %d remain. Will try starting some more.' % (len(self.dict['Running']),checkn,len(self.list['NewBatch'])))
                    for j in range(1,self.nprocs()):
                        if j not in self.dict['Running']: self.nextJob(j)    # Skip first node
                    big_problem = True
                else:
                    self.printLog('#ERR','Only %d of %d jobs running but %d remain. Funny 5 minutes?' % (len(self.dict['Running']),checkn,len(self.list['NewBatch'])))
                    funny5min = True
                    time.sleep(300)
                    continue
            else: funny5min = False; big_problem = False
            time.sleep(self.stat['SubSleep'])
    #####################################################################################################################
    def cleanUp(self):  ### Tries to clean up random run files that are still hanging about
        '''Tries to clean up random run files that are still hanging about.'''
        try:### ~ [1] Remove jobs without files from self.list['JRan'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            precleanx = len(self.list['JRan'])
            for jran in self.list['JRan'][0:]:
                jlog = '%s%s.log' % (self.info['RunPath'],jran)
                jcsv = '%s%s.csv' % (self.info['RunPath'],jran)
                if not os.path.exists(jlog) and not os.path.exists(jcsv): self.list['JRan'].remove(jran)
            self.printLog('#CLEAN','%s of %s job runs have files left' % (rje.iLen(self.list['JRan']),rje.iStr(precleanx)))
            if not self.list['JRan']: return
            ### ~ [2] Try cleaning up the files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tranx = 0; pretranx = len(self.list['JRan'])
            for jran in self.list['JRan'][0:]:
                jlog = '%s%s.log' % (self.info['RunPath'],jran)
                jcsv = '%s%s.csv' % (self.info['RunPath'],jran)
                jocc = '%s%s.occ.csv' % (self.info['RunPath'],jran)
                socc = '%s.occ.csv' % self.info['ResFile'][:-4]
                if os.path.exists(jlog):
                    jloglines = open(jlog,'r').readlines()
                    if jloglines:
                        if rje.split(jloglines[-1])[0] == '#LOG' and 'End:' in rje.split(jloglines[-1]):
                            self.printLog('#END','Random job %s looks finished' % jran)
                            if os.path.exists(jcsv):
                                if os.path.exists(self.info['ResFile']): open(self.info['ResFile'],'a').writelines(open(jcsv,'r').readlines()[1:])
                                else: open(self.info['ResFile'],'w').writelines(open(jcsv,'r').readlines())     # Headers too
                            if os.path.exists(jocc):
                                if os.path.exists(socc): open(socc,'a').writelines(open(jocc,'r').readlines()[1:])
                                else: open(socc,'w').writelines(open(jocc,'r').readlines())     # Headers too
                            open(self.log.info['LogFile'],'a').writelines(open(jlog,'r').readlines()[5:-1])
                            self.printLog('#END','%s results & log content transferred' % jran); tranx += 1
                        else: self.printLog('#ERR','%s only reached: %s' % (jran,jloglines[-1]))
                    os.unlink(jlog)
                    self.printLog('#DEL','%s deleted in cleanup' % jlog)
                if os.path.exists(jcsv):
                    os.unlink(jcsv)
                    self.printLog('#DEL','%s deleted in cleanup' % jcsv)
            self.printLog('#CLEAN','%s of %s job runs had results & logs transferred' % (rje.iStr(tranx),rje.iLen(self.list['JRan'])))
        except: self.errorLog('Cleanup error. %s runs not checked.' % rje.iLen(self.list['JRan']))
    #####################################################################################################################
    def endJob(self,host_id):   ### Replaced endJob to set a new SLiMFinder going
        '''Replaced endJob to set a new Gopher going.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            waiting = False
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                if 'Log' in jdict:
                    if os.path.exists(jdict['Res']):
                        if os.path.exists(self.info['ResFile']): open(self.info['ResFile'],'a').writelines(open(jdict['Res'],'r').readlines()[1:])
                        else: open(self.info['ResFile'],'w').writelines(open(jdict['Res'],'r').readlines())     # Headers too
                        os.unlink(jdict['Res'])
                    jfile = '%s.occ.csv' % jdict['Res'][:-4]
                    sfile = '%s.occ.csv' % self.info['ResFile'][:-4]
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
                        time.sleep(self.stat['SubSleep'])
                    self.printLog('#END',etxt)
                    self.printLog('#~~#','#~~#',timeout=False)
                #x#if 'DAT' in jdict: os.system('rm %s*' % jdict['DAT'])
                elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': waiting = True
                else: self.printLog('#END','Job on processor %d ended.' % host_id)
        except IOError:
            if self.stat['IOError'] == 1: self.errorLog('iSLiMFinder.endJob IOError limit reached'); raise
            else: self.stat['IOError'] -= 1; self.errorLog('iSLiMFinder.endJob')
        except: self.errorLog('iSLiMFinder.endJob error')
        self.nextJob(host_id)   # Carry on regardless
        if waiting:
            try:
                if rje.split('%s' % self.dict['Running'][host_id]['PID'])[0] == 'WAIT':
                    if self.opt['Test']: self.printLog('#WAIT','Resumed after waiting for memory. %d sec sleep!' % self.stat['SubSleep'])
                    time.sleep(self.stat['SubSleep'])
            except:
                if self.opt['Test']: self.printLog('#WAIT','Resumed after waiting for memory but no new PID!')
    #####################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            node = self.list['Hosts'][host_id]
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.stat['MemFree'] > 0.0: freemem = self.freeMem(node)
            if self.stat['MemFree'] > 0.0 and freemem < self.stat['MemFree']:
                jdict['PID'] = 'WAIT - %.1f%% %s mem (host %d)' % (freemem*100.0,node,host_id)
                if self.opt['Test']: self.printLog('#JOB',jdict['PID'])
                return
            if self.list['NewBatch']:
                if host_id < self.stat['PPN']: next = self.list['NewBatch'].pop(-1)     # Keep head node running small jobs
                elif self.opt['SortRun'] and self.opt['LoadBalance'] and not len(self.list['NewBatch']) % 2: next = self.list['NewBatch'].pop(-1)
                else: next = self.list['NewBatch'].pop(0)   # SLiMFinder this dataset
            else: return                                                    # Out of datasets: stop
            jran = 'i_%s_%s' % (host_id,rje.randomString(6))
            while jran in self.list['JRan'] or os.path.exists('%s%s.log' % (self.info['RunPath'],jran)):
                jran = 'i_%s_%s' % (host_id,rje.randomString(6))
            self.list['JRan'].append(jran)
            jdict['JRan'] = jran; jdict['Next'] = next
            jdict['Log'] = '%s%s.log' % (self.info['RunPath'],jran)
            jdict['Res'] = '%s%s.csv' % (self.info['RunPath'],jran)
            jdict['Start'] = time.time()
            try: open(jdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            #x#initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as %s on `hostname` ; setenv IUPred_PATH /home/re1u06/Bioware/iupred/ ;' % (next,jran)
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as %s on `hostname` ;' % (next,jran)  # bash: setenv: command not found
            if next[-3:] == 'acc':
                jdict['DAT'] = '%sdat' % next[:-3]
                while os.path.exists(jdict['DAT']):   # No need to run - skip onto the next one
                    if self.list['NewBatch']: next = self.list['NewBatch'].pop(0)
                    else: return
                    jdict['DAT'] = '%sdat' % next[:-3]
                    initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as %s on `hostname` ; ' % (next,jran)
                job = 'python %srje_uniprot.py extract=%s datout=%s i=-1 v=-1 log=%s' % (self.info['PyPath'],next,jdict['DAT'],jdict['Log'])
                #initial_cmds = '%s python %srje_uniprot.py extract=%s datout=%s i=-1 v=-1 log=%s ; ' % (initial_cmds,self.info['PyPath'],next,jdict['DAT'],jdict['Log'])
                #next = jdict['DAT']
            else:
                if self.opt['SLiMCore']: job = 'python %srje_slimcore.py' % (self.info['PyPath'])
                elif self.opt['SLiMSearch']: job = 'python %sslimsearch.py' % (self.info['PyPath'])
                elif self.opt['QSLiMFinder']: job = 'python %sqslimfinder.py' % (self.info['PyPath'])
                else: job = 'python %sslimfinder.py' % (self.info['PyPath'])
                if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
                job = '%s pickup=F batch= seqin=%s i=-1 v=-1 runid=%s' % (job,next,self.getStr('RunID'))
                job = '%s resfile=%s log=%s' % (job,jdict['Res'],jdict['Log'])
            rsh = "rsh %s '%s %s'" % (self.list['Hosts'][host_id],initial_cmds,job)
            if self.getBool('RSH'): self.printLog('#RSH',rsh)
            else: self.printLog('#FORK',job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                remaining = len(self.list['NewBatch']); runnum = self.stat['BatchNum'] - remaining
                if self.stat['MemFree'] > 0.0: freemem = '; %.1f%% mem free.' % (freemem*100.0)
                else: freemem = '.'
                self.printLog('#BATCH','Running %s as %s (pid:%s) [%d::%s]: %d run; %d remain%s' % (os.path.basename(next),jran,cpid,host_id,self.list['Hosts'][host_id],runnum,remaining,freemem))
            else:                   # child process
                if self.getBool('RSH'): os.system(rsh)
                else: os.system(job)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('iSLiMFinder.nextJob error')
#########################################################################################################################
class iSLiMSearch(IRIDIS):  ### Farms a SLiMSearch run out over hosts
    '''Farms a SLiMSearch run out over hosts.'''
    #####################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)      # General Options ###
                self._cmdReadList(cmd,'file',['ResFile'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self._iridisCmdList()               # Allows easy incorporation when inheriting class
    #####################################################################################################################
    def run(self):
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list+['reverse=F','wildscram=F'])
            self.obj['SlimList'].loadMotifs()
            self.list['Slims'] = self.obj['SlimList'].slims()
            self.runJobs()
            return True        
        except SystemExit: raise    # Child
        except: self.errorLog('Problem with iSLiMSearch.run()')
        return False
    #####################################################################################################################
    def endJob(self,host_id):   ### Replaced endJob to set a new Gopher going
        '''Replaced endJob to set a new Gopher going.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                if 'Log' in jdict:
                    if os.path.exists(jdict['Res']):
                        if os.path.exists(self.info['ResFile']): open(self.info['ResFile'],'a').writelines(open(jdict['Res'],'r').readlines()[1:])
                        else: open(self.info['ResFile'],'w').writelines(open(jdict['Res'],'r').readlines())     # Headers too
                        os.unlink(jdict['Res'])
                    jfile = '%s.summary.csv' % jdict['Res'][:-4]
                    sfile = '%s.summary.csv' % self.info['ResFile'][:-4]
                    if os.path.exists(jfile):
                        if os.path.exists(sfile): open(sfile,'a').writelines(open(jfile,'r').readlines()[1:])
                        else: open(sfile,'w').writelines(open(jfile,'r').readlines())     # Headers too
                        os.unlink(jfile)                        
                    open(self.log.info['LogFile'],'a').writelines(open(jdict['Log'],'r').readlines()[5:-1])
                    self.printLog('#END','Job on processor %d ended: results & log content transferred' % host_id)
                    self.printLog('#~~#','#~~#',timeout=False)
                    os.unlink(jdict['Log'])
                elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': pass
                else: self.printLog('#END','Job on processor %d ended.' % host_id)
                if 'Motif' in jdict: os.unlink(jdict['Motif'])
        except IOError:
            if self.stat['IOError'] == 1: self.errorLog('iSLiMSearch.endJob IOError limit reached'); raise
            else: self.stat['IOError'] -= 1; self.errorLog('iSLiMSearch.endJob')
        except: self.errorLog('iSLiMSearch.endJob error')
        self.nextJob(host_id)   # Carry on regardless        
    #####################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.list['Slims']:
                next = self.list['Slims'].pop(0)    # SLiMSearch this slim
                name = next.info['Name']
                if next.info['Description'] != name: name = '%s %s' % (name,next.info['Description'])
                slimlist = rje_slimlist.SLiMList(self.log,self.cmd_list+['reverse=F','wildscram=F'])
                slimlist._addMotif(name,next.pattern(),reverse=False,check=False,logrem=False,wildscram=False)
                slimlist._addMotif(name,next.pattern(),reverse=True,check=False,logrem=False,wildscram=False)
                slimlist._addMotif(name,next.pattern(),reverse=False,check=False,logrem=False,wildscram=True)
            else: return                                                    # Out of datasets: stop
            jran = 'i_%s' % rje.randomString(6)
            jdict['Log'] = '%s%s.log' % (self.info['RunPath'],jran)
            jdict['Motif'] = '%s%s.motif' % (self.info['RunPath'],jran)
            jdict['Res'] = '%s%s.csv' % (self.info['RunPath'],jran)
            try: open(jdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            slimlist.motifOut(filename=jdict['Motif'])
            #x#initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as %s on `hostname` ; setenv IUPred_PATH /home/re1u06/Bioware/iupred/ ;' % (next,jran)
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s as %s on `hostname` ;' % (next,jran)  # bash: setenv: command not found
            job = '%s python %sslimsearch.py motifs=%s i=-1 v=-1 resfile=%s log=%s' % (initial_cmds,self.info['PyPath'],jdict['Motif'],jdict['Res'],jdict['Log'])
            if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
            job = '%s reverse=F wildscram=F' % job
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            self.log.printLog('#RSH',rsh)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#BATCH','Running %s as %s [%d::%s]: %d remain' % (next.info['Name'],cpid,host_id,self.list['Hosts'][host_id],len(self.list['Slims'])))
            else:                   # child process
                os.system(rsh)
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('iSLiMSearch.nextJob error')
#########################################################################################################################
class iGopher(IRIDIS):  ### Farms a GOPHER run out over hosts
    '''Farms a GOPHER run out over hosts.'''
    #####################################################################################################################
    def run(self):
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqcmd = self.cmd_list+['autoload=T','accnr=F','seqnr=F','autofilter=F']
            try: gopher.Gopher(log=self.log,cmd_list=self.cmd_list[0:]).setupBlast()
            except: self.errorLog('Problem with Gopher.setupBlast()')
            if self.info['StartFrom'].lower() == 'none': self.info['StartFrom'] = ''
            self.list['Seq'] = rje_seq.SeqList(self.log,seqcmd).seq[0:]
            self.printLog('#SEQN','%s sequences for iGopher run' % rje.integerString(len(self.list['Seq'])))
            while self.info['StartFrom'] and self.list['Seq']:
                if self.list['Seq'][0].shortName() != self.info['StartFrom']: self.list['Seq'] = self.list['Seq'][1:]
                else: self.info['StartFrom'] = ''
            #x#open('igopher.ini','w').write(rje.join(self.cmd_list+['i=-1','v=-1'],'\n'))
            self.runJobs()
            return True        
        except SystemExit: raise    # Child
        except: self.errorLog('Problem with main run()')
        return False
    #####################################################################################################################
    def endJob(self,host_id):   ### Replaced endJob to set a new Gopher going
        '''Replaced endJob to set a new Gopher going.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                if 'Log' in jdict:
                    loglines = open(jdict['Log'],'r').readlines()[5:-1]
                    if len(loglines) > 2:
                        open(self.log.info['LogFile'],'a').writelines(loglines)
                        self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                        self.printLog('#~~#','#~~#',timeout=False)
                    else: self.printLog('#END','Job on processor %d ended: GOPHER already complete' % host_id)
                    os.unlink(jdict['Log'])
                    os.system('rm %s*' % jdict['Qry'])
                elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': pass
                else: self.printLog('#END','Job on processor %d ended.' % host_id)
        except IOError:
            if self.stat['IOError'] == 1: self.errorLog('iGOPHER.endJob IOError limit reached'); raise
            else: self.stat['IOError'] -= 1; self.errorLog('iGOPHER.endJob')
        except: self.errorLog('iGOPHER.endJob error')
        self.nextJob(host_id)   # Carry on regardless        
    #####################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.list['Seq']: seq = self.list['Seq'].pop(0)  # GOPHER this sequence
            else: return                                        # Out of sequences: stop
            while os.path.exists('ALN/%s.orthaln.fas' % seq.info['AccNum']):   # No need to run - skip onto the next one
                self.printLog('#SKIP','Skipping %s - orthaln found' % seq.shortName())
                if self.list['Seq']: seq = self.list['Seq'].pop(0)  # GOPHER this sequence
                else: return                                        # Out of sequences: stop
            jran = 'i_%s' % rje.randomString(6)
            jdict['Log'] = '%s%s.log' % (self.info['RunPath'],jran)
            jdict['Qry'] = '%s%s.qry' % (self.info['RunPath'],jran)
            try: open(jdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            open(jdict['Qry'],'w').write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s on `hostname` ; ' % (seq.shortName())
            job = '%s python %sgopher_V2.py seqin=%s i=-1 v=-1 log=%s' % (initial_cmds,self.info['PyPath'],jdict['Qry'],jdict['Log'])
            if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
            job = '%s ; echo Finishing on `hostname`' % job
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: cpid = os.fork()        # Fork child process
            except OSError:
                self.errorLog('IRIDIS.nextJob error. Will sleep for 10 mins then continue.')
                time.sleep(600)
                self.printLog('#YAWN','Re-awakening IRIDIS nextJob. Fingers crossed.')
                jdict = {'Error':'IRIDIS.nextJob OSError.'}     # Make sure this node is retried.
                return
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                #x#print '::%d:: %s' % (host_id,rsh)
                #x#self.printLog('#HOSTS','%s' % self.list['Hosts'])
                self.printLog('#SEQ','Running %s as %s [%d::%s]: %d remain' % (seq.shortName(),cpid,host_id,self.list['Hosts'][host_id],len(self.list['Seq'])))
                #x#self.printLog('#RSH',rsh)
            else:                   # child process
                #x#print '::%d:: %s' % (host_id,rsh)
                os.system(rsh)
                #x#sys.exit()
                os._exit(0) 
        except SystemExit: raise    # Child
        except: self.errorLog('IRIDIS.nextJob error')
#########################################################################################################################
class iUniFake(IRIDIS):  ### Farms a UniFake run out over hosts
    '''Farms a UniFake run out over hosts.'''
    #####################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self.info['ResFile'] = 'unifake.dat'
        self.stat['Counter'] = 0
        self.list['TMCrap'] = []
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)      # General Options ###
                self._cmdReadList(cmd,'file',['ResFile'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self._iridisCmdList()               # Allows easy incorporation when inheriting class
    #####################################################################################################################
    def run(self):
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['StartFrom'].lower() == 'none': self.info['StartFrom'] = ''
            self.list['Seq'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','accnr=F','seqnr=F']).seq[0:]
            while self.info['StartFrom'] and self.list['Seq']:
                if self.list['Seq'][0].shortName() != self.info['StartFrom']: self.list['Seq'] = self.list['Seq'][1:]
                else: self.info['StartFrom'] = ''
            self.runJobs()
            ### ~ [3] ~ Cleanup crap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for tmp in glob.glob('TMHMM*'):
                if os.path.isdir(tmp): os.rmdir(tmp)
            for tmp in glob.glob('i_*dat'): os.unlink(tmp)
            ### ~ [4] ~ Make index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.stat['Interactive'] = -1
            self.info['UniPath'] = rje.makePath(os.path.split(self.info['ResFile'])[0])
            rje_uniprot.processUniProt(self,makeindex=True,makespec=False,makefas=False)
            return True        
        except SystemExit: raise    # Child
        except: self.errorLog('Problem with main run()')
        return False
    #####################################################################################################################
    def endJob(self,host_id):   ### Replaced endJob to set a new Gopher going
        '''Replaced endJob to set a new UniFake going.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if host_id in self.dict['Running']:
                jdict = self.dict['Running'].pop(host_id)
                if 'Log' in jdict:
                    if os.path.exists(jdict['Res']):
                        open(self.info['ResFile'],'a').writelines(open(jdict['Res'],'r').readlines())     
                        os.unlink(jdict['Res'])
                    if os.path.exists(jdict['PRes']):
                        presfile = self.info['ResFile'][:-3] + 'pfam.tdt'
                        if os.path.exists(presfile): open(presfile,'a').writelines(open(jdict['PRes'],'r').readlines()[1:])
                        else: open(presfile,'w').writelines(open(jdict['PRes'],'r').readlines())     # Headers too
                        os.unlink(jdict['PRes'])
                    loglines = open(jdict['Log'],'r').readlines()[5:-1]
                    if len(loglines) > 2:
                        open(self.log.info['LogFile'],'a').writelines(loglines)
                        self.printLog('#END','Job on processor %d ended: log content transferred' % host_id)
                        self.printLog('#~~#','#~~#',timeout=False)
                    else: self.printLog('#END','Job on processor %d ended: GOPHER already complete' % host_id)
                    os.unlink(jdict['Log'])
                    os.system('rm %s*' % jdict['Qry'])
                elif 'PID' in jdict and rje.split('%s' % jdict['PID'])[0] == 'WAIT': pass
                else: self.printLog('#END','Job on processor %d ended.' % host_id)
            ### ~ [2] ~ Cleanup TMHMM Crap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.stat['Counter'] += 1
                if self.stat['Counter'] >= 500:
                    for tmp in self.list['TMCrap']:
                        if os.path.isdir(tmp): os.rmdir(tmp)
                    self.list['TMCrap'] = glob.glob('TMHMM*')
        except IOError:
            if self.stat['IOError'] == 1: self.errorLog('iUniFake.endJob IOError limit reached'); raise
            else: self.stat['IOError'] -= 1; self.errorLog('iUniFake.endJob')
        except: self.errorLog('iUniFake.endJob error')
        self.nextJob(host_id)   # Carry on regardless        
    #####################################################################################################################
    def nextJob(self,host_id):  ### Sets an new job running on host with given index
        '''Sets an new job running on host with given index.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            jdict = self.dict['Running'][host_id] = {}          # Setup empty dictionary to fill, if jobs available
            if self.list['Seq']: seq = self.list['Seq'].pop(0)  # UniFake this sequence
            else: return                                        # Out of sequences: stop
            jran = 'i_%s' % rje.randomString(6)
            jdict['Log'] = '%s%s.log' % (self.info['RunPath'],jran)
            jdict['Qry'] = '%s%s.qry' % (self.info['RunPath'],jran)
            jdict['Res'] = '%s%s.dat' % (self.info['RunPath'],jran)
            jdict['PRes'] = '%s%s.pfam.tdt' % (self.info['RunPath'],jran)
            try: open(jdict['Log'],'w')
            except: self.errorLog('Log problem. Aborting %s job.' % host_id); return self.endJob(host_id)
            open(jdict['Qry'],'w').write('>%s\n%s\n' % (seq.info['Name'],seq.info['Sequence']))
            initial_cmds = 'cd ' + self.info['RunPath'] + ' ; echo %s on `hostname` as %s ; ' % (seq.shortName(),jran)
            job = '%s python %sunifake.py' % (initial_cmds,self.info['PyPath'])
            if self.info['iIni']: job = '%s ini=%s' % (job,self.info['iIni'])
            job = '%s datout=%s seqin=%s i=-1 v=-1 log=%s cleanup=F makeindex=F' % (job,jdict['Res'],jdict['Qry'],jdict['Log'])
            job = '%s ; echo Finishing on `hostname`' % job
            rsh = "rsh %s '%s'" % (self.list['Hosts'][host_id],job)
            ### ~ [2] ~ Add Job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: cpid = os.fork()        # Fork child process
            except OSError:
                self.errorLog('IRIDIS.nextJob error. Will sleep for 10 mins then continue.')
                time.sleep(600)
                self.printLog('#YAWN','Re-awakening IRIDIS nextJob. Fingers crossed.')
                jdict = {'Error':'IRIDIS.nextJob OSError.'}     # Make sure this node is retried.
                return
            if cpid:                # parent process records pid of child rsh process
                jdict['PID'] = cpid
                self.printLog('#SEQ','Running %s as %s [%d::%s]: %d remain' % (seq.shortName(),cpid,host_id,self.list['Hosts'][host_id],len(self.list['Seq'])))
            else:                   # child process
                os.system(rsh)
                os._exit(0) 
        except SystemExit: raise    # Child
        except: self.errorLog('IRIDIS.nextJob error')
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
    try: IRIDIS(mainlog,cmd_list).iridisRun()
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: os._exit(0)  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: rje.printf('Cataclysmic run error: {0}'.format(sys.exc_info()[0]))
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
