#!/usr/bin/python

# See below for name and description
# Copyright (C) 2013 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_forker
Description:  Generic RJE Forking Module
Version:      0.1.0
Last Edit:    10/12/20
Copyright (C) 2013  Richard J. Edwards - See source code for GNU License Notice

Function:
    The primary function of this module is to provide a generic Forker class that can be used by other objects. This is
    loosely based on the IRIDIS Class of rje_iridis but with the difference that it does not fork out processes to other
    host nodes.

Commandline:

    ### ~ FORK CONTROL OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    forkdir=PATH    : Alternative directory for forking output. [`./`]
    tofork=LIST     : List of system commands to fork out (standalone Forker only) []
    iolimit=X       : Limit of number of IOErrors before termination [50]
    memfree=X       : Min. proportion of node memory to be free before spawning new fork [0.1]
    noforks=T/F     : Whether to avoid forks [False]
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]
    killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [True]
    forksleep=X     : Sleep time (seconds) between cycles of forking out more process [0]
    rjepy=T/F       : Whether forked commands are rje Python commands [False]
    logfork=T/F     : Whether to log forking in main log [True]
    resfile=LIST    : List of results files (BASEFILE.X) that will need transferring []
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
See also rje.py generic commandline options.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.0.1 - Tweaked the logfork=F log forking output.
    # 0.0.2 - Fixed formatting for Python 2.6 back compatibility for servers.
    # 0.1.0 - Added killmain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [True]
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add MemFree capability.
    # [ ] : Add standalone function to fork out list of commandline commands from file.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('rje_forker', '0.1.0', 'December 2020', '2013')
    description = 'Generic RJE Forking Module'
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
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### SECTION II: Forker Class                                                                                            #
#########################################################################################################################
class Forker(rje_obj.RJE_Object):
    '''
    Forker Class. Author: Rich Edwards (2013).

    Str:str
    - ForkDir = Alternative directory for forking output. [`./`]
    
    Bool:boolean
    - LogFork = Whether to log forking (if quick, might simply be a counter in main script)
    - PIDCheck = Whether to generate a `*.pid` file to check PID status.
    - RjePy = Whether the jobs being run are rje_python jobs and therefore need log files.

    Int:integer
    - IOLimit = Limit of number of IOErrors before termination [50]

    Num:float
    - ForkSleep = Sleep time (seconds) between cycles of forking out more process [0]
    - KillForks = Time to monitor for killing of hanging forks [36000]
    - KillMain=T/F    : Whether to kill main thread rather than individual forks when killforks reached. [True]
    - KillTime = Monitors start of time period to be assessed for hanging forks [time.time()]
    - MemFree = Minimum % memory that should be free to start new fork. (Not yet implemented.)

    List:list
    - Forked = List of dictionaries containing data regarding forked processes.
    - ResFile = List of results files (BASEFILE.X) that will need transferring []
    - ToFork = List of data to be given to runFork() method. (runFork() usually replaced by inheriting class)

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = []
        self.boollist = ['PIDCheck','RjePy','LogFork','KillMain']
        self.intlist = ['IOLimit']
        self.numlist = ['MemFree','ForkSleep','KillForks','KillTime']
        self.listlist = ['ToFork','Forked','ResFile']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({})
        self.setBool({'RjePy':False,'LogFork':True,'KillMain':True})
        self.setInt({'IOLimit':50})
        self.setNum({'MemFree':0.0,'ForkSleep':0.0,'KillForks':36000,'KillTime':time.time()})
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
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                self._cmdReadList(cmd,'path',['ForkDir'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['KillMain','LogFork','PIDCheck','RjePy'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['IOLimit'])   # Integers
                self._cmdReadList(cmd,'float',['MemFree','ForkSleep','KillForks']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ToFork'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.parent(): self.baseFile(self.parent().baseFile()) # Set basefile to be the same as parent object
        if self.getNum('MemFree') > 0.0 and self.getBool('Win32'): self.warnLog('Cannot use memfree=X in win32 mode')
        if self.getNum('MemFree') > 0.0 and self.getBool('OSX'): self.warnLog('Cannot use memfree=X in win32 mode')
        self.int['IOError'] = self.int['IOLimit']
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            forkx = len(self.list['ToFork'])
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.forking()
            self.printLog('#FORK','Forking of %s jobs completed.' % (rje.iStr(forkx)),log=self.getBool('LogFork'))
        except:  self.errorLog('Forker.run() Error')
        if self.list['Forked']:
            self.warnLog('%s fork jobs remain unforked.' % rje.iLen(self.list['Forked']))
            return False
        return True
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while len(self.list['Forked']) < self.getNum('Forks') and self.list['ToFork']: self.nextFork()
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def forking(self):  ### Keeps forking out and processing jobs until no more jobs in self.list['Forked'].
        '''Keeps forking out and processing jobs until no more jobs in self.list['Forked'].'''
        ### ~ [1] ~ Start first set of jobs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('PIDCheck') or self.dev(): pidcheck = '%s.pid' % rje.baseFile(self.log.info['LogFile'])    # Set *.pid object to match log
        else: pidcheck = None
        #self.deBug(pidcheck)
        ### ~ [2] ~ Monitor jobs and set next one running as they finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while self.list['Forked']:
            if not self.getBool('LogFork'):
                self.progLog('\r#FORK','Forking jobs: {0} running; {1} remain.'.format(len(self.list['Forked']),rje.iLen(self.list['ToFork'])))
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
                if (self.i() < 0 and self.getBool('KillMain')) or rje.yesNo('Kill Main Thread?'):
                    raise ValueError('%d seconds of main thread inactivity. %d forks still active!' % (self.getNum('KillForks'),len(self.list['Forked'])))
                elif self.i() < 0 or rje.yesNo('Kill hanging forks?'):
                    self.printLog('#KILL','KillForks=%d seconds walltime reached.' % (self.getNum('KillForks')))
                    for fdict in self.list['Forked']:
                        self.printLog('#KILL','Killing Fork %s, PID %d.' % (fdict['ID'],fdict['PID']))
                        os.system('kill %d' % fdict['PID'])
                else: self.setNum({'KillTime':time.time()})
            ## ~ [2b] Sleep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            time.sleep(self.getNum('ForkSleep'))
#########################################################################################################################
    def nextFork(self):  ### Sets an new fork running.
        '''Sets an new fork running.'''
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
            fdict['cmd'] = self.list['ToFork'].pop(0)
            fdict['ID'] = 'Fork -%d' % (len(self.list['ToFork'])+1)
            fdict['FID'] = 'f_%s' % rje.randomString(6)
            if self.getBool('RjePy'):
                fdict['Log'] = '%s%s.log' % (self.getStr('RunPath'),fdict['FID'])
                fdict['cmd'] += ' basefile=%s' % (fdict['Log'])
                fdict['ResFile'] = self.list['ResFile'][0:]
                try: open(fdict['Log'],'w')
                except: self.errorLog('Log problem. Aborting fork.'); return self.endJob(fdict)
            ### ~ [2] ~ Add Fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setNum({'KillTime':time.time()})
            cpid = os.fork()        # Fork child process
            if cpid:                # parent process records pid of child rsh process
                fdict['PID'] = cpid
                self.printLog('\r#FORK','Forking cmd as %s: %d remain; %.1f%% mem free' % (cpid,len(self.list['ToFork']),fdict['Mem']),log=self.getBool('LogFork'),screen=self.getBool('LogFork') or self.v() > 1)
                self.printLog('#FORK','%s cmd: %s' % (cpid,fdict['cmd']),log=self.getBool('LogFork'),screen=self.getBool('LogFork') or self.v() > 1)
            else:                   # child process
                os.system(fdict['cmd'])
                os._exit(0)
        except SystemExit: raise    # Child
        except: self.errorLog('Forker.startFork error')
#########################################################################################################################
    def endFork(self,fdict):   ### Ends fork, tidies and sets new one running
        '''Ends fork, tidies and sets new one running.'''
        try:### ~ [1] ~ End and tidy current job ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'ResFile' in fdict:
                for resfile in fdict['ResFile']:
                    fromfile = '%s.%s' % (fdict['FID'],resfile)
                    if not rje.exists(fromfile): continue #self.warnLog('Results file %s missing!' % fromfile); continue
                    tofile = '%s.%s' % (self.baseFile(),resfile)
                    if rje.exists(tofile): open(tofile,'a').writelines(open(fromfile,'r').readlines()[1:])
                    else: rje.fileTransfer(fromfile,tofile)
            if 'Log' in fdict:
                if 'cmd' in fdict:
                    open(self.log.info['LogFile'],'a').writelines(open(fdict['Log'],'r').readlines()[5:-1])
                    os.unlink(fdict['Log'])
                else: rje.fileTransfer(fdict['Log'],self.log.info['LogFile'])
                if self.getBool('LogFork'):
                    self.printLog('#END','Fork %s ended: log content transferred' % fdict['PID'])
                    self.printLog('#~~#','#~~#',timeout=False)
                #if self.dev(): self.deBug(fdict['Log'])
                #if self.dev(): self.deBug(rje.exists(fdict['Log']))
            elif 'PID' in fdict and string.split('%s' % fdict['PID'])[0] == 'WAIT': pass
            else: self.printLog('#END','Fork %s ended.' % fdict['PID'],log=self.getBool('LogFork'),screen=self.getBool('LogFork') or self.v() > 1)
        except IOError:
            if self.getInt('IOError') == 1: self.errorLog('Forker.endFork IOError limit reached'); raise
            else: self.int['IOError'] -= 1; self.errorLog('Forker.endFork')
        except: self.errorLog('Forker.endFork error')
        self.nextFork()   # Carry on regardless
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def checkMemory(self): ### Checks free memory and returns True/False if OK to run                               #V0.0
        '''Checks free memory and returns True/False if OK to run.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getNum('MemFree') <= 0.0: return True
            ### ~ [2] Check details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.freeMem() < self.getNum('MemFree'): return False
            return True
        except: self.errorLog('Problem with checkMemory(%s)' % node); return False
#########################################################################################################################
    def freeMem(self):  ### Returns proportion of free system memory                                                #V0.0
        if self.getBool('Win32') or self.getBool('OSX'): return 0.0
        try: memdata = string.split(os.popen('free').readlines()[1])
        except: self.errorLog('Unable to read free memory',printerror=False); return 0.0
        return string.atof(memdata[3]) / string.atof(memdata[1])
#########################################################################################################################
### End of SECTION II: Forker Class                                                                                     #
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
    try:#NewClass(mainlog,cmd_list).run()
        print rje_obj.zen(), '\n\n *** No standalone functionality! *** \n\n'

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
