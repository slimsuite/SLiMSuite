#!/usr/bin/python

# rje.py - RJE General object module
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
#
# This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Author contact: <redwards@cabbagesofdoom.co.uk> / School of Biological Sciences, University of Southampton, UK.

"""
Module:       rje
Description:  Contains General Objects for all my (Rich's) scripts
Version:      4.15.0
Last Edit:    27/07/15
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    General module containing Classes used by all my scripts plus a number of miscellaneous methods.
    - Output to Screen, Commandline parameters and Log Files

    Commandline options are all in the form X=Y. Where Y is to include spaces, use X="Y".

General Commandline:
    v=X             : Sets verbosity (-1 for silent) [0]
    i=X             : Sets interactivity (-1 for full auto) [0]
    log=FILE        : Redirect log to FILE [Default = calling_program.log]
    newlog=T/F      : Create new log file. [Default = False: append log file]
    silent=T/F      : If set to True will not write to screen or log. [False]
    errorlog=FILE   : If given, will write errors to an additional error file. [None]
    help            : Print help to screen

Program-Specific Commands: (Some programs only)
    basefile=FILE   : This will set the 'root' filename for output files (FILE.*), including the log
    outfile=FILE    : This will set the 'root' filename for output files (FILE.*), excluding the log
    delimit=X       : Sets standard delimiter for results output files [\t]
    mysql=T/F       : MySQL output
    append=T/F      : Append to results files rather than overwrite [False]
    force=T/F       : Force to regenerate data rather than keep old results [False]
    backups=T/F     : Whether to generate backup files (True) or just overwrite without asking (False) [True]
    maxbin=X        : Maximum number of trials for using binomial (else use Poisson) [-]
        
System Commandline:
    win32=T/F       : Run in Win32 Mode [False]
    osx=T/F         : Run in MacOSX Mode [False]
    pwin            : Run in PythonWin (** Must be 'commandline', not in ini file! **)
    cerberus        : Run on Cerberus cluster at RCSI
    memsaver=T/F    : Some modules will have a memsaver option to save memory usage [False]
    runpath=PATH    : Run program from given path (log files and some programs only) [path called from]
    rpath=PATH      : Path to installation of R ['R']
    webserver=T/F   : Trigger webserver run and output [False]
    soaplab=T/F     : Implement special options/defaults for SoapLab implementations [False]
    rest=X          : Variable that sets the output to be returned by REST services [None]

Forking Commandline:
    noforks=T/F     : Whether to avoid forks [False]
    forks=X         : Number of parallel sequences to process at once [0]
    killforks=X     : Number of seconds of no activity before killing all remaining forks. [36000]

Development Commandline:
    debug=T/F       : Turn on additional debugging prints and prompts [False]
    warn=T/F        : Turn on program integrity check warnings (unless silent) [True]
    test=T/F        : Run additional testing methods and/or produce additional test outputs [False]
    dev=T/F         : Run development-specific code. (Added to keep main coding working during dev) [False]

Classes:
    RJE_Object(log=None,cmd_list=[]):
        - Metclass for inheritance by other classes.
        >> log:Log = rje.Log object
        >> cmd_list:List = List of commandline variables
        On intiation, this object:
        - sets the Log object (if any)
        - sets verbosity and interactive attributes
        - calls the _setAttributes() method to setup class attributes
        - calls the _cmdList() method to process relevant Commandline Parameters   
    Log(itime=time.time(),cmd_list=[]):
        - Handles log output; printing to log file and error reporting
        >> itime:float = initiation time
        >> cmd_list:list of commandline variables
    Info(prog='Unknown',vers='X',edit='??/??/??',desc='Python script',author='Unknown',ptime=None):
        - Stores intro information for a program.
        >> prog:str = program name
        >> vers:str = version number
        >> edit:str = last edit date
        >> desc:str = program description
        >> author:str = author name
        >> ptime:float = starting time of program, time.time()
    Out(cmd=[]):
        - Handles basic generic output to screen based on Verbosity and Interactivity for modules without classes.
        >> cmd:list = list of command-line arguments

Uses general modules: glob, math, os, random, re, resource, string, sys, time, traceback
"""
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on RJE_General01.plx. Simplified Class Names
    # 0.1 - Added comments and changed capitilisation etc.
    # 0.2 - Added RJE_Object class and 'No Log' concept
    # 0.3 - Updated RJE_Object to not need Out object (restrict Out Object to Log object?)
    # 1.0 - Better Documentation to go with GASP V:1.2
    # 1.1 - Added make Path and sorted ErrorLog
    # 1.1 - def errorLog(self, text='Missing text for errorLog() call!',quitchoice=False,printerror=True):
    # 1.2 - Added Norman's error-tracking and 'basefile/outfile' options as Object defaults
    # 2.0 - Major Overhaul of Module. Out, Log and RJE_Object will now inherit RJE_Object_Shell (see also changes below)
    # 2.1 - Added X="Y" options
    # 2.2 - Added confirm=T/F option to choice()
    # 2.3 - Added more delimited text functions
    # 2.4 - Added listFromCommand() function
    # 2.5 - Added Force opt
    # 3.0 - Added self.list and self.dict dictionaries to RJE_Object. Added 'list', 'clist' and 'glist' to cmdRead.
    # 3.1 - Added object save and load methods
    # 3.2 - Added more refined path options
    # 3.3 - Added general delimitedFileOut() method based on slim_pickings.py compileOut()
    # 3.4 - Added cmdReadList for handling lists of "standard" options
    # 3.5 - Added getFileName for interactive input of file names with confirmation/existence options and checkInputFiles()
    # 3.6 - Added dataDict() method for extracting data from delimited file into dictionary
    # 3.7 - Added extra dictionary methods for storage and retrieval of results in a dict['Data'] dictionary
    # 3.8 - Gave module a general tidy up.
    # 3.9 - Added RunPath to control where program is run.
    # 3.10- Added errorlog=FILE option to redirect errors to an additional error file.
    # 3.11- Added dictionary ranking method.
    # 3.12- Added more list methods.
    # 4.0 - Added RJE_ObjectLite for data storage objects with minimal generic methods and attributes.
    # 4.1 - Added a bioware_server switch for reading iniServer/ ini files on bioware.
    # 4.2 - Modified INI reading across the board to look in ../settings/ and look for defaults.ini as well as rje.ini.
    # 4.2 - Enabled handing on -ini FILE in addition to ini=FILE.
    # 4.3 - Added ilist and nlist types to cmdRead for objects. (Lists of integers and floats). Add ratio() function.
    # 4.4 - Added lineFromIndex(target,file,re_index='^(\S+)\s',sortunique=False,xreplace=True).
    # 4.5 - Modified randomString() and added stringShuffle() methods.
    # 4.6 - Added dev and warn options. Fixed -h lack of help.
    # 4.7 - Added self.warn list and self.warnLog() functions to Log object. Modified i=-1 quitchoice to raise not quit.
    # 4.8 - Added perc cmdtype = float that is multiplied by 100.0 if < 1.0. Removed server option from iniCmds().
    # 4.9 - Added rje.slimsuite, which determines the slimsuite home directory from rje.py file path.
    # 4.10- Added osx=T/F option for Mac-specific running options.
    # 4.11- Enabled '\t#' comments in ini files. Modified getStrLC to return '' for 'none' by default. Added listMax().
    # 4.11- Added self.name() to basic object class.
    # 4.12- Added 'bool' and 'str' to _cmdRead() to ease switchover to new RJE_Objects.
    # 4.13.0 - Added new built-in attributes/options for REST services.
    # 4.13.1 - Fixed MemSaver typo in WarnLog output. Modified mkDir() to avoid clashes raising errors.
    # 4.13.2 - Removed excess REST HTML methods.
    # 4.13.3 - Added uselower=False to dataDict() method.
    # 4.13.4 - Added maxrep=X to listCombos() method.
    # 4.14.0 - Added listToDict() method.
    # 4.14.1 - Fixed matchExp method to be able to handline multilines. (Shame re.DOTALL doesn't work!)
    # 4.14.2 - Modified integer commands to read/convert floats.
    # 4.15.0 - Added intList() and numList() functions.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Split general functions into groups, like delimited text functions
    '''
#########################################################################################################################
import glob, math, os, pickle, random, re, string, sys, time, traceback, urllib2
try:
   set
except NameError:
   from sets import Set as set
try:
    import resource
except: pass    # Not UNIX
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
ini_dir = os.path.join(slimsuitepath,'settings/')
docs_dir = os.path.join(slimsuitepath,'docs/')
#ini_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),'../settings/')
#########################################################################################################################

#########################################################################################################################
###  RJE_Object_Shell Class:                                                                                            #
#########################################################################################################################
class RJE_Object_Shell(object):     ### Metaclass for inheritance by other classes
    '''Metaclass for inheritance by other classes. Forms shell for RJE_Object, Out and Log objects.'''
    ### ~ [1] General Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    cmd_list = []   # List of commandline parameters from which options are parsed.
    log = None      # Log Object (if any) for handling errors etc.
    ### ~ [2] Basic Attribute Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    info = {}       # Stores information variables, such as names and decriptions
    infolist = ['Name','Basefile','Path','Delimit','RunPath','RPath','ErrorLog','Rest']
    stat = {}       # Stores numeric variables
    statlist = ['Verbose','Interactive']
    opt = {}        # Stores boolean variables
    optlist = ['DeBug','Win32','PWin','Memsaver','Append','MySQL','Force','Pickle','SoapLab','Test','Backups','Silent',
               'Webserver','ProgLog','Dev','Warn','OSX']
    obj = {}        # Stores a dictionary of other RJE_Objects 'owned' by object 
    objlist = []
    list = {}       # Dictionary of list-type attributes for object
    listlist = []
    dict = {}       # Dictionary of dictionary attributes for object
    dictlist = []   
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def __init__(self,log=None,cmd_list=[],parent=None):  
        '''
        RJE_Object:
        > log:Log = rje.Log object
        > cmd_list:List = List of commandline variables

        On intiation, this object:
        - sets the Log object (can be None)
        - sets verbosity and interactive attributes
        - calls the _setAttributes() method to setup class attributes
        - calls the _cmdList() method to process relevant Commandline Parameters       
        '''
        ### ~ [1] Setup Log Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.cmd_list = cmd_list
        if not log: self.log = Log(cmd_list=cmd_list)
        else: self.log = log
        ### ~ [2] Setup Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setGeneralAttributes()        # Sets general attributes used (potentially) by all objects
        self.obj['Parent'] = parent
        self._setAttributes()               # Specific method in other classes that set own attributes
        ### ~ [3] Commandline Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._cmdList()                     # Read attribute settings from user-defined parameters
#########################################################################################################################
    def _setGeneralAttributes(self):    ### Sets general attributes for use in all classes
        '''Sets general attributes for use in all classes.'''
        self.info = {'Name':'None','Basefile':'None','Delimit':getDelimit(self.cmd_list),'Rest':'None',
                     'RunPath':makePath(os.path.abspath(os.curdir)),'ErrorLog':'None',
                     'RPath':'R'}
        self.info['Path'] = makePath(os.path.abspath(string.join(string.split(sys.argv[0],os.sep)[:-1]+[''],os.sep)))
        self.stat = {'Verbose':1,'Interactive':0}
        self.opt = {'DeBug':False,'Win32':False,'PWin':False,'MemSaver':False,'Append':False,'MySQL':False,'Force':False,
                    'Pickle':True,'SoapLab':False,'Test':False,'Backups':True,'Silent':False,'Webserver':False,
                    'ProgLog':True,'Warn':True,'Dev':False,'OSX':False}
        self.list = {}
        self.dict = {'Output':{}}
        self.obj = {'DB':None}
#########################################################################################################################
    def _setForkAttributes(self):       ### Sets forking attributes for use in all classes
        '''Sets general forking attributes for use in all classes.'''
        self.statlist += ['Forks','KillForks']
        self.stat['Forks'] = 0
        self.stat['KillForks'] = 36000
        self.optlist += ['NoForks']
        self.opt['NoForks'] = False
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object. Should be replaced in other objects. (See RJE_Object.)
        '''
        Sets Attributes of Object:
        - Info, Stat, Opt and Obj
        '''
        self.infolist = []
        self.statlist = []
        self.optlist = []
        self.objlist = []
        self.listlist = []
        self.dictlist = []
#########################################################################################################################
    def interactive(self): return self.stat['Interactive']
    def i(self): return self.stat['Interactive']
    def v(self): return self.stat['Verbose']
    def force(self): return self.opt['Force']
    def test(self): return self.opt['Test']
    def server(self): return self.opt['Webserver']
    def warn(self): return self.opt['Warn'] and not self.opt['Silent']
    def dev(self): return self.opt['Dev']
    def debugging(self): return self.opt['DeBug']
    def parent(self): return self.obj['Parent']
    def me(self): return string.split('%s' % self)[0][1:]
#########################################################################################################################
    def yesNo(self,text='',default='Y',confirm=False,i=0):
        if self.i() < i: return {'Y':True,'N':False}[default.upper()]
        else: return yesNo(text,default,confirm=False)
#########################################################################################################################
    def db(self,table=None):    ### Return rje_dbase Database object, if one associated with Object
        try:
            db = self.obj['DB']
            if table: return db.getTable(table)
            else: return db
        except: return None
#########################################################################################################################
    def data(self,table,strict=False):
        try: return self.obj['DB'].getTable(table).data()
        except:
            if strict: self.errorLog('No DB table "%s"?' % table); raise
            else: return {}
#########################################################################################################################
    def warnChecks(self):   ### Checks certain input paths etc. and warns of anomalies
        '''Checks certain input paths etc. and warns of anomalies.'''
        if self.warn(): pass
#########################################################################################################################
    ### <2> ### Command-line parameters                                                                                 #
#########################################################################################################################
    def _generalCmd(self,cmd=''):     ### Sets General Attributes from commandline
        '''
        Sets general attributes according to commandline parameters:
        - see rje.__doc__ or run with 'help' option
        '''
        try:
            self._cmdRead(cmd,type='int',att='Verbose',arg='v')
            self._cmdRead(cmd,type='int',att='Interactive',arg='i')
            self._cmdReadList(cmd,type='int',attlist=['Verbose','Interactive'])
            self._cmdReadList(cmd,type='stat',attlist=['MaxBin'])
            self._cmdRead(cmd,type='file',att='Basefile',arg='outfile')
            self._cmdRead(cmd,type='info',att='Rest',arg='outfmt')
            self._cmdReadList(cmd,'info',['Rest'])
            self._cmdReadList(cmd,'file',['Basefile','RPath','ErrorLog'])
            self._cmdReadList(cmd,'abspath',['Path','RunPath'])
            self._cmdRead(cmd,type='opt',att='Win32',arg='pwin')
            self._cmdReadList(cmd,'opt',['DeBug','Win32','PWin','MemSaver','Append','Force','MySQL','Pickle','Test',
                                         'SoapLab','Backups','Webserver','ProgLog','Dev','Warn','OSX'])
        except:
            self.deBug(self.cmd_list)
            if self.log: self.log.errorLog('Problem with %s.cmd:%s' % (self.me(),cmd))
            else: print 'ERROR! Problem with %s.cmd:%s' % (self.me(),cmd)
#########################################################################################################################
    def _cmdReadList(self,cmd=None,type='info',attlist=[]):     ### Sets self.type[att] from commandline command cmd
        '''
        Sets self.type[att] from commandline command cmd.
        >> cmd:str = commandline command
        >> type:str = type of attribute (info,path,opt,int,float,min,max,list,clist,glist,file)
        >> att:list of attributes (key of dictionary) where commandline argument is att.lower()  []
        '''
        for att in attlist: self._cmdRead(cmd,type,att)
#########################################################################################################################
    def _cmdRead(self,cmd=None,type='info',att=None,arg=None):     ### Sets self.type[att] from commandline command cmd
        '''
        Sets self.type[att] from commandline command cmd.
        >> type:str = type of attribute (info,path,opt,int,float,min,max,list,clist,glist,ilist,nlist,file)
        >> att:str = attribute (key of dictionary)
        >> arg:str = commandline argument[att.lower()]
        >> cmd:str = commandline command
        '''
        ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if arg == None: arg = att.lower()
        cmdarg = string.split(cmd,'=')[0].lower()
        if cmdarg not in [arg,'-%s' % arg]: return
        value = cmd[len('%s=' % arg):]
        value = string.replace(value,'#DATE',dateTime(dateonly=True))
        ### ~ [1] Basic commandline types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if type in ['opt','bool']:
            if value[:1].lower() in ['f','0']: self.opt[att] = False
            else: self.opt[att] = True
        elif type in ['info','str']: self.info[att] = value
        elif type == 'path': self.info[att] = makePath(value)
        elif type == 'abspath': self.info[att] = makePath(os.path.abspath(value))
        elif type in ['fullpath','file']: self.info[att] = makePath(value,wholepath=True)
        elif type == 'int':
            try: self.stat[att] = string.atoi(value)
            except:
                self.stat[att] = int(string.atof(value))
                self.log.warnLog('%s=%s needs integer -> %s=%d' % (arg,value,arg,self.stat[att]))
        elif type in ['float','stat','num']: self.stat[att] = string.atof(value)
        elif type == 'max' and matchExp('^\d+,(\d+)',value): self.stat[att] = string.atoi(matchExp('^\d+,(\d+)',value)[0])
        elif type in ['min','max'] and matchExp('^(\d+)',value): self.stat[att] = string.atoi(matchExp('^(\d+)',value)[0])
        elif type == 'list': self.list[att] = listFromCommand(value)
        elif type == 'lclist': self.list[att] = listLower(listFromCommand(value))
        elif type == 'uclist': self.list[att] = listUpper(listFromCommand(value))
        ### ~ [2] Special types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        elif type == 'perc':
            self.stat[att] = string.atof(value)
            if self.stat[att] < 1.0: self.stat[att] *= 100.0
        elif type in ['ilist','nlist']:
            nlist = listFromCommand(value)
            self.list[att] = []
            for nvalue in nlist:
                if type == 'ilist': self.list[att].append(string.atoi(nvalue))
                else: self.list[att].append(string.atof(nvalue))
        elif type == 'clist':   # Returns a *string* of a CSV list
            self.info[att] = string.join(listFromCommand(value),',')
        elif type == 'glist':   # 'Glob' List - returns a list of files using wildcards & glob
            globlist = string.split(value,',')
            self.list[att] = []
            for g in globlist:
                newg = glob.glob(g); newg.sort()
                self.list[att] += newg
        elif type in ['cdict','cdictlist']:   # Converts a CSV list into a dictionary
            self.dict[att] = {}
            if type == 'cdictlist': self.list[att] = []
            clist = listFromCommand(value)
            for c in clist:
                data = string.split(c,':')
                if len(data) == 2:
                    self.dict[att][data[0]] = data[1]
                    if type == 'cdictlist': self.list[att].append(data[0])
        else: raise ValueError
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try: self._generalCmd(cmd)
            except: self.log.errorLog('Problem with %s.cmd:%s' % (self.me(),cmd))
#########################################################################################################################
    ### <3> ### Input/Output                                                                                            #
#########################################################################################################################
    def vPrint(self,text,v=1): return self.verbose(v,text=text)  
#########################################################################################################################
    def verbose(self,v=0,i=None,text='',newline=1):     ### Prints text to screen and maybe pauses program.
        '''
        Prints text to screen if verbosity high enough. Pauses program if interactivity high enough.
        >> v:int = verbosity cut-off for statement
        >> i:int = interactivity cut-off for pause
        >> text:str = text to be printed
        >> newline:int = no. of newlines to follow text. (Pause counts as 1 newline)
        '''
        if i == None: i = self.stat['Interactive'] + 1
        if not self.opt['Silent'] and (self.stat['Verbose'] >= v or self.stat['Interactive'] >= i):
            print text,
            if self.stat['Interactive'] >= i:
                raw_input(" <ENTER> to continue.")
                if 'pwin' not in sys.argv + self.cmd_list: newline -= 1
            while newline > 0:
                print
                newline -= 1
            try: sys.stdout.flush()
            except: pass
#########################################################################################################################
    def debug(self,text): self.deBug(text)
#########################################################################################################################
    def bugPrint(self,text): self.deBug(text,pause=False)
    def devPrint(self,text):
        if self.dev(): self.deBug(text,pause=False)
#########################################################################################################################
    def deBug(self,text,pause=True):   ### Prints text to screen if self.debug 
        '''
        Prints text to screen if self.opt['DeBug'].
        >> text:str = Debugging text to print.
        '''
        pause = pause and self.i() >= 0
        try:
            if self.opt['DeBug'] and pause: self.verbose(self.v(),self.i(),text,1)
            elif self.opt['DeBug']: self.verbose(self.v(),self.i()+1,text,1)
        except KeyboardInterrupt:
            if yesNo('Interrupt program?'): raise
            if yesNo('Switch off Debugging?'): self.opt['DeBug'] = False
#########################################################################################################################
    def pickleMe(self,basefile=None,gzip=True,replace=True):   ### Saves self object to pickle and zips
        '''
        Saves self object to pickle and zips.
        >> basefile:str [None] = if none, will use self.info['Basefile']
        >> gzip:bool [True] = whether to GZIP (win32=F only)
        >> replace:bool [True] = whether to replace existing Pickle
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Pickle']: self.printLog('#PICKLE','%s pickling disabled. (pickle=F)' % self.prog()); return False
            if not basefile or basefile.lower() == 'none': basefile = self.info['Basefile']
            pfile = '%s.pickle' % basefile
            if not replace and (os.path.exists(pfile) or os.path.exists('%s.gz' % pfile)):
                self.printLog('#PICKLE','No pickling - pickle already exists (no replace)!'); return False
            ### ~ [2] ~ Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#SAVE','Attempting to save %s to %s.' % (self.prog(),pfile))
            pickle.dump(self,open(pfile,'w'))
            self.printLog('\r#SAVE','%s Intermediate saved as %s (Python pickle).' % (self.prog(),pfile))
            ### ~ [3] ~ GZip and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['Win32'] and gzip:
                try:
                    if os.path.exists('%s.gz' % pfile): os.unlink('%s.gz' % pfile)
                    os.system('gzip %s' % pfile)
                    self.printLog('\r#GZIP','%s zipped.' % pfile)
                except: self.errorLog('Cannot gzip %s' % pfile)
            return True
        except: self.errorLog('Problem during %s.pickleMe()' % self); return False
#########################################################################################################################
    def unpickleMe(self,basefile=None,process=True): ### (Unzips and) loads pickle object and returns.
        '''(Unzips and) loads pickle object and returns.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basefile or basefile.lower() == 'none': basefile = self.info['Basefile']
            pfile = '%s.pickle' % basefile
            gzfile = '%s.gz' % pfile
            if not (os.path.exists(pfile) or os.path.exists(gzfile)): return None
            if not self.opt['Pickle']: self.printLog('#PICKLE','No loading of pickle. (pickle=F)'); return None
            newme = None    ## New object, loaded from Pickle
            ### ~ [2] ~ Check for file and load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(gzfile) and isYounger(gzfile,pfile) != pfile:
                if self.opt['Win32']: self.errorLog('Cannot unzip %s (win32=T)' % (gzfile)); return None
                else:
                    if os.path.exists(pfile): os.unlink(pfile)
                    try:
                        self.progLog('\r#GUNZIP','Unzipping %s...' % gzfile)
                        gu = os.popen('gunzip %s' % gzfile).read()
                        self.printLog('\r#GUNZIP','Pickle %s unzipped.' % gzfile)
                    except: self.errorLog('Cannot unzip %s' % (gzfile)); return None
            if os.path.exists(pfile):
                self.progLog('\r#LOAD','Attempting to load %s.' % pfile)
                newme = pickle.load(open(pfile,'r'))
                self.printLog('\r#LOAD','%s Intermediate loaded: %s.' % (self.prog(),pfile))
                if not self.opt['Win32']:
                    try:
                        if os.path.exists(gzfile): os.unlink(gzfile)
                        os.system('gzip %s' % pfile)
                        self.printLog('\r#GZIP','%s Pickle %s (re)zipped.' % (self.prog(),pfile))
                    except: self.errorLog('Cannot (re)gzip %s' % pfile)
            if not newme: return None
            ### ~ [3] ~ Process and End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if process: return self.processPickle(newme)
            else: return newme
        except: self.errorLog('Problem during %s.unpickleMe()' % self); return None
#########################################################################################################################
    def processPickle(self,newme):  ### Changes attributes accordingly
        '''Changes attributes accordingly. Replace this method in subclasses.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newme.setLog(self.log)      # Replaces old log with new log.
            if 'Silent' not in newme.opt: newme.opt['Silent'] = False
            for obj in self.obj.values():
                if not obj: continue
                if obj and 'Silent' not in obj.opt: obj.opt['Silent'] = False
                for obj2 in obj.obj.values():
                    if not obj2: continue
                    if obj2 and 'Silent' not in obj2.opt: obj2.opt['Silent'] = False
                    for obj3 in obj2.obj.values():
                        if obj3 and 'Silent' not in obj3.opt: obj3.opt['Silent'] = False
            #newme.cmd_list = self.cmd_list
            #newme.setInfo(self.info)
            #newme.setStat(self.stat)
            #newme.setOpt(self.opt)
            #newme.list = self.list
            #newme.dict = self.dict
            return newme
        except: self.errorLog('Problem during %s.processPickle()' % self); return None
#########################################################################################################################
    def replaceMe(self,newme):  ### Replaces all my attributes with those of newme
        '''Replaces all my attributes with those of newme.'''
        self.info = newme.info
        self.opt = newme.opt
        self.stat = newme.stat
        self.list = newme.list
        self.dict = newme.dict
        self.obj = newme.obj
#########################################################################################################################
### End of RJE_Object_Shell Class                                                                                       #
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
###  RJE_Object Class:                                                                                                  # 
#########################################################################################################################
class RJE_Object(RJE_Object_Shell):     ### Metaclass for inheritance by other classes
    '''Metaclass for inheritance by other classes.'''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
    ### => Now mostly handled in inherited RJE_Object_Shell Class                                                       #
#########################################################################################################################
    def prog(self): return self.log.info['Name']
    def name(self): return self.getStr('Name',default='%s' % self)
#########################################################################################################################
    def baseFile(self,newbase=None,runpath=False): return self.basefile(newbase,runpath)
    def basefile(self,newbase=None,runpath=False):
        if newbase: return self.setBasefile(newbase)
        if self.getStr('Basefile').lower() in ['','none']:
            if runpath: return self.getStr('RunPath')
            else: return ''
        fullpath = os.path.abspath(self.getStr('Basefile')) == self.getStr('Basefile')
        if runpath and not fullpath and not self.getStr('Basefile').startswith(self.getStr('RunPath')):
            return '%s%s' % (self.getStr('RunPath'),self.getStr('Basefile'))    #!# Check that this does not break stuff!
        else: return self.getStr('Basefile')
    def setBasefile(self,basefile=None,cascade=True): ### Sets basefile and cascades to daughter objects
        '''Sets basefile and cascades to daughter objects.'''
        if basefile: self.setStr({'Basefile':basefile})
        else: basefile = self.basefile(runpath=False)
        for obj in self.obj.values():
            try:
                if obj and cascade and obj.basefile() != basefile: obj.setBasefile(basefile)
            except: pass
#########################################################################################################################
    def _setDefaults(self,info='None',opt=False,stat=0.0,obj=None,setlist=False,setdict=False):     ### Default defaults!
        '''
        Sets default defaults.
        >> info:str = default info setting
        >> opt:bool = default opt setting
        >> stat:float = default stat setting
        >> obj:Object = default obj setting
        >> setlist:boolean = whether to set list attributes [False]
        >> setdict:boolean = whether to set dictionary attributes [False]
        '''
        for i in self.infolist: self.info[i] = info
        for o in self.optlist: self.opt[o] = opt
        for s in self.statlist: self.stat[s] = stat
        for j in self.objlist: self.obj[j] = obj
        if setlist:
            for l in self.listlist: self.list[l] = []
        if setdict:
            for d in self.dictlist: self.dict[d] = {}
#########################################################################################################################
    def progLog(self, id='#ERR', text='Log Text Missing!',screen=True,rand=0.0,clear=0):
        if self.v() < 1: return
        if rand > 0 and random.random() > rand: return
        if 'ProgLog' in self.opt and not self.opt['ProgLog']: return False
        return self.printLog('\r%s' % id,text,screen=screen,log=False,newline=False,clear=clear)
    def printLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True,clear=0):
        return self.log.printLog(id,text,timeout,screen and not self.opt['Silent'],log and not self.opt['Silent'],newline,clear=clear)
    def errorLog(self, text='Missing text for errorLog() call!',quitchoice=False,printerror=True,nextline=True,log=True,errorlog=True,warnlist=True):
        #try:
        return self.log.errorLog(text,quitchoice,printerror,nextline,log,errorlog,warnlist)
        #except: return self.log.errorLog(text,quitchoice,False,nextline,log,errorlog)
    def warnLog(self,message,warntype=None,quitchoice=False,suppress=False,dev=False,screen=True):
        return self.log.warnLog(message,warntype,quitchoice,suppress,dev,screen)
#########################################################################################################################
    ### <2> ### Command-line parameters                                                                                 #
    ### => Now mostly handled in inherited RJE_Object_Shell Class                                                       #
#########################################################################################################################
    def _forkCmd(self,cmd=''):     ### Sets General Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see rje.__doc__ or run with 'help' option
        '''
        try:
            self._cmdReadList(cmd,'int',['Forks','KillForks'])
            self._cmdRead(cmd,type='opt',att='NoForks')
            self._cmdRead(cmd,type='opt',att='NoForks',arg='nofork')
        except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <3> ### Object Attributes                                                                                       #
#########################################################################################################################
    def details(self):  ### Prints Details to screen
        '''Returns object details as text.'''
        try:
            ### <a> ### Summary
            if self.info.has_key('Type'): details = '%s (%s)\n' % (self.info['Name'], self.info['Type'])
            else: details = '%s\n' % self.info['Name']
            ### <b> ### Info
            for info in self.infolist:
                if info != 'Name' and info != 'Type' and self.info[info] != '':
                    details += '%s: %s\n' % (info,self.info[info])
            ### <c> ### Stats
            if len(self.statlist) > 0:
                details += 'Stats: [ '
                for stat in self.statlist: details += '%s: %s; ' % (stat,self.stat[stat])
                details += ']\n'
            ### <d> ### Options
            if len(self.optlist) > 0:
                details += 'Options: [ '
                for option in self.optlist: details += '%s: %s; ' % (option,self.opt[option])
                details += ']\n'
            ### <e> ### Lists
            if len(self.listlist) > 0:
                details += 'Lists: [\n'
                for _list in self.listlist: details += '%s: %s;\n' % (_list,self.list[_list])
                details += ']\n'
            ### <f> ### Dictionaries
            if len(self.dictlist) > 0:
                details += 'Dictionaries: [\n'
                for dict in self.dictlist: details += '%s: %s;\n' % (dict,self.dict[dict])
                details += ']\n'
            ### <g> ### Objects
            if len(self.objlist) > 0:
                for object in self.objlist:
                    try:
                        if not self.obj[object]: details += '- %s: None\n' % object
                        else:
                            details += '- %s: %s (%s)\n' % (object, self.obj[object].info['Name'],self.obj[object]) #self.obj[object].details())
                            #!# Add option to add details? - Be careful of circular additions. Store objlist? #!#
                    except: self.log.errorLog('%s (%s) not an acceptable object.' % (self.obj[object],object),False,False)
            details += '\n'
            ### <f> ### Return
            return details
        except:
            self.log.errorLog('Major problem with details()')
            return 'Major problem with details(): %s\n' % sys.exc_info()[0]
#########################################################################################################################
    def attDetails(self,types=['All'],printblanks=True):     ### Prints Details to screen
        '''Returns object details as text.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info.has_key('Type'): details = '%s (%s)\n' % (self.info['Name'], self.info['Type'])
            else: details = '%s\n' % self.info['Name']
            ### ~ [1] ~ Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Info' in types or 'All' in types:
                for info in sortKeys(self.info):
                    if info in ['Name','Type']: continue
                    if not self.info[info] and not printblanks: continue
                    details += '%s: %s\n' % (info,self.info[info])
            ### ~ [2] ~ Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Stat' in types or 'All' in types:
                if self.stat:
                    details += 'Stats: [ '
                    for stat in sortKeys(self.stat): details += '%s: %s; ' % (stat,self.stat[stat])
                    details += ']\n'
            ### ~ [3] ~ Options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Opt' in types or 'All' in types:
                if self.opt:
                    details += 'Options: [ '
                    for opt in sortKeys(self.opt): details += '%s: %s; ' % (opt,self.opt[opt])
                    details += ']\n'
            ### ~ [4] ~ Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'List' in types or 'All' in types:
                for _list in sortKeys(self.list):
                    if not self.list[_list] and not printblanks: continue
                    details += '%s: %s\n' % (_list,self.list[_list])
            ### ~ [5] ~ Dicts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Dict' in types or 'All' in types:
                for _dict in sortKeys(self.dict):
                    if not self.dict[_dict] and not printblanks: continue
                    details += '%s: %s\n' % (_dict,self.dict[_dict])
            ### ~ [6] ~ Obj __~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ('Obj' in types or 'All' in types) and self.obj:
                for object in sortKeys(self.obj):
                    try:
                        if not self.obj[object] and printblanks: details += '- %s: None\n' % object
                        else:
                            details += '- %s: %s (%s)\n' % (object, self.obj[object].info['Name'],self.obj[object]) #self.obj[object].details())
                            #!# Add option to add details? - Be careful of circular additions. Store objlist? #!#
                    except: self.log.errorLog('%s (%s) not an acceptable object.' % (self.obj[object],object),False,False)
            details += '\n'
            ### ~ [7] ~ Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return details
        except:
            self.errorLog('Major problem with attDetails()')
            return 'Major problem with attDetails(): %s\n' % sys.exc_info()[0]
#########################################################################################################################
    def edit(self):     ### Options to change details
        '''
        Options to change all object details (Info, Stat, Opt) and associated objects' details too.
        #!# Lists and Dictionaries not included. #!#
        '''
        try:
            print '\nEdit %s Attributes.\nEnter new values or leave Blank to retain.\n' % self.info['Name']
            for info in self.infolist: self.info[info] = self._editChoice(info,self.info[info])
            for stat in self.statlist: self.stat[stat] = self._editChoice(stat,self.stat[stat],numeric=True)
            for opt in self.optlist: self.opt[opt] = self._editChoice(opt,self.opt[opt],boolean=True)
            for object in self.objlist:
                if yesNo('Edit %s object?' % object): self.obj[object].edit()
            if self.info['Name'] != '' and yesNo('%s\nOK?: ' % self.details()): return
            self.edit()
        except: self.log.errorLog('Major problem with Object edit()')
#########################################################################################################################
    def _editChoice(self,text,value,numeric=False,boolean=False):  ### Returns Current or New Value as chosen
        '''
        Returns Current or New Value as chosen.
        >> text:str = Text to dislay for choice
        >> value:str/int/float/bool = Existing value
        >> numeric:bool = whether a numeric value is wanted
        >> boolean:bool = whether a True/False value is wanted
        '''
        try:
            if numeric: return getFloat(text=text,default=value)
            elif boolean: return getBool(text=text,default=value)
            else: return choice(text=text,default=value)
        except: self.log.errorLog('Major problem with _editChoice()')
#########################################################################################################################
    def getAtt(self,type,key,default=None): return self.getAttribute(type,key,default)
    def getAttribute(self,type,key,default=None):    ### Gets object information of correct type
        '''
        Gets object information of correct type.
        >> type:str = 'info','list','dict','opt','stat' or 'obj'
        >> key:str = key for given attribute dictionary
        >> default:anything = what to return if type is wrong or key not found in type
        '''
        ### Setup ###
        att = {}
        if type in ['info','str']: att = self.info
        elif type in ['opt','bool']: att = self.opt
        elif type in ['int','num','stat']: att = self.stat
        elif type == 'list': att = self.list
        elif type == 'dict': att = self.dict
        elif type == 'obj': att = self.obj
        ### Return ###
        if att.has_key(key): return att[key]
        else: return default
#########################################################################################################################
    def setAttribute(self,type,key,newvalue):    ### Sets object information of correct type from string
        '''
        Gets object information of correct type.
        >> type:str = 'info','list','dict','opt','stat' or 'obj'
        >> key:str = key for given attribute dictionary
        >> newvalue:str = string version of new value
        '''
        self._cmdRead('%s=%s' % (key.lower(),newvalue),type=type,att=key)
#########################################################################################################################
    def getStrLC(self,ikey=None,return_none=False):
        istr =  self.getInfo(ikey).lower()
        if istr != 'none' or return_none: return istr
        return ''
    def getStrUC(self,ikey=None,return_none=False):
        istr =  self.getInfo(ikey).upper()
        if istr != 'none' or return_none: return istr
        return ''
    def getStr(self,ikey=None,default='None',checkdata=True): return self.getInfo(ikey,default,checkdata)
#########################################################################################################################
    def getInfo(self,ikey=None,default='None',checkdata=True):  ### Gets object information or returns default if not found
        '''
        Gets object information or returns default if not found.
        >> ikey:str = key for self.info dictionary
        >> default:str = values returned if ikey not found.
        >> checkdata:boolean = whether to check self.dict['Data'] if missing
        '''
        if not ikey: return self.info
        if self.info.has_key(ikey): return self.info[ikey]
        elif checkdata and self.dict.has_key('Data'):
            return getFromDict(self.dict['Data'],ikey,returnkey=False,case=False,default=default)
        return default
#########################################################################################################################
    def getBool(self,ikey=None,default=False): return self.getOpt(ikey,default)
#########################################################################################################################
    def getOpt(self,okey=None,default=False):    ### Gets object opt or returns default
        '''Gets object opt or returns default.'''
        if not okey: return self.opt
        if self.opt.has_key(okey): return self.opt[okey]
        return default
#########################################################################################################################
    def getInt(self,ikey=None,default=0,checkdata=True):
        if ikey: return int(self.getStat(ikey,default,checkdata))
        else: return self.stat
#########################################################################################################################
    def getNum(self,ikey=None,default=0.0,checkdata=True): return self.getStat(ikey,default,checkdata)
#########################################################################################################################
    def getStat(self,skey=None,default=0,checkdata=True):  ### Gets object stat or returns default if not found
        '''
        Gets object stat or returns default if not found.
        >> skey:str = key for self.stat dictionary
        >> default:num = values returned if ikey not found.
        >> checkdata:boolean = whether to check self.dict['Data'] if missing
        '''
        if not skey: return self.stat
        if self.stat.has_key(skey): return self.stat[skey]
        elif checkdata and self.dict.has_key('Data'):
            val =  getFromDict(self.dict['Data'],skey,returnkey=False,case=False,default=default)
            if val != default: return string.atof(val)
        return default
#########################################################################################################################
    def getDict(self,dkey=None,dkeykey=None,default=None):    ### Returns value of self.dict[dkey][dkeykey] else default
        '''Returns value of self.dict[dkey][dkeykey] else default.'''
        if not dkey: return self.dict
        if not dkeykey: return self.dict[dkey]
        if self.dict.has_key(dkey) and self.dict[dkey].has_key(dkeykey): return self.dict[dkey][dkeykey]
        return default
#########################################################################################################################
    def setInfo(self,infodic,addtolist=True):  ### Sets object information
        '''
        Sets object information.
        >> infodic:dictionary with keys corresponding to self.info keys
        >> addtolist:boolean = whether to add to self.infolist if missing
        '''
        try:
            for key in infodic.keys():
                self.info[key] = infodic[key]
                if addtolist and key not in self.infolist: self.infolist.append(key)
        except: self.log.errorLog('Problem with setInfo()',True)
#########################################################################################################################
    def setStat(self,statdic,addtolist=True):  ### Sets object Parameters
        '''
        Sets object information.
        >> statdic:dictionary with keys corresponding to self.stat keys
        >> addtolist:boolean = whether to add to self.statlist if missing
        '''
        try:
            for key in statdic.keys():
                self.stat[key] = statdic[key]
                if addtolist and key not in self.statlist: self.statlist.append(key)
        except: self.log.errorLog('Problem with setStat()',True)
#########################################################################################################################
    def setOpt(self,optdic,addtolist=True):  ### Sets object Options
        '''
        Sets object information.
        >> optdic:dictionary with keys corresponding to self.opt keys
        >> addtolist:boolean = whether to add to self.optlist if missing
        '''
        try:
            for key in optdic.keys():
                self.opt[key] = optdic[key]
                if addtolist and key not in self.optlist: self.optlist.append(key)
        except: self.log.errorLog('Problem with setOpt()',True)
#########################################################################################################################
    def setStr(self,attdic,addtolist=True): self.setInfo(attdic,addtolist)
    def setInt(self,attdic,addtolist=True): self.setStat(attdic,addtolist)
    def setNum(self,attdic,addtolist=True): self.setStat(attdic,addtolist)
    def setBool(self,attdic,addtolist=True): self.setOpt(attdic,addtolist)
#########################################################################################################################
    def setList(self,listdic,addtolist=True):  ### Sets object Lists
        '''
        Sets object list attributes.
        >> listdic:dictionary with keys corresponding to self.list keys
        >> addtolist:boolean = whether to add to self.listlist if missing
        '''
        try:
            for key in listdic.keys():
                self.list[key] = listdic[key]
                if addtolist and key not in self.listlist: self.listlist.append(key)
        except: self.log.errorLog('Problem with setList()',True)
#########################################################################################################################
    def setDict(self,dictdic,addtolist=True):  ### Sets object Dictionaries
        '''
        Sets object dictionary attributes.
        >> dictdic:dictionary with keys corresponding to self.dict keys
        >> addtolist:boolean = whether to add to self.dictlist if missing
        '''
        try:
            for key in dictdic.keys():
                self.dict[key] = dictdic[key]
                if addtolist and key not in self.dictlist: self.dictlist.append(key)
        except: self.log.errorLog('Problem with setDict()',True)
#########################################################################################################################
    def setObj(self,objdic,addtolist=True):  ### Sets object Objects
        '''
        Sets object information.
        >> objdic:dictionary with keys corresponding to self.obj keys
        >> addtolist:boolean = whether to add to self.objlist if missing
        '''
        try:
            for key in objdic.keys():
                self.obj[key] = objdic[key]
                if addtolist and key not in self.objlist: self.objlist.append(key)
        except: self.log.errorLog('Problem with setObj()',True)
#########################################################################################################################
    def setDictData(self,datadict,dictkey='Data'):  ### Sets object dict values
        '''
        Sets object Data dict values.
        >> datadict = Dictionary of values to add to self.dict[dictkey]
        '''
        try:
            if not self.dict.has_key(dictkey): self.dict[dictkey] = {}
            for key in datadict.keys(): self.dict[dictkey][key] = datadict[key]
        except: self.log.errorLog('Problem with setDictData()',True)
#########################################################################################################################
    def getData(self,dkey,dlist=['stat','info','opt'],case=False,str=True,default=None,dp=-1):  ### Gets data from dict['Data']
        '''
        Gets data from dict['Data'] if it exists and has key, else tries dlist dictionaries.
        >> dkey:str = Key for dictionaries
        >> dlist:list [self.stat,self.info,self.opt] = list of dictionaries to try after dict['Data']
        >> case:bool [False] = whether to match case for dkey
        >> str:bool [True] = whether to return all values as a string
        >> default [None] = what to return if no dictionary has key
        >> dp:int [-1] = Number of decimal places to use for stats if str=True (-1 = return as is)
        << returns value from data or stat/info/opt
        '''
        try:
            ### Setup ###
            dictlist = []
            if self.dict.has_key('Data'): dictlist = [self.dict['Data']]
            ddict = {'stat':self.stat,'info':self.info,'opt':self.opt}
            for dict in dlist:
                if ddict.has_key(dict): dictlist.append(ddict[dict])
                else: dictlist.append(dict)
            ### Look in dictionaries ###
            data = default
            for dict in dictlist:
                if data == default: data = getFromDict(dict,dkey,returnkey=False,case=case,default=default)
                if data != default and dict == self.stat:
                    if dp > 0: data = int(data * (10 ** dp) + 0.5) / float(10 ** dp)
                    elif dp == 0: data = int(data + 0.5)
            ### Return ###
            if str: return '%s' % data
            else: return data
        except: self.log.errorLog('Problem with %s.getData(%s)' % (self,dkey))
        return default
#########################################################################################################################
    def setLog(self,log,cascade=True):  ### Sets given log as log object and cascades through self.obj
        '''Sets given log as log object and cascades through self.obj.'''
        self.log = log
        for obj in self.obj.values():
            if cascade and obj and log != obj.log: obj.setLog(log)
#########################################################################################################################
    ### <4> ### Input/Output
#########################################################################################################################
    def checkInputFiles(self,ilist,ask=True):    ### Checks for existence of Input Files and asks for them if missing
        '''
        Checks for existence of Input Files and asks for them if missing. Raises error if not interactive and missing.
        >> ilist:list of keys for self.info that should be sequence files
        >> ask:boolean = whether to ask for file if missing and i >= 0
        '''
        for seqfile in ilist:
            if not os.path.exists(self.info[seqfile]):
                if not ask or self.stat['Interactive'] < 0:
                    self.log.errorLog('%s file "%s" does not exist!' % (seqfile,self.info[seqfile]),printerror=False)
                    raise ValueError
                else: self.info[seqfile] = getFileName('%s File?' % seqfile)
#########################################################################################################################
    def loadFromFile(self,filename=None,v=0,checkpath=True,chomplines=False):   ### Loads all data from file and returns readlines() list
        '''
        Loads all data from file and returns readlines() list. Will look in same directory and self.info['Path']
        >> filename:str = Name of file
        >> v:int = verbosity setting for loading message
        >> checkpath:boolean [True] = whether to check in self.info['Path'] if filename missing.
        >> chomp:boolean [False] = whether to remove \\r and \\n from lines
        << filelines:list of lines from file
        '''
        try:
            lookinpath = True
            altfile = makePath(self.info['Path']) + filename
            alt2file = makePath(self.info['Path']+'../libraries/') + filename
            alt3file = makePath(self.info['Path']+'../data/') + filename
            if filename.lower() in ['none','']: return []
            if checkForFile(filename): openfile = filename
            elif checkpath and checkForFile(altfile): openfile = altfile
            elif checkpath and checkForFile(alt2file): openfile = alt2file
            elif checkpath and checkForFile(alt3file): openfile = alt3file
            elif checkpath:
                filetxt = filename
                if altfile != filename: filetxt += '(and %s)' % altfile
                if alt2file not in [filename,altfile]: filetxt += '(and %s)' % alt2file
                self.errorLog('File %s not found!' % (filetxt),printerror=False)
                return []
            else:
                self.log.errorLog('File %s not found!' % filename,printerror=False)
                return []
            if v <= self.stat['Verbose']:
                self.log.printLog('#LOAD','Loading data from %s ...' % openfile,newline=False,log=False)
            file_lines = open(openfile, 'r').readlines()
            if len(file_lines) == 1: file_lines = string.split(file_lines[0],'\r')
            if chomplines: file_lines = string.split(chomp(string.join(file_lines,'!#ENDOFLINE#!')),'!#ENDOFLINE#!')
            if v <= self.stat['Verbose']:
                self.log.printLog('\r#LOAD','Loading data from %s complete: %s lines.' % (openfile,integerString(len(file_lines))),log=False)
            return file_lines
        except:
            self.log.errorLog('Major problem in loadFromFile(%s)' % filename)
            return []
#########################################################################################################################
    def saveSelf(self,filename=None,append=False):  ### Saves a dump of own basic data (info,stat,opt) into file
        '''
        Saves a dump of own basic data (info,stat,opt) into file. Superceded by pickling.
        >> filename:str = name of file (else self name.object.txt)
        >> append:boolean = whether to append filename or write anew
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not filename: filename = '%s.object.txt' % self.info['Name']
            if append: OUTFILE = open(filename,'a')
            else: OUTFILE = open(filename,'w')
            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.saveData(OUTFILE)
            OUTFILE.close()
            return True
        except: self.log.errorLog('Major problem in saveSelf(%s)' % filename)
        return False
#########################################################################################################################
    def saveData(self,OUTFILE,skiplist=[]):     ### Saves a dump of basic data into OUTFILE handle
        '''
        Saves a dump of basic data into OUTFILE handle.
        >> OUTFILE:open file handle for writing
        >> skiplist:list of things to skip output for.
        - 'info','opt','stat','list','dict','obj'
        '''
        OUTFILE.write('>#:%s:#\n' % self)
        if 'info' not in skiplist:
            for info in self.infolist: OUTFILE.write('Info:::%s:::%s\n' % (info,self.info[info]))
        if 'stat' not in skiplist:
            for stat in self.statlist: OUTFILE.write('Stat:::%s:::%s\n' % (stat,self.stat[stat]))
        if 'opt' not in skiplist:
            for opt in self.optlist: OUTFILE.write('Opt:::%s:::%s\n' % (opt,self.opt[opt]))
            OUTFILE.write('# Lists #\n')
        if 'list' not in skiplist:
            for obj in self.list.keys():
                for o in self.list[obj]: OUTFILE.write('%s=%s\n' % (obj,o))
            OUTFILE.write('# Dictionaries #\n')
        if 'dict' not in skiplist:
            for obj in self.dict.keys():
                for o in self.dict[obj].keys(): OUTFILE.write('%s=%s:::%s\n' % (obj,o,self.dict[obj][o]))
            OUTFILE.write('# Objects #\n')
        if 'obj' not in skiplist:
            for obj in self.obj.keys(): OUTFILE.write('%s=%s\n' % (obj,self.obj[obj]))
        OUTFILE.write('//\n\n')
#########################################################################################################################
    def loadSelf(self,filename=None):  ### Loads a dump of own basic data (info,stat,opt) into file
        '''
        Loads a saved dump of own basic data (info,stat,opt) into file. Superceded by pickling.
        >> filename:str = name of file
        '''
        try:
            ### Setup ###
            if not os.path.exists(filename): raise IOError
            loadlines = self.loadFromFile(filename,1)
            selfcode = None

            ### <1> ### Build Object Dictionary ###
            obj_dict = {'None':None,'True':True,'False':False}   # Code:Object
            for line in loadlines:
                line = chomp(line)
                if matchExp('^>#:(<(\S+) object \S.+>):#',line):   # New Object
                    (ocode,otype) = matchExp('^>#:(<(\S+) object \S.+>):#',line)
                    if selfcode: obj_dict[ocode] = self.addLoadObj(otype,self.log,self.cmd_list)
                    else:
                        obj_dict[ocode] = self
                        selfcode = ocode

            ### <2> ### Load Details ###
            i = 0
            while i < len(loadlines):
                line = chomp(loadlines[i])
                if matchExp('^>#:(<(\S+) object \S.+>):#',line):   # New Object
                    (ocode,otype) = matchExp('^>#:(<(\S+) object \S.+>):#',line)
                    obj = obj_dict[ocode]
                    if not obj:
                        self.log.errorLog('No Object for %s!' % ocode,False,False)
                        i += 1
                        continue
                    ## Basic details ##
                    i += 1
                    line = chomp(loadlines[i])
                    basics = []
                    while matchExp('^(\S+):::(.+):::(.+)$',line):   # Matched line
                        lmatch = matchExp('^(\S+):::(.+):::(.+)$',line)
                        if lmatch == 'Special': obj.loadSpecial(line)
                        elif lmatch[0] != '#':  # Ignore comment
                            obj.loadData(type=lmatch[0],key=lmatch[1],data=lmatch[2])
                        i += 1
                        line = chomp(loadlines[i])
                    ## Extra ##
                    etype = 'None'  # Lists/Dictionaries/Objects
                    while line != '//': # End!
                        if matchExp('^# (\S+) #',line): etype = matchExp('^# (\S+) #',line)[0]
                        elif line.find('#') == 0:   # Skip comment
                            etype = etype
                        elif matchExp('^(\S+)=(\S.*)$',line):
                            (key,val) = matchExp('^(\S+)=(\S.*)$',line)
                            if obj_dict.has_key(val): val = obj_dict[val]
                            elif matchExp('^<(\S+) object \S.+>',val):
                                self.log.errorLog('Object %s missing from %s.' % (val,filename),False,False)
                                val = None
                            if etype == 'Lists':    # Lists
                                if key in obj.list.keys(): obj.list[key].append(val)
                                else: obj.list[key] = [val]
                                if key not in obj.listlist: obj.listlist.append(key)
                            if etype == 'Objects':    # Objects
                                obj.obj[key] = val
                                if key not in obj.objlist: obj.objlist.append(key)
                            if etype == 'Dictionaries':    # Dictionaries
                                if matchExp('^(\S.*):::(\S*.*)$',val):
                                    (newkey,newval) = matchExp('^(\S.*):::(\S*.*)$',val)
                                    if obj_dict.has_key(newkey): newkey = obj_dict[newkey]
                                    if obj_dict.has_key(newval): newval = obj_dict[newval]
                                    if key in obj.dict.keys(): obj.dict[key][newkey] = newval
                                    else: obj.dict[key] = {newkey:newval}
                                    if key not in obj.dictlist: obj.dictlist.append(key)
                                else: self.log.errorLog('Cannot find key:value for %s:\n%s' % (key,val),False,False)
                                    
                        i += 1
                        if i == len(loadlines): break
                        line = chomp(loadlines[i])
                    obj.selfLoadTidy(obj_dict)  # Special attributes that need more work
                ## Loop ##
                i += 1
            return True
        except: self.log.errorLog('Major problem in loadSelf(%s)' % filename)
        return False
#########################################################################################################################
    def addLoadObj(self,object_type=None,log=None,cmd_list=[]):  ### Returns a new object of the right type
        '''
        Returns a new object of the right type.
        >> object_type:str = Code from object used to identify correct object to create.
        >> log:Log object to feed new object
        >> cmd_list:List of commands for new object
        '''
        try: return None     #!# Edit this subroutine for the class #!#
        except: self.log.errorLog('Oh, the shame! Something has gone wrong in addLoadObj()')
        return None
#########################################################################################################################
    def selfLoadTidy(self,obj_dict):    ### Tidies up unusual data types, such as dictionaries
        '''
        Tidies up unusual data types, such as dictionaries.
        >> obj_dict:Dictionary of object codes and objects
        '''
        try: return True    #!# Add class-specific code #!#
        except: self.log.errorLog('Oh dear! Major problem in %s.selfLoadTidy()' % self)
        return False
#########################################################################################################################
    def loadData(self,type=None,key=None,data=None):  ### Loads a single piece of data 
        '''
        Loads a single piece of data.
        >> type:str = Info/Opt/Stat
        >> key:Key for relevant data type
        >> data:Actual data
        '''
        try:
            ### Input ###
            if type == 'Info': self.info[key] = data
            elif type == 'Stat':
                self.stat[key] = string.atof(data)  #!# Remember to convert to integer later if needed
                if (self.stat[key] - int(self.stat[key])) == 0: self.stat[key] = int(self.stat[key])
            elif type == 'Opt' and data == 'True': self.opt[key] = True
            elif type == 'Opt' and data == 'False': self.opt[key] = False
            elif type[0] != '#':    # Ignore comments
                return self.loadSpecial('%s:::%s:::%s' % (type,key,data))
            return True
        except: self.log.errorLog('Major problem in loadSelf(%s)' % filename)
        return False
#########################################################################################################################
    def loadSpecial(self,line=''):  ### Loads a single piece of data 
        '''
        Loads a single piece of data.
        >> line:str = Special data to load
        '''
        try:
            ### Setup ###
            #!# Replace this method with the appropriate one for 'Special:::Data' input #!#
            self.log.errorLog('%s cannot handle Special input: %s.' % (self.info['Name'],line),False,False)
            return False
        except: self.log.errorLog('Major problem in loadSpecial(%s)' % line)
        return False
#########################################################################################################################
    def gUnzip(self,file,log=True):  ### Unzips a given file
        '''Unzips a given file.'''
        if log: self.log.printLog('#UNZIP','Unzipping %s ...' % file,newline=False,log=False)
        try:
            if not os.path.exists(file): raise IOError
            if self.opt['Win32']: raise ValueError
            os.system('gunzip %s' % file)
            if log: self.log.printLog('\r#UNZIP','Unzipped %s successfully' % file)
        except: self.log.errorLog('Problem unzipping %s' % file)
#########################################################################################################################
    def gZip(self,filename,log=True,unlink=True):    ### Zips file, if appropriate
        '''Zips file, if appropriate.'''
        if self.opt['Win32']: self.printLog('#WIN32','Cannot GZIP %s (Win32=T)' % filename); return False
        if log: self.progLog('#GZIP','Zip %s ...' % filename)
        if unlink and os.path.exists('%s.gz' % filename): os.unlink('%s.gz' % filename)
        if os.path.exists(filename): 
            os.system('gzip %s' % filename)
            self.printLog('\r#GZIP','%s zipped.' % filename); return True
        else: self.printLog('\r#ERR','%s missing!' % filename); return False
#########################################################################################################################
    def needToRemake(self,checkfile,parentfile,checkdate=None,checkforce=True): ### Checks whether checkfile needs remake
        '''
        Checks whether checkfile needs remake.
        >> checkfile:str = File name of file that may need remaking.
        >> parentfile:str = File name of file that was used to make checkfile. (Should be older)
        >> checkdate:bool = whether to bother checking the comparative dates.
        >> checkforce:bool = whether to use self.force() to identify whether remake should be forced
        '''
        if checkforce and self.force(): return True
        if not os.path.exists(checkfile): return True
        if checkdate == None and self.getBool('IgnoreDate',default=False) or checkdate == False: return True
        if isYounger(checkfile,parentfile) != parentfile: return False
        return True
#########################################################################################################################
    ### <4a> ### REST Output                                                                                            #
#########################################################################################################################
    def restOutputError(self,errormsg): ### Returns full error message.
        '''Returns full error message..'''
        try:
            print errormsg
            ### Setup error variables
            error_type = str(sys.exc_info()[0])         # Error Type       : exceptions.IOError
            error_type = string.replace(error_type,'exceptions.','')
            error_value = str(sys.exc_info()[1])        # Error Value      : [Errno 9] Bad file descriptor
            error_traceback = traceback.extract_tb(sys.exc_info()[2])
            error_file = str(error_traceback[-1][0])    # File             : C:\Documents and Settings\normdavey\Desktop\Python\BLAST\Main.py
            error_method = str(error_traceback[-1][2])  # Method           : readFile
            error_line = str(error_traceback[-1][1])    # Line             : 15
            error_error = str(error_traceback[-1][3])   # Error            : for lines in fileIn.readlines():
            ### Output error
            error_txt = []
            error_txt.append('%s\n' % dateTime())
            error_txt.append('%s\n' % errormsg)
            error_txt.append('%s %s\n' % (sys.exc_info()[0],sys.exc_info()[1]))
            error_txt.append('%s: %s\n' % (error_type, error_value))
            error_txt.append('File: %s\n' % error_file)
            error_txt.append('Method: %s (line %s)\n' % (error_method, error_line))
            error_txt.append('Error: %s\n' % error_error)
            error_txt.append('\nContact webmaster for more details.\n')
            return string.join(error_txt,'')
        except: os._exit(0)
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        There is currently no specific help available on REST output for this program. Run with &rest=help for general
        options. Run with &rest=full to get full server output. Individual outputs can be identified/parsed:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ...

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.info key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutput(self,outfmt=None,maxfilesize=0):   ### Controls what is returned by REST server output
        '''Controls what is returned by REST server output.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not outfmt: outfmt = self.getStrLC('Rest')
            if outfmt == 'default': outfmt = self.restOutputDefault()
            if not outfmt: return self  # If no REST output (rest=None), the object itself is returned.
            ### ~ [1] RESTful output to return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfmt == 'log': return open(self.log.info['LogFile'],'r').read()
            if outfmt in ['full','text']: return self.restFullOutput()
            if outfmt == 'version': return Out(self.log,['v=-1']).printIntro(self.log.obj['Info'])
            if outfmt == 'outfmt': return self.restSetup.__doc__
            if outfmt == 'warnings': return string.join(self.log.list['WarnLog'],'\n')  # List of log warning messages.
            if outfmt == 'errors': return string.join(self.log.list['ErrorLog'],'\n')   # List of log error messages.
            if outfmt in self.dict['Output']:
                outdata = self.dict['Output'][outfmt]
                if outdata in self.info: outdata = self.info[outdata]
                if exists(outdata):
                    nbytes = os.path.getsize(outdata)
                    if nbytes > maxfilesize > 0:
                        otext = '%s is too large to return (> %s)' % (outdata,humanByteSize(nbytes))
                        if outfmt == self.getStrLC('Rest'): return 'ERROR: %s.' % otext
                        return '%s in full output. Try retrieve&rest=%s.' % (otext,outfmt)
                    return open(outdata,'r').read()
                elif outdata: return '%s\n' % outdata
                else: return 'No output generated.\n'
            elif self.db(outfmt):
                ext = delimitExt(self.getStr('Delimit'))
                dbfile = '%s.%s.%s' % (self.baseFile(runpath=True),outfmt,ext)
                if not exists(dbfile): self.db(outfmt).saveToFile(dbfile)
                return os.path.abspath(dbfile)
            # summary = Reduced output
            # Default: return full text parsed into tabs
            return '# SERVER ERROR: "&rest=%s" not a recognised output\n%s' % (outfmt,self.restFullOutput())
        except: return self.restOutputError('%s.restOutput() error.' % self.prog())
#########################################################################################################################
    def restOutputOrder(self): return sortKeys(self.dict['Output'])
    def restOutputDefault(self): return 'full'
#########################################################################################################################
    def restFullOutput(self,outfile=None):  ### Full REST server output method. Should be replaced in individual Classes to customise.
        '''Full REST server output method. Should be replaced in individual Classes to customise.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outputs = self.restOutputOrder()
            if not outputs: return 'No full REST output for %s' % self.prog()
            info = self.log.obj['Info']
            outtxt = '# Output for %s V%s: %s\n' % (info.program,info.version,time.asctime(time.localtime(info.start_time)))
            warntxt = '# %s warnings; %s errors' % (iStr(self.log.warnx),iLen(self.log.list['ErrorLog']))
            if self.log.warnx or self.log.list['ErrorLog']: warntxt += ': see log output (below) for details'
            outtxt += '%s.\n' % warntxt
            if 'jobid' in self.dict['Output']: outtxt += '# JobID: %s\n' % self.dict['Output']['jobid']
            if 'jobid' in outputs: outputs.remove('jobid')
            outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            ### ~ [1] RESTful output to return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'status' in outputs: outputs.remove('status')
            outtxt += '# status:\n%s\n' % self.restOutput('status')
            outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            if 'version' not in outputs:
                outtxt += '# version:\n%s\n' % info.version
                outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            ## ~ [1a] Output dictionary contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for outkey in outputs:
                if outkey in ['help','check']: continue
                if outkey not in self.dict['Output'] and self.db(outkey):
                    ext = delimitExt(self.getStr('Delimit'))
                    dbfile = '%s.%s.%s' % (self.baseFile(runpath=True),outkey,ext)
                    if not exists(dbfile): self.db(outkey).saveToFile(dbfile)
                    self.dict['Output'][outkey] = os.path.abspath(dbfile)
                if outkey not in self.dict['Output']:
                    outtxt += '# %s:\nNo output generated.\n' % outkey
                    outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
                    continue
                outdata = self.dict['Output'][outkey]
                try:
                    if self.getStrLC(outdata) and os.path.exists(self.getStr(outdata)): outtxt += '# %s: %s\n' % (outkey,self.getStr(outdata))
                    elif os.path.exists(outdata): outtxt += '# %s: %s\n' % (outkey,outdata)
                    else: outtxt += '# %s:\n' % outkey
                except: outtxt += '# %s:\n' % outkey
                outtxt += self.restOutput(outkey)
                outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            ## ~ [1b] Log output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outtxt += '# log: %s\n%s' % (self.log.info['LogFile'],open(self.log.info['LogFile'],'r').read())
            outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            outtxt += '# warnings:\n'
            if self.log.warnx: outtxt += '%s\n' % self.restOutput('warnings')
            else: outtxt += '0 warning messages.\n'
            #outtxt += '%s warning messages - check log or: retrieve&jobid=%s&rest=warnings\n' % (rje.iStr(self.log.warnx),self.dict['Output']['jobid'])
            if self.log.list['ErrorLog']:
                outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
                outtxt += '# errors:\n'
                outtxt += self.restOutput('errors')
            ## ~ [1c] Return output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if outfile:
                open(outfile,'w').write(outtxt)
                self.printLog('#OUT','Full REST output: %s' % outfile)
            return outtxt
        except:
            self.errorLog('%s.restFullOutput() error.' % self.prog(),warnlist=False)
            return self.restOutputError('%s.restFullOutput() error.' % self.prog())
#########################################################################################################################
    ### <5> ### Forks
#########################################################################################################################
    def _activeForks(self,pidlist=[]):   ### Checks Process IDs of list and returns list of those still running.
        '''
        Checks Process IDs of list and returns list of those still running.
        >> pidlist:list of integers = Process IDs
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            oldpids = pidlist[0:]
            ### ~ [1] Check PIDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pid in oldpids[0:]:
                (newpid,exitcode) = os.waitpid(pid,os.WNOHANG)
                if newpid == pid and exitcode in [0,256]: oldpids.remove(pid)
                elif newpid == pid:
                    oldpids.remove(pid)
                    self.errorLog('WARNING!: PID %d returned with exit code %d.' % (pid,exitcode),printerror=False)
            return oldpids
        except:
            self.log.errorLog('Error in _activeForks(%s)' % (pidlist))
            raise   
#########################################################################################################################
### End of RJE_Object Class                                                                                             #
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
###  RJE_Object Class:                                                                                                  # 
#########################################################################################################################
class RJE_ObjectLite(RJE_Object_Shell):     ### Metclass for inheritance by other classes
    '''Metaclass for inheritance by other classes.'''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
    ### => Now mostly handled in inherited RJE_Object_Shell Class                                                       #
#########################################################################################################################
    def _setGeneralAttributes(self):    ### Sets general attributes for use in all classes
        '''Sets general attributes for use in all classes.'''
        self.info = {'Name':'None'}
        self.stat = {}  # Remember to add 'Verbose' and 'Interactive' if needed to be different from parent
        self.opt = {}
        self.obj = {}
        self.list = {}
        self.dict = {}
#########################################################################################################################
    def prog(self): return self.log.info['Name']
    def name(self): return self.getStr('Name',default='%s' % self)
#########################################################################################################################
    def _setDefaults(self,info='None',opt=False,stat=0.0,obj=None,setlist=False,setdict=False):     ### Default defaults!
        '''
        Sets default defaults.
        >> info:str = default info setting
        >> opt:bool = default opt setting
        >> stat:float = default stat setting
        >> obj:Object = default obj setting
        >> setlist:boolean = whether to set list attributes [False]
        >> setdict:boolean = whether to set dictionary attributes [False]
        '''
        for i in self.infolist: self.info[i] = info
        for o in self.optlist: self.opt[o] = opt
        for s in self.statlist: self.stat[s] = stat
        for j in self.objlist: self.obj[j] = obj
        if setlist:
            for l in self.listlist: self.list[l] = []
        if setdict:
            for d in self.dictlist: self.dict[d] = {}
#########################################################################################################################
    def progLog(self, id='#ERR', text='Log Text Missing!',screen=True,rand=0.0,clear=0):
        if self.v() < 1: return
        if rand > 0 and random.random() > rand: return
        if 'ProgLog' in self.opt and not self.opt['ProgLog']: return False
        return self.printLog('\r%s' % id,text,screen=screen,log=False,newline=False,clear=clear)
    def printLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True,clear=0):
        return self.log.printLog(id,text,timeout,screen and not self.getAttribute('opt','Silent',False),log and not self.getAttribute('opt','Silent',False),newline,clear=clear)
    def errorLog(self, text='Missing text for errorLog() call!',quitchoice=False,printerror=True,nextline=True,log=True,errorlog=True,warnlist=True):
        #try:
        return self.log.errorLog(text,quitchoice,printerror,nextline,log,errorlog,warnlist)
        #except: return self.log.errorLog(text,quitchoice,False,nextline,log,errorlog)
    def warnLog(self,message,warntype=None,quitchoice=False,suppress=False,dev=False,screen=True):
        return self.log.warnLog(message,warntype,quitchoice,suppress,dev,screen)
#########################################################################################################################
    ### <2> ### Object Attributes                                                                                       #
#########################################################################################################################
    def details(self): return attDetails(printblanks=False)  ### Prints Details to screen
    def attDetails(self,types=['All'],printblanks=True):     ### Prints Details to screen
        '''Returns object details as text.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info.has_key('Type'): details = '%s (%s)\n' % (self.info['Name'], self.info['Type'])
            else: details = '%s\n' % self.info['Name']
            ### ~ [1] ~ Info ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Info' in types or 'All' in types:
                for info in sortKeys(self.info):
                    if info in ['Name','Type']: continue
                    if not self.info[info] and not printblanks: continue
                    details += '%s: %s\n' % (info,self.info[info])
            ### ~ [2] ~ Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Stat' in types or 'All' in types:
                if self.stat:
                    details += 'Stats: [ '
                    for stat in sortKeys(self.stat): details += '%s: %s; ' % (stat,self.stat[stat])
                    details += ']\n'
            ### ~ [3] ~ Options~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Opt' in types or 'All' in types:
                if self.opt:
                    details += 'Options: [ '
                    for opt in sortKeys(self.opt): details += '%s: %s; ' % (opt,self.opt[opt])
                    details += ']\n'
            ### ~ [4] ~ Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'List' in types or 'All' in types:
                for _list in sortKeys(self.list):
                    if not self.list[_list] and not printblanks: continue
                    details += '%s: %s\n' % (_list,self.list[_list])
            ### ~ [5] ~ Dicts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Dict' in types or 'All' in types:
                for _dict in sortKeys(self.dict):
                    if not self.dict[_dict] and not printblanks: continue
                    details += '%s: %s\n' % (_dict,self.dict[_dict])
            ### ~ [6] ~ Obj __~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if ('Obj' in types or 'All' in types) and self.obj:
                for object in sortKeys(self.obj):
                    try:
                        if not self.obj[object] and printblanks: details += '- %s: None\n' % object
                        else:
                            details += '- %s: %s (%s)\n' % (object, self.obj[object].info['Name'],self.obj[object]) #self.obj[object].details())
                            #!# Add option to add details? - Be careful of circular additions. Store objlist? #!#
                    except: self.log.errorLog('%s (%s) not an acceptable object.' % (self.obj[object],object),False,False)
            details += '\n'
            ### ~ [7] ~ Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return details
        except:
            self.errorLog('Major problem with attDetails()')
            return 'Major problem with attDetails(): %s\n' % sys.exc_info()[0]
#########################################################################################################################
    def edit(self):     ### Options to change details
        '''
        Options to change all object details (Info, Stat, Opt) and associated objects' details too.
        #!# Lists and Dictionaries not included. #!#
        '''
        try:
            print '\nEdit %s Attributes.\nEnter new values or leave Blank to retain.\n' % self.info['Name']
            for info in self.infolist: self.info[info] = self._editChoice(info,self.info[info])
            for stat in self.statlist: self.stat[stat] = self._editChoice(stat,self.stat[stat],numeric=True)
            for opt in self.optlist: self.opt[opt] = self._editChoice(opt,self.opt[opt],boolean=True)
            for object in self.objlist:
                if yesNo('Edit %s object?' % object): self.obj[object].edit()
            if self.info['Name'] != '' and yesNo('%s\nOK?: ' % self.details()): return
            self.edit()
        except: self.log.errorLog('Major problem with Object edit()')
#########################################################################################################################
    def _editChoice(self,text,value,numeric=False,boolean=False):  ### Returns Current or New Value as chosen
        '''
        Returns Current or New Value as chosen.
        >> text:str = Text to dislay for choice
        >> value:str/int/float/bool = Existing value
        >> numeric:bool = whether a numeric value is wanted
        >> boolean:bool = whether a True/False value is wanted
        '''
        try:
            if numeric: return getFloat(text=text,default=value)
            elif boolean: return getBool(text=text,default=value)
            else: return choice(text=text,default=value)
        except: self.log.errorLog('Major problem with _editChoice()')
#########################################################################################################################
    def getAtt(self,type,key,default=None): return self.getAttribute(type,key,default)
    def getAttribute(self,type,key,default=None):    ### Gets object information of correct type
        '''
        Gets object information of correct type.
        >> type:str = 'info','list','dict','opt','stat' or 'obj'
        >> key:str = key for given attribute dictionary
        >> default:anything = what to return if type is wrong or key not found in type
        '''
        ### Setup ###
        att = {}
        if type == 'info': att = self.info
        elif type == 'opt': att = self.opt
        elif type == 'stat': att = self.stat
        elif type == 'list': att = self.list
        elif type == 'dict': att = self.dict
        elif type == 'obj': att = self.obj
        elif type == 'data' and self.dict.has_key('Data'): att = self.dict['Data']
        ### Return ###
        if att.has_key(key): return att[key]
        elif 'Parent' in self.obj and self.obj['Parent']: return self.obj['Parent'].getAttribute(type,key,default)
        else: return default
#########################################################################################################################
    def getStr(self,ikey,default='None',checkdata=True): return self.getInfo(ikey,default,checkdata)
#########################################################################################################################
    def getInfo(self,ikey,default='None',checkdata=True):  ### Gets object information or returns default if not found
        '''
        Gets object information or returns default if not found.
        >> ikey:str = key for self.info dictionary
        >> default:str = values returned if ikey not found.
        >> checkdata:boolean = whether to check self.dict['Data'] if missing
        '''
        if self.info.has_key(ikey): return self.info[ikey]
        elif checkdata and self.dict.has_key('Data'):
            return getFromDict(self.dict['Data'],ikey,returnkey=False,case=False,default=default)
        return default
#########################################################################################################################
    def getOpt(self,okey,default=False):    ### Gets object opt or returns default
        '''Gets object opt or returns default.'''
        if self.opt.has_key(okey): return self.opt[okey]
        return default
#########################################################################################################################
    def getStat(self,skey,default=0,checkdata=True):  ### Gets object stat or returns default if not found
        '''
        Gets object stat or returns default if not found.
        >> skey:str = key for self.stat dictionary
        >> default:num = values returned if ikey not found.
        >> checkdata:boolean = whether to check self.dict['Data'] if missing
        '''
        if self.stat.has_key(skey): return self.stat[skey]
        elif checkdata and self.dict.has_key('Data'):
            val =  getFromDict(self.dict['Data'],skey,returnkey=False,case=False,default=default)
            if val != default: return string.atof(val)
        return default
#########################################################################################################################
    def getDict(self,dkey,dkeykey,default=None):    ### Returns value of self.dict[dkey][dkeykey] else default
        '''Returns value of self.dict[dkey][dkeykey] else default.'''
        if self.dict.has_key(dkey) and self.dict[dkey].has_key(dkeykey): return self.dict[dkey][dkeykey]
        return default
#########################################################################################################################
    def setAttribute(self,type,key,newvalue):    ### Sets object information of correct type from string
        '''
        Gets object information of correct type.
        >> type:str = 'info','list','dict','opt','stat' or 'obj'
        >> key:str = key for given attribute dictionary
        >> newvalue:str = string version of new value
        '''
        self._cmdRead('%s=%s' % (key.lower(),newvalue),type=type,att=key)
#########################################################################################################################
    def setInfo(self,infodic,addtolist=True):  ### Sets object information
        '''
        Sets object information.
        >> infodic:dictionary with keys corresponding to self.info keys
        >> addtolist:boolean = whether to add to self.infolist if missing
        '''
        try:
            for key in infodic.keys():
                self.info[key] = infodic[key]
                if addtolist and key not in self.infolist: self.infolist.append(key)
        except: self.log.errorLog('Problem with setInfo()',True)
#########################################################################################################################
    def setStat(self,statdic,addtolist=True):  ### Sets object Parameters
        '''
        Sets object information.
        >> statdic:dictionary with keys corresponding to self.stat keys
        >> addtolist:boolean = whether to add to self.statlist if missing
        '''
        try:
            for key in statdic.keys():
                self.stat[key] = statdic[key]
                if addtolist and key not in self.statlist: self.statlist.append(key)
        except: self.log.errorLog('Problem with setStat()',True)
#########################################################################################################################
    def setOpt(self,optdic,addtolist=True):  ### Sets object Options
        '''
        Sets object information.
        >> optdic:dictionary with keys corresponding to self.opt keys
        >> addtolist:boolean = whether to add to self.optlist if missing
        '''
        try:
            for key in optdic.keys():
                self.opt[key] = optdic[key]
                if addtolist and key not in self.optlist: self.optlist.append(key)
        except: self.log.errorLog('Problem with setOpt()',True)
#########################################################################################################################
    def setStr(self,attdic,addtolist=True): self.setInfo(attdic,addtolist)
    def setInt(self,attdic,addtolist=True): self.setStat(attdic,addtolist)
    def setNum(self,attdic,addtolist=True): self.setStat(attdic,addtolist)
    def setBool(self,attdic,addtolist=True): self.setOpt(attdic,addtolist)
#########################################################################################################################
    def setList(self,listdic,addtolist=True):  ### Sets object Lists
        '''
        Sets object list attributes.
        >> listdic:dictionary with keys corresponding to self.list keys
        >> addtolist:boolean = whether to add to self.listlist if missing
        '''
        try:
            for key in listdic.keys():
                self.list[key] = listdic[key]
                if addtolist and key not in self.listlist: self.listlist.append(key)
        except: self.log.errorLog('Problem with setList()',True)
#########################################################################################################################
    def setDict(self,dictdic,addtolist=True):  ### Sets object Dictionaries
        '''
        Sets object dictionary attributes.
        >> dictdic:dictionary with keys corresponding to self.dict keys
        >> addtolist:boolean = whether to add to self.dictlist if missing
        '''
        try:
            for key in dictdic.keys():
                self.dict[key] = dictdic[key]
                if addtolist and key not in self.dictlist: self.dictlist.append(key)
        except: self.log.errorLog('Problem with setDict()',True)
#########################################################################################################################
    def setObj(self,objdic,addtolist=True):  ### Sets object Objects
        '''
        Sets object information.
        >> objdic:dictionary with keys corresponding to self.obj keys
        >> addtolist:boolean = whether to add to self.objlist if missing
        '''
        try:
            for key in objdic.keys():
                self.obj[key] = objdic[key]
                if addtolist and key not in self.objlist: self.objlist.append(key)
        except: self.log.errorLog('Problem with setObj()',True)
#########################################################################################################################
    def setDictData(self,datadict,dictkey='Data'):  ### Sets object dict values
        '''
        Sets object Data dict values.
        >> datadict = Dictionary of values to add to self.dict[dictkey]
        '''
        try:
            if not self.dict.has_key(dictkey): self.dict[dictkey] = {}
            for key in datadict.keys(): self.dict[dictkey][key] = datadict[key]
        except: self.log.errorLog('Problem with setDictData()',True)
#########################################################################################################################
    def getData(self,dkey,dlist=['stat','info','opt'],case=False,str=True,default=None,dp=-1):  ### Gets data from dict['Data']
        '''
        Gets data from dict['Data'] if it exists and has key, else tries dlist dictionaries.
        >> dkey:str = Key for dictionaries
        >> dlist:list [self.stat,self.info,self.opt] = list of dictionaries to try after dict['Data']
        >> case:bool [False] = whether to match case for dkey
        >> str:bool [True] = whether to return all values as a string
        >> default [None] = what to return if no dictionary has key
        >> dp:int [-1] = Number of decimal places to use for stats if str=True (-1 = return as is)
        << returns value from data or stat/info/opt
        '''
        try:
            ### Setup ###
            dictlist = []
            if self.dict.has_key('Data'): dictlist = [self.dict['Data']]
            ddict = {'stat':self.stat,'info':self.info,'opt':self.opt}
            for dict in dlist:
                if ddict.has_key(dict): dictlist.append(ddict[dict])
                else: dictlist.append(dict)
            ### Look in dictionaries ###
            data = default
            for dict in dictlist:
                if data == default: data = getFromDict(dict,dkey,returnkey=False,case=case,default=default)
                if data != default and dict == self.stat:
                    if dp > 0: data = int(data * (10 ** dp) + 0.5) / float(10 ** dp)
                    elif dp == 0: data = int(data + 0.5)
            ### Return ###
            if str: return '%s' % data
            else: return data
        except: self.log.errorLog('Problem with %s.getData(%s)' % (self,dkey))
        return default
#########################################################################################################################
    def setLog(self,log,cascade=True):  ### Sets given log as log object and cascades through self.obj
        '''Sets given log as log object and cascades through self.obj.'''
        self.log = log
        for obj in self.obj.values():
            if cascade and obj and log != obj.log: obj.setLog(log)
#########################################################################################################################
### End of RJE_ObjectLite Class                                                                                         #
#########################################################################################################################


                                                    ### ~ ### ~ ###
#########################################################################################################################
###  Log Class: Error and Activity Logs                                                                                 #
#########################################################################################################################
class Log(RJE_Object_Shell):
    '''Class to handle log output: printing to log file and error reporting.'''
    # info['Name'] = Calling program
    # info['LogFile'] : log_file = None # General File Name
    # opt['NewLog'] : newlog = False  # whether to create new or just append 
    # stat['StartTime'] : start_time = 0  # Time at start of program
    # stat['HourMod'] : hourmod = 0     # Hour modifier to make time print right!
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes on basis of command-line parameters
#########################################################################################################################
    def __init__(self,itime=time.time(),cmd_list=[]):
        '''
        Handles log output; printing to log file and error reporting.
        >> itime:float = initiation time
        >> cmd_list:list of commandline variables
        '''
        try:### ~ [1] Setup Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.warnings = []              # Stores a list of warning types that have already been triggered.
            self.warnx = 0                  # Count of warning messages.
            self.no_suppression = []        # Stores a list of warning types that could be suppressed but should not be.
            self.cmd_list = cmd_list
            self._setGeneralAttributes()
            self.list['Log'] = []           # List of log messages. Stored for rest services when silent=True.
            self.list['WarnLog'] = []       # List of log warning messages. Stored when silent=True or memsaver=False.
            self.list['ErrorLog'] = []      # List of log error messages.
            self.stat['StartTime'] = itime
            self.info['LogFile'] = None
            self.opt['NewLog'] = False
            self.stat['HourMod'] = 0
            self.opt['Silent'] = False      # When True, will not write to screen or log
            self._cmdList()
            ### ~ [2] Log File details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['LogFile']:
                if self.info['LogFile'][-4:] != '.log': self.info['LogFile'] += '.log'
                if self.info['LogFile'][:1] == '/': self.info['LogFile'] = makePath(self.info['LogFile'],wholepath=True)
                elif self.info['LogFile'][1:2] == ':':
                    self.info['LogFile'] = makePath(self.info['LogFile'],wholepath=True)
                    if not self.opt['Win32']: self.printLog('#WIN32','Log path looks like Windows but Win32=F. Check!')
                else:
                    try:
                        mkDir(self,self.info['RunPath'])
                        self.info['LogFile'] = makePath('%s/%s' % (os.path.abspath(self.info['RunPath']),self.info['LogFile']),wholepath=True)
                    except: self.info['LogFile'] = makePath('%s/%s' % (os.path.abspath(os.curdir),self.info['LogFile']),wholepath=True)

                if '#DATE' in self.info['LogFile']: self.info['LogFile'] = string.replace(self.info['LogFile'],'#DATE',dateTime(dateonly=True))
                if self.opt['NewLog']:
                    self.verbose(0,2,'Make new file: %s' % self.info['LogFile'],2)
                    if checkForFile(self.info['LogFile']): os.unlink(self.info['LogFile'])
                else: self.verbose(0,2,'Append file: %s' % self.info['LogFile'],2)
                self.info['LogFile'] = os.path.abspath(self.info['LogFile'])
            ### ~ [3] ErrorLog File details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['ErrorLog'].lower() not in ['','none']:
                if self.opt['NewLog']:
                    self.verbose(0,2,'Make new error file: %s' % self.info['ErrorLog'],2)
                    if checkForFile(self.info['ErrorLog']): os.unlink(self.info['ErrorLog'])
                else: self.verbose(0,2,'Append error file: %s' % self.info['ErrorLog'],2)
                self.info['ErrorLog'] = os.path.abspath(self.info['ErrorLog'])
            ### ~ [4] Update files to full paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.info['Name'] = self.info['LogFile']
        except:
            self.errorLog('Problem with rje.Log (%s) initiation!' % self.info['LogFile'])
            raise
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Log
                self._cmdRead(cmd,type='info',att='LogFile',arg='log')
                self._cmdRead(cmd,type='info',att='LogFile',arg='logfile')
                self._cmdRead(cmd,type='info',att='LogFile',arg='basefile')
                self._cmdReadList(cmd,'opt',['NewLog','Silent'])
            except: self.errorLog('Problem with cmd: %s' % cmd)
        if self.opt['Silent']: self.stat['Verbose'] = -1
#########################################################################################################################
    ### <2> ### Log Output                                                                                              #
#########################################################################################################################
    def printLog(self,id='#ERR',text='',timeout=True,screen=True,log=True,newline=True,error=False,clear=0,warnlist=True):
        '''
        Prints text to log with or without run time.
        >> id:str = identifier for type of information
        >> text:str = log text
        >> timeout:boolean = whether to print run time
        >> screen:boolean = whether to print to screen (v>=0)
        >> log:boolean = whether to print to log file
        >> newline:boolean = whether to add newline if missing [True]
        >> error:bool [False] = whether to print to errorlog file
        >> clear:int [0] = Number of characters to clear on screen prior to printing
        >> warnlist:bool [True] = Whether to add the text to the stored list of warnings and errors.
        '''
        try:### ~ [1] ~ Setup text ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clear > 0: self.printLog('\r'+id,' ' * clear,timeout,screen,False,False)
            t = '#~+~~+~#'
            if timeout: t = self.myRunTime(time.time() - self.stat['StartTime'])
            text = string.join(['\r'+id,t,string.strip(text,'\r\n')],'\t')
            if text[-1:] <> '\n' and newline: text += '\n'
            ### ~ [2] ~ Storage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Webserver']: self.list['Log'].append(string.strip(text,'\r\n'))
            if warnlist and newline:
                warntext = string.strip(text,'\r\n')
                if id.endswith('#WARN') and (self.opt['Webserver'] or not self.opt['MemSaver']):
                    if warntext not in self.list['WarnLog'][-1:]: self.list['WarnLog'].append(warntext); self.warnx += 1
                if id.endswith('#ERR') and warntext not in self.list['ErrorLog'][-1:]: self.list['ErrorLog'].append(warntext)
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Silent'] and id != '#ERR': return text      # Do not write to log or screen
            if error and self.info['ErrorLog'].lower() not in ['','none']: logfile = self.info['ErrorLog']
            else: logfile = self.info['LogFile']
            if logfile and log: open(logfile,'a').write('%s\n' % string.strip(text,'\r\n'))
            if screen and error: print text
            elif screen and (log or self.opt['ProgLog']): self.verbose(0,4,text,0)
            return text
        except:
            if id == '#ERR': os._exit(1)
            self.errorLog('printLog() problem')
#########################################################################################################################
    def errorLog(self, text='Missing text for errorLog() call!',quitchoice=False,printerror=True,nextline=True,log=True,errorlog=True,warnlist=True):
        '''
        Raises text as error and prints to log.
        >> text:str = Error Description Text to print to log.
        >> quitchoice:bool = whether to give user choice to terminate program prematurely.
        >> printerror:bool = whether to print the system error. (Only if no return given.)
        >> nextline:bool [True] = whether to print error on next line
        '''
        ### Handle SystemExit and KeyboardInterrupt from previous sources
        try:
            new = {True:'\n',False:'\r'}[nextline]
            if not printerror or not sys.exc_info()[0]: return self.printLog('%s#ERR' % new,text,log=log,warnlist=warnlist)
            try: raise
            except SystemExit: os._exit(1)
            except KeyboardInterrupt: quitchoice = True
            except: sys.stderr = text
            if string.split(str(sys.exc_info()[0]),'.')[1] == 'SystemExit': os._exit(1)
        except: return self.printLog('%s#ERR' % new,text,log=log,warnlist=warnlist)

        try:
            ### Setup error variables
            error_type = str(sys.exc_info()[0])         # Error Type       : exceptions.IOError
            error_type = string.replace(error_type,'exceptions.','')
            error_value = str(sys.exc_info()[1])        # Error Value      : [Errno 9] Bad file descriptor
            error_traceback = traceback.extract_tb(sys.exc_info()[2])
            error_file = str(error_traceback[-1][0])    # File             : C:\Documents and Settings\normdavey\Desktop\Python\BLAST\Main.py
            error_method = str(error_traceback[-1][2])  # Method           : readFile
            error_line = str(error_traceback[-1][1])    # Line             : 15
            error_error = str(error_traceback[-1][3])   # Error            : for lines in fileIn.readlines():
            
            ### Log Error
            if text[-1:] <> '\n':
                if printerror:
                    if text[-1:] != ':' and text[-2:-1] != ':': text = '%s: ' % text
                    text = '%s%s (%s line %s) %s' % (text,error_type,error_method,error_line,error_value)
                    #text = '%s%s' % (text,error_type)
                text = '%s  \n' % text
            myreturn = self.printLog('%s#ERR' % new,text,log=log,warnlist=warnlist)
            if errorlog: self.printLog('\r#ERR',text,screen=not log,log=False,error=True,warnlist=warnlist)

            ### Quit Option and Additional Error Info
            logprog = self.info['Name']
            if logprog == 'None': logprog = 'Program'
            if (self.opt['DeBug'] or quitchoice) and (self.i() >= 0 and yesNo('Quit %s?' % logprog)):
                if error_type not in ['KeyboardInterrupt']:
                    self.verbose(-1,4,'%s: %s\nFile: %s\nMethod: %s (line %s)\nError: %s' % (error_type, error_value, error_file, error_method, error_line, error_error),2)
                os._exit(1)
            elif printerror and error_type not in ['KeyboardInterrupt']:
                self.verbose(1,1,'%s: %s\nFile: %s\nMethod: %s (line %s)\nError: %s' % (error_type, error_value, error_file, error_method, error_line, error_error),2)
                if self.dev(): self.printLog('#DEV','%s: %s;| File: %s;| Method: %s (line %s);| Error: %s' % (error_type, error_value, error_file, error_method, error_line, error_error),log=log)
        except SystemExit: os._exit(1)      #!# Why is this no longer working? #!#
        except KeyboardInterrupt: os._exit(0)
        except MemoryError: raise
        except:
            print 'That\'s not right: Error in errorLog()!!:',
            try:
                print '%s %s' % (sys.exc_info()[0],sys.exc_info()[1])
                print '%s: %s' % (error_type, error_value)
                print 'File: %s' % error_file
                print 'Method: %s (line %s)' % (error_method, error_line)
                print 'Error: %s' % error_error
            except: print 'This is really not right!! (Check Excel etc. isn\'t preventing writing to log)'
            os._exit(1)
        if quitchoice and self.i() < 0: raise
        return myreturn
#########################################################################################################################
    def warnLog(self,message,warntype=None,quitchoice=False,suppress=False,dev=False,screen=True,warnlist=True):   ### Prints warning text if appropriate
        '''
        Prints warning text if appropriate.
        >> message:str = Warning message to print
        >> warntype:str [None] = Warning type to check and/or store in self.warnings.
        >> quitchoice:str [False] = Whether to give quit option if warning given.
        >> suppress:bool [False] = Whether choice given to suppress warnings. Will show all otherwise if False.
        >> dev:bool [False] = Whether to give choice to not suppress regular warnings.
        '''
        try:
            if not self.warn() or warntype in self.warnings: return
            self.printLog('\r#WARN',message,clear=min(120,len(message)+5),screen=screen,warnlist=warnlist)
            if dev: self.deBug('>> Dev Warning! <<')
            logprog = self.info['Name']
            if logprog == 'None': logprog = 'Program'
            if quitchoice and self.i() >= 0 and warntype not in self.warnings + self.no_suppression and yesNo('Quit %s following warning?' % logprog):
                self.printLog('\r#EXIT','%s terminated by user following warning.' % logprog)
                os._exit(1)
            elif quitchoice and self.i() < 0 and warntype == 'fatal':
                self.printLog('\r#EXIT','%s terminated after warning.' % logprog)
                os._exit(1)
            if not warntype or warntype in self.no_suppression: return
            if suppress:
                if self.i() < 0 or not yesNo('Suppress future "%s" warnings?' % warntype): self.no_suppression.append(warntype)
                else:
                    self.warnings.append(warntype)
                    self.printLog('\r#WARN','Note: future "%s" warnings have been supressed' % warntype,warnlist=False)
        except: self.errorLog('Log.warnLog(%s) error' % warntype,quitchoice=quitchoice)
#########################################################################################################################
    def myRunTime(self, secs):  ### Converts time in seconds into time for easy comprehension
        '''
        Converts time in seconds into time for easy comprehension
        >> secs:int = time in seconds
        << str = hours:minutes:seconds
        '''
        (day,hour,min,sec) = (0,0,0,int(secs))
        while sec > 59: min += 1; sec -= 60
        while min > 59: hour += 1; min -= 60
        return '%s:%s:%s' % (preZero(hour,24),preZero(min,60),preZero(sec,60))
#########################################################################################################################
    def name(self): return self.info['Name']
#########################################################################################################################
    ### <3> ### Summarise/Display Warnings/ErrorLog?                                                                    #
#########################################################################################################################
    def endLog(self,proginfo=None): ### Method to run at end of program and summarise warnings/errors.
        '''
        Method to run at end of program and summarise warnings/errors.
        >> proginfo:Info object [None] = Will generate Log end statement if given.
        '''
        ### Summarise warnings and errors in log
        warnings = self.list['WarnLog'][0:]    # List of log warning messages. Stored when silent=True or memsaver=False.
        errors = self.list['ErrorLog'][0:]     # List of log error messages.
        if self.warnx: self.warnLog('%s warning messages: check log for details.' % iStr(self.warnx),warnlist=False)
        if errors: self.warnLog('%s error messages! Check log for details.' % iLen(errors),warnlist=False)
        ### End of Program Log entry
        if proginfo: self.printLog('#LOG', '%s V%s End: %s\n' % (proginfo.program,proginfo.version,time.asctime(time.localtime(time.time()))))
        ### Optional repeat of warnings and error messages
        if warnings and self.i() >= 0 and yesNo('Repeat %s warnings?' % iLen(warnings),default='N'): print string.join(warnings+[''],'\n')
        if errors and self.i() >= 0 and yesNo('Repeat %s error messages?' % iLen(errors)): print string.join(errors,'\n')
#########################################################################################################################
### End of Log                                                                                                          #
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
### Info Class: Basic Program Information                                                                               #
#########################################################################################################################
class Info(object):    ### Stores intro information for a program
    '''
    Stores intro information for a program.
    >> program:str = program name
    >> version:str = version number
    >> last_edit:str = last edit date
    >> description:str = program description
    >> author:str = author name
    >> start_time:float = starting time of program, time.time()
    '''
    program = None
    version = None
    last_edit = None
    description = None
    author = None
    start_time = None
    copyright = None
    comments = []
#########################################################################################################################
    def __init__(self,prog='Unknown',vers='X',edit='??/??/??',desc='Python script',author='Unknown',ptime=None,copyright='2007',comments=[]):
        '''
        Stores intro information for a program.
        >> prog:str = program name
        >> vers:str = version number
        >> edit:str = last edit date
        >> desc:str = program description
        >> author:str = author name
        >> ptime:float = starting time of program, time.time()
        >> copyright:str = year of copyright
        '''
        self.program = prog
        self.version = vers
        self.last_edit = edit
        self.description = desc
        self.author = author
        if ptime == None: self.start_time = time.time()
        else: self.start_time = ptime
        self.copyright = copyright
        self.comments = comments 
#########################################################################################################################
## End of Info
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
##  Out Class: verbosity & interactivity                                                                                #
#########################################################################################################################
class Out(RJE_Object_Shell):
    '''Class to handle basic generic output to screen based on Verbosity and Interactivity outside of Classes.''' 
#########################################################################################################################
    ### <2> ### General Class Methods
#########################################################################################################################
    def printIntro(self,info):  ### Prints introductory program information to screen.
        '''
        Prints introductory program information to screen.
        >> info:Info Object
        '''
        ## Setup lines ##
        self.verbose(0,3,'\n%s V%s run: %s' % (info.program,info.version,time.asctime(time.localtime(info.start_time))),2)
        line = ['  #############################################']
        line.append('  #|#> %s version %s : %s' % (info.program,info.version,info.description))
        line.append('  #|#> Copyright (c) %s %s <#~#> Last Modified: %s' % (info.copyright,info.author,info.last_edit))
        line.append('  #|#> Disclaimer: %s comes with ABSOLUTELY NO WARRANTY;' % info.program)
        line.append('  #|#>   This is free software, and you are welcome to redistribute it;')
        line.append('  #|#>   For details see attached license file (gnu_general_public_license.txt)')
        line.append('  #|#>   or http://www.gnu.org/copyleft/gpl.html.')
        for comment in info.comments: line.append('  #|#> %s' % comment)
        ## Find longest line ##
        maxlen = 0
        for li in line: maxlen = max(len(li),maxlen)
        ## Sort out spacer ##
        while len(line[0]) < maxlen: line[0] = line[0] + '#'
        line[0] = line[0] + '#####\n'
        ## Sort out end of lines
        while len(line[2]) < len(line[1]): line[2] = string.replace(line[2],'<#~','<#~~')
        for i in range(1,len(line)):
            while len(line[i]) < maxlen: line[i] = line[i] + ' '
            line[i] = line[i] + ' <#|#\n'
        ## Print out ##
        spacer = line[0:1]
        spacer = ['  #|#>' + '~' * (maxlen - 5) + '<#|#\n']
        if info.comments: introtext = string.join(line[0:2]+spacer+line[2:3]+spacer+line[3:7]+spacer+line[7:]+line[0:1],'')
        else: introtext = string.join(line[0:2]+spacer+line[2:3]+spacer+line[3:]+spacer+line[0:1],'')
        self.verbose(0,2,introtext,2)
        return introtext
#########################################################################################################################
    def SafeprintIntro(self,info):  ### Prints introductory program information to screen.
        '''
        Prints introductory program information to screen.
        >> info:Info Object
        '''
        self.verbose(0,3,"\n" + info.program + " run: " + time.asctime(time.localtime(info.start_time)),1)
        line=["  #############################################","",""]
        line[1]="  # #> " + info.program + " version " + info.version + " : " + info.description
        line[2]="  # #> Copyright (c) 2005 " + info.author + " :: Last Modified: " + info.last_edit
        if len(line[1])>len(line[2]): l=len(line[1])
        else: l=len(line[2])
        while len(line[0])<l: line[0]=line[0] + "#"
        line[0]=line[0]+"#####\n"
        while len(line[1])<l: line[1]=line[1] + " "
        while len(line[2])<l: line[2]=line[2] + " "
        self.verbose(0,4,"\n"+line[0]+line[1]+" <# #",1)
        self.verbose(0,2,line[0]+line[2]+" <# #\n"+line[0],1)

        gnu_license = '%s comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it;\n' % info.program
        gnu_license += 'For details see attached license file (gnu_general_public_license.txt) or http://www.gnu.org/copyleft/gpl.html.'
        self.verbose(0,1,gnu_license,2)
        self.verbose(0,4,line[0],2)
#########################################################################################################################
## End of Out
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
##  UserExit Exception Class: used to propogate user-selected program exit in controlled fashion                        #
#########################################################################################################################
class UserExit(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
#########################################################################################################################
## End of Out
#########################################################################################################################


                                                    ### ~ ### ~ ###


#########################################################################################################################
##  General Module Functions                                                                                           ##
#########################################################################################################################
if __name__ == "__main__":
    print "This module is not for standalone running!"
    os._exit(0)
#########################################################################################################################
def objtype(obj): ### Returns the type of obj as string
    typestr = '%s' % type(obj)
    return string.split(typestr,"'")[1]
#########################################################################################################################
def dateTime(t=(),yymmdd=False,dateonly=False): ### Returns date time string given time tuple t (makes if empty)
    '''Returns date time string given time tuple t (makes if empty).'''
    if not t: t = time.localtime(time.time())
    try:
        if yymmdd: return '%s%s%s' % (str(t[0])[-2:],preZero(t[1],12),preZero(t[2],31))
        if dateonly: return '%s-%s-%s' % (str(t[0]),preZero(t[1],12),preZero(t[2],31))
        return '%s-%s-%s %s:%s:%s' % (str(t[0]),preZero(t[1],12),preZero(t[2],31),preZero(t[3],24),preZero(t[4],60),preZero(t[5],60))
    except: return 'dateTime Error'
#########################################################################################################################
def regExp(re_object, text):     ### Returns matched groups
    '''
    Returns matched groups.
    >> re_object:re Object = re.compiled regular expression
    >> text:str = string to match
    '''
    m_pattern = re_object.search(text)
    if m_pattern: return m_pattern.groups()
    else: return None
#########################################################################################################################
def matchExp(re_pattern, text, flags=0):     ### Returns matched groups or None if no match.
    '''
    Returns matched groups or None if no match.
    >> re_pattern:str = regular expression
    >> text:str = string to match
    >> flags:int = regex flags, e.g. re.DOTALL (see https://docs.python.org/2/library/re.html)
    '''
    re_object = re.compile(re_pattern,flags)    #!# Should improve efficiency by optionally storing this?
    if re.search(re_pattern,text): return regExp(re_object,text)
    else: return None
#########################################################################################################################
def perCount(num,max,steps,prints):     ### Returns string to print (if any)
    '''
    Returns string to print (if any).
    >> num:int = current number
    >> max:int = max number
    >> steps:int = number of counts at which to print '.'
    >> prints:int = number of counts at which to print 'X.X%'
    '''
    text = ''
    if int(num) % int(prints) == 0: text += ' %.1f%% ' % (num * 100.0 / max)
    elif int(num) % int(steps) == 0: text += '.'
    return text
#########################################################################################################################
def perCounter(countlist):  ### Counter for functions
    '''
    Counter for functions. Adds one to 
    >> countlist = list of counter parameters:
    .. (sloop,ploop,total,max,steps,prints)       
        > sloop:int = current number of step loop
        > ploop:int = current number of print loop
        > total:int = total count
        > max:int = max number
        > steps:int = steps at which to print '.'
        > prints:int = steps at which to print 'X%'
        > verbosity:int = verbosity level of object
    Setup: perc = [0,0,0,max,steps,prints,verbosity], e.g. [0,0,0,x,100,1000,verbosity]
    .. use rje.setPerc(max,steps,prints,verbosity)
    << countlist
    '''
    ### <a> ### Increment
    if len(countlist) != 7: return  # Error!
    [sloop, ploop, total, max, steps, prints, verbosity] = countlist
    sloop += 1
    ploop += 1
    total += 1
    if sloop == steps: sloop = 0
    if ploop == prints:
        ploop = 0
        sloop = 0

    ### <b> ### Print if appropriate
    if verbosity >= 0 and (ploop == 0 or total == max):
        perc = float(total) / float(max)
        print '%d%%.' % (perc*100),
    elif verbosity >= 0 and sloop == 0: print '.',
    return [sloop, ploop, total, max, steps, prints, verbosity]
#########################################################################################################################
def setPerc(max,steps,prints,verbosity): ### Set up countlist for perCounter
    '''
    Set up countlist for perCounter.
    >> max:int = max number
    >> steps:int = steps at which to print '.'
    >> prints:int = steps at which to print 'X%'
    >> verbosity:int = verbosity level of object
    '''
    return [0,0,0,max,steps,prints,verbosity]
#########################################################################################################################
def progressPrint(callobj,x,dotx=1000,numx=10000,v=0):  ### Prints a dot or a number to screen dependent on numbers
    '''
    Prints a dot or a number to screen dependent on numbers.
    >> callobj:RJE_Object = object controlling verbosity
    >> x:int = current count
    >> dotx:int = count for a printed dot
    >> numx:int = count for a printed number
    >> v:int = verbosity level for output
    '''
    if x/numx == x/float(numx): callobj.verbose(v,4,integerString(x),0)
    elif x/dotx == x/float(dotx): callobj.verbose(v,4,'.',0)
#########################################################################################################################
def hhmmss(secs):  ### Converts time in seconds into time for easy comprehension
        '''
        Converts time in seconds into time for easy comprehension
        >> secs:int = time in seconds
        << str = hours:minutes:seconds
        '''
        while secs < 0: secs += 3600
        (day,hour,min,sec) = time.localtime(secs)[2:6]
        return '%s:%s:%s' % (preZero(24*(day-1)+hour,24*day),preZero(min,60),preZero(sec,60))
#########################################################################################################################
def memoryUse(who=None): ### Returns memory usage in kb
    '''
    Returns memory usage in kb. (Maybe b)
    Index	Field	Resource
    0	ru_utime	time in user mode (float)
    1	ru_stime	time in system mode (float)
    2	ru_maxrss	maximum resident set size
    3	ru_ixrss	shared memory size
    4	ru_idrss	unshared memory size
    5	ru_isrss	unshared stack size
    6	ru_minflt	page faults not requiring I/O
    7	ru_majflt	page faults requiring I/O
    8	ru_nswap	number of swap outs
    9	ru_inblock	block input operations
    10	ru_oublock	block output operations
    11	ru_msgsnd	messages sent
    12	ru_msgrcv	messages received
    13	ru_nsignals	signals received
    14	ru_nvcsw	voluntary context switches
    15	ru_nivcsw	involuntary context switches
    '''
    try:### ~ [0] This will only work on UNIX and MAC OSX. Look into Windows equivalent? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not who: who = resource.RUSAGE_SELF
        return resource.getrusage(who).ru_maxrss 
    except: return 0.0
#########################################################################################################################
### End of General Functions                                                                                            #
#########################################################################################################################


#########################################################################################################################
### String Functions                                                                                                    #
#########################################################################################################################
def chomp(text):
    if text: return text.replace('\n','').replace('\r','')
    else: return text
#########################################################################################################################
def preZero(num,max):   ### Adds leading zeros to number and returns as string
    '''
    Adds trailing zeros to integer values for output of equal length
    >> num:int = number to be extended
    >> max:int = maximum number in set (defines 'length' of output number)
    << prezero:str = string of num with leading zeros
    '''
    neg = num < 0
    if neg: num = -num
    prezero = []
    while (len(str(num))+len(prezero)) < len(str(max)): prezero.append('0')
    prezero.append(str(num))
    prezero = string.join(prezero,'')
    if neg: return '-%s' % prezero
    return prezero
#########################################################################################################################
def strSentence(instr,allwords=False,pure=False):     # Changes to sentence case
    '''
    Changes to sentence case. i.e. capitalisation of first letter.
    >> instr:str = string to convert
    >> allwords:bool [False] = whether to split on whitespace and convert all words.
    >> pure:bool [False] = whether to enforce lower case on rest of letters.
    '''
    if allwords:
        newstr = []
        for word in string.split(instr): newstr.append(strSentence(word,pure=pure))
        return string.join(newstr)
    if pure: instr[:1].upper() + instr[1:].lower()
    return instr[:1].upper() + instr[1:]
#########################################################################################################################
def strSub(instr,start,end,sub):   ### Replaces part of string and returns. Start and end are inclusive.
    '''
    Replaces part of string and returns. Start and end are inclusive.
    >> instr:str = full string
    >> start:int
    >> end:int
    >> sub:str = substitution
    '''
    end+=1
    newstr = instr[:start] + sub + instr[end:]
    return newstr
#########################################################################################################################
def strReplace(strtext,deltext,newtext='',case_sens=True,allocc=False,max=1): ### Replaces deltext with newtext in strtext
    '''
    Replaces first (or all) occurrences of deltext within strtext with newtext. Note that this does not use regular
    expressions and looks for exact matches only. Use re.sub() for regular expression replacements.
    >> strtext:str = text in which replacements to take place
    >> deltext:str = text to be replaced
    >> newtext:str = replacement text ['']
    >> case_sens:bool = whether to replace in case-sensitive manner [True]
    >> allocc:bool = whether to replace all occurrences [False]
    >> max:int = maximum number of replacement (if allocc=False) [1]
    << returns a tuple of (strtext,count), where strtext has replacements made and count is number of replacements.
    '''
    try:
        _stage = '<0> Setup'
        count = 0
        if deltext == newtext or newtext.find(deltext) >= 0:
            if deltext.find('#') >= 0:
                print 'Error with strReplace Setup. deltext in newtext - cannot replace like with like.'
                raise ValueError
            strtext = strReplace(strtext,deltext,'#',case_sens,allocc,max)[0]
            deltext = '#'

        _stage = '<1> Find and Replace'
        while allocc or count < max:
            if case_sens: findtext = strtext.find(deltext)
            else: findtext = strtext.lower().find(deltext.lower())
            if findtext >= 0:
                strtext = strtext[:findtext] + newtext + strtext[(findtext+len(deltext)):]
                count += 1
            else: break
    except: print 'Error with strReplace(%s): %s.' % (_stage,sys.exc_info()[0])
    return (strtext,count)        
#########################################################################################################################
def strEscape(text,charlist):   ### Adds escape \ characters in front of charlist
    '''Adds escape \ characters in front of charlist.'''
    for x in charlist: text = string.replace(text,x,'\\%s' % x)
    return text
#########################################################################################################################
def strReverse(text):   ### Returns reversed string
    '''Returns reversed string.'''
    return text[::-1]
#########################################################################################################################
def strSort(text,unique=False):   ### Returns sorted string
    '''Returns sorted string.'''
    letters = strList(text,unique)
    letters.sort()    
    return string.join(letters,sep='')
#########################################################################################################################
def strList(text,unique=False):   ### Returns string as list
    '''Returns string as list.'''
    if not unique: return string.split(string.join(text))
    letters = []
    for x in text:
        if not unique or x not in letters: letters.append(x)
    return letters
#########################################################################################################################
def strRearrange(text): ### Returns all possible orders of letters as list
    '''Returns all possible orders of letters as list.'''
    letters = string.split(string.join(text))
    bases = ['']; variants = []
    for i in range(len(text)):
        variants = []
        for base in bases:
            available = letters[0:]
            for b in range(len(base)): available.remove(base[b])
            for a in range(len(available)):
                newvar = base + available[a]
                if newvar not in variants: variants.append(newvar)
        bases = variants[0:]
    return variants
#########################################################################################################################
def iStr(number): return integerString(number)
#########################################################################################################################
def integerString(number):  ### Returns a string with commas for long integer, e.g. 1,000
    '''
    Returns a string with commas for long integer, e.g. 1,000.
    >> number:integer to be returned as string
    << intstring
    '''
    try: intstring = '%d' % int(number)
    except:
        try: intstring = str(int(number))
        except: intstring = str(number)
    intlist = []
    while len(intstring) > 3:
        intlist = [intstring[-3:]] + intlist
        intstring = intstring[:-3]
    intlist = [intstring] + intlist
    if intlist[0] == '-': return '-%s' % string.join(intlist[1:],',')
    return string.join(intlist,',')
#########################################################################################################################
def randomString(length,choices=''):   ### Returns a random string of given length
    '''
    Returns a random string of given length.
    >> length:int = length of string to return
    >> choices:str = choices of characters to randomly pick from (with replacement).
    '''
    if not choices: choices = 'ABCDEFGHIJKLMNOPQRSTUVWYXZabcdefghijklmnopqrstuvwxyz0123456789'
    rstring = ''
    for r in range(length): rstring += choices[-random.randint(1,len(choices))]
    return rstring
#########################################################################################################################
def shuffleString(instr):   ### Returns a randomly shuffled version of input string.
    '''Returns a randomly shuffled version of input string.'''
    rstring = strList(instr)
    random.shuffle(rstring)
    return string.join(rstring,'')
#########################################################################################################################
def stringStrip(instr,striplist):   ### Returns string with striplist strings removed
    '''Returns string with striplist strings removed.'''
    for x in striplist: instr = string.replace(instr,x,'')
    return instr
#########################################################################################################################
def fileSafeString(instr):  ### Returns a string that is safe for file names
    '''Returns a string that is safe for file names.'''
    for x in ['/','\\','?','%','*',':','|','"','>','<']: instr = string.replace(instr,x,'')
    return instr
#########################################################################################################################
### End of String Functions                                                                                             #
#########################################################################################################################


#########################################################################################################################
### User Input Functions                                                                                                #
#########################################################################################################################
def getInt(text='Integer Value?:',blank0=False,default='0',confirm=False):     ### Asks for a choice and returns integer
    '''
    Asks for a choice and returns integer.
    >> text:str = Prompt Text
    >> default:str = Default value as string
    '''
    try:
        if blank0: default = '0'
        c = choice(text,default,confirm)
        integer = string.atoi(c)
        return integer
    except KeyboardInterrupt: raise
    except:
        print 'Must be a number! (Integer) [Or blank for %s]' % default
        return getInt(text,blank0,default,confirm)
#########################################################################################################################
def getFloat(text='Numerical Value?:',default='0.0',confirm=False):     ### Asks for a choice and returns integer
    '''
    Asks for a choice and returns float.
    >> text:str = Prompt Text
    >> default:str = Default value as string
    '''
    try: return string.atof(choice(text,default,confirm))
    except KeyboardInterrupt: raise
    except:
        print 'Must be a number! (Float) [Or blank for %s]' % default
        return getFloat(text,default,confirm)
#########################################################################################################################
def getBool(text='Numerical Value?:',default=False,confirm=False):     ### Asks for a choice and returns boolean
    '''
    Asks for a choice and returns boolean.
    >> text:str = Prompt Text
    >> default:bool = Default value
    '''
    bool = choice(text,default,confirm)
    if bool.lower().find('t') == 0: return True
    elif bool.lower().find('f') == 0: return False
    else:
        print 'Must be True or False! [Or blank for %s]' % default
        return getBool(text,default,confirm)
#########################################################################################################################
def getFileName(text='File Name?',default='',mustexist=True,confirm=False): ### Asks for a filename
    '''
    Asks for a filename as part of interactive menus etc.
    >> text:str = Text to print as prompt for file name
    >> default:str = Default file name
    >> mustexist:boolean = Whether the file must exist (ask again if not)
    >> confirm:boolean = Whether the user should confirm the new filename
    << filename:str = file name entered by user
    '''
    while 1:
        filename = choice(text,default,confirm)
        if mustexist and not os.path.exists(filename): print 'File "%s" not found!' % filename
        else: return filename
#########################################################################################################################
def choice(text='?: ',default='',confirm=False): ### Asks for a choice and returns input
    '''
    Asks for a choice and returns input.
    >> default:str = default value given for blank entry ['']
    '''
    while text[-1] == ' ': text = text[:-1]
    if text[-1] == ':': text = text[:-1]
    if default: text = '%s [default=%s]' % (text,default)
    print '%s: ' % text,
    if 'pwin' in sys.argv: mychoice = raw_input('%s: ' % text)
    else: mychoice = raw_input()
    if mychoice == '': mychoice = '%s' % default
    if confirm and mychoice and yesNo('=> New value = "%s"?' % mychoice) == False: return choice(text,default,confirm)
    if 'pwin' in sys.argv: print mychoice   ### Restate if running in PythonWin
    if mychoice == '\\t': mychoice = '\t'
    return mychoice
#########################################################################################################################
def yesNo(text='',default='Y',confirm=False):    ### Asks for yes or no and returns True or False
    '''Asks for yes or no and returns True or False.'''
    try:
        answer = choice('%s (y/n)' % text,default,confirm).upper()[:1]
        if answer == 'Y': return True
        elif answer == 'N': return False
        else: return yesNo(text,default,confirm)
    except KeyboardInterrupt: raise
    except: return False
#########################################################################################################################
def errorMsg(): ### Prints error message to screen and asks for death
    '''Prints error message to screen and asks for death.'''
    ### Setup error variables
    error_type = str(sys.exc_info()[0])         # Error Type       : exceptions.IOError
    error_type = string.replace(error_type,'exceptions.','')
    error_value = str(sys.exc_info()[1])        # Error Value      : [Errno 9] Bad file descriptor
    error_traceback = traceback.extract_tb(sys.exc_info()[2])
    error_file = str(error_traceback[-1][0])    # File             : C:\Documents and Settings\normdavey\Desktop\Python\BLAST\Main.py
    error_method = str(error_traceback[-1][2])  # Method           : readFile
    error_line = str(error_traceback[-1][1])    # Line             : 15
    error_error = str(error_traceback[-1][3])   # Error            : for lines in fileIn.readlines():
    ### Print Error ###
    print '%s %s' % (sys.exc_info()[0],sys.exc_info()[1])
    print '%s: %s' % (error_type, error_value)
    print 'File: %s' % error_file
    print 'Method: %s (line %s)' % (error_method, error_line)
    print 'Error: %s' % error_error
    if yesNo('Kill me now?'): os.exit(1)
#########################################################################################################################
### End of User Input Functions                                                                                         #
#########################################################################################################################


#########################################################################################################################
### Maths Functions                                                                                                     #
#########################################################################################################################
def meanse(numlist):    ### Returns (mean,standard error) for a list of numbers
    '''Returns (mean,standard error) for a list of numbers.'''
    if not numlist: return (0.0,0.0)
    n = len(numlist)
    mean = float(sum(numlist)) / n
    variance = 0.0
    for num in numlist: variance += (num - mean) * (num - mean)
    variance /= n
    sd = math.sqrt(variance)
    se = sd / math.sqrt(n)
    return (mean,se)
#########################################################################################################################
def meansd(numlist):    ### Returns (mean,standard deviation) for a list of numbers
    '''Returns (mean,standard deviation) for a list of numbers.'''
    if not numlist: return (0.0,0.0)
    n = len(numlist)
    mean = float(sum(numlist)) / n
    variance = 0.0
    for num in numlist: variance += (num - mean) * (num - mean)
    variance /= n
    sd = math.sqrt(variance)
    return (mean,sd)
#########################################################################################################################
def modulus(num):   ### Returns modulus of number
    '''Returns modulus of number.'''
    if num < 0: return -num
    else: return num
#########################################################################################################################
def objFactorial(callobj,m): ### Returns the factorial of the number m
    '''Returns the factorial of the number m.'''
    if not callobj.list.has_key('Factorial'): callobj.list['Factorial'] = [1,1]
    if m < len(callobj.list['Factorial']): return callobj.list['Factorial'][m]
    while len(callobj.list['Factorial']) <= m:
        callobj.list['Factorial'].append(callobj.list['Factorial'][-1] * len(callobj.list['Factorial']))
    return callobj.list['Factorial'][m]
#########################################################################################################################
def factorial(m,callobj=None): ### Returns the factorial of the number m
    '''Returns the factorial of the number m.'''
    if callobj: return objFactorial(callobj,m)
    value = 1
    if m != 0:
        while m != 1:
            value = value*m
            m = m - 1
    return value
#########################################################################################################################
def isEven(num): return not isOdd(num)
def isOdd(num):     ### Returns True if Odd or False if Even
    '''Returns True if Odd or False if Even.'''
    if (float(num) / 2) == (int(num) / 2): return False
    return True
#########################################################################################################################
def geoMean(numlist=[]):    ### Returns geometric mean of numbers
    '''Returns geometric mean of numbers in list.'''
    geomean = 0.0
    for num in numlist: geomean += math.log(num)
    geomean /= len(numlist)
    geomean = math.exp(geomean)
    return geomean
#########################################################################################################################
def geoMeanSE(numlist=[]):    ### Returns geometric mean of numbers
    '''Returns geometric mean of numbers in list.'''
    geomean = []
    for num in numlist: geomean.append(math.log(num))
    (geomean,geose) = meanse(geomean)
    return (math.exp(geomean),math.exp(geose))
#########################################################################################################################
def formula(callobj=None,formula='',data={},varlist=[],operators=[],check=False,calculate=True):  ### Calculates formula
    '''
    Calculates formula using data dictionary, restricting to varlist and operators if desired. This calculation is
    executed in lower case, so varlist, data and formula need not match case. However, this obviously means that case-
    sensitive variables cannot be used.
    >> callobj:Object [None] = calling object, used for error messages if given.
    >> formula:str [''] = Formula as a string. Will be split on operators. Can have variables in data or numbers.
    >> data:dict {} = Dictionary of {variable:value}, where values should be numbers.
    >> varlist:list [] = List of data.keys() to include in calculation (in case some are strings etc.)
    >> operators:list [] = List of restricted operators. Will use ()+-*/^ if none given. (Other brackets replaced.)
    >> check:boolean [False] = Whether to check for any variables missing from data.keys()/varlist
    >> calculate:boolean [True] = Whether to calculate result of formula and return
    << value:float = results of calculation
    '''
    try:
        ### Setup operator list ###
        oplist = [')','(','+','-','*','/','^']
        for op in oplist[0:]:
            if operators and op not in operators:
                oplist.remove(op)
                if check and callobj: callobj.log.errorLog('Formula operator "%s" not recognised!' % op,printerror=False)
                if not calculate: return False

        ### Setup vardict ###
        vardict = {}
        for key in data.keys(): vardict[key.lower()] = data[key]
                    
        ### Setup varlist ###
        if varlist:
            varlist = string.split(string.join(varlist,',').lower(),',')
            if vardict:
                for var in varlist[0:]:
                    if not vardict.has_key(var): varlist.remove(var)
        else: varlist = vardict.keys()
                
        ### Setup formula list ###
        for left in ['[','{']: formula = string.replace(formula,left,'(')
        for right in [']','}']: formula = string.replace(formula,right,')')
        for op in oplist: formula = string.replace(formula,op,',%s,' % op).lower()
        formula = string.split('(,%s,)' % formula,',')
        while formula.count('') > 0: formula.remove('')
        if formula.count('(') != formula.count(')'):
            if callobj: callobj.log.errorLog('Formula brackets do not balance!',printerror=False)
            else: print 'Formula brackets do not balance!'
            if not calculate: return False
            else: return '!ERR!'
        for part in formula:
            if part not in oplist and part not in varlist:
                try: float(part)
                except:
                    if check and callobj: callobj.log.errorLog('Variable "%s" not recognised!' % part,printerror=False)
                    if not calculate: return False
                    else: return '!ERR!'
        if not calculate: return True  ## Assume check=True ... end of checking
            
        ### Process Formula ###
        while formula:
            ## Find single pair of brackets ##
            (x,y) = (0,0)
            while formula[x:(y+1)].count('(') != 1 or formula[x:(y+1)].count(')') != 1:
                while formula[x] == '(' and y <= x: y += 1
                while formula[y] not in [')','(']: y += 1
                if formula[y] == '(': x = y
            ## Identify stuff between brackets ##
            part = formula[(x+1):y]
            ## Peform calculation ##
            value = float(getFromDict(vardict,part.pop(0)))
            while part:
                op = part.pop(0)
                if op == '+': value = value + float(getFromDict(vardict,part.pop(0)))
                if op == '-': value = value - float(getFromDict(vardict,part.pop(0)))
                if op == '/':
                    try: value = value / float(getFromDict(vardict,part.pop(0)))
                    except ZeroDivisionError: value = 0.0
                    except: raise
                if op == '*': value = value * float(getFromDict(vardict,part.pop(0)))
                if op == '^': value = value ** float(getFromDict(vardict,part.pop(0)))
            if x == 0: return value
            else: formula = formula[:x] + [value] + formula[(y+1):]
        print '!ERR!'
        return value                    
    except:
        if callobj: callobj.log.errorLog('Major problem with rje.formula()')
        else: Log().errorLog('Major problem with rje.formula()')
        if calculate: return '!ERR!'
        else: return False
#########################################################################################################################
def poisson(observed,expected,exact=False,callobj=None,uselog=True): ### Returns the poisson probability of observed+ occurrences, given expected
    '''Returns the poisson probability of observed+ occurrences, given expected. (Or of exactly observed occurrences if exact=True.'''
    ### Exact ###
    if uselog:
        try: return logPoisson(observed,expected,exact,callobj)
        except: pass 
    expected = float(expected)
    if exact:
        try: return (math.exp(-expected) * pow(expected,observed) / factorial(observed,callobj))
        except:
            if callobj: callobj.errorLog('Poisson error: will return p = 0.0')
            return 0.0
    ### Cumulative ###
    prob = 1.0
    for x in range(0,observed):
        try:        #!# Fudge for OverflowError: long int too large to convert to float
            prob -= (math.exp(-expected) * pow(expected,x) / float(factorial(x,callobj)))
            #X#print x, prob
        except KeyboardInterrupt: raise
        except:
            if callobj: callobj.errorLog('Poisson error: will return p = %.3f' % (prob))
            #X#raise
            break
    if prob >= 0: return prob
    else: return 0.0
#########################################################################################################################
def logPoisson(observed,expected,exact=False,callobj=None): ### Returns the poisson probability of observed+ occurrences, given expected
    '''Returns the poisson probability of observed+ occurrences, given expected. (Or of exactly observed occurrences
    if exact=True.'''
    ### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if observed <= 0 and (expected == 0 or not exact): return 1.0
    if expected < 0:
        if callobj: callobj.errorLog('Warning: LogPoisson expected < 0',printerror=False)
        expected = 0
    expected = float(expected)
    if expected == 0: return 0.0  # Already handled observed=0 case
    ### ~ [1] Exact probability of observed given expected ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if exact:
        try: return math.exp((math.log(expected)*observed) - expected - logFactorial(observed,callobj))
        except:
            if callobj: callobj.errorLog('LogPoisson error: will return p = 0.0')
            return 0.0
    ### ~ [2] Cumulative probability of observed+ given expected ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    prob = 1.0
    for x in range(0,observed):
        try: prob -= math.exp((math.log(expected)*x) - expected - logFactorial(x,callobj))
        except KeyboardInterrupt: raise
        except:
            if callobj: callobj.errorLog('LogPoisson error: will return p = %.3f' % (prob))
            break
    if prob >= 0: return prob
    else: return 0.0
#########################################################################################################################
def logFactorial(m,callobj=None): ### Returns the log factorial of the number m
    '''Returns the factorial of the number m.'''
    try: return callobj.dict['LogFactorial'][m]
    except:
        if callobj and not callobj.dict.has_key('LogFactorial'): callobj.dict['LogFactorial'] = {}
    if callobj: callobj.dict['LogFactorial'][0] = callobj.dict['LogFactorial'][1] = 0.0
    try:
        x = sortKeys(callobj.dict['LogFactorial'])[-1]
        value = callobj.dict['LogFactorial'][x]
        x += 1
    except: (value,x) = (0.0,1)
    for i in range(x,m+1):
        value += math.log(i)
        if callobj: callobj.dict['LogFactorial'][i] = value
    return value
#########################################################################################################################
def logBinomial(observed,trials,prob,exact=False,callobj=None): ### Returns the binomial probability of observed+ occurrences
    '''
    Returns the binomial probability of observed+ occurrences.
    >> observed:int = number of successes (k)
    >> trials:int = number of trials (n)
    >> prob:int = probability of success of each trial
    >> usepoisson:bool = whether to use Poisson if Binomial fails [True]
    '''
    try:
        if getFromDict(callobj.stat,'MaxBin',returnkey=False,case=True,default=trials) < trials: return logPoisson(observed,trials*prob,exact,callobj)
    except: pass
    ### Exact Binomial Probability ###
    if exact:
        #x#print observed, trials, prob,
        logp = logFactorial(trials,callobj) - (logFactorial(observed,callobj) + logFactorial(trials - observed,callobj))    #+ math.log(pow(prob,observed)) + math.log(pow(1-prob,trials-observed))
        if logp > 0: p = math.exp(logp) * pow(prob,observed) * pow(1-prob,trials-observed)
        else: p = pow(prob,observed) * pow(1-prob,trials-observed)
        #x#print p
        return p
    ### Cumulative Probability ###
    if observed == 0: return 1.0
    elif observed == 1:
        p = 1 - pow((1 - prob),trials)
        if p > 0.0: return p
    p = 0.0
    (observed, trials) = (int(observed),int(trials+0.99))
    for k in range(observed,trials+1):
        bp = logBinomial(k,trials,prob,exact=True,callobj=callobj)
        p += bp
        if p > 0 and bp == 0: break
    if p > 1: return 1.0
    return p
#########################################################################################################################
def binomial(observed,trials,prob,exact=False,usepoisson=True,callobj=None): ### Returns the binomial probability of observed+ occurrences
    '''
    Returns the binomial probability of observed+ occurrences.
    >> observed:int = number of successes (k)
    >> trials:int = number of trials (n)
    >> prob:int = probability of success of each trial
    >> usepoisson:bool = whether to use Poisson if Binomial fails [True]
    '''
    if not prob or not trials:
        if observed: return 0.0
        else: return 1.0
    try: return logBinomial(observed,trials,prob,exact,callobj)
    except: pass #X#if callobj: callobj.log.errorLog('Cannot use logBinomial! Will try old binomial with poisson backup.')
    try:
        if getFromDict(callobj.stat,'MaxBin',returnkey=False,case=True,default=trials) < trials: return poisson(observed,trials*prob,exact,callobj)
    except: pass
    ### Exact Binomial Probability ###
    if exact:
        try:
            return (float(factorial(trials,callobj))/(factorial(observed,callobj)*factorial(trials-observed,callobj))) * pow(prob,observed) * pow(1-prob,trials-observed)
        except:
            if usepoisson: return poisson(observed,trials*prob,exact,callobj)
            raise
    ### Cumulative Probability ###
    if observed == 0: return 1.0
    elif observed == 1:
        p = 1 - pow((1 - prob),trials)
        if p > 0.0: return p
    p = 0.0
    (observed, trials) = (int(observed),int(trials+0.99))
    for k in range(observed,trials+1):
        try:        
            p += (1.0 * factorial(trials,callobj)/(factorial(k,callobj)*factorial(trials-k,callobj))) * pow(prob,k) * pow(1-prob,trials-k)
        except:
            if usepoisson or p == 0.0: return poisson(observed,trials*prob,callobj=callobj)
            break
    if p > 1: return 1.0
    return p
#########################################################################################################################
def OLDbinomial(observed,trials,prob,exact=False,usepoisson=True,callobj=None): ### Returns the binomial probability of observed+ occurrences
    '''
    Returns the binomial probability of observed+ occurrences.
    >> observed:int = number of successes (k)
    >> trials:int = number of trials (n)
    >> prob:int = probability of success of each trial
    >> usepoisson:bool = whether to use Poisson if Binomial fails [True]
    '''
    #x#print exact
    if exact:
        try:
            return (float(factorial(trials,callobj))/(factorial(observed,callobj)*factorial(trials-observed,callobj))) * pow(prob,observed) * pow(1-prob,trials-observed)
        except:
            if usepoisson:
                return poisson(observed,trials*prob,exact,callobj)
            raise
    p = 1.0
    for k in range(0,observed):
        try:        
            p -= (1.0 * factorial(trials,callobj)/(factorial(k,callobj)*factorial(trials-k,callobj))) * pow(prob,k) * pow(1-prob,trials-k)
            #X#print (1.0 * factorial(trials,callobj)/(factorial(k,callobj)*factorial(trials-k,callobj))) * pow(prob,k) * pow(1-prob,trials-k)
            #x#print k, p
        except KeyboardInterrupt:
            raise
        except:
            if usepoisson:
                return poisson(observed,trials*prob,callobj=callobj)
            raise
    if p < 0 and usepoisson:
        return poisson(observed,trials*prob,callobj=callobj)
    return p
#########################################################################################################################
def dp(data,dp): ### Returns number rounded to X dp
    '''Returns number rounded to X dp.'''
    #if dp == 1: data = string.atof('%.1f' % data)
    if dp > 0: data = int(data * (10 ** dp) + 0.5) / float(10 ** dp)
    elif dp == 0: data = int(data + 0.5)
    return data
#########################################################################################################################
def sf(data,sf=3): ### Returns number rounded to X sf
    '''Returns number rounded to X sf.'''
    neg = data < 0
    if neg: data = -data
    if not data: return dp(data,sf-1)
    if data > 1:    # Make smaller, round, increase
        x = 0
        while data >= 1:
            data /= 10.0
            x += 1
        data = dp(data,sf)
        for i in range(x): data *= 10.0
    else:   # Make bigger, round, decrease
        x = 0
        while (data*10) < 1:
            data *= 10.0
            x += 1
        data = dp(data,sf)
        for i in range(x): data /= 10.0
    if neg: data = -data
    if int(data) == data: data = int(data)
    return data
#########################################################################################################################
def eStr(_expect,strict=True): return expectString(_expect,strict)
def expectString(_expect,strict=True):  ### Returns formatted string for _expect value
    '''Returns formatted string for _expect value.'''
    try:
        _expect = float(_expect)
        if _expect < 0.0: return '-%s' % expectString(-_expect,strict)
        if not _expect: return '0.000'
        if _expect >= 10: return '%.1f' % _expect
        elif _expect >= 1: return '%.2f' % _expect
        elif _expect >= 0.1: return '%.3f' % _expect
        elif _expect >= 0.001: return '%.4f' % _expect
        else: return '%.2e' % _expect
    except:
        if not strict: return '%s' % _expect
        print _expect
        raise
#########################################################################################################################
def ratio(num,denom,dividezero=0.0):   ### Returns num/denom unless denom=0.0 - returns dividezero unless dividezero=None
    '''Returns num/denom unless denom=0.0 - returns dividezero unless dividezero=None.'''
    try: return num / float(denom)
    except:
        if dividezero == None: raise
        return dividezero
#########################################################################################################################
### End of Maths Functions                                                                                              #
#########################################################################################################################


#########################################################################################################################
### Dictionary Functions                                                                                                #
#########################################################################################################################
def sortKeys(dic,revsort=False):  ### Returns sorted keys of dictionary as list
    '''
    Returns sorted keys of dictionary as list.
    >> dic:dictionary object
    >> revsort:boolean = whether to reverse list before returning
    '''
    dkeys = dic.keys()
    dkeys.sort()
    if revsort: dkeys.reverse()
    return dkeys
#########################################################################################################################
def dictValues(dict,key,valtype='list'):   ### Returns dict values or empty list if dict does not have key
    '''Returns dict values or empty list if dict does not have key.'''
    if dict.has_key(key): return dict[key]
    if valtype == 'list': return []
    elif valtype == 'dict': return {}
    elif valtype == 'num': return 0
    elif valtype == 'str': return ''
    return None
#########################################################################################################################
def scaledict(dict={},scale=1.0):   ### Scales all values by scale and returns new dictionary
    '''Scales all values by scale and returns new dictionary.'''
    newdict = {}
    for key in dict.keys(): newdict[key] = dict[key] * scale
    return newdict
#########################################################################################################################
def getFromDict(dict,key,returnkey=True,case=True,default=None):   ### Returns dict value if it has key, else default/key 
    '''
    Returns dict value if it has key, else given key (or None).
    >> dict:dictionary from which to get value
    >> key:dictionary key
    >> returnkey:bool [True] = whether to return the key itself if not found in dictionary
    >> case:bool [True] = whether to use case-sensitive key matching
    >> default [None] = what to return if no entry
    '''
    if dict.has_key(key): return dict[key]
    elif not case:
        for dkey in dict.keys():
            if dkey.lower() == key.lower(): return dict[dkey]
    if returnkey: return key
    return default
#########################################################################################################################
def dictFreq(dict,total=True,newdict=False):  ### Normalises values of dict by total. Adds 'Total' key if desired. Must be numeric!
    '''Normalises values of dict by total. Adds 'Total' key if desired (total=True). Must be numeric values!'''
    dsum = 0.0
    if total and dict.has_key('Total'): dsum = float(dict.pop('Total'))
    if not dsum: dsum = float(sum(dict.values()))
    if newdict: fdict = {}
    else: fdict = dict
    for dkey in dict.keys():
        if dsum: fdict[dkey] = dict[dkey] / dsum
        else: fdict[dkey] = False
    if total: fdict['Total'] = dsum
    return fdict
#########################################################################################################################
def entropyDict(data,ikeys=[],fillblanks=True):    ### Calculate entropy of input dictionary
    '''
    Calculate entropy of input dictionary.
    >> data:dict = input data dictionary {key:numeric}
    >> ikeys:list = keys for entropy calculation. Uses data.keys() if []
    >> fillblanks:bool [True] = whether to fill in zero values for missing keys
    '''
    if not ikeys: ikeys = sortKeys(data)
    n = len(ikeys)
    if n == 1: return 1.0
    freq = {}
    for k in ikeys:
        try: freq[k] = data[k]
        except:
            if fillblanks: freq[k] = 0
            else: raise
    dictFreq(freq,total=False)
    entropy = 1.0
    for k in ikeys:
        if freq[k]: entropy += (freq[k] * math.log(freq[k],n))
    return max(0.0,entropy)
#########################################################################################################################
def combineDict(targetdict,sourcedict,overwrite=True,replaceblanks=True,copyblanks=False):  ### Adds data from sourcedict to targetdict (targetdict changes)
    '''
    Adds data from sourcedict to targetdict (targetdict changes).
    >> targetdict:dictionary that will be altered
    >> sourcedict:dictionary containing data to add to targetdict
    >> overwrite:bool [True] = whether keys from sourcedict will overwrite same data in targetdict
    >> replaceblanks:bool [True] = whether to replace existing but empty targetdict entries
    >> copyblanks:bool [False] = whether to copy blank entries over existing target entries
    '''
    for key in sourcedict:
        if key not in targetdict or (replaceblanks and not targetdict[key]): targetdict[key] = sourcedict[key]
        elif overwrite and (copyblanks or sourcedict[key]): targetdict[key] = sourcedict[key]
    return targetdict
#########################################################################################################################
def rankDict(data,rev=False,absolute=False,lowest=False):   ### Returns rank of values as new dictionary
    '''
    Returns rank of values as new dictionary.
    >> data:dict = input dictionary of scores.
    >> rev:Boolean = if True will return 0 for Highest
    >> absolute:boolean [False] = return 1 to n, rather than 0 to 1
    >> lowest:boolean [False] = returns lowest rank rather mean rank in case of ties
    << ranklist:list of ranks (0 = Lowest, 1 = Highest)
    '''
    ranked = {}
    scorelist = []
    for rkey in sortKeys(data): scorelist.append(data[rkey])
    scorelist = rankList(scorelist,rev,absolute,lowest)
    for rkey in sortKeys(data): ranked[rkey] = scorelist.pop(0)
    return ranked
#########################################################################################################################
def valueSortedKeys(data,rev=False):    ### Returns list of keys, sorted by values
    '''Returns list of keys, sorted by values.'''
    sortdict = combineDict({},data)
    sorter = data.values()
    sorter.sort()
    if rev: sorter.reverse()
    valsorted = []
    for val in sorter:
        for d in sortdict.keys()[0:]:
            if data[d] == val: valsorted.append(d); sortdict.pop(d)
    return valsorted
#########################################################################################################################
def dictKeysSortedByValues(dict,revsort=False):     ### Returns dictionary keys sorted according to corresponding values
    '''Returns dictionary keys sorted according to corresponding values.'''
    sortdict = {}
    for key in sortKeys(dict):
        if dict[key] not in sortdict: sortdict[dict[key]] = []
        sortdict[dict[key]].append(key)
    sortedkeys = []
    for v in sortKeys(sortdict,revsort): sortedkeys += sortdict[v]
    return sortedkeys
#########################################################################################################################
def dictValueList(dict,keys):   ### Returns dictionary (subset) values in same order as keys
    '''Returns dictionary (subset) values in same order as keys.'''
    vals = []
    for key in keys: vals.append(dict[key])
    return vals
#########################################################################################################################
### End of Dictionary Functions                                                                                         #
#########################################################################################################################


#########################################################################################################################
###  List Manipulation Functions                                                                                        #
#########################################################################################################################
def collapseTupleList(inlist):  ### Collapses a list of (x,y) values by merging/removing overlaps where required.
    '''Collapses a list of (x,y) values by merging/removing overlaps where required.'''
    inlist = inlist[0:] # Create a copy
    inlist.sort()
    i = 1
    while i < len(inlist):
        if inlist[i][0] <= inlist[i-1][1]+1:
            inlist[i-1] = (inlist[i-1][0],max(inlist[i-1][1],inlist[i][1]))
            inlist.pop(i)
        else: i += 1
    return inlist
#########################################################################################################################
def intList(inlist): ### Converts inlist to integers and returns
    '''Converts inlist to integers and returns.'''
    outlist = []
    for x in inlist: outlist.append(int(x))
    return outlist
#########################################################################################################################
def numList(inlist): ### Converts inlist to floats and returns
    '''Converts inlist to integers and returns.'''
    outlist = []
    for x in inlist: outlist.append(float(x))
    return outlist
#########################################################################################################################
def iLen(inlist): return integerString(len(inlist))
#########################################################################################################################
def asList(inlist): ### Checks whether inlist is list or string and returns as list
    tlist = inlist[0:]
    try: tlist.sort(); return inlist[0:]
    except: return [inlist]
#########################################################################################################################
def rankList(scorelist=[],rev=False,absolute=False,lowest=False,unique=False):  ### Returns rank of scores as list
    '''
    Returns rank of scores as list.
    >> scorelist = list of scores
    >> rev:Boolean = if True will return 0 for Highest
    >> absolute:boolean [False] = return 1 to n, rather than 0 to 1
    >> lowest:boolean [False] = returns lowest rank rather mean rank in case of ties
    >> unique:boolean [False] = give each element a unique rank (ties rank in order of entry)
    << ranklist:list of ranks (0 = Lowest, 1 = Highest)
    '''
    rsorted = scorelist[0:]
    rsorted.sort()
    if rev: rsorted.reverse()
    slen = len(rsorted) - 1
    ranklist = []
    for score in scorelist:
        if lowest: rank = rsorted.index(score)
        else: rank = (rsorted.index(score) * 2.0 + rsorted.count(score) - 1) / 2.0
        if absolute: ranklist.append(rank+1)
        else: ranklist.append(rank/slen)
        if unique: rsorted[rsorted.index(score)] = None
    return ranklist
#########################################################################################################################
def randomList(inlist): ### Returns inlist in randomised order
    '''
    Returns inlist in randomised order.
    >> inlist:List object
    << ranlist:randomised order list
    '''
    ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ordlist = inlist[0:]    # Don't mess up original
    ranlist = []            # Random list
    ### ~ [2] ~ Randomise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    while ordlist:
        r = random.randint(0,len(ordlist)-1)
        ranlist.append(ordlist.pop(r))
    if len(ranlist) == len(inlist): return ranlist
    else: raise ValueError
#########################################################################################################################
def sortUnique(inlist,xreplace=True,num=False): ### Returns sorted unique list: *** Case-insensitive collation ***
    '''
    Returns sorted unique list. Uses a dictionary to allow modification of sorting whilst keep list contents unchanged.
    >> inlist:List to be sorted
    >> xreplace:boolean [True] = whether to replace dots with xs for determining sort
    >> num:boolean [False] = whether the sortlist is numbers rather than strings
    '''
    ### ~ [1] Try, assuming that xreplace and num settings are right ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        if num: raise ValueError    # Trigger simple sorting
        repdict = {}
        for i in inlist:
            j = string.replace(i,'_','')     # Underscores are invisible to UNIX sort
            if xreplace:  # Special sort with . replaced by x
                repdict[string.replace(j,'.','x').upper()] = i
            else: repdict[j.upper()] = i
    except:
        repdict = {}
        for i in inlist: repdict[i] = i
    ## ~ [2] Return unique sorted list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    outlist = []
    for i in sortKeys(repdict):
        o = repdict[i]
        if o not in outlist[-1:]: outlist.append(o)
    return outlist
#########################################################################################################################
def listCombos(inlist,maxrep=0): ### Returns a list of all possible combinations of inlist entries
    '''
    Returns a list of all possible combinations of inlist entries.
    E.g. [AB,D,EF] will return [ [A,D,E], [A,D,F], [B,D,E], [B,D,F] ]
    >> maxrep:int [0] = Max. no. repetitions of an individual element
    '''
    bases = [[]]
    for el in inlist:
        variants = []
        for base in bases:
            for var in el:
                if base.count(el) >= maxrep > 0: continue
                variants.append(base + [var])
        bases = variants[0:]
    return bases
#########################################################################################################################
def listRearrange(inlist): ### Returns all possible orders of inlist
    '''Returns all possible orders of letters as list.'''
    bases = [[]]; variants = []
    for i in range(len(inlist)):
        variants = []
        for base in bases:
            available = inlist[0:]
            for b in base: available.remove(b)
            for a in available:
                newvar = base + [a]
                if newvar not in variants: variants.append(newvar)
        bases = variants[0:]
    return variants
#########################################################################################################################
#def listIntersect(list1,list2): return list(sets.Set(list1).intersection(sets.Set(list2)))
def listIntersect(list1,list2): return list(set(list1).intersection(set(list2)))
#########################################################################################################################
#def listUnion(list1,list2): return list(sets.Set(list1).union(sets.Set(list2)))
def listUnion(list1,list2): return list(set(list1).union(set(list2)))
#########################################################################################################################
#def listDifference(list1,list2): return list(sets.Set(list1).difference(sets.Set(list2)))
def listDifference(list1,list2):    ### Returns the elements of list1 that are not found in list 2.
    return list(set(list1).difference(set(list2)))
#########################################################################################################################
def listMax(inlist,jointmax=True):    # returns the list element(s) with max counts
    listmax = inlist[:1]
    for el in inlist[1:]:
        if inlist.count(el) > inlist.count(listmax[0]): listmax = [el]
        elif inlist.count(el) == inlist.count(listmax[0]): listmax.append(el)
    if not jointmax and len(listmax) > 1: listmax = []
    return listmax
#########################################################################################################################
def sortListsByLen(listoflists,rev=False):
    sortdict = {}
    for slist in listoflists:
        L = len(slist)
        if L not in sortdict: sortdict[L] = []
        sortdict[L].append(slist)
    sortlist = []
    for L in sortKeys(sortdict,rev): sortlist += sortdict[L]
    return sortlist
#########################################################################################################################
def listLower(inlist,convert=True,strict=False): ### Returns lower case version of list
    '''
    Returns lower case version of list.
    >> list = List object to convert
    >> convert:bool [True] = whether to convert non-strings into strings first ('%s' % s)
    >> strict:bool [False] = whether to raise error if convert=False and non-string found
    '''
    newlist = []
    for s in inlist:
        if convert: s = '%s' % s
        try: newlist.append(s.lower())
        except:
            if strict: raise
            newlist.append(s)
    return newlist
#########################################################################################################################
def listUpper(inlist,convert=True,strict=False): ### Returns upper case version of list
    '''
    Returns upper case version of list.
    >> list = List object to convert
    >> convert:bool [True] = whether to convert non-strings into strings first ('%s' % s)
    >> strict:bool [False] = whether to raise error if convert=False and non-string found
    '''
    newlist = []
    for s in inlist:
        if convert: s = '%s' % s
        try: newlist.append(s.upper())
        except:
            if strict: raise
            newlist.append(s)
    return newlist
#########################################################################################################################
def listJoin(inlist,sortunique=False): # Joins a list of lists into a single list
    '''Joins a list of lists into a single list.'''
    newlist = []
    for partlist in inlist: newlist += partlist
    if sortunique: return sortUnique(newlist)
    return newlist
#########################################################################################################################
def listToDict(inlist,keylist): ### Returns list as a dictionary, using keylist as keys
    '''Returns list as a dictionary, using keylist as keys.'''
    if not inlist: return {}
    if len(inlist) != len(keylist): raise ValueError('Mismatch of listToDict inlist and keylist lengths.')
    listdict = {}
    inlist = inlist[0:]; keylist = keylist[0:]
    while inlist: listdict[keylist.pop(0)] = inlist.pop(0)
    return listdict
#########################################################################################################################
### End of List Functions                                                                                               #
#########################################################################################################################


#########################################################################################################################
###  File Manipulation Functions                                                                                        #
#########################################################################################################################
def isYounger(file1,file2):     ### Returns younger file or None if either does not exist
    '''
    Compares age of files. file2 should be desired new file. file1 is returned in the case of a tie.
    Returns:
    - younger file, or
    - None if either does not exist, or
    - First file if of same age
    '''
    if os.access(file1, os.F_OK) == False or os.access(file2, os.F_OK) == False: return None
    birth = (os.stat(file1)[8],os.stat(file2)[8])
    if birth[0] >= birth[1]:     # file 1 was modified after file 2: file1 is younger
        return file1
    elif birth[1] > birth[0]:   # file 2 was modified after file 1: file2 is younger
        return file2
    else: return None
#########################################################################################################################
def exists(file,none_response=False):     ### Returns True if file exists or False if not
    '''Returns True if file exists or False if not.'''
    if not file or file.lower() in ['','none']: return none_response
    return os.path.exists(file)
#########################################################################################################################
def checkForFile(file):     ### Returns True if file exists or False if not 
    '''Returns True if file exists or False if not.'''
    return file and os.path.exists(file)
#########################################################################################################################
def makePath(path='',wholepath=False,return_blank=True):  ### Returns path that can be used for calling programs etc.
    '''
    Returns path that can be used for calling programs etc.
    >> path:str = Given path with directory separators as '/'
    >> wholepath:boolean = whether path includes the program call (True) or not (False)
    >> return_blank:boolean [True] = whether to return '' if given or replace with '.'
    << os_path:str = Returned path with appropriate separators
    '''
    if path == '':
        if return_blank: return ''
        else: path = '.'
    os_path = path.split('/')
    if os_path[-1] != '' and wholepath == False: os_path.append('')
    if wholepath and os_path[-1] == '': os_path = os_path[:-1]
    return string.join(os_path,os.sep)
#########################################################################################################################
def fileTransfer(fromfile=None,tofile=None,deletefrom=True,append=True):    ### Appends fromfile to tofile and deletes fromfile
    '''
    Appends fromfile to tofile and deletes fromfile.
    >> fromfile:str = name of file to be copied and deleted
    >> tofile:str = name of file to be deleted
    >> deletefrom:bool = whether fromfile to be deleted [True]
    >> append:bool = whether to append tofile [True]
    '''
    if not fromfile or not tofile or not checkForFile(fromfile): return False
    #if append: open(tofile,'a').writelines(open(fromfile,'r').readlines())
    #else: open(tofile,'w').writelines(open(fromfile,'r').readlines())
    if append: open(tofile,'a').write(open(fromfile,'r').read())
    else: open(tofile,'w').write(open(fromfile,'r').read())
    if deletefrom: os.unlink(fromfile)
    return True
#########################################################################################################################
def subDir(pathname,exclude=[]):   ### Returns the subdirectories given by glob.glob(), i.e. no *.* returned
    '''
    Returns the subdirectories given by glob.glob()
    >> pathname:pathname for glob.glob()
    >> exclude:list of directories to leave out of list
    '''
    dirlist = glob.glob('%s*' % makePath(pathname))
    subdir = []
    for element in dirlist[0:]:
        if os.path.isdir(element):      #X#.find('.') < 0:
            subdir.append(string.split(element,os.sep)[-1])
    for ex in exclude:
        if ex in subdir: subdir.remove(ex)
    return subdir
#########################################################################################################################
def stripPath(path): return os.path.basename(path)
def baseFile(filename,strip_path=False,extlist=[],keepext=False):   ### Returns file without extension, with or without path
    '''
    Returns file without extension, with or without path.
    >> filename:str = file to reduce to basefile
    >> strip_path:bool = whether to strip any path information [False]
    >> extlist:list of acceptable file extensions to remove []
    >> keepext:bool = default keep extension behaviour in no extlist
    << basefile:str = returned filename base
    '''
    (basefile,_ext) = os.path.splitext(filename)
    if extlist and _ext not in extlist and _ext[1:] not in extlist: basefile = filename  # Check _ext with(out) leading .
    elif keepext and not extlist: basefile = filename
    if strip_path: basefile = os.path.basename(basefile)
    return basefile
#########################################################################################################################
def listDir(callobj=None,folder=os.getcwd(),subfolders=True,folders=True,files=True,summary=True,asksub=False,dircut=0,dirdepth=-1):
    '''
    Returns a full list of files and/or folders using os.listdir().
    >> callobj:RJE_Object = object used for verbosity etc. (if any) [None]
    >> folder:str = folder to start looking in (for os.listdir(folder)) [Current directory]
    >> subfolders:bool = whether to also look in subfolders
    >> folders:bool [True] = whether to return folders in list
    >> files:bool [True] = whether to return files in list
    >> summary:bool [True] = whether to print output summary
    >> asksub:bool [False] = whether to ask for confirmation before scanning subdir
    >> dircut:int [0] = Maximum number of subdirectories to delve into [<1 = all]
    >> dirdepth:int [-1] = Maximum depth of subdirectories to delve into [<0 = all]
    << pathlist:list of strings = paths to files from current directory
    '''
    try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not folders and not files: return []
        asksub = asksub and (not callobj or callobj.i() >= 0)
        if folders and files: sumtxt = 'files/folders'
        elif folders: sumtxt = 'folders'
        else: sumtxt = 'files'
        sumtxt += ' read from %s' % os.path.basename(folder)
        pathlist = []; pathx = 0   # List of files to return
        ## ~ [0a] Check Subfolders ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        allpaths = os.listdir(folder)
        subfolders = subfolders and dirdepth
        if subfolders:
            dirx = 0
            for path in allpaths:
                fullpath = os.path.join(folder,path)
                if os.path.isdir(fullpath): dirx += 1
            if (dirx > dircut > 0):
                subfolders = False
                if callobj: callobj.printLog('\r#LSDIR','%s subdir in %s will be skipped (dircut=%s)' % (iStr(dirx),folder,iStr(dircut)))
            else: sumtxt += ' (and subdir)'
        ### ~ [1] Build directory listing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for path in allpaths:
            if summary and callobj: callobj.progLog('\r#LSDIR','%s %s' % (iStr(pathx),sumtxt))
            fullpath = os.path.join(folder,path)
            if os.path.isdir(fullpath):
                if folders: pathlist.append(fullpath); pathx += 1
                if subfolders and (not asksub or yesNo('Scan %s?' % fullpath)): 
                    pathlist += listDir(callobj,fullpath,subfolders,folders,files,False,asksub,dircut,dirdepth-1)
                    pathx = len(pathlist)
            elif files: pathlist.append(fullpath); pathx += 1
        if summary and callobj: callobj.printLog('\r#LSDIR','%s %s' % (iStr(pathx),sumtxt))
    except:
        if callobj: callobj.errorLog('Error in rje.listDir(%s,%s)' % (folder,subfolders),printerror=True,quitchoice=False)
    return pathlist
#########################################################################################################################
def getFileList(callobj=None,folder=os.getcwd(),filelist=['*'],subfolders=True,summary=True,filecount=0,asksub=False):   ### Returns a list of files 
    '''
    Returns a list of files with appropriate filenames.
    >> callobj:RJE_Object = object used for verbosity etc. (if any) [None]
    >> folder:str = folder to start looking in (for os.listdir(folder)) [Current directory]
    >> filelist:list of files ['*']
    >> subfolders:bool = whether to also look in subfolders
    >> summary:bool [True] = whether to print output summary
    >> filecount:int [0] = running total of files so far for progress reporting
    >> asksub:bool [False] = whether to ask for confirmation before scanning subdir
    << globlist:list of strings = paths to files from current directory
    '''
    try:
        ### Setup ###
        asksub = asksub and (not callobj or callobj.i() >= 0)
        globlist = []   # List of files to return
        ### Get files for this directory
        for file in filelist:     # Each file in turn
            for newfile in glob.glob(os.path.abspath(os.path.join(folder,file))):    # Looks in folder/file
                if newfile not in globlist: globlist.append(newfile)
        if callobj and summary:
            callobj.progLog('\r#FILES','Getting files: %5s' % integerString(filecount+len(globlist)))
        ### Look in subdirectories? ###
        if subfolders:
            for file in os.listdir(folder):     # Get full list of directory contents
                fullfile = os.path.join(folder,file)
                #print file, os.path.isdir(os.path.join(folder,file))
                if os.path.isdir(fullfile) and (not asksub or yesNo('Scan %s?' % fullfile)):     # Found a subdirectory
                    for newfile in getFileList(callobj,fullfile,filelist,subfolders,summary,filecount=(filecount+len(globlist)),asksub=asksub):    # Perform iterative calling of method
                        if newfile not in globlist: globlist.append(newfile)
                    if callobj and summary: callobj.progLog('\r#FILES','Getting files: %5s' % integerString(filecount+len(globlist)))
                        
        ### Return data ###
        return globlist
    except:
        if callobj: callobj.log.errorLog('Error in rje.getFileList(%s,%s,%s)' % (folder,filelist,subfolders),printerror=True,quitchoice=False)
        return globlist
#########################################################################################################################
def nextLine(FILEHANDLE=None,strip=True):   ### Returns next line or None if no line
    '''
    Returns next line or None if end of file.
    >> FILEHANDLE:File object
    >> strip:boolean = whether to strip \r and \n from line [True]
    '''
    if not FILEHANDLE: return None
    line = FILEHANDLE.readline()
    if not line: return None
    if strip: return line.strip('\r\n')
    else: return line
#########################################################################################################################
def fileLineFromSeek(FILE=None,pos=0,reseek=False,next=False):  ### Returns full line & new pos for seek position pos
    '''
    Returns full line & new pos for seek position pos.
    >> FILE:Open file object for reading
    >> pos:int [0] = position to be within line
    >> reseek:boolean [False] = seek to new position before returning
    >> next:boolean [False] = return next line, rather than line containing pos
    << returns (line(str),pos(long))
    '''
    ### Setup ###
    preline = 'ni'
    line = 'i'
    if pos == 0:
        FILE.seek(pos)
        line = FILE.readline()
    else: pos += 1
    ### Find beginning of line ###
    while preline.find(line) >= 0 and pos > 0:
        pos -= 1
        FILE.seek(pos)
        line = FILE.readline()        
        FILE.seek(pos-1)
        preline = FILE.readline()
    ### Return ###
    if next: return linePosFromPos(FILE,FILE.tell(),reseek)
    if reseek: FILE.seek(pos)
    return (line,pos)
#########################################################################################################################
def lineFromIndex(target,file,re_index='^(\S+)\s',sortunique=False,xreplace=True):    ### Returns line matching target
    '''
    Returns line matching target from sorted file.
    >> target:str = String matching the first group returned by re_index.
    >> file:str = Name of file in which to locate line.
    >> re_index:str ['^(\S+)\s'] = RegExp matching sorted element to match.
    >> sortunique:boolean [False] = whether to attempt to match the UNIX sort as best possible (dodgy)
    >> xreplace:boolean [True] = whether to replace dots with xs for determining sort 
    << returns line if found, else '' if not found.
    '''
    FILE = open(file,'r')
    ipos = posFromIndex(target,FILE,re_index=re_index,sortunique=sortunique,xreplace=xreplace)
    if ipos >= 0: iline = fileLineFromSeek(FILE,ipos)[0]
    else: iline = ''
    FILE.close()
    return iline
#########################################################################################################################
def endPos(FILE=None,filename=None):    ### Returns end position from file
    '''Returns end position from file.'''
    if FILE: fpos = FILE.tell()
    else: FILE = open(filename,'r'); fpos = -1
    FILE.seek(0,2)
    fend = FILE.tell()
    if fpos < 0: FILE.close()
    else: FILE.seek(fpos)
    return fend
#########################################################################################################################
def posFromIndex(target,INDEX,start_pos=0,end_pos=-1,re_index='^(\S+)=',sortunique=False,xreplace=True):    ### Returns position to extract target
    '''
    Returns position in file from which to extract target using fileLineFromSeek(INDEX,pos).
    NB. Index should be sorted and every line from start_pos on should match re_index.
    >> target:str = target string from INDEX using re_index.
    >> INDEX:filehandle open for reading.
    >> start_pos:int [0] = position in file to start looking
    >> end_pos:int [0] = position at end of file (seek(0,2)->tell())
    >> re_index:str ['^(\S+)='] = regular expression to use to identify match to target
    >> sortunique:boolean [False] = whether to attempt to match the UNIX sort as best possible (dodgy)
    >> xreplace:boolean [True] = whether to replace dots with xs for determining sort 
    << returns ipos for use with fileLineFromSeek, else -1 if not found.
    '''
    ### Setup ###
    if end_pos < start_pos:
        INDEX.seek(0,2)
        end_pos = INDEX.tell()
    (ipos,jpos) = (start_pos,end_pos)   # Positions used in search to home in on target
    check_me = None   # attribute to compare to target
    ## Special RE ##
    leader = ''
    if re_index[:2] != '^(':
        leader = string.split(re_index,'(')[0]
        if leader[:1] == '^':
            leader = leader[1:]

    ### Search ###
    while 'I am a looping wonder':
        ## Check entry at current ipos ##
        line = fileLineFromSeek(INDEX,ipos)[0]  # String of INDEX line
        if line and matchExp(re_index,line):    # Match
            check_me = matchExp(re_index,line)[0]
            assess_me = leader + check_me
        elif line and matchExp('^(\S+)=(\S.+)$',line):    # Backup Match
            check_me = matchExp('^(\S+)=(\S.+)$',line)[0]
            assess_me = check_me
        else:   # Cannot match - problem with INDEX or re_index
            check_me = None
            break
        if check_me == target:  # Found!
            break
        ## Assess whether target earlier or later in file ##
        if sortunique:
            assessment = sortUnique([leader+target,assess_me],xreplace=xreplace)
        else:
            assessment = [leader+target,assess_me]
            assessment.sort()
        if assessment[0] == assess_me:   ## This line is before the target: move on
            start_pos = ipos
        else:   ## This line is after the accession number: move back
            jpos = ipos
        ## Redefine ipos and repeat ##
        ipos = start_pos + ((jpos - start_pos) / 2)
        if ipos in [start_pos,jpos]:   # target not found in Index
            ipos = -1
            break

    ### Return ###
    if check_me: return ipos
    else: return -1
#########################################################################################################################
def urlToFile(sourceurl,filename,callobj,appendable=True,backupfile=True,log=True):  ### Downloads URL to file
    '''
    Downloads URL to file with backup and append options taken from callobj.
    >> sourceurl,filename,callobj,appendable=True,backupfile=True,log=True
    << True/False
    '''
    try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        mkDir(callobj,filename)
        if backupfile and callobj: backup(callobj,filename,appendable=appendable)
        ### ~ [1] Download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if log and callobj: callobj.progLog('#URL','Downloading %s' % sourceurl)
        open(filename,'a').write(urllib2.urlopen(sourceurl).read())
        if log and callobj: callobj.printLog('\r#URL','Downloaded %s -> %s' % (sourceurl,filename),log=log)
    except:
        if callobj: callobj.errorLog('urlToFile error!'); return False
        raise
#########################################################################################################################
def backup(callobj,filename,unlink=True,appendable=True,autobackup=False):  ### Checks for existence of file and gives backup options
    '''
    Checks for existence of file and gives backup options.
    >> callobj:Object calling the method (use for interactive/append)
    >> filename:str = filename to backup
    >> unlink:boolean [True] = whether to delete file if found but not backed up (backup will delete anyway)
    >> autobackup:boolean [False] = whether to backup file (and overwrite backup) if i<0.
    '''
    if (appendable and callobj.getBool('Append')) or not os.path.exists(filename): return False      # No need
    backupfile = callobj.getAtt('opt','Backups',default=True)
    if backupfile and callobj.i() < 0 and autobackup: backupfile = '%s.bak' % filename
    elif backupfile and callobj.i() >= 0 and yesNo('%s already exists. Backup?' % filename):
        backupfile = '%s.bak' % filename
        while os.path.exists(backupfile) and not yesNo('%s already exists. Overwrite?' % backupfile):
            backupfile = ''
            while not backupfile: backupfile = choice('New name for backup file?')
    else: backupfile = False
    if backupfile:
        if os.path.exists(backupfile): os.unlink(backupfile)
        os.rename(filename,backupfile)
        callobj.printLog('#BAK','%s backed up as %s' % (filename,backupfile))
        return True
    elif unlink: os.unlink(filename)
    return False
#########################################################################################################################
def deleteDir(callobj,deldir,contentsonly=True,confirm=True,report=True):   ### Deletes directory contents
    '''
    Deletes directory contents, including subdirectories. Use with care!!
    >> callobj:Object calling the method (use for interactive/append)
    >> deldir:str = path to directory to delete
    >> contentsonly:boolean [True] = whether to delete files within directory only or directory too
    >> confirm:boolean [True] = ask for confirmation before deleting files
    >> report:boolean [True] = whether to summarise deletion in log if callobj given
    << True if deleted, False if not, KeyboardInterrupt error if mind changed
    '''
    ### Setup ###
    if not os.path.exists(deldir): return False    # No need to delete
    files = getFileList(callobj,deldir,filelist=['*'],subfolders=True,summary=False)
    if contentsonly and not files: return False
    dtxt = '%s files from %s' % (len(files),deldir)
    if not contentsonly: dtxt += ' (and %s)' % deldir
    ### Confirm ###
    if confirm and callobj and (callobj.i() >= 0 or callobj.getBool('DeBug')):
        print '\n%d files in %s:\n - %s\n' % (len(files),deldir,string.join(files,'\n - '))
        if not yesNo('Delete %s?' % dtxt): raise KeyboardInterrupt
    ### Delete files ###
    for f in files:
        if not os.path.isdir(f): os.unlink(f)
    if not contentsonly:
        for f in files:
            if os.path.isdir(f): os.removedirs(f)
        os.removedirs(deldir)
    if callobj and report: callobj.printLog('#DEL','Deleted %s' % dtxt)
    return True
#########################################################################################################################
def mkDir(callobj,newdir,log=False):     ### Makes directory and necessary parent directories
    '''Makes directory and necessary parent directories.'''
    try:
        ### Setup ###
        addback = []
        makedir = os.path.abspath(os.path.split(newdir)[0])    # Full path of directory part of newdir.
        ### Find missing directories ###
        while not os.path.exists(makedir):
            (makedir,missing) = os.path.split(makedir)
            addback.insert(0,missing)
        ### Make missing directories ###
        while addback:
            makedir = os.path.join(makedir,addback.pop(0))
            try:
                if log: callobj.printLog('#MKDIR','mkdir %s' % makedir); os.mkdir(makedir)
                else: os.mkdir(makedir)
            except:
                if not os.path.exists(makedir): raise   # Might just be 2+ programs making same directory at same time.
    except:
        try: callobj.errorLog('Problem with rje.mkDir(%s)' % newdir,quitchoice=True)
        except: raise
#########################################################################################################################
def cleanDir(callobj=None,keepfiles=[],cleandir='',log=True): ### Cleanup directory by removing files not in list
    '''
    Cleanup directory by removing files not in list.
    >> callobj:Object = calling RJE Object
    >> keepfiles:list = Files to ignore and not delete
    >> cleandir:str [] = Directory to cleanup
    '''
    try:### ~ [1] Glob current files and try deleting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        remx = 0
        for file in glob.glob('%s*' % cleandir):
            if file in keepfiles or os.path.isdir(file): continue
            os.unlink(file); remx += 1
        if callobj and log: callobj.printLog('#DEL','%s files deleted from %s' % (integerString(remx),makePath(cleandir,return_blank=False)))
    except:
        try: callobj.errorLog('Problem with rje.cleanDir()')
        except: raise
#########################################################################################################################
def humanFileSize(filename): return humanByteSize(os.path.getsize(filename))
#########################################################################################################################
def humanByteSize(nbytes):
    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    if nbytes == 0: return '0 B'
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.0
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])
#########################################################################################################################
###  End of File Manipulation Functions                                                                                 #
#########################################################################################################################


#########################################################################################################################
###  Binary Counting Functions                                                                                          #
#########################################################################################################################
def binaryCount(binlist):   ### Adds one in binary and returns binlist
    '''
    Adds one in binary and returns binlist.
    >> binlist:list of 0s and 1s
    '''
    i = len(binlist) - 1
    while binlist[i] == 1 and i >= 0:
        binlist[i] = 0
        i -= 1
    if binlist[i] == 0: binlist[i] = 1
    return binlist
#########################################################################################################################
def binComb(positions,cmin,cmax): ### Returns the number of binary combinations within certain ranges
    '''
    Returns the number of binary combinations within certain ranges.
    >> positions:int = the number of positions that can be 1 or 0
    >> min:int = the min. no. of 1s to have
    >> max:int = the max. no. of 1s to have
    '''
    bincomb = 0
    binlist = [0] * positions
    while sum(binlist) < len(binlist):
        if sum(binlist) <= cmax and sum(binlist) >= cmin: bincomb += 1
        binlist = binaryCount(binlist)
    if sum(binlist) <= cmax and sum(binlist) >= cmin: bincomb += 1
    return bincomb
#########################################################################################################################
###  End of Binary Counting Functions                                                                                   #
#########################################################################################################################



#########################################################################################################################
###  Delimited Text Functions                                                                                           #
#########################################################################################################################
def getDelimit(cmd_list=[],default='\t'):   ### Returns delimit from command list
    '''Returns delimit from command list.'''
    delimit = default
    for cmd in cmd_list:
        if cmd.find('delimit=') == 0:
            delimit = string.replace(cmd[len('delimit='):],'\\t','\t')
    if delimit.lower() == 'tab': delimit = '\t'
    return delimit
#########################################################################################################################
def delimitExt(delimit):    ### Returns file extension for text file with given delimiter
    '''Returns file extension for text file with given delimiter.'''
    if delimit == ',': return 'csv'
    elif delimit == '\t': return 'tdt'
    else: return 'txt'
#########################################################################################################################
def delimitFromExt(ext='',filename='',write=False):    ### Returns default delimiter for given file extension
    '''Returns default delimiter for given file extension.'''
    if not ext and filename: ext = os.path.splitext(filename)[1]
    if ext in ['.csv','csv']: return ','
    elif ext in ['.tsv','tsv','.tdt','tdt','.tab','tab']: return '\t'
    elif write and ext not in ['txt','.txt']: return '\t'
    else: return ' '
#########################################################################################################################
def writeDelimit(OUTFILE=None,outlist=[],delimit='\t',outfile=None): ### Writes given list of strings to file (if given) or just returns
    '''
    Writes given list of strings to file (if given) or just returns string.
    >> OUTFILE:file handle = output file handle [None]
    >> outlist:list of string objects to write to file []
    >> delimit:text limiter for file [tab]
    >> outfile:string = name of file to append (use if OUTFILE not given)
    << returns the output string
    '''
    if outfile:
        try: OUTFILE = open(outfile,'a')
        except: OUTFILE = None
    writelist = []
    for element in outlist:
        element = '%s' % element
        if element.find(delimit) >= 0:
            element.replace('"','\"')
            element = '"%s"' % element
        writelist.append(element)
    if OUTFILE:
        OUTFILE.write('%s\n' % string.join(writelist,delimit))
        if outfile: OUTFILE.close()
    return string.join(writelist,delimit)
#########################################################################################################################
def readDelimit(line='',delimit='\t'):  ### Returns list of strings from file line, removing "" where necessary
    '''
    Reads a line into a list of strings.
    >> line:string object line from input file 
    >> delimit:text limiter for file [tab]
    << returns list of strings
    '''
    line = chomp(line)
    if not line: return []
    if delimit in [' ','\s+']: splitlist = line.split()
    else: splitlist = line.split(delimit)
    readlist = []
    s = 0
    while s < len(splitlist):
        if splitlist[s] == '': readlist.append(splitlist[s])
        elif splitlist[s][0] == '"' and splitlist[s][-1] == '"': readlist.append(splitlist[s][1:-1])
        elif splitlist[s][0] == '"':    # Long, split string to handle
            minilist = []               # This stores the split entry to be recombined
            minilist.append(splitlist[s][1:])
            s += 1
            while s < len(splitlist):   # Now need to look for the end of the " "
                if not splitlist[s]:
                    s += 1
                    continue
                elif splitlist[s][-1] == '"':
                    minilist.append(splitlist[s][:-1])
                    break
                elif splitlist[s][-1] == ' ':
                    splitlist[s] = splitlist[s][:-1]
                    continue
                minilist.append(splitlist[s])
                s += 1
            readlist.append(string.join(minilist,delimit))
        else: readlist.append(splitlist[s])
        s += 1
    return readlist
#########################################################################################################################
def delimitedFileOutput(callobj,filename,headers,delimit=None,datadict={},rje_backup=False):   ### Outputs a single line of compiled data to file
    '''
    Outputs a single line of compiled data to file.
    >> callobj:Object = calling object
    >> filename:str = Name of output file or open file handle
    >> headers:list of field headers
    >> delimit:str = text delimiter. If None, will try to work out from filename and callobj
    >> datadict:dictionary of {header:str} = data to output. If none, will output headers themselves.
    >> rje_backup:boolean [False] = if no datadict and true, will call rje.backup(callobj,filename,unlink=True)
    '''
    ### Delimiter ###
    if not delimit: delimit = getDelimit(callobj.cmd_list,default=delimitFromExt(filename=filename))
    ### Output Data ###
    if datadict:
        outlist = []
        for h in headers:
            if datadict.has_key(h): outlist.append(str(datadict[h]))
            else: outlist.append('')
    else:
        if rje_backup: backup(callobj,filename,unlink=True)
        if callobj.getOpt('Append') and os.path.exists(filename): return  # Don't output headers
        if callobj.getOpt('MySQL'): outlist = string.split(string.join(headers,',').lower(),',')
        else: outlist = headers[0:]
    ### Write to file ###
    if string.split('%s' % filename)[0] == '<open': writeDelimit(filename,outlist=outlist,delimit=delimit)
    else: writeDelimit(outlist=outlist,delimit=delimit,outfile=filename)
#########################################################################################################################
def delimitedObjDataOutput(callobj,filename,headers,delimit=None,dpdict={}):    ### Outputs object data as single delimited line
    '''
    Outputs object data as single delimited line.
    >> callobj:Object containing data in obj.dict['Data'], obj.stat, obj.opt and obj.info
    >> filename:str = Name of output file
    >> headers:list of field headers. Should correspond to obj.getData() keys
    >> delimit:str = text delimiter. If None, will try to work out from filename and callobj
    >> dpdict:dictionary of dp settings (integers) for obj.getData() and corresponding lists of headers. (All else -1)
    '''
    ### Data Dictionary ###
    datadict = {}
    for dp in sortKeys(dpdict):
        for stat in dpdict[dp]:
            datadict[stat] = callobj.getData(stat,default='',dp=dp)
    for stat in headers:
        if not datadict.has_key(stat):
            datadict[stat] = callobj.getData(stat,default='')
    delimitedFileOutput(callobj,filename,headers,delimit,datadict)
#########################################################################################################################
def dataDict(callobj,filename,mainkeys=[],datakeys=[],delimit=None,headers=[],getheaders=False,ignore=[],lists=False,debug=False,enforce=False,uselower=False): ### Extracts data from delimited file into dictionary
    '''
    Extracts data from delimited file into dictionary.
    >> callobj:RJE_Object controlling logs and error-handling
    >> filename:str = file to read from
    >> mainkeys:list = List of headers to be used as key for returned dictionary. If None, will use first header.
    >> datakeys = List of headers to be used as keys for data returned for each mainkey (all headers if [])
    >> delimit = string delimiter. If None, will identify from filename
    >> headers = List of headers to use instead of reading from first line or use datakeys
    >> getheaders:bool [False] = whether to add an extra 'Headers' dictionary value containing list of headers
    >> ignore:list = Leading strings for lines to ignore (e.g. #)
    >> lists:bool [False] = whether to return values as lists. (Otherwise, later entries will overwrite earlier ones)
    >> enforce:bool [False] = whether to enforce correct length of delimited lines [True] or [False] truncate/extend
    >> uselower:bool [False] = Whether to convert headers into lower case.
    '''
    try:
        ### ~ Setup method variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if filename.lower() in ['','none'] or not os.path.exists(filename): return {}
        ## Delimiter ##
        datadict = {}; warnx = 0; warn10 = []
        if not delimit: delimit = delimitFromExt(filename=filename)
        FILE = open(filename,'r')
        ## Get headers, if not given ##
        while not headers:
            fline = FILE.readline()
            if ignoreLine(fline,ignore): continue
            if uselower: fline = fline.lower()
            headers = readDelimit(fline,delimit)
            while headers and not headers[-1]: headers.pop(-1)
        ## Establish keys for data dictionary. Mainkeys form dictionary keys. Datakeys determine which data is returned.
        autoid = mainkeys in ['All','all','auto','#']
        if not mainkeys: mainkeys = headers[:1]
        if not datakeys:
            datakeys = headers[0:]
            for key in mainkeys: datakeys.remove(key)
        if datakeys in ['All','all']: datakeys = headers[0:]
        keylen = 0
        #callobj.debug('%s <- %s' % (mainkeys,headers))
        if autoid: keylen = len(headers) - 1
        else: 
            for key in mainkeys:
                try: keylen = max(keylen,headers.index(key))
                except:
                    if callobj: callobj.errorLog('Key "%s" not found in %s headers: %s' % (key,filename,headers))
                    raise ValueError
        
        ### ~ Read data from file into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        fline = FILE.readline(); ix = 0
        while fline:
            try:
                #x#callobj.deBug('%s -> %d' % (fline,len(datadict)))
                ## Check whether line is to be ignored ##
                if ignoreLine(fline,ignore):
                    fline = FILE.readline()
                    continue
                ## Convert to data list and check for headers ##
                data = readDelimit(fline,delimit)
                if (len(data) != len(headers) and enforce) or len(data) < keylen:
                    fline = FILE.readline()
                    continue
                linedata = {}
                for h in range(len(headers)):
                    try: linedata[headers[h]] = data[h]
                    except: linedata[headers[h]] = ''
                ## Main Key ##
                ix += 1
                if autoid: mainkey = ix
                else:
                    mainkey = []
                    for key in mainkeys: mainkey.append(linedata[key])
                    if not string.join(mainkey,''): fline = FILE.readline(); continue
                    mainkey = string.join(mainkey,delimit)
                    #x#callobj.deBug('%s:%s' % (mainkey,datadict.has_key(mainkey)))
                if not datadict.has_key(mainkey): datadict[mainkey] = {}
                elif not lists: 
                    warnx += 1
                    if warnx <= 10: warn10.append(string.replace(mainkey,'\t',','))
                    if debug: callobj.deBug('Dup: %s' % mainkey)
                ## Other Data ##
                for key in datakeys:
                    if lists:
                        if key not in datadict[mainkey]: datadict[mainkey][key] = []
                        if linedata[key] and linedata[key] not in datadict[mainkey][key]: datadict[mainkey][key].append(linedata[key])
                    else: datadict[mainkey][key] = linedata[key]
                fline = FILE.readline()
            except: callobj.deBug(fline); raise
            if debug: callobj.deBug(fline); callobj.deBug(datadict)

        ### ~ Finish and return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        FILE.close()
        if warnx: 
            callobj.printLog('\r#WARN','Warning: %s entries overwritten due to common key (lists=False)' % integerString(warnx))
            if warnx > 10: callobj.printLog('\r#WARN','Dups: %s ...' % string.join(warn10,' | '))
            else: callobj.printLog('\r#WARN','Dups: %s.' % string.join(warn10,' | '))
        if getheaders: datadict['Headers'] = headers[0:]
        return datadict            
    except:
        callobj.errorLog('Problem with rje.dataDict(%s)' % filename,quitchoice=False); raise
#########################################################################################################################
def ignoreLine(fileline,ignorelist):    ### Returns whether line should be ignored
    '''Returns whether line should be ignored.'''
    for i in ignorelist:
        if fileline.find(i) == 0: return True
    return False
#########################################################################################################################
### End of Delimited Text Functions                                                                                     #
#########################################################################################################################



#########################################################################################################################
###  Input Command Functions                                                                                            #
#########################################################################################################################
def longCmd(cmd_list):  ### Extracts long command from command list and returns altered command list
    '''
    Extracts long command from command list and returns altered command list.
    cmd_list:list of commandline arguments
    '''
    longcmd = []
    i = 0
    while i < len(cmd_list):
        cmd = cmd_list[i]
        if re.search('^(\S+=)"(.+)$',cmd):
            thiscmd = [string.join(matchExp('^(\S+=)"(.+)$',cmd),'')]
            while cmd_list[i][-1] != '"' and i < (len(cmd_list) - 1):
                i += 1
                thiscmd.append(cmd_list[i])
            longcmd.append(string.join(thiscmd,' ')[:-1])
        elif matchExp('^"(.+)$',cmd):
            thiscmd = [matchExp('^"(.+)$',cmd)[0]]
            while cmd_list[i][-1] != '"' and i < (len(cmd_list) - 1):
                i += 1
                thiscmd.append(cmd_list[i])
            longcmd.append(string.join(thiscmd,' ')[:-1])
        else: longcmd.append(cmd)
        if len(longcmd) > 1 and longcmd[-2][:1] == '-':
            try: float(longcmd[-1]); numval = True
            except: numval = False
            if '=' in longcmd[-1] and not longcmd[-1].startswith('http:'): longcmd[-2] = longcmd[-2][1:]
            elif longcmd[-1][:1] == '-' and not numval: longcmd[-2] = longcmd[-2][1:]
            else:
                longcmd[-2] = '%s=%s' % (longcmd[-2][1:],longcmd[-1])
                longcmd.pop(-1)
        i += 1
    if longcmd:
        try: float(longcmd[-1]); numval = True
        except: numval = False
        if longcmd[-1][:1] == '-' and not numval and '=' not in longcmd[-1]: longcmd[-1] = longcmd[-1][1:]
    if 'h' in longcmd: longcmd.append('help')
    return longcmd
#########################################################################################################################
def iniCmds(ini_path,ini_file,iowarning=True,altpaths=True):  ### Reads a list of commands from an inifile.         #V4.2
    '''
    Reads a list of commands from an inifile
    >> ini_path:str = filepath. 
    >> ini_file:str = filename. Add ini_path if not found in current directory
    >> iowarning:boolean [True] = whether to warn and kill if *.ini missing.
    >> altpaths:bool [True] = whether to look in alternative paths for ini_file
    << inicmds:list = list of commands
    '''
    ### ~ [0] ~ Setup INI file path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    if not ini_file: return []
    path_ini = makePath(ini_path) + ini_file    # First place to look after current directory is source directory
    alt_ini = makePath(ini_dir) + ini_file      # Next place to look is settings/ directory for overall defaults
    ## ~ [0a] ~ Check for INI File and warn if missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
    if checkForFile(ini_file): openfile = ini_file
    elif altpaths and checkForFile(path_ini): openfile = path_ini
    elif altpaths and checkForFile(alt_ini): openfile = alt_ini
    elif iowarning:
        print '#ERR\tINI file %s not found! (Also checked %s and %s)' % (ini_file,ini_path,ini_dir)
        if not yesNo('Continue without %s arguments?' % ini_file): os._exit(0)
        return []
    else: return []
    ### ~ [1] ~ Load arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    ilines = open(openfile,'r').readlines()
    inicmds = []
    for iniline in ilines:
        if iniline.find('#') == 0: continue  # Comment line
        if iniline.find(' #') >= 0: # Comment in line
            comment = iniline.find(' #')
            iniline = iniline[:comment]
        elif iniline.find('\t#') >= 0: # Comment in line
            comment = iniline.find('\t#')
            iniline = iniline[:comment]
        else: iniline = chomp(iniline)
        inicmds += iniline.split(None,)
    return inicmds
#########################################################################################################################
def getCmdList(argcmd,info=None):   ### Converts arguments into list of commands, reading from ini file as appropriate.
    '''
    Converts arguments into list of commands, reading from ini file as appropriate.
    >> argcmd:list = list of commands from commandline separated by whitespace (typically sys.argv[1:])
    >> info:Info Object [None] = if given, will try to read defaults from 'info.program.ini'
    << cmd_list:list = list of commands!
    '''
    cmd_list=[]
    path = ''
    server = 'bioware_server' in argcmd
    ### <a> ### Defaults - hierarchy of rje.ini, defaults.ini, global program.ini and local program.ini
    try:
        if info != None:
            defaults = info.program.lower() + '.ini'
            rjepath = makePath(os.path.split(sys.argv[0])[0]); path = rjepath
            if server: argcmd = longCmd(iniCmds(rjepath,'defaults.ini',iowarning=False)) + longCmd(iniCmds(rjepath,defaults,iowarning=False)) + argcmd
            else: argcmd = longCmd(iniCmds(rjepath,'defaults.ini',iowarning=False)) + longCmd(iniCmds(rjepath,'rje.ini',iowarning=False)) + longCmd(iniCmds(rjepath,defaults,iowarning=False)) + argcmd
    except:
        print '#ERR\tMajor error in getCmdList() looking for defaults:', sys.exc_info()[0]
        errorMsg()
        if yesNo('Proceed?') == False: os._exit(0)
    ### <b> ### System Arguments
    try:
        while argcmd:
            nextcmd = argcmd.pop(0)
            #X#print '#', nextcmd, '#'
            cmd_list.append(nextcmd)
            if nextcmd[:5].lower() in ['path=']: path = makePath(path=nextcmd[5:],wholepath=False)
            if argcmd and nextcmd[:5].lower() in ['-path']: path = makePath(path=argcmd.pop(0),wholepath=False)
            if nextcmd[:4].lower() in ['ini=']: argcmd = longCmd(iniCmds(path,nextcmd[4:])) + argcmd
            if argcmd and nextcmd[:4].lower() in ['-ini']: argcmd = longCmd(iniCmds(path,argcmd.pop(0))) + argcmd
        return longCmd(cmd_list)
    except:
        print '#ERR\tMajor error in getCmdList() processing system arguments:', sys.exc_info()[0]
        errorMsg()
        if yesNo('Proceed?') == False: os._exit(0)
        return longCmd(cmd_list)
#########################################################################################################################
def inputCmds(cmd_out,cmd_list):    ### Reads extra commands from user prompt
    '''
    Reads extra commands from user prompt
    >> cmd_out:Out 
    >> cmd_list:list = current list of commands
    << newcmd:list = new list of commands
    '''
    cmd_out.verbose(0,3,"Arguments: " + str(cmd_list),1)
    newcmd = longCmd(choice('Additional Arguments?: ').split(None,))
    cmd_out.verbose(0,3,"New Arguments: " + str(newcmd),1)
    return getCmdList(newcmd)
#########################################################################################################################
def setLog(info,out,cmd_list,printlog=True,fullcmd=True):  ### Makes Log Object and outputs general program run info to start of file. Returns log.
    '''
    Makes Log Object and outputs general program run info to start of file. Returns log.
    >> info:rje.Info object containing program info
    >> out:rje.Out object controlling output to screen
    >> cmd_list:list = full list of commandline options
    >> printlog:bool [True] = whether to print summaries to the log
    >> fullcmd:bool [True] = whether full command output should be included in log output
    << log:Log object
    '''
    ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    out.verbose(2,3,'Commands: %s' % cmd_list,2)
    log_file = 'log=%s.log' % info.program.lower()
    argcmd = longCmd(sys.argv[1:])
    cmd_list = [log_file] + cmd_list[0:]
    ### ~ [2] ~ Make and return Log ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        log = Log(info.start_time, cmd_list)
        log.obj = {'Info':info}
        log.info['Name'] = info.program
        if printlog:
            log.printLog('#~~#','#~~#',timeout=False,screen=False)
            log.printLog('#LOG','Activity Log for %s V%s: %s' % (info.program,info.version,time.asctime(time.localtime(info.start_time))),screen=False)
            log.printLog('#DIR','Run from directory: %s' % os.path.abspath(os.curdir),screen=False)
            log.printLog('#ARG','Commandline arguments: %s' % string.join(argcmd),screen=False)
            if fullcmd: log.printLog('#CMD','Full Command List: %s' % argString(tidyArgs(cmd_list)),screen=False)
            log.printLog('#VIO','Verbosity: %d; Interactivity: %d.' % (log.v(),log.i()))
            if log.dev() or log.test(): log.printLog('#DEV','Development mode: %s; Testing mode: %s.' % (log.dev(),log.test()))
            if log.info['ErrorLog'].lower() not in ['','none']:
                log.printLog('#~~#','#~~#',timeout=False,screen=False,error=True)
                log.printLog('#LOG','Error Log for %s %s: %s' % (info.program,info.version,time.asctime(time.localtime(info.start_time))),screen=False,error=True)
                log.printLog('#DIR','Run from directory: %s' % os.path.abspath(os.curdir),screen=False,error=True)
                log.printLog('#ARG','Commandline arguments: %s' % str(argcmd),screen=False,error=True)
                if fullcmd: log.printLog('#CMD','Full Command List: %s' % argString(tidyArgs(cmd_list)),screen=False,error=True)
                if log.dev() or log.test(): log.printLog('#DEV','Development mode: %s; Testing mode: %s.' % (log.dev(),log.test()))
                elif not log.warn(): log.warnLog('Runtime warnings switched off (warn=F).')
        return log
    except:
        print 'Log problem'
        raise
#########################################################################################################################
def listFromCommand(command,checkfile=True):   ### Returns a list object from a given command string
    '''
    Returns a list object from a given command string. Reads list from file if found. If not, will split on commas.
    >> command:string
    >> checkfile:boolean = whether to check for presence of file and read list from it
    << comlist:list
    '''
    try:
        ### Empty list if None ###
        if command.lower() == 'none': return []
        ### Check for file ###
        if checkfile and os.path.exists(command) and not os.path.isdir(command):
            comlist = open(command,'r').readlines()
            comlist = string.split(string.join(comlist,''),'\r')
            comlist = string.split(string.join(comlist,''),'\n')
            while '' in comlist: comlist.remove('')
            comlist = readDelimit(string.join(comlist,','),',')
        ### Split list ###
        else: comlist = readDelimit(command,',')
        ### Tidy and return ###
        while '' in comlist: comlist.remove('')
        return comlist
    except:
        print 'listFromCommand problem!'
        raise
#########################################################################################################################
def tidyArgs(argcmd,skiplist=['ini'],nopath=False,purgelist=[]):  ### Tidies commandline arguments to remove redundancy
    '''
    Tidies commandline arguments to remove redundancy.
    >> argcmd:list = List of commandline arguments.
    >> skiplist:list = List of commands to skip (duplicates OK).
    >> nopath:bool [False] = Whether to screen out any *path or *dir commands.
    >> purgelist:list [] = List of additional commands to screen out.
    '''
    ### ~ [1] ~ Replace with file aliases and remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    jobargs = []
    jobopts = []
    argcmd = argcmd[0:]
    while argcmd:
        while argcmd.count(argcmd[0]) > 1: argcmd.pop(0)
        #cmd = self.cmdAlias(argcmd.pop(0)) #!# Add general cmdAlias function from slimsuiteREST at some point?
        cmd = argcmd.pop(0)
        jobargs.append(cmd)
        jobopts.append(string.split(cmd,'=')[0])
    ### ~ [2] ~ Reduce to latest version of commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    argnr = []
    while jobargs:
        opt = jobopts.pop(0); cmd = jobargs.pop(0)
        if nopath and (opt.endswith('path') or opt.endswith('dir')): continue
        elif opt in purgelist: continue
        if opt in skiplist or opt not in jobopts:
            if nopath and '=/' in cmd:
                (arg,val) = string.split(cmd,'=',1)
                val = os.path.basename(val)
                argnr.append('%s=%s' % (arg,val))
            else: argnr.append(cmd)
    return argnr
#########################################################################################################################
def argString(arglist):    ### Returns correctly formatted string of commandline arguments
    '''Returns correctly formatted string of commandline arguments.'''
    argstr = []
    for cmd in arglist:
        if len(string.split(cmd)) > 1 and len(string.split(cmd,'=',1)) == 2:
            (opt,val) = string.split(cmd,'=',1)
            argstr.append('%s="%s"' % (opt,val))
        else: argstr.append(cmd)
    return string.join(argstr)
#########################################################################################################################
###  End of Input Command Functions                                                                                     #
#########################################################################################################################

#########################################################################################################################
###  Generic Error Handling Functions                                                                                   #
#########################################################################################################################
def errorReport(text='',quitchoice=False,killme=True):
    '''
    Raises text as error and prints to log.
    >> text:str = Optional error description text to print.
    >> quitchoice:bool = whether to give user choice to terminate program prematurely.
    '''
    ### Handle SystemExit and KeyboardInterrupt from previous sources
    try: raise
    except SystemExit: os._exit(1)
    except KeyboardInterrupt: quitchoice = True
    except: pass
    if string.split(str(sys.exc_info()[0]),'.')[1] == 'SystemExit': os._exit(1)

    try:
        ### Setup error variables
        error_type = str(sys.exc_info()[0])         # Error Type       : exceptions.IOError
        error_type = string.replace(error_type,'exceptions.','')
        error_value = str(sys.exc_info()[1])        # Error Value      : [Errno 9] Bad file descriptor
        error_traceback = traceback.extract_tb(sys.exc_info()[2])
        error_file = str(error_traceback[-1][0])    # File             : C:\Documents and Settings\normdavey\Desktop\Python\BLAST\Main.py
        error_method = str(error_traceback[-1][2])  # Method           : readFile
        error_line = str(error_traceback[-1][1])    # Line             : 15
        error_error = str(error_traceback[-1][3])   # Error            : for lines in fileIn.readlines():
        ### Output
        print '### ----- ERROR! ----- ###'
        if text: print text
        print '%s %s' % (sys.exc_info()[0],sys.exc_info()[1])
        print '%s: %s' % (error_type, error_value)
        print 'File: %s' % error_file
        print 'Method: %s (line %s)' % (error_method, error_line)
        print 'Error: %s' % error_error
        print '### ------------------ ###'
        ### Kill program
        if quitchoice:
            print 'Quit program (y/n)?',
            if raw_input().upper().startswith('Y'): os._exit(1)
        elif killme: os._exit(1)
    except:
        print 'That\'s not right: Error in errorReport()!!:'
#########################################################################################################################
###  End of Generic Error Handling Functions                                                                            #
#########################################################################################################################

