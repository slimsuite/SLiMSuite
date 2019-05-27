#!/usr/bin/python

# rje_obj.py - RJE General object module
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_obj
Description:  Contains revised General Object templates for Rich Edwards scripts and bioinformatics programs
Version:      2.4.1
Last Edit:    22/05/19
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains updated Classes for use by all RJE scripts and bioinformatics programs. With time, all
    programs should be migrated over to the new objects. Miscellaneous methods will, for the most part, be kept in
    the rje module.

    Commandline options are all in the form X=Y. Where Y is to include spaces, use X="Y". Commands in the form "-X Y"
    will also be recognised.

General Commandline:
    v=X             : Sets verbosity (-1 for silent, 0 for no progress counters) [1]
    i=X             : Sets interactivity (-1 for full auto) [0]
    log=FILE        : Redirect log to FILE [Default = calling_program.log]
    newlog=T/F      : Create new log file. [Default = False: append log file]
    silent=T/F      : If set to True will not write to screen or log. [False]
    quiet=T/F       : If set to True will not write to screen or log except for errors. [False]
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
    memsaver=T/F    : Some modules will have a memsaver option to save memory usage [False]

System Commandline:
    win32=T/F       : Run in Win32 Mode [False]
    osx=T/F         : Run in MacOSX Mode [False]
    pwin            : Run in PythonWin (** Must be 'commandline', not in ini file! **)
    runpath=PATH    : Run program from given path (log files and some programs only) [path called from]
    rpath=PATH      : Path to installation of R ['R']
    webserver=T/F   : Trigger webserver run and output [False]
    soaplab=T/F     : Implement special options/defaults for SoapLab implementations [False]
    rest=X          : Variable that sets the output to be returned by REST services [None]
    screenwrap=X    : Maximum width for some screen outputs [200]

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

Uses general modules: glob, math, os, random, re, string, sys, time, traceback
"""
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on rje.py objects.
    # 1.0 - Fully working version, so upgraded to 1.0. Added dev and warn options.
    # 1.1 - Added rje_zen import and self.zen() to call rje_zen.Zen().wisdom().
    # 1.2 - Added warnLog functions.
    # 1.3 - Added perc cmdtype = float that is multiplied by 100.0 if < 1.0. Also added cmdtype = date for YYYY-MM-DD.
    # 1.4 - Added sourceDataFile() method from SLiMBench for wider use.
    # 1.5 - Added 'basefile' cmdlist types.
    # 1.6 - Added osx=T/F option for Mac-specific running options.
    # 1.7 - Added self.name() to basic object class.
    # 1.8 - Cleaned up some erroeneous opt, stat and info references.
    # 2.0 - Added self.file dictionary and methods for handling file handles with matching self.str filenames.
    # 2.1.0 - Added new built-in attributes/options for REST services.
    # 2.1.1 - Removed excess REST HTML methods.
    # 2.1.2 - Tweaked glist cmdRead warnings.
    # 2.1.3 - Modified integer commands to read/convert floats.
    # 2.2.0 - Added screenwrap=X.
    # 2.2.1 - Improved handling of integer parameters when given bad commands.
    # 2.2.2 - Updated error handling for full REST output.
    # 2.3.0 - Added quiet mode to object and stderr output.
    # 2.4.0 - Added vLog() and bugLog() methods.
    # 2.4.1 - Fixed bug where silent=T wasn't running silent.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Get working objects with new str, int, num, bool, dict, list and obj.
    # [ ] : Add typeCheck() method to cycle through each attribute dictionary and check values are the correct types.
    '''
#########################################################################################################################
import glob, os, pickle, random, string, sys, time, traceback
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
import rje, rje_html, rje_zen
#########################################################################################################################

#########################################################################################################################
###  RJE_Object Class:                                                                                                  #
#########################################################################################################################
class RJE_Object(object):     ### Metaclass for inheritance by other classes
    '''Metaclass for inheritance by other classes.'''
    ### ~ [1] General Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    #cmd_list = []   # List of commandline parameters from which options are parsed.
    #log = None      # Log Object (if any) for handling errors etc.
    ### ~ [2] Basic Attribute Lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    #str = {}        # Stores string variables, such as names and decriptions
    #int = {}        # Stores numeric variables, constrained to integers
    #num = {}        # Stores numeric variables that can be floats (integers allowed)
    #bool = {}       # Stores boolean variables
    #file = {}       # Stores file handle objects with matching self.str filenames
    #obj = {}        # Stores a dictionary of other RJE_Objects 'owned' by object
    #list = {}       # Dictionary of list-type attributes for object
    #dict = {}       # Dictionary of dictionary attributes for object
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def __init__(self,log=None,cmd_list=[],parent=None,newstyle=True):  
        '''
        RJE_Object:
        > log:Log = rje.Log object
        > cmd_list:List = List of commandline variables
        > parent:Object [None] = Optional parent object. Can be used for propagating changes back "up".
        > newstyle:bool [True] = Whether to use the new-style attribute dictionaries.

        On intiation, this object:
        - sets the Log object (can be None)
        - sets verbosity and interactive attributes
        - calls the _setAttributes() method to setup class attributes
        - calls the _cmdList() method to process relevant Commandline Parameters       
        '''
        ### ~ [1] Setup Log Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.warnings = []              # Stores a list of warnings that have already been triggered.
        self.no_suppression = []        # Stores a list of warnings that could be suppressed but should not be.
        self.cmd_list = cmd_list
        if not log: self.log = Log(cmd_list=cmd_list)
        else: self.log = log
        ### ~ [2] Setup Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.old = False
        self._setGeneralAttributes(parent)        # Sets general attributes used (potentially) by all objects
        self._setAttributes()               # Specific method in other classes that set own attributes
        ### ~ [3] Commandline Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._cmdList()                     # Read attribute settings from user-defined parameters
        if not newstyle: self.backConvert() 
        if self.getBool('Silent') or self.getBool('Quiet'):
            self._cmdRead('v=-1',type='int',att='Verbose',arg='v')
            self._cmdRead('i=-1',type='int',att='Interactive',arg='i')
#########################################################################################################################
    def _setGeneralAttributes(self,parent=None):    ### Sets general attributes for use in all classes
        '''Sets general attributes for use in all classes.'''
        self.str = {'Name':'None'}
        self.int = {}       # Stores numeric variables, constrained to integers
        self.num = {}       # Stores numeric variables that can be floats (integers allowed)
        self.bool = {}      # Stores boolean variables
        self.file = {}      # Stores open file handles with matching self.str filenames
        self.obj = {'Parent':parent}    # Stores a dictionary of other RJE_Objects 'owned' by object
        self.list = {}      # Dictionary of list-type attributes for object
        self.dict = {}      # Dictionary of list-type attributes for object
        if parent: return   ### ObjectLite settings - no need for object Lite this way!
        self.str = {'Name':'None','Basefile':'None','Delimit':rje.getDelimit(self.cmd_list),'ErrorLog':'None',
                    'RunPath':rje.makePath(os.path.abspath(os.curdir)),'RPath':'R','Rest':'None'}
        self.str['Path'] = rje.makePath(os.path.abspath(string.join(string.split(sys.argv[0],os.sep)[:-1]+[''],os.sep)))
        self.int = {'Verbose':1,'Interactive':0,'ScreenWrap':200}
        self.bool = {'DeBug':False,'Win32':False,'PWin':False,'MemSaver':False,'Append':False,'MySQL':False,
                     'Force':False,'Pickle':True,'SoapLab':False,'Test':False,'Backups':True,'Silent':False,'Quiet':False,
                     'Webserver':False,'ProgLog':True,'Warn':True,'Dev':False,'Setup':False,'OSX':False}
        self.dict = {'Output':{}}
        self.obj['DB'] = None
#########################################################################################################################
    def _setForkAttributes(self):       ### Sets forking attributes for use in all classes
        '''Sets general forking attributes for use in all classes.'''
        self.int['Forks'] = 0
        self.num['KillForks'] = 36000
        self.bool['NoForks'] = False
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object. Should be replaced in other objects. (See RJE_Object.)
        '''
        Sets Attributes of Object:
        - Str, Int, Num, Bool, List, Dict and Obj
        '''
        self.strlist = []; self.intlist = []; self.numlist = []; self.boollist = []
        self.objlist = []; self.listlist = []; self.dictlist = []; self.filelist = []
#########################################################################################################################
    def _setDefaults(self,str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=False):  ### Default defaults!
        '''Sets default defaults.'''
        for i in self.strlist: self.str[i] = str
        for o in self.boollist: self.bool[o] = bool
        for s in self.intlist: self.int[s] = int
        for s in self.numlist: self.num[s] = num; 
        for j in self.objlist: self.obj[j] = obj
        if setlist:
            for l in self.listlist: self.list[l] = []
        if setdict:
            for d in self.dictlist: self.dict[d] = {}
        if setfile:
            for f in self.filelist: self.file[f] = None
#########################################################################################################################
    def setLog(self,log,cascade=True):  ### Sets given log as log object and cascades through self.obj
        '''Sets given log as log object and cascades through self.obj.'''
        self.log = log
        for obj in self.obj.values():
            if cascade and obj and log != obj.log: obj.setLog(log)
#########################################################################################################################
    def warnChecks(self):   ### Checks certain input paths etc. and warns of anomalies
        '''Checks certain input paths etc. and warns of anomalies.'''
        if self.warn(): pass
#########################################################################################################################
    def warnLog(self,message,warntype=None,quitchoice=False,suppress=False,dev=False,screen=True):
        return self.log.warnLog(message,warntype,quitchoice,suppress,dev,screen)
    def infoLog(self,message):
        if self.v() < 1 or self.getBool('Silent'): return message
        return self.log.printLog('#INFO',message)
#########################################################################################################################
    ### <2> ### General Attribute Get/Set methods                                                                       #
#########################################################################################################################
    def prog(self): return self.log.name()
    def type(self): return string.split('%s' % self)[0][1:]
    def name(self): return self.getStr('Name',default='%s' % self)
#########################################################################################################################
    def parent(self): return self.obj['Parent']
    def baseFile(self,newbase=None,return_none='None',runpath=False,strip_path=False): return self.basefile(newbase,return_none,runpath,strip_path)
    def basefile(self,newbase=None,return_none='None',runpath=False,strip_path=False):
        if newbase: return self.setBasefile(newbase)
        elif not self.getStrLC('Basefile'):
            if runpath: return self.getStr('RunPath')
            return return_none
        fullpath = os.path.abspath(self.getStr('Basefile')) == self.getStr('Basefile')
        if runpath and not fullpath and not self.getStr('Basefile').startswith(self.getStr('RunPath')):
            return '%s%s' % (self.getStr('RunPath'),self.getStr('Basefile'))    #!# Check that this does not break stuff!
        if strip_path: return rje.stripPath(self.getStr('Basefile'))
        else: return self.getStr('Basefile')
    def setBaseFile(self,basefile=None,cascade=True): return self.setBasefile(basefile,cascade)
    def setBasefile(self,basefile=None,cascade=True): ### Sets basefile and cascades to daughter objects
        '''Sets basefile and cascades to daughter objects.'''
        if basefile: self.setStr({'Basefile':basefile})
        else: basefile = self.basefile(runpath=False)
        for obj in self.obj.values():
            try:
                if obj and cascade and obj.basefile(runpath=False) != basefile: obj.setBasefile(basefile)
            except: pass
#########################################################################################################################
    def interactive(self): return self.getInt('Interactive')
    def i(self): return self.getInt('Interactive')
    def v(self): return self.getInt('Verbose')
    def force(self): return self.getBool('Force')
    def test(self): return self.getBool('Test')
    def server(self): return self.getBool('Webserver')
    def warn(self): return self.getBool('Warn') and not self.getBool('Silent')
    def dev(self): return self.getBool('Dev')
    def zen(self): return rje_zen.Zen().wisdom()
    def debugging(self): return self.getBool('DeBug')
    def win32(self): return self.getBool('Win32')
    def osx(self): return self.getBool('OSX')
    def silence(self): return self.log.silence()
    def quiet(self): return self.log.quiet()
    def talk(self): return self.log.talk()
#########################################################################################################################
    def getStrLC(self,key,return_none=False):
        if not key: return ''
        getstr = self.getStr(key).lower()
        if getstr == 'none' and not return_none: return ''
        return getstr
#########################################################################################################################
    def getStrUC(self,key,return_none=False):
        if not key: return ''
        getstr = self.getStr(key).upper()
        if getstr == 'NONE' and not return_none: return ''
        return getstr
#########################################################################################################################
    def getStr(self,key=None,default='',checkdata=False):    ### Returns string attribute
        '''Returns string attribute.'''
        try:
            if not key: return self.str
            if key in self.str: return self.str[key]
        except:
            if not key: return self.info
            if key in self.info: return self.info[key]
        if checkdata and 'Data' in self.dict and key in self.dict['Data']: return self.dict['Data'][key]
        elif self.obj['Parent']: return self.obj['Parent'].getStr(key,default,checkdata)
        else: return default
#########################################################################################################################
    def getStrBase(self,key=None,default='',checkdata=False,strip_path=False,extlist=[]):   ### Returns file without extension, with or without path):    ### Returns string attribute
        '''Returns basefile for string attribute.'''
        return rje.baseFile(self.getStr(key,default,checkdata),strip_path,extlist)
#########################################################################################################################
    def getInfo(self,key=None,default='',checkdata=False,warn=True):   ### Returns string attribute for backwards compatibility
        '''Returns string attribute for backwards compatibility.'''
        if self.dev() and warn: self.printLog('#DEV','Code calling getInfo() for new %s object' % self)
        return self.getStr(key,default,checkdata)
#########################################################################################################################
    def getInt(self,key=None,default=0,checkdata=False):    ### Returns integer attribute
        '''Returns integer attribute.'''
        try:
            if not key: return self.int
            if key in self.int: return self.int[key]
            elif key in self.num: return int(self.num[key])
        except:
            if not key: return self.stat
            if key in self.stat: return self.stat[key]
        if checkdata and 'Data' in self.dict and key in self.dict['Data']: return int(self.dict['Data'][key])
        elif self.obj['Parent']: return self.obj['Parent'].getInt(key,default,checkdata)
        else: return default
#########################################################################################################################
    def getPerc(self,key=None,default=0.0,checkdata=False): return self.getNum(key,default,checkdata)/100.0
    def getNum(self,key=None,default=0.0,checkdata=False):    ### Returns float attribute
        '''Returns float attribute.'''
        try:
            if not key: return self.num
            if key in self.num: return self.num[key]
            elif key in self.int: return self.int[key]
        except:
            if not key: return self.stat
            if key in self.stat: return self.stat[key]
        if checkdata and 'Data' in self.dict and key in self.dict['Data']: return float(self.dict['Data'][key])
        elif self.obj['Parent']: return self.obj['Parent'].getNum(key,default,checkdata)
        else: return default
#########################################################################################################################
    def getStat(self,key=None,default=0.0,checkdata=False,warn=True):   ### Returns numerical attribute for backwards compatibility
        '''Returns numerical attribute for backwards compatibility.'''
        if self.dev() and warn: self.printLog('#DEV','Code calling getStat() for new %s object' % self)
        if not key:
            try: return rje.combineDict(rje.combineDict({},self.num),self.int)
            except: return self.stat
        return self.getNum(key,default,checkdata)
#########################################################################################################################
    def getBool(self,key=None,default=False,checkdata=False):    ### Returns boolean attribute
        '''Returns boolean attribute.'''
        try:
            if not key: return self.bool
            if key in self.bool: return self.bool[key]
        except:
            if not key: return self.opt
            if key in self.opt: return self.opt[key]
        if checkdata and 'Data' in self.dict and key in self.dict['Data']:
            if self.dict['Data'][key] == 'False': return False
            if self.dict['Data'][key]: return True
            else: return False
        elif self.obj['Parent']:
            try: return self.obj['Parent'].getBool(key,default,checkdata)
            except: return self.obj['Parent'].getBool(key,default)
        else: return default
#########################################################################################################################
    def getOpt(self,key=None,default=False,checkdata=False,warn=True):   ### Returns boolean attribute for backwards compatibility
        '''Returns boolean attribute for backwards compatibility.'''
        if self.dev() and warn: self.printLog('#DEV','Code calling getOpt() for new %s object' % self)
        return self.getBool(key,default,checkdata)
#########################################################################################################################
    def yesNo(self,text='',default='Y',confirm=False,i=0,log=True):
        if default in [True,False]: default = {True:'Y',False:'N'}
        if self.i() < i:
            if log: self.printLog('#AUTO','%s: %s' % (text,default.upper()))
            return {'Y':True,'N':False}[default.upper()]
        else: return rje.yesNo(text,default,confirm)
#########################################################################################################################
    def choice(self,text='',default='Y',confirm=False,i=0):
        if self.i() < i:
            self.printLog('#AUTO','%s: %s' % (text,default))
            return default
        else: return rje.choice(text,default,confirm)
#########################################################################################################################
    def db(self,table=None,add=False,forcecheck=True,mainkeys=[],uselower=False):   ### Return rje_dbase Database object or Table
        '''
        Return rje_dbase Database object or Table, if one associated with Object.
        >> table:str or Table [None] = Table to look for and return from Database Object 
        >> add:bool [False] = Whether to try and add table with minimal info (first field as key) if missing
        >> forcecheck:bool [True] = Whether to only add missing table if not self.force()
        >> mainkeys:list [] = Keys to use if adding table.
        >> uselower:bool [False] = Whether to convert headers into lower case.
        '''
        try:### ~ [1] ~ Return Database object if not table given ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB']
            if not table: return db
            ### ~ [2] ~ Return Table if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tdb = db.getTable(table)
            if tdb or not add: return tdb
            ### ~ [3] ~ Add Table if missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if forcecheck and self.force():
                if rje.exists(db.dbFileName(table)): self.infoLog('Ignoring %s: force=T' % db.dbFileName(table))
                return None
            #self.printLog('#TAB','Table "%s" not found' % table)
            return db.addTable(mainkeys=mainkeys,name=table,expect=False,uselower=uselower)
        except: return None
#########################################################################################################################
    def data(self,table,strict=False):
        try: return self.obj['DB'].getTable(table).data()
        except:
            if strict: self.errorLog('No DB table "%s"?' % table); raise
            else: return {}
#########################################################################################################################
    def getAtt(self,type,key,default=None): return self.getAttribute(type,key,default)
    def getAttribute(self,type,key,default=None):    ### Gets object information of correct type
        '''
        Gets object information of correct type.
        >> type:str = 'info','list','dict','opt','stat' or 'obj'
        >> key:str = key for given attribute dictionary
        >> default:anything = what to return if type is wrong or key not found in type
        '''
        att = {}
        if type in ['info','str']: return self.getStr(key,default)
        elif type in ['opt','bool']: return self.getBool(key,default)
        elif type == 'stat': return self.getStat(key,default)
        elif type == 'int': return self.getInt(key,default)
        elif type == 'num': return self.getNum(key,default)
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
     ### <3> ### Command-line parameters                                                                                 #
#########################################################################################################################
    def setStr(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if self.old: return self.setInfo(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.str[key] = attdict[key]
                if addtolist and key not in self.strlist: self.strlist.append(key)
        except: self.errorLog('Problem with %s.setStr()' % self,True)
#########################################################################################################################
    def setBool(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if self.old: return self.setOpt(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.bool[key] = attdict[key]
                if addtolist and key not in self.boollist: self.boollist.append(key)
        except: self.errorLog('Problem with %s.setBool()' % self,True)
#########################################################################################################################
    def setInt(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if self.old: return self.setStat(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.int[key] = attdict[key]
                if addtolist and key not in self.intlist: self.intlist.append(key)
        except: self.errorLog('Problem with %s.setInt()' % self,True)
#########################################################################################################################
    def setNum(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if self.old: return self.setStat(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.num[key] = attdict[key]
                if addtolist and key not in self.numlist: self.numlist.append(key)
        except: self.errorLog('Problem with %s.setNum()' % self,True)
#########################################################################################################################
    def setObj(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        try:
            for key in attdict.keys():
                self.obj[key] = attdict[key]
                if addtolist and key not in self.objlist: self.objlist.append(key)
        except: self.errorLog('Problem with %s.setObj()' % self,True)
#########################################################################################################################
    def setList(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        try:
            for key in attdict.keys():
                self.list[key] = attdict[key]
                if addtolist and key not in self.listlist: self.listlist.append(key)
        except: self.errorLog('Problem with %s.setList()' % self,True)
#########################################################################################################################
    def setDict(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        try:
            for key in attdict.keys():
                self.dict[key] = attdict[key]
                if addtolist and key not in self.dictlist: self.dictlist.append(key)
        except: self.errorLog('Problem with %s.setDict()' % self,True)
#########################################################################################################################
    def _generalCmd(self,cmd=''):     ### Sets General Attributes from commandline
        '''
        Sets general attributes according to commandline parameters:
        - see rje.__doc__ or run with 'help' option
        '''
        try:
            self._cmdRead(cmd,type='int',att='Verbose',arg='v')
            self._cmdRead(cmd,type='int',att='Interactive',arg='i')
            self._cmdReadList(cmd,type='int',attlist=['Verbose','Interactive','ScreenWrap'])
            self._cmdReadList(cmd,type='num',attlist=['MaxBin'])
            self._cmdRead(cmd,type='file',att='Basefile',arg='outfile')
            self._cmdRead(cmd,type='str',att='Rest',arg='outfmt')
            self._cmdReadList(cmd,'str',['Rest'])
            self._cmdReadList(cmd,'file',['Basefile','RPath','ErrorLog'])
            self._cmdReadList(cmd,'abspath',['Path','RunPath'])
            self._cmdRead(cmd,type='bool',att='Win32',arg='pwin')
            self._cmdReadList(cmd,'bool',['DeBug','Win32','PWin','MemSaver','Append','Force','MySQL','Pickle','Test',
                                         'SoapLab','Backups','Webserver','ProgLog','Dev','Warn','OSX'])
        except:
            self.deBug(self.cmd_list)
            self.errorLog('Problem with %s.cmd:%s' % (self,cmd))
#########################################################################################################################
    def _cmdReadList(self,cmd=None,type='str',attlist=[]):     ### Sets self.type[att] from commandline command cmd
        '''
        Sets self.type[att] from commandline command cmd.
        >> cmd:str = commandline command
        >> type:str = type of attribute (info,path,opt,int,float,min,max,list,clist,glist,file)
        >> att:list of attributes (key of dictionary) where commandline argument is att.lower()  []
        '''
        for att in attlist:
            if type.startswith('base'): self._cmdRead(cmd,type[4:],att,'basefile')
            else: self._cmdRead(cmd,type,att)
#########################################################################################################################
    def _cmdRead(self,cmd=None,type='str',att=None,arg=None):     ### Sets self.type[att] from commandline command cmd
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
        if cmdarg != arg: return
        value = cmd[len('%s=' % arg):]
        value = string.replace(value,'#DATE',rje.dateTime(dateonly=True))
        ### ~ [1] Basic commandline types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if type in ['opt','bool']:
            if value[:1].lower() in ['f','0']: self.bool[att] = False
            else: self.bool[att] = True
        elif type in ['info','str']:
            if arg == 'basefile': self.str[att] = rje.baseFile(value,True)
            else: self.str[att] = value
        elif type == 'path': self.str[att] = rje.makePath(os.path.expanduser(value))
        elif type == 'abspath': self.str[att] = rje.makePath(os.path.abspath(os.path.expanduser(value)))
        elif type in ['fullpath','file']: self.str[att] = rje.makePath(os.path.expanduser(value),wholepath=True)
        elif type == 'int':
            try: self.int[att] = string.atoi(value)
            except:
                if rje.matchExp('^(\d+)',value): self.int[att] = string.atoi(rje.matchExp('^(\d+)',value)[0])
                else: self.int[att] = int(string.atof(value))
                self.warnLog('%s=%s needs integer -> %s=%d' % (arg,value,arg,self.int[att]))
        elif type in ['float','stat','num']: self.num[att] = string.atof(value)
        elif type == 'max' and rje.matchExp('^\d+,(\d+)',value): self.int[att] = string.atoi(rje.matchExp('^\d+,(\d+)',value)[0])
        elif type in ['min','max'] and rje.matchExp('^(\d+)',value): self.int[att] = string.atoi(rje.matchExp('^(\d+)',value)[0])
        elif type == 'fmax' and rje.matchExp('^[\.\d]+,([\.\d]+)',value): self.num[att] = string.atof(rje.matchExp('^[\.\d]+,([\.\d]+)',value)[0])
        elif type in ['fmin','fmax'] and rje.matchExp('^([\.\d]+)',value): self.num[att] = string.atof(rje.matchExp('^([\.\d]+)',value)[0])
        elif type == 'date' and rje.matchExp('^(\d\d\d\d)-?(\d\d)-?(\d\d)$',value):
            self.str[att] = '%s-%s-%s' % rje.matchExp('^(\d\d\d\d)-?(\d\d)-?(\d\d)$',value)
        elif type == 'list': self.list[att] = rje.listFromCommand(value)
        elif type == 'lclist': self.list[att] = rje.listLower(rje.listFromCommand(value))
        elif type == 'uclist': self.list[att] = rje.listUpper(rje.listFromCommand(value))
        ### ~ [2] Special types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        elif type == 'perc':
            self.num[att] = string.atof(value)
            if self.num[att] < 1.0: self.num[att] *= 100.0
        elif type in ['ilist','nlist']:
            nlist = rje.listFromCommand(value)
            self.list[att] = []
            for nvalue in nlist:
                if type == 'ilist': self.list[att].append(string.atoi(nvalue))
                else: self.list[att].append(string.atof(nvalue))
        elif type == 'clist':   # Returns a *string* of a CSV list
            self.str[att] = string.join(rje.listFromCommand(value),',')
        elif type == 'glist':   # 'Glob' List - returns a list of files using wildcards & glob
            globlist = string.split(value,',')
            if len(globlist) == 1 and globlist[0].endswith('.fofn'):  # File of file names
                globlist = rje.listFromCommand(value)
            self.list[att] = []
            for g in globlist:
                if not g: continue
                newg = glob.glob(g); newg.sort()
                if not newg:
                    if not '*' in g: self.warnLog('No "%s" %s files found' % (g,arg),warntype='noglob',suppress=True,quitchoice=True)
                    else: self.warnLog('No "%s" %s files found' % (g,arg),warntype='noglob')
                self.list[att] += newg
        elif type in ['cdict','cdictlist']:   # Converts a CSV list into a dictionary
            self.dict[att] = {}
            if type == 'cdictlist': self.list[att] = []
            clist = rje.listFromCommand(value)
            for c in clist:
                data = string.split(c,':')
                if ':' in c:
                    if len(data) > 2: data = [data[0],string.join(data[1:],':')]
                    self.dict[att][data[0]] = data[1]
                    if type == 'cdictlist': self.list[att].append(data[0])
        else:
            self.deBug('%s=%s' % (type,value))
            raise ValueError
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try: self._generalCmd(cmd)
            except: self.errorLog('Problem with %s.cmd:%s' % (self,cmd))
#########################################################################################################################
    def _forkCmd(self,cmd=''):     ### Sets General Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see rje.__doc__ or run with 'help' option
        '''
        try:
            self._cmdReadList(cmd,'int',['Forks','KillForks'])
            self._cmdRead(cmd,type='bool',att='NoForks')
            self._cmdRead(cmd,type='bool',att='NoForks',arg='nofork')
        except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
     ### <4> ### Input/Output                                                                                            #
#########################################################################################################################
    def bugProg(self, id='#ERR', text='Log Text Missing!',screen=True,rand=0.0,clear=0):
        return self.progLog(id, text,screen and self.debugging(),rand,clear)
    def progLog(self, id='#ERR', text='Log Text Missing!',screen=True,rand=0.0,clear=0):
        if self.v() < 1: return
        if rand > 0 and random.random() > rand: return
        if self.getBool('ProgLog',default=True): return self.printLog('\r%s' % id,text,screen=screen,log=False,newline=False,clear=clear)
        else: return False
    def printLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True, clear=0):
        return self.log.printLog(id,text,timeout,screen and not self.getBool('Silent'),log and not self.getBool('Silent'),newline,clear=clear)
    def bugLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True, clear=0):
        return self.log.printLog(id,text,timeout,screen and not self.getBool('Silent') and self.debugging(),log and not self.getBool('Silent') and self.debugging(),newline,clear=clear)
    def vLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True, clear=0,v=1):
        return self.log.printLog(id,text,timeout,screen and not self.getBool('Silent') and self.v() >= v,log and not self.getBool('Silent') and self.v() >= v,newline,clear=clear)
    def errorLog(self, text='Missing text for errorLog() call!',quitchoice=False,printerror=True,nextline=True,log=True,errorlog=True,warnlist=True):
        return self.log.errorLog(text,quitchoice,printerror,nextline,log,errorlog,warnlist)
    def devLog(self, lid='#DEV', text='Log Text Missing!', debug=True):
        if not self.dev(): return
        self.printLog(lid,text)
        if debug: self.debug('')
#########################################################################################################################
    def headLog(self,text,line='~',hash='#',width=0,minside=4):  # Generates a header-style printLog command
        '''
        Generates a header-style printLog command.
        >> text: header text
        >> line: str for lines
        >> hash:str [#] = borders for main text
        >> width: total width for the printed message
        >> minside: shortest
        '''
        if not width and line in '~-=': width = {'=':90,'-':75,'~':60}[line]
        strlist = [hash,line*minside,text,line*minside,hash]
        while len(string.join(strlist)) < width:
            strlist[1] += line
            strlist[3] += line
        if len(string.join(strlist)) == width + 1 and len(strlist[1]) > minside: strlist[1] = strlist[1][:-1]
        self.printLog('%s%s%s%s' % (hash,line,line,hash),string.join(strlist))
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
        if i == None: i = self.getInt('Interactive') + 1
        if not self.getBool('Silent') and (self.getInt('Verbose') >= v or self.getInt('Interactive') >= i):
            print text,
            if self.getInt('Interactive') >= i:
                raw_input(" <ENTER> to continue.")
                if 'pwin' not in sys.argv + self.cmd_list: newline -= 1
            while newline > 0: print; newline -= 1
            try: sys.stdout.flush()
            except: pass
#########################################################################################################################
    def screenWrap(self,intext,screenwrap=0,prefix='',suffix='',fixwidth=False,wordsplit=True):   ### Returns wrapped text
        '''
        Returns wrapped text, framed by prefix and suffix.
        @param intext:str Text to wrap.
        @param screenwrap:int [0] width at which to wrap. Will use self.getInt('ScreenWrap') if < 1.
        @param prefix:str [''] Text to place at start of each line.
        @param suffix:str  [''] Text to place at end of each line.
        @param fixwidth:bool [False] Whether to pad wraptext with spaces to always hit fixwidth
        @param wordsplit:bool [True] Whether to wrap at whitespace [True] or at any point [False]
        @return: wrapped text
        '''
        if screenwrap < 1: screenwrap = self.getInt('ScreenWrap')
        wrapped = []
        if wordsplit: wraptext = string.split(intext)
        else: wraptext = rje.strList(intext)
        splitlen = screenwrap - len(prefix) - len(suffix)
        while wraptext:
            addtext = wraptext.pop(0)
            while wraptext and (len(addtext) + len(wraptext[0])) <= splitlen:
                if wordsplit: addtext += ' '
                addtext += wraptext.pop(0)
            if fixwidth and len(addtext) < splitlen: addtext += ' ' * (splitlen - len(addtext))
            wrapped.append('%s%s%s' % (prefix,addtext,suffix))
        while intext.endswith('\n') and len(wrapped) < 2: wrapped.append('')
        return string.join(wrapped,'\n')
#########################################################################################################################
    def debug(self,text): self.deBug(text)
#########################################################################################################################
    def bugPrint(self,text): self.deBug(text,pause=False)
#########################################################################################################################
    def deBug(self,text=None,pause=True):   ### Prints text to screen if self.debug 
        '''
        Prints text to screen if self.bool['DeBug'].
        >> text:str = Debugging text to print.
        '''
        if text == None: return self.getBool('DeBug')
        pause = pause and self.i() >= 0
        try:
            if self.getBool('DeBug') and pause: self.verbose(self.v(),self.i(),text,1)
            elif self.getBool('DeBug'): self.verbose(self.v(),self.i()+1,text,1)
        except KeyboardInterrupt:
            if self.yesNo('Interrupt program?'): raise
            if self.yesNo('Switch off Debugging?'): self.bool['DeBug'] = False; self.cmd_list.append('debug=F')
#########################################################################################################################
    def pickleMe(self,basefile=None,gzip=True,replace=True):   ### Saves self object to pickle and zips
        '''
        Saves self object to pickle and zips.
        >> basefile:str [None] = if none, will use self.info['Basefile']
        >> gzip:bool [True] = whether to GZIP (win32=F only)
        >> replace:bool [True] = whether to replace existing Pickle
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Pickle'): self.printLog('#PICKLE','%s pickling disabled. (pickle=F)' % self.prog()); return False
            if not basefile or basefile.lower() == 'none': basefile = self.getStr('Basefile')
            pfile = '%s.pickle' % basefile
            if not replace and (os.path.exists(pfile) or os.path.exists('%s.gz' % pfile)):
                self.printLog('#PICKLE','No pickling - pickle already exists (no replace)!'); return False
            self.close()    # Cannot pickle file handles.
            ### ~ [2] ~ Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#SAVE','Attempting to save %s to %s.' % (self.prog(),pfile))
            pickle.dump(self,open(pfile,'w'))
            self.printLog('\r#SAVE','%s Intermediate saved as %s (Python pickle).' % (self.prog(),pfile))
            ### ~ [3] ~ GZip and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Win32') and gzip:
                try:
                    if os.path.exists('%s.gz' % pfile): os.unlink('%s.gz' % pfile)
                    os.system('gzip %s' % pfile)
                    self.printLog('#GZIP','%s zipped.' % pfile)
                except: self.errorLog('Cannot gzip %s' % pfile)
            return True
        except: self.errorLog('Problem during %s.pickleMe()' % self); return False
#########################################################################################################################
    def unpickleMe(self,basefile=None,process=True): ### (Unzips and) loads pickle object and returns.
        '''(Unzips and) loads pickle object and returns.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basefile or basefile.lower() == 'none': basefile = self.getStr('Basefile')
            pfile = '%s.pickle' % basefile
            gzfile = '%s.gz' % pfile
            if not (os.path.exists(pfile) or os.path.exists(gzfile)): return None
            if not self.getBool('Pickle'): self.printLog('#PICKLE','No loading of pickle. (pickle=F)'); return None
            newme = None    ## New object, loaded from Pickle
            ### ~ [2] ~ Check for file and load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(gzfile) and rje.isYounger(gzfile,pfile) != pfile:
                if self.getBool('Win32'): self.errorLog('Cannot unzip %s (win32=T)' % (gzfile)); return None
                else:
                    if os.path.exists(pfile): os.unlink(pfile)
                    try:
                        self.progLog('\r#GUNZIP','Unzipping %s...' % gzfile)
                        gu = os.popen('gunzip %s' % gzfile).read()
                        self.printLog('\r#GUNZIP','Pickle %s unzipped.' % gzfile)
                    except: self.errorLog('Cannot unzip %s' % (gzfile)); return None
            if os.path.exists(pfile):
                self.printLog('\r#LOAD','Attempting to load %s.' % pfile,log=False)
                newme = pickle.load(open(pfile,'r'))
                self.printLog('\r#LOAD','%s Intermediate loaded: %s.' % (self.prog(),pfile))
                if not self.getBool('Win32'):
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
            if 'Silent' not in newme.bool: newme.bool['Silent'] = False
            for obj in self.obj.values():
                if not obj: continue
                if obj and 'Silent' not in obj.bool: obj.bool['Silent'] = False
                for obj2 in obj.obj.values():
                    if not obj2: continue
                    if obj2 and 'Silent' not in obj2.bool: obj2.bool['Silent'] = False
                    for obj3 in obj2.obj.values():
                        if obj3 and 'Silent' not in obj3.bool: obj3.bool['Silent'] = False
            return newme
        except: self.errorLog('Problem during %s.processPickle()' % self); return None
#########################################################################################################################
    def replaceMe(self,newme):  ### Replaces all my attributes with those of newme
        '''Replaces all my attributes with those of newme.'''
        self.str = newme.str; self.bool = newme.bool; self.int = newme.int; self.num = newme.num
        self.list = newme.list; self.dict = newme.dict; self.obj = newme.obj
#########################################################################################################################
    def checkInputFiles(self,ilist,ask=True):    ### Checks for existence of Input Files and asks for them if missing
        '''
        Checks for existence of Input Files and asks for them if missing. Raises error if not interactive and missing.
        >> ilist:list of keys for self.str that should be sequence files
        >> ask:boolean = whether to ask for file if missing and i >= 0
        '''
        for ifile in ilist:
            while not (self.getStrLC(ifile) and os.path.exists(self.getStr(ifile))):
                if self.getStrLC(ifile):
                    self.errorLog('%s file "%s" does not exist!' % (ifile,self.getStr(ifile)),printerror=False)
                    if not ask or self.i() < 0: raise IOError('%s not found!' % self.getStr(ifile))
                else:
                    self.errorLog('%s=FILE required!' % (ifile),printerror=False)
                    if not ask or self.i() < 0: raise ValueError('%s not given!' % self.getStr(ifile))
                self.str[ifile] = rje.getFileName('%s File?' % ifile)
#########################################################################################################################
    def loadFromFile(self,filename=None,v=0,checkpath=True,chomplines=False):   ### Loads all data from file and returns readlines() list
        '''
        Loads all data from file and returns readlines() list. Will look in same directory and self.str['Path']
        >> filename:str = Name of file
        >> v:int = verbosity setting for loading message
        >> checkpath:boolean [True] = whether to check in self.str['Path'] if filename missing.
        >> chomp:boolean [False] = whether to remove \\r and \\n from lines
        << filelines:list of lines from file
        '''
        try:
            lookinpath = True
            altfile = rje.makePath(self.getStr('Path')) + filename
            if filename.lower() in ['none','']: return []
            if rje.checkForFile(filename): openfile = filename
            elif checkpath and altfile != filename:
                openfile = altfile
                if not rje.checkForFile(openfile):
                    self.errorLog('File %s (and %s) not found!' % (filename,openfile),printerror=False)
                    return []
            else:
                self.errorLog('File %s not found!' % filename,printerror=False)
                return []
            if v <= self.v():
                self.printLog('#LOAD','Loading data from %s ...' % openfile,newline=False,log=False)
            file_lines = open(openfile, 'r').readlines()
            if len(file_lines) == 1: file_lines = string.split(file_lines[0],'\r')
            if chomplines: file_lines = string.split(rje.chomp(string.join(file_lines,'!#ENDOFLINE#!')),'!#ENDOFLINE#!')
            if v <= self.v():
                self.printLog('\r#LOAD','Loading data from %s complete: %s lines.' % (openfile,rje.iLen(file_lines)),log=False)
            return file_lines
        except:
            self.errorLog('Major problem in loadFromFile(%s)' % filename)
            return []
#########################################################################################################################
    def gUnzip(self,file,log=True):  ### Unzips a given file
        '''Unzips a given file.'''
        if log: self.printLog('#UNZIP','Unzipping %s ...' % file,newline=False,log=False)
        try:
            if not os.path.exists(file): raise IOError
            if self.getBool('Win32'): raise ValueError
            os.system('gunzip %s' % file)
            if log: self.printLog('\r#UNZIP','Unzipped %s successfully' % file)
        except: self.errorLog('Problem unzipping %s' % file)
#########################################################################################################################
    def gZip(self,filename,log=True,unlink=True):    ### Zips file, if appropriate
        '''Zips file, if appropriate.'''
        if self.getBool('Win32'): self.printLog('#WIN32','Cannot GZIP %s (Win32=T)' % filename); return False
        if log: self.printLog('#GZIP','Zip %s ...' % filename,newline=False,log=False)
        if unlink and os.path.exists('%s.gz' % filename): os.unlink('%s.gz' % filename)
        if os.path.exists(filename): 
            os.system('gzip %s' % filename)
            self.printLog('#GZIP','%s zipped.' % filename); return True
        else: self.printLog('#ERR','%s missing!' % filename); return False
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
        if rje.isYounger(checkfile,parentfile) != parentfile: return False
        return True
#########################################################################################################################
    def sourceDataFile(self,str,ask=True,expect=True,check=True,force=False,download=None,sourceurl=None):   ### Returns source data file.        #V2.0
        '''
        Returns source data file. This will be stored in the variable self.str[str] and also added to the Source table
        in self.db().
        >> str = Key for self.str dictionary corresponding to the source data file being sought.
        >> ask:bool [True] = Whether to ask for file if not found.
        >> expect:bool [True] = Whether to download/ask for file name if missing
        >> check:bool [True] = Whether to check for file presence. If False will return the desired file download name.
        >> force:bool [True] = Whether to use self.force() to govern file regeneration.
        >> download:bool [None] = Whether to download. If None, will use self.bool['Download'] if possible.
        >> sourceurl:str [None] = String of source URL for download. Will try self.dict['SourceURL'] and 'URL' field
         of self.db('Source') if missing.
        << returns filename if present or None if missing.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if download == None:
                try: download = self.getBool('Download')
                except: self.errorLog('No Download boolean for %s' % self); download = False
            force = force and self.force()
            ask = ask and self.i() >= 0
            if not self.db('Source',add=False): self.db().addEmptyTable('Source',['Name','File','Status','Entries','URL'],keys=['Name'],log=False)   # Store Source info
            sentry = self.db('Source').data(str)
            if not sentry: sentry = self.db('Source').addEntry({'Name':str})
            ## ~ [0a] Filename lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getStrLC('SourceDate'): self.setStr({'SourceDate':''})
            filename = os.path.split(self.getStr(str))[-1]
            #self.deBug(filename)
            sourcefile = self.getStr('SourcePath') + filename
            fileparts = os.path.splitext(filename)
            if fileparts[0].endswith(self.getStr('SourceDate')):
                datefile = sourcefile
            else:
                if rje.matchExp('^(.+)\.(\d\d\d\d-\d\d-\d\d)$',fileparts[0]):
                    fileparts = (rje.matchExp('^(.+)\.(\d\d\d\d-\d\d-\d\d)$',fileparts[0])[0],fileparts[1])
                datefile = '%s%s.%s%s' % (self.getStr('SourcePath'),fileparts[0],self.getStr('SourceDate'),fileparts[1])
            nowfile = '%s%s.%s%s' % (self.getStr('SourcePath'),fileparts[0],rje.dateTime(dateonly=True),fileparts[1])
            if not check: return nowfile
            lastfile = None
            for gfile in glob.glob('%s%s.*%s' % (self.getStr('SourcePath'),fileparts[0],fileparts[1])):
                if rje.matchExp('^%s%s\.(\d\d\d\d)-(\d\d)-(\d\d)%s$' % (self.getStr('SourcePath'),fileparts[0],fileparts[1]),gfile):
                    lastfile = gfile    # Possibly use the latest dated version
            if lastfile in [self.getStr(str), datefile, sourcefile, nowfile]: lastfile = None   # Only interested if different

            ### ~ [1] Return file if it exists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sentry['Status'] = 'Not found'
            checked = [None]
            for checkfile in [self.getStr(str), datefile, sourcefile, nowfile, lastfile]:
                if checkfile in checked: continue
                checked.append(checkfile)
                if not rje.exists(checkfile): continue
                if checkfile != nowfile and force: self.printLog('#FORCE','Ignoring %s. (force=T)' % checkfile); continue
                if checkfile == lastfile and not self.yesNo('%s not found. Use %s?' % (sourcefile,lastfile)): continue
                if checkfile not in [self.getStr(str),datefile,nowfile]:
                    if not self.getStr('SourceDate'):
                        newdate = rje.matchExp('^(\d\d\d\d-\d\d-\d\d)$',lastfile.split('.')[-2])[0]
                        if self.yesNo('Set sourcedate=%s?' % newdate):
                            self.setStr({'SourceDate':newdate})
                            self.printLog('#DATE','sourcedate=%s' % newdate)    # Can parse out to add to ini
                        else:
                            self.setStr({'SourceDate':rje.dateTime(dateonly=True)})
                            self.warnLog('Using %s rather than %s (sourcedate=%s)' % (checkfile,datefile,self.getStr('SourceDate')))
                    else: self.warnLog('Using %s rather than %s (sourcedate=%s)' % (checkfile,datefile,self.getStr('SourceDate')))
                if force and self.yesNo('%s found but force=T. Regenerate?' % checkfile): self.printLog('#FORCE','Ignoring %s. (force=T)' % checkfile); continue
                sentry['Status'] = 'Found'
                sentry['File'] = checkfile
                if checkfile != self.getStr(str): self.printLog('#SOURCE','Set %s=%s.' % (str.lower(),checkfile))
                self.setStr({str:checkfile})
                return checkfile
            if not expect: return None

            ### ~ [2] Optional download ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not sourceurl:
                if 'SourceURL' in self.dict and str in self.dict['SourceURL']:
                    sourceurl = self.dict['SourceURL'][str]
                    sentry['URL'] = sourceurl
                elif 'URL' in sentry: sourceurl = sentry['URL']
            #!# Replace this section with a method that can be edited within the class? #!#
            if download and str in ['HINT.DROME','HINT.CAEEL']:
                self.warnLog('Trying to correct for dodgy HINT.DROME/HINT.CAEEL download w/o headers')
                ## Special fudge
                open(nowfile,'w').write('%s\n' % string.join('Id_A Id_B Gene_A Gene_B Pubmedid,EvidenceCode,HT'.split(),'\t'))
                sourcefile = nowfile
                rje.urlToFile(sourceurl,nowfile,self,backupfile=False)
                sentry['Status'] = 'Downloaded'
            elif download and str == 'TaxMap':
                try:
                    sourcefile = nowfile
                    datecode = rje.dateTime(dateonly=True)
                    taxdump = '%staxdump.%s.tar.gz' % (self.getStr('SourcePath'),datecode)
                    #rje.urlToFile(sourceurl,taxdump,self)
                    #self.debug(sourceurl)
                    if self.force() or not rje.exists(taxdump):
                        self.printLog('#FTP','Downloading %s...' % sourceurl,log=False)
                        if self.getBool('OSX'):
                            os.system("curl -O %s" % (sourceurl))
                            os.rename('taxdump.tar.gz',taxdump)
                        elif self.getBool('Win32'): self.warnLog('Cannnot use wget with Win32=T'); raise ValueError
                        else: os.system("wget -O %s %s" % (taxdump,sourceurl))
                    sentry['Status'] = 'Downloaded TarGZ'
                    predump = glob.glob('*.*')
                    self.printLog('#TARGZ','tar -xzf %s' % taxdump)
                    os.system('tar -xzf %s' % taxdump)
                    for dfile in glob.glob('*.*')[0:]:
                        if dfile in predump: continue
                        fileparts = os.path.splitext(dfile)
                        datefile = '%s%s.%s%s' % (self.getStr('SourcePath'),fileparts[0],datecode,fileparts[1])
                        os.rename(dfile,datefile)
                    sentry['Status'] = 'Downloaded'
                except: self.errorLog('Problem processing NCBI Taxa download')
            elif download and sourceurl:
                sourcefile = nowfile
                rje.urlToFile(sourceurl,nowfile,self)
                sentry['Status'] = 'Downloaded'
            elif self.getBool('Download') and expect: self.warnLog('Downloads for %s not supported' % str)

            ### ~ [3] Ask for replacement file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while sourcefile and ask and not rje.exists(sourcefile):
                btext = {True:'exit',False:'ignore'}[expect and self.getBool('Integrity')]
                sourcefile = rje.choice('%s not found. Enter %s file path (blank to %s)' % (sourcefile,str,btext))
                sentry['Status'] = 'Manual'
            sentry['File'] = sourcefile
            if rje.exists(sourcefile):
                if sourcefile != self.getStr(str): self.printLog('#SET','Set %s=%s.' % (str.lower(),sourcefile))
                self.setStr({str:sourcefile})
                return sourcefile
            elif expect and sourcefile: raise IOError
            elif expect and self.getBool('Integrity'): raise KeyboardInterrupt
            self.warnLog('Problem locating %s source data file' % str,quitchoice=self.getBool('Integrity'))
            sentry['Status'] = 'Missing'
            return ''
        except KeyboardInterrupt: raise
        except: self.errorLog('Problem locating %s source data file' % str,quitchoice=self.getBool('Integrity')); return ''
#########################################################################################################################
    ### <4b> ### REST server output methods. (Should be replaced in mature Classes)                                     #
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
            error_txt.append('%s\n' % rje.dateTime())
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
    def restOutputDefault(self): return 'full'
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
            if outfmt == 'version': return rje.Out(self.log,['v=-1']).printIntro(self.log.obj['Info'])
            if outfmt == 'outfmt': return self.restSetup.__doc__
            if outfmt == 'warnings': return string.join(self.log.list['WarnLog'],'\n')  # List of log warning messages.
            if outfmt == 'errors': return string.join(self.log.list['ErrorLog'],'\n')   # List of log error messages.
            if outfmt in self.dict['Output']:
                outdata = self.dict['Output'][outfmt]
                if outdata in self.str: outdata = self.str[outdata]
                if rje.exists(outdata):
                    nbytes = os.path.getsize(outdata)
                    if nbytes > maxfilesize > 0:
                        otext = '%s is too large to return (> %s)' % (outdata,rje.humanByteSize(nbytes))
                        if outfmt == self.getStrLC('Rest'): return 'ERROR: %s.' % otext
                        return '%s in full output. Try retrieve&rest=%s.' % (otext,outfmt)
                    return rje.fixASCII(open(outdata,'r').read())
                elif outdata: return '%s\n' % outdata
                else: return 'No output generated.\n'
            elif self.db(outfmt):
                ext = rje.delimitExt(self.getStr('Delimit'))
                dbfile = '%s.%s.%s' % (self.baseFile(runpath=True),outfmt,ext)
                if not rje.exists(dbfile): self.db(outfmt).saveToFile(dbfile)
                return os.path.abspath(dbfile)
            # summary = Reduced output
            # Default: return full text
            return '# SERVER ERROR: "&rest=%s" not a recognised output\n%s' % (outfmt,self.restFullOutput())
        except: return self.restOutputError('%s.restOutput() error.' % self.prog())
#########################################################################################################################
    def restFullOutput(self,outfile=None):  ### Full REST server output method. Should be replaced in individual Classes to customise.
        '''
        Full REST server output method. Should be replaced in individual Classes to customise.
        >> outfile:str [None] = Name of output file for REST text. If given, will save full file paths rather than file
        contents for subsequent loading and parsing by SLiMParser.py
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outputs = self.restOutputOrder()
            if not outputs: return 'No full REST output for %s' % self.prog()
            info = self.log.obj['Info']
            outtxt = '# Output for %s V%s: %s\n' % (info.program,info.version,time.asctime(time.localtime(info.start_time)))
            warntxt = '# %s warnings; %s errors' % (rje.iStr(self.log.warnx),rje.iLen(self.log.list['ErrorLog']))
            if self.log.warnx or self.log.list['ErrorLog']: warntxt += ': see log output (below) for details'
            outtxt += '%s.\n' % warntxt
            if 'jobid' in self.dict['Output']: outtxt += '# JobID: %s\n' % self.dict['Output']['jobid']
            if 'jobid' in outputs: outputs.remove('jobid')
            outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            ### ~ [1] RESTful output to return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'status' in outputs: outputs.remove('status')
            outtxt += '# status:\n%s' % self.restOutput('status')
            outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            if 'version' not in outputs:
                outtxt += '# version:\n%s\n' % info.version
                outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
            ## ~ [1a] Output dictionary contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for outkey in outputs:
                try:
                    if outkey in ['help','check']: continue
                    if outkey not in self.dict['Output'] and self.db(outkey):
                        ext = rje.delimitExt(self.getStr('Delimit'))
                        dbfile = '%s.%s.%s' % (self.baseFile(runpath=True),outkey,ext)
                        if not rje.exists(dbfile): self.db(outkey).saveToFile(dbfile)
                        self.dict['Output'][outkey] = os.path.abspath(dbfile)
                    if outkey not in self.dict['Output']:
                        outtxt += '# %s:\nNo output generated.\n' % outkey
                        outtxt += '###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###\n'
                        continue
                    outdata = self.dict['Output'][outkey]
                    keyfile = None
                    try:
                        if self.getStrLC(outdata) and os.path.exists(self.getStr(outdata)):
                            keyfile = self.getStr(outdata)
                            outtxt += '# %s: %s\n' % (outkey,self.getStr(outdata))
                        elif os.path.exists(outdata):
                            keyfile = outdata
                            outtxt += '# %s: %s\n' % (outkey,outdata)
                        else: outtxt += '# %s:\n' % outkey
                    except: outtxt += '# %s:\n' % outkey
                    if outfile and keyfile: outtxt += '%s\n' % os.path.abspath(keyfile)
                    else: outtxt += self.restOutput(outkey)
                except: outtxt += 'ERROR: %s\n' % self.errorLog(outkey)
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
    ### <5> ### Forks                                                                                                   #
#########################################################################################################################
    def _activeForks(self,pidlist=[],nowarn=[0]):   ### Checks Process IDs of list and returns list of those still running.
        '''
        Checks Process IDs of list and returns list of those still running.
        >> pidlist:list of integers = Process IDs
        >> nowarn:list [[0]] = List of exit codes that will not trigger a warning
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            oldpids = pidlist[0:]
            ### ~ [1] Check PIDs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pid in oldpids[0:]:
                (newpid,exitcode) = os.waitpid(pid,os.WNOHANG)
                while exitcode > 255: exitcode -= 256
                if newpid == pid and exitcode in nowarn: oldpids.remove(pid)
                elif newpid == pid:
                    oldpids.remove(pid)
                    self.warnLog('PID %d returned with exit code %d.' % (pid,exitcode))
            return oldpids
        except:
            self.errorLog('Error in _activeForks(%s)' % (pidlist))
            raise   
#########################################################################################################################
    ### <6> ### Data dictionary methods                                                                                 #
#########################################################################################################################
    def setDictData(self,datadict,dictkey='Data'):  ### Sets object dict values
        '''
        Sets object Data dict values.
        >> datadict = Dictionary of values to add to self.dict[dictkey]
        '''
        try:
            if not self.dict.has_key(dictkey): self.dict[dictkey] = {}
            for key in datadict.keys(): self.dict[dictkey][key] = datadict[key]
        except: self.errorLog('Problem with setDictData()',True)
#########################################################################################################################
    def getData(self,dkey,dlist=['str','int','num','bool'],case=False,str=True,default=None,dp=-1):  ### Gets data from dict['Data']
        '''
        Gets data from dict['Data'] if it exists and has key, else tries dlist dictionaries.
        >> dkey:str = Key for dictionaries
        >> dlist:list [self.int,self.num,self.str,self.bool] = list of dictionaries to try after dict['Data']
        >> case:bool [False] = whether to match case for dkey
        >> str:bool [True] = whether to return all values as a string
        >> default [None] = what to return if no dictionary has key
        >> dp:int [-1] = Number of decimal places to use for stats if str=True (-1 = return as is)
        << returns value from data or int/num/str/bool
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dictlist = []
            if self.dict.has_key('Data'): dictlist = [self.dict['Data']]
            ddict = {'str':self.str,'int':self.int,'num':self.num,'bool':self.bool}
            for dict in dlist:
                if ddict.has_key(dict): dictlist.append(ddict[dict])
                else: dictlist.append(dict)
            ### ~ [1] ~ Look in dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = default
            for dict in dictlist:
                if data == default: data = rje.getFromDict(dict,dkey,returnkey=False,case=case,default=default)
                if data != default and dict in [self.int,self.num]:
                    if dp > 0: data = int(data * (10 ** dp) + 0.5) / float(10 ** dp)
                    elif dp == 0: data = int(data + 0.5)
            ### ~ [2] ~ Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if str: return '%s' % data
            else: return data
        except: self.log.errorLog('Problem with %s.getData(%s)' % (self,dkey))
        return default
#########################################################################################################################
    ### <7> ### Backwards compatibility methods                                                                         #
#########################################################################################################################
    def backConvert(self):  ### Back-converts object into old style with info, opt and stat dictionaries.
        '''Back-converts object into old style with info, opt and stat dictionaries.'''
        self.old = True
        self.info = self.str; self.infolist = self.strlist; self.str = {}; self.strlist = []
        self.opt = self.bool; self.optlist = self.boollist; self.bool = {}; self.boollist = []
        self.stat = rje.combineDict(self.int,self.num); self.statlist = self.intlist + self.numlist
        self.int = {}; self.num = {}; self.intlist = []; self.numlist = []
#########################################################################################################################
    def setInfo(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if not self.old: return self.setStr(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.info[key] = attdict[key]
                if addtolist and key not in self.infolist: self.infolist.append(key)
        except: self.errorLog('Problem with %s.setInfo()' % self,True)
#########################################################################################################################
    def setOpt(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if not self.old: return self.setBool(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.opt[key] = attdict[key]
                if addtolist and key not in self.optlist: self.optlist.append(key)
        except: self.errorLog('Problem with %s.setOpt()' % self,True)
#########################################################################################################################
    def setStat(self,attdict,addtolist=True):  ### Sets object information and optionally add keys to list
        '''Sets object information and optionally add keys to list.'''
        if not self.old: return self.setNum(attdict,addtolist)
        try:
            for key in attdict.keys():
                self.stat[key] = attdict[key]
                if addtolist and key not in self.statlist: self.statlist.append(key)
        except: self.errorLog('Problem with %s.setNum()' % self,True)
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
            if numeric: return rje.getFloat(text=text,default=value)
            elif boolean: return rje.getBool(text=text,default=value)
            else: return rje.choice(text=text,default=value)
        except: self.log.errorLog('Major problem with _editChoice()')
#########################################################################################################################
    ### <8> ### Methods for open file handles in self.file                                                              #
#########################################################################################################################
    def FILE(self,key): ### Returns file handle for key. Opens if needed.
        '''Returns file handle for key. Opens if needed.'''
        if key in self.file and self.file[key]: return self.file[key]
        return self.open(key)
#########################################################################################################################
    def open(self,key,filename=None,warnings=True,action='r'):     ### Opens file handle. If filename will populate self.str[key]
        '''
        Opens file handle. If filename will populate self.str[key].
        >> key:str = Key for self.file and self.str dictionaries.
        >> filename:str [None] = Will place into self.str or take from here if None.
        >> warnings:bool [True] = Whether to warn about things being overwritten etc.
        >> action:str ['r'] = Action for file handle (r/w/a)
        << returns open file handle object.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Close open file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.close(key)     # self.file[key] will now be None
            ## ~ [0b] Get filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if filename and filename.lower() != 'none':
                if self.getStrLC(key) and self.getStr(key) != filename and warnings:
                    self.warnLog('Over-writing %s filename (%s -> %s)' % (key,self.getStr(key),filename),'filename_overwrite',suppress=True)
            else: filename = self.getStr(key)
            self.setStr({key:filename})
            ### ~ [1] Open file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile(filename): self.file[key] = open(filename,action)
            elif warnings: self.warnLog('File "%s" does not exist.' % filename,'file_missing')
        except:
            self.errorLog('Problem with %s.open(%s)' % (self,key),True)
            self.close(key)
        return self.file[key]
#########################################################################################################################
    def close(self,key=None,expect=False):   ### Closes file handle if open. If no key given (or all/*), will close all.
        '''
        Closes file handle if open. If no key given (or all/*), will close all.
        >> key:str [None] = matches self.file key.
        >> expect:bool [False] = whether to warn if file not found
        '''
        try:### ~ [1] Close all if no key ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not key or (key not in self.file and key.lower() in ['all','*']):
                for fkey in self.file: self.close(fkey)
                return
            elif key not in self.file:
                if expect: self.warnLog('Trying to close non-existent "%s" filehandle.' % key,'No_filehandle')
                self.file[key] = None
                return
            ### ~ [2] Close file handle if open ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.file[key]: self.file[key].close(); self.file[key] = None
        except: self.errorLog('Problem with %s.close(%s)' % (self,key),True)
        self.file[key] = None
#########################################################################################################################
    def endPos(self,key): return rje.endPos(self.FILE(key))
    def fileProg(self,key,endpos): return '%.2f%%' % (100.0 * self.FILE(key).tell() / endpos)
#########################################################################################################################
    def readline(self,key,chomp=True): return rje.nextLine(self.FILE(key),strip=chomp)  ### Returns next line or None.
    def readDelimit(self,key,delimit='\t'): return rje.readDelimit(self.readline(key),delimit)  ### Returns list or None.
#########################################################################################################################
    def findline(self,key,word,wrap=True,chomp=True):  ### Scan through file to find line starting with word.
        '''
        Scan through file to find line starting with word.
        >> key:str = self.file key.
        >> word:str = word to find at beginning of line using startswith()
        >> wrap:bool [True] = Whether to scan whole file, jumping to start if end reached.
        >> chomp:bool [True] = Whether to strip /r and /n from end of line.
        << returns line as string or empty string if not found
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] File handle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            FILE = self.FILE(key)
            if not FILE: raise IOError('Cannot open "%s" file' % key)
            ## ~ [0b] Current file position & first line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            startpos = FILE.tell()   # Current file position
            wrapped = False          # Whether reading has wrapped around to start of file
            fline =  rje.fileLineFromSeek(FILE,startpos,reseek=False,next=False)[0]    # Current line
            if not fline and wrap:
                FILE.seek(0)
                wrapped = True
                fline = FILE.readline()
            ### ~ [1] Scan through file for word ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while fline:
                #self.debug('%s in %s? %s' % (word,fline,fline.startswith(word)))
                if fline.startswith(word):  # Found word!
                    if chomp: return rje.chomp(fline)
                    else: return fline
                if wrapped and FILE.tell() >= startpos: break   # Gone full cycle
                fline = FILE.readline()
                if not fline and wrap and not wrapped: FILE.seek(0); wrapped = True; fline = FILE.readline()
        except: self.errorLog('Problem with %s.findline(%s)' % (self,key),True)
        return ''
#########################################################################################################################
    def findlines(self,key,wordlist,asdict=True,wrap=True,chomp=True,warnings=True,nextonly=False):  ### Scan through file to find lines starting with word.
        '''
        Scan through file to find line starting with word.
        >> key:str = self.file key.
        >> wordlist:str = list of words to find at beginning of line using string.split()
        >> asdict:bool [True] = return {word:line} dictionary (blank if missing). !!! Assumes unique line per word !!!
        >> wrap:bool [True] = Whether to scan whole file, jumping to start if end reached.
        >> chomp:bool [True] = Whether to strip /r and /n from end of line.
        >> warnings:bool [True] = Whether to warn if dictionary values are over-written.
        >> nextonly:bool [False] = Over-rules asdict and will return the next line found or '' if EOF reached.
        << returns lines as list or {word:line} dictionary.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Return dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if asdict:
                found = {}
                for word in wordlist: found[word] = {}
            else: found = []
            ## ~ [0b] File handle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            FILE = self.FILE(key)
            if not FILE: raise IOError('Cannot open "%s" file' % key)
            ## ~ [0c] Current file position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            startpos = FILE.tell()  # Current file position
            wrapped = False         # Whether reading has wrapped around to start of file
            #fline =  rje.fileLineFromSeek(FILE,startpos,reseek=False,next=False)[0]    # Current line
            fline = FILE.readline()
            if not fline and wrap:  # Return to start of file and read line
                FILE.seek(0)
                wrapped = True
                fline = FILE.readline()
            ### ~ [1] Scan through file for word ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            while fline:
                word = string.split(fline)[0]
                if word in wordlist:  # Found word!
                    if chomp: fline = rje.chomp(fline)
                    if nextonly: return fline
                    if asdict:
                        if warnings and found[word]: self.warnLog('%s file line for %s being overwritten!' % (key,word),'findlines_overwrite',suppress=True)
                        found[word] = fline
                    else: found.append(fline)   #!# Why is first line being repeated if first line encountered?
                if wrapped and FILE.tell() >= startpos: break   # Gone full cycle
                fline = FILE.readline()
                if not fline and wrap and not wrapped: FILE.seek(0); wrapped = True; fline = FILE.readline()
            if nextonly: return ''
        except: self.errorLog('Problem with %s.findlines(%s)' % (self,key),True)
        return found
#########################################################################################################################
### End of RJE_Object Class                                                                                             #
#########################################################################################################################

#########################################################################################################################
class Log(rje.Log):
    '''Class to handle log output: printing to log file and error reporting.'''
#########################################################################################################################
class Out(rje.Out):
    '''Class to handle basic generic output to screen based on Verbosity and Interactivity outside of Classes.'''
#########################################################################################################################
def zen(): return rje_zen.Zen().wisdom()
#########################################################################################################################
