#!/usr/local/bin/python

# See below for name and description
# Copyright (C) 2005 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 25 Bassett Court, Bassett Avenue, Southampton SO16 7DR, UK.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       file_monster
Description:  Goes through directories etc and collects information on files etc. 
Version:      2.2
Last Edit:    05/12/13
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Functions:
    <1> File Monster:   [scavenge=T, organise=X and/or cleanup=T]
    Goes through directories etc and collects information on files etc. "scavenge=T" will output a file containing
    relevant file names and locations, ages and sizes. "monster=T" will compare for identical files and farm them off
    into another directory for possible deletion. 

    <2> DirSum: [dirsum=T]
    Goes through directories and subdirectories and summarises the number of files and subdirectories they contain. If a
    directory contains less than X subdirectories [dircut=X] then the subdirectories of that directory will also be
    summarised. Output is in the form: PATH,files,subdir

    <3> Rename: [rename=T]
    Renames all the chosen files [filelist] with the given prefix [prefix=X] into outdir. If usedate=T, dates will be
    added to the prefix (e.g. outdir/prefixdate_num) otherwise new names are just outdir/prefix_num.

    <4> Fix line endings [fixendings=LIST]
    Replace Mac \\r with \\n line endings in place with option to backup old file (unless backups=F).

Output:
    File = File (or directory) name (no path)
    Parent = Parent directory  from DirList
    Folder = Path containing file (or directory) 
    Type = File extension (or "DIR")
    Size = Size in bytes
    Date = Age converted to human readable string
    Age = Age of file (MTime) in seconds
    CTime = String representation of Creation Time
    MTime = String representation of Modified Time
    ATime = String representation of ATime
    FilePath = Full path to file

General Commandline:
    filelist=X,Y,..,Z   : List of files of interest. Can have wildcards. [*]
    dirlist=LIST        : List of directories to look in, in order of preference good -> bad for duplicates. [./]
    subfolders=T/F      : Whether to look in subfolders [True]
    stripnum=T/F        : Whether files may have -XXX numerical suffix from renaming, which should be stripped [False]

File Monster Commandline:
    oldmonster=T/F      : Whether to run old File Monster (V1.x). Will be retired once update complete. [False]
    outdir=PATH         : Output directory for renamed/reorganised files [./Organised/]
    dumpdir=PATH        : Directory in which to dump redundant files (don't move if "None") [./Redundant/]
    cleanup=T/F         : Whether to delete empty directories (and move/delete stripnum) [False]
    cleanfiles=LIST     : List of hidden files to delete during cleanup if directory seems empty ['.picasa.ini','.DS_Store']
    sizematch=X         : Size % similarity threshold to count as match [99.9]
    skiplist=LIST       : List of filenames to skip [thumbs.db]
    organise=X          : File reorganisation mode (none/date/month/compile) [None]
    orgprefix=X         : Prefix for organised outdir subdirectories ['']
    redundancy=T/F      : Whether to check/rate redundancy for scavenge etc. [True]
    scavenge=T/F        : Whether to perform collation of file information [False]
    searchid=X          : ID for search - allows multiple searches to be compared easily [None]

DirSum Commandline:    
    dirsum=T/F          : Whether to perform summary of directory contents [False]
    dircut=X            : Max number of subdirectories to have and still delve into them [50]
    dirdepth=X          : Max depth of subdirectories to delve into. Negative = all. [-1]
    extlist=LIST        : List of file extenstions to report individual stats for ['']

LineEndings Commandline:
    fixendings=FILELIST : Replace Mac with UNIX line endings for FILELIST (wildcards allowed) []

OLD COMMANDS:

File Monster Commandline:
    monster=T/F         : Whether to perform monster cleanup of redundant files [False]
    gooddir=LIST        : List of "good" directories to be automatically kept if i<1 (including subdirs) []
    baddir=LIST         : List of "bad" directories to be automatically screened if i<1 (including subdirs) []
    keepnew=T/F         : Preferentially keep newer files of same size if good/bad status equal [True]
    purgelist=LIST      : List of filenames (allowing wildcards) to purge (move/delete) [WS_FTP.LOG]

Rename Commandline:
    rename=T/F          : Whether to rename files [False]
    sortby=X            : Whether to sort by date or name [date]
    prefix=X            : Text prefix for new file names []
    usedate=T/F         : Whether to use date in new name [False]

Uses general modules: copy, glob, os, stat, re, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import stat, re, traceback
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial Working version
    # 1.1 - Broadened away from strict extension-based scavenging to whole file names with wildcards
    # 1.2 - Added DirSum function and updated FileMonster slightly.
    # 1.3 - Added redundant file cleanup
    # 1.4 - Added skiplist and purgelist
    # 1.5 - Added rename function (to replace rename.pl Perl module)
    # 1.6 - Minor bug fix.
    # 2.0 - Major reworking with new object making use of rje_db tables etc. Old functions to be ported with time.
    # 2.1 - Added dirsum function.
    # 2.2 - Added fixendings=FILELIST to convert Mac \\r into UNIX \\n
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Removal of empty directories during dirsum
    # [ ] : BadDir and GoodDir lists for automated filtering (i=0)
    # [ ] : More options in yes/no - add to baddir/gooddir lists. Move previous file instead of later.
    # [ ] : Check it still works and add du -h function.
    # [ ] : Use Database objects for better handling and analysis of directory contents etc.
    # [ ] : Convert dirsum=T function into creating a summary table (also parsing du -h) rather than log lines only.
    # [ ] : Add a rename function that will replace part of file name.
    # [ ] : Convert all old functions to new Class.
    # [Y] : Add reorganise function and organise by date.
    # [Y] : Allow wildcards for dirlist selection.
    # [Y] : Use os.listdir() instead of glob tp retrieve files lists
    # [ ] : Add abspath option to use absolute paths. (Currently relative)
    # [ ] : Consider integrating the fix endings function into scavenge
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('FileMonster', '2.2', 'December 2013', '2007')
    description = 'File/directory summary/manipulation program'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: FileMonster Class                                                                                       #
#########################################################################################################################
class FileMonster(rje_obj.RJE_Object):     
    '''
    FileMonster Class. Author: Rich Edwards (2012).

General Commandline:
    filelist=X,Y,..,Z   : List of files of interest. Can have wildcards. [*]
    dirlist=LIST        : List of directories to look in, in order of preference good -> bad for duplicates. Over-rules folder=PATH []
    subfolders=T/F      : Whether to look in subfolders [True]

File Monster Commandline:
    outfile=FILE        : File for output [file_monster.csv]
    searchid=X          : ID for search - allows multiple searches to be compared easily [None]
    scavenge=T/F        : Whether to perform collation of file information [False]
    monster=T/F         : Whether to perform monster cleanup of redundant files [False]
    dumpdir=PATH        : Directory in which to dump redundant files [Redundant]
    gooddir=LIST        : List of "good" directories to be automatically kept if i<1 (including subdirs) []
    baddir=LIST         : List of "bad" directories to be automatically screened if i<1 (including subdirs) []
    keepnew=T/F         : Preferentially keep newer files of same size if good/bad status equal [True]
    skiplist=LIST       : List of filenames to skip [thumbs.db]
    purgelist=LIST      : List of filenames (allowing wildcards) to purge (move/delete) [WS_FTP.LOG]

DirSum Commandline:    
    dirsum=T/F          : Whether to perform summary of directory contents [False]
    dircut=X            : Max number of subdirectories to have and still delve into them [20]

Rename Commandline:
    rename=T/F          : Whether to rename files [False]
    outdir=PATH         : Output directory for renamed files [./]
    sortby=X            : Whether to sort by date or name [date]
    prefix=X            : Text prefix for new file names []
    usedate=T/F         : Whether to use date in new name [False]

Obselete:
    folder=PATH         : Folder to analyse [./]

    - Name = ID for search - allows multiple searches to be compared easily [None]
    - Folder = Folder to analyse [./]
    - OrgPrefix = Prefix for organised outdir subdirectories ['']
    - OutFile = File for output [file_monster.txt]
    - SortBy = Whether to sort by date or name [date]
    - PreFix = Text prefix for new file names []
    
    Opt:boolean
    - Scavenge = Whether to perform collation of file information [False]
    - CleanUp = Whether to delete empty directories (and move/delete stripnum) [False]
    - Monster = Whether to perform monster cleanup of redundant files [False]
    - KeepNew = Preferentially keep newer files of same size if good/bad status equal [True]
    - Rename = Whether to rename files [False]
    - StripNum = Whether files may have -XXX numerical suffix from renaming, which should be stripped [False]
    - UseDate = Whether to use date in new name [False]

    Stat:numeric

    List:list
    - GoodDir = List of "good" directories to be given priority in decision making (Can be partial paths) []
    - BadDir = List of "bad" directories to be automatically screened if i<1 (Can be partial paths) []
    - ExtList = List of file extenstions to report individual stats for ['']
    - PurgeList = List of filenames (allowing wildcards) to purge (move/delete) [WS_FTP.LOG]

    <<<< >>>>

    Str:str
    - DumpDir = Directory in which to dump redundant files [./Redundant/]
    - Organise = File reorganisation mode (none/date) [None]
    - OutDir = Output directory for renamed/reorganised files [./]
    
    Bool:boolean
    - DirSum = Whether to perform summary of directory contents [False]
    - OldMonster = Whether to run old File Monster (V1.x). Will be retired once update complete. [False]
    - Redundancy = Whether to check/rate redundancy for scavenge etc. [True]
    - SubFolders = Whether to look in subfolders [True]

    Int:integer
    - DirCut = Max number of subdirectories to have and still delve into them [50]
    - DirDepth = Max depth of subdirectories to delve into. Negative = all. [-1]

    Num:float
    - SizeMatch = Size % similarity threshold to count as match [99.95]
    
    List:list
    - CleanFiles = List of hidden files to delete during cleanup if directory seems empty ['.picasa.ini']
    - DirList = List of directories to look in, in order of preference good -> bad. Over-rules folder=PATH []
    - FileList = List of files of interest (X,Y,..,Z)  [*]
    - FixEndings = Replace Mac with UNIX line endings for FILELIST (wildcards allowed) []
    - SkipList = List of filenames to skip [thumbs.db]
    
    Dict:dictionary    

    Obj:RJE_Objects
    - DB = Database Object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['DumpDir','Organise','OrgPrefix','OutDir','SearchID']
        self.boollist = ['CleanUp','DirSum','OldMonster','Redundancy','SubFolders','Scavenge','StripNum']
        self.intlist = ['DirCut','DirDepth']
        self.numlist = ['SizeMatch']
        self.listlist = ['CleanFiles','DirList','ExtList','FileList','FixEndings','SkipList']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'DumpDir':rje.makePath('./Redundant/'),'OutDir':rje.makePath('./Organised/'),'OrgPrefix':''})
        self.setBool({'Redundancy':True,'SubFolders':True})
        self.setInt({'DirCut':50,'DirDepth':-1})
        self.setNum({'SizeMatch':0.999})
        self.list['DirList'] = ['./']
        self.list['CleanFiles'] = ['.picasa.ini','.DS_Store']
        self.list['FileList'] = ['*']
        self.list['SkipList'] = ['thumbs.db','Thumbs.db']
        self.baseFile('filemonster')
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
                self._cmdReadList(cmd,'str',['Organise','SearchID','OrgPrefix'])   # Normal strings
                self._cmdReadList(cmd,'path',['DumpDir','OutDir'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['Att'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['CleanUp','DirSum','OldMonster','Redundancy','SubFolders','Scavenge','StripNum'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['DirCut','DirDepth'])   
                self._cmdReadList(cmd,'float',['SizeMatch']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CleanFiles','ExtList','FileList','SkipList'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['DirList','FixEndings']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.deBug(self.boollist); self.deBug(self.bool)
        while self.getNum('SizeMatch') > 100.0: self.num['SizeMatch'] /= 100.0
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('OldMonster'): return FileScout(self.log,self.cmd_list).run()
            ### ~ [4] Fix Line Endings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['FixEndings']:
                self.printLog('#FIX','Attempting to fix line endings for %s files.' % rje.iStr(len(self.list['FixEndings'])))
                for file in self.list['FixEndings']: self.fixEnding(file)
                return True
            elif not self.setup(): return self.printLog('#END','No files/folders found.')    # This now reads in all files and directories
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('CleanUp'): self.cleanup()
            ## Invoke skiplist here
            if self.list['SkipList']:
                self.printLog('#SKIP','%s SkipList files/types to skip.' % (rje.iLen(self.list['SkipList'])))            
                self.db('Files').dropEntriesDirect('File',self.list['SkipList'],inverse=False,log=True)
                self.db('Files').dropEntriesDirect('Type',self.list['SkipList'],inverse=False,log=True)
            if self.getBool('Scavenge') and self.getBool('Redundancy'): self.redundancy()
            if self.organise() and self.getBool('CleanUp'):
                if self.getBool('Scavenge') and rje.yesNo('Output pre-tidy file scavenge?'):
                    fdb = self.db('Files')
                    if self.getStr('SearchID').lower() not in ['none','']: fdb.addField('SearchID',evalue=self.getStr('SearchID'))
                    fdb.saveToFile(append=self.getBool('Append'))
                    self.setBool({'Scavenge':False})
                if rje.yesNo('Cleanup remaining directories etc.?') and self.setup(): self.cleanup(stripnum=False)
                cleandir = []
                if rje.yesNo('Cleanup post-organise %s?' % self.getStr('OutDir')): cleandir.append(self.getStr('OutDir'))
                if os.path.exists(self.getStr('DumpDir')) and rje.yesNo('Cleanup post-organise %s?' % self.getStr('DumpDir')): cleandir.append(self.getStr('DumpDir'))
                if cleandir:
                    self.list['DirList'] = cleandir
                    if self.setup(): self.cleanup()
            if self.getBool('Scavenge'):
                fdb = self.db('Files')
                if self.getStr('SearchID').lower() not in ['none','']: fdb.addField('SearchID',evalue=self.getStr('SearchID'))
                fdb.saveToFile(append=self.getBool('Append'))
            ### ~ [3] ~ Directory Summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DirSum'):
                ddb = self.db().copyTable(self.db('Files'),'DirSum')
                ddb.keepFields(['File','Parent','Folder','Type','Size','Age','FilePath'])
                dtot = ddb.entryNum(); dx = 0.0
                for entry in ddb.entries():
                    self.progLog('\r#DIR','Adding entries for DirSum compression: %.2f%%' % (dx/dtot)); dx += 100.0
                    if entry['Type'] not in ['DIR'] + self.list['ExtList']: entry['Type'] = 'file'
                    elif entry['Type'] != ['DIR']: ddb.addEntry(rje.combineDict({'Type':'file','FilePath':entry['FilePath']+'filemonstercopy'},entry,False))
                self.printLog('\r#DIR','Added entries for DirSum compression: %s "file" entries.' % rje.iStr(ddb.entryNum()-dtot))
                ddb.addField('N',evalue=1)
                folders = ddb.indexDataList('Type','DIR','File')
                ddb.compress(['Folder','Type'],{'N':'sum','Age':'max','Size':'sum'})
                ddb.dropFields(['File','FilePath'])
                ddb.reshapeWide('Type',['N','Age','Size'])
                #!# Add up within parents and filter DirCut folders #!#
                ddb.saveToFile()                

            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            db.baseFile(self.baseFile())
            inx = len(self.list['DirList'])
            for folder in self.list['DirList'][0:]:
                if not os.path.isdir(folder): self.list['DirList'].remove(folder)
            self.printLog('#DIR','%s folders from %s possible dirlist paths' % (rje.iLen(self.list['DirList']),rje.iStr(inx)))
            ### ~ [2] Read Files into Database Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = db.addEmptyTable('Files',['File','Parent','Folder','Type','Size','Date','Age','CTime','MTime','ATime','FilePath'],['FilePath'])
            ## ~ [2a] Read full file list into table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #skipx = 0
            self.progLog('\r#FILES','Scavenging from %d top level directories...' % (len(self.list['DirList'])))
            self.deBug(self.progLog)
            for folder in self.list['DirList']:
                folderfirst = True
                for gfile in [folder] + rje.listDir(self,folder,self.getBool('SubFolders'),dircut=self.getInt('DirCut'),dirdepth=self.getInt('DirDepth')):
                    #self.deBug(file)
                    details = {'Parent':folder,'Size':0,'ATime':'???','MTime':'???','CTime':'???'}
                    details['File'] = os.path.basename(gfile)                            # Basic file name without path
                    details['Type'] = os.path.splitext(gfile)[1].lower()[1:]                         # File type (extension)
                    if os.path.isdir(gfile): details['Type'] = 'DIR'
                    #if details['File'] in self.list['SkipList']: skipx += 1; continue
                    #if details['Type'] in self.list['SkipList']: skipx += 1; continue
                    if details['Type'] != 'DIR' and self.getBool('StripNum') and rje.matchExp('^(.+)-\d+$',os.path.splitext(details['File'])[0]):
                        details['File'] = '%s.%s' % (rje.matchExp('^(.+)-\d+$',os.path.splitext(details['File'])[0])[0],details['Type'])
                    details['FilePath'] = gfile #os.path.abspath(file)                         # Full path to file
                    details['Path'] = os.path.dirname(details['FilePath'])              # Full path to folder (no file)
                    if folderfirst: details['Parent'] = details['Folder'] = ''; folderfirst = False
                    else: details['Folder'] = string.split(details['Path'])[-1]               # Immediate folder containing file
                    if details['Type'] != 'DIR': details['File'] = '%s.%s' % (os.path.splitext(details['File'])[0],details['Type'])
                    try:    # File details - size and age
                        #size = os.path.getsize(file)
                        details['Age'] = os.path.getmtime(gfile)
                        details['Date'] = string.split(rje.dateTime(time.localtime(details['Age'])))[0]
                        details['Size'] = int(os.path.getsize(gfile))
                        details['CTime'] = time.ctime(os.path.getctime(gfile))
                        details['MTime'] = time.ctime(os.path.getmtime(gfile))
                        details['ATime'] = time.ctime(os.path.getatime(gfile))
                    except: self.errorLog('Error in scavenge %s size/date.' % gfile,printerror=True,quitchoice=False)
                    #if details['MTime'] != details['CTime']: self.deBug('%s\nCTime: %s\nMTime: %s\n' % (file,details['CTime'],details['MTime']))
                    fdb.addEntry(details)
            self.printLog('\r#FILES','%s files scavenged from %d top level directories.' % (rje.iStr(fdb.entryNum()),len(self.list['DirList'])))            
            #self.printLog('#SKIP','%s SkipList files skipped.' % (rje.iStr(skipx)))            
            return fdb.entryNum()     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def cleanup(self,stripnum=True):  ### Deletes empty directories and move/delete stripnum copies
        '''Deletes empty directories and move/delete stripnum copies.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('Files')
            ### ~ [1] Delete PurgeList Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add removal of purgelist
            #fdb.dropEntriesDirect('File',self.list['PurgeList'],log=True)
            ### ~ [2] Delete/Move Duplicate Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if stripnum and self.getBool('StripNum'):
                deldup = rje.yesNo('Delete StripNum Duplicates?',default='N')
                movedup = not deldup and rje.yesNo('Move StripNum Duplicates to %s?' % self.getStr('DumpDir'),default='Y')
                confirm = (deldup or movedup) and rje.yesNo('Confirm?')
                checkskip = rje.yesNo('Manually confirm skipping of copies < %.2f%% size similarity?' % (100.0 * self.getNum('SizeMatch')))
                redx = 0
                for ekey in fdb.dataKeys():
                    if not deldup and not movedup: break
                    entry = fdb.data(ekey)
                    self.deBug('%s (%s)' % (entry['FilePath'],entry['File']))
                    if os.path.splitext(os.path.basename(entry['FilePath']))[0] != os.path.splitext(entry['File'])[0]:
                        origfile = self.getFilePath(os.path.join(entry['Folder'],entry['File']))
                        self.deBug(origfile)
                        if origfile in fdb.index('FilePath'): 
                            try:
                                ### Check Size ###
                                if fdb.data(origfile)['Size'] != entry['Size']:
                                    if fdb.data(origfile)['Size'] > entry['Size']: sizecomp = entry['Size'] / fdb.data(origfile)['Size']
                                    else: sizecomp = fdb.data(origfile)['Size'] / entry['Size']
                                    self.printLog('#WARN','%s and %s are different sizes! (%.2f%%)' % (entry['FilePath'],origfile,100*sizecomp))
                                    if sizecomp < self.getNum('SizeMatch'):
                                        if self.i() < 0 or not checkskip or rje.yesNo('Skip %s?' % entry['FilePath']): continue
                                ### Move or Delete ###
                                redx += 1
                                if deldup and (not confirm or rje.yesNo('Delete %s?' % entry['FilePath'])):
                                    os.unlink(entry['FilePath'])
                                    fdb.dropEntry(entry)
                                elif movedup:
                                    newpath = self.getStr('DumpDir') + os.path.basename(entry['FilePath'])
                                    dupx = 0
                                    while os.path.exists(newpath):
                                        dupx += 1
                                        newpath = '%s-%s.%s' % (os.path.splitext(self.getStr('DumpDir') + os.path.basename(entry['File']))[0],rje.preZero(dupx,1000),entry['Type'])
                                    if not confirm or rje.yesNo('Move %s to %s?' % (entry['FilePath'],newpath)):
                                        os.rename(entry['FilePath'],newpath)
                                        entry['FilePath'] = newpath
                            except KeyboardInterrupt: break
                            except: raise
                self.printLog('#RED','%s redundant stripnum files' % rje.iStr(redx))
            ### ~ [3] Delete Empty Directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            deldir = rje.yesNo('Delete Empty Directories?',default='N')
            confirm = deldir and rje.yesNo('Confirm?')
            dirdata = fdb.subset('Type','DIR')
            dirpaths = rje.sortKeys(dirdata)
            #for fdir in dirpaths[0:]:   ### BAD FUDGE. NEED TO REWORK FROM START SOME TIME USING os.listdir(folder)
            #    for hpath in glob.glob('%s/.*' % fdir):
            #        if os.path.isdir(hpath) and hpath not in dirpaths: 
            #            dirpaths.insert(0,hpath)
            #            dirdata[hpath] = fdb.addEntry({'Folder':fdir,'FilePath':hpath})
            dirx = 0; delx = 0
            remhid = {}
            for hfile in self.list['CleanFiles']: remhid[hfile] = 0
            cfiles = rje.sortKeys(remhid)
            self.deBug(cfiles)
            self.deBug(dirpaths)
            while dirpaths and deldir:
                if len(dirpaths) == dirx: break
                self.progLog('\r#RMDIR','Deleting %s empty directories' % rje.iStr(delx))
                dirx = len(dirpaths)
                for fdir in dirpaths[0:]:
                    self.deBug(fdir)
                    ffiles = os.listdir(fdir)
                    if ffiles:
                        subdir = False
                        for ffile in ffiles:
                            subdir = subdir or os.path.isdir(os.path.join(fdir,ffile))
                            if subdir: break
                        if subdir: continue     # Wait until subdirs dealth with
                        ffiles.sort()
                        self.deBug(ffiles)
                        cleanup = rje.listIntersect(cfiles,ffiles)
                        cleanup.sort()
                        self.deBug(cleanup)
                        if cleanup == ffiles:   # All files are cleanup files!
                            for ffile in ffiles:
                                os.unlink(os.path.join(fdir,ffile))
                                remhid[ffile] += 1
                            fdb.dropEntriesDirect('Folder',[fdir],log=False)
                    if fdir not in fdb.index('Folder'):     # Should be empty
                        try:
                            if not confirm or rje.yesNo('Delete empty directory: %s?' % fdir):
                                os.rmdir(fdir)
                                fdb.dropEntry(dirdata.pop(fdir))
                                dirpaths.remove(fdir)
                                delx += 1
                        except KeyboardInterrupt:
                            if rje.yesNo('Switch off confirmation?',default='N'):
                                confirm = False; dirx = 0
                            else: dirx = len(dirpaths); break
                        except: 
                            hidden = glob.glob('%s/.*' % fdir)
                            if hidden: self.printLog('#HID','Hidden files in %s: %s' % (fdir,hidden))
                            else: self.errorLog('Error!',quitchoice=False)
            self.printLog('\r#RMDIR','Deleted %s empty directories.' % rje.iStr(delx))
            for hfile in rje.sortKeys(remhid):
                if remhid[hfile]: self.printLog('#DEL','Deleted %s %s cleanup files' % (rje.iStr(remhid[hfile]),hfile))
                
        except: self.errorLog('%s.cleanup() error' % self); raise   # Delete this if method error not terrible
#########################################################################################################################
    def getFilePath(self,filepath): ### Returns the relevant FilePath key from the Files table, else None
        '''Returns the relevant FilePath key from the Files table, else None.'''
        try:### ~ [1] Easy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('Files')
            if fdb.data(filepath): return filepath
            ## ~ [1a] Try Full path ~~~~~~~~~~~~~~ ##
            #filepath = os.path.abspath(filepath)  
            #if  fdb.data(filepath): return filepath
            ## ~ [1b] Try switching case ext ~~~~~~~~~~~~~~ ##
            file = os.path.splitext(filepath)
            if file[1] != file[1].lower(): filepath = '%s%s' % (file[0],file[1].lower())
            else: filepath = '%s%s' % (file[0],file[1].upper())
            if  fdb.data(filepath): return filepath
            if file[1].lower() in ['.jpg','.png'] and fdb.data(file[0]+'.jpg'): return file[0]+'.jpg'
            ## ~ [1c] Failure ~~~~~~~~~~~~~~~~~~~ ##
            return None
        except: self.errorLog('%s.getFilePath() error' % self); raise   # Delete this if method error not terrible
#########################################################################################################################
    def makeDate(self,ftime):   ### Returns YYYY-MM-DD from CTime or MTime etc.
        '''Returns YYYY-MM-DD from CTime or MTime etc.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            noday = self.getStr('Organise').lower() == 'month'
            #Wed Dec 19 00:00:00 2012
            fdata = string.split(ftime)
            year = fdata[-1]
            month = rje.preZero(months.index(fdata[1])+1,12)
            day = rje.preZero(int(fdata[2]),31)
            if noday: return '%s-%s' % (year, month)
            return '%s-%s-%s' % (year, month, day)
        except: self.errorLog('%s.makeDate() error' % self); raise   # Delete this if method error not terrible
#########################################################################################################################
    def organise(self):   ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            orgmethod = self.getStr('Organise').lower()
            if orgmethod in ['','none']: return False
            fdb = self.db('Files')
            fdb.dropEntriesDirect('Type','DIR')
            fdb.addField('NewPath',evalue='')
            if 'Redundancy' in fdb.fields(): fdb.dropField('Redundancy')
            fdb.addField('Redundancy',evalue='')
            newfiles = {}   # Dictionary of {New Path:entry}
            ### ~ [2] Organisation by date ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if orgmethod in ['date','month','ctime','mtime','atime']:
                if orgmethod in ['month','date']: orgmethod = 'mtime'
                orgmethod = orgmethod[:2].upper() + orgmethod[2:]
                self.printLog('#ORG','Organised files by %s' % orgmethod)
                ## ~ [2a] First path for special connections, e.g. Picassa ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries():
                    if '.picasaoriginals/' not in entry['FilePath']: continue
                    entry['NewPath'] = '%s%s%s/.picasaoriginals/%s' % (self.getStr('OutDir'),self.getStr('OrgPrefix'),self.makeDate(entry[orgmethod]),entry['File'])
                    #self.deBug('%s (%s) -> %s' % (entry['FilePath'],entry[orgmethod],entry['NewPath']))
                    editpath = self.getFilePath(string.replace(entry['FilePath'],'.picasaoriginals/',''))
                    if editpath:
                        editentry = fdb.data(editpath)
                        editentry['NewPath'] = '%s%s%s/%s' % (self.getStr('OutDir'),self.getStr('OrgPrefix'),self.makeDate(entry[orgmethod]),editentry['File'])
                        #self.deBug('%s (%s) -> %s' % (editentry['FilePath'],entry[orgmethod],editentry['NewPath']))
                    else: self.errorLog('Cannot find Picassa edit for %s!' % entry['FilePath'],printerror=False)
                ## ~ [2b] Rest of files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries():
                    if entry['NewPath']: continue
                    entry['NewPath'] = '%s%s%s/%s' % (self.getStr('OutDir'),self.getStr('OrgPrefix'),self.makeDate(entry[orgmethod]),entry['File'])
                    #self.deBug('%s (%s) -> %s' % (entry['FilePath'],entry[orgmethod],entry['NewPath']))
            elif orgmethod == 'compile':
                ## ~ [2a] First path for special connections, e.g. Picassa ~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries():
                    if '.picasaoriginals/' not in entry['FilePath']: continue
                    entry['NewPath'] = '%s.picasaoriginals/%s' % (self.getStr('OutDir'),entry['File'])
                    #self.deBug('%s (%s) -> %s' % (entry['FilePath'],entry[orgmethod],entry['NewPath']))
                    editpath = self.getFilePath(string.replace(entry['FilePath'],'.picasaoriginals/',''))
                    if editpath:
                        editentry = fdb.data(editpath)
                        editentry['NewPath'] = '%s%s' % (self.getStr('OutDir'),editentry['File'])
                        #self.deBug('%s (%s) -> %s' % (editentry['FilePath'],entry[orgmethod],editentry['NewPath']))
                    else: self.errorLog('Cannot find Picassa edit for %s!' % entry['FilePath'],printerror=False)
                ## ~ [2b] Rest of files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for entry in fdb.entries():
                    if entry['NewPath']: continue
                    entry['NewPath'] = '%s%s' % (self.getStr('OutDir'),entry['File'])
                    #self.deBug('%s (%s) -> %s' % (entry['FilePath'],entry[orgmethod],entry['NewPath']))
                
            ### ~ [3] Check New Paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Build lists of files for each new path ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for entry in fdb.entries():
                newpath = entry['NewPath']
                if newpath not in newfiles: newfiles[newpath] = [entry]
                else:
                    nentry = newfiles[newpath][0]
                    if nentry['Age'] < entry['Age']: 
                        newfiles[newpath].append(entry)
                        if not nentry['Redundancy']: 
                            nentry['Redundancy'] = 'Oldest'; entry['Redundancy'] = 'Copy?'
                    elif entry['Age'] < nentry['Age']: 
                        newfiles[newpath].insert(0,entry)
                        entry['Redundancy'] = 'Oldest'
                        nentry['Redundancy'] = 'Copy?'
                    elif self.list['DirList'].index(nentry['Parent']) < self.list['DirList'].index(entry['Parent']): 
                        newfiles[newpath].append(entry)
                        if not nentry['Redundancy']: nentry['Redundancy'] = entry['Parent']
                    elif self.list['DirList'].index(entry['Parent']) < self.list['DirList'].index(nentry['Parent']): 
                        newfiles[newpath].insert(0,entry)
                        entry['Redundancy'] = entry['Parent']
                    elif entry['Size'] < nentry['Size']: 
                        newfiles[newpath].append(entry)
                        if not nentry['Redundancy']: nentry['Redundancy'] = 'Biggest'
                    elif entry['Size'] > nentry['Size']: 
                        newfiles[newpath].insert(0,entry)
                        entry['Redundancy'] = 'Biggest'
                    else: newfiles[newpath].append(entry)
            ## ~ [3b] Update newfile dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodump = self.getStr('DumpDir').lower() == rje.makePath('none')
            renumber = rje.yesNo('Renumber duplicates with same name into %s?' % self.getStr('OutDir'))
            #!# Add option to rename with numerical suffix in DumpDir #!#
            if renumber: newpathlist = newfiles.keys()[0:]
            for newpath in newfiles:
                if len(newfiles[newpath]) == 1: newfiles[newpath][0]['Redundancy'] = 'Unique'
                for entry in newfiles[newpath][1:]:
                    entry['Redundancy'] = 'Redundant'
                    if renumber:
                        if orgmethod == 'compile': dumppath = self.getStr('OutDir')
                        else: dumppath = '%s%s%s/' % (self.getStr('OutDir'),self.getStr('OrgPrefix'),self.makeDate(entry[orgmethod]))
                        entry['NewPath'] = '%s%s' % (dumppath,entry['File'])
                        rx = 0
                        while entry['NewPath'] in newpathlist:
                            rx += 1
                            (root,ext) = os.path.splitext(entry['File'])
                            entry['NewPath'] = '%s%s-%s%s' % (dumppath,root,rje.preZero(rx,1000),ext)
                            self.deBug(entry['NewPath'])
                        newpathlist.append(entry['NewPath'])
                    elif nodump: entry['NewPath'] = entry['FilePath']
                    elif entry == newfiles[newpath][1]:
                        entry['NewPath'] = '%s%s%s/%s' % (self.getStr('DumpDir'),self.getStr('OrgPrefix'),self.makeDate(entry[orgmethod]),entry['File'])
                        newpathlist.append(entry['NewPath'])
                    else:
                        mainpath = rje.matchExp('^%s(.+)' % entry['Parent'],entry['FilePath'])[0]
                        entry['NewPath'] = rje.makePath(self.getStr('DumpDir')+mainpath,wholepath=True)
                newfiles[newpath] = newfiles[newpath][0]
            ## ~ [3c] Move appropriate files, creating directories along the way ~~~~~~~~~~~~~~~~~~ ##
            mkdir = rje.yesNo('Make Directories?')
            mvfiles = mkdir and rje.yesNo('Move files?')
            if not mkdir: return
            mx = 0
            for ekey in fdb.dataKeys():
                entry = fdb.data(ekey)
                if mkdir: rje.mkDir(self,entry['NewPath'],log=True)
                if entry['FilePath'] != entry['NewPath']:
                    if mvfiles: os.rename(entry['FilePath'],entry['NewPath']); mx += 1
                    else: self.printLog('#MV','%s -> %s' % (entry['FilePath'],entry['NewPath']))
            self.printLog('#MOVE','%s files moved' % rje.iStr(mx))
            return mx

            #x#self.redundancy('NewPath')
                    
        except: self.errorLog('%s.organise() error' % self); raise   # Delete this if method error not terrible
#########################################################################################################################
    def redundancy(self,redfield='File',redatt=True):   ### Identifies duplicate files and adds Redundancy to Table
        '''
        Identifies duplicate files and adds Redundancy to Table.
        >> redfield:str ['File'] = Main field of Files table to use to identify redundancy.
        >> redatt:bool [True] = Whether to also look for redundant files based on attributes.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdb = self.db('Files')  #!# Deal with different case extensions at some point!
            if 'Redundancy' in fdb.fields(): fdb.dropField('Redundancy')
            fdb.addField('Redundancy',evalue='')
            if 'Att' not in fdb.fields(): fdb.makeField('#MTime#|#Size#',fieldname='Att')
            redlist = ['Conflict','Duplicate','Renamed','Backup','Copy?','Copy','Old Version','New Version','Copied','Original','']
            ### ~ [2] Redundancy by redfield ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(fdb.index(redfield)) == fdb.entryNum():
                self.printLog('#NR','All files unique for %s' % redfield)
            else:
                self.printLog('#RED','Only %s unique %s from %s files' % (rje.iLen(fdb.index(redfield)),redfield,rje.iStr(fdb.entryNum())))
                for unique in fdb.index(redfield):
                    if len(fdb.index(redfield)[unique]) == 1: continue
                    redfiles = fdb.indexEntries(redfield,unique)
                    for i in range(len(redfiles)-1):
                        ientry = redfiles[i]
                        irenamed = self.getBool('StripNum') and os.path.splitext(os.path.basename(ientry['FilePath']))[0] != os.path.splitext(ientry['File'])[0]
                        for j in range(i+1,len(redfiles)):
                            jentry = redfiles[j]
                            jrenamed = self.getBool('StripNum') and os.path.splitext(os.path.basename(jentry['FilePath']))[0] != os.path.splitext(jentry['File'])[0]
                            ## Total Match ##
                            if '.picasaoriginals/' in jentry['FilePath'] and '.picasaoriginals/' not in ientry['FilePath']:
                                if redlist.index('Original') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Original'
                                if redlist.index('Backup') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Backup'
                            elif '.picasaoriginals/' in ientry['FilePath'] and '.picasaoriginals/' not in jentry['FilePath']:
                                if redlist.index('Original') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Original'
                                if redlist.index('Backup') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Backup'
                            elif irenamed and not jrenamed:
                                if redlist.index('Renamed') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Renamed'
                                if redlist.index('Copied') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Copied'
                            elif jrenamed and not irenamed:
                                if redlist.index('Renamed') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Renamed'
                                if redlist.index('Copied') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Copied'
                            elif ientry['Att'] == jentry['Att']:
                                if redlist.index('Copied') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Duplicate'
                                if redlist.index('Copy') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Duplicate'
                            elif ientry['Size'] == jentry['Size'] and ientry['Age'] < jentry['Age']:
                                if redlist.index('Original') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Original'
                                if redlist.index('Copy') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Copy?'
                            elif ientry['Size'] == jentry['Size'] and jentry['Age'] < ientry['Age']:
                                if redlist.index('Original') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Original'
                                if redlist.index('Copy') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Copy?'
                            elif ientry['Age'] > jentry['Age']:
                                if redlist.index('Old Version') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Old Version'
                                if redlist.index('New Version') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'New Version'
                            elif ientry['Age'] < jentry['Age']:
                                if redlist.index('Old Version') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Old Version'
                                if redlist.index('New Version') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'New Version'
                            else:
                                if redlist.index('Conflict') < redlist.index(ientry['Redundancy']): ientry['Redundancy'] = 'Conflict'
                                if redlist.index('Conflict') < redlist.index(jentry['Redundancy']): jentry['Redundancy'] = 'Conflict'
                            continue
                            #!# Need to think about samesize for photos etc. #!#

                            #else:
                            #    ientry['Redundancy'] = jentry['Redundancy'] = 'Same Name'
                                
                for redtype in fdb.index('Redundancy'):
                    if redtype: self.printLog('#RED','%s: %s files' % (redtype,rje.iLen(fdb.index('Redundancy')[redtype])))
                    else: self.printLog('#RED','Unique: %s files' % (rje.iLen(fdb.index('Redundancy')[redtype])))                                             
                return
                if self.opt['Monster']:
                    if file not in skiplist and filedict.has_key(filename):  # Potential redundancy!
                        # 1. Total redundancy #
                        if filedict[filename].has_key(details['FilePath']):     
                            if filedict[filename].has_key[details['FilePath']]['Age'] != age:
                                self.log.errorLog('File "%s" appears twice with different ages!' % details['FilePath'],printerror=False)
                            if filedict[filename].has_key[details['FilePath']]['Size'] != size:
                                self.log.errorLog('File "%s" appears twice with different sizes!' % details['FilePath'],printerror=False)
                            overlap += 1
                            #continue
                        # 2. Duplicate files in different locations? #
                        details['Redundant'] = ''
                        moved = False
                        for prevfile in filedict[filename].keys()[0:]:
                            # Setup comparison #
                            reversed = False    # Whether files are "reversed" due to good vs. bad comparison
                            prevage = filedict[filename][prevfile]['Age']
                            prevsize = filedict[filename][prevfile]['Size']
                            good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                            prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                            if (prevgood == 'bad' and good != 'bad') or (good == 'good' and prevgood != 'good'):
                                reversed = True
                            elif self.opt['KeepNew'] and prevage > age:   # Previous file newer - reverse
                                reversed = True
                                    
                            # ID Type of redundancy #
                            if prevage == age and prevsize == size:
                                details['Redundant'] = 'Duplicate'
                            elif prevage < age and prevsize == size:
                                details['Redundant'] = 'Newer copy'
                            elif prevage > age and prevsize == size:
                                details['Redundant'] = 'Older copy'
                            elif prevage > age and prevsize > size:
                                details['Redundant'] = 'Older version? (Smaller)'
                            elif prevage < age and prevsize < size:
                                details['Redundant'] = 'Newer version? (Bigger)'

                            # Choose default option #
                            options = '<0> Keep both; Move <1> or <2>; Add <G>oodDir or <B>adDir; Add to <S>kipList; <P>urge; <Q>uit'
                            autochoose = False
                            default = self.moveDefault(good,prevgood,details['Redundant'])
                            if prevgood == 'good' and good == 'good':
                                autochoose = True
                            elif details['Redundant'] in ['Duplicate','Newer copy','Older copy'] and good == 'bad':
                                autochoose = True
                            elif details['Redundant'] in ['Duplicate','Newer copy','Older copy'] and prevgood == 'bad':
                                autochoose = True
                            
                            # Give and execute options #
                            if details['Redundant']:
                                choicestr = (self.reducedPath(details['FilePath']),good,self.reducedPath(prevfile),prevgood,good,details['Redundant'],prevgood,options)
                                choicetxt = '\n\n<1> %s [%s]\n<2> %s [%s]\n\n\t=> <1> [%s] is %s of <2> [%s]\n\n%s' % choicestr
                                if self.stat['Interactive'] < 0 or (self.stat['Interactive'] < 1 and autochoose):
                                    choice = default
                                else:
                                    choice = 'X'
                                    while choice not in ['0','1','2','Q','S']:
                                        choice = rje.choice(choicetxt,default).upper()
                                        while choice == 'G':    # Add GoodDir
                                            print '\nAdd directory as "Good":\n<0> Cancel'
                                            print '<1>', os.path.dirname(self.reducedPath(details['FilePath']))
                                            print '<2>', os.path.dirname(self.reducedPath(prevfile))
                                            choice = rje.choice('Select dir [0/1/2]','0')
                                            if choice == '0':
                                                choice = 'X'    # Exit this but not full choice
                                            elif choice == '1':
                                                self.addGoodDir(os.path.dirname(details['FilePath']))
                                                good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                                                prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                                                default = self.moveDefault(good,prevgood,details['Redundant'])
                                                choice = 'X'
                                            elif choice == '2':
                                                self.addGoodDir(os.path.dirname(prevfile))
                                                good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                                                prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                                                default = self.moveDefault(good,prevgood,details['Redundant'])
                                                choice = 'X'
                                            else:
                                                choice == 'G'   # Unallowed choice!
                                        while choice == 'B':    # Add BadDir
                                            print '\nAdd directory as "Bad":\n<0> Cancel'
                                            print '<1>', os.path.dirname(self.reducedPath(details['FilePath']))
                                            print '<2>', os.path.dirname(self.reducedPath(prevfile))
                                            choice = rje.choice('Select dir [0/1/2]','0')
                                            if choice == '0':
                                                choice = 'X'    # Exit this but not full choice
                                            elif choice == '1': # Keep for original menu - move!
                                                self.addBadDir(os.path.dirname(details['FilePath']))
                                                if self.stat['Interactive'] >= 1:
                                                    default = '1'
                                                    choice = 'X'
                                            elif choice == '2': # Keep for original menu - move!
                                                self.addBadDir(os.path.dirname(prevfile))
                                                if self.stat['Interactive'] >= 1:
                                                    default = '2'
                                                    choice = 'X'
                                            else:
                                                choice == 'B'   # Unallowed choice!
                                        if choice == 'P':
                                            if self.purgeList(filename):
                                                choice = '0'
                                        if choice == 'S':   #!# Does not work yet! #!#
                                            addskip = filename
                                            if rje.yesNo('Skip all *%s files?' % os.path.splitext(filename)[1],default='N'):
                                                addskip = os.path.splitext(filename)[1]
                                            self.log.printLog('#SKIP','Added %s to skiplist' % filename)
                                            for folder in self.list['DirList']:
                                                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                                                self.verbose(0,3,'\nScavenging %s files to skip:' % folder,1)
                                                skiplist += rje.getFileList(self,folder,addskip,self.opt['SubFolders'])
                                            self.log.printLog('\n#SKIP','%s files to skip from %d top level directories.' % (rje.integerString(len(skiplist)),len(self.list['DirList'])))
                                            choice = '0'
                                            
                                # Execute choice #
                                if choice == '1':   # Moving new file
                                    moved = self.monsterMove(details['FilePath'])
                                    if moved:
                                        movex += 1
                                        details['Redundant'] = details['Redundant'] + ' (Moved)'
                                        break
                                elif choice == '2': # Moving old file
                                    moved = self.monsterMove(prevfile)
                                    if moved:
                                        movex += 1
                                    filedict[filename].pop(prevfile)
                                elif choice == 'S':
                                    break
                                elif choice == 'Q' and rje.yesNo('Quit Monster?'):
                                    return

                        # 3. Same name but different size #
                        if not details['Redundant']:
                            details['Redundant'] = 'Size conflict!'
                        if os.path.exists(details['FilePath']):
                            filedict[filename][details['FilePath']] = {'Age':age,'Size':size}
                    else:                        
                        filedict[filename] = {details['FilePath']:{'Age':age,'Size':size}}
            
        except: self.errorLog('%s.redundancy error' % self)
#########################################################################################################################
    def fixEnding(self,file):    ### Fixes Mac -> UNIX line endings for file
        '''Fixes Mac -> UNIX line endings for file.'''
        try:
            fixed = string.replace(open(file,'r').read(),'\r','\n')
            rje.backup(self,file)
            open(file,'w').write(fixed)
            return True
        except: return False
#########################################################################################################################
### End of SECTION II: FileMonster Class                                                                                #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: FileScout Class:                                                                                        #
#########################################################################################################################
class FileScout(rje.RJE_Object):    ### Searches for files
    '''
    FileScout Class. Author: Rich Edwards (2005).

    Info:str
    - Name = ID for search - allows multiple searches to be compared easily [None]
    - Folder = Folder to analyse [./]
    - OutFile = File for output [file_monster.txt]
    - DumpDir = Directory in which to dump redundant files [Redundant]
    - OutDir = Output directory for renamed files [./]
    - SortBy = Whether to sort by date or name [date]
    - PreFix = Text prefix for new file names []
    
    Opt:boolean
    - Append = Whether to append to output file [False]
    - SubFolders = Whether to look in subfolders [True]
    - Scavenge = Whether to perform collation of file information [False]
    - DirSum = Whether to perform summary of directory contents [False]
    - Monster = Whether to perform monster cleanup of redundant files [False]
    - KeepNew = Preferentially keep newer files of same size if good/bad status equal [True]
    - Rename = Whether to rename files [False]
    - UseDate = Whether to use date in new name [False]

    Stat:numeric
    - DirCut = Max number of subdirectories to have and still delve into them [20]

    List:list
    - FileList = List of file extensions of interest (X,Y,..,Z)  [*]
    - DirList = List of directories to look in, in order of preference good -> bad. Over-rules folder=PATH []
    - GoodDir = List of "good" directories to be given priority in decision making (Can be partial paths) []
    - BadDir = List of "bad" directories to be automatically screened if i<1 (Can be partial paths) []
    - SkipList = List of filenames to skip [thumbs.db]
    - PurgeList = List of filenames (allowing wildcards) to purge (move/delete) [WS_FTP.LOG]
    
    Dict:dictionary    

    Obj:RJE_Objects
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Folder','OutFile','DumpDir','OutDir','SortBy','PreFix']
        self.optlist = ['SubFolders','Append','Monster','DirSum','Scavenge','KeepNew','Rename','UseDate']
        self.statlist = ['DirCut']
        self.listlist = ['FileList','DirList','GoodDir','BadDir','SkipList','PurgeList']
        self.dictlist = []
        self.objlist = []
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Folder':rje.makePath('./',return_blank=False),'OutDir':rje.makePath('',return_blank=False),
                      'OutFile':'file_monster.csv','DumpDir':rje.makePath('Redundant'),'SortBy':'date'})
        self.setOpt({'SubFolders':True,'KeepNew':True})
        self.setStat({'DirCut':20})
        self.list['FileList'] = ['*']
        self.list['DirList'] = ['']
        self.list['SkipList'] = ['Thumbs.db']
        self.list['PurgeList'] = ['WS_FTP.LOG']
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdRead(cmd,type='info',att='Name',arg='searchid')  
                self._cmdReadList(cmd,'path',['Folder','DumpDir'])  
                self._cmdReadList(cmd,'info',['Name','PreFix','SortBy'])  
                self._cmdReadList(cmd,'file',['OutFile'])  
                self._cmdReadList(cmd,'opt',['SubFolders','Append','Monster','DirSum','Scavenge','KeepNew','Rename','UseDate'])  
                self._cmdRead(cmd,type='int',att='DirCut')
                self._cmdReadList(cmd,'list',['FileList','DirList','GoodDir','BadDir','SkipList','PurgeList'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def run(self):      ### Main run control method
        '''Main run control method.'''
        if self.opt['DirSum']: self.dirSum()
        if self.opt['Scavenge'] or self.opt['Monster']: self.monster()
        if self.opt['Rename']: self.rename()
#########################################################################################################################
    def monster(self):      ### Goes through directory (and subdirectories) and scavenges files. Yum.
        '''
        Goes through directory (and subdirectories) and scavenges files. Yum.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            ## General setup ##
            delimit = rje.getDelimit(self.cmd_list,',')
            if not self.list['DirList']:
                self.list['DirList'] = [self.info['Folder']]
            filelist = []
            if self.info['OutFile'].lower() in ['none','']:
                self.info['OutFile'] = ''
            if self.opt['Monster']:
                if self.info['DumpDir'].lower() in ['none','']:
                    self.info['DumpDir'] = ''
                elif not os.path.exists(self.info['DumpDir']):
                    os.mkdir(self.info['DumpDir'])
                if not os.path.exists(self.info['DumpDir']):
                    self.info['DumpDir'] = ''       # Mark redundant in scavenge only
            headers = ['Search','File','Folder','Path','Type','Size','Accessed','Modified']
            if self.opt['Monster']:
                headers.append('Redundant')
            ## File Dictionary ##
            filedict = {}   # Dictionary of files: {file:{fullpath:{age,size}}}
            overlap = 0     # Number of files duplicated in filelist due to overlapping input directories
            movex = 0       # Number of redundant files moved to the DumpDir

            ### PurgeList ###
            self.purgeList()

            ### Scavenging Time! ###
            _stage = 'Scavenge'
            for folder in self.list['DirList']:
                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                self.verbose(0,3,'\nScavenging %s files:' % folder,1)
                filelist += rje.getFileList(self,folder,self.list['FileList'],self.opt['SubFolders'])
            self.log.printLog('\n#FILES','%s files scavenged from %d top level directories.' % (rje.integerString(len(filelist)),len(self.list['DirList'])))

            ### SkipList ###
            skiplist = []
            for folder in self.list['DirList']:
                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                self.verbose(0,3,'\nScavenging %s files to skip:' % folder,1)
                skiplist += rje.getFileList(self,folder,self.list['SkipList'],self.opt['SubFolders'])
            self.log.printLog('\n#SKIP','%s files to skip from %d top level directories.' % (rje.integerString(len(skiplist)),len(self.list['DirList'])))
            
            ### Setup Output ###
            _stage = 'Setup output'
            if self.opt['Append'] and self.info['OutFile']:
                self.verbose(0,3,'Appending %s file details to %s...' % (rje.integerString(len(filelist)),self.info['OutFile']),1)
                OUTFILE = open(self.info['OutFile'],'a')
            elif self.info['OutFile']:
                rje.backup(self,self.info['OutFile'])
                self.verbose(0,3,'Writing %s file details to %s...' % (rje.integerString(len(filelist)),self.info['OutFile']),1)
                rje.delimitedFileOutput(self,self.info['OutFile'],headers,delimit)

            ### Process Files ###
            _stage = 'Process'
            fx = 0.0    # Progress counter
            dirx = 0    # Number of directories
            for file in filelist:
                ## Log ##
                fx += 100.00
                logtxt = 'Checking %s files: %.2f%% (%d dir); %s moved (%s overlaps)' % (rje.integerString(len(filelist)),fx/len(filelist),dirx,rje.integerString(movex),rje.integerString(overlap))
                self.log.printLog('\r#CLEAN',logtxt,log=False,newline=False)
                ## Skip directory ##
                if os.path.isdir(file):
                    dirx += 1
                    continue
                ## Get file details ##
                details = {'Search':self.info['Name'],'Size':0,'Accessed':'???','Modified':'???'}
                details['File'] = os.path.basename(file)                            # Basic file name without path
                filename = details['File']
                details['FilePath'] = os.path.abspath(file)                         # Full path to file
                details['Path'] = os.path.dirname(details['FilePath'])              # Full path to folder (no file)
                details['Folder'] = string.split(details['Path'])[-1]               # Immediate folder containing file
                details['Type'] = os.path.splitext(file)[1]                         # File type (extension)
                try:    # File details - size and age
                    size = os.path.getsize(file)
                    age = os.path.getmtime(file)
                    details['Size'] = '%d' % float(os.path.getsize(file))
                    details['Accessed'] = time.ctime(os.path.getatime(file))
                    details['Modified'] = time.ctime(os.path.getmtime(file))
                except:
                    self.log.errorLog('Error in scavenge %s size/date.' % file,printerror=True,quitchoice=False)

                ## Add to dictionary and check redundancy ##
                if self.opt['Monster']:
                    if file not in skiplist and filedict.has_key(filename):  # Potential redundancy!
                        # 1. Total redundancy #
                        if filedict[filename].has_key(details['FilePath']):     
                            if filedict[filename].has_key[details['FilePath']]['Age'] != age:
                                self.log.errorLog('File "%s" appears twice with different ages!' % details['FilePath'],printerror=False)
                            if filedict[filename].has_key[details['FilePath']]['Size'] != size:
                                self.log.errorLog('File "%s" appears twice with different sizes!' % details['FilePath'],printerror=False)
                            overlap += 1
                            continue
                        # 2. Duplicate files in different locations? #
                        details['Redundant'] = ''
                        moved = False
                        for prevfile in filedict[filename].keys()[0:]:
                            # Setup comparison #
                            reversed = False    # Whether files are "reversed" due to good vs. bad comparison
                            prevage = filedict[filename][prevfile]['Age']
                            prevsize = filedict[filename][prevfile]['Size']
                            good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                            prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                            if (prevgood == 'bad' and good != 'bad') or (good == 'good' and prevgood != 'good'):
                                reversed = True
                            elif self.opt['KeepNew'] and prevage > age:   # Previous file newer - reverse
                                reversed = True
                                    
                            # ID Type of redundancy #
                            if prevage == age and prevsize == size:
                                details['Redundant'] = 'Duplicate'
                            elif prevage < age and prevsize == size:
                                details['Redundant'] = 'Newer copy'
                            elif prevage > age and prevsize == size:
                                details['Redundant'] = 'Older copy'
                            elif prevage > age and prevsize > size:
                                details['Redundant'] = 'Older version? (Smaller)'
                            elif prevage < age and prevsize < size:
                                details['Redundant'] = 'Newer version? (Bigger)'

                            # Choose default option #
                            options = '<0> Keep both; Move <1> or <2>; Add <G>oodDir or <B>adDir; Add to <S>kipList; <P>urge; <Q>uit'
                            autochoose = False
                            default = self.moveDefault(good,prevgood,details['Redundant'])
                            if prevgood == 'good' and good == 'good':
                                autochoose = True
                            elif details['Redundant'] in ['Duplicate','Newer copy','Older copy'] and good == 'bad':
                                autochoose = True
                            elif details['Redundant'] in ['Duplicate','Newer copy','Older copy'] and prevgood == 'bad':
                                autochoose = True
                            
                            # Give and execute options #
                            if details['Redundant']:
                                choicestr = (self.reducedPath(details['FilePath']),good,self.reducedPath(prevfile),prevgood,good,details['Redundant'],prevgood,options)
                                choicetxt = '\n\n<1> %s [%s]\n<2> %s [%s]\n\n\t=> <1> [%s] is %s of <2> [%s]\n\n%s' % choicestr
                                if self.stat['Interactive'] < 0 or (self.stat['Interactive'] < 1 and autochoose):
                                    choice = default
                                else:
                                    choice = 'X'
                                    while choice not in ['0','1','2','Q','S']:
                                        choice = rje.choice(choicetxt,default).upper()
                                        while choice == 'G':    # Add GoodDir
                                            print '\nAdd directory as "Good":\n<0> Cancel'
                                            print '<1>', os.path.dirname(self.reducedPath(details['FilePath']))
                                            print '<2>', os.path.dirname(self.reducedPath(prevfile))
                                            choice = rje.choice('Select dir [0/1/2]','0')
                                            if choice == '0':
                                                choice = 'X'    # Exit this but not full choice
                                            elif choice == '1':
                                                self.addGoodDir(os.path.dirname(details['FilePath']))
                                                good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                                                prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                                                default = self.moveDefault(good,prevgood,details['Redundant'])
                                                choice = 'X'
                                            elif choice == '2':
                                                self.addGoodDir(os.path.dirname(prevfile))
                                                good = self.goodStatus(details['FilePath'])     # Whether current file is in good/bad/neutral dir
                                                prevgood = self.goodStatus(prevfile)            # Whether previous file is in good/bad/neutral dir
                                                default = self.moveDefault(good,prevgood,details['Redundant'])
                                                choice = 'X'
                                            else:
                                                choice == 'G'   # Unallowed choice!
                                        while choice == 'B':    # Add BadDir
                                            print '\nAdd directory as "Bad":\n<0> Cancel'
                                            print '<1>', os.path.dirname(self.reducedPath(details['FilePath']))
                                            print '<2>', os.path.dirname(self.reducedPath(prevfile))
                                            choice = rje.choice('Select dir [0/1/2]','0')
                                            if choice == '0':
                                                choice = 'X'    # Exit this but not full choice
                                            elif choice == '1': # Keep for original menu - move!
                                                self.addBadDir(os.path.dirname(details['FilePath']))
                                                if self.stat['Interactive'] >= 1:
                                                    default = '1'
                                                    choice = 'X'
                                            elif choice == '2': # Keep for original menu - move!
                                                self.addBadDir(os.path.dirname(prevfile))
                                                if self.stat['Interactive'] >= 1:
                                                    default = '2'
                                                    choice = 'X'
                                            else:
                                                choice == 'B'   # Unallowed choice!
                                        if choice == 'P':
                                            if self.purgeList(filename):
                                                choice = '0'
                                        if choice == 'S':   #!# Does not work yet! #!#
                                            addskip = filename
                                            if rje.yesNo('Skip all *%s files?' % os.path.splitext(filename)[1],default='N'):
                                                addskip = os.path.splitext(filename)[1]
                                            self.log.printLog('#SKIP','Added %s to skiplist' % filename)
                                            for folder in self.list['DirList']:
                                                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                                                #??#
                                                self.verbose(0,3,'\nScavenging %s files to skip:' % folder,1)
                                                skiplist += rje.getFileList(self,folder,addskip,self.opt['SubFolders'])
                                            self.log.printLog('\n#SKIP','%s files to skip from %d top level directories.' % (rje.integerString(len(skiplist)),len(self.list['DirList'])))
                                            choice = '0'
                                            
                                # Execute choice #
                                if choice == '1':   # Moving new file
                                    moved = self.monsterMove(details['FilePath'])
                                    if moved:
                                        movex += 1
                                        details['Redundant'] = details['Redundant'] + ' (Moved)'
                                        break
                                elif choice == '2': # Moving old file
                                    moved = self.monsterMove(prevfile)
                                    if moved:
                                        movex += 1
                                    filedict[filename].pop(prevfile)
                                elif choice == 'S':
                                    break
                                elif choice == 'Q' and rje.yesNo('Quit Monster?'):
                                    return

                        # 3. Same name but different size #
                        if not details['Redundant']:
                            details['Redundant'] = 'Size conflict!'
                        if os.path.exists(details['FilePath']):
                            filedict[filename][details['FilePath']] = {'Age':age,'Size':size}
                    else:                        
                        filedict[filename] = {details['FilePath']:{'Age':age,'Size':size}}

                ## Output ##                
                #X#rje.writeDelimit(OUTFILE,[self.info['Name'],outfile,outfolder.lower(),outpath,outtype,outsize,outatime,outmtime],delimit)
                rje.delimitedFileOutput(self,self.info['OutFile'],headers,delimit,details)
            logtxt = 'Checking %s files complete: %s moved (%s overlaps)' % (rje.integerString(len(filelist)),rje.integerString(movex),rje.integerString(overlap))
            self.log.printLog('\r#CLEAN',logtxt)
                
            
        except:
            self.log.errorLog('Error in monster()',printerror=True,quitchoice=True)
#########################################################################################################################
    def moveDefault(self,good,prevgood,redundant):  ### Sets the default choice for file move options
        '''
        Sets the default choice for file move options.
        >> good:str = "Good" status of the current file being considered (file 1)
        >> prevgood:str = "Good" status of previous file being compared to (file 2)
        >> redundant:str = Redundancy status of file 1 vs. file 2
        '''
        try:
            ### Setup ###
            default = '0'
            if prevgood == 'good' and good == 'good':
                default = '0'
            elif redundant in ['Duplicate','Newer copy','Older copy','Older version? (Smaller)'] and good == 'bad':
                default = '1'
            elif redundant in ['Duplicate','Newer copy','Older copy','Newer version? (Bigger)'] and prevgood == 'bad':
                default = '2'
            elif redundant in ['Duplicate','Newer copy','Older copy','Older version? (Smaller)'] and prevgood == 'good':
                default = '1'
            elif redundant in ['Duplicate','Newer copy','Older copy','Newer version? (Bigger)'] and good == 'good':
                default = '2'
            elif redundant in ['Duplicate','Newer copy']:
                default = '1'
            elif redundant == 'Older copy' and self.opt['KeepNew']:
                default = '2'
            elif redundant == 'Older copy':
                default = '1'
            return default
        except:
            self.log.errorLog('Problem with moveDefault()')
            return default
#########################################################################################################################
    def purgeList(self,filename=None,clear=True):      ### Purges files from directories
        '''
        Purges files from directories.
        >> filename:str = file to consider
        >> clear:boolean = whether to clear PurgeList after purge [True]
        '''
        try:
            ### Add new PurgeList ###
            purgedfile = False
            if not self.list['PurgeList']:     # Menu to add
                if not filename: return False
                addlist = ['Cancel',filename,'*.%s' % os.path.splitext(filename)[1],'Other','Purge']
                default = 0
                choice = ''
                while choice != 'Purge':
                    print '\nAdd file(s) to "PurgeList":'
                    for i in range(len(addlist)):
                        print '<%d> %s' % (i, addlist[i])
                    x = rje.getInt('Choice?',default=default)
                    if x < len(addlist):
                        choice = addlist[i]
                    if choice == 'Cancel':
                        return False
                    elif choice == 'Other':
                        choice = rje.choice('Filename to add to PurgeList:',default='',confirm=True)
                    if choice and choice != 'Purge':
                        self.list['PurgeList'].append(choice)
                        if choice in addlist:
                            addlist.remove(choice)
                            purgedfile = True

            ### Purge/Delete Option ###
            delete = False
            if self.stat['Interactive'] >= 0 and rje.yesNo('Delete purged files? (No = move only.)',default='N',confirm=True):
                delete = True
                
            ### Purge files ###
            purgelist = []
            for folder in self.list['DirList']:
                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                self.verbose(0,3,'\nScavenging %s files:' % folder,1)
                purgelist += rje.getFileList(self,folder,self.list['PurgeList'],self.opt['SubFolders'])
            self.log.printLog('\n#FILES','%s files scavenged from %d top level directories.' % (rje.integerString(len(purgelist)),len(self.list['DirList'])))
            if delete and (self.stat['Interactive'] < 0 or rje.yesNo('Delete %s files. Are you sure?' % rje.integerString(len(purgelist)),default='Y')):
                for file in purgelist:
                    os.unlink(file)
            else:
                for file in purgelist:
                    self.monsterMove(os.path.abspath(file))
            if clear:
                self.list['PurgeList'] = []
                                               
            ### End ###
            return purgedfile
        except:
            self.log.errorLog('Error in purgeList()')
            return False
#########################################################################################################################
    def addBadDir(self,path):   ### Adds path to BadDir
        '''Adds path to BadDir.'''
        el = string.split(path, os.sep)
        i = len(el)
        while el:
            print 'Add BadDir: %s' % rje.makePath(string.join(el[:i],'/'))
            choice = rje.choice('<> to lengthen/shorten or Enter to accept?',default='')
            if choice == '>' and i < len(el):
                i += 1
            elif choice == '<' and i > 1:
                i -= 1
            elif choice == '':
                self.list['BadDir'].append(rje.makePath(string.join(el[:i],'/')))
                return
#########################################################################################################################
    def addGoodDir(self,path):   ### Adds path to GoodDir
        '''Adds path to GoodDir.'''
        el = string.split(path, os.sep)
        i = len(el)
        while el:
            print 'Add GoodDir: %s' % rje.makePath(string.join(el[:i],'/'))
            choice = rje.choice('<> to lengthen/shorten or Enter to accept?',default='')
            if choice == '>' and i < len(el):
                i += 1
            elif choice == '<' and i > 1:
                i -= 1
            elif choice == '':
                self.list['GoodDir'].append(rje.makePath(string.join(el[:i],'/')))
                return
#########################################################################################################################
    def goodStatus(self,filepath):      ### Returns good/bad/neutral
        '''Returns good/bad/neutral.'''
        for dir in self.list['GoodDir']:
            if filepath.find(dir) >= 0:
                return 'good'
        for dir in self.list['BadDir']:
            if filepath.find(dir) >= 0:
                return 'bad'
        return 'neutral'
#########################################################################################################################
    def reducedPath(self,filepath):     ### Returns a reduced filepath on basis of self.list['DirList']
        '''Returns a reduced filepath on basis of self.list['DirList']'''
        try:
            filedirs = string.split(filepath,os.sep)
            for folder in self.list['DirList']:
                fulldirs = string.split(os.path.abspath(rje.makePath(folder)),os.sep)
                dirmatch = True
                for dir in fulldirs:
                    if dir not in filedirs:
                        dirmatch = False
                #X#print fulldirs, filedirs, dirmatch
                if dirmatch:
                    filedirs = string.split(folder,'/') + filedirs[len(fulldirs):]
                    break
            return rje.makePath(string.join(filedirs,'/'),wholepath=True)
        except:
            self.log.errorLog('Error with reducedPath()')
            return filepath
#########################################################################################################################
    def monsterMove(self,filepath):     ### Moves file into self.info['DumpDir']
        '''Moves file into self.info['DumpDir'].'''
        try:
            ### Setup ###
            if not self.info['DumpDir']:
                return False
            filename = os.path.basename(filepath)

            ### Try straight ###
            newfile = rje.makePath(self.info['DumpDir'] + filename,wholepath=True)
            if not os.path.exists(newfile):
                os.rename(filepath,newfile)
                self.log.printLog('#MOVE','%s moved to %s' % (filepath,newfile),screen=False)
                return True

            ### Try subdirs ###            
            extra = 1
            while extra:
                dump = rje.makePath('%s%d' % (self.info['DumpDir'],extra))
                if not os.path.exists(dump):
                    os.mkdir(dump)
                newfile = rje.makePath(dump + filename,wholepath=True)
                if not os.path.exists(newfile):
                    os.rename(filepath,newfile)
                    self.log.printLog('#MOVE','%s moved to %s' % (filepath,newfile),screen=False)
                    return True
                extra += 1
                if extra == 100:
                    raise ValueError
        except:
            self.log.errorLog('Problem with monsterMove(%s)' % filepath,quitchoice=True)
            return False
#########################################################################################################################
    def dirSum(self):   ### Summarises number of files and subdirectories
        '''Summarises number of files and subdirectories.'''
        try:

            ### DirSum ###
            my_dirs = [self.info['Folder']]     # Add subfolders here to look in
            d = 0
            while d < len(my_dirs):     # Continue looping
                folder = my_dirs[d]
                d += 1
                ## Count subdir ##
                subdir = []
                try:
                    for file in os.listdir(folder):     # Get full list of directory contents
                        fullfile = os.path.join(folder,file)
                        if os.path.isdir(fullfile):     # Found a subdirectory
                            subdir.append(fullfile)
                except: self.errorLog('Problem with folder "%s" during FileScout.dirSum()' % folder); continue
                    
                ## Count Files ##
                filex = ''
                for file in self.list['FileList']:
                    if '*' in file: filex += ' %5s %s;' % (rje.integerString(len(glob.glob(os.path.join(folder,file)))),file)
                    else: filex += ' %5s %s;' % (rje.integerString(len(glob.glob(os.path.join(folder,'*.%s' % file)))),file)
                ## Output ##
                self.log.printLog('#DIR','%s\t%s %5s dir.' % (folder,filex,rje.integerString(len(subdir))))
                if len(subdir) <= self.stat['DirCut']:
                    my_dirs = my_dirs[:d] + subdir + my_dirs[d:]
            
        except:
            self.log.errorLog('Major error in FileScout.dirSum()',quitchoice=True)
            return False
#########################################################################################################################
    def rename(self):   ### Scavenges list of files and then renames them into a new directory
        '''Scavenges list of files and then renames them into a new directory.'''
        try:
            ### ~ [1] Get list of files to rename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            filelist = []
            for folder in self.list['DirList']:
                folder = rje.makePath(folder,wholepath=True,return_blank=False)
                self.verbose(0,3,'\nScavenging "%s" files:' % folder,1)
                print folder,self.list['FileList'],self.opt['SubFolders']
                filelist += rje.getFileList(self,folder,self.list['FileList'],self.opt['SubFolders'])
            self.log.printLog('\n#FILES','%s files scavenged from %d top level directories.' % (rje.integerString(len(filelist)),len(self.list['DirList'])))
            filelist.sort()
            ## ~ [1a] Sort according to name/date ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sorter = {}
            #x#self.deBug(filelist[:10])
            for file in filelist:
                #x#print file, os.path.getmtime(file), rje.baseFile(file,strip_path=True)
                if self.info['SortBy'] == 'date':
                    skey = os.path.getmtime(file)
                    while skey in sorter: skey += 0.0001     # Resolve ties
                else:
                    skey = rje.baseFile(file,strip_path=True)
                    s = 0
                    while skey in sorter: (skey,s) = ('%s%d' % (rje.baseFile(file,strip_path=True),s), s + 1)
                sorter[skey] = file
                #self.deBug(sorter)

            ### ~ [2] Output new file names ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Make dictionary of new prefixes and old files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfiles = {}       # Dictionary of {new prefix:[old files]}
            for skey in rje.sortKeys(sorter):
                file = sorter[skey]
                prefix = self.info['PreFix']
                if self.opt['UseDate']:
                    age = os.path.getmtime(file)
                    t = time.gmtime(age)
                    prefix = '%s%s%s%s' % (prefix,str(t[0]),rje.preZero(t[1],12),rje.preZero(t[2],31))
                if prefix not in outfiles: outfiles[prefix] = []
                self.deBug('%s :: %s :: %s' % (skey,file,prefix))
                #self.deBug(outfiles[prefix])
                outfiles[prefix].append(file)
                #self.deBug(outfiles[prefix])
            ## ~ [2b] Rename files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.mkDir(self,self.info['OutDir'])
            for prefix in rje.sortKeys(outfiles):
                prex = len(outfiles[prefix])
                if self.stat['Interactive'] >= 0 and not rje.yesNo('Rename %s %s files to %s?' % (rje.integerString(prex),prefix,self.info['OutDir'])): continue
                for f in range(prex):
                    file = outfiles[prefix][f]
                    ext = os.path.splitext(file)[-1]
                    newfile = '%s%s_%s%s' % (self.info['OutDir'],prefix,rje.preZero(f+1,prex),ext)
                    print '%s -> %s' % (file,newfile)
                    os.rename(file,newfile)
                self.log.printLog('#FILE','%s %s files renamed to %s' % (rje.integerString(prex),prefix,self.info['OutDir']))

        except: self.log.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
### End of SECTION II: FileScout                                                                                        #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
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
    try: FileMonster(mainlog,cmd_list).run()

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
