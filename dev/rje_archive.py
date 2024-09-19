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
Module:       rje_archive
Description:  KDM Archive Manager
Version:      0.7.3
Last Edit:    28/04/20
Copyright (C) 2017  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is for backing up data to the UNSW Research Data Store (RDS) and reporting on the status of backups.
    Details to be added.

Development notes:
    This module is designed to run in the following modes:

    1. Backup. This archives a specific directory (or list) to a given data archive project.

    2. Remove. This checks that directories are in the archive and then deletes them according to some criteria.

    `du testing` -> parse `^(\d+)\s_\S.+$' -> size, directory

    ```
    ls -lt | head -2
    total 60
    drwxr-xr-x.  3 z3452659 unsw 4096 Mar 30 16:40 tools
    ```

    -> Store this for each directory: should not change between upload check.

    `ls -1p dev | grep "/$" -cv` -> Number of files
    `ls -1p dev | grep "/$" -c` -> Number of directories

    ```
    module add unswdataarchive/2015-09-10
    upload.sh /home/z3452659/bioinf/redwards/projects/ManefieldPacBio-Sep15 "/UNSW_RDS/D0234444/"
    ```

    NOTE: Will need to get the system set up so that it does not ask for a password. This can be achieved by requesting
    a token from IT and then updating the config file.

    The first time it runs, it will report X imported file(s). If the files are already uploaded, it will report 0 files
    imported. All files are listed in "Consume" lines, though. This should (presumably) match the number of files listed
    by ls with the -a flag, ignoring '.' and '..'? Or could try matching against:
    rje.listDir(callobj=None,folder=os.getcwd(),subfolders=True,folders=True,files=True,summary=True,asksub=False,dircut=0,dirdepth=-1)

    => Make checking file numbers a toggle. (checknum=T/F)

    $ upload.sh /srv/scratch/z3452659/CaneToad-May15/analysis/2016-11-22.Tyr "/UNSW_RDS/D0234445/CaneToad-May15/analysis"
    Picked up _JAVA_OPTIONS: -Xmx1g
    Password:
    Consume: 2016-11-22.Tyr/BLASTFAS/F7CL37.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/canetoad.20161122A.tyr_XENTR.vs.canetoad.20161122A.est_hits.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/fiesta.ini [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/fiesta.log [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/gablam.log [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nhr [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nin [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nog [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsd [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsi [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsq [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.gablam.tdt [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.hitsum.tdt [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.local.tdt [construct: null, logical: null, encapsulation: null]
    live: imported 15 file(s)

    $ upload.sh /srv/scratch/z3452659/CaneToad-May15/analysis/2016-11-22.Tyr "/UNSW_RDS/D0234445/CaneToad-May15/analysis"
    Picked up _JAVA_OPTIONS: -Xmx1g
    Password:
    Consume: 2016-11-22.Tyr/canetoad.20161122A.tyr_XENTR.vs.canetoad.20161122A.est_hits.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/BLASTFAS/F7CL37.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/fiesta.ini [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/fiesta.log [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/gablam.log [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nhr [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nin [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nog [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsd [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsi [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.fas.nsq [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.gablam.tdt [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.hitsum.tdt [construct: null, logical: null, encapsulation: null]
    Consume: 2016-11-22.Tyr/tyr_XENTR.vs.canetoad.20161122A.local.tdt [construct: null, logical: null, encapsulation: null]
    live: imported 0 file(s)

    #!# Should generate a log file that has the time, total number of files and total number imported.

    .name' failed: The namespace '/UNSW_RDS/D0234445/CaneToad-May15/delete' does not exist or is not accessible

Output:
    Add details of backups.tdt and archived.tdt here: (dir, project, rds, date, files, imports)

Commandline:
    ### ~ Main Archive Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    rds=X           : UNSW_RDS ResData ID project code to use (e.g. D0234444) []
    uploadsh=X      : Full path for runnining upload.sh script ['/home/z3452659/unswdataarchive/upload.sh']
    homedir=PATH    : Home directory from which the archive script will be run ['~']
    projects=FILE   : Delimited file of `Project` and `RDS` code. If provided, will not use rds=X ['projects.tdt']
    strict=T/F      : Restrict processing to projects found in Projects file and add no new ones. [False]
    backupdirs=LIST : List of directories to backup (should be project subdirectory full paths) []
    archivedirs=LIST: List of directories to check archive and tar/delete (should be project subdirectory full paths) []
    rmdirs=T/F      : Delete archived directories. (Will ask if i>0) [False]
    targz=T/F       : Whether to tar and zip directories to be deleted [True]
    checknum=T/F    : Whether to check numbers of files consumed by upload.sh versus directory contents [True]
    tryparent=T/F   : Whether to try to run backup parent directory in case of failure [False]
    basefile=FILE   : This will set the 'root' filename for output files (FILE.*), including the log ['rds']
    backupdb=FILE   : File to output backup summaries into  ['BASEFILE.backups.tdt']
    archived=FILE   : File to output archive summaries into  ['BASEFILE.archived.tdt']
    cleanup=T/F     : Whether to perform post-upload cleanup of backups and archived files [True]
    quiet=X         : Min number of days of inactivity before a directory gets rates as quiet [1]
    skipquiet=T/F   : Whether to skip uploads for quiet directories [True]
    checkarchive=T/F: Whether to run upload.sh on on directories that have been uploaded in last run [False]
    dormancy=X      : Min number of days of inactivity before a directory gets rated as dormant (0=no dormancy) [30]
    skipdormant=T/F : Whether to skip uploads for dormant directories [True]
    dormant=FILE    : File to output dormant directories into  ['BASEFILE.dormant.tdt']
    maxfiles=INT    : Maximum number of files in a directory to generate backups (0 = no limit) [10000]
    maxdirsize=INT  : Maximum directory size in bytes to generate backup (0 = no limit) [1e11 (~100Gb)]
    archivetgz=T/F  : Whether to tar and zip then backup directories exceeding maxfiles=INT cutoff [False]
    useqsub=T/F     : Whether to use QSub for tarballing and then archiving the tarball [False]
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
import rje, rje_db, rje_obj
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Initial (partially) functional version.
    # 0.2.0 - Updated functions with additional status measures and backups/archives division.
    # 0.2.1 - Fixed float division error.
    # 0.3.0 - Modified default homepath. Renamed backups to backupdb.
    # 0.3.1 - Fixed backups/backupdb bug.
    # 0.4.0 - Replaced module with uploadsh=X Full path for runnining upload.sh script ['/share/apps/unswdataarchive/2015-09-10/']
    # 0.4.1 - Added tryparent=T/F : Whether to try to run backup parent directory in case of failure [True]
    # 0.5.0 - Updated to run on Mac with OSX=T.
    # 0.6.0 - Added toggle to skip quiet created/updated/modified.
    # 0.7.0 - Added maxfiles cap and option to targz directories exceeding threshold.
    # 0.7.1 - Set checkarchive=F tryparent=F by default to make standard running quicker.
    # 0.7.2 - Added maxdirsize=INT  : Maximum directory size in bytes to generate backup [1e11 (~100Gb)]
    # 0.7.3 - Python 2.6 compatibility.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [ ] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add to SLiMSuite or SeqSuite.
    # [Y] : Add strict=T/F : Restrict processing to projects found in Projects file and add no new ones.
    # [?] : append=T/F      : Whether to append new data to archived file (True) or replace (False) [True]
    # [Y] : Fix the size measurement: du return kb not bytes.
    # [Y] : Not sure whether ftot counts should include directories? Might only consume directory if not there already.
    # [ ] : Add fullcheck=T/F : Special mode that reads all directories from backups.tdt and checks their status.
    # [ ] : Add projmode=X setting to recognise projects using different criteria?
    # [ ] : Update docstring to current running version.
    # [ ] : Add generation of the backup directory list when backing up whole project?
    # [ ] : Use tree for better documentation/checking of contents?
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('Archive', '0.7.3', 'April 2020', '2017')
    description = 'KDM Archive Manager'
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
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: Archive Class                                                                                           #
#########################################################################################################################
class Archive(rje_obj.RJE_Object):
    '''
    Archive Class. Author: Rich Edwards (2017).

    Str:str
    - Archived=FILE   : File to output archive summaries into (dir, project, rds, date, files, imports) ['archived.tdt']
    - BackupDB=FILE   : File to output backup summaries into  ['backups.tdt']
    - Dormant=FILE    : File to output dormant directories into  ['dormant.tdt']
    - HomeDir=PATH    : Home directory from which the archive script will be run ['/home/z3452659']
    - Projects=FILE   : Delimited file of Project and RDS code. If provided, will not use rds=X ['projects.tdt']
    - RDS=X           : UNSW_RDS project code to use (e.g. D0234444) []
    - UploadSH=X      : Full path for runnining upload.sh script ['/share/apps/unswdataarchive/2015-09-10/']

    Bool:boolean
    - ArchiveTGZ=T/F  : Whether to tar and zip then backup directories exceeding maxfiles=INT cutoff [False]
    - CheckArchive=T/F: Whether to run upload.sh on on directories that have been uploaded in last run [False]
    - CheckNum=T/F    : Whether to check numbers of files consumed by upload.sh versus directory contents [True]
    - Cleanup=T/F     : Whether to perform post-upload cleanup of backups and archived files [True]
    - RmDirs=T/F      : Delete archived directories. (Will ask if i>0) [False]
    - SkipDormant=T/F : Whether to skip uploads for dormant directories [True]
    - SkipQuiet=T/F   : Whether to skip uploads for quiet directories [True]
    - Strict=T/F      : Restrict processing to projects found in Projects file and add no new ones. [False]
    - TarGZ=T/F       : Whether to tar and zip directories to be deleted [True]
    - TryParent=T/F   : Whether to try to run backup parent directory in case of failure [False]
    - UseQSub=T/F     : Whether to use QSub for tarballing and then archiving the tarball [False]

    Int:integer
    - Dormancy=X      : Min number of days of inactivity before a directory gets rated as dormant [30]
    - MaxDirSize=INT  : Maximum directory size in bytes to generate backup (0 = no limit) [1e11 (~100Gb)]
    - MaxFiles=INT    : Maximum number of files in a directory to generate backups [10000]
    - Quiet=X     : Min number of days of inactivity before a directory will be deleted [1]

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - ArchiveDirs=LIST : List of directories to archive (should be project subdirectory full paths) []
    - BackupDirs=LIST : List of directories to backup (should be project subdirectory full paths) []
    - RmDirs=LIST     : List of directories to check archive status and then delete []

    Dict:dictionary

    Obj:RJE_Objects
    - DB:Database object =
         - Projects = [Project:str], RDS:str (loaded from Projects=FILE)

    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Archived','BackupDB','Dormant','HomeDir','Projects','RDS','UploadSH']
        self.boollist = ['CheckArchive','CheckNum','Cleanup','RmDirs','SkipDormant','SkipQuiet','Strict','TarGZ','TryParent','ArchiveTGZ','UseQSub']
        self.intlist = ['Dormancy','MaxDirSize','MaxFiles','Quiet']
        self.numlist = []
        self.filelist = []
        self.listlist = ['ArchiveDirs','BackupDirs']
        self.dictlist = []
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'HomeDir':os.path.expanduser('~'),'UploadSH':'/home/z3452659/unswdataarchive/upload.sh','Projects':'projects.tdt'})
        self.setBool({'CheckArchive':False,'CheckNum':True,'Cleanup':True,'RmDirs':False,'SkipDormant':True,'SkipQuiet':True,'Strict':False,'TarGZ':True,'TryParent':False,'ArchiveTGZ':False,'UseQSub':False})
        self.setInt({'Dormancy':30,'MaxDirSize':1e11,'MaxFiles':10000,'Quiet':1})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.baseFile('rds')
        self._setForkAttributes()   # Delete if no forking
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
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
                self._cmdReadList(cmd,'str',['RDS'])   # Normal strings
                self._cmdReadList(cmd,'path',['HomeDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['Archived','BackupDB','Dormant','Projects','UploadSH'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['CheckArchive','CheckNum','Cleanup','RmDirs','SkipDormant','SkipQuiet','Strict','TarGZ','TryParent','ArchiveTGZ','UseQSub'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['Dormancy','MaxDirSize','MaxFiles','Quiet'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['ArchiveDirs','BackupDirs'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        for skey in ['Archived','BackupDB','Dormant']:
            if skey == 'BackupDB': sfile = 'backups'
            else: sfile = skey.lower()
            if not self.getStrLC(skey): self.setStr({skey:'%s.%s.tdt' % (self.baseFile(),sfile)})
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): raise ValueError('Archive.setup() failed.')
            bdb = self.db('Backups')
            afields = bdb.fields()
            #!# Add a new method for selecting directories based on latest bdb.index('status')
            skipstatus = ['archived','dormant']
            if not self.getBool('CheckArchive'): skipstatus += ['modified','updated','created']
            ## ~ [1a] Setup up rmdirs list for removing redundancy from backupdirs ~~~~~~~~~~~~~~~~ ##
            rmdirs = []
            for rdir in self.list['ArchiveDirs']:
                rpath = rje.split(rdir)[0]
                if rpath not in rmdirs: rmdirs.append(rpath)
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Backup directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# backupdirs=LIST: List of directories to archive (should be project subdirectory full paths) []
            backups = []
            if self.list['BackupDirs']:
                self.printLog('#BACKUP','%s directories to backup' % rje.iLen(self.list['BackupDirs']))
                self.printLog('#BACKUP','Details will be output to %s' % self.getStr('BackupDB'))
                for adir in self.list['BackupDirs']:
                    apath = rje.split(adir)[0]
                    self.headLog(apath,line='=')
                    if not rje.exists(apath):
                        self.warnLog('%s missing: cannot backup' % apath)
                        if apath in bdb.index('dir') and bdb.indexEntries('dir',apath)[-1]['status'] != 'deleted':
                            project = self.getProject(apath)
                            rds = self.getRDS(project,add=not self.getBool('Strict'))
                            dirdata = {'dir':apath,'project':project,'rds':rds,'status':'deleted',
                                       'size':0,'fnum':0,'dnum':0,'ftot':0,'dtot':0,'files':0,'imports':0}
                            rje.delimitedFileOutput(self,self.getStr('BackupDB'),afields,datadict=dirdata,rje_backup=False)
                            bdb.addEntry(dirdata)
                        continue
                    if apath in backups:
                        self.printLog('#SKIP','Skipping duplicate: %s' % apath)
                        continue
                    #!#elif apath in rmdirs:
                    #!#    self.printLog('#SKIP','Skipping ArchiveDir: %s' % apath)
                    if (self.getBool('SkipDormant') or self.getBool('SkipQuiet')) and apath in bdb.index('dir'):
                        lastentry = bdb.indexEntries('dir',apath)[-1]
                        if lastentry['status'] in skipstatus:
                            postdata = self.dirDetails(apath)
                            postdata['date'] = rje.dateTime()
                            editdays = self.dayDif(self.lastEdit(apath,postdata),postdata['date'])
                            if self.getBool('SkipDormant') and editdays > self.getInt('Dormancy') > 0:
                                self.printLog('#DORMANT','Skipped %s: no change in %.1f (> %d) days' % (apath,editdays,self.getInt('Dormancy')))
                                if lastentry['status'] == 'archived':
                                    postdata['files'] = 0            # No. files consumed
                                    postdata['imports'] = 0          # No. files imported
                                    postdata['status'] = 'dormant'   # Whether directory looks the same post-import
                                    if not rje.exists(self.getStr('BackupDB')):
                                        rje.delimitedFileOutput(self,self.getStr('BackupDB'),afields,rje_backup=False)
                                    rje.delimitedFileOutput(self,self.getStr('BackupDB'),afields,datadict=postdata,rje_backup=False)
                                    bdb.addEntry(postdata)
                            elif self.getBool('SkipQuiet') and editdays > self.getInt('Quiet') > 0:
                                self.printLog('#QUIET','Skipped %s: no change in %.1f (> %d) days' % (apath,editdays,self.getInt('Quiet')))
                            elif self.archiveDirectory(apath,postdata)['status'] not in ['failed','denied','targz']:
                                backups.append(apath)
                            continue
                    if self.archiveDirectory(apath)['status'] not in ['failed','denied','targz']:
                        backups.append(apath)

            ## ~ [2b] Remove directories ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# archivedirs=LIST: List of directories to archive (should be project subdirectory full paths) []
            #quiet=X     : Min number of days of inactivity before a directory will be deleted [0]
            if rmdirs:
                self.warnLog('ArchiveDirs option not yet implemented!')
                #self.getBool('TarGZ')
                #self.getBool('RmDirs')


            ### ~ [3] Finish run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] Save updated projects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Cleanup'):
                if bdb.entryNum():
                    bdb.saveToFile(self.getStr('BackupDB'))

                self.warnLog('Cleanup option not yet implemented!')

                if self.dev():

                    #!# Sort out Backups, Dormant and Archived
                    if rje.exists(self.getStr('Archived')):
                        adb = self.db().addTable(self.getStr('Archived'),mainkeys=['dir','date'],name='Archived',expect=True,replace=True)
                        if adb.entryNum():

                            #!# Go through and check for reappearance of directory
                            #!# If exists, generate new stats and add to archive with status=restored

                            self.db('Archived').saveToFile(self.getStr('Archived'))



                    if rje.exists(self.getStr('Dormant')):
                        ddb = self.db().addTable(self.getStr('Dormant'),mainkeys=['dir','date'],name='Dormant',expect=True,replace=True)
                        ddb.dataFormat({'size':'int','fnum':'int','dnum':'int','ftot':'int','dtot':'int','files':'int','imports':'int'})
                        if ddb.entryNum():
                            self.db('Dormant').saveToFile(self.getStr('Dormant'))
            else: bdb.saveToFile(self.getStr('BackupDB'))
            ## ~ [3b] Save updated projects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = self.db('Projects')
            if self.getInt('ProjNum') < pdb.entryNum():
                self.printLog('#PROJ','%s Project-RDS links added' % rje.iStr(pdb.entryNum() - self.getInt('ProjNum')))
                if self.i() < 1 or rje.yesNo('Save projects table?'): pdb.saveToFile(self.getStr('Projects'))
            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] DB Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [1a] Projects Table setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = None
            self.setInt({'ProjNum':0})  # Number of projects in the projects table
            if rje.exists(self.getStr('Projects')):
                pdb = db.addTable(self.getStr('Projects'),mainkeys=['Project'],name='Projects',expect=True,replace=True)
                # Check file OK
                if not pdb or 'RDS' not in pdb.fields():
                    raise ValueError('Problem with Projects file "%s": needs Project and RDS fields.')
                self.setInt({'ProjNum':pdb.entryNum()})
            elif self.getBool('Strict') and self.getStrLC('Projects'): raise IOError('Projects file "%s" not found!')
            else: pdb = db.addEmptyTable('Projects',['Project','RDS'],['Project'],log=False)
            self.printLog('#PROJ','%d Project-RDS links loaded.' % self.getInt('ProjNum'))
            if not self.getInt('ProjNum') and not self.getStrLC('RDS'):
                raise ValueError('Need to have Project-RDS links or RDS given with rds=X.')
            ## ~ [1b] Backups table setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.exists(self.getStr('BackupDB')):
                bdb = db.addTable(self.getStr('BackupDB'),mainkeys=['dir','date'],name='Backups',expect=True,replace=True)
                bdb.dataFormat({'size':'int','fnum':'int','dnum':'int','ftot':'int','dtot':'int','files':'int','imports':'int'})
            else:
                afields = rje.split('dir project rds date files imports size fnum dnum ftot dtot status')
                db.addEmptyTable('Backups',afields,['dir','date'],log=True)
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### General Project Methods                                                                                 #
#########################################################################################################################
    def getRDS(self,project,add=True):   ### Returns the RDS to be used for a given project. If no RDS, returns None.
        '''
        Returns the RDS to be used for a given project. If no RDS, returns None. If add=True, add to Projects Table.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not project: return None
            pdb = self.db('Projects')
            ### ~ [2] Return set RDS if no Project File loaded ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if project not in pdb.data() and self.getStrLC('RDS'):
                if add or not pdb.entryNum():
                    if add: pdb.addEntry({'Project':project,'RDS':self.getStr('RDS')})
                    return self.getStr('RDS')
                else: self.warnLog('Project "%s" missing from Project table.' % project)
            ### ~ [3] Return RDS from Project File loaded ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return pdb.data(project)['RDS']
        except: self.errorLog('%s.method error' % self.prog())
#########################################################################################################################
    def getProject(self,path):  ### Parses the project from the path. Raises ValueError if fails.
        '''
        Parses the project from the path. Raises ValueError if fails.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not path.startswith('/'): raise ValueError('"%s" not a full path' % path)
            psplit = rje.split(path,os.sep)
            project = None
            i = 0
            while not project and i < len(psplit):
                if rje.matchExp('^(\S+-[A-Z][a-z][a-z]\d+)$',psplit[i]): project = psplit[i]
                i += 1
            if not project: raise ValueError('Cannot parse ProjectID from path: "%s"' % path)
            self.printLog('#PROJ','%s -> Project:%s' % (path,project))
            return project
        except ValueError: return None
        except: self.errorLog('%s.getProject error' % self.prog())
#########################################################################################################################
    def projectSuffix(self,path,project):  ### Parses the project suffix that must be appended to RDS from the path.
        '''
        Parses the project from the path. Raises ValueError if fails.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not path.startswith('/'): raise ValueError('"%s" not a full path' % path)
            if path.endswith('/'): path = path[:-1]
            psplit = rje.split(path,os.sep)
            if project not in psplit: raise ValueError('Cannot find ProjectID %s in path: "%s"' % (project,path))
            while psplit and psplit[0] != project: psplit.pop(0)
            if psplit: psplit = psplit[:-1]     # Drop the directory being added.
            return rje.join(psplit,os.sep)
        except: self.errorLog('%s.projectSuffix error' % self.prog()); raise
#########################################################################################################################
    ### <4> ### Directory information parsing methods                                                                   #
#########################################################################################################################
    def dirDetails(self,path,log=True):     ### Parses details of directory contents and returns as dictionary.
        '''Parses details of directory contents and returns as dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dpath = os.path.expanduser(path)
            ### ~ [2] Parse details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Total size of directory contents
            if self.getBool('OSX'): #i# Will be in 512 byte blocks
                size = int(rje.split(os.popen('du -d 1 %s' % dpath).readlines()[-1])[0])
            else:
                size = int(rje.split(os.popen('du %s --max-depth 1 --bytes' % dpath).readlines()[-1])[0])
            fnum = int(rje.split(os.popen('ls -1p %s | grep "/$" -cv' % dpath).readlines()[-1])[0])
            dnum = int(rje.split(os.popen('ls -1p %s | grep "/$" -c' % dpath).readlines()[-1])[0])
            project = self.getProject(dpath)
            rds = self.getRDS(project,add=not self.getBool('Strict'))
            #if self.getBool('CheckNum'): ftot = len(rje.getFileList(callobj=self,folder=dpath))
            if self.getBool('CheckNum') or self.getInt('MaxFiles') > 0:
                ftot = len(rje.listDir(callobj=self,folder=dpath,folders=False))
                dtot = len(rje.listDir(callobj=self,folder=dpath,files=False))
            else: ftot = dtot = 0
            ### ~ [3] Return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dirdata = {'dir':dpath,'project':project,'rds':rds,'size':size,'fnum':fnum,'dnum':dnum,'ftot':ftot,'dtot':dtot}
            if log:
                for dkey in rje.split('dir project rds'):
                    self.printLog('#%s' % dkey[:4].upper(),dirdata[dkey])
                for dkey in rje.split('size fnum dnum ftot dtot'):
                    self.printLog('#%s' % dkey[:4].upper(),rje.iStr(dirdata[dkey]))
            return {'dir':dpath,'project':project,'rds':rds,'size':size,'fnum':fnum,'dnum':dnum,'ftot':ftot,'dtot':dtot}
        except: self.errorLog('%s.dirDetails error' % self.prog()); raise
#########################################################################################################################
    ### <5> ### Data Archive Script Methods                                                                             #
#########################################################################################################################
    def archiveDirectory(self,directory,predata=None,parent=False):     ### Archives directory contents and returns details as dictionary.
        '''Archives directory contents and returns details as dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.headLog('Backup: %s' % directory,line='~')
            #self.printLog('#~~#','## ~~~~~ Backup: %s ~~~~~ ##' % directory)
            adb = self.db('Backups')    # dir project rds date files imports size fnum dnum status
            afields = adb.fields()
            if not predata: predata = self.dirDetails(directory)    # Directory details pre-import
            date = rje.dateTime()
            ## ~ [1a] Check size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getInt('MaxDirSize') > 0 and predata['size'] > self.getInt('MaxDirSize'):
                predata['date'] = date
                predata['files'] = 0            # No. files consumed
                predata['imports'] = 0          # No. files imported
                predata['status'] = 'denied'     # Whether directory looks the same post-import
                self.printLog('#STATUS','%s exceeds maxdirsize=%s threshold' % (predata['dir'],rje.iStr(self.getInt('MaxDirSize'))))
                self.warnLog('%s exceeds maxdirsize=%s threshold. Does this really need backup? If so, consider deleting/compressing some files first.' % (predata['dir'],rje.iStr(self.getInt('MaxDirSize'))))
                adb.addEntry(predata)
                return predata
            ## ~ [1b] Check filecount ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getInt('MaxFiles') > 0 and predata['ftot'] > self.getInt('MaxFiles'):
                predata['date'] = date
                predata['files'] = 0            # No. files consumed
                predata['imports'] = 0          # No. files imported
                predata['status'] = 'denied'     # Whether directory looks the same post-import
                self.printLog('#STATUS','%s exceeds maxfiles=%s threshold' % (predata['dir'],rje.iStr(self.getInt('MaxFiles'))))
                #!# Perform tgz and update directory then continue archiving
                if self.getBool('ArchiveTGZ') and not parent:
                    while directory.endswith('/'): directory = directory[:-1]
                    (parentdir, tardir) = os.path.split(directory)
                    tgzfile = '{0}.tgz'.format(directory)
                    if rje.exists(tgzfile) and not self.force() and rje.isYounger(tgzfile,directory) != directory:
                        self.printLog('\r#TARGZ','Tarball {0} found (force=F)'.format(tgzfile))
                        predata['status'] = 'targz'
                        adb.addEntry(predata)
                        return predata
                    elif rje.exists(tgzfile) and rje.isYounger(tgzfile,directory) == directory:
                        self.printLog('\r#TARGZ','Tarball {0} found but {1} younger (ignoredate=F)'.format(tgzfile,directory))
                    cwd = os.getcwd()
                    os.chdir(parentdir)
                    TARGZ = os.popen('tar -cvzf {0}.tgz {1}'.format(tardir,tardir))
                    line = TARGZ.readline()
                    filex = 0
                    while line:
                        filex += 1
                        if self.v() < 2:
                            self.progLog('\r#TARGZ','Tarballing %s files: %.f%%' % (rje.iStr(predata['ftot']),100.0*filex/predata['ftot']))
                        else: self.vPrint(rje.chomp(line),v=2)
                        line = TARGZ.readline()
                    TARGZ.close()
                    os.chdir(cwd)
                    self.printLog('\r#TARGZ','Tarballed %s files -> %s/%s.tgz' % (rje.iStr(predata['ftot']),parentdir,tardir))
                    predata['status'] = 'targz'
                    adb.addEntry(predata)
                    predata['dir'] = '%s/%s.tgz' % (parentdir,tardir)
                    predata['ftot'] = 1
                #predata['status'] = 'targz'     # Whether directory looks the same post-import
                else:
                    self.warnLog('%s exceeds maxfiles=%s threshold (archivetgz=F)' % (predata['dir'],rje.iStr(self.getInt('MaxFiles'))))
                    adb.addEntry(predata)
                    return predata
            filex = 0          # No. files consumed
            importx = -1       # No. files imported
            skipx = 0          # No. files skipped
            ignorex = 0        # No. files ignored
            ### ~ [2] Run Archive and parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cwd = os.getcwd()
            suffix = self.projectSuffix(predata['dir'],predata['project'])
            ucmd = '%s %s "/UNSW_RDS/%s/%s"' % (self.getStr('UploadSH'),predata['dir'],predata['rds'],suffix)
            self.printLog('#HOME','Running from %s' % self.getStr('HomeDir'))
            self.printLog('#UPLOAD',ucmd)
            if self.getBool('Test'):
                self.printLog('#TEST','No files consumed/imported')
            else:
                os.chdir(self.getStr('HomeDir'))
                if self.i() < 1 or self.yesNo('Run command?'):
                    #for line in os.popen(ucmd).readlines():
                    UPLOAD = os.popen(ucmd)
                    line = UPLOAD.readline()
                    while line:
                        if line.startswith('Consume:'): filex += 1
                        if line.startswith('live: imported'): importx = int(rje.split(line)[2])
                        elif rje.matchExp('^live:\s+imported\s+(\d+)',line): importx = int(rje.matchExp('^live:\s+imported\s+(\d+)',line)[0])
                        if rje.matchExp('^\s+skipped\s+(\d+)',line): skipx = int(rje.matchExp('^\s+skipped\s+(\d+)',line)[0])
                        if rje.matchExp('^\s+ignored\s+(\d+)',line): ignorex = int(rje.matchExp('^\s+ignored\s+(\d+)',line)[0])
                        if self.v() < 2:
                            if predata['ftot']:
                                self.progLog('\r#UPLOAD','Consuming %s files: %.f%%' % (rje.iStr(predata['ftot']),100.0*filex/predata['ftot']))
                        else: self.vPrint(rje.chomp(line),v=2)
                        line = UPLOAD.readline()
                    UPLOAD.close()
                    self.progLog('\r#UPLOAD','upload.sh consumed %s of %s files' % (rje.iStr(filex),rje.iStr(predata['ftot'])))
                os.chdir(cwd)
                if self.getBool('CheckNum'):
                    self.printLog('#FILES','%s of %s files consumed by upload.sh' % (rje.iStr(filex),rje.iStr(predata['ftot'])))
                    if filex < predata['ftot']:
                        self.warnLog('File count vs file consumption mismatch: %s' % predata['dir'])
                    elif filex > predata['ftot']:
                        self.warnLog('File count vs file consumption mismatch - possible directory creation: %s' % predata['dir'])
                else: self.printLog('#FILES','%s files consumed by upload.sh' % rje.iStr(filex))
                self.printLog('#IMPORT','%s files imported to %s' % (rje.iStr(importx),predata['rds']))
                if skipx: self.printLog('#SKIP','%s files skipped.' % (rje.iStr(skipx)))
                if ignorex: self.printLog('#IGNORE','%s files ignored.' % (rje.iStr(ignorex)))
            ### ~ [3] Check for changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            postdata = self.dirDetails(directory)    # Directory details pre-import
            static = self.dirStatic(predata,postdata)
            if importx < 0:
                status = 'failed'
                self.printLog('#STATUS','%s backup to %s failed!' % (predata['dir'],predata['rds']))
                if self.getBool('TryParent'):  # Option to try parent directory
                    project = predata['project']
                    psplit = rje.split(directory, os.sep)
                    psplit = psplit[:-1]
                    if project in psplit:   # Can go up a level
                        directory = rje.join(psplit, os.sep)
                        self.printLog('#PARENT','Trying to backup parent directory: %s' % directory)
                        if self.archiveDirectory(directory,parent=True)['status'] != 'failed':
                            status = 'parent'
                    else: self.printLog('#STATUS','Total backup failure for %s' % project)
            elif importx > 0:
                if filex == importx:
                    status = 'created'
                    self.printLog('#STATUS','%s created in %s' % (predata['dir'],predata['rds']))
                else:
                    status = 'updated'
                    self.printLog('#STATUS','%s updated in %s' % (predata['dir'],predata['rds']))
            elif static:
                #!# Add dormancy assessment
                status = 'archived'
                self.printLog('#STATUS','%s archived in %s' % (predata['dir'],predata['rds']))
            else:
                status = 'modified'
                self.printLog('#STATUS','%s modified since update to %s' % (predata['dir'],predata['rds']))
            ### ~ [4] Add information to archive table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            predata['date'] = date
            predata['files'] = filex         # No. files consumed
            predata['imports'] = importx     # No. files imported
            predata['status'] = status       # Whether directory looks the same post-import
            if not rje.exists(self.getStr('Backups')):
                rje.delimitedFileOutput(self,self.getStr('BackupDB'),afields,rje_backup=False)
            rje.delimitedFileOutput(self,self.getStr('BackupDB'),afields,datadict=predata,rje_backup=False)
            adb.addEntry(predata)
            return predata
        except: self.errorLog('%s.archiveDirectory error' % self.prog()); raise
#########################################################################################################################
    def lastEdit(self,directory,postdata={}):   ### Returns last date and time of edit.
        '''
        Returns last date and time of edit.
        >> directory:str = directory to assess using the contents of backups.tdt.
        >> postdata = data from self.dirDetails(directory) to compare to entry
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bdb = self.db('Backups')
            blist = bdb.indexEntries('dir',directory)[0:]
            if not postdata: postdata = self.dirDetails(directory)
            edate = postdata['date']
            while blist:
                predata = blist.pop(-1)
                if self.dirStatic(predata,postdata): edate = predata['date']
                else: break     # Changed, even if it was earlier the same
                #x# if entry['status'] in ['created','updated','modified']
            return edate
        except: self.errorLog('%s.lastEdit error' % self.prog()); raise
#########################################################################################################################
    def dayDif(self,date1,date2):  ### Calculates and returns the difference between date1 and date2 in days
        '''Calculates and returns the difference between date1 and date2 in days.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            time1 = rje.matchExp('(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+)',date1)
            time2 = rje.matchExp('(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+)',date2)
            ### ~ [2] Calculate difference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            diff = []
            ## ~ [2a] Start in regular units ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in range(6):
                diff.append(int(time2[i])-int(time1[i]))
            ## ~ [2b] Convert to days ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            days = diff[0:]
            months = [31,28,31,30,31,30,31,31,30,31,30,31]
            days[0] *= 365.0
            days[1] = sum(months[:int(time2[1])]) - sum(months[:int(time1[1])])
            days[3] /= 24.0
            days[4] /= (24.0*60)
            days[5] /= (24.0*60*60)
            ### ~ [3] Return sum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return sum(days)
        except: self.errorLog('%s.dayDif error' % self.prog()); raise
#########################################################################################################################
    def isDormant(self,entry,postdata={}):  ### Returns True/False whether directory is dormant
        '''
        Returns True/False whether directory is dormant.
        >> entry = backups.tdt entry (or data from self.dirDetails(directory))
        >> postdata = data from self.dirDetails(directory) to compare to entry
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not postdata: postdata = self.dirDetails(entry['dir'])
            if self.dirStatic(entry,postdata) and self.getInt('Dormancy') > 0:
                return self.dayDif(self.lastEdit(entry['dir'],postdata),postdata['date']) > self.getInt('Dormancy')
            return False
        except: self.errorLog('%s.isDormant error' % self.prog()); raise
#########################################################################################################################
    def isQuiet(self,entry,postdata={}):  ### Returns True/False whether directory is dormant
        '''
        Returns True/False whether directory is dormant.
        >> entry = backups.tdt entry (or data from self.dirDetails(directory))
        >> postdata = data from self.dirDetails(directory) to compare to entry
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not postdata: postdata = self.dirDetails(entry['dir'])
            if self.dirStatic(entry,postdata) and self.getInt('Quiet') > -1:
                return self.dayDif(self.lastEdit(entry['dir'],postdata),postdata['date']) > self.getInt('Quiet')
            return False
        except: self.errorLog('%s.isQuiet error' % self.prog()); raise
#########################################################################################################################
    def dirStatic(self,predata,postdata):  ### Returns True/False whether directory data has changed
        '''
        Returns True/False whether directory is dormant.
        >> entry = backups.tdt entry (or data from self.dirDetails(directory))
        >> postdata = data from self.dirDetails(directory) to compare to entry
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            static = True
            for d in ['size','fnum','dnum','ftot','dtot']: static = static and predata[d] == postdata[d]
            return static
        except: self.errorLog('%s.dirStatic error' % self.prog()); raise
#########################################################################################################################
### End of SECTION II: Archive Class                                                                                    #
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
    try: Archive(mainlog,cmd_list).run()

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
