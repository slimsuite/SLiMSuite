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
Module:       rje_kat
Description:  KAT wrapper and parser
Version:      0.1.0
Last Edit:    27/09/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module has no standalone functions. The KAT object is designed to be inherited by other RJE program objects.

Commandline:
    ### ~ Main KAT wrapper run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE          : Input sequence assembly [None]
    basefile=FILE       : Root of output file names [$SEQIN basefile]
    kmerreads=FILELIST  : File of high quality reads for KAT kmer analysis []
    10xtrim=T/F         : Whether to trim 16bp 10x barcodes from Read 1 of Kmer Reads data for KAT analysis [False]
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
    # 0.1.0 - Fixed _setKatAttributes bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Populate Module Docstring with basic info.
    # [ ] : Populate makeInfo() method with basic info.
    # [ ] : Add full description of program to module docstring.
    # [ ] : Create initial working version of program.
    # [ ] : Add REST outputs to restSetup() and restOutputOrder()
    # [ ] : Add to SLiMSuite or SeqSuite.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_KAT', '0.0.0', 'September 2021', '2021')
    description = 'KAT wrapper and parser'
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
### SECTION II: KAT Class                                                                                               #
#########################################################################################################################
class KAT(rje_obj.RJE_Object):
    '''
    KAT Class. Author: Rich Edwards (2021).

    This class is designed for being inherited. In this case:
    - add self._setKatAttributes() to the _setAttributes() method
    - add self._katCmd(cmd) to parse kat options

    Str:str
    - SeqIn=FILE      : Input sequence assembly []

    Bool:boolean
    - 10xTrim=T/F     : Whether to trim 16bp 10x barcodes from Read 1 of Kmer reads data [False]

    Int:integer

    Num:float

    File:file handles with matching str filenames
    
    List:list
    - KmerReads=FILELIST   : File of reads for KAT kmer analysis []

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
        self.strlist = ['SeqIn']
        self.boollist = ['10xTrim']
        self.intlist = []
        self.numlist = []
        self.filelist = []
        self.listlist = ['KmerReads']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'10xTrim':False})
        self.setInt({})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        self._setKatAttributes()
#########################################################################################################################
    def _setKatAttributes(self):       ### Sets forking attributes for use in all classes
        '''Sets general forking attributes for use in all classes.'''
        self.str['SeqIn'] = 'None'
        self.list['KmerReads'] = []
        self.bool['10xTrim'] = False
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
                self._katCmd(cmd)   # Parse kat options
                ### Class Options (No need for arg if arg = att.lower()) ###
                #self._cmdRead(cmd,type='str',att='Att',arg='Cmd')  # No need for arg if arg = att.lower()
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                #self._cmdReadList(cmd,'file',['SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['10xTrim'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['KmerReads']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    def _katCmd(self,cmd):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self._cmdReadList(cmd,'file',['SeqIn'])  # String representing file path
        self._cmdReadList(cmd,'bool',['10xTrim'])  # True/False Booleans
        self._cmdReadList(cmd,'glist',['KmerReads']) # List of files using wildcards and glob
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.katKmers()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if not self.baseFile(return_none=''):
                if self.getStrLC('SeqIn'): self.baseFile(rje.baseFile(self.getStr('SeqIn'),strip_path=True))
            self.printLog('#BASE','Output file basename: %s' % self.baseFile())
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### KAT run Methods                                                                                         #
#########################################################################################################################
    def checkForKat(self,report=True):
        return self.checkForProgram('kat',report,needed=False)
#########################################################################################################################
    def katKmers(self,assembly=None,kmerfiles=None,basefile=None,force=None,trim10x=True):    ### Performs read kmer kat sect analysis
        '''
        Performs read kmer kat sect analysis. Generates:
        - '{0}.kat-stats.tsv'.format(basefile) = kmer summary per sequence
        - '{1}.kat-counts.cvg'.format(basefile) = kmer counts per position (CVG format)
        >> assembly:str [None] = Assembly file. Will use self.getStr('SeqIn') if None
        >> kmerfiles:list [None] = files for setting kmers to count (self.list['KmerReads'] if None)
        >> basefile:bool [None] = output file prefix (self.baseFile() if None)
        >> force:bool [None] = whether to overwrite existing files (self.force() if None)
        >> trim10x:bool [True] = Whether to check 10xtrim setting.
        << katfile or None if failed
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.checkForKat(report=True): return None
            if not assembly: assembly = self.getStr('SeqIn')
            seqin = assembly
            if kmerfiles:
                if type(kmerfiles) == type('str'):
                    kmerfiles = [kmerfiles]
            else:
                if not self.list['KmerReads']:
                    self.printLog('#KAT','Cannot use KAT kmer analysis without KmerReads data')
                    return None
                kmerfiles = self.list['KmerReads']
            rje.checkForFiles(filelist=[seqin]+kmerfiles,basename='',log=self.log,cutshort=False,ioerror=True,missingtext='Not found: aborting KAT run.')
            if not basefile: basefile = self.baseFile(return_none=None)
            if not basefile: self.baseFile(rje.baseFile(assembly,strip_path=True))
            if force == None: force = self.force()
            katfile = '{}.kat-stats.tsv'.format(basefile)
            # seq_name        median  mean    gc%     seq_length      kmers_in_seq    invalid_kmers   %_invalid       non_zero_kmers  %_non_zero      %_non_zero_corrected
            katcvg =  '{}.kat-counts.cvg'.format(basefile)
            #i# Check for files
            if not force and rje.checkForFiles(filelist=[katfile,katcvg],basename='',log=self.log,cutshort=False,ioerror=False,missingtext='Not found: will generate.'):
                return katfile
            self.backup(katfile,appendable=False)
            self.backup(katcvg,appendable=False)
            ### ~ [2] Run KAT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            katcall = 'kat sect -t {} -o {}.kat {} {}'.format(self.threads(),basefile,seqin,' '.join(kmerfiles))
            if trim10x and self.getBool('10xTrim'):
                trim5 = ['16'] + ['0'] * (len(self.list['KmerReads']) - 1)
                trim5 = ','.join(trim5)
                katcall = 'kat sect -t {} --5ptrim {} -o {}.kat {} {}'.format(self.threads(),trim5,basefile,seqin,' '.join(kmerfiles))
            self.printLog('#SYS',katcall)
            #i# Catching completion in case KAT hangs after running
            KAT = os.popen(katcall)
            while not KAT.readline().startswith('Total runtime'): continue
            KAT.close()
            ### ~ [3] Check for outputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFiles(filelist=[katfile,katcvg],basename='',log=self.log,cutshort=False,ioerror=True,missingtext='Not found: KAT failed?'):
                return katfile
        except: self.errorLog('%s.katKmers error' % self.prog()); return None
#########################################################################################################################
    def selfKat(self,assembly=None): ### Performs self KAT search
        '''
        Performs self KAT search.
        >> assembly:str [None] = Assembly file. Will use self.getStr('SeqIn') if None
        '''
        try:### ~ [1] Assembly versus self ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not assembly: assembly = self.getStr('SeqIn')
            basefile = self.baseFile(return_none=None)
            if not basefile: self.baseFile(rje.baseFile(assembly, strip_path=True))
            basefile = '{0}.self'.format(self.basefile())
            return self.katKmers(assembly=assembly,basefile=basefile,kmerfiles=[assembly],trim10x=False)
        except: self.errorLog('%s.selfKat error' % self.prog())
#########################################################################################################################
    ### <4> ### KAT Parse Methods                                                                                       #
#########################################################################################################################
    def parseKat(self,katfile=None,name='kat',diploidocus=False):   ### Parse main KAT file into 'kat' database table and returns.
        '''
        Parse main KAT file into 'kat' database table and returns.
        >> katfile:str [$BASEFILE.kat-stats.tsv] = KAT summary file to read into table.
        >> name:str ['kat'] = Name for database table
        >> diploidocus:bool [False] = Whether to rename fields in line with Diploidocus
        << katdb
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            if not katfile: katfile = '{}.kat-stats.tsv'.format(self.baseFile())
            ### ~ [2] Parse ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # seq_name median mean gc% seq_length kmers_in_seq invalid_kmers %_invalid non_zero_kmers %_non_zero %_non_zero_corrected
            katdb = db.addTable(katfile, mainkeys=['seq_name'], expect=True, name=name)
            ## ~ [2a] Diploidocus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # SeqName MedK AvgK SeqGC KPerc
            if katdb and diploidocus:
                katdb.renameField('seq_name','SeqName')
                katdb.renameField('median','MedK')
                katdb.renameField('mean','AvgK')
                katdb.renameField('gc%','SeqGC')
                katdb.renameField('%_non_zero_corrected','KPerc')
                katdb.setFields(['SeqName','MedK','AvgK','SeqGC','KPerc'])
                for entry in katdb.entries(): entry['SeqName'] = entry['SeqName'].split()[0]
                katdb.remakeKeys()
            return katdb
        except: self.errorLog('%s.parseKat error' % self.prog()); return None
#########################################################################################################################
### End of SECTION II: KAT Class                                                                                        #
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
    try: KAT(mainlog,cmd_list).run()

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
