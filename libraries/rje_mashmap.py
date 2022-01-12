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
Module:       rje_mashmap
Description:  MashMap control and parsing module
Version:      0.0.0
Last Edit:    19/11/21
Copyright (C) 2021  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module will run MashMap on two assemblies and have the code for additional parsing and filtering methods where
    appropriate.

Commandline:
    ### ~ Main MashMap run options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Genome assembly used as the query (also recognises assembly=FILE) []
    reference=FILE  : Genome assembly used as the searchdb (also refgenome=FILE or searchdb=FILE) []
    basefile=X      : Prefix for output files [$SEQIN.$REFERENCE]
    threads=INT     : Number of threads to use for MashMap runs []
    seglength=INT   : MashMap segment length (--seqLength <value>) [1000]
    mashpid=NUM     : MashMap threshold for identity [85]
    mashfilt=X      : MashMap filtering mode (map/both/none; both=one-to-one) [both]
    ### ~ MashMap filtering options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    selfhit=T/F     : Whether to include self hits (matching name, pos, strand) [True]
    minloclen=INT   : Minimum length for aligned chunk to be kept (shortest hit length in bp) [1000]
    minlocid=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]
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
import rje, rje_obj
import rje_db
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
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
    # [ ] : Add rapid R-based filtering?
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_MASHMAP', '0.0.0', 'November 2021', '2021')
    description = 'MashMap control and parsing module'
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
### SECTION II: MashMap Class                                                                                           #
#########################################################################################################################
class MashMap(rje_obj.RJE_Object):
    '''
    MashMap Class. Author: Rich Edwards (2021).

    Str:str
    - MashFilt=X      : MashMap filtering mode (map/both/none; both=one-to-one) [both]
    - Reference=FILE  : Genome assembly used as the searchdb (also refgenome=FILE or searchdb=FILE) []
    - SeqIn=FILE      : Genome assembly used as the query (also recognises assembly=FILE) []

    Bool:boolean
    - SelfHit=T/F     : Whether to include self hits (matching name, pos, strand) [True]

    Int:integer
    - MinLocLen=INT   : Minimum length for aligned chunk to be kept (shortest hit length in bp) [1000]
    - SegLength=INT   : MashMap segment length (--seqLength <value>) [1000]

    Num:float
    - MashPID=PERC     : MashMap threshold for identity [85]
    - MinLocID=PERC   : Minimum percentage identity for aligned chunk to be kept (local %identity) [50]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = RJE_DB Database object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['MashFilt','Reference','SeqIn']
        self.boollist = ['SelfHit']
        self.intlist = ['MinLocLen','SegLength']
        self.numlist = ['MashPID','MinLocID']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = ['DB']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'MashFilt':'one-to-one'})
        self.setBool({'SelfHit':False})
        self.setInt({'MinLocLen':1000,'SegLength':1000})
        self.setNum({'MashPID':80,'MinLocID':50})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                self._cmdReadList(cmd,'str',['MashFilt'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Reference','SeqIn'])  # String representing file path
                self._cmdRead(cmd, type='file', att='Reference', arg='refgenome')
                self._cmdRead(cmd, type='file', att='Reference', arg='searchdb')
                self._cmdRead(cmd, type='file', att='SeqIn', arg='assembly')
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['SelfHit'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['MinLocLen','SegLength'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                self._cmdReadList(cmd,'perc',['MashPID','MinLocID'])  # Percentage, converts to 1-100 scale.
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'list',['Att'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                #self._cmdReadList(cmd,'glist',['Att']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''
        Main run method.

        MashMap output is space-delimited with each line consisting of:
         query name, length, 0-based start, end, strand, target name, length, start, end and mapping nucleotide identity.

        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.runMashMap()
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.checkForProgram('mashmap',report=True,warn=True,needed=True):
                return False
            if not rje.checkForFiles([self.getStr('SeqIn'),self.getStr('Reference')],basename='',log=self.log,cutshort=True,ioerror=True,missingtext='Not found.'):
                return False
            if not self.getStrLC('Basefile'):
                self.baseFile('{0}.{1}'.format(rje.baseFile(self.getStr('SeqIn'),strip_path=True),rje.baseFile(self.getStr('Reference'),strip_path=True)))
            if self.getStrLC('MashFilt') in ['both','one-to-one']:
                self.setStr({'MashFilt':'one-to-one'})
            elif self.getStrLC('MashFilt') in ['map','query']:
                self.setStr({'MashFilt':'map'})
            else: self.setStr({'MashFilt':'none'})
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def runMashMap(self):      ### Call mashmap and return output file, or raise error if failed.
        '''
        Call mashmap and return output file, or raise error if failed.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.infoLog('MashMap filter mode: {0}'.format(self.getStr('MashFilt')))
            mashout = '{0}.s{1}.pi{2}.out'.format(self.baseFile(),self.getInt('SeqLength'),self.getNum('MashPID'))
            if rje.checkForFiles([mashout],log=self.log):
                if self.force():
                    self.printLog('#MASH','MashMap output found but force=T: regenerating')
                    rje.backup(self,mashout)
                else:
                    return mashout
            mashcmd = 'mashmap -q {0} -r {1} -f {2} -t {3} -s {4} --pi {5} -o {6}'.format(self.getStr('SeqIn'),self.getStr('Reference'),self.getStr('MashFilt'),self.threads(),self.getInt('SegLength'),self.getNum('MashPID'),mashout)
            self.loggedSysCall(mashcmd)
            rje.checkForFiles([mashout],basename='',log=self.log,ioerror=True)
            return mashout
        except: self.errorLog('%s.runMashMap error' % self.prog())
#########################################################################################################################
    def parseOut(self,filename,table='mashmap'):    ### Parses output file into db object
        '''
        Parses output file into db object.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mhead = ['Qry','QryLen','QryStart','QryEnd','Strand','Hit','HitLen','HitStart','HitEnd','PercID']
            mkeys = ['Qry','QryStart','QryEnd','Strand','Hit','HitStart','HitEnd']
            db = self.db()
            ### ~ [2] Load table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            delimit = db.info['Delimit']
            mashdb = db.addTable(filename,mainkeys=mkeys,datakeys='All',delimit=' ',headers=mhead,ignore=['#'],name=table,expect=True,replace=True)
            mashdb.info['Delimit'] = delimit
            mashdb.dataFormat({'QryLen':'int','QryStart':'int','QryEnd':'int','HitLen':'int','HitStart':'int','HitEnd':'int','PercID':'num'})
            ### ~ [3] Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mashdb.addField('AlnNum',evalue=0); ix = 1
            mashdb.addField('Length',evalue=0)
            mashdb.addField('Identity',evalue=0)
            ex = 0.0; etot = mashdb.entryNum()
            for entry in mashdb.entries():
                self.progLog('Formatting mashmap output: {0:.2f}%'.format(ex/etot)); ex += 100.0
                for f in ['QryStart','QryEnd','HitStart','HitEnd']: entry[f] += 1
                entry['AlnNum'] = ix; ix += 1
                entry['Length'] = min((entry['HitEnd']-entry['HitStart']),(entry['QryEnd']-entry['QryStart'])) + 1
                entry['Identity'] = int(0.5 + (entry['Length'] * entry['PercID'] / 100.0))
            self.printLog('Formatting mashmap output complete.')
            ### ~ [4] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mashdb.dropEntries(['Length<%d' % self.getInt('MinLocLen')],logtxt='MinLocLen',log=True)
            mashdb.dropEntries(['PercID<%f' % self.getNum('MinLocID')],logtxt='MinLocID',log=True)
            return mashdb
        except: self.errorLog('%s.parseOut error' % self.prog())
#########################################################################################################################
### End of SECTION II: MashMap Class                                                                                        #
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
    try:#NewClass(mainlog,cmd_list).run()
        rje.printf('\n\n{0}\n\n *** No standalone functionality! *** \n\n'.format(rje_obj.zen()))

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
