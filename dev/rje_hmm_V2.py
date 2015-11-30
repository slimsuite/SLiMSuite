#!/usr/bin/python

# See below for name and description
# Copyright (C) 2012 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_hmm
Description:  HMMer Control Module
Version:      2.0
Last Edit:    23/01/13
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to perform basic HMM functions using the HMMer program. Currently, there are three functions
    that may be performed, separately or consecutively:
    * 1. Use hmmbuild to construct HMMs from input sequence files
    * 2. Search a sequence database with HMMs files
    * 3. Convert HMMer output into a delimited text file of results.

    Version 2.x onwards are redeveloped to work with HMMER3.0. Please used rje_hmm_V1.py for legacy runs.

Commandline:
    ## Build Options ##
    makehmm=LIST        : Sequence file(s). Can include wildcards [None]
    hmmcalibrate=T/F    : Whether to calibrate HMM files once made [True]

    ## Search Options ##    
    hmm=LIST        : HMM file(s). Can include wildcards. [*.hmm]
    searchdb=FILE   : Fasta file to search with HMMs [None]
    hmmoptions=LIST : List or file of additional HMMer search options (joined by whitespace) []
    hmmpfam=T/F     : Performs standard HMMer PFam search (--cut_ga) (or processes if present) [False]
    hmmout=FILE     : Pipe results of HMM searches into FILE [None]
    hmmres=LIST     : List of HMM search results files to convert (wildcards allowed) []
    hmmtab=FILE     : Delimited table of results ('None' to skip) [searchdb.tdt]
    cleanres=T/F    : Option to reduce size of HMM results file by removing no-hit sequences [True]

    ## System Parameters ##
    hmmerpath=PATH  : Path for hmmer files [/home/richard/Bioware/hmmer-2.3.2/src/] 
    force=T/F       : Whether to force regeneration of new HMMer results if already existing [False]
    gzip=T/F        : Whether to gzip (and gunzip) HMMer results files (not Windows) [True]
    
Classes:
    HMMRun Object = Full HMM run
    rje_db.Database Object = Stores detailed hit information as with the rje_blast V2.x module.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_zen
import rje_blast_V1 as rje_blast #?#
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 2.0 - Initial HMMER3.0 version based on Version 1.2 and RJE_BLAST 2.0.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Create initial working version of program based on V1.2 and rje_blast 2.0.
    # [ ] : Add stats and more HMM options
    # [ ] : Tidy up classes in line with other RJE modules
    # [ ] : Add a cleanres=T option to reduce size of HMM results file by removing no-hit sequences.
    # [ ] : Update to work with HMMER3.0
    # - [ ] : hmmpress to format Pfam HMM downloads into binary database
    # - [ ] : hmmscan to replace hmmpfam --cut_ga
    # - [ ] : phammer to replace BLASTP searches
    # - [ ] : jackhammer to replace PSI-BLAST searches
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_HMM', '2.0', 'January 2013', '2005')
    description = 'RJE HMM Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is under development and may contain bugs!',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
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
### SECTION II: HMMRun Class                                                                                            #
#########################################################################################################################
class HMMRun(rje_obj.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2012).

    Str:str
    - SearchDB = Fasta file to search with HMMs [None]
    - HMMTab = Delimited table of results ('None' to skip) [searchdb.hmmer.tdt]
    - HMMOut = Pipe results of HMM searches into FILE [None]
    - HMMerPath = path for hmmer files [c:/bioware/hmmer/] *Use fwd slashes
    
    Bool:boolean
    - CleanRes = Option to reduce size of HMM results file by removing no-hit sequences [True]
    - GZip = Whether to gzip (and gunzip) HMMer results files (not Windows) [True]
    - HMMCalibrate = Whether to calibrate HMM files once made [True]
    - HMMPFam = Performs standard HMMer PFam search (--cut_ga) (or processes if present) [False]

    Int:integer

    Num:float
    
    List:list
    - MakeHMM = Sequence file(s). Can include wildcards [None]
    - HMM = HMM file(s). Can include wildcards. [*.hmm]
    - HMMOptions = List or file of additional HMMer search options (joined by whitespace) []

    Dict:dictionary    

    Obj:RJE_Objects
    - DB = rje_db.Database object for storing settings and results
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['SearchDB','HMMOut','HMMTab','HMMerPath']
        self.boollist = ['HMMCalibrate','HMMPFam','GZip','CleanRes']
        self.intlist = []
        self.numlist = []
        self.listlist = ['MakeHMM','HMM','HMMOptions']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'HMMerPath':rje.makePath('/home/richard/Bioware/hmmer-2.3.2/src/'),'HMMOut':'','HMMTab':''})
        self.setBool({'HMMCalibrate':True,'GZip':True,'CleanRes':True})
        self.setInt({})
        self.setNum({})
        self._cmdRead(cmd='hmm=*.hmm',type='glist',att='HMM')
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        ### Other Attributes ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
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
                self._cmdReadList(cmd,'path',['HMMerPath'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['SearchDB','HMMTab','HMMOut'])  # String representing file path 
                self._cmdReadList(cmd,'bool',['HMMCalibrate','HMMPFam','GZip','CleanRes'])  # True/False Booleans
                #self._cmdReadList(cmd,'int',['Att'])   # Integers
                #self._cmdReadList(cmd,'float',['Att']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['HMMOptions'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['HMM','MakeHMM','HMMRes']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def _run(self):     ### Controls main Class functions                                                           #V2.0
        '''
        Controls main Class functions:
        * 1. Use hmmbuild to construct HMMs from input sequence files
        * 2. Search a sequence database with HMMs files
        * 3. Convert HMMer output into a delimited text file of results.
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.setup(): return False
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Build ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for seqfile in self.list['MakeHMM']:
                hmmfile = self.buildHMM(seqfile)
                if hmmfile: self.list['HMM'].append(hmmfile)
            ## ~ [2b] Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.deBug(self.list['HMM'])
            if self.list['HMM'] and os.path.exists(self.str['SearchDB']) and self.str['HMMOut'].lower() not in ['','none']:
                rje.backup(self,self.str['HMMOut'],unlink=True)
            for hmm in self.list['HMM']: self.hmmSearch(hmm,outfile=self.str['HMMOut'])
            ## ~ [2c] Read in results and tabulate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if rje.checkForFile(self.getStr('HMMOut')): self.readHMM(self.getStr('HMMOut'))
            if self.getStr('HMMTab').lower() not in ['','none']: self.db('Hit').saveToFile(self.getStr('HMMTab'),append=self.getBool('Append'))
            return True                   
        except:
            self.errorLog('Fatal Error during rje_hmm._run()',quitchoice=True)
            return False
#########################################################################################################################
    def setup(self):    ### Main class setup method.                                                                #V2.0
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            ## ~ [1a] HMM Runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db.addEmptyTable('Run',['Name','Type','HMMFile','SeqFile','HMMCmd','HMMNum','SeqNum'],['Name'])
            ## ~ [1b] HMM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #db.addEmptyTable('Search',['Run','HMM','Length','Hits'],['Run','HMM'])  #?#
            ## ~ [1c] HMM Hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db.addEmptyTable('Hit',['Run','HMM','Hit','Description','Score','Eval','Start','End','Hit_Start','Hit_End'],['Run','HMM','Hit'])
            ## ~ [1d] Local Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #db.addEmptyTable('Local',['Run','Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],['Run','Query','Hit','AlnID'])
            ## ~ [1e] GABLAM Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #db.addEmptyTable('GABLAM',['Run','Query','Hit','QryHit','GABLAM Start','GABLAM End','GABLAM Dirn','GABLAM Len','GABLAM Sim','GABLAM ID','GABLAM Frag','GABLAMO Start','GABLAMO End','GABLAMO Dirn','GABLAMO Len','GABLAMO Sim','GABLAMO ID','GABLAMO Frag'],['Run','Query','Hit','QryHit'])
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### General Class Data Methods                                                                              #
#########################################################################################################################
    def clear(self): self.setup()   # Currently, this will generate a new Database object full of empty tables.     #V2.0
#########################################################################################################################
    def addDBRun(self):     # Adds DBRun using current option settings                                              #V2.0
        '''Adds DBRun using current option settings.'''
        dbrun = {}
        for add_dict in [self.str,self.bool,self.int,self.num]: rje.combineDict(dbrun,add_dict)
        self.db('Run').addEntry(dbrun)
#########################################################################################################################
#    def searchNum(self,run=None):   # Returns number of Searches for specific run or all if no run given            #V2.0
#        '''Returns number of Searches for specific run or all if no run given.'''
#        if run: return len(self.db('Search').indexEntries('Run',run))
#        else: return self.db('Search').entryNum()
#########################################################################################################################
    def hitNum(self,run=None):  # Returns number of Hits for specific run or all if no run given                    #V2.0
        '''Returns number of Hits for specific run or all if no run given.'''
        if run: return len(self.db('Search').indexEntries('Run',run))
        else: return self.db('Hit').entryNum()
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def _method(self):      ### Generic method
        '''
        Generic method. Add description here (and arguments.)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return
        except: self.errorLog('%s.method error' % self)
#########################################################################################################################


#########################################################################################################################
    ### <2> ### OLD Methods                                                                                        #
#########################################################################################################################
    ### <3> ### HMM Build Methods                                                                                       #
#########################################################################################################################
    def buildHMM(self,seqfile,hmmfile=None):    ### Makes an HMM from a sequence alignment file
        '''
        Makes an HMM from a sequence alignment file.
        >> seqfile:str = Name of sequence file
        >> hmmfile:str = Name of HMM file [*.hmm]
        << hmmfile if made, None if failed.
        '''
        try:
            ### Setup ###
            _hmmpath = self.info['HMMerPath']
            if not hmmfile:
                hmmfile = '%s.hmm' % rje.baseFile(seqfile)

            ### Build HMM ##
            os.system('%shmmbuild %s %s' % (_hmmpath,hmmfile,seqfile))
            if self.opt['HMMCalibrate']:
                os.system('%shmmcalibrate %s' % (_hmmpath,hmmfile))
            return hmmfile      #!# Add error catching during build/calibrate (How?!) #!#
        except:
            self.log.errorLog('Oh my, what a calamity during buildHMM(%s)!' % seqfile)
            return None
#########################################################################################################################
    ### <4> ### HMM Search Methods                                                                                      #
#########################################################################################################################
    def hmmSearch(self,hmm,dbase=None,outfile=None,wait=True):    ### Performs HMMer Search using object attributes
        '''
        Performs HMMer Search using object attributes.
        >> hmm:str = Name of HMM file 
        >> dbase:str = Name of DBase file [self.info['SearchDB']]
        >> outfile:str = Name of Output file file [self.info['HMMOut']]
        >> wait:boolean  = whether to wait for HMMer. [True]
        << returns outfile or None if fails
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Input files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not rje.checkForFile(hmm): self.printLog('#ERR','HMM file %s is missing!' % hmm); return None
            if not dbase: dbase = self.info['SearchDB']
            if not rje.checkForFile(dbase): self.printLog('#ERR','Database file "%s" is missing!' % dbase); return None
            ## ~ [1b] ~ Output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not outfile or outfile.lower() in ['','none']:       # Make an outfile per search
                outfile = '%s.%s.hmmer' % (rje.baseFile(hmm,True),rje.baseFile(dbase,True))
                resfile = outfile
                if not os.path.exists(outfile) and self.opt['GZip'] and os.path.exists('%s.gz' % outfile) and not self.opt['Force']:
                    resfile = '%s.gz' % outfile
                if not self.opt['Force'] and rje.isYounger(resfile,hmm) == resfile and rje.isYounger(resfile,dbase) == resfile:
                    self.printLog('#HMM','HMM results file "%s" exists.' % resfile)
                    return outfile      # Already exists
                else: rje.backup(self,outfile,unlink=True)
            ### ~ [2] ~ HMM Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['HMMPFam']:
                _command = 'hmmpfam --cut_ga %s %s %s > %s' % (string.join(self.list['HMMOptions']),hmm,dbase,outfile)
            else: _command = 'hmmsearch %s %s %s > %s' % (string.join(self.list['HMMOptions']),hmm,dbase,outfile)
            self.log.printLog('#HMM',_command)
            if not wait: os.system(self.info['HMMerPath'] + _command + ' &')
            elif not os.path.exists(outfile) or self.opt['Force']: open(outfile,'a').write(os.popen(self.info['HMMerPath'] + _command).read())
            self.printLog('#HMM','Outfile produced for %s: %s.' % (hmm,outfile))
            if self.opt['GZip']:
                rje.backup(self,'%s.gz' % outfile,unlink=True)
                os.system('gzip %s' % outfile)
                self.printLog('#GZIP','%s gzipped to save space' % outfile)
            return outfile
        except:
            self.log.errorLog('Fatal Error during hmmSearch(%s)' % hmm)
            return None
#########################################################################################################################
    ### <5> ### HMM Search Results Processing and Output                                                                #
#########################################################################################################################
    def hmmTable(self,outfile='',append=False,delimit=None):    ### Outputs results table
        '''
        Outputs results table.
        >> outfile:str = Name of output file
        >> append:boolean = whether to append file
        >> delimit:str = Delimiter to use [\t]
        '''
        try:
            ### Setup ###
            if not outfile: outfile = self.info['HMMTab']
            if outfile.lower() == 'none':
                self.log.printLog('#TAB','HMMTab = "None": No table output')
                return False
            if not delimit: delimit = rje.getDelimit(self.cmd_list,'\t')
            if not outfile: outfile = '%s.hmmer.%s' % (rje.baseFile(self.info['SearchDB'],True),rje.delimitExt(delimit))
            self.readResults()
            self.log.printLog('#TAB','Tabulating results for %s searches into %s' % (len(self.search),outfile),log=False)

            ### Setup Resfile ###
            if self.opt['MySQL']: headers = ['HMM','Hit','Hit_Start','Hit_End','Eval','Score']
            else: headers = ['Type','Name','Start','End','Eval','Score']
            if not append or not os.path.exists(outfile): rje.delimitedFileOutput(self,outfile,headers,delimit,rje_backup=True)
            
            ### Output Search details ###
            for search in self.search:
                for hit in search.hit:
                    for aln in hit.aln:
                        out = {'HMM':search.info['Name'],'Type':search.info['Name'],
                               'Name':hit.info['Name'],'Hit':hit.info['Name'],
                               'Start':'%d' % aln.stat['SbjStart'], 'End':'%d' % aln.stat['SbjEnd'],
                               'Hit_Start':'%d' % aln.stat['SbjStart'], 'Hit_End':'%d' % aln.stat['SbjEnd'],
                               'Eval':'%.2e' % aln.stat['Expect'],'Score':'%.1f' % aln.stat['BitScore']}
                        rje.delimitedFileOutput(self,outfile,headers,delimit,out)
            self.log.printLog('#OUT','Results for %s searches output to %s.' % (len(self.search),outfile))
        except:
            self.log.errorLog('Fatal Error during hmmTable(%s).' % outfile)
            raise
#########################################################################################################################
    def readResults(self,clear=True,readaln=False):  ### Reads results from self.list['HMMRes'] into objects
        '''
        Reads results from self.list['HMMRes'] into objects.
        >> clear:boolean = whether to clear self.search before reading [True]
        >> readaln:boolean = whether to bother reading Alignments into objects [False]
        '''
        try:
            if clear: self.search = []
            for resfile in rje.sortUnique(self.list['HMMRes'],xreplace=False):
                if not os.path.exists(resfile) and self.opt['GZip'] and os.path.exists('%s.gz' % resfile):
                    os.system('gunzip %s.gz' % resfile)
                    self.printLog('#GUNZIP','Gunzipped %s.gz' % resfile)
                if self.opt['HMMPFam']: self.readHMMPFamSearch(resfile,readaln)
                else: self.readHMMSearch(resfile,readaln)
                if self.opt['GZip'] and os.path.exists(resfile):
                    rje.backup(self,'%s.gz' % resfile,unlink=True)
                    os.system('gzip %s' % resfile)
                    self.printLog('#GZIP','%s gzipped to save space' % resfile)
        except:
            self.log.errorLog('Hmm indeed. rje_hmm.readResults() gone awry!',quitchoice=True)
            return False
#########################################################################################################################
    def readHMMPFamSearch(self,resfile=None,readaln=False):  ### Reads HMM PFam Search Results into objects    
        '''
        Reads HMM Search Results into objects.
        >> resfile:str = Results File (set as self.info['OutFile'])
        >> readaln:boolean = whether to bother reading Alignments into objects [False] !!! Currently always False !!!
        '''
        try:
            ### Setup ###
            if not resfile or not os.path.exists(resfile):
                self.log.errorLog('Results file "%s" missing!' % resfile,printerror=False)
                return False
            ## Make RegExp for starting next alignment ##
            re_hit = string.join(['^(\S+):','domain','(\d+)','of','(\d+),','from','(\d+)','to','(\d+):','score','(\S+),','E','=','(\S+)'],'\s+')
            ## Search dictionary as results come back per sequence, not per HMM! ##
            pfam = {}   # Dictionary of {PFam name:search}
            hitx = 0    # Total number of hits
            hitlist = []        # List of sequences processed from file (may or may not include zero hit sequences)
            ### Read in Search results ###
            if open(resfile,'r').readline().find('hmmpfam') != 0:
                self.errorLog('File "%s" does not appear to be an hmmpfam results file' % resfile,printerror=False)
                if rje.yesNo('Delete incorrect results file? (Check that hmmpfam=T is right!)',default='N'):
                    os.unlink(resfile)
                    self.printLog('#DEL','Dodgy results file "%s" deleted.' % resfile)
                return False
            hitname = None
            i = 0; hx = 0; seqx = 0
            RESFILE = open(resfile,'r')
            #x#resline = self.loadFromFile(resfile,chomplines=True)
            #x#while i < len(resline):
            line = RESFILE.readline()
            newres = [rje.chomp(line)]; newresout = True; newresfile = '%s.partial' % resfile
            if os.path.exists(newresfile): os.unlink(newresfile)
            while line:
                self.progLog('\r#RES','Reading %s: %s Seqs; %s Domains; %s Hits' % (resfile,rje.integerString(hx),rje.integerString(len(pfam)),rje.integerString(hitx)))
                line = rje.chomp(line)
                #print line
                ## New Sequence ##
                if rje.matchExp('^Query sequence:\s+(\S+)',line):
                    if newres and newresout and self.opt['CleanRes']: open(newresfile,'a').write(string.join(newres,'\n'))
                    newres = ['',line]; newresout = False
                    hitname = rje.matchExp('^Query sequence:\s+(\S+)',line)[0]; hx += 1
                    #x#if hitname not in hitlist: hitlist.append(hitname)
                ## One Line Data for hits ##
                elif line.find('Parsed for domains:') == 0:
                    #x#i += 3      # Skip two complete lines
                    newres += [line,rje.chomp(RESFILE.readline()),rje.chomp(RESFILE.readline())]
                    line = rje.chomp(RESFILE.readline()); newres.append(line)
                    #Model           Domain  seq-f seq-t    hmm-f hmm-t      score  E-value
                    #--------        ------- ----- -----    ----- -----      -----  -------
                    #Lep_receptor_Ig   1/1      24   114 ..     1   103 []   158.4  1.7e-44
                    # ... else ...
                    #         [no hits above thresholds]
                    while rje.matchExp(string.join(['^(\S+)','\S+','(\d+)','(\d+)\D.+','(\S+)','(\S+)\s*$'],'\s+'),line):
                        newresout = True
                        (dom,start,end,score,eval) = rje.matchExp(string.join(['^(\S+)','\S+','(\d+)','(\d+)\D.+','(\S+)','(\S+)\s*$'],'\s+'),line)
                        if not pfam.has_key(dom):
                            pfam[dom] = self._addSearch()
                            pfam[dom].info['Name'] = dom
                        hit = pfam[dom]._addHit()
                        hit.info['Name'] = hitname
                        aln = hit._addAln()
                        aln.setStat({'SbjStart':string.atoi(start),'SbjEnd':string.atoi(end),'Expect':string.atof(eval),'BitScore':string.atof(score)})
                        hitx += 1
                        self.progLog('\r#RES','Reading %s: %s Seqs; %s Domains; %s Hits' % (resfile,rje.integerString(hx),rje.integerString(len(pfam)),rje.integerString(hitx)))
                        line = rje.chomp(RESFILE.readline()); newres.append(line)
                ## End of Protein ##
                elif line[:2] == '//': hitname = None; newres.append(line)
                elif rje.matchExp('End of rje_hmm reduced results file: (%d) sequences in original',line):
                    seqx = string.atoi(rje.matchExp('End of rje_hmm reduced results file: (\d+) sequences in original',line)[0])
                elif newres: newres.append(line)
                #x#i += 1
                line = RESFILE.readline()
            if newres and newresout and self.opt['CleanRes']: open(newresfile,'a').write(string.join(newres,'\n'))
            if not seqx: seqx = hx
            if self.opt['CleanRes']:
                open(newresfile,'a').write(string.join(['','End of rje_hmm reduced results file: %d sequences in original' % seqx],'\n'))
                os.unlink(resfile)
                os.rename(newresfile,resfile)
                self.printLog('\r#RED','Results file %s replaced with reduced version (%s Hits only)' % (resfile,rje.integerString(hitx)))
            self.printLog('\r#RES','Reading %s complete: %s Seqs; %s Domains; %s Hits' % (resfile,rje.integerString(seqx),rje.integerString(len(pfam)),rje.integerString(hitx)))
            return True
        except:
            self.log.errorLog('Calamity during readHMMSearch(%s)' % (resfile))
            return False
#########################################################################################################################
    def readHMMSearch(self,resfile=None,readaln=False):  ### Reads HMM Search Results into objects    #!# Needs tidying! #!#
        '''
        Reads HMM Search Results into objects.
        >> resfile:str = Results File (set as self.info['OutFile'])
        >> readaln:boolean = whether to bother reading Alignments into objects [False] (!!!currently always True!!!)
        '''
        try:
            ### <a> ### Setup
            _stage = '<a> Setup'
            #print resfile
            if not resfile or not os.path.exists(resfile):
                self.log.errorLog('Results file (%) missing!' % resfile,False,False)
                raise IOError
            _hit_elements = ['^(\S+):','domain','(\d+)','of','(\d+),','from','(\d+)','to','(\d+):','score','(\S+),','E','=','(\S+)']
            _hit_re = string.join(_hit_elements,'\s+')

            ### <b> ### Read in Search results
            _stage = '<b> Read Results'
            self.verbose(0,4,'Reading %s HMMer search results' % resfile,0)
            RESFILE = open(resfile, 'r')
            lines = RESFILE.readlines()
            RESFILE.close()
            resline = []
            for line in lines:
                resline.append(re.sub('\n','',line))
            search = None
            i = 0
            hitaln = 0
            if resline[i].find('hmmsearch') != 0:
                self.log.errorLog("File %s does not appear to be an hmmsearch results file" % resfile)
                raise

            while i < len(resline):
                line = resline[i]
                #print line
                ## <i> ## Basic Search Info
                _stage = '<b-i> Basic Search Info'
                if line.find('HMM file:') == 0:
                    search = self._addSearch()
                    search.info['Name'] = rje.matchExp('HMM file:\s+(\S+)',line)[0]
                    self.verbose(0,4,'.',0)
                    self.verbose(1,3,'\n%s' % search.info['Name'],0)
                elif line.find('Sequence database:') == 0:
                    search.info['SearchDB'] = rje.matchExp('Sequence database:\s+(\S+)', line)[0]
                elif line.find('Total sequences searched:') == 0:
                    dbnum = rje.matchExp('Total sequences searched:\s+(\d\S*)', line)[0]
                    dbnum = re.sub('\D','',dbnum)
                    search.stat['DBNum'] = string.atoi(dbnum)
                ## <ii> ## One-line hit data (BLASTHit)
                elif line.find('Scores for complete sequences') == 0: # One-line hits
                    _stage = '<b-ii> One-line hits'
                    i += 3  # Skip two lines
                    while re.search('^(\S+)\s.+\s(\S*\d)\s+(\S*\d)\s+(\d+)\s*$',resline[i]):
                        match = rje.matchExp('^(\S+)\s.+\s(\S*\d)\s+(\S*\d)\s+\d+\s*$',resline[i])
                        self.verbose(2,3,'\n - %s (%s, %s)' % match,0)
                        hit = search._addHit()
                        hit.info['Name'] = match[0]
                        hit.stat['BitScore'] = string.atof(match[1])
                        #print hit.stat['BitScore'], resline[i], match
                        eval = match[2]
                        if eval.find('e') == 0:
                            eval = '1' + eval
                        hit.stat['E-Value'] = string.atof(eval)
                        i += 1
                    line = resline[i]   # End of one-lines (blank line)
                    self.verbose(1,3,'=> %d Hits' % search.hitNum(),1)
                    hitaln = 0
                #!# Make new No hits pattern match
                elif line.find('***** No hits found ******') >= 0:  # No Hits
                    search.hit = []
                    self.verbose(1,3,'=> %d Hits' % search.hitNum(),1)
                    hitaln = 0
                ## <iii> ## Aln Hit data (PWAln)
                #!# Consider reading in the 'parsed for domains' section instead/as well
                elif re.search(_hit_re,line):   # New aln hit
                    _stage = '<b-iii> Aln Hit Info'
                    # Identify hit object
                    _hit_detail = rje.matchExp(_hit_re,line)
                    #print _hit_detail
                    hitname = _hit_detail[0]
                    #hitaln += 1 - string.atoi(_hit_detail[1])
                    #print hitname
                    try:
                        #if hitname != search.hit[hitaln].info['Name']:
                        for hit in search.hit:
                            if hit.info['Name'] == hitname:
                                hitaln = search.hit.index(hit)
                        if hitname != search.hit[hitaln].info['Name']:
                            self.log.errorLog('Problem with HMM results %s - %s single-line hits and alignments do not match' % (hitname,search.info['Name']),printerror=False,quitchoice=True)
                            i += 1
                            continue
                    except:
                        self.log.errorLog('Problem with HMM results reconciling %s - %s single-line hits and alignments.' % (hitname,search.info['Name']),True,True)
                        i += 1
                        continue                        
                    hit = search.hit[hitaln]
                    #print hit
                    hitaln += 1
                    # Add details
                    _stage = '<b-iii> Add Aln Hit Info'
                    aln = hit._addAln()
                    aln.stat['SbjStart'] = string.atoi(_hit_detail[3])
                    aln.stat['SbjEnd'] = string.atoi(_hit_detail[4])
                    aln.stat['BitScore'] = string.atof(_hit_detail[5])
                    aln.stat['Expect'] = string.atof(_hit_detail[6])
                    ## <iv> ## Alignments
                    readaln = True
                    i += 1
                    while readaln:
                        _stage = '<b-iv> Read alignments'
                        line = resline[i]
                        #print line
                        block = rje.matchExp('^(\s+)(\S+)',line)
                        #print block
                        if block:
                            # Query Line
                            leadlen = len(block[0])
                            seqblock = block[1]
                            #print block, leadlen, (leadlen+len(seqblock))
                            if block[1][:3] == '*->':    # Start
                                leadlen += 3
                                #print seqblock[3:]
                                seqblock = seqblock[3:]
                            if block[1][-3:] == '<-*':  # End
                                #print seqblock[:-3]
                                seqblock = seqblock[:-3]
                                readaln = False
                            #print block, leadlen, (leadlen+len(seqblock))
                            aln.info['QrySeq'] += seqblock
                            # Alignment Line
                            i += 1
                            aln.info['AlnSeq'] += resline[i][leadlen:(leadlen+len(seqblock))]
                            # Subject Line
                            i += 1
                            aln.info['SbjSeq'] += resline[i][leadlen:(leadlen+len(seqblock))]
                            # Skip Blank line
                            i += 2
                        else:
                            #print 'This should be a block!:\n', line
                            i += 1

                i += 1
            #print self.search
            #print self.search[0].hit
            #print self.search[0].hit[0].aln
            self.verbose(0,1,'Reading of %s HMM results complete! (%d Searches)' % (resfile,len(self.search)),2)
            return True
        except:
            self.log.errorLog('Calamity during readHMMSearch(%s) %s.' % (resfile,_stage))
            return False
#########################################################################################################################
    def _addSearch(self):   ### Adds and returns a new search object
        '''Adds and returns a new search object.'''
        newsearch = HMMSearch(log=self.log)
        self.search.append(newsearch)
        return newsearch
#########################################################################################################################
### END OF SECTION I: HMMRun                                                                                            #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: HMMSearch Class                                                                                        #
#########################################################################################################################
class HMMSearch(rje_blast.BLASTSearch):     
    '''
    HMM Search Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Name of HMM Search File
    - DBase = Name of database File Searched
    
    Opt:boolean

    Stat:numeric
    - DBNum = Number of sequences in database searched

    Obj:RJE_Objects

    Other:
    hit : list of HMMHit Objects
    '''
    ### Attributes
    hit = []
    def hitNum(self): return len(self.hit)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','DBase']
        - Stats:float ['DBNum']
        - Opt:boolean []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','DBase']
        for info in self.infolist:
            self.info[info] = 'None'
        self.statlist = ['DBNum']
        for stat in self.statlist:
            self.stat[stat] = 0.0
        self.optlist = []
        for opt in self.optlist:
            self.opt[opt] = False
        self.objlist = []
        for obj in self.objlist:
            self.obj[obj] = None
        ### <b> ### Defaults
        ### <c> ### Other Attributes
        self.hit = []
#########################################################################################################################
    def _addHit(self):   ### Adds and returns a new hit object
        '''Adds and returns a new hit object.'''
        newhit = HMMHit(log=self.log)
        self.hit.append(newhit)
        return newhit
#########################################################################################################################
### Methods inherited from BLASTSearch:
#    def hitSeq(self,seqlist):   ### Returns dictionary of hits as sequences from seqlist (or None if missing)
#    def saveHitIDs(self,outfile=None):   ### Saves IDs of hits in file (e.g. for later fastacmd extraction)
#########################################################################################################################
### END OF SECTION III: HMMSearch                                                                                       #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: HMMHit Class                                                                                            #
#########################################################################################################################
class HMMHit(rje.RJE_Object):     
    '''
    BLAST Hit Class. Author: Rich Edwards (2005).

    Info:str
    - Name = Short Name of Hit Sequence (First Word of description line)
    
    Opt:boolean

    Stat:numeric
    - BitScore
    - E-Value
    - Length

    Obj:RJE_Objects

    Other:
    - aln:list of PWAln Objects
    '''
    ### Attributes
    aln = []    # List of PWAln Objects
    def alnNum(self): return len(self.aln)
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''
        Sets Attributes of Object:
        - Info:str ['Name','Type','Description']
        - Stats:float ['BitScore','E-Value','Length']
        - Opt:boolean []
        - Obj:RJE_Object []
        '''
        ### <a> ### Basics 
        self.infolist = ['Name','Type','Description']
        for info in self.infolist:
            self.info[info] = 'None'
        self.statlist = ['BitScore','E-Value','Length']
        for stat in self.statlist:
            self.stat[stat] = 0.0
        self.optlist = []
        for opt in self.optlist:
            self.opt[opt] = False
        self.objlist = []
        for obj in self.objlist:
            self.obj[obj] = None
        ### <b> ### Defaults
        ### <c> ### Other Attributes
        self.aln = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### <a> ### General Options
                self._generalCmd(cmd)
                ### <b> ### Class Options
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
        return
#########################################################################################################################
    def _addAln(self):  ### Adds and returns a PWAln Object
        '''
        Adds and returns a PWAln Object.
        '''
        newaln = rje_blast.PWAln(log=self.log)
        self.aln.append(newaln)
        return newaln
#########################################################################################################################
### END OF SECTION IV: BLASTHit                                                                                         #
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
    try: HMMRun(mainlog,cmd_list).run()

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
