#!/usr/bin/python

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
Module:       rje_markov
Description:  RJE Module for protein sequence Markov Chain related gubbins
Version:      2.2
Last Edit:    17/01/13
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to handle a number of things to do with amino acid frequencies and the like, including:
    - Observed and expected 1mer .. Xmer frequencies

    Some old functions including 1mer .. Xmer counts per protein are now only found in V1.2.

Commandline:
    ## General ##
    seqin=FILE          : File with input sequences [None]
    alphabet=X,Y,..Z    : List of letters in alphabet of interest
    split=X             : Splits file into numbered files of X sequences and recombines at end.
    autoload=T/F        : Whether to load sequences automatically. If False, will try memory-efficient seqlist handling [False]
    aafreq=FILE         : Generate expected 1mer frequencies from FILE of aafreq [None]
    xmerfile=FILE       : File from which to load Xmer counts and build suffix tree [None]
    direction=X         : Direction to read chains (fwd/bwd/both) [fwd]
    sorted=T/F          : Whether to use sorted xmer fragments to reduce memory requirements [False] 
    scap=T/F            : Whether to use special SCAP sorting of xmers (sorts all but last aa) [False]
    negvpos=X           : Perform Negatives versus Positives analysis on X.* files [None]

    ## Output ##
    markov=T/F      : Whether to perform Markov Chain Analysis of observed vs Expected for Xmers [True]
    xmers=X[,Y]     : Deal with Xmers [upto Ymers] [1]

Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje, rje_seq
Other modules needed: rje_blast, rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Basic working model, with suffix tree
    # 1.1 - Added splitting. Committed to suffix tree. Removed autoload and suftree options.
    # 1.2 - Removed obs/exp calculation from output - can calculate it from frequencies. Removed redundant methods.
    # 2.0 - Added reverse reading of markov chains and general tidy up. Removed XCount method - use V1.2 for this.
    # 2.1 - Temp fixed suftree problem. Module needs a major overhaul! Still broken in places.
    # 2.2 - Minor modifications. Main methods still need checking/fixing.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Reinstate AACount using suffix trees.
    # [ ] : Add reading of reverse markovs (prob given following sequences)
    # [ ] : Finish tidying with self.dict and self.list organisation.
    # [ ] : Sorted=T option is all messed up - add properly later
    # [ ] : Needs a major overhaul!
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_MARKOV', '2.1', 'June 2011', '2005')
    description = 'RJE Markov Chain Gubbins Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
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
### SECTION II: Markov Class                                                                                            #
#########################################################################################################################
class Markov(rje.RJE_Object):     
    '''
    Markov Sequence Class. Author: Rich Edwards (2005).

    Info:str
    - Name = File from which Markov chains last generated
    - AAFreqFile = File of AA Frequencies from which to generate expected 1mer frequencies
    - Direction = Direction to read chains (fwd/bwd/both)
    - XmerFile = File from which to load Xmer counts and build suffix tree
    - NegVPos = Perform Negatives versus Positives analysis on X.* files [None]
    
    Opt:boolean
    - AutoLoad = Whether to load sequences automatically. If False, will try memory-efficient seqlist handling [False]
    - Markov = Whether to perform Markov Chain Analysis of observed vs Expected for Xmers [True]
    - SCAP = Whether to use special SCAP sorting of xmers (sorts all but last aa) [False] 
    - Silent = Whether to cancel log output [False]
    - Sorted = Whether to use sorted xmer fragments to reduce memory requirements [False]

    Stat:numeric
    - MinXmer = Minimum length Xmers for 
    - MaxXmer = Maximum length Xmers for 
    - Split = Splits file into numbered files of X sequences and recombines at end.

    List:list
    - Alphabet = Alphabet allowed in sequences [standard 1 letter codes]

    Dict:dictionary
    - AAFreq = Dictionary of aa frequencies read from input sequences
    - SufTree = Suffix tree dictionary. Keys: AAs or '=' for counts, '+' for sum / xmers lengths with '=' and 'X' for no. dif xmers
    - PreTree = Prefix tree dictionary (Reversed Suffix Tree)
    
    Obj:RJE_Objects
    - DB = rje_db.Database object
    - SeqList = rje_seq.SeqList object
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['AAFreqFile','XmerFile','Direction','NegVPos']
        self.optlist = ['AutoLoad','Markov','Silent','Sorted','SCAP']
        self.statlist = ['MinXmer','MaxXmer','Split']
        self.listlist = ['Alphabet']
        self.dictlist = ['AAFreq','SufTree','PreTree']
        self.objlist = ['SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=True,stat=1.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'Direction':'fwd'})
        self.setOpt({'AutoLoad':False,'Silent':False,'Sorted':False,'SCAP':False})
        self.stat['Split'] = 0
        self.obj['SeqList'] = rje_seq.SeqList(log=self.log,cmd_list=['autoload=F']+self.cmd_list+['autofilter=F'])
        if self.obj['SeqList'].opt['DNA']: self.list['Alphabet'] = rje_seq.alph_dna[:-1]
        else: self.list['Alphabet'] = rje_seq.alph_protx[:-1]
        self.cmd_list = ['delimit=,'] + self.cmd_list
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
                ### <a> ### General Options
                self._generalCmd(cmd)
                self._forkCmd(cmd)  # Delete if no forking
                ### <b> ### Class Options
                self._cmdRead(cmd,type='info',att='AAFreqFile',arg='aafreq')
                self._cmdReadList(cmd,'info',['XmerFile','Direction','NegVPos'])
                self._cmdRead(cmd,type='list',att='Alphabet')
                self._cmdReadList(cmd,'opt',['Markov','Silent','AutoLoad','Sorted','SCAP'])
                self._cmdRead(cmd,type='int',att='Split')
                self._cmdRead(cmd,type='min',att='MinXmer',arg='xmers')                
                self._cmdRead(cmd,type='max',att='MaxXmer',arg='xmers')                
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Run Methods                                                                                        #
#########################################################################################################################
    #def verbose(self,a,b,c,d): return
    def suftree(self): return self.dict['SufTree']
    def pretree(self): return self.dict['PreTree']
    def sufOpt(self): return self.info['Direction'].lower() in ['suf','fwd','both','suffix','forward','forwards']
    def preOpt(self): return self.info['Direction'].lower() in ['pre','bwd','both','prefix','backward','backwards']
    def _xmerCountFromSufTree(self,xmer,prefix=False): return self.xmerCount(xmer)   ### Returns count for given xmer from self.dict tree
#########################################################################################################################
    def run(self):      ### Main Run Method
        '''Main Run Method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['NegVPos'].lower() not in ['','none']: return self.negVPos()
            self.setup()
            ### ~ [2] ~ Whole dataset analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['Markov'] and self.stat['Split'] > 0: self.splitFileMarkov()
            else:
                self.buildMarkov(self.obj['SeqList'])
                if self.opt['Markov']: self.markovAnalysis()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def setup(self):      ### Main Run setup Method
        '''Main Run setup Method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if self.opt['Sorted']:
            #    self.opt['Sorted'] = False
            #    self.printLog('#SORTED','Sorted Markov option currently unavailable!')
            self.obj['SeqList'].makeBaseFile()
            if not self.opt['AutoLoad']: self.opt['MemSaver'] = True
            if rje.checkForFile(self.info['AAFreqFile']):
                self.dict['AAFreq'] = self.obj['SeqList'].aaFreq(alphabet=self.list['Alphabet'],loadfile=self.info['AAFreqFile'])
            self.stat['MinXmer'] = int(self.stat['MinXmer'])
            self.stat['MaxXmer'] = max(self.stat['MinXmer'],int(self.stat['MaxXmer']))
            self.info['Basefile'] = '%s.%dmer' % (self.obj['SeqList'].info['Basefile'],self.stat['MaxXmer'])
            if self.opt['Sorted']:
                self.stat['MinXmer'] = self.stat['MaxXmer']      # Cannot look at lower Xmers
                self.info['Basefile'] = '%s.sorted' % self.info['Basefile'] 
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <3> ### Build Suffix Tree Methods                                                                               #
#########################################################################################################################
    def _buildSuffixTree(self,seqlist,partial=False): return self.buildMarkov(seqlist,partial)
    def buildMarkov(self,seqlist,partial=False):    ### Uses self.attributes to make the self.dict tree dictionaries
        '''
        Uses self.attributes to make the self.dict tree dictionaries.
        >> seqlist:SeqList object on which to perform analysis
        >> partial:boolean = Whether this is a partial analysis to be recombined later
        << returns name of output file. (Useful if data split)
        '''
        try:### ~ [1] Setup & Checks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Some duplication of self._setup() removed #i#
            ## ~ [1a] Input Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = self.obj['SeqList']
            if self.opt['MemSaver']:
                if not os.path.exists(seqlist.info['Name']):
                    self.log.errorLog('File <<%s>> does not exist!' % seqlist.info['Name'],False,False)
                    raise IOError
                SEQFILE = open(seqlist.info['Name'], 'r')
                lastline = 'Start'
            elif seqlist.seqNum() < 1:
                seqlist.loadSeqs()                               
                seqlist.degapSeq()
            self.info['Name'] = seqlist.info['Name']
            ## ~ [1b] ~ Global variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for tree in ['SufTree','PreTree']:
                self.dict[tree] = {'=':0}
                for x in range(self.stat['MaxXmer']): self.dict[tree][x+1] = {'=':0,'X':0}            
            ## ~ [1c] ~ Method variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            max_xmer = int(self.stat['MaxXmer'])
            
            ### ~ [2] ~ Process Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0
            if self.opt['Sorted']: xtxt = 'Building Trees for Sorted %dmers' % max_xmer
            else: xtxt = 'Building Trees for <= %dmers' % max_xmer
            while self.opt['MemSaver'] or sx < seqlist.seqNum():
                self.progLog('\r#XMERS', '%s in %s (%s seq)' % (xtxt,seqlist.info['Name'],rje.integerString(sx)))
                ## ~ [2a] ~ Get Seq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ## 
                if self.opt['MemSaver']:
                    (seq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
                    if seq: seqlist.seq = [seq]
                    else: break
                else: seq = seqlist.seq[sx]    
                ## ~ [2b] ~ Build suftree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self._addToTrees(seq.info['Sequence'],max_xmer)
                sx += 1
            self.printLog('\r#XMERS', '%s in %s (%s seq) complete!' % (xtxt,seqlist.info['Name'],rje.integerString(sx)))
            if self.opt['MemSaver']: SEQFILE.close()

            ### ~ [3] ~ Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for x in range(max_xmer):
                self.printLog('#SUFF','%s counts of %s different Suffix %dmers.' % (rje.integerString(self.suftree()[x+1]['=']),rje.integerString(self.suftree()[x+1]['X']),x+1))
                self.printLog('#PREF','%s counts of %s different Prefix %dmers.' % (rje.integerString(self.pretree()[x+1]['=']),rje.integerString(self.pretree()[x+1]['X']),x+1))
        except: self.errorLog('Error in Markov.buildMarkov()',printerror=True,quitchoice=False); raise
#########################################################################################################################
    def markovAnalysis(self,partial=False): ### Calculates Observed vs Expected etc.
        '''Calculates Observed vs Expected etc.'''
        try:#!# Not modified in V2 #!#
            ## ~ [1c] ~ Method variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqlist = self.obj['SeqList']
            alphabet = self.list['Alphabet']
            max_xmer = int(self.stat['MaxXmer'])
            ## ~ [1d] Output File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            delimit = rje.getDelimit(self.cmd_list)
            filename = '%s.markov.%s' % (seqlist.info['Basefile'],rje.delimitExt(delimit))
            ### Markov Chain Analysis ###
            self.verbose(0,4,'Markov Chain Analysis for <= %dmers...' % max_xmer,0)
            if partial:
                MARKOV = open(filename,'w')
                MARKOV.write('aacount=%d\n' % self.suftree()['='])
                MARKOV.write('seqnum=%d\n' % sx)
                outlist = ['xmer','obs','pref_obs','suff_obs','comm_obs']
                MARKOV.write('%s\n' % string.join(outlist,delimit))
                MARKOV.close()
            else:
                for x in range(max_xmer):
                    xfile = '%s.%dmer%s' % (rje.baseFile(filename),(x+1),os.path.splitext(filename)[1])
                    MARKOV = open(xfile,'w')
                    outlist = ['xmer','obs','f_obs','f_exp']
                    MARKOV.write('%s\n' % string.join(outlist,delimit))
                    MARKOV.close()
            self._markovSufTree(self.suftree(),'',filename,delimit,seqlist.seqNum(),partial)  # sufdic,xmer,filename,delimit,seqnum,partial,typex=0)
            self.verbose(0,1,'Done!',1)
            if not self.opt['Silent']: self.log.printLog('#INF','Markov Chain Analysis of %s for <= %dmers complete!' % (seqlist.info['Name'],max_xmer))
            return filename
        except:
            self.errorLog('Error in _buildSuffixTree',printerror=True,quitchoice=True)
            raise
#########################################################################################################################
    def _addToTrees(self,sequence,max_xmer): ### Adds markov chains from sequence into self.suftree
        '''
        Adds markov chains from sequence into self.suftree.
        >> sequence:str = sequence to add
        >> max_xmer:int = maximum size of Xmer to add
        '''
        for r in range(len(sequence)):   # Progressive starting position
            if self.sufOpt(): self._addXmerToSufTree(sequence[r:(r+max_xmer)])
            if self.preOpt(): self._addXmerToSufTree(rje.strReverse(sequence[max(0,(r+1-max_xmer)):r+1]),prefix=True)
#########################################################################################################################
    def _addXmerToSufTree(self,xmer,count=1,prefix=False):  ### Adds a single Xmer to self.dict trees
        '''
        Adds a single Xmer to self.dict trees.
        >> xmer:str = Xmerto add
        >> count:int = Number of occurrences of Xmer
        >> prefix:bool [False] = Add to prefix tree rather than suffix tree
        '''
        ### ~ [1] ~ Setup dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['Sorted']:
            if len(xmer) < self.stat['MaxXmer']: return
            if self.opt['SCAP']: xmer = rje.strSort(xmer[:-1]) + xmer[-1]
            else: xmer = rje.strSort(xmer)
        if not xmer or xmer[0].upper() not in self.list['Alphabet']: return
        if prefix: _subdic = _sufdic = self.pretree()
        else: _subdic = _sufdic = self.suftree()
        ### ~ [2] ~ Iterate through ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        _sufdic['='] += count
        for a in range(len(xmer)):   # Expanding Xmer
            aa = xmer[a].upper()
            if aa not in self.list['Alphabet']: return
            if _sufdic.has_key(a+1): _sufdic[a+1]['='] += count
            else: _sufdic[a+1] = {'=':count,'X':0}
            if _subdic.has_key(aa):     # Xmer exists - add count
                _subdic[aa]['='] += count
            else:                       # New xmer of length a+1
                _subdic[aa] = {'=':count}
                _sufdic[a+1]['X'] += 1
            _subdic = _subdic[aa]
#########################################################################################################################
    def readXmers(self,xfile=None,clear=True):    ### Uses xmer file content to make the self.suftree dictionary
        '''
        Uses xmer file content to make the self.suftree dictionary.
        >> xfile:str = Name of input file
        >> clear:boolean = Whether to clear self.suftree before adding xmers
        '''
        try:#!# Needs to be updated #!#
            ### Setup & Checks ###
            ## Suffix Tree ##
            if clear:
                self.dict['SufTree'] = {'=':0}
            ## File ##
            if not xfile:
                xfile = self.info['XmerFile']
            if not os.path.exists(xfile):
                self.log.errorLog('File %s missing!' % xfile,False,False)
                raise IOError
            delimit = rje.getDelimit(self.cmd_list)

            ### Process Sequences ###
            self.verbose(0,4,'\nBuilding suffix tree from %s Xmers...' % xfile,0)
            self.info['Name'] = xfile
            XFILE = open(xfile,'r')
            header = rje.readDelimit(XFILE.readline(),delimit)
            if header[0] != 'xmer':
                self.log.errorLog('File %s in wrong format? Header missing or incorrect delimiter given!' % xfile,False,False)
                self.log.errorLog('Header info read: %s' % header,False,True)
                raise ValueError
            line = XFILE.readline()
            sx = 0
            while line not in [None,'']:
                readlist = rje.readDelimit(line,delimit)
                line = XFILE.readline()
                xmer = readlist[0]
                count = string.atoi(readlist[1])
                ## Build suftree ##
                self._addXmerToSufTree(xmer,count)
                ## Next seq
                sx += 1
                rje.progressPrint(self,sx,dotx=1000,numx=10000)
            self.verbose(0,2,'Done!',1)
            XFILE.close()

            ### Summarise ###
            max_xmer = self.getMaxXmer(self.suftree())
            for x in range(max_xmer):
                if not self.opt['Silent']: self.log.printLog('#XMERS','%s counts of %s different %dmers from %s.' % (rje.integerString(self.suftree()[x+1]['=']),rje.integerString(self.suftree()[x+1]['X']),(x+1),xfile))
            if not self.opt['Silent']: self.log.printLog('#NB','Xmers < %dmers will be under-represented as sequence end information missing.' % max_xmer)
            if clear or self.stat['MaxXmer'] < max_xmer:
                self.stat['MaxXmer'] = max_xmer
        except:
            self.log.errorLog('Error in readXmers',printerror=True,quitchoice=True)
            raise
#########################################################################################################################
    def _aaFreqToSufTree(self): ### Uses self.dict['AAFreq'] to build self.suftree
        '''
        Uses self.dict['AAFreq'] to build self.suftree (and self.probtree).
        '''
        #!# Needs to be updated #!#
        self.dict['SufTree'] = {'=':0,'X':len(self.list['Alphabet'])}
        total = 0.0
        for aa in self.list['Alphabet']:
            if self.dict['AAFreq'].has_key(aa):
                total += float(self.dict['AAFreq'][aa])
            else:
                self.dict['AAFreq'][aa] = 0.0
        self.suftree()['='] = 1.0
        for aa in self.list['Alphabet']:
            self.suftree()[aa] = {'=':self.dict['AAFreq'][aa] / total}
        #?!#self._buildProbTree()
        self.info['Name'] = self.info['AAFreqFile']
#########################################################################################################################
    ### <4> ### Removal of Xmers from Suffix Tree Methods                                                               #
#########################################################################################################################
    def _removeFromSufTree(self,sequence,max_xmer): ### Removes markov chains from sequence from self.suftree
        '''
        Removes markov chains from sequence from self.suftree.
        >> sequence:str = sequence to remove
        >> max_xmer:int = maximum size of Xmer to remove
        '''
        #!# Needs to be updated #!#
        for r in range(len(sequence)):   # Progressive starting position
            self._removeXmerFromSufTree(sequence[r:(r+max_xmer)])
#########################################################################################################################
    def _removeXmerFromSufTree(self,xmer): ### Removes a single Xmer from self.suftree
        '''
        Removes a single Xmer from self.suftree.
        >> xmer:str = Xmer to remove
        '''
        #!# Needs to be updated #!#
        try:
            _sufdic = self.suftree()
            _sufdic['='] -= 1
            for a in range(len(xmer)):   # Expanding Xmer
                aa = xmer[a].upper()
                self.suftree()[(a+1)]['='] -= 1
                #print aa, 'in', _sufdic.keys(), '?'
                if _sufdic.has_key(aa): # Xmer exists - reduce count
                    #print aa, _sufdic[aa]['='], '=>',
                    _sufdic[aa]['='] -= 1
                    #print _sufdic[aa]['=']
                else:   # New xmer of length a+1. Shouldn't happen when removing!
                    _sufdic[aa] = {'=':-1}
                    self.log.errorLog('Removing non-existent Xmer %s!' % xmer,False,False)
                _sufdic = _sufdic[aa]
        except:
            self.deBug('%s: %s' % (xmer,self.suftree()))
            raise
#########################################################################################################################
    def stripToX(self,sufdic,maxmer,xlen=0): ### Strips a suffix tree down to max Xmer length maxmer
        '''
        Strips a suffix tree down to max Xmer length maxmer.
        >> sufdic:suffix tree dictionary
        >> maxmer:int max length of xmer
        >> xlen:int = length of xmer
        '''
        #!# Needs to be updated #!#
        if xlen == 0:
            for r in range(maxmer,self.getMaxXmer(sufdic)):
                sufdic.pop(r+1)
        if xlen == maxmer:
            sufdic = {'=':sufdic['=']}
        else:
            for aa in self.list['Alphabet']:
                if sufdic.has_key(aa):
                    self.stripToX(sufdic[aa],maxmer,xlen+1)
        #print xlen,sufdic
#########################################################################################################################
    ### <5> ### Interrogate Suffix Tree Methods                                                                         #
#########################################################################################################################
    def getMaxXmer(self,suftree):   ### Returns max Xmer size from suftree
        '''
        Returns max Xmer size from suftree.
        >> suftree:suffix tree dictionary
        '''
        max = 0
        while suftree.has_key(max+1): max += 1
        return max
#########################################################################################################################
    def xmerCount(self,xmer,prefix=False,prob=False):    ### Returns count for given xmer from self.dict tree
        '''
        Returns count for given xmer from self.dict tree.
        >> xmer:str = Xmer of interest
        >> prefix:bool [False] = Use Prefix tree rather than suffix tree
        >> prob:bool [False] = Return probability rather than count
        '''
        ### ~ [1] ~ Choose tree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.opt['Sorted'] and len(xmer) < self.stat['MaxXmer']: return {True:1.0,False:0}[prob]
        if prefix: _sufdic = self.pretree(); xmer = rje.strReverse(xmer)
        else: _sufdic = self.suftree()
        if self.opt['Sorted']:
            if self.opt['SCAP']: xmer = rje.strSort(xmer[:-1]) + xmer[-1]
            else: xmer = rje.strSort(xmer)
        prex = 0
        self.deBug('%s :: %s' % (xmer,rje.sortKeys(_sufdic)))
        ### ~ [2] ~ Find subtree ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for x in range(len(xmer)):
            if xmer[x] in _sufdic.keys():
                prex = _sufdic['=']; _sufdic = _sufdic[xmer[x]]
            else: break
            self.deBug('%s %d [%s] :: %d = %s' % (xmer,x,xmer[x],prex,_sufdic))
        if prob:
            if prex: return _sufdic['='] / float(prex)
            else: return 1.0
        else: return _sufdic['=']
#########################################################################################################################
    def _markovSufTree(self,sufdic,xmer,filename,delimit,seqnum,partial,typex=0):  ### Performs Observed vs Expected analysis for xmer
        '''
        Performs Observed vs Expected analysis for xmer.
        >> sufdic:suftree level
        >> xmer:str = basal xmer to extend with alphabet
        >> filename:str = output file to append
        >> delimit:str = delimiter for file
        >> partial:boolean = Whether this is a partial analysis to be recombined later
        << returns typex = the number of xmer types
        '''
        #!# Not modified yet #!#
        if sufdic['='] == 0:     # No entries at this level
            return typex
        for aa in self.list['Alphabet']:
            ## This xmer ##
            if aa not in sufdic.keys():
                continue
            _xmer = xmer + aa
            (pref,suff,comm) = (_xmer[:-1],_xmer[1:],_xmer[1:-1])
            num = {}
            ## Counts ##
            for _str in [_xmer,pref,suff,comm]:
                num[_str] = self._xmerCountFromSufTree(_str)
            ## Output ##
            if partial:
                MARKOV = open(filename,'a')
                outlist = [_xmer]
                for _str in [_xmer,pref,suff,comm]:
                    outlist.append('%d' % num[_str])
                MARKOV.write('%s\n' % string.join(outlist,delimit))
                MARKOV.close()
            else:
                xfile = '%s.%dmer%s' % (rje.baseFile(filename),len(_xmer),os.path.splitext(filename)[1])
                self._markovOutput(xfile,_xmer,num,self.suftree()['='],seqnum,delimit)
            typex += 1
            if typex/50000 == typex/50000.0:
                self.verbose(0,4,rje.integerString(typex),0)
            elif typex/1000 ==typex/1000.0:
                self.verbose(0,4,'.',0)
            ## Extend further ##
            typex = self._markovSufTree(sufdic[aa],_xmer,filename,delimit,seqnum,partial,typex)
        return typex
#########################################################################################################################
    def _markovOutput(self,filename,xmer,numdic,aacount,seqnum,delimit):    ### Outputs observed vs expected to file
        '''
        Outputs observed vs expected to file.
        >> filename:str = output file name
        >> xmer:str = Xmer for output
        >> numdic:dictionary of observed counts for xmer and substrings
        >> aacount:int = Total aa count of dataset
        >> seqnum:int = Number of sequences in dataset
        >> delimit:str = delimiter for file
        '''
        #!# Not modified yet #!#
        ### Setup ###
        (pref,suff,comm) = (xmer[:-1],xmer[1:],xmer[1:-1])
        freq = {}
        ### Frequency Data ###
        for _str in [xmer,pref,suff,comm]:
            xcount = float(aacount - ((len(_str) - 1) * seqnum))
            if _str == '' or xcount == 0:
                freq[_str] = 1.0
            else:
                freq[_str] = numdic[_str] / xcount
        ### Expected ###
        if len(xmer) > 1:
            f_expected = freq[pref] * freq[suff] / freq[comm]
        else:
            if self.dict['AAFreq'].has_key(xmer):
                f_expected = self.dict['AAFreq'][xmer]
            else:
                f_expected = freq[xmer] #!# Add expected AA freq option later
        ### Output ###
        MARKOV = open(filename,'a')
        outlist = [xmer]
        outlist.append('%d' % numdic[xmer])
        outlist.append('%.3e' % freq[xmer])
        outlist.append('%.3e' % f_expected)
        #outlist.append('%f' % (freq[xmer] / f_expected))
        MARKOV.write('%s\n' % string.join(outlist,delimit))
        MARKOV.close()
#########################################################################################################################
    ### <6> ### Suffix Tree => Markov Probability Tree Methods                                                          #
#########################################################################################################################
    #!# >>> This section removed - can calculate directly from '=' data - now in self._xmerCountFromSufTree()
#########################################################################################################################
    ### <7> ### Markov Probability Tree => Suffix Tree Methods                                                          #
#########################################################################################################################
    def _sufTreeFromProbTree(self,total_count):    ### Populates sufdic using probdic and given total number of letters
        '''
        Populates sufdic using probdic and given total number of letters. Uses most common letters to fill in remainders,
        selecting a random one in the case of a tie.
        >> total_count:int = Total Count of AAs.
        '''
        #!# Yet to be modified #!#
        try:
            ### Setup ###
            self.dict['SufTree'] = {'=':total_count,'X':0}
            ### Build tree ###
            self.verbose(0,4,'Building suffice tree from markov extension probabilities...',0)
            typex = self.sufFromProb(self.suftree(),self.probtree)
            self.verbose(0,4,'Counting Xmer length occurences...',0)
            self._xmerCountPerLen(self.suftree(),0)
            self.verbose(0,1,'Done!',1)
            self.deBug(self.suftree())
            if not self.opt['Silent']: self.log.printLog('#INF','Suffix Tree populated using Markov extension probabilities for %d xmers.' % typex)
            for x in range(self.getMaxXmer(self.suftree())):
                if not self.opt['Silent']: self.log.printLog('#XMERS','%s counts of %s different %dmers.' % (rje.integerString(self.suftree()[x+1]['=']),rje.integerString(self.suftree()[x+1]['X']),(x+1)))
        except:
            self.log.errorLog('Error in _sufTreeFromProbTree',printerror=True,quitchoice=True)
            raise
#########################################################################################################################
    def _xmerCountPerLen(self,sufdic,xlen=0):  ### Populates self.suftree xmer length counts using sufdic
        '''
        Populates self.suftree() xmer length counts using sufdic.
        >> sufdic:suftree() level
        >> xlen:int = length of current Xmer (suftree() level)
        '''
        #!# Yet to be modified #!#
        if sufdic['='] == 0:     # No entries at this level
            return 
        xlen += 1
        if xlen not in self.suftree().keys():
            self.suftree()[xlen] = {'=':0,'X':0}
        for a in self.list['Alphabet']:
            if sufdic.has_key(a):
                self.suftree()[xlen]['X'] += 1
                self.suftree()[xlen]['='] += sufdic[a]['=']
                self._xmerCountPerLen(sufdic[a],xlen)
#########################################################################################################################
    def sufFromProb(self,sufdic,probdic,typex=0):  ### Populates sufdic using probdic
        '''
        Populates sufdic using probdic.
        >> sufdic:suftree level
        >> probdic:probtree level
        << returns typex = the number of xmer types processed
        '''
        #!# Yet to be modified #!#
        try:
            ### Setup ###
            if sufdic['='] == 0 or probdic.keys() == ['=']:     # No entries at this level
                return typex
            ### Generate exact counts and remainders ###
            count = 0   # Total count of generated Xmers
            exact = {}  # Dictionary of exact numbers
            rounded = {}    # Dictionary of integer occurrances
            remainders= {}  # Dictionary of AAs per given remainder
            for a in self.list['Alphabet']:
                if probdic.has_key(a):
                    exact[a] = probdic[a]['='] * sufdic['=']
                else:
                    exact[a] = 0.0
                rounded[a] = int(exact[a])
                remainder = exact[a] - rounded[a]
                if remainders.has_key(remainder):
                    remainders[remainder].append(a)
                else:
                    remainders[remainder] = [a]
                count += rounded[a]
            ### Sort remainders and add additional instances of chains ###
            while count < sufdic['=']:
                if remainders == {}:  # Something has gone wrong!
                    self.log.errorLog('Something has gone wrong! No remainders left but only %d Xmers where %d needed!' % (count,sufdic['=']),False,False)
                    raise ValueError
                elif remainders.keys() == ['0.0']:
                    self.log.errorLog('Something has gone wrong! Only 0.0 remainders left but only %d Xmers where %d needed!' % (count,sufdic['=']),False,False)
                    raise ValueError                    
                #self.deBug(remainders)
                best_rem = rje.sortKeys(remainders,revsort=True)[0]
                if len(remainders[best_rem]) == 1:
                    add_aa = remainders.pop(best_rem)[0]
                else:
                    ra = random.randint(1,len(remainders[best_rem]))
                    add_aa = remainders[best_rem].pop(ra-1)
                count += 1
                rounded[add_aa] += 1
                #self.deBug(remainders)
            ### Update sufdic
            for aa in self.list['Alphabet']:
                if rounded[aa] > 0:
                    sufdic[aa] = {'=':rounded[aa]}
                    typex += 1
                    if typex/50000 == typex/50000.0:
                        self.verbose(0,4,rje.integerString(typex),0)
                    elif typex/1000 ==typex/1000.0:
                        self.verbose(0,4,'.',0)
                    ## Extend further ##
                    if probdic.has_key(aa): # Longer chains to look at!
                        typex = self.sufFromProb(sufdic[aa],probdic[aa],typex)
            return typex
        except:
            self.deBug('SufDic: %s\n' % sufdic)
            self.deBug('ProbDic: %s\n' % probdic)
            self.deBug('Exact: %s\n' % exact)
            self.deBug('Rounded: %s\n' % rounded)
            self.log.errorLog('Error in _sufTreeFromProbTree',printerror=True,quitchoice=True)
            raise
#########################################################################################################################
    ### <8> ### Split file Markov generation methods                                                                    #
#########################################################################################################################
    def splitFileMarkov(self):  ### Generates markov chains for split file and recombines (saves memory) 
        '''Generates markov chains for split file and recombines (saves memory).'''
        try:
            #!# Add forking here for split and check functionality - not modified in V2 yet. #!#
            ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList']
            ### Deal with split ###
            seqlist.loadSeqs()
            if self.stat['Split'] < seqlist.seqNum():
                filelist = seqlist.splitSeq(self.stat['Split'])
            else:
                return self._buildSuffixTree(self.obj['SeqList'])
                
            ## Perform split ##
            resfiles = []
            for splitfile in filelist:
                self._buildSuffixTree(rje_seq.SeqList(log=self.log,cmd_list=['autoload=F']+self.cmd_list+['autofilter=F','seqin=%s' % splitfile]),partial=True)
                resfiles.append(self.markovAnalysis(partial=True))
            self.dict['SufTree'] = {}
            ## Combine ##
            if not self.opt['Silent']: self.log.printLog('#INF','Combining <= %dmer details from %s files.' % (max_xmer,rje.integerString(len(filelist))))
            SPLITFILES = [] # List of file handles
            lastline = {}   # Dictionary of last line read from SPLITFILES
            aacount = 0
            seqnum = 0
            # Setup Input files #
            for splitfile in resfiles:
                SPLITFILES.append(open(splitfile,'r'))
                aacount += string.atoi(rje.matchExp('^aacount=(\d+)',SPLITFILES[-1].readline())[0])
                seqnum += string.atoi(rje.matchExp('^seqnum=(\d+)',SPLITFILES[-1].readline())[0])
                SPLITFILES[-1].readline()   # Remaining header line
                lastline[SPLITFILES[-1]] = SPLITFILES[-1].readline()    # First xmer line
            # Setup Output file #
            delimit = rje.getDelimit(self.cmd_list)
            filebase = '%s.markov' % (self.obj['SeqList'].info['Basefile'])
            for x in range(min_xmer,(max_xmer+1)):
                filename = '%s.%dmer.%s' % (filebase,x,rje.delimitExt(delimit))
                MARKOV = open(filename,'w')
                outlist = ['xmer','obs','f_obs','f_exp']
                MARKOV.write('%s\n' % string.join(outlist,delimit))
                MARKOV.close()
            # Work through possible Xmers #
            types = self._combineSplits(filebase,'',lastline,aacount,seqnum,delimit)
            if not self.opt['Silent']: self.verbose(0,2,'Done!',2)
            if not self.opt['Silent']: self.log.printLog('#INF','Data Combined for <= %dmer from %s files.' % (max_xmer,rje.integerString(len(lastline))))

            ## Clean up ##
            for SPLITFILE in SPLITFILES:
                SPLITFILE.close()
            for splitfile in resfiles:  #!# Why were files missing?!
                os.unlink(splitfile)
            if not self.opt['Silent']: self.log.printLog('#INF','%s partial results files deleted.' % (rje.integerString(len(resfiles))))
            if rje.yesNo('Delete split fasta files?'):
                for splitfile in filelist:
                    os.unlink(splitfile)
                if not self.opt['Silent']: self.log.printLog('#INF','%s partial fasta files deleted.' % (rje.integerString(len(filelist))))

        except:
            self.errorLog('Error in Markov.run()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _combineSplits(self,filebase,xmer,lastline,aacount,seqnum,delimit,typex=0):    ### Combines and outputs observed vs expected 
        '''
        Combines and outputs observed vs expected to file.
        >> filebase:str = output file name base - add .xmer.ext
        >> xmer:str = Xmer for output
        >> lastline:dictionary of the last line read in from each split file
        >> aacount:int = Total aa count of dataset
        >> seqnum:int = Number of sequences in dataset
        >> delimit:str = Text delimiter for output file
        << typex:int = count of different Xmer types
        '''
        if len(xmer) > self.stat['MaxXmer']:
            return typex
        for aa in self.list['Alphabet']:
            ## This xmer ##
            _xmer = xmer + aa
            (pref,suff,comm) = (_xmer[:-1],_xmer[1:],_xmer[1:-1])
            num = {}
            ## Counts ##
            for _str in [_xmer,pref,suff,comm]:
                num[_str] = 0
            ## Split Files ##
            for SPLIT in lastline.keys():
                line = lastline[SPLIT].split(delimit)
                if line[0] == _xmer:    # Found it!
                    num[_xmer] += string.atoi(line[1])
                    num[pref] += string.atoi(line[2])
                    num[suff] += string.atoi(line[3])
                    num[comm] += string.atoi(line[4])
                    lastline[SPLIT] = SPLIT.readline()  # Advance to next line
            ## Counts ##
            if num[_xmer] == 0:     # This xmer and all descendants absent from dataset
                continue
            ## Output ##
            typex += 1  # Add this xmer to type count
            filename = '%s.%dmer.%s' % (filebase,len(_xmer),rje.delimitExt(delimit))
            self._markovOutput(filename,_xmer,num,aacount,seqnum,delimit)
            if typex/500000 == typex/500000.0:
                self.verbose(0,4,rje.integerString(typex),0)
            elif typex/10000 == typex/10000.0:
                self.verbose(0,4,'.',0)
            ## Extend further ##
            typex = self._combineSplits(filebase,_xmer,lastline,aacount,seqnum,delimit,typex)            
        return typex
#########################################################################################################################
    ### <8> ### Split file Markov generation methods                                                                    #
#########################################################################################################################
    def negVPos(self):  ### Compares Negative versus Positive sets and introduces stats too                         !V2.1
        '''Compares Negative versus Positive sets and introduces stats too.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['basefile=%s' % self.info['NegVPos']])
            if rje.exists('%s.markov.tdt' % self.info['NegVPos']) and not self.opt['Force']:
                db = self.db().addTable('%s.markov.tdt' % self.info['NegVPos'],['xmer','type'],name='markov')
                db.dataFormat({'obs':'int','f_obs':'num','f_exp':'num','obs/exp':'num','p_obs':'num','sig':'num','x':'int'})
            else:
                db = None
                for t in ['negatives','positives']:
                    base = '%s.%s' % (self.info['NegVPos'],t)
                    for x in range(self.stat['MinXmer'],self.stat['MaxXmer']+1):
                        newdb = self.db().addTable('%s.markov.%dmer.csv' % (base,x),['xmer'],name='%s%d' % (t[:3],x))
                        newdb.addField('type','xmer',evalue=t[:3])
                        newdb.addField('x','type',evalue=x)
                        newdb.list['Keys'] = ['xmer','type']
                        newdb.remakeKeys()
                        if db: self.db().mergeTables(db,newdb,overwrite=True,matchfields=True)
                        else: db = newdb
            ### ~ [1] ~ Calculate basic stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                db.info['Name'] = 'markov'; db.info['Delimit'] = '\t'; db.remakeKeys()
                db.dataFormat({'obs':'int','f_obs':'num','f_exp':'num'})
                db.addField('obs/exp'); db.addField('p_obs'); db.addField('sig')
                ex = 0.0; etot = db.entryNum()
                for entry in db.entries():
                    self.progLog('\r#PROB','Calculating enrichment probabilities (%.2f%%)' % (ex/etot)); ex += 100.0
                    entry['obs/exp'] = entry['f_obs'] / entry['f_exp']
                    if entry['obs/exp'] < 1.0: entry['p_obs'] = 1.0 - rje.logPoisson(entry['obs']+1,entry['obs']/entry['obs/exp'],exact=False,callobj=self)
                    else: entry['p_obs'] = rje.logPoisson(entry['obs'],entry['obs']/entry['obs/exp'],exact=False,callobj=self)
                    entry['sig'] = rje.logBinomial(1,20**entry['x'],entry['p_obs'],callobj=self)
                    self.deBug(entry)
                self.printLog('\r#PROB','Calculating enrichment probabilities (%.2f%%)' % (ex/etot))
                db.saveToFile()
            ### ~ [2] ~ Do complicated calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db.info['Name'] = 'negvpos'
            db.reshapeWide('type',['obs','f_obs','f_exp','obs/exp','p_obs','sig'])
            db.fillBlanks(blank=0.0,fields=['f_obs|neg','f_exp|neg','obs/exp|neg','p_obs|neg'],fillempty=True); db.fillBlanks(0,['obs|neg'],True); db.fillBlanks(1.0,['sig|neg'],True)
            db.fillBlanks(blank=0.0,fields=['f_obs|pos','f_exp|pos','obs/exp|pos','p_obs|pos'],fillempty=True); db.fillBlanks(0,['obs|pos'],True); db.fillBlanks(1.0,['sig|pos'],True)
            db.addField('obs_neg/pos'); db.addField('p_neg/pos'); db.addField('obs/exp_neg/pos'); db.addField('enr'); db.addField('p_enr')
            ex = 0.0; etot = db.entryNum()
            for entry in db.entries():
                self.progLog('\r#PROB','Calculating enrichment probabilities (%.2f%%)' % (ex/etot)); ex += 100.0
                try: entry['obs_neg/pos'] = entry['f_obs|neg'] / entry['f_obs|pos']
                except: entry['obs_neg/pos'] = entry['obs|neg']
                if entry['obs_neg/pos'] < 1.0:
                    if entry['obs_neg/pos']:
                        entry['p_neg/pos'] = 1.0 - rje.logPoisson(entry['obs|neg']+1,entry['obs|neg']/entry['obs_neg/pos'],exact=False,callobj=self)
                    else: entry['p_neg/pos'] = 1.0
                else:
                    if entry['obs_neg/pos']: entry['p_neg/pos'] = rje.logPoisson(entry['obs|neg'],entry['obs|neg']/entry['obs_neg/pos'],exact=False,callobj=self)
                    else: entry['p_neg/pos'] = 0.0
                try: entry['p_neg/pos'] = rje.logBinomial(1,20**entry['x'],entry['p_neg/pos'],callobj=self)
                except: self.deBug(entry)
                try: entry['obs/exp_neg/pos'] = entry['obs/exp|neg'] / entry['obs/exp|pos']
                except: entry['obs/exp_neg/pos'] = 0.0
                if entry['obs/exp|neg'] == 1.0 or entry['obs/exp|pos'] == 1.0: entry['enr'] = 'none'
                elif entry['obs/exp|neg'] > 1.0 and entry['obs/exp|pos'] > 1.0: entry['enr'] = 'up'
                elif entry['obs/exp|neg'] < 1.0 and entry['obs/exp|pos'] < 1.0: entry['enr'] = 'down'
                elif entry['obs/exp|neg'] > 1.0 and entry['obs/exp|pos'] < 1.0: entry['enr'] = 'neg'
                elif entry['obs/exp|neg'] < 1.0 and entry['obs/exp|pos'] > 1.0: entry['enr'] = 'pos'
                entry['p_enr'] = entry['sig|neg'] * entry['sig|pos']
            self.printLog('\r#PROB','Calculating enrichment probabilities (%.2f%%)' % (ex/etot))
            db.saveToFile()
        except: self.errorLog(rje_zen.Zen().wisdom())            
#########################################################################################################################
### End of SECTION II: Markov Class                                                                                     #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def xmersFromYmers(callobj=None,ymers=[],alphabet=None):     ### Makes all possible Xmers from shorter Ymers
    '''
    Makes all possible Xmers from shorter Ymers.
    >> ymers:list of str = possible Ymers from which to construct Xmers (Should be length x-1)
    >> alphabet:list of str = alphabet from which to append ymers to make xmers
    << xmers:list of str = return Xmers
    '''
    try:
        if not alphabet: return []
        xmers = []
        for ymer in ymers:
            for a in alphabet:
                if '%s%s' % (ymer[1:],a) in ymers:  # End of new Xmer allowed
                    xmers.append('%s%s' % (ymer,a))
        return xmers
    except:
        if callobj: callobj.errorLog('Error in rje_markov.xmersFromYmers()',printerror=True,quitchoice=True)
    return []
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
    try: Markov(mainlog,cmd_list).run()

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