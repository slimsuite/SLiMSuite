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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 29 Kingsland Parade, Portobello, Dublin 8, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Module:       rje_seqgen
Description:  Random Sequence Generator Module
Version:      1.7
Last Edit:    17/01/13
Copyright (C) 2006  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to generate a number a random sequences based on input AA or Xmer frequencies and the desired
    order of markov chain from which to draw the amino acid (or nucleotide) probabilites.

    If poolgen=T, then the amino acid frequencies will be used to make a finite pool of amino acids from which sequences
    will be built. This will ensure that the total dataset has the correct amino acid frequencies. Because of the
    potential of this to get 'stuck' in impossible sequence space - especially if screenx > 0 - an additional parameter
    poolcyc=X determines how many times to retry the generation of sequences. If poolgen=F, then generation of sequences
    will be faster but the resulting dataset may have Xmer frequencies that differ greatly from the input frequencies,
    depending on how many (and which) redundant and/or screened Xmer-containing sequences are removed. (If seqin contains
    only one peptide then each random peptide will be a scramble of that peptide.)

    !!!NEW!!! Verson 1.3 has new scramble function, which takes in a list of peptides and tries to construct scrambled
    versions of them. In this case, screenx=X sets the length of common Xmers between the scrambled peptide and the
    original peptide at which a scrambled peptide will be rejected. This should set > 1, else all peptides will be
    rejected. (If left at the default of zero, no peptides will be rejected.) In this mode, outfile=FILE will set the
    name of a delimited output file containing two columns: peptide & scramble. (Default filename = scramble.tdt)

    !!!NEW!!! Version 1.5 has a new BLAST-centred method for making a random dataset from an input dataset, retaining the
    approximate evolutionary relationships as defined by BLAST homology, which should result in similar GABLAM statistics
    for the randomised dataset. For this, a random sequence is created first. Any BLAST hits between this and other
    sequences are then mapped, keeping the required percentage identity (and using different amino acids drawn from the
    frequency pool for the rest). The next sequence is taken, completed and then the same process followed, until all
    sequences have been made. Improvements to make: (a) incorporate similarity too; (b) adjust aa frequencies after BLAST
    mapping. This method is activated by the blastgen=T option and has limited options as yet.
    NB. The input dataset will *not* be subject to rje_seq filtering.

    !!!NEW!!! Version 1.6 has an EST randomiser. This will go through each sequence in turn and generate a new sequence of
    the same length using the NT frequencies (or markov chain frequencies) of just that sequence. Updated in V1.7 to make
    this work for proteins too.

Commandline:
    ## Generation options ##
    seqnum=X        : Number of random sequences to generate [24]
    seqlen=X,Y      : Range of lengths for random sequences [10]
    markovx=X       : Order of markov chain to use for sequence construction [1]
    aafreq=FILE     : File from which to read AA Freqs [None]
    xmerfile=FILE   : File from which to read Xmer frequencies for sequence generation [None]
    xmerseq=FILE    : Sequence file from which to calculate Xmer frequencies [None]
    * xmerseq is overwridden by xmerfile and aafreq. aafreq only works if markovx=1 and is over-ridden by xmerfile *
    nrgen=T/F       : Whether to generate a non-redundant sequence list (whole-sequence redundancies only) [True]
    poolgen=T/F     : Whether to build sequences using a fixed AA pool (exact freqs) or probabilities only [False]
    poolcyc=X       : Number of times to retry making sequences if rules are broken [1]
    maxhyd=X        : Maximum mean hydrophobicity score [10]

    ## Output & Naming ##
    outfile=FILE: Output file name [randseq.fas]
    randname=X  : Name 'leader' for output fasta file [randseq]
    randdesc=T/F: Whether to include construction details in description line of output file [True]
    idmin=X     : Starting numerical ID for randseq (allows appending) [1]
    idmax=X     : Max number for randseq ID. If < seqnum, will use seqnum. If <0, no zero-prefixing of IDs. [0]
    append=T/F  : Whether to append to outfile [False]

    ## Other Xmers of Interest ##
    screenfile=FILE : File of Xmers to screen in generated sequences [None]
    xmerocc=T/F     : Whether to output occurrences of screened Xmers [True]
    screenx=X       : Reject generated sequences containing screened Xmers >= X [0]
    screenrev=T/F   : Whether to screen reverse Xmers too [False]

    ## Peptide Scrambling Parameters ##
    # Uses seqnum=X, randdesc, idmin/max=X and randname=X but for each input peptide (oldname_randnameID)
    scramble=T/F    : Run peptide scrambler [False]
    fullscramble=T/F: Generate all possible scrambles for each peptide in TDT [False]
    scramblecyc=X   : Number of attempts to try each scramble before giving up [10000]
    seqin=FILE      : Sequence file containing peptides to scramble [None]
    peptides=LIST   : Alternative peptide sequence input for scrambling []
    outfile=FILE    : Output delimited file of scrambled peptides or peptide and scrambled sequence. [scramble.tdt]
    teiresias=X		    : Length of patterns to be screened by additional TEIRESIAS search on scrambled vs original [0]
    teiresiaspath=PATH	: Path to TEIRESIAS ['c:/bioware/Teiresias/teiresias_char.exe'] * Use forward slashes (/)

    ### BLAST-based dataset randomiser (uses some of the Output options listed) ###
    blastgen=T/F    : Activate the BLASTGen method [False]
    seqin=FILE      : Input sequence file to randomise [None]
    keepnames=T/F   : Whether to keep same input names in outfile [False]

    ### EST Randomiser ###
    estgen=T/F      : Whether to run EST randomiser method [False]
    
Uses general modules: copy, os, re, string, sys, time
Uses RJE modules: rje, rje_markov, rje_seq, rje_blast
Other modules needed: rje_dismatrix, rje_pam, rje_sequence, rje_uniprot
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, random, re, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje, rje_markov, rje_seq, rje_sequence, rje_zen
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Initial Working version
    # 1.1 - Added max hydrophobicity
    # 1.3 - Added peptide scrambler
    # 1.4 - Separated Xmer screen and Teiresias pattern screen
    # 1.5 - Added BLASTGen Method
    # 1.6 - Checked function with DNA. Added EST randomiser function.
    # 1.7 - Modified/fixed ESTgen function to work for protein sequences.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Fix counting of xmer types in read
    # [ ] : Allow poolgen=T when markovx > 1
    # [Y] : Append=T option
    # [ ] : Additional xmer screeing (e.g. chemistry synthesis rules)
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyright) = ('RJE_SEQGEN', '1.7', 'January 2013', '2006')
    description = 'Random Sequence Generator'
    author = 'Dr Richard J. Edwards.'
    comments = [rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyright,comments)
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
### SECTION II: SeqGen Class                                                                                            #
#########################################################################################################################
class SeqGen(rje.RJE_Object):     
    '''
    Sequence Generator Class. Author: Rich Edwards (2005).

    Info:str
    - AAFreqFile = File from which to read AA Freqs [None]
    - XmerFile = File from which to read Xmer frequencies for sequence generation [None]
    - XmerSeq = Sequence file from which to generate build markov chains [None]
    - OutFile = Output file name [randseq.fas]
    - RandName = Name 'leader' for output fasta file [randseq]
    - ScreenFile = File of Xmers to screen in generated sequences [None]
    - SeqIn = Sequence file containing peptides to scramble [None]
    - TeiresiasPath = Path to TEIRESIAS ['c:/bioware/Teiresias/teiresias_char.exe'] * Use forward slashes (/)

    Opt:boolean
    - NRGen = Whether to generate a non-redundant sequence list [True]
    - PoolGen = Whether to build sequences using a fixed AA pool (exact freqs) or probabilities only [False]
    - RandDesc = Whether to include construction details in description line of output file [True]
    - XmerOcc = Whether to output occurrences of screened Xmers [True]
    - ScreenRev = Whether to screen reverse Xmers too [False]
    - Scramble = Run peptide scrambler [False]
    - BLASTGen = Activate the BLASTGen method [False]
    - KeepNames = Whether to keep same names in outfile [False]
    - FullScramble = Generate all possible scrambles for each peptide in TDT [False]

    Stat:numeric
    - SeqNum = Number of random sequences to generate [24]
    - MinSeqLen & MaxSeqLen = Range of lengths for random sequences [10]
    - MarkovX = Order of markov chain to use for sequence construction [1]
    - PoolCyc = Number of times to retry making sequences if rules are broken [1]
    - IDMin = Starting numerical ID for randseq (allows appending) [1]
    - IDMax = Max number for randseq ID. If < seqnum, will use seqnum. If <0, no zero-prefixing of IDs. [0]
    - ScreenX = Reject generated sequences containing screened Xmers >= X [0]
    - MaxHyd = Maximum mean hydrophobicity score [10]
    - Teiresias = Length of patterns to be screened by additional TEIRESIAS search on scrambled peptides vs original peptide [0]
    - ScrambleCyc = Number of attempts to try each scramble before giving up [10000]

    List:list
    - Peptides = Alternative peptide sequence input for scrambling []

    Dict:dictionary

    Obj:RJE_Objects
    - BuildMarkov = rje_markov.Markov : Markov chain probabilities from which to make sequences
    - ScreenMarkov = rje_markov.Markov : Xmers to screen from sequences
    - SeqList = rje_seq.SeqList : random sequences
    '''
    ### Attributes
    lenlist = []    # List of sequence lengths, where first entry is total combined length
    debuglenlist = []
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AAFreqFile','XmerFile','XmerSeq','OutFile','RandName','ScreenFile','SeqIn','TeiresiasPath']
        self.optlist = ['NRGen','PoolGen','RandDesc','XmerOcc','ScreenRev','Scramble','FullScramble','Teiresias',
                        'BLASTGen','KeepNames','ESTGen']
        self.statlist = ['SeqNum','MinSeqLen','MaxSeqLen','MarkovX','PoolCyc','IDMin','IDMax','ScreenX','MaxHyd',
                         'Teiresias','ScrambleCyc']
        self.listlist = ['Peptides']
        self.dictlist = []
        self.objlist = ['BuildMarkov','ScreenMarkov','SeqList']
        ### Defaults ###
        self._setDefaults(info='None',opt=True,stat=0.0,obj=None,setlist=True,setdict=True)
        self.setInfo({'OutFile':'randseq.fas','RandName':'randseq'})
        self.setStat({'SeqNum':24,'MinSeqLen':10,'MaxSeqLen':10,'MarkovX':1,'PoolCyc':1,'IDMin':1,'ScreenX':0,
                      'MaxHyd':10,'ScrambleCyc':10000})
        self.setOpt({'PoolGen':False,'ScreenRev':False,'Scramble':False,'BLASTGen':False,'KeepNames':False,
                     'FullScramble':False,'ESTGen':False})
        ### Other Attributes ###
        self.lenlist = []
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
                self._cmdRead(cmd,type='info',att='AAFreqFile',arg='aafreq')
                self._cmdReadList(cmd,'info',['XmerFile','XmerSeq','OutFile','RandName','ScreenFile'])
                self._cmdRead(cmd,type='stat',att='MaxHyd')
                self._cmdReadList(cmd,'opt',['NRGen','PoolGen','RandDesc','XmerOcc','ScreenRev','ESTGen'])
                self._cmdRead(cmd,type='min',att='MinSeqLen',arg='seqlen')
                self._cmdRead(cmd,type='max',att='MaxSeqLen',arg='seqlen')
                self._cmdReadList(cmd,'int',['SeqNum','MarkovX','PoolCyc','IDMin','IDMax','ScreenX'])
                ## Peptide Scramble Options ##
                self._cmdReadList(cmd,'opt',['Scramble','FullScramble'])
                self._cmdRead(cmd,type='list',att='Peptides')
                self._cmdRead(cmd,type='info',att='SeqIn')
                self._cmdReadList(cmd,'int',['Teiresias','ScrambleCyc'])
                self._cmdRead(cmd,type='fullpath',att='TeiresiasPath')
                ## BLASTGen Options ##
                self._cmdReadList(cmd,'opt',['BLASTGen','KeepNames'])
            except:
                self.log.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Methods                                                                                      #
#########################################################################################################################
    def seqGen(self,setobj=True,save=True):  ### Master run method for random sequence generation
        '''
        Master run method for random sequence generation.
        >> setobj:bool = whether to clear and reset all objects
        >> save:bool = whether to save sequences to file
        << returns rje_seq.SeqList object with sequences
        '''
        try:
            ### <0> ### Setup
            _stage = '<0> Setup'
            ## General variables ##
            if self.stat['SeqNum'] < 1:
                self.log.errorLog('Trying to generate %d sequences! Reset to 1.' % self.stat['SeqNum'],False,False)
                self.stat['SeqNum'] = 1
            if self.stat['MarkovX'] < 1:
                self.log.errorLog('Trying to use %d Order Markov Chain! Reset to 1st order.' % self.stat['MarkovX'],False,False)
                self.stat['MarkovX'] = 1
            if self.opt['PoolGen'] and self.stat['MarkovX'] > 1:
                self.log.errorLog('Cannot generate sequences using AA Pool when MarkovX > 1.',False,False)
                if rje.yesNo('Use raw frequencies rather than AA pool?'):
                    self.opt['PoolGen'] = False
                else:
                    self.stat['MarkovX'] = 1
                    self.log.printLog('#CMD','MarkovX reduced to 1.')
            if self.stat['IDMax'] >= 0 and self.stat['IDMax'] < self.stat['IDMin']: self.stat['IDMax'] = self.stat['SeqNum']
            if self.stat['MaxSeqLen'] < self.stat['MinSeqLen']:
                self.log.errorLog('MaxSeqLen < MinSeqLen! Increased to MinSeqLen (%daa).' % self.stat['MinSeqLen'],False,False)
                self.stat['MaxSeqLen'] = self.stat['MinSeqLen']
            if self.stat['PoolCyc'] < 0: self.stat['PoolCyc'] = 0

            ## Setup objects ##
            # Markov Build Object #
            markovbuild = self.obj['BuildMarkov']
            if setobj or not markovbuild:
                markovbuild = rje_markov.Markov(log=self.log,cmd_list=self.cmd_list+['seqin=None'])
                self.obj['BuildMarkov'] = markovbuild
                #markovbuild.stripToX(markovbuild.suftree,self.stat['MarkovX'])
            if not markovbuild.probtree:    ### Need to generate markov probabilities
                if not markovbuild.suftree(): ### Need to generate markov counts
                    if self.info['XmerFile'] != 'None':
                        markovbuild.readXmers(xfile=self.info['XmerFile'],clear=True)
                        self.info['Name'] = self.info['XmerFile']
                    elif self.info['AAFreqFile'] != 'None' and self.stat['MarkovX'] == 1:
                        markovbuild.info['AAFreqFile'] = self.info['AAFreqFile']
                        self.info['Name'] = self.info['AAFreqFile']
                        markovbuild.aafreq = markovbuild.obj['SeqList'].aaFreq(loadfile=self.info['AAFreqFile'])
                        markovbuild._aaFreqToSufTree()
                    elif self.info['XmerSeq']:
                        markovbuild.obj['SeqList'].info['Name'] = self.info['XmerSeq']
                        self.info['Name'] = self.info['XmerSeq']
                        markovbuild._buildSuffixTree(markovbuild.obj['SeqList'])
                markovbuild.stripToX(markovbuild.suftree(),self.stat['MarkovX'])
                #?!#markovbuild._buildProbTree()
            self.deBug(markovbuild.suftree())

            print(' '.join(['Xmers', self.stat['MarkovX'], 'vs loaded', markovbuild.getMaxXmer(markovbuild.suftree())]))
            self.deBug(markovbuild.suftree())
            # Markov Screen Object #
            markovscreen = self.obj['ScreenMarkov']
            if setobj or not markovscreen:
                markovscreen = rje_markov.Markov(log=self.log,cmd_list=self.cmd_list+['seqin=None'])
                self.obj['ScreenMarkov'] = markovscreen
            if self.info['ScreenFile'] != 'None' and not markovscreen.suftree(): ### Need to generate markov counts
                markovscreen.readXmers(xfile=self.info['ScreenFile'],clear=True)
            # SeqList object #
            if setobj or not self.obj['SeqList']:
                self.obj['SeqList'] = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None'])
            seqlist = self.obj['SeqList']
            seqlist.info['Name'] = self.info['OutFile']
            
            ### Build List of sequence lengths ###
            _stage = 'LenList'
            self.lenlist = [0]
            if self.stat['MaxSeqLen'] == self.stat['MinSeqLen']:
                self.lenlist = [self.stat['SeqNum'] * self.stat['MinSeqLen']] + [self.stat['MinSeqLen']] * self.stat['SeqNum']
            else:
                for s in range(self.stat['SeqNum']):
                    rlen = random.randint(self.stat['MinSeqLen'],self.stat['MaxSeqLen'])
                    self.lenlist[0] += rlen
                    self.lenlist.append(rlen)

            ### Generate Sequences ###
            self.log.printLog('#SEQ','Generating %s seqs, Total Len = %s.' % (rje.integerString(self.stat['SeqNum']),rje.integerString(self.lenlist[0])))
            g = 0
            while g < self.stat['PoolCyc']:
                g += 1
                if self._genSeqs(retry=self.stat['PoolCyc']):    # Sequences made OK
                    break
                if g < self.stat['PoolCyc']:
                    self.log.printLog('#CYC','Attempt %d of %d to generate sequences failed. Trying Again.' % (g,self.stat['PoolCyc']))
                else:
                    self.log.printLog('#CYC','Attempt %d of %d to generate sequences failed. Exiting.' % (g,self.stat['PoolCyc']))
                    return
            self.log.printLog('#SEQ','%s seqs generated by Attempt %d.' % (rje.integerString(seqlist.seqNum()),g))
            seqlist.saveFasta(append=self.opt['Append'])

            ### Screen for Xmers ###
            if self.info['ScreenFile'] != 'None' and self.opt['XmerOcc']:
                ## Setup file ##
                maxmer = int(markovscreen.stat['MaxXmer'])
                delimit = rje.getDelimit(self.cmd_list)
                occfile = '%s.%s.xmer_occ.%s' % (rje.baseFile(seqlist.info['Name']),rje.baseFile(self.info['ScreenFile']),rje.delimitExt(delimit))
                XOCC = open(occfile,'w')
                headlist = ['randseq']
                for x in range(maxmer):
                    headlist.append('%dmer' % (x+1))
                rje.writeDelimit(XOCC,headlist,delimit)
                ## Calculations ##
                for seq in seqlist.seq:
                    xocc = {}
                    for x in range(maxmer):
                        xocc[x+1] = 0
                    outlist = [seq.shortName()] 
                    for r in range(seq.seqLen()):
                        subseq = seq.info['Sequence'][r:(r+maxmer)]
                        while subseq != '':
                            if markovscreen._xmerCountFromSufTree(subseq) > 0:
                                xocc[len(subseq)] += 1
                            subseq = subseq[:-1]
                    for x in range(maxmer):
                        outlist.append(rje.integerString(xocc[x+1]))
                    rje.writeDelimit(XOCC,outlist,delimit)
                XOCC.close()
                self.log.printLog('#OCC','%s <=%dmer occurrence for %s seqs output to %s.' % (markovscreen.info['Name'],maxmer,rje.integerString(seqlist.seqNum()),occfile))
                
        except:
            self.log.errorLog('Error in seqGen (%s)' % _stage,printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _genSeqs(self,retry=1): ### Builds a set of seqences using self.attributes
        '''
        Builds a set of seqences using self.attributes.
        >> retry:int = The number of times to try generating sequences.
        << True/False depending on success
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            lenlist = self.lenlist[0:]
            self.debuglenlist = lenlist
            #self.deBug(lenlist)
            buildmarkov = self.obj['BuildMarkov']
            if self.opt['PoolGen']:
                buildmarkov._sufTreeFromProbTree(lenlist[0])
                self.log.printLog('#GEN','Suffix Tree "Pool" generated for dataset, Total Len = %s.' % rje.integerString(lenlist[0]))
            seqlist = self.obj['SeqList']
            seqlist.seq = []
            seqcyc = self.stat['PoolCyc']
            stuckseq = 0    # Number of last sequence to get 'stuck' during generation
            nrlist = []     # List full sequences (or crc64 values) in NR mode
            seqx = 0

            ### Sequence Creation ###
            _stage = 'Sequence Creation'
            while len(lenlist) > 1: # Still sequences to make
                seqx = len(nrlist) + 1
                seqlen = lenlist.pop(-1)
                #self.verbose(0,3,'Generating sequence %d (%d aa)...' % (seqx,seqlen),0)
                self.log.printLog('#GEN','Generating sequence %d (%d aa).' % (seqx,seqlen))
                seqname = self._randSeqName(seqx,seqlen)
                sequence = self._makeSeq(seqlen)
                ## Check for dodgy sequence ##
                _stage = 'Checking Dodgy Sequence'
                seqstate = 'ok'
                if len(sequence) < seqlen:
                    seqstate = 'bad'
                elif (rje_sequence.eisenbergHydropathy(sequence) / len(sequence)) > self.stat['MaxHyd']:
                    seqstate = 'hydrophobic'
                #print sequence, rje_sequence.eisenbergHydropathy(sequence)
                if self.opt['NRGen'] and sequence in nrlist:
                    seqstate = 'redundant'
                if seqstate != 'ok':
                    ## Put sequence 'back' ##
                    lenlist.append(seqlen)
                    self._addToBuildMarkov(sequence)
                    ## Loop again for this sequence? ##
                    if seqcyc > 0:  
                        self.log.printLog('#ERR','%s = %s sequence. Regenerating.' % (sequence,seqstate))
                        seqcyc -= 1
                        continue
                    ## Try losing several sequences and trying again? ##
                    elif seqx > stuckseq and self.stat['PoolCyc'] > 1 and seqx > self.stat['PoolCyc']:
                        self.log.printLog('#ERR','Stuck generating sequence %d. Regenerating last %d sequences.' % (seqx,self.stat['PoolCyc']))
                        stuckseq = seqx
                        nrlist = nrlist[:-self.stat['PoolCyc']]
                        for x in range(self.stat['PoolCyc']):
                            seq = seqlist.seq.pop(-1)
                            lenlist.append(seq.seqLen())
                            self._addToBuildMarkov(seq.info['Sequence'])
                        continue
                    ## Clear all sequences and try again! ##
                    else:
                        if self.opt['PoolGen']:
                            for seq in seqlist.seq:
                                self._addToBuildMarkov(seq.info['Sequence'])
                        if retry > 1 or (self.stat['Interactive'] > 0  and rje.yesNo('Failed to generate sequence %d after all attempts! Try again?' % seqx)):
                            self.log.printLog('#ERR','Really stuck generating sequence %d. Try again!' % seqx)
                            return self._genSeqs(retry=retry-1)
                        else:
                            self.log.printLog('#ERR','Failed to generate sequence %d after all attempts!' % seqx)
                            return False
                ## Sequence OK! ##
                if self.opt['NRGen']:
                    nrlist.append(sequence)
                else:
                    nrlist.append(seqx)
                seqlist._addSeq(seqname,sequence)
                seqcyc = self.stat['PoolCyc']
            ### Success! ###
            self.log.printLog('#GEN','Sequence Generation Successful! (%s seqs, Total Len = %s)' % (rje.integerString(seqx),rje.integerString(lenlist[0])))
            return True
        except:
            self.log.errorLog('Error in _genSeqs (%s)' % _stage,printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def _makeSeq(self,seqlen):    ### Generates and returns random sequence using self.attributes
        '''
        Generates and returns random sequence using self.attributes.
        >> seqlen:int = length of sequence to make
        << sequence:str = generated sequence
        '''
        try:
            ### Setup ###
            self._debugTree(seqlen,'')
            sequence = ''
            chopcyc = self.stat['PoolCyc']
            buildmarkov = self.obj['BuildMarkov']
            maxcyc = self.stat['PoolCyc'] * seqlen
            ### Build Sequence ###
            while len(sequence) < seqlen:   # More AAs needed
                ## Extend ##
                badseq = False  # Whether something wrong with sequence extension
                extseq = self._addAA(sequence)
                self._debugTree(seqlen - len(extseq),extseq)
                ## Check extension ##
                if extseq == sequence: # Something wrong in _addAA
                    badseq = True
                elif len(sequence) >= self.stat['ScreenX'] and self.stat['ScreenX'] > 0: # Check for bad xmers
                    last_xmer = extseq[-self.stat['ScreenX']:]
                    if self.obj['ScreenMarkov']._xmerCountFromSufTree(last_xmer) > 0:
                        badseq = True
                    elif self.opt['ScreenRev'] and self.obj['ScreenMarkov']._xmerCountFromSufTree(rje.strReverse(last_xmer)) > 0:
                        badseq = True
                ## Retry if something wrong ##
                if badseq:
                    if chopcyc > 1 and maxcyc > 0: # Try again
                        self._addToBuildMarkov(extseq)
                        sequence = extseq[:-1]
                        if self.opt['PoolGen']:
                            buildmarkov._removeFromSufTree(sequence,self.stat['MarkovX'])
                        chopcyc -= 1
                        maxcyc -= 1
                        continue
                    else:   # Give up!
                        return extseq
                ## Continue OK ##
                chopcyc = self.stat['PoolCyc']
                sequence = extseq[0:]
            return sequence
        except:
            self.log.errorLog('Error in _makeSeq',quitchoice=True)
            return sequence
#########################################################################################################################
    def _addAA(self,sequence):  ### Returns extended sequence using self.attributes
        '''
        Returns extended sequence using self.attributes.
        >> sequence:str = input sequence
        << sequence:str = extended sequence
        '''
        try:
            ### Setup ###
            buildmarkov = self.obj['BuildMarkov']
            sufdic = buildmarkov.probtree   # Default = use probabilities
            if self.opt['PoolGen']:
                sufdic = buildmarkov.suftree()
            if self.stat['MarkovX'] > 1:
                pref = sequence[(1-self.stat['MarkovX']):]
                while pref:
                    sufdic = sufdic[pref[0]]
                    pref = pref[1:]
            ### Get random number ###
            if self.opt['PoolGen']:
                if sufdic['='] < 1:
                    self.log.errorLog('Problem in _addAA to %s. Sufdic[=] = %d but need to add AA!' % (sequence,sufdic['=']),True,False)
                    return sequence
                p = random.randint(1,sufdic['='])   # Use <=
            else:
                p = random.random()
            ### Get next letter ###
            for aa in buildmarkov.alphabet:
                #x#print p, aa
                if sufdic.has_key(aa) and sufdic[aa]['='] > 0:
                    #x#print p, aa, 'vs', sufdic[aa]['=']
                    if p <= sufdic[aa]['=']:  # Found our letter!
                        extseq = '%s%s' % (sequence,aa)
                        if self.opt['PoolGen'] and len(extseq) >= self.stat['MarkovX']:
                            extmer = extseq[-self.stat['MarkovX']:]
                            #x#print extseq, 'removing', extmer
                            buildmarkov._removeXmerFromSufTree(extmer)
                        return extseq
                    else:
                        p -= sufdic[aa]['=']
            ### Something wrong! No letter added. ###
            self.log.errorLog('Error in _addAA (%s): Could not add AA (p=%f). Sufdic: %s' % (sequence,p,sufdic),False,False)
            return sequence            
        except:
            self.log.errorLog('Error in _addAA')
            return sequence
#########################################################################################################################
    def _addToBuildMarkov(self,sequence):   ### Adds sequence to BuildMarkov using PoolGen and MarkovX parameters
        '''
        Adds sequence to BuildMarkov using PoolGen and MarkovX parameters.
        >> sequence:str = sequence to be added to suftree
        '''
        if self.opt['PoolGen']:
            self.obj['BuildMarkov']._addToSufTree(sequence,self.stat['MarkovX'])
#########################################################################################################################
    def _debugTree(self,remlen,extseq):    ### Checks integrity of sufdic counts etc.
        '''
        Checks integrity of sufdic counts etc. Has adding/removing from pool screwed things up?
        >> remlen:int = Remaining length for current sequence
        '''
        ### Tree ###
        sufdic= self.obj['BuildMarkov'].suftree()
        self._checkSufDic(sufdic,extseq)

        ### Pool ###
        if not self.opt['DeBug']:
            return
        if self.opt['PoolGen']:
            for slen in self.debuglenlist[1:]:
                remlen += slen
        if remlen <= self.obj['BuildMarkov'].suftree()['=']:
            return
        self.log.printLog('#DEBUG','Suftree has %d at root. Remaining length = %d.' % (self.obj['BuildMarkov'].suftree()['='],remlen))
        max_xmer = self.obj['BuildMarkov'].getMaxXmer(self.obj['BuildMarkov'].suftree())
        for x in range(max_xmer):
            self.log.printLog('#DEBUG','%s counts of %s different %dmers.' % (rje.integerString(self.obj['BuildMarkov'].suftree()[x+1]['=']),rje.integerString(self.obj['BuildMarkov'].suftree()[x+1]['X']),(x+1)))
        self.deBug('')
#########################################################################################################################
    def _checkSufDic(self,sufdic,extseq):
        if sufdic.keys() == ['=']: return
        count = sufdic['=']
        acount = 0
        for a in self.obj['BuildMarkov'].alphabet:
            if sufdic.has_key(a):
                acount += sufdic[a]['=']
                self._checkSufDic(sufdic[a],extseq)
        if count != acount:
            self.deBug('%s ... Total count %d != Summed count %d.\nSufdic: %s' % (extseq,count,acount,sufdic))
#########################################################################################################################
    ### <3> ### Sequence Details Methods                                                                                #
#########################################################################################################################
    def _randSeqName(self,seqx,seqlen):    ### Generates random sequence name/description using self.attributes
        '''
        Generates random sequence name/description using self.attributes.
        >> seqx:int = Sequence number
        >> seqlen:int= Sequence length
        << name:str = sequence name
        '''
        name = '%s%s' % (self.info['RandName'],rje.preZero(seqx+self.stat['IDMin']-1,self.stat['IDMax']))
        if self.opt['RandDesc']:
            if self.opt['Scramble']:
                desc = 'RJE_SEQGEN scrambled sequence'
                if self.stat['ScreenX'] > 0:
                    desc += ' screening %dmer linear motifs' % self.stat['ScreenX']
                if self.stat['Teiresias'] > 0:
                    desc += ' (screening %dmer TEIRESIAS motifs)' % self.stat['Teiresias']
            else:                
                desc = 'Random sequence using %d Order Markov Chains from %s.' % (self.stat['MarkovX'],self.info['Name'])
                if self.info['ScreenFile'] != 'None' and self.stat['ScreenX'] > 0:
                    desc = '%s Constrained to exclude %d Order Xmers from %s, ScreenRev:%s.' % (desc,self.stat['ScreenX'],self.info['ScreenFile'],self.opt['ScreenRev'])
            name = '%s %s (%d letters)' % (name,desc,seqlen)
        return name
#########################################################################################################################
    ### <4> ### Peptide Scramble Methods                                                                                #
#########################################################################################################################
    def scramble(self):     ### Generates scramble peptides
        '''Generates scramble peptides.'''
        try:
            ### Setup Output ###
            outfile = self.info['OutFile']
            if self.opt['FullScramble'] and outfile.lower() in ['','none','randseq.fas']:
                outfile = 'scramble.tdt'
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=outfile))
            rje.backup(self,outfile)     # Checks append and existence
            if self.opt['FullScramble'] and not self.opt['Append']:
                rje.writeDelimit(outlist=['peptide','scramble'],delimit=delimit,outfile=outfile)

            ### Setup Peptides ###
            if not self.list['Peptides']:
                pepseq = rje_seq.SeqList(self.log,self.cmd_list)
                for seq in pepseq.seq:
                    self.list['Peptides'].append(seq.info['Sequence'])

            ### Setup Teiresias ###
            teiresias = rje.makePath(self.info['TeiresiasPath'],wholepath=True)

            ### Setup simple scramble ###
            if self.stat['SeqNum'] < 1:
                self.log.errorLog('Trying to generate %d sequences! Reset to 1.' % self.stat['SeqNum'],False,False)
            if self.stat['IDMax'] >= 0 and self.stat['IDMax'] < self.stat['IDMin']:
                self.stat['IDMax'] = self.stat['SeqNum']
            if self.stat['ScrambleCyc'] < 1:
                self.stat['ScrambleCyc'] = 1
            randname = self.info['RandName']

            ### Scramble! ###
            if self.opt['FullScramble']:
                for peptide in self.list['Peptides']:
                    aalist = []
                    for a in peptide:
                        aalist.append(a)
                    varlist = self.grow('',aalist,peptide)
                    varx = len(varlist)
                    while varlist:
                        self.log.printLog('\r#SCRAMBLE','Output of %s scrambled peptides for %s to %s          ' % (rje.integerString(len(varlist)),peptide,outfile),log=False,newline=False)
                        scramble = varlist.pop(0)
                        if self.stat['Teiresias'] > 0:
                            tmpscr = 'tmp%s' % rje.randomString(5)
                            SCRAM = open('%s.fas' % tmpscr,'w')
                            SCRAM.write('>peptide 1\n%s\n>scramble 2\n%s\n' % (peptide,scramble))
                            SCRAM.close()
                            command = teiresias + ' -i%s.fas -o%s.out -l%d -w10 -c1 -k2 -p' % (tmpscr,tmpscr,self.stat['Teiresias'])
                            #X#os.system(command)
                            os.popen(command).read()
                            pattern = False
                            if os.path.exists('%s.out' % tmpscr):
                                SCRAM = open('%s.out' % tmpscr,'r')
                                scrampat = SCRAM.readlines()
                                SCRAM.close()
                                os.unlink('%s.out' % tmpscr)
                                for line in scrampat:
                                    if rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',line):
                                        pattern = True
                                        break
                            if os.path.exists('%s.fas' % tmpscr):
                                os.unlink('%s.fas' % tmpscr)
                            if pattern:     # TEIRESIAS match: Skip!
                                varx -= 1
                                continue
                        rje.writeDelimit(outlist=[peptide,scramble],delimit=delimit,outfile=outfile)
                    self.log.printLog('\r#SCRAMBLE','Output of %s scrambled peptides for %s to %s complete!' % (rje.integerString(varx),peptide,outfile))

            ### Simpler scrambling ###
            else:
                nrlist = []
                for seq in pepseq.seq:
                    peptide = seq.info['Sequence']
                    self.info['RandName'] = seq.shortName() + '_' + randname
                    for s in range(self.stat['SeqNum']):
                        cyc = 0
                        newpep = ''
                        while cyc < self.stat['ScrambleCyc']:
                            cyc += 1
                            scramble = rje.join(rje.randomList(rje.strList(peptide)),'')
                            self.log.printLog('\r#SCR','Scrambling %s %4s: %7s' % (seq.shortName(),s+1,cyc),newline=False,log=False)
                            ## Check ScreenX violation ##
                            r = 0
                            while len(scramble[r:(r+self.stat['ScreenX'])]) == self.stat['ScreenX']:
                                if peptide.find(scramble[r:(r+self.stat['ScreenX'])]) >= 0:     # Bad!
                                    r = -1
                                    break
                                r += 1
                            if r < 0:
                                continue
                            ## Check Teiresias violation ##
                            if self.stat['Teiresias'] > 0:
                                tmpscr = 'tmp%s' % rje.randomString(5)
                                SCRAM = open('%s.fas' % tmpscr,'w')
                                SCRAM.write('>peptide 1\n%s\n>scramble 2\n%s\n' % (peptide,scramble))
                                SCRAM.close()
                                command = teiresias + ' -i%s.fas -o%s.out -l%d -w10 -c1 -k2 -p' % (tmpscr,tmpscr,self.stat['Teiresias'])
                                #X#os.system(command)
                                os.popen(command).read()
                                pattern = False
                                if os.path.exists('%s.out' % tmpscr):
                                    SCRAM = open('%s.out' % tmpscr,'r')
                                    scrampat = SCRAM.readlines()
                                    SCRAM.close()
                                    os.unlink('%s.out' % tmpscr)
                                    for line in scrampat:
                                        if rje.matchExp('^(\d+)\s+(\d+)\s+(\S+)',line):
                                            pattern = True
                                            break
                                if os.path.exists('%s.fas' % tmpscr):
                                    os.unlink('%s.fas' % tmpscr)
                                if pattern:     # TEIRESIAS match: Skip!
                                    continue
                            ## Check redundancy ##
                            if self.opt['NRGen'] and scramble in nrlist:
                                continue
                            ## Make newpep ##
                            newpep = scramble
                            newname = self._randSeqName(s+1,len(newpep))
                            open(outfile,'a').write('>%s\n%s\n' % (newname,newpep))
                            nrlist.append(newpep)
                            break
                        ## Check cyc ##
                        if not newpep:
                            self.log.printLog('\r#CYC','%s cycles exceeded for %s peptide %d. Giving up.' % (rje.integerString(self.stat['ScrambleCyc']),seq.shortName(),s))
                            break
                self.log.printLog('\r#SCR','Scrambling of peptides complete!')

        except:
            self.log.errorLog('SeqGen.scramble() doesn\'t fecking work. :o(',quitchoice=True)
#########################################################################################################################
    def grow(self,base,aalist,peptide): ### Grows peptide from base and returns full list of scrambles
        '''
        Grows peptide from base and returns full list of scrambles.
        >> base:str = base sequence for peptide growth
        >> aalist:list = list of aas left to add to peptide
        >> peptide:str = original peptide sequence
        << varlist:list = list of scrambled peptides
        '''
        try:
            ### Check base for ScreenX violation ###
            if len(base) >= self.stat['ScreenX'] and peptide.find(base[-self.stat['ScreenX']:]) >= 0:   # Screen!
                return []
            self.log.printLog('\r#SCRAMBLE','Scrambling %s: %s%s' % (peptide,base,' ' * (len(peptide) - len(base))),log=False,newline=False)

            ### Full Length? ###
            if not aalist:  # No more aas to add!
                return [base]
            
            ### Grow ###
            varlist = []
            for a in rje.sortUnique(aalist,False):
                passlist = aalist[0:]
                passlist.remove(a)
                varlist += self.grow(base+a,passlist,peptide)
            return varlist

        except:
            self.log.errorLog('SeqGen.grow() doesn\'t fecking work. :o(',quitchoice=True)
#########################################################################################################################
    ### <5> ### BLASTGen Random Dataset Generator                                                                       #
#########################################################################################################################
    def blastGen(self):     ### Generates random dataset using BLAST similarities and AA frequencies
        '''Generates random dataset using BLAST similarities and AA frequencies.'''
        try:
            ### Input Dataset/BLAST ###
            seqin = rje_seq.SeqList(self.log,self.cmd_list+['autoload=T','autofilter=F','accnr=F','seqnr=F'])
            self.info['Name'] = seqin.info['Name']
            blast = rje_blast.BLASTRun(self.log,self.cmd_list+['blastf=F'])
            blastres = '%s.self.blast' % seqin.info['Name']
            blast.setInfo({'Name':blastres,'DBase':seqin.info['Name'],'InFile':seqin.info['Name']})
            blast.setStat({'OneLine':seqin.seqNum(),'HitAln':seqin.seqNum()})
            blast.formatDB(force=False)
            blast.blast(cleandb=False,use_existing=True)    #!# Will leave BLAST DB Files #!#
            blast.readBLAST(clear=True,gablam=False,unlink=False)   #!# Add Unlink later #!#

            ### Output Dataset Setup ###
            #X#self.opt['RandDesc'] = False
            seqout = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','autoload=F'])
            for seq in seqin.seq:
                seqout.seq.append(rje_sequence.Sequence(log=self.log))
                #!# Add option to keep shortName() or not - will need for mapping? #!#
                if self.opt['KeepNames']:
                    seqout.seq[-1].info['Name'] = seq.shortName() + ' ' + self._randSeqName(seqout.seqNum(),seq.aaLen())
                else:
                    seqout.seq[-1].info['Name'] = self._randSeqName(seqout.seqNum(),seq.aaLen())
                seqout.seq[-1].info['Sequence'] = 'Z' * seq.aaLen()

            ### Make new sequences ###
            markovbuild = rje_markov.Markov(self.log,self.cmd_list+['xcount=F','markov=F','xmers=1'])
            markovbuild.run()
            markovbuild._buildSuffixTree(markovbuild.obj['SeqList'])
            #?!#markovbuild._buildProbTree()
            self.obj['BuildMarkov']= markovbuild
            seqi = rje.randomList(range(seqin.seqNum()))     ### Indices of input sequences in random order
            for i in seqi:
                seq = seqin.seq[i]
                rseq = seqout.seq[i]
                ## Complete newseq generation ##
                sequence = ''
                while len(sequence) < seq.aaLen():
                    if rseq.info['Sequence'][len(sequence)] == 'Z':
                        sequence = self._addAA(sequence)
                    else:
                        sequence = sequence + rseq.info['Sequence'][len(sequence)]
                rseq.info['Sequence'] = sequence
                self.log.printLog('#RAND','Random generation of %s complete.' % rseq.info['Name'])
                ## BLAST Map ##
                search = blast.search[i]
                hitseq = search.hitSeq(seqin)
                for hit in search.hit:
                    h = seqin.seq.index(hitseq[hit])
                    if seqi.index(h) <= seqi.index(i):      ## Same or earlier sequence
                        continue
                    hseq = seqout.seq[h]
                    self.log.printLog('#MAP','Mapping %d BLAST alignments from %s to %s' % (len(hit.aln),search.info['Name'],hit.info['Name']))
                    for aln in hit.aln:
                        qseq = aln.info['QrySeq']
                        sseq = aln.info['SbjSeq']
                        aseq = aln.info['AlnSeq']
                        q = aln.stat['QryStart'] - 2
                        s = aln.stat['SbjStart'] - 2
                        for a in range(len(qseq)):
                            if qseq[a] != '-':
                                q += 1
                            if sseq[a] != '-':
                                s += 1
                            if qseq[a] == sseq[a]:  # AA and ID
                                hseq.info['Sequence'] = rje.strSub(hseq.info['Sequence'],s,s,rseq.info['Sequence'][q])
                        #!# Add check?    - QryEnd    - SbjEnd
                
            ### Output ###
            seqout.info['Name'] = self.info['OutFile']
            seqout.saveFasta(append=self.opt['Append'])
            
            ### Finish ###
            #option to run rje_blast.cleanupDB(blast,blast.info['DBase'])
        except:
            self.log.errorLog('Major error with SeqGen.BLASTGen()',quitchoice=True)
            return False
#########################################################################################################################
    ### <6> ### ESTGen Random Dataset Generator                                                                         #
#########################################################################################################################
    def estGen(self):     ### Generates random dataset using individual sequence details.
        '''Generates random dataset using individual sequence details.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Setup Markov chain object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.stat['MarkovX'] < 1:
                self.log.errorLog('Trying to use %d Order Markov Chain! Reset to 1st order.' % self.stat['MarkovX'],False,False)
                self.stat['MarkovX'] = 1
            self.obj['BuildMarkov'] = markovbuild = rje_markov.Markov(log=self.log,cmd_list=self.cmd_list+['seqin=None'])
            markovbuild.setStat({'MinXmer':1,'MaxXmer':self.stat['MarkovX']})
            markovbuild.opt['Silent'] = True
            ## ~ [0b] Setup Sequence List objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            inseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['autoload=T','autofilter=F','accnr=F','seqnr=F','gnspacc=T'])     # Given input file
            if not inseq.seq: return self.log.errorLog('No sequences loaded from "%s"' % inseq.info['Name'],printerror=False)
            self.info['Name'] = inseq.info['Name']
            if inseq.info['Type'].lower() == 'dna': markovbuild.alphabet = ['A','C','G','T','N']
            if inseq.info['Type'].lower() == 'rna': markovbuild.alphabet = ['A','C','G','U','N']
            buildseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None'])  # Will hold singles sequences
            buildseq.info['Name'] = 'seqgen.tmp.%s' % rje.randomString(8)
            markovbuild.obj['SeqList'] = buildseq
            outseq = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+['seqin=None'])    # Will hold output sequences
            outseq.info['Name'] = self.info['OutFile']
            ## ~ [0c] ~ Setup additional self attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.stat['IDMin'] = max(1,self.stat['IDMin']) + 1
            self.stat['SeqNum'] = inseq.seqNum() 
            self.stat['IDMax'] = self.stat['IDMin'] + self.stat['SeqNum'] - 1
            self.opt['PoolGen'] = False

            ### ~ [1] Generate new sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sx,snum) = (0.0,self.stat['SeqNum'])
            for seq in inseq.seq:
                self.progLog('\r#GEN','Generating %s seqs: %.2f%%' % (rje.integerString(snum),(sx/snum)))
                sx += 100.0
                ## ~ [1a] Build Markov chain for sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sequence = seq.info['Full'] = seq.info['Sequence']
                t = 0   # Look for run of 5' Ts (i.e. rev complement of poly-A tail
                if inseq.dna():
                    while sequence[t] == 'T': t += 1
                #a = -1
                #while sequence[a] == 'A': a -= 1
                #if t-a > 10: print 'T -> %d; A -> %d' % (t,(-a - 1))
                #continue
                seq.info['Sequence'] = sequence[t:]
                markovbuild.obj['SeqList'].seq = [seq]
                markovbuild._buildSuffixTree(markovbuild.obj['SeqList'])
                self.deBug(markovbuild.suftree())
                markovbuild.stripToX(markovbuild.suftree(),self.stat['MarkovX'])
                self.deBug(markovbuild.suftree())
                #?!#markovbuild._buildProbTree()
                self.deBug(markovbuild.suftree())
                ## ~ [1b] Generate new sequence using length and markov chain of seq ~~~~~~~~~~~~~~ ##
                newseq = 'T' * t + self._makeSeq(seq.aaLen())
                newname = self._randSeqName(outseq.seqNum(),len(newseq))
                outseq._addSeq(newname,newseq)
                seq.info['Sequence'] = seq.info['Full']
                self.deBug(newname)
            #return
            self.printLog('\r#GEN','Generating %s seqs complete!' % (rje.integerString(snum)))
            outseq.saveFasta(append=self.opt['Append'])

            ### ~ [2] ~ Clean up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mfile in glob.glob('%s*' % buildseq.info['Name']): os.unlink(mfile)

        except: self.log.errorLog(rje_zen.Zen().wisdom())    
#########################################################################################################################
## End of SECTION II: SeqGen Class                                                                                      #
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
    except: rje.printf('Unexpected error during program setup:', sys.exc_info()[0]); return

    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try:
        seqgen = SeqGen(mainlog,cmd_list)
        if seqgen.opt['Scramble']: seqgen.scramble()
        elif seqgen.opt['BLASTGen']: seqgen.blastGen()
        elif seqgen.opt['ESTGen']: seqgen.estGen()
        else: seqgen.seqGen()
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
