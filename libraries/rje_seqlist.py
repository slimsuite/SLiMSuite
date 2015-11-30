#!/usr/bin/python

# See below for name and description
# Copyright (C) 2011 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_seqlist
Description:  RJE Nucleotide and Protein Sequence List Object (Revised)
Version:      1.15.3
Last Edit:    16/11/15
Copyright (C) 2011  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is designed to replace rje_seq. The scale of projects has grown substantially, and rje_seq cannot deal
    well with large datasets. An important feature of rje_seqlist.SeqList objects, therefore, is to offer different
    sequence modes for different applications. To simplify matters, rje_seqlist will now only cope with single format
    sequences, which includes a single naming format. 

    This version of the SeqList object therefore has several distinct modes that determine how the sequences are stored.
    - full = Full loading into Sequence Objects.
    - list = Lists of (name,sequence) tuples only.
    - file = List of file positions.
    - index = No loading of sequences. Use index file to find sequences on the fly.
    - db = Store sequence data in database object.

SeqShuffle:
    Version 1.2 introduced the seqshuffle function for randomising input sequences. This generates a set of biologically
    unrealistic sequences by randomly shuffling each input sequence without replacement, such that the output sequences
    have the same primary monomer composition as the input but any dimer/trimer biases etc. are removed. This is executed
    by the shuffleSeq() method, which can also generate sequences shuffled with replacement, i.e. based on frequencies.

Sampler:
    Version 1.5 introduced a sequence sampling function for pulling out a random selection of input sequences into one or
    more output files. This is controlled by `sampler=N(,X)` where the X setting is optional. Random selections of N
    sequences will be output into a file named according to the `seqout=FILE` option (or the input file appended with
    `.nN` if none given). X defines the number of replicate datasets to generate and will be set to 1 if not given.
    If X>1 then the output filenames will be appended with `.rx` for each replicate, where x is 1 to X. If 0.0 < N < 1.0
    then a proportion of the input sequences (rounding to the nearest integer) will be selected.

SortSeq:
    In Version 1.8, the `sizesort=T/F` function is replaced with `sortseq=X` (or `seqsort=X`), where X is a choice of:
    - size = Sort sequences by size small -> big
    - accnum = Alphabetical by accession number
    - name = Alphabetical by name
    - seq[X] = Alphabetical by sequence with option to use first X aa/nt only (to save memory)
    - species = Alphabetical by species code
    - desc = Alphabetical by description
    - invsize = Sort by size big -> small re-output prior to loading/filtering (old sizesort - still sets sortseq)
    - invX / revX (Note adding `inv` or `rev` in front of any selection will reverse sort.)

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Sequence input file name. [None]
    seqmode=X       : Sequence mode, determining method of sequence storage (full/list/file/index/db). [file]
    seqdb=FILE      : Sequence file from which to extract sequences (fastacmd/index formats) [None]
    seqindex=T/F    : Whether to save (and load) sequence index file in file mode. [True]
    seqformat=X     : Expected format of sequence file [None]
    seqtype=X       : Sequence type (prot(ein)/dna/rna/mix(ed)) [None]
    mixed=T/F       : Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    dna=T/F         : Alternative option to indicate dealing with nucleotide sequences [False]
    autoload=T/F    : Whether to automatically load sequences upon initialisation. [True]
    autofilter=T/F  : Whether to automatically apply sequence filtering. [True]

    ### ~ SEQUENCE FORMATTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    reformat=X      : Output format for sequence files (fasta/short/acc/acclist/speclist/index/dna2prot/peptides/(q)region) [fasta]
    rename=T/F      : Whether to rename sequences [False]
    spcode=X        : Species code for non-gnspacc format sequences [None]
    newacc=X        : New base for sequence accession numbers - will rename sequences [None]
    newgene=X       : New gene for renamed sequences (if blank will use newacc or 'seq') [None]
    concatenate=T   : Concenate sequences into single output sequence named after file [False]
    split=X         : String to be inserted between each concatenated sequence [''].
    seqshuffle=T/F  : Randomly shuffle each sequence without replacement (maintains monomer composition) [False]
    region=X,Y      : Alignment/Query region to use for peptides/(q)region reformatting of fasta alignment (1-L) [1,-1]

    ### ~ DNA TRANSLATIONS (reformat=dna2prot) ~~~~~~~~~~~~ ###
    minorf=X        # Min. ORF length for translated sequences output. -1 for single translation inc stop codons [-1]
    terminorf=X     # Min. length for terminal ORFs, only if no minorf=X ORFs found (good for short sequences) [-1]
    orfmet=T/F      # Whether ORFs must start with a methionine (before minorf cutoff) [True]
    rftran=X        # No. reading frames (RF) into which to translate (1,3,6) [1]

    ### ~ FILTERING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqnr=T/F       : Whether to check for redundancy on loading. (Will remove, save and reload if found) [False]
    revcompnr=T/F   : Whether to check reverse complement for redundancy too [False]
    goodX=LIST      : Inclusive filtering, only retaining sequences matching list []
    badX=LIST       : Exclusive filtering, removing sequences matching list []
    - where X is 'Acc', Accession number; 'Seq', Sequence name; 'Spec', Species code; 'Desc', part of name;
    minlen=X        : Minimum sequence length [0]
    maxlen=X	    : Maximum length of sequences (<=0 = No maximum) [0]

    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqout=FILE     : Whether to output sequences to new file after loading and filtering [None]
    usecase=T/F     : Whether to return sequences in the same case as input (True), or convert to Upper (False) [False]
    sortseq=X       : Whether to sort sequences prior to output (size/invsize/accnum/name/seq/species/desc) [None]
    sampler=N(,X)   : Generate (X) file(s) sampling a random N sequences from input into seqout.N.X.fas [0]
    summarise=T/F   : Generate some summary statistics in log file for sequence data after loading [False]
    splitseq=X      : Split output sequence file according to X (gene/species) [None]

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, re, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_sequence, rje_uniprot, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation. Based on rje_seq 3.10.
    # 0.1 - Added basic species filtering and sequence output.
    # 0.2 - Added upper case filtering.
    # 0.3 - Added accnum filtering and sequence renaming.
    # 0.4 - Added sequence redundancy filtering.
    # 0.5 - Added newgene=X for sequence renaming (newgene_spcode__newaccXXX). NewAcc no longer fixed Upper Case.
    # 1.0 - Upgraded to "ready" Version 1.0. Added concatenate=T and split=X options for sequence concatenation.
    # 1.0 - Added reading of sequence type from rje_seq.py and mixed=T/F.
    # 1.1 - Added shortName() and modified SeqDict.
    # 1.2 - Added seqshuffle option for randomising sequences.
    # 1.3 - Modified use of index file (appends, not replaces, file extension)
    # 1.4 - Added dna2prot reformat function.
    # 1.5 - Added sampler=N(,X)   : Generate (X) file(s) sampling a random N sequences from input into seqout.N.X.fas [0]
    # 1.6 - Modified currSeq() and nextSeq() slightly to fix index mode breakage. Look out for other programs breaking.
    # 1.6 - Add sequence fragment extraction.
    # 1.7 - Added code to create rje_sequence.Sequence objects.
    # 1.8 - Added sortseq=X : Whether to sort sequences prior to output (size/invsize/accnum/name/seq/species/desc) [None]
    # 1.9.0 - Added extra functions for returning sequence AccNum, ID or Species code.
    # 1.10.0 - Added extraction of uniprot IDs for seqin.
    # 1.11.0 - Added more dna2prot reformatting options.
    # 1.12.0 - Added peptides/qregion reformatting and region=X,Y.
    # 1.13.0 - Added summarise=T option for generating some summary statistics for sequence data. Added minlen & maxlen.
    # 1.14.0 - Added splitseq=X split output sequence file according to X (gene/species) [None]
    # 1.15.0 - Added names() method.
    # 1.15.1 - Fixed bug with storage and return of summary stats.
    # 1.15.2 - Fixed dna2prot reformatting.
    # 1.15.3 - Fixed summarise bug (n=1).
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Make sure index date is checked versus sequence file.
    # [Y] : Add seqnr redundancy filter.
    # [Y] : Add reverse complement redundancy filter?
    # [ ] : Check revcompnr filter is working?
    # [ ] : Add sequence masking based on rje_slimcore: mask and add second SeqList object containing masked sequences.
    # [ ] : Add method/options to return old-style rje_seq.SeqList object.
    # [Y] : Add additional sorting methods using sort=X: size/name/sequence/revsize/accnum
    # [ ] : Fix dna2prot reformatting output using minlen.
    # [ ] : Add assemble=FILE mode with reformat=assemble: list sequences on each line to join/revcomp into single seqs.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SeqList', '1.15.3', 'November 2015', '2011')
    description = 'RJE Nucleotide and Protein Sequence List Object (Revised)'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copy_right,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        ### ~ [2] ~ Look for help commands and print options if found ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
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
### SECTION II: SeqList Class                                                                                           #
#########################################################################################################################
class SeqList(rje_obj.RJE_Object):     
    '''
    SeqList Class. Author: Rich Edwards (2011).

    Str:str
    - Name = Sequence file name - specifies output. [None]
    - NewAcc = New base for sequence accession numbers - will rename sequences [None]
    - NewGene = New gene for renamed sequences (if blank will use newacc) [None]
    - Region = Query region to use for peptides/qregion reformatting of fasta alignment [0,-1]
    - ReFormat = Output format for sequence files (fasta/short/acc/acclist/speclist) [fasta]
    - SeqDB = Sequence file from which to extra sequences (fastacmd/index formats) [None]
    - SeqDictType = String identifier of the type of sequence dictionary made (accnum/short/name/max) [None]
    - SeqFormat = Expected format of sequence file [None]
    - SeqIn = Sequence input file name. [None]
    - SeqMode = Sequence mode, determining method of sequence storage (full/list/file/index/db). [list]
    - SeqOut = Whether to output sequences to new file after loading and filtering [None]
    - SeqType = Sequence type (prot(ein)/dna/rna/mix(ed)) [None]
    - SortSeq = Whether to sort sequences prior to output (size/invsize/accnum/name/seq/species/desc) [None]
    - SpCode = Species code for non-gnpacc format sequences [None]
    - Split = String to be inserted between each concatenated sequence [''].
    - SplitSeq = Split output sequence file according to X (gene/species) [None]

    Bool:boolean
    - AutoFilter = Whether to automatically apply sequence filtering. [True]
    - AutoLoad = Whether to automatically load sequences upon initialisation. [True]
    - Concatenate = Concatenate sequences into single output sequence named after file [False]
    - DNA = Alternative option to indicate dealing with nucleotide sequences [False]
    - Mixed = Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    - ORFMet = Whether ORFs must start with a methionine (before minorf cutoff) [True]
    - ReName = Whether to rename sequences (will need newacc and spcode) [False]
    - RevCompNR = Whether to check reverse complement for redundancy too [False]
    - SeqIndex = Whether to save (and load) sequence index file in file mode. [True]
    - SeqNR = Whether to check for redundancy on loading. (Will remove, save and reload if found) [False]
    - SeqShuffle = Randomly shuffle each sequence (cannot use file or index mode) [False]
    - SizeSort = Sort sequences by size big -> small re-output prior to loading/filtering [False]
    - Summarise = Generate some summary statistics in log file for sequence data after loading [False]
    - UseCase = Whether to return sequences in the same case as input (True), or convert to Upper (False) [False]

    Int:integer
    - MinLen = Minimum sequence length [0]
    - MaxLen = Maximum length of sequences (<=0 = No maximum) [0]
    - MinORF = Min. ORF length for translated sequences output. -1 for single translation inc stop codons [-1]
    - RFTran = No. reading frames (RF) into which to translate (1,3,6) [1]
    - TerMinORF = Min. length for terminal ORFs, only if no minorf=X ORFs found (good for short sequences) [-1]

    Num:float
    
    List:list
    - Sampler = N(,X) = Generate (X) file(s) sampling a random N sequences from input into seqout.N.X.fas [0]
    - Seq = List of sequences - nature depends on SeqMode []

    Dict:dictionary
    - SeqDict = Dictionary of {key:seq object}

    Obj:RJE_Objects
    - Current = Current sequence "object" from list. [None]
    - CurrSeq = Formatted current sequence "object" from list. [None]
    - DB = Database object. Used to store sequence data and index data. [None]
    - SEQFILE = Stores open file handle for reading sequences in file mode [None]
    - INDEX = Stores open index file handle for reading sequences in index mode [None]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Name','NameFormat','NewAcc','Region',
                        'SeqDB','SeqDictType','SeqFormat','SeqIn','SeqMode','SeqType','SeqOut',
                        'Reformat','SpCode','SeqNR','NewGene','Split','SortSeq','SplitSeq']
        self.boollist = ['AutoFilter','AutoLoad','Concatenate','DNA','Mixed','ORFMet','ReName','RevCompNR','SizeSort',
                         'SeqIndex','SeqShuffle','Summarise','UseCase']
        self.intlist = ['MinLen','MaxLen','MinORF','RFTran','TerMinORF']
        self.numlist = []
        self.listlist = ['Sampler','Seq']
        self.dictlist = ['Filter','SeqDict']
        self.objlist = ['Current','CurrSeq','DB','SEQFILE','INDEX']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'SeqMode':'file','ReFormat':'fasta','Region':'1,-1'})
        self.setBool({'AutoFilter':True,'AutoLoad':True,'ORFMet':True,'SeqIndex':True})
        self.setInt({'MinORF':-1,'RFTran':1,'TerMinORF':-1})
        self.list['Sampler'] = [0,1]
        #self.setInt({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### ~ [1] ~ Commandline options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for cmd in self.cmd_list:
            try:
                self._generalCmd(cmd)   ### General Options ### 
                self._forkCmd(cmd)  # Delete if no forking
                ### Class Options ### 
                self._cmdRead(cmd,type='file',att='SeqDB',arg='fasdb')  # No need for arg if arg = att.lower()
                self._cmdRead(cmd,type='str',att='SortSeq',arg='seqsort')  # No need for arg if arg = att.lower()
                self._cmdReadList(cmd,'str',['NewAcc','NewGene','Region','SeqFormat','SeqMode','ReFormat','SpCode','SeqType','Split','SortSeq','SplitSeq'])
                self._cmdReadList(cmd,'file',['SeqDB','SeqIn','SeqOut'])
                self._cmdReadList(cmd,'int',['MinLen','MaxLen','MinORF','RFTran','TerMinORF'])
                self._cmdReadList(cmd,'nlist',['Sampler'])
                self._cmdReadList(cmd,'bool',['Align','AutoFilter','AutoLoad','Concatenate','DNA','Mixed','ORFMet','ReName','RevCompNR','SizeSort','SeqIndex','SeqNR','SeqShuffle','Summarise','UseCase'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ## ~ [1a] ~ Tidy Commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.getStrLC('SeqMode') == 'tuple': self.setStr({'SeqMode':'list'})
        if not self.list['Sampler']: self.list['Sampler'] = [0,0]
        elif len(self.list['Sampler']) < 2: self.list['Sampler'].append(1)
        elif len(self.list['Sampler']) > 2:
            self.warnLog('Too many values given to sampler=N,X. Will use sampler=%s,%s' % self.list['Sampler'][:2])
            self.list['Sampler'] = self.list['Sampler'][:2]
        if not self.getStrLC('SortSeq') and self.getBool('SizeSort'): self.setStr({'SortSeq':'size'})
        else: self.setStr({'SortSeq':self.getStrLC('SortSeq')})
        if self.getInt('RFTran') not in [1,3,6]:
            self.warnLog('rftran=%d not recognised: will use rftran=1' % self.getInt('RFTran'))
            self.setInt({'RFTran':1})
        ## ~ [1b] ~ REST Command setup/adjustment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.getStrLC('Rest') in string.split('fasta/short/acc/acclist/speclist/index/dna2prot/peptides/qregion/region','/'):
            self.setStr({'ReFormat':self.getStr('Rest')})
            self.dict['Output'][self.getStrLC('Rest')] = 'SeqOut'
        if self.getStrLC('ReFormat') == 'dna2prot' and self.getStrLC('Rest'): self.setBool({'DNA':True})
        ### ~ [2] ~ AutoLoad option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('DNA') and not self.nt(): self.str['SeqType'] = 'dna'
        self._filterCmd()
        if self.getBool('AutoLoad'):
            if self.getBool('SeqShuffle'): self.shuffleSeq()
            #if self.getBool('SizeSort'): self.sizeSort()
            if self.getStr('SortSeq'): self.seqSort()
            if self.getBool('SeqNR'): self.seqNR(twopass=not self.getStrLC('SortSeq'))
            if self.getBool('Concatenate'): self.concatenate()
            if self.getBool('ReName'): self.rename()
            self.loadSeq()
            if self.getBool('AutoFilter'): self.filterSeqs(screen=self.v()>0)
            if self.getBool('Summarise'): self.summarise()
            if self.list['Sampler'][0] > 0: self.sampler()
            elif self.getStrLC('SplitSeq'): self.splitSeq()
            elif self.getStr('SeqOut').lower() not in ['','none']: self.saveSeq()
            elif self.getStrLC('ReFormat') and self.getStrLC('ReFormat') != 'fasta':
                if self.getStrLC('Basefile'):  seqout = '%s.%s.fas' % (self.baseFile(),self.getStrLC('ReFormat'))
                else:  seqout = '%s.%s.fas' % (rje.baseFile(self.getStr('SeqIn')),self.getStrLC('ReFormat'))
                if self.i() < 0 or rje.yesNo('Reformat (%s) and save to %s?' % (self.getStrLC('ReFormat'),seqout)):
                    self.saveSeq(seqfile=seqout)
#########################################################################################################################
    def _filterCmd(self,cmd_list=None,clear=True):   ### Reads filter commands into attributes
        '''
        Reads filter commands into attributes.
        >> cmd_list:list of commands from which to get filter options [None = self.cmd_list]
        '''
        #!# This method needs to be edited to introduce full filtering #!#
        ### ~ [1] ~ Setup filter attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if 'Filter' not in self.dict: self.dict['Filter'] = {}
        self.dict['Filter']['MinLen'] = 0; self.dict['Filter']['MaxLen'] = 0
        for filt in ['Acc','Seq','Spec','DB','Desc']:
            self.dict['Filter']['Good%s' % filt] = 0
            if clear: self.list['Good%s' % filt] = []
            self.dict['Filter']['Bad%s' % filt] = 0
            if clear: self.list['Bad%s' % filt] = []

        ### ~ [2] ~ Commandline options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if cmd_list == None: cmd_list = self.cmd_list
        for cmd in cmd_list:
            try:
                #self._cmdRead(cmd,type='str',att='FilterOut')
                #self._cmdReadList(cmd,'int',['MinLen','MaxLen'])
                #self._cmdReadList(cmd,'stat',['MaxGap','MaxX','MaxGlob'])
                #self._cmdRead(cmd,type='stat',att='NR ID',arg='nrid')
                #self._cmdRead(cmd,type='stat',att='NR Sim',arg='nrsim')
                #self._cmdRead(cmd,type='opt',att='NR Align',arg='nralign')
                #self._cmdRead(cmd,type='clist',att='DBList')
                #self._cmdReadList(cmd,'opt',['LogRem','DBOnly','UnkSpec','AccNR','SeqNR','SpecNR','QueryNR'])
                for filt in ['Acc','Seq','Spec','DB','Desc']:
                    self._cmdRead(cmd,type='list',att='Good%s' % filt)
                    self._cmdRead(cmd,type='list',att='Bad%s' % filt)
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ### ~ [3] ~ Check Filtering commands against sequence mode and warn if any will be ignored ~~~~~~~~~~~~~~~~~~ ###
        #!# Add this code
#########################################################################################################################
    ### <2> ### Main attribute retrieval                                                                                #
#########################################################################################################################
    def mode(self): return self.getStr('SeqMode').lower()
    def nt(self): return self.getStr('SeqType').lower()[-2:] == 'na'
    def dna(self): return self.getStr('SeqType').lower() == 'dna'
    def rna(self): return self.getStr('SeqType').lower()[-3:] == 'rna'
    def protein(self): return self.getStr('SeqType').lower()[:4] == 'prot'
    def seqNum(self): return len(self.list['Seq'])
    def progress(self): return '%.2f%%' % (100.0 * self.list['Seq'].index(self.obj['Current']) / self.seqNum())
#########################################################################################################################
    ## <2a> ## Sequence retrieval                                                                                       #
#########################################################################################################################
    def seqs(self,copy=False):
        if copy: return self.list['Seq'][0:]
        else: return self.list['Seq']
    def current(self): return self.currSeq()
    def currSeq(self):  ### Returns the current sequence of interest
        if self.obj['Current'] in self.list['Seq']:
            #!# This cannot work for tuples #!# self.obj['CurrSeq'] = self.getSeq(self.obj['Current'])
            return self.getSeq(self.obj['Current']) #self.obj['CurrSeq']
        else: return None
#########################################################################################################################
    def nextSeq(self):  ### Returns next sequence in list and updates current
        try:
            if self.obj['Current'] in self.list['Seq']:
                try: self.obj['Current'] = self.list['Seq'][self.list['Seq'].index(self.obj['Current'])+1]
                except: return None    # Reached end of sequence list
            else: self.obj['Current'] = self.list['Seq'][0]
            #self.obj['CurrSeq'] = self.currSeq()   #!# Replacing with getSeq messes everything up!
            return self.currSeq()   #self.obj['CurrSeq']
        except:
            if self.dev(): self.errorLog('Ugg')
            return None
#########################################################################################################################
    def prevSeq(self):  ### Returns previous sequence in list and updates current
        try:
            if self.obj['Current'] in self.list['Seq']: 
                try: self.obj['Current'] = self.list['Seq'][self.list['Seq'].index(self.obj['Current'])-1]
                except: return None    # Reached start of sequence list
            else: self.obj['Current'] = self.list['Seq'][-1]
            self.obj['CurrSeq'] = self.currSeq()
            return self.obj['CurrSeq']
        except: return None    
#########################################################################################################################
    def SEQFILE(self):  ### Returns open Sequence File handle
        '''Returns open Sequence File handle.'''
        if not self.obj['SEQFILE']: self.obj['SEQFILE'] = open(self.getStr('SeqIn'),'r')
        return self.obj['SEQFILE']           
#########################################################################################################################
    def INDEX(self):  ### Returns open Sequence File handle
        '''Returns open Sequence File handle.'''
        if not self.obj['INDEX']: self.obj['INDEX'] = open('%s.index' % self.getStr('SeqIn'),'r')
        return self.obj['INDEX']           
#########################################################################################################################
    def seqNameDic(self):
        if self.dict['SeqDict']: return self.dict['SeqDict']
        else: return self.makeSeqNameDic()
#########################################################################################################################
    def makeSeqNameDic(self,keytype=None,clear=True,warnings=True):  ### Make SeqDict sequence name dictionary.
        '''
        >> keytype:str [None] = Type of data to use as key for dictionary (accnum/short/name/max)
        >> clear:bool [True] = whether to clear self.dict['SeqDict'] before filling
        >> warnings:bool [True] = whether to warn when keys are getting over-written
        << returns self.dict['SeqDict']
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if keytype: self.setStr({'SeqDictType':keytype.lower()})
            if not self.getStrLC('SeqDictType'): self.setStr({'SeqDictType':'short'})
            keytype = self.getStrLC('SeqDictType')
            if keytype not in ['name','max','short','acc','accnum','id']: raise ValueError('SeqNameDic keytype "%s" not recognised!' % keytype)
            if clear: self.dict['SeqDict'] = {}
            ### ~ [1] Populate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                skeys = []
                name = self.getSeq(seq,'tuple')[0]
                if keytype in ['name','max']: skeys.append(name)
                if keytype in ['short','max']: skeys.append(string.split(name)[0])
                if keytype in ['acc','accnum','max']: skeys.append(string.split(string.split(name)[0],'__')[-1])
                if keytype in ['id','max']: skeys.append(string.split(string.split(name)[0],'__')[0])
                skeys = rje.sortUnique(skeys)
                for skey in skeys:
                    if warnings and skey in self.dict['SeqDict']: self.warnLog('Sequence "%s" already in SeqDict.' % skey,suppress=True)
                    self.dict['SeqDict'][skey] = seq
            return  self.dict['SeqDict']
        except:
            self.errorLog('SeqNameDic problem')
            raise ValueError
#########################################################################################################################
    def shortName(self,seq=None):    ### Returns short name (first word) of given sequence
        '''Returns short name (first word) of given sequence.'''
        if seq == None: seq = self.obj['Current']
        return string.split(self.getSeq(seq,'tuple')[0])[0]
#########################################################################################################################
    def seqLen(self,seq=None):       ### Returns length of given sequence
        '''Returns length of given sequence.'''
        if seq == None: seq = self.obj['Current']
        return len(self.getSeq(seq,'tuple')[1])
#########################################################################################################################
    def seqNonX(self,seq=None):  ### Returns number of resolved positons
        '''Returns number of resolved positons.'''
        if seq == None: seq = self.obj['Current']
        sequence = self.getSeq(seq)[1]
        if self.nt(): return len(sequence) - string.count(sequence.upper(),'N')
        else: return len(sequence) - string.count(sequence.upper(),'X')
#########################################################################################################################
    def aaLen(self,seq=None):  ### Returns number of resolved positons
        '''Returns number of resolved positons.'''
        if seq == None: seq = self.obj['Current']
        sequence = self.getSeq(seq)[1]
        return len(sequence) - sequence.count('-')
#########################################################################################################################
    def seqAcc(self,seq=None): return string.split(self.shortName(seq),'__')[-1]
#########################################################################################################################
    def seqID(self,seq=None): return string.split(self.shortName(seq),'__')[0]
#########################################################################################################################
    def seqGene(self,seq=None): return string.split(self.shortName(seq),'_')[0]
#########################################################################################################################
    def seqSpec(self,seq=None): return string.split(self.seqID(seq),'_')[-1]
#########################################################################################################################
    def seqName(self,seq=None):
        if seq == None: seq = self.obj['Current']
        return self.getSeq(seq,'tuple')[0]
#########################################################################################################################
    def seqDesc(self,seq=None): return string.join(string.split(self.seqName(seq))[1:])
#########################################################################################################################
    def isSwiss(self,seq=None):  ### Returns whether sequence appears to be SwissProt
        '''Returns whether sequence appears to be SwissProt.'''
        sname = self.shortName(seq)
        sid = string.split(sname,'_')[0]
        if sid.upper() != sid: return False
        sacc = string.split(sname,'__')[-1]
        if sid == sacc: return False
        return True
#########################################################################################################################
    def names(self,short=True): ### Returns list of sequence (short) names
        '''Returns list of sequence (short) names.'''
        names = []
        for seq in self.seqs():
            if short: names.append(self.shortName(seq))
            else: names.append(self.seqName(seq))
        return names
#########################################################################################################################
    def getDictSeq(self,dkey,format=None,mode=None,case=None,errors=True):   ### Returns sequence from dictionary key
        '''Returns sequence from dictionary key with option to raise error or return None if missing.'''
        try: return self.getSeq(self.dict['SeqDict'][dkey],format,mode,case)
        except:
            if errors:
                self.deBug('%s not in %s' % (dkey,rje.sortKeys(self.dict['SeqDict'])))
                raise
            else: return None
#########################################################################################################################
    def getSeqObj(self,name,sequence):  ### Returns an rje_sequence.Sequence object for given (name,sequence)
        '''Returns an rje_sequence.Sequence object for given (name,sequence).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newseq = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list,parent=self)
            newseq.setStr({'Name':name,'Type':self.getStr('SeqType')})
            newseq.addSequence(sequence)
            newseq.extractDetails(gnspacc=True)
            return newseq
        except: self.errorLog('SeqList.getSeqObj() error'); return None
#########################################################################################################################
    def getSeq(self,seq=None,format=None,mode=None,case=None):   ### Returns sequence as (name,seqence) tuple, db entry or object as appropriate
        '''
        Returns sequence as (name,seqence) tuple, db entry or object as appropriate.
        >> seq:various = The type of variable will depend on self.mode()
        >> format:str [None] = This can be used to over-ride self.mode() and return a specific format (tuple/entry/obj/short)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq == None: seq = self.obj['Current']
            if not mode: mode = self.mode()
            if case == None: case = self.getBool('UseCase')
            #self.deBug('%s: %s (%s)' % (seq,format,mode))
            ### ~ [1] ~ Return Sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if mode == 'full':   # Full loading into Sequence Objects.
                if format == 'tuple': return (seq.getStr('Name'),seq.getSequence(case))
                elif format == 'entry': return seq.info #?#
                elif format == 'short': return seq.shortName()
                else: return seq
            elif mode == 'list':     # Lists of (name,sequence) tuples only.
                (name,sequence) = seq
                if format == 'obj': return self.getSeqObj(name,sequence)
                elif format == 'entry':
                    if case: return {'Name':name,'Sequence':sequence}
                    else: return {'Name':name,'Sequence':sequence.upper()}
                elif format == 'short': return string.split(name)[0]
                else: return seq
            elif mode == 'db':       # Store sequence data in database object.
                (name,sequence) = (seq['Name'],seq['Sequence'])
                if format == 'obj': return self.getSeqObj(name,sequence)
                elif format == 'tuple':
                    if case: return (name,sequence)
                    else: return (name,sequence.upper())
                elif format == 'short': return string.split(name)[0]
                else: return seq
            elif mode == 'file':     # List of file positions.
                SEQFILE = self.SEQFILE()
                SEQFILE.seek(seq)
                name = rje.chomp(SEQFILE.readline())
                if name[:1] != '>':
                    self.deBug(name)
                    self.errorLog('Given file position that is not name line. May have SeqList objects mixed up?',printerror=False)
                    raise ValueError    # Must be Fasta!
                name = name[1:]
                sequence = ''; line = rje.chomp(SEQFILE.readline())
                while line and line[:1] != '>':
                    sequence += line
                    line = rje.chomp(SEQFILE.readline())
                if format == 'obj': return self.getSeqObj(name,sequence)
                elif format == 'entry': 
                    if case: return {'Name':name,'Sequence':sequence}
                    else: return {'Name':name,'Sequence':sequence.upper()}
                elif format == 'short': return string.split(name)[0]
                else:
                    if case: return (name,sequence)
                    else: return (name,sequence.upper())
            elif mode == 'index':    # No loading of sequences. Use index file to find sequences on the fly.
                INDEX = self.INDEX()
                ipos = rje.posFromIndex(seq,INDEX,re_index='^(\S+)\s')
                iline = rje.fileLineFromSeek(INDEX,ipos)
                #self.deBug(iline)
                fpos = string.atol(string.split(rje.chomp(iline[0]))[1])
                #self.deBug('%s: %d' % (seq,fpos))
                return self.getSeq(fpos,format,mode='file')
        except: self.errorLog('%s.getSeq(%s) error' % (self,type(seq))); return None
#########################################################################################################################
    def getSeqFrag(self,seq=None,fragstart=0,fragend=0,case=None):    ### Returns sequence fragment as (name,seqence) tuple with modified name.
        '''
        Returns sequence fragment as (name,seqence) tuple with modified name.
        >> seq:various = The type of variable will depend on self.mode()
        >> fragstart:int [0] = Fragment start position 1-L. Will be sequence start if <=0
        >> fragend:int [0] = Fragment start position 1-L. Will be sequence start if <=0
        >> case:bool [None] = Whether to use sequence case.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq == None: seq = self.obj['Current']
            if case == None: case = self.getBool('UseCase')
            if fragstart <= 0: fragstart = 1
            ### ~ [1] ~ Crude case not using index ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#if self.mode() not in ['file','index']:

            (seqname,fullseq) = self.getSeq(seq,format='tuple',case=case)
            seqlen = len(fullseq)
            if fragend <= 0 or fragend > seqlen: fragend = seqlen
            seqname = string.split(seqname)
            #i# Note that positions are 1 to L and not 0 to (L-1)
            seqname[0] = '%s-%s.%s' % (seqname[0],rje.preZero(fragstart,seqlen),rje.preZero(fragend,seqlen))
            seqname.insert(1,'(Pos %s - %s)' % (rje.iStr(fragstart),rje.iStr(fragend)))
            #seqname.insert(1,'GABLAM Fragment %d of %d (%s - %s)' % (fx+1,ftot,rje.iStr(fragstart+1),rje.iStr(fragend+1)))
            seqname = string.join(seqname)
            sequence = fullseq[fragstart-1:fragend]
            return (seqname,sequence)

            #!# Delete when safe ...
            ### ~ [2] ~ Using index or file mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlen = 0
            ## ~ [2a] ~ Get file position for sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.mode() == 'index':
                INDEX = self.INDEX()
                ipos = rje.posFromIndex(seq,INDEX,re_index='^(\S+)\s')
                iline = rje.fileLineFromSeek(INDEX,ipos)
                fpos = string.atol(string.split(rje.chomp(iline[0]))[1])
            else: fpos = seq
            ## ~ [2b] ~ Read sequence name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            SEQFILE = self.SEQFILE()
            SEQFILE.seek(seq)
            name = rje.chomp(SEQFILE.readline())
            if name[:1] != '>':
                self.deBug(name)
                self.errorLog('Given file position that is not name line. May have SeqList objects mixed up?',printerror=False)
                raise ValueError    # Must be Fasta!
            seqname = name[1:]
            sequence = ''; line = rje.chomp(SEQFILE.readline())
            while line and line[:1] != '>':
                prelen = seqlen
                seqlen += len(line)
                if seqlen >= fragstart:
                    if prelen >= fragstart: sequence += line
                sequence += line
                line = rje.chomp(SEQFILE.readline())

            if not case: sequence = sequence.upper()
            #!# Delete when safe ^^^

        except: self.errorLog('%s.getSeqFrag(%s) error' % (self,type(seq))); raise
#########################################################################################################################
    def readSeqType(self,log=True): ### Calculates sequence format from sequences
        '''Calculates sequence format from sequences.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prot = dna = rna = 0
            if self.getBool('Mixed'): tseq = self.seqs()
            else: tseq = self.seqs()[:1]
            ### ~ [1] ~ Read Sequence Types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in tseq:
                (name,sequence) = self.getSeq(seq,format='tuple')
                seqtype = seqType(sequence)
                if seqtype == 'DNA': dna = 1
                elif seqtype == 'Protein': prot = 2
                elif seqtype == 'RNA': rna = 4
            total = dna + rna + prot
            if total == 1: self.str['SeqType'] = 'DNA'
            elif total == 2: self.str['SeqType'] = 'Protein'
            elif total == 4: self.str['SeqType'] = 'RNA'
            elif total == 5: self.str['SeqType'] = 'NA'
            else: self.str['SeqType'] = 'Mixed'
            if log: self.printLog('#TYPE','Sequence type: %s' % self.str['SeqType'])
            return self.str['SeqType']
        except: self.errorLog('Error in "%s" readSeqType' % self.str['Name']); raise
#########################################################################################################################
    ### <3> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.tidy()
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] AutoLoad ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def tidy(self):     ### Shuts open file handles etc. ready for closure
        '''Shuts open file handles etc. ready for closure.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['SEQFILE']: self.obj['SEQFILE'].close(); self.obj['SEQFILE'] = None
            if self.obj['INDEX']: self.obj['INDEX'].close(); self.obj['INDEX'] = None
            return True
        except: self.errorLog('Problem during %s tidy.' % self); return False   # Tidy failed
#########################################################################################################################
    def summarise(self,seqs=None,basename=None):    ### Generates summary statistics for sequences
        '''Generates summary statistics for sequences.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqbase = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if not self.getStrLC('SeqIn'): seqbase = self.getStr('Basefile')
            if not basename: basename = seqbase
            seqdata = {}    # SeqNum, TotLength, MinLength, MaxLength, MeanLength, MedLength, N50Length
            self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SEQUENCE SUMMARY FOR %s ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % basename)
            seqlen = []
            if not seqs: seqs = self.seqs()
            # Total number of sequences
            for seq in seqs:
                self.progLog('\r#SUM','Total number of sequences: %s' % rje.iLen(seqlen),rand=0.01)
                seqlen.append(self.seqLen(seq))
            self.printLog('#SUM','Total number of sequences: %s' % rje.iLen(seqlen))
            if not seqs: return {}
            # Total sequence length
            sumlen = sum(seqlen)
            self.printLog('#SUM','Total length of sequences: %s' % rje.iStr(sumlen))
            seqdata['SeqNum'] = len(seqlen)
            seqdata['TotLength'] = sumlen
            # Min, Max
            seqlen.sort()
            self.printLog('#SUM','Min. length of sequences: %s' % rje.iStr(seqlen[0]))
            self.printLog('#SUM','Max. length of sequences: %s' % rje.iStr(seqlen[-1]))
            seqdata['MinLength'] = seqlen[0]
            seqdata['MaxLength'] = seqlen[-1]
            # Mean & Median sequence lengths
            meanlen = float(sumlen)/len(seqlen)
            meansplit = string.split('%.2f' % meanlen,'.')
            self.printLog('#SUM','Mean length of sequences: %s.%s' % (rje.iStr(meansplit[0]),meansplit[1]))
            seqdata['MeanLength'] = meanlen
            if rje.isOdd(len(seqlen)): median = seqlen[len(seqlen)/2]
            else: median = sum(seqlen[len(seqlen)/2:][:2]) / 2.0
            self.printLog('#SUM','Median length of sequences: %s' % (rje.iStr(median)))
            seqdata['MedLength'] = median
            ## N50 calculation
            n50len = sumlen / 2.0
            n50 = seqlen[0:]
            while n50len > 0 and n50: n50len -= n50.pop(-1)
            if n50:
                self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(n50[-1]))
                seqdata['N50Length'] = n50[-1]
            else:
                self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(seqlen[-1]))
                seqdata['N50Length'] = seqlen[-1]
            return seqdata
        except: self.errorLog('Problem during %s summarise.' % self.prog()); return {}   # Summarise failed
#########################################################################################################################
    ### <4> ### Class Loading Methods                                                                                   #
#########################################################################################################################
    def shuffleSeq(self,replacement=False):   ### Randomly shuffles sequences and saves them to file before reloading as normal
        '''
        Randomly shuffles sequences and saves them to file before reloading as normal.
        >> replacement:bool [False] = Whether to shuffle with or without replacment. (e.g. Use freq versus pure shuffle).
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']: outfile = '%s.shuffled.fas' % rje.baseFile(self.getStr('SeqIn'))
            else: self.setStr({'SeqOut':'None'})
            rje.backup(self,outfile)
            ### ~ [1] Shuffle Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN = open(self.getStr('SeqIn'),'r')
            OUT = open(outfile,'a')
            ## ~ [1a] Count sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0
            iline = IN.readline(); name = None; seq = None
            acc = self.getStr('NewAcc').lower()
            while iline:
                if iline[:1] == '>':
                    if name:
                        if replacement: seq = rje.randomString(len(seq),seq)
                        else: seq = rje.shuffleString(seq)
                        OUT.write('>%s\n%s\n' % (name,seq))
                    sx += 1; self.progLog('\r#SHUF','Shuffling %s sequences' % rje.iStr(sx))
                    name = rje.chomp(iline[1:])
                    seq = ''
                else: seq += rje.chomp(iline)
                iline = IN.readline()
            if name:
                if replacement: seq = rje.randomString(len(seq),seq)
                else: seq = rje.shuffleString(seq)
                OUT.write('>%s\n%s\n' % (name,seq))
            self.printLog('\r#SHUF','Shuffled %s sequences -> %s' % (rje.iStr(sx),outfile))
            IN.close(); OUT.close()
            self.setStr({'SeqIn':outfile})
        except: self.errorLog('Problem during %s shuffleSeq.' % self); return False   # Tidy failed
#########################################################################################################################
    def concatenate(self,outfile=''):   ### Concatenates all sequences into a single sequence and saves to *.cat.fas
        '''Concatenates all sequences into a single sequence and saves to *.cat.fas or outfile.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if outfile == self.getStr('SeqIn'): outfile = ''
            if outfile.lower() in ['','none']: outfile = '%s.cat.fas' % rje.baseFile(self.getStr('SeqIn'))
            rje.backup(self,outfile)
            ### ~ [1] Read in Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN = open(self.getStr('SeqIn'),'r')
            seqx = 0
            name = rje.baseFile(self.getStr('SeqIn'),True)
            seq = ''
            iline = IN.readline()
            while iline:
                if iline[:1] == '>':
                    seqx += 1
                    self.progLog('\r#CAT','Concatenating %s sequences' % rje.iStr(seqx))
                    if seqx > 1: seq += self.getStr('Split')
                else: seq += rje.chomp(iline)
                iline = IN.readline()
            self.printLog('\r#CAT','Concatenated %s sequences: length = %s.' % (rje.iStr(seqx),rje.iLen(seq)))
            IN.close()
            ### ~ [2] Save and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            open(outfile,'a').write('>%s\n%s\n' % (name,seq))
            self.printLog('\r#SAVE','Saved concatenated sequences to %s.' % (outfile))
            self.setStr({'SeqIn':outfile})
            return True
        except: self.errorLog('Problem during %s concatenate.' % self); return False   # Failed
#########################################################################################################################
    def rename(self,keepsprotgene=False):   ### Renames sequences and saves them to file before reloading as normal
        '''
        Renames sequences and saves them to file.
        >> keepsprotgene:bool[False] = Whether to keep SwissProt gene if recognised.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#if self.getStr('NewAcc').lower() in ['','none']: return False
            if self.getStr('NewGene').lower() in ['','none']: self.str['NewGene'] = self.getStr('NewAcc').lower()
            if self.getStr('NewGene').lower() in ['','none']: self.str['NewGene'] = 'seq'
            if self.getStr('SpCode').lower() in ['','none']:
                self.setStr({'SpCode':'UNK'})
                self.printLog('#SPEC','spcode=X not specified. Will use "UNK" unless species read.')
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']: outfile = '%s.new.fas' % rje.baseFile(self.getStr('SeqIn'))
            else: self.setStr({'SeqOut':'None'})
            rje.backup(self,outfile)
            ### ~ [1] Rename Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN = open(self.getStr('SeqIn'),'r')
            OUT = open(outfile,'a')
            ## ~ [1a] Count sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqx = 0
            iline = IN.readline()
            while iline:
                if iline[:1] == '>': seqx += 1; self.progLog('\r#NAME','Counting %s sequences' % rje.iStr(seqx))
                iline = IN.readline()
            IN.seek(0); sx = 0
            iline = IN.readline(); name = None; seq = None
            acc = self.getStr('NewAcc').lower()
            while iline:
                if iline[:1] == '>':
                    if name: OUT.write('>%s\n%s\n' % (name,seq))
                    sx += 1; self.progLog('\r#NAME','Renaming %6s sequences' % rje.iStr(seqx-sx))
                    try:
                        (gene,spcode,accnum) = rje.matchExp('^(\S+)_(\S+)__(\S+)\s*',iline[1:])
                        if gene.upper() != gene or not keepsprotgene: gene = self.getStr('NewGene')
                        if acc in ['','none']: name = '%s_%s__%s %s' % (gene,spcode,accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                        else: name = '%s_%s__%s%s %s %s' % (gene,spcode,self.getStr('NewAcc'),rje.preZero(sx,seqx),accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                    except:
                        accnum = string.split(rje.chomp(iline)[1:])[0]
                        spcode = self.getStr('SpCode').upper()
                        gene = self.getStr('NewGene')
                        if acc in ['','none']: name = '%s_%s__%s %s' % (gene,spcode,accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                        else: name = '%s_%s__%s%s %s' % (gene,spcode,self.getStr('NewAcc'),rje.preZero(sx,seqx),rje.chomp(iline[1:]))
                    seq = ''
                else: seq += rje.chomp(iline)
                iline = IN.readline()
            if name: OUT.write('>%s\n%s\n' % (name,seq))
            self.printLog('\r#NAME','Renamed %s sequences -> %s' % (rje.iStr(sx),outfile))
            IN.close(); OUT.close()
            self.setStr({'SeqIn':outfile})
        except: self.errorLog('Problem during %s rename.' % self); return False   # Tidy failed
#########################################################################################################################
    def sizeSort(self,nodup=True):     ### Sort sequences by size and removes duplicates then saves them to file before reloading as normal
        '''Removes redundancy from sequences and saves them to file before reloading as normal.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']:
                outfile = '%s.sort.fas' % rje.baseFile(self.getStr('SeqIn'))
                if rje.isYounger(outfile,self.getStr('SeqIn')) == outfile and not self.force():
                    self.printLog('#SORT','Will use existing sorted sequences (%s)' % outfile)
                    self.setStr({'SeqIn':outfile}); return
            else: self.setStr({'SeqOut':'None'})
            dupcheck = []; self.dict['Filter']['Duplicate'] = 0
            ### ~ [1] ~ Read into size dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sizedict = {}
            SEQ = self.SEQFILE(); fpos = 0; SEQ.seek(0,2); fend = SEQ.tell(); SEQ.seek(0); sx = 0
            sequence = None; name = None; fprev = 0
            while fpos < fend:
                self.progLog('\r#SIZE','Reading sequences sizes: %.2f%%' % (100.0*fpos/fend),rand=0.2)
                line = rje.chomp(SEQ.readline())
                if line.find('>') == 0:                                     # New Sequence
                    if sequence and name:
                        if nodup and name in dupcheck: self.dict['Filter']['Duplicate'] += 1
                        else:
                            slen = len(sequence)
                            if slen not in sizedict: sizedict[slen] = []
                            sizedict[slen].append(fprev)
                            if nodup: dupcheck.append(name)
                    name = string.split(line[1:])[0]; sequence = ''; fprev = fpos; sx += 1
                else: sequence += line[0:]
                fpos = SEQ.tell()
            if sequence and name:
                if nodup and name in dupcheck: self.dict['Filter']['Duplicate'] += 1
                else:
                    slen = len(sequence)
                    if slen not in sizedict: sizedict[slen] = []
                    sizedict[slen].append(fprev)
                    if nodup: dupcheck.append(name)
            self.printLog('\r#SIZE','Read sizes for %s sequences for sorting.' % rje.iStr(sx))
            if self.dict['Filter']['Duplicate']: self.printLog('\r#DUP','%s duplicate sequences ignored' % (rje.iStr(self.dict['Filter']['Duplicate'])))
            ### ~ [2] ~ Sort and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqx = sx; sx = 0.0
            if nodup: seqx = len(dupcheck)
            rje.backup(self,outfile)
            OUT = open(outfile,'a')
            for slen in rje.sortKeys(sizedict,revsort=True):
                for fpos in sizedict[slen]:
                    (name,seq) = self.getSeq(seq=fpos,format='tuple',mode='file',case=None)
                    OUT.write('>%s\n%s\n' % (name,seq))
                    self.progLog('\r#SORT','Saving size-sorted sequences: %.2f%%' % (sx/seqx),rand=0.1); sx += 100.0
            self.printLog('\r#SORT','Saved %s size-sorted sequences -> %s' % (rje.iStr(seqx),outfile))
            OUT.close(); SEQ.close()
            self.setStr({'SeqIn':outfile})
        except: self.errorLog('Problem during %s sizesort.' % self); return False   # Tidy failed
#########################################################################################################################
    def seqSort(self,nodup=True,sortmode=None):   ### Sort sequences and removes duplicates then saves them to file before reloading as normal
        '''Removes redundancy from sequences and saves them to file before reloading as normal.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not sortmode: sortmode = self.getStrLC('SortSeq')
            else: sortmode = sortmode.lower()
            invert = sortmode[:3] in ['inv','rev']
            if invert: sortmode = sortmode[3:]
            self.debug(sortmode)
            if sortmode[:3] in ['siz','acc','nam','seq','spe','des']:
                sortdesc = {'siz':'size','acc':'accnum','nam':'name','seq':'sequence','spe':'species','des':'description'}[sortmode[:3]]
            else: raise ValueError('Unrecognised sort method: "%s"' % sortmode)
            if rje.matchExp('^seq(\d+)',sortmode):
                sortseqlen = string.atoi(rje.matchExp('^seq(\d+)',sortmode)[0])
                sortdesc = 'first %d positions in sequence' % sortseqlen
            else: sortseqlen = 0
            if invert: sortdesc += '(reverse sorted)'
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']:
                outfile = '%s.sort.fas' % rje.baseFile(self.getStr('SeqIn'))
                rje.backup(self,outfile)
            self.setStr({'SeqOut':'None'})  # Prevent over-writing by later method
            dupcheck = []; self.dict['Filter']['Duplicate'] = 0
            ### ~ [1] ~ Read into size dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sortdict = {}
            SEQ = self.SEQFILE(); fpos = 0; SEQ.seek(0,2); fend = SEQ.tell(); SEQ.seek(0); sx = 0
            sequence = None; name = None; fprev = 0
            while fpos < fend or name:
                self.progLog('\r#SORT','Reading sequences: %.2f%%' % (100.0*fpos/fend),rand=0.2)
                if fpos < fend: line = rje.chomp(SEQ.readline())
                if line.find('>') == 0 or fpos >= fend:   # New Sequence or end of file
                    if sequence and name:
                        if nodup and name in dupcheck: self.dict['Filter']['Duplicate'] += 1
                        else:
                            if sortmode.startswith('siz'): skey = len(sequence)
                            elif sortmode.startswith('seq'):
                                if sortseqlen: skey = sequence[:sortseqlen]
                                else: skey = sequence
                            else:
                                ndata = string.split(name)
                                if sortmode.startswith('des'): skey = string.join(ndata[1:])
                                elif sortmode.startswith('nam'): skey = ndata[0]
                                elif sortmode.startswith('acc'): skey = string.split(ndata[0],'__')[-1]
                                elif sortmode.startswith('spe'): skey = string.split(ndata[0],'_')[1]
                                else: raise ValueError('Unrecognised sort method: "%s"' % sortmode)
                            if skey not in sortdict: sortdict[skey] = []
                            sortdict[skey].append(fprev)
                            if nodup: dupcheck.append(name)
                    if fpos < fend: name = string.split(line[1:])[0]; sequence = ''; fprev = fpos; sx += 1
                    else: name = ''
                else: sequence += line[0:]
                fpos = SEQ.tell()
            self.printLog('\r#SORT','Read %s for %s sequences for sorting.' % (sortdesc,rje.iStr(sx)))
            if self.dict['Filter']['Duplicate']: self.printLog('\r#DUP','%s duplicate sequences ignored' % (rje.iStr(self.dict['Filter']['Duplicate'])))
            ### ~ [2] ~ Sort and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqx = sx; sx = 0.0
            if nodup: seqx = len(dupcheck)
            rje.backup(self,outfile)
            OUT = open(outfile,'a')
            for skey in rje.sortKeys(sortdict,revsort=invert):
                for fpos in sortdict[skey]:
                    (name,seq) = self.getSeq(seq=fpos,format='tuple',mode='file',case=None)
                    OUT.write('>%s\n%s\n' % (name,seq))
                    self.progLog('\r#SORT','Saving %s-sorted sequences: %.2f%%' % (sortdesc,sx/seqx),rand=0.1); sx += 100.0
            self.printLog('\r#SORT','Saved %s %s-sorted sequences -> %s' % (rje.iStr(seqx),sortdesc,outfile))
            OUT.close(); SEQ.close()
            self.setStr({'SeqIn':outfile})  # Additional processing will modify, not over-write, this file name.
        except: self.errorLog('Problem during %s seqsort.' % self); raise
#########################################################################################################################
    def seqNR(self,twopass=False,grepnr=True):  ### Removes redundancy from sequences and saves them to file before reloading as normal
        '''
        Removes redundancy from sequences and saves them to file before reloading as normal.
        >> twopass:bool [False] = Whether to use two-pass NR removal (if sequences not sorted)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']: outfile = '%s.nr.fas' % rje.baseFile(self.getStr('SeqIn'))
            else: self.setStr({'SeqOut':'None'})
            self.dict['Filter']['NR'] = 0
            seqdict = {}    # Dictionary of sequence: name for output if needed
            sequences = ''  # String that gets built for NR check
            ## ~ [0a] Count sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #x#self.deBug(os.popen("grep '>' %s -c" % self.getStr('SeqIn')).readlines())
            grepnr = grepnr and not twopass #and not self.getBool('RevCompNR')
            if grepnr:
                try:
                    seqx = string.atoi(rje.chomp(os.popen("grep -c '>' %s" % self.getStr('SeqIn')).readlines()[0])); grepnr = [True]
                    self.printLog('#GREP','Identified %s sequences using grep' % rje.iStr(seqx))
                except: self.printLog('#GREP','grep failure: will use python NR mode'); grepnr = False
            IN = open(self.getStr('SeqIn'),'r')
            if not grepnr:
                seqx = 0
                iline = IN.readline()
                while iline:
                    if iline[:1] == '>': seqx += 1; self.progLog('\r#NAME','Counting %s sequences' % rje.iStr(seqx))
                    iline = IN.readline()
            ### ~ [1] Forward Pass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN.seek(0); sx = 0.0; goodnames = []
            iline = IN.readline(); name = None; seq = None
            while seqx:
                if iline[:1] == '>' or not iline:
                    if grepnr and name and name not in grepnr:
                        goodrev = goodname = ''
                        foundme = False
                        #self.deBug("grep %s %s -B 1 | grep '>'" % (seq,self.getStr('SeqIn')))
                        glines = os.popen("grep %s %s -B 1 | grep '>'" % (seq,self.getStr('SeqIn'))).readlines()
                        #self.deBug(name)
                        #self.deBug(glines)
                        for line in glines:
                            gname = string.split(rje.chomp(line))[0][1:]
                            if gname == name: foundme = True
                            if not goodname: goodname = gname
                            elif gname not in grepnr:
                                self.printLog('\r#SEQNR','%s redundant with %s' % (gname, goodname))
                                grepnr.append(gname)
                        if not goodname and not goodrev:
                            self.deBug(name); self.deBug(seq); self.deBug(glines)
                        if not foundme:
                            self.errorLog('Problem: %s does not find itself! Will try memory-hungry method.' % name,printerror=False)
                            self.bugPrint("grep -B 1 %s %s" % (seq,self.getStr('SeqIn')))
                            self.deBug(os.popen("grep -B 1 %s %s" % (seq,self.getStr('SeqIn'))).read())
                            IN.close()
                            if self.getStr('SeqOut') == 'None': self.setStr({'SeqOut':'None'})
                            return self.seqNR(grepnr=False)
                        if self.getBool('RevCompNR') and self.nt():
                            revseq = rje_sequence.reverseComplement(seq,rna=self.rna())
                            glines = os.popen("grep -B 1 %s %s | grep '>'" % (revseq,self.getStr('SeqIn'))).readlines()
                            for line in glines:
                                gname = string.split(rje.chomp(line))[0][1:]
                                if not goodrev: goodrev = gname
                                elif gname not in grepnr:
                                    self.printLog('\r#SEQNR','Revcomp %s redundant with %s' % (gname, goodrev))
                                    grepnr.append(gname)
                        if not goodname and not goodrev: self.deBug(glines)
                        elif goodname and not goodrev:
                            if goodname not in goodnames: goodnames.append(goodname)
                        elif goodrev not in goodnames:
                            if goodname not in goodnames: goodnames.append(goodname)
                            self.printLog('\r#SEQNR','Revcomp %s redundant with %s' % (goodrev, goodname))
                            grepnr.append(goodrev)
                        else:
                            self.printLog('\r#SEQNR','Revcomp %s redundant with %s' % (goodname, goodrev))
                            grepnr.append(goodname)
                    elif name and not grepnr:
                        if self.getBool('RevCompNR') and self.nt(): revseq = rje_sequence.reverseComplement(seq,rna=self.rna())
                        else: revseq = ''
                        if seq in seqdict:
                            self.printLog('\r#SEQNR','%s removed: 100%% identical to %s' % (string.split(name)[0],string.split(seqdict[seq])[0]))
                            self.dict['Filter']['NR'] += 1
                        elif revseq and revseq in seqdict:
                            self.printLog('\r#SEQNR','%s removed: reverse complement of %s' % (string.split(name)[0],string.split(seqdict[revseq])[0]))
                            self.dict['Filter']['NR'] += 1
                        elif seq in sequences:
                            matched = rje.matchExp('\s(\S*%s\S*)\s' % seq,sequences)[0]
                            self.printLog('\r#SEQNR','%s removed: contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                            self.dict['Filter']['NR'] += 1
                        elif revseq and revseq in sequences:
                            matched = rje.matchExp('\s(\S*%s\S*)\s' % revseq,sequences)[0]
                            self.printLog('\r#SEQNR','%s removed: reverse complement contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                            self.dict['Filter']['NR'] += 1
                        else:
                            sequences += ' %s ' % seq
                            seqdict[seq] = fullname
                    if twopass: self.progLog('\r#SEQNR','Removing redundancy (Pass I): %.3f%%' % (sx/seqx)); sx += 100.0
                    elif grepnr: self.progLog('\r#SEQNR','Removing redundancy (%s:%s|%s): %.3f%%' % (rje.iStr(int(sx)/100), rje.iLen(goodnames),rje.iLen(grepnr[1:]),sx/seqx)); sx += 100.0
                    else: self.progLog('\r#SEQNR','Removing redundancy (%s-%s): %.3f%%' % (rje.iStr(int(sx)/100), rje.iStr(self.dict['Filter']['NR']),sx/seqx)); sx += 100.0
                    if not iline: break
                    fullname = rje.chomp(iline[1:])
                    name = string.split(fullname)[0]
                    seq = ''
                else: seq += rje.chomp(iline).upper()
                iline = IN.readline()
            IN.close(); 
            if grepnr: self.printLog('\r#SEQNR','Removing redundancy (%s:%s|%s): ready to filter.' % (rje.iStr(int(sx)/100), rje.iLen(goodnames),rje.iLen(grepnr[1:])))
            ### ~ [2] Reverse pass ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = string.split(sequences); seqx = len(seqlist); sx = 0.0
            if twopass:
                goodseq = []
                while seqlist:
                    check = seqlist.pop(0)
                    sequences = string.join(seqlist)
                    if self.getBool('RevCompNR') and self.nt(): revseq = rje_sequence.reverseComplement(check,rna=self.rna())
                    else: revseq = ''
                    if check in sequences:
                        matched = rje.matchExp('\s(\S*%s\S*)\s' % seq,' %s ' % sequences)[0]
                        name = seqdict.pop(check)
                        self.printLog('\r#SEQNR','%s removed: contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                        self.dict['Filter']['NR'] += 1
                    elif revseq and revseq in sequences:
                        name = seqdict.pop(check)
                        matched = rje.matchExp('\s(\S*%s\S*)\s' % revseq,' %s ' % sequences)[0]
                        self.printLog('\r#SEQNR','%s removed: reverse complement contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                        self.dict['Filter']['NR'] += 1
                    else: goodseq.append(check)
                self.progLog('\r#SEQNR','Removing redundancy (Pass II): %.2f%%  ' % (sx/seqx)); sx += 100.0
            else: goodseq = seqlist
            ### ~ [3] Save if filtering done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if grepnr: grepnr = grepnr[1:]
            if self.dict['Filter']['NR']:
                sx = 0.0; seqx = len(goodseq)
                rje.backup(self,outfile)
                OUT = open(outfile,'a')
                while goodseq:
                    seq = goodseq.pop(0)
                    name = seqdict.pop(seq)
                    OUT.write('>%s\n%s\n' % (name,seq))
                    self.progLog('\r#SEQNR','Saving %s NR sequences: %.2f%%' % (rje.iStr(seqx),sx/seqx)); sx += 100.0
                self.printLog('\r#SEQNR','Saved %s NR sequences (%s filtered) -> %s' % (rje.iStr(seqx),rje.iStr(self.dict['Filter']['NR']),outfile))
                OUT.close()
                self.setStr({'SeqIn':outfile})
            elif grepnr:
                open('%s.redundant.txt' % rje.baseFile(outfile),'w').write(string.join(grepnr,'\n'))
                self.printLog('\r#SEQNR','%s redundant sequences flagged for removal' % (rje.iLen(grepnr)))
                self._filterCmd(clear=True)
                self.list['BadSeq'] = grepnr
                self.loadSeq(seqfile=self.getStr('SeqIn'),nodup=True,clearseq=True,mode='file')
                self.filterSeqs()
                self.dict['Filter']['NR'] = self.dict['Filter']['BadSeq']
                self.dict['Filter']['BadSeq'] = 0
                self.list['BadSeq'] = []
                self._filterCmd(self.cmd_list,clear=False)
                self.saveSeq(seqfile=outfile)
                self.printLog('\r#SEQNR','Saved %s NR sequences (%s filtered) -> %s' % (rje.iStr(self.seqNum()),rje.iStr(self.dict['Filter']['NR']),outfile))
                self.setStr({'SeqIn':outfile})
            else: self.printLog('\r#SEQNR','No sequence redundancy found: %s unique sequences.' % rje.iStr(seqx))
        except: self.errorLog('Problem during %s rename.' % self); return False   # Tidy failed
#########################################################################################################################
    def loadSeq(self,seqfile=None,filetype=None,seqtype=None,nodup=True,clearseq=True,mode=None,screen=True):     ### Loads sequences from file
        '''
        Loads sequences from file.
        >> seqfile:str = file name
        >> filetype:str = format of sequence file
        - 'fas' = fasta, 'phy' = phylip, 'aln' = clustalW alignment
        >> seqtype:str = type of sequence in file (dna/protein/rna/mixed)
        >> nodup:Boolean = whether to check for (and remove) duplicate sequences.
        >> clearseq:Boolean = whether to clear existing sequences prior to loading [True]
        >> mode:str [None] = Sequence Mode to use and over-ride self.mode() [None]
        >> screen:bool [True] = Whether to output log messages to the screen.
        '''
        try:### ~ [0] ~ SetUp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not mode: mode = self.mode()
            if clearseq: self.list['Seq'] = []; self.tidy(); self.dict['SeqDict'] = {}
            startx = self.seqNum()
            #!# Note. Sequence alignment has been removed for the current time. Will be reinstated from rje_seq. #!#
            if nodup: self.dict['Filter']['Duplicate'] = 0
            #!# Note. Nodup currently only works for list and file (with index) seqmodes #!#
            if seqtype: self.setStr({'SeqType':seqtype})
            ## ~ [0a] ~ Check File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not seqfile: seqfile = self.getStr('SeqIn')
            if seqfile.lower() in ['','none']:
                if self.i() >= 0:
                    seqfile = rje.choice(text='\nInput filename? (Blank to exit.)')
                    if seqfile == '': sys.exit()
                else:
                    self.errorLog('No file name given: cannot load sequences!')
                    return False
            while not os.path.exists(seqfile):
                if ',' in seqfile or '.' not in seqfile: # Interpret as uniprot extraction list
                    uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
                    uniprot._extractProteinsFromURL(string.split(seqfile,','))
                    self.setStr({'SeqMode':'list','SeqType':'protein'})
                    sx = 0
                    for uentry in uniprot.entries(): self._addSeq(uentry.seqname(),uentry.sequence()); sx += 1
                    if sx:
                        self.printLog('\r#SEQ','%s of %s sequences loaded from Uniprot download.' % (rje.iStr(self.seqNum()),rje.iStr(sx)))
                        return True
                    else: self.printLog('\r#FAIL','%s not recognised as Uniprot accnum.' % seqfile)
                if self.i() >= 0:
                    seqfile = rje.choice(text='Input file "%s" not found. Input filename? (Blank to exit.)' % seqfile)
                    if seqfile == '': sys.exit()
                else:
                    self.errorLog('File %s not found. Cannot load sequences!' % seqfile)
                    return False
            if not rje.matchExp('(\S)',open(seqfile,'r').readline()):
                self.printLog('#ERR','Cannot load sequences from %s. Might be empty?!' % seqfile)
                return False
            self.printLog('#LOAD','Load sequences from %s' % seqfile,log=False,screen=screen)
            ## ~ [0b] ~ Check and report filtering options? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('AutoFilter'):
                for filt in ['Acc','Seq','Spec','DB','Desc']:
                    if self.list['Good%s' % filt]: self.printLog('#FILT','Filter on Good%s: %s %s' % (filt,rje.iLen(self.list['Good%s' % filt]),filt))
                    if self.list['Bad%s' % filt]: self.printLog('#FILT','Filter on Bad%s: %s %s' % (filt,rje.iLen(self.list['Bad%s' % filt]),filt))
            ## ~ [0b] ~ Index file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            makeindex = False
            if self.getBool('SeqIndex') and mode in ['file','index']:
                indexfile = seqfile + '.index'
                if indexfile != seqfile and rje.isYounger(indexfile,seqfile) == indexfile:  # Check index file is newer
                    self.setStr({'SeqDB':seqfile})
                    seqfile = indexfile
                elif indexfile != seqfile:
                    makeindex = True
                    self.setStr({'SeqDictType':'short'})
            #self.deBug('%s: %s' % (mode,makeindex))
            ## ~ [0c] ~ Check file type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            filetype = self.checkFileType(seqfile,filetype)
            if not filetype or filetype == 'unknown': raise ValueError("Unknown sequence input file type")
            ## ~ [0d] ~ Database object and Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setStr({'Name':rje.baseFile(seqfile)})
            if mode == 'db':
                if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
                self.obj['DB'].addEmptyTable(self.getStr('Name'),['name','sequence','fpos'],['name'])

            ### ~ [1] ~ Load sequences into self.list['Seq'] ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SEQ = open(seqfile,'r'); fpos = 0; SEQ.seek(0,2); fend = SEQ.tell(); SEQ.seek(0); sx = 0
            self.obj['SEQFILE'] = SEQ
            logseqfile = seqfile
            if len(logseqfile) > 64: logseqfile = '...%s' % seqfile[-64:]
            ## ~ [1a] ~ Fasta format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if filetype == 'fas':
                sequence = None; name = None; fprev = 0
                if mode == 'index' and not makeindex: fpos = fend
                while fpos < fend:
                    self.progLog('\r#SEQ','Loading seq from %s: %.2f%%' % (logseqfile,(100.0*fpos/fend)),rand=0.2,screen=screen)
                    line = rje.chomp(SEQ.readline())
                    if line.find('>') == 0:                                     # New Sequence
                        if sequence and name: self._addSeq(name,sequence,fprev,makeindex,nodup=nodup)  # Previous Sequence to Add
                        name = line[1:]; sequence = ''; fprev = fpos; sx += 1
                    else: sequence += line[0:]
                    fpos = SEQ.tell()
                if sequence and name: self._addSeq(name,sequence,fprev,makeindex,nodup=nodup)   # Previous Sequence to Add
            ## ~ [1b] ~ Index format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'index':
                if self.mode() == 'index': self.obj['INDEX'] = SEQ; self.obj['SEQFILE'] = None
                else:
                    self.setStr({'SeqDictType':'short'})
                    while fpos < fend:
                        self.progLog('\r#SEQ','Loading seq from %s: %.2f%%' % (logseqfile,(100.0*fpos/fend)),rand=0.01,screen=screen)
                        line = rje.matchExp('^(\S+)\s+(\d+)',rje.chomp(SEQ.readline()))
                        if line:
                            self.list['Seq'].append(string.atol(line[1])); sx += 1
                            self.dict['SeqDict'][line[0]] = self.list['Seq'][-1]
                        fpos = SEQ.tell()
                    self.list['Seq'].sort()
                    SEQ.close(); self.obj['INDEX'] = None
                    self.obj['SEQFILE'] = open(self.getStr('SeqDB'),'r')
            ##  ~ [2b] ~ Phylip Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'phy':
                self.printLog('#SEQ','Phylip format not currently supported. Please use rje_seq to reformat first.')
                raise ValueError
            ##  ~ [2c] ~ Aln Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'aln':
                self.printLog('#SEQ','Alignment format not currently supported. Please use rje_seq to reformat first.')
                raise ValueError
            ## ~ [2d] ~ Fasta Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'mafft':
                self.printLog('#SEQ','MAFFT format not currently supported. Please use rje_seq to reformat first.')
                raise ValueError
            ## ~ [2e] ~ UniProt DAT file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'uniprot':
                self.printLog('#SEQ','UniProt format not currently supported. Please use rje_seq to reformat first.')
                raise ValueError
            ## ~ [2f] ~ UniProt DAT file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'fastacmd':
                self.printLog('#SEQ','Fastacmd format not currently supported. Please use rje_seq to reformat first.')
                raise ValueError
            if mode == 'index' and not makeindex: self.printLog('\r#SEQ','No sequences loaded from %s (Format: %s) - index mode.' % (seqfile,filetype))
            elif startx: self.printLog('\r#SEQ','%s of %s sequences loaded from %s (Format: %s): %s total' % (rje.iStr(self.seqNum()-startx),rje.iStr(sx),seqfile,filetype,rje.iStr(self.seqNum())),screen=screen)
            else: self.printLog('\r#SEQ','%s of %s sequences loaded from %s (Format: %s).' % (rje.iStr(self.seqNum()),rje.iStr(sx),seqfile,filetype),screen=screen)

            ### ~ [3] ~ Create Index for file mode if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if makeindex:
                INDEX = open(indexfile,'w')
                for name in rje.sortKeys(self.dict['SeqDict']): INDEX.write('%s %d\n' % (name,self.dict['SeqDict'][name]))
                INDEX.close()
                self.printLog('#INDEX','Index file %s made' % indexfile,screen=screen)

            ### ~ [4] ~ Summarise ignored duplicate sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if nodup and self.dict['Filter']['Duplicate']: self.printLog('\r#DUP','%s duplicate sequences ignored' % (rje.iStr(self.dict['Filter']['Duplicate'])),screen=screen)
            if 'Sequence unavailable' in self.dict['Filter']: self.printLog('\r#NOSEQ','%s unavailable sequences ignored' % (rje.iStr(self.dict['Filter']['Sequence unavailable'])),screen=screen)
            #!# Add rest of code #!#
            return True
        except: self.errorLog('%s.loadSeq error' % self.prog()); return False
#########################################################################################################################
    def checkFileType(self,seqfile,filetype=None):  ### Checks file format if type given. Returns filetype, or False if wrong.
        '''
        Checks file format if type given. Returns filetype, or False if wrong.
        >> seqfile:str = Sequence file name.
        >> filetype:str [None] = sequence format that the file should have
        << filetype:str = identified file type (from first few lines), or False if desired format not found.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not rje.checkForFile(seqfile): raise IOError
            if filetype and filetype.lower() in ['','none','unknown']: filetype = None
            type_re = [('fas','^>(\S+)'),               # Fasta format
                       ('phy','^\s+(\d+)\s+(\d+)'),     # Phylip
                       ('aln','^(CLUSTAL)'),            # ClustalW alignment
                       ('uniprot','^ID\s+(\S+)'),       # UniProt
                       ('mafft','^##### (atgcfreq)'),   # MAFFT alignment
                       ('fastacmd','^(\S+)\s*$'),       # List of names for fastacmd retrieval - needs seqdb
                       ('index','^(\S+)\s+(\d+)'),      # Sequence position index file - needs seqdb
                       ('unknown','(\S+)')]
            readtype = None
            ### ~ [1] ~ Look for file format type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            firstline = open(seqfile, 'r').readline()
            for (sformat,regexp) in type_re:
                if rje.matchExp(regexp,firstline): readtype = sformat; break
            ## ~ [1a] ~ Check other options if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if filetype and filetype != readtype:
                self.errorLog('"%s" format required but "%s" format read.' % (filetype,readtype),printerror=False)
                return False
            if readtype in ['fastacmd','index'] and not rje.checkForFile(self.getStr('SeqDB')):
                self.errorLog('"%s" format requires SeqDB file - "%s" missing.' % (readtype,self.getStr('SeqDB')),printerror=False)
                return False
            return readtype
        except: self.errorLog('%s.checkFileType error' % self)
#########################################################################################################################
    def _addSeq(self,name,sequence,fpos=-1,makeindex=False,filter=None,nodup=False):   ### Adds sequence to object according to mode
        '''
        Adds sequence to object according to mode.
        - full = Full loading into Sequence Objects.
        - list = Lists of (name,sequence) tuples only.
        - file = List of file positions.
        - index = No loading of sequences. Use index file to find sequences on the fly.
        - db = Store sequence data in database object.
        >> name:str = sequence name line (inc. description)
        >> sequence:str = sequence
        >> fpos:int = position in sequence file [-1]
        >> makeindex:bool = whether needing to make index (therefore create SeqDict)
        >> filter:bool [None] = whether to filter sequences (if None, will use self.getBool['AutoFilter']
        >> nodup:bool [False] = Whether to filter duplicate sequences based on names
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.mode() == 'index' and not makeindex: return   # No need!
            if sequence == 'Sequence unavailable':
                try: self.dict['Filter']['Sequence unavailable'] += 1
                except: self.dict['Filter']['Sequence unavailable'] = 1
                return   # Don't add this!
            if filter == None: filter = self.getBool('AutoFilter')
            if nodup and 'Duplicate' not in self.dict['Filter']: self.dict['Filter']['Duplicate'] = 0
            #gnspacc = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
            ### ~ [1] ~ Pre-Processing filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #if gnspacc and filter:
            #    (gene,spec,acc) = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
            #    if self.list['GoodSpec'] and spec not in self.list['GoodSpec']: self.dict['Filter']['GoodSpec'] += 1; return
            #    if self.list['BadSpec'] and spec in self.list['BadSpec']: self.dict['Filter']['GoodSpec'] += 1; return
            ### ~ [2] ~ Processing-free modes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.mode() == 'file':     # List of file positions.
                fkey = string.split(name)[0]
                if makeindex:
                    if fkey in self.dict['SeqDict']:
                        if nodup: self.dict['Filter']['Duplicate'] += 1; return
                        else: self.errorLog('Sequence name "%s" occurs 2+ times!' % fkey)
                    else: self.dict['SeqDict'][fkey] = fpos
                    self.list['Seq'].append(fpos)
                else: self.list['Seq'].append(fpos)     #!# No duplicate checking if index not used #!#
                return
            elif self.mode() == 'index':    # No loading of sequences. Use index file to find sequences on the fly.
                fkey = string.split(name)[0]
                if fkey in self.dict['SeqDict']:
                    if nodup: self.dict['Filter']['Duplicate'] += 1; return
                    else: self.errorLog('Sequence name "%s" occurs 2+ times!' % fkey)
                else: self.dict['SeqDict'][fkey] = fpos
                return
            ### ~ [3] ~ Sequence Processing and filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# At some point, add processing and filtering based on name and sequence #!#
            ### ~ [4] ~ Full ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.mode() == 'full':   # Full loading into Sequence Objects.
                return self.OLD_addSeq(name,sequence)     # Temporary method
            ### ~ [5] ~ Reduced ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elif self.mode() in ['list','tuple']:     # Lists of (name,sequence) tuples only.
                if nodup and (name,sequence) in self.list['Seq']: self.dict['Filter']['Duplicate'] += 1
                else: self.list['Seq'].append((name,sequence))
            elif self.mode() == 'db':       # Store sequence data in database object.
                db = self.db(self.getStr('Name'))
                entry = {'Name':name,'Sequence':sequence,'FPos':fpos}
                db.addEntry(entry)  #!# Later, add method to expand/extract data
                self.list['Seq'].append(entry)
        except: self.errorLog('%s._addSeq error (%s)' % (self,name))            
#########################################################################################################################
    def OLD_addSeq(self,name,sequence):    ### Adds a new Sequence Object to list                               !TEMP!
        '''
        Adds a new Sequence Object to list.
        >> name:str = sequence name line (inc. description)
        >> sequence:str = sequence
        '''
        try:### ~ [1] ~ Setup new sequence object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newseq = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list,parent=self)
            newseq.info['Name'] = name
            #i# Add uppercase sequence and stores case as tuples in newseq.dict['Case']
            if sequence == 'Sequence unavailable': return   # Don't add this!
            newseq.addSequence(sequence,caselist=self.list['Case'],stripnum=self.opt['ReplaceChar'])    
            newseq.info['Type'] = self.info['Type']
            newseq.extractDetails(gnspacc=self.opt['GeneSpAcc'])
            ### ~ [2] ~ Exclude sequence if appropriate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['DBOnly'] and newseq.info['DBase'].lower() not in string.split(self.info['DBList'].lower(),','):
                self.printLog('\r#REM','Sequence %s excluded as not from given database list.' % newseq.shortName())
                return
            elif not self.opt['UnkSpec'] and newseq.info['SpecCode'] == 'UNK':
                if self.stat['Interactive'] >= 1:
                    newseq.info['SpecCode'] = rje.choice('Enter Species Code for %s. (Blank to Exclude.)' % newseq.info['Name'],default='')
                if newseq.info['SpecCode'] in ['UNK','']:
                    self.printLog('\r#REM','Sequence %s excluded as Species Unknown and unkspec=F.' % newseq.shortName())
                    return
            ### ~ [3] ~ Replace characters in sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Termination * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['ReplaceChar'] and newseq.info['Sequence'][-1:] == '*':
                newseq.info['Sequence'] = newseq.info['Sequence'][:-1]
                #self.printLog('#SEQ','Removed termination signal (*) from end of %s.' % newseq.shortName(),screen=False)
            ## ~ [3c] ~ Bad sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['ReplaceChar']:
                sequence = newseq.info['Sequence']
                for r in range(len(sequence)):
                    if sequence[r] not in self.getAlphabet():
                        if newseq.info['Type'] in ['RNA','DNA']: sequence = rje.strSub(sequence,r,r,'N')
                        elif sequence[r] == 'U':
                            sequence = rje.strSub(sequence,r,r,'C')
                            self.printLog('#SEQ','Replaced assumed selenocysteine %s U%d with C.' % (newseq.shortName(),r+1),screen=False)
                        else: sequence = rje.strSub(sequence,r,r,'X')
                if sequence != newseq.info['Sequence']:
                    if self.list['GoodSpec'] and newseq.info['SpecCode'] not in self.list['GoodSpec'] and newseq.info['Species'] not in self.list['GoodSpec']:
                        remseq = 'GoodSpec'     #!# Don't report replacement
                    elif self.list['BadSpec'] and (newseq.info['SpecCode'] in self.list['BadSpec'] or newseq.info['Species'] in self.list['BadSpec']):
                        remseq = 'BadSpec'  #!# Don't report replacement
                    elif newseq.info['Type'] in ['RNA','DNA']: self.printLog('#SEQ','Replaced non-standard characters in %s with Nss.' % newseq.shortName(),screen=False)
                    else: self.printLog('#SEQ','Replaced non-standard characters in %s with Xs.' % newseq.shortName(),screen=False)
                    newseq.info['Sequence'] = sequence
            ## ~ [3c] ~ NTrim ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nchar = {True:'N',False:'X'}[newseq.dna()]
            if self.stat['NTrim'] and nchar in newseq.info['Sequence']:
                nprop = self.stat['NTrim']; 
                nsplit = string.split(newseq.info['Sequence'],nchar)
                while len(nsplit) > 1:
                    nx = len(nsplit) - 1
                    if float(nx) / len(string.join(['']+nsplit[1:],nchar)) >= nprop: nsplit = nsplit[:1]
                    else: nsplit = [string.join(nsplit[:2],nchar)] + nsplit[2:]
                sequence = string.join(nsplit,nchar)
                if sequence != newseq.info['Sequence']:
                    self.printLog('#NTRIM','Trimmed %d trailing %s-rich characters from %s' % (len(newseq.info['Sequence'])-len(sequence),nchar,newseq.shortName()),screen=False)#self.opt['DeBug'])
                    newseq.info['Sequence'] = sequence
                nsplit = string.split(newseq.info['Sequence'],nchar)
                while len(nsplit) > 1:
                    nx = len(nsplit) - 1
                    if float(nx) / len(string.join(nsplit[:-1]+[''],nchar)) >= nprop: nsplit = nsplit[-1:]
                    else: nsplit = nsplit[:-2] + [string.join(nsplit[-2:],nchar)]
                sequence = string.join(nsplit,nchar)
                if sequence != newseq.info['Sequence']:
                    self.printLog('#NTRIM','Trimmed %d leading %s-rich characters from %s' % (len(newseq.info['Sequence'])-len(sequence),nchar,newseq.shortName()),screen=False)#self.opt['DeBug'])
                    newseq.info['Sequence'] = sequence
            self.seq.append(newseq)
            return newseq
        except:
            self.errorLog('Major error during _addSeq(%s)' % name)
            raise
#########################################################################################################################
    def filterSeqs(self,log=True,screen=None):  # Filters sequences using filter command options
        '''Filters sequences using filter command options.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodseq = []
            #!# Add seqfilter option for switching on sequence-based (not name-based) filtering? Or auto-assess?)
            minlen = max(0,self.getInt('MinLen')); maxlen = max(0,self.getInt('MaxLen'))
            #self.debug('Min: %s; Max: %s' % (minlen,maxlen))
            seqfilter = minlen + maxlen
            ### ~ [1] ~ Simple Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            simplefilter = self.mode() == 'file' and self.dict['SeqDict'] and not seqfilter and self.getStr('SeqDictType') in ['short','name']
            if simplefilter:
                sx = 0.0; stot = len(self.dict['SeqDict'])
                for name in self.dict['SeqDict'].keys()[0:]:
                    if screen: self.progLog('\r#FILT','Filtering sequences %.2f%%' % (sx/stot)); sx += 100
                    ok = True
                    if self.list['GoodSeq'] and string.split(name)[0] not in self.list['GoodSeq']: self.dict['Filter']['GoodSeq'] += 1; ok = False
                    if ok and self.list['BadSeq'] and string.split(name)[0] in self.list['BadSeq']: self.dict['Filter']['BadSeq'] += 1; ok = False
                    if ok and rje.matchExp('^(\S+)_(\S+)__(\S+)',name):
                        (gene,spec,acc) = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
                        if self.list['GoodSpec'] and spec not in self.list['GoodSpec']: self.dict['Filter']['GoodSpec'] += 1; ok = False
                        if self.list['BadSpec'] and spec in self.list['BadSpec']: self.dict['Filter']['BadSpec'] += 1; ok = False
                        if self.list['GoodAcc'] and acc not in self.list['GoodAcc']: self.dict['Filter']['GoodAcc'] += 1; ok = False
                        if self.list['BadAcc'] and acc in self.list['BadAcc']: self.dict['Filter']['BadAcc'] += 1; ok = False
                    elif ok and rje.matchExp('^gi\|(\d+)\|(\S+)\|(\S+)\|',name):
                        (gi,db,acc) = rje.matchExp('^gi\|(\d+)\|(\S+)\|(\S+)\|',name)
                        if self.list['GoodAcc'] and acc not in self.list['GoodAcc']: self.dict['Filter']['GoodAcc'] += 1; ok = False
                        if self.list['BadAcc'] and acc in self.list['BadAcc']: self.dict['Filter']['BadAcc'] += 1; ok = False
                    elif ok and self.list['GoodAcc']: self.errorLog('Wrong format for GoodAcc filter! Switching off. (Try as GoodDesc or GoodSeq?)',printerror=False); self.list['GoodAcc'] = []
                    elif ok and self.list['BadAcc']: self.errorLog('Wrong format for BadAcc filter! Switching off. (Try as GoodDesc or GoodSeq?)',printerror=False); self.list['BadAcc'] = []
                    if ok and (self.list['GoodDesc'] or self.list['BadDesc']):
                        fullname = self.getDictSeq(name,format='tuple')[0]
                        if self.list['GoodDesc']:
                            ok = False
                            for desc in self.list['GoodDesc']:
                                if desc.lower() in fullname.lower(): ok = True#; self.deBug('%s Good!' % desc); break
                        if not ok: self.dict['Filter']['GoodDesc'] += 1
                        elif self.list['BadDesc']:
                            for desc in self.list['BadDesc']:
                                if desc.lower() in fullname.lower(): ok = False; self.dict['Filter']['BadDesc'] += 1; break
                    if ok:
                        goodseq.append(self.dict['SeqDict'][name])
                        #if (self.list['GoodDesc'] or self.list['BadDesc']): self.deBug('%s Good!' % fullname)
                    else: self.dict['SeqDict'].pop(name)
            ### ~ [2] ~ Full Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                sx = 0.0; stot = len(self.seqs())
                for seq in self.seqs():
                    if screen: self.progLog('\r#FILT','Filtering sequences %.2f%%' % (sx/stot)); sx += 100
                    (name,sequence) = self.getSeq(seq,format='tuple')
                    ok = False
                    if minlen and len(sequence) < minlen: self.dict['Filter']['MinLen'] += 1
                    elif maxlen and len(sequence) > maxlen: self.dict['Filter']['MaxLen'] += 1
                    else: ok = True
                    #self.debug('%s -> %s vs %s (%s)' % (seq,self.seqLen(seq),len(sequence),ok))
                    if ok and self.list['GoodSeq'] and string.split(name)[0] not in self.list['GoodSeq']: self.dict['Filter']['GoodSeq'] += 1; ok = False
                    if ok and self.list['BadSeq'] and string.split(name)[0] in self.list['BadSeq']: self.dict['Filter']['BadSeq'] += 1; ok = False
                    if ok and rje.matchExp('^(\S+)_(\S+)__(\S+)',name):
                        (gene,spec,acc) = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
                        if self.list['GoodSpec'] and spec not in self.list['GoodSpec']: self.dict['Filter']['GoodSpec'] += 1; ok = False
                        if self.list['BadSpec'] and spec in self.list['BadSpec']: self.dict['Filter']['BadSpec'] += 1; ok = False
                        if self.list['GoodAcc'] and acc not in self.list['GoodAcc']: self.dict['Filter']['GoodAcc'] += 1; ok = False
                        if self.list['BadAcc'] and acc in self.list['BadAcc']: self.dict['Filter']['BadAcc'] += 1; ok = False
                    elif ok and self.list['GoodAcc']: self.errorLog('Wrong format for GoodAcc filter! Switching off. (Try as GoodDesc or GoodSeq?)',printerror=False); self.list['GoodAcc'] = []
                    elif ok and self.list['BadAcc']: self.errorLog('Wrong format for BadAcc filter! Switching off. (Try as GoodDesc or GoodSeq?)',printerror=False); self.list['BadAcc'] = []
                    if ok: goodseq.append(seq)
                    elif name in self.dict['SeqDict']: self.dict['SeqDict'].pop(name)
                    if not ok: self.bugPrint('Filter %s = %s' % (name,seq))
            ### ~ [3] ~ Report filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#FILT','%s of %s sequences retained.' % (rje.iLen(goodseq),rje.iLen(self.list['Seq'])),log=log,screen=screen)
            filtered = False
            for filtertype in rje.sortKeys(self.dict['Filter']):
                if self.dict['Filter'][filtertype]:
                    self.printLog('\r#FILT','%s sequences filtered on %s' % (rje.iStr(self.dict['Filter'][filtertype]),filtertype))
                    filtered = True
            if filtered:
                if simplefilter:
                    goodseq = self.dict['SeqDict'].values()
                    goodseq.sort()
                self.list['Seq'] = goodseq
                #self.debug(goodseq)
                #self.debug(len(self.seqs()))
                self.obj['Current'] = None
        except: self.errorLog('Major error during filterSeqs()')
#########################################################################################################################
     ### <5> ### Class Saving Methods                                                                                   #
#########################################################################################################################
    def saveSeq(self,seqs=None,seqfile=None,reformat=None,append=None,log=True,screen=None,backup=True):  ### Saves sequences in fasta format
        '''
        Saves sequences in SeqList object in fasta format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> reformat:str = Type of output format (fasta/short/acc/acclist/speclist/peptide/qregion/region)
        >> append:boolean [None] = append, do not overwrite, file. Use self.getBool('Append') if None
        >> log:boolean [True] = Whether to log output
        >> screen:bool [None] = Whether to print log output to screen (None will use log setting)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqs == None: seqs = self.seqs()
            if not seqfile: seqfile = self.getStr('SeqOut')
            if not seqs: return self.errorLog('Cannot output to "%s": No sequences.' % seqfile,printerror=False); return
            if screen == None: screen = log
            if append == None: append = self.getBool('Append')
            if reformat == None: reformat = self.getStr('ReFormat')
            if reformat == 'speclist': speclist = []
            if reformat in ['dna2prot','rna2prot','translate','nt2prot'] and not self.nt():
                self.readSeqType()
                if not self.nt(): self.warnLog('Trying to translate possible non-nucleotide sequence (dna=F).')
            ifile = '%s.index' % rje.baseFile(seqfile)
            if rje.exists(ifile):
                self.warnLog('%s deleted.' % ifile)
                os.unlink(ifile)
            ## ~ [0a] ~ Setup QRegion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            qregion = False
            qstart = 0; qend = -1
            if reformat in ['peptides','qregion','region']:
                qregion = True
                (qstart,qend) = string.split(self.getStr('Region'),',')
                (qstart,qend) = (int(qstart),int(qend))
                self.printLog('#REGION','Reformatting to %s %d -> %d' % (reformat,qstart,qend))
                if qstart > 0: qstart -= 1  # qstart is on a 0<L scale
                #qend -= 1
                # Check alignment if query region being use
                if reformat == 'qregion':
                    qseq = self.getSeq(seqs[0],format='tuple')[1]
                    qlen = len(qseq)
                    for seq in seqs[1:]:
                        if reformat == 'qregion' and len(self.getSeq(seq,format='tuple')[1]) != len(qseq):
                            raise ValueError('Sequences different lengths. Cannot generate qregion output.')
                    # Identify alignment region
                    if qstart < 0: qstart = max(0,qlen + qstart)
                    if qend < 0: qend = qlen + qend + 1
                    ax = -1; ix = -1
                    while ax < qstart and ix < qlen:
                        ix += 1
                        if qseq[ix] != '-': ax += 1
                    qstart = ix
                    while ax < (qend-1) and ix < qlen:
                        ix += 1
                        if qseq[ix] != '-': ax += 1
                    qend = ix + 1
            ## ~ [0b] ~ Open output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outlog = 'output to'
            if append and rje.exists(seqfile):
                outlog = 'appended to'
                #self.printLog('#OUT','Append sequences to "%s"' % seqfile,log=log,screen=screen)
                SEQOUT = open(seqfile,'a')
            else:
                if rje.exists(seqfile):
                    if backup: rje.backup(self,seqfile)
                    outlog = 'output overwriting'
                    #self.printLog('#OUT','Overwrite file: "%s"' % seqfile,log=log,screen=screen)
                #else: self.printLog('#OUT','Create new file: "%s"' % seqfile,log=log,screen=screen)
                SEQOUT = open(seqfile,'w')
            ### ~ [1] ~ Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = len(seqs); fasx = stot
            for seq in seqs:
                if screen: self.progLog('\r#OUT','Sequence output: %.2f%%' % (sx/stot)); sx += 100.0
                (name,sequence) = self.getSeq(seq,format='tuple')
                slen = len(sequence); startchop = 0; endchop = 0
                if qregion:
                    startchop = qstart - sequence[:qstart].count('-')   # Positions lost
                    if qend != -1:
                        endchop = slen - len(sequence[:qend]) - sequence[qend:].count('-')   # Positions lost
                        sequence = sequence[qstart:qend]
                    else: sequence = sequence[qstart:]
                if reformat[:3] in ['fas']: SEQOUT.write('>%s\n%s\n' % (name,sequence))
                elif reformat == 'peptides':
                    sequence = string.replace(sequence,'-','')
                    if sequence: SEQOUT.write('%s\n' % (sequence))
                    else: self.warnLog('No peptide for %s (%d -> %d): 100%% gaps' % (string.split(name)[0],qstart+1,qend))
                elif reformat == 'short': SEQOUT.write('>%s\n%s\n' % (string.split(name)[0],sequence))
                elif reformat.endswith('region'):
                    sstart = startchop + 1
                    send = self.aaLen(seq) - endchop
                    SEQOUT.write('>%s.%s-%s\n%s\n' % (string.split(name)[0],rje.preZero(sstart,slen),rje.preZero(send,slen),sequence))
                elif reformat[:3] in ['acc','spe']:
                    try: (gene,spec,acc) = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
                    except: acc = string.split(name)[0]; gene = 'seq'; spec = 'UNKSP'
                    if reformat in ['acc','accfas']: SEQOUT.write('>%s\n%s\n' % (acc,sequence))
                    elif reformat == 'acclist': SEQOUT.write('%s\n' % (acc))
                    elif reformat == 'speclist':
                        if spec not in speclist: speclist.append(spec)
                elif reformat in ['dna2prot','rna2prot','translate','nt2prot']:
                    #SEQOUT.write('>%s\n%s\n' % (name,rje_sequence.dna2prot(sequence)))
                    pfasta = self.dna2protFasta(name,sequence)
                    fasx = fasx + (pfasta.count('\n')/2) - 1
                    SEQOUT.write(pfasta)
                else:
                    self.errorLog('Cannot recognise reformat type "%s": changing to Fasta' % reformat)
                    reformat = 'fasta'
                    SEQOUT.write('>%s\n%s\n' % (name,sequence))
            ### ~ [3] ~ Close files and tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if reformat == 'speclist': SEQOUT.write(string.join(speclist+[''],'\n'))
            SEQOUT.close()
            self.printLog('\r#OUT','%s Sequences %s %s' % (rje.iStr(fasx),outlog,seqfile),log=log,screen=screen)
            
        except(IOError): self.errorLog("Cannot create %s" % seqfile); raise
        except: self.errorLog("Problem saving sequences to %s" % seqfile); raise
#########################################################################################################################
    def splitSeq(self,basename=None,splitseq=None):   ### Outputs sequence sets into separate files
        '''
        Outputs sequence sets into separate files.
        >> basename:str [None] = basefile for output sequences (basefile or seqin by default).
        >> splitseq:str [None] = type of splitting. (self.getStr('SplitSeq') by default).
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not basename:
                if self.getStrLC('Basefile'): basename = self.getStr('Basefile')
                else: basename = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if not splitseq: splitseq = self.getStrLC('SplitSeq')
            if not splitseq in ['gene','species','spcode','spec']:
                raise ValueError('Cannot divide sequences based on "%s". (Gene/Species only.)' % splitseq)
            ### ~ [1] ~ Output sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = self.seqNum()
            splits = {}     # Dictionary of splitx seqs
            for seq in self.seqs():
                self.progLog('\r#SPLIT','Splitting sequences by %s: %.1f%%' % (splitseq,sx/stot)); sx += 100.0
                if splitseq in ['gene']: splitx = self.seqGene(seq)
                else: splitx = self.seqSpec(seq)
                if splitx not in splits: splits[splitx] = [seq]
                else: splits[splitx].append(seq)
            self.printLog('\r#SPLIT','Split sequences by %s: %s files to generate.' % (splitseq,rje.iLen(splits)))
            for splitx in rje.sortKeys(splits):
                sfile = '%s.%s.fas' % (basename,splitx)
                rje.backup(self,sfile)
                open(sfile,'w').write(self.fasta(splits[splitx]))
                self.printLog('#OUT','%s sequences output to %s.' % (rje.iLen(splits[splitx]),sfile))
        except: self.errorLog("Problem saving %s-split sequences" % splitseq); raise
#########################################################################################################################
    def dna2protFasta(self,name,sequence): ### Returns fasta text for translated DNA sequence, using self.attributes.
        '''Returns fasta text for translated DNA sequence, using self.attributes.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rf = self.getInt('RFTran')
            tranfas = ''
            namesplit = string.split(name)
            ### ~ [1] ~ Translate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rf == 1: rfseq = {1:rje_sequence.dna2prot(sequence)}
            elif rf == 3: rfseq = rje_sequence.threeFrameTranslation(sequence)
            elif rf == 6: rfseq = rje_sequence.sixFrameTranslation(sequence)
            else: raise ValueError('rftran=%d not recognised!' % rf)
            ### ~ [2] ~ Full translation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('MinORF') < 1:   # Return full translated sequences, including STOP *
                if rf == 1: return '>%s\n%s\n' % (name,rfseq[1])
                for frame in rje.sortKeys(rfseq):
                    tranfas += '>%s\n%s\n' % (string.join([namesplit[0]+'.RF%d' % frame]+namesplit[1:]),rfseq[frame])
                return tranfas
            ### ~ [3] ~ Selected ORFs only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            minorf = self.getInt('MinORF'); terminorf = self.getInt('TerMinORF')
            if terminorf < 0: terminorf = minorf
            rforfs = {}; longorf = False
            for frame in rje.sortKeys(rfseq):
                self.progLog('\r#ORF','Sequence RF%s ORFs...   ' % frame)
                osplit = string.split(string.replace(rfseq[frame].upper(),'*','*|'),'|')
                orfs = osplit[:1]
                for orf in osplit[1:-1]:
                    if self.getBool('ORFMet'):  # Trim to Nterminal Met
                        if 'M' in orf: orf = orf[orf.find('M'):]
                        else: continue
                    if len(orf) < (minorf + 1): continue
                    orfs.append(orf); longorf = True    # More than termini - internal ORFs meet length requirement
                if len(osplit) > 1: orfs.append(osplit[-1])
                rforfs[frame] = orfs
            for frame in rje.sortKeys(rfseq):
                self.progLog('\r#ORF','Sequence ORF Fasta...   ')
                orfs = rforfs[frame]
                if longorf:
                    if len(orfs[-1]) < minorf: orfs = orfs[:-1]
                    if len(orfs[0]) < (minorf + 1): orfs = orfs[1:]
                else:
                    if len(orfs[-1]) < terminorf: orfs = orfs[:-1]
                    if orfs and len(orfs[0]) < (terminorf + 1): orfs = orfs[1:]
                # Remaining ORFs should meet length requirements
                for i in range(len(orfs)):
                    tranfas += '>%s\n%s\n' % (string.join([namesplit[0]+'.RF%d.ORF%d' % (frame,i+1)]+namesplit[1:]+['Length=%d' % len(orfs[i])]),orfs[i])
            return tranfas
        except: self.errorLog("Problem with dna2protFasta()"); raise
#########################################################################################################################
    def OLDdna2protFasta(self,name,sequence): ### Returns fasta text for translated DNA sequence, using self.attributes.
        '''Returns fasta text for translated DNA sequence, using self.attributes.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rf = self.getInt('RFTran')
            tranfas = ''
            namesplit = string.split(name)
            ### ~ [1] ~ Translate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rf == 1: rfseq = {1:rje_sequence.dna2prot(sequence)}
            elif rf == 3: rfseq = rje_sequence.threeFrameTranslation(sequence)
            elif rf == 6: rfseq = rje_sequence.sixFrameTranslation(sequence)
            else: raise ValueError('rftran=%d not recognised!' % rf)
            ### ~ [2] ~ Full translation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('MinORF') < 1:   # Return full translated sequences, including STOP *
                if rf == 1: return '>%s\n%s\n' % (name,rfseq[1])
                for frame in rje.sortKeys(rfseq):
                    tranfas += '>%s\n%s\n' % (string.join([namesplit[0]+'.RF%d' % frame]+namesplit[1:]),rfseq[frame])
                return tranfas
            ### ~ [3] ~ Selected ORFs only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            minorf = self.getInt('MinORF'); terminorf = self.getInt('TerMinORF')
            if terminorf < 0: terminorf = minorf
            rforfs = {}; longorf = False
            for frame in rje.sortKeys(rfseq):
                orfs = string.split(string.replace(rfseq[frame].upper(),'*','*|'),'|')
                if self.getBool('ORFMet'):  # Trim to Nterminal Met
                    for i in range(len(orfs)-1,0,-1):
                        if 'M' in orfs[i]: orfs[i] = orfs[i][orfs[i].find('M'):]
                        else: orfs[i] = ''
                for i in range(len(orfs)-2,0,-1):
                    if len(orfs[i]) < (minorf + 1): orfs.pop(i)
                if len(orfs) > 2: longorf = True    # More than termini - internal ORFs meet length requirement
                rforfs[frame] = orfs
            for frame in rje.sortKeys(rfseq):
                orfs = rforfs[frame]
                if longorf:
                    if len(orfs[-1]) < minorf: orfs = orfs[:-1]
                    if len(orfs[0]) < (minorf + 1): orfs = orfs[1:]
                else:
                    if len(orfs[-1]) < terminorf: orfs = orfs[:-1]
                    if orfs and len(orfs[0]) < (terminorf + 1): orfs = orfs[1:]
                # Remaining ORFs should meet length requirements
                for i in range(len(orfs)):
                    tranfas += '>%s\n%s\n' % (string.join([namesplit[0]+'.RF%d.ORF%d' % (frame,i+1)]+namesplit[1:]+['Length=%d' % len(orfs[i])]),orfs[i])
            return tranfas
        except: self.errorLog("Problem with dna2protFasta()"); raise
#########################################################################################################################
    def fasta(self,seqs=[]):    ### Returns text of sequences in fasta format
        '''Returns text of sequences in fasta format.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not type(seqs) == list: seqs = [seqs]
            if not seqs: seqs = self.seqs()
            fastxt = ''
            for seq in seqs: fastxt += '>%s\n%s\n' % self.getSeq(seq,format='tuple')
            return fastxt
        except: self.errorLog("Problem with fasta()"); raise
#########################################################################################################################
    def sampler(self):   ### Randomly samples sequences and outputs to files.
        '''Randomly samples sequences and outputs to files.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['Sampler'][1] = int(self.list['Sampler'][1] + 0.5)
            self.list['Sampler'][0] = max(self.list['Sampler'][0],0.0)
            if self.list['Sampler'][0] <= 0.0: return False
            if self.list['Sampler'][0] > self.seqNum():
                self.warnLog('Cannot sample %s sequences: only %s sequences in total!' % (rje.iStr(self.list['Sampler'][0]), rje.iStr(self.seqNum())))
                return False
            elif self.list['Sampler'][0] < 1:
                sperc = self.list['Sampler'][0] * 100.0
                self.list['Sampler'][0] = int(self.list['Sampler'][0] * self.seqNum() + 0.5)
                self.printLog('#SAMPLE','Sampling %d x %s (%.1f%%) sequences.' % (self.list['Sampler'][1],rje.iStr(self.list['Sampler'][0]),sperc))
            else:
                self.list['Sampler'][0] = int(self.list['Sampler'][0] + 0.5)
                self.printLog('#SAMPLE','Sampling %d x %s sequences.' % (self.list['Sampler'][1],rje.iStr(self.list['Sampler'][0])))
            ## ~ [0a] ~ SeqOut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('SeqOut'): seqout = rje.baseFile(self.getStr('SeqOut'))
            else: seqout = '%s.n%d' % (rje.baseFile(self.getStr('SeqIn')),self.list['Sampler'][0])
            ### ~ [1] ~ Sample and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for x in range(self.list['Sampler'][1]):
                if self.list['Sampler'][1] > 1: rfile = '%s.r%s.fas' % (seqout,rje.preZero(x+1,self.list['Sampler'][1]))
                else: rfile = '%s.fas' % seqout
                rseqs = rje.randomList(self.seqs())
                self.saveSeq(seqs=rseqs[:self.list['Sampler'][0]],seqfile=rfile)
            if self.list['Sampler'][1] > 1:
                self.printLog('#SAMPLE','%s file(s) of %s sequences output to %s.r*.fas' % (self.list['Sampler'][1],rje.iStr(self.list['Sampler'][0]),seqout))
            else: self.printLog('#SAMPLE','%s sequences output to %s.' % (rje.iStr(self.list['Sampler'][0]),rfile))
        except: self.errorLog("Problem with SeqList.sampler()"); raise
#########################################################################################################################
### End of SECTION II: SeqList Class                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def seqType(sequence):  ### Returns (and possible guesses) Sequence Type - Protein/DNA/RNA/NA
    '''
    Returns (and possible guesses) Sequence Type
    - Protein if non-ATGCUN
    - DNA if ATGCN only
    - RNA if AUGCN only
    - NA if AUTGCN only
    '''
    #!# Look for B, J etc. as unrecognised? #!#
    if re.search('[DEFHIKLMPQRSVWY]',sequence.upper()): return 'Protein'
    elif re.search('U',sequence.upper()): return 'RNA'
    else: return 'DNA'
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
    try: SeqList(mainlog,cmd_list).run()

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
