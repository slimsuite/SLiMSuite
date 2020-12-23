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
Module:       rje_jellyfish
Description:  RJE Jellyfish processing module
Version:      0.0.0
Last Edit:    01/10/19
Copyright (C) 2019  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

    The `regions=FILE` contains a set of regions to focus on for generating a kmer frequency table. This should contain
    `Locus` (or `Contig`), `Start` and `End` if a delimited file, or can be a GFF3 file. If `Strand` is given, it will be
    used. Otherwise, it will be added and reverse strand identified as `Start` > `End`. BUSCO full results files will
    also be recognised.

    The sampling frequency is set by `sampling=X`:

    * `random` : Sample random kmers until maxsample=X or walltime=X is met
    * `1/N` : Sample every N bases (N must be an integer)
    * `ends` : Sample each end of the regions given by regions=FILE
    * `full` : Sample every position in input/regions


CompJF Mode:
    This samples kmers from the `seqin=FASFILE` sequences (or given `regions=FILE` regions) in each jellyfish file given
    by `compjf=FILES`, where `FILES` are `*.jf` files generated using the same `k=X` setting.

    NOTE: Jellyfish must be run with reverse complementarity switched on.


Commandline:
    ### ~ Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Source DNA sequences for Jellyfish kmer analysis in fasta format (basename for jellyfish files) []
    jffile=FILE     : Jellyfish file that is associated with seqin=FILE. [$SEQIN.kXX.jf]
    k=INT           : Kmer length used in analysis [21]
    regions=FILE    : Delimited file with Locus, Start, End for subset of input to analyse []
    gfftype=LIST    : List of feature types to use if regions=FILE is a GFF3 file. (All if empty.) []
    jfdir=PATH      : Directory for jellyfish files ("source" for source directory of seqin=FASFILE) [./]
    sampling=X      : Kmer sampling frequency (random,1/N,ends,full) [full]
    uniquekmer=T/F  : Whether to restrict output to unique kmers [False]
    maxsample=X     : Maximum number of kmers to sample (if > 0) [10000]
    walltime=X      : Sampling walltime - will quit execution after X hours [0.0]
    ### ~ Frequency comparison options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    compjf=FILES    : List of jellyfish files (with same k=X) for pulling out frequencies []
    compnames=LIST  : List of field names to be used for files (must be same length as compfj=FILES) to use []
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
import rje, rje_obj, rje_seqlist, rje_db
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
    # [ ] : Add options to toggle output of kmer, locus and position (outkmer=T/F, outloc=T/F, outpos=T/F)
    # [ ] : Add sampling option to give a list of pre-defined kmers (sampling=FILE)
    # [ ] : Swap to use KMC (and rename rje_kmer!)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('RJE_Jellyfish', '0.0.0', 'September 2019', '2019')
    description = 'RJE Jellyfish processing module'
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
### SECTION II: Jellyfish Class                                                                                         #
#########################################################################################################################
class Jellyfish(rje_obj.RJE_Object):
    '''
    Jellyfish Class. Author: Rich Edwards (2019).

    Str:str
    - JFDir=PATH      : Directory for jellyfish files ("source" for source directory of seqin=FASFILE) [./]
    - JFFile=FILE     : Jellyfish file that is associated with seqin=FILE. [$SEQIN.jf]
    - Regions=FILE    : Delimited file with Locus, Start, End for subset of input to analyse []
    - Sampling=X      : Kmer sampling frequency (random,1/N,ends,full) [full]

    Bool:boolean
    - UniqueKmer=T/F : Whether to restrict output to unique kmers [False]

    Int:integer
    - K=INT           : Kmer length used in analysis [21]
    - MaxSample=X     : Maximum number of kmers to sample (if > 0) [10000]

    Num:float
    - WallTime=X      : Sampling walltime - will quit execution after X hours [0.0]

    File:file handles with matching str filenames

    List:list
    - CompJF=FILES    : List of jellyfish files (with same k=X) for pulling out frequencies []
    - CompNames=LIST  : List of field names to be used for files (must be same length as compfj=FILES) to use []
    - GFFType=LIST    : List of feature types to use if regions=FILE is a GFF3 file. (All if empty.) []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['JFDir','JFFile','Regions','Sampling']
        self.boollist = ['UniqueKmer']
        self.intlist = ['K','MaxSample']
        self.numlist = ['WallTime']
        self.filelist = []
        self.listlist = ['CompJF','CompNames','GFFType']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({'JFDir':rje.makePath('./'),'Regions':'None','Sampling':'full'})
        self.setBool({})
        self.setInt({'K':21,'MaxSample':10000})
        self.setNum({})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
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
                self._cmdReadList(cmd,'str',['Sampling'])   # Normal strings
                self._cmdReadList(cmd,'path',['JFDir'])  # String representing directory path
                self._cmdReadList(cmd,'file',['JFFile','Regions'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                self._cmdReadList(cmd,'bool',['UniqueKmer'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['K','MaxSample'])   # Integers
                self._cmdReadList(cmd,'float',['WallTime']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CompNames','GFFType'])  # List of strings (split on commas or file lines)
                #self._cmdReadList(cmd,'clist',['Att']) # Comma separated list as a *string* (self.str)
                self._cmdReadList(cmd,'glist',['CompJF']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.compJF()
        except ValueError: raise
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def k(self): return self.getInt('K')
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''
        Main class setup method.

        First, the input sequences given by `seqin=FILE` are loaded and mapped onto short names and accession numbers for
        matching with loaded `Locus` information from `regions=FILE`.

        Next, regions for kmer sampling are (optionally) loaded from `regions=FILE`.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            ## ~ [1a] Input sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['SeqList'] = rje_seqlist.SeqList(self.log,['dna=T']+self.cmd_list+['autoload=T','seqmode=file'])
            #i# Dictionary matching Locus via short name or accnum
            self.obj['SeqList'].makeSeqNameDic('loci')
            #i# Jellyfish file
            if not self.getStrLC('JFFile'): self.setStr({'JFFile':'%s.k%d.jf' % (self.obj['SeqList'].getStr('SeqIn'),self.k())})
            ## ~ [1b] Regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('Regions'):
                if not rje.exists(self.getStr('Regions')): raise IOError('Cannot find regions=%s' % self.getStr('Regions'))
                #i# Check type and load Regions
                firstline = open(self.getStr('Regions'),'r').readline()
                if firstline.startswith('# BUSCO'):
                    headers = string.split('BuscoID	Status	Locus	Start	End	Score	Length')
                    regdb = db.addTable(self.getStr('Regions'),mainkeys='auto',datakeys='All',headers=headers,ignore=['#'],lists=False,name='regions',expect=True)
                    regdb.dropEntriesDirect('Status',['Complete'],inverse=True)
                else:
                    regdb = db.addTable(self.getStr('Regions'),mainkeys='auto',datakeys='All',headers=headers,ignore=['#'],lists=False,name='regions',expect=True)
                if not regdb.entryNum(): raise ValueError('Failed to load any entries from regions=%s' % self.getStr('Regions'))
                if 'Locus' not in regdb.fields():
                    if 'Contig' in regdb.fields(): regdb.renameField('Contig','Locus')
                regdb.dataFormat({'Start':'int','End':'int'})
                if 'Strand' not in regdb.fields():
                    regdb.addField('Strand',evalue='+')
                #i# Make sure Start < End and record Strand. (Not used at present as kmers symmetrical)
                for entry in regdb.entries():
                    if entry['Start'] > entry['End']:
                        entry['Strand'] = '-'
                        (entry['Start'],entry['End']) = (entry['End'],entry['Start'])
                regdb.newKey(['Locus','Start','End'])
            ## ~ [1c] JF files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for jf in self.list['CompJF']:
                if not rje.exists(jf): raise IOError('CompJF Jellyfish file "%s" not found!' % jf)

            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
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
#########################################################################################################################
    ### <3> ### Additional Class Methods                                                                                #
#########################################################################################################################
    def compJF(self):   ### Jellyfish Kmer Frequency comparison analysis
        '''
        Function:
        The function of this module will be added here.

        The `regions=FILE` contains a set of regions to focus on for generating a kmer frequency table. This should contain
        `Locus` (or `Contig`), `Start` and `End` if a delimited file, or can be a GFF3 file. If `Strand` is given, it will be
        used. Otherwise, it will be added and reverse strand identified as `Start` > `End`. BUSCO full results files will
        also be recognised.

        The sampling frequency is set by `sampling=X`:

        * `random` : Sample random kmers until maxsample=X or walltime=X is met
        * `1/N` : Sample every N bases (N must be an integer)
        * `ends` : Sample each end of the regions given by regions=FILE
        * `full` : Sample every position in input/regions


        CompJF Mode:
        This samples kmers from the `seqin=FASFILE` sequences (or given `regions=FILE` regions) in each jellyfish file given
        by `compjf=FILES`, where `FILES` are `*.jf` files generated using the same `k=X` setting.

        NOTE: Jellyfish must be run with reverse complementarity switched on.


        Commandline:
        ### ~ Input/Output Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        seqin=FASFILE   : Source DNA sequences for Jellyfish kmer analysis in fasta format (basename for jellyfish files) []
        k=INT           : Kmer length used in analysis [21]
        regions=FILE    : Delimited file with Locus, Start, End for subset of input to analyse []
        gfftype=LIST    : List of feature types to use if regions=FILE is a GFF3 file. (All if empty.) []
        jfdir=PATH      : Directory for jellyfish files ("source" for source directory of seqin=FASFILE) [./]
        sampling=X      : Kmer sampling frequency (random,1/N,ends,full) [full]
        uniquekmer=T/F  : Whether to restrict output to unique kmers [False]
        maxsample=X     : Maximum number of kmers to sample (if > 0) [10000]
        walltime=X      : Sampling walltime - will quit execution after X hours [0.0]
        ### ~ Frequency comparison options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        compjf=FILES    : List of jellyfish files (with same k=X) for pulling out frequencies []
        compnames=LIST  : List of field names to be used for files (must be same length as compfj=FILES) to use []
        ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Set up objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Sequences and regions are loaded in setup
            seqlist = self.obj['SeqList']
            seqdict = seqlist.seqNameDic()
            regdb = self.db('regions')
            #i# If no regions loaded, make the full sequences into regions!
            #!# Add an endshift value for sampling=ends when no regions given? Or sampling=POSLIST (with -ves allowed)
            if not regdb:
                regdb = self.db().addEmptyTable('regions',['Locus','Start','End','Strand'],['Locus','Start','End'],log=self.debugging())
                for seq in seqlist.seqs():
                    entry = {'Locus':seqlist.shortName(seq),'Start':1,'End':seqlist.seqLen(seq),'Strand':'+'}
                    regdb.addEntry(entry)
            ## ~ [1b] Set up variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# maxsample and walltime are set at the commandline
            samplex = 0             # Number of kmers sampled
            samplet = time.time()   # Time at start of sampling
            self.list['Sampled'] = []              # List of processed kmers. Only used if uniquekmer=T.
            k = self.k()
            uniquek = self.getBool('UniqueKmer')
            ## ~ [1c] Set up output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.list['CompJF']:
                if not rje.exists(self.getStr('JFFile')): raise IOError('Jellyfish file "%s" not found!' % self.getStr('JFFile'))
                self.list['CompJF'] = [self.getStr('JFFile')]
                self.list['CompNames'] = ['N']
            if not self.list['CompNames']: self.list['CompNames'] = self.list['CompJF']
            if not len(self.list['CompNames']) == len(self.list['CompJF']):
                raise ValueError('%d CompJF files but %d CompNames!' % (len(self.list['CompJF']), len(self.list['CompNames'])))
            #i# Setup output file with kmer,NAMES = `basefile.compjf.tsv`
            outhead = ['#','kmer'] + self.list['CompNames']
            outfile = '%s.compjf.tsv' % self.baseFile()
            rje.backup(self,outfile)
            outmode = {True:'a',False:'w'}[self.getBool('Append')]
            OUT = open(outfile,outmode)
            if not rje.exists(outfile) or not self.getBool('Append'):
                OUT.write('%s\n' % string.join(outhead,'\t'))
            ## ~ [1d] Set up sampling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Random needs to know total kmer number available
            totalk = 0
            for entry in regdb.entries():
                totalk += max(0,(entry['End']-entry['Start']+1)-k+1)
            if not totalk: raise ValueError('No possible k=%d kmers read from %s' % (k,seqlist.getStr('SeqIn')))
            if self.getStrLC('Sampling') == 'random' and self.getInt('MaxSample') < 1 and self.getInt('WallTime') <= 0:
                raise ValueError('Cannot use sampling=random without maxsample=X and/or walltime=X above zero!')
            current = None  # This is either a sequence or a regdb entry
            kstart = 0      # This is the position of the last kmer start (1 to L)

            ### ~ [2] Sample Kmers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            finished = True
            rx = 0.0; rtot = regdb.entryNum()
            sampling = self.getStrLC('Sampling')
            if self.getStrLC('Sampling') == 'full' or not self.getStrLC('Sampling'): sampling = '1/1'

            # * `random` : Sample random kmers until maxsample=X or walltime=X is met
            # * `1/N` : Sample every N bases (N must be an integer)
            # * `ends` : Sample each end of the regions given by regions=FILE
            # * `full` : Sample every position in input/regions
            #!# Add forking
            # Optionally Check kmer against list, then extract from compjf files

            ## ~ [2a] Random sampling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if sampling == 'random': raise ValueError('sampling=random not yet implemented!')
            #            while samplex < self.getInt('MaxSample') or self.getInt('MaxSample') < 1:

            ## ~ [2b] 1/N sampling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif sampling.startswith('1/'):
                try:
                    N = string.atoi(sampling[2:])
                except:
                    raise ValueError('sampling=%s - 1/N must be integer!' % sampling)
                if N == 1: sampling = 'Full'
                for entry in regdb.entries():
                    stepx = 100.0 / (entry['End'] - entry['Start'] + 1)
                    px = 0.0
                    kpos = entry['Start'] - 1
                    while kpos < entry['End']:
                        if self.v() < 2: self.progLog('\r#KMERS','%s-sampling %s kmers: %.2f%%' % (sampling,rje.iStr(samplex),(rx+px)/rtot)); px += (N * stepx)
                        sequence = seqlist.seqSequence(seqdict[entry['Locus']])
                        kmer = sequence[kpos:][:k]
                        if len(kmer) == k:
                            samplex += self.sampleOut(kmer,samplex,OUT,forked=False)
                            if self.checkFinish(samplex,samplet): finished = False; break
                        kpos += N
                    rx += 100.0
                if uniquek:
                    self.printLog('\r#KMERS','%s-sampling %s unique kmers from %s regions completed' % (sampling,rje.iStr(samplex),rje.iStr(rtot)))
                else:
                    self.printLog('\r#KMERS','%s-sampling %s kmers from %s regions completed' % (sampling,rje.iStr(samplex),rje.iStr(rtot)))

            ## ~ [2c] end sampling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif sampling == 'ends':
                for entry in regdb.entries():
                    if self.v() < 2: self.progLog('\r#KMERS','End-sampling %s kmers: %.2f%%' % (rje.iStr(samplex),rx/rtot)); rx += 100.0
                    sequence = seqlist.seqSequence(seqdict[entry['Locus']])
                    # Start
                    kmer = sequence[entry['Start']-1:][:k]
                    if len(kmer) == k:
                        samplex += self.sampleOut(kmer,samplex,OUT,forked=False)
                        if self.checkFinish(samplex,samplet): finished = False; break
                    # End
                    kmer = sequence[:entry['End']][-k:]
                    if len(kmer) == k:
                        samplex += self.sampleOut(kmer,samplex,OUT,forked=False)
                        if self.checkFinish(samplex,samplet): finished = False; break
                if uniquek:
                    self.printLog('\r#KMERS','End-sampling %s unique kmers from %s regions completed' % (rje.iStr(samplex),rje.iStr(rtot)))
                else:
                    self.printLog('\r#KMERS','End-sampling %s kmers from %s regions completed' % (rje.iStr(samplex),rje.iStr(rtot)))

            ## ~ [2d] unknown sampling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                raise ValueError('sampling=%s not recognised!' % self.getStrLC('Sampling'))


            #!# Add uniquekmer=T/F : Whether to restrict output to unique kmers [False]

            # Update output file

            # Check against walltime and maxsample


            #i# Returns True if sampling completed, False if limit reached, or raises error if error
            self.printLog('#OUT','%s sampled kmers output to %s' % (rje.iStr(samplex),outfile))
            OUT.close()
            return finished
        except: self.errorLog('%s.compJF error' % self.prog()); raise
#########################################################################################################################
    def sampleOut(self,kmer,samplex,OUT,forked=False):  ### Samples kmer from self.list['JFFiles'] and saves to OUT
        '''
        Samples kmer from self.list['JFFiles'] and saves to OUT. Returns 1 if kmer sampled, or 0 if not.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('UniqueKmer') and kmer in self.list['Sampled']: return 0
            out = ['%d' % (samplex+1),kmer]
            ksample = self.sampleKmer(kmer,forked)
            for kname in self.list['CompNames']: out.append(ksample[kname])
            OUT.write('%s\n' % string.join(out,'\t'))
            if self.getBool('UniqueKmer'): self.list['Sampled'].append(kmer)
            return 1
        except ValueError: raise
        except: self.errorLog('%s.sampleOut(%s) error' % (self.prog(),kmer)); raise
#########################################################################################################################
    def sampleKmer(self,kmer,forked=False):  ### Samples kmer from self.list['JFFiles'] and returns dictionary
        '''Samples kmer from self.list['JFFiles'] and returns dictionary.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ksample = {'kmer':kmer}
            knames = self.list['CompNames'][0:]
            ### ~ [2] Sample files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for jf in self.list['CompJF']:
                jcmd = 'jellyfish query %s %s' % (jf,kmer)
                if self.v() > 1 and not forked: self.printLog('#JF',jcmd)
                jdata = string.split(os.popen(jcmd).readline())
                try:
                    if len(jdata[0]) != len(kmer): raise ValueError()
                    #kcount = string.atoi(jdata[1])
                except:
                    raise ValueError('Problem with: %s' % jcmd)
                #ksample[knames.pop(0)] = kcount
                ksample[knames.pop(0)] = jdata[1]
            return ksample
        except: self.errorLog('%s.sampleKmer(%s) error' % (self.prog(),kmer)); raise
#########################################################################################################################
    def checkFinish(self,samplex,samplet):  # Return whether to stop sampling now
        maxnum = (samplex >= self.getInt('MaxSample') and self.getInt('MaxSample') > 0)
        if maxnum:
            self.printLog('#MAX','Kmer sampling limit maxsample=%d reached!' % self.getInt('MaxSample'))
            return True
        maxwall = (self.getNum('WallTime') > 0 and (time.time() - samplet) < (self.getNum('WallTime')*3600))
        if maxwall:
            self.printLog('#WALL','Walltime (%.2f hr) reached!' % self.getNum('WallTime'))
            return True
        return False
#########################################################################################################################
### End of SECTION II: Jellyfish Class                                                                                  #
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
    try: Jellyfish(mainlog,['basefile=jellyfish','append=F']+cmd_list).run()

    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.endLog(info)
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
