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
Module:       SAMPhaser
Description:  Diploid chromosome phasing from SAMTools Pileup format.
Version:      0.1.0
Last Edit:    30/08/16
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input genome to phase variants in []
    pileup=FILE     : Pileup file of reads against input genome. []
    basefile=FILE   : Root of output file names (same as Pileup input file by default) []
    mincut=X        : Minimum read count for minor allele (proportion if <1) [0.25]
    absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [5]
    ### ~ Processing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    snperr=X        : Probability of an incorrect (biallelic) SNP call for individual read nucleotides [0.05]
    snpcalc=X       : Max number of SNPs to use for read probability calculations (fewer = quicker) [10]
    trackprob=X     : Min probability for assigning a read/SNP to Track A/B [0.95]
    minsnp=X        : Min number of SNPs per phased haplotype block [5]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, string, sys, time, math
slimsuitepath = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../')) + os.path.sep
sys.path.append(os.path.join(slimsuitepath,'libraries/'))
sys.path.append(os.path.join(slimsuitepath,'tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_obj, rje_db, rje_samtools, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Updated SAMPhaser to be more memory efficient.
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
    (program, version, last_edit, copy_right) = ('SAMPhaser', '0.1.0', 'August 2016', '2016')
    description = 'Diploid chromosome phasing from SAMTools Pileup format'
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
### SECTION II: SAMPhaser Class                                                                                         #
#########################################################################################################################
class SAMPhaser(rje_obj.RJE_Object):
    '''
    SAMPhaser Class. Author: Rich Edwards (2015).

    Str:str
    - Pileup=FILE     : Pileup file of reads against input genome. []
    - Seqin=FASFILE   : Input genome to phase variants in []

    Bool:boolean

    Int:integer
    - AbsMinCut=X     : Absolute minimum read count for minor allele (used if mincut<1) [3]
    - MinSNP=X        : Min number of SNPs per phased haplotype block [5]
    - SNPCalc=X       : Max number of SNPs to use for read probability calculations (fewer = quicker) [10]

    Num:float
    - MinCut=X        : Minimum read count for minor allele (proportion if <1) [0.25]
    - SNPErr=X        : Probability of an incorrect (biallelic) SNP call for individual read nucleotides [0.05]
    - TrackProb=X     : Min probability for assigning a read/SNP to Track A/B [0.95]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Pileup','SeqIn']
        self.boollist = []
        self.intlist = ['AbsMinCut','MinSNP','SNPCalc']
        self.numlist = ['MinCut','SNPErr','TrackProb']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({})
        self.setInt({'AbsMinCut':5,'MinSNP':5,'SNPCalc':10})
        self.setNum({'MinCut':0.25,'SNPErr':0.05,'TrackProb':0.95})
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
                #self._cmdReadList(cmd,'str',['Att'])   # Normal strings
                #self._cmdReadList(cmd,'path',['Att'])  # String representing directory path 
                self._cmdReadList(cmd,'file',['Pileup','SeqIn'])  # String representing file path
                #self._cmdReadList(cmd,'date',['Att'])  # String representing date YYYY-MM-DD
                #self._cmdReadList(cmd,'bool',['Att'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['AbsMinCut','MinSNP','SNPCalc'])   # Integers
                self._cmdReadList(cmd,'float',['MinCut','SNPErr','TrackProb']) # Floats
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
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sam = self.obj['SAMTools']
            sam.parsePileup(self.getStr('Pileup'))
            #==> S288C-ISH.chrI.Q30.10.tdt <==
            #Locus	Pos	Ref	N	QN	Seq	Dep	RID
            #chrI_YEAST__BK006935	26	A	148	148	A:107|C:2	0.74	A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            #
            #==> S288C-ISH.chrI.QC.tdt <==
            #Qual	Count
            #1	0
            #
            #==> S288C-ISH.chrI.rid.tdt <==
            #RID	Locus	Start	End
            #chrI_YEAST__BK006935	176	61	351

            outfile = '%s.Q%d.%d.tdt' % (rje.baseFile(self.getStr('Pileup'),strip_path=True),sam.getInt('QCut'),sam.getInt('MinQN'))
            snpdb = self.db().openTable(outfile,mainkeys=['Locus','Pos'],name='snp',expect=True)
            # Filter based on allele frequency
            snpx = 0
            snpend = rje.endPos(snpdb.obj['File'])
            self.progLog('#SNP','%s SNP reduced to %s biallelic SNP meeting requirements.' % (rje.iStr(snpx),rje.iStr(snpdb.entryNum())))
            while snpdb:
                entry = snpdb.readEntry(add=False,close=True)
                if not entry: break
                snpx += 1
                alleles = {}
                for aseq in string.split(entry['Seq'],'|'):
                    [a,n] = string.split(aseq,':')
                    n = int(n)
                    if self.getNum('MinCut') < 1.0:
                        if n / float(entry['QN']) < self.getNum('MinCut'): continue
                        if n < self.getInt('AbsMinCut'): continue
                    elif n < self.getInt('MinCut'): continue
                    alleles[a] = n
                if len(alleles) != 2: continue   # Not biallelic: Need to filter better!
                entry['Seq'] = []
                for a in rje.sortKeys(alleles): entry['Seq'].append('%s:%d' % (a,alleles[a]))
                entry['Seq'] = string.join(entry['Seq'],'|')
                entry['RID'] = string.split(entry['RID'],'|')
                for rid in entry['RID'][0:]:
                    if string.split(rid,':')[0] not in alleles: entry['RID'].remove(rid)
                entry['RID'] = string.join(entry['RID'],'|')
                snpdb.addEntry(entry)
                self.progLog('#SNP','Parsing SNPs: %.2f%%; %s Pos -> %s biallelic SNP.' % (100.0 * snpdb.obj['File'].tell() / snpend,rje.iStr(snpx),rje.iStr(snpdb.entryNum())))
            self.printLog('#SNP','%s positions reduced to %s biallelic SNP meeting mincut requirements.' % (rje.iStr(snpx),rje.iStr(snpdb.entryNum())))

            #prex = snpdb.entryNum()
            #for entry in snpdb.entries():
            #    try: (a1,n1,a2,n2) = string.split(string.replace(entry['Seq'],'|',':'),':')
            #    except: self.warnLog('SNP entry format issues'); snpdb.dropEntry(entry); continue
            #    n1 = int(n1)
            #    n2 = int(n2)
            #    if min(n1,n2) < self.getInt('AbsMinCut'): snpdb.dropEntry(entry)
            #    elif self.getNum('MinCut') < 1.0:
            #        if float(n2) / (n1 + n2) < self.getNum('MinCut'): snpdb.dropEntry(entry)
            #    elif n2 < self.getInt('MinCut'): snpdb.dropEntry(entry)
            #self.printLog('#SNP','%s SNP reduced to %s biallelic SNP meeting requirements.' % (rje.iStr(prex),rje.iStr(snpdb.entryNum())))

            ridfile = '%s.rid.tdt' % (rje.baseFile(self.getStr('Pileup'),strip_path=True))
            self.db().addTable(ridfile,mainkeys=['RID'],datakeys='All',name='rid',expect=True)

            self.phase()

            #self.makeBlocks()

            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            samcmd = ['biallelic=T','rid=T','snponly=T','indels=F']
            sam = self.obj['SAMTools'] = rje_samtools.SAMtools(self.log,self.cmd_list+samcmd)
            sam.obj['DB'] = self.obj['DB']
            if not self.baseFile(return_none=None): self.baseFile(rje.baseFile(self.getStr('Pileup')))
            self.printLog('#BASE',self.baseFile(return_none=None))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT.
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
    def snpAX(self,locus,pos,allele):   ### Returns the probability allele in (locus,pos) being trackA
        '''Returns the probability allele in (locus,pos) being trackA.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db('snp')  #Locus Pos Ref N QN Seq Dep RID V1 R1 A1 V2 R2 A2 Block
            sentry = snpdb.data((locus,pos))
            if allele == sentry['V1']: return sentry['A1']
            elif allele == sentry['V2']: return sentry['A2']
            else: raise ValueError('(%s,%s,%s)' % (locus,pos,allele))
        except: self.errorLog('%s.snpA() error' % self.prog()); raise
#########################################################################################################################
    def snpA(self,locus,pos,phaserid,erate):    ### Sets the probability of each allele in (locus,pos) being trackA
        '''Sets the probability of each allele in (locus,pos) being trackA using phaserid RID.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpdb = self.db('snp')  #Locus Pos Ref N QN Seq Dep RID V1 R1 A1 V2 R2 A2 Block
            sentry = snpdb.data((locus,pos))
            riddb = self.db('rid')  #RID Locus Start End A B SNP
            ### ~ [2] Create Raw TrackA probability for each allele if new block ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not phaserid:
                if sentry['Ref'] == sentry['V2']:
                    sentry['A1'] = 0.0
                    sentry['A2'] = 1.0
                else:
                    sentry['A2'] = 0.0
                    sentry['A1'] = 1.0
                return 1    # A new block has been reached (no active phaserid)
            ### ~ [3] Calculate Raw TrackA probability for each allele ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            snpa = []
            for a in '12':
                pa = 0.0    # Log Prob of allele in Track A
                pb = 0.0    # Log Prob of allele in Track B
                for rid in sentry['R%s' % a]:
                    if rid not in phaserid: continue    # RID that have had previous SNPs
                    rp = riddb.data(rid)['A']
                    #if rp <= 0.0 or rp > 1.0: raise ValueError('RID: %s' % riddb.data(rid))
                    pa += math.log(rp*(1-erate) + (1-rp)*erate)
                    pb += math.log((1-rp)*(1-erate) + rp*erate)
                #    self.bugPrint('%s: %s => logA = %s; logB = %s' % (rid,rp,pa,pb))
                #self.bugPrint('%s => %s' % (pa,math.exp(pa)))
                pa = math.exp(pa)
                #self.bugPrint('%s => %s' % (pb,math.exp(pb)))
                pb = math.exp(pb)
                ptot = pa + pb
                #self.bugPrint(ptot)
                if ptot: snpa.append(pa/ptot)     # Relative Likelihood of Allele in Track A based on A RID
                else: snpa.append(0.5)
                #self.bugPrint('= %s\n' % snpa)
            try:
                if snpa[0] == snpa[1]: sentry['A1'] = 0.5   #?# Add warning for 0.0 values?
                else: sentry['A1'] = snpa[0] / (snpa[0]+snpa[1])
                sentry['A2'] = 1.0 - sentry['A1']
            except:
                self.debug('%s' % snpa)
                self.errorLog('%s' % (snpa[0]+snpa[1]))
                raise
            return 0    # Not a new block
        except: self.errorLog('%s.snpA() error' % self.prog()); self.debug(sentry); raise
#########################################################################################################################
    def readA(self,rid,erate,snpx=0):    ### Calculates the probability of read RID being trackA
        '''
        Calculates the probability of read RID being trackA.
        >> rid:str = RID identifier
        >> erate:float = Allele error rate
        >> snpx:int = Number of SNPs to use for probability calculation. (<1 = all)
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            riddb = self.db('rid')      #RID Locus Start End A B SNP
            rentry = riddb.data(rid)
            locus = rentry['Locus']
            snplist = rentry['SNP']     # List of (pos,allele)
            if snpx > 0: snplist = snplist[-snpx:]  # (pos,sentry['V1'])
            rentry['A'] = 0.0   # Log probability
            rentry['B'] = 0.0   # Log probability
            for (pos,var) in snplist:
                pa = self.snpAX(locus,pos,var)
                rentry['A'] += math.log(pa*(1-erate) + (1-pa)*erate)
                rentry['B'] += math.log((1-pa)*(1-erate) + pa*erate)
                #self.bugPrint(rentry)
            # Convert from log
            rentry['A'] = math.exp(rentry['A'])
            rentry['B'] = math.exp(rentry['B'])
            # Convert to relative likelihood
            ptot = rentry['A'] + rentry['B']
            if ptot:
                rentry['A'] /= ptot
                rentry['B'] /= ptot
            else:
                rentry['A'] = rentry['B'] = 0.5
            if max(ptot,rentry['A'],rentry['B']) > 1 or min(ptot,rentry['A'],rentry['B']) < 0:
                self.debug(rentry)
        except: self.errorLog('%s.readA() error' % self.prog()); raise
#########################################################################################################################
    def phase(self):    ### Main SNP phasing method.
        '''
        Main SNP phasing method.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            erate = self.getNum('SNPErr')    # Error rate (probability of incorrect base call)
            if erate <= 0: raise ValueError('Error rate (%s) <= 0!' % erate)
            snpx = self.getInt('SNPCalc')      # Number of SNPs to use for readA and to recalculate each round
            ## ~ [1a] Load and format SNPs for phasing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb = self.db('snp')  #Locus	Pos	Ref	N	QN	Seq	Dep	RID
            #chrI_YEAST__BK006935	26	A	148	148	A:107|C:2	0.74	A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            snpdb.dataFormat({'Pos':'int','N':'int','QN':'int','Dep':'num','RID':'str'})
            snpdb.index('Locus')
            snpdb.addFields(['V1','R1','A1','V2','R2','A2','Block'])
            for sentry in snpdb.entries():
                alldata = []
                for snpdata in string.split(sentry['RID'],'|'):
                    alldata += string.split(snpdata,':')
                sentry['V1'] = alldata[0]
                sentry['R1'] = string.split(alldata[1],',')
                sentry['A1'] = 0.0
                sentry['V2'] = alldata[2]
                sentry['R2'] = string.split(alldata[3],',')
                sentry['A2'] = 0.0
                sentry['Block'] = 0
            #?# Add filter based on read counts?
            snpdb.dropField('RID')
            ## ~ [1b] Load and format RID for phasing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            riddb = self.db('rid')  #RID	Locus	Start	End
            riddb.dataFormat({'RID':'str','Start':'int','End':'int'})
            #chrI_YEAST__BK006935	176	61	351
            riddb.addField('A')     # Prob of trackA
            riddb.addField('B')     # Prob of trackB
            riddb.addField('SNP')   # List of SNPs (pos,allele)
            riddb.addField('Block') # Haplotype Block
            for rentry in riddb.entries():
                rentry['A'] = 0.5
                rentry['B'] = 0.5
                rentry['SNP'] = []
                rentry['Block'] = 0

            ### ~ [2] Phase ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            blockx = 0          # Number of phased blocks
            blockpos = []       # List of SNP positions for each block (0<N for blocks 1-N)
            for locus in snpdb.indexKeys('Locus'):
                phasing = []    # List of reads currently being phased
                #trackA = []     # List of (pos,allele) assigned to track A
                #readA = {}      # Dictionary of {read:prob in A}
                #readB = {}      # Dictionary of {read:prob in B}
                #readV = {}      # Dictionary of {read:[(pos,allele)]}
                #snpA = {}       # Dictionary of {(pos,allele):prob in A} [Should always have 2 alleles per Locus/pos]
                px = 0.0; ptot = len(snpdb.index('Locus')[locus]); pstr = rje.iStr(ptot)
                for skey in snpdb.index('Locus')[locus]:
                    self.progLog('\r#PHASE','Phasing %s %s SNPs: %.2f%% => %d blocks' % (pstr,locus,px/ptot,blockx)); px += 100.0
                    sentry = snpdb.data(skey)
                    pos = sentry['Pos']
                    ## ~ [2a] Drop readA entries for finished reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for rid in phasing[0:]:
                        if riddb.data(rid)['End'] < pos: phasing.remove(rid)
                    ## ~ [2b] Calculate snpA statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    blockx += self.snpA(locus,pos,phasing,erate)
                    sentry['Block'] = blockx
                    if blockx > len(blockpos): blockpos.append([pos])
                    else: blockpos[blockx-1].append(pos)
                    ## ~ [2c] Update read SNP lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    for a in '12':
                        for rid in sentry['R%s' % a]:
                            riddb.data(rid)['SNP'].append((pos,sentry['V%s' % a]))
                            if rid not in phasing: phasing.append(rid); riddb.data(rid)['Block'] = blockx
                    ## ~ [2c] Calculate readA statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    # Make this a method using last x SNPs #
                    for rid in phasing: self.readA(rid,erate,snpx)
                    ## ~ [2d] Recalculate X snpA statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    # Keep a list of SNP positions in block and call self.snpA for last X using all RID
                    if snpx >= 0: recalc = blockpos[blockx-1][-snpx:]
                    else: recalc = blockpos[blockx-1][0:]
                    for pos in recalc:
                        snprid = snpdb.data((locus,pos))['R1'] + snpdb.data((locus,pos))['R2']
                        self.snpA(locus,pos,snprid,erate)
                self.printLog('\r#PHASE','Phased %s %s SNPs: %d phased total haplotype blocks' % (pstr,locus,blockx)); px += 100.0

            ### ~ [3] Resolve Blocks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #badblock = []
            #for bi in range(blockx):
            #    if len(blockpos[bi]) < self.getInt('MinSNP'): badblock.append(bi+1)
            #snpdb.dropEntriesDirect('Block',badblock)
            #self.printLog('\r#PHASE','%d of %d phased haplotype blocks have %d+ SNPs' % (blockx-len(badblock),blockx,self.getInt('MinSNP')))
            if not snpdb.entryNum(): self.printLog('#PHASE','No SNPs left in haplotype blocks!'); return False
            #for rentry in riddb.entries():
            #    if rentry['Block'] in badblock: rentry['Block'] = 0
            ## ~ [3a] Dev output for checking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.dev():
                for entry in snpdb.entries():
                    entry['R1'] = string.join(entry['R1'],',')
                    entry['R2'] = string.join(entry['R2'],',')
                snpdb.saveToFile()
                riddb.saveToFile(filename='%s.devrid.tdt' % self.baseFile())
            ## ~ [3b] Assign RID to tracks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            riddb.addField('Track')
            riddb.addField('pTrack')
            for rentry in riddb.entries():
                if rentry['A'] >= rentry['B']: rentry['Track'] = 'A'
                else: rentry['Track'] = 'B'
                rentry['pTrack'] = rentry[rentry['Track']]
                rentry['SNP'] = self.snpListToStr(rentry['SNP'])
            riddb.setFields(['RID','Locus','Start','End','Block','Track','pTrack','SNP'])
            riddb.dropEntries(['pTrack<%f' % self.getNum('TrackProb')])
            riddb.saveToFile(filename='%s.haprid.tdt' % self.baseFile())
            ## ~ [3c] Create Blocks Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb = self.db().copyTable(riddb,'blocks',replace=True,add=True)
            bdb.dropFields(['pTrack','SNP'])
            bdb.compress(['Locus','Block','Track'],rules={'RID':'list','Start':'min','End':'max'})
            bdb.setFields(['Locus','Block','Track','Start','End','HapStart','HapEnd','SNP'])
            for bentry in bdb.entries(): bentry['SNP'] = []   # Populate from snpdb
            ## ~ [3b] Assign SNPs to tracks  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb.addFields(['A','B','pA'])
            for sentry in snpdb.entries():
                if sentry['A1'] >= sentry['A2']:
                    sentry['A'] = sentry['V1']
                    sentry['B'] = sentry['V2']
                    sentry['pA'] = sentry['A1']
                else:
                    sentry['A'] = sentry['V2']
                    sentry['B'] = sentry['V1']
                    sentry['pA'] = sentry['A2']
            snpdb.dropEntries(['pA<%f' % self.getNum('TrackProb')])
            snpdb.setFields(['Locus','Block','Pos','Ref','A','B','pA'])
            # Add SNPs to Block Tables
            for skey in snpdb.dataKeys():
                sentry = snpdb.data(skey)
                for track in 'AB':
                    bdb.data((sentry['Locus'],sentry['Block'],track))['SNP'].append((sentry['Pos'],sentry[track]))
            for bentry in bdb.entries():
                if bentry['Block'] == 0:
                    bentry['HapStart'] = 0
                    bentry['HapEnd'] = 0
                    bentry['SNP'] = ''
                else:
                    bentry['HapStart'] = bentry['SNP'][0][0]
                    bentry['HapEnd'] = bentry['SNP'][-1][0]
                    bentry['SNP'] = self.snpListToStr(bentry['SNP'])
            # Add Block Starts and Ends
            for bentry in bdb.entries():
                if bentry['Block'] == 0: continue
                sentry = rje.combineDict({'A':'.','B':'.','Ref':'.','pA':1.0},bentry)
                sentry['Pos'] = sentry.pop('Start'); sentry.pop('End'); sentry.pop('SNP')
                track = sentry.pop('Track')
                sentry[track] = '^'
                if snpdb.makeKey(sentry) in snpdb.dataKeys() and snpdb.data(snpdb.makeKey(sentry))[track] == '.':
                    snpdb.data(snpdb.makeKey(sentry))[track] = '^'    # Tied start point
                elif not snpdb.addEntry(sentry,overwrite=False): self.warnLog('%s:%s%s Start has existing SNP/Block end! (Not added to SNP table.)' % (sentry['Locus'],sentry['Block'],bentry['Track']))
                sentry = rje.combineDict({'A':'.','B':'.','Ref':'.','pA':1.0},bentry)
                sentry['Pos'] = sentry.pop('End'); sentry.pop('Start'); sentry.pop('SNP')
                track = sentry.pop('Track')
                sentry[track] = '$'
                if snpdb.makeKey(sentry) in snpdb.dataKeys() and snpdb.data(snpdb.makeKey(sentry))[track] == '.':
                    snpdb.data(snpdb.makeKey(sentry))[track] = '$'    # Tied start point
                elif not snpdb.addEntry(sentry,overwrite=False): self.warnLog('%s:%s%s End has existing SNP/Block end! (Not added to SNP table.)' % (sentry['Locus'],sentry['Block'],bentry['Track']))
            snpdb.saveToFile(filename='%s.hapsnp.tdt' % self.baseFile())
            bdb.saveToFile()

        except: self.errorLog('%s.phase() error' % self.prog())
#########################################################################################################################
    def snpListToStr(self,snplist,snpjoin='|'):
        snpstr = []
        for snp in snplist: snpstr.append('%d:%s' % snp)
        return string.join(snpstr,snpjoin)
#########################################################################################################################
    def unzipSeq(self):     ### Unzip sequences based on SNP phasing
        '''Unzip sequences based on SNP phasing.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
            minx = self.getInt('MinSNP')
            #!# To be added!
        except: self.errorLog('%s.unzipSeq() error' % self.prog())
#########################################################################################################################
### End of SECTION II: SAMPhaser Class                                                                                  #
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
    try: SAMPhaser(mainlog,cmd_list).run()

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
