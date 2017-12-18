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
Version:      0.5.0
Last Edit:    22/11/17
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice

Function:
    The function of this module will be added here.

Commandline:
    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input genome to phase variants in []
    pileup=FILE     : Pileup file of reads against input genome. []
    basefile=FILE   : Root of output file names (same as Pileup input file by default) []
    mincut=X        : Minimum read count for minor allele (proportion of QN if <1) for pileup parsing [0.1]
    absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    indels=T/F      : Whether to include indels in "SNP" parsing [True]
    snptableout=T/F : Output filtered alleles to SNP Table [False]
    ### ~ Phasing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    phasecut=X      : Minimum read count for minor allele for phasing (proportion of QN if <1) [0.25]
    absphasecut=X   : Absolute minimum read count for phasecut (used if phasecut<1) [5]
    phaseindels=T/F : Whether to include indels in "SNP" phasing [False]
    snperr=X        : Probability of an incorrect (biallelic) SNP call for individual read nucleotides [0.05]
    snpcalc=X       : Max number of SNPs to use for read probability calculations (fewer = quicker) [10]
    trackprob=X     : Min probability for assigning a read/SNP to Track A/B [0.95]
    ### ~ Unzipping options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minsnp=X        : Min number of SNPs per phased haplotype block [5]
    endmargin=X     : Extend block ends within X nucleotides of sequence ends [10]
    unzipcut=X      : Minimum read count for allele for unzipping haplotigs (proportion of QN if <1) [0.1]
    absunzipcut=X   : Absolute minimum read count for unzipcut (used if unzipcut<1) [3]
    minhapx=X       : Minimum mean coverage for haplotig [5]
    halfhap=T/F     : Whether to allow "half haplotigs" where one halpotig in a pair is removed by minhapx [True]
    splitzero=X     : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100] (dev only)
    rgraphics=T/F   : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
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
import rje, rje_html, rje_obj, rje_db, rje_samtools, rje_seqlist
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Updated SAMPhaser to be more memory efficient.
    # 0.2.0 - Added reading of sequence and generation of SNP-altered haplotype blocks.
    # 0.2.1 - Fixed bug in which zero-phasing sequences were being excluded from blocks output.
    # 0.3.0 - Made a new unzip process.
    # 0.4.0 - Added RGraphics for unzip.
    # 0.4.1 - Fixed MeanX bug in devUnzip.
    # 0.4.2 - Made phaseindels=F by default: mononucleotide indel errors will probably add phasing noise. Fixed basefile R bug.
    # 0.4.3 - Fixed bug introduced by adding depthplot code. Fixed phaseindels bug. (Wasn't working!)
    # 0.4.4 - Modified mincut=X to adjust for samtools V1.12.0.
    # 0.4.5 - Updated for modified RJE_SAMTools output.
    # 0.4.6 - splitzero=X : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100] (dev only)
    # 0.5.0 - snptable=T/F    : Output filtered alleles to SNP Table [False]
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
    # [ ] : Add forking to unzip to do one locus per fork.
    # [ ] : Add rje.sf() output for the probabilites. This should be part of db.saveToFile().
    # [ ] : Add depthplot output for the original input and then generate stacked plots like PAGSAT assembly plots.
    # [ ] : -- Plot haplotigs against input Locus and mark SNP regions. (And phased SNPs or SNP density?)
    # [ ] : Give code a clean and tidy! (Lose old unzip.)
    # [ ] : Option to ignore mononucleotide repeat indels.
    # [ ] : splitzero=X : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100] (dev only)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SAMPhaser', '0.5.0', 'November 2017', '2016')
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
    - HalfHap=T/F     : Whether to allow "half haplotigs" where one halpotig in a pair is removed by minhapx [True]
    - PhaseIndels=T/F : Whether to include indels in "SNP" phasing [False]
    - RGraphics=T/F   : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
    - SNPTableOut=T/F    : Output filtered alleles to SNP Table [False]

    Int:integer
    - AbsPhaseCut=X   : Absolute minimum read count for phasecut (used if phasecut<1) [5]
    - AbsUnzipCut=X   : Absolute minimum read count for unzipcut (used if unzipcut<1) [3]
    - EndMargin=X     : Extend block ends within X nucleotides of sequence ends [10]
    - MinSNP=X        : Min number of SNPs per phased haplotype block [5]
    - SNPCalc=X       : Max number of SNPs to use for read probability calculations (fewer = quicker) [10]

    Num:float
    - MinHapX=X       : Minimum mean coverage for haplotig [5]
    - PhaseCut=X      : Minimum read count for minor allele for phasing (proportion if <1) [0.25]
    - SNPErr=X        : Probability of an incorrect (biallelic) SNP call for individual read nucleotides [0.05]
    - SplitZero=X     : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100] (dev only)
    - TrackProb=X     : Min probability for assigning a read/SNP to Track A/B [0.95]
    - UnzipCut=X      : Minimum read count for allele for unzipping haplotigs (proportion if <1) [0.1]

    File:file handles with matching str filenames
    
    List:list

    Dict:dictionary    

    Obj:RJE_Objects
    - SAMTools = rje_samtools.SAMtools()

    SAMTools settings used by SAMPhaser:
    - AbsMinCut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    - MinCut=X        : Minimum read count for minor allele (proportion if <1) [0.25]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['Pileup','SeqIn']
        self.boollist = ['HalfHap','PhaseIndels','RGraphics','SNPTableOut']
        self.intlist = ['AbsPhaseCut','AbsUnzipCut','EndMargin','MinSNP','SNPCalc']
        self.numlist = ['MinHapX','PhaseCut','SNPErr','SplitZero','TrackProb','UnzipCut']
        self.filelist = []
        self.listlist = []
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True,setfile=True)
        self.setStr({})
        self.setBool({'HalfHap':True,'PhaseIndels':False,'RGraphics':True,'SNPTableOut':False})
        self.setInt({'AbsPhaseCut':5,'AbsUnzipCut':3,'EndMargin':10,'MinSNP':5,'SNPCalc':10})
        self.setNum({'MinHapX':5.0,'PhaseCut':0.25,'SNPErr':0.05,'SplitZero':100,'TrackProb':0.95,'UnzipCut':0.1})
        ### ~ Other Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setForkAttributes()   # Delete if no forking
        #SeqList loaded in self.unzipSeq
        #seqcmd = ['seqmode=file','autoload=T']
        #self.obj['SeqList'] = rje_seqlist.SeqList(self.log,self.cmd_list+seqcmd)
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
                self._cmdReadList(cmd,'bool',['HalfHap','PhaseIndels','RGraphics','SNPTableOut'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['EndMargin','MinSNP','SNPCalc','SplitZero'])   # Integers
                self._cmdReadList(cmd,'float',['MinHapX','PhaseCut','SNPErr','TrackProb','UnzipCut']) # Floats
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
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            #i# Removed biallelic=T from default as might want to use non-biallelic positions unzipping.
            #i# indels=T is now an option
            samcmd = ['rid=T','snponly=T']
            sam = self.obj['SAMTools'] = rje_samtools.SAMtools(self.log,['mincut=0.1']+self.cmd_list+samcmd)
            sam.obj['DB'] = self.obj['DB']
            if not self.baseFile(return_none=None): self.baseFile(rje.baseFile(self.getStr('Pileup')))
            self.printLog('#BASE',self.baseFile(return_none=None))
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self.prog()); return False  # Setup failed
#########################################################################################################################
    def getSNPFile(self):   ### Returns the name of the parsed SNP File from RJE_SAMTools
        '''Returns the name of the parsed SNP File from RJE_SAMTools.'''
        sam = self.obj['SAMTools']
        qtxt = 'Q%d.%d' % (sam.getInt('QCut'),sam.getInt('MinQN'))
        outbase = '%s.%s' % (rje.baseFile(self.getStr('Pileup'),strip_path=True),qtxt)
        if sam.getBool('SNPOnly'):
            suffix = 'b%si%s' % (str(sam.getBool('Biallelic'))[:1],str(sam.getBool('Indels'))[:1])
            outfile = '%s.snponly.%s.tdt' % (outbase,suffix)
        else: outfile = '%s.tdt' % (outbase)
        return outfile
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()
            if self.db('blocks',add=True,mainkeys=['Locus','Block','Track']):
                self.unzipSeq()
                if self.getBool('RGraphics'): self.report()

                #if self.dev() and self.getInt('SplitZero') > 0:
                #    self.devZeroSplit()

                return True
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Parse SNPs with SAMTools object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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

            #!# Add indels filter? #!#



            ## ~ [2b] Filter based on allele frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb = self.db().openTable(self.getSNPFile(),mainkeys=['Locus','Pos'],name='snp',expect=True)
            #i# At this point we have filtered SNPs
            if self.getBool('SNPTableOut'):
                snptabdb = self.db().addEmptyTable('snptable',['Locus','Pos','REF','ALT','N','Freq'],['Locus','Pos','REF','ALT'],log=True)
            snpx = 0  # Number of SNPs read from SNP table
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
                    if self.getNum('PhaseCut') < 1.0:
                        if n / float(entry['QN']) < self.getNum('PhaseCut'): continue
                        if n < self.getInt('AbsPhaseCut'): continue
                    elif n < self.getInt('PhaseCut'): continue
                    if (a[0] == '-' or len(a) > 1) and not self.getBool('PhaseIndels'): continue
                    alleles[a] = n
                #?# Add biallelic toggle?
                if len(alleles) != 2: continue   # Not biallelic: Need to filter better!
                entry['Seq'] = []
                for a in rje.sortKeys(alleles):
                    entry['Seq'].append('%s:%d' % (a,alleles[a]))
                    if self.getBool('SNPTableOut') and entry['Ref'] != a:
                        allentry = {'Locus':entry['Locus'],'Pos':int(entry['Pos']),'REF':entry['Ref'],'ALT':a,'N':alleles[a],
                                    'Freq':float(alleles[a])/sum(alleles.values())}
                        snptabdb.addEntry(allentry)
                entry['Seq'] = string.join(entry['Seq'],'|')
                entry['RID'] = string.split(entry['RID'],'|')
                for rid in entry['RID'][0:]:
                    if string.split(rid,':')[0] not in alleles: entry['RID'].remove(rid)
                entry['RID'] = string.join(entry['RID'],'|')
                snpdb.addEntry(entry)
                self.progLog('#SNP','Parsing SNPs: %.2f%%; %s Pos -> %s biallelic SNP.' % (100.0 * snpdb.obj['File'].tell() / snpend,rje.iStr(snpx),rje.iStr(snpdb.entryNum())))
            self.printLog('#SNP','%s positions reduced to %s biallelic SNP meeting PhaseCut requirements (indels=%s).' % (rje.iStr(snpx),rje.iStr(snpdb.entryNum()),self.getBool('PhaseIndels')))
            if self.getBool('SNPTableOut'):
                tabfile = '%s.snptable.tdt' % rje.baseFile()
                self.printLog('#SNP','%s alleles parsed from %s biallelic SNP meeting PhaseCut requirements (indels=%s).' % (rje.iStr(snptabdb.entryNum()),rje.iStr(snpdb.entryNum()),self.getBool('PhaseIndels')))
                snptabdb.saveToFile(tabfile)



            ridfile = '%s.rid.tdt' % (rje.baseFile(self.getStr('Pileup'),strip_path=True))
            riddb = self.db().addTable(ridfile,mainkeys=['RID'],datakeys='All',name='rid',expect=True)
            riddb.dataFormat({'RID':'str','Start':'int','End':'int'})


            self.phase()

            # The unzipSeq() method needs blocks.tdt output
            # This needs to use the full (not biallelic) data!
            self.unzipSeq()

            return
        except:
            self.errorLog(self.zen())
            raise   # Delete this if method error not terrible
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
        except: self.errorLog('%s.snpAX() error' % self.prog()); raise
#########################################################################################################################
    def snpA(self,locus,pos,phaserid,erate):    ### Sets the probability of each allele in (locus,pos) being trackA
        '''
        Sets the probability of each allele in (locus,pos) being trackA using phaserid RID.
        >> locus:str = Locus from pileup file
        >> pos:int = SNP position under consideration.
        >> phaserid:list = List of RIDs being currently phased.
        >> erate:float = Probability of an incorrect (biallelic) SNP call for individual read nucleotides (snperr=X)
        '''
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
            if not self.force() and self.db('blocks',add=True,mainkeys=['Locus','Block','Track']):
                self.printLog('#PHASE','%s.blocks.tdt found: phasing skipped' % self.baseFile())
                return True
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
                sentry['V1'] = alldata[0]                       # Allele 1
                sentry['R1'] = string.split(alldata[1],',')     # List of RID for allele 1
                sentry['A1'] = 0.0                              # Probability of Track A for Allele 1
                sentry['V2'] = alldata[2]                       # Allele 2
                sentry['R2'] = string.split(alldata[3],',')     # List of RID for allele 2
                sentry['A2'] = 0.0                              # Probability of Track A for Allele 2
                sentry['Block'] = 0                             # Assigned haplotype block
            #?# Add filter based on read counts?
            snpdb.dropField('RID')
            ## ~ [1b] Load and format RID for phasing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            riddb = self.db('rid')  #RID	Locus	Start	End
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
                            rentry = riddb.data(rid)
                            if not rentry:
                                self.warnLog('RID %s not found in RIDdb' % rid)
                                continue
                            rentry['SNP'].append((pos,sentry['V%s' % a]))
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
            for entry in snpdb.entries():
                entry['R1'] = string.join(entry['R1'],',')
                entry['R2'] = string.join(entry['R2'],',')
            if self.dev():
                snpdb.saveToFile(sfdict={'A1':3,'A2':3})
                riddb.saveToFile(filename='%s.devrid.tdt' % self.baseFile(),sfdict={'pTrack':3})
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
            riddb.saveToFile(filename='%s.haprid.tdt' % self.baseFile(),sfdict={'pTrack':3})
            riddb.rename('haprid')
            ## ~ [3c] Create Blocks Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bdb = self.db().copyTable(riddb,'blocks',replace=True,add=True)
            bdb.dropFields(['pTrack','SNP'])
            bdb.compress(['Locus','Block','Track'],rules={'RID':'list','Start':'min','End':'max'})
            bdb.setFields(['Locus','Block','Track','Start','End','HapStart','HapEnd','SNP'])
            for bentry in bdb.entries(): bentry['SNP'] = []   # Populate from snpdb
            ## ~ [3b] Assign SNPs to tracks  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb.addFields(['A','B','nA','nB','pA'])
            for sentry in snpdb.entries():
                if sentry['A1'] >= sentry['A2']:
                    sentry['A'] = sentry['V1']
                    sentry['B'] = sentry['V2']
                    sentry['pA'] = sentry['A1']
                    sentry['nA'] = len(string.split(sentry['R1'],','))
                    sentry['nB'] = len(string.split(sentry['R2'],','))
                else:
                    sentry['A'] = sentry['V2']
                    sentry['B'] = sentry['V1']
                    sentry['pA'] = sentry['A2']
                    sentry['nA'] = len(string.split(sentry['R2'],','))
                    sentry['nB'] = len(string.split(sentry['R1'],','))
            snpdb.dropEntries(['pA<%f' % self.getNum('TrackProb')])
            snpdb.setFields(['Locus','Block','Pos','Ref','A','B','nA','nB','pA'])
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
                sentry = rje.combineDict({'A':'.','B':'.','Ref':'.','pA':1.0,'nA':0,'nB':0},bentry)
                sentry['Pos'] = sentry.pop('Start'); sentry.pop('End'); sentry.pop('SNP')
                track = sentry.pop('Track')
                sentry[track] = '^'
                if snpdb.makeKey(sentry) in snpdb.dataKeys() and snpdb.data(snpdb.makeKey(sentry))[track] == '.':
                    snpdb.data(snpdb.makeKey(sentry))[track] = '^'    # Tied start point
                elif not snpdb.addEntry(sentry,overwrite=False): self.warnLog('%s:%s%s Start has existing SNP/Block end! (Not added to SNP table.)' % (sentry['Locus'],sentry['Block'],bentry['Track']))
                sentry = rje.combineDict({'A':'.','B':'.','Ref':'.','pA':1.0,'nA':0,'nB':0},bentry)
                sentry['Pos'] = sentry.pop('End'); sentry.pop('Start'); sentry.pop('SNP')
                track = sentry.pop('Track')
                sentry[track] = '$'
                if snpdb.makeKey(sentry) in snpdb.dataKeys() and snpdb.data(snpdb.makeKey(sentry))[track] == '.':
                    snpdb.data(snpdb.makeKey(sentry))[track] = '$'    # Tied start point
                elif not snpdb.addEntry(sentry,overwrite=False): self.warnLog('%s:%s%s End has existing SNP/Block end! (Not added to SNP table.)' % (sentry['Locus'],sentry['Block'],bentry['Track']))
            snpdb.saveToFile(filename='%s.hapsnp.tdt' % self.baseFile(),sfdict={'pA':3})
            bdb.saveToFile()

        except: self.errorLog('%s.phase() error' % self.prog())
#########################################################################################################################
    def snpListToStr(self,snplist,snpjoin='|'):
        snpstr = []
        for snp in snplist: snpstr.append('%d:%s' % snp)
        return string.join(snpstr,snpjoin)
#########################################################################################################################
    def devUnzipSeq(self):     ### Unzip sequences based on SNP phasing
        '''Unzip sequences based on SNP phasing.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            hapfile = '%s.haplotigs.tdt' % self.baseFile()
            if not self.force() and rje.exists(hapfile):
                return self.printLog('#UNZIP','%s found (force=F): unzip cancelled.' % hapfile)
            self.printLog('#UNZIP','DEVELOPMENTAL UNZIP METHOD')
            sam = self.obj['SAMTools']
            sam.obj['DB'] = self.obj['DB']
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.makeSeqNameDic('max')
            minx = self.getInt('MinSNP')    # Min number of SNPs to unzip block


            #unzipcut=X      : Minimum read count for allele for unzipping haplotigs (proportion if <1) [0.1]
            #absunzipcut=X   : Absolute minimum read count for unzipcut (used if unzipcut<1) [3]
            #minhapx=X       : Minimum mean coverage for haplotig [5]
            #halfhap=T/F     : Whether to allow "half haplotigs" where one halpotig in a pair is removed by minhapx [True]



            #!# Add pickup at some point but need to save new blocks when doing so #!#
            #pickupfile = '%s.samphaser.pickup' % self.baseFile()
            #pickuplist = []
            #if rje.exists(pickupfile):
            #    if self.force(): os.unlink(pickupfile)
            #    else: pickuplist = self.loadFromFile(pickupfile,chomplines=True)
            #if pickuplist and pickuplist[-1] == '>>':
            #    self.warnLog('Warning! May have incomplete sequence output from previous run: check for duplicates')
            #while '>>' in pickuplist: pickuplist.remove('>>')
            #if pickuplist: self.printLog('#PICKUP','%s loci read from %s to skip' % (rje.iLen(pickuplist),pickupfile))

            #i# This will use the main QX.Y.tdt file to adjust/construct the sequence for each block from its reads
            #># Locus   Pos     Ref     N       QN      Seq     Dep     RID
            #># chrI_YEAST__BK006935    26      A       148     148     A:107|C:2       0.612   A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            #># chrI_YEAST__BK006935    29      A       151     151     A:120|C:2       0.624   A:1,2,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,26,28,29,31,32,33,34,35,36,38,39,40,41,42,43,45,46,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,68,72,73,74,75,77,78,80,81,82,83,84,85,86,88,89,90,91,92,94,96,97,98,99,100,101,102,105,106,107,110,112,113,114,115,116,117,119,120,121,122,124,127,128,129,130,132,133,134,135,136,137,138,139,140,141,142,144,145,146,147,149,150,151|C:30,93
            #i# Can load a chromosome at a time using table.readSet(['locus'],entry=None,clear=True)

            #i# Compress into RID list for each block
            #i# Also get a RID list for unassigned RIDs
            #i# Use simple winner-takes-all consensus generation
            #i# First identify the changes to make as (pos,sub) tuples and then work backwards to replace (in case of indels)
            #i# Also identify start/end positions

            #i# The "collapsed" C regions will be numbered separately.

            ### ~ [1] Open the data files for extraction of data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Use RID assignments from *.haprid.tdt
            #># RID     Locus   Start   End     Block   Track   pTrack  SNP
            #># 1       chrI_YEAST__BK006935    1       12483   15      B       0.997237569061  119:G|530:G|929:C|1263:T|1645:A|1779:T|2068:A|2137:A|2207:G|2295:C|2377:A|3049:T|3279:A|3518:A|4508:C|4696:A|4907:G|5446:T|5811:C|6645:A|6785:G|6972:A|11676:A|11732:A|11835:G|11997:A
            hapdb = self.db('haprid',add=True,mainkeys=['RID'])
            if not hapdb: raise IOError('%s.haprid.tdt not found!' % self.baseFile())
            hapdb.dataFormat({'Block':'int','Start':'int','End':'int','pTrack':'num'})
            hapdb.makeField('#Block##Track#','Haplotig')
            hapdb.indexReport('Locus')
            hapdb.indexReport('Haplotig')
            #i# Locus is being kept but will be replaced with new sequence AccNum for coverage plot.
            hapdb.addField('HapAcc')
            hapdb.dropFields(['Locus','Block','Track','SNP'])
            hapdb.rename('haplotigs.rid')

            # Load/use data from the *.blocks.tdt output
            #Locus	Block	Track	Start	End	HapStart	HapEnd	SNP
            #chrXII_YEAST__BK006945	1	A	1	468925	176	450984	176:G|11572:G|11993:C|13147:G|13381:G|13486:A|13571:...
            # Block table: will add to this with new blocks
            bdb = self.db('blocks',add=True,mainkeys=['Locus','Block','Track'])
            bdb.dataFormat({'Block':'int','Start':'int','End':'int','HapStart':'int','HapEnd':'int'})
            # Set number of blocks: will add to with "Collapsed" haplotigs
            blockx = max(bdb.dataList(bdb.entries(),'Block',sortunique=False,empties=False))
            # Generate haplotig IDs
            bdb.makeField('#Block##Track#','Haplotig')
            bdb.addField('HapAcc')
            bdb.newKey(['Haplotig'])
            bdb.addField('nSNP')
            for entry in bdb.entries(): entry['nSNP'] = len(string.split(entry['SNP'],'|'))
            bdb.dropField('SNP')
            bdb.dropEntries(['nSNP<%d' % minx],logtxt='nSNP<%d phased SNP filter' % minx)
            #self.printLog('#MINSNP','%s blocks entries with nSNP >= %d' % (rje.iStr(bdb.entryNum()),minx))
            if not bdb.entries(): return False
            bdb.rename('haplotigs') # Save at end
            #!# Save and append during output for pickup #!#
            #!# Add mimimum number of reads supporting (a) variant (3?), and (b) haplotig (10?)

            #i# Filter haplotigs based on coverage
            bdb.addField('nRID',evalue=0)
            bdb.addField('MeanX',evalue=0)
            for bentry in bdb.entries():
                hap = bentry['Haplotig']
                haplen = bentry['End'] - bentry['Start'] + 1.0
                for hentry in hapdb.indexEntries('Haplotig',hap):
                    bentry['MeanX'] += (hentry['End'] - hentry['Start'] + 1.0) / haplen
            bdb.dropEntries(['MeanX<%d' % self.getNum('MinHapX')],logtxt='MeanX<%d X Coverage filter' % self.getNum('MinHapX'))

            #i# Check for orphan haplotigs
            if not self.getBool('HalfHap'):
                drophap = []
                for hap in bdb.index('Block'):
                    if len(bdb.index('Block')[hap]) < 2: drophap.append(hap)
                    elif len(bdb.index('Block')[hap]) > 2: self.warnLog('Block %d has %d haplotigs!: %s' % (hap,len(bdb.index('Block')[hap]),string.join(bdb.index('Block')[hap],'; ')))
                bdb.dropEntries(drophap,keylist=True,logtxt='HalfHap=False filtering of solo haplotigs')

            #i# Extend haplotigs to ends of sequence if close
            for bkey in bdb.dataKeys():
                bentry = bdb.data(bkey)
                seqlen = seqlist.seqLen(seqdict[bentry['Locus']])
                if bentry['Start'] <= self.getInt('EndMargin'): bentry['Start'] = 1
                if bentry['End'] >= seqlen - self.getInt('EndMargin'): bentry['End'] = seqlen

                #!# Add splitting on zero-coverage here. Need to update bdb and hapdb.
                if self.dev() and self.getInt('SplitZero') >= 0:
                    hap = bentry['Haplotig'] # 1A etc.
                    self.progLog('\r#SPLIT','Zero Splitting %s...     ' % hap)
                    ridpos = []
                    for hentry in hapdb.indexEntries('Haplotig',hap):
                        ridpos.append((hentry['Start'],hentry['End']))
                    ridpos = rje.collapseTupleList(ridpos,joindistance=self.getInt('SplitZero'),overlaps=False)
                    self.bugPrint('%s: %s' % (hap,ridpos))
                    #i# No split
                    if len(ridpos) == 1: continue
                    #i# ridpos should now contain all the regions with coverage these are the new fragments!
                    if ridpos[0][0] <= self.getInt('EndMargin'): ridpos[0] = (1,ridpos[0][1])
                    if ridpos[-1][1] >= seqlen - self.getInt('EndMargin'): (ridpos[-1][0],seqlen)
                    #i# Clean up small chunks
                    splitfrag = []
                    for frag in ridpos:
                        if (frag[1] - frag[0] + 1) >= self.getInt('SplitZero'): splitfrag.append(frag)
                        else:
                            self.printLog('\r#SPLIT','Short %s fragment (%d-%d) removed.' % (hap,frag[0],frag[1]))
                    #i# Split
                    self.printLog('\r#SPLIT','%s -> %d fragments' % (hap,len(splitfrag)))
                    spliti = 0
                    while splitfrag:
                        spliti += 1
                        # Generate new block entry
                        sentry = rje.combineDict({},bentry)
                        #i# Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	HapAcc	nSNP	nRID	MeanX
                        #i# NOTE: nSNP will be from ORIGINAL Haplotig, not the split version!
                        #i# NOTE: nRID is calculated later
                        for field in ['Block','Haplotig','HapAcc']:
                            sentry[field] = '%s.%s' % (sentry[field],rje.preZero(spliti,len(splitfrag)))
                        (begpos,endpos) = splitfrag.pop(0)
                        self.printLog('\r#SPLIT','%s (%s-%s) -> %s' % (hap,rje.iStr(begpos),rje.iStr(endpos),sentry['Haplotig']))
                        sentry['Start'] = begpos
                        sentry['End'] = endpos
                        sentry['MeanX'] = 0
                        # Upate RID and RID table
                        haplen = endpos - begpos + 1.0
                        for hentry in hapdb.indexEntries('Haplotig',hap):
                            if hentry['End'] < begpos or hentry['Start'] > endpos: continue
                            hentry['Haplotig'] = sentry['Haplotig']
                            hentry['HapAcc'] = sentry['HapAcc']
                            sentry['MeanX'] += (hentry['End'] - hentry['Start'] + 1.0) / haplen
                        bdb.addEntry(sentry)
                    bdb.dropEntry(bentry)
            bdb.indexReport('Haplotig','HAPTIG',force=True)
            hapdb.indexReport('Haplotig','HAPTIG',force=True)

            # Parsed pileup data
            snpdb = self.db().openTable(self.getSNPFile(),mainkeys=['Locus','Pos'],name='pileup',expect=True)
            # Raw rids
            #># RID     Locus   Start   End
            #># 176     chrI_YEAST__BK006935    61      351
            #># 74      chrI_YEAST__BK006935    2       464
            #i# riddb is already loaded and formatted in main run method prior to phase()
            riddb = self.db('rid',add=True)
            ridfile = '%s.rid.tdt' % (rje.baseFile(self.getStr('Pileup'),strip_path=True))  #!# Check this!
            if not riddb:
                riddb = self.db().addTable(ridfile,mainkeys=['RID'],name='rid',expect=True)
            riddb.dataFormat({'RID':'str','Start':'int','End':'int'})
            riddb.rename('allrid')
            # Want coverage plot ready for unzipping
            sam.obj['DB'] = self.obj['DB']
            sam.db().baseFile(rje.baseFile(self.getStr('Pileup')))
            sam.setBool({'RGraphics':False})
            sam.coverageFromRID(ridfile=ridfile,depthplot=True)
            sam.db().baseFile(self.baseFile())

            # Set up output fasta file
            bfile = '%s.haplotigs.fas' % self.baseFile()
            #if pickuplist: HAPFAS = open(bfile,'a'); bx = 0; cx = 0; hx = 0
            #else:
            rje.backup(self,bfile)
            HAPFAS = open(bfile,'w'); bx = 0; cx = 0; hx = 0

            # Filter based on allele frequency
            locx = 0.0; loctot = seqlist.seqNum()
            self.progLog('#UNZIP','Unzipping SNP into haplotigs: %.f%%' % (locx/loctot))
            # Next collapsed block
            blockx += 1; colhap = '%dC' % blockx
            colentry = {'Locus':'None','Haplotig':colhap,'Block':blockx,'Track':'C','Start':0,'End':0,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
            while snpdb:
                locus = snpdb.readSet(['Locus'],entry=None,clear=True)
                if not locus: break
                locus = locus[0]    # readSet() returns a list
                #if locus in pickuplist: self.printLog('#SKIP','Locus %s in %s (force=F)' % (locus,pickupfile)); continue
                #else:
                self.printLog('#LOCUS',locus)
                snpdb.dataFormat({'Pos':'int'})
                (seqname,sequence) = seqlist.getSeq(seqdict[locus])
                #ridloc = riddb.readSet(['Locus'],entry={'Locus':locus},clear=True)[0]
                #ridloc = riddb.readSet(['Locus'],clear=True)[0]
                #if ridloc != locus: raise ValueError()
                #riddb.dataFormat({'Start':'int','End':'int'})
                # Generate a list of all RID in locus - will be reduced as blocks determined and used for remaining Collapsed region
                #locrid = riddb.dataKeys()           # NOTE: will be string numbers
                locrid = riddb.indexDataKeys('Locus',locus)
                # Generate list of haplotigs for this locus
                lochap = bdb.indexDataList('Locus',locus,'Haplotig')
                # Identify subset of reads in collapsed sections
                colrid = rje.listDifference(locrid,hapdb.dataKeys())    # Returns the elements of list1 that are not found in list 2.
                # Removed colrid that are wholly within phased region
                #i# Rejected by: riddb.dropEntries(['pTrack<%f' % self.getNum('TrackProb')])
                rejx = 0
                for rid in colrid[0:]:
                    colstart = riddb.data(rid)['Start']
                    colend = riddb.data(rid)['End']
                    for hentry in bdb.indexEntries('Locus',locus):
                        #if hentry['HapStart'] <= colstart and hentry['HapEnd'] >= colend:
                        if hentry['Start'] <= colstart and hentry['End'] >= colend:
                            colrid.remove(rid)
                            rejx += 1
                            break
                self.printLog('#REJECT','%s %s reads rejected as wholly within phased haplotig but poor pTrack' % (rje.iStr(rejx),locus))
                # Generate new collapsed block if needed (newxt locus)
                if colentry['Start']:
                    blockx += 1; colhap = '%dC' % blockx
                    colentry = {'Locus':'None','Haplotig':colhap,'Block':blockx,'Track':'C','Start':0,'End':0,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
                # Identify all RIDs in haplotig blocks
                haprid = hapdb.index('Haplotig')    # Dictionary of {haplotig:[RID list]}
                # Generate list of variants for making new sequences.
                locvar = {}     # Dictionary of {haplotig:[(pos,variant)]} for locus
                for hap in lochap: locvar[hap] = []

                # Work through variant positions and establish variants, along with Collapsed blocks
                cutx = 0
                snpx = 100.0 / snpdb.entryNum()
                for (locid,pos) in snpdb.dataKeys():
                    self.progLog('#UNZIP','Unzipping %s SNP into haplotigs: %.2f%%' % (locus,locx/loctot)); locx += snpx
                    if locid != locus: raise ValueError
                    entry = snpdb.data((locid,pos))
                    hapall = {}     # Allele counts per haplotig
                    posrid = []     # List of RID with parsed allele data for this position
                    allrid=  {}     # Dictionary of {alleles:ridlist} for alleles that have met unzip cutoffs
                    allcut = self.getNum('UnzipCut')
                    if allcut < 1: allcut = max(self.getNum('UnzipCut')*int(entry['QN']),self.getNum('AbsUnzipCut'))
                    #self.bugPrint(entry)
                    #self.bugPrint(allcut)
                    for alldata in string.split(entry['RID'],'|'):
                        (allele,ridlist) = string.split(alldata,':')
                        allrid[allele] = string.split(ridlist,',')
                        if len(allrid[allele]) >= allcut: posrid += allrid[allele]
                        else: allrid.pop(allele)
                    #self.debug(allrid.keys())
                    if allrid.keys() in [ [], [entry['Ref']] ]: cutx += 1; continue
                    colvar = rje.listIntersect(posrid,colrid)
                    # Check for new collapsed block and/or update current collapsed block
                    if colvar:
                        colrentries = riddb.entries(colvar)
                        colstart = min(riddb.dataList(colrentries,'Start'))
                        colend = max(riddb.dataList(colrentries,'End'))
                        if colentry['Start'] == 0:  # First one of locus
                            colentry['Locus'] = locus
                            colentry['Start'] = colstart
                            colentry['End'] = colend
                            colentry = bdb.addEntry(colentry)
                            lochap.append(colhap)
                            haprid[colhap] = colvar
                            locvar[colhap] = []
                        elif colstart > colentry['End']:    # Add another block
                            blockx += 1; colhap = '%dC' % blockx
                            colentry = {'Locus':locus,'Haplotig':colhap,'Block':blockx,'Track':'C',
                                        'Start':colstart,'End':colend,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
                            colentry = bdb.addEntry(colentry)
                            lochap.append(colhap)
                            haprid[colhap] = colvar
                            locvar[colhap] = []
                        else:
                            haprid[colhap] = rje.listUnion(haprid[colhap],colvar)
                            colentry['Start'] = min(colstart,colentry['Start'])
                            colentry['End'] = max(colend,colentry['End'])
                        for rid in colvar:
                            if hapdb.data(rid): continue
                            hapdb.addEntry(rje.combineDict({'Haplotig':colhap},riddb.data(rid)))

                    # Read in variants
                    for allele in allrid:
                        hapvar = {}
                        for hap in lochap: hapvar[hap] = rje.listIntersect(allrid[allele],haprid[hap])
                        for hap in hapvar:
                            if not hapvar[hap]: continue
                            if hap not in hapall: hapall[hap] = []
                            # hapall is a list of (number, allele) that can be reverse sorted to get dominant allele
                            hapall[hap].append((len(hapvar[hap]),allele))
                    for hap in hapall:
                        #i# hapall is a list of (number, allele) that can be reverse sorted to get dominant allele
                        hapall[hap].sort(reverse=True)
                        locvar[hap].append((pos,hapall[hap][0][1]))
                self.printLog('#UNZIP','Unzipping %s SNP into haplotigs: %s of %s positions failed to meet UnzipCut variant criteria' % (locus,rje.iStr(cutx),rje.iStr(snpdb.entryNum())))

                # Generate variant sequences for locus
                lseq = seqdict[locus]
                (seqname,sequence) = seqlist.getSeq(lseq,format='tuple')
                seqlen = len(sequence)
                sacc = seqlist.seqAcc(lseq)
                spec = seqlist.seqSpec(lseq)
                sgene = seqlist.seqGene(lseq)
                #open(pickupfile,'a').write('>>\n')
                for hap in lochap:
                    bentry = bdb.data(hap)
                    bentry['nRID'] = len(haprid[hap])
                    locvar[hap].sort(reverse=True)  # Make changes from end
                    #self.deBug(locvar[hap])
                    # Construct sequence name
                    bname = '%s%s%s_%s__%s.%s%s' % (sgene,bentry['Track'],bentry['Block'],spec,sacc,bentry['Track'],bentry['Block'])
                    bname += ' %s %s-%s of %s (%s block SNP' % (locus,rje.iStr(bentry['Start']),rje.iStr(bentry['End']),rje.iStr(seqlen),bentry['nSNP'])
                    if bentry['nSNP']:
                        bname = '%s; %s-%s)' % (bname,bentry['HapStart'],bentry['HapEnd']); hx += 1
                    else: bname = '%s; collapsed)' % (bname); cx += 1
                    # Construct sequence
                    bseq = sequence[0:bentry['End']]
                    for (pos,nt) in locvar[hap]:
                        pos -= 1                # Adjust for 0<L numbering
                        if nt == '-': nt = ''   # Do not add gaps!
                        elif nt[:1] in '-+': self.warnLog('Cannot process variant: %s' % nt); continue
                        if bseq[pos] != nt:
                            #before = bseq[pos-10:pos+10]
                            bseq = bseq[:pos] + nt + bseq[pos+1:]
                            #self.debug('%s %s >>\n %s ->\n %s\n' % (locus,snp,before,bseq[pos-10:pos+10]))
                    bseq = bseq[bentry['Start']-1:]
                    self.printLog('#SEQ',bname)
                    HAPFAS.write('>%s\n%s\n' % (bname,bseq)); bx +=1
                    bentry['HapAcc'] = '%s.%s' % (sacc,hap)
                    for hentry in hapdb.indexEntries('Haplotig',hap):
                        hentry['Start'] -= (bentry['Start']-1)
                        hentry['End'] -= (bentry['Start']-1)
                        hentry['HapAcc'] = '%s.%s' % (sacc,hap)
                        if min(hentry['Start'],hentry['End']) < 1: self.warnLog('%s %s RID%s starts/ends before haplotig!' % (locus,hap,hentry['RID']))
                #open(pickupfile,'a').write('%s\n' % locus)
            HAPFAS.close()
            self.printLog('#UNZIP','Unzipped SNP into %s haplotigs (%s phased) for %s loci.' % (rje.iStr(bx),rje.iStr(hx),rje.iStr(seqlist.seqNum())))
            self.printLog('#FAS','%s gapped phased and %s unphased blocks output to %s' % (rje.iStr(hx),rje.iStr(cx),bfile))
            self.db().deleteTable('allrid')
            self.db().deleteTable('pileup')

            # Save haplotigs file
            bdb.fillBlanks(blank=0,fields=['MeanX'],fillempty=True)
            for bentry in bdb.indexEntries('Track','C'):
                hap = bentry['Haplotig']
                haplen = bentry['End'] - bentry['Start'] + 1.0
                for hentry in hapdb.indexEntries('Haplotig',hap):
                    bentry['MeanX'] += (hentry['End'] - hentry['Start'] + 1.0) / haplen
            bdb.saveToFile(sfdict={'MeanX':3})
            hapdb.fillBlanks(blank=0.0,fields=['pTrack'],fillempty=True)
            hapdb.saveToFile(sfdict={'pTrack':4})
            #!# NOTE: There is an output in *.haplotigs.tdt without a HapAcc!

            # Read coverage plot
            sam.baseFile('%s.haplotigs' % self.baseFile())
            self.db().baseFile('%s.haplotigs' % self.baseFile())
            hapdb.rename('rid')
            hapdb.renameField('HapAcc','Locus')
            #x#sam.list['CheckFields'] = []
            sam.coverageFromRID(depthplot=True)
            sam.baseFile(self.baseFile())

            # Generate SAMPhase stacked graphics
            if self.getBool('RGraphics'):
                sam.setBool({'RGraphics':True})
                sam.rGraphics('samunzip','"sambase=%s"' % rje.baseFile(self.getStr('Pileup')))

        except: self.errorLog('%s.unzipSeq() error' % self.prog())
#########################################################################################################################
    def devZeroSplit(self): ###
        '''Unzip sequences based on SNP phasing.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#UNZIP','DEVELOPMENTAL ZERO-SPLIT METHOD')
            #!# Add force check for existing output
            #?# Add a table tracking the splitting?
            hapbase = '%s.haplotigs' % self.baseFile()
            db = self.db()
            depdb = db.addTable('%s.depthplot.tdt' % hapbase,mainkeys=['Locus','Pos'],name='depthplot',expect=True)      # ['Locus', 'Pos', 'X']
            depdb.dataFormat({'Pos':'int','X':'int'})
            depdb.index('Locus')
            bfile = '%s.fas' % hapbase
            splitfile = '%s.split.fas' % hapbase
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file','seqin=%s' % bfile])
            #i# seqdict keys will match depdb 'Locus' field
            seqdict = seqlist.makeSeqNameDic('short')

            #i# Loci will be split into N fragments and (if N>1) renamed Locus.1 - Locus.N
            locx = 0; splitx = 0; fragx = 0
            seqfrag = []    # List of (name,sequence) fragmented sequences
            for locus in rje.sortKeys(seqdict):
                locx += 1
                (sname,sequence) = seqlist.getSeq(seqdict[locus])
                #i# >sgdIA6_YEAST__BK006935.A6 sgdI_YEAST__BK006935 1-14,881 of 230,218 (109 SNP; 72-13302)
                desc = string.split(sname)[1:]  # desc[1] should be the sequence positions of the original sequence
                beghap = int(string.split(desc[1],'-')[0].replace(',',''))
                endhap = int(string.split(desc[1], '-')[1].replace(',', ''))
                splitn = []
                begpos = 1
                endpos = 0
                #!# Convert locus: need to standardise this later #!#
                #!# Modify the unsplit gapped sequences to match the depthplot loci
                #!# Then switch the ordering here to make it clear that sequence lengths are different!
                dlocus = string.split(string.split(locus,'_')[-1],'.')
                if rje.matchExp('^(\D+)(\d+)',dlocus[1]): dlocus[1] = dlocus[1][1:] + dlocus[1][:1]
                dlocus = string.join(dlocus,'.')
                self.bugPrint('%s -> %s: %s' % (locus,dlocus,len(depdb.indexEntries('Locus',dlocus))))
                self.deBug('%d == %d?' % (len(sequence),max(depdb.indexDataList('Locus',dlocus,'Pos'))))
                for entry in depdb.indexEntries('Locus',dlocus):
                    if entry['X']: continue  # Coverage
                    self.bugPrint(entry)
                    # Identify start of zero-coverage region
                    if not endpos:
                        endpos = entry['Pos']
                    # End of zero-coverage but too short
                    elif (entry['Pos'] - endpos + 1) < self.getInt('SplitZero'):
                        self.bugPrint('%d->%d < %d' % (endpos,entry['Pos'],self.getInt('SplitZero')))
                        endpos = 0
                    # End of zero-coverage for split
                    else:
                        splitn.append((begpos,endpos))
                        self.bugPrint(splitn)
                        begpos = entry['Pos'] + 1
                        endpos = 0
                #i# Final position
                endpos = len(sequence)
                splitn.append((begpos, endpos))
                self.debug(splitn)
                #i# Clean up small chunks
                splitfrag = []
                for frag in splitn:
                    if (frag[1] - frag[0] + 1) >= self.getInt('SplitZero'): splitfrag.append(frag)
                    else:
                        self.printLog('#SPLIT','Short %s fragment (%d-%d) removed.' % (locus,frag[0],frag[1]))
                #i# Split if N > 1
                if len(splitn) > 1:
                    self.printLog('#SPLIT','%s -> %d fragments' % (sname,len(splitfrag)))
                    splitx += 1
                    fragx += len(splitfrag)
                    # Cycle through splits and make new sequence
                    spliti = 0
                    while splitfrag:
                        spliti += 1
                        (begpos,endpos) = splitfrag.pop(0)
                        newshort = '%s.%s' % (locus,rje.preZero(spliti,len(splitfrag)))
                        # Desc - update the start and end positions
                        newdesc = desc[:1] + ['%s-%s' % (rje.iStr(begpos-beghap+1),rje.iStr(endpos-beghap))] + desc[2:]
                        newdesc[-1] = '%s; ZeroSplit %d-%d)' % (desc[-1][:-1],beghap,endhap)
                        newname = string.join([newshort]+newdesc)
                        # Truncate sequence
                        newseq = sequence[begpos-1:endpos-1]
                        newseq = newseq.replace('-','')
                        seqfrag.append((newname,newseq))    # Add fragment
                        self.printLog('#SPLIT',newname)
                else:
                    fragx += 1
                    sequence = sequence.replace('-','')
                    seqfrag.append((sname,sequence))    # Not changed
            self.printLog('#SPLIT','%s haplotigs split on zero coverage: %s -> %s haplotigs.' % (rje.iStr(splitx),rje.iStr(locx),rje.iStr(fragx)))
            seqlist.saveSeq(seqfrag,seqfile=splitfile,seqtuples=True)

        except: self.errorLog('%s.devZeroSplit() error' % self.prog())
#########################################################################################################################
    def unzipSeq(self):     ### Unzip sequences based on SNP phasing
        '''Unzip sequences based on SNP phasing.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dev():
                return self.devUnzipSeq()
            hapfile = '%s.haplotigs.tdt' % self.baseFile()
            if not self.force() and rje.exists(hapfile):
                return self.printLog('#UNZIP','%s found (force=F): unzip cancelled.' % hapfile)
            self.printLog('#UNZIP','SAMPHASER HAPLOTIG UNZIP')
            sam = self.obj['SAMTools']
            sam.obj['DB'] = self.obj['DB']
            seqlist = rje_seqlist.SeqList(self.log,self.cmd_list+['autoload=T','seqmode=file'])
            seqdict = seqlist.makeSeqNameDic('max')
            minx = self.getInt('MinSNP')    # Min number of SNPs to unzip block


            #unzipcut=X      : Minimum read count for allele for unzipping haplotigs (proportion if <1) [0.1]
            #absunzipcut=X   : Absolute minimum read count for unzipcut (used if unzipcut<1) [3]
            #minhapx=X       : Minimum mean coverage for haplotig [5]
            #halfhap=T/F     : Whether to allow "half haplotigs" where one halpotig in a pair is removed by minhapx [True]



            #!# Add pickup at some point but need to save new blocks when doing so #!#
            #pickupfile = '%s.samphaser.pickup' % self.baseFile()
            #pickuplist = []
            #if rje.exists(pickupfile):
            #    if self.force(): os.unlink(pickupfile)
            #    else: pickuplist = self.loadFromFile(pickupfile,chomplines=True)
            #if pickuplist and pickuplist[-1] == '>>':
            #    self.warnLog('Warning! May have incomplete sequence output from previous run: check for duplicates')
            #while '>>' in pickuplist: pickuplist.remove('>>')
            #if pickuplist: self.printLog('#PICKUP','%s loci read from %s to skip' % (rje.iLen(pickuplist),pickupfile))

            #i# This will use the main QX.Y.tdt file to adjust/construct the sequence for each block from its reads
            #># Locus   Pos     Ref     N       QN      Seq     Dep     RID
            #># chrI_YEAST__BK006935    26      A       148     148     A:107|C:2       0.612   A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            #># chrI_YEAST__BK006935    29      A       151     151     A:120|C:2       0.624   A:1,2,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,22,23,24,25,26,28,29,31,32,33,34,35,36,38,39,40,41,42,43,45,46,48,49,50,51,52,53,54,55,56,57,58,59,61,62,63,64,65,66,68,72,73,74,75,77,78,80,81,82,83,84,85,86,88,89,90,91,92,94,96,97,98,99,100,101,102,105,106,107,110,112,113,114,115,116,117,119,120,121,122,124,127,128,129,130,132,133,134,135,136,137,138,139,140,141,142,144,145,146,147,149,150,151|C:30,93
            #i# Can load a chromosome at a time using table.readSet(['locus'],entry=None,clear=True)

            #i# Compress into RID list for each block
            #i# Also get a RID list for unassigned RIDs
            #i# Use simple winner-takes-all consensus generation
            #i# First identify the changes to make as (pos,sub) tuples and then work backwards to replace (in case of indels)
            #i# Also identify start/end positions

            #i# The "collapsed" C regions will be numbered separately.

            ### ~ [1] Open the data files for extraction of data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Use RID assignments from *.haprid.tdt
            #># RID     Locus   Start   End     Block   Track   pTrack  SNP
            #># 1       chrI_YEAST__BK006935    1       12483   15      B       0.997237569061  119:G|530:G|929:C|1263:T|1645:A|1779:T|2068:A|2137:A|2207:G|2295:C|2377:A|3049:T|3279:A|3518:A|4508:C|4696:A|4907:G|5446:T|5811:C|6645:A|6785:G|6972:A|11676:A|11732:A|11835:G|11997:A
            hapdb = self.db('haprid',add=True,mainkeys=['RID'])
            if not hapdb: raise IOError('%s.haprid.tdt not found!' % self.baseFile())
            hapdb.dataFormat({'Block':'int','Start':'int','End':'int','pTrack':'num'})
            hapdb.makeField('#Block##Track#','Haplotig')
            hapdb.indexReport('Locus')
            hapdb.indexReport('Haplotig')
            #i# Locus is being kept but will be replaced with new sequence AccNum for coverage plot.
            hapdb.addField('HapAcc')
            hapdb.dropFields(['Locus','Block','Track','SNP'])
            hapdb.rename('haplotigs.rid')

            # Load/use data from the *.blocks.tdt output
            #Locus	Block	Track	Start	End	HapStart	HapEnd	SNP
            #chrXII_YEAST__BK006945	1	A	1	468925	176	450984	176:G|11572:G|11993:C|13147:G|13381:G|13486:A|13571:...
            # Block table: will add to this with new blocks
            bdb = self.db('blocks',add=True,mainkeys=['Locus','Block','Track'])
            bdb.dataFormat({'Block':'int','Start':'int','End':'int','HapStart':'int','HapEnd':'int'})
            # Set number of blocks: will add to with "Collapsed" haplotigs
            blockx = max(bdb.dataList(bdb.entries(),'Block',sortunique=False,empties=False))
            # Generate haplotig IDs
            bdb.makeField('#Block##Track#','Haplotig')
            bdb.addField('HapAcc')
            bdb.newKey(['Haplotig'])
            bdb.addField('nSNP')
            for entry in bdb.entries(): entry['nSNP'] = len(string.split(entry['SNP'],'|'))
            bdb.dropField('SNP')
            bdb.dropEntries(['nSNP<%d' % minx],logtxt='nSNP<%d phased SNP filter' % minx)
            #self.printLog('#MINSNP','%s blocks entries with nSNP >= %d' % (rje.iStr(bdb.entryNum()),minx))
            if not bdb.entries(): return False
            bdb.rename('haplotigs') # Save at end
            #!# Save and append during output for pickup #!#
            #!# Add mimimum number of reads supporting (a) variant (3?), and (b) haplotig (10?)

            #i# Filter haplotigs based on coverage
            bdb.addField('nRID',evalue=0)
            bdb.addField('MeanX',evalue=0)
            for bentry in bdb.entries():
                hap = bentry['Haplotig']
                haplen = bentry['End'] - bentry['Start'] + 1.0
                for hentry in hapdb.indexEntries('Haplotig',hap):
                    bentry['MeanX'] += (hentry['End'] - hentry['Start'] + 1.0) / haplen
            bdb.dropEntries(['MeanX<%d' % self.getNum('MinHapX')],logtxt='MeanX<%d X Coverage filter' % self.getNum('MinHapX'))

            #i# Check for orphan haplotigs
            if not self.getBool('HalfHap'):
                drophap = []
                for hap in bdb.index('Block'):
                    if len(bdb.index('Block')[hap]) < 2: drophap.append(hap)
                    elif len(bdb.index('Block')[hap]) > 2: self.warnLog('Block %d has %d haplotigs!: %s' % (hap,len(bdb.index('Block')[hap]),string.join(bdb.index('Block')[hap],'; ')))
                bdb.dropEntries(drophap,keylist=True,logtxt='HalfHap=False filtering of solo haplotigs')

            #i# Extend haplotigs to ends of sequence if close
            for bentry in bdb.entries():
                seqlen = seqlist.seqLen(seqdict[bentry['Locus']])
                if bentry['Start'] <= self.getInt('EndMargin'): bentry['Start'] = 1
                if bentry['End'] >= seqlen - self.getInt('EndMargin'): bentry['End'] = seqlen


            # Parsed pileup data
            snpdb = self.db().openTable(self.getSNPFile(),mainkeys=['Locus','Pos'],name='pileup',expect=True)
            # Raw rids
            #># RID     Locus   Start   End
            #># 176     chrI_YEAST__BK006935    61      351
            #># 74      chrI_YEAST__BK006935    2       464
            #i# riddb is already loaded and formatted in main run method prior to phase()
            riddb = self.db('rid',add=True)
            ridfile = '%s.rid.tdt' % (rje.baseFile(self.getStr('Pileup'),strip_path=True))  #!# Check this!
            if not riddb:
                riddb = self.db().addTable(ridfile,mainkeys=['RID'],name='rid',expect=True)
            riddb.dataFormat({'RID':'str','Start':'int','End':'int'})
            riddb.rename('allrid')
            # Want coverage plot ready for unzipping
            sam.obj['DB'] = self.obj['DB']
            sam.db().baseFile(rje.baseFile(self.getStr('Pileup')))
            sam.setBool({'RGraphics':False})
            sam.coverageFromRID(ridfile=ridfile,depthplot=True)
            sam.db().baseFile(self.baseFile())

            # Set up output fasta file
            bfile = '%s.haplotigs.fas' % self.baseFile()
            #if pickuplist: HAPFAS = open(bfile,'a'); bx = 0; cx = 0; hx = 0
            #else:
            rje.backup(self,bfile)
            HAPFAS = open(bfile,'w'); bx = 0; cx = 0; hx = 0

            # Filter based on allele frequency
            locx = 0.0; loctot = seqlist.seqNum()
            self.progLog('#UNZIP','Unzipping SNP into haplotigs: %.f%%' % (locx/loctot))
            # Next collapsed block
            blockx += 1; colhap = '%dC' % blockx
            colentry = {'Locus':'None','Haplotig':colhap,'Block':blockx,'Track':'C','Start':0,'End':0,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
            while snpdb:
                locus = snpdb.readSet(['Locus'],entry=None,clear=True)
                if not locus: break
                locus = locus[0]    # readSet() returns a list
                #if locus in pickuplist: self.printLog('#SKIP','Locus %s in %s (force=F)' % (locus,pickupfile)); continue
                #else:
                self.printLog('#LOCUS',locus)
                snpdb.dataFormat({'Pos':'int'})
                (seqname,sequence) = seqlist.getSeq(seqdict[locus])
                #ridloc = riddb.readSet(['Locus'],entry={'Locus':locus},clear=True)[0]
                #ridloc = riddb.readSet(['Locus'],clear=True)[0]
                #if ridloc != locus: raise ValueError()
                #riddb.dataFormat({'Start':'int','End':'int'})
                # Generate a list of all RID in locus - will be reduced as blocks determined and used for remaining Collapsed region
                #locrid = riddb.dataKeys()           # NOTE: will be string numbers
                locrid = riddb.indexDataKeys('Locus',locus)
                # Generate list of haplotigs for this locus
                lochap = bdb.indexDataList('Locus',locus,'Haplotig')
                # Identify subset of reads in collapsed sections
                colrid = rje.listDifference(locrid,hapdb.dataKeys())    # Returns the elements of list1 that are not found in list 2.
                # Removed colrid that are wholly within phased region
                #i# Rejected by: riddb.dropEntries(['pTrack<%f' % self.getNum('TrackProb')])
                rejx = 0
                for rid in colrid[0:]:
                    colstart = riddb.data(rid)['Start']
                    colend = riddb.data(rid)['End']
                    for hentry in bdb.indexEntries('Locus',locus):
                        #if hentry['HapStart'] <= colstart and hentry['HapEnd'] >= colend:
                        if hentry['Start'] <= colstart and hentry['End'] >= colend:
                            colrid.remove(rid)
                            rejx += 1
                            break
                self.printLog('#REJECT','%s %s reads rejected as wholly within phased haplotig but poor pTrack' % (rje.iStr(rejx),locus))
                # Generate new collapsed block if needed (newxt locus)
                if colentry['Start']:
                    blockx += 1; colhap = '%dC' % blockx
                    colentry = {'Locus':'None','Haplotig':colhap,'Block':blockx,'Track':'C','Start':0,'End':0,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
                # Identify all RIDs in haplotig blocks
                haprid = hapdb.index('Haplotig')    # Dictionary of {haplotig:[RID list]}
                # Generate list of variants for making new sequences.
                locvar = {}     # Dictionary of {haplotig:[(pos,variant)]} for locus
                for hap in lochap: locvar[hap] = []

                # Work through variant positions
                cutx = 0
                snpx = 100.0 / snpdb.entryNum()
                for (locid,pos) in snpdb.dataKeys():
                    self.progLog('#UNZIP','Unzipping SNP into haplotigs: %.2f%%' % (locx/loctot)); locx += snpx
                    if locid != locus: raise ValueError
                    entry = snpdb.data((locid,pos))
                    hapall = {}     # Allele counts per haplotig
                    posrid = []     # List of RID with parsed allele data for this position
                    allrid=  {}     # Dictionary of {alleles:ridlist} for alleles that have met unzip cutoffs
                    allcut = self.getNum('UnzipCut')
                    if allcut < 1: allcut = max(self.getNum('UnzipCut')*int(entry['QN']),self.getNum('AbsUnzipCut'))
                    #self.bugPrint(entry)
                    #self.bugPrint(allcut)
                    for alldata in string.split(entry['RID'],'|'):
                        (allele,ridlist) = string.split(alldata,':')
                        allrid[allele] = string.split(ridlist,',')
                        if len(allrid[allele]) >= allcut: posrid += allrid[allele]
                        else: allrid.pop(allele)
                    #self.debug(allrid.keys())
                    if allrid.keys() in [ [], [entry['Ref']] ]: cutx += 1; continue
                    colvar = rje.listIntersect(posrid,colrid)
                    # Check for new collapsed block and/or update current collapsed block
                    if colvar:
                        colrentries = riddb.entries(colvar)
                        colstart = min(riddb.dataList(colrentries,'Start'))
                        colend = max(riddb.dataList(colrentries,'End'))
                        if colentry['Start'] == 0:  # First one of locus
                            colentry['Locus'] = locus
                            colentry['Start'] = colstart
                            colentry['End'] = colend
                            colentry = bdb.addEntry(colentry)
                            lochap.append(colhap)
                            haprid[colhap] = colvar
                            locvar[colhap] = []
                        elif colstart > colentry['End']:    # Add another block
                            blockx += 1; colhap = '%dC' % blockx
                            colentry = {'Locus':locus,'Haplotig':colhap,'Block':blockx,'Track':'C',
                                        'Start':colstart,'End':colend,'HapStart':0,'HapEnd':0,'SNP':0,'nSNP':0}
                            colentry = bdb.addEntry(colentry)
                            lochap.append(colhap)
                            haprid[colhap] = colvar
                            locvar[colhap] = []
                        else:
                            haprid[colhap] = rje.listUnion(haprid[colhap],colvar)
                            colentry['Start'] = min(colstart,colentry['Start'])
                            colentry['End'] = max(colend,colentry['End'])
                        for rid in colvar:
                            if hapdb.data(rid): continue
                            hapdb.addEntry(rje.combineDict({'Haplotig':colhap},riddb.data(rid)))

                    # Read in variants
                    for allele in allrid:
                        hapvar = {}
                        for hap in lochap: hapvar[hap] = rje.listIntersect(allrid[allele],haprid[hap])
                        for hap in hapvar:
                            if not hapvar[hap]: continue
                            if hap not in hapall: hapall[hap] = []
                            # hapall is a list of (number, allele) that can be reverse sorted to get dominant allele
                            hapall[hap].append((len(hapvar[hap]),allele))
                    for hap in hapall:
                        #i# hapall is a list of (number, allele) that can be reverse sorted to get dominant allele
                        hapall[hap].sort(reverse=True)
                        locvar[hap].append((pos,hapall[hap][0][1]))
                self.printLog('#UNZIP','Unzipping SNP into haplotigs: %s of %s %s positions failed to meet UnzipCut variant criteria' % (rje.iStr(cutx),rje.iStr(snpdb.entryNum()),locus))

                # Generate variant sequences for locus
                lseq = seqdict[locus]
                (seqname,sequence) = seqlist.getSeq(lseq,format='tuple')
                seqlen = len(sequence)
                sacc = seqlist.seqAcc(lseq)
                spec = seqlist.seqSpec(lseq)
                sgene = seqlist.seqGene(lseq)
                #open(pickupfile,'a').write('>>\n')
                for hap in lochap:
                    bentry = bdb.data(hap)
                    bentry['nRID'] = len(haprid[hap])
                    locvar[hap].sort(reverse=True)  # Make changes from end
                    # Construct sequence name
                    bname = '%s%s%s_%s__%s.%s%s' % (sgene,bentry['Track'],bentry['Block'],spec,sacc,bentry['Track'],bentry['Block'])
                    bname += ' %s %s-%s of %s (%s SNP' % (locus,rje.iStr(bentry['Start']),rje.iStr(bentry['End']),rje.iStr(seqlen),bentry['nSNP'])
                    if bentry['nSNP']:
                        bname = '%s; %s-%s)' % (bname,bentry['HapStart'],bentry['HapEnd']); hx += 1
                    else: bname = '%s; collapsed)' % (bname); cx += 1
                    # Construct sequence
                    bseq = sequence[0:bentry['End']]
                    for (pos,nt) in locvar[hap]:
                        pos -= 1                # Adjust for 0<L numbering
                        if nt == '-': nt = ''   # Do not add gaps!
                        elif nt[:1] in '-+': self.warnLog('Cannot process variant: %s' % nt); continue
                        if bseq[pos] != nt:
                            #before = bseq[pos-10:pos+10]
                            bseq = bseq[:pos] + nt + bseq[pos+1:]
                            #self.debug('%s %s >>\n %s ->\n %s\n' % (locus,snp,before,bseq[pos-10:pos+10]))
                    bseq = bseq[bentry['Start']-1:]
                    self.printLog('#SEQ',bname)
                    HAPFAS.write('>%s\n%s\n' % (bname,bseq)); bx +=1
                    bentry['HapAcc'] = '%s.%s' % (sacc,hap)
                    for hentry in hapdb.indexEntries('Haplotig',hap):
                        hentry['Start'] -= (bentry['Start']-1)
                        hentry['End'] -= (bentry['Start']-1)
                        hentry['HapAcc'] = '%s.%s' % (sacc,hap)
                        if min(hentry['Start'],hentry['End']) < 1: self.warnLog('%s %s RID%s starts/ends before haplotig!' % (locus,hap,hentry['RID']))
                #open(pickupfile,'a').write('%s\n' % locus)
            HAPFAS.close()
            self.printLog('#UNZIP','Unzipped SNP into %s haplotigs (%s phased) for %s loci.' % (rje.iStr(bx),rje.iStr(hx),rje.iStr(seqlist.seqNum())))
            self.printLog('#FAS','%s phased and %s unphased blocks output to %s' % (rje.iStr(hx),rje.iStr(cx),bfile))
            self.db().deleteTable('allrid')
            self.db().deleteTable('pileup')

            # Save haplotigs file
            bdb.fillBlanks(blank=0,fields=['MeanX'],fillempty=True)
            for bentry in bdb.indexEntries('Track','C'):
                hap = bentry['Haplotig']
                haplen = bentry['End'] - bentry['Start'] + 1.0
                for hentry in hapdb.indexEntries('Haplotig',hap):
                    bentry['MeanX'] += (hentry['End'] - hentry['Start'] + 1.0) / haplen
            bdb.saveToFile(sfdict={'MeanX':3})
            hapdb.fillBlanks(blank=0.0,fields=['pTrack'],fillempty=True)
            hapdb.saveToFile(sfdict={'pTrack':4})
            #!# NOTE: There is an output in *.haplotigs.tdt without a HapAcc!

            # Read coverage plot
            sam.baseFile('%s.haplotigs' % self.baseFile())
            self.db().baseFile('%s.haplotigs' % self.baseFile())
            hapdb.rename('rid')
            hapdb.renameField('HapAcc','Locus')
            #x#sam.list['CheckFields'] = []
            sam.coverageFromRID(depthplot=True)
            sam.baseFile(self.baseFile())

            # Generate SAMPhase stacked graphics
            if self.getBool('RGraphics'):
                sam.setBool({'RGraphics':True})
                sam.rGraphics('samunzip','"sambase=%s"' % rje.baseFile(self.getStr('Pileup')))

        except: self.errorLog('%s.unzipSeq() error' % self.prog())
#########################################################################################################################
    def report(self):   ### Generates HTML reports of SAMPhaser haplotigs
        '''Generates HTML reports of SAMPhaser haplotigs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.devLog('#RUN','report',debug=False)
            db = self.db()
            self.printLog('#~~#','## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SAMPhaser Report HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##')
            ## ~ [0a] Setup file paths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            hfile = '%s.report.html' % self.baseFile()
            if rje.exists(hfile) and not self.force():
                self.printLog('#HTML','HTML report found: %s' % hfile)
                return True
            # plotbase will point to files within the PAGSAT plots folder
            basename = self.baseFile(strip_path=True)
            plotbase = '%s.haplotigs.Plots/%s' % (self.baseFile(),basename)
            ## ~ [0b] Check for required files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Not all files are listed but these should indicate that things have run correctly.
            wanted = ['haplotigs.tdt']  # Summary delimited text file
            for wext in wanted:
                wfile = '%s.%s' % (self.baseFile(),wext)
                if not os.path.exists(wfile): raise IOError('Cannot find %s!' % wfile)

            #!# Edited to here #!#

            ## ~ [0c] Setup HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Consider making an HTML table where each entry is a section, with fields:
            #!# |-- Section, Order, Link, LinkDesc, HTML
            sectdata = {}   # Start with a secdata{Link, LinkDesc, HTML} dictionary that could become table.dict['Data']

            #i# Each section should have a description in the desc dictionary and a corresponding <a name=""> xref.
            sections = ['Contents','Summary Table','Haplotig Plots']
            #i# links is a dictionary to link the TEXT description to ANAME when the latter is not just the first word of TEXT
            links = {}
            #i# desc is a dictionary of ANAME identifiers and the corresponding mouseover description (title="X")
            #i# These are now laid out in order of appearance:
            desc = {'contents':'Report contents and quick links',
                    'summary':'Summary table of phased haplotigs',
                    'haplotig':'Plots of haplotig read depth and SNPs against parent sequence'
                    }
            #i# hlink is HTML that can be repeated in each section.
            #i# The link identifier should be used with <a name=""> at the start of the section.
            hlink = ['<p>','Quick Links:']
            for section in sections:
                if section in links: link = links[section]
                else: link = string.split(section)[0].lower()
                hlink.append('~ <a href="#%s" title="%s">%s</a>' % (link,desc[link],section))
                sectdata[section] = {'Link':link, 'LinkDesc':desc[link],
                                     'LinkHTML':'<a href="#%s" title="%s">%s</a>' % (link,desc[link],section),
                                     'HTML':[]}  # This will need to be joined at the end.
            hlink += ['</p>','']
            hlink = string.join(hlink,'\n')
            #i# Elsewhere, a link to the top can be provided:
            toplink = '[<a href="#head" title="Return to top of page">Top</a>]'

            ## ~ [0d] Sort out sequence links and PNG files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# Each locus gets its own <A NAME="Locus"> link subsection with its PNG file.
            #i# These should be linked from the summary tables.
            hapdb = db.copyTable(self.db('haplotigs',add=True,mainkeys=['Block','Track']),'HapSummary')
            hapdb.dataFormat({'Block':'int','Start':'int','End':'int','HapStart':'int','HapEnd':'int','nSNP':'int','nRID':'int','MeanX':'num'})
            hapdb.remakeKeys()
            hapfields = hapdb.fields()
            #># Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	HapAcc	nSNP	nRID	MeanX
            tdtitle = {'Locus':'Sequence/Locus Accession Number',
                       'HapAcc':'Sequence name'}
            # Replace certain field values with <a href="#xxx"> links to PNGs
            # Use # for sorting the order
            hapdb.addFields(['PNG','Link'])
            for entry in hapdb.entries():
                entry['PNG'] = '%s.haplotigs.%s.png' % (plotbase,entry['Locus'])
                entry['Link'] = '<a href="#%s" title="%s">%s</a>' % (entry['Locus'],entry['Locus'],entry['Locus'])

            ### ~ [1] ~ Contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Report contents and quick links to sections
            section = 'Contents'
            if section in sections:
                hbody = sectdata[section]['HTML']
                hbody += ['<a name="%s"></a><h2 title="%s">%s %s</h2>' % (sectdata[section]['Link'],sectdata[section]['LinkDesc'],basename,section),hlink,'']

            ### ~ [2] ~ Summary Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Summary table of assembly against reference chromosomes
            section = 'Summary Table'
            if section in sections:
                hbody = sectdata[section]['HTML']
                hbody += ['<a name="%s"></a><h2 title="%s">%s %s</h2>' % (sectdata[section]['Link'],sectdata[section]['LinkDesc'],basename,section),hlink,'']
                htxt = '''
                <p>SAMPhaser haplotig phasing summary. Locus names link to read depth/SNP plots.</p>
                '''
                hbody.append(htxt)
                for entry in hapdb.entries(): [entry['Locus'],entry['Link']] = [entry['Link'],entry['Locus']]
                hbody.append(rje_html.dbTableToHTML(hapdb,hapfields,tdtitles=tdtitle))  # Could modify CSS for tabid formatting?
                for entry in hapdb.entries(): [entry['Locus'],entry['Link']] = [entry['Link'],entry['Locus']]

            ### ~ [3] ~ Haplotig Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #i# Summary plot of assembly against reference chromosomes
            section = 'Haplotig Plots'
            plotlink = '%s [<a href="#haplotig" title="Haplotig Plots">Haplotig Plots</a>]' % toplink
            if section in sections:
                hbody = sectdata[section]['HTML']
                hbody += ['<a name="%s"></a><h2 title="%s">%s %s</h2>' % (sectdata[section]['Link'],sectdata[section]['LinkDesc'],basename,section),hlink,'']
                loci = []
                for ekey in hapdb.dataKeys():
                    entry = hapdb.data(ekey)
                    if entry['Locus'] in loci: continue
                    loci.append(entry['Locus'])
                    lkeys = hapdb.index('Locus')[entry['Locus']]
                    hbody.append('<hr width="80%%"><a name="%s"></a><h3 title="%s">%s ~ %s</h3>' % (entry['Locus'],entry['Locus'],entry['Locus'],plotlink))
                    pngfile = entry['PNG']
                    hbody.append(rje_html.dbTableToHTML(hapdb,hapfields,datakeys=lkeys,tabwidth=1050,tdtitles=tdtitle))
                    hbody.append('<p><code>%s</code></p>' % pngfile)
                    hbody.append('<a href="%s"><img src="%s" width="100%%" title="%s"></a>' % (pngfile,pngfile,'%s haplotigs read depth plot' % entry['Locus']))

            ### ~ [X] Output HTML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            html = rje_html.HTML(self.log,self.cmd_list)    #!# Add default css here. (Put in SLiMSuite? Or Plots/?)
            if rje.exists(hfile) and not self.force() and not assembly:
                self.printLog('#HTML','HTML report found: %s' % hfile)
            else:
                HTML = open(hfile,'w')
                HTML.write(html.htmlHead(title=basename,tabber=False,frontpage=True,keywords=[],redirect='',refresh=0))
                HTML.write('<a name="head"><h1>%s SAMPhaser Report</h1></a>\n\n' % basename)
                HTML.write(rje_html.progStartHTML(self))
                for section in sections: HTML.write(string.join(sectdata[section]['HTML'],'\n'))
                HTML.write(html.htmlTail(tabber=False))
                HTML.close()
                self.printLog('#HTML','HTML report output: %s' % hfile)
            return sectdata
        except: self.errorLog('%s.report error' % self.prog()); return False
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
