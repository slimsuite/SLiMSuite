    #!/usr/bin/python

# See below for name and description
# Copyright (C) 2009 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Module:       rje_aic
Description:  Alternative Initiation Codon Module
Version:      0.8
Last Edit:    30/01/14
Copyright (C) 2009  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module contains a bunch of miscellaneous methods etc. for exploring AIC. set by mode=X. Appropriate modes are
    described in more detail below.

Aug11:
    This mode is initially experimenting with the new rje_seqlist module and whether flanking regions can be extracted
    directly from chromosomal DNA data.

Feb10:
    This mode constitutes a revised Kozak-style analysis of human TIC (Translation Initiation Codon). Data is read in
    from four sources:
    1. Full length EnsEMBL cDNAs (reference only)
    2. 1000nt 5' of TSS. This should not include 5' UTR where annotation is good.
    3. 5'UTR of transcripts (where available).
    4. Coding sequences of each transcript.

    Sequences are compiled and cross-checked with each other. Sequence <3> should match the start of sequence <1> for a
    given transcript and be followed by sequence <4>. Where <3> is missing, the full length cDNA and CDS are used and the
    TIC is assumed to be the first codon of the CDS. In terms of TIC context, this gives different types of sequence
    ('etype'):
    - AL. 5'UTR annotated and longer than kozwin.
    - AS. 5'UTR annotated but shorter than kozwin.
    - US. 5'UTR inferred from cDNA but shorter than kozwin.
    - These are also combined for additional combos (XL,XS,AX,XX) for calculations
    - NL. 5'UTR annotated and longer than kozwin but non-AUG TIC. Not used for calculations (but still rated).
    - NS. 5'UTR annotated but shorter than kozwin, plus non-AUG TIC. Not used for calculations (but still rated).

    For each transcript, the following sequences are then extracted:
    - The "TICwin"  window around TIC (as set by kozwin). [*.kozwin.fas]
    - A control "Conwin" window around the first AUG encountered in seq <2>. [*.conwin.fas]

    Each transcript is then classified according to its core sequence:
    - Good = "best" Kozak gccaugg [set by GoodIC list]
    - augG = Not "Good" but has G at +4
    - Mid (25-70% activity) where the -3/+4 combination is G/A, A/A, G/C, A/C, U/C, G/U, A/U, C/U
    - Weak (<25% activity) where -3/+4 combination is C/C, C/A, U/A, U/U

    For each gene, transcripts are checked for redundancy and only one representative of each different TICwin sequence
    is kept. Regional GC content is estimated from the first 500nt of seq <2>.

    The following calculations are then performed for both TICwin and Conwin sequences:
    - Nucleotide counts for each position of the window, within each eType 
    - Frequency and Ranking for "core" sequences (-3+4) within each eType
    - Frequency "f" scores for each window using sums of observed frequencies at each position

    These calculations are performed for all sequences and then repeated for each CoreType independently. This is output
    in the "CoreSplit" field of each output (All/Good/augG/Rxx/Poor).
    
    The following data is recorded for each transcript:
    - CoreSplit = All/Good/augG/Rxx/Poor.
    - Transcript ID
    - Gene ID
    - TICwin = Window around TIC (as set by kozwin)
    - Core = Core Kozak window -3 to + 4 (XXXAUGX)
    - CoreType = Good/augG/Rxx/Poor
    - WegCore = Wegrzyn Core sequence XXnnXnnAUGX
    - Conwin = Control window around 5' AUG (as set by kozwin)
    - Redundancy = No. filtered identical TIC windows for this gene
    - AIC = No. of *different* TIC windows for this gene.
    - e5utr = UTR length from input seq <3>
    - c5utr = UTR length from input seq <1>
    - eType = AL/AS/US (see above)
    - regG = G count from first 500nt of seq <2>
    - regC = G count from first 500nt of seq <2>
    - regA = A count from first 500nt of seq <2>
    - regT = T count from first 500nt of seq <2>
    - tfX = TICwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
    - cfX = Conwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
    - trX = TICwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
    - crX = Conwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
    - twX = TICwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
    - cwX = Conwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.

    The following data is recorded for each Core:
    - CoreSplit = All/Good/augG/Rxx/Poor.
    - Core = Core sequence
    - tfX = Mean frequency scores (excluding Core) for TICwin with this core
    - tiX = Information content (excluding Core) for TICwin with this core
    - trX = Rank for TICwin with this score
    - cfX = Mean frequency scores (excluding Core) for Conwin with this core
    - ciX = Information content (excluding Core) for Conwin with this core
    - crX = Rank for Conwin with this core

    The following data is recorded for each position in the window:
    - eType = AL/AS/US/XL/XS/AX/XX.
    - CoreSplit = All/Good/augG/Rxx/Poor.
    - G = Count of G
    - A = Count of G
    - T = Count of G
    - C = Count of G
    - Info = Information Content of position.

BIOL3050
    This method analyses the genome of choice for good (and bad) candidates for BIOL3050 project genes. Initiation codons
    are divided into:
    - Weak = -3[UC]NNXXX[UCA]
    - Strong = -3[AG]NNXXXG, where X is AUG/CUG/ACG
    - Mid = All others
    Genes are then classified according to the following criteria for experimentation:
    - Good = Weak annotated IC and Strong eORF or tORF, 40aa <= eLen <= 200aa. No RECuts in 5'UTR to IC+300nt.
    - Cut = Weak annotated IC and Strong eORF or tORF, 40aa <= eLen <= 200aa but RECuts in 5'UTR to IC+300nt.
    - Len = Weak annotated IC and Strong eORF or tORF, 40aa > eLen > 200aa.
    - Poor = All others.

PRIDE
    This is a very basic parser of PRIDE data that extracts the peptides identified from each XML file.

Commandline:
    ### ~ GENERAL ~ ###
    mode=LIST       : Run mode (refcheck/kozak/feb10/pride/jan14) [pride]
    track=LIST      : List of genes or IDs to track through analysis []
    enspath=PATH    : Path to EnsEMBL files [./EnsEMBL/]
    ### ~ REFCHECK ~ ###
    refseq=FILE     : File containing RefSeq download in GenBank format []
    biomart=FILE    : File containing BioMart download []
    shortutr=X      : 5' UTR <X bp will be marked as "Short" [10]
    ### ~ KOZAK ~ ###
    kozak=X         : Basename for Kozak analysis input ['Homo_sapiens.GRCh37.56']
    kozwin=X        : No. of nucleotides either side of start codon [21]
    noutr=T/F       : Whether to include sequences without a 5' UTR [True]
    flank=X         : Length for 5' flanking sequence [1000]
    output=LIST     : List of outputs to generate with this run [all]
    rep=X           : Number of random replicates [1000]
    coreal=LIST     : List of Sequence eTypes for CoreAL analysis ['AL']
    nrgene=X        : Which gene type to use for redundancy removal (Gene/EnsG) ['Gene']
    nonaug=LIST     : List of non-AUG codons to consider in good context [CTG,GTG]
    ### ~ BIOL3050 ~ ###
    recut=LIST      : List of recognition sequences for restriction enzymes ['CTCGAG','GCTAGC']
    ### ~ PRIDE ~ ###
    pridefile=FILE  : Delimited file containing PRIDE IDs and Protein IDs []
    pridepath=PATH  : Path to XML downloads from PRIDE ['./PRIDE_XML/']
       
See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_zen
Other modules needed: None
"""
#    goodic=LIST     : List of Core sequence to be considered "good" ['gccaugg']
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_seq, rje_seqlist, rje_sequence, rje_uniprot, rje_zen
import rje_blast_V1 as rje_blast
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Added run modes and Kozak analysis.
    # 0.2 - Added Feb '10 Kozak+ analysis.
    # 0.3 - Changed to use new downloads.
    # 0.4 - Added *very* basic PRIDE mode.
    # 0.5 - Added BIOL3050 experimental design mode.
    # 0.6 - Added Aug11 mode for clean analysis.
    # 0.7 - Updated the Weak and Mid definitions in the light of more recent experimental data. Added nonaug=LIST.
    # 0.8 - Added jan14 mode for specific AIC paper analysis. In general need of a remake and tidy!
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : General tidy up and documentation. See also GATIC and PATIS!
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyright) = ('RJE_AIC', '0.7', 'November 2013', '2009')
    description = 'Misc Alt Initiation Codon Module'
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
### SECTION II: AIC Class                                                                                               #
#########################################################################################################################
class AIC(rje.RJE_Object):     
    '''
    AIC Class. Author: Rich Edwards (2009).

    Info:str
    - EnsPath = Path to EnsEMBL files [./EnsEMBL/]
    - NRGene = Which gene type to use for redundancy removal (Gene/EnsG) ['Gene']
    - Kozak = Basename for Kozak analysis input (*.flank5.fasta, *.utr5.fasta, *.cds.fasta)
    - PrideFile = Delimited file containing PRIDE IDs and Protein IDs []
    - PridePath = Path to XML downloads from PRIDE ['./PRIDE_XML/']
    - RefSeq = File containing RefSeq download in GenBank format []
    
    Opt:boolean
    - NoUTR = Whether to include sequences without a 5' UTR [True]

    Stat:numeric
    - Flank = Length for 5' flanking sequence [1000]
    - KozWin = No. of nucleotides either side of start codon [9]
    - Rep = Number of random replicates [1000]
    - ShortUTR = 5' UTR <X bp will be marked as "Short" [10]

    List:list
    - CoreAL = List of Sequence eTypes for CoreAL analysis ['AL']
    - GoodIC = List of Core sequence to be considered "good" ['gccaugg']
    - IgnoreFT = List of feature types to ignore when reading in GenBank
    - Mode = Run mode (refcheck/kozak) []
    - NonAUG = List of non-AUG codons to consider in good context [CUG,GUG]
    - Output = List of outputs to generate with this run [enst,core,wegcore,pos]
    - RECut = List of recognition sequences for restriction enzymes ['CTCGAG','GCTAGC']
    - Track = List of genes or IDs to track through analysis []
    
    Dict:dictionary
    - Data = Dictionary of data to store along with each Transcript
    - Freq = Dictionary of Frequencies {Pos:{nt:freq}}
    - RefSeq = Dictionary of {ID:GenBankEntry Object}

    Obj:RJE_Objects
    - SeqList = SeqList object with assembled sequences, start +/- 200.
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.infolist = ['Kozak','RefSeq','BioMart','NRGene','PrideFile','PridePath','EnsPath']
        self.optlist = ['NoUTR']
        self.statlist = ['Flank','KozWin','ShortUTR','Rep']
        self.listlist = ['IgnoreFT','Mode','GoodIC','Output','CoreAL','Track','RECut','NonAUG']
        self.dictlist = ['RefSeq','Data']
        self.objlist = ['SeqList']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.info['RefSeq'] = 'D:\\Projects - Proteomics\\2008-12 AIC Proteomics\\2009-04-06 UTR Assessment\\human_refseq_mrna.gb'
        self.setInfo({'BioMart':'human_utr.csv','Kozak':'ens_HUMAN',#'Basefile':'KozakPlus',
                      'Kozak':'Homo_sapiens.GRCh37.63','NRGene':'Gene','PridePath':rje.makePath('./PRIDE_XML/'),
                      'EnsPath':rje.makePath('./EnsEMBL/')})
        self.setStat({'ShortUTR':10,'KozWin':21,'Flank':1000,'Rep':1000})
        self.setOpt({'NoUTR':True})
        self.list['IgnoreFT'] = ['misc_feature']
        self.list['Mode'] = ['pride']
        self.list['GoodIC'] = ['gccaugg']
        self.list['Output'] = ['all']
        self.list['CoreAL'] = ['AL']
        self.list['RECut'] = ['CTCGAG','GCTAGC']
        self.list['NonAUG'] = ['GTG','CTG']
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
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['BioMart','RefSeq','Kozak','PrideFile'])  
                self._cmdReadList(cmd,'path',['PridePath','EnsPath'])  
                self._cmdReadList(cmd,'info',['NRGene'])
                self._cmdReadList(cmd,'int',['KozWin','ShortUTR','Flank','Rep'])
                self._cmdReadList(cmd,'list',['IgnoreFT','Mode','GoodIC','Output','CoreAL','Track','RECut','NonAUG'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.list['Mode'] = string.split(string.join(self.list['Mode']).lower())
        self.list['NonAUG'] = string.split(string.replace(string.join(self.list['NonAUG']).upper(),'U','T'))
        if 'all' in self.list['Output']: self.list['Output'] = ['enst','core','wegcore','pos','alcore','eorf','types']
        if self.info['Basefile'].lower() in ['','none']: self.info['Basefile'] = self.info['Kozak']
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'jan14' in self.list['Mode']: self.jan14()
            if 'aug11' in self.list['Mode']: self.aug11()
            if 'refcheck' in self.list['Mode']: self.refCheck()
            if 'kozak' in self.list['Mode']: self.kozak()
            if 'feb10' in self.list['Mode']: self.kozak10()
            if 'pride' in self.list['Mode']: self.pride()
            if 'jun10' in self.list['Mode']:
                self.jun10(plus=True); self.jun10('ATG',True); self.jun10();
                self.jun10('ATG')
                for nonaug in self.list['NonAUG']: self.jun10(nonaug)
            if 'biol3050' in self.list['Mode']: self.biol3050()
        except: self.errorLog(rje_zen.Zen().wisdom()); raise
#########################################################################################################################
    ### <3> ### RefCheck Methods                                                                                        #
#########################################################################################################################
    def refCheck(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setupRefCheck()
            ### ~ [2] ~ Generate Max UTR dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            max_utr = {}    # Dictionary of AccNum:Max5'UTR
            (gx,gtot) = (0.0,len(self.dict['RefSeq']))
            for id in rje.sortKeys(self.dict['RefSeq']):
                self.progLog('\r#UTR','Calculating max 5\' UTR: %.2f%%' % (gx/gtot)); gx += 100.0
                gb = self.dict['RefSeq'][id]
                max_utr[id] = -1
                for ft in gb.list['Feature']:
                    if ft['Type'].upper() == 'CDS': max_utr[id] = max(max_utr[id],ft['Start']-1)
            self.printLog('\r#UTR','Calculated max 5\' UTR for %s entries.' % rje.integerString(gtot))
            ### ~ [3] ~ Report ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            utr = {-1:0,1:0,'Short':0,'OK':0}
            ufile = rje.baseFile(self.info['RefSeq']) + '.utr.tdt'
            uhead = ['SeqID','MaxUTR']
            rje.delimitedFileOutput(self,ufile,uhead,rje_backup=True)
            for id in rje.sortKeys(max_utr):
                rje.delimitedFileOutput(self,ufile,uhead,datadict={'SeqID':id,'MaxUTR':max_utr[id]})
                if max_utr[id] == -1: utr[-1] += 1
                else:
                    if max_utr[id] == 0: utr[1] += 1
                    if max_utr[id] < self.stat['ShortUTR']: utr['Short'] += 1
                    else: utr['OK'] += 1
            self.printLog('#TDT','Max 5\' UTR lengths saved to %s' % ufile)
            self.printLog('#UTR','%s entries missing CDS features.' % rje.integerString(utr[-1]))
            self.printLog('#UTR','%s entries have no 5\' UTR.' % rje.integerString(utr[1]))
            self.printLog('#UTR','%s entries have short 5\' UTR (< %d bp).' % (rje.integerString(utr['Short']),self.stat['ShortUTR']))
            self.printLog('#UTR','%s entries have OK 5\' UTR.' % rje.integerString(utr['OK']))
            ### ~ [4] ~ Generate Max UTR from BioMart too ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ensg_list = []      # List of EnsEMBL genes - added to max_utr
            ## ~ [4a] ~ Concatenate UTR for transcripts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            (bx,btot) = (0.0,len(self.dict['BioMart']))
            for bkey in rje.sortKeys(self.dict['BioMart']):
                self.progLog('\r#MART','Converting BioMart: %.2f%%' % (bx/btot)); bx += 50.0
                bdat = self.dict['BioMart'].pop(bkey)
                tran = bdat['Ensembl Transcript ID']
                try:
                    ulen = string.atoi(bdat["5' UTR End"]) - string.atoi(bdat["5' UTR Start"]) + 1
                    if tran in self.dict['BioMart']: self.dict['BioMart'][tran]['ULen'] += ulen
                    else: self.dict['BioMart'][tran] = {'Gene':bdat['Ensembl Gene ID'], 'ULen':ulen}
                except: 
                    if tran not in self.dict['BioMart']: self.dict['BioMart'][tran] = {'Gene':bdat['Ensembl Gene ID'], 'ULen':0}
            ## ~ [4b] ~ Convert transcripts to Genes (Max UTR per gene) ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            btot += len(self.dict['BioMart'])
            for tran in rje.sortKeys(self.dict['BioMart']):
                self.progLog('\r#MART','Converting BioMart: %.2f%%' % (bx/btot)); bx += 50.0
                gene = self.dict['BioMart'][tran]['Gene']
                if gene not in ensg_list: ensg_list.append(gene); max_utr[gene] = self.dict['BioMart'][tran]['ULen']
                else: max_utr[gene] = max(max_utr[gene],self.dict['BioMart'][tran]['ULen'])
            self.printLog('\r#MART','Converted BioMart for %s genes.' % rje.integerString(len(ensg_list)))
            ensg_list.sort()
            ### ~ [4c] ~ BioMart UTR output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            utr = {1:0,'Short':0,'OK':0}
            ufile = rje.baseFile(self.info['BioMart']) + '.utr.tdt'
            uhead = ['SeqID','MaxUTR']
            rje.delimitedFileOutput(self,ufile,uhead,rje_backup=True)
            for id in ensg_list:
                rje.delimitedFileOutput(self,ufile,uhead,datadict={'SeqID':id,'MaxUTR':max_utr[id]})
                if max_utr[id] == 1: utr[1] += 1
                if max_utr[id] < self.stat['ShortUTR']: utr['Short'] += 1
                else: utr['OK'] += 1
            self.printLog('#TDT','Max 5\' UTR lengths saved to %s' % ufile)
            self.printLog('#UTR','%s entries have no 5\' UTR.' % rje.integerString(utr[1]))
            self.printLog('#UTR','%s entries have short 5\' UTR (< %d bp).' % (rje.integerString(utr['Short']),self.stat['ShortUTR']))
            self.printLog('#UTR','%s entries have OK 5\' UTR.' % rje.integerString(utr['OK']))
            return
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setupRefCheck(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Read GenBank ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.exists(self.info['RefSeq']):
                REF = open(self.info['RefSeq'],'r')
                rline = REF.readline()
                rlines = [rje.chomp(rline)]
                while rline:
                    rline = REF.readline()
                    rlines.append(rje.chomp(rline))
                    if rlines[-1][:2] == '//':      # Complete entry
                        self.addGenBank(rlines)
                        rlines = []
                    self.progLog('\r#REF','Reading RefSeq GenBank data: %s entries' % (rje.integerString(len(self.dict['RefSeq']))))
                self.printLog('\r#REF','Reading RefSeq GenBank data complete: %s entries' % (rje.integerString(len(self.dict['RefSeq']))))
            else: self.printLog('\r#REF','Cannot find RefSeq file "%s"' % self.info['RefSeq'])
            ### ~ [2] ~ Read BioMart Ensembl Gene ID,Ensembl Transcript ID,5' UTR Start,5' UTR End ~~~~~~~~~~~~~~~~~~ ###
            if rje.exists(self.info['BioMart']):
                bkeys = ["Ensembl Gene ID","Ensembl Transcript ID","5' UTR Start","5' UTR End"]
                self.dict['BioMart'] = rje.dataDict(self,self.info['BioMart'],bkeys,bkeys)
            else: self.printLog('\r#ENS','Cannot find BioMart file "%s"' % self.info['BioMart'])
            return True     # Setup successful
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    def addGenBank(self,rlines): ### Add a GenBank entry
        '''Add a GenBank entry.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gb = GenBankEntry(self.log,self.cmd_list)
            gb.process(rlines,ft_only=['cds'])
            id = gb.info['Name']
            self.dict['RefSeq'][id] = gb
            return True            
        except: self.errorLog('Problem during %s addGenBank.' % self); return False  
#########################################################################################################################
    ### <4> ### Kozak analysis prep methods                                                                             #
#########################################################################################################################
    def kozak(self):    ### Generate sequences for Kozak style analysis
        '''Generate sequences for Kozak style analysis.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            kwin = self.stat['KozWin']; ksize = (2 * kwin) + 3
            ## ~ [1a] Load Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqs = {}
            seqcmd = self.cmd_list + ['autoload=T','accnr=F','seqnr=F','gnspacc=F']
            seqs['utr5'] = rje_seq.SeqList(self.log,seqcmd+['seqin=%s.utr5.fasta' % self.info['Kozak']])
            seqs['flank'] = rje_seq.SeqList(self.log,seqcmd+['seqin=%s.flank5.fasta' % self.info['Kozak']])
            seqs['cds'] = rje_seq.SeqList(self.log,seqcmd+['seqin=%s.cds.fasta' % self.info['Kozak']])
            ## ~ [1b] Dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqdict = {}
            for key in ['utr5','flank','cds']:
                sx = 0.0; stot = seqs[key].seqNum()
                for seq in seqs[key].seqs():
                    self.progLog('\r#%s' % key,'Parsing %s %s sequences: %.1f%%' % (rje.integerString(stot),key,sx/stot)); sx += 100.0
                    seqid = string.split(string.split(seq.info['FullName'])[0],'|')
                    if len(seqid) > 1: desc = string.join(seqid[-1:]+string.split(seq.info['FullName'])[1:])
                    else: desc = string.join(string.split(seq.info['FullName'])[1:])
                    for id in seqid[:2]:
                        if id not in seqdict: seqdict[id] = {'desc':desc}
                        if key not in seqdict[id]: seqdict[id][key] = seq.info['Sequence'].lower()
                self.printLog('\r#%s' % key,'Parsed %s %s sequences: %s seq ID' % (rje.integerString(stot),key,rje.integerString(len(seqdict))))
            ### ~ [2] Check and Generate Kozak sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#KOZAK','Generating Kozak %dnt test sequences (win=%d)' % (ksize,kwin),log=False)
            kozak = {}
            ix = 0; itot = len(seqdict); bx = 0; d = self.opt['DeBug']
            for id in rje.sortKeys(seqdict):
                ## ~ [2a] Check for all information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#KOZAK','Processing %s IDs %s >> %s kozak test IDs' % (rje.integerString(itot),rje.integerString(ix),rje.integerString(len(kozak)))); ix += 1
                skip = False
                #self.deBug('%s: %s' % (id, rje.sortKeys(seqdict[id])))
                for key in ['utr5','flank','cds']:
                    if key not in seqdict[id]: skip = True
                if skip: continue
                self.opt['DeBug'] = d
                #for key in ['flank','utr5','cds']: self.deBug('%s = %s' % (key,seqdict[id][key]))
                ## ~ [2b] Check UTR vs CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if seqdict[id]['cds'].find(seqdict[id]['utr5']): continue            # CDS should start with UTR
                ## ~ [2c] Remove UTR from CDS and add UTR3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqdict[id]['cds'] = seqdict[id]['cds'][len(seqdict[id]['utr5']):]   # Remove UTR from CDS
                key = 'cds'; self.deBug('%s = %s' % (key,seqdict[id][key]))
                if seqdict[id]['cds'][:3] != 'atg': continue                        # CDS should start with ATG
                i = 0
                while seqdict[id]['cds'][i:i+3] not in ['tga','taa','tag','']: i += 3
                i +=3
                seqdict[id]['utr3'] = seqdict[id]['cds'][i:]
                seqdict[id]['cds'] = seqdict[id]['cds'][:i]
                key = 'cds'; self.deBug('%s = %s' % (key,seqdict[id][key]))
                key = 'utr3'; self.deBug('%s = %s' % (key,seqdict[id][key]))
                ### ~ [2d] Generate Kozak test sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if len(seqdict[id]['utr5']) < kwin or len(seqdict[id]['cds']) < (kwin + 3): continue
                kozak[id] = {'start':seqdict[id]['utr5'][-kwin:]+seqdict[id]['cds'][:kwin+3]}
                i = 3
                while i < len(seqdict[id]['cds']) and seqdict[id]['cds'][-i:][:3] != '':
                    if seqdict[id]['cds'][-i:][:3] != 'atg': i += 3; continue
                    testseq = seqdict[id]['cds'][-i-kwin:][:ksize]
                    if len(testseq) == ksize:
                        kozak[id]['cds'] = testseq
                        break
                    i += 3
                for type in ['flank','utr5','utr3']:
                    i = 0
                    while seqdict[id][type][i:i+3] != '':
                        if seqdict[id][type][i:i+3] != 'atg': i += 1; continue
                        testseq = seqdict[id][type][i-kwin:i+3+kwin]
                        if len(testseq) == ksize:
                            kozak[id][type] = testseq
                            break
                        i += 1
                if len(kozak[id]) == 5 or (len(kozak[id]) == 4 and 'utr5' not in kozak[id]): bx += 1
            self.printLog('\r#KOZAK','Processed %s IDs >> %s kozak test IDs; %s "balanced"' % (rje.integerString(itot),rje.integerString(len(kozak)),rje.integerString(bx)))
            ### ~ [3] Generate output files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for type in ['start','cds','flank','utr5','utr3']:
                open('%s.ktest.%d.%s.fas' % (self.info['Kozak'],kwin,type),'w')
                open('%s.ktest.%d.%s.txt' % (self.info['Kozak'],kwin,type),'w')
                open('%s.ktest.%d.balanced.%s.fas' % (self.info['Kozak'],kwin,type),'w')
                open('%s.ktest.%d.balanced.%s.txt' % (self.info['Kozak'],kwin,type),'w')
            ix = 0.0; itot = len(kozak); bx = 0
            for id in rje.sortKeys(kozak):
                self.progLog('\r#OUT','Ouput to files: %.1f%% (%s balanced)' % (ix/itot,rje.integerString(bx))); ix += 100.0
                balanced = True
                for type in ['start','cds','flank','utr5','utr3']:
                    if type in kozak[id]:
                        open('%s.ktest.%d.%s.fas' % (self.info['Kozak'],kwin,type),'a').write('>%s %s\n%s\n' % (id,seqdict[id]['desc'],kozak[id][type]))
                        open('%s.ktest.%d.%s.txt' % (self.info['Kozak'],kwin,type),'a').write('%s\n' % kozak[id][type])
                    elif type != 'utr5': balanced = False
                if balanced:
                    bx += 1
                    for type in ['start','cds','flank','utr3']:
                        open('%s.ktest.%d.balanced.%s.fas' % (self.info['Kozak'],kwin,type),'a').write('>%s %s\n%s\n' % (id,seqdict[id]['desc'],kozak[id][type]))
                        open('%s.ktest.%d.balanced.%s.txt' % (self.info['Kozak'],kwin,type),'a').write('%s\n' % kozak[id][type])
            self.printLog('\r#OUT','Ouput to files %s ID (%s balanced)' % (rje.integerString(itot),rje.integerString(bx)))
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    ### <5> ### Feb 10 Kozak analysis methods                                                                           #
#########################################################################################################################
    def kozak10_pickle(self):   ### Check for and load pickles
        '''Check for and load pickles.'''
        try:
            self.info['Pickle'] = 'none'
            for p in ['score','counts','loaded']:
                self.printLog('\r#PICK','Checking for %s pickle.' % p,log=False)
                try: pickle = self.unpickleMe('%s.%s' % (self.info['Basefile'],p))
                except: pickle = None
                if pickle:
                    self.obj['SeqList'] = pickle.obj['SeqList']
                    self.dict['Data'] = pickle.dict['Data']
                    self.dict['AIC'] = pickle.dict['AIC']
                    self.info['Pickle'] = p
                    if p == 'score': self.dict['CoreData'] = pickle.dict['CoreData']
                    return pickle
        except: pass
        return None
#########################################################################################################################
    def kozak10(self):    ### Perform revised Feb '10 Kozak style analysis
        '''This mode constitutes a revised Kozak-style analysis of human TIC (Translation Initiation Codon).'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['GoodIC'] = string.split(string.replace(string.join(self.list['GoodIC']).upper(),'U','T'))
            #unpickleMe('%s.loaded' % self.info['Basefile'])
            #if pickle: seqlist = pickle.obj['SeqList']
            #else: seqlist = self._kozak10_loadSeq()
            #self.obj['SeqList'] = seqlist
            pickle = self.kozak10_pickle()
            if not self.obj['SeqList']: self._kozak10_loadSeq()
            if self.list['Output']: self._kozak10_headers()
            self._kozak10_outputTypes()
            ### ~ [2] ~ Perform Calculations for both TICwin and Conwin sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Nucleotide and Core sequence counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #pickle = self.unpickleMe('%s.counts' % self.info['Basefile'])
            #if pickle: self.dict['Data'] = pickle.dict['Data']
            #else: self.dict['Data'] = self._kozak10_counts()
            if not self.dict['Data']: self._kozak10_counts()
            ## ~ [2b] ~ Rank Cores within each eType & cType & score each sequence ~~~~~~~~~~~~~~~~ ##
            #pickle = self.unpickleMe('%s.score' % self.info['Basefile'])
            #if pickle:
            #    self.obj['SeqList'] = pickle.obj['SeqList']
            #    self.dict['Data'] = pickle.dict['Data']
            #else: self.obj['SeqList'] = self._kozak10_scoring()
            if self.info['Pickle'] not in ['score']: self._kozak10_scoring()
            ## ~ [2c] ~ eORF analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'eorf' in self.list['Output']: self._kozak10_eORFs()
            ## ~ [2d] ~ AL core analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'alcore' in self.list['Output']: self._kozak10_randCoreAL()
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'enst' in self.list['Output']: self._kozak10_outputENST()
            if 'core' in self.list['Output']: self._kozak10_outputCore('core')
            if 'wegcore' in self.list['Output']: self._kozak10_outputCore('wegcore')
            if 'alcore' in self.list['Output']: self._kozak10_outputALCore()
            if 'pos' in self.list['Output']: self._kozak10_outputPos()

        except: self.errorLog('Feb-2010-Kozak Error')
#########################################################################################################################
    def _kozak10_eORFs(self):   ### Extended 5' ORF analysis
        '''Extended 5' ORF analysis.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            maindata = self.dict['Data']['All']['tic']
            etdt = '%s.eorf.tdt' % self.info['Basefile']
            ehead = ['eORF','eLen','eUTR','Score','Known','ENST','eType','CoreType','NewCoreType']
            for h in ['Core','NewCore','WegCore','NewWeg']: ehead.append(h); ehead.append('%sRank' % h)
            ehead += ['EnsG','Gene','Description']
            rje.delimitedFileOutput(self,etdt,ehead,rje_backup=True)
            efile = '%s.eorf.fas' % self.info['Basefile']
            exfile = '%s.eorf-extra.fas' % self.info['Basefile']
            open(efile,'w'); open(exfile,'w')
            edata = {}
            seqlist = self.obj['SeqList']
            minext = 16     # Minimum extension for ORF
            ## ~ [1a] ~ Make known ORF dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            known = {}; sx = 0.0; stot = seqlist.seqNum(); kx = 0
            for seq in seqlist.seq:
                self.progLog('\r#ORF','Assembling known ORFs: %.2f%%' % (sx/stot)); sx += 100.0
                gene = seq.info['Gene']
                orf = rje_sequence.dna2prot(seq.info['CDS'])
                if gene not in known: known[gene] = []
                if orf not in known[gene]: known[gene].append(orf); kx += 1
                if gene in self.list['Track']:
                    self.printLog('#TRACK','%s -> %s orfs' % (gene,len(known[gene])))
                    if seq.info['AccNum'] not in self.list['Track']: self.list['Track'].append(seq.info['AccNum'])
            self.printLog('\r#ORF','Assembled %s known ORFs for %s genes.' % (rje.integerString(kx),rje.integerString(len(known))))
            ### ~ [2] ~ Check extended ORF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = seqlist.seqNum(); ex = 0; tx = 0; exx = 0; badcore = []
            for seq in seqlist.seq:
                acc = seq.info['AccNum']; gene = seq.info['Gene']
                eorf = 0
                self.progLog('\r#EORF','Extending ORFs 5\': %.2f%%' % (sx/stot)); sx += 100.0
                if gene in self.list['Track']: self.deBug('\n >>> eORF | %s <<<' % seq.info['Gene'])
                upstream = seq.info['Flank'] + seq.info['5utr']
                i = 1
                while ((i*3)+7) < len(upstream) and upstream[-3*i:][:3] not in ['TAA','TGA','TAG','']:
                    if acc in self.list['Track']: self.bugPrint('%s\n%s\n' % (upstream[-3*i:],rje_sequence.dna2prot(upstream[-3*i:])))
                    if 'N' in upstream[-3*i:][:3]: break
                    newstart = False
                    if upstream[-3*i:][:3] == 'ATG': newstart = True
                    if upstream[-3*(i+1):][:1] in ['G','A'] and upstream[-3*i:][:4] in self.list['NonAUG']: newstart = True
                    if newstart:
                        eorf += 1
                        prot = rje_sequence.dna2prot(upstream[-3*i:]+seq.info['CDS'])
                        newcore = upstream[(-3*i)-3:][:7]
                        newweg = upstream[(-3*i)-7:][:11]
                        if i == 1: newcore += seq.info['CDS'][0]; newweg += seq.info['CDS'][0]
                        newweg = rje.strSub(newweg,2,3,'nn')
                        newweg = rje.strSub(newweg,5,6,'nn')
                        edata = {'eORF':'%s-e%d' % (acc,eorf),'ENST':acc,'eLen':i,'eUTR':len(seq.info['5utr'])-(i*3)}
                        edata['NewCore'] = newcore
                        edata['NewCoreType'] = self._coreType(newcore)
                        edata['NewWeg'] = newweg
                        edata['Score'] = 0
                        seq.stat['eORF'] = eorf
                        desc = '%s|%s|%s|e-%s|%s' % (seq.info['Gene'],seq.info['eType'],seq.info['CoreType'],edata['NewCoreType'],seq.info['Description'])
                        if prot in known[seq.info['Gene']]: edata['Known'] = 'Known'
                        else:
                            edata['Known'] = 'Novel'
                            newcds = rje_sequence.dna2prot(upstream[-3*i:]+seq.info['CDS'][:21])
                            newcds50 = rje_sequence.dna2prot(upstream[-3*i:]+seq.info['CDS'][:21])[:50]
                            for cds in known[seq.info['Gene']]:
                                found = cds.find(newcds)
                                found50 = cds.find(newcds50)
                                if found == 0: edata['Known'] = 'Trunc'; break
                                elif found > 0: edata['Known'] = 'Internal'
                                elif found50 == 0 and edata['Known'] not in ['Internal']: edata['Known'] = 'Trunc50'
                                elif found50 > 0 and edata['Known'] not in ['Internal','Trunc50']: edata['Known'] = 'Int50'
                        try: edata['CoreRank'] = maindata['kozrank']['AL'][self.coreXXX(seq.info['Core'],'ATG')]
                        except:
                            if self.coreXXX(seq.info['Core'],'ATG') not in badcore:
                                badcore.append(self.coreXXX(seq.info['Core'],'ATG'))
                                self.printLog('\r#BAD','Cannot rank core "%s"' % self.coreXXX(seq.info['Core'],'ATG'))
                                #self.deBug(edata)
                        try: edata['WegCoreRank'] = maindata['wegrank']['AL'][self.coreXXX(seq.info['WegCore'],'ATG')]
                        except:
                            if self.coreXXX(seq.info['WegCore'],'ATG') not in badcore:
                                badcore.append(self.coreXXX(seq.info['WegCore'],'ATG'))
                                self.printLog('\r#BAD','Cannot rank core "%s"' % self.coreXXX(seq.info['WegCore'],'ATG'))
                                #self.deBug(edata)
                        try: edata['NewCoreRank'] = maindata['kozrank']['AL'][newcore]
                        except:
                            if newcore not in badcore:
                                badcore.append(newcore)
                                self.printLog('\r#NEW','Cannot rank new core "%s"' % newcore)
                                #self.deBug(edata)
                        try: edata['NewWegRank'] = maindata['wegrank']['AL'][newweg]
                        except:
                            if newweg not in badcore:
                                badcore.append(newweg)
                                self.printLog('\r#NEW','Cannot rank new core "%s"' % newweg)
                                #self.deBug(edata)
                        try: edata['Score'] += 500 - (edata['NewWegRank'] - edata['WegCoreRank'])
                        except: pass
                        try: edata['Score'] += 500000 - (1000 * (edata['NewCoreRank'] - edata['CoreRank']))
                        except: pass
                        if (i*3) < len(seq.info['5utr']):
                            name = '%s-e%d (+%daa UTR %s) %s' % (acc,eorf,i,edata['Known'],desc)
                            if i >= minext and edata['Known'] == 'Novel':
                                open(exfile,'a').write('>%s\n%s%s\n' % (name,prot[:i].lower(),prot[i:])); exx += 1
                        else:
                            name = '%s-e%d (+%daa UTR+Flank %s) %s' % (acc,eorf,i,edata['Known'],desc)
                            edata['Score'] = 0 - edata['Score']
                        open(efile,'a').write('>%s\n%s%s\n' % (name,prot[:i].lower(),prot[i:])); ex += 1
                        rje.delimitedFileOutput(self,etdt,ehead,datadict=rje.combineDict(edata,seq.info))
                        if acc in self.list['Track']: self.deBug('\n >>> %s' % name)
                    i += 1
                if acc in self.list['Track']: self.deBug('!%s\n!%s\n' % (upstream[-3*i:],rje_sequence.dna2prot(upstream[-3*i:])))
                if eorf:
                    desc = '%s|%s|%s|%s' % (seq.info['Gene'],seq.info['eType'],seq.info['CoreType'],seq.info['Description'])
                    cds = rje_sequence.dna2prot(seq.info['CDS'])
                    open(efile,'a').write('>%s %s\n%s\n' % (acc,desc,cds))
                if acc in self.list['Track']:
                    self.deBug('\n >>> eORF | %s -> %d eORF <<<' % (seq.info['Gene'],eorf))
                    #self.deBug(seq.attDetails())
            ### ~ [3] ~ Check truncated ORF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                cds = seq.info['CDS']
                i = 1; torf = 0
                while ((i*3)+self.stat['KozWin']) < len(cds) and cds[3*i:][:3] not in ['TAA','TGA','TAG','']:
                    newstart = False
                    if cds[3*(i-1):][:1] in ['G','A'] and cds[3*(i+1):][:1] == 'G' and cds[3*i:][:3] in ['ATG'] + self.list['NonAUG']: newstart = True
                    if newstart:
                        torf += 1
                        prot = rje_sequence.dna2prot(cds[3*i:])
                        newcore = cds[(3*i)-3:][:7]
                        if i > 2: newweg = cds[(3*i)-7:][:11]
                        else:
                            newweg = cds[:(3*i)+4]
                            if (len(newweg)) < 11: newweg = upstream[-(11-len(newweg)):] + newweg
                        newweg = rje.strSub(newweg,2,3,'nn')
                        newweg = rje.strSub(newweg,5,6,'nn')
                        edata = {'eORF':'%s-t%d' % (acc,torf),'ENST':acc,'eLen':i,'eUTR':len(seq.info['5utr'])+(i*3)}
                        edata['NewCore'] = newcore
                        edata['NewCoreType'] = self._coreType(newcore)
                        edata['NewWeg'] = newweg
                        edata['Score'] = 0
                        seq.stat['eORF'] = torf
                        desc = '%s|%s|%s|e-%s|%s' % (seq.info['Gene'],seq.info['eType'],seq.info['CoreType'],edata['NewCoreType'],seq.info['Description'])
                        if prot in known[seq.info['Gene']]: edata['Known'] = 'Known'
                        else:
                            edata['Known'] = 'Internal'
                            newcds50 = prot[:50]
                            for gcds in known[seq.info['Gene']]:
                                found50 = gcds.find(newcds50)
                                if found50 == 0: edata['Known'] = 'Trunc'
                        try: edata['CoreRank'] = maindata['kozrank']['AL'][self.coreXXX(seq.info['Core'],'ATG')]
                        except:
                            if self.coreXXX(seq.info['Core'],'ATG') not in badcore:
                                badcore.append(self.coreXXX(seq.info['Core'],'ATG'))
                                self.printLog('\r#BAD','Cannot rank core "%s"' % self.coreXXX(seq.info['Core'],'ATG'))
                                #self.deBug(edata)
                        try: edata['WegCoreRank'] = maindata['wegrank']['AL'][self.coreXXX(seq.info['WegCore'],'ATG')]
                        except:
                            if self.coreXXX(seq.info['WegCore'],'ATG') not in badcore:
                                badcore.append(self.coreXXX(seq.info['WegCore'],'ATG'))
                                self.printLog('\r#BAD','Cannot rank core "%s"' % self.coreXXX(seq.info['WegCore'],'ATG'))
                                #self.deBug(edata)
                        try: edata['NewCoreRank'] = maindata['kozrank']['AL'][newcore]
                        except:
                            if newcore not in badcore:
                                badcore.append(newcore)
                                self.printLog('\r#NEW','Cannot rank trunc core "%s"' % newcore)
                                #self.deBug(edata)
                        try: edata['NewWegRank'] = maindata['wegrank']['AL'][newweg]
                        except:
                            if newweg not in badcore:
                                badcore.append(newweg)
                                self.printLog('\r#NEW','Cannot rank trunc core "%s"' % newweg)
                                #self.deBug(edata)
                        try: edata['Score'] += 500 - (edata['NewWegRank'] - edata['WegCoreRank'])
                        except: pass
                        try: edata['Score'] += 500000 - (1000 * (edata['NewCoreRank'] - edata['CoreRank']))
                        except: pass
                        prot = rje_sequence.dna2prot(cds)
                        name = '%s-t%d (-%daa CDS %s) %s' % (acc,torf,i,edata['Known'],desc)
                        if i >= minext: open(exfile,'a').write('>%s\n%s%s\n' % (name,prot[:i].lower(),prot[i:])); exx += 1
                        open(efile,'a').write('>%s\n%s%s\n' % (name,prot[:i].lower(),prot[i:])); tx += 1
                        rje.delimitedFileOutput(self,etdt,ehead,datadict=rje.combineDict(edata,seq.info))
                        if acc in self.list['Track']: self.deBug('\n >>> %s' % name)
                    i += 1
                if acc in self.list['Track']: self.deBug('!%s\n!%s\n' % (upstream[-3*i:],rje_sequence.dna2prot(upstream[-3*i:])))
                if torf or eorf:
                    desc = '%s|%s|%s|%s' % (seq.info['Gene'],seq.info['eType'],seq.info['CoreType'],seq.info['Description'])
                    cds = rje_sequence.dna2prot(seq.info['CDS'])
                    open(efile,'a').write('>%s %s\n%s\n' % (acc,desc,cds))
                if acc in self.list['Track']:
                    self.deBug('\n >>> eORF | %s -> %d eORF <<<' % (seq.info['Gene'],eorf))
                    #self.deBug(seq.attDetails())
            self.printLog('\r#EORF','Extending ORFs 5\' complete: %d e-ORF; %d t-ORF & %d e-ORF-extra.' % (ex,tx,exx))
        except: self.errorLog('Feb-2010-Kozak eORF Error')
#########################################################################################################################
    def _kozak10_eORF_Summary(self):   ### Extended 5' ORF analysis summary
        '''Extended 5' ORF analysis summary.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            etdt = '%s.eorf.tdt' % self.info['Basefile']
            eorf = rje.dataDict(self,etdt)
            esum = {}
            for c1 in ['Good','augG','Mid','Weak']:
                esum[c1] = {}
                for c2 in ['Good','augG','Mid','Weak']:
                    esum[c1][c2] = []

            #!# Think about this and finish it off #!#            

            ehead = ['eORF','eLen','eUTR','ENST','eType','CoreType']
            for h in ['Core','NewCore','WegCore','NewWeg']: ehead.append(h); ehead.append('%sRank' % h)
            ehead += ['EnsG','Gene','Description']
        except: self.errorLog('Feb-2010-Kozak eORF Summary Error')
#########################################################################################################################
    def _kozak10_randCoreAL(self):  ### Calculates probabilites of ranked cores for AL sequences only
        '''Calculates probabilites of ranked cores for AL sequences only.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            d = self.dict['Data']
            k = self.stat['KozWin']
            seqlist = self.obj['SeqList']
            cdata = self.dict['CoreData']
            rep = self.stat['Rep']  #!# Add commandline option? #!#
            ## ~ [1a] ~ Setup sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            mingatc = 100   # Min. number of non-N nucleotides #        #!# Add commandline option? #!#
            alseq = []; goodal = 0; badal = 0
            sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seq:
                sx += 100.0
                if seq.info['eType'] not in self.list['CoreAL']: continue   #!= 'AL': continue
                self.progLog('\r#AL','Setting up AL sequences: %.2f%%' % (sx/stot))
                seq.dict['GATC'] = {}
                for n in 'GATC':
                    try: seq.dict['GATC'][n] = seq.dict['GC'][n]
                    except: seq.dict['GATC'][n] = 0
                rje.dictFreq(seq.dict['GATC'],total=True)
                alseq.append(seq)
                if seq.dict['GATC']['Total'] < mingatc: badal += 1
                else: goodal += 1
            self.printLog('\r#AL','Setup of %s AL sequences complete.' % rje.integerString(len(alseq)))
            self.printLog('\r#ALNT','%s AL sequences with sufficient G/A/T/C nt. %s without.' % (rje.integerString(goodal),rje.integerString(badal)))
            if not goodal: raise ValueError
            ## ~ [1b] ~ Check for Pickle of CoreAL data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pbase = '%s.alcore.%s.%s' % (self.info['Basefile'],goodal,rep)
            alpickle = self.unpickleMe(pbase)
            if alpickle: self.dict['CoreAL'] = alpickle.dict['CoreAL']; return
            ## ~ [1c] ~ Setup AL core data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            al = self.dict['CoreAL'] = {}       # Dictionary of data for AL cores only
            alists = {}
            ax = 0.0; atot = 5 * 5
            for c in self._coreType():
                self.progLog('\r#AL','Setting up AL cores: %.2f%%' % (ax/atot)); ax += 100.0
                al[c] = {}
                alists[c] = {}
                for x in ['tk','ck','tw','cw']: alists[c][x] = []
                for x1 in 'GATC':
                    for x2 in 'GATC':
                        for x3 in 'GATC':
                            for x4 in 'GATC':
                                kozcore = '%s%s%sATG%s' % (x1,x2,x3,x4)
                                wegcore = '%s%snn%snnATG%s' % (x1,x2,x3,x4)
                                al[c][kozcore] = cdata[c][kozcore]['AL']
                                al[c][wegcore] = cdata[c][wegcore]['AL']
                                alists[c]['tk'].append(al[c][kozcore]['tn'])
                                alists[c]['ck'].append(al[c][kozcore]['cn'])
                                alists[c]['tw'].append(al[c][wegcore]['tn'])
                                alists[c]['cw'].append(al[c][wegcore]['cn'])
                for x in ['tk','ck','tw','cw']:
                    self.progLog('\r#AL','Setting up AL cores: %.2f%%' % (ax/atot)); ax += 100.0
                    alists[c][x].sort()
                    alists[c][x].reverse()
                    for i in range(len(alists[c][x])):
                        try: alists[c][x][i] = (float(alists[c][x][i])/float(sum(alists[c][x][i:])))
                        except: alists[c][x][i] = 0.0
            self.printLog('\r#AL','Setting up of AL cores complete.')

            ### ~ [2] ~ Generate random core count lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            randranks = []; random.seed(rep)
            sx = 0.0; stot = rep * len(alseq)
            for r in range(rep):   # Generate 1000 random samples
                rcore = {}      # Random Core data
                for x1 in 'GATC':
                    for x2 in 'GATC':
                        for x3 in 'GATC':
                            for x4 in 'GATC':
                                rcore['%s%s%s%s' % (x1,x2,x3,x4)] = 0
                ## ~ [2a] ~ Generate random cores based on nt freqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for seq in alseq:
                    sx += 100.0
                    while seq.dict['GATC']['Total'] < mingatc: seq = random.choice(alseq)   # Pick random sequence if too many Ns
                    if seq.info['eType'] != 'AL': continue
                    self.progLog('\r#RAND','Generating random core data (%d): %.2f%%' % (r+1,sx/stot))
                    core = ''
                    for x in range(4):
                        p = random.random()
                        for n in 'GATCG':
                            if p <= seq.dict['GATC'][n]: core += n; break
                            p -= seq.dict['GATC'][n]
                    if len(core) < 4:
                        self.printLog('\r#ERR','Problem generating random core for %s: "%s"' % (seq.info['AccNum'],core),screen=False)
                        #self.deBug(seq.info)
                        #self.deBug(seq.dict['GATC'])
                        continue
                    rcore[core] += 1
                ## ~ [2b] ~ Rank counts and add to randranks list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rrank = rcore.values()
                rrank.sort()
                rrank.reverse()
                randranks.append(rrank[0:])
            self.printLog('\r#RAND','%s x Random core data generation complete.' % rep)
            ### ~ [3] ~ Assess actual core ranks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [3a] ~ Convert to frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            randfreqs = []
            for rrank in randranks:
                rfreq = []
                for i in range(len(rrank)):
                    try: rfreq.append(float(rrank[i])/float(sum(rrank[i:])))
                    except: rfreq.append(0.0)
                randfreqs.append(rfreq[0:])
            ## ~ [3b] ~ Assess actual cores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cx = 0.0; ctot = 5 * len(al['All'])
            for c in self._coreType():
                rje.valueSortedKeys(al[c],rev=True)  #!#
                for core in al[c]:      # (nerps alread)    # o = over-representation; f = over-freq
                    self.progLog('\r#CORE','Assessment of core ranks using randomised data: %.2f%%' % (cx/ctot)); cx += 100.0
                    al[c][core]['to'] = 0.0
                    al[c][core]['tf'] = 0.0
                    al[c][core]['co'] = 0.0
                    al[c][core]['cf'] = 0.0
                    al[c][core]['CoreSplit'] = c
                    al[c][core]['Core'] = core
                    for r in range(rep):   # 1000 random samples
                        try:
                            if randranks[r][al[c][core]['tr']-1] >= al[c][core]['tn']: al[c][core]['to'] += 1
                            if randranks[r][al[c][core]['cr']-1] >= al[c][core]['cn']: al[c][core]['co'] += 1
                        except:
                            #self.deBug('%s|%s' %(r,core))
                            #self.deBug(randranks[r])
                            #self.deBug(al[c][core])
                            continue
                        if len(core) == 7:
                            if randfreqs[r][al[c][core]['tr']-1] >= alists[c]['tk'][al[c][core]['tr']-1]: al[c][core]['tf'] += 1
                            if randfreqs[r][al[c][core]['cr']-1] >= alists[c]['ck'][al[c][core]['cr']-1]: al[c][core]['cf'] += 1
                            al[c][core]['CoreType'] = 'Kozak'
                        else:
                            if randfreqs[r][al[c][core]['tr']-1] >= alists[c]['tw'][al[c][core]['tr']-1]: al[c][core]['tf'] += 1
                            if randfreqs[r][al[c][core]['cr']-1] >= alists[c]['cw'][al[c][core]['cr']-1]: al[c][core]['cf'] += 1
                            al[c][core]['CoreType'] = 'Weg'
                    al[c][core]['to'] /= rep
                    al[c][core]['tf'] /= rep
                    al[c][core]['co'] /= rep
                    al[c][core]['cf'] /= rep
            self.printLog('\r#CORE','Assessment of core ranks using randomised data complete.')
            ## ~ [3c] ~ Save ALCore dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            alpickle = AIC(self.log,self.cmd_list)
            alpickle.dict['CoreAL'] = self.dict['CoreAL']
            alpickle.pickleMe(pbase)

        except: self.errorLog('Feb-2010-Kozak Random AL Core Scoring Error')
#########################################################################################################################
    def _kozak10_scoring(self):     ### Ranks cores and the score each sequence
        '''Ranks cores and the score each sequence.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            d = self.dict['Data']
            k = self.stat['KozWin']
            seqlist = self.obj['SeqList']
            db = self.opt['DeBug']
            ## ~ [1a] ~ Convert nucleotide counts to frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('#FREQ','Converting nucleotide counts to frequencies ...')
            for c in self._coreType():
                for t in ['tic','con']:
                    d[c][t]['freq'] = {}
                    for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                        d[c][t]['freq'][e] = {}
                        for i in range(self.stat['KozWin']):
                            d[c][t]['freq'][e][i-k] = rje.dictFreq(d[c][t]['pos'][e][i-k],total=True,newdict=True)
                            d[c][t]['freq'][e][1+i] = rje.dictFreq(d[c][t]['pos'][e][1+i],total=True,newdict=True)
                        #if c in ['Good']:
                            #self.deBug('%s %s %s' % (c,t,e))
                            #self.deBug(d[c][t]['pos'][e])
                            #self.deBug(d[c][t]['freq'][e])
            self.printLog('\r#FREQ','Converting nucleotide counts to frequencies complete.')
            self.opt['DeBug'] = db
            ## ~ [1b] ~ Calculate Information contents ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.progLog('#INFO','Calculating Positional Information Content...')
            for c in self._coreType():
                for t in ['tic','con']:
                    d[c][t]['info'] = {}
                    for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                        d[c][t]['info'][e] = {}
                        for i in range(self.stat['KozWin']):
                            d[c][t]['info'][e][i-k] = rje.entropyDict(d[c][t]['pos'][e][i-k],ikeys=['G','A','T','C'])
                            d[c][t]['info'][e][1+i] = rje.entropyDict(d[c][t]['pos'][e][i+1],ikeys=['G','A','T','C'])
                        #self.deBug('%s %s %s' % (c,t,e))
                        #self.deBug(d[c][t]['pos'][e])
                        #self.deBug(d[c][t]['freq'][e])
            self.printLog('\r#INFO','Calculating Positional Information Content complete.')
            ### ~ [2] ~ Rank cores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.progLog('#RANK','Ranking core sequences ...')
            for coresplit in self._coreType():
                for t in ['tic','con']:
                    for s in ['kozcore','wegcore']:
                        r = string.replace(s,'core','rank')
                        for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                            d[coresplit][t][r][e] = rje.rankDict(d[coresplit][t][s][e],rev=True,absolute=True,lowest=True)
            #self.deBug(d['All']['tic']['kozrank']['AL'])
            self.printLog('\r#RANK','Ranking core sequences complete.')
            ### ~ [3] ~ Score each sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            '''
            - Ranking for "core" sequences within each eType
            - Frequency "f" scores for each window using sums of observed frequencies at each position
            '''
            sx = 0.0; stot = seqlist.seqNum(); 
            for seq in seqlist.seqs():
                self.progLog('\r#SCORE','Scoring sequences: %.2f%%' % (sx/stot)); sx += 100.0
                if len(seq.info['TICWin']) < (2 * self.stat['KozWin']): continue                    
                seq.dict['Score'] = {}
                for coresplit in self._coreType():
                    sdata = seq.dict['Score'][coresplit] = {}
                    for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                        try:
                            # - tfX = TICwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
                            score = 0
                            for i in range(self.stat['KozWin']):
                                try: score += d[coresplit]['tic']['freq'][e][i-k][seq.info['TICWin'][i]]
                                except: pass    #self.deBug(d[coresplit]['tic']['freq'][e][i-k])
                                try: score += d[coresplit]['tic']['freq'][e][1+i][seq.info['TICWin'][i+self.stat['KozWin']]]
                                except: pass    #self.deBug(d[coresplit]['tic']['freq'][e][1+i])
                            sdata['tf%s' % e] = score
                            # - cfX = Conwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
                            score = 0
                            for i in range(self.stat['KozWin']):
                                try:
                                    if seq.info['TICWin'][i] not in 'GATC': score += 0.25
                                    else: score += d[coresplit]['con']['freq'][e][i-k][seq.info['TICWin'][i]]
                                    if seq.info['TICWin'][i+self.stat['KozWin']] not in 'GATC': score += 0.25
                                    else: score += d[coresplit]['con']['freq'][e][1+i][seq.info['TICWin'][i+self.stat['KozWin']]]
                                except:
                                    self.errorLog('!')
                                    #self.deBug(d[coresplit]['con']['freq'][e])
                                    #self.deBug(d[coresplit]['con']['pos'][e])
                            sdata['cf%s' % e] = score
                            kozcore = rje.strSub(seq.info['Core'],-4,-2,'ATG')
                            if 'N' not in kozcore:
                                # - trX = TICwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
                                sdata['tr%s' % e] = d[coresplit]['tic']['kozrank'][e][kozcore]
                                # - crX = Conwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
                                sdata['cr%s' % e] = d[coresplit]['con']['kozrank'][e][kozcore]
                            wegcore = rje.strSub(seq.info['WegCore'],-4,-2,'ATG')
                            if 'N' not in wegcore:
                                # - twX = TICwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
                                sdata['tw%s' % e] = d[coresplit]['tic']['wegrank'][e][wegcore]
                                # - cwX = Conwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
                                sdata['cw%s' % e] = d[coresplit]['con']['wegrank'][e][wegcore]
                        except:
                            self.errorLog('%s %s %s scoring error' % (seq.info['AccNum'],coresplit,e))
                            #self.deBug(seq.info)
            self.printLog('\r#SCORE','Scoring of sequences complete.')
            ### ~ [4] ~ Perform Core-specifc calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            '''
            The following data is recorded for each Core and WegCore: (*.core.tdt, *.wegcore.tdt)
            - tfX = Mean frequency scores (excluding Core) for TICwin with this core
            - tiX = Information content (excluding Core) for TICwin with this core
            - cfX = Mean frequency scores (excluding Core) for Conwin with this core
            - ciX = Information content (excluding Core) for Conwin with this core
            '''
            ## ~ [4a] ~ Setup GC content dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0.0; stot = seqlist.seqNum(); gatc = {}; seqn = {}
            for c in self._coreType():
                gatc[c] = {}; seqn[c] = {}
                for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                    gatc[c][e] = {}; seqn[c][e] = 0
                    for n in 'GATC': gatc[c][e][n] = 0.0
            for seq in seqlist.seqs():
                self.progLog('\r#GATC','Calculating full SeqList GATC frequency: %.2f%%' % (sx/stot)); sx += 100.0
                e = seq.info['eType']
                etypes = [e,'X%s' % e[-1],'%sX' % e[0],'XX']
                for c in ['All',seq.info['CoreType']]:
                    for e in etypes:
                        if e == 'UX': continue
                        seqn[c][e] += 1
                        for n in 'GATC':
                            try: gatc[c][e][n] += seq.dict['GC'][n]
                            except: pass
            for c in self._coreType():
                for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'): rje.dictFreq(gatc[c][e])
            self.printLog('\r#GATC','Calculating full SeqList GATC frequency complete.')

            ## ~ [4b] ~ Setup CoreData dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cx = 0.0; ctot = 5 * 256
            cdata = self.dict['CoreData'] = {}
            for c in self._coreType():
                cdata[c] = {}
                for x1 in 'GATC':
                    for x2 in 'GATC':
                        for x3 in 'GATC':
                            for x4 in 'GATC':
                                self.progLog('\r#CORE','Calculating CoreData: %.2f%%' % (cx/ctot)); cx += 100.0
                                kozcore = '%s%s%sATG%s' % (x1,x2,x3,x4)
                                wegcore = '%s%snn%snnATG%s' % (x1,x2,x3,x4)
                                cdata[c][kozcore] = {}
                                cdata[c][wegcore] = {}
                                for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                                    ex = seqn[c][e]
                                    for n in (x1,x2,x3,x4): ex *= gatc[c][e][n]
                                    cdata[c][kozcore][e] = {'te':ex,'ce':ex,
                                                            'tn':d[c]['tic']['kozcore'][e][kozcore],
                                                            'tr':d[c]['tic']['kozrank'][e][kozcore],
                                                            'cn':d[c]['con']['kozcore'][e][kozcore],
                                                            'cr':d[c]['con']['kozrank'][e][kozcore]}
                                    cdata[c][kozcore][e]['tp'] = rje.logPoisson(d[c]['tic']['kozcore'][e][kozcore],ex,callobj=self)
                                    cdata[c][kozcore][e]['ts'] = rje.binomial(1,256,cdata[c][kozcore][e]['tp'],callobj=self)
                                    cdata[c][kozcore][e]['cp'] = rje.logPoisson(d[c]['con']['kozcore'][e][kozcore],ex,callobj=self)
                                    cdata[c][kozcore][e]['cs'] = rje.binomial(1,256,cdata[c][kozcore][e]['cp'],callobj=self)
                                    cdata[c][wegcore][e] = {'te':ex,'ce':ex,
                                                            'tn':d[c]['tic']['wegcore'][e][wegcore],
                                                            'tr':d[c]['tic']['wegrank'][e][wegcore],
                                                            'cn':d[c]['con']['wegcore'][e][wegcore],
                                                            'cr':d[c]['con']['wegrank'][e][wegcore]}
                                    cdata[c][wegcore][e]['tp'] = rje.logPoisson(d[c]['tic']['wegcore'][e][wegcore],ex,callobj=self)
                                    cdata[c][wegcore][e]['ts'] = rje.binomial(1,256,cdata[c][wegcore][e]['tp'],callobj=self)
                                    cdata[c][wegcore][e]['cp'] = rje.logPoisson(d[c]['con']['wegcore'][e][wegcore],ex,callobj=self)
                                    cdata[c][wegcore][e]['cs'] = rje.binomial(1,256,cdata[c][wegcore][e]['cp'],callobj=self)
            self.printLog('\r#CORE','Calculating CoreData complete.')
            #for bugcore in ['AAAATGA','AAAATGC','AAAATGG','AAAATGT']:
            #    self.deBug(d['All']['tic']['kozcore']['AL'][bugcore])
            #    self.deBug(cdata['All'][bugcore])
            ## ~ [4c] ~ Calculate scores for cores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #sx = 0.0; stot = seqlist.seqNum()
            #for seq in seqlist.seqs():
            #    self.progLog('\r#SCORE','Scoring cores: %.2f%%' % (sx/stot)); sx += 100.0
            #    if len(seq.info['TICWin']) < (2 * self.stat['KozWin']): continue
            #    kozcore = rje.strSub(seq.info['Core'],-4,-2,'ATG')
            #    wegcore = rje.strSub(seq.info['WegCore'],-4,-2,'ATG')
            #    for c in ['All',seq.info['CoreType']]:
            #        e = seq.info['eType']
            #        etypes = [e,'X%s' % e[-1],'%sX' % e[0],'XX']
            #        for e in etypes:
            #            if e == 'UX': continue
            #                cdata[c][kozcore][e]['tf'].append(0.0)
            #                cdata[c][kozcore][e]['ti'].append(0.0)
            #                cdata[c][kozcore][e]['cf'].append(0.0)
            #                cdata[c][kozcore][e]['ci'].append(0.0)
            #self.printLog('\r#SCORE','Scoring of cores complete.')
                      
            ### ~ [5] ~ End and return sequence lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            #x#self.pickleMe('%s.score' % self.info['Basefile'])
            return seqlist
        except: self.errorLog('Feb-2010-Kozak Scoring Error')
#########################################################################################################################
    def _kozak10_outputENST(self):  ### Output transcript data
        '''Output transcript data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = self.dict['Headers']['enst']
            seqlist = self.obj['SeqList']
            ### ~ [2] ~ Output data from each Transcript ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seqs():
                self.progLog('\r#OUT','Transcript output: %.2f%%' % (sx/stot)); sx += 100.0
                for coresplit in self._coreType():
                    sdata = rje.combineDict({'CoreSplit':coresplit},seq.dict['Score'][coresplit])
                    sdata = rje.combineDict(sdata,seq.info)
                    sdata['Transcript ID'] = sdata['AccNum']
                    sdata['Gene ID'] = sdata['EnsG']
                    sdata = rje.combineDict(sdata,seq.stat,overwrite=False)
                    for n in 'GATC':
                        try: sdata['reg%s' %n] = seq.dict['GC'][n]
                        except: sdata['reg%s' %n] = 0.0
                    rje.delimitedFileOutput(self,'%s.enst.tdt' % (self.info['Basefile']),headers,datadict=sdata)
            self.printLog('\r#OUT','Transcript output complete.')
            return
        except: self.errorLog('Feb-2010-Kozak Transcript Output Error')
#########################################################################################################################
    def _kozak10_outputALCore(self,ctype='alcore'):  ### Output core data
        '''Output core data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = self.dict['Headers'][ctype]
            d = self.dict['CoreAL']
            dx = 0.0; dtot = 5
            ### ~ [2] ~ Output data from each Core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            for c in self._coreType():
                self.progLog('\r#OUT','ALCore data output: %.2f%%' % (dx/dtot)); dx += 100.0
                for core in rje.sortKeys(d[c]):
                    rje.delimitedFileOutput(self,'%s.%s.tdt' % (self.info['Basefile'],ctype),headers,datadict=d[c][core])    
            self.printLog('\r#OUT','ALCore data output complete.')
        except: self.errorLog('Feb-2010-Kozak Core Output Error')
#########################################################################################################################
    def _kozak10_outputCore(self,ctype='core'):  ### Output core data
        '''Output core data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = self.dict['Headers'][ctype]
            d = self.dict['CoreData']
            dx = 0.0; dtot = 5 * len(d['All'])
            ### ~ [2] ~ Output data from each Core ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            for c in self._coreType():
                for core in rje.sortKeys(d[c]):
                    for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                        for x in d[c][core][e]: d[c][core]['%s%s' % (x,e)] = d[c][core][e][x]
                    self.progLog('\r#OUT','Core (%s) data output: %.2f%%' % (ctype,(dx/dtot))); dx += 100.0
                    if ctype == 'core' and len(core) != 7: continue
                    elif ctype == 'wegcore' and len(core) == 7: continue
                    cdata = rje.combineDict({'CoreSplit':c,'Core':core},d[c][core])
                    rje.delimitedFileOutput(self,'%s.%s.tdt' % (self.info['Basefile'],ctype),headers,datadict=cdata)    
            self.printLog('\r#OUT','Core (%s) data output complete.' % (ctype))
        except: self.errorLog('Feb-2010-Kozak Core Output Error')
#########################################################################################################################
    def _kozak10_outputPos(self,ctype='core'):  ### Output position data
        '''Output position data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = self.dict['Headers'][ctype]
            d = self.dict['CoreData']
            ### ~ [2] ~ Output data from each Transcript ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!#
            return
            '''
            The following data is recorded for each position in the window: (*.pos.tdt)
            - eType = AL/AS/US/XL/XS/AX/XX.
            - CoreSplit = All/Good/augG/Rxx/Poor.
            - G = Count of G
            - A = Count of G
            - T = Count of G
            - C = Count of G
            - Info = Information Content of position.
            '''
            self.dict['Headers']['pos'] = ['CoreSplit','eType','Pos','G','A','T','C','Info']
            for out in self.dict['Headers']:
                rje.delimitedFileOutput(self,'%s.%s.tdt' % (self.info['Basefile'],out),self.dict['Headers'][out],rje_backup=True)

        except: self.errorLog('Feb-2010-Kozak Headers Error')
#########################################################################################################################
    def _coreType(self,core=None):   ### 
        '''Returns core type.'''
        if not core: return ['All','Good','augG','Mid','Weak','Poor','Alt']
        if core in self.list['GoodIC']: return 'Good'
        if len(core) != 7: return 'Error'
        if core[0] in ['G','A'] and core[-4:-1] in self.list['NonAUG'] and core[-1] == 'G': return 'Alt'
        elif core[-4:-1] != 'ATG': return 'Poor'    #'Bad'
        if core[-1] == 'G' and core[0] in ['G','A']: return 'Good'
        elif core[-1] == 'G': return 'augG'
        elif core[0] == 'C' and core[-1] == 'C': return 'Weak'
        elif core[0] == 'C' and core[-1] == 'A': return 'Weak'
        elif core[0] == 'U' and core[-1] == 'A': return 'Weak'
        elif core[0] == 'U' and core[-1] == 'U': return 'Weak'
        return 'Mid'
#########################################################################################################################
    def _kozak10_outputTypes(self):     ### Output position data
        '''Output position data.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            headers = self.dict['Headers']['types']
            tfile = '%s.types.tdt' % self.info['Basefile']
            tdata = {}
            for coresplit in self._coreType():   #['All','Good','augG','Rxx','Poor']:
                tdata[coresplit] = {'CoreType':coresplit}
                for h3 in string.split('AL/AS/US/NL/NS/NX/XL/XS/AX/XX','/'):
                    tdata[coresplit][h3] = 0
            sx = 0.0; stot = self.obj['SeqList'].seqNum()
            for seq in self.obj['SeqList'].seqs():
                self.progLog('\r#TYPE','Generating Type data: %.2f%%' % (sx/stot)); sx += 100.0
                seq.info['CoreType'] = self._coreType(seq.info['Core'])
                ctypes = ['All',seq.info['CoreType']]
                e = seq.info['eType']
                etypes = [e,'X%s' % e[-1],'%sX' % e[0],'XX']
                for c in ctypes:
                    for e in etypes:
                        if e in ['UX']: continue
                        tdata[c][e] += 1    
            self.printLog('\r#TYPE','Generating Type data complete.')
            ### ~ [2] ~ Output data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'types' in self.list['Output']:
                for coresplit in self._coreType():
                    rje.delimitedFileOutput(self,tfile,headers,datadict=tdata[coresplit])

        except: self.errorLog('Feb-2010-Kozak Headers Error')
#########################################################################################################################
    def _kozak10_headers(self):     ### Setup outfiles
        '''Setup outfiles.'''
        try:## ~ [1b] ~ Setup outfiles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['Headers'] = {}
            '''
            The following data is recorded for each transcript: (*.enst.tdt)
            - CoreSplit = All/Good/augG/Rxx/Poor.
            - Transcript ID
            - Gene ID
            - TICwin = Window around TIC (as set by kozwin)
            - Core = Core Kozak window -3 to + 4 (XXXAUGX)
            - CoreType = Good/augG/Rxx/Poor
            - Conwin = Control window around 5' AUG (as set by kozwin)
            - WegCore = Wegrzyn Core sequence XXnnXnnAUGX
            - Redundancy = No. filtered identical TIC windows for this gene
            - AIC = No. of *different* TIC windows for this gene.
            - e5utr = UTR length from input seq <3>
            - c5utr = UTR length from input seq <1>
            - eType = AL/AS/US (see above)
            - regG = G count from first 500nt of seq <2>
            - regC = G count from first 500nt of seq <2>
            - regA = A count from first 500nt of seq <2>
            - regT = T count from first 500nt of seq <2>
            - tfX = TICwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
            - cfX = Conwin Frequency scores, where X is AL/AS/US/XL/XS/AX/XX.
            - trX = TICwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
            - crX = Conwin core ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
            - twX = TICwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
            - cwX = Conwin WegCore ranking scores, where X is AL/AS/US/XL/XS/AX/XX.
            '''
            self.dict['Headers']['enst'] = ['CoreSplit','Transcript ID','Gene ID','TICwin',
                                            'Core','CoreType','Conwin','WegCore','Redundancy','AIC',
                                            'e5utr','c5utr','eType','regG','regC','regA','regT']
            for h1 in 'tc':
                for h2 in 'frw':
                    for h3 in string.split('AL/AS/US/XL/XS/AX/XX','/'):
                        self.dict['Headers']['enst'].append('%s%s%s' % (h1,h2,h3))
            '''
            The following data is recorded for each Core and WegCore: (*.core.tdt, *.wegcore.tdt)
            - CoreSplit = All/Good/augG/Rxx/Poor.
            - Core = Core sequence
            - tnX = No. of TICwin with this core
            - trX = Rank for TICwin with this score
            - cnX = No. of Conwin with this core
            - crX = Rank for Conwin with this core
            '''
            self.dict['Headers']['core'] = ['CoreSplit','Core']
            for h1 in 'tc':
                for h3 in string.split('AL/AS/US/XL/XS/AX/XX','/'):
                    for h2 in 'nreps':  # 'fir':    #!# Instead of this crap, output probability of observation
                        self.dict['Headers']['core'].append('%s%s%s' % (h1,h2,h3))
            self.dict['Headers']['wegcore'] = self.dict['Headers']['core'][0:]
            self.dict['Headers']['alcore'] = ['CoreType','CoreSplit','Core']
            for h1 in 'tc':
                #for h3 in string.split('AL/AS/US/XL/XS/AX/XX','/'):
                for h2 in 'nrepsof':  
                    self.dict['Headers']['alcore'].append('%s%s' % (h1,h2))
            '''
            The following data is recorded for each position in the window: (*.pos.tdt)
            - eType = AL/AS/US/XL/XS/AX/XX.
            - CoreSplit = All/Good/augG/Rxx/Poor.
            - G = Count of G
            - A = Count of G
            - T = Count of G
            - C = Count of G
            - Info = Information Content of position.
            '''
            self.dict['Headers']['pos'] = ['CoreSplit','eType','Pos','G','A','T','C','Info']
            
            self.dict['Headers']['types'] = ['CoreType'] + string.split('AL/AS/AX/US/NL/NS/NX/XL/XS/XX','/')           

            for out in self.dict['Headers']:
                if out not in self.list['Output']: continue
                rje.delimitedFileOutput(self,'%s.%s.tdt' % (self.info['Basefile'],out),self.dict['Headers'][out],rje_backup=True)

        except: self.errorLog('Feb-2010-Kozak Headers Error')
#########################################################################################################################
    def _kozak10_counts(self):     ### Perform Calculations for both TICwin and Conwin sequences 
        '''Perform Calculations for both TICwin and Conwin sequences.'''
        try:### ~ [1] ~ Setup dictionaries to store count data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = self.obj['SeqList']
            # These calculations are performed for all sequences and then repeated for each CoreType independently.
            # This is output in the "CoreSplit" field of each output (All/Good/augG/Rxx/Poor).
            self.dict['Data'] = d = {}  # Main data dictionary
            cx = {}; ex = {}            # Counts of different types
            for coresplit in self._coreType():
                cx[coresplit] = 0
                d[coresplit] = {}   # Data dictionary for each split
                for t in ['tic','con']:
                    d[coresplit][t] = {}
                    for s in ['pos','kozcore','wegcore','kozrank','wegrank']:
                        d[coresplit][t][s] = {}
                        for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
                            ex[e] = 0
                            d[coresplit][t][s][e] = {}
                            if s == 'pos':
                                for i in range(self.stat['KozWin']):
                                    for p in [i-self.stat['KozWin'],i+1]:
                                        d[coresplit][t][s][e][p] = {}
                                        for n in 'GATC': d[coresplit][t][s][e][p][n] = 0
                            else:
                                for x1 in 'GATC':
                                    for x2 in 'GATC':
                                        for x3 in 'GATC':
                                            for x4 in 'GATC':
                                                if s[:3] == 'koz': x = '%s%s%sATG%s' % (x1,x2,x3,x4)
                                                else: x = '%s%snn%snnATG%s' % (x1,x2,x3,x4)
                                                d[coresplit][t][s][e][x] = 0
            ### ~ [2] ~ Perform Calculations for both TICwin and Conwin sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = seqlist.seqNum(); k = self.stat['KozWin']
            db = self.opt['DeBug']; badcore = []
            for seq in seqlist.seqs():
                self.progLog('\r#POS','Generating position-specific & core frequency data: %.2f%%' % (sx/stot)); sx += 100.0
                if len(seq.info['TICWin']) < (2 * self.stat['KozWin']): continue                    
                ctypes = ['All',seq.info['CoreType']]
                for c in ctypes: cx[c] += 1
                e = seq.info['eType']
                ex[e] += 1
                etypes = [e,'X%s' % e[-1],'%sX' % e[0],'XX']
                kozcore = seq.info['Core']
                wegcore = seq.info['WegCore']
                if seq.info['ConWin']:
                    concore = seq.info['ConCore']
                    conweg = seq.info['ConWeg']
                #x#self.deBug(seq.info)
                ## ~ [2a] ~ Core counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for c in ctypes:
                    for e in etypes:
                        if e == 'UX': continue
                        # Frequency for "core" sequences within each eType (GATC only)
                        if kozcore in d[c]['tic']['kozcore'][e]: d[c]['tic']['kozcore'][e][kozcore] += 1 #= 0
                        if wegcore in d[c]['tic']['wegcore'][e]: d[c]['tic']['wegcore'][e][wegcore] += 1 #= 0
                        if seq.info['ConWin']:
                            if concore in d[c]['con']['kozcore'][e]: d[c]['con']['kozcore'][e][concore] += 1 #= 0
                            else:
                                if concore not in badcore:
                                    self.printLog('\r#BAD','Cannot use control core "%s" for core counts' % concore)
                                    badcore.append(concore)
                                    #x#self.deBug(d[c]['con']['kozcore'][e])
                            if conweg in d[c]['con']['wegcore'][e]: d[c]['con']['wegcore'][e][conweg] += 1 #= 0
                            else:
                                if conweg not in badcore:
                                    self.printLog('\r#BAD','Cannot use control core "%s" for core counts' % conweg)
                                    badcore.append(conweg)
                                    #x#self.deBug(d[c]['con']['wegcore'][e])
                ## ~ [2b] ~ Nucleotide counts for each position of the window, within each eType ~~ ##
                for i in range(self.stat['KozWin']):
                    t5 = seq.info['TICWin'][i]
                    t3 = seq.info['TICWin'][i+self.stat['KozWin']]
                    if seq.info['ConWin']:
                        c5 = seq.info['ConWin'][i]
                        c3 = seq.info['ConWin'][i+self.stat['KozWin']]
                    for c in ctypes:
                        for e in etypes:
                            if e == 'UX': continue
                            # TIC #
                            if t5 not in d[c]['tic']['pos'][e][i-k]: d[c]['tic']['pos'][e][i-k][t5] = 0
                            d[c]['tic']['pos'][e][i-k][t5] += 1
                            if t3 not in d[c]['tic']['pos'][e][1+i]: d[c]['tic']['pos'][e][1+i][t3] = 0
                            d[c]['tic']['pos'][e][1+i][t3] += 1
                            # Control #
                            if seq.info['ConWin']:
                                if c5 not in d[c]['con']['pos'][e][i-k]: d[c]['con']['pos'][e][i-k][c5] = 0
                                d[c]['con']['pos'][e][i-k][c5] += 1
                                if c3 not in d[c]['con']['pos'][e][1+i]: d[c]['con']['pos'][e][1+i][c3] = 0
                                d[c]['con']['pos'][e][1+i][c3] += 1
            self.printLog('\r#POS','Generation of position-specific & core frequency data complete.')
            for c in self._coreType():
                self.printLog('#%s' % c,'"%s" sequence count: %s' % (c,rje.integerString(cx[c])))
            for e in rje.sortKeys(ex):
                if e[0] != 'X' and e[1] != 'X': self.printLog('#%s' % e,'"%s" sequence count: %s' % (e,rje.integerString(ex[e])))
            #x#self.pickleMe('%s.counts' % self.info['Basefile'])
            self.opt['DeBug'] = db
            #for c in self._coreType():
            #    for e in string.split('AL/AS/US/NL/NS/XL/XS/AX/NX/XX','/'):
            #       for s in ['kozcore','wegcore','kozrank','wegrank']:
            #           #print c, e, s
            #           self.deBug(d[c]['tic'][s][e])
            #           self.deBug(d[c]['con'][s][e])
            #self.opt['DeBug'] = db
            return d
        except: self.errorLog('Feb-2010-Kozak Counts Error')
#########################################################################################################################
    def coreXXX(self,core,xxx='XXX'): return rje.strSub(core,-4,-2,xxx)
#########################################################################################################################
    def _kozak10_cleanSeq(self,seqlist,type='input',goodacc=[],minlen=0):     ### Cleans up sequence and reduces by goodacc
        '''Cleans up sequence and reduces by goodacc.'''
        sx = 0; stot = seqlist.seqNum()
        for seq in seqlist.seqs()[0:]:
            try:
                self.progLog('\r#INFO','Extracting %s sequence info: %.2f%%' % (type,sx/stot)); sx += 100.0
                details = string.split(seq.info['FullName'],'|')
                seq.info['AccNum'] = details[0]
                if goodacc and seq.info['AccNum'] not in goodacc: seqlist.seq.remove(seq); continue
                seq.info['Description'] = string.join(details[1:],'|')
                try: seq.info['Gene'] = details[2]
                except: seq.info['Gene'] = ''
                seq.info['DBase'] = 'EnsEMBL'
                if seq.info['Sequence'][:8] == 'SEQUENCE': seq.info['Sequence'] = ''
                if minlen > 0 and seq.seqLen() < minlen: seqlist.seq.remove(seq); continue
            except:
                self.errorLog('Problem extracting info for "%s"' % seq.info['FullName'])
                seq.info['AccNum'] = '!ERR!'
                seqlist.seq.remove(seq); continue
        self.printLog('\r#INFO','Extracted %s sequence info: %s -> %s retained.' % (type,rje.integerString(stot),rje.integerString(seqlist.seqNum())))
#########################################################################################################################
    def _kozak10_cleanLoadSeq(self,type='input',goodacc=[],minlen=0):     ### Cleans up sequence and reduces by goodacc
        '''Cleans up sequence and reduces by goodacc.'''
        sx = 0; sj = 0; lastline = 'START'; seqdict = {}
        LOAD = open('%s.%s.fasta' % (self.info['Kozak'],type),'r')
        while lastline:
            (rawseq,lastline) = self.obj['SeqList'].nextFasSeq(LOAD,lastline,raw=True)
            if not rawseq: break
            accnum = string.split(rawseq[0],'|')[0]
            sequence = rawseq[1].upper()
            sx += 1; sj += 1
            if sj == 100: self.progLog('\r#LOAD','Loading %s sequences: %s' % (type,rje.integerString(sx))); sj = 0
            if goodacc and accnum not in goodacc: continue
            if sequence[:8] == 'SEQUENCE': continue
            if minlen > 0 and len(sequence) < minlen: continue
            seqdict[accnum] = sequence
        LOAD.close()
        self.printLog('\r#LOAD','Loaded %s %s sequences: %s retained.' % (type,rje.integerString(sx),rje.integerString(len(seqdict))))
        return seqdict
#########################################################################################################################
    def _kozak10_loadSeq(self):     ### Loads sequences and peforms initial processing before pickling
        '''Loads sequences and peforms initial processing before pickling.'''
        try:### ~ [1] ~ Load Sequences for processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Data is read in from four sources:
            # 1. Full length EnsEMBL cDNAs
            # 2. 5' flank of TSS. This should not include 5' UTR where annotation is good.
            # 3. 5'UTR sequences.
            # 4. CDS.
            seqcmd = self.cmd_list + ['autoload','accnr=F','seqnr=F','gnspacc=F','seqin=None']
            self.obj['SeqList'] = seqlist = rje_seq.SeqList(self.log,seqcmd)  # This object will store all the final sequence data
            ## ~ [1a] ~ CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqcmd[-1] = 'seqin=%s.cds.fasta' % self.info['Kozak']
            seqcds = rje_seq.SeqList(self.log,seqcmd)
            self._kozak10_cleanSeq(seqcds,type='CDS',goodacc=[])
            seqcds = seqcds.seqNameDic('AccNum')
            if '!ERR!' in seqcds: seqcds.pop('!ERR!')
            goodacc = rje.sortKeys(seqcds)
            open('%s.cds.enst' % self.info['Kozak'],'w').write(string.join(goodacc,'\n'))
            self.printLog('#ACC','%s CDS enst output for future reduction' % rje.integerString(len(goodacc)))
            ## ~ [1b] ~ 5' Flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #seqcmd[-1] = 'seqin=%s.5flank.fasta' % self.info['Kozak']
            #seq5flank = rje_seq.SeqList(self.log,seqcmd)
            #self._kozak10_cleanSeq(seq5flank,'5\' flank',goodacc,minlen=self.stat['Flank'])
            #seq5flank = seq5flank.seqNameDic('AccNum')
            #if '!ERR!' in seq5flank: seq5flank.pop('!ERR!')
            seq5flank = self._kozak10_cleanLoadSeq('5flank',goodacc,minlen=self.stat['Flank'])
            ## ~ [1c] ~ 5' UTR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #seqcmd[-1] = 'seqin=%s.5utr.fasta' % self.info['Kozak']
            #seq5utr = rje_seq.SeqList(self.log,seqcmd)
            #self._kozak10_cleanSeq(seq5utr,'5\' UTR',goodacc)
            #seq5utr = seq5utr.seqNameDic('AccNum')
            #if '!ERR!' in seq5utr: seq5utr.pop('!ERR!')
            seq5utr = self._kozak10_cleanLoadSeq('5utr',goodacc)
            ## ~ [1d] ~ cDNA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            g2t = {}    # Dictionary of {gene:[transcripts]}
            t2g = {}    # Dictionary of {trans:gene}
            seqcmd[-1] = 'seqin=%s.cdna.fasta' % self.info['Kozak']
            #x#seqcmd += ['goodacc=%s.cds.enst' % self.info['Kozak']]
            seqcdna = rje_seq.SeqList(self.log,seqcmd)
            sx = 0; stot = seqcdna.seqNum()
            for seq in seqcdna.seqs()[0:]:
                try:
                    self.progLog('\r#INFO','Extracting cDNA sequence info: %.2f%%' % (sx/stot)); sx += 100.0
                    if goodacc and seq.info['AccNum'] not in goodacc: seqcdna.seq.remove(seq); continue
                    seq.info['Gene'] = rje.matchExp('gene:(\S+)',seq.info['Name'])[0]
                    seq.info['DBase'] = 'EnsEMBL'
                    t2g[seq.info['AccNum']] = seq.info['Gene']
                    if seq.info['Gene'] not in g2t: g2t[seq.info['Gene']] = []
                    g2t[seq.info['Gene']].append(seq.info['AccNum'])
                    if seq.info['Gene'] in self.list['Track'] and seq.info['AccNum'] not in self.list['Track']: self.list['Track'].append(seq.info['AccNum'])
                    if seq.info['AccNum'] in self.list['Track'] and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])
                except:
                    self.errorLog('Problem extracting info for "%s"' % seq.info['FullName'])
                    seq.info['AccNum'] = '!ERR!'
                    seqcdna.seq.remove(seq); continue
            self.printLog('\r#INFO','Extracted cDNA sequence info: %s -> %s retained.' % (rje.integerString(stot),rje.integerString(seqcdna.seqNum())))
            seqcdna = seqcdna.seqNameDic('AccNum')
            if '!ERR!' in seqcdna: seqcdna.pop('!ERR!')
            ## ~ [1e] ~ Summarise input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','#~~# %s Input #~~#' % self.info['Kozak'],timeout=False,log=False)
            self.printLog('#CDS','%s CDS sequences' % rje.integerString(len(seqcds)))
            self.printLog('#FLANK','%s 5\' flank sequences' % rje.integerString(len(seq5flank)))
            self.printLog('#UTR','%s 5\' UTR sequences' % rje.integerString(len(seq5utr)))
            self.printLog('#CDNA','%s cDNA sequences' % rje.integerString(len(seqcdna)))
            self.printLog('#ENST','%s ENST IDs' % rje.integerString(len(t2g)))
            self.printLog('#GENE','%s genes' % rje.integerString(len(g2t)))
            self.printLog('#~~#','#~~#',timeout=False,log=False)

            ### ~ [2] ~ Compile Sequences for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Sequences are compiled and cross-checked with each other. Where sequence <3> is over 200nt, it is assumed
            # that the TIC is at position -200 (i.e. 200 from end). Sequence <3> should also match the start of sequence
            # <1> for a given transcript. Where <3> is <200nt, the full length cDNA is used and the TIC is assumed to be
            # the first AUG. In terms of TIC context, this gives different types of sequence ('etype'):
            # - AL. 5'UTR annotated and longer than kozwin.
            # - AS. 5'UTR annotated but shorter than kozwin.
            # - US. 5'UTR inferred from cDNA but shorter than kozwin.
            # - These are also combined for additional combos (XL,XS,AX,XX) for calculations
            # - NL. 5'UTR annotated and longer than kozwin but non-AUG TIC. Not used for calculations (but still rated).
            # - NS. 5'UTR annotated but shorter than kozwin, plus non-AUG TIC. Not used for calculations (but still rated).
            genewin = {}; gx = 0.0; gtot = len(g2t); missx = 0; no_atg = []; atgx = 0
            self.dict['AIC'] = {}       # Dictionary of non-ATG codons and their AccNum
            seqlist.seq = []
            no_cds = 0; errx = 0
            for gene in rje.sortKeys(g2t):
                self.progLog('\r#GENE','Processing genes: %.2f%%' % (gx/gtot)); gx += 100.0
                for enst in g2t[gene]:
                    ## Get all details for sequence into one sequence object ##
                    seq = seqcdna.pop(enst)
                    seq.info['cDNA'] = seq.info['Sequence']
                    seq.info['EnsG'] = seq.info['Gene']
                    try:
                        cds = seqcds.pop(enst)
                        if not cds.info['Sequence']: raise ValueError
                        seq.info['Gene'] = cds.info['Gene']
                        seq.info['CDS'] = cds.info['Sequence']
                        if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']) and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])                    
                    except: no_cds += 1; continue    # Not a CDS
                    if self.info['NRGene'].lower() == 'ensg': nrgene = seq.info['EnsG']
                    else: nrgene = seq.info['Gene']
                    if nrgene not in genewin: genewin[nrgene] = {}      
                    try:
                        us = seq5flank.pop(enst)
                        #seq.info['Gene'] = us.info['Gene']
                        seq.info['Flank'] = us#.info['Sequence']
                    except:
                        self.printLog('\r#NOFLANK','No 5\' flank sequence found for %s' % enst,screen=False)
                        try: utr = seq5utr.pop(enst)
                        except: self.printLog('\r#NOUTR','No 5\'utr sequence found for %s' % enst,screen=False)
                        missx += 1; continue
                    try:
                        utr = seq5utr.pop(enst)
                        seq.info['5utr'] = utr  #.info['Sequence']
                        if seq.info['5utr']:
                            if seq.info['Sequence'][:len(seq.info['5utr'])] != seq.info['5utr']:
                                self.printLog('\r#ERR','%s 5\' UTR sequence does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['Sequence'][:len(seq.info['5utr'])])
                                raise ValueError
                            no3 = seq.info['5utr'] + seq.info['CDS']
                            if seq.info['Sequence'].find(no3) != 0:
                                self.printLog('\r#ERR','%s 5\' UTR + CDS does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['CDS'])
                                self.deBug(no3)
                                self.deBug(seq.info['Sequence'])
                                raise ValueError
                        else: self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    except ValueError: errx += 1; continue
                    except:
                        seq.info['5utr'] = ''
                        self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    ## Assess sequence and dictate eType ##
                    seq.info['utr5'] = seq.info['5utr']
                    if seq.info['5utr']:
                        etype = 'A'
                        seq.stat['e5utr'] = len(seq.info['utr5'])
                    else:
                        etype = 'U'
                        seq.stat['e5utr'] = 0
                    kozwin = seq.info['CDS'][:self.stat['KozWin']]
                    if kozwin[:3] != 'ATG':
                        self.printLog('\r#AIC','%s AIC found for %s' % (kozwin[:3],enst),screen=False)
                        etype = 'N'
                        if kozwin[:3] not in self.dict['AIC']: self.dict['AIC'][kozwin[:3]] = []
                        self.dict['AIC'][kozwin[:3]].append(enst)
                    if len(seq.info['utr5']) < self.stat['KozWin']: etype += 'S'
                    else: etype += 'L'
                    seq.info['eType'] = etype
                    # For each transcript, the following sequences are then extracted:
                    # - The "TICwin"  window around TIC (as set by kozwin). [*.kozwin.fas]
                    # - A control "Conwin" window around the first AUG encountered in seq <2>. [*.conwin.fas]
                    flank = seq.info['Flank'] + seq.info['utr5']
                    seq.info['TICWin'] = flank[-self.stat['KozWin']:] + kozwin
                    if len(seq.info['TICWin']) < (2 * self.stat['KozWin']): missx += 1; continue
                    if seq.info['TICWin'] not in genewin[nrgene]: genewin[nrgene][seq.info['TICWin']] = []
                    genewin[nrgene][seq.info['TICWin']].append(seq)
                    seq.info['Core'] = flank[-3:] + kozwin[:4]
                    seq.info['WegCore'] = flank[-7:-5] + 'nn' + flank[-3] + 'nn' + kozwin[:4] # XXnnXnnAUGX
                    ## Add ConWin ##
                    flank = seq.info['Flank'][-self.stat['Flank']:]
                    conatg = flank[self.stat['KozWin']:-self.stat['KozWin']].find('ATG')
                    if conatg >= 0:
                        seq.info['ConWin'] = flank[conatg:][:2*self.stat['KozWin']]
                        seq.info['ConCore'] = seq.info['ConWin'][self.stat['KozWin']-3:][:7]
                        seq.info['ConWeg'] = seq.info['ConWin'][self.stat['KozWin']-7:][:2] + 'nn' + seq.info['ConWin'][self.stat['KozWin']-3] + 'nn' + seq.info['ConWin'][self.stat['KozWin']:][:4] # XXnnXnnAUGX
                    else:
                        seq.info['ConWin'] = ''
                        self.printLog('\r#CON','No upstream control ATG window found for %s' % (enst),screen=False)
                    # Each transcript is then classified according to its core sequence:
                    # - Good = "best" Kozak gccaugg [set by GoodIC list]
                    # - augG = Not "Good" but has G at +4
                    # - Rxx = Not "Good" or augG but has purine at -3
                    # - Poor = All the rest
                    seq.info['CoreType'] = self._coreType(seq.info['Core'])
                    #if seq.info['Core'] in self.list['GoodIC']: seq.info['CoreType'] = 'Good'
                    #elif seq.info['Core'][-1] == 'G': seq.info['CoreType'] = 'augG'
                    #elif seq.info['Core'][0] in ['G','A'] and seq.info['Core'][-1] == 'A': seq.info['CoreType'] = 'Rxx'     #!# Modified #!#
                    #elif seq.info['Core'][0] in ['G','A'] and seq.info['Core'][-1] == 'G' and seq.info['Core'][-4:-1] in ['CTG','GTG','ACG']: seq.info['CoreType'] = 'Alt'     #!# Modified #!#
                    #else: seq.info['CoreType'] = 'Poor'
                    #if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']): self.deBug(seq.details())
            self.printLog('\r#GENE','Processing %s genes complete: %s Genes with TICwin' % (rje.integerString(gtot),rje.integerString(len(genewin))))
            self.printLog('#MISS','%s cDNA transcripts excluded due to missing/short 5\' flanking sequence' % rje.integerString(missx))
            self.printLog('#NOCDS','%s cDNA transcripts excluded due to no CDS' % rje.integerString(no_cds))
            self.printLog('#ERRX','%s cDNA transcripts excluded due to sequence errors' % rje.integerString(errx))
            #x#self.printLog('#ATG','%s cDNA transcripts excluded due to missing 5\' UTR and no ATG in cDNA' % rje.integerString(len(no_atg)))
            open('%s.no_atg.txt' % self.info['Basefile'],'w').write(string.join(no_atg,'\n'))
            aicx = 0
            for aic in rje.sortKeys(self.dict['AIC']):
                self.printLog('#AIC','%s start -> %s transcripts.' % (aic,rje.integerString(len(self.dict['AIC'][aic]))),screen=False)
                aicx += len(self.dict['AIC'][aic])
            self.printLog('#AIC','Total AIC (non-ATG) transcripts = %s.' % (rje.integerString(aicx)))

            ### ~ [3] ~ Remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # For each NR gene, transcripts are checked for redundancy and only one representative of each
            # different TICwin sequence is kept.
            gx = 0.0; gtot = len(genewin); notx = 0
            for gene in rje.sortKeys(genewin):
                self.progLog('\r#GENE','Removing redundancy genes: %.2f%%' % (gx/gtot)); gx += 100.0
                if not genewin[gene]: notx += 1; continue
                for ticwin in genewin[gene]:
                    for seq in genewin[gene][ticwin]:
                        seq.stat['Redundancy'] = len(genewin[gene][ticwin]) - 1
                        seq.stat['AIC'] = len(genewin[gene])
                    for type in ['AL','AS','US','NL','NS']:    # Keep the "best"
                        if len(genewin[gene][ticwin]) == 1: break
                        for seq in genewin[gene][ticwin]:
                            if seq.info['eType'] == type: genewin[gene][ticwin] = [seq]; break
                    seqlist.seq.append(genewin[gene][ticwin][0])     # Add to sequence object
                    if seqlist.seq[-1].info['Core'][-4:-1] == 'ATG': atgx += 1
                    if gene in self.list['Track']: self.bugPrint('\rNR: %s' % seqlist.seq[-1].info)
                if gene in self.list['Track']: self.deBug('NR for gene %s complete' % gene)
            self.printLog('\r#GENE','Processing %s genes complete: %s NR sequences; %s without sequences.' % (rje.integerString(gtot),rje.integerString(seqlist.seqNum()),rje.integerString(notx)))
            self.printLog('\r#ATG','%s of %s NR sequences have ATG start codon.' % (rje.integerString(atgx),rje.integerString(seqlist.seqNum())))

            ### ~ [4] ~ Regional GC calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Regional GC content is estimated from the end of flanking seq <2>.
            sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seq:
                self.progLog('\r#GC','Calculating GC content from 5\' flanking sequence: %2.f%%' % (sx/stot)); sx += 100.0
                seq.dict['GC'] = {}
                for n in seq.info['Flank'][-self.stat['Flank']:]:
                    if n not in seq.dict['GC']: seq.dict['GC'][n] = 0
                    seq.dict['GC'][n] += 1
            self.printLog('\r#GC','Calculation of GC content from 5\' flanking sequence complete.')

            ### ~ [5] ~ Output unprocessed sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FLANK','%s 5\' flanking sequence entries remain unprocessed.' % rje.integerString(len(seq5flank)))
            if seq5flank: open('%s.seq5flank.left' % self.info['Basefile'],'w').write(string.join(rje.sortKeys(seq5flank),'\n'))
            self.printLog('#UTR','%s 5\' UTR sequence entries remain unprocessed.' % rje.integerString(len(seq5utr)))
            if seq5utr: open('%s.seq5utr.left' % self.info['Basefile'],'w').write(string.join(rje.sortKeys(seq5utr),'\n'))

            ### ~ [6] ~ Save SeqList object and pickle for later time-saving ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = seqlist
            self.pickleMe('%s.loaded' % self.info['Basefile'])
            return seqlist

        except: self.errorLog('Feb-2010-Kozak Loading Sequence Error')
#########################################################################################################################
    def _biol3050_loadSeq(self):    ### Loads sequences and peforms initial processing before pickling
        '''Loads sequences and peforms initial processing before pickling.'''
        try:### ~ [1] ~ Load Sequences for processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            pickle = self.kozak10_pickle()
            if pickle and pickle.obj['SeqList']:
                self.deBug(pickle.obj['SeqList'].seqNum())
                return pickle.obj['SeqList']
            # Data is read in from four sources:
            # 1. Full length EnsEMBL cDNAs
            # 2. 5' flank of TSS. This should not include 5' UTR where annotation is good.
            # 3. 5'UTR sequences.
            # 4. CDS.
            seqcmd = self.cmd_list + ['autoload','accnr=F','seqnr=F','gnspacc=F','seqin=None']
            self.obj['SeqList'] = seqlist = rje_seq.SeqList(self.log,seqcmd)  # This object will store all the final sequence data
            ## ~ [1a] ~ CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqcmd[-1] = 'seqin=%s.cds.fasta' % self.info['Kozak']
            seqcds = rje_seq.SeqList(self.log,seqcmd)
            self._kozak10_cleanSeq(seqcds,type='CDS',goodacc=[])
            seqcds = seqcds.seqNameDic('AccNum')
            if '!ERR!' in seqcds: seqcds.pop('!ERR!')
            goodacc = rje.sortKeys(seqcds)
            open('%s.cds.enst' % self.info['Kozak'],'w').write(string.join(goodacc,'\n'))
            self.printLog('#ACC','%s CDS enst output for future reduction' % rje.integerString(len(goodacc)))
            ## ~ [1c] ~ 5' UTR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seq5utr = self._kozak10_cleanLoadSeq('5utr',goodacc)
            ## ~ [1d] ~ cDNA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            g2t = {}    # Dictionary of {gene:[transcripts]}
            t2g = {}    # Dictionary of {trans:gene}
            seqcmd[-1] = 'seqin=%s.cdna.fasta' % self.info['Kozak']
            seqcdna = rje_seq.SeqList(self.log,seqcmd)
            sx = 0; stot = seqcdna.seqNum()
            for seq in seqcdna.seqs()[0:]:
                try:
                    self.progLog('\r#INFO','Extracting cDNA sequence info: %.2f%%' % (sx/stot)); sx += 100.0
                    if goodacc and seq.info['AccNum'] not in goodacc: seqcdna.seq.remove(seq); continue
                    seq.info['Gene'] = rje.matchExp('gene:(\S+)',seq.info['Name'])[0]
                    seq.info['DBase'] = 'EnsEMBL'
                    t2g[seq.info['AccNum']] = seq.info['Gene']
                    if seq.info['Gene'] not in g2t: g2t[seq.info['Gene']] = []
                    g2t[seq.info['Gene']].append(seq.info['AccNum'])
                    if seq.info['Gene'] in self.list['Track'] and seq.info['AccNum'] not in self.list['Track']: self.list['Track'].append(seq.info['AccNum'])
                    if seq.info['AccNum'] in self.list['Track'] and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])
                except:
                    self.errorLog('Problem extracting info for "%s"' % seq.info['FullName'])
                    seq.info['AccNum'] = '!ERR!'
                    seqcdna.seq.remove(seq); continue
            self.printLog('\r#INFO','Extracted cDNA sequence info: %s -> %s retained.' % (rje.integerString(stot),rje.integerString(seqcdna.seqNum())))
            seqcdna = seqcdna.seqNameDic('AccNum')
            if '!ERR!' in seqcdna: seqcdna.pop('!ERR!')
            ## ~ [1e] ~ Summarise input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','#~~# %s Input #~~#' % self.info['Kozak'],timeout=False,log=False)
            self.printLog('#CDS','%s CDS sequences' % rje.integerString(len(seqcds)))
            self.printLog('#UTR','%s 5\' UTR sequences' % rje.integerString(len(seq5utr)))
            self.printLog('#CDNA','%s cDNA sequences' % rje.integerString(len(seqcdna)))
            self.printLog('#ENST','%s ENST IDs' % rje.integerString(len(t2g)))
            self.printLog('#GENE','%s genes' % rje.integerString(len(g2t)))
            self.printLog('#~~#','#~~#',timeout=False,log=False)

            ### ~ [2] ~ Compile Sequences for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            genewin = {}; gx = 0.0; gtot = len(g2t); missx = 0; no_atg = []; atgx = 0
            self.dict['AIC'] = {}       # Dictionary of non-ATG codons and their AccNum
            seqlist.seq = []
            no_cds = 0; errx = 0
            for gene in rje.sortKeys(g2t):
                self.progLog('\r#GENE','Processing genes: %.2f%%' % (gx/gtot)); gx += 100.0
                for enst in g2t[gene]:
                    ## Get all details for sequence into one sequence object ##
                    seq = seqcdna.pop(enst)
                    seq.info['cDNA'] = seq.info['Sequence']
                    seq.info['EnsG'] = seq.info['Gene']
                    try:
                        cds = seqcds.pop(enst)
                        if not cds.info['Sequence']: raise ValueError
                        seq.info['Gene'] = cds.info['Gene']
                        seq.info['CDS'] = cds.info['Sequence']
                        if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']) and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])                    
                    except: no_cds += 1; continue    # Not a CDS
                    if self.info['NRGene'].lower() == 'ensg': nrgene = seq.info['EnsG']
                    else: nrgene = seq.info['Gene']
                    try:
                        utr = seq5utr.pop(enst)
                        seq.info['5utr'] = utr  #.info['Sequence']
                        if seq.info['5utr']:
                            if seq.info['Sequence'][:len(seq.info['5utr'])] != seq.info['5utr']:
                                self.printLog('\r#ERR','%s 5\' UTR sequence does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['Sequence'][:len(seq.info['5utr'])])
                                raise ValueError
                            no3 = seq.info['5utr'] + seq.info['CDS']
                            if seq.info['Sequence'].find(no3) != 0:
                                self.printLog('\r#ERR','%s 5\' UTR + CDS does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['CDS'])
                                self.deBug(no3)
                                self.deBug(seq.info['Sequence'])
                                raise ValueError
                        else: self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    except ValueError: errx += 1; continue
                    except:
                        seq.info['5utr'] = ''
                        self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    ## Assess sequence and dictate eType ##
                    seq.info['utr5'] = seq.info['5utr']
                    #if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']): self.deBug(seq.details())
                    seqlist.seq.append(seq)
            self.printLog('\r#GENE','Processed %s genes: %s sequence objects' % (rje.iStr(gtot),rje.iStr(seqlist.seqNum())))
            self.printLog('#NOCDS','%s cDNA transcripts excluded due to no CDS' % rje.integerString(no_cds))
            self.printLog('#ERRX','%s cDNA transcripts excluded due to sequence errors' % rje.integerString(errx))

            ### ~ [3] ~ Save SeqList object and pickle for later time-saving ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = seqlist
            self.pickleMe('%s.loaded' % self.info['Basefile'])
            self.deBug(seqlist.seqNum())
            return seqlist
        except: self.errorLog('Feb-2010-Kozak Loading Sequence Error')
#########################################################################################################################
    def biol3050(self):    ### June 2010 analysis for poster
        '''
        This method analyses the genome of choice for good (and bad) candidates for BIOL3050 project genes. Initiation
        codons are divided into:
        - Weak = -3[UC]NNXXX[UCA]
        - Strong = -3[AG]NNXXXG, where X is AUG/CUG/ACG
        - Mid = All others
        Genes are then classified according to the following criteria for experimentation:
        - Good = Weak annotated IC and Strong eORF or tORF, 40aa <= eLen <= 200aa and not "UTR" or "Cut". 
        - UTR = Weak annotated IC and Strong eORF or tORF, 40aa <= eLen <= 200aa but eORF not within UTR
        - Cut = Weak annotated IC and Strong eORF or tORF, 40aa <= eLen <= 200aa but RECuts in 5'UTR to IC+300nt.
        - Len = Weak annotated IC and Strong eORF or tORF, 40aa > eLen > 200aa.
        - Poor = All others.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            if os.path.exists('%s.CandidateGene.tdt' % (self.info['Basefile'])):
                gdb = db.addTable('%s.CandidateGene.tdt' % (self.info['Basefile']),mainkeys=['Gene'],name='Gene')
                gdb.dropFields(['Description'])
                xdb = db.addTable('../../../../Databases/DBase_100203/HGNC/genemap_100312.data.tdt',mainkeys=['Gene'],name='XRef')
                xdb.dropFields(['Gene','Entrez','OMIM','UniProt','EnsEMBL','HPRD','EnsDesc'],inverse=True)
                xdb.renameField('EnsDesc','Description')
                xdb.dropEntriesDirect('Gene',gdb.dataKeys(),inverse=True)
                cdb = db.joinTables('CandidateGeneXRef',join=[(gdb,'Gene'),(xdb,'Gene')],newkey=['Gene'],cleanup=True,delimit='\t',empties=True,check=True)
                cdb.fillBlanks()
                ## ~ [0a] ~ Annotate for BIOL3050 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Good x 2	Genes with one Good eORF or tORF.
                # Biol x 4	Genes with no Good eORFs but UTR/Len or Cuts.
                # Weak x 2	Genes with weak starts but no Strong alternatives.
                # Poor x 2	Genes with reasonable start codons. Possible make these poorly characterised genes - bonus work can be to try and work out what they are?

                # Revised Rules for Candidate Selection	
                # Poor/Weak x 4	Too few genes with weak starts and no strong alternatives
                # Total number of transcripts and eORFs should be limited for each student.	
                # PLUS	Genes with 2+ Good eORFs/tORFs AND 2+ transcripts. (Can have 2+ in same transcript.)
                cdb.addField('BIOL3050','Poor'); cdb.addField('Notes','BIOL3050')
                cdb.addField('eORF','EnsG')
                for entry in cdb.entries(): entry['eORF'] = 0
                reformat = {}
                for f in string.split('ENST	EnsG	eORF	Good	UTR	Len	Cut	Weak	Poor'): reformat[f] = 'int'
                cdb.dataFormat(reformat)
                for entry in cdb.entries():
                    for f in string.split('Good	UTR	Len	Cut	Weak	Poor'): entry['eORF'] += entry[f]
                    if entry['Good'] > 1: entry['BIOL3050'] = 'Good'
                    elif (entry['UTR'] + entry['Len']) > 1: entry['BIOL3050'] = 'Biol'
                    elif entry['Weak'] > 1: entry['BIOL3050'] = 'Weak'
                    else: entry['BIOL3050'] = 'Poor'
                    if entry['eORF'] > 10: entry['Notes'] = 'eORF > 10'
                    elif entry['ENST'] > 5: entry['Notes'] = 'Transcript > 5'
                    elif not entry['Description']: entry['Notes'] = 'Annotation'
                    else: entry['Notes'] = 'Candidate'
                cdb.saveToFile('%s.%s.tdt' % (self.info['Basefile'],cdb.info['Name']))
                cdb.addField('Class','Notes')
                for entry in cdb.entries():
                    if entry['Notes'] == 'Annotation': entry['Class'] = 'Reject'
                    else: entry['Class'] = entry['BIOL3050']
                cdb.saveToFile('%s.biol3050.tdt' % (self.info['Basefile']))
                return
            self.obj['SeqList'] = seqlist = self._biol3050_loadSeq()
            edb = db.addTable('%s.eorf.tdt' % self.info['Kozak'],mainkeys=['eORF'],name='eorf')
            fieldformat = {}    #'Dataset','RunID','Masking','Build','RunTime','Pattern',
            for f in ['eLen','eUTR','Score']: fieldformat[f] = 'int'
            edb.dataFormat(fieldformat)

            ### ~ [1] ~ Reclassify ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edb.addField('Dirn',after='eORF'); edb.addField('Class',after='eORF')
            edb.addField('Codon',after='CoreType'); edb.addField('NewCodon',after='NewCoreType')
            ex = 0.0; etot = edb.entryNum(); x3 = 0
            for entry in edb.entries():
                self.progLog('\r#EDB','Updating eORF data: %.2f%%' % (ex/etot)); ex += 100.0
                for c in ['','New']:
                    if entry['%sCoreType' % c] in ['Good','Alt']: entry['%sCoreType' % c] = 'Strong'
                    #elif entry['%sCoreType' % c] == 'Rxx' and entry['%sCore' % c][-1] == 'C': entry['%sCoreType' % c] = 'Weak'
                    elif entry['%sCore' % c][0] in ['U','C'] and entry['%sCore' % c][-1] != 'G': entry['%sCoreType' % c] = 'Weak'
                    else: entry['%sCoreType' % c] = 'Mid'
                    entry['%sCodon' % c] = entry['%sCore' % c][3:-1]
                if string.split(entry['eORF'],'-')[1][0] == 't':
                    entry['Dirn'] = '3'; x3 += 1
                    entry['eUTR'] = entry['eUTR'] + 2 * entry['eLen']
                else: entry['Dirn'] = '5'
                entry['Class'] = 'Poor'
                if entry['CoreType'] == 'Weak' and entry['NewCoreType'] == 'Strong':
                    if entry['eUTR'] < 0: entry['Class'] = 'UTR'
                    elif 40 <= entry['eLen'] <= 200: entry['Class'] = 'Good'
                    else: entry['Class'] = 'Len'   #!# Need to check RECut
                elif entry['CoreType'] == 'Weak': entry['Class'] = 'Weak'
            self.printLog('\r#EDB','Updating eORF data complete. %s 3\' AIC' % rje.integerString(x3))

            ## ~ [2] ~ Restriction Enzyme Cutting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqdict = seqlist.seqNameDic(key='AccNum',proglog=True)
            self.deBug(rje.sortKeys(seqdict)[:20])
            ex = 0.0; etot = edb.entryNum(); rx = 0
            for entry in edb.entries():
                self.progLog('\r#RE','Check restriction sites: %.2f%%' % (ex/etot)); ex += 100.0
                if entry['Class'] != 'Good': continue
                enst = entry['ENST']
                try: seq = seqdict[enst]
                except: entry['Class'] = 'ERR'; continue
                if entry['Dirn'] == '5': aic3 = entry['eUTR'] + entry['eLen'] * 3
                else: aic3 = entry['eUTR']
                clone = seq.info['Sequence'][:aic3 + 303]
                for recut in self.list['RECut']:
                    if recut.upper() in clone.upper(): entry['Class'] = 'Cut'; rx += 1; break
            self.printLog('\r#RE','Check restriction sites complete: %s cut.' % (rx))
        
            ### ~ [3] ~ Split and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            edb.saveToFile('%s.Full.tdt' % (self.info['Basefile']))
            splits = self.db().splitTable(edb,'Class',asdict=True)
            for split in splits: splits[split].saveToFile('%s.%s.tdt' % (self.info['Basefile'],split))

            ### ~ [4] ~ Output full Alt Initiation table for Good Candidates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            classes = ['Good','UTR','Len','Cut','Weak','Poor']
            cfields = ['ENST','EnsG','Gene','Description'] 
            tdb = self.db().addEmptyTable('CandidateENST',cfields+classes,['ENST'])
            gdb = self.db().addEmptyTable('CandidateENSG',cfields[1:]+cfields[:1]+classes,['EnsG'])
            hdb = self.db().addEmptyTable('CandidateGene',cfields[2:]+cfields[:2]+classes,['Gene'])
            ex = 0.0; etot = edb.entryNum(); 
            for entry in edb.entries():
                self.progLog('\r#CAND','Establishing Candidate Genes/Transcripts: %.2f%%' % (ex/etot)); ex += 100.0
                ensg = entry['EnsG']; enst = entry['ENST']; gene = entry['Gene']
                if enst not in tdb.data():
                    tdb.data()[enst] = {}
                    for f in cfields: tdb.data()[enst][f] = entry[f]
                    for c in classes: tdb.data()[enst][c] = 0
                tdb.data()[enst][entry['Class']] += 1
                if ensg not in gdb.data():
                    gdb.data()[ensg] = {'ENST': len(edb.indexDataList('EnsG',ensg,'ENST',sortunique=True))}
                    for f in cfields[1:]: gdb.data()[ensg][f] = entry[f]
                    for c in classes: gdb.data()[ensg][c] = 0
                gdb.data()[ensg][entry['Class']] += 1
                if gene not in hdb.data():
                    hdb.data()[gene] = {'ENST': len(edb.indexDataList('Gene',gene,'ENST',sortunique=True)),
                                        'EnsG': len(edb.indexDataList('Gene',gene,'EnsG',sortunique=True))}
                    for f in cfields[2:]: hdb.data()[gene][f] = entry[f]
                    for c in classes: hdb.data()[gene][c] = 0
                hdb.data()[gene][entry['Class']] += 1
            tx = len(edb.indexDataList('Class','Good','ENST'))
            gx = len(edb.indexDataList('Class','Good','EnsG'))
            hx = len(edb.indexDataList('Class','Good','Gene'))
            self.printLog('\r#CAND','Established Candidate Genes/Transcripts: %s enst; %s ensg; %s genes.' % (rje.iStr(tx),rje.iStr(gx),rje.iStr(hx)))            
            for table in [tdb,gdb,hdb]: table.saveToFile('%s.%s.tdt' % (self.info['Basefile'],table.info['Name']))
            edb.dropEntriesDirect('Class','Poor')
            edb.saveToFile('%s.Candidates.tdt' % (self.info['Basefile']))
        except: self.errorLog('Ugg')
#########################################################################################################################
    def jun10(self,ptype=None,plus=False):    ### June 2010 analysis for poster   #!# NB. had errors - re-run! #!#
        '''June 2010 analysis for poster.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodlist = ['Good','Alt']
            if plus: goodlist = ['Good','Alt','augG','Rxx']
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            edb = db.addTable('Homo_sapiens.GRCh37.56.eorf.tdt',mainkeys=['eORF'],name='eorf')
            fieldformat = {}    #'Dataset','RunID','Masking','Build','RunTime','Pattern',
            for f in ['eLen','eUTR','Score']: fieldformat[f] = 'int'
            edb.dataFormat(fieldformat)
            ## ~ [0a] ~ Reclassify ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb.addField('Dirn'); edb.addField('Keep')
            ex = 0.0; etot = edb.entryNum(); x3 = 0
            for entry in edb.entries():
                self.progLog('\r#EDB','Updating eORF data: %.2f%%' % (ex/etot)); ex += 100.0
                if entry['CoreType'] == 'Rxx' and entry['Core'][-1] == 'C': entry['CoreType'] = 'Poor' 
                if entry['NewCoreType'] == 'Rxx' and entry['NewCore'][-1] == 'C': entry['NewCoreType'] = 'Poor'
                if string.split(entry['eORF'],'-')[1][0] == 't':
                    entry['Dirn'] = '3'; x3 += 1
                    entry['eUTR'] = entry['eUTR'] + 2 * entry['eLen']
                else: entry['Dirn'] = '5'
                entry['Keep'] = True
            self.printLog('\r#EDB','Updating eORF data complete. %s 3\' AIC' % rje.integerString(x3))
            tdb = db.addTable('Homo_sapiens.GRCh37.56.types.tdt',mainkeys=['CoreType'],name='type')
            for field in tdb.fields():
                if field != 'CoreType': fieldformat[field] = 'int'
            tdb.dataFormat(fieldformat)
            for entry in edb.entries():
                entry['Keep'] = not ptype
                entry['Keep'] = entry['Keep'] or entry['NewCore'][-4:-1] == ptype 
            edb.dropEntriesDirect('Keep',[False],log=True)
            ## ~ [0b] ~ Reduce data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb.dropEntriesDirect('Dirn',['3'])
            #x#edb.dropEntriesDirect('CoreType',['Poor'],inverse=True,log=True)
            edb.dropEntriesDirect('NewCoreType',goodlist,inverse=True,log=True)
            ## ~ [0c] ~ Type Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tdata = tdb.data()['All']
            ### ~ [1] ~ Calculate transcript data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = len(edb.index('ENST'))
            aic = {0:{'Len':0,'eAIC':0,'uAIC':0,'eORF':0}}
            for enst in edb.index('ENST'):
                self.progLog('\r#AIC','Generating AIC summary data: %.2f%%' % (ex/etot)); ex += 100.0
                best = [1e6,1e6,0]  # e / u / e 
                for entry in edb.indexEntries('ENST',enst):
                    self.deBug(entry)
                    if entry['Dirn'] == '5': best[0] = min(entry['eLen'],best[0])
                    if entry['Dirn'] == '5' and entry['eUTR'] > 0:
                        best[1] = min(entry['eLen'],best[1])
                        best[2] = max(entry['eLen'],best[2])
                self.deBug('%s: %s' % (enst,best))
                for i in best:
                    if i not in aic: aic[i] = {'Len':i,'eAIC':0,'uAIC':0,'eORF':0}
                aic[best[0]]['eAIC'] += 1
                aic[best[1]]['uAIC'] += 1
                aic[best[2]]['eORF'] += 1
            self.printLog('\r#AIC','Generation of AIC summary data complete.')
            ## ~ [1a] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ahead = ['Len','eAIC','uAIC','eORF','peAIC','puAIC','peORF']
            if not ptype: afile = '%s.ok_eorf.tdt' % self.info['Basefile']
            else: afile = '%s.ok_eorf.%s.tdt' % (self.info['Basefile'],ptype)
            if ptype or plus: db.list['Tables'].remove(edb)
            if plus: afile = '%s.plus.tdt' % afile[:-4]
            rje.delimitedFileOutput(self,afile,ahead,rje_backup=True)
            p = {}
            for e in ['eAIC','uAIC','eORF']: p['p%s' % e] = 0.0
            #if aic.has_key(1e6): aic[0] = aic.pop(1e6)
            for i in rje.sortKeys(aic,revsort=True):
                p['peORF'] += (float(aic[i]['eORF']) / tdata['XX'])
                aic[i]['peORF'] = p['peORF']
            aic[0]['peORF'] = aic[rje.sortKeys(aic)[1]]['peORF'] 
            for i in rje.sortKeys(aic):
                for e in ['eAIC','uAIC']:
                    p['p%s' % e] += (float(aic[i][e]) / tdata['XX'])
                    aic[i]['p%s' % e] = p['p%s' % e]
                #aic[i]["p5'"] = float(aic[i]["5'"]) / tdata['AX']
                #aic[i]["p3'"] = float(aic[i]["3'"]) / tdata['XX']
                #aic[i]["pAIC"] = float(aic[i]["AIC"]) / tdata['XX']
                rje.delimitedFileOutput(self,afile,ahead,datadict=aic[i])
            self.printLog('\r#AIC','AIC summary data output to %s.' % afile)
        except: self.errorLog('Ugg')
#########################################################################################################################
    def OLDjun10(self,ptype=None,plus=False):    ### June 2010 analysis for poster
        '''June 2010 analysis for poster.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            goodlist = ['Good','Alt']
            if plus: goodlist = ['Good','Alt','augG','Rxx']
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            edb = db.addTable('Homo_sapiens.GRCh37.56.eorf.tdt',mainkeys=['eORF'],name='eorf')
            fieldformat = {}    #'Dataset','RunID','Masking','Build','RunTime','Pattern',
            for f in ['eLen','eUTR','Score']: fieldformat[f] = 'int'
            edb.dataFormat(fieldformat)
            ## ~ [0a] ~ Reclassify ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb.addField('Dirn'); edb.addField('Keep')
            ex = 0.0; etot = edb.entryNum(); x3 = 0
            for entry in edb.entries():
                self.progLog('\r#EDB','Updating eORF data: %.2f%%' % (ex/etot)); ex += 100.0
                if entry['CoreType'] == 'Rxx' and entry['Core'][-1] == 'C': entry['CoreType'] == 'Poor' 
                if entry['NewCoreType'] == 'Rxx' and entry['NewCore'][-1] == 'C': entry['NewCoreType'] == 'Poor'
                if string.split(entry['eORF'],'-')[1][0] == 't':
                    entry['Dirn'] = '3'; x3 += 1
                    entry['eUTR'] = entry['eUTR'] + 2 * entry['eLen']
                else: entry['Dirn'] = '5'
                entry['Keep'] = True
            self.printLog('\r#EDB','Updating eORF data complete. %s 3\' AIC' % rje.integerString(x3))
            tdb = db.addTable('Homo_sapiens.GRCh37.56.types.tdt',mainkeys=['CoreType'],name='type')
            for field in tdb.fields():
                if field != 'CoreType': fieldformat[field] = 'int'
            tdb.dataFormat(fieldformat)
            for entry in edb.entries():
                entry['Keep'] = not ptype
                entry['Keep'] = entry['Keep'] or entry['NewCore'][-4:-1] == ptype 
            edb.dropEntriesDirect('Keep',[False],log=True)
            ## ~ [0b] ~ Reduce data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edb.index('Dirn'); self.deBug(edb.index('Dirn').keys())
            edb.index('CoreType'); self.deBug(edb.index('CoreType').keys())
            edb.dropEntriesDirect('CoreType',['Poor'],inverse=True,log=True)
            edb.index('CoreType'); self.deBug(edb.index('CoreType').keys())
            edb.index('Dirn'); self.deBug(edb.index('Dirn').keys())
            edb.index('NewCoreType'); self.deBug(edb.index('NewCoreType').keys())
            edb.dropEntriesDirect('NewCoreType',goodlist,inverse=True,log=True)
            edb.index('NewCoreType'); self.deBug(edb.index('NewCoreType').keys())
            edb.index('Dirn'); self.deBug(edb.index('Dirn').keys())
            ## ~ [0c] ~ Type Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tdata = tdb.data()['Poor']
            ### ~ [1] ~ Calculate transcript data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ex = 0.0; etot = len(edb.index('ENST'))
            aic = {0:{'Len':0,'eAIC':0,'uAIC':0,'tAIC':0,'AIC':0}}
            for enst in edb.index('ENST'):
                self.progLog('\r#AIC','Generating AIC summary data: %.2f%%' % (ex/etot)); ex += 100.0
                best = [1e6,1e6,1e6]  # e / u / t 
                for entry in edb.indexEntries('ENST',enst):
                    self.deBug(entry)
                    if entry['Dirn'] == '5': best[0] = min(entry['eLen'],best[0])
                    if entry['Dirn'] == '5' and entry['eUTR'] > 0: best[1] = min(entry['eLen'],best[1])
                    if entry['Dirn'] == '3': best[2] = min(entry['eLen'],best[2])
                self.deBug('%s: %s' % (enst,best))
                for i in best:
                    if i not in aic: aic[i] = {'Len':i,'eAIC':0,'uAIC':0,'tAIC':0,'AIC':0}
                aic[best[0]]['eAIC'] += 1
                aic[best[1]]['uAIC'] += 1
                aic[best[2]]['tAIC'] += 1
                aic[min(best)]['AIC'] += 1
            self.printLog('\r#AIC','Generation of AIC summary data complete.')
            ## ~ [1a] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ahead = ['Len','eAIC','uAIC','tAIC','AIC','peAIC','puAIC','ptAIC','pAIC']
            if not ptype: afile = '%s.poor_aic.tdt' % self.info['Basefile']
            else: afile = '%s.poor_aic.%s.tdt' % (self.info['Basefile'],ptype)
            if ptype or plus: db.list['Tables'].remove(edb)
            if plus: afile = '%s.plus.tdt' % afile[:-4]
            rje.delimitedFileOutput(self,afile,ahead,rje_backup=True)
            p = {}
            for e in ['eAIC','uAIC','tAIC','AIC']: p['p%s' % e] = 0.0
            #if aic.has_key(1e6): aic[0] = aic.pop(1e6)
            for i in rje.sortKeys(aic):
                for e in ['eAIC','uAIC','tAIC','AIC']:
                    p['p%s' % e] += (float(aic[i][e]) / tdata['XX'])
                    aic[i]['p%s' % e] = p['p%s' % e]
                #aic[i]["p5'"] = float(aic[i]["5'"]) / tdata['AX']
                #aic[i]["p3'"] = float(aic[i]["3'"]) / tdata['XX']
                #aic[i]["pAIC"] = float(aic[i]["AIC"]) / tdata['XX']
                rje.delimitedFileOutput(self,afile,ahead,datadict=aic[i])
            self.printLog('\r#AIC','AIC summary data output to %s.' % afile)
        except: self.errorLog('Ugg')
#########################################################################################################################
    ### <5> ### PRIDE analysis methods                                                                                  #
#########################################################################################################################
    def pride(self):    ### Performs basic PRIDE analysis. (Make more complex with time.)
        '''Performs basic PRIDE analysis. (Make more complex with time.).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['PRIDE'] = rje.dataDict(self,self.info['PrideFile'],['PRIDE Experiment Accession'],['Ensembl Human'],lists=True)
            #self.downloadPRIDE()  wget 'ftp://ftp.ebi.ac.uk/pub/databases/pride/PRIDE_Exp_Complete_Ac_%s.xml.gz' % id
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            ### ~ [1] ~ Extract Peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            px = 0.0; ptot = len(self.dict['PRIDE']); success = 0; pepx = 0; genex = 0
            self.printLog('\r#PRIDE','%d PRIDE files to parse.' % (ptot))
            for pid in rje.sortKeys(self.dict['PRIDE']):
                self.progLog('\r#PARSE','Parsing Peptides from PRIDE files: %.2f%%' % (px/ptot)); px += 100.0
                self.dict['PRIDE'][pid]['Peptides'] = []
                if self.opt['Test'] and genex > 100 and pepx > 100 and px >= 1500: continue
                gfile = '%sPRIDE_%s.grep' % (self.info['PridePath'],pid)
                if not os.path.exists(gfile):
                    pfile = '%sPRIDE_Exp_Complete_Ac_%s.xml' % (self.info['PridePath'],pid)
                    if os.path.exists('%s.gz' % pfile): os.system('gunzip %s.gz' % pfile)
                    if not os.path.exists('%s' % pfile): continue
                    os.system('grep Sequence %s > %s' % (pfile,gfile))
                    os.system('gzip %s' % pfile)
                for xline in open(gfile,'r').readlines():
                    try: pep = rje.matchExp('<Sequence>(\S+)</Sequence>',xline)[0]
                    except: continue
                    #if pep not in self.dict['PRIDE'][pid]['Peptides']:  #?# Do we want this? #?#
                    self.dict['PRIDE'][pid]['Peptides'].append(pep)
                self.dict['PRIDE'][pid]['Peptides'].sort()
                if self.dict['PRIDE'][pid]['Peptides']: success += 1; pepx += len(self.dict['PRIDE'][pid]['Peptides'])
                genex += len(self.dict['PRIDE'][pid]['Ensembl Human'])
            self.printLog('\r#PARSE','Parsed %s Peptides from %d of %d PRIDE files.' % (rje.integerString(pepx),success,ptot))
            ## ~ [1a] ~ Make peptide DBase table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = db.addEmptyTable('Peptides',['#','PID','Peptide'],keys=['#'])
            gdb = db.addEmptyTable('PRIDE',['#','PID','EnsHUMAN'],keys=['#'])
            px = 0; gx = 0
            for pid in rje.sortKeys(self.dict['PRIDE']):
                for pep in self.dict['PRIDE'][pid]['Peptides']:
                    px += 1
                    entry = {'#':rje.preZero(px,pepx),'PID':pid,'Peptide':pep}
                    pdb.data()[entry['#']] = entry
                for ens in self.dict['PRIDE'][pid]['Ensembl Human']:
                    gx += 1
                    entry = {'#':rje.preZero(gx,genex),'PID':pid,'EnsHUMAN':ens}
                    gdb.data()[entry['#']] = entry
            pdb.index('PID'); gdb.index('PID'); gdb.index('EnsHUMAN')
            #pdb.saveToFile('%s.peptides.tdt' % self.info['Basefile'])
            #gdb.saveToFile('%s.pride.tdt' % self.info['Basefile'])
                    
            # Load mappings from ...EnsEMBL/ens_HUMAN.map.100505.tdt
            mdb = db.addTable('/scratch/RJE_Filestore/SBSBINF/Databases/DBase_100505/EnsEMBL/ens_HUMAN.map.100505.tdt',mainkeys=['Ensembl Protein ID'],name='Map')
            self.deBug(mdb.fields())
            mdb.renameField('HGNC symbol','Gene')
            mdb.renameField('Ensembl Protein ID','EnsHUMAN')    # Key
            mdb.index('Gene')
            # Load "known" eORFs from *.eorf.tdt
            etdt = '../2010-02-11-Initiation_Codons/Homo_sapiens.GRCh37.56.eorf.tdt'
            etdt = '../2010-02-11-Initiation_Codons/2010-06-02_tORF/Homo_sapiens.GRCh37.56.eorf.tdt'
            edb = db.addTable(etdt,mainkeys=['eORF'],name='eORF')
            edb.index('Gene')
            # Load eORF sequences from *.eorf.fas
            try:
                eseqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=%s.fas' % etdt[:-4],'seqnr=F','accnr=F,usecase=T'])
                edict = eseqlist.seqNameDic()
                if not edict: raise ValueError
            except:
                edict = {}
                eseqlist = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None','seqnr=F','accnr=F,usecase=T'])
                blast = rje_blast.BLASTRun(self.log,self.cmd_list)
                blast.formatDB('%s.fas' % etdt[:-4])

            # Map known eORFs onto Genes on Peptides onto PRIDE experiments
            edb.addField('ePepN'); edb.addField('ePepX'); edb.addField('iPepN'); edb.addField('iPepX');
            edb.addField('eN~iN'); edb.addField('eX~iX'); edb.addField('PRIDE')
            ex = 0.0; etot = edb.entryNum()
            for eorf in edb.datakeys():
                self.progLog('\r#EPEP','Mapping eORF peptides: %.2f%%' % (ex/etot)); ex += 100.0
                try:
                    entry = edb.data()[eorf]; entry['PRIDE'] = 'na'
                    ## ~ Setup peptides ~ ##
                    try:
                        if edict: eseq = edict[eorf].getSequence(case=True)
                        else: eseq = eseqlist.seqFromFastaCmd(eorf,'%s.fas' % etdt[:-4]).getSequence(case=True)
                    except: continue
                    for a in 'krKR': eseq = string.replace(eseq,a,'%s ' % a)
                    peplist = string.split(eseq)
                    if len(peplist) < 2: continue
                    if peplist[1][0] == peplist[1].lower()[0]: epep = peplist[1]
                    else: continue
                    ipep = 'na'
                    for pep in peplist[2:]:
                        if pep[0] in ['M',pep[0].lower()]: continue
                        ipep = pep; break
                    if min(len(ipep),len(epep)) < 5: continue
                    ## ~ Setup PRIDE list ~ ##
                    gene = entry['Gene']
                    pidlist = []
                    self.deBug(gene)
                    self.deBug(mdb.dataList(mdb.indexEntries('Gene',gene),'EnsHUMAN',sortunique=True))
                    for ensp in mdb.dataList(mdb.indexEntries('Gene',gene),'EnsHUMAN',sortunique=True):
                        self.deBug(gdb.indexEntries('EnsHUMAN',ensp))
                        for gentry in gdb.indexEntries('EnsHUMAN',ensp):
                            if gentry['PID'] not in pidlist: pidlist.append(gentry['PID'])
                    entry['PRIDE'] = len(pidlist)
                    if not pidlist: continue
                    ## ~ Map peptides to PRIDE ~ ##
                    epepn = []; epepx = 0; ipepn = []; ipepx = 0
                    for pid in pidlist:
                        for pep in self.dict['PRIDE'][pid]['Peptides']:
                            if pep.upper() == epep.upper():
                                if pid not in epepn: epepn.append(pid)
                                epepx += 1
                            if pep.upper() == ipep.upper():
                                if pid not in ipepn: ipepn.append(pid)
                                ipepx += 1
                    entry['ePepX'] = epepx
                    entry['ePepN'] = len(epepn)
                    entry['iPepX'] = ipepx
                    entry['iPepN'] = len(ipepn)
                    entry['eX~iX'] = 'na'
                    entry['eN~iN'] = 'na'
                    if max(ipepx,epepx) == 0: continue
                    entry['eX~iX'] = float(epepx) / (ipepx + epepx)
                    entry['eN~iN'] = float(len(epepn)) / (len(ipepn) + len(epepn))                
                except: self.errorLog('%s problem' % eorf)
            self.printLog('\r#EPEP','Mapping eORF peptides complete.')
            edb.saveToFile('%s.eorf-pride.tdt' % self.info['Basefile'])
            
            # Look for eORF peptides versus ORF peptides in PRIDE peptide list
            # => % experiments returning AUG -1 and AUG +1
            # eORF  ExpNum eFrag aFrag eRatio(e/a)
            # => Bad eORF starts should have lower e/a on average
            return
        except: self.errorLog('RJE_AIC.PRIDE error')
#########################################################################################################################
    ### <5> ### PRIDE analysis methods                                                                                  #
#########################################################################################################################
    def aug11(self):    ### Temp experimental method
        '''Temp experimental method.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#ZEN',rje_zen.Zen().wisdom())
            db = self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
            inbase = self.getStr('EnsPath') + self.getStr('Kozak')
            seqcmd = self.cmd_list + ['seqmode=file','seqindex=T','autoload=T']
            cdna = rje_seqlist.SeqList(self.log,seqcmd+['seqin=%s.cdna.all.fa' % inbase])
            pep = rje_seqlist.SeqList(self.log,seqcmd+['seqin=%s.pep.all.fa' % inbase])
            chrom = rje_seqlist.SeqList(self.log,seqcmd+['seqin=%s.chromosome.all.fa' % inbase])
            ### ~ [1] ~ Generate Flanking Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdbfile = self.getStr('Kozak') + '.cdna.tdt'
            ## ~ [1a] ~ Check for  dictionary of sequences on each chromosome ~~~~~~~~~~~~~~~~~~~~~ ##
            if os.path.exists(cdbfile) and not self.getBool('Force'):
                cdb = db.addTable(cdbfile,mainkeys=['ENST'],name='cdna')
            ## ~ [1b] ~ Make a dictionary of sequences on each chromosome ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                cdb = db.addEmptyTable('cdna',['ENST','Type','Csome','Pos5','Pos3','Strand','ESNG','Len'],keys=['ENST'])
                csomedict = {}
                cdna.setStr({'SeqDictType':'short'})
                while cdna.nextSeq():
                    self.progLog('\r#CSOME','Assigning sequences to chromosomes: %s' % cdna.progress())
                    seq = cdna.currSeq()
                    self.deBug(seq[0])
                    data = string.split(seq[0])
                    locus = string.split(data[2],':')
                    entry = {'ENST':data[0],'Type':string.split(data[1],':')[1],'Csome':locus[2],'Pos5':locus[3],
                             'Pos3':locus[4],'Strand':locus[5],'ESNG':string.split(data[3],':')[1],'Len':len(seq[1])}
                    cdb.addEntry(entry)
                    cdna.dict['SeqDict'][data[0]] = cdna.obj['Current']
                self.printLog('\r#CSOME','Assignment of %s sequences to chromosomes complete' % rje.iStr(cdb.entryNum()))
                cdb.saveToFile(cdbfile)
            ## ~ [1c] ~ Format Table for flank generation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            cdb.dataFormat({'ENST':'str','Type':'str','Csome':'str','Pos5':'int','Pos3':'int','Strand':'str','ESNG':'str','Len':'int'})
            cdb.index('Csome')
            ## ~ [1d] ~ Generate flanks, one chromosome at a time ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            FLANK = open('%s.5flank.fas' % self.getStr('Kozak'),'w')
            UNSPLICED = open('%s.unspliced.fas' % self.getStr('Kozak'),'w')
            while chrom.nextSeq():
                self.progLog('\r#FLANK','Making flanking sequences: %s' % chrom.progress())
                cseq = chrom.currSeq()
                csome = string.split(cseq[0])[0]
                for entry in cdb.indexEntries('Csome',csome):
                    name = cdna.getDictSeq(entry['ENST'])[0]
                    pos5 = entry['Pos5']
                    pos3 = entry['Pos3']
                    strand = entry['Strand']
                    if strand == '1':
                        flank = cseq[1][max(0,pos5-self.getInt('Flank')):pos5]
                        unspliced =  cseq[1][pos5+1:pos3]
                    else:
                        flank = rje_sequence.reverseComplement(cseq[1][pos3+1:pos3+self.getInt('Flank')+1])
                        unspliced =  rje_sequence.reverseComplement(cseq[1][pos5+1:pos3])
                    FLANK.write('>%s\n%s\n' % (name,flank))
                    UNSPLICED.write('>%s\n%s\n' % (name,unspliced))
            FLANK.close(); UNSPLICED.close()
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    ### <6> ### Jan14 analysis methods                                                                                  #
#########################################################################################################################
    def jan14(self):    ### Perform revised Feb '10 Kozak style analysis
        '''This mode constitutes a revised Kozak-style prediction of human TIC (Translation Initiation Codon).'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['GoodIC'] = string.split(string.replace(string.join(self.list['GoodIC']).upper(),'U','T'))
            pickle = self.kozak10_pickle()
            if not self.obj['SeqList']: self._kozak14_loadSeq()



            if self.list['Output']: self._kozak10_headers()
            self._kozak10_outputTypes()
            ### ~ [2] ~ Perform Calculations for both TICwin and Conwin sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Nucleotide and Core sequence counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #pickle = self.unpickleMe('%s.counts' % self.info['Basefile'])
            #if pickle: self.dict['Data'] = pickle.dict['Data']
            #else: self.dict['Data'] = self._kozak10_counts()
            if not self.dict['Data']: self._kozak10_counts()
            ## ~ [2b] ~ Rank Cores within each eType & cType & score each sequence ~~~~~~~~~~~~~~~~ ##
            #pickle = self.unpickleMe('%s.score' % self.info['Basefile'])
            #if pickle:
            #    self.obj['SeqList'] = pickle.obj['SeqList']
            #    self.dict['Data'] = pickle.dict['Data']
            #else: self.obj['SeqList'] = self._kozak10_scoring()
            if self.info['Pickle'] not in ['score']: self._kozak10_scoring()
            ## ~ [2c] ~ eORF analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'eorf' in self.list['Output']: self._kozak10_eORFs()
            ## ~ [2d] ~ AL core analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if 'alcore' in self.list['Output']: self._kozak10_randCoreAL()
            ### ~ [3] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'enst' in self.list['Output']: self._kozak10_outputENST()
            if 'core' in self.list['Output']: self._kozak10_outputCore('core')
            if 'wegcore' in self.list['Output']: self._kozak10_outputCore('wegcore')
            if 'alcore' in self.list['Output']: self._kozak10_outputALCore()
            if 'pos' in self.list['Output']: self._kozak10_outputPos()

        except: self.errorLog('Feb-2010-Kozak Error')
#########################################################################################################################
    def _kozak14_loadSeq(self):     ### Loads sequences and peforms initial processing before pickling
        '''Loads sequences and peforms initial processing before pickling.'''
        try:### ~ [1] ~ Load Sequences for processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Data is read in from four sources:
            # 1. Full length EnsEMBL cDNAs
            # 2. 5' flank of TSS. This should not include 5' UTR where annotation is good.
            # 3. 5'UTR sequences.
            # 4. CDS.
            seqcmd = self.cmd_list + ['autoload','accnr=F','seqnr=F','gnspacc=F','seqin=None']
            self.obj['SeqList'] = seqlist = rje_seq.SeqList(self.log,seqcmd)  # This object will store all the final sequence data
            ## ~ [1a] ~ CDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqcmd[-1] = 'seqin=%s.cds.fasta' % self.info['Kozak']
            seqcds = rje_seq.SeqList(self.log,seqcmd)
            self._kozak10_cleanSeq(seqcds,type='CDS',goodacc=[])
            seqcds = seqcds.seqNameDic('AccNum')
            if '!ERR!' in seqcds: seqcds.pop('!ERR!')
            goodacc = rje.sortKeys(seqcds)
            open('%s.cds.enst' % self.info['Kozak'],'w').write(string.join(goodacc,'\n'))
            self.printLog('#ACC','%s CDS enst output for future reduction' % rje.integerString(len(goodacc)))
            ## ~ [1b] ~ 5' Flanks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seq5flank = self._kozak10_cleanLoadSeq('5flank',goodacc,minlen=self.stat['Flank'])
            ## ~ [1c] ~ 5' UTR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seq5utr = self._kozak10_cleanLoadSeq('5utr',goodacc)
            ## ~ [1d] ~ cDNA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            g2t = {}    # Dictionary of {gene:[transcripts]}
            t2g = {}    # Dictionary of {trans:gene}
            seqcmd[-1] = 'seqin=%s.cdna.fasta' % self.info['Kozak']
            seqcdna = rje_seq.SeqList(self.log,seqcmd)
            sx = 0; stot = seqcdna.seqNum()
            for seq in seqcdna.seqs()[0:]:
                try:
                    self.progLog('\r#INFO','Extracting cDNA sequence info: %.2f%%' % (sx/stot)); sx += 100.0
                    if goodacc and seq.info['AccNum'] not in goodacc: seqcdna.seq.remove(seq); continue
                    seq.info['Gene'] = rje.matchExp('gene:(\S+)',seq.info['Name'])[0]
                    seq.info['DBase'] = 'EnsEMBL'
                    t2g[seq.info['AccNum']] = seq.info['Gene']
                    if seq.info['Gene'] not in g2t: g2t[seq.info['Gene']] = []
                    g2t[seq.info['Gene']].append(seq.info['AccNum'])
                    if seq.info['Gene'] in self.list['Track'] and seq.info['AccNum'] not in self.list['Track']: self.list['Track'].append(seq.info['AccNum'])
                    if seq.info['AccNum'] in self.list['Track'] and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])
                except:
                    self.errorLog('Problem extracting info for "%s"' % seq.info['FullName'])
                    seq.info['AccNum'] = '!ERR!'
                    seqcdna.seq.remove(seq); continue
            self.printLog('\r#INFO','Extracted cDNA sequence info: %s -> %s retained.' % (rje.integerString(stot),rje.integerString(seqcdna.seqNum())))
            seqcdna = seqcdna.seqNameDic('AccNum')
            if '!ERR!' in seqcdna: seqcdna.pop('!ERR!')
            ## ~ [1e] ~ Summarise input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#~~#','#~~# %s Input #~~#' % self.info['Kozak'],timeout=False,log=False)
            self.printLog('#CDS','%s CDS sequences' % rje.integerString(len(seqcds)))
            self.printLog('#FLANK','%s 5\' flank sequences' % rje.integerString(len(seq5flank)))
            self.printLog('#UTR','%s 5\' UTR sequences' % rje.integerString(len(seq5utr)))
            self.printLog('#CDNA','%s cDNA sequences' % rje.integerString(len(seqcdna)))
            self.printLog('#ENST','%s ENST IDs' % rje.integerString(len(t2g)))
            self.printLog('#GENE','%s genes' % rje.integerString(len(g2t)))
            self.printLog('#~~#','#~~#',timeout=False,log=False)

            ### ~ [2] ~ Compile Sequences for analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Sequences are compiled and cross-checked with each other. Where sequence <3> is over 200nt, it is assumed
            # that the TIC is at position -200 (i.e. 200 from end). Sequence <3> should also match the start of sequence
            # <1> for a given transcript. Where <3> is <200nt, the full length cDNA is used and the TIC is assumed to be
            # the first AUG. In terms of TIC context, this gives different types of sequence ('etype'):
            # - AL. 5'UTR annotated and longer than kozwin.
            # - AS. 5'UTR annotated but shorter than kozwin.
            # - US. 5'UTR inferred from cDNA but shorter than kozwin.
            # - These are also combined for additional combos (XL,XS,AX,XX) for calculations
            # - NL. 5'UTR annotated and longer than kozwin but non-AUG TIC. Not used for calculations (but still rated).
            # - NS. 5'UTR annotated but shorter than kozwin, plus non-AUG TIC. Not used for calculations (but still rated).
            genewin = {}; gx = 0.0; gtot = len(g2t); missx = 0; no_atg = []; atgx = 0
            self.dict['AIC'] = {}       # Dictionary of non-ATG codons and their AccNum
            seqlist.seq = []
            no_cds = 0; errx = 0
            for gene in rje.sortKeys(g2t):
                self.progLog('\r#GENE','Processing genes: %.2f%%' % (gx/gtot)); gx += 100.0
                for enst in g2t[gene]:
                    ## Get all details for sequence into one sequence object ##
                    seq = seqcdna.pop(enst)
                    seq.info['cDNA'] = seq.info['Sequence']
                    seq.info['EnsG'] = seq.info['Gene']
                    try:
                        cds = seqcds.pop(enst)
                        if not cds.info['Sequence']: raise ValueError
                        seq.info['Gene'] = cds.info['Gene']
                        seq.info['CDS'] = cds.info['Sequence']
                        if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']) and seq.info['Gene'] not in self.list['Track']: self.list['Track'].append(seq.info['Gene'])
                    except: no_cds += 1; continue    # Not a CDS
                    if self.info['NRGene'].lower() == 'ensg': nrgene = seq.info['EnsG']
                    else: nrgene = seq.info['Gene']
                    if nrgene not in genewin: genewin[nrgene] = {}
                    try:
                        us = seq5flank.pop(enst)
                        #seq.info['Gene'] = us.info['Gene']
                        seq.info['Flank'] = us#.info['Sequence']
                    except:
                        self.printLog('\r#NOFLANK','No 5\' flank sequence found for %s' % enst,screen=False)
                        try: utr = seq5utr.pop(enst)
                        except: self.printLog('\r#NOUTR','No 5\'utr sequence found for %s' % enst,screen=False)
                        missx += 1; continue
                    try:
                        utr = seq5utr.pop(enst)
                        seq.info['5utr'] = utr  #.info['Sequence']
                        if seq.info['5utr']:
                            if seq.info['Sequence'][:len(seq.info['5utr'])] != seq.info['5utr']:
                                self.printLog('\r#ERR','%s 5\' UTR sequence does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['Sequence'][:len(seq.info['5utr'])])
                                raise ValueError
                            no3 = seq.info['5utr'] + seq.info['CDS']
                            if seq.info['Sequence'].find(no3) != 0:
                                self.printLog('\r#ERR','%s 5\' UTR + CDS does not match cDNA sequence!' % enst)
                                self.deBug(seq.info['5utr'])
                                self.deBug(seq.info['CDS'])
                                self.deBug(no3)
                                self.deBug(seq.info['Sequence'])
                                raise ValueError
                        else: self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    except ValueError: errx += 1; continue
                    except:
                        seq.info['5utr'] = ''
                        self.printLog('\r#NOUTR','No 5\' UTR sequence for %s' % enst,screen=False)
                    ## Assess sequence and dictate eType ##
                    seq.info['utr5'] = seq.info['5utr']
                    if seq.info['5utr']:
                        etype = 'A'
                        seq.stat['e5utr'] = len(seq.info['utr5'])
                    else:
                        etype = 'U'
                        seq.stat['e5utr'] = 0
                    kozwin = seq.info['CDS'][:self.stat['KozWin']]
                    if kozwin[:3] != 'ATG':
                        self.printLog('\r#AIC','%s AIC found for %s' % (kozwin[:3],enst),screen=False)
                        etype = 'N'
                        if kozwin[:3] not in self.dict['AIC']: self.dict['AIC'][kozwin[:3]] = []
                        self.dict['AIC'][kozwin[:3]].append(enst)
                    if len(seq.info['utr5']) < self.stat['KozWin']: etype += 'S'
                    else: etype += 'L'
                    seq.info['eType'] = etype
                    # For each transcript, the following sequences are then extracted:
                    # - The "TICwin"  window around TIC (as set by kozwin). [*.kozwin.fas]
                    # - A control "Conwin" window around the first AUG encountered in seq <2>. [*.conwin.fas]
                    flank = seq.info['Flank'] + seq.info['utr5']
                    seq.info['TICWin'] = flank[-self.stat['KozWin']:] + kozwin
                    if len(seq.info['TICWin']) < (2 * self.stat['KozWin']): missx += 1; continue
                    if seq.info['TICWin'] not in genewin[nrgene]: genewin[nrgene][seq.info['TICWin']] = []
                    genewin[nrgene][seq.info['TICWin']].append(seq)
                    seq.info['Core'] = flank[-3:] + kozwin[:4]
                    seq.info['WegCore'] = flank[-7:-5] + 'nn' + flank[-3] + 'nn' + kozwin[:4] # XXnnXnnAUGX
                    ## Add ConWin ##
                    flank = seq.info['Flank'][-self.stat['Flank']:]
                    conatg = flank[self.stat['KozWin']:-self.stat['KozWin']].find('ATG')
                    if conatg >= 0:
                        seq.info['ConWin'] = flank[conatg:][:2*self.stat['KozWin']]
                        seq.info['ConCore'] = seq.info['ConWin'][self.stat['KozWin']-3:][:7]
                        seq.info['ConWeg'] = seq.info['ConWin'][self.stat['KozWin']-7:][:2] + 'nn' + seq.info['ConWin'][self.stat['KozWin']-3] + 'nn' + seq.info['ConWin'][self.stat['KozWin']:][:4] # XXnnXnnAUGX
                    else:
                        seq.info['ConWin'] = ''
                        self.printLog('\r#CON','No upstream control ATG window found for %s' % (enst),screen=False)
                    # Each transcript is then classified according to its core sequence:
                    # - Good = "best" Kozak gccaugg [set by GoodIC list]
                    # - augG = Not "Good" but has G at +4
                    # - Rxx = Not "Good" or augG but has purine at -3
                    # - Poor = All the rest
                    seq.info['CoreType'] = self._coreType(seq.info['Core'])
                    #if seq.info['Core'] in self.list['GoodIC']: seq.info['CoreType'] = 'Good'
                    #elif seq.info['Core'][-1] == 'G': seq.info['CoreType'] = 'augG'
                    #elif seq.info['Core'][0] in ['G','A'] and seq.info['Core'][-1] == 'A': seq.info['CoreType'] = 'Rxx'     #!# Modified #!#
                    #elif seq.info['Core'][0] in ['G','A'] and seq.info['Core'][-1] == 'G' and seq.info['Core'][-4:-1] in ['CTG','GTG','ACG']: seq.info['CoreType'] = 'Alt'     #!# Modified #!#
                    #else: seq.info['CoreType'] = 'Poor'
                    #if (seq.info['AccNum'] in self.list['Track'] or seq.info['EnsG'] in self.list['Track']): self.deBug(seq.details())
            self.printLog('\r#GENE','Processing %s genes complete: %s Genes with TICwin' % (rje.integerString(gtot),rje.integerString(len(genewin))))
            self.printLog('#MISS','%s cDNA transcripts excluded due to missing/short 5\' flanking sequence' % rje.integerString(missx))
            self.printLog('#NOCDS','%s cDNA transcripts excluded due to no CDS' % rje.integerString(no_cds))
            self.printLog('#ERRX','%s cDNA transcripts excluded due to sequence errors' % rje.integerString(errx))
            #x#self.printLog('#ATG','%s cDNA transcripts excluded due to missing 5\' UTR and no ATG in cDNA' % rje.integerString(len(no_atg)))
            open('%s.no_atg.txt' % self.info['Basefile'],'w').write(string.join(no_atg,'\n'))
            aicx = 0
            for aic in rje.sortKeys(self.dict['AIC']):
                self.printLog('#AIC','%s start -> %s transcripts.' % (aic,rje.integerString(len(self.dict['AIC'][aic]))),screen=False)
                aicx += len(self.dict['AIC'][aic])
            self.printLog('#AIC','Total AIC (non-ATG) transcripts = %s.' % (rje.integerString(aicx)))

            ### ~ [3] ~ Remove redundancy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # For each NR gene, transcripts are checked for redundancy and only one representative of each
            # different TICwin sequence is kept.
            gx = 0.0; gtot = len(genewin); notx = 0
            for gene in rje.sortKeys(genewin):
                self.progLog('\r#GENE','Removing redundancy genes: %.2f%%' % (gx/gtot)); gx += 100.0
                if not genewin[gene]: notx += 1; continue
                for ticwin in genewin[gene]:
                    for seq in genewin[gene][ticwin]:
                        seq.stat['Redundancy'] = len(genewin[gene][ticwin]) - 1
                        seq.stat['AIC'] = len(genewin[gene])
                    for type in ['AL','AS','US','NL','NS']:    # Keep the "best"
                        if len(genewin[gene][ticwin]) == 1: break
                        for seq in genewin[gene][ticwin]:
                            if seq.info['eType'] == type: genewin[gene][ticwin] = [seq]; break
                    seqlist.seq.append(genewin[gene][ticwin][0])     # Add to sequence object
                    if seqlist.seq[-1].info['Core'][-4:-1] == 'ATG': atgx += 1
                    if gene in self.list['Track']: self.bugPrint('\rNR: %s' % seqlist.seq[-1].info)
                if gene in self.list['Track']: self.deBug('NR for gene %s complete' % gene)
            self.printLog('\r#GENE','Processing %s genes complete: %s NR sequences; %s without sequences.' % (rje.integerString(gtot),rje.integerString(seqlist.seqNum()),rje.integerString(notx)))
            self.printLog('\r#ATG','%s of %s NR sequences have ATG start codon.' % (rje.integerString(atgx),rje.integerString(seqlist.seqNum())))

            ### ~ [4] ~ Regional GC calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Regional GC content is estimated from the end of flanking seq <2>.
            sx = 0.0; stot = seqlist.seqNum()
            for seq in seqlist.seq:
                self.progLog('\r#GC','Calculating GC content from 5\' flanking sequence: %2.f%%' % (sx/stot)); sx += 100.0
                seq.dict['GC'] = {}
                for n in seq.info['Flank'][-self.stat['Flank']:]:
                    if n not in seq.dict['GC']: seq.dict['GC'][n] = 0
                    seq.dict['GC'][n] += 1
            self.printLog('\r#GC','Calculation of GC content from 5\' flanking sequence complete.')

            ### ~ [5] ~ Output unprocessed sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FLANK','%s 5\' flanking sequence entries remain unprocessed.' % rje.integerString(len(seq5flank)))
            if seq5flank: open('%s.seq5flank.left' % self.info['Basefile'],'w').write(string.join(rje.sortKeys(seq5flank),'\n'))
            self.printLog('#UTR','%s 5\' UTR sequence entries remain unprocessed.' % rje.integerString(len(seq5utr)))
            if seq5utr: open('%s.seq5utr.left' % self.info['Basefile'],'w').write(string.join(rje.sortKeys(seq5utr),'\n'))

            ### ~ [6] ~ Save SeqList object and pickle for later time-saving ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'] = seqlist
            self.pickleMe('%s.loaded' % self.info['Basefile'])
            return seqlist

        except: self.errorLog('Jan-2014-Kozak Loading Sequence Error')
#########################################################################################################################
### End of SECTION II: AIC Class                                                                                        #
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
### SECTION IV: EMBLEntry Class                                                                                         # 
#########################################################################################################################
class GenBankEntry(rje_uniprot.UniProtEntry):     
    '''
    GenBank Entry Class. Author: Rich Edwards (2009).

    Info:str
    - Name = GenBank ID of Entry
    - Type = Molecule type
    - FullText = Full Text of Entry
    
    Opt:boolean
    - InvMask = Whether to invert the masking and only retain maskft features [False]    

    Stat:numeric
    - Length = Length of Sequence as annotated

    List:list
    - CaseFT = List of Features to make upper case with rest of sequence lower case []
    - Feature = List of feature dictionaries: [Type,Start,End,Desc]
    - FTSkip = List of feature details to skip ['transl_table','translation']
    - MaskFT = List of Features to mask out []
    - PubMed = List of PubMed IDs (as strings)
    - Keywords = List of UniProt Keywords
    - Tissues = List of UniProt Tissues
    - Synonyms = List of Gene synonyms
    
    Dict:dictionary

    Obj:RJE_Objects
    - Sequence = rje_sequence.Sequence object
    '''
    ### Attributes
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['Name','Type','FullText']
        self.statlist = ['Length']
        self.optlist = ['InvMask']
        self.listlist = ['Feature','PubMed','Keywords','Tissues','Synonyms','MaskFT','CaseFT','FTSkip']
        self.dictlist = []
        self.objlist = ['Sequence']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.obj['Sequence'] = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list)
        self.info['FullText'] = ''
        self.list['Feature'] = []   # List of features = {'Type':str,'Start':int,'End':int,'Desc':str}
        self.list['MaskFT'] = []
        #?#self.list['FTSkip'] = ftskip[0:]
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:### General Options ###
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'opt',['CC2FT','InvMask','TMConvert'])
                self._cmdReadList(cmd,'list',['CaseFT','MaskFT','FTSkip'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Attribute Processing                                                                                    #
#########################################################################################################################
    def process(self,rlines,ft_only=[]):   ### Limited processing of GenBank Entry
        '''
        Limited processing of GenBank Entry.
        >> rlines:list of lines containing Entry data
        >> ft_only:bool [] = List of features to exclusively read in.
        '''
        try:### ~ [1] ~ Setup and read through lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.deBug(string.join(rlines,'\n'))
            i = -1; atype = ''; x = 12
            while (i+1) < len(rlines):
                ## ~ [1a] ~ Next line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                i += 1
                iline = rlines[i]
                itype = string.replace(iline[:x],' ','')   # Annotation type
                if atype == 'FEATURES' and iline[:1] == ' ':
                    itype = ''
                    ftlines.append(iline)
                if itype:
                    if atype == 'FEATURES': self.addFeatures(ftlines,ft_only)
                    atype = itype
                iline = iline[x:]
                ## ~ [1b] ~ Cont detail ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not itype:
                    if atype == 'DEFINITION': self.info['Description'] = string.replace('%s %s' % (self.info['Description'],iline),'  ',' ')
                    if atype == 'ORGANISM': self.info['Taxa'] = string.replace('%s %s' % (self.info['Taxa'],iline),'  ',' ')
                    continue
                ## ~ [1c] ~ New detail ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if atype == 'LOCUS':
                    self.info['Name'] = string.split(iline)[0]
                    if rje.matchExp('(\S+)\s+(\d+) bp\s+(\S+)',iline):  # NM_001159380            3002 bp    mRNA    linear   PRI 16-APR-2009
                        (name,length,type) = rje.matchExp('(\S+)\s+(\d+) bp\s+(\S+)',iline)
                        self.stat['Length'] = string.atoi(length)
                        self.info['Type'] = type
                if atype == 'DEFINITION': self.info['Description'] = iline
                if atype == 'ACCESSION': self.info['AccNum'] = iline
                if atype == 'VERSION':
                    self.info['Version'] = string.split(iline)[0]
                    try: self.info['GI'] = rje.matchExp('GI:(\S+)',string.split(iline)[1])[0]
                    except: pass    # Not...  VERSION     NM_001159380.1  GI:226823325
                if atype == 'ORGANISM': self.info['Species'] = iline; self.info['Taxa'] = ''
                ## ~ [1d] ~ Features ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if atype == 'FEATURES': ftlines = []
            self.deBug(self.info)
            self.deBug(self.list)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
    def addFeatures(self,ftlines,ft_only=[],x=21):   ### Limited processing of features from self.info['FTText']
        '''Limited processing of features from self.list['FTText'].'''
        try:### ~ [1] ~ Setup and read through lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            i = -1; ftype = ''; ft = {}; fdet = ''
            self.deBug(string.join(ftlines,'\n'))
            while (i+1) < len(ftlines):
                ## ~ [1a] ~ Next line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                i += 1
                iline = ftlines[i]
                itype = string.replace(iline[:x],' ','')   # Annotation type
                if itype: ftype = itype
                iline = iline[x:]
                ## ~ [1b] ~ Cont feature ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not itype:
                    if ft:
                        ft['Desc'] = string.replace('%s %s' % (ft['Desc'],iline),'  ',' ')
                        if rje.matchExp('/(\S+)=(\S.+)',iline):
                            (fdet,dtxt) = rje.matchExp('/(\S+)=(\S.+)',iline)
                            ft[fdet] = string.strip(dtxt,'"')
                        else: ft[fdet] = string.strip(string.replace('%s %s' % (ft[fdet],iline),'  ',' '),'"')
                    continue
                ## ~ [1c] ~ New detail ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.deBug(itype)
                self.deBug(iline)
                if ft: self.list['Feature'].append(ft); self.deBug(ft)
                if ft_only and ftype.lower() not in ft_only: ft = {}; continue
                ft = {'Type':ftype,'Desc':''}
                if rje.matchExp('join\((\d+).+(\d+)\)',iline): (start,end) = rje.matchExp('join\((\d+).+(\d+)\)',iline)
                else:
                    try: [start,end] = string.split(string.replace(string.replace(iline,'<',''),'>',''),'..')
                    except: start = end = string.split(iline,'..')[0]
                try:
                    ft['Start'] = string.atoi(start)
                    ft['End'] = string.atoi(end)
                except: self.errorLog('Problem with FT: %s' % ftlines[i]); ft = {}
            if ft: self.list['Feature'].append(ft); self.deBug(ft)
        except: self.errorLog(rje_zen.Zen().wisdom())
#########################################################################################################################
## End of SECTION III : GenBankEntry Class                                                                              #
#########################################################################################################################



#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### ~ [1] ~ Basic Setup of Program  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: (info,out,mainlog,cmd_list) = setupProgram()
    except SystemExit: return  
    except: print 'Unexpected error during program setup:', sys.exc_info()[0]; return
    
    ### ~ [2] ~ Rest of Functionality... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    try: AIC(mainlog,cmd_list).run()
        
    ### ~ [3] ~ End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except:mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program,info.version,time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
