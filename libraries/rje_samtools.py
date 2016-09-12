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
Module:       rje_samtools
Description:  RJE SAMtools parser and processor
Version:      1.8.1
Last Edit:    29/08/16
Copyright (C) 2013  Richard J. Edwards - See source code for GNU License Notice

Function:
    The initial function of this program is for calling/assessing genetic changes following MPileup mapping of multiple
    short read datasets against the same reference genome (provided as BAM files). The MPileup files should be generated
    by piping the output of the following into a file (*.mpileup):

    samtools mpileup -BQ0 -d10000000 -f <Ref Genome Fasta> <BAM file>

    Initial parsing and filtering converts the MPileup format into a delimited text file with quality score and read
    depth filtering, converting mapped read data into an allele list in the form: allele1:count|allele2:count|... These
    will be sorted by frequency, so that allele1 is the "major allele". Output of the *.QX.tdt file will have the fields:
    * Locus = Reference locus (contig/chromosome)
    * Pos = Position in reference (1-L)
    * Ref = Reference sequence
    * N = Total read depth at that position
    * QN = Read depth after quality filtering
    * Seq = Mapped (filtered) allele list in the form: allele1:count|allele2:count|...

SNP Frequency Calculations:
    A second function of this tool is to compare the SNP frequencies of two populations/datasets and identify major
    alleles in the "Treatment" that have significantly increased in frequency compared to the "Control". This mode takes
    two pileup files as control=FILE and treatment=FILE. These file names (minus extension) will be output fields. For
    something more user-friendly, use `labels=X,Y` to give them better labels (where X is control and Y is treatment).
    These file names will also be used to set the output files `CONTROL.vs.TREATMENT` unless `basefile=X` is set.

    Parsed pileup files (see above) are read in and combined. Only locus positions with entries in both files (i.e. both
    meet the `qcut=X` and `minqn=X` criteria) are kept and base calls combined. Any alleles failing to meet the minimum
    count criteria (mincut=X) are removed. If `mincut` < 1.0, it is treated as a proportion of all reads, with a minimum
    value set by `absmincut=X`. The exception to this allele removal is the Reference allele, unless `ignoreref=T`. By
    default, `N` alleles are removed and do not contribute to overall read counts for allele frequency calculations. This
    can be changed by setting `ignoren=F`. Filtered SNPs are output to `*.snp.tdt`. Unless `ignoreref=T`, fixed
    differences to the Reference (i.e. 100% frequency in both control and treatment) are also output.

    The focus of analysis is the "major allele" for the treatment population, which is the allele with highest frequency.
    (When tied, alleles will be arbitrarily selected in alphabetical order). The frequency of  the major allele in the
    treatment (`MajFreq`) is output along with the difference in frequency of this allele relative to the control
    (`MajDiff`). A positive value indicates that the allele is at higher frequency in the treatment, whereas a negative
    value means higher frequency in the control. Finally, the probability of the treatment frequency is calculated using
    a binomial distribution, where: p is probability of a read being the major allele based on the control frequency (if
    zero, this is conservatively set to 1/(N+1), i.e. assuming that the next control read would be that allele); k is the
    number of major allele treatment reads and n is the total number of treatment reads.

    The next stage is to filter these SNPs and calculate FDR statistics. By default, only positions with non-reference
    major treatment alleles are considered. This can be switched off with `majmut=F` or reversed with `majmut=F` and
    `majref=T` in combination. Following this filtering, the total
    number of SNPs is used for an FDR correction of the `MajProb` probabilities, using the Benjamini & Hochberg approach,
    where FDR = pN/n. To speed this up and reduce output, a probability cutoff can be applied with `sigcut=X` (0.05 by
    default). To only keep positions where the treatment major allele is different to the control major allele, use
    `majdif=T`. `majcut=X` will add an additional min. frequency criterion. These are applied after FDR correction.
    Output is then further reduced to those SNPs where `MajDiff` > 0 and filtered with a final FDR cutoff (`fdrcut=X`).
    Remaining SNPs are output to the `*.fdr.tdt` table.

    Optionally, a table of SNP annotation (snptable=FILE) can be merged with the `*.fdr.tdt` output. This table must have
    `Locus` and `Pos` fields. If it has other fields that are required to identify unique entries, these should be set
    with `snptabkeys=LIST`.

Read Coverage Analysis:
    Version 1.8.0 added read coverage analysis, which will calculate average depth of coverage for each Locus based on a
    SAM or mpileup file, or a Read ID file (`*.rid.tdt`) previously generated by `rje_samtools`. This file is given by
    `readcheck=FILE`. If the source sequence file is also given using `seqin=FASFILE` then the true sequence lengths are
    used for the calculation. Otherwise, the last position covered by a read in the `readcheck` file is used. Read depths
    are output to `*.coverage.tdt`.

    A table of regions to be checked can be provided using `checkpos=FILE`. This should have `Locus`, `Start` and `End`
    fields for the regions to be checked. Reads from `readcheck` will be scanned to (a) identify reads completely
    spanning the region, and (b) calculate the mean depth of coverage for that region. Reads spanning the region must
    also cover a number of nucleotides flanking the region, as set by `checkflanks=LIST`. Each flanking distance in the
    list will have its own `SpanX` output field. For example, the default settings check 0 flanks (just the region), plus
    flanks of 100, 500 and 1000 nt. This generates output fields `Span0`, `Span100`, `Span500` and `Span1000`.

    Alternative `Locus`, `Start` and `End` fields for `checkpos=FILE` can be given with `checkfields=LIST`. NOTE: This is
    case-sensitive and needs all three elements, even if only one or two fields are being changed.

Commandline:
    ### ~ MPileup Parsing Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    batch=FILELIST  : List of MPileup/SAM files to be parsed and filtered (e.g. *.pileup) []
    qcut=X          : Min. quality score for a call to include [30]
    minqn=X         : Min. number of reads meeting qcut (QN) for output [10]
    rid=T/F         : Whether to include Read ID (number) lists for each allele [True]
    snponly=T/F     : Whether to restrict parsing output to SNP positions (will use mincut settings below) [False]
    indels=T/F      : Whether to include indels in "SNP" parsing [True]
    ### ~ SNP Frequency Options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    control=FILE    : MPileup or processed TDT file to be used for control SNP frequencies []
    treatment=FILE  : MPileup or processed TDT file to be used for treatment SNP frequencies []
    labels=X,Y      : Optional labels for Control and Treatment fields in output (other file basename used) []
    mincut=X        : Minimum read count for minor allele (proportion if <1) [3]
    absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    biallelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    ignoren=T/F     : Whether to exclude "N" calls from alleles [True]
    ignoreref=T/F   : If False will always keep Reference allele and call fixed change as SNP [False]
    basefile=X      : Basename for frequency comparison output [<CONTROLBASE>.v.<TREATMENTBASE>]
    majfocus=T/F    : Whether the focus is on Major Alleles (True) or Mutant/Reference Alleles (False) [True]
    majdif=T/F      : Whether to restrict output and stats to positions with Major Allele differences in sample [False]
    majmut=T/F      : Whether to restrict output and stats to positions with non-reference Major Allele [True]
    majref=T/F      : Whether to restrict output and stats to positions with reference Major Allele (if majmut=F) [False]
    majcut=X        : Frequency cutoff for Major allele [0.0]
    sigcut=X        : Significance cutoff for enriched treatment SNPs [0.05]
    fdrcut=X        : Additional FDR cutoff for enriched treatment SNPs [1.0]
    snptable=FILE   : Table of SNPs of cross-reference with FDR SNP output []
    snptabkeys=LIST : Fields that make unique key entries for snptable (with Locus, Pos) []
    snptabmap=X,Y   : Optional SNPTable fields to replace for mapping onto FDR Locus,Pos fields (Locus,Pos) []
    ### ~ Double Genome Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    altcontrol=FILE     : MPileup or processed TDT file to be used for control SNP frequencies in Alt genome []
    alttreatment=FILE   : MPileup or processed TDT file to be used for treatment SNP frequencies in Alt genome []
    ### ~ Read Coverage Analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    readcheck=FILE      : SAM/Pileup/RID file with read mappings [None]
    seqin=FASFILE       : Sequence file for loci in MPileup/SAM files (e.g. matching all relevant RID files) [None]
    checkpos=TDTFILE    : File of Locus, Start, End positions for read coverage checking [None]
    checkfields=LIST    : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    checkflanks=LIST    : List of lengths flanking check regions that must also be spanned by reads [0,100,500,1000]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_db, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import glob, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_obj, rje_seqlist, rje_sequence, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1.0 - Modified version to handle multiple loci per file. (Original was for single bacterial chromosomes.)
    # 0.2.0 - Added majmut=T/F : Whether to restrict output and stats to positions with non-reference Major Allele [False]
    # 1.0.0 - Major reworking. Old version frozen as rje_samtools_V0.
    # 1.1.0 - Added snptabmap=X,Y alternative SNPTable mapping and read_depth statistics []. Added majref=T/F.
    # 1.2.0 - Added developmental combining of read mapping onto two different genomes.
    # 1.3.0 - Major debugging and code clean up.
    # 1.4.0 - Added parsing of read number (to link SNPs) and fixed deletion error at same time. Added rid=T/F and snponly=T/F.
    # 1.5.0 - Added biallelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    # 1.5.1 - Fixed REF/Ref ALT/Alt bug.
    # 1.6.0 - Added majfocus=T/F : Whether the focus is on Major Alleles (True) or Mutant/Reference Alleles (False) [True]
    # 1.7.0 - Added parsing of *.sam files for generating RID table.
    # 1.8.0 - Added read coverage summary/checks.
    # 1.8.1 - Fixed issue when RID file not generated by pileup parsing. Set RID=True by default to avoid issues.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [Y] : Populate Module Docstring with basic info.
    # [Y] : Populate makeInfo() method with basic info.
    # [Y] : Add full description of program to module docstring.
    # [Y] : Create initial working version of program.
    # [ ] : Add option to limit output to different Major Alleles only. (And/or Allele Freq cutoffs?)
    # [ ] : Add option to output data for specific subset of positions (e.g. from different analysis)
    # [ ] : Add option to compare outputs with different QCs?
    # [ ] : Add options/warnings for low QN counts?
    # [ ] : Add Locus to SNP table.
    # [ ] : Add Locus to mpileup reading.
    # [ ] : Add pvalue threshold and p-value cutoff summaries.
    # [ ] : Add a minimum count for WT and Mut for output. (Stop reporting errors and warnings unless flag set.)
    # [ ] : Add output of major allele frequency distributions? (Or just make R code and output PNG?)
    # [ ] : Check formatting of numbers in tables.
    # [ ] : Calculate SNP densities?
    # [ ] : Add feature/region-focused reporting. Input table of regions/tables for XRef of data instead of SNP Table.
    # [Y] : Add indels=T/F option to only retain SNPs and not indels.
    # [ ] : Need to deal with insertions properly as will all be mapped onto same reference position!
    # [ ] : Consider whether to put the QX.N.tdt and rid.tdt files in the pileup/sam file directory.
    # [ ] : Tidy up error handling and messages with rid.tdt file missing.
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copyyear) = ('rje_samtools', '1.8.1', 'August 2016', '2013')
    description = 'RJE SAMtools parser and processor'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is still in development and has not been published.',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
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
### SECTION II: SAMtools Class                                                                                          #
#########################################################################################################################
class SAMtools(rje_obj.RJE_Object):     
    '''
    Class. Author: Rich Edwards (2012).

    Str:str
    - AltControl=FILE    : MPileup or processed TDT file to be used for control SNP frequencies in Alt genome []
    - AltTreatment=FILE  : MPileup or processed TDT file to be used for treatment SNP frequencies in Alt genome []
    - Basefile=X      : Basename for frequency comparison output [<CONTROLBASE>.v.<TREATMENTBASE>]
    - CheckPos=TDTFILE    : File of Locus, Start, End positions for read coverage checking [None]
    - Control=FILE    : MPileup file to be used for control SNP frequencies []
    - ReadCheck=FILE      : SAM/Pileup/RID file with read mappings [None]
    - Seqin=FASFILE       : Sequence file for loci in MPileup/SAM files (e.g. matching all relevant RID files) [None]
    - SNPTable        : File with existing SNPs to map
    - Treatment=FILE  : MPileup file to be used for treatment SNP frequencies []

    Bool:boolean
    - Biallelic=T/F   : Whether to restrict SNPs to pure biallelic SNPs (two alleles meeting mincut) [False]
    - IgnoreN = Whether to exclude "N" calls for major/minor alleles [True]
    - IgnoreRef=T/F : If False will always keep Reference allele and call fixed change as SNP [False]
    - Indels=T/F      : Whether to include indels in "SNP" parsing [True]
    - MajDif = Whether to restrict output and stats to positions with Major Allele differences [True]
    - MajFocus=T/F    : Whether the focus is on Major Alleles (True) or Mutant/Reference Alleles (False) [True]
    - MajMut = Whether to restrict output and stats to positions with non-reference Major Allele [True]
    - MajRef = Whether to restrict output and stats to positions with non-reference Major Allele [False]
    - RID=T/F         : Whether to include Read ID (number) lists for each allele [False]
    - SNPOnly=T/F     : Whether to restrict parsing output to SNP positions (will use mincut settings below) [False]

    Int:integer
    - AbsMinCut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    - FDRCut=X        : Additional FDR cutoff for enriched treatment SNPs [1.0]
    - MajCut=X        : Frequency cutoff for Treatment major allele [0.0]
    - MinCut=X        : Minimum read count for minor allele (proportion if <1) [2]
    - MinQN = Min. number of reads meeting qcut (QN) for output [10]
    - QCut = Min. quality score for a call to include [30]
    - SigCut=X        : Significance cutoff for enriched treatment SNPs [0.05]

    Num:float
    - MinFreq = Minor allele(s) frequency correction for zero counts (e.g. Sequencing error) [0.01]
    
    List:list
    - Batch=FILELIST  : List of MPileup files to be parsed and filtered []
    - CheckFields=LIST    : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    - CheckFlanks=LIST       : List of lengths flanking check regions that must also be spanned by reads [0,100,500,1000]
    - Labels=X,Y      : Optional labels for Control and Treatment fields in output (other file basename used) []
    - SNPTabKeys=LIST : Fields that make unique key entries for snptable (with Locus, Pos) []
    - SNPTabMap=X,Y   : Optional SNPTable fields to replace for mapping onto FDR Locus,Pos fields (Locus,Pos) []

    Dict:dictionary    

    Obj:RJE_Objects
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### ~ Basics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.strlist = ['AltControl','AltTreatment','CheckPos','Control','ReadCheck','SNPTable','Treatment']
        self.boollist = ['Biallelic','IgnoreN','IgnoreRef','Indels','MajDif','MajFocus','MajMut','MajRef','RID','SNPOnly']
        self.intlist = ['AbsMinCut','CheckFlanks','MinQN','QCut']
        self.numlist = ['FDRCut','MinCut','MinFreq','SigCut']
        self.listlist = ['Batch','CheckFields','CheckFlanks','Labels','SNPTabKeys','SNPTabMap']
        self.dictlist = []
        self.objlist = []
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'RefSeq':''})
        self.setBool({'Biallelic':False,'IgnoreN':True,'IgnoreRef':False,'Indels':True,'MajDif':False,'MajFocus':True,'MajMut':True,'MajRef':False,'RID':True})
        self.setInt({'QCut':30,'MinQN':10,'AbsMinCut':2})
        self.setNum({'MajCut':0.0,'MinFreq':0.001,'MinCut':2.0,'SigCut':0.05,'FDRCut':1.0})
        self.list['Batch'] = [] #x# glob.glob('*.pileup')
        self.list['CheckFlanks'] = [0,100,500,1000]
        self.list['CheckFields'] = ['Locus','Start','End']
        self.list['SNPTabMap'] = ['Locus','Pos']
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
                self._cmdReadList(cmd,'file',['AltControl','AltTreatment','CheckPos','Control','ReadCheck','SNPTable','Treatment'])  # String representing file path
                self._cmdReadList(cmd,'bool',['Biallelic','IgnoreN','IgnoreRef','Indels','MajDif','MajFocus','MajMut','MajRef','RID','SNPOnly'])  # True/False Booleans
                self._cmdReadList(cmd,'int',['AbsMinCut','QCut','MinQN'])   # Integers
                self._cmdReadList(cmd,'float',['MinCut','FDRCut','MajCut','MinFreq','SigCut']) # Floats
                #self._cmdReadList(cmd,'min',['Att'])   # Integer value part of min,max command
                #self._cmdReadList(cmd,'max',['Att'])   # Integer value part of min,max command
                self._cmdReadList(cmd,'list',['CheckFields','Labels','SNPTabKeys','SNPTabMap'])  # List of strings (split on commas or file lines)
                self._cmdReadList(cmd,'ilist',['CheckFlanks']) # Comma separated list as integers
                self._cmdReadList(cmd,'glist',['Batch']) # List of files using wildcards and glob
                #self._cmdReadList(cmd,'cdict',['Att']) # Splits comma separated X:Y pairs into dictionary
                #self._cmdReadList(cmd,'cdictlist',['Att']) # As cdict but also enters keys into list
            except: self.errorLog('Problem with cmd:%s' % cmd)
        if self.getBool('MajMut') and self.getBool('MajRef'):
            if self.i() < 0 or rje.yesNo('Cannot have both MajMut=T and MajRef=T. Switch MajRef=F?'):
                self.setBool({'MajRef':False})
                self.printLog('#MAJREF','Cannot have both MajMut=T and MajRef=T: MajRef=F')
            else:
                self.setBool({'MajMut':False})
                self.printLog('#MAJMUT','Cannot have both MajMut=T and MajRef=T: MajMut=F')
        if self.getBool('MajFocus') and not (self.getBool('MajMut') or self.getBool('MajRef')):
            self.warnLog('Cannot have majfocus=T without majmut=T or majref=T.')
            self.setBool({'MajFocus':False})
            #!# Add option to change settings
#########################################################################################################################
    ### <2> ### Main Class Backbone                                                                                     #
#########################################################################################################################
    def run(self):  ### Main run method
        '''Main run method.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setup()    # This includes basic parsing of pileup files.
            ### ~ [2] ~ Add main run code here ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Control') and self.getStrLC('Treatment'):
                if self.getStrLC('AltControl') and self.getStrLC('AltTreatment'):
                    snpfile = self.combinePileUpStats()
                else: snpfile = self.pileUpStats()
                fdrdb = self.rateSNPs(snpfile)
                self.combineSNPs(fdrdb)
            if self.getStrLC('ReadCheck'):
                ridfile = self.getStr('ReadCheck')
                ridbase = rje.baseFile(ridfile,strip_path=True)
                if not self.baseFile(return_none=None): self.baseFile(ridbase)
                if not ridfile.endswith('.rid.tdt'):
                    self.setBool({'RID':True})
                    self.parsePileup(ridfile)
                    ridfile = '%s.rid.tdt' % (ridbase)  #?# Should this strip path? Probably not!
                self.coverageFromRID(ridfile)
        except:
            self.errorLog(rje_zen.Zen().wisdom())
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def setup(self):    ### Main class setup method.
        '''Main class setup method.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            if self.list['Labels'] and len(self.list['Labels']) != 2:
                self.warnLog('Labels not given exactly two values (Control,Treatment): cannot use',quitchoice=True)
                self.list['Labels'] = []
            addlabels = not self.list['Labels']
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStrLC(fstr): self.list['Batch'] = []
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStrLC(fstr):
                    fbase = rje.baseFile(self.getStr(fstr),strip_path=True)
                    if addlabels: self.list['Labels'].append(fbase)
                    fcuts = rje.matchExp('\.Q(\d+).(\d+)$',fbase)
                    if fcuts:
                        if int(fcuts[0]) != self.getInt('QCut'): self.warnLog('%s QCut mismatch (%s vs qcut=%d)' % (fstr,fcuts[0],self.getInt('QCut')))
                        if int(fcuts[1]) != self.getInt('MinQN'): self.warnLog('%s QCut mismatch (%s vs minqn=%d)' % (fstr,fcuts[0],self.getInt('MinQN')))
                    else:
                        self.list['Batch'].append(self.getStr(fstr))
                        self.printLog('#FILES','%s file %s added for pileup parsing' % (fstr,self.getStr(fstr)))
            if not self.getStrLC('Basefile') and self.getStrLC('Control') and self.getStrLC('Treatment'):
                self.baseFile('%s.vs.%s' % (rje.baseFile(self.getStr('Control'),strip_path=True),rje.baseFile(self.getStr('Treatment'),strip_path=True)))
            ### ~ [2] Parse Pileup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#FILES','%d pileup files identified for parsing' % len(self.list['Batch']))
            for pfile in self.list['Batch']: self.parsePileup(pfile)
            return True
        except: self.errorLog('Problem during %s setup.' % self); return False  # Setup failed
#########################################################################################################################
    ### <3> ### Pileup parsing Methods                                                                                  #
#########################################################################################################################
    def parseSAM(self,filename):  ### Extracts read data from SAM file
        '''
        Extracts read data from SAM file.
        >> filename:str = Pileup file name
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ridfile = '%s.rid.tdt' % (rje.baseFile(filename,strip_path=True))
            if rje.exists(ridfile) and not self.force():
                self.printLog('#SKIP','%s found! (force=F)' % ridfile)
                return True
            self.printLog('#~~#','## ~~~~~ Parsing SAM File: %s ~~~~~ ##' % filename)
            RIDOUT = open(ridfile,'w')
            rje.writeDelimit(RIDOUT,outlist=['RID','Locus','Start','End'],delimit='\t')
            SAM = open(filename,'r')
            rid = 0          # Read counter (ID counter)
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in SAM:
                ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if line.startswith('@'): continue
                samdata = string.split(line)
                if len(samdata) < 11: continue
                self.progLog('\r#PARSE','Parsing %s: %s reads...' % (filename,rje.iStr(rid)),rand=0.01)
                rid += 1
                locus = samdata[2]
                rpos = int(samdata[3])
                cigstr = samdata[5]
                cigdata = parseCigar(cigstr)
                rje.writeDelimit(RIDOUT,outlist=[rid,locus,rpos,rpos+cigarAlnLen(cigdata)-1],delimit='\t')
            self.printLog('#RID','Parsed %s: %s read start/end positions output to %s' % (filename,rje.iStr(rid),ridfile))
            RIDOUT.close()
            return True
        except: self.errorLog('%s.parseSAM() error' % (self.prog())); return False
#########################################################################################################################
    def parsePileup(self,filename,depth=True):  ### Extracts, filters and processes PileUp data
        '''
        Extracts, filters and processes PileUp data.
        >> filename:str = Pileup file name
        >> depth:bool [True] = Whether to generate read depth calculation.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if filename.endswith('.sam'): return self.parseSAM(filename)
            outfile = '%s.Q%d.%d.tdt' % (rje.baseFile(filename,strip_path=True),self.getInt('QCut'),self.getInt('MinQN'))
            ridfile = '%s.rid.tdt' % (rje.baseFile(filename,strip_path=True))
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStr(fstr) == filename: self.setStr({fstr:outfile})
            skiprun = not self.force()
            if self.getBool('RID'):
                skiprun = skiprun and rje.exists(ridfile)
            skiprun = skiprun and rje.exists(outfile)
            if skiprun:
                #!# Should we add something to check integrity/complete? #!#
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return True
            self.printLog('#~~#','## ~~~~~ Parsing PileUp File: %s ~~~~~ ##' % filename)
            ## ~ [1a] Calculate mean read depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Currently this is based on unfiltered 'N' values.
            PILEUP = open(filename,'r'); px = 0; rx = 0
            for line in PILEUP:
                data = string.split(rje.chomp(line))
                if not data: break
                self.progLog('\r#DEPTH','Parsing %s depth: %s pos; %s read bases' % (filename,rje.iStr(px),rje.iStr(rx)),rand=0.01); px += 1
                rx += int(data[3])
            meandepth = float(rx) / px
            self.printLog('\r#DEPTH','Parsed %s depth: %s pos; %s read bases -> mean = %.1f' % (filename,rje.iStr(px),rje.iStr(rx),meandepth))
            PILEUP.close()
            ## ~ [1b] Setup files for main processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,outfile)
            outfields = ['Locus','Pos','Ref','N','QN','Seq']    # Fields for output file
            if depth: outfields.append('Dep')
            RIDOUT = None
            if self.getBool('RID'):
                outfields.append('RID')
                RIDOUT = open(ridfile,'w')
                rje.writeDelimit(RIDOUT,outlist=['RID','Locus','Start','End'],delimit='\t')
            PILEUP = open(filename,'r'); px = 0; ex = 0
            PILEOUT = open(outfile,'w')
            rje.writeDelimit(PILEOUT,outlist=outfields,delimit='\t')
            qc = [0]        # List of counts at different quality scores
            ri = 0          # Read counter (ID counter)
            ridlist = []    # List of read IDs
            rdel = {}       # Dictionary of {rid:current deletion}
            rstart = {}     # Dictionary of {rid:start point}
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for line in PILEUP:
                ## ~ [2a] Parse pileup data into dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = self.parsePileupLine(line,ri,ridlist,rdel,qc)
                if not entry: break
                if ridlist: ri = max(ri,max(ridlist))
                else: self.debug(entry)
                self.progLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),rand=0.01); px += 1
                entry['Dep'] = rje.sf(entry['N']/meandepth)
                for rid in entry.pop('Start'): rstart[rid] = entry['Pos']
                for rid in entry.pop('End'):
                    startpos = rstart.pop(rid)
                    if self.getBool('RID'):
                        rje.writeDelimit(RIDOUT,outlist=[rid,entry['Locus'],startpos,entry['Pos']],delimit='\t')
                ## ~ [2b] Filter low quality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                # Remove (from back) any reads than do not meet QV cutoff
                entry['QN'] = len(entry['Reads'])
                if entry['QN'] < self.getInt('MinQN'): continue     # Not enough data for output
                ## ~ [2c] Alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if not self.pileupAlleles(entry): continue
                #table.addEntry(entry)
                outlist = []
                for field in outfields: outlist.append(entry[field])
                rje.writeDelimit(PILEOUT,outlist,delimit='\t'); ex += 1
            self.printLog('\r#PARSE','Parsed %s: %s entries from %s lines: %s reads.' % (filename,rje.iStr(ex),rje.iStr(px),rje.iStr(ri)))
            PILEOUT.close()
            PILEUP.close()
            if self.getBool('RID'):
                self.printLog('#RID','%s read start/end positions output to %s' % (rje.iStr(ri),ridfile))
                RIDOUT.close()
            ### ~ [3] Save QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qfile = '%s.QC.tdt' % rje.baseFile(filename,strip_path=True)
            if rje.exists(qfile) and not self.force(): self.printLog('#QC','%s already exists (force=F).' % qfile)
            else:
                QC = open(qfile,'w')
                QC.write('Qual\tCount\n')
                for q in range(len(qc)):
                    try: QC.write('%d\t%d\n' % (q+1,qc[q]))
                    except: self.errorLog('!')
                QC.close()
                self.printLog('#QC','Quality score counts output to: %s' % qfile)
            return True
        except: self.errorLog('%s.parsePileup() error' % (self.prog())); return False
#########################################################################################################################
    def pileupAlleles(self,entry):  ### Parses allele data from pileup data dictionary
        '''
        Parses allele data from pileup data dictionary.
        >> entry:dict = dictionary from self.parsePileupLine() [** Edited in place **]
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qread = entry['Reads']  # Dictionary of {rid:allele} passing qscore
            alleles = {}            # Dictionary of {nt:count}
            allrid = {}             # Dictionary of {nt:[ridlist]}
            if self.getBool('SNPOnly'):         # Parsing allele mincut values for comparison
                acut = self.getNum('MinCut')    # Minimum allele count allowed
                if acut < 1: acut = max(self.getInt('AbsMinCut'),len(qread)*acut)
            else: acut = 0
            ### ~ [1] Parse Alleles (and corresponding Read IDs) from entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for rid in rje.sortKeys(qread):
                read = qread[rid]
                if read in alleles: alleles[read] += 1; allrid[read].append(str(rid))
                else: alleles[read] = 1; allrid[read] = [str(rid)]
            asort = []  # List of tuples for sorting
            for allele in alleles: asort.append((alleles[allele],allele))
            asort.sort()
            asort.reverse()     # This should now be in order of allele frequency
            acheckx = 0         # Count check
            aseq = []           # List of (filtered?) alleles and counts
            akept = []          # List of kept alleles
            for i in range(len(asort)):
                acheckx += asort[i][0]
                if asort[i][1] == 'N' and self.getBool('SNPOnly') and self.getBool('IgnoreN'): continue
                if (asort[i][1] == '-' or len(asort[i][1]) > 1) and self.getBool('SNPOnly') and not self.getBool('Indels'): continue
                if asort[i][0] >= acut:
                    aseq.append('%s:%d' % (asort[i][1],asort[i][0]))
                    akept.append(asort[i][1])
            if not aseq: return False   # No alleles kept
            if self.getBool('SNPOnly'):
                if len(aseq) == 1: return False         # Not a SNP
                elif self.getBool('Biallelic') and len(aseq) > 2: return False     # Not a biallelic SNP
            entry['Seq'] = string.join(aseq,'|')
            entry['RID'] = []
            akept.sort()
            for allele in akept:    #rje.sortKeys(allrid):
                entry['RID'].append('%s:%s' % (allele,string.join(allrid[allele],',')))
            entry['RID'] = string.join(entry['RID'],'|')
            if acheckx != entry['QN']:  #!# Convert these to test statements?
                self.errorLog('Allele versus Quality count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
            return True
        except: self.errorLog('%s.pileupAlleles(%s) error' % (self.prog(),entry)); return False
#########################################################################################################################
    def parsePileupLine(self,line,ri=-1,ridlist=[],rdel={},qc=[]):  ### Parses pileUp data from single line
        '''
        Parses pileUp data from single line.
        >> line:str = Line from pileup file.
        >> ri:int [-1] = Read counter (ID counter). If -1 will not use/track read IDs.
        >> ridlist:list [] = List of current read IDs from previous line. (**Updated in place**)
        >> rdel:dict {} = Dictionary of {rid:current deletion} for sequence checking. (**Updated in place**)
        >> qc:list [] = Optional QC score value counter. (**Updated in place**)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            data = string.split(rje.chomp(line))
            if not data: return False
            ### ~ [1] Extract Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            entry = {'Locus':data[0],'Pos':int(data[1]),'Ref':data[2].upper(),'N':int(data[3]),'QN':0}
            if ri > -1: entry['Start'] = []; entry['End'] = []
            rseq = data[4]  # Read sequences from pileup
            reads = []      # List of read alleles, extracted from rseq
            # Note that ends of reads are [element]$
            # Beginning of reads are ^[QV][element]
            rid = -1        # Current Read ID
            ridend = []     # List of reads ending in this row to remove after processing. (Won't be in next row.)
            while rseq:
                try:
                    # Establish current read ID (Needed for deletions)
                    if rseq[:1] not in '^-+$' and ridlist: rid = ridlist[len(reads)]
                    # Update read data
                    if rseq[:1] in ['.',',']:   # Matches reference (+/- strand)
                        reads.append(entry['Ref']); rseq = rseq[1:]   # Matches reference
                    elif rseq[:1] == '*':   # Check for existing deletions
                        if ri > -1:
                            if rdel[rid][:1] != entry['Ref']: self.warnLog('Deletion sequence mismatch @ %s:%s (%s vs %s)' % (data[0],data[1],rdel[rid][:1],entry['Ref']),warntype='del_mismatch',quitchoice=True,suppress=True,dev=True)
                            rdel[rid] = rdel[rid][1:]
                            if not rdel[rid]: rdel.pop(rid) # End of deletion
                        reads.append('-'); rseq = rseq[1:]
                    elif rseq[:1] == '^':
                        rseq = rseq[2:]   # Indicates a new read
                        if ri > -1:
                            ri += 1; ridlist.insert(len(reads),ri)
                            entry['Start'].append(ri)
                            #self.bugPrint('Read %d start @ %s:%s' % (ri,data[0],data[1]))
                    elif rseq[:1] in ['-','+']:
                        ilen = string.atoi(rje.matchExp('^(\d+)',rseq[1:])[0])
                        indel = rseq[len('%s' % ilen)+1:][:ilen]    # Just the indel sequences
                        if rseq[:1] == '-':     # Deletion
                            if ri > -1:     # Otherwise will just trust the * characters in later lines.
                                rid = ridlist[len(reads)-1]
                                if rid in rdel: self.warnLog('Conflicting deletions for RID %s @ %s:%s' % (rid,data[0],data[1]))
                                rdel[rid] = indel.upper()
                                #self.debug(rdel)
                        else:   # Insertion
                            reads[-1] += indel.upper()  # Insertion just has whole sequence AFTER the position
                        rseq = rseq[len('%s' % ilen)+ilen+1:]
                    elif rseq[:1] in ['$']:
                        rseq = rseq[1:]     # Indicated end of a read
                        if ri > -1:
                            ridend.append(rid)
                            entry['End'].append(rid)
                            #self.deBug('Read %d end @ %s:%s' % (rid,data[0],data[1]))
                    else:
                        if rseq[0].upper() not in 'ATGCN': print ' ???', rseq[0].upper(), '???'
                        reads.append(rseq[0].upper()); rseq = rseq[1:]
                except:
                    self.bugPrint(reads)
                    self.bugPrint(ridlist)
                    self.errorLog('!')
                    self.deBug(rseq)
                    raise ValueError
            self.bugPrint('%s:%s : %s'% (entry['Locus'],entry['Pos'],ridlist))
            ## ~ [1a] Check data integrity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if ri > -1 and len(reads) != len(ridlist):
                raise ValueError('Read allele versus Read ID count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']))
            if len(reads) != entry['N']:    # Formerly erroenous: (entry['N'] + delx):
                self.deBug('%s = %d' % (data[4],entry['N']))
                self.deBug('%s = %d' % (reads,len(reads)))
                self.errorLog('Read versus Read Count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                raise ValueError
            ### ~ [2] Assess Quality Scores and generate alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qual = data[5]
            if len(reads) != len(qual):
                self.deBug('%s = %d' % (reads,len(reads)))
                self.deBug('%s = %d' % (qual,len(qual)))
                self.deBug(data)
                self.errorLog('Read versus Quality length mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                raise ValueError
            qread = {}    # Dictionary of {rid:allele} passing qscore
            for r in range(len(qual)):
                qrid = ridlist[r]
                q = ord(qual[r]) - 33
                if qc:
                    qc += [0] * (q - len(qc)); qc[q-1] += 1
                if q >= self.getInt('QCut'): qread[qrid] = reads[r]
            while ridend: ridlist.remove(ridend.pop(0))     # Remove ended reads
            entry['Reads'] = qread
            return entry
        except: self.errorLog('%s.parsePileupLine() error' % (self.prog())); return False
#########################################################################################################################
    def NEWparsePileup(self,filename,depth=True):  ### Extracts, filters and processes PileUp data
        '''
        Extracts, filters and processes PileUp data.
        >> filename:str = Pileup file name
        >> depth:bool [True] = Whether to generate read depth calculation.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.Q%d.%d.tdt' % (rje.baseFile(filename,strip_path=True),self.getInt('QCut'),self.getInt('MinQN'))
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStr(fstr) == filename: self.setStr({fstr:outfile})
            if rje.exists(outfile) and not self.force():
                #!# Should we add something to check integrity/complete? #!#
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return True
            self.printLog('#~~#','## ~~~~~ Parsing PileUp File: %s ~~~~~ ##' % filename)
            ## ~ [1a] Calculate mean read depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Currently this is based on unfiltered 'N' values.
            PILEUP = open(filename,'r'); px = 0; rx = 0
            for line in PILEUP:
                data = string.split(rje.chomp(line))
                if not data: break
                self.progLog('\r#DEPTH','Parsing %s depth: %s pos; %s read bases' % (filename,rje.iStr(px),rje.iStr(rx)),rand=0.01); px += 1
                rx += int(data[3])
            meandepth = float(rx) / px
            self.printLog('\r#DEPTH','Parsed %s depth: %s pos; %s read bases -> mean = %.1f' % (filename,rje.iStr(px),rje.iStr(rx),meandepth))
            PILEUP.close()
            ## ~ [1b] Setup files for main processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,outfile)
            outfields = ['Locus','Pos','Ref','N','QN','Seq']    # Fields for output file
            if depth: outfields.append('Dep')
            PILEUP = open(filename,'r'); px = 0; ex = 0
            PILEOUT = open(outfile,'w')
            rje.writeDelimit(PILEOUT,outlist=outfields,delimit='\t')
            qc = [] # List of counts at different quality scores
            ri = 0      # Read counter (ID counter)
            ridlist = []    # List of read IDs
            rdel = {}   # Dictionary of {rid:current deletion}
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Replace this with a method to parse a single line and return data or entry
            for line in PILEUP:
                # Split line up into data. Should be: locus, position, reference, no. reads, read data, qualscores
                data = string.split(rje.chomp(line))
                if not data: break
                self.progLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),rand=0.01); px += 1
                ## ~ [2a] Extract Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = {'Locus':data[0],'Pos':int(data[1]),'Ref':data[2].upper(),'N':int(data[3]),'QN':0}
                entry['Dep'] = rje.sf(entry['N']/meandepth)
                rseq = data[4]  # Read sequences from pileup
                reads = []      # List of read alleles, extracted from rseq
                delx = 0        # Number of deletions: these are not in the count depth (N) nor have a QV
                #!# How should/do people cope with deletions? #!#
                #!# Make it a toggle to include/exclude? #!#
                # Note that ends of reads are [element]$
                # Beginning of reads are ^[QV][element]
                rid = -1
                #while '$' in ridlist: ridlist.remove('$')   # Get rid of reads that ended in previous row.
                ridend = []     # List of reads ending in this row
                while rseq:
                    try:
                        # Establish current read ID (Needed for deletions)
                        if rseq[:1] not in '^-+$' and ridlist: rid = ridlist[len(reads)]
                        # Update read data
                        if rseq[:1] in ['.',',']:   # Matches reference (+/- strand)
                            reads.append(entry['Ref']); rseq = rseq[1:]   # Matches reference
                        elif rseq[:1] == '*':   # Check for existing deletions
                            #x#if rid in rdel:
                            if rdel[rid][:1] != entry['Ref']: self.warnLog('Deletion sequence mismatch @ %s:%s (%s vs %s)' % (data[0],data[1],rdel[rid][:1],entry['Ref']),warntype='del_mismatch',quitchoice=True,suppress=True,dev=True)
                            rdel[rid] = rdel[rid][1:]
                            if not rdel[rid]: rdel.pop(rid) # End of deletion
                            reads.append('-'); rseq = rseq[1:]
                        elif rseq[:1] == '^':
                            rseq = rseq[2:]   # Indicates a new read
                            ri += 1; ridlist.insert(len(reads),ri)
                            self.bugPrint('Read %d start @ %s:%s' % (ri,data[0],data[1]))
                        elif rseq[:1] in ['-','+']:
                            ilen = string.atoi(rje.matchExp('^(\d+)',rseq[1:])[0])
                            indel = rseq[len('%s' % ilen)+1:][:ilen]    # Just the indel sequences
                            if rseq[:1] == '-':     # Deletion
                                delx += 1           # These do not have QV scores
                                #X wrong! X# reads.append(rseq[:len('%s' % ilen)+ilen+1].upper())    # Deletion includes -n part.
                                rid = ridlist[len(reads)-1]
                                if rid in rdel: self.warnLog('Conflicting deletions for RID %s @ %s:%s' % (rid,data[0],data[1]))
                                rdel[rid] = indel.upper()
                                self.debug(rdel)
                            else:   # Insertion
                                reads[-1] += indel.upper()  # Insertion just has whole sequence AFTER the position
                            rseq = rseq[len('%s' % ilen)+ilen+1:]
                        elif rseq[:1] in ['$']:
                            rseq = rseq[1:]     # Indicated end of a read
                            #rfin = ridlist[len(reads)-1]
                            #ridlist[len(reads)-1] = '$'
                            ridend.append(rid)
                            self.deBug('Read %d end @ %s:%s' % (rid,data[0],data[1]))
                        else:
                            if rseq[0].upper() not in 'ATGCN': print ' ???', rseq[0].upper(), '???'
                            # NB. * = single bp deletion, handled below
                            reads.append(rseq[0].upper()); rseq = rseq[1:]
                    except:
                        self.bugPrint(reads)
                        self.bugPrint(ridlist)
                        self.errorLog('!')
                        self.deBug(rseq)
                        raise ValueError
                self.bugPrint('%s:%s : %s'% (entry['Locus'],entry['Pos'],ridlist))
                if len(reads) != len(ridlist):
                    raise ValueError('Read allele versus Read ID count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']))
                if len(reads) != entry['N']:    # Formerly erroenous: (entry['N'] + delx):
                    self.deBug('%s = %d' % (data[4],entry['N']))
                    self.deBug('%s = %d' % (reads,len(reads)))
                    self.errorLog('Read versus Read Count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                    raise ValueError
                #!# Need to replace reads and qual with dictionaries of {rid:value} ?
                #!# Or convert the finished rid to $ in ridlist AFTER processing of alleles #!# - THIS #!#
                ## ~ [2b] Assess Quality Scores and generate alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qual = data[5]
                if len(reads) != len(qual):
                    self.deBug('%s = %d' % (reads,len(reads)))
                    self.deBug('%s = %d' % (qual,len(qual)))
                    self.deBug(data)
                    self.errorLog('Read versus Quality length mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                    raise ValueError
                qread = {}    # Dictionary of {rid:allele} passing qscore
                for r in range(len(qual)):
                    qrid = ridlist[r]
                    q = ord(qual[r]) - 33
                    qc += [0] * (q - len(qc)); qc[q-1] += 1
                    if q >= self.getInt('QCut'): qread[qrid] = reads[r]
                while ridend: ridlist.remove(ridend.pop(0)) # Remove ended reads
                #while len(qual) < len(reads) and reads[len(qual)][0] == '-': qual.append(self.getInt('QCut'))
                #while '*' in reads: reads[reads.index('*')] = '-'   #'-1%s' % entry['Seq'].upper()
                ## ~ [2c] Filter low quality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #if entry['Pos'] in [190359]:    #100,98901,183697,169284,
                #    self.deBug(qual)
                #    self.deBug(reads)
                #    self.deBug(qc)
                # Remove (from back) any reads than do not meet QV cutoff
                entry['QN'] = len(qread)
                if entry['QN'] < self.getInt('MinQN'): continue     # Not enough data for output
                ## ~ [2d] Alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                alleles = {}    # Dictionary of {nt:count}
                for read in qread.values():
                    if read in alleles: alleles[read] += 1
                    else: alleles[read] = 1
                asort = []  # List of tuples for sorting
                for allele in alleles: asort.append((alleles[allele],allele))
                asort.sort()
                asort.reverse()     # This should now be in order of allele frequency
                acheckx = 0         # Count check
                for i in range(len(asort)):
                    acheckx += asort[i][0]
                    asort[i] = '%s:%d' % (asort[i][1],asort[i][0])
                entry['Seq'] = string.join(asort,'|')
                #if entry['Pos'] in [190359]:    #100,98901,183697,169284,
                #    self.deBug(qual)
                #    self.deBug(reads)
                #    self.deBug(alleles)
                #    self.deBug(entry)
                #    self.deBug(line)
                if acheckx != entry['QN']:  #!# Convert these to test statements?
                    self.errorLog('Allele versus Quality count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                #table.addEntry(entry)
                outlist = []
                for field in outfields: outlist.append(entry[field])
                rje.writeDelimit(PILEOUT,outlist,delimit='\t'); ex += 1
            self.printLog('\r#PARSE','Parsed %s: %s entries from %s lines: %s reads.' % (filename,rje.iStr(ex),rje.iStr(px),rje.iStr(ri)))
            PILEOUT.close()
            PILEUP.close()
            ### ~ [3] Save QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qfile = '%s.QC.tdt' % rje.baseFile(filename,strip_path=True)
            if rje.exists(qfile) and not self.force(): self.printLog('#QC','%s already exists (force=F).' % qfile)
            else:
                QC = open(qfile,'w')
                QC.write('Qual\tCount\n')
                for q in range(len(qc)):
                    try: QC.write('%d\t%d\n' % (q+1,qc[q]))
                    except: self.errorLog('!')
                QC.close()
                self.printLog('#QC','Quality score counts output to: %s' % qfile)
            return True
        except: self.errorLog('%s.parsePileupLine() error' % (self.prog())); return False
#########################################################################################################################
    def OLDparsePileup(self,filename,depth=True):  ### Extracts, filters and processes PileUp data
        '''
        Extracts, filters and processes PileUp data.
        >> filename:str = Pileup file name
        >> depth:bool [True] = Whether to generate read depth calculation.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfile = '%s.Q%d.%d.tdt' % (rje.baseFile(filename,strip_path=True),self.getInt('QCut'),self.getInt('MinQN'))
            for fstr in ['Control','Treatment','AltControl','AltTreatment']:
                if self.getStr(fstr) == filename: self.setStr({fstr:outfile})
            if rje.exists(outfile) and not self.force():
                #!# Should we add something to check integrity/complete? #!#
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return True
            self.printLog('#~~#','## ~~~~~ Parsing PileUp File: %s ~~~~~ ##' % filename)
            ## ~ [1a] Calculate mean read depth ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # Currently this is based on unfiltered 'N' values.
            PILEUP = open(filename,'r'); px = 0; rx = 0
            for line in PILEUP:
                data = string.split(rje.chomp(line))
                if not data: break
                self.progLog('\r#DEPTH','Parsing %s depth: %s pos; %s read bases' % (filename,rje.iStr(px),rje.iStr(rx)),rand=0.01); px += 1
                rx += int(data[3])
            meandepth = float(rx) / px
            self.printLog('\r#DEPTH','Parsed %s depth: %s pos; %s read bases -> mean = %.1f' % (filename,rje.iStr(px),rje.iStr(rx),meandepth))
            PILEUP.close()
            ## ~ [1b] Setup files for main processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.backup(self,outfile)
            outfields = ['Locus','Pos','Ref','N','QN','Seq']    # Fields for output file
            if depth: outfields.append('Dep')
            PILEUP = open(filename,'r'); px = 0; ex = 0
            PILEOUT = open(outfile,'w')
            rje.writeDelimit(PILEOUT,outlist=outfields,delimit='\t')
            qc = [] # List of counts at different quality scores
            ri = 0      # Read counter (ID counter)
            ridlist = []    # List of read IDs
            rdel = {}   # Dictionary of {rid:current deletion}
            ### ~ [2] Process each entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Replace this with a method to parse a single line and return data or entry
            for line in PILEUP:
                # Split line up into data. Should be: locus, position, reference, no. reads, read data, qualscores
                data = string.split(rje.chomp(line))
                if not data: break
                self.progLog('\r#PARSE','Parsing %s: %s pos...' % (filename,rje.iStr(px)),rand=0.01); px += 1
                ## ~ [2a] Extract Read Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                entry = {'Locus':data[0],'Pos':int(data[1]),'Ref':data[2],'N':int(data[3]),'QN':0}
                entry['Dep'] = rje.sf(entry['N']/meandepth)
                rseq = data[4]  # Read sequences from pileup
                reads = []      # List of read alleles, extracted from rseq
                delx = 0        # Number of deletions: these are not in the count depth (N) nor have a QV
                #!# How should/do people cope with deletions? #!#
                #!# Make it a toggle to include/exclude? #!#
                # Note that ends of reads are [element]$
                # Beginning of reads are ^[QV][element]
                rid = -1
                while '$' in ridlist: ridlist.remove('$')   # Get rid of reads that ended in previous row.
                while rseq:
                    try:
                        # Establish current read ID (Needed for deletions)
                        if rseq[:1] not in '^-+$' and ridlist: rid = ridlist[len(reads)]
                        # Update read data
                        if rseq[:1] in ['.',',']:   # Matches reference (+/- strand)
                            reads.append(entry['Ref']); rseq = rseq[1:]   # Matches reference
                        elif rseq[:1] == '*':   # Check for existing deletions
                            #x#if rid in rdel:
                            if rdel[rid][:1] != entry['Ref']: self.warnLog('Deletion sequence mismatch @ %s:%s' % (data[0],data[1]))
                            rdel[rid] = rdel[rid][1:]
                            if not rdel[rid]: rdel.pop(rid) # End of deletion
                            reads.append('-'); rseq = rseq[1:]
                        elif rseq[:1] == '^':
                            rseq = rseq[2:]   # Indicates a new read
                            ri += 1; ridlist.insert(len(reads),ri)
                            self.bugPrint('Read %d start @ %s:%s' % (ri,data[0],data[1]))
                        elif rseq[:1] in ['-','+']:
                            ilen = string.atoi(rje.matchExp('^(\d+)',rseq[1:])[0])
                            indel = rseq[len('%s' % ilen)+1:][:ilen]    # Just the indel sequences
                            if rseq[:1] == '-':     # Deletion
                                delx += 1           # These do not have QV scores
                                #X wrong! X# reads.append(rseq[:len('%s' % ilen)+ilen+1].upper())    # Deletion includes -n part.
                                rid = ridlist[len(reads)-1]
                                if rid in rdel: self.warnLog('Conflicting deletions for RID %s @ %s:%s' % (rid,data[0],data[1]))
                                rdel[rid] = indel.upper()
                                self.debug(rdel)
                            else:   # Insertion
                                reads[-1] += indel.upper()  # Insertion just has whole sequence AFTER the position
                            rseq = rseq[len('%s' % ilen)+ilen+1:]
                        elif rseq[:1] in ['$']:
                            rseq = rseq[1:]     # Indicated end of a read
                            rfin = ridlist[len(reads)-1]
                            ridlist[len(reads)-1] = '$'
                            self.deBug('Read %d end @ %s:%s' % (rfin,data[0],data[1]))
                        else:
                            if rseq[0].upper() not in 'ATGCN': print ' ???', rseq[0].upper(), '???'
                            # NB. * = single bp deletion, handled below
                            reads.append(rseq[0].upper()); rseq = rseq[1:]
                    except:
                        self.bugPrint(reads)
                        self.bugPrint(ridlist)
                        self.errorLog('!')
                        self.deBug(rseq)
                        raise ValueError
                self.bugPrint('%s:%s : %s'% (entry['Locus'],entry['Pos'],ridlist))
                if len(reads) != len(ridlist):
                    raise ValueError('Read allele versus Read ID count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']))
                if len(reads) != entry['N']:    # Formerly erroenous: (entry['N'] + delx):
                    self.deBug('%s = %d' % (data[4],entry['N']))
                    self.deBug('%s = %d' % (reads,len(reads)))
                    self.errorLog('Read versus Read Count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                    raise ValueError
                ## ~ [2b] Convert Quality Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                qual = []
                for q in data[5]:
                    # Gaps do not have a quality score, so fill these in first - Actually, they do! They are from the row(s) above
                    #!# Make inclusion of gaps optional? Currently assuming quality is good enough.
                    #X# WRONG: #X# while len(qual) < len(reads) and reads[len(qual)][0] == '-': qual.append(self.getInt('QCut'))
                    # Then append actual qv
                    qual.append(ord(q) - 33)
                    qc += [0] * (qual[-1] - len(qc)); qc[qual[-1]-1] += 1
                #while len(qual) < len(reads) and reads[len(qual)][0] == '-': qual.append(self.getInt('QCut'))
                #while '*' in reads: reads[reads.index('*')] = '-'   #'-1%s' % entry['Seq'].upper()
                if len(reads) != len(qual):
                    self.deBug('%s = %d' % (reads,len(reads)))
                    self.deBug('%s = %d' % (qual,len(qual)))
                    self.deBug(data)
                    self.errorLog('Read versus Quality length mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                    raise ValueError
                ## ~ [2c] Filter low quality ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #if entry['Pos'] in [190359]:    #100,98901,183697,169284,
                #    self.deBug(qual)
                #    self.deBug(reads)
                #    self.deBug(qc)
                # Remove (from back) any reads than do not meet QV cutoff
                for r in range(len(qual)-1,-1,-1):
                    if qual[r] < self.getInt('QCut'): qual.pop(r); reads.pop(r)
                entry['QN'] = len(reads)
                if entry['QN'] < self.getInt('MinQN'): continue     # Not enough data for output
                ## ~ [2d] Alleles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                alleles = {}    # Dictionary of {nt:count}
                for read in reads:
                    if read in alleles: alleles[read] += 1
                    else: alleles[read] = 1
                asort = []  # List of tuples for sorting
                for allele in alleles: asort.append((alleles[allele],allele))
                asort.sort()
                asort.reverse()     # This should now be in order of allele frequency
                acheckx = 0         # Count check
                for i in range(len(asort)):
                    acheckx += asort[i][0]
                    asort[i] = '%s:%d' % (asort[i][1],asort[i][0])
                entry['Seq'] = string.join(asort,'|')
                #if entry['Pos'] in [190359]:    #100,98901,183697,169284,
                #    self.deBug(qual)
                #    self.deBug(reads)
                #    self.deBug(alleles)
                #    self.deBug(entry)
                #    self.deBug(line)
                if acheckx != entry['QN']:  #!# Convert these to test statements?
                    self.errorLog('Allele versus Quality count mismatch for %s Pos %s' % (entry['Locus'],entry['Pos']),printerror=False)
                #table.addEntry(entry)
                outlist = []
                for field in outfields: outlist.append(entry[field])
                rje.writeDelimit(PILEOUT,outlist,delimit='\t'); ex += 1
            self.printLog('\r#PARSE','Parsed %s: %s entries from %s lines: %s reads.' % (filename,rje.iStr(ex),rje.iStr(px),rje.iStr(ri)))
            PILEOUT.close()
            PILEUP.close()
            ### ~ [3] Save QC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            qfile = '%s.QC.tdt' % rje.baseFile(filename,strip_path=True)
            if rje.exists(qfile) and not self.force(): self.printLog('#QC','%s already exists (force=F).' % qfile)
            else:
                QC = open(qfile,'w')
                QC.write('Qual\tCount\n')
                for q in range(len(qc)):
                    try: QC.write('%d\t%d\n' % (q+1,qc[q]))
                    except: self.errorLog('!')
                QC.close()
                self.printLog('#QC','Quality score counts output to: %s' % qfile)
            return True
        except: self.errorLog('%s.parsePileup(%s) error' % (self,filename)); return False
#########################################################################################################################
    ### <4> ### Pileup SNP stats methods                                                                                #
#########################################################################################################################
    def pileUpStats(self,snpdb=None,alt=False,locfmt='full'):  ### Calculates statistics of genetic differences from parsed PileUp Tables
        '''Calculates statistics of genetic differences from parsed PileUp Tables.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alt:
                self.printLog('#~~#','## ~~~~~ Parsing & Combining Alt PileUp Tables (%s format) ~~~~~ ##' % locfmt)
                self.printLog('#ALT','AltControl and AltTreatment')
            else:
                self.printLog('#~~#','## ~~~~~ Parsing & Combining PileUp Tables (%s format) ~~~~~ ##' % locfmt)
                self.printLog('#BASE','Output basefile: %s' % self.baseFile())
            if snpdb:
                if alt: outfile = '%s.altcomb.tdt' % self.baseFile()
                else: outfile = '%s.snpcomb.tdt' % self.baseFile()
            else: outfile = '%s.snp.tdt' % self.baseFile()
            if not self.force() and os.path.exists(outfile):
                self.printLog('#SKIP','%s found! (force=F)' % outfile)
                return outfile
            rje.backup(self,outfile)
            outfields = ['Locus','Pos','Ref']
            for fstr in self.list['Labels']:
                for field in ['N','Dep','QN','AN','Seq']:
                    outfields.append('%s|%s' % (field,fstr))
            outfields += ['MajFreq','MajDiff','MajProb']
            ## ~ [0a] Open file handles ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            infields = ['Locus','Pos','Ref','N','QN','Seq','Dep']
            if alt: self.file['Control'] = open(self.getStr('AltControl'),'r')
            else: self.file['Control'] = open(self.getStr('Control'),'r')
            cfields = self.readDelimit('Control')
            if not cfields == infields: raise ValueError('Control data field error! Expect: %s; Read: %s' % (string.join(infields),string.join(cfields)))
            if alt: self.file['Treatment'] = open(self.getStr('AltTreatment'),'r')
            else: self.file['Treatment'] = open(self.getStr('Treatment'),'r')
            tfields = self.readDelimit('Treatment')
            if not tfields == infields: raise ValueError('Treatment data field error! Expect: %s; Read: %s' % (string.join(infields),string.join(tfields)))
            OUTFILE = open(outfile,'w')
            rje.writeDelimit(OUTFILE,outlist=outfields,delimit='\t')

            ### ~ [1] Cycle through data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            locus = None
            cdata = self.readDelimit('Control'); cx = 1
            tdata = None; tx = 0
            bx = 0  # Covered by both
            sx = 0  # SNPs covered by both
            cutx = 0    # Alleles dropped due to low frequency
            difx = 0    # Fixed difference from reference
            refx = 0    # Fixed reference sequence
            nosnpx = 0  # Positions not found in SNP Table
            outx = 0    # Counter of output data lines
            notbix = 0  # Number of SNPs rejected as not biallelic (self.getBool('Biallelic'))
            indx = 0    # Number of indels rejected if not self.getBool('Indels')
            self.progLog('\r#SNP','Combining SNP data: %s Control & %s Treatment => %s common; %s SNP (%s low freq alleles < %s rejected)' % (rje.iStr(cx),rje.iStr(tx),rje.iStr(bx),rje.iStr(sx),rje.iStr(cutx),self.getNum('MinCut')))
            while cdata:
                if cdata[0] != locus:   # New locus
                    locus = cdata[0]
                    tline = self.findline('Treatment',locus)
                    if not tline:   # Locus missing from Treatment
                        self.warnLog('No %s reads in %s!' % (locus,self.getStr('Treatment')))
                        while cdata and cdata[0] == locus: cdata = self.readDelimit('Control'); cx += 1
                        continue
                    tdata = string.split(tline); tx += 1
                    #self.bugPrint('TData: %s' % tdata)
                    #?# Why is the next readDelimit the same?!
                    #self.deBug(self.readDelimit('Treatment'))
                elif not tdata: # Reached end of treatment
                    bx = 0
                    while cdata and cdata[0] == locus: cdata = self.readDelimit('Control'); cx += 1; bx += 1
                    if self.dev(): self.bugPrint('Skipped %d Control lines for %s' % (bx,locus))
                    continue
                elif tdata[0] != locus:   # Overshot locus
                    self.bugPrint('CData: %s' % cdata)
                    self.bugPrint('TData: %s' % tdata)
                    self.debug('What is this?! Should never happen?')
                    while cdata and cdata[0] == locus: cdata = self.readDelimit('Control'); cx += 1
                    continue
                # Should now have cdata and tdata
                if cdata[0] != tdata[0]: raise ValueError('Locus mismatch!')
                if int(cdata[1]) < int(tdata[1]): cdata = self.readDelimit('Control'); cx += 1; continue
                if int(cdata[1]) > int(tdata[1]): tdata = self.readDelimit('Treatment'); tx += 1; continue
                bx += 1
                # Check against SNPDB
                if snpdb:
                    acc = self.mapLocus(locus,locfmt)
                    pos = int(cdata[1])
                    ikey = '%s|%d' % (acc,pos)
                    self.debug(ikey)
                    #?# Add warning for multiple alleles in Pileup data #?#
                    if alt and ikey not in snpdb.index('#AltLocus#|#AltPos#'):  cdata = self.readDelimit('Control'); cx += 1; nosnpx += 1; continue
                    if not alt and ikey not in snpdb.index('#Locus#|#Pos#'):  cdata = self.readDelimit('Control'); cx += 1; nosnpx += 1; continue
                # Convert to dictionaries
                # infields = ['Locus','Pos','Ref','N','QN','Seq','Dep']
                #self.bugPrint('C: %s' % cdata)
                #self.bugPrint('T: %s' % tdata)
                cdict = rje.list2dict(cdata,infields)
                tdict = rje.list2dict(tdata,infields)
                #self.bugPrint('C: %s' % cdict)
                #self.bugPrint('T: %s' % tdict)
                # Combine and assess
                refseq = cdict['Ref']
                totx = 0.0
                alleles = {}
                if not self.getBool('IgnoreRef'): alleles[refseq] = 0
                for allele in string.split('%s|%s' % (cdict['Seq'],tdict['Seq']),'|'):
                    adata = string.split(allele,':')
                    if self.getBool('IgnoreN') and adata[0] == 'N': continue    # Skip N alleles
                    ax = int(adata[1])
                    totx += ax
                    if adata[0] not in alleles: alleles[adata[0]] = 0
                    alleles[adata[0]] += ax
                #self.debug(alleles)
                # Check min allele frequency
                acut = self.getNum('MinCut')
                if acut < 1: acut = max(self.getInt('AbsMinCut'),totx*acut)
                for aseq in rje.sortKeys(alleles):
                    if aseq == refseq and not self.getBool('IgnoreRef'): continue   # Keep Reference allele
                    if alleles[aseq] < acut:    # Too low frequency: dump
                        cutx += 1
                        alleles.pop(aseq)
                if not alleles: raise ValueError('No alleles survived mincut!')
                # Check for SNP
                if len(alleles) == 1:
                    if refseq not in alleles: difx += 1 # Fixed difference from reference
                    else: refx += 1
                    cdata = self.readDelimit('Control'); cx += 1
                    if not snpdb: cdata = self.readDelimit('Control'); cx += 1; continue  # Not a SNP!
                elif len(alleles) == 2 and refseq in alleles and not alleles[refseq]: difx += 1
                elif self.getBool('Biallelic') and len(alleles) != 2: notbix += 1; cdata = self.readDelimit('Control'); cx += 1; continue
                elif '-' in alleles and not self.getBool('Indels'): indx += 1; cdata = self.readDelimit('Control'); cx += 1; continue
                sx += 1
                # Basic Output data with new allele frequencies
                #!# Need to fix this for depth! (Make more versatile)
                major = string.split(tdict['Seq'],':')[0]
                cfreq = 0.0; ctot = 0
                #odata = cdata[:-1] + [0,[]] # Add AN and Allele data
                odata = rje.dict2list(cdict,['Locus','Pos','Ref','N','Dep','QN'])
                odata += [0,[]] # Add AN and Allele data
                for callele in string.split(cdict['Seq'],'|'):
                    if string.split(callele,':')[0] in alleles:
                        odata[-1].append(callele)
                        ctot += int(string.split(callele,':')[1])
                        if string.split(callele,':')[0] == major: cfreq += int(string.split(callele,':')[1])
                odata[-1] = string.join(odata[-1],'|') # Revised alleles
                odata[-2] = ctot    # AN
                tfreq = 0.0; ttot = 0
                #odata += tdata[-3:-1] + [0,[]]
                odata += rje.dict2list(tdict,['N','Dep','QN'])
                odata += [0,[]]   # Add AN and Allele data
                for tallele in string.split(tdict['Seq'],'|'):
                    if string.split(tallele,':')[0] in alleles:
                        odata[-1].append(tallele)
                        ttot += int(string.split(tallele,':')[1])
                        if string.split(tallele,':')[0] == major: tfreq += int(string.split(tallele,':')[1])
                odata[-1] = string.join(odata[-1],'|')
                odata[-2] = ttot
                # Allele frequencies
                cfreq /= ctot
                majx = int(tfreq)
                tfreq /= ttot
                odata += [rje.dp(tfreq,3),rje.dp(tfreq-cfreq,3)]
                if cfreq == 0.0: cfreq = 1.0 / (ctot + 1)   # Assume next read would be Treatment allele
                odata.append(rje.eStr(rje.binomial(majx,ttot,cfreq,exact=False,callobj=self)))
                rje.writeDelimit(OUTFILE,outlist=odata,delimit='\t'); outx += 1
                self.progLog('\r#SNP','Combining SNP data: %s Control & %s Treatment => %s common; %s SNP (%s low freq alleles < %s rejected)' % (rje.iStr(cx),rje.iStr(tx),rje.iStr(bx),rje.iStr(sx),rje.iStr(cutx),self.getNum('MinCut')))
                cdata = self.readDelimit('Control'); cx += 1
            self.printLog('\r#SNP','Combined SNP data: %s Control & %s Treatment => %s common; %s SNP (%s low freq alleles < %s rejected).' % (rje.iStr(cx),rje.iStr(tx),rje.iStr(bx),rje.iStr(sx),rje.iStr(cutx),self.getNum('MinCut')))
            if self.getBool('Biallelic'): self.printLog('#SNP','%s SNPs rejected as not biallelic.' % rje.iStr(notbix))
            if not self.getBool('Indels'): self.printLog('#INDEL','%s Indels rejected (indels=F).' % rje.iStr(indx))
            if snpdb: self.printLog('#FIX','%s combined positions: %s fixed reference; %s fixed differences; %s rejected due to lack of SNP Table data.' % (rje.iStr(bx),rje.iStr(refx),rje.iStr(difx),rje.iStr(nosnpx)))
            else: self.printLog('#FIX','%s combined positions: %s fixed reference; %s fixed differences.' % (rje.iStr(bx),rje.iStr(refx),rje.iStr(difx)))
            self.printLog('#OUT','%s combined SNP data output to %s' % (rje.iStr(outx),outfile))
            return outfile
        except: self.errorLog('%s.pileUpStats() error' % (self)); return None
#########################################################################################################################
    def rateSNPs(self,snpfile):  ### Calculates statistics of genetic differences from parsed PileUp Tables
        '''Calculates statistics of genetic differences from parsed PileUp Tables.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not snpfile: self.printLog('#FDR','No SNP Rating and FDR output.'); return None
            self.printLog('#~~#','## ~~~~~ Rating compiled SNPs ~~~~~ ##')
            alt = self.getStrLC('AltControl') and self.getStrLC('AltTreatment')
            fdrkeys = ['Locus','Pos']
            if alt: fdrkeys += ['AltLocus','AltPos']
            db = self.db()
            fdrfile = '%s.fdr.tdt' % self.baseFile()
            if not self.force() and os.path.exists(fdrfile):
                return db.addTable(fdrfile,mainkeys=fdrkeys,datakeys='All',name='fdr',expect=True)
            snpdb = db.addTable(snpfile,mainkeys=fdrkeys,datakeys='All',name='fdr',expect=True)
            snpdb.dataFormat({'Pos':'int','MajProb':'num','MajDiff':'num'})
            clabel = self.list['Labels'][0]; self.printLog('#LABEL','Control: %s' % clabel)
            tlabel = self.list['Labels'][1]; self.printLog('#LABEL','Treatment: %s' % tlabel)
            ### ~ [1] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MajMut'):
                mx = 0
                for entry in snpdb.entries():
                    if string.split(entry['Seq|%s' % tlabel],':')[0] == entry['Ref']: snpdb.dropEntry(entry); mx += 1
                self.printLog('#DROP','Dropped %s entries with Major %s allele matching Ref' % (rje.iStr(mx),tlabel))
            elif self.getBool('MajRef'):
                mx = 0
                for entry in snpdb.entries():
                    if string.split(entry['Seq|%s' % tlabel],':')[0] != entry['Ref']: snpdb.dropEntry(entry); mx += 1
                self.printLog('#DROP','Dropped %s entries with Major %s allele NOT matching Ref' % (rje.iStr(mx),tlabel))
            ### ~ [2] Calculate FDR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdrx = snpdb.entryNum()     # Total for FDR correction.
            snpdb.dropEntries(['MajProb>%s' % self.getNum('SigCut')],inverse=False,logtxt='SigCut > %s' % self.getNum('SigCut'))
            totx = snpdb.entryNum()     # Total for FDR correction.
            self.printLog('#RANK','Ranking by MajProb for FDR calculation')
            snpdb.rankField('MajProb',newfield='FDR',rev=True,absolute=True,lowest=True,unique=False)
            # FDR field now has number of entries with higher Prob (+1)
            for entry in snpdb.entries():
                self.bugPrint(entry)
                self.debug('%s vs %s' % (entry['MajProb'] * fdrx,(totx - entry['FDR'] + 1)))
                try: entry['FDR'] = entry['MajProb'] * fdrx / (totx - entry['FDR'] + 1)
                except:
                    self.warnLog('FDR Error: %s vs %s' % (entry['MajProb'] * fdrx,(totx - entry['FDR'] + 1)))
                    entry['FDR'] = 1.0  #???
            ### ~ [3] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MajDif'):
                mx = 0
                for entry in snpdb.entries():
                    if string.split(entry['Seq|%s' % tlabel],':')[0] == string.split(entry['Seq|%s' % clabel],':')[0]: snpdb.dropEntry(entry); mx += 1
                self.printLog('#DROP','Dropped %s entries with Major %s allele matching %s' % (rje.iStr(mx),tlabel,clabel))
            snpdb.dropEntries(['MajFreq<%s' % self.getNum('MajCut')],inverse=False,logtxt='MajFreq < majcut=%s' % self.getNum('MajCut'))
            snpdb.dropEntries(['MajDiff<0'],inverse=False,logtxt='%s > %s' % (tlabel,clabel))
            snpdb.dropEntries(['FDR>%s' % self.getNum('FDRCut')],inverse=False,logtxt='FDRCut > %s' % self.getNum('FDRCut'))
            for entry in snpdb.entries(): entry['FDR'] = rje.eStr(entry['FDR'])
            snpdb.saveToFile()
            return snpdb
        except: self.errorLog('%s.pileUpFDR() error' % (self)); return None
#########################################################################################################################
    def loadSNPTable(self,name='SNP'):     ### Loads and returns a SNP table for the genome(s) being analysed
        '''Loads and returns a SNP table for the genome(s) being analysed.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SNPTable'): self.printLog('\r#SNP','No SNP table to add.'); return None
            snptable = self.getStr('SNPTable')
            snpkeys = self.list['SNPTabKeys']
            for field in ['Pos','Locus']:
                if field in snpkeys: snpkeys.remove(field)
                snpkeys.insert(0,field)
            self.printLog('#KEYS','Loading %s with keys: %s' % (snptable,string.join(snpkeys,'; ')))
            snpdb = self.db().addTable(snptable,name=name,expect=True,mainkeys=snpkeys)
            snpdb.dataFormat({'Pos':'int','AltPos':'int'})
            if not self.getBool('Indels'):
                indx = snpdb.entryNum()
                if 'REF' in snpdb.fields(): snpdb.dropEntriesDirect('REF',['-'])
                else: snpdb.dropEntriesDirect('Ref',['-'])
                if 'ALT' in snpdb.fields(): snpdb.dropEntriesDirect('ALT',['-'])
                else: snpdb.dropEntriesDirect('Alt',['-'])
                indx -= snpdb.entryNum()
                self.printLog('#INDEL','%s Indels dropped (indels=F).' % rje.iStr(indx))
            return snpdb
        except: self.errorLog('%s.loadSNPTable() error' % (self)); return None
#########################################################################################################################
    def mapLocus(self,locus,locformat='full'):   ### Reduce locus to accession number for mapping
        '''Reduce locus to accession number for mapping.'''
        acc = locus
        if locformat == 'fullnum': return acc
        if rje.matchExp('^(\S+_\S+__\S+)\.\d+$',acc): acc = rje.matchExp('^(\S+_\S+__\S+)\.\d+$',acc)[0]
        if locformat == 'full': return acc
        if rje.matchExp('^\S+_\S+__(\S+)',acc): acc = rje.matchExp('^\S+_\S+__(\S+)',acc)[0]
        if locformat == 'acc': return acc
        if rje.matchExp('^(\S+)\.\d+$',acc): acc = rje.matchExp('^(\S+)\.\d+$',acc)[0]
        return acc
#########################################################################################################################
    def locusFormat(self,locus):   ### Returns locus format for later mapping.
        '''Returns locus format for later mapping.'''
        acc = locus
        if rje.matchExp('^\S+_\S+__(\S+)\.\d+$',acc): return 'fullnum'
        if rje.matchExp('^\S+_\S+__(\S+)',acc): return 'full'
        if rje.matchExp('^(\S+)\.\d+$',acc): return 'acc'
        return 'locus'
#########################################################################################################################
    def combineSNPs(self,fdb):  ### Calculates statistics of genetic differences from parsed PileUp Tables
        '''Calculates statistics of genetic differences from parsed PileUp Tables.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SNPTable'): self.printLog('\r#SNP','No SNP table to add.'); return fdb
            self.printLog('#~~#','## ~~~~~ Combining PileUp SNPs with SNPTable data ~~~~~ ##')
            snptable = self.getStr('SNPTable')
            snpkeys = self.list['SNPTabKeys']
            for field in ['Pos','Locus']:
                if field in snpkeys: snpkeys.remove(field)
                snpkeys.insert(0,field)
            snpdb = self.db().addTable(snptable,name='SNP',expect=True,mainkeys=snpkeys)
            if not fdb: fdb = self.db().addTable(name='fdr',expect=True,mainkeys=['Locus','Pos'])
            fdb.remakeKeys()   #!# Delete once tuple thing OK
            fdbkeys = fdb.dataKeys()
            if len(self.list['SNPTabMap']) > 2:
                self.warnLog('#SNPMAP','SNPTabMap %s reduced to %s.' % (string.join(self.list['SNPTabMap'],','),string.join(self.list['SNPTabMap'][:2],',')))
                self.list['SNPTabMap'] = self.list['SNPTabMap'][:2]
            elif not self.list['SNPTabMap']: self.list['SNPTabMap'] = ['Locus','Pos']
            elif len(self.list['SNPTabMap']) == 1:
                raise ValueError('SNPTabMap (%s) needs 2 elements!' % self.list['SNPTabMap'][0])
            for jfield in self.list['SNPTabMap']:
                if jfield not in snpdb.fields(): raise ValueError('SNPTabMap field (%s) not in SNPTable!' % jfield)
            snpjoin = self.list['SNPTabMap']
            ## ~ [0a] Check Loci naming ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            shared = rje.listIntersect(fdb.indexKeys('Locus'),snpdb.indexKeys(snpjoin[0]))
            self.printLog('#LOCUS','%s shared FDR and SNP Table loci.' % rje.iLen(shared))
            if shared:
                fdb.makeField('#Locus#|#Pos#')
                snpjoinfield = '#%s#|#%s#' % (snpjoin[0],snpjoin[1])
                snpdb.makeField(snpjoinfield)
                if snpjoinfield != '#Locus#|#Pos#': snpdb.renameField(snpjoinfield,'#Locus#|#Pos#')
            else:
                self.printLog('#LOCUS','Will try pulling out accession numbers for matching.')
                for table in [fdb,snpdb]:
                    table.addField('#Locus#|#Pos#')
                    for entry in table.entries():
                        if table == snpdb: acc = entry[snpjoin[0]]
                        else: acc = entry['Locus']
                        if rje.matchExp('^\S+_\S+__(\S+)',acc): acc = rje.matchExp('^\S+_\S+__(\S+)',acc)[0]
                        if rje.matchExp('^(\S+)\.\d+$',acc): acc = rje.matchExp('^(\S+)\.\d+$',acc)[0]
                        if table == snpdb: entry['#Locus#|#Pos#'] = '%s|%s' % (acc, entry[snpjoin[1]])
                        else: entry['#Locus#|#Pos#'] = '%s|%s' % (acc, entry['Pos'])
            ### ~ [1] Join Tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newsnpkeys = snpkeys[0:]
            if self.list['SNPTabMap'] != ['Locus','Pos']:  # Mapping onto a different field
                for field in self.list['SNPTabMap']:
                    newsnpkeys.remove(field)
                newsnpkeys += ['SNP_Locus','SNP_Pos']
                self.printLog('#KEYS','Join table keys: %s -> %s' % (snpkeys,newsnpkeys))
            joindb = self.db().joinTables(name='snpmap',join=[(fdb,'#Locus#|#Pos#'),(snpdb,'#Locus#|#Pos#')],newkey=newsnpkeys,keeptable=True)
            self.printLog('#SNP','Added SNPs from %s' % snptable)
            joindb.dropField('#Locus#|#Pos#')
            for field in ['SNP_Locus','SNP_Pos']:
                if field not in newsnpkeys: joindb.dropField(field)
            joindb.saveToFile()
            return joindb
        except: self.errorLog('%s.pileUpStats() error' % (self)); return None
#########################################################################################################################
    ### <5> ### Combined Pileup SNP stats methods                                                                       #
#########################################################################################################################
    def combinePileUpStats(self):  ### Combines statistics of genetic differences from parsed PileUp of 2 genomes
        '''Combines statistics of genetic differences from parsed PileUp Tables of 2 genomes.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ Combining PileUp SNP Statistics ~~~~~ ##')
            db = self.db()
            revcomp = {'A':'T','C':'G','G':'C','T':'A','-':'-'}
            revx = 0    # Number of sites reverse complemented
            ## ~ [0a] Load SNP Mapping table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb = self.loadSNPTable('combined')
            if not snpdb: raise IOError('Cannot perform AltControl/AltTreatment analysis without SNPTable.')
            self.debug(snpdb.keys())
            mapfields = ['Locus','Pos','AltLocus','AltPos','#Locus#|#Pos#','#AltLocus#|#AltPos#']
            for field in mapfields[0:4]:
                if field not in snpdb.fields(): raise ValueError('Cannot perform AltControl/AltTreatment analysis without SNPTable field "%s".' % field)
            ## ~ [0b] Drop Indels: cannot match positions via short reads ~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb.dropEntriesDirect('REF',['-'])
            snpdb.dropEntriesDirect('ALT',['-'])
            ## ~ [0c] Sort out keys and locus mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            snpdb.addField('#Locus#|#Pos#')
            snpdb.addField('#AltLocus#|#AltPos#')
            for entry in snpdb.entries():
                #entry['#Locus#|#Pos#'] = '%s|%s' % (self.mapLocus(entry['Locus']),entry['Pos'])
                #entry['#AltLocus#|#AltPos#'] = '%s|%s' % (self.mapLocus(entry['AltLocus']),entry['AltPos'])
                # Now assuming that the pileup files with have fuller versions of identifiers.
                entry['#Locus#|#Pos#'] = '%s|%s' % (entry['Locus'],entry['Pos'])
                entry['#AltLocus#|#AltPos#'] = '%s|%s' % (entry['AltLocus'],entry['AltPos'])
            snpdb.compress(mapfields[:4])
            snpdb.keepFields(mapfields+['REF','ALT'])
            snpdb.index('#Locus#|#Pos#')
            snpdb.index('#AltLocus#|#AltPos#')
            snpdb.remakeKeys()
            self.debug(snpdb.datakeys()[:10])
            sentry = snpdb.entries()[0]
            ### ~ [1] Generate SNP tables for the two genomes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            refdb = db.addTable(self.pileUpStats(snpdb,alt=False,locfmt=self.locusFormat(sentry['Locus'])),mainkeys=['Locus','Pos'],datakeys='All',name='ref',expect=True)
            refdb.dataFormat({'Pos':'int','MajProb':'num','MajDiff':'num'})
            altdb = db.addTable(self.pileUpStats(snpdb,alt=True,locfmt=self.locusFormat(sentry['AltLocus'])),mainkeys=['Locus','Pos'],datakeys='All',name='alt',expect=True)
            altdb.dataFormat({'Pos':'int','MajProb':'num','MajDiff':'num'})
            locfmt = {refdb:self.locusFormat(sentry['Locus']),altdb:self.locusFormat(sentry['AltLocus'])}
            for samdb in [refdb,altdb]:
                samdb.addField('#Locus#|#Pos#')
                for entry in samdb.entries(): entry['#Locus#|#Pos#'] = '%s|%s' % (self.mapLocus(entry['Locus'],locfmt[samdb]),entry['Pos'])
                samdb.index('#Locus#|#Pos#')
            ### ~ [2] Combine SNP tables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#~~#','## ~~~~~ Combining SNP Tables ~~~~~ ##')
            if not self.getBool('MajFocus') and self.getBool('MajMut'): self.printLog('#ALTSNP','Generating output for all Alt ("Mutant") alleles.')
            if not self.getBool('MajFocus') and self.getBool('MajRef'): self.printLog('#REFSNP','Generating output for all Alt ("Mutant") alleles.')
            clabel = self.list['Labels'][0]
            tlabel = self.list['Labels'][1]
            combfields = ['Locus','Pos','AltLocus','AltPos','Ref','Alt','AN|%s' % clabel,'Seq|%s' % clabel,'AN|%s' % tlabel,'Seq|%s' % tlabel,'MajFreq','MajDiff','MajProb']
            snpdb.list['Fields'] += combfields[4:]
            ex = 0.0; etot = snpdb.entryNum(); snpx = 0; nosnpx = 0; nofocusx = 0
            for entry in snpdb.entries():
                self.progLog('\r#SNP','Combining %s Ref and Alt SNPs: %.1f%%' % (rje.iStr(snpx),ex/etot)); ex += 100.0

                #!# This needs to be fixed if it wanted to keep entries that only appear in one or the other!
                if len(refdb.indexEntries('#Locus#|#Pos#',entry['#Locus#|#Pos#'])) < 1: entry['Ref'] = 'X'; continue
                if len(altdb.indexEntries('#Locus#|#Pos#',entry['#AltLocus#|#AltPos#'])) < 1: entry['Alt'] = 'X'; continue
                refentry = refdb.indexEntries('#Locus#|#Pos#',entry['#Locus#|#Pos#'])[0]
                altentry = altdb.indexEntries('#Locus#|#Pos#',entry['#AltLocus#|#AltPos#'])[0]
                if len(refdb.indexEntries('#Locus#|#Pos#',entry['#Locus#|#Pos#'])) > 1:
                    self.warnLog('%s "%s" #Locus#|#Pos# Ref SAM entries!' % (len(refdb.indexEntries('#Locus#|#Pos#',entry['#Locus#|#Pos#'])),entry['#Locus#|#Pos#']))
                if len(altdb.indexEntries('#Locus#|#Pos#',entry['#AltLocus#|#AltPos#'])) > 1:
                    self.warnLog('%s "%s" #Locus#|#Pos# Alt SAM entries!' % (len(altdb.indexEntries('#Locus#|#Pos#',entry['#AltLocus#|#AltPos#'])),entry['#AltLocus#|#AltPos#']))
                entry['Ref'] = refentry['Ref']

                #!# Check for reverse complement sequence and adjust alleles if required! #!#
                if entry['REF'] == entry['ALT']:    # Not a SNP! (Snapper included feature start/end positions)
                    entry['Alt'] = entry['Ref'] = '.'
                    nosnpx += 1
                    continue
                elif entry['ALT'] not in revcomp: #Fix later for GATT etc.
                    self.warnLog('Odd Alt Sequence: %s' % entry)
                elif altentry['Ref'] == rje_sequence.reverseComplement(entry['ALT']):    # Hit on opposite strand!
                    try: altentry['Ref'] = rje_sequence.reverseComplement(altentry['Ref'])
                    except: self.warnLog('Cannot RevComp "%s"' % altentry['Ref'])
                    for ct in [clabel,tlabel]:
                        revcompalleles = []
                        for allele in string.split(altentry['Seq|%s' % ct],'|'):
                            [nt,ntx] = string.split(allele,':')
                            nt = rje_sequence.reverseComplement(nt)
                            revcompalleles.append(string.join([nt,ntx],':'))
                        altentry['Seq|%s' % ct] = string.join(revcompalleles,'|')
                    revx += 1
                elif entry['ALT'] == '-': continue
                elif altentry['Ref'] != entry['ALT']: self.warnLog('Alt Sequence mismatch: %s' % entry)
                entry['Alt'] = altentry['Ref']

                # Combine alleles
                for ct in [tlabel,clabel]:
                    asort = []  # List of tuples for sorting
                    mutallele = None
                    refallele = None
                    for aseq in string.split(refentry['Seq|%s' % ct],'|') + string.split(altentry['Seq|%s' % ct],'|'):
                        adata = string.split(aseq,':')  # Allele:count
                        asort.append((int(adata[1]),adata[0]))
                        if adata[0] == entry['Alt']: mutallele = asort[-1]
                        if adata[0] == entry['Ref']: refallele = asort[-1]
                    asort.sort()
                    asort.reverse()     # This should now be in order of allele frequency
                    # Move Mutant or Reference allele to "Major" position if majdif=F and majmut=T or majref=T
                    #!# Might need some extra checks for this that mutalle and refallele are present!
                    if not self.getBool('MajFocus'):
                        if self.getBool('MajMut'):
                            if mutallele: asort.remove(mutallele)
                            else: mutallele = (0,entry['Alt'])
                            asort.insert(0,mutallele)
                        elif self.getBool('MajRef'):
                            if refallele: asort.remove(refallele)
                            else: refallele = (0,entry['Ref'])
                            asort.insert(0,refallele)
                        else:
                            self.warnLog('Cannot set "Major" allele for %s' % entry,'majfocus',suppress=True)
                    alleles = {}; aseq = []
                    acount = 0
                    for i in range(len(asort)):
                        if asort[i][1] in alleles: continue    # Partial hit from other mapping
                        #i# NOTE: This explicitly assumes that each allele will 100% hit either ref or alt
                        acount += asort[i][0]
                        alleles[asort[i][1]] = asort[i][0]
                        aseq.append('%s:%d' % (asort[i][1],asort[i][0]))
                    entry['Seq|%s' % ct] = string.join(aseq,'|')
                    entry['AN|%s' % ct] = acount
                    # Calculate 'MajFreq','MajDiff','MajProb'
                    # The Major Allele is the main one present in the Treatment group!
                    major = string.split(entry['Seq|%s' % tlabel],':')[0]
                    if ct == clabel:
                        cfreq = 0.0
                        if major in alleles: cfreq = float(alleles[major])/acount
                        ctot = acount
                    else:
                        ttot = acount
                        tfreq = float(alleles[major])/acount
                        majx = alleles[major]   # Observed number in treatment
                entry['MajFreq'] = rje.dp(tfreq,3)
                entry['MajDiff'] = rje.dp(tfreq-cfreq,3)
                if cfreq == 0.0: cfreq = 1.0 / (ctot + 1)   # Assume next read would be Treatment allele
                entry['MajProb'] = rje.eStr(rje.binomial(majx,ttot,cfreq,exact=False,callobj=self))
                snpx += 1
            self.printLog('\r#SNP','Combined %s of %s Ref and Alt SNPs. %s RevComp sites. %s not SNPs in SNP Table.'  % (rje.iStr(snpx),rje.iStr(etot),rje.iStr(revx),rje.iStr(nosnpx)))
            outfile = '%s.combinedsnp.tdt' % self.baseFile()
            ## Check REF and ALT match Ref and Alt
            snpdb.dropEntriesDirect('Ref',['X','.'])
            snpdb.dropEntriesDirect('Alt',['X','.'])
            snpdb.fillBlanks(fields=['Ref','Alt'])
            for entry in snpdb.entries():
                if entry['ALT'] != entry['Alt'] and entry['Alt'] != 'X':
                    self.warnLog('ALT/Alt Sequence mismatch: %s' % entry)   #!# Make these better
                    entry['Alt'] = 'X'
                if entry['REF'] != entry['Ref'] and entry['Ref'] != 'X':
                    self.warnLog('REF/Ref Sequence mismatch: %s' % entry)   #!# Make these better
                    entry['Ref'] = 'X'
            self.debug(snpdb.keys())
            snpdb.dropEntriesDirect('Ref',['X'])
            snpdb.dropEntriesDirect('Alt',['X'])
            snpdb.newKey(['Locus','Pos','AltLocus','AltPos'])
            snpdb.dropFields(['#Locus#|#Pos#','#AltLocus#|#AltPos#','ALT','REF'])
            snpdb.saveToFile(outfile)
            return outfile
        except:
            self.errorLog('%s.combinePileUpStats() error' % (self))
            return None
#########################################################################################################################
    def pileUpFDR(self):  ### Calculates statistics of genetic differences from parsed PileUp Tables
        '''Calculates statistics of genetic differences from parsed PileUp Tables.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fdrfile = '%s.fdr.tdt' % self.baseFile()
            if not self.force() and os.path.exists(fdrfile): return 
            sigpval = {}    # pval:[fpos]
            npos = 0; nx = 0
            for locus in rje.sortKeys(self.dict['RefSeq']):
                npos += len(self.dict['RefSeq'][locus]) - self.dict['RefSeq'][locus].count('?')
            ### ~ [1] Parse out stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SAMSIG = open('%s.pdiff.tdt' % self.baseFile(),'r')
            headers = string.split(SAMSIG.readline()) + ['p.FDR']
            fpos = SAMSIG.tell(); fline = SAMSIG.readline(); px = 0
            while fline:
                self.progLog('\r#SIG','Reading Pvalues: %s p <= 0.05...' % rje.iStr(px))
                try: pval = float(string.split(fline)[-1])
                except: break
                if pval <= 0.05:
                    if pval not in sigpval: sigpval[pval] = []
                    sigpval[pval].append(fpos); px += 1
                fpos = SAMSIG.tell(); fline = SAMSIG.readline()
            self.printLog('\r#SIG','Reading Pvalues complete: %s p <= 0.05.' % rje.iStr(px))
            ### ~ [2] Calculate FDR and output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SAMFDR = open(fdrfile,'w')
            rje.writeDelimit(SAMFDR, headers)
            px = 0; sx = 0.0; stot = len(sigpval)
            for pval in rje.sortKeys(sigpval):
                self.progLog('\r#FDR','Calculating FDR: %.2f%%' % (sx/stot)); sx += 100.0
                px += len(sigpval[pval])
                if pval: fdr = (pval * npos) / px
                else: fdr = 0.0
                for fpos in sigpval[pval]:
                    SAMSIG.seek(fpos)
                    rje.writeDelimit(SAMFDR,rje.readDelimit(SAMSIG.readline())+[rje.expectString(fdr)])
            SAMSIG.close()
            SAMFDR.close()
            self.printLog('\r#FDR','%s FDR lines output to %s' % (rje.iStr(px),fdrfile))
        except: self.errorLog('%s.pileUpFDR() error' % (self)); return None
#########################################################################################################################
    ### <6> ### Read coverage methods                                                                                   #
#########################################################################################################################
    def coverageFromRID(self,ridfile=None):  ### Extracts read data from RID file and summarises read coverage
        '''
        Extracts read data from SAM file.
        >> filename:str = Pileup file name
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            db = self.db()
            covflanks = self.list['CheckFlanks']
            if not len(self.list['CheckFields']) == 3:
                raise ValueError('checkfields=LIST must have exactly 3 elements: Locus, Start, End. %d found!' % len(self.list['CheckFields']))
            [locusfield,startfield,endfield] = self.list['CheckFields']
            ## ~ [1a] Setup Check Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # This will be set by self.getStr('ReadCheck') but might need parsing of SAM/Pileup
            # ['RID','Locus','Start','End']
            if ridfile: rdb = db.addTable(ridfile,mainkeys=['RID'],name='rid',expect=True)
            else: rdb = self.db('rid',add=True,mainkeys=['RID'])
            if not rdb: raise ValueError('Cannot perform coverage analysis without RID table')
            rdb.dataFormat({'Start':'int','End':'int'})
            self.printLog('#~~#','## ~~~~~ Calculating read coverage ~~~~~ ##')
            ## ~ [1b] Setup Check Table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # This is self.getStr('CheckPos')
            if self.getStrLC('CheckPos'):
                cdb = db.addTable(self.getStr('CheckPos'),mainkeys=self.list['CheckFields'],name='check',expect=True)
                cdb.dataFormat({'Start':'int','End':'int'})
                for covx in covflanks: cdb.addField('Span%d' % covx,evalue=0)
                cdb.addField('MeanX',evalue=0.0)
                shared = rje.listIntersect(rdb.indexKeys('Locus'),cdb.indexKeys(locusfield))
                if not shared:
                    self.warnLog('No overlapping Loci in RID file and CheckPos file. Will attempt AccNum conversion.')
                    for entry in rdb.entries(): entry['Locus'] = string.split(entry['Locus'],'__')[-1]
                    rdb.remakeKeys()
                    for entry in cdb.entries(): entry[locusfield] = string.split(entry[locusfield],'__')[-1]
                    cdb.remakeKeys()
                    shared = rje.listIntersect(rdb.indexKeys('Locus'),cdb.indexKeys(locusfield))
                self.printLog('#LOCI','%s shared loci between RID and CheckPos files.' % rje.iLen(shared))
                cdict = cdb.indexEntries(locusfield,asdict=True)
            else: cdb = None; cdict = {}    # No position checking
            ## ~ [1c] Setup locus sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            # This is from seqin=FILE so just needs SeqList object
            #if self.getStrLC('SeqIn'):
            #   seqlist = rje_seqlist.SeqList(self.log,self.cmd_list)
            #   seqdict = seqlist.makeSeqNameDic('max')
            #else:
            seqdict = {}
            #calculate mean depth for check table too so could give start/end positions here if desired

            ### ~ [2] Calculate read coverage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            covdb = db.addEmptyTable('coverage',['Locus','Length','MeanX'],['Locus'])
            rx = 0.0; rtot = rdb.entryNum(); lx = 0; ltot = len(rdb.index('Locus'))
            for locus in rdb.indexKeys('Locus'):
                lx += 1
                self.progLog('\r#COV','Calculating read coverage: %.2f%% (Locus %d of %d)' % (rx/rtot,lx,ltot),rand=0.01); rx +=1
                centry = {'Locus':locus,'Length':0,'MeanX':0.0}
                if locus in seqdict: centry['Length'] = seqlist.seqLen(seqdict[locus])
                else: centry['Length'] = max(rdb.indexDataList('Locus',locus,'End',sortunique=False))
                for rentry in rdb.indexEntries('Locus',locus):
                    self.progLog('\r#COV','Calculating read coverage: %.2f%% (Locus %d of %d)' % (rx/rtot,lx,ltot),rand=0.01); rx +=1
                    centry['MeanX'] += rentry['End'] - rentry['Start'] + 1
                    if locus in cdict:
                        # Cycle through CheckPos entries
                        for lentry in cdict[locus]:
                            # Skip if no overlap
                            if lentry[startfield] > rentry['End'] or lentry[endfield] < rentry['Start']: continue
                            # Check for span
                            for covx in covflanks:
                                if rentry['Start'] <= max(1,lentry[startfield]-covx) and rentry['End'] >= min(centry['Length'],lentry[endfield]+covx):
                                    lentry['Span%d' % covx] += 1
                            # Update X coverage as proportion of checkpos entry covered by read
                            spany = min(rentry['End'],lentry[endfield])
                            spanx = max(rentry['Start'],lentry[startfield])
                            lentry['MeanX'] += (spany - spanx + 1.0) / (lentry[endfield] - lentry[startfield] + 1.0)
                centry['MeanX'] /= centry['Length']
                covdb.addEntry(centry)
            self.printLog('\r#COV','Calculating read coverage: complete (%d of %d loci)' % (lx,ltot))
            if cdb:
                #!# Remerge the cdict tables for output!
                cdb.saveToFile()
            covdb.saveToFile()
            return True
        except: self.errorLog('%s.coverageFromRID() error' % (self.prog())); return False
#########################################################################################################################
### End of SECTION II: SAMtools Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MODULE METHODS                                                                                         #
#########################################################################################################################
def parseCigar(cigarstr): ### Returns dictionary of parsed data from cigar string
    '''
    Returns dictionary of parsed data from cigar string:
    M	 	Match (alignment column containing two letters). This could contain two different letters (mismatch) or two identical letters. USEARCH generates CIGAR strings containing Ms rather than X's and ='s (see below).
    D	 	Deletion (gap in the target sequence).
    I	 	Insertion (gap in the query sequence).
    S	 	Segment of the query sequence that does not appear in the alignment. This is used with soft clipping, where the full-length query sequence is given (field 10 in the SAM record). In this case, S operations specify segments at the start and/or end of the query that do not appear in a local alignment.
    H	 	Segment of the query sequence that does not appear in the alignment. This is used with hard clipping, where only the aligned segment of the query sequences is given (field 10 in the SAM record). In this case, H operations specify segments at the start and/or end of the query that do not appear in the SAM record.
    =	 	Alignment column containing two identical letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
    X	 	Alignment column containing a mismatch, i.e. two different letters. USEARCH can read CIGAR strings using this operation, but does not generate them.
    Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    '''
    cigdata = {}
    cigstr = cigarstr
    while cigstr:
        (clen,ctype) = rje.matchExp('^(\d+)(\D)',cigstr)
        clen = int(clen)
        if ctype in cigdata: cigdata[ctype] += clen
        else: cigdata[ctype] = clen
        cigstr = cigstr[len('%d%s' % (clen,ctype)):]
    return cigdata
#########################################################################################################################
def cigarAlnLen(cigdata):   ### Returns the total length of the aligned region in the reference sequence
    '''Returns the total length of the aligned region in the reference sequence.'''
    alnlen = 0
    for ctype in 'MD=X':
        if ctype in cigdata: alnlen += cigdata[ctype]
    return alnlen
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
    try: SAMtools(mainlog,cmd_list).run()

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
