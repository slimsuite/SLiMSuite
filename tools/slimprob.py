#!/usr/bin/python

# See below for name and description
# Copyright (C) 2007 Richard J. Edwards <redwards@cabbagesofdoom.co.uk>
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
Program:      SLiMProb
Description:  Short Linear Motif Probability tool
Version:      2.2.5
Last Edit:    07/12/15
Citation:     Davey, Haslam, Shields & Edwards (2010), Lecture Notes in Bioinformatics 6282: 50-61. 
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMProb is a tool for finding pre-defined SLiMs (Short Linear Motifs) in a protein sequence database. SLiMProb
    can make use of corrections for evolutionary relationships and a variation of the SLiMChance alogrithm from
    SLiMFinder to assess motifs for statistical over- and under-representation. SLiMProb is replace for the original
    SLiMSearch, which itself was a replacement for PRESTO. The basic architecture is the same but it was felt that having
    two different "SLiMSearch" servers was confusing. 

    Benefits of SLiMProb that make it more useful than a lot of existing tools include:
    * searching with mismatches rather than restricting hits to perfect matches.
    * optional equivalency files for searching with specific allowed mismatched (e.g. charge conservation)
    * generation or reading of alignment files from which to calculate conservation statistics for motif occurrences.
    * additional statistics, including protein disorder, surface accessibility and hydrophobicity predictions
    * recognition of "n of m" motif elements in the form <X:n:m>, where X is one or more amino acids that must occur n+
    times across which m positions. E.g. <IL:3:5> must have 3+ Is and/or Ls in a 5aa stretch.

    Main output for SLiMProb is a delimited file of motif/peptide occurrences but the motifaln=T and proteinaln=T also
    allow output of alignments of motifs and their occurrences. The primary outputs are named *.occ.csv for the occurrence
    data and *.csv for the summary data for each motif/dataset pair. (This is a change since SLiMSearch.)

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ###
    motifs=FILE     : File of input motifs/peptides (also motif=X) [None]
                      Single line per motif format = 'Name Sequence #Comments' (Comments are optional and ignored)
                      Alternative formats include fasta, SLiMDisc output and raw motif lists.
    seqin=SEQFILE   : Sequence file to search. Over-rules batch=FILE and uniprotid=LIST [None]
    batch=FILELIST  : List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    uniprotid=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.dat and use as seqin=FILE. []
    maxseq=X        : Maximum number of sequences to process [0]
    maxsize=X       : Maximum dataset size to process in AA (or NT) [100,000]
    maxocc=X        : Filter out Motifs with more than maximum number of occurrences [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resfile=FILE    : Main SLiMProb results table (*.csv and *.occ.csv) [slimprob.csv]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMProb/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMProb/]
    force=T/F       : Force re-running of BLAST, UPC generation and search [False]
    dna=T/F         : Whether the sequences files are DNA rather than protein [False]
    alphabet=LIST   : List of characters to include in search (e.g. AAs or NTs) [default AA or NT codes]
    megaslim=FILE   : Make/use precomputed results for a proteome (FILE) in fasta format [None]
    megablam=T/F    : Whether to create and use all-by-all GABLAM results for (gablamdis) UPC generation [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SearchDB Options I: Input Protein Sequence Masking ###
    masking=T/F     : Master control switch to turn off all masking if False [False]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out (True=EM,DOMAIN,TRANSMEM) []
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [None]
    casemask=X      : Mask Upper or Lower case [None]
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    metmask=T/F     : Masks the N-terminal M [False]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas  []
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []

    ### SearchDB Options II: Evolutionary Filtering  ###
    efilter=T/F     : Whether to use evolutionary filter [True]
    blastf=T/F      : Use BLAST Complexity filter when determining relationships [True]
    blaste=X        : BLAST e-value threshold for determining relationships [1e=4]
    altdis=FILE     : Alternative all by all distance matrix for relationships [None]
    gablamdis=FILE  : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
    occupc=T/F      : Whether to output the UPC ID number in the occurrence output file [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SLiMChance Options ###
    maskfreq=T/F    : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [True]
    aafreq=FILE     : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    aadimerfreq=FILE: Use empirical dimer frequencies from FILE (fasta or *.aadimer.tdt) [None]
    negatives=FILE  : Multiply raw probabilities by under-representation in FILE [None]
    background=FILE : Use observed support in background file for over-representation calculations [None]
    smearfreq=T/F   : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    seqocc=X        : Restrict to sequences with X+ occurrences (adjust for high frequency SLiMs) [1]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Output Options ###
    extras=X        : Whether to generate additional output files (alignments etc.) [2]
                        - 0 = No output beyond main results file
                        - 1 = Saved masked input sequences [*.masked.fas]
                        - 2 = Generate additional outputs (alignments etc.)
                        - 3 = Additional distance matrices for input sequences
    pickle=T/F      : Whether to save/use pickles [True]
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc and *.pickle files
                        - 2 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)
    * See also rje_slimcalc options for occurrence-based calculations and filtering *
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import pickle, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_sequence, rje_scoring, rje_slim, rje_slimcore, rje_slimcalc, rje_slimlist, rje_zen
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 1.0 - SLiMProb 1.0 based on SLiMSearch 1.7. Altered output files to be *.csv and *.occ.csv.
    # 1.1 - Tidied import commands.
    # 1.2 - Increased extras=X levels. Adjusted maxsize=X assessment to be post-masking.
    # 1.3 - Consolidating output file naming for consistency across SLiMSuite. (SLiMBuild = Motif input)
    # 1.4 - Preparation for SLiMProb V2.0 & SLiMCore V2.0 using newer RJE_Object.
    # 2.0 - Converted to use rje_obj.RJE_Object as base. Version 1.4 moved to legacy/.
    # 2.1 - Modified output of N-terminal motifs to correctly start at position 1.
    # 2.2.0 - Added basic REST functionality.
    # 2.2.1 - Updated REST output.
    # 2.2.2 - Modified input to allow motif=X in addition to motifs=X.
    # 2.2.3 - Tweaked basefile setting and citation.
    # 2.2.4 - Improved slimcalc output (s.f.).
    # 2.2.5 - Fixed FTMask=T/F bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Finish implementation with current methods/options!
    # [ ] : Add occurrence and SLiM filtering.
    # [ ] : Reinstate the expcut=X option for filtering motifs based on expected occurrences.
    # [Y] : Remove ORTHALN and make mapping like SLiMFinder
    # [ ] : Add "startfrom" option.
    # [Y] : Incorporate and test new BLAST module (blast+)
    # [ ] : Test/Check TEIRESIAS output
    # [ ] : Add numpy etc.
    # [ ] : Make sure that Force option does not use old pickles.
    # [ ] : Generally tidy up use of force/debug/test/warn/dev/v/i options. Improve error handling and checking.
    # [Y] : Adjust the maxsize=X setting to be enforced AFTER masking!
    # [ ] : Make slimsearch variant (prog=slimsearch efilter=F)
    # [ ] : ResFile default is currently over-ruling basefile=X. Change and test that REST server is still OK.
    # [ ] : Add pickup=T/F option, as with SLiMFinder.
    # [ ] : Add maxupc=X        : Maximum UPC size of dataset to process [0]
    # [ ] : Review masking defaults.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('SLiMProb', '2.2.5', 'December 2015', '2007')
    description = 'Short Linear Motif Probability tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please cite: Davey, Haslam, Shields & Edwards (2010), Lecture Notes in Bioinformatics 6282: 50-61.']
    return rje.Info(program,version,last_edit,description,author,time.time(),copyyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        helpx = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if helpx > 0:
            print '\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time)))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show SLiMList commandline options?'): out.verbose(-1,4,text=rje_slimlist.__doc__)
            if rje.yesNo('Show SLiMCalc commandline options?'): out.verbose(-1,4,text=rje_slimcalc.__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()
            cmd_list += rje.inputCmds(out,cmd_list)
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print 'Major Problem with cmdHelp()'
#########################################################################################################################
def setupProgram(): ### Basic Setup of Program
    '''
    Basic setup of Program:
    - Reads sys.argv and augments if appropriate
    - Makes Info, Out and Log objects
    - Returns [info,out,log,cmd_list]
    '''
    try:
        ### Initial Command Setup & Info ###
        info = makeInfo()
        cmd_list = rje.getCmdList(sys.argv[1:],info=info)      ### Load defaults from program.ini
        ### Out object ###
        out = rje.Out(cmd_list=cmd_list)
        out.verbose(2,2,cmd_list,1)
        out.printIntro(info)
        ### Additional commands ###
        cmd_list = cmdHelp(info,out,cmd_list)
        ### Log ###
        log = rje.setLog(info=info,out=out,cmd_list=cmd_list)
        return [info,out,log,cmd_list]
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except:
        print 'Problem during initial setup.'
        raise
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMProb Class                                                                                          #
#########################################################################################################################
class SLiMProb(rje_slimcore.SLiMCore):     
    '''
    Short Linear Motif Regular Expression Search Tool Class. Author: Rich Edwards (2007). This module inherits the
    SLiMFinder class from SLiMFinder, which handles the masking of datasets and correcting for evolutionary relationships
    (if this desired).

    Str:str
    - AAFreq = Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    - AltDis = Alternative all by all distance matrix for relationships [None]
    - Background = Use observed support in background file for over-representation calculations [None]
    - BuildPath = Alternative path to look for existing intermediate files [SLiMProb/]
    - CaseMask = Mask Upper or Lower case [None]
    - CompMask = Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    - GablamDis = Alternative GABLAM results file [None]
    - Input = Original name (and path) of input file
    - ResDir = Redirect individual output files to specified directory [SLiMProb/]
    - ResFile = If FILE is given, will also produce a table of results in resfile [slimprob.csv]
    - RunID = Run ID for resfile (allows multiple runs on same data) [DATE:TIME]
    
    Bool:boolean
    - DisMask = Whether to mask ordered regions (see rje_disorder for options) [False]
    - Force = whether to force recreation of key files [False]
    - EFilter = Whether to use evolutionary filter [True]
    - LogMask = Whether to log the masking of individual sequences [True]
    - Masked = Whether dataset has been masked [False]
    - Masking = Master control switch to turn off all masking if False [True]
    - MaskM = Masks the N-terminal M (can be useful if termini=T) [False]
    - MaskFreq = Whether to mask input before any analysis, or after frequency calculations [True]
    - OccUPC = Whether to output the UPC ID number in the occurrence output file [False]
    - Pickle = Whether to save/use pickles [True]
    - SmearFreq = Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    - TarGZ = Whether to tar and zip dataset result files (UNIX only) [False]
    - Teiresias = Replace TEIRESIAS only, making *.out and *.mask.fas files [False]

    Int:numeric
    - MaxOcc = Filter out Motifs with more than maximum number of occurrences [0]
    - MaxSeq = Maximum number of sequences to process [500]
    - MaxSize = Maximum dataset size to process in AA (or NT) [1e5]
    - SaveSpace = Delete "unneccessary" files following run (see Manual for details) [0]
    - SeqOcc = Restrict to sequences with X+ occurrences (adjust for high frequency SLiMs) [1]

    Num:numeric
    - MST = MST corrected size for whole dataset
    - StartTime = Starting time in seconds (for output when using shared log file)
    - WallTime = Time in hours before program will terminate [1.0]

    List:list
    - AAMask = Masks list of AAs from all sequences (reduces alphabet) []
    - Alphabet = List of characters to include in search (e.g. AAs or NTs) 
    - Batch = List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    - FTMask = UniProt features to mask out [EM,DOMAIN,TRANSMEM]
    - Headers = Headers for main output table
    - IMask = UniProt features to inversely ("inclusively") mask [IM]
    - NoHits = List of sequences without any motif hits []
    - OccHeaders = Headers for main occurrence output table
    - UP = List of UP cluster tuples
     
    Dict:dictionary
    - AAFreq = AA frequency dictionary for each seq / UPC
    - MaskPos = Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    - MST = MST corrected size for UPC {UPC:MST}
    - P1 = Probabilities of 1+ occurrences {SLiM:{Seq/UPC:p}}

    Obj:RJE_Objects
    - Background = Background SLiMProb object containing data for background occurrences
    - SeqList = main SeqList object containing dataset to be searched
    - SlimList = MotifList object handling motif stats and filtering options
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.strlist = ['AAFreq','AltDis','BuildPath','CaseMask','CompMask','GablamDis','Input','ResDir','ResFile',
                         'RunID','Background']
        self.boollist = ['DisMask','Force','EFilter','LogMask','Masked','Masking','MaskM','MaskFreq','SmearFreq',
                        'TarGZ','Teiresias','ConsMask','Pickle','OccUPC']
        self.numlist = ['MST','StartTime','WallTime']
        self.intlist = ['MaxSeq','SaveSpace','HomCut','SeqOcc','MaxSize','MaxOcc']
        self.listlist = ['AAMask','Alphabet','Batch','FTMask','Headers','IMask','NewScore','OccFilter','OccHeaders','OccStats',
                         'UP','NoHits']
        self.dictlist = ['AAFreq','MST','NewScore','OccFilter','StatFilter']
        self.objlist = ['Background','SeqList','SlimList']
        ### Defaults ###
        self._setDefaults(str='None',bool=True,num=0.0,int=0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
        self.setStr({'BuildPath':rje.makePath('SLiMProb/'),'CompMask':'None',
                      'ResDir':rje.makePath('SLiMProb/'),'ResFile':'slimprob.csv'})
        self.setInt({'SeqOcc':1,'MaxSeq':0,'MaxSize':1e5,'Extras':2})
        self.setBool({'ConsMask':False,'EFilter':True,'MaskM':False,'DisMask':False,'Masking':True})
        self.dict['MaskPos'] = {}
        t = time.localtime(time.time())
        self.setStr({'Date': '%s%s%s-%s:%s' % (str(t[0])[-2:],rje.preZero(t[1],12),rje.preZero(t[2],31),rje.preZero(t[3],24),rje.preZero(t[4],60))})
        self.setStr({'RunID': self.getStr('Date')})
        ### Object setup ###
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        self.obj['SlimList'].loadMotifs()
        self.dict['Output']['motifs'] = self.obj['SlimList'].motifOut(pattern=True)
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self.coreCmd()
        for cmd in self.cmd_list:
            try:### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ### 
                self._cmdReadList(cmd,'file',['AAFreq','AltDis','GablamDis','ResFile','Background'])
                self._cmdReadList(cmd,'path',['ResDir'])
                self._cmdReadList(cmd,'str',['CaseMask','CompMask','RunID'])
                self._cmdReadList(cmd,'bool',['DisMask','Force','EFilter','LogMask','Masked','Masking','MaskM',
                                             'MaskFreq','SmearFreq','TarGZ','Teiresias','OccUPC'])
                self._cmdReadList(cmd,'int',['MaxSeq','MaxSize','SaveSpace','HomCut','SeqOcc','MaxOcc','Extras'])
                self._cmdReadList(cmd,'num',['MST','WallTime'])
                self._cmdRead(cmd,'int','MaxSize','maxaa')
                self._cmdReadList(cmd,'list',['Batch'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
#########################################################################################################################
    ### <2> ### Simple Stats Methods                                                                                    #
#########################################################################################################################
    def slimNum(self): return self.obj['SlimList'].slimNum()
    def slims(self): return self.obj['SlimList'].slims()
    def slimOccNum(self,slim): return slim.occNum()
    def slimSeqNum(self,slim): return slim.seqNum()
    def slimIC(self,slim): return slim.stat['IC'] 
#########################################################################################################################
    def slimUP(self,slim):  ### Returns UP Num SLiM.
        '''Returns UP Num SLiM.'''
        uplist = []
        for seq in slim.seqs():
            upc = self.getUP(seq)
            #if upc: self.deBug('%s >> %s (%d seq)' % (seq.shortName(),self.list['UP'].index(upc),len(upc)))
            #else: self.deBug('%s !!>> %s' % (seq.shortName(),self.list['UP']))
            if upc not in uplist: uplist.append(upc)
        return len(uplist)
#########################################################################################################################
    ### <3> ### General Run Methods                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main SLiMProb Run Method
        '''
        Main SLiMProb Run Method:
        1. Input:
            - Read sequences into SeqList
            - or - Identify appropriate Batch datasets and rerun each with batch=True
        2. Generate UPC:
            - Check for existing UPC load if found. 
            - or - Save sequences as fasta and mask sequences in SeqList
            -  Perform BLAST and generate UPC based on saved fasta.
            - Calculate AAFreq for each sequence and UPC.
        3. Perform SLiM Search.
        4. Output results and tidy files.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:### ~ [1] Batch/SeqIn input ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Webserver') and self.obj['SlimList'].motifNum() < 1:
                self.serverEnd('NoMotif',exit=False)
                return False
            #seqcmd = ['gnspacc=T','usecase=T','accnr=F','seqnr=F'] + self.cmd_list + ['autoload=T','query=None']
            #self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            self.setupSeqIn()
            if self.background(): self.printLog('#BACK','Background data from "%s" successful' % self.getStr('Background'))
            self.setupBasefile()
            ## ~ [1a] Batch Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch and self.seqNum() < 1:   # No sequences loaded - use batch mode
                batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                if not batchfiles: self.errorLog('No input files found!',printerror=False)
                else:
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        self.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.run(batch=True)
                        self.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.setBool({'Append':True})
                        bx += 1
                        self.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.getBool('Win32') and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return
            self.setupResults()                 
            self.backupOrCreateResFile()
            ## ~ [1b] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~~~ ##
            if 0 < self.getInt('MaxSeq') < self.seqNum():
                self.printLog('#MAX','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.iStr(self.seqNum()),rje.iStr(self.getInt('MaxSeq'))))
                self.serverEnd('MaxSeq',exit=False)
                return False
            self.setNum({'StartTime':time.time()})

            ### ~ [2] Read/Generate UPC, mask input & calculate/read AAFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.makeUPC():
                self.errorLog('Error during makeUPC(). Abandoning %s run' % self.dataset(),printerror=False)
                self.serverEnd('Crash','makeUPC()')
                return False
            pickled = self.searchPickleMe(load=self.getBool('Pickle'))  # Returns appropriate pickled SLiMProb Object, else None
            if pickled: self = pickled  ## Replace me with my pickle!
            else: self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
            if self.getBool('MaskFreq'): self.makeAAFreq()
            else: 
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                self.makeAAFreq()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
            self.adjustAATotals()
            if 0 < self.getInt('MaxSize') < self.dict['AAFreq']['Dataset']['Total']:
                self.printLog('#MAX','%s = %s %s > Max %s %s. Analysis terminated.' % (self.dataset(),rje.iStr(self.aaNum()),self.units(),rje.iStr(self.getInt('MaxSize')),self.units()))
                self.serverEnd('MaxSize',exit=False)
                return False

            ### ~ [3] Perform basic search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                self.searchDB()     #!# Or read in results?! #!#
                self.searchPickleMe()     # Generates pickling for speedy re-running
                self.wallTime()

            ### ~ [4] Generate extra stats, filter and output etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Setup output & filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['SlimList'].setupFilters(occheaders=self.list['OccHeaders'])
            #x#self.dict['OccFilter'] = rje_scoring.setupStatFilter(self,self.list['OccHeaders'],self.list['OccFilter'])
            ## ~ [4b] Calculate extra occurrence attributes & output occurrences ~~~~~~~~~~~~~~~~~~ ##
            self.calculateOccAttributes()
            self.processSeqOccs()       # Converts for output, filtering as desired
            ## ~ [4c] Combine occurrences and generate summary, with statistics ~~~~~~~~~~~~~~~~~~~ ##
            self.combMotifOccStats()
            if self.getInt('Extras'): self.extraOutput()   # MotifList Outputs
            self.tarZipSaveSpace()      # Tarring, Zipping and Saving Space 

            ### ~ [5] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            resfilestr = string.split(self.getStr('SummaryFile'),'.')
            resfilestr[-1] = '(occ.)' + resfilestr[-1]
            resfilestr = string.join(resfilestr,'.')
            self.printLog('#RES','SLiMProb results output to %s and %s.*' % (resfilestr,self.seqBaseFile()))
            targz = '%s.tgz' % self.runBase()
            if self.getBool('TarGZ') and not self.getBool('Win32') and rje.exists(targz):
                self.printLog('#TGZ','Additional dataset results tarred to %s' % targz)
            if self.getBool('Win32') and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.getNum('WallTime') <= 0 or (time.time() - self.getNum('StartTime')) < (self.getNum('WallTime')*3600): raise
            return False # Walltime reached
        except:
            self.errorLog('Error in %s.run()' % self.prog(),printerror=True,quitchoice=False)
            self.serverEnd(endcause='Crash',details='main run')
            return False
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMProb object for new batch run
        '''Returns SLiMFinder object for new batch run.'''
        return SLiMProb(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.getBool('Append'),'basefile=None'])
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def background(self):   ### Sets up and/or returns Background SLiMProb object
        '''Sets up and/or returns Background SLiMProb object.'''
        try:### ~ [0] ~ Simple case first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['Background']: return self.obj['Background']
            if not rje.checkForFile(self.getStr('Background')): return None

            ### ~ [1] ~ Setup Background SLiMProb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bg = SLiMProb(self.log,self.cmd_list)
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','seqin=%s' % self.getStr('Background')]
            bg.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            bg.setupBasefile()
            ## ~ [1a] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~~~ ##
            if 0 < self.getInt('MaxSeq') < bg.seqNum():
                self.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (bg.dataset(),rje.iStr(bg.obj['SeqList'].seqNum()),rje.iStr(self.getInt('MaxSeq'))))
                raise ValueError
            self.setNum({'StartTime':time.time()})

            ### ~ [2] Read/Generate UPC, mask input & calculate/read AAFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not bg.makeUPC(): self.errorLog('Error during makeUPC(). Abandoning %s run' % bg.dataset(),printerror=False); raise ValueError
            pickled = bg.searchPickleMe(load=True)  # Returns appropriate pickled SLiMFinder Object, else None
            if pickled: bg = pickled  ## Replace me with my pickle!
            else: bg.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']

            ### ~ [3] Perform basic search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                bg.searchDB()     #!# Or read in results?! #!#
                bg.searchPickleMe()     # Generates pickling for speedy re-running
            self.obj['Background'] = bg
            return bg

        except ValueError: self.setStr({'Background':'None'}); return None
        except: self.errorLog('SLiMProb.background error!'); return None
#########################################################################################################################
    def OLDsetupBasefile(self):    ### Sets up self.info['Basefile'] and self.info['Input']
        '''Sets up self.info['Basefile'].'''
        ### ~ [0] ~ Store original input filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.getStrLC('Input'): self.setStr({'Input': self.obj['SeqList'].info['Name']})
        ### ~ [1] ~ Basefile and ResDir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        rje.mkDir(self,self.getStr('ResDir'))
        if self.dev():
            if not self.baseFile(): self.baseFile(self.prog().lower())
            return
        #!# Delete this method once dev tested - in SLiMCore #!#
        if not self.getStrLC('Basefile'):
            self.baseFile(rje.baseFile(self.obj['SeqList'].info['Name']))
            if self.getStrLC('ResDir'):
                rje.mkDir(self,self.getStr('ResDir'))
                slimbase = os.path.split(rje.baseFile(self.obj['SlimList'].info['Name']))[1]
                if slimbase.lower() in ['','none']: self.baseFile(self.getStr('ResDir') + '%s' % (self.dataset()))
                else: self.baseFile(self.getStr('ResDir') + '%s.%s' % (self.dataset(),slimbase))
#########################################################################################################################
    def buildText(self): return rje.baseFile(self.obj['SlimList'].name(),strip_path=True)
#########################################################################################################################
#CORE# def maskInput(self):    ### Masks input sequences, replacing masked regions with Xs
#CORE# def makeUPC(self):  ### Generates UP Clusters from self.obj['SeqList'] using BLAST
#CORE# def readUPC(self):  ### Generates UP Clusters from self.obj['SeqList'] using BLAST
#CORE# makeMST(self,gablam=None):   ### Makes UPC dictionary from GABLAM
#########################################################################################################################
    ### <5> ### SearchDB Methods                                                                                        #
#########################################################################################################################
    def searchDB(self):      ### Main method for searching database with motifs. 
        '''Main method for searching database with motifs. This just does basic matches. Stats etc. dealt with later.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.slimNum(): return self.errorLog('No SLiMs for SLiMProb.searchDB()',printerror=False)
            (ox,sx,stot) = (0,0.0,self.seqNum()*self.slimNum())    # Counter for Log output
            stxt = 'Searching %s sequences with %d SLiMs: ' % (rje.integerString(self.seqNum()),self.slimNum())
            for seq in self.seqs():
                ## ~ [1a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.progLog('\r#SEARCH','%s %.1f%% (%s occ)' % (stxt,(sx/stot),rje.integerString(ox)))
                seq.deGap()
                ## ~ [1b] Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for slim in self.slims():
                    sx += 100.0
                    #!# Note that MegaSLiMProb would have to replace this search function here.
                    ox += len(slim.searchSequence(seq))    # {Pos,Variant,Match,ID,MisMatch}
                    self.progLog('\r#SEARCH','%s %.1f%% (%s occ)' % (stxt,(sx/stot),rje.integerString(ox)))
                    self.wallTime()

            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
            self.setBool({'Searched':True})
            self.printLog('\n#SEARCH','%s sequences searched for %s motifs: %s occ.' % (rje.integerString(self.seqNum()),rje.integerString(self.slimNum()),rje.integerString(ox)))
        except:
            self.errorLog('Error in SLiMProb.searchDB()')
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def calculateOccAttributes(self):   ### Executes rje_slimcalc calculations via rje_slimlist object
        '''Executes rje_slimcalc calculations via rje_slimlist object.'''
        self.obj['SlimList'].calculateOccAttributes(wallobj=self)   ### Executes rje_slimcalc calculations via rje_slimlist object        
#########################################################################################################################
    def processSeqOccs(self):   ### Processes Occurrences after search - outputs to results file
        '''Processes Occurrences after search - output to results file.'''
        try:### ~ [1] Process and Output hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (sb,sx,stot,ox) = (0.0,0.0,self.slimNum(),0)
            for slim in self.slims():
                if slim.dict['Occ']: sj = 100.0 / len(slim.dict['Occ'])
                for seq in slim.dict['Occ'].keys():
                    self.progLog('\r#OCC','Processing occurrences: %.1f%%' % ((sb+sx)/stot))
                    sx += sj
                    ## ~ [1a] Update dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #i# [Dataset,RunID,Masking,Motif,RunID,Seq,Start_Pos,End_Pos,Prot_Len,Match,Variant,MisMatch,Desc]
                    results = {}    # Dictionary of results dictionaries to output
                    occlist = slim.dict['Occ'][seq][0:]
                    for occ in occlist[0:]:
                        ox += 1
                        occ['Dataset'] = self.dataset()
                        occ['Motif'] = slim.info['Name']
                        occ['RunID'] = self.getStr('RunID')
                        occ['Masking'] = self.maskText()
                        occ['Seq'] = seq.shortName()
                        occ['Start_Pos'] = max(1,occ['Pos'])
                        occ['End_Pos'] = occ['Start_Pos'] + len(occ['Match']) - 1
                        occ['Prot_Len'] = seq.aaLen()
                        occ['Desc'] = seq.info['Description']
                        occ['Pattern'] = slim.pattern()
                        if self.getBool('OccUPC'): occ['UPC'] = self.list['UP'].index(self.getUP(seq)) + 1
                        for occstat in ['Cons','SA','Hyd','Fold','IUP','Comp']:
                            if occstat in occ: occ[occstat] = rje_slim.expectString(occ[occstat])
                        results['%s-%s-%s' % (rje.preZero(occ['Start_Pos'],occ['Prot_Len']),rje.preZero(occ['End_Pos'],occ['Prot_Len']),occ['Variant'])] = occ
                    ## ~ [1b] Filtering & Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    results = rje_scoring.statFilter(self,results,self.dict['OccFilter'])  #!# Add restrict/exclude #!#
                    for occ in occlist[0:]:
                        okey = '%s-%s-%s' % (rje.preZero(occ['Start_Pos'],occ['Prot_Len']),rje.preZero(occ['End_Pos'],occ['Prot_Len']),occ['Variant'])
                        if okey in results: rje.delimitedFileOutput(self,self.getStr('ResFile'),self.resHead('OccHeaders'),datadict=results[okey])
                        else: slim.dict['Occ'][seq].remove(occ)
                    ## ~ [1c] SeqOcc Reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if self.getInt('SeqOcc') > 1 and self.getInt('SeqOcc') > len(slim.dict['Occ'][seq]): slim.dict['Occ'].pop(seq)
                    self.wallTime()
                (sx,sb) = (0.0,sb+100.0)
            self.printLog('\r#OCC','Processing %s occurrences complete.' % rje.integerString(ox))
        except: self.errorLog('Error in processSeqOccs()')
#########################################################################################################################
    def combMotifOccStats(self):    ### Combined occurrence stats for each Motif and summary output
        '''Combined occurrence stats for each Motif and summary output.'''
        try:### ~ [1] Generate missing values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.background(): ptxt = 'Calculating SLiM Probabilities using Background frequencies'
            else: ptxt = 'Calculating SLiM Probabilities'
            (sx,stot) = (0.0,len(self.slims()))
            for slim in self.slims():
                self.progLog('\r#PROB','%s: %.1f%%' % (ptxt,sx/stot)); sx+=100;
                self.slimProb(slim)
                self.wallTime()
            self.printLog('\r#PROB','%s complete.' % ptxt)
            self.obj['SlimList'].combMotifOccStats()

            ### ~ [2] Filtering & Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] General output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sec = int(time.time() - self.getNum('StartTime') + 0.5)
            (hour,min) = (0,0)
            while sec >= 60: (min,sec) = (min+1,sec-60)
            while min >= 60: (hour,min) = (hour+1,min-60)
            self.setStr({'RunTime':'%s:%s:%s' % (rje.preZero(hour,24),rje.preZero(min,60),rje.preZero(sec,60))})
            totalaa = 0     # Total number of AA in dataset
            for seq in self.seqs():
                if self.getBool('Masked'): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                totalaa += seq.nonX()
            general = {'Dataset':self.dataset(),'RunID':self.getStr('RunID'),'Masking':self.maskText(),
                       'RunTime':self.getStr('RunTime'),'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa}
            ## ~ [2b] SLiM Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for slim in self.slims():
                datadict = rje.combineDict({'Motif':slim.info['Name'],'Pattern':slim.pattern()},general)
                datadict = rje.combineDict(datadict,slim.stat)
                for dkey in datadict.keys():
                    if string.split(dkey,'_')[0] in ['E','p','pUnd']: datadict[dkey] = rje_slim.expectString(datadict[dkey])
                    elif string.split(dkey,'_')[-1] in ['mean']: datadict[dkey] = rje_slim.expectString(datadict[dkey])
                rje.delimitedFileOutput(self,self.getStr('SummaryFile'),self.resHead(),datadict=datadict)
            self.printLog('#OUT','Summary data for %s saved to %s' % (self.dataset(),self.getStr('SummaryFile')))
        except: self.errorLog('Error in combMotifOccStats()')
#########################################################################################################################
    ### <6> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
#CORE# def makeAAFreq(self):   ### Makes an initial AAFreq dictionary containing AA counts only (including Xs)
#########################################################################################################################
    def adjustAATotals(self):   ### Adjusts AA Totals following masking     #!# Changed from Core #!#
        '''
        Adjusts AA Totals following masking and makes appropriate frequencies. If maskfreq=T then the amino acid counts
        will be exactly as they were. If maskfreq=F, however, frequencies will need to be adjust for the new number of
        masked and non-masked amino acids.
        '''
        try:
            ### Simple PreMasking Procedure ###
            if self.getBool('MaskFreq'):
                for upc in self.list['UP']:
                    if self.dict['AAFreq'][upc].has_key('X'): self.dict['AAFreq'][upc].pop('X') # Ignore Xs
                    self.dict['AAFreq'][upc].pop('Total')  #!# Total remade by dictFreq #!#
                    self.dict['AAFreq'][upc]['^'] = 0
                    self.dict['AAFreq'][upc]['$'] = 0
                    rje.dictFreq(self.dict['AAFreq'][upc])
                    ## Make MST adjustments for UPC ##
                    self.dict['AAFreq'][upc]['Total'] = int(0.5+(self.dict['AAFreq'][upc]['Total']*self.dict['MST'][upc]))
                    if self.dict['AAFreq'][upc]['Total']:
                        self.dict['AAFreq'][upc]['^'] = self.dict['MST'][upc] / float(self.dict['AAFreq'][upc]['Total'])
                        self.dict['AAFreq'][upc]['$'] = self.dict['MST'][upc] / float(self.dict['AAFreq'][upc]['Total'])
                    else:
                        self.dict['AAFreq'][upc]['^'] = 0.5
                        self.dict['AAFreq'][upc]['$'] = 0.5
                if self.getBool('SmearFreq'): self.smearAAFreq()
                return

            ### More complicated PostMasking Procedure ###
            (prex,postx) = (0,0)
            for upc in self.list['UP']:
                x = 0
                if self.dict['AAFreq'][upc].has_key('X'): x = self.dict['AAFreq'][upc].pop('X') # Ignore Xs
                totalaa = self.dict['AAFreq'][upc].pop('Total')
                preaa = totalaa - x   # Want to calculate new total
                prex += preaa
                nonx = 0.0
                for seq in upc: nonx += seq.nonX()
                ## Termini ##
                self.dict['AAFreq'][upc]['^'] = 0
                self.dict['AAFreq'][upc]['$'] = 0
                rje.dictFreq(self.dict['AAFreq'][upc])
                ## Make MST adjustments for UPC ##
                self.dict['AAFreq'][upc]['Total'] = int(0.5+(nonx*self.dict['MST'][upc])) 
                postx += self.dict['AAFreq'][upc]['Total']
                ## Termini ##
                if nonx:
                    self.dict['AAFreq'][upc]['^'] = self.dict['MST'][upc] / nonx
                    self.dict['AAFreq'][upc]['$'] = self.dict['MST'][upc] / nonx
                else:
                    self.dict['AAFreq'][upc]['^'] = 0.5
                    self.dict['AAFreq'][upc]['$'] = 0.5
                    
            ### Finish ###
            if self.getStrLC('AAFreq'): prex = self.dict['AAFreq']['Dataset']['Total']
            if prex != postx: self.printLog('#ADJ','Effective dataset size reduced from %s %s to %s %s' % (rje.integerString(prex),self.units(),rje.integerString(postx),self.units()))
            if self.getBool('SmearFreq'): self.smearAAFreq()
        except:
            self.errorLog('Problem during SLiMProb.adjustAATotals()')
            raise
#########################################################################################################################
#CORE# def smearAAFreq(self):  ### Equalises AAFreq across UPC. Leaves Totals unchanged.
#########################################################################################################################
    def slimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM. Modified from slimcore for SLiMProb.'''
        try:### ~ [1] Setup SLiMProb Motif stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slim.stat['N_Seq'] = slim.seqNum()      # also 'E_Occ','E_Seq','E_UPC','p_Occ','p_Seq','p_UPC'
            slim.stat['N_UPC'] = self.slimUP(slim)
            slim.stat['N_Occ'] = slim.occNum()
            slim.stat['E_Occ'] = 0.0
            if self.background(): return self.slimProbBG(slim)
            ## ~ [1a] Setup pattern and variable-lenght multiplier ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            poslist = []    # List of AA positions in SLiM
            wildlist = []   # List of wildcard lengths in SLiM
            wild = False    # Whether next part is a wildcard length
            mult = 1        # Variable-length multiplier
            for part in string.split(slim.slim(),'-'):      # Split SLiM code in components
                ## Update lists ##
                if wild: wildlist.append(part)
                else: poslist.append(part)
                ## Calculate multiplier ##
                if wild:
                    (minx,maxx) = (int(part[0]),0)
                    for x in part:
                        minx = min(minx,int(x))
                        maxx = max(maxx,int(x))
                    mult *= (int(maxx) - int(minx) + 1)
                wild = not wild

            ### ~ [2] Calculate prob of 1+ occ for each UPC/Seq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            p1 = {'UPC':{},'Seq':{}}         # Dictionary of {upc:chance of 1+ occ in upc/seq}
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.list['UP']:
                ## Setup  parameters for binomial ##
                N = (self.dict['AAFreq'][upc]['Total'] - ((slim.slimMinLen()-1) * self.dict['MST'][upc])) * mult       # Each length variant is effectively another position the SLiM could occur
                N = max(0,N)
                p = 1.0                                 # Probability of SLiM at each position
                k = 1                                   # Number of successful trials (occurrences)
                if self.getInt('SeqOcc') > 1: k = self.getInt('SeqOcc')
                ## Calculate p and N from AAFreq ##
                for pos in poslist:     # AA position
                    posfreq = 0.0
                    nofm = rje.matchExp('<(\D+):(\d+):(\d+)>',pos)
                    nofmorb = rje.matchExp('<(\D+):(\d+):(\d+):(\D+)>',pos)
                    if pos.find('^') > 0: raise ValueError
                    if nofm:        # Special N of M position
                        combo = rje.binComb(int(nofm[2]),int(nofm[1]),int(nofm[1]))
                        for aa in nofm[0]: posfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        posfreq = min(1.0,posfreq) ** int(nofm[1])  # Cannot have more that a 100% chance of matching a position
                        posfreq = min(1.0, posfreq * combo)         # Cannot have more that a 100% chance of matching a position
                    elif nofmorb:   # Special N of M or B position
                        combo = rje.binComb(int(nofmorb[2]),int(nofmorb[1]),int(nofmorb[1]))
                        mfreq = 0.0; bfreq = 0.0
                        for aa in nofmorb[0]: mfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        for aa in nofmorb[3]: bfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        posfreq = (min(1.0,mfreq) ** int(nofmorb[1])) * (min(1.0,bfreq) ** (int(nofmorb[2]) - int(nofmorb[1]))) # Cannot have more that a 100% chance of matching a position
                        posfreq = min(1.0, posfreq * combo)         # Cannot have more that a 100% chance of matching a position
                    #!# Add (XY|AB) formats too #!#                    
                    else:   # Normal ambiguity
                        if pos.count('X') > 0: posfreq = 1.0
                        else:
                            for aa in pos: posfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        posfreq = min(1.0,posfreq) # Cannot have more that a 100% chance of matching a position
                    p *= posfreq   
                if p > 1: p = 1.0           # Cannot in reality have p > 1!
                ## Calculate binomial ##
                try: p1['UPC'][upc] = rje.binomial(k,N,p,usepoisson=False,callobj=self)      # logBinomial gets maths range errors sometimes!
                except: p1['UPC'][upc] = 0.0; print k, N , p
                #if slim.pattern() == 'KLY': open('kly.tmp','a').write('%s::SS| k = %d; p = %s; N = %d; p1+ = %s\n' % (upc[0].shortName(),k,p,N,p1['UPC'][upc]))                        
                ## Sequence-specific stats ##
                for seq in upc:
                    N = (seq.nonX() - (slim.slimMinLen()-1)) * mult 
                    N = max(0,N)
                    slim.stat['E_Occ'] += (N * p)
                    if self.getInt('SeqOcc') > 1: k = self.getInt('SeqOcc')
                    ## Calculate binomial ##
                    try: p1['Seq'][seq] = rje.binomial(k,N,p,usepoisson=False,callobj=self)
                    except: print '!', k , N, p
            ## Extra verbosity. Remove at some point? ##
            self.verbose(2,3,'%s: %s' % (slim.pattern(),p1),1)

            ### ~ [3] Calculate overall probability of observed support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#print '>>>\n', slim.stat
            ## ~ [3a] All observed occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if slim.stat['N_Occ'] <= 0: slim.stat['p_Occ'] = 1.0
            elif slim.stat['E_Occ'] > 0.0: slim.stat['p_Occ'] = rje.logPoisson(slim.stat['N_Occ'],slim.stat['E_Occ'],callobj=self)     #!# logPoisson() #!#
            else: slim.stat['p_Occ'] = 0.0
            ## ~ [3b] UPC occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slim.stat['E_UPC'] = sum(p1['UPC'].values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimUP(slim), self.UPNum(), slim.stat['E_UPC']/self.UPNum()) # Use mean p1+
            if k <= 0: slim.stat['p_UPC'] = 1.0
            else:
                try: slim.stat['p_UPC'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
                except: print '!!', k , N, p
            ## ~ [3c] Seq occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slim.stat['E_Seq'] = sum(p1['Seq'].values())    # Expected number of observed UPCs
            (k,n,p) = (slim.seqNum(), self.seqNum(), slim.stat['E_Seq']/self.seqNum()) # Use mean p1+
            if k <= 0: slim.stat['p_Seq'] = 1.0
            else:
                try: slim.stat['p_Seq'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
                except: print '!!!', k , N, p
            ## ~ [3d] Under-representation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for lvl in ['Occ','UPC','Seq']:
                k = slim.stat['N_%s' % lvl]
                expected = slim.stat['E_%s' % lvl]
                po = rje.logPoisson(k,expected,exact=True,callobj=self)       # Exact probability of observed occ
                #self.bugPrint('k=%s, exp=%s, pk=%s' % (k,expected,po))
                #self.bugPrint('pOver=%s' % slim.stat['p_%s' % lvl])
                slim.stat['pUnd_%s' % lvl] = (1 - slim.stat['p_%s' % lvl]) + po
                #self.bugPrint('(1 - %s) + %s = %s !' % (slim.stat['p_%s' % lvl],po,slim.stat['pUnd_%s' % lvl]))
                slim.stat['pUnd_%s' % lvl] = 0.0
                for x in range(0,k+1): slim.stat['pUnd_%s' % lvl] += rje.logPoisson(x,expected,exact=True,callobj=self)
                #self.deBug('Recalc => %s' % rje_slim.expectString(slim.stat['pUnd_%s' % lvl]))
                slim.stat['pUnd_%s' % lvl] = min(1.0,slim.stat['pUnd_%s' % lvl])    # Cannot exceed 1 (rounding error?)
            #self.deBug('%s >> %s' % (slim.pattern(),slim.stat))
        except:
            self.errorLog('Error with slimProb(%s)' % slim.info['Sequence'])
            slim.stat['p_Occ'] = slim.stat['p_Seq'] = slim.stat['p_UPC'] = 1.0
            #self.deBug(slim)
#########################################################################################################################
    def slimProbBG(self,slim):  ### Calculate Probabilities for given SLiM using Background Support
        '''Calculate Probabilities for given SLiM using Background Support.'''
        try:### ~ [1] Setup SLiMProb Motif stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bg = self.background()
            bgslim = bg.obj['SlimList'].mapPattern(slim.info['Sequence'],update=False)

            ### ~ [2] Calculate prob of 1+ occ for each UPC/Seq~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            p1 = {'UPC':{},'Seq':{}}         # Dictionary of {upc:chance of 1+ occ in upc/seq}
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            bgseqx = bgslim.seqNum()      
            bgupcx = bg.slimUP(bgslim)
            bgoccx = bgslim.occNum()
            for upc in bg.list['UP']:
                p1['UPC'][upc] = bgupcx / float(bg.UPNum())
                ## Sequence-specific stats ##
                for seq in upc:
                    slim.stat['E_Occ'] += (bgoccx / float(bg.seqNum()))
                    p1['Seq'][seq] = bgseqx / float(bg.seqNum())
            ## Extra verbosity. Remove at some point? ##
            self.verbose(2,3,'%s: %s' % (slim.pattern(),p1),1)

            ### ~ [3] Calculate overall probability of observed support ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #x#print '>>>\n', slim.stat
            ## ~ [3a] All observed occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if slim.stat['N_Occ'] <= 0: slim.stat['p_Occ'] = 1.0
            elif slim.stat['E_Occ'] > 0.0: slim.stat['p_Occ'] = rje.logPoisson(slim.stat['N_Occ'],slim.stat['E_Occ'],callobj=self)     #!# logPoisson() #!#
            else: slim.stat['p_Occ'] = 0.0
            ## ~ [3b] UPC occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slim.stat['E_UPC'] = sum(p1['UPC'].values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimUP(slim), self.UPNum(), slim.stat['E_UPC']/self.UPNum()) # Use mean p1+
            if k <= 0: slim.stat['p_UPC'] = 1.0
            else: slim.stat['p_UPC'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            ## ~ [3c] Seq occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slim.stat['E_Seq'] = sum(p1['Seq'].values())    # Expected number of observed UPCs
            (k,n,p) = (slim.seqNum(), self.seqNum(), slim.stat['E_Seq']/self.seqNum()) # Use mean p1+
            if k <= 0: slim.stat['p_Seq'] = 1.0
            else: slim.stat['p_Seq'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            #self.deBug('%s >> %s' % (slim.pattern(),slim.stat))
            ## ~ [3d] Under-representation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for lvl in ['Occ','UPC','Seq']:
                k = slim.stat['N_%s' % lvl]
                expected = slim.stat['E_%s' % lvl]
                po = rje.logPoisson(k,expected,exact=True,callobj=self)       # Exact probability of observed occ
                #self.bugPrint('k=%s, exp=%s, pk=%s' % (k,expected,po))
                #self.bugPrint('pOver=%s' % slim.stat['p_%s' % lvl])
                slim.stat['pUnd_%s' % lvl] = (1 - slim.stat['p_%s' % lvl]) + po
                #self.bugPrint('(1 - %s) + %s = %s !' % (slim.stat['p_%s' % lvl],po,slim.stat['pUnd_%s' % lvl]))
                slim.stat['pUnd_%s' % lvl] = 0.0
                for x in range(0,k+1): slim.stat['pUnd_%s' % lvl] += rje.logPoisson(x,expected,exact=True,callobj=self)
                #self.deBug('Recalc => %s' % rje_slim.expectString(slim.stat['pUnd_%s' % lvl]))
                slim.stat['pUnd_%s' % lvl] = min(1.0,slim.stat['pUnd_%s' % lvl])    # Cannot exceed 1 (rounding error?)
        except:
            self.errorLog('Error with slimProb(%s)' % slim.info['Sequence'])
            slim.stat['p_Occ'] = slim.stat['p_Seq'] = slim.stat['p_UPC'] = 1.0
            #self.deBug(slim)
#########################################################################################################################
    ### <8> ### SLiMBuild Pickle Methods                                                                                #
#########################################################################################################################
    def processPickle(self,newme):  ### Changes attributes accordingly
        '''Changes attributes accordingly. Replace this method in subclasses.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Check other pertinent attributes - masking and additional filtering ##
            if self.obj['SlimList'].nameList() != newme.obj['SlimList'].nameList():
                self.printLog('#PICKLE','Motif attributes changed since pickle: will not use'); return None
            ## Note that MustHave and OccFilter filtering currently occur *after* SLiMBuild only ##
            changes = []
            for var in ['CompMask','CaseMask','MotifMask']:         # Info
                if self.getStr(var) != newme.getStr(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
            for var in ['Masking','DisMask','MaskM','ConsMask']:   # Opt
                if self.getBool(var) != newme.getBool(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
            for var in ['FTMask','IMask']:  #x#,'Equiv']:      # List
                slist = self.list[var][0:]
                nlist = newme.list[var][0:]
                slist.sort()
                nlist.sort()
                if slist != nlist:
                    changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
            ## Recreate or use pickle but add new commands ##
            if changes and (self.i() < 0 or rje.yesNo('%d SLiMBuild parameter mismatches with pickle. Create new pickle?' % len(changes))):
                self.printLog('#PICKLE','Parameters changed. Making new pickle.')
                return None
            newme.cmd_list = self.cmd_list
            newme.setStr(self.str)
            newme.setInt(self.int)
            newme.setNum(self.num)
            self.setBool({'Masked':newme.getBool('Masked'),'DNA':newme.getBool('DNA')})
            newme.setBool(self.bool)
            newme.setLog(self.log)
            self.list['UP'] = newme.list['UP']
            self.dict['AAFreq'] = newme.dict['AAFreq']
            self.dict['MST'] = newme.dict['MST']
            newme.list = self.list
            newme.dict = self.dict
            return newme
        except: self.errorLog('Problem during %s.processPickle()' % self); return None
#########################################################################################################################
    ### <9> ### Results Output Methods                                                                                  #
#########################################################################################################################
    def setupResults(self):     ### Sets up Main Results File as well as StatFilters etc.
        '''Sets up Main Results File as well as StatFilters etc.'''
        try:### ~ [1] Setup Initial Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slist = self.obj['SlimList']
            self.setupBasefile()
            self.list['Headers'] = ['Motif','Pattern','IC']
            for lvl in ['Occ','Seq','UPC']:
                for s in ['N','E','p','pUnd']: self.list['Headers'].append('%s_%s' % (s,lvl))
            self.list['OccHeaders'] = ['Motif','Seq','Start_Pos','End_Pos','Prot_Len','Pattern','Match','Variant','MisMatch','Desc']
            if self.getBool('OccUPC'): self.list['OccHeaders'].append('UPC')
            if slist.opt['Peptides']: self.list['OccHeaders'] += ['PepSeq','PepDesign']

            ### ~ [2] Special Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for h in slist.obj['SLiMCalc'].list['Headers']:
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            for h in slist.obj['SLiMCalc'].list['OccHeaders']:
                if h not in self.list['OccHeaders']: self.list['OccHeaders'].append(h)

            ### ~ [3] Custom Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            (newheads,self.list['NewScore'],self.dict['NewScore']) = rje_scoring.setupCustomScores(self,self.list['OccHeaders']+self.list['OccStats'],self.list['NewScore'],self.dict['NewScore'])
            for new in newheads:
                if new not in self.list['OccHeaders'] and new not in self.list['OccStats']:
                    if rje.formula(self,self.dict['NewScore'][new],varlist=self.list['OccHeaders'],check=False,calculate=False):
                        self.list['OccHeaders'].append(new)
                    else:
                        self.list['OccStats'].append(new)   #!# Made of OccStats, not means
                        self.list['OccHeaders'].append(self.list['OccStats'][-1])
                        self.list['Headers'].append('%s_mean' % self.list['OccStats'][-1])

            ### ~ [4] ResFile(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('ResFile') and self.getStr('ResFile').find('.') < 0: self.setStr({'ResFile':self.getStr('ResFile') + '.csv'})
            self.setStr({'SummaryFile': self.getStr('ResFile')})
            spex = os.path.splitext(self.getStr('ResFile'))
            if spex[0].endswith('.occ'):    # Old-style resfile naming of occ file
                self.setStr({'SummaryFile': os.path.splitext(spex[0])[0] + spex[1]})
            else: self.setStr({'ResFile': spex[0] + '.occ' + spex[1]})
            self.dict['Output']['main'] = self.getStr('SummaryFile')
            self.dict['Output']['occ'] = self.getStr('ResFile')
        except:
            self.errorLog('Problem with SLiMProb.setupResults()')
            raise
#########################################################################################################################
    def backupOrCreateResFile(self):    ### Backups up and/or creates main results file
        '''Backups up and/or creates main results file.'''
        if self.getStrLC('ResFile'):
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.getStr('ResFile'),write=True))
            rje.delimitedFileOutput(self,self.getStr('ResFile'),self.resHead('OccHeaders'),delimit,rje_backup=True)
        if self.getStrLC('SummaryFile'):
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.getStr('SummaryFile'),write=True))
            rje.delimitedFileOutput(self,self.getStr('SummaryFile'),self.resHead(),delimit,rje_backup=True)
#########################################################################################################################
    def resHead(self,htype='Headers'):  ### Returns main Output headers
        '''Returns main Output headers.'''
        if htype == 'Headers': heads = ['Dataset','RunID','Masking','RunTime','SeqNum','UPNum','AANum']
        else: heads = ['Dataset','RunID','Masking']
        return heads + self.list[htype]
#########################################################################################################################
    def extraOutput(self):  ### Method controlling additional outputs (primarily MotifList alignments)
        '''
        Method controlling additional outputs (primarily MotifList alignments):
        - Full protein alignments (in subdirectory) with orthlogues (no masking)
        - Protein Alignments with Motifs and masking marked
        - Motif alignments with and without masking
        '''
        try:
            ### ~ [1] Masked Sequence output (fasta) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#EXTRA','Generating additional outputs...',log=False)
            if self.getBool('Masked'):
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq']
                self.obj['SeqList'].saveFasta(seqfile='%s.%s.masked.fas' % (self.seqBaseFile(),self.maskText()))
                self.dict['Output']['masked'] = '%s.%s.masked.fas' % (self.seqBaseFile(),self.maskText())
            if self.getInt('Extra') < 2: return
            
            ### ~ [2] Full Protein Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #alndir = rje.makePath('%s_ORTHALN' % self.info['Basefile'])     #!# Modify to match SLiMFinder #!#
            #self.obj['SlimList'].proteinAlignments(alndir=alndir)

            ### ~ [3] Protein Alignments with Motifs and Masking marked ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alncmd = ['seqin=None','accnr=F','seqnr=F','autofilter=F','align=F','gnspacc=F','unkspec=T'] 
            paln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
            paln.info['Name'] = '%s.mapping.fas' % self.runBase()
            self.dict['Output']['mapping'] = paln.info['Name']
            if self.getBool('Webserver'): waln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd+['replacechar=F'])
            ## ~ [3a] Generate Alignment Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#MAP','Generating motif mapping alignments')
            motif_occ = self.obj['SlimList'].motifOcc(byseq=True)
            for seq in self.seqs():
                if self.getBool('Masked'): seq.info['Sequence'] = seq.info['PreMask']   #x# seq.info['MaskSeq']
                occlist = []
                if motif_occ.has_key(seq):
                    for Motif in motif_occ[seq].keys():
                        occlist += motif_occ[seq][Motif]
                if not occlist: continue
                #self.bugPrint(seq); self.bugPrint(occlist)
                #self.deBug(self.obj['SlimList'].obj['SLiMCalc'].info)
                #self.deBug(self.obj['SlimList'].obj['SLiMCalc'].opt)
                #self.deBug(self.obj['SlimList'].obj['SLiMCalc'].stat)
                saln = self.obj['SlimList'].obj['SLiMCalc'].singleProteinAlignment(seq,occlist,usegopher=False,savefasta=False,wintuple=30)
                saln.seq[0].info['Name'] = '%s-%s Motifs' % (seq.shortName(),self.dataset())
                #self.deBug('%d:%s' % (saln.seqNum(),saln.accList()))
                if self.getBool('Masked'):
                    alnseq = saln.seq[1].info['Sequence']
                    #x# saln.seq[1].info['Sequence'] = rje_sequence.mapGaps(seq.info['PreMask'],alnseq,self)
                    saln._addSeq('%s-masked' % seq.shortName(),rje_sequence.mapGaps(seq.info['MaskSeq'],alnseq,self))
                    saln.seq = saln.seq[:1] + saln.seq[-1:] + saln.seq[1:-1]
                paln.seq += saln.seq[0:]
                ## Special webserver output ##
                if self.getBool('Webserver') and occlist:
                    #self.deBug(occlist)
                    waln.info['Name'] = '%s.webserver.%s.fas' % (self.getStr('Basefile'),seq.shortName())
                    waln.seq = []
                    for wseq in saln.seq:
                        makeme = []
                        for i in range(len(occlist)):
                            Occ = occlist[i]
                            if i > occlist.index(Occ): continue     # Stupid duplicates!
                            (start,end,v,lshift) = saln.list['WinTuple'][i]
                            #X#v = Occ['Motif'].pattern()
                            makeme.append(wseq.info['Sequence'][start:end])
                            if wseq == saln.seq[0]:
                                makeme[-1] = '-' * len(makeme[-1][:lshift]) + v + '-' * len(makeme[-1][lshift+len(v):])
                                wpos = Occ['Pos']
                            else: wpos = waln.seqAlnPos(wseq,start,next=True)
                            makeme[-1] = '%s-' % rje.preZero(wpos,wseq.seqLen()) + makeme[-1]
                        waln._addSeq(wseq.info['Name'],string.join(makeme,'-XXXXXXXXXX-'))    #!#Add positions at some point
                        #x#self.deBug(waln.seq[-1].info)   #x#string.join(makeme,'-XXXXXXXXXX-'))
                    waln.saveFasta()


            ## ~ [3b] OLD SLiMSearch Code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not 'new coding working':
                motif_occ = self.obj['SlimList'].motifOcc(byseq=True)
                for seq in self.seqs():
                    if self.getBool('Masked'): seq.info['Sequence'] = seq.info['MaskSeq']
                    occlist = []
                    if motif_occ.has_key(seq):
                        for Motif in motif_occ[seq].keys():
                            occlist += motif_occ[seq][Motif]
                    saln = self.obj['SlimList'].obj['SLiMCalc'].singleProteinAlignment(seq,occlist,usegopher=False,savefasta=False)
                    saln.seq[0].info['Name'] = '%s-%s Motifs' % (seq.shortName(),self.dataset())
                    if self.getBool('Masked'):
                        saln.seq[1].info['Sequence'] = seq.info['PreMask'] 
                        saln._addSeq('%s-masked' % seq.shortName(),seq.info['MaskSeq'])
                    paln.seq += saln.seq[0:]
            ## ~ [3c] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            paln.saveFasta()
            
            ### ~ [4] Motif alignments with and without masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] No Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['SlimList'].info['Name'] = '%s SLiMProb Motifs' % self.dataset()
            self.obj['SlimList'].motifAlignments('%s.motifaln.fas' % self.runBase())
            self.dict['Output']['motifaln'] = '%s.motifaln.fas' % self.runBase()
            ## ~ [4b] With Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Masked'):
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq']
                self.obj['SlimList'].motifAlignments('%s.maskaln.fas' % self.runBase())
                self.dict['Output']['maskaln'] = '%s.maskaln.fas' % self.runBase()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask']

        except: self.errorLog('Problem with additional SLiMProb output for %s' % self.dataset())
#########################################################################################################################
    def teiresias(self):     ### Output in TEIRESIAS format along with masked sequence (input) fasta file.
        '''Output in TEIRESIAS format along with masked sequence (input) fasta file.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('Basefile'):
                self.baseFile(rje.baseFile(self.obj['SeqList'].info['Name']))
                if self.getStrLC('ResDir'):
                    rje.mkDir(self,self.getStr('ResDir'))
                    self.baseFile(self.getStr('ResDir') + self.dataset())

            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Rank File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfile = self.runBase() + '.out'
            OUT = open(outfile,'w')
            ## ~ [2b] Header ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            OUT.write(string.join(['##########################################################',
                                   '#                                                        #',
                                   '#                       FINAL RESULTS                    #',
                                   '#                                                        #',
                                   '##########################################################',''],'\n'))
            ## ~ [2c] Motifs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slimocc = self.obj['SLiMList'].motifOcc(nested=False)
            sx = len(slimocc)
            for slim in slimocc:    ### These have been pre-filtered using self.filterSLiMs()
                if not slimocc[slim]:
                    sx -= 1
                    continue
                OUT.write('%d\t%d\t%s' % (slim.occNum(),slim.seqNum(),slim.pattern()))
                for occ in slimocc[slim]: OUT.write(' %d %d' % (self.seqs().index(slim['Seq']),slim['Pos']))
                OUT.write('\n')
            OUT.close()
            self.printLog('#OUT','%s shared patterns output to %s' % (rje.integerString(sx),outfile))
                        
        except:
            self.errorLog('Major disaster with SLiMProb.teiresias()')
            raise
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output.
        Individual outputs can be identified/parsed:

        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        # OUTFMT:
        ...

        Outputs available:
            main = main results file (extras=-1)
            motifs = Input motifs for searching (extras=-1)
            seqin = Input file (extras=-1)
            occ = occurrence file (extras=0)
            upc = UPC file (extras=0)
            slimdb = Fasta file used for UPC generation etc. (extras=0)
            masked = masked.fas (extras=1)
            mapping = mapping.fas file (extras=2)
            motifaln = motif alignments (extras=2)
            maskaln = masked motif alignments (extras=2)
            dismatrix = *.dis.tdt file (extras=3)

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfmt = self.getStrLC('Rest')
            if not self.extras(0) and outfmt in ['main','occ','upc','slimdb','seqin','motifs']: self.setInt({'Extras':0})
            if not self.extras(1) and outfmt in ['masked']: self.setInt({'Extras':1})
            if not self.extras(2) and outfmt in ['mapping','motifaln','maskaln']: self.setInt({'Extras':2})
            if not self.extras(3) and outfmt in ['dismatrix']: self.setInt({'Extras':3})
            if outfmt == 'default' and not self.extras(2): self.setInt({'Extras':2})
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        output = ['main','motifs','seqin']
        if self.extras(0): output = ['main','occ','upc','motifs']
        if self.extras(1): output += ['masked']
        if self.extras(2): output += ['mapping','motifaln','maskaln']
        if self.extras(3): output += ['dismatrix']
        if self.extras(0): output += ['slimdb','seqin']
        return output
#########################################################################################################################
### End of SECTION II: SLiMProb Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: MAIN PROGRAM                                                                                           #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return 
        
    ### Rest of Functionality... ###
    try: SLiMProb(mainlog,cmd_list).run()
        
    ### End ###
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
### END OF SECTION III                                                                                                  #
#########################################################################################################################
