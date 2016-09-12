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
Program:      SLiMSearch
Description:  Short Linear Motif Search tool
Version:      1.7.1
Last Edit:    03/12/15
Citation:     Davey, Haslam, Shields & Edwards (2010), Lecture Notes in Bioinformatics 6282: 50-61. 
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    SLiMSearch is a tool for finding pre-defined SLiMs (Short Linear Motifs) in a protein sequence database. SLiMSearch
    can make use of corrections for evolutionary relationships and a variation of the SLiMChance alogrithm from
    SLiMFinder to assess motifs for statistical over- and under-representation. SLiMSearch is a replacement for PRESTO
    and uses many of the same underlying modules.

    Benefits of SLiMSearch that make it more useful than a lot of existing tools include:
    * searching with mismatches rather than restricting hits to perfect matches.
    * optional equivalency files for searching with specific allowed mismatched (e.g. charge conservation)
    * generation or reading of alignment files from which to calculate conservation statistics for motif occurrences.
    * additional statistics, including protein disorder, surface accessibility and hydrophobicity predictions
    * recognition of "n of m" motif elements in the form <X:n:m>, where X is one or more amino acids that must occur n+
    times across which m positions. E.g. <IL:3:5> must have 3+ Is and/or Ls in a 5aa stretch.

    Main output for SLiMSearch is a delimited file of motif/peptide occurrences but the motifaln=T and proteinaln=T also
    allow output of alignments of motifs and their occurrences. The primary outputs are named *.csv for the occurrence
    data and *.summary.csv for the summary data for each motif/dataset pair. 
    
    NOTE: SLiMSearch has now been largely superseded by SLiMProb for motif statistics.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ###
    motifs=FILE     : File of input motifs/peptides [None]
                      Single line per motif format = 'Name Sequence #Comments' (Comments are optional and ignored)
                      Alternative formats include fasta, SLiMDisc output and raw motif lists.
    seqin=FILE      : Sequence file to search [None]
    batch=LIST      : List of sequence files for batch input (wildcard * permitted) []
    maxseq=X        : Maximum number of sequences to process [0]
    maxsize=X       : Maximum dataset size to process in AA (or NT) [100,000]
    maxocc=X        : Filter out Motifs with more than maximum number of occurrences [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resfile=FILE    : Main SLiMSearch results table [slimsearch.csv]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMSearch/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMSearch/]
    force=T/F       : Force re-running of BLAST, UPC generation and search [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SearchDB Options I: Input Protein Sequence Masking ###
    masking=T/F     : Master control switch to turn off all masking if False [False]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out [EM,DOMAIN,TRANSMEM]
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    casemask=X      : Mask Upper or Lower case [None]
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    metmask=T/F     : Masks the N-terminal M [False]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []

    ### SearchDB Options II: Evolutionary Filtering  ###
    efilter=T/F     : Whether to use evolutionary filter [False]
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
    extras=X        : Whether to generate additional output files (alignments etc.) [1]
                        - 0 = No output beyond main results file
                        - 1 = Generate additional outputs (alignments etc.)
    pickle=T/F      : Whether to save/use pickles [True]
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc and *.pickle files
                        - 2 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)
    * See also rje_slimcalc options for occurrence-based calculations and filtering *
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import pickle, os, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),"../libraries/"))
import rje, rje_seq, rje_sequence, rje_scoring, rje_slim, rje_slimcore, rje_slimcalc, rje_slimlist, rje_zen
#import rje_dismatrix_V2 as rje_dismatrix
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Standardised masking options. Still not fully tested.
    # 1.1 - Added background=FILE option for determing mean(p1+) for SLiMs based on background file.
    # 1.2 - Added maxsize option.
    # 1.3 - Add aamask option (and alphabet)
    # 1.4 - Fixed zero-size UPC bug.
    # 1.5 - Add MaxOcc setting.
    # 1.6 - Minor tweaks to Log output. Add option for UPC number in occ output.
    # 1.7 - Modified to work with GOPHER V3.0.
    # 1.7.1 - Minor modification to docstring. Preparation for update to SLiMSearch 2.0 optimised for proteome searches.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Finish implementation with current methods/options!
    # [ ] : Add occurrence and SLiM filtering.
    # [ ] : Reinstate the expcut=X option for filtering motifs based on expected occurrences.
    # [ ] : Remove ORTHALN and make mapping like SLiMFinder
    # [ ] : Add "startfrom" option.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('SLiMSearch', '1.7', 'December 2012', '2007')
    description = 'Short Linear Motif Regular Expression Search Tool'
    author = 'Dr Richard J. Edwards.'
    comments = ['This program is under development and may contain bugs!',rje_zen.Zen().wisdom()]
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
#########################################################################################################################
def cmdHelp(info=None,out=None,cmd_list=[]):   ### Prints *.__doc__ and asks for more sys.argv commands
    '''Prints *.__doc__ and asks for more sys.argv commands.'''
    try:
        if not info: info = makeInfo()
        if not out: out = rje.Out()
        help = cmd_list.count('help') + cmd_list.count('-help') + cmd_list.count('-h')
        if help > 0:
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
### SECTION II: SLiMSearch Class                                                                                        #
#########################################################################################################################
class SLiMSearch(rje_slimcore.SLiMCore):     
    '''
    Short Linear Motif Regular Expression Search Tool Class. Author: Rich Edwards (2007). This module inherits the
    SLiMFinder class from SLiMFinder, which handles the masking of datasets and correcting for evolutionary relationships
    (if this desired).

    Info:str
    - AAFreq = Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    - AltDis = Alternative all by all distance matrix for relationships [None]
    - Background = Use observed support in background file for over-representation calculations [None]
    - BuildPath = Alternative path to look for existing intermediate files [SLiMSearch/]
    - CaseMask = Mask Upper or Lower case [None]
    - CompMask = Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    - GablamDis = Alternative GABLAM results file [None]
    - Input = Original name (and path) of input file
    - ResDir = Redirect individual output files to specified directory [SLiMSearch/]
    - ResFile = If FILE is given, will also produce a table of results in resfile [slimsearch.csv]
    - RunID = Run ID for resfile (allows multiple runs on same data) [DATE:TIME]
    
    Opt:boolean
    - DisMask = Whether to mask ordered regions (see rje_disorder for options) [False]
    - Force = whether to force recreation of key files [False]
    - EFilter = Whether to use evolutionary filter [True]
    - Extras = Whether to generate additional output files (alignments etc.) [False]
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

    Stat:numeric
    - MaxOcc = Filter out Motifs with more than maximum number of occurrences [0]
    - MaxSeq = Maximum number of sequences to process [500]
    - MaxSize = Maximum dataset size to process in AA (or NT) [1e5]
    - MST = MST corrected size for whole dataset
    - SaveSpace = Delete "unneccessary" files following run (see Manual for details) [0]
    - SeqOcc = Restrict to sequences with X+ occurrences (adjust for high frequency SLiMs) [1]
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
    - Background = Background SLiMSearch object containing data for background occurrences
    - SeqList = main SeqList object containing dataset to be searched
    - SlimList = MotifList object handling motif stats and filtering options
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.infolist = ['AAFreq','AltDis','BuildPath','CaseMask','CompMask','GablamDis','Input','ResDir','ResFile',
                         'RunID','Background']
        self.optlist = ['DisMask','Force','EFilter','Extras','LogMask','Masked','Masking','MaskM','MaskFreq','SmearFreq',
                        'TarGZ','Teiresias','ConsMask','Pickle','OccUPC']
        self.statlist = ['MaxSeq','MST','SaveSpace','StartTime','WallTime','HomCut','SeqOcc','MaxSize','MaxOcc']
        self.listlist = ['AAMask','Alphabet','Batch','FTMask','Headers','IMask','NewScore','OccFilter','OccHeaders','OccStats',
                         'UP','NoHits']
        self.dictlist = ['AAFreq','MST','NewScore','OccFilter','StatFilter']
        self.objlist = ['Background','SeqList','SlimList']
        ### Defaults ###
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
        self.setInfo({'BuildPath':rje.makePath('SLiMSearch/'),'CompMask':'None',
                      'ResDir':rje.makePath('SLiMSearch/'),'ResFile':'slimsearch.csv'})
        self.setStat({'SeqOcc':1,'MaxSeq':0,'MaxSize':1e5})
        self.setOpt({'ConsMask':False,'EFilter':False,'MaskM':False,'DisMask':False,'Masking':True})
        self.dict['MaskPos'] = {}
        t = time.localtime(time.time())
        self.info['RunID'] = self.info['Date'] = '%s%s%s-%s:%s' % (str(t[0])[-2:],rje.preZero(t[1],12),rje.preZero(t[2],31),rje.preZero(t[3],24),rje.preZero(t[4],60)) #x#time.ctime(self.stat['StartTime'])
        ### Object setup ###
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        self.obj['SlimList'].loadMotifs()
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
                self._cmdReadList(cmd,'info',['CaseMask','CompMask','RunID'])
                self._cmdReadList(cmd,'opt',['DisMask','Force','EFilter','Extras','LogMask','Masked','Masking','MaskM',
                                             'MaskFreq','SmearFreq','TarGZ','Teiresias','OccUPC'])
                self._cmdReadList(cmd,'int',['MaxSeq','SaveSpace','HomCut','SeqOcc','MaxOcc'])
                self._cmdReadList(cmd,'stat',['MST','WallTime','MaxSize'])
                self._cmdReadList(cmd,'list',['Batch','FTMask','IMask'])
            except: self.log.errorLog('Problem with cmd:%s' % cmd)
        self.deBug(self.obj['SlimList'].obj['SLiMCalc'].info)
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
    def run(self,batch=False):  ### Main SLiMSearch Run Method
        '''
        Main SLiMSearch Run Method:
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
            if self.opt['Webserver'] and self.obj['SlimList'].motifNum() < 1:
                self.serverEnd('NoMotif',exit=False)
                return False
            seqcmd = ['gnspacc=T','usecase=T','accnr=F','seqnr=F'] + self.cmd_list + ['autoload=T','query=None']
            self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            if self.background(): self.printLog('#BACK','Background data from "%s" successful' % self.info['Background'])
            self.setupBasefile()
            ## ~ [1a] Batch Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch and self.seqNum() < 1:   # No sequences loaded - use batch mode
                batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                if not batchfiles: self.log.errorLog('No input files found!',printerror=False)
                else:
                    mycmd = self.cmd_list[0:]
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        self.log.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.run(batch=True)
                        self.log.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.opt['Append'] = True
                        bx += 1
                        self.log.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.opt['Win32'] and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return
            ## ~ [1b] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~~~ ##
            if self.stat['MaxSeq'] > 0 and self.stat['MaxSeq'] < self.seqNum():
                self.log.printLog('#MAX','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.integerString(self.obj['SeqList'].seqNum()),rje.integerString(self.stat['MaxSeq'])))
                self.serverEnd('MaxSeq',exit=False)
                return False
            if self.stat['MaxSize'] > 0 and self.stat['MaxSize'] < self.aaNum():
                self.log.printLog('#MAX','%s = %s %s > Max %s %s. Analysis terminated.' % (self.dataset(),rje.integerString(self.aaNum()),self.units(),rje.integerString(self.stat['MaxSize']),self.units()))
                self.serverEnd('MaxSize',exit=False)
                return False
            self.stat['StartTime'] = time.time()

            ### ~ [2] Read/Generate UPC, mask input & calculate/read AAFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.makeUPC():
                self.errorLog('Error during makeUPC(). Abandoning %s run' % self.dataset(),printerror=False)
                self.serverEnd('Crash','makeUPC()')
                return False
            pickled = self.pickleMe(load=self.opt['Pickle'])  # Returns appropriate pickled SLiMFinder Object, else None
            #x#pickled = None  #!# Not reloading SlimList properly. #!#
            if pickled: self = pickled  ## Replace me with my pickle!
            else: self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
            if self.opt['MaskFreq']: self.makeAAFreq()
            else: 
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                self.makeAAFreq()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
            self.adjustAATotals()

            ### ~ [3] Perform basic search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                self.searchDB()     #!# Or read in results?! #!#
                self.pickleMe()     # Generates pickling for speedy re-running
                self.wallTime()
            #self.deBug(self.obj['SeqList'].seq[0].info)

            ### ~ [4] Generate extra stats, filter and output etc. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Setup output & filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupResults()                 
            self.backupOrCreateResFile()
            self.obj['SlimList'].setupFilters(occheaders=self.list['OccHeaders'])
            #x#self.dict['OccFilter'] = rje_scoring.setupStatFilter(self,self.list['OccHeaders'],self.list['OccFilter'])
            ## ~ [4b] Calculate extra occurrence attributes & output occurrences ~~~~~~~~~~~~~~~~~~ ##
            self.calculateOccAttributes()
            self.processSeqOccs()       # Converts for output, filtering as desired
            ## ~ [4c] Combine occurrences and generate summary, with statistics ~~~~~~~~~~~~~~~~~~~ ##
            self.combMotifOccStats()
            if self.opt['Extras']: self.extraOutput()   # MotifList Outputs 
            self.tarZipSaveSpace()      # Tarring, Zipping and Saving Space 

            ### ~ [5] End ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.log.printLog('#RES','SLiMSearch results output to %s and %s*' % (self.info['ResFile'],self.info['Basefile']))
            if self.opt['Win32'] and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.stat['WallTime'] <= 0 or (time.time() - self.stat['StartTime']) < (self.stat['WallTime']*3600): raise
            return False # Walltime reached
        except:
            self.log.errorLog('Error in SLiMSearch.run()',printerror=True,quitchoice=False)
            raise   # Delete this if method error not terrible
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMSearch object for new batch run
        '''Returns SLiMFinder object for new batch run.'''
        return SLiMSearch(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.opt['Append']])
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def background(self):   ### Sets up and/or returns Background SLiMSearch object
        '''Sets up and/or returns Background SLiMSearch object.'''
        try:### ~ [0] ~ Simple case first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['Background']: return self.obj['Background']
            if not rje.checkForFile(self.info['Background']): return None

            ### ~ [1] ~ Setup Background SLiMSearch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            bg = SLiMSearch(self.log,self.cmd_list)
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','seqin=%s' % self.info['Background']]
            bg.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            bg.setupBasefile()
            ## ~ [1a] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~~~ ##
            if self.stat['MaxSeq'] > 0 and self.stat['MaxSeq'] < bg.seqNum():
                self.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (bg.dataset(),rje.integerString(bg.obj['SeqList'].seqNum()),rje.integerString(self.stat['MaxSeq'])))
                raise ValueError
            self.stat['StartTime'] = time.time()

            ### ~ [2] Read/Generate UPC, mask input & calculate/read AAFreq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not bg.makeUPC(): self.errorLog('Error during makeUPC(). Abandoning %s run' % bg.dataset(),printerror=False); raise ValueError
            pickled = bg.pickleMe(load=True)  # Returns appropriate pickled SLiMFinder Object, else None
            if pickled: bg = pickled  ## Replace me with my pickle!
            else: bg.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']

            ### ~ [3] Perform basic search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                bg.searchDB()     #!# Or read in results?! #!#
                bg.pickleMe()     # Generates pickling for speedy re-running
            self.obj['Background'] = bg
            return bg

        except ValueError: self.info['Background'] = 'None'; return None
        except: self.errorLog('SLiMSearch.background error!'); return None
#########################################################################################################################
    def setupBasefile(self):    ### Sets up self.info['Basefile'] and self.info['Input']
        '''Sets up self.info['Basefile'].'''
        ### Basefile and ResDir ###
        if self.info['Basefile'].lower() in ['','none']:
            self.info['Basefile'] = rje.baseFile(self.obj['SeqList'].info['Name'])
            if self.info['ResDir'].lower() not in ['','none']:
                rje.mkDir(self,self.info['ResDir'])
                slimbase = os.path.split(rje.baseFile(self.obj['SlimList'].info['Name']))[1]
                if slimbase.lower() in ['','none']: self.info['Basefile'] = self.info['ResDir'] + '%s' % (self.dataset())
                else: self.info['Basefile'] = self.info['ResDir'] + '%s.%s' % (self.dataset(),slimbase)
        ### Input ###
        if not self.info['Input']: self.info['Input'] = self.obj['SeqList'].info['Name']
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
            if not self.slimNum(): return self.log.errorLog('No SLiMs for SLiMSearch.searchDB()',printerror=False)
            (ox,sx,stot) = (0,0.0,self.seqNum()*self.slimNum())    # Counter for Log output
            stxt = 'Searching %s sequences with %d SLiMs: ' % (rje.integerString(self.seqNum()),self.slimNum())
            for seq in self.seqs():
                ## ~ [1a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.log.printLog('\r#SEARCH','%s %.1f%% (%s occ)' % (stxt,(sx/stot),rje.integerString(ox)),log=False,newline=False)
                seq.deGap()
                ## ~ [1b] Search ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for slim in self.slims():
                    sx += 100.0
                    ox += len(slim.searchSequence(seq))    # {Pos,Variant,Match,ID,MisMatch}
                    self.log.printLog('\r#SEARCH','%s %.1f%% (%s occ)' % (stxt,(sx/stot),rje.integerString(ox)),log=False,newline=False)
                    self.wallTime()

            ### ~ [2] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
            self.opt['Searched'] = True
            self.printLog('\n#SEARCH','%s sequences searched for %s motifs: %s occ.' % (rje.integerString(self.seqNum()),rje.integerString(self.slimNum()),rje.integerString(ox)))
        except:
            self.log.errorLog('Error in SLiMSearch.searchDB()')
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
                    self.log.printLog('\r#OCC','Processing occurrences: %.1f%%' % ((sb+sx)/stot),newline=False,log=False)
                    sx += sj
                    ## ~ [1a] Update dictionaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    #i# [Dataset,RunID,Masking,Motif,RunID,Seq,Start_Pos,End_Pos,Prot_Len,Match,Variant,MisMatch,Desc]
                    results = {}    # Dictionary of results dictionaries to output
                    occlist = slim.dict['Occ'][seq][0:]
                    for occ in occlist[0:]:
                        ox += 1
                        occ['Dataset'] = self.dataset()
                        occ['Motif'] = slim.info['Name']
                        occ['RunID'] = self.info['RunID']
                        occ['Masking'] = self.maskText()
                        occ['Seq'] = seq.shortName()
                        occ['Start_Pos'] = occ['Pos']
                        occ['End_Pos'] = occ['Pos'] + len(occ['Match']) - 1
                        occ['Prot_Len'] = seq.aaLen()
                        occ['Desc'] = seq.info['Description']
                        occ['Pattern'] = slim.pattern()
                        if self.opt['OccUPC']: occ['UPC'] = self.list['UP'].index(self.getUP(seq)) + 1
                        results['%s-%s-%s' % (rje.preZero(occ['Start_Pos'],occ['Prot_Len']),rje.preZero(occ['End_Pos'],occ['Prot_Len']),occ['Variant'])] = occ
                    ## ~ [1b] Filtering & Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    results = rje_scoring.statFilter(self,results,self.dict['OccFilter'])  #!# Add restrict/exclude #!#
                    for occ in occlist[0:]:
                        okey = '%s-%s-%s' % (rje.preZero(occ['Start_Pos'],occ['Prot_Len']),rje.preZero(occ['End_Pos'],occ['Prot_Len']),occ['Variant'])
                        if okey in results: rje.delimitedFileOutput(self,self.info['ResFile'],self.resHead('OccHeaders'),datadict=results[okey])
                        else: slim.dict['Occ'][seq].remove(occ)
                    ## ~ [1c] SeqOcc Reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    if self.stat['SeqOcc'] > 1 and self.stat['SeqOcc'] > len(slim.dict['Occ'][seq]): slim.dict['Occ'].pop(seq)
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
            sec = int(time.time() - self.stat['StartTime'] + 0.5)
            (hour,min) = (0,0)
            while sec >= 60: (min,sec) = (min+1,sec-60)
            while min >= 60: (hour,min) = (hour+1,min-60)
            self.info['RunTime'] = '%s:%s:%s' % (rje.preZero(hour,24),rje.preZero(min,60),rje.preZero(sec,60))
            totalaa = 0     # Total number of AA in dataset
            for seq in self.seqs():
                if self.opt['Masked']: seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                totalaa += seq.nonX()
            general = {'Dataset':self.dataset(),'RunID':self.info['RunID'],'Masking':self.maskText(),
                       'RunTime':self.info['RunTime'],'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa}
            ## ~ [2b] SLiM Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for slim in self.slims():
                datadict = rje.combineDict({'Motif':slim.info['Name'],'Pattern':slim.pattern()},general)
                datadict = rje.combineDict(datadict,slim.stat)
                for dkey in datadict.keys():
                    if string.split(dkey,'_')[0] in ['E','p','pUnd']: datadict[dkey] = rje_slim.expectString(datadict[dkey])

                rje.delimitedFileOutput(self,self.info['SummaryFile'],self.resHead(),datadict=datadict)
            self.printLog('#OUT','Summary data for %s saved to %s' % (self.dataset(),self.info['SummaryFile']))
        except: self.log.errorLog('Error in combMotifOccStats()')
#########################################################################################################################
    ### <6> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
#CORE# def makeAAFreq(self):   ### Makes an initial AAFreq dictionary containing AA counts only (including Xs)
#########################################################################################################################
    def adjustAATotals(self):   ### Adjusts AA Totals following masking     #!# Changed from Core #!#
        '''Adjusts AA Totals following masking and makes appropriate frequencies. If maskfreq=T then the amino acid counts
        will be exactly as they were. If maskfreq=F, however, frequencies will need to be adjust for the new number of
        masked and non-masked amino acids.'''
        try:
            ### Simple PreMasking Procedure ###
            if self.opt['MaskFreq']:
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
                if self.opt['SmearFreq']: self.smearAAFreq()
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
            if self.info['AAFreq'].lower() not in ['','none']: prex = self.dict['AAFreq']['Dataset']['Total']
            if prex != postx: self.printLog('#ADJ','Effective dataset size reduced from %s %s to %s %s' % (rje.integerString(prex),self.units(),rje.integerString(postx),self.units()))
            if self.opt['SmearFreq']: self.smearAAFreq()
        except:
            self.log.errorLog('Problem during SLiMSearch.adjustAATotals()')
            raise
#########################################################################################################################
#CORE# def smearAAFreq(self):  ### Equalises AAFreq across UPC. Leaves Totals unchanged.
#########################################################################################################################
    def slimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM. Modified from slimcore for SLiMSearch.'''
        try:### ~ [1] Setup SLiMSearch Motif stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
                if self.stat['SeqOcc'] > 1: k = self.stat['SeqOcc']
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
                    if self.stat['SeqOcc'] > 1: k = self.stat['SeqOcc']
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
            self.log.errorLog('Error with slimProb(%s)' % slim.info['Sequence'])
            slim.stat['p_Occ'] = slim.stat['p_Seq'] = slim.stat['p_UPC'] = 1.0
            #self.deBug(slim)
#########################################################################################################################
    def slimProbBG(self,slim):  ### Calculate Probabilities for given SLiM using Background Support
        '''Calculate Probabilities for given SLiM using Background Support.'''
        try:### ~ [1] Setup SLiMSearch Motif stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
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
            self.log.errorLog('Error with slimProb(%s)' % slim.info['Sequence'])
            slim.stat['p_Occ'] = slim.stat['p_Seq'] = slim.stat['p_UPC'] = 1.0
            #self.deBug(slim)
#########################################################################################################################
    ### <8> ### SLiMBuild Pickle Methods                                                                                #
#########################################################################################################################
    def myPickle(self): return self.maskText('')   ### Returns pickle identifier
#########################################################################################################################
    def pickleMe(self,load=False):  ### Loads existing pickle, or saves pickle for later!
        '''Saves pickle for later!.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load and self.opt['Force']: return None      # Re-run search!
            elif not load and not self.opt['Pickle']: return None  # Do not save pickle
            mypickle = self.myPickle()  ## Pickle name ##
            
            ### ~ [2] Load existing pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load:
                newme = None    ## New SLiMSearch object, loaded from Pickle
                ## Check for file and load ##
                for pfile in [self.info['Basefile'],'%s%s' % (self.info['BuildPath'],self.dataset())]:
                    if not self.opt['Win32'] and os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)):
                        if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            if rje.isYounger('%s.%s.pickle.gz' % (pfile,mypickle),'%s.%s.pickle' % (pfile,mypickle)) == '%s.%s.pickle.gz' % (pfile,mypickle):
                                os.unlink('%s.%s.pickle' % (pfile,mypickle))
                        if not os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            try: os.system('gunzip %s.%s.pickle.gz' % (pfile,mypickle))
                            except: self.log.errorLog('Cannot unzip %s.%s.pickle.gz' % (pfile,mypickle))
                    if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                        self.log.printLog('#LOAD','Attempting to load SLiMSearch pickle.',log=False)
                        newme = pickle.load(open('%s.%s.pickle' % (pfile,mypickle),'r'))
                        self.log.printLog('#LOAD','SLiMSearch intermediate loaded: %s.%s.pickle.' % (pfile,mypickle))
                        if not self.opt['Win32']:
                            try:
                                if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                                os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                                self.log.printLog('#GZIP','SLiMSearch %s.%s.pickle zipped.' % (pfile,mypickle))
                            except: self.log.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
                        break
                if not newme: return None
                ## Check other pertinent attributes - masking and additional filtering ##
                if self.obj['SlimList'].nameList() != newme.obj['SlimList'].nameList():
                    self.printLog('#PICKLE','Motif attributes changed since pickle: will not use'); return None
                ## Note that MustHave and OccFilter filtering currently occur *after* SLiMBuild only ##
                changes = []
                for var in ['CompMask','CaseMask','MotifMask']:         # Info
                    if self.info[var] != newme.info[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['Masking','DisMask','MaskM','ConsMask']:   # Opt
                    if self.opt[var] != newme.opt[var]: changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['FTMask','IMask']:  #x#,'Equiv']:      # List   
                    slist = self.list[var][0:]
                    nlist = newme.list[var][0:]
                    slist.sort()
                    nlist.sort()
                    if slist != nlist:
                        changes.append(self.log.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                ## Recreate or use pickle but add new commands ##
                if changes and (self.stat['Interactive'] < 0 or rje.yesNo('%d SLiMBuild parameter mismatches with pickle. Create new pickle?' % len(changes))):
                    self.log.printLog('#PICKLE','Parameters changed. Making new pickle.')
                    return None
                newme.cmd_list = self.cmd_list
                newme.setInfo(self.info)
                newme.setStat(self.stat)
                self.opt['Masked'] = newme.opt['Masked']
                newme.setOpt(self.opt)
                newme.info['ResFile'] = self.info['ResFile']
                newme.info['ResDir'] = self.info['ResDir']
                newme.info['BuildPath'] = self.info['BuildPath']
                newme.stat['StartTime'] = self.stat['StartTime']
                newme.setLog(self.log)
                self.list['UP'] = newme.list['UP']
                self.dict['AAFreq'] = newme.dict['AAFreq']
                self.dict['MST'] = newme.dict['MST']
                newme.list = self.list
                newme.dict = self.dict
                return newme
                
            ### ~ [3] Save new pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if 'Pickle' in self.opt and not self.opt['Pickle']:
                self.log.printLog('#PICKLE','SLiMSearch pickling disabled with pickle=F.',log=True)
                return None
            self.log.printLog('#SAVE','Attempting to save SLiMSearch with pickle.',log=False)
            pickle.dump(self,open('%s.%s.pickle' % (self.info['Basefile'],mypickle),'w'))
            self.log.printLog('#SAVE','SLiMSearch intermediate saved as %s.%s.pickle (Python pickle).' % (self.info['Basefile'],mypickle))
            if not self.opt['Win32']:
                try:
                    pfile = self.info['Basefile']
                    if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                    os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                    self.log.printLog('#GZIP','SLiMSearch %s.%s.pickle zipped.' % (pfile,mypickle))
                except: self.log.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
            return None
        
        except:
            self.log.errorLog('Major problem with SLiMSearch pickling!')
            return None
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
            if self.opt['OccUPC']: self.list['OccHeaders'].append('UPC')
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
            if self.info['ResFile'].lower() not in ['','none'] and self.info['ResFile'].find('.') < 0: self.info['ResFile'] += '.csv'
            spex = os.path.splitext(self.info['ResFile'])
            self.info['SummaryFile'] = spex[0] + '.summary' + spex[1]
        except:
            self.log.errorLog('Problem with SLiMSearch.setupResults()')
            raise
#########################################################################################################################
    def backupOrCreateResFile(self):    ### Backups up and/or creates main results file
        '''Backups up and/or creates main results file.'''
        if self.info['ResFile'].lower() not in ['','none']:
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.info['ResFile'],write=True))
            rje.delimitedFileOutput(self,self.info['ResFile'],self.resHead('OccHeaders'),delimit,rje_backup=True)
        if self.info['SummaryFile'].lower() not in ['','none']:
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.info['SummaryFile'],write=True))
            rje.delimitedFileOutput(self,self.info['SummaryFile'],self.resHead(),delimit,rje_backup=True)
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
            if self.opt['Masked']:
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq']
                self.obj['SeqList'].saveFasta(seqfile='%s.masked.fas' % self.info['Basefile'])
            
            ### ~ [2] Full Protein Alignments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #alndir = rje.makePath('%s_ORTHALN' % self.info['Basefile'])     #!# Modify to match SLiMFinder #!#
            #self.obj['SlimList'].proteinAlignments(alndir=alndir)

            ### ~ [3] Protein Alignments with Motifs and Masking marked ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alncmd = ['seqin=None','accnr=F','seqnr=F','autofilter=F','align=F','gnspacc=F','unkspec=T'] 
            paln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
            paln.info['Name'] = '%s.mapping.fas' % self.info['Basefile']
            if self.opt['Webserver']: waln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd+['replacechar=F'])
            ## ~ [3a] Generate Alignment Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('#MAP','Generating motif mapping alignments')
            motif_occ = self.obj['SlimList'].motifOcc(byseq=True)
            for seq in self.seqs():
                if self.opt['Masked']: seq.info['Sequence'] = seq.info['PreMask']   #x# seq.info['MaskSeq']
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
                if self.opt['Masked']:
                    alnseq = saln.seq[1].info['Sequence']
                    #x# saln.seq[1].info['Sequence'] = rje_sequence.mapGaps(seq.info['PreMask'],alnseq,self)
                    saln._addSeq('%s-masked' % seq.shortName(),rje_sequence.mapGaps(seq.info['MaskSeq'],alnseq,self))
                    saln.seq = saln.seq[:1] + saln.seq[-1:] + saln.seq[1:-1]
                paln.seq += saln.seq[0:]
                ## Special webserver output ##
                if self.opt['Webserver'] and occlist:
                    #self.deBug(occlist)
                    waln.info['Name'] = '%s.webserver.%s.fas' % (self.info['Basefile'],seq.shortName())
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
                    if self.opt['Masked']: seq.info['Sequence'] = seq.info['MaskSeq']
                    occlist = []
                    if motif_occ.has_key(seq):
                        for Motif in motif_occ[seq].keys():
                            occlist += motif_occ[seq][Motif]
                    saln = self.obj['SlimList'].obj['SLiMCalc'].singleProteinAlignment(seq,occlist,usegopher=False,savefasta=False)
                    saln.seq[0].info['Name'] = '%s-%s Motifs' % (seq.shortName(),self.dataset())
                    if self.opt['Masked']:
                        saln.seq[1].info['Sequence'] = seq.info['PreMask'] 
                        saln._addSeq('%s-masked' % seq.shortName(),seq.info['MaskSeq'])
                    paln.seq += saln.seq[0:]
            ## ~ [3c] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            paln.saveFasta()
            
            ### ~ [4] Motif alignments with and without masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] No Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.obj['SlimList'].info['Name'] = '%s SLiMSearch Motifs' % self.dataset()
            self.obj['SlimList'].motifAlignments('%s.motifaln.fas' % self.info['Basefile'])
            ## ~ [4b] With Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['Masked']:
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq']
                self.obj['SlimList'].motifAlignments('%s.maskaln.fas' % self.info['Basefile'])
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask']

        except: self.log.errorLog('Problem with additional SLiMSearch output for %s' % self.dataset())
#########################################################################################################################
    def teiresias(self):     ### Output in TEIRESIAS format along with masked sequence (input) fasta file.
        '''Output in TEIRESIAS format along with masked sequence (input) fasta file.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Basefile'].lower() in ['','none']:
                self.info['Basefile'] = rje.baseFile(self.obj['SeqList'].info['Name'])
                if self.info['ResDir'].lower() not in ['','none']:
                    rje.mkDir(self,self.info['ResDir'])
                    self.info['Basefile'] = self.info['ResDir'] + self.dataset()

            ### ~ [2] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] Rank File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfile = self.info['Basefile'] + '.out'
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
            self.log.printLog('#OUT','%s shared patterns output to %s' % (rje.integerString(sx),outfile))
                        
        except:
            self.log.errorLog('Major disaster with SLiMSearch.teiresias()')
            raise
#########################################################################################################################
### End of SECTION II: SLiMSearch Class                                                                                 #
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
    try: SLiMSearch(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################
