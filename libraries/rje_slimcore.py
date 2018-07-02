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
Module:       SLiMCore
Description:  SLiMSuite core module
Version:      2.9.0
Last Edit:    19/10/17
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    This module is primarily to contain core dataset processing methods for both SLiMFinder and SLiMProv to inherit and
    use. This primarily consists of the options and methods for masking datasets and generating UPC. This module can
    therefore be run in standalone mode to generate UPC files for SLiMFinder or SLiMProb. It can also be used to generate
    "MegaSLiM" files of precomputed scores that can be used for subsequent subdata.

    In addition, the secondary MotifSeq and Randomise functions are handled here.
    
Secondary Functions:
    The "MotifSeq" option will output fasta files for a list of X:Y, where X is a motif pattern and Y is the output file.

    The "Randomise" function will take a set of input datasets (as in Batch Mode) and regenerate a set of new datasets
    by shuffling the UPC among datasets. Note that, at this stage, this is quite crude and may result in the final
    datasets having fewer UPC due to common sequences and/or relationships between UPC clusters in different datasets.

    The "EquivMaker" function will read in a BLOSUM matrix and, for a particular score cut-off, generate all equivalence
    groups for which all members have pairwise BLOSUM scores that equal or exceed the cut-off. This equivalence file can
    then be used as input for TEIRESIAS or SLiMFinder.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options ### 
    seqin=FILE      : Sequence file to search. Over-rules batch mode (and uniprotid=LIST). [None]
    batch=LIST      : List of files to search, wildcards allowed. [*.dat,*.fas]
    uniprotid=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.dat and use as seqin=FILE. []
    maxseq=X        : Maximum number of sequences to process [0]
    maxupc=X        : Maximum UPC size of dataset to process [0]
    sizesort=X      : Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMFinder/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMFinder/]
    force=T/F       : Force re-running of BLAST, UPC generation and SLiMBuild [False]
    dna=T/F         : Whether the sequences files are DNA rather than protein [False]
    alphabet=LIST   : List of characters to include in search (e.g. AAs or NTs) [default AA or NT codes]
    megaslim=FILE   : Make/use precomputed results for a proteome (FILE) in fasta format [None]
    megablam=T/F    : Whether to create and use all-by-all GABLAM results for (gablamdis) UPC generation [False]
    megaslimfix=T/F : Whether to run megaslim in "fix" mode to tidy/repair existing files [False]
    ptmlist=LIST    : List of PTM letters to add to alphabet for analysis and restrict PTM data []
    ptmdata=DSVFILE : File containing PTM data, including AccNum, ModType, ModPos, ModAA, ModCode
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Evolutionary Filtering Options ###
    efilter=T/F     : Whether to use evolutionary filter [True]
    blastf=T/F      : Use BLAST Complexity filter when determining relationships [True]
    blaste=X        : BLAST e-value threshold for determining relationships [1e=4]
    altdis=FILE     : Alternative all by all distance matrix for relationships [None]
    gablamdis=FILE  : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
    fupc=T/F        : Whether to use experimental "Fragment UPC" approach for UPC of large proteomes [False]
    domtable=FILE   : Domain table containing domain ("Type") and sequence ("Name") pairings for additional UPC [None]
    homcut=X        : Max number of homologues to allow (to reduce large multi-domain families) [0]
    extras=T/F      : Whether to generate additional output files (distance matrices etc.) [True]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Input Masking and AA Frequency Options ###
    masking=T/F     : Master control switch to turn off all masking if False [True]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out (True=EM,DOMAIN,TRANSMEM) []
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    fakemask=T/F    : Whether to invoke UniFake to generate additional features for masking [False]
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    casemask=X      : Mask Upper or Lower case [None]
    metmask=T/F     : Masks the N-terminal M (can be useful if SLiMFinder termini=T) *Also maskm=T/F* [True]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas *Also maskpos=LIST* [2:A]
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    logmask=T/F     : Whether to output the log messages for masking of individual sequences to screen [False]
    masktext=X      : Text ID to over-ride automated masking text and identify specific masking settings [None]
    maskpickle=T/F  : Whether to save/load pickle of masked input data, independent of main pickling [False]
    maskfreq=T/F    : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [True]
    aafreq=FILE     : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]    
    smearfreq=T/F   : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    qregion=X,Y     : Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
    megaslimdp=X    : Accuracy (d.p.) for MegaSLiM masking tool raw scores [4]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Output Options ###
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    pickle=T/F      : Whether to save/use pickles [True]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc, *.pickle (Pickle excluded from tar.gz with this setting)
                        - 2 = Delete all bar *.upc files (Pickle included in tar.gz with this setting)
                        - 3 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)
    iuscoredir=PATH : Path in which to save protein acc.iupred.txt score files for megaslim analysis
    protscores=T/F  : Whether to save individual protein rlc.txt files in alignment directory [False]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Additional Functions I: MotifSeq ###    
    motifseq=LIST   : Outputs fasta files for a list of X:Y, where X is the pattern and Y is the output file []
    slimbuild=T/F   : Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]

    ### Additional Functions II: Randomised datasets ###
    randomise=T/F   : Randomise UPC within batch files and output new datasets [False]
    randir=PATH     : Output path for creation of randomised datasets [Random/]
    randbase=X      : Base for random dataset name [rand]
    randsource=FILE : Source for new sequences for random datasets (replaces UPCs) [None]

    ### Additional Functions III: Equivalence File Maker ###
    blosumfile=FILE : BLOSUM file from which to make equivalence file
    equivout=FILE   : File for equivalence list output [equiv.txt]
    equivcut=X      : BLOSUM score cut-off for equivalence groups [1]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import copy, glob, os, pickle, random, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_disorder, rje_obj, rje_seq, rje_sequence, rje_seqlist, rje_slim, rje_slimcalc, rje_slimlist
import rje_uniprot
#import rje_motif_V3 as rje_motif
import rje_dismatrix_V2 as rje_dismatrix
import rje_dismatrix_V3
#import comparimotif_V3 as comparimotif
import rje_blast_V2 as rje_blast
import gablam, unifake
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation based on SLiMFinder 3.1.
    # 0.1 - Tidied with respect to SLiMFinder and SLiMSearch.
    # 0.2 - Added DNA mode.
    # 0.3 - Added relative conservation masking.
    # 0.4 - Altered TarGZ to *not* include *.pickle.gz
    # 1.0 - Standardised masking options. Added motifmask and motifcull.
    # 1.1 - Checked/updated Randomise option.
    # 1.2 - Added generation of UniFake file to maskInput() method. (fakemask=T)
    # 1.3 - Added aamask effort
    # 1.4 - Added DomTable option
    # 1.5 - Added masktext and maskpickle options to accelerate runs with large masked datasets.
    # 1.6 - Fixed occurrence table bugs.
    # 1.7 - Added SizeSort and NewUPC. Add server #END statements.
    # 1.8 - Added BLOSUM2equiv method for making equivalent lists from a BLOSUM matrix.
    # 1.9 - Minor modifications to Log output. Updated motifSeq() function to output unmasked sequences.
    # 1.10- Bypass UPC generation for single sequences.
    # 1.11- Tidied some of the module imports.
    # 1.12- Upgraded BLAST to BLAST+. Can use old BLAST with oldblast=T.
    # 1.13- Modified the savespace settings to reduce numbers of files. targz file now uses RunID not Build Info.
    # 1.14- Started adding code for Fragmented UPC (FUPC) clustering.
    # 1.15- Added pre-running GOPHER if no alndir and usegopher=T. Updated dataset() to use Input not Basefile.
    # 1.16- Preparation for SLiMCore V2.0 using newer RJE_Object.
    # 2.0 - Converted to use rje_obj.RJE_Object as base. Version 1.16 moved to legacy/.
    # 2.1 - Added megaslim=FILE option to make/use precomputed results for a proteome. Upgraded MotifSeq method.
    # 2.2 - Modified aa frequency calculations to use alphabet to generate 0.0 frequencies (rather than missing aa).
    # 2.3 - Docstring edits. Minor tweak to walltime() to close open files.
    # 2.4 - Added megaslimfix=T/F : Whether to run megaslim in "fix" mode to tidy/repair existing files [False]
    # 2.5 - Added (hidden) slimmutant=T/F : Whether to ignore '.p.\D\d+\D' at end of accnum. Made default append=True.
    # 2.6.0 - Added uniprotid=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.dat and use as seqin=FILE. []
    # 2.6.1 - Removed the maxseq default setting.
    # 2.7.0 - Updating MegaSLiM function to work with REST server. Allow megaslim=seqin. Added iuscoredir=PATH and protscores=T/F.
    # 2.7.1 - Modified iuscoredir=PATH and protscores=T/F to work without megaslim. Fixed UPC/SLiMdb issue for GOPHER.
    # 2.7.2 - Fixed iuscoredir=PATH to stop raising errors when file not previously made.
    # 2.7.3 - Fixed serverend message error.
    # 2.7.4 - Fixed walltime server bug.
    # 2.7.5 - Fixed feature masking.
    # 2.7.6 - Added feature masking log info or warning.
    # 2.7.7 - Switched feature masking OFF by default to give consistent Uniprot versus FASTA behaviour.
    # 2.7.8 - Fixed batch=FILE error for single input files.
    # 2.8.0 - Added map and failed output to REST servers and standalone uniprotid=LIST input runs.
    # 2.8.1 - Updated resfile to be set by basefile if no resfile=X setting given
    # 2.9.0 - Added separate IUPred long suffix for reusing predictions
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Reduce methods and attributes of SLiMCore object to only those necessary for function.
    # [ ] : Reduce commandline options accordingly - proper use of coreCmd() and add coreDefaults().
    # [Y] : Add proper use of seqbase and basefile for SLiMSearch results output.
    # [ ] : Update to use new rje_seqlist module?
    # [ ] : Sort out masking log output, especially for conservation masking.
    # [Y] : Update rje_blast used for UPC.
    # [ ] : Replace final rje_motiflist functions?
    # [Y] : Double-check the BLAST V2 ID clustering and then disable the test() option in makeUPC()
    # [Y] : Consider running GOPHER fully before masking if usegopher=T otherwise not sure that forks will be used.
    # [ ] : Check the basefile test()/dev() option and then delete excess code.
    # [Y] : Tidy in preparation for update to new RJE_Object.
    # [Y] : Update to new RJE_Object.
    # [Y] : Add megaslim=FILE option to precompute results for a proteome.
    # [ ] : Remove slimbuild=T/F? Appears to do nothing.
    # [ ] : Add datasets=FILE : Uses seqlist=FILE hubfield=X spokefield=X and resdir=PATH to make resdir/hub.fas datasets
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('SLiMCore', '2.9.0', 'October 2017', '2007')
    description = 'SLiMSuite Core Module'
    author = 'Richard J. Edwards'
    comments = ['Please report bugs to Richard.Edwards@UNSW.edu.au']
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
            if rje.yesNo('Show general commandline options?',default='N'): out.verbose(-1,4,text=rje.__doc__)
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
### CONSTANTS ###                                                                                                     
wildcards = ['.','X','x']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMCore Class                                                                                          #
#########################################################################################################################
class SLiMCore(rje_obj.RJE_Object):
    '''
    SLiMFinder Class. Author: Rich Edwards (2007).

    Str:str
    - AAFreq = Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    - AltDis = Alternative all by all distance matrix for relationships [None]
    - BlosumFile = BLOSUM file from which to make equivalence file
    - Build = String given summary of key SLiMBuild options
    - BuildPath = Alternative path to look for existing intermediate files [SLiMFinder/]
    - CaseMask = Mask Upper or Lower case [None]
    - CompMask = Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    - DomTable = Domain table containing domain ("Type") and sequence ("Name") pairings for additional UPC [None]
    - EquivOut = File for equivalence list output [equiv.txt]
    - GablamDis = Alternative GABLAM results file [None]
    - Input = Original name (and path) of input file
    - IUScoreDir = Path in which to save protein acc.iupred.txt score files for megaslim analysis
    - MaskText = Text ID to over-ride automated masking text and identify specific masking settings [None]
    - MegaSLiM = Make/use precomputed results for a proteome (FILE) in fasta format [None]
    - MotifMask = List (or file) of motifs to mask from input sequences []
    - Prog = Program calling SLiMCore functions. Could be provided by rest server, for example.
    - RanDir = Output path for creation of randomised datasets [./]
    - Randbase = Base for random dataset name [rand_]
    - RandSource = Source for new sequences for random datasets (replaces UPCs) [None]
    - ResDir = Redirect individual output files to specified directory [SLiMFinder/]
    - RunID = Run ID for resfile (allows multiple runs on same data) [DATE]
    - SLiMDB = Path to *.slimdb file corresponding to UPC file used. []
    
    Bool:boolean
    - ConsMask = Whether to use relative conservation masking [False]
    - DisMask = Whether to mask ordered regions (see rje_disorder for options) [False]
    - DNA = Whether the sequences files are DNA rather than protein [False]
    - Force = whether to force recreation of key files [False]
    - EFilter = Whether to use evolutionary filter [True]
    - Extras = Whether to generate additional output files (alignments etc.) [False]
    - FakeMask = Whether to invoke UniFake to generate additional features for masking [False]
    - FUPC = Whether to use experimental "Fragment UPC" approach for UPC of large proteomes [False]
    - LogMask = Whether to log the masking of individual sequences [True]
    - Masked = Whether dataset has been masked [False]
    - Masking = Master control switch to turn off all masking if False [True]
    - MaskM = Masks the N-terminal M (can be useful if termini=T) [False]
    - MaskPickle = Whether to save/load pickle of masked input data, independent of main pickling [False]
    - MegaGABLAM = Whether to create and use all-by-all GABLAM results for (gablamdis) UPC generation [False]
    - MegaSLiMFix = Whether to run megaslim in "fix" mode to tidy/repair existing files [False]
    - OccStatsCalculated = Whether OccStats have been calculated for all occurrence [False]
    - MaskFreq = Whether to mask input before any analysis, or after frequency calculations [True]
    - Pickle = Whether to save/use pickles [True]
    - ProtScores = Whether to save individual protein rlc.txt files in alignment directory [False]
    - PTMData = File containing PTM data, including AccNum, ModType, ModPos, ModAA, ModCode
    - PTMod = Whether any of the sequences have been modified by PTMs
    - Randomise = Randomise UPC within batch files [False]
    - SlimBuild = Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]
    - SLiMMutant=T/F : Whether to ignore '.p.\D\d+\D' at end of accnum.
    - SmearFreq = Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    - TarGZ = Whether to tar and zip dataset result files (UNIX only) [False]
    - Test = Special Test parameter for experimentation with code [False]
    - Webserver = Generate additional webserver-specific output [False]

    Int:numeric
    - EquivCut = BLOSUM score cut-off for equivalence groups [1]
    - HomCut = Max number of homologues to allow (to reduce large multi-domain families) [0]
    - MaxSeq = Maximum number of sequences to process [500]
    - MaxUPC = Maximum UPC size of dataset to process [0]
    - MegaSLiMdp = Accuracy (d.p.) for MegaSLiM masking tool raw scores [4]
    - SaveSpace = Delete "unneccessary" files following run (see Manual for details) [0]
    - SizeSort = Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]

    Num:numeric
    - MST = MST corrected size for whole dataset
    - WallTime = Time in hours before program will terminate [1.0]

    List:list
    - AAMask = Masks list of AAs from all sequences (reduces alphabet) []
    - Alphabet = List of characters to include in search (e.g. AAs or NTs) 
    - Batch = List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    - FTMask = UniProt features to mask out [EM,DOMAIN,TRANSMEM]
    - Headers = Headers for main SLiMFinder output table
    - IMask = UniProt features to inversely ("inclusively") mask [IM]
    - PTMList = List of PTM letters to add to alphabet for analysis and restrict PTM data []
    - QRegion = Mask all but the region of the query from (and including) residue X to residue Y [0,-1]
    - SigSlim = List of significant SLiMs - matches keys to self.dict['Slim(Freq)'] - *in rank order*
    - UniprotID=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.dat and use as seqin=FILE. []
    - UP = List of UP cluster tuples
    - Warning = List of text (log) warnings to reproduce at end of run
     
    Dict:dictionary
    - AAFreq = AA frequency dictionary for each seq / UPC
    - DimFreq = Frequency of dimers of each X length per upc {upc:[freqs]}
    - Dimers = main nested dictionary for SLiMFinder {Ai:{X:{Aj:{'UP':[UPC],'Occ':[(Seq,Pos)]}}}}
    - ElementIC = dictionary of {Motif Element:Information Content}
    - Extremf. = Dictionary of {length:extremferroni correction}
    - MaskPos = Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    - MotifSeq = Dictionary of {pattern:output file for sequences}
    - MST = MST corrected size for UPC {UPC:MST}
    - Slim = main dictionary containing SLiMs with enough support {Slim:{'UPC':[UPC],'Occ':[Seq,Pos]}}
    - SeqOcc = dictionary of {Slim:{Seq:Count}}
    - SmearFreq = Smeared AA frequency dictionary for IC calculations

    Obj:RJE_Objects
    - DB = Database object. Used for storing PTMData. Might expand to use for other data output.
    - SeqList = main SeqList object containing dataset to be searched
    - SLiMCalc = SLiMCalc object used for single sequence masking
    - SlimList = MotifList object handling motif stats and filtering options

    File:FILE handles
    - MegaSLiM = Make/use precomputed results for a proteome (FILE) in fasta format [None]
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.strlist = ['AAFreq','Input','CaseMask','CompMask','ResDir','RanDir','Randbase','Build','BuildPath',
                         'AltDis','GablamDis','RandSource','DomTable','SeqIn','MaskText','EquivOut','BlosumFile',
                         'MegaSLiM','SLiMDB','RunID','IUScoreDir','Prog','PTMData']
        self.boollist = ['BuildProb','DisMask','Force','MaskFreq','Test','LogMask','MaskM','Masked','Masking','EFilter',
                        'SlimBuild','Randomise','Extras','TarGZ','SmearFreq','OccStatsCalculated','DNA','ConsMask',
                        'FakeMask','MaskPickle','Webserver','FUPC','MegaGABLAM','MegaSLiMFix','SLiMMutant','ProtScores',
                        'PTMod']
        self.intlist = ['MaxSeq','SaveSpace','MaxUPC','HomCut','Extras','SizeSort','EquivCut','MegaSLiMdp']
        self.numlist = ['StartTime','MST','WallTime']
        self.listlist = ['AAMask','Alphabet','Batch','FTMask','IMask','SigSlim','UP','Headers','Warning',
                         'PTMList','QRegion','UniprotID']
        self.dictlist = ['AAFreq','DimFreq','Dimers','Slim','SeqOcc','Extremf.','ElementIC','MST','MotifSeq','MaskPos',
                         'SmearFreq']
        self.objlist = ['DB','SeqList','SLiMCalc','SlimList']
        ### Defaults ###
        self._setDefaults(str='None',bool=True,num=0.0,int=0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
#########################################################################################################################
    def coreDefaults(self): ### Sets the default parameters used by SLiMCore functions (for easy inheritance).
        '''### Sets the default parameters used by SLiMCore functions (for easy inheritance).'''
        self.setStr({'CompMask':'5,8','ResDir':rje.makePath('SLiMCore/'),'Input':'','Basefile':'','MotifMask':'',
                      'RanDir':rje.makePath('Random/'),'Randbase':'rand','BuildPath':rje.makePath('SLiMCore/'),
                     'MaskText':'None','EquivOut':'equiv.txt','IUScoreDir':'','Prog':'','PTMData':''})
        self.setInt({'MaxSeq':0,'SaveSpace':0,'MaxUPC':0,'HomCut':0,'Extras':1,'SizeSort':0,'EquivCut':1,'MegaSLiMdp':4})
        self.setNum({'StartTime':time.time(),'WallTime':1.0})
        self.setBool({'MaskFreq':True,'Extras':True,'Masked':False,'MaskM':True,'Randomise':False,'TarGZ':False,
                     'Force':False,'SmearFreq':False,'DisMask':False,'Test':False,'OccStatsCalculated':False,'Append':True,
                     'LogMask':True,'DNA':False,'ConsMask':False,'Pickle':True,'FakeMask':False,'MaskPickle':False,
                     'Webserver':False,'FUPC':False,'MegaGABLAM':False,'MegaSLiMFix':False,'SLiMMutant':False,
                     'ProtScores':False,'PTMod':False})
        for lkey in ['AAMask','Alphabet','Batch','FTMask','IMask','SigSlim','UP','Headers','Warning','QRegion']: self.list[lkey] = []
        self.list['FTMask'] = []
        self.list['Batch'] = ['*.dat','*.fas']
        self.list['PTMList'] = []
        self.list['UniprotID'] = []
        self.dict['MaskPos'] = {2:'A'}
        self.dict['SmearFreq'] = {}
        ### Other Attributes ###
        self.obj['DB'] = rje_db.Database(self.log,self.cmd_list)
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        self.obj['SlimList'].obj['SLiMCalc'].setupFilters(slimheaders=[],occheaders=[])    ### Sets up SLiM/Occ Filters
#########################################################################################################################
    def _cmdList(self):
        for cmd in self.cmd_list:
            try:    ### SLiMCore-specific functions, not used by parents
                self._cmdReadList(cmd,'str',['Prog'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.coreCmd()
#########################################################################################################################
    def coreCmd(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'file',['AAFreq','AltDis','GablamDis','RandSource','DomTable','SeqIn','EquivOut',
                                              'MegaSLiM','BlosumFile','PTMData'])
                self._cmdReadList(cmd,'str',['CaseMask','CompMask','Randbase','MotifMask','MaskText','RunID'])
                self._cmdReadList(cmd,'path',['ResDir','RanDir','BuildPath','IUScoreDir'])
                self._cmdReadList(cmd,'int',['SaveSpace','MaxSeq','MaxUPC','HomCut','Extras','SizeSort','EquivCut','MegaSLiMdp'])
                self._cmdReadList(cmd,'num',['WallTime'])
                self._cmdReadList(cmd,'bool',['DisMask','Force','MaskFreq','Extras','Test','LogMask','Randomise','MaskM',
                                             'EFilter','TarGZ','Masking','SlimBuild','SmearFreq','DNA','ConsMask',
                                             'Pickle','FakeMask','MaskPickle','Webserver','FUPC','MegaGABLAM',
                                             'MegaSLiMFix','SLiMMutant','ProtScores'])
                self._cmdRead(cmd,'bool','MaskM','metmask')
                self._cmdRead(cmd,'bool','MegaGABLAM','megablam')
                self._cmdRead(cmd,'list','FTMask','featuremask')
                self._cmdReadList(cmd,'list',['FTMask','IMask','AAMask','Alphabet','QRegion','PTMList','UniprotID'])
                self._cmdReadList(cmd,'glist',['Batch'])
                self._cmdReadList(cmd,'cdict',['MotifSeq','MaskPos'])
                self._cmdRead(cmd,'cdict','MaskPos','posmask')
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self.setAlphabet()
        self.setupPTMs()
        self.setQRegion()
        if self.getStrLC('ResDir') == 'none': self.setStr({'ResDir':''})
        if self.getStrLC('MegaSLiM') == 'seqin' and self.getStrLC('SeqIn'): self.setStr({'MegaSLiM':self.getStr('SeqIn')})
        if self.getStrLC('MegaSLiM') and self.getBool('MegaGABLAM'):
            self.setStr({'GablamDis':"%s.gablam.tdt" % rje.baseFile(self.getStr('MegaSLiM'))})
            self.printLog('#GDIS','MegaSLiM: set GABLAMDis=%s for UPC generation.' % self.getStr('GablamDis'))
        self.list['FTMask'] = rje.listUpper(self.list['FTMask'])
        if self.list['FTMask'] == ['F'] or self.list['FTMask'] == ['FALSE']:
            self.list['FTMask'] = []
            self.warnLog('Command ftmask=F interpreted as empty FTMask list.')
        elif self.list['FTMask'] == ['T'] or self.list['FTMask'] == ['TRUE']:
            self.list['FTMask'] = ['EM','DOMAIN','TRANSMEM']
            self.warnLog('Command ftmask=T converted to default ftmask=EM,DOMAIN,TRANSMEM.')
        #else: self.printLog('#FTCMD',string.join(self.list['FTMask']))
#########################################################################################################################
    def prog(self):
        if self.getStrLC('Prog'): return self.getStrLC('Prog')
        else: return self.log.name()
#########################################################################################################################
    ### <2> ### Simple Stats Methods                                                                                    #
#########################################################################################################################
    def occFilter(self): return self.obj['SlimList'].obj['SLiMCalc'].dict['OccFilter']
    def statFilter(self): return self.obj['SlimList'].obj['SLiMCalc'].dict['SlimFilter']
    def occStats(self): return self.obj['SlimList'].obj['SLiMCalc'].list['SLiMCalc']
    def slims(self): return self.obj['SlimList'].slims()
    def dataset(self): return rje.baseFile(self.getStr('Input'),strip_path=True)
    def seqs(self): return self.obj['SeqList'].seq[0:]
    def seqNum(self): return self.obj['SeqList'].seqNum()
    def units(self): return self.obj['SeqList'].units()
    def aaNum(self): return self.obj['SeqList'].aaTotal()
    def nonX(self): return self.obj['SeqList'].aaTotal(nonx=True)
#########################################################################################################################
    def slimNum(self): return len(self.dict['Slim'])
#########################################################################################################################
    def slimUP(self,slim):  ### Returns UP Num for SLiM if in dictionary, else 0.
        if self.dict['Slim'].has_key(slim): return len(self.dict['Slim'][slim]['UP'])
        return 0
#########################################################################################################################
    def slimOccNum(self,slim,upc=None):    ### Returns number of occ of Slim in given UPC
        '''Returns number of occ of Slim in given UPC.'''
        ## No UPC - return total! ##
        if not upc: return len(self.dict['Slim'][slim]['Occ'])
        ## See if SLiM occurs at all! ##
        if not self.dict['SeqOcc'].has_key(slim): return 0
        ## Mean over UPC ##
        sx = 0
        for seq in upc:
            if self.dict['SeqOcc'][slim].has_key(seq) and self.dict['SeqOcc'][slim][seq] > sx:
                sx = self.dict['SeqOcc'][slim][seq]
        return sx
#########################################################################################################################
    def slimSeqNum(self,slim):  ### Returns the number of sequences SLiM occurs in.
        '''Returns the number of sequences SLiM occurs in.'''
        if not self.dict['SeqOcc'].has_key(slim): return 0      ## See if SLiM occurs at all! ##
        return len(self.dict['SeqOcc'][slim])
#########################################################################################################################
    def UPNum(self): return len(self.list['UP'])
#########################################################################################################################
    def getUP(self,seq):
        for upc in self.list['UP']:
            if seq in upc: return upc
        return None
#########################################################################################################################
    def slimIC(self,slim,usefreq=False):   ### Returns IC of slim  #!# Add aafreq & wildcard pen etc? 
        '''Returns IC of slim. Does not account for variable length wildcards!'''
        ic = 0.0
        slist = string.split(slim,'-')
        if usefreq and ('SmearFreq' not in self.dict or not self.dict['SmearFreq']): self.dict['SmearFreq'] = self.smearAAFreq(False)
        i = 0
        while i < (len(slist)):     #!# Update this for new calculation and DNA?! #!#
            el = slist[i]
            if not self.dict['ElementIC'].has_key(el):
                #x#self.dict['ElementIC'][el] = rje_motif.elementIC(el,callobj=self)
                if usefreq: self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,self.dict['SmearFreq'])
                elif self.getBool('DNA'): self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,rje_slim.even_ntfreq)
                else: self.dict['ElementIC'][el] = rje_slim.weightedElementIC(el,rje_slim.even_aafreq)
            ic += self.dict['ElementIC'][el]
            i += 2
        return ic
#########################################################################################################################
    ### <3> ### General Run Methods                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main Run Method. SLiMFinder, SLiMSearch etc. should have their own methods.
        '''
        Main Run Method. SLiMFinder, SLiMSearch etc. should have their own methods.
        1. Check for randomise function and execute if appropriate
        2. Input:
            - Read sequences into SeqList
            - or - Identify appropriate Batch datasets and rerun each with batch=True
        3. Masking and UPC construction.
        4. MotifSeq option if desired.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:### ~ [1] Special Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Randomise Function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('Randomise') and not batch: return self.randomise()
            ## ~ [1b] Equivalence Builder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getStrLC('BlosumFile'):
                equiv = self.BLOSUM2equiv(self.getStr('BlosumFile'),self.getNum('EquivCut'))    # getInt?
                self.printLog('#EQUIV',string.join(equiv,', '))
                rje.backup(self,self.getStr('EquivOut'))
                open(self.getStr('EquivOut'),'w').write(string.join(equiv,'\n'))
                self.printLog('#SAVE','Equivalence list saved to %s' % self.getStr('EquivOut'))
                return equiv
            ## ~ [1c] Setup MegaSLiM files for GABLAM and masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setupSeqIn()
            self.setupBasefile()    # Includes restSetup
            if not batch and self.getStrLC('MegaSLiM'): self.megaSLiM()

            ### ~ [2] Input/Batch handling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            #self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            ## ~ [2a] Batch Mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not batch and self.seqNum() < 1:   # No sequences loaded - use batch mode
                #batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
                #self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
                batchfiles = self.batchFiles()
                if not batchfiles and self.getStrLC('MegaSLiM'): self.printLog('#MEGA','MegaSLiM run finished.')
                elif not batchfiles: self.errorLog('No input files found!',printerror=False)
                else:
                    mycmd = self.cmd_list[0:]
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        bx += 1
                        self.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.run(batch=True)
                        self.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.setBool({'Append':True})
                        self.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.iStr(bx),rje.iStr(len(batchfiles)-bx)),log=False)
                if self.getBool('Win32') and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return not not batchfiles
            ## ~ [2b] Check whether to bother running dataset at all - Check Input versus Min and Max Seq ~ ##
            if 0 < self.getInt('MaxSeq') < self.obj['SeqList'].seqNum():
                self.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.iStr(self.obj['SeqList'].seqNum()),rje.iStr(self.getInt('MaxSeq'))))
                return False

            ### ~ [3] UPC and Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.setNum({'StartTime':time.time()})
            ## ~ [3a] UPC and MinOcc settings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.makeUPC(): return self.errorLog('Error during %s.makeUPC(). Abandoning %s run' % (self.prog(),self.dataset()),printerror=False)
            if 0 < self.getInt('MaxUPC') < self.UPNum():
                self.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). Run aborted.' % (self.UPNum(),self.getInt('MaxUPC')))
                return
            ## ~ [3b] AA Frequency Calculations made early as needed superficially in SLiMBuild ~~~ ##
            if self.dict['MotifSeq']:
                self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
                if self.getBool('Masked'): self.obj['SeqList'].saveFasta(seqfile='%s.%s.masked.fas' % (self.seqBaseFile(),self.maskText()))

            ### ~ [4] Special MotifSeq Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                self.motifSeq()
    
            ### End ###
                self.printLog('#CORE','SLiMCore results output to %s*' % (self.seqBaseFile()))  #['Basefile']))
            if self.getBool('Win32') and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            return True
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.getNum('WallTime') <= 0 or (time.time() - self.getNum('StartTime')) < (self.getNum('WallTime')*3600): raise
            if self.list['Headers']: self.results(aborted=True)
            return False # Walltime reached
        except:
            self.serverEnd(endcause='Crash',details=self.errorLog('Error in %s.run()' % self.prog(),printerror=True,quitchoice=False))
            return False
#########################################################################################################################
    def batchFiles(self,pickup=[]):   ### Returns batch file list
        '''Returns batch file list.'''
        ### ~ [0] ~ Standard Batch file list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
        self.printLog('\r#FILES','Getting files: %5s files for batch run' % rje.integerString(len(batchfiles)))
        if not batchfiles and self.list['Batch'] and not self.getStrLC('MegaSLiM'):
            self.errorLog('No batch files identified from %s!' % self.list['Batch'],printerror=False)
            return []
        try:
            if not self.getInt('SizeSort'): return batchfiles
        except: return batchfiles
        ### ~ [1] SizeSort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            newbatch = []       # New list of batch files in size order, starting with the biggest
            sizedict = {}
            for file in batchfiles:
                if os.path.split(rje.baseFile(file))[1] in pickup:
                    self.printLog('#PICKUP','Skipping %s batch file %s' % (self.getStr('RunID'),file))
                    continue
                fsize = os.path.getsize(file)
                if fsize not in sizedict: sizedict[fsize] = []
                sizedict[fsize].append(file)
            sizes = rje.sortKeys(sizedict,revsort=self.getInt('SizeSort')<0)
            for fsize in sizes: newbatch += sizedict.pop(fsize)
        except: self.errorLog('Problem with %s.batchFiles() SizeSort' % self.prog()); return batchfiles
        return newbatch
#########################################################################################################################
    def wallTime(self):     ### Exits if walltime has been reached
        '''Exits if walltime has been reached.'''
        if self.getNum('WallTime') <= 0 or (time.time() - self.getNum('StartTime')) < (self.getNum('WallTime')*3600): return
        self.errorLog('%s Walltime (%.2f hours) reached! Try increasing WallTime=X.' % (self.prog(),self.getNum('WallTime')),printerror=False)
        self.close()
        self.serverEnd('Walltime',exit=True)
        sys.exit()
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMCore object for new batch run
        '''Returns SLiMCore object for new batch run.'''
        return SLiMCore(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.getBool('Append')])
#########################################################################################################################
    def serverEnd(self,endcause,details=None,exit=True):  ### Adds #END statement for webserver
        '''
        Adds #END statement for webserver.
        >> endcause:str = Identifier for ending program run.
        >> details:str [None] = Extra information for certain end causes only.
        >> exit:bool [True] = whether to terminate program after #END statement.
        '''
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        try:
            dtxt = 'Download install a local copy of SLiMFinder <a href="http://www.southampton.ac.uk/~re1u06/software/slimsuite/">here</a> for larger runs.'
            if not self.getBool('Webserver') and not self.getStrLC('Rest'):
                if endcause == 'Crash': sys.exit(1)
                else: return
            ### ~ [2] END statements ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if endcause == 'Walltime': end = '%s server Walltime (%.2f hours) reached! Suggestions: reduce dataset size, increase masking, reduce slimlen and/or maxwild, switch ambiguity off.' % (self.prog(),self.getNum('WallTime'))
            elif endcause == 'MaxSeq': end = 'Server limit on unrelated sequence (maxseq=%s) exceeded: %s input seq. Analysis terminated.' % (rje.iStr(self.getInt('MaxSeq')),rje.integerString(self.seqNum()))
            elif endcause == 'MaxSize': end = 'Server limit on size of dataset (maxsize=%s %s) exceeded: %s post-masking AA. Analysis terminated.' % (rje.integerString(self.getInt('MaxSize')),self.units(),rje.integerString(self.aaNum()))
            elif endcause == 'FewSeq': end = 'Insufficient sequences for minimum SLiM occurrence setting (minocc=%d): %d input sequences; %d+ unrelated sequences needed. Run aborted.' % (self.getInt('MinOcc'),self.seqNum(),self.getInt('MinOcc'))
            elif endcause == 'MaxUPC': end = 'Server limit on unrelated sequence (maxupc=%s) exceeded: %s unrelated input sequences (UPC). Analysis terminated.' % (rje.integerString(self.getInt('MaxUPC')),rje.integerString(self.UPNum()))
            elif endcause == 'FewUPC': end = 'Insufficient <i>unrelated</i> sequences for minimum SLiM occurrence setting (minocc=%d): %d unrelated input sequences (UPC); %d+ unrelated sequences needed. Run aborted.' % (self.getInt('MinOcc'),self.UPNum(),self.getInt('MinOcc'))
            elif endcause == 'NoSLiM': end = 'No SLiMs for output. Try relaxing probability threshold (probcut) and setting number of SLiMs to return (topranks) > 0.'
            elif endcause == 'NoMotif': end = 'No SLiMs retained for searching. Check format of input file and information content of motifs.'
            elif endcause == 'Crash':
                if details: end = 'Sorry - %s crashed during %s! Please check log or webmaster with job ID.' % (self.prog(),details)
                else: end = 'Sorry - %s crashed! Please check log or contact webmaster with job ID.' % self.prog()
            ### ~ [3] Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else: end = '%s server run ended (%s): see log file for details.' % (self.prog(),endcause)
            if endcause in ['Walltime','MaxSeq','MaxUPC']: end = '%s %s' % (end,dtxt)
            self.printLog('#END',end)
            #self.dict['Output']['status'] = end
            if exit:
                if self.getStrLC('Rest'): return
                elif endcause == 'Crash': sys.exit(1)
                else: sys.exit(0)
        except SystemExit: raise
        except: self.errorLog('ServerEnd error.')
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def setupSeqIn(self):   ### Sets up SeqIn, using UniprotID list if appropriate
        '''Sets up SeqIn, using UniprotID list if appropriate.'''
        seqcmd = ['gnspacc=T','usecase=T','accnr=F','seqnr=F'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
        if self.getStrLC('Basefile'): datfile = '%s.dat' % self.baseFile(runpath=True)
        else: datfile = '%s%s.dat' % (self.getStr('RunPath'),self.prog().lower())
        if self.list['UniprotID'] and not rje.exists(self.getStr('SeqIn')):
            unicmd = self.cmd_list + ['datout=%s' % datfile]
            uni = rje_uniprot.UniProt(self.log,unicmd)
            uni.run()
            self.dict['Output']['map'] = rje.baseFile(datfile) + '.map.tdt'
            mfile = rje.baseFile(datfile) + '.missing.acc'
            if rje.exists(mfile): self.dict['Output']['failed'] = mfile
            self.setStr({'SeqIn':datfile})
            seqcmd += ['seqin=%s' % datfile]
        #!# Add proteome=X input (into buildpath: warn about this output location #!#
        seqlist = self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
        if rje.exists(seqlist.name()): self.dict['Output']['seqin'] = seqlist.name()
        else: self.dict['Output']['seqin'] = seqlist.fasta()
        if self.getStrLC('MegaSLiM') == 'seqin':
            seqfile = self.getStr('SeqIn')
            if not open(seqfile,'r').readline()[:1] == '>':
                self.printLog('#CONV','Converting %s to fasta format for MegaSLiM compatibility' % seqfile)
                newfile = rje.baseFile(seqfile) + '.fas'
                seqlist.saveFasta(seqfile=newfile)
                seqfile = newfile
            self.setStr({'MegaSLiM':seqfile})
            if self.getBool('MegaGABLAM'): self.setStr({'GablamDis':"%s.gablam.tdt" % rje.baseFile(self.getStr('MegaSLiM'))})
#########################################################################################################################
    def setupBasefile(self):    ### Sets up self.info['Basefile'] and self.info['Input']
        '''Sets up self.info['Basefile'].'''
        ### ~ [0] ~ Store original input filename ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.getStr('Input'): self.setStr({'Input': self.obj['SeqList'].name()})
        ## ~ [0a] ~ Default RunID: used for pickup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        t = time.localtime(self.log.stat['StartTime'])
        self.setStr({'Date': '%s-%s-%s' % (str(t[0])[-4:],rje.preZero(t[1],12),rje.preZero(t[2],31))})
        if not self.getStrLC('RunID'): self.setStr({'RunID': self.getStr('Date')})
        ### ~ [1] ~ Basefile and ResDir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.seqNum(): rje.mkDir(self,self.getStr('ResDir'))
        if not self.baseFile(): self.baseFile(self.prog().lower())
        if not self.getStrLC('ResFile'): self.setStr({'ResFile':'%s.csv' % self.baseFile()})
        self.restSetup()
#########################################################################################################################
    def seqBaseFile(self,resdir='ResDir'):  ### Returns the results directory basefile as specified by sequence alone (for SLiMSearch)
        '''Returns the results directory basefile as specified by sequence alone (e.g. for SLiMSearch)'''
        return self.getStr(resdir) + self.dataset()
#########################################################################################################################
    def setupPTMs(self):    ### Sets up PTM data and modifies Alphabet.
        '''Sets up PTM data and modifies Alphabet.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #self.debug('SetupPTMs')
            #if not self.dev() and (self.getStrLC('PTMData') or self.list['PTMList']): return self.printLog('#DEV','PTMList/PTMData development only (set dev=T)')
            if not (self.getStrLC('PTMData') or self.list['PTMList']): return
            ## ~ [0a] ~ Load PTMData ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pdb = None
            if self.getStrLC('PTMData'):
                pdb = self.db().addTable(self.getStr('PTMData'),mainkeys=['accnum','modtype','modpos'],datakeys=['accnum','modtype','modpos','modaa','modcode'],name='PTMData',uselower=True)
            ### ~ [1] ~ Process PTMData ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ If PTMData loaded, restrict using PTMList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if pdb:
                ex = pdb.entryNum()
                if self.list['PTMList']:
                    self.printLog('#PTM','Filter PTMData (type/code) using PTMList: %s' % string.join(rje.sortUnique(self.list['PTMList']),'|'))
                for entry in pdb.entries()[0:]:
                    if not entry['modtype'] or not entry['modcode'] or not entry['modaa'] or not entry['modpos']: pdb.dropEntry(entry)
                    elif self.list['PTMList']:
                        if entry['modtype'] not in self.list['PTMList'] and entry['modcode'] not in self.list['PTMList']: pdb.dropEntry(entry)
                self.printLog('#PTM','%s PTMData entries reduced to %s.' % (rje.iStr(ex),rje.iStr(pdb.entryNum())))
                pdb.index('accnum')
                pdb.indexReport('modtype','#PTM'); pdb.indexReport('modcode','#PTMX')
            ## ~ [1b] ~ Use PTMData to update PTMList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ptmcodes = []
            for ptm in self.list['PTMList']:
                if pdb and ptm in pdb.index('modtype'):
                    for code in pdb.indexDataList('modtype',ptm,'modcode'): ptmcodes.append(code)
                else: ptmcodes.append(ptm)
            ptmcodes = rje.sortUnique(ptmcodes)
            if self.list['PTMList']: self.printLog('#PTM','%s: %s' % (string.join(self.list['PTMList'],'|'),string.join(ptmcodes)))
            self.list['PTMList'] = ptmcodes
            ## ~ [1c] ~ Use PTMList to update Alphabet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ptx = 0
            for ptm in ptmcodes:
                if ptm not in self.list['Alphabet']: self.list['Alphabet'].append(ptm); ptx += 1
            self.list['Alphabet'] = rje.sortUnique(self.list['Alphabet'])
            if ptx: self.printLog('#ALPH','PTM-Adjusted Alphabet: %s' % string.join(self.list['Alphabet']))
        except: self.errorLog('%s.setupPTMs failure' % self.prog())
#########################################################################################################################
    def setAlphabet(self):  ### Sets up self.list['Alphabet']
        '''Sets up self.list['Alphabet'].'''
        ### ~ [1] ~ Setup basic alphabet ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if len(self.list['Alphabet']) == 1: self.list['Alphabet'] = string.split(string.join(self.list['Alphabet'][0]).upper())
        elif self.list['Alphabet']: self.list['Alphabet'] = string.split(string.join(self.list['Alphabet']).upper())
        else:
            if self.getBool('DNA'): self.list['Alphabet'] = rje_seq.alph_dna[:-1]
            else: self.list['Alphabet'] = rje_seq.alph_protx[:-1]
        self.list['Alphabet'] = rje.sortUnique(self.list['Alphabet'])
        self.printLog('#ALPH','Alphabet: %s' % string.join(self.list['Alphabet']))
        ### ~ [2] ~ Remove masked AAs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for aa in self.list['AAMask']:
            if aa in self.list['Alphabet']: self.list['Alphabet'].remove(aa)
            self.printLog('#ALPH','AA masked alphabet: %s' % string.join(self.list['Alphabet']))
        self.dict['Output']['alphabet'] = string.join(self.list['Alphabet'])
#########################################################################################################################
    def setQRegion(self):   ### Sets up Query Region masking
        '''Sets up Query Region masking.'''
        if 'QRegion' not in self.list: self.list['QRegion'] = []
        try:
            qreg = []; pairs = True
            for x in self.list['QRegion']:
                r = string.atoi(x)
                qreg.append(r)
                pairs = not pairs
            if not pairs: qreg.append(-1)
            self.list['QRegion'] = qreg
        except:
            self.errorLog('Problem with QRegion masking %s - will not use' % self.list['QRegion'])
            self.list['QRegion'] = []
        if self.list['QRegion']: self.printLog('#QREG','%d Query region(s)' % int((len(self.list['QRegion'])+0.5)/2))
#########################################################################################################################
    def maskInput(self):    ### Masks input sequences, replacing masked regions with Xs
        '''Masks input sequences, replacing masked regions with Xs. Creates seq.info['PreMask' & 'MaskSeq']'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('FUPC'): self.warnLog('Masking not yet implemented for FUPC',quitchoice=True)
            if self.maskPickle(load=True): return True
            masked = False
            if not self.getBool('Masking'):
                #for seq in self.seqs(): seq.info['PreMask'] = seq.info['MaskSeq'] = seq.info['Sequence'][0:]
                for seq in self.seqs(): seq.info['PreMask'] = seq.info['Sequence'][0:]
                self.setBool({'Masked':False})
                return True
            ## ~ [0a] ~ UniFake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('FakeMask'):
                fake = unifake.UniFake(self.log,self.cmd_list)
                fake.setup()
                fake.obj['SeqList'] = self.obj['SeqList']
                fake.info['DatOut'] = '%s.dat' % self.seqBaseFile() #['Basefile']
                fake.uniFake(store=True)
                try: self.printLog('#FAKE','%d UniFake entries added for masking' % len(self.obj['SeqList'].obj['UniProt'].list['Entry']))
                except: self.errorLog('UniFake Masking error')
                self.wallTime()
            ### ~ [1] ~ Case Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                seq.info['PreMask'] = seq.info['Sequence'][0:]
                if seq.maskCase(case=self.getStr('CaseMask'),mask='X',log=self.getBool('LogMask')): masked = True
                self.wallTime()
            ### ~ [2] ~ Disorder (Inclusive Filtering) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('DisMask'):
                if self.getStrLC('MegaSLiM'): self.megaSLiMMaskDisorder()
                elif self.getStrLC('IUScoreDir'):
                    self.setStr({'MegaSLiM':self.getStr('SLiMDB')})     # Generate/reuse IUScores for all proteins
                    self.megaSLiMMaskDisorder()
                else:
                    for seq in self.seqs():
                        seq.maskDisorder(inverse=True,log=self.getBool('LogMask'),ikey='PreMask')
                        #seq.maskDisorder(inverse=True,log=self.getBool('LogMask'))
                        self.wallTime()
                masked = True
                self.wallTime()
            ### ~ [3] ~ Relative conservation masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('ConsMask'):
                if self.getStrLC('MegaSLiM'): self.megaSLiMMaskConservation()
                elif self.getStrLC('IUScoreDir') or self.getBool('ProtScores'):
                    self.setStr({'MegaSLiM':self.getStr('SLiMDB')})     # Generate/reuse IUScores for all proteins
                    self.megaSLiMMaskConservation()
                else:
                    slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc=','usealn=T'])
                    slimcalc.loadBLOSUM()
                    if slimcalc.getBool('UseGopher') and not slimcalc.getStr('AlnDir') and slimcalc.obj['Gopher'].getInt('Forks') > 1 and not slimcalc.obj['Gopher'].getBool('NoForks'):
                        slimcalc.obj['Gopher'].setStr({'Name':self.getStr('SLiMDB')})
                        slimcalc.obj['Gopher'].cmd_list.append('orthaln')
                        slimcalc.obj['Gopher'].run()
                    cmx = 0
                    for seq in self.seqs():
                        (seq.info['MaskSeq'],seq.info['Sequence']) = (seq.info['Sequence'],seq.info['PreMask'])
                        seq.deGap()
                        if self.getBool('SLiMMutant') and rje.matchExp('^(\S+)\.p\.\D\d+\D$',seq.info['AccNum']):
                            oldacc = seq.info['AccNum']
                            seq.info['AccNum'] = rje.matchExp('^(\S+)\.p\.\D\d+\D$',seq.info['AccNum'])[0]
                            relcon = slimcalc.relConListFromSeq(seq,window=30)
                            seq.info['AccNum'] = oldacc
                        else: relcon = slimcalc.relConListFromSeq(seq,window=30)
                        newseq = []
                        try:
                            for i in range(seq.seqLen()):
                                if relcon[i] < 0: newseq.append('X')
                                else: newseq.append(seq.info['MaskSeq'][i])
                            maskx = newseq.count('X') - seq.info['MaskSeq'].count('X'); cmx += maskx
                            seq.info['Sequence'] = string.join(newseq,'')
                            if maskx > 0 and self.getBool('LogMask'): self.printLog('#MASK','Masked %s low relative conservation. (%d X added to %daa seq.)' % (seq.shortName(),maskx,seq.aaLen()))
                        except:
                            self.errorLog('Problem with relative conservation masking for %s' % seq.shortName())
                            self.printLog('#REL','%daa (%d) => %d conscores' % (seq.seqLen(),seq.aaLen(),len(relcon)))
                            seq.info['Sequence'] = seq.info['MaskSeq']
                        self.wallTime()
                    self.printLog('#REL','%s aa masked using relative conservation' % rje.integerString(cmx))
                masked = True
            ### ~ [4] ~ PTMData Modifications ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.db('PTMData'):   # Table of ['accnum','modtype','modpos','modaa','modcode']
                pdb = self.db('PTMData')
                pdb.dataFormat({'modpos':'int'})
                moddict = {}    # {(sname,pos):mod} dictionary
                for seq in self.seqs():
                    accnum = seq.info['AccNum']; sname = seq.shortName()
                    if accnum not in pdb.index('accnum'): self.printLog('#PTM','No PTM modifications for %s.' % sname)
                    for entry in pdb.indexEntries('accnum',accnum):
                        pos = entry['modpos'] - 1
                        aa = seq.info['PreMask'][pos].upper()
                        if (sname,pos) in moddict:
                            self.warnLog('Conflicting %s (%s) at %s Pos %s; Already modified: %s' % (entry['modtype'],entry['modcode'],sname,pos+1,moddict[(sname,pos)]))
                        if seq.info['PreMask'][pos] != entry['modaa']:
                            self.warnLog('%s position %d should be modified %s but %s found in sequence.' % (sname,entry['modpos'],entry['modaa'],aa),'bad modpos',suppress=True)
                            continue
                        moddict[(sname,pos)] = '%s (%s%s%s)' % (entry['modtype'],entry['modaa'],entry['modpos'],entry['modcode'])
                        self.setBool({'PTMod':True})
                        #self.debug(moddict[(sname,pos)])
                        for stype in ['PreMask','Sequence']:    # No need to adjust 'MaskSeq'
                            #self.bugPrint('%s (%s aa)' % (stype,len(seq.info[stype])))
                            if seq.info[stype][pos].upper() == aa: seq.info[stype] = seq.info[stype][:pos] + entry['modcode'] + seq.info[stype][pos+1:]
                            #self.debug(seq.info[stype])
                for (sname,pos) in rje.sortKeys(moddict): self.printLog('#PTMOD','%s %s' % (sname,moddict[(sname,pos)]))
            ### ~ [5] ~ UniProt Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['FTMask']: ftmasktxt = 'FTMask: %s' % string.join(self.list['FTMask'],',')
            else: ftmasktxt = 'FTMask: None'
            if self.list['IMask']: ftmasktxt += '; Inclusive Masking: %s' % string.join(self.list['IMask'],',')
            if self.obj['SeqList'].obj['UniProt']:
                self.printLog('#FTMASK',ftmasktxt)
                for entry in self.obj['SeqList'].obj['UniProt'].list['Entry']:
                    if self.list['IMask']:
                        entry.maskFT(types=self.list['IMask'],inverse=True,mask='X',log=self.getBool('LogMask'))
                        masked = True
                    if self.list['FTMask']:
                        entry.maskFT(types=self.list['FTMask'],inverse=False,mask='X',log=self.getBool('LogMask'))
                        masked = True
                    self.wallTime()
            elif self.list['FTMask'] or self.list['IMask']: self.warnLog('No feature based masking for non-Uniprot input (%s)' % ftmasktxt)
            ### ~ [6] ~ Low Complexity Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.matchExp('(\d+),(\d+)',self.getStr('CompMask')):
                (x,y) = rje.matchExp('(\d+),(\d+)',self.getStr('CompMask'))
                for seq in self.seqs(): seq.maskLowComplexity(lowfreq=int(x),winsize=int(y),log=self.getBool('LogMask'))
                masked = True
                self.wallTime()
            elif self.getStr('CompMask').lower()[:1] not in ['f','n']:
                self.errorLog('CompMask "%s" wrong format for complexity masking' % self.getStr('CompMask'),printerror=False)
            ### ~ [7] ~ N-terminal Met ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MaskM'):
                mx = 0
                for seq in self.seqs():
                    if seq.info['Sequence'][0] == 'M':
                        seq.info['Sequence'] = 'X' + seq.info['Sequence'][1:]
                        mx += 1
                self.printLog('#MASK','MetMask: %d N-terminal Ms masked' % mx)
                self.wallTime()
            ### ~ [8] ~ Position-specific masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.dict['MaskPos']:
                for pos in rje.sortKeys(self.dict['MaskPos']):
                    try: (x,y) = (int(pos),rje.strList(self.dict['MaskPos'][pos]))
                    except: self.errorLog('Problem with MaskPos entry "%s:%s"' % (pos,self.dict['MaskPos'].pop(pos)))
                for seq in self.seqs(): seq.maskPosAA(self.dict['MaskPos'],mask='X',log=self.getBool('LogMask'))
                masked = True; self.wallTime()
            ### ~ [9] ~ Motif masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStr('MotifMask'):
                maskslims = rje_slimlist.SLiMList(self.log,self.cmd_list)
                maskslims.loadMotifs(self.getStr('MotifMask'))
                mx = 0
                for seq in self.seqs():
                    sx = -seq.info['Sequence'].count('X')
                    for slim in maskslims.slims():
                        for hit in slim.searchSequence(sequence=seq.info['PreMask']):    # {Pos,Variant,Match,ID,MisMatch}
                            #self.deBug(hit)
                            try:
                                r = hit['Pos'] - 1
                                for i in range(len(hit['Variant'])):
                                    if hit['Variant'][i] != '.': seq.info['Sequence'] = rje.strSub(seq.info['Sequence'],r+i,r+i,'X')
                            except: self.errorLog('Problem with motifmask %s' % hit)
                    sx += seq.info['Sequence'].count('X')
                    if sx > 0:
                        mx += sx
                        if self.getBool('LogMask'): self.printLog('#MASK','Motif-Masked %s. (%d X added to %daa seq.)' % (seq.shortName(),sx,seq.aaLen()))
                    self.wallTime()
                self.printLog('#MASK','MotifMask: %d motif match residues masked' % mx)
            ### ~ [10] ~ AA Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['AAMask']:
                mx = 0
                for seq in self.seqs(): mx += seq.maskAA(self.list['AAMask'],log=self.getBool('LogMask'))
                self.printLog('#MASK','AA-Masked %s: %s X added.' % (string.join(self.list['AAMask'],';'),rje.integerString(mx)))
                self.wallTime()
            ### ~ [10] ~ QRegion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #try:
            #    if self.dev():
            #        self.bugPrint('Development debugging of QRegion...')
            #        self.deBug(self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus'])
            #        self.deBug(self.list['QRegion'])
            #        self.deBug('Focus' in self.dict)
            #        self.deBug('Query' in self.dict['Focus'])
            #except: pass
            if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
                mx = 0
                for seq in self.dict['Focus']['Query']: mx += seq.maskRegion(self.list['QRegion'],inverse=True)
                self.printLog('#MASK','QRegion %s: %s X added.' % (self.list['QRegion'],rje.integerString(mx)))
                self.wallTime()
            ### ~ [11] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs(): seq.info['MaskSeq'] = seq.info['Sequence'][0:]
            if masked:
                self.printLog('#MASK','%s: Masking of input sequences complete.' % self.dataset())
                self.setBool({'Masked':True})
            #self.deBug(self.obj['SeqList'].seq[0].info)
            self.maskPickle()
            return True
        except: self.errorLog('%s.maskInput error!' % self.prog()); raise
#########################################################################################################################
    def maskPickle(self,load=False):    ### Loads or Saves post-masking pickle
        '''Loads or Saves post-masking pickle.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('MaskPickle'): return False
            #maskpickle = '%s.mask.%s' % (rje.baseFile(self.info['Input'],strip_path=True),self.maskText(freq=False))  ## Pickle name ##
            maskpickle = '%s.mask.%s' % (self.dataset(),self.maskText(freq=False))  ## Pickle name ##

            ### ~ [1] Load existing pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load:
                if self.force(): return None
                newme = None    ## New object, loaded from Pickle
                ## Check for file and load ##
                for pdir in [self.getStr('ResDir'),self.getStr('BuildPath')]:
                    pfile = '%s%s' % (pdir,maskpickle)
                    newme = self.unpickleMe(basefile=pfile,process=False)
                    if newme: break
                if not newme: return None
                ## Check for masking differences ##
                changes = []
                for var in ['CompMask','CaseMask','MotifMask']:         # Info
                    if self.getStr(var) != newme.getStr(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['Masking','DisMask','MaskM','ConsMask']:   # Opt
                    if self.getBool(var) != newme.getBool(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['FTMask','IMask']:                          # List   
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
                newme.setInfo(self.str)    #!# Will need to change these when switching to new Object template
                newme.setInt(self.int)
                newme.setNum(self.num)
                self.setBool({'Masked':newme.getBool('Masked')})
                newme.setBool(self.bool)
                newme.setLog(self.log)
                for o in self.obj:
                    if o != 'SeqList': newme.obj[o] = self.obj[o]
                self.list['UP'] = newme.list['UP']
                self.dict['AAFreq'] = newme.dict['AAFreq']
                self.dict['MST'] = newme.dict['MST']
                newme.list = self.list
                newme.dict = self.dict
                self.replaceMe(newme)
                return True
                
            ### ~ [2] Save new pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.maskPickleMe(basefile='%s%s' % (self.getStr('ResDir'),maskpickle),gzip=True,replace=True)
            return None
        except: self.errorLog('Problem with %s.maskPickle()' % self.prog()); return False
#########################################################################################################################
    def maskPickleMe(self,basefile=None,gzip=True,replace=True):   ### Saves self object to pickle and zips
        '''
        Saves self object to pickle and zips.
        >> basefile:str [None] = if none, will use self.getStr('Basefile')
        >> gzip:bool [True] = whether to GZIP (win32=F only)
        >> replace:bool [True] = whether to replace existing Pickle
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Pickle'): self.printLog('#PICKLE','No pickling. (pickle=F)'); return False
            if not basefile or basefile.lower() == 'none': basefile = self.basefile()
            pfile = '%s.pickle' % basefile
            if self.getBool('Webserver'): pfile = '%s%s' % (self.getStr('BuildPath'),os.path.basename(pfile))
            if not replace and (os.path.exists(pfile) or os.path.exists('%s.gz' % pfile)):
                self.printLog('#PICKLE','No pickling - pickle already exists (no replace)!'); return False
            ### ~ [2] ~ Pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('#SAVE','Attempting to save to %s.' % pfile,log=False)
            pickle.dump(self,open(pfile,'w'))
            self.printLog('#SAVE','Intermediate saved as %s (Python pickle).' % pfile)
            ### ~ [3] ~ GZip and finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getBool('Win32') and gzip:
                try:
                    if os.path.exists('%s.gz' % pfile): os.unlink('%s.gz' % pfile)
                    os.system('gzip %s' % pfile)
                    self.printLog('#GZIP','%s zipped.' % pfile)
                except: self.errorLog('Cannot gzip %s' % pfile)
            return True
        except: self.errorLog('Problem during %s.maskPickleMe()' % self); return False
#########################################################################################################################
    def searchPickleMe(self,load=False):    ### Loads existing search pickle or save pickle for later
        '''Loads existing search pickle or save pickle for later.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load and self.force(): return None      # Re-run search!
            ### ~ [2] Load existing pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if load:    # Subsequent processing of pickle now handled by processPickle()
                newme = self.unpickleMe(self.searchBase(resdir='ResDir'))
                if not newme: newme = self.unpickleMe(self.searchBase(resdir='BuildPath'))
                return newme
            ### ~ [3] Save new pickle ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.pickleMe(self.searchBase())
            return None
        except:
            self.errorLog('Major problem with SLiMProb pickling!')
            return None
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
            newme.setStr(self.str)        #!# Will need to update when switching to new Object
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
    def makeUPC(self):  ### Generates UP Clusters from self.obj['SeqList'] using BLAST
        '''Generates UP Clusters from self.obj['SeqList'] using BLAST.'''
        try:### ~ [1] Check whether BLAST and MST is necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            upcbase = '%s%s' % (self.getStr('ResDir'),rje.baseFile(self.getStr('Input'),strip_path=True))
            seqfile = '%s.slimdb' % upcbase
            ## ~ [1a] No evolutionary filtering (not a good idea generally!) ~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.getBool('EFilter'):
                self.list['UP'] = []
                for seq in self.obj['SeqList'].seq: self.list['UP'].append((seq,))
                self.makeMST()
                if self.seqNum() > 1: self.printLog('#EVOL','WARNING! No evolutionary filtering (efilter=F)')
                self.setStr({'SLiMDB':seqfile})
                return True
            ## ~ [1b] Direct reading of *.upc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.force() and self.readUPC(): return True
            if self.getBool('FUPC'): return self.makeFUPC()
            self.list['UP'] = []
            self.dict['MST'] = {}
            self.dict['Output']['upc'] = '%s.upc' % upcbase
            ## ~ [1c] Single sequence datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.seqNum() < 2:
                for seq in self.obj['SeqList'].seq: self.list['UP'].append((seq,))
                self.makeMST(); gentxt = 'Single sequence'
                UFILE = open('%s.upc' % upcbase,'w')
                UFILE.write('#%s# %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))
                rje.writeDelimit(UFILE,['UP','N','MST','Seqs'],'\t')
                for u in range(self.UPNum()):
                    upc = self.list['UP'][u]
                    seqs = upc[0].shortName()
                    for seq in upc[1:]: seqs += ' %s' % seq.shortName()
                    rje.writeDelimit(UFILE,[u+1,len(upc),'%.3f' % (len(upc)*self.dict['MST'][upc]),seqs],'\t')
                UFILE.close()
                self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))
                return True
            ### ~ [2] Setup BLAST etc - necessary for UPC and/or Focal Group MST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #for seq in self.obj['SeqList'].seq: self.debug(seq.info)
            gentxt = ''
            self.dict['Output']['slimdb'] = '%s.slimdb' % upcbase
            self.obj['SeqList'].setStr({'Name':seqfile})
            self.setStr({'SLiMDB':seqfile})
            if self.force() or not rje.exists(seqfile): self.obj['SeqList'].saveFasta()
            seqdict = self.obj['SeqList'].seqNameDic()
            gablam = rje_dismatrix_V3.DisMatrix(self.log,self.cmd_list)
            ## ~ [2a] Read distance matrix if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('MegaSLiM') and self.getBool('MegaGABLAM'): self.megaSLiMGABLAM(seqfile)
            readdis = False
            gablamdis = False
            mfiles = ['%s%s.dis.tdt' % (self.getStr('ResDir'),rje.baseFile(self.getStr('Input'),strip_path=True)),
                      '%s%s.dis.tdt' % (self.getStr('BuildPath'),rje.baseFile(self.getStr('Input'),strip_path=True))]
            if self.force(): mfiles = []
            if self.getStrLC('AltDis'): mfiles.append(self.getStr('AltDis'))
            if self.getStrLC('GablamDis'): mfiles.append(self.getStr('GablamDis'))
            for dismatrix in mfiles:
                if os.path.exists(dismatrix):
                    if dismatrix == self.getStr('GablamDis'):
                        gablam.loadFromDataTable(dismatrix,objlist=seqdict.keys())
                        gentxt = 'GABLAM DisMat %s' % dismatrix
                    else:
                        gablam.loadMatrix(dismatrix,checksym=False,default=1.0)
                        gentxt = 'DisMat %s' % dismatrix
                    self.printLog('#UPC','Loaded distance matrix from %s' % dismatrix)
                    if dismatrix == self.getStr('GablamDis'): gablamdis = True
                    self.dict['Output']['dismatrix'] = dismatrix
                    break
            if gablam.dict['Matrix']:   # Something read... convert
                updict = {}
                readdis = True
                newmatrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)
                (cx,ctot) = (0.0,self.seqNum())
                for seq1 in self.seqs():
                    self.progLog('\r#UPC','Converting distance matrix: %.1f%%' % (cx/ctot))
                    cx += 100.0
                    updict[seq1] = []
                    if seq1.shortName() not in gablam.dict['Matrix']: continue
                    for upseq in gablam.dict['Matrix'][seq1.shortName()]:   #x#for seq2 in self.seqs():
                        if upseq not in seqdict: continue
                        seq2 = seqdict[upseq]
                        try: newmatrix.addDis(seq1,seq2,gablam.dict['Matrix'][seq1.shortName()][seq2.shortName()])
                        except:
                            if not gablamdis: 
                                self.errorLog('Error converting distance matrix')
                                readdis = False
                                break
                        if newmatrix.getDis(seq1,seq2,default=1.0) < 1 or newmatrix.getDis(seq2,seq1,default=1.0) < 1: updict[seq1].append(seq2)
                    if not readdis:
                        del newmatrix
                        break
                self.printLog('\r#UPC','Convertion of distance matrix complete.',log=False)

            blast = rje_blast.blastObj(self.log,cmd_list=['blastcf=F','blastf=T','blaste=1e-4']+self.cmd_list,type='New')
            savedis = True
            if readdis: gablam.dict['Matrix'] = newmatrix.dict['Matrix']; savedis = False
            ## ~ [2b] BLAST generate matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif not blast.getBool('OldBLAST'): # and not self.test(): #self.dev():    # Try new BLAST+ clustering and dismatrix UPC generation
                # Make GABLAM Distance Matrix from BLAST 
                gablam = blast.blastIDMatrix(seqfile,seqdict,dna=self.getBool('DNA'),keepblast=self.extras(2))
                gentxt = 'blaste=%s, blastcf=%s, blastf=%s' % (rje_slim.expectString(blast.getNum('E-Value')),str(blast.getBool('Composition Statistics'))[0],str(blast.getBool('Complexity Filter'))[0])
                # UPC from BLAST Clusters
                updict = {}
                for bclust in gablam.cluster(maxdis=1.0,singletons=True):
                    for seq in bclust:
                        if seq not in updict: updict[seq] = []
                        for seq2 in bclust:
                            if seq2 not in updict[seq]: updict[seq].append(seq2)
                for seq in self.seqs():
                    if seq not in updict: self.errorLog('No sequence for hit "%s"!' % seq.shortName(),printerror=False); updict[seq] = [seq]
            else:
                self.warnLog('Using old BLAST (oldblast=T). Consider using BLAST+ (oldblast=F)')
                blast = rje_blast.blastObj(self.log,cmd_list=['blastcf=F','blastf=T','blaste=1e-4']+self.cmd_list,type='Old')
                #self.deBug(blast)
                blast.setInfo({'InFile':seqfile,'DBase':seqfile,'Name':'%s.self.blast' % upcbase})
                blast.setStat({'OneLine':self.seqNum(),'HitAln':self.seqNum()})
                if self.getBool('DNA'): blast.info['Type'] = 'blastn'
                blast.formatDB(fasfile=seqfile,force=self.force(),protein=not self.getBool('DNA'))
                if not rje_blast.checkForDB(dbfile=seqfile,checkage=False,log=self.log,protein=not self.getBool('DNA'),oldblast=True):
                    self.errorLog('FormatDB failed for unknown reasons.',printerror=False)
                    raise IOError
                ### Perform and read BLAST ###
                blast.blast(cleandb=True,use_existing=not self.force())
                if not blast.readBLAST(gablam=True,unlink=not self.extras(2)):
                    self.errorLog('Major problem with BLAST for unknown reasons.')
                    raise IOError
                gentxt = 'blaste=%s, blastcf=%s, blastf=%s' % (rje_slim.expectString(blast.stat['E-Value']),str(blast.opt['Composition Statistics'])[0],str(blast.opt['Complexity Filter'])[0])
                ### Make GABLAM Distance Matrix from BLAST ###
                updict = {}
                for seq in self.seqs(): gablam.addDis(seq,seq,0.0)      # Always self-distance of zero, even without BLAST hits
                for search in blast.search:
                    seq = seqdict[search.info['Name']]
                    updict[seq] = []
                    for hit in search.hit:
                        if seqdict.has_key(hit.info['Name']):
                            updict[seq].append(seqdict[hit.info['Name']])
                            if seqdict[hit.info['Name']] != seq:
                                gablam.addDis(seq,seqdict[hit.info['Name']],1 - (hit.dict['GABLAM']['Query']['GABLAM ID'] / float(seq.aaLen())))
                        else: self.errorLog('No sequence for hit "%s"!' % hit.info['Name'],printerror=False)
            ## ~ [2c] Check and enforce symmetry of UP information ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gablam.forceSymmetry(method='min',missing=1.0)
            for seq in rje.sortKeys(updict):
                for hom in updict[seq][0:]:
                    if hom not in updict:
                        self.deBug(hom)
                        self.warnLog('%s not found in UP dictionary.' % hom.shortName())
                        updict[hom] = [hom]
                    while updict[seq].count(hom) > 1: updict[seq].remove(hom) 
                    if seq not in updict[hom]: updict[hom].append(seq)

            ##!# New DomTable supplement of GABLAM #!##
            try: self.getStr('DomTable')
            except: self.setStr({'DomTable': ''})
            if rje.checkForFile(self.getStr('DomTable')):
                domdata = {}
                seqdom = rje.dataDict(self,self.getStr('DomTable'),['Name'],['Type'],lists=True)
                self.printLog('#DOM','Using domain table for additional UPC grouping.')
                for seq in self.seqs():
                    if seq.shortName() in seqdom:
                        for dom in seqdom[seq.shortName()]['Type']:
                            if dom not in domdata: domdata[dom] = []
                            domdata[dom].append(seq)
                #self.deBug(domdata)
                for dom in domdata:
                    for seq1 in domdata[dom]:
                        for seq2 in domdata[dom]:
                            if seq1 != seq2 and seq1 not in updict[seq2]: updict[seq2].append(seq1)
            elif self.getStrLC('DomTable'): self.errorLog('Domain table "%s" not found' % self.getStr('DomTable'),printerror=False)

            ##!# New HomCut reduction of sequences based on BLAST #!##
            if self.getInt('HomCut'):
                (px,prex) = (0.0,self.seqNum())
                ## Remove hub proteins ##
                for seqname in rje.sortKeys(seqdict):
                    self.progLog('\r#HOM','Checking HomCut (%d) %.1f%%' % (self.getInt('HomCut'),px/prex))
                    px += 100.0
                    seq = seqdict[seqname]
                    if len(updict[seq]) > self.getInt('HomCut'):
                        updict.pop(seq)
                        seqdict.pop(seqname)
                        self.obj['SeqList'].removeSeq(text='Exceeds HomCut (%d)' % self.getInt('HomCut'),seq=seq)
                self.printLog('\r#HOM','Checked HomCut (%d). %s of %s seqs remain.' % (self.getInt('HomCut'),rje.iStr(self.seqNum()),rje.iStr(prex)))
                gentxt += '; HomCut=%d' % self.getInt('HomCut')
                ## Clean up ##
                if self.seqNum() < prex:
                    (px,ptot) = (0.0,self.seqNum())
                    for seq in updict:
                        self.progLog('\r#HOM','Cleaning up after HomCut %.1f%%' % (px/ptot))
                        px += 100.0
                        for hit in updict[seq][0:]:
                            if hit not in self.seqs(): updict[seq].remove(hit)
                    self.verbose(0,5,'\r',0)
                    self.obj['SeqList'].saveFasta()
                        
            ### Make UP Clusters ##
            self.printLog('#UP','UP clusters',log=False,newline=False)
            seqdict = self.obj['SeqList'].seqNameDic()
            sortedseq = rje.sortKeys(seqdict)
            while updict:
                sname = sortedseq.pop(0); seq = seqdict[sname]
                while seq not in updict: sname = sortedseq.pop(0); seq = seqdict[sname]
                self.list['UP'].append((seq,))
                uplen = 0
                while uplen < len(self.list['UP'][-1]):
                    self.progLog('\r#UP','%d UP clusters (%d seq remaining)     ' % (self.UPNum(),len(updict)))
                    uplen = len(self.list['UP'][-1])
                    for seq in self.list['UP'][-1]:
                        if updict.has_key(seq):
                            for p in updict.pop(seq):
                                if p not in self.list['UP'][-1]:
                                    self.list['UP'][-1] += (p,)
            self.printLog('\r#UP','%d UP clusters generated from %s sequences.' % (self.UPNum(),rje.iStr(self.seqNum())),log=False)

            ### Make MST ###
            if gablam and gablamdis: self.meanDisMST(gablam)
            else: self.makeMST(gablam)

            ### Output *.upc ###
            UFILE = open('%s.upc' % upcbase,'w')
            UFILE.write('#%s# %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))
            rje.writeDelimit(UFILE,['UP','N','MST','Seqs'],'\t')
            for u in range(self.UPNum()):
                upc = self.list['UP'][u]
                seqs = upc[0].shortName()
                for seq in upc[1:]: seqs += ' %s' % seq.shortName()
                rje.writeDelimit(UFILE,[u+1,len(upc),'%.3f' % (len(upc)*self.dict['MST'][upc]),seqs],'\t')
            UFILE.close()
            self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))

            ### Output DisMatrix ###
            gablam.setStr({'Name':'%s GABLAM' % self.dataset()})
            if 0 < self.getInt('MaxUPC') < self.UPNum():
                self.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). No *.dis.tdt.' % (self.UPNum(),self.getInt('MaxUPC')))
            else:
                if savedis or self.extras(3):
                    gablam.saveMatrix(self.seqs(),filename='%s.dis.tdt' % upcbase,delimit='\t',format='text',log=True)
                    self.dict['Output']['dismatrix'] = '%s.dis.tdt' % upcbase
                if self.extras(3):
                    gablam.saveMatrix(self.seqs(),filename='%s.phydis.txt' % upcbase,format='phylip',log=True)
            return True                
                
        except: self.errorLog('Fatal Error in %s.makeUPC(). Check for old results files with same name but different sequences.' % self.prog())
        return False
#########################################################################################################################
    def readUPC(self):  ### Looks for UPC file and loads details from file.
        '''Looks for UPC file and loads details from file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('FUPC'): self.printLog('#DEV','Need to add reading of FUPC sequence file.'); return False
            #!#NB. seqbase and upcbase are probably the same. Should check and unify at some point
            upcbase = '%s%s' % (self.getStr('ResDir'),rje.baseFile(self.getStr('Input'),strip_path=True))
            seqbase = self.getStr('ResDir') + rje.baseFile(self.obj['SeqList'].name(),strip_path=True)
            ## ~ [0a] ~ List of possible UPC file basenames to check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            bfiles = []
            for bfile in [upcbase,seqbase,'%s%s' % (self.getStr('BuildPath'),os.path.split(seqbase)[1]),self.getStr('Basefile'),'%s%s' % (self.getStr('BuildPath'),self.dataset())]:
                if bfile not in bfiles: bfiles.append(bfile)
            ### ~ [1] ~ Check UPC files and load if possible ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for pfile in bfiles:
                if self.getBool('FUPC'): ufile = '%s.fupc' % pfile
                else: ufile = '%s.upc' % pfile
                if os.path.exists(ufile):
                    if self.loadUPCFromFile(ufile):
                        seqfile = '%s.slimdb' % seqbase
                        if not os.path.exists(seqfile): seqfile = '%s.slimdb' % upcbase
                        self.dict['Output']['upc'] = ufile
                        self.dict['Output']['slimdb'] = seqfile
                        self.obj['SeqList'].setStr({'Name':seqfile})
                        self.setStr({'SLiMDB':seqfile})
                        if not rje.exists(seqfile): self.obj['SeqList'].saveFasta()
                        return True
                    else: self.warnLog('UPC file "%s" found but load failed.' % ufile)
                elif self.dev(): self.printLog('#UPC','Cannot read: UPC file "%s" not found.' % ufile)
            self.printLog('#UPC','No UPC file found.')
            return False
        except SystemExit: raise
        except ValueError: self.errorLog('Problem with *.upc formatting in %s.readUPC().' % self.prog())
        except: self.errorLog('Error in %s.readUPC().' % self.prog())
        return False
#########################################################################################################################
    def loadUPCFromFile(self,ufile):    ### Loads UPC details from specific file
        '''Loads UPC details from specific file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(ufile):
                self.printLog('#UPC','Loading UPC from "%s".' % ufile)
                ulines = self.loadFromFile(ufile,chomplines=True)
            else:
                self.printLog('#UPC','Cannot load: UPC file "%s" not found.' % ufile)
                return False
            seqdict = self.obj['SeqList'].seqNameDic()
            ### ~ [1] ~ Read UPC File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['UP'] = []
            self.dict['MST'] = {}
            ## ~ [1a] ~ Dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.printLog('\r#UPC',ulines[0])
            data = rje.matchExp('#(\S.*)# (\d+) Seq; (\d+) UPC; (\S+) MST',ulines[0])
            base = data[0]
            seqx = int(data[1])
            if seqx != self.seqNum(): self.printLog('#ERR','Wrong number of sequences in UPC file!'); return False
            upx = int(data[2])
            self.setNum({'MST':float(data[3])})
            ## ~ [1b] ~ UPCs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for u in range(upx):
                data = string.split(ulines[u+2],'\t')
                if int(data[0]) != (u+1): self.printLog('#ERR','UP identifier number wrong in file!'); raise ValueError
                useqs = []
                for seq in string.split(data[-1]): useqs.append(seqdict[seq])
                if len(useqs) != int(data[1]): self.printLog('#ERR','UP SeqNum mismatch!'); raise ValueError
                upc = (useqs[0],)
                for seq in useqs[1:]: upc += (seq,)
                if len(upc) != int(data[1]): self.printLog('#ERR','No. seqs in UP does not match file!'); raise ValueError
                seqx -= len(upc)
                self.dict['MST'][upc] = float(data[2]) / len(upc)
                self.list['UP'].append(upc)
            if upx != self.UPNum(): self.printLog('#ERR','Wrong number of UPC in UPC file!'); return False
            if seqx != 0: self.printLog('#ERR','Wrong number of sequences read from UPC file!'); return False
            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST')))
            return True
        except KeyError:
            if self.i() > 0:
                if not rje.yesNo('Possible change in sequence names. Force file regeneration and continue (Y) or quit (N)?'):
                    sys.exit(1)
            self.errorLog('Possible change in sequence names. Forcing file regeneration [force=T]',printerror=False)
            self.setBool({'Force':True})
        except ValueError: self.errorLog('Problem with *.upc formatting or content in %s.loadUPCFromFile().' % self.prog())
        except: self.errorLog('Error in %s.loadUPCFromFile().' % self.prog())
        return False
#########################################################################################################################
    def makeMST(self,gablam=None):   ### Makes UPC dictionary from GABLAM
        '''Makes UPC dictionary from GABLAM.'''
        self.dict['MST'] = {}
        if gablam:
            gablam.checkSymmetry(force=True)
            self.setNum({'MST':gablam.MST()})
            ux = 0.0; utot = len(self.list['UP'])
            for upc in self.list['UP']:
                self.progLog('\r#MST','Calculating MST: %.1f%%' % (ux/utot)); ux += 100.0
                uplist = []
                for seq in upc: uplist.append(seq)
                self.dict['MST'][upc] = gablam.MST(uplist) / len(upc)   # Mean Total * MST Value
            self.printLog('\r#MST','MST for %s sequences = %.3f.' % (rje.iStr(self.seqNum()),self.getNum('MST')))
        else:
            self.setNum({'MST':self.seqNum()})
            for upc in self.list['UP']: self.dict['MST'][upc] = 1.0
#########################################################################################################################
    def meanDisMST(self,gablam=None):   ### Makes UPC dictionary from GABLAM using mean distance instead of MST
        '''Makes UPC dictionary from GABLAM using mean distance instead of MST.'''
        self.dict['MST'] = {}
        self.setNum({'MST': 0.0})
        self.printLog('#SYM','Checking symmetry of GABLAM DisMatrix')
        gablam.checkSymmetry(force=True)
        (ux,upx) = (0.0,100.0/len(self.list['UP']))
        for upc in self.list['UP']:
            self.progLog('\r#MST','Calculating MeanDis MST replacement %.1f%%' % ux)
            ux += upx
            uplist = []
            if len(upc) == 1:
                self.dict['MST'][upc] = 1.0
                self.setNum({'MST':self.getNum('MST') + 1.0})
                continue
            self.dict['MST'][upc] = float((len(upc) * (len(upc)-1)))
            for s1 in upc:
                if s1 not in gablam.dict['Matrix']: continue
                for s2 in gablam.dict['Matrix'][s1]:
                    if s1 == s2 or s2 not in upc: continue
                    self.dict['MST'][upc] -= gablam.dict['Matrix'][s1][s2]
            self.dict['MST'][upc] = self.dict['MST'][upc] / (len(upc) * (len(upc)-1))   # Mean Total * MST Value
            self.setNum({'MST':self.getNum('MST') + self.dict['MST'][upc] * len(upc)})
        self.printLog('\r#MST','Calculation of MeanDis MST replacement complete.')
#########################################################################################################################
    def makeFUPC(self):  ### Generates Fragmented UP Clusters from self.obj['SeqList'] using BLAST
        '''Generates Fragmented UP Clusters from self.obj['SeqList'] using BLAST.'''
        try:### ~ [0] General set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['UP'] = []
            self.dict['MST'] = {}
            self.dict['FragSeq'] = fragseq = {}    # Dictionary of {seq:[blast sections]}
            upcbase = '%s%s' % (self.getStr('ResDir'),rje.baseFile(self.getStr('Input'),strip_path=True))
            seqlist = self.obj['SeqList']
            seqfile = '%s.slimdb' % upcbase
            seqlist.setStr({'Name':seqfile})
            if self.force() or not rje.exists(seqfile): seqlist.saveFasta()
            seqdict = seqlist.seqNameDic()
            fupcfas = '%s.fupc.fas' % upcbase
            fupcdis = '%s.fupc.dis.tdt' % upcbase
            fupcloc = '%s.fupc.Local.tdt' % upcbase     #?# Add looking for this at appropriate time    #?#
            bcmd = ['blastcf=F','blastf=T','blaste=1e-4'] + self.cmd_list
            bcmd += ['blasto=%s.self.blast' % upcbase,'restab=Search,Hit,Local']
            blast = rje_blast.blastObj(self.log,bcmd,type='New')

            ### ~ [1] All-by-all BLAST and read local hits table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # RJE_BLASTV2 should read in local results if present.
            # - only delete Local table if (not test/dev and) right output made. (Unless self.extras(2))
            #?# Add a method to delete results tables to BLAST object?
            blast.setStr({'InFile':seqfile,'DBase':seqfile})
            blast.setInt({'OneLine':self.seqNum(),'HitAln':self.seqNum()})
            blast.run(format=True,gablam=False)

            ### ~ [2] From local hit table, split into fragseq sections ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ldb = blast.db('Local') # ['Query','Hit','AlnID']
            # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives','QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq'],
            ldb.dropFields(['BitScore','Expect'])
            ldb.list['Fields'].insert(ldb.list['Fields'].index('Query')+1,'QParent')
            ldb.list['Fields'].insert(ldb.list['Fields'].index('QryEnd')+1,'QryPar')
            ldb.list['Fields'].insert(ldb.list['Fields'].index('Hit')+1,'HParent')
            ldb.list['Fields'].insert(ldb.list['Fields'].index('SbjEnd')+1,'SbjPar')
            self.printLog('\r#FIELD','Added fields "QryPar" and "SbjPar" to table "Local".')
            selfhitx = 0
            for entry in ldb.entries()[0:]:
                if entry['Query'] == entry['Hit']: ldb.dropEntry(entry); selfhitx += 1
                else:
                    entry['QryPar'] = entry['QryStart']; entry['SbjPar'] = entry['SbjStart'];
                    entry['QParent'] = entry['Query']; entry['HParent'] = entry['Hit']
            self.printLog('#SELF','%s self hits removed from Local hit table.' % rje.iStr(selfhitx))
            fragmenting = True  # Whether BLAST local hits are still being fragmented
            skiplist = []       # List of proteins that are completely fragmented (all local hits for a frag are different Hits)
            db = blast.db()
            fragdb = db.addEmptyTable('SeqFrag',['Name','Start','End'],['Name'],log=True)
            ldb.index('Query'); ldb.index('Hit')
            fragseq = {}

            while fragmenting:
                fragmenting = False
                seqx = seqlist.seqNum(); sx = 0.0; fx = 0
                if self.dev(): self.setBool({'DeBug':True})
                if skiplist: self.printLog('#SKIP','%s Parents to skip; %s Sequences with FragSeq data' % (rje.iLen(skiplist),rje.iLen(fragseq))); self.deBug('...')
                minfrag = 5
                fragdb.clear()
                for sname in seqdict:
                    self.progLog('\r#FRAG','Fragmenting proteins on Local BLAST searches: %.2f%%' % (sx/seqx)); sx += 100.0
                    seq = seqdict[sname]
                    #!# fragseq[seq] created here. If keeping unchanged fragments, will need to check (and skip)
                    # >> existing entries here. Can also use to make a list of sequences to analyse later.
                    #if (seq in fragseq) != (sname in skiplist): raise ValueError('Skiplist/Fragseq mismatch!')
                    if seq in fragseq: continue     # Skipped from earlier analysis
                    ## Build list of BLAST fragment start and end sites
                    frag = ['%s(' % rje.preZero(1,seq.seqLen()),'%d)' % seq.seqLen()]
                    for (nfield,sfield,efield) in [('Query','QryStart','QryEnd'),('Hit','SbjStart','SbjEnd')]:
                        for entry in ldb.indexEntries(nfield,sname):
                            frag.append('%s(' % rje.preZero(entry[sfield],seq.seqLen()))
                            frag.append('%s)' % rje.preZero(entry[efield],seq.seqLen()))
                    ## Collapse lists using minfrag
                    #self.setBool({'DeBug':sname in ['APC_HUMAN__P25054']})
                    frag.sort()     # ['001(', '006(', '006)', '007(', '010)', '180)']
                    #self.bugPrint(frag)
                    collapse = True
                    ## i. Replace 'x(', 'x)' with 'x-1)','x+1(' and collapse identities
                    while collapse:
                        collapse = False; i = 0
                        while i < (len(frag) - 1):  # Need to compare current position with next
                            if frag[i] == frag[i+1]: frag.pop(i+1); collapse = True
                            elif frag[i][:-1] == frag[i+1][:-1]:
                                if frag[i][-1] != '(' or frag[i+1][-1] != ')': raise ValueError
                                ipos = int(frag[i][:-1]); jpos = int(frag[i+1][:-1])
                                frag[i] = '%s)' % rje.preZero(ipos-1,seq.seqLen())
                                frag[i+1] = '%s(' % rje.preZero(jpos+1,seq.seqLen())
                                i += 1; collapse = True
                            else: i += 1
                        frag.sort()
                    ## ii. Collapse opens and closures using minfrag
                    collapse = True
                    while collapse:
                        collapse = False; i = 0
                        while i < (len(frag) - 1):  # Need to compare current position with next
                            if '(' in [frag[i][-1], frag[i+1][-1]]: i += 1;  continue
                            ipos = int(frag[i][:-1]); jpos = int(frag[i+1][:-1])
                            if jpos - ipos < minfrag: frag.pop(i); collapse = True
                            else: i += 1
                        i = len(frag) - 1
                        while i > 0:
                            if ')' in [frag[i][-1], frag[i-1][-1]]: i -= 1;  continue
                            jpos = int(frag[i][:-1]); ipos = int(frag[i-1][:-1])
                            if jpos - ipos < minfrag: frag.pop(i); collapse = True
                            i -= 1
                    #self.debug(frag)
                    ## iii. Convert into tuples
                    fragseq[seq] = []
                    fragtmp = frag[0:]
                    while frag:
                        if frag[0][-1] == ')': raise ValueError('%s' % frag)
                        if fragseq[seq] and int(frag[0][:-1]) != fragseq[seq][-1][1] + 1:
                            fragseq[seq].append((fragseq[seq][-1][1] + 1,int(frag[0][:-1])-1))
                        if frag[1][-1] == '(':
                            fragseq[seq].append((int(frag.pop(0)[:-1]),int(frag[0][:-1])-1))
                        else:
                            fragseq[seq].append((int(frag.pop(0)[:-1]),int(frag.pop(0)[:-1])))
                            if frag and frag[0][-1] == ')': frag.insert(0,'%s(' % rje.preZero(fragseq[seq][-1][1]+1,seq.seqLen()))
                    ### iv. Check fragmentation is complete
                    for i in range(len(fragseq[seq])-1):
                        try:
                            if fragseq[seq][i][1] != fragseq[seq][i+1][0] - 1:
                                self.debug(fragtmp)
                                raise ValueError('%s: %d != %d; %s' % (sname,fragseq[seq][i][1],fragseq[seq][i+1][0] + 1,fragseq[seq]))
                        except:
                            self.errorLog('%s (%d): %s' % (sname,i,fragseq[seq]))
                    fx += len(fragseq[seq])
                self.printLog('\r#FRAG','Fragmented %s proteins on Local BLAST searches into %s fragments.' % (rje.iStr(seqx),rje.iStr(fx)))

                ### ~ [3] Split proteins and generate fragseq GABLAM matrix *.fupc.dis.tdt &  and save *.fupc.fas ~~~~~~~ ###
                if self.dev(): self.setBool({'DeBug':True})
                ## ~ [3a] Save new sequences and update Local BLAST table with frag names ~~~~~~~~~~~~~ ##
                #x# Moved further up for skiplist efficiency: fragdb = db.addEmptyTable('SeqFrag',['Name','Start','End'],['Name'],log=True)
                sx = 0.0; fx = 0
                for sname in seqdict:
                    self.progLog('\r#FRAG','Regenerating BLAST hits from fragmented proteins: %.2f%%' % (sx/seqx)); sx += 100.0
                    seq = seqdict[sname]
                    #!# fragdb entries are added here. If keeping unchanged elements, will need to assess for existence here
                    # >> remake those that are missing.
                    if sname in skiplist: continue
                    #!# Can we add a fragdb check?! raise ValueError('Skiplist/Fragseq mismatch!')
                    if len(fragseq[seq]) > 1:
                        #self.bugPrint(fragseq[seq])
                        for i in range(len(fragseq[seq])):
                            (x,y) = fragseq[seq][i]
                            fragdb.addEntry({'Name':'%s.%d' % (sname,i+1),'Start':x,'End':y})
                    else:
                        fragdb.addEntry({'Name':sname,'Start':1,'End':seq.seqLen()})
                    for (nfield,sfield,efield,alt) in [('Query','QryStart','QryEnd','Sbj'),('Hit','SbjStart','SbjEnd','Qry')]:
                        if sname not in ldb.index(nfield): continue
                        for entry in ldb.indexEntries(nfield,sname)[0:]:
                            dropentry = False   # Whether to drop the entry due to splitting

                            for i in range(len(fragseq[seq])):
                                (nt,ct) = fragseq[seq][i]
                                if nt > entry[efield] or ct < entry[sfield]: continue       # Outside of hit
                                elif nt <= entry[sfield] and ct >= entry[efield]: continue  # Local hit contained within fragment
                                else:     # Need to fragment local hit at fragment end
                                    dropentry = True; fragmenting = True
                                    fentry = rje.combineDict({nfield:sname},entry,False)
                                    # ['Query','Hit','AlnID','BitScore','Expect','Length','Identity','Positives',
                                    # 'QryStart','QryEnd','SbjStart','SbjEnd','QrySeq','SbjSeq','AlnSeq']
                                    fragaln = entry[sfield[:3] + 'Seq']
                                    alnct = alnnt = 0; alnpos = entry[sfield]
                                    while fragaln[alnnt] == '-': alnnt += 1
                                    while alnpos < nt:
                                        alnnt += 1
                                        try:
                                            while fragaln[alnnt] == '-': alnnt += 1
                                        except:
                                            self.bugPrint(fentry)
                                            self.bugPrint('nt: %s; ct: %s; alnnt: %s; alnct: %s; alnpos: %s' % (nt,ct,alnnt,alnct,alnpos))
                                            self.debug(fragaln)
                                            self.errorLog('Fuck')
                                            break
                                        alnpos += 1
                                    if alnpos < nt:
                                        self.warnLog('Fragment alignment issue: %s' % fentry)
                                        continue    # Problem with alignment
                                    # alnpos is the starting position for the fragment; alnnt is the corresponding position in the alignments
                                    alnct = alnnt
                                    while alnct < len(fragaln) and fragaln[alnct] == '-': alnct += 1
                                    while alnpos < (min(ct,entry[efield])):
                                        alnct += 1
                                        while alnct < len(fragaln) and fragaln[alnct] == '-': alnct += 1
                                        if alnct < len(fragaln): alnpos += 1
                                        else: break
                                    # Should now have alnnt and alnct marking the alignment start and end positions
                                    # Truncate alignments & recalculate ends
                                    truncx = 0
                                    while '-' in [fentry['QrySeq'][alnnt], fentry['SbjSeq'][alnnt]]: alnnt += 1; truncx += 1
                                    while '-' in [fentry['QrySeq'][alnct], fentry['SbjSeq'][alnct]]: alnct -= 1; truncx += 1
                                    for afield in ['QrySeq','SbjSeq','AlnSeq']: fentry[afield] = fentry[afield][alnnt:alnct+1]
                                    for prot in ['Qry','Sbj']:
                                        fentry[prot+'Par'] += (alnnt - string.count(entry[prot+'Seq'][:alnnt],'-'))
                                        fentry[prot+'Start'] += (alnnt - string.count(entry[prot+'Seq'][:alnnt],'-'))
                                        fentry[prot+'End'] = fentry[prot+'Start'] + len(fentry[prot+'Seq']) - string.count(fentry[prot+'Seq'],'-') - 1
                                    fentry['AlnID'] = '%s.%s' % (fentry['AlnID'],fentry[sfield])
                                    # Recalculate Identity and Positives
                                    if nfield == 'Hit': fentry['Length'] = fragseq[seq][i][1] - fragseq[seq][i][0] + 1
                                    fentry['Positives'] = len(fentry['AlnSeq']) - string.count(fentry['AlnSeq'],' ')
                                    fentry['Identity'] = fentry['Positives'] - string.count(fentry['AlnSeq'],'+')
                                    if min(fentry['Identity'],fentry['Positives']) < 0: raise ValueError('Identity/Positives calculation error!\nAlnSeq: %s\n' % fentry['AlnSeq'])
                                    if fentry['Identity']: ldb.addEntry(fentry)
                            if dropentry:
                                ldb.dropEntry(entry)
                self.printLog('\r#FRAG','Regenerated %s local BLAST hits from fragmented proteins.' % rje.iStr(ldb.entryNum()))
                if self.test() or self.dev():
                    ldb.saveToFile('%s.fupc.Local.tdt' % upcbase)
                    fragdb.saveToFile() # For debugging.

            ### ~ Rename and renumber the local BLAST hits ~ ###
            sx = 0.0
            for sname in seqdict:
                self.progLog('\r#FRAG','Renumbering/renaming BLAST hits from fragmented proteins: %.2f%%' % (sx/seqx)); sx += 100.0
                seq = seqdict[sname]
                for (nfield,sfield,efield,alt) in [('Query','QryStart','QryEnd','Sbj'),('Hit','SbjStart','SbjEnd','Qry')]:
                    if sname not in ldb.index(nfield): continue
                    for entry in ldb.indexEntries(nfield,sname)[0:]:
                        for i in range(len(fragseq[seq])):
                            (nt,ct) = fragseq[seq][i]
                            if nt > entry[efield] or ct < entry[sfield]: continue       # Outside of hit
                            elif nt <= entry[sfield] and ct >= entry[efield]:  # Local hit contained within fragment
                                if len(fragseq[seq]) > 1: entry[nfield] = '%s.%d' % (sname,i+1)
                                #?# Not convinced that start/end positions need changing #?#
                                #entry[sfield] -= (fragseq[seq][i][0]-1)
                                #entry[efield] -= (fragseq[seq][i][0]-1)
                            else: raise ValueError('Still have %s local hits overlapping fragment ends!' % sname)
            self.printLog('\r#FRAG','Renumbered/renamed BLAST hits from fragmented proteins: %.2f%%' % (sx/seqx)); sx += 100.0
            ldb.remakeKeys(); ldb.dict['Index'] = {}
            ldb.saveToFile('%s.fupc.Local.tdt' % upcbase)

            if self.dev() or self.test():   # Integrity check
                locfrag = rje.sortUnique(ldb.index('Query').keys() + ldb.index('Hit').keys()); mx = 0
                for sname in locfrag:
                    if sname not in fragdb.data(): self.warnLog('Local fragment %s not found in fragdb.' % sname); mx += 1
                if mx: raise ValueError('%s local hits missing protein fragments.' % rje.iStr(mx))
                self.printLog('#FRAG','%s local hits missing protein fragments.' % rje.iStr(mx)); mx = 0
                for sname in fragdb.dataKeys():
                    if sname not in locfrag and sname not in seqdict: mx += 1
                self.printLog('#FRAG','%s fragments w/o local BLAST hits.' % rje.iStr(mx))

            if self.dev(): self.setBool({'DeBug':True})
            #self.deBug('')

            #fragdb.saveToFile()
            #?# rje_seqlist.SeqList(self.log,self.cmd_list+['seqin=%s' % fupcfas,'seqmode=list'])
            ## ~ [3b] Regenerate Search & Hit Table from Local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            db.deleteTable('Search')
            sdb = db.addEmptyTable('Search',['Query','Length','Hits'],['Query'],log=True)
            db.deleteTable('Hit')
            hdb = db.addEmptyTable('Hit',['Query','Hit','Length','Aln','GablamFrag','LocalCut','GABLAM','HitAln'],['Query','Hit','HitAln'],log=True)
            for entry in ldb.entries():
                hdb.addEntry({'Query':entry['Query'],'Hit':entry['Hit'],'Aln':1,'HitAln':entry['AlnID']})
            hdb.compress(['Query','Hit'],rules={'Aln':'sum','HitAln':'str'})
            hdb.renameField('HitAln','Length')
            for sname in fragdb.dataKeys():
                for entry in hdb.indexEntries('Hit',sname): entry['Length'] = fragdb.data(sname)['End'] - fragdb.data(sname)['Start'] + 1
                sdb.addEntry({'Query':sname,'Hits':len(hdb.indexEntries('Hit',sname)),
                              'Length':fragdb.data(sname)['End'] - fragdb.data(sname)['Start'] + 1})
            if self.test() or self.dev():
                sdb.saveToFile('%s.fupc.Search.tdt' % upcbase)
                hdb.saveToFile('%s.fupc.Hit.tdt' % upcbase)
                ldb.saveToFile('%s.fupc.Local.tdt' % upcbase)

            if self.dev(): self.setBool({'DeBug':True})
            #self.deBug('')

            ## ~ [3c] Generate GABLAM matrix for fragmented proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# hitaln[hit['Hit']].append(aln); hit['Aln'] += 1
            #i# if local: self.db('Local').addEntry(aln)
            #for search in sdb.entries():
            #    for hit in hdb.indexEntries('Query',search['Query']):
            #        hitaln = ldb.indexEntries('#Query#|#Hit#','%s|%s' % (search['Query'],hit['Hit']))
            #        gablam = self.globalFromLocal(search,hit,hitaln)     # This will also clear the alignment data.
            #        for qh in ['Query','Hit']:
            #            blast.db('GABLAM').addEntry(rje.combineDict({'Query':search['Query'],'Hit':hit['Hit'],'QryHit':qh},gablam[qh]))

            ## ~ [3d] Convert BLAST hits into hit matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            minfrag = 20
            self.progLog('#DIS','Generating FUPC distance matrix...')
            idmatrix = rje_dismatrix.DisMatrix(self.log,self.cmd_list)      # blast.resultsIDMatrix()
            for entry in hdb.entries():
                # Don't add any fragments with len < minfrag
                if fragdb.data(entry['Query'])['End'] - fragdb.data(entry['Query'])['Start'] + 1 < minfrag: continue
                if fragdb.data(entry['Hit'])['End'] - fragdb.data(entry['Hit'])['Start'] + 1 < minfrag: continue
                idmatrix.addDis(entry['Query'],entry['Hit'],0.0)
            self.printLog('\r#DIS','Generated FUPC distance matrix for clustering.')

            ## ~ [3e] Generate UPC from BLAST Clusters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# NB. Currently clustering sequence names not sequence objects: will need to be replaced later.
            updict = {}; seqx = fragdb.entryNum(); sx = 0.0
            for bclust in idmatrix.cluster(maxdis=1.0,singletons=True):
                for seq in bclust:
                    self.progLog('\r#CLUST','Converting clusters for FUPC: %.2f%%' % (sx/seqx)); sx += 100.0
                    if seq not in updict: updict[seq] = []
                    for seq2 in bclust:
                        if seq2 not in updict[seq]: updict[seq].append(seq2)
            for seq in fragdb.dataKeys():
                if seq not in updict: updict[seq] = [seq]
            self.printLog('\r#UPDIC','Basic UP Dictionary made for %s fragments.' % fragdb.entryNum(),log=False)
            ## ~ Check and enforce symmetry of UP information ~ ##
            seqx = len(updict); sx = 0.0
            for seq in rje.sortKeys(updict):
                self.progLog('\r#UPC','Tidying FUP dictionary: %.2f%%' % (sx/seqx)); sx += 100.0
                for hom in updict[seq][0:]:
                    if hom not in updict:
                        self.deBug(hom)
                        self.warnLog('%s not found in FUP dictionary.' % hom.shortName())
                        updict[hom] = [hom]
                    while updict[seq].count(hom) > 1: updict[seq].remove(hom)
                    if seq not in updict[hom]: updict[hom].append(seq)
            self.printLog('\r#UPDIC','FUP symmetry check complete.',log=False,clear=35)
            ## ~ Make UP Clusters ~ ##
            fragupc = {}    # {i: List of seqname}
            sequpc = {}     # {parent:[fragupc ids]}    - use to combine fragupc
            upcseq = {}     # {i: List of parents}
            self.progLog('\r#FUP','FUP clusters')
            sortedseq = rje.sortKeys(updict)
            while updict:
                sname = sortedseq.pop(0)
                if sname not in updict: continue    # Already incorporated
                i = len(fragupc)
                fragupc[i] = [sname]
                uplen = 0
                while uplen < len(fragupc[i]):
                    self.progLog('\r#FUP','%s FUP clusters (%s seq remaining)     ' % (rje.iLen(fragupc),rje.iLen(updict)))
                    uplen = len(fragupc[i])
                    for seq in fragupc[i]:
                        if updict.has_key(seq):
                            for p in updict.pop(seq):
                                if p not in fragupc[i]: fragupc[i].append(p)
                for sname in fragupc[i]:
                    if sname not in seqdict:
                        parent = rje.matchExp('^(\S+)\.(\d+)$',sname)[0]
                        if parent not in sequpc: sequpc[parent] = []
                        if i not in sequpc[parent]: sequpc[parent].append(i)
            self.printLog('\r#FUP','%s FUP clusters generated from %s sequence fragments.' % (rje.iLen(fragupc),rje.iStr(fragdb.entryNum())))
            ## ~ Combine UPC with same sequence composition ~ ##
            self.progLog('#FUP','Combining FUPC with same sequence composition...')
            for parent in rje.sortKeys(sequpc):
                if len(sequpc[parent]) == 1: sequpc.pop(parent); continue
                for i in sequpc[parent]:
                    if i not in upcseq: upcseq[i] = [parent]
                    else: upcseq[i].append(parent)  # This should now have all parents in 2+ UPC in alphabetical order
            for i in rje.sortKeys(upcseq):
                if i not in upcseq: continue    # Absorbed
                for j in rje.sortKeys(upcseq,revsort=True):
                    if j <= i: break
                    if upcseq[i] == upcseq[j]:  # Combine
                        fragupc[i] += fragupc.pop(j)
                        upcseq.pop(j)
            self.printLog('\r#FUP','%s tidied FUPC from %s sequence fragments.' % (rje.iLen(fragupc),rje.iStr(fragdb.entryNum())))

            ### ~ [4] Review UPCs and combine X.i and X.(i+1) that fall into the same UPC ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['FullSeq'] = rje_seq.SeqList(self.log,self.cmd_list+['autoload=F'])
            self.obj['FullSeq'].seq = seqlist.seq[0:]
            seqlist.seq = []
            #FFAS = open(fupcfas,'w')    #?# Do we need this yet?    #?# Can save from seqlist later
            newseqs = []    # list of (seqname,sequence) tuples
            for i in fragupc:
                upc = fragupc[i]
                seqfrag = {}    # Parent: fraglist
                upmade = False
                for sname in upc:
                    if sname in seqdict:    # Full-length protein
                        seqlist.seq.append(seqdict[sname])
                        if upmade: self.list['UP'][-1] += (seq,)
                        else: self.list['UP'].append((seq,)); upmade = True
                    else:
                        (parent,frag) = rje.matchExp('^(\S+)\.(\d+)$',sname)
                        try: seqfrag[parent].append(int(frag))
                        except: seqfrag[parent] = [int(frag)]
                ## ~ Merge fragments and expand seqdict with new sequences ~ ##
                for parent in seqfrag:
                    seq = seqdict[parent]
                    seqfrag[parent].sort()
                    while seqfrag[parent]:
                        i = j = seqfrag[parent].pop(0)
                        while seqfrag[parent] and seqfrag[parent][0] == (j+1): j = seqfrag[parent].pop(0)
                        x = fragdb.data('%s.%d' % (parent,i))['Start']
                        y = fragdb.data('%s.%d' % (parent,j))['End']
                        if i == j: sname = '%s.%d %d-%d' % (parent,i,x,y)
                        else:
                            sname = '%s.%d-%d %d-%d' % (parent,i,j,x,y)
                            fragdb.addEntry({'Name':'%s.%d-%d' % (parent,i,j),'Start':x,'End':y})
                            for oldi in range(i,j+1):
                                oldfrag = '%s.%d' % (parent,oldi)
                                fragdb.data().pop(oldfrag)
                                for entry in ldb.indexEntries('Query',oldfrag): entry['Query'] = string.split(sname)[0]
                                for entry in ldb.indexEntries('Hit',oldfrag): entry['Hit'] = string.split(sname)[0]
                            #!# Update blast tables #!#
                            # sdb ['Query','Length','Hits'],['Query']
                            # hdb ['Query','Hit','Length','Aln','GablamFrag','LocalCut','GABLAM'],['Query','Hit']
                        sequence = '%s%s' % ('x' * (x-1),seq.info['Sequence'][x-1:y])
                        #FFAS.write('>%s\n%s\n' % (sname,sequence))
                        #self.debug('>%s\n%s\n' % (sname,sequence))
                        fseq = seqlist._addSeq(sname,sequence)
                        if upmade: self.list['UP'][-1] += (fseq,)
                        else: self.list['UP'].append((fseq,)); upmade = True
                        seqdict[string.split(sname)[0]] = fseq
            fragdb.saveToFile()
            #FFAS.close()

            ### ~ [5] Save revised sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist.saveFasta(seqfile=fupcfas)
            seqdict = seqlist.seqNameDic()

            ### ~ [6] Calculate MST for each UPC and save UPC File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [6a] Regenerate Search & Hit Table from Local hits ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            ldb.remakeKeys()
            db.deleteTable('Search')
            sdb = db.addEmptyTable('Search',['Query','Length','Hits'],['Query'],log=True)
            db.deleteTable('Hit')
            hdb = db.addEmptyTable('Hit',['Query','Hit','Length','Aln','GablamFrag','LocalCut','GABLAM','HitAln'],['Query','Hit','HitAln'],log=True)
            for entry in ldb.entries():
                hdb.addEntry({'Query':entry['Query'],'Hit':entry['Hit'],'Aln':1,'HitAln':entry['AlnID']})
            hdb.compress(['Query','Hit'],rules={'Aln':'sum','HitAln':'max'})
            hdb.renameField('HitAln','Length')
            #for sname in seqdict:
            for sname in rje.sortKeys(fragdb.data()):
                self.debug(sname)
                #i# The shorter fragment lengths need to be used in order to make the correct % calculations in GABLAM
                for entry in hdb.indexEntries('Hit',sname):
                    entry['Length'] = fragdb.data(sname)['End'] - fragdb.data(sname)['Start'] + 1
                sdb.addEntry({'Query':sname,'Hits':len(hdb.indexEntries('Hit',sname)),
                              'Length':fragdb.data(sname)['End'] - fragdb.data(sname)['Start'] + 1})
            #i# Local hit positions will need to be re-jigged to match shortened lengths
            for entry in ldb.entries():
                for (nfield,sfield,efield,alt) in [('Query','QryStart','QryEnd','Sbj'),('Hit','SbjStart','SbjEnd','Qry')]:
                    fentry = fragdb.data(entry[nfield])
                    entry[sfield] = entry[sfield] - fentry['Start'] + 1
                    entry[efield] = entry[efield] - fentry['Start'] + 1

            if self.test() or self.dev():
                sdb.saveToFile('%s.fupc.Search.tdt' % upcbase)
                hdb.saveToFile('%s.fupc.Hit.tdt' % upcbase)
                ldb.saveToFile('%s.fupc.Local.tdt' % upcbase)
            ## ~ [6b] Generate GABLAM matrix for fragmented proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #i# hitaln[hit['Hit']].append(aln); hit['Aln'] += 1
            #i# if local: self.db('Local').addEntry(aln)
            if self.dev(): self.setBool({'DeBug':True})
            ldb.makeField('#Query#|#Hit#')
            for search in sdb.entries():
                for hit in hdb.indexEntries('Query',search['Query']):
                    hitaln = ldb.indexEntries('#Query#|#Hit#','%s|%s' % (search['Query'],hit['Hit']))
                    self.debug('%s\n%s\n%s\n' % (search,hit,hitaln))
                    gablam = blast.globalFromLocal(search,hit,hitaln)     # This will also clear the alignment data.
                    for qh in ['Query','Hit']:
                        blast.db('GABLAM').addEntry(rje.combineDict({'Query':search['Query'],'Hit':hit['Hit'],'QryHit':qh},gablam[qh]))
            ## ~ [6c] Calculate MST from GABLAM ID Matrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            idmatrix = blast.resultsIDMatrix()
            idmatrix.forceSymmetry(method='min',missing=1.0)
            self.makeMST(idmatrix)
            ## ~ [6d] Output *.fupc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            gentxt = 'blaste=%s, blastcf=%s, blastf=%s' % (rje_slim.expectString(blast.getNum('E-Value')),str(blast.getBool('Composition Statistics'))[0],str(blast.getBool('Complexity Filter'))[0])
            UFILE = open('%s.fupc' % upcbase,'w')
            UFILE.write('#%s# %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))
            rje.writeDelimit(UFILE,['UP','N','MST','Seqs'],'\t')
            for u in range(self.UPNum()):
                upc = self.list['UP'][u]
                seqs = upc[0].shortName()
                for seq in upc[1:]: seqs += ' %s' % seq.shortName()
                rje.writeDelimit(UFILE,[u+1,len(upc),'%.3f' % (len(upc)*self.dict['MST'][upc]),seqs],'\t')
            UFILE.close()
            self.printLog('\r#UP','%s: %d Seq; %d UPC; %.3f MST; %s\n' % (self.dataset(),self.seqNum(),self.UPNum(),self.getNum('MST'),gentxt))
            ## ~ [6e] Output DisMatrix ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #?# Not sure if this is desirable and/or useful. #?#
            idmatrix.setStr({'Name': '%s GABLAM' % self.dataset()})
            if 0 < self.getInt('MaxUPC') < self.UPNum():
                self.printLog('#FUPC','Too many FUPC (%d) for MaxUPC setting (%d). No *.fupc.dis.tdt.' % (self.UPNum(),self.getInt('MaxUPC')))
            elif self.extras(1):
                idmatrix.saveMatrix(self.seqs(),filename='%s.fupc.dis.tdt' % upcbase,delimit='\t',format='text',log=True)
                if self.extras(2):
                    idmatrix.saveMatrix(self.seqs(),filename='%s.fupc.phydis.txt' % upcbase,format='phylip',log=True)
            else:
                self.printLog('#FUPC','No *.fupc.dis.tdt. (extras<1)')
            return True

        except: self.errorLog('Fatal Error in %s.makeFUPC(). Check for old results files with same name but different sequences.' % self.prog())
        return False
#########################################################################################################################
    ### <5> ### SLiM Calculation method                                                                                 #
#########################################################################################################################
    def addSLiMToList(self,slim):   ### Add slims from self.dict['Slim'] to self.obj['SlimList']
        '''Add slims from self.dict['Slim'] to self.obj['SlimList'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = rje.dictValues(self.dict['Slim'][slim],'SLiM','SLiM')
            if Motif: return Motif
            pattern = patternFromCode(slim)
            Motif = self.obj['SlimList']._addMotif(slim,pattern,reverse=False,check=False,logrem=True)
            self.dict['Slim'][slim]['SLiM'] = Motif

            ### ~ [2] Add Calculate Statistics and make MotifOcc objects ###
            for (Seq,pos) in self.dict['Slim'][slim]['Occ'][0:]:
                if Seq not in Motif.dict['Occ']: Motif.dict['Occ'][Seq] = []
                if pos < 0: pos = 0     ## Re-find motif in sequence to get proper object data ##
                sequence = Seq.info['PreMask'][pos:pos+Motif.slimLen()]     #x#stat['FullLength']]
                #x#self.deBug('%s %s : %s\n%s' % (pos, slim, sequence,Seq.info['PreMask'][:pos+Motif.slimLen()]))
                for occ in Motif.searchSequence(sequence=sequence):    # {Pos,Variant,Match,ID,MisMatch}
                    try:
                        occ['Start_Pos'] = occ['Pos'] + pos
                        if slim[0] == '^': occ['Start_Pos'] += 1
                        occ['End_Pos'] = occ['Start_Pos'] + len(occ['Match']) - 1
                        occ['Expect'] = self.dict['Slim'][slim]['ExpUP']
                    except:
                        print Seq.info['Sequence']
                        print Seq.info['PreMask']
                        print slim, pos, Seq.shortName(), sequence[:10], Motif.dict['Search']
                        occ = {'Match':'!ERR!','SearchVar':'!ERR!','Variant':'!ERR!'}
                    occ['Pos'] = pos+1
                    occ['Prot_Len'] = Seq.aaLen()
                    occ['Seq'] = Seq
                    occ['Motif'] = Motif
                    #self.deBug(occ)
                    Motif.dict['Occ'][Seq].append(occ)
            return Motif
        except:
            self.errorLog('Problem with %s.addSLiMToList(%s)' % (self.prog(),slim))
            return None
#########################################################################################################################
    def OLDaddSLiMToList(self,slim):   ### Add slims from self.dict['Slim'] to self.obj['SlimList']
        '''Add slims from self.dict['Slim'] to self.obj['SlimList'].'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = rje.dictValues(self.dict['Slim'][slim],'SLiM','SLiM')
            if Motif: return Motif
            pattern = patternFromCode(slim)
            Motif = self.obj['SlimList']._addMotif(slim,pattern,reverse=False,check=False,logrem=True)
            self.dict['Slim'][slim]['SLiM'] = Motif

            ### ~ [2] Add Calculate Statistics and make MotifOcc objects ###
            for (Seq,pos) in self.dict['Slim'][slim]['Occ'][0:]:
                if Seq not in Motif.dict['Occ']: Motif.dict['Occ'][Seq] = []
                if pos < 0: pos = 0     ## Re-find motif in sequence to get proper object data ##
                sequence = Seq.info['PreMask'][pos:pos+Motif.slimLen()]     #x#stat['FullLength']]
                #x#self.deBug('%s %s : %s\n%s' % (pos, slim, sequence,Seq.info['PreMask'][:pos+Motif.slimLen()]))
                try:
                    occ = Motif.searchSequence(sequence=sequence)[0]    # {Pos,Variant,Match,ID,MisMatch}
                    occ['Start_Pos'] = occ['Pos'] + pos
                    if slim[0] == '^': occ['Start_Pos'] += 1
                    occ['End_Pos'] = occ['Start_Pos'] + len(occ['Match']) - 1
                    occ['Expect'] = self.dict['Slim'][slim]['ExpUP']
                except:
                    print Seq.info['Sequence']
                    print Seq.info['PreMask']
                    print slim, pos, Seq.shortName(), sequence[:10], Motif.dict['Search']
                    occ = {'Match':'!ERR!','SearchVar':'!ERR!','Variant':'!ERR!'}
                occ['Pos'] = pos+1
                occ['Prot_Len'] = Seq.aaLen()
                occ['Seq'] = Seq
                occ['Motif'] = Motif
                #self.deBug(occ)
                Motif.dict['Occ'][Seq].append(occ)
            return Motif
        except:
            self.errorLog('Problem with %s.addSLiMToList(%s)' % (self.prog(),slim))
            return None
#########################################################################################################################
    def calculateSLiMOccStats(self):   ### Makes entries to SLiMList object and calculates attributes with slimcalc
        '''Makes entries to SLiMList object and calculates attributes with slimcalc.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for slim in self.dict['Slim']: self.addSLiMToList(slim)
            self.obj['SlimList'].calculateOccAttributes(silent=False)
            self.setBool({'OccStatsCalculated':True})
        except: self.errorLog('Problem with %s.calculateSLiMOccStats()' % self.prog())
#########################################################################################################################
    ### <6> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
    def makeAAFreq(self):   ### Makes an initial AAFreq dictionary containing AA counts only (including Xs)
        '''Makes an initial AAFreq dictionary containing AA counts only (including Xs).'''
        try:### ~ [1] ~ Load AA Frequencies for file and use for whole dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('AAFreq'):
                self.setBool({'MaskFreq':False})        # Over-riding AA Frequencies invalidates premasking
                self.setBool({'SmearFreq':False})       # Amino acids are already constant across UPC
                ## Read Freqs ##
                aafile = self.getStr('AAFreq')
                if os.path.exists(aafile) and open(aafile,'r').read()[:1] == '>':   # Fasta
                    aafreq = self.obj['SeqList'].aaFreq(fromfile=aafile,total=True)    # Returns proportions and total
                else: aafreq = self.obj['SeqList'].aaFreq(loadfile=aafile,total=True)    # Returns proportions and total
                ## Adjust to Counts ##
                for aa in aafreq.keys():
                    if aa != 'Total': aafreq[aa] = int(aafreq[aa]*aafreq['Total']+0.5)
                ## Copy to seqs and UPCs ##
                for dkey in ['Dataset'] + self.list['UP']: self.dict['AAFreq'][dkey] = copy.copy(aafreq)
                self.dict['AAFreq']['Dataset']['Total'] = self.obj['SeqList'].aaTotal(nonx=True)  ### Returns total number of AA in dataset
                self.printLog('#AAFREQ','Using aa frequencies from %s (no mask/smear freq)' % self.getStr('AAFreq'))
            ### ~ [2] ~ Calculate from Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                self.dict['AAFreq']['Dataset'] = {}
                for aa in self.list['Alphabet']: self.dict['AAFreq']['Dataset'][aa] = 0.0
                for upc in self.list['UP']:
                    self.dict['AAFreq'][upc] = {}
                    for aa in self.list['Alphabet']: self.dict['AAFreq'][upc][aa] = 0.0
                    for seq in upc:
                        seq.aaFreq(aafreq=self.dict['AAFreq'][upc])
                        seq.aaFreq(aafreq=self.dict['AAFreq']['Dataset'])
                ## Add Totals to AAFreq dict ## 
                for dkey in ['Dataset'] + self.list['UP']:    #!# Don't need totals? #!#
                    self.dict['AAFreq'][dkey]['Total'] = sum(self.dict['AAFreq'][dkey].values())
            ### ~ [3] ~ Convert Dataset into Freqs - not bothered with masking for this ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            rje.dictFreq(self.dict['AAFreq']['Dataset'])
        except: self.errorLog('Error in %s.makeAAFreq().' % self.prog());  raise
#########################################################################################################################
    def smearAAFreq(self,update=True):  ### Equalises AAFreq across UPC. Leaves Totals unchanged.
        '''Equalises AAFreq across UPC. Leaves Totals unchanged. Updates or returns smearfreq.'''
        try:
            ###~Setup dictionary containing mean AAFreq across UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            smearfreq = {}
            letters = {False:'ACDEFGHIKLMNPQRSTVWY^$',True:'ACGT^$'}[self.getBool('DNA')]
            for aa in letters:
                smearfreq[aa] = 0.0
                for upc in self.list['UP']: smearfreq[aa] += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)
            rje.dictFreq(smearfreq,total=False)
            if not update: return smearfreq
            ###~Copy frequencies to all AAFreq dictionaries~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for grp in self.dict['AAFreq']:
                for aa in letters: self.dict['AAFreq'][grp][aa] = smearfreq[aa]
            if self.getBool('DNA'): self.printLog('#NTFREQ','NT frequencies "smeared" over %d UPC' % self.UPNum())
            else: self.printLog('#AAFREQ','AA frequencies "smeared" over %d UPC' % self.UPNum())
        except:
            self.errorLog('Major problem during %s.smearAAFreq()' % self.prog())
            raise
#########################################################################################################################
    ### <7> ### Results Output Methods                                                                                  #
#########################################################################################################################
    def extras(self,level=1):     ### Returns whether extras settings permit output at this level
        '''Returns whether extras settings permit output at this level.'''
        try: return self.getInt('Extras') >= level
        except: return self.getBool('Extras')
#########################################################################################################################
    def searchBase(self,resdir='ResDir'): return '%s.%s.%s' % (self.seqBaseFile(resdir),self.buildText(),self.maskText())
#########################################################################################################################
    def runBase(self): return '%s.%s.%s' % (self.seqBaseFile(),self.buildText(),self.getStr('RunID'))
#########################################################################################################################
    def buildText(self): return self.prog()
#########################################################################################################################
    def maskText(self,joiner='',freq=True): ### Returns masking text
        '''Returns masking text.'''
        if self.getStrLC('MaskText'): return self.getStr('MaskText')
        if not 'MaskText' in self.list: 
            masking = []
            if self.getBool('Masking'):
                if self.getBool('ConsMask'): masking.append('Cons')
                if self.getBool('DisMask'): masking.append('Dis')
                if rje.matchExp('(\d+),(\d+)',self.getStr('CompMask')) and self.getStr('CompMask') != '0,0':
                    masking.append('Comp-%s-%s' % rje.matchExp('(\d+),(\d+)',self.getStr('CompMask')))
                if self.list['FTMask']:
                    if self.getBool('Webserver'):
                        ftxt = []
                        for f in self.list['FTMask']: ftxt.append(f[:1])
                        ftxt.sort()
                        masking.append('FT-%s' % string.join(ftxt,''))
                    else: masking.append('FT')
                if self.list['IMask']:
                    if self.getBool('Webserver'):
                        ftxt = []
                        for f in self.list['IMask']: ftxt.append(f[:1])
                        ftxt.sort()
                        masking.append('Inc-%s' % string.join(ftxt,''))
                    else: masking.append('Inc')
                if self.getStrLC('CaseMask'): masking.append('%s%s' % (self.getStr('CaseMask')[:1].upper(),self.getStr('CaseMask')[1:3].lower()))
                if self.getStrLC('MotifMask'): masking.append('Mot')
                if self.list['AAMask']: masking.append('AA')
                if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']: masking.append('QReg')
            if not masking: masking = ['NoMask']
            self.list['MaskText'] = masking
        if freq and self.getBool('MaskFreq') and 'Freq' not in self.list['MaskText']: self.list['MaskText'].insert(0,'Freq')
        if not freq and 'Freq' in self.list['MaskText']: self.list['MaskText'].remove('Freq')
        return string.join(self.list['MaskText'],joiner)
#########################################################################################################################
    def tidyMotifObjects(self): ### This should not be necessary but somehow is!
        '''This should not be necessary but somehow is!'''
        if self.dev(): self.warnLog('tidyMotifObjects() is being skipped and should not be required!')
        return #!# Try without this! #!#

        #>>>

        ### Tidy MotifList Object ###
        for Motif in self.obj['MotifList'].motifs()[0:]:
            slim = Motif.slimCode() 
            if not slim or slim not in self.list['SigSlim']: self.obj['MotifList'].removeMotif(Motif)
        ### Check completeness of MotifObject ###
        for slim in self.list['SigSlim']:   #!# Add SigNum? #!#
            pattern = patternFromCode(slim)
            Motif = self.obj['MotifList'].mapPattern(pattern,update=False)
            if not Motif in self.obj['MotifList'].motifs():
                self.printLog('#MOT','Motif "%s" missing from Motif Occurrence dictionary.' % pattern)
                try:
                    Motif = self.obj['MotifList'].mapPattern(pattern,update=True)
                    Motif = self.motifOccStats(slim)
                    Motif.stat['OccNum'] = self.slimOccNum(slim)
                    Motif.stat['OccSeq'] = self.slimUP(slim)
                    Motif.dict['Expect'] = {self.basefile():self.dict['Slim'][slim]['ExpUP']}
                    self.obj['MotifList'].combMotifOccStats(statlist=self.list['OccStats'],revlist=['Hyd'],log=False,motiflist=[Motif])
                    self.printLog('#MOT','Successfully added Motif Object for "%s" (%s)!' % (slim,pattern))
                except: self.errorLog('Cannot add/process Motif Object for "%s" (%s)!' % (slim,pattern))
#########################################################################################################################
    def abNprob(self,a,b,N,overlap=0,dirn='less'):     ### Returns probabilities of overlap in a given b/N and in b given a/N
        '''Returns probabilities of no overlap in a given b/N and in b given a/N.'''
        p = 0.0
        if dirn == 'less':
            for k in range(overlap+1):
                p += rje.binomial(k,a,float(b)/N,exact=True,callobj=self)
                p += rje.binomial(k,b,float(a)/N,exact=True,callobj=self)
        else: p += rje.binomial(overlap,b,float(a)/N,callobj=self) + rje.binomial(overlap,a,float(b)/N,callobj=self)
        return p / 2.0
#########################################################################################################################
#    def myPickle(self):  ### Returns pickle identifier, also used for Outputs "Build" column (self.info['Build'])
#        '''Returns pickle identifier, also used for Outputs "Build" column.'''
#        ## Pickle name ##
#        return '%s.%s' % (self.buildText(),self.maskText())
#########################################################################################################################
    def tarZipSaveSpace(self):  ### Tars and Zips output and/or deletes extra files as appropriate
        '''Tars and Zips output and/or deletes extra files as appropriate.'''
        try:### ~ [0] Setup file lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #targz = '%s.%s.tgz' % (self.info['Basefile'],self.info['RunID'])    #self.myPickle())
            targz = '%s.tgz' % self.runBase()
            ## ~ [0a] Pickle for this run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #mypickle = '%s.%s.pickle' % (self.info['Basefile'],self.myPickle())
            mypickle = '%s.pickle' % self.searchBase()
            gzpickle = '%s.gz' % mypickle
            if rje.isYounger(mypickle,gzpickle) == mypickle: os.unlink(gzpickle)
            if rje.isYounger(gzpickle,mypickle) == gzpickle: os.unlink(mypickle)
            if os.path.exists(mypickle) and not self.getBool('Win32'):
                try:
                    os.system('gzip %s' % (mypickle))
                    self.printLog('#GZIP','%s %s zipped.' % (self.prog(),mypickle))
                except: self.errorLog('Cannot gzip %s' % (mypickle))
            ## ~ [0b] All dataset files in directory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            files = glob.glob('%s*' % self.seqBaseFile())   #info['Basefile'])
            #!# Need to rationalise this now that different basefiles are used at different points
            ### ~ [1] Tar and Zip Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('TarGZ'):
                ## ~ [1a] Check Windows ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getBool('Win32'):
                    self.printLog('#TAR','Sorry! TarGZ option not available for Windows yet!')
                    if self.getInt('SaveSpace') > 0: self.printLog('#DEL','No files deleted to save space. Control output with extras=X.')
                    return
                ## ~ [1b] Delete existing archive ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if os.path.exists(targz):
                    os.unlink(targz)
                    self.printLog('#TAR','Deleted old TAR archive.')
                ## ~ [1c] Make archive ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #tcmd = 'tar -czf %s %s*' % (targz,self.info['Basefile'])
                tarfiles = []
                for tfile in files[0:]:
                    if tfile.endswith('tar.gz') or tfile.endswith('.tgz'): continue
                    if tfile.find('.pickle') > 0: continue
                    tarfiles.append(tfile)
                if self.getInt('SaveSpace') > 1 and os.path.exists(mypickle): tarfiles.append(mypickle)
                if self.getInt('SaveSpace') > 1 and os.path.exists(gzpickle): tarfiles.append(gzpickle)
                if tarfiles:
                    tcmd = 'tar -czf %s %s' % (targz,string.join(tarfiles))
                    self.printLog('#TAR',tcmd)
                    if os.system(tcmd):
                        self.printLog('#ERR','Problem executing TarGZ option for %s!' % self.dataset())
                        if self.getInt('SaveSpace') > 0: self.printLog('#DEL','No files deleted to save space.')
                        return
            ### ~ [2] Delete files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getInt('SaveSpace') > 0:
                keeplist = ['tar.gz','tgz']
                if self.getInt('SaveSpace') < 2:  keeplist += ['pickle*']
                if self.getInt('SaveSpace') < 3:  keeplist += ['upc']
                ## ~ [2a] Setup File List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for keepext in keeplist:
                    keepers = glob.glob('%s*%s' % (self.seqBaseFile(),keepext))
                    for k in keepers:
                        if k in files: files.remove(k)
                self.printLog('#DEL','Deleting of %d files/dirs to save space...' % len(files),newline=False,log=False)
                ## ~ [2b] Delete Files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for dfile in files:
                    if os.path.isdir(dfile): rje.deleteDir(self,dfile,contentsonly=False,confirm=False,report=False)
                    else: os.remove(dfile)
                self.printLog('\r#DEL','Deletion of %d files/dirs to save space complete.' % len(files))
        except: self.errorLog('Problem with %s.tarZipSaveSpace()' % self.prog())
#########################################################################################################################
    ### <10> ### MegaSLiM Functions                                                                                     #
#########################################################################################################################
    def megaSLiM(self): ### Sets up MegaSLiM data (GABLAM and masking) for subsequent use
        '''
        Sets up MegaSLiM data (GABLAM and masking) for subsequent use. Will run GABLAM all-by-all if needed (or force=T)
        and save masking program/method raw scores (disorder and conservation) as FILE.X.tdt using sequence shortname and
        simple space delimited scores. Scores will be saved to an accuracy determined by self.int['MegaSLiMdp']
        '''
        try:### ~ [0] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing input file!
            scmd = self.cmd_list + ['seqin=%s' % seqfile,'autoload=T','seqmode=file','seqindex=T']
            seqlist = rje_seqlist.SeqList(self.log,scmd)
            ### ~ [1] GABLAM All-by-all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('MegaGABLAM'):
                gcmd = ['local=F','dismat=F','distrees=F','disgraph=F','clusters=F'] + scmd + ['searchdb=%s' % seqfile,'fullres=T','hitsum=T','qryacc=F','basefile=%s' % basefile]
                if self.force() or not rje.exists('%s.hitsum.tdt' % basefile):  # Run full GABLAM
                    if not gablam.GABLAM(self.log,gcmd+['append=F']).gablam(): raise ValueError('GABLAM run failure!')
                elif not self.megaSLiMGABLAM(): raise ValueError('GABLAM run failure!')
            ### ~ [2] Masking raw scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.prog() in ['rlc','disorder']:
                #self.printLog('#CMD','IUScoreDir: %s; ProtScores: %s' % (self.getStr('IUScoreDir'),self.getBool('ProtScores')))
                pass  # Effectively, masking=T
            elif not self.getBool('Masking'): return True
            ## ~ [2a] Disorder Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('DisMask') or self.getBool('ConsMask'): self.megaSLiMDisorder(seqlist)
            ## ~ [2b] Conservation Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.getBool('ConsMask'): self.megaSLiMConservation(seqlist)
            ### ~ [3] Special REST functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Add these to output
            return True
        except: self.errorLog('MegaSLiM failure.'); return False
#########################################################################################################################
    def megaSLiMGABLAM(self,seqfile=None):  ### Creates/updates GABLAM for seqlist sequences.
        '''
        Creates/updates GABLAM for seqlist sequences.
        >> seqlist:obj [None] = rje_seqlist.SeqList object with sequences requiring GABLAM. Uses MegaSLiM if None.
        '''
        try:### ~ [0] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqfile: seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing input file!
            scmd = self.cmd_list + ['seqin=%s' % seqfile,'autoload=T','seqmode=file','seqindex=T']
            seqlist = rje_seqlist.SeqList(self.log,scmd)
            ### ~ [1] GABLAM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gcmd = ['local=F','dismat=F','distrees=F','disgraph=F','clusters=F'] + scmd + ['searchdb=%s' % self.getStr('MegaSLiM'),'fullres=T','hitsum=T','qryacc=F','basefile=%s' % basefile]
            missing = []
            db = rje_db.Database(self.log,self.cmd_list)
            hdb = db.addTable('%s.hitsum.tdt' % basefile,['Qry'])
            while seqlist.nextSeq():
                sname = string.split(seqlist.currSeq()[0])[0]   # Shortname
                if not hdb.data(sname): missing.append(sname)
            if missing: return gablam.GABLAM(self.log,gcmd+['append=T','missing=%s' % string.join(missing,',')]).gablam()
            self.printLog('#GABLAM','All input sequences found in %s.hitsum.tdt' % basefile)
            return True
        except: self.errorLog('%s.megaSLiMGABLAM() failure.' % self.prog()); return False
#########################################################################################################################
    def megaSLiMDisorder(self,seqlist=None):    ### Creates/updates disorder score lists for seqlist sequences.
        '''
        Creates/updates disorder score lists for seqlist sequences, saved to an accuracy determined by MegaSLiMdp.
        >> seqlist:obj [None] = rje_seqlist.SeqList object with sequences requiring scores. Uses MegaSLiM if None.
        '''
        try:### ~ [0] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing input file!
            self.dict['Output']['disorder'] = ''
            if not seqlist:
                scmd = self.cmd_list + ['seqin=%s' % seqfile,'autoload=T','seqmode=file','seqindex=T']
                seqlist = rje_seqlist.SeqList(self.log,scmd)
            ### ~ [1] Disorder Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disobj = rje_disorder.Disorder(self.log,self.cmd_list)
            dkey = disobj.getStrLC('Disorder'); wkey = '%s-write' % dkey
            if dkey == 'iupred' and disobj.getStrLC('IUMethod') == 'long': dkey = 'iulong'
            readfile = writefile = '%s.%s.txt' % (basefile,dkey)
            if self.getBool('MegaSLiMFix'):
                readfile = '%s.%s.fix.txt' % (basefile,dkey)
                if rje.exists(writefile): os.rename(writefile,readfile)
            self.setStr({dkey:readfile,wkey:writefile})
            force = self.force() or not rje.exists(readfile)
            if force or not rje.exists(writefile): RFILE = open(writefile,'w')
            else:
                RFILE = open(writefile,'a')
                if open(writefile,'r').read()[:-1] != '\n': RFILE.write('\n')   # Try to fix odd crash concatenation issue
            self.file[wkey] = RFILE
            ## ~ [2b] Scan and update from sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dp = self.getInt('MegaSLiMdp'); dx = 0; ox = 0; sx = 0; stot = seqlist.seqNum()
            seqlist.obj['Current'] = None
            while seqlist.nextSeq():
                self.progLog('\r#DIS','Scanning %s for %s sequences: %.1f%% (%s missing/run)' % (readfile,rje.iStr(stot),sx/stot,rje.iStr(dx))); sx += 100.0
                sname = string.split(seqlist.currSeq()[0])[0]   # Shortname
                acc = seqlist.seqAcc()
                sequence = seqlist.currSeq()[1]
                if force: fline = ''
                else: fline = self.findline(dkey,sname)
                disobj.list['ResidueDisorder'] = []
                if fline:
                    for dstr in string.split(fline)[1:]: disobj.list['ResidueDisorder'].append(string.atof(dstr))
                    if len(disobj.list['ResidueDisorder']) != len(sequence):
                        self.errorLog('Disorder score length mismatch (%d score vs %d pos)' % (len(disobj.list['ResidueDisorder']),len(sequence)),printerror=False)
                        fline = ''; disobj.list['ResidueDisorder'] = []
                elif not self.getBool('MegaSLiMFix'): force = True  # Reached end of reading. Assumes same order.

                if not fline and self.getStr('IUScoreDir'):    # Look for individual result
                    ifile = '%s%s.%s.txt' % (self.getStr('IUScoreDir'),acc,dkey)
                    if os.path.exists(ifile): fline = open(ifile,'r').readline()
                    for dstr in string.split(fline)[1:]: disobj.list['ResidueDisorder'].append(string.atof(dstr))
                    if len(disobj.list['ResidueDisorder']) != len(sequence):
                        self.errorLog('%s Disorder score length mismatch (%d score vs %d pos)' % (sname,len(disobj.list['ResidueDisorder']),len(sequence)),printerror=False)
                        fline = ''; disobj.list['ResidueDisorder'] = []

                if not fline:   # Calculate
                    disobj.disorder(sequence=sequence,name=sname); dx += 1

                if self.getBool('MegaSLiMFix') or not fline:    # Append
                    if disobj.list['ResidueDisorder']:  # Should have scores
                        dlist = []
                        for x in disobj.list['ResidueDisorder']:
                            if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                            else: dlist.append('%f' % x)
                        RFILE.write('%s\t%s\n' % (sname,string.join(dlist))); ox += 1
                        if self.getStr('IUScoreDir'):
                            rje.mkDir(self,ifile)
                            open(ifile,'w').write('%s\t%s\n' % (sname,string.join(dlist)))
                            fline = '%s\t%s\n' % (sname,string.join(dlist))
                    else: raise ValueError('Disorder scoring for %s failed!' % sname)
                if self.getStr('IUScoreDir'): self.dict['Output']['disorder'] += fline
            self.printLog('\r#DIS','Scanned %s for %s sequences: %.1f%% (%s missing/run).' % (readfile,rje.iStr(stot),sx/stot,rje.iStr(dx)))
            if ox or self.getBool('MegaSLiMFix'): self.printLog('#%s' % dkey.upper(),'%s disorder score lists output to %s for %s sequences.' % (dkey,writefile,rje.iStr(ox)))
            RFILE.close()
            if self.getBool('MegaSLiMFix') and rje.exists(readfile): os.unlink(readfile)
        except: self.errorLog('%s.megaSLiMDisorder() failure.' % self.prog()); raise
#########################################################################################################################
    def megaSLiMSeqDisorderObj(self,sname,sequence):    ### Returns a disorder object for given (name,sequence).
        '''
        Returns a disorder object for given (sname,sequence).
        '''
        try:### ~ [0] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing MegaSLiM file!
            ## ~ [0a] Disorder Object/Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            disobj = rje_disorder.Disorder(self.log,self.cmd_list)
            dkey = disobj.getStrLC('Disorder')
            if dkey == 'iupred' and disobj.getStrLC('IUMethod') == 'long': dkey = 'iulong'
            dfile = '%s.%s.txt' % (basefile,dkey)
            self.setStr({dkey:dfile})
            ### ~ [1] Scan and update disorder & mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dp = self.getInt('MegaSLiMdp')
            ## ~ [1a] Update MegaSLiM disorder file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            fline = self.findline(dkey,sname)
            if fline:
                disobj.list['ResidueDisorder'] = []
                for dstr in string.split(fline)[1:]: disobj.list['ResidueDisorder'].append(string.atof(dstr))
                if len(disobj.list['ResidueDisorder']) != len(sequence):
                    self.debug('%s: %s' % (sname,sequence))
                    self.errorLog('%s disorder score length mismatch (%d score vs %d pos)! Will regenerate.' % (sname,len(disobj.list['ResidueDisorder']),len(sequence)))
                    fline = ''; disobj.list['ResidueDisorder'] = []
            if not fline:   # Append!
                disobj.disorder(sequence=sequence,name=sname)
                if disobj.list['ResidueDisorder']:  # Should have scores
                    if len(disobj.list['ResidueDisorder']) != len(sequence): raise ValueError('%s disorder score length mismatch! (%d score vs %d pos)!' % (sname,len(disobj.list['ResidueDisorder']),len(sequence)))
                    dlist = []
                    for x in disobj.list['ResidueDisorder']:
                        if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                        else: dlist.append('%f' % x)
                    open(dfile,'a').write('%s\t%s\n' % (sname,string.join(dlist)))
                else: raise ValueError('Disorder scoring for %s failed!' % sname)
            ## ~ [1b] Make RegionDisorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            disobj.makeRegions(sequence)
            return disobj
        except: self.errorLog('%s.megaSLiMSeqDisorderObj() failure.' % self.prog()); raise
#########################################################################################################################
    def megaSLiMMaskDisorder(self): ### Uses disorder score lists for sequence masking, updating if required.
        '''
        Uses disorder score lists for sequence masking, updating if required.
        '''
        try:### ~ [0] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing MegaSLiM file!
            self.dict['Output']['disorder'] = ''
            ## ~ [0a] Disorder Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            disobj = rje_disorder.Disorder(self.log,self.cmd_list)
            dkey = disobj.getStrLC('Disorder')
            if dkey == 'iupred' and disobj.getStrLC('IUMethod') == 'long': dkey = 'iulong'
            dfile = '%s.%s.txt' % (basefile,dkey)
            self.printLog('#DIS','Disorder score file: %s' % dfile)
            self.setStr({dkey:dfile})
            force = self.force() or not rje.exists(dfile)
            if force: DFILE = open(dfile,'w')
            else: DFILE = open(dfile,'a')
            ### ~ [1] Scan and update disorder & mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dp = self.getInt('MegaSLiMdp'); dx = 0
            for seq in self.seqs():
                ## ~ [1a] Update MegaSLiM disorder file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sname = seq.shortName()
                acc = seq.getStr('AccNum')
                if force: fline = ''
                else: fline = self.findline(dkey,sname)
                disobj.list['ResidueDisorder'] = []
                if fline:
                    for dstr in string.split(fline)[1:]: disobj.list['ResidueDisorder'].append(string.atof(dstr))
                    if len(disobj.list['ResidueDisorder']) != seq.aaLen():
                        self.errorLog('%s Disorder score length mismatch (%d score vs %d pos)' % (sname,len(disobj.list['ResidueDisorder']),seq.aaLen()),printerror=False)
                        fline = ''; disobj.list['ResidueDisorder'] = []
                if not fline and self.getStr('IUScoreDir'):    # Look for individual result
                    ifile = '%s%s.%s.txt' % (self.getStr('IUScoreDir'),acc,dkey)
                    self.bugPrint(ifile)
                    if os.path.exists(ifile):
                        fline = open(ifile,'r').readline()
                        for dstr in string.split(fline)[1:]: disobj.list['ResidueDisorder'].append(string.atof(dstr))
                        if len(disobj.list['ResidueDisorder']) != seq.aaLen():
                            self.errorLog('%s Disorder score length mismatch (%d score vs %d pos)' % (sname,len(disobj.list['ResidueDisorder']),seq.aaLen()),printerror=False)
                            fline = ''; disobj.list['ResidueDisorder'] = []
                if not fline:   # Calculate & Append
                    disobj.disorder(sequence=seq.info['PreMask'],name=sname)
                    if len(disobj.list['ResidueDisorder']) == seq.aaLen():  # Should have scores
                        dlist = []
                        for x in disobj.list['ResidueDisorder']:
                            if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                            else: dlist.append('%f' % x)
                        DFILE.write('%s\t%s\n' % (sname,string.join(dlist))); dx += 1
                        if self.getStr('IUScoreDir'):
                            rje.mkDir(self,ifile)
                            open(ifile,'w').write('%s\t%s\n' % (sname,string.join(dlist)))
                            fline = '%s\t%s\n' % (sname,string.join(dlist))
                    else: raise ValueError('Disorder scoring for %s failed!' % sname)
                #?# Not sure if this is the best way to work out with REST output is required
                if self.getStr('IUScoreDir'): self.dict['Output']['disorder'] += fline
                ## ~ [1b] Make RegionDisorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try: disobj.makeRegions(seq.info['PreMask'])
                except: raise ValueError('Disorder scoring for %s failed!' % sname)
                ## ~ [1c] Mask sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                mask='X'
                oldseq = seq.info['Sequence'][0:]
                newseq = mask * len(oldseq)
                prex = oldseq.count(mask)
                (newseq,oldseq) = (oldseq,newseq)
                for region in disobj.list['RegionDisorder']: oldseq = oldseq[:region[0]-1] + newseq[region[0]-1:region[1]] + oldseq[region[1]:]
                seq.info['Sequence'] = oldseq
                maskx = oldseq.count(mask) - prex
                if maskx > 0:
                    mtxt = '%s masked %d ' % (seq.shortName(),len(disobj.list['RegionFold']))
                    self.printLog('#MASK','%sordered regions. (%d %s added.)' % (mtxt,maskx,mask),screen=self.getBool('LogMask'))
            if dx: self.printLog('#%s' % dkey.upper(),'%s disorder score lists output to %s for %s sequences.' % (dkey,dfile,rje.iStr(dx)))
            DFILE.close()
        except: self.errorLog('%s.megaSLiMMaskDisorder() failure.' % self.prog()); raise
#########################################################################################################################
    def megaSLiMMaskConservation(self):    ### Uses RLC score lists for sequence masking, updating if required.
        '''
        Uses RLC score lists for sequence masking, updating if required.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
             ## ~ [0a] Setup MegaSLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing input file!
            self.dict['Output']['rlc'] = ''
            ## ~ [0b] Setup SLiMCalc object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc=','usealn=T'])
            slimcalc.loadBLOSUM()
            # Note. MegaSLiM masking will not check full GOPHER using forking.
            ## ~ [0c] RLC Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rfile = '%s.rlc.txt' % (basefile); dkey = 'rlc'; wkey = '%s-write' % dkey
            self.setStr({dkey:rfile,wkey:rfile})
            force = self.force() or not rje.exists(rfile)
            if force or not rje.exists(rfile): RFILE = open(rfile,'w')
            else:
                RFILE = open(rfile,'a')
                if open(rfile,'r').read()[:-1] != '\n': RFILE.write('\n')   # Try to fix odd crash concatenation issue
            self.file['%s-write' % dkey] = RFILE
            ### ~ [1] Scan and update from sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dp = self.getInt('MegaSLiMdp'); dx = 0
            cmx = 0
            for seq in self.seqs():     #!# Could return sorted seqlist to accelerate process?
                ## ~ [1a] Update MegaSLiM RLC score file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                (seq.info['MaskSeq'],seq.info['Sequence']) = (seq.info['Sequence'],seq.info['PreMask'])
                sname = seq.shortName()
                if self.getBool('SLiMMutant') and rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname): sname = rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname)[0]
                if force: fline = ''
                else: fline = self.findline(dkey,sname)
                relcon = []
                if fline:
                    for dstr in string.split(fline)[1:]: relcon.append(string.atof(dstr))
                    if len(relcon) != seq.aaLen():
                        self.errorLog('MegaSLiM %s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),seq.aaLen()),printerror=False)
                        fline = ''; relcon = []
                if not fline and self.getBool('ProtScores'):
                    sfile = '%s.rlc.txt' % (rje.baseFile(slimcalc.loadOrthAln(seq,alnfile=True)))
                    if os.path.exists(sfile):
                        fline = open(sfile,'r').readline()
                        for dstr in string.split(fline)[1:]: relcon.append(string.atof(dstr))
                        if len(relcon) != seq.aaLen():
                            self.errorLog('ProtScore %s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),seq.aaLen()),printerror=False)
                            fline = ''; relcon = []
                if not fline:
                    seq.obj['Disorder'] = self.megaSLiMSeqDisorderObj(sname,seq.info['PreMask'])
                    relcon = slimcalc.relConListFromSeq(seq,window=30,screen=self.v()>1)  #!# Need to use MegaSLiM disorder here too! #!#
                    if relcon:
                        dlist = []
                        for x in relcon:
                            if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                            else: dlist.append('%f' % x)
                        RFILE.write('%s\t%s\n' % (sname,string.join(dlist))); dx += 1
                        if self.getBool('ProtScores'):
                            open(sfile,'w').write('%s\t%s\n' % (sname,string.join(dlist)))
                            fline = '%s\t%s\n' % (sname,string.join(dlist))
                    else: raise ValueError('RLC scoring for %s failed!' % sname)
                #?# Not sure if this is the best way to work out with REST output is required
                if self.getBool('ProtScores'): self.dict['Output']['rlc'] += fline
                ## ~ [1b] Mask sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                newseq = []
                try:
                    for i in range(seq.seqLen()):
                        if relcon[i] < 0: newseq.append('X')
                        else: newseq.append(seq.info['MaskSeq'][i])
                    maskx = newseq.count('X') - seq.info['MaskSeq'].count('X'); cmx += maskx
                    seq.info['Sequence'] = string.join(newseq,'')
                    if maskx > 0 and self.getBool('LogMask'): self.printLog('#MASK','Masked %s low relative conservation. (%d X added to %daa seq.)' % (seq.shortName(),maskx,seq.aaLen()))
                except:
                    self.errorLog('Problem with relative conservation masking for %s' % seq.shortName())
                    self.printLog('#REL','%daa (%d) => %d conscores' % (seq.seqLen(),seq.aaLen(),len(relcon)))
                    seq.info['Sequence'] = seq.info['MaskSeq']
                self.wallTime()
            self.printLog('#REL','%s aa masked using relative conservation' % rje.integerString(cmx))
            if dx: self.printLog('#RLC','%s RLC score lists output to %s for %s sequences.' % (dkey,rfile,rje.iStr(dx)))
            RFILE.close()
        except: self.errorLog('%s.megaSLiMConservation() failure.' % self.prog()); raise
#########################################################################################################################
    def megaSLiMConservation(self,seqlist=None):    ### Creates/updates RLC score lists for seqlist sequences.
        '''
        Creates/updates disorder score lists for seqlist sequences, saved to an accuracy determined by MegaSLiMdp.
        >> seqlist:obj [None] = rje_seqlist.SeqList object with sequences requiring scores. Uses MegaSLiM if None.
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] Setup SeqList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqfile = self.getStr('MegaSLiM')
            basefile = rje.baseFile(seqfile)    # NOTE: Output will go to path containing input file!
            self.dict['Output']['rlc'] = ''
            if not seqlist:
                scmd = self.cmd_list + ['seqin=%s' % seqfile,'autoload=T','seqmode=file','seqindex=T']
                seqlist = rje_seqlist.SeqList(self.log,scmd)
            ## ~ [0b] Setup SLiMCalc object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc=','usealn=T'])
            slimcalc.loadBLOSUM()
            ## ~ [0c] Pre-run GOPHER if forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if slimcalc.getBool('UseGopher') and not slimcalc.getStr('AlnDir') and slimcalc.obj['Gopher'].getInt('Forks') > 1 and not slimcalc.obj['Gopher'].getBool('NoForks'):
                slimcalc.obj['Gopher'].setStr({'Name':seqfile})
                slimcalc.obj['Gopher'].cmd_list.append('orthaln')
                slimcalc.obj['Gopher'].run()
            ### ~ [1] RLC Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            readfile = writefile = '%s.rlc.txt' % (basefile); dkey = 'rlc'; wkey = '%s-write' % dkey
            if self.getBool('MegaSLiMFix'):
                readfile = '%s.rlc.fix.txt' % (basefile)
                if rje.exists(writefile): os.rename(writefile,readfile)
            self.setStr({dkey:readfile,wkey:writefile})
            force = self.force() or not rje.exists(readfile)
            if force or not rje.exists(writefile): RFILE = open(writefile,'w')
            else:
                RFILE = open(writefile,'a')
                if open(writefile,'r').read()[:-1] != '\n': RFILE.write('\n')   # Try to fix odd crash concatenation issue
            self.file['%s-write' % dkey] = RFILE
            ## ~ [1a] Scan and update from sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dp = self.getInt('MegaSLiMdp'); dx = 0; ox = 0; sx = 0; stot = seqlist.seqNum()
            seqlist.obj['Current'] = None
            while seqlist.nextSeq():    #!# Could use sorted seqlist to accelerate process?
                self.progLog('\r#RLC','Scanning %s for %s sequences: %.1f%% (%s missing/run)' % (readfile,rje.iStr(stot),sx/stot,rje.iStr(dx))); sx += 100.0
                (sname,sequence) = seqlist.currSeq()
                sname = string.split(sname)[0]   # Shortname
                if self.getBool('SLiMMutant') and rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname): sname = rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname)[0]

                if force: fline = ''
                else: fline = self.findline(dkey,sname)
                relcon = []
                if fline:
                    for dstr in string.split(fline)[1:]: relcon.append(string.atof(dstr))
                    if len(relcon) != len(sequence):
                        self.errorLog('%s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),len(sequence)),printerror=False)
                        self.debug('%s: %s' % (sname,sequence))
                        fline = ''
                if not fline: seq = seqlist.getSeqObj(sname,sequence)

                if not fline and self.getBool('ProtScores'):
                    sfile = '%s.rlc.txt' % (rje.baseFile(slimcalc.loadOrthAln(seq,alnfile=True)))
                    if os.path.exists(sfile):
                        fline = open(sfile,'r').readline()
                        for dstr in string.split(fline)[1:]: relcon.append(string.atof(dstr))
                        if len(relcon) != seq.aaLen():
                            self.errorLog('ProtScore %s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),seq.aaLen()),printerror=False)
                            fline = ''; relcon = []

                if not fline:
                    seq.obj['Disorder'] = self.megaSLiMSeqDisorderObj(sname,sequence)
                    relcon = slimcalc.relConListFromSeq(seq,window=30,screen=self.v()>1); dx += 1

                if self.getBool('MegaSLiMFix') or not fline:
                    if relcon:
                        dlist = []
                        for x in relcon:
                            if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                            else: dlist.append('%f' % x)
                        RFILE.write('%s\t%s\n' % (sname,string.join(dlist))); ox += 1
                        if self.getBool('ProtScores'):
                            open(sfile,'w').write('%s\t%s\n' % (sname,string.join(dlist)))
                            fline = '%s\t%s\n' % (sname,string.join(dlist))
                    else: raise ValueError('RLC scoring for %s failed!' % sname)
                if self.getBool('ProtScores'): self.dict['Output']['rlc'] += fline

            self.printLog('\r#RLC','Scanned %s for %s sequences: %.1f%% (%s missing/run).' % (readfile,rje.iStr(stot),sx/stot,rje.iStr(dx)))
            if ox or self.getBool('MegaSLiMFix'): self.printLog('#RLC','%s RLC score lists output to %s for %s sequences.' % (dkey,writefile,rje.iStr(ox)))
            RFILE.close()
            if self.getBool('MegaSLiMFix') and rje.exists(readfile): os.unlink(readfile)
        except: self.errorLog('%s.megaSLiMConservation() failure.' % self.prog()); raise
#########################################################################################################################
    ### <11> ### Single sequence masking                                                                                #
#########################################################################################################################
    def maskSequence(self,sequence,name=None,log=True,screen=True,mask='X'):   ### Masks single input sequence and returns, replacing masked regions with Xs
        '''Masks single input sequence and returns, replacing masked regions with Xs.'''
        ### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not name: name = 'sequence'
        if not self.getBool('Masking'): return sequence
        sname = string.split(name)[0]
        if self.getBool('SLiMMutant') and rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname): sname = rje.matchExp('^(\S+)\.p\.\D\d+\D$',sname)[0]
        nomask = sequence[0:]   # Used for disorder prediction etc.
        masked = []     # List of 'Type: %d X' to be joined.
        ## ~ [0a] ~ UniFake ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.getBool('FakeMask'):
            self.warnLog('UniFake masking not permitted for single sequences. (Switched off.)')
            self.setBool({'FakeMask':False})
        ### ~ [1] ~ Case Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getStrLC('CaseMask'):
            premask = sequence[0:]
            case = self.getStrLC('CaseMask')[:1]
            if case not in 'ul': raise ValueError('Do not recognise casemask=%s' % self.getStr('CaseMask'))
            sequence = ''; i = 0
            while i < len(premask):
                if case == 'u' and premask[i] == premask[i].upper(): sequence += mask
                elif case == 'l' and premask[i] == premask[i].lower(): sequence += mask
                else: sequence += premask[i].upper()
            maskx = sequence.count(mask) - premask.upper().count(mask)
            if maskx: masked.append('CaseMask (%sC): %d %s' % (case.upper(),maskx,mask))
        ### ~ [2] ~ Disorder (Inclusive Filtering) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        disobj = None   # Could be used for RLC and disorder masking
        if self.getBool('DisMask'):
            if self.getStrLC('MegaSLiM'): disobj = self.megaSLiMSeqDisorderObj(sname,nomask)
            else: disobj = self.seqDisorderObj(sname,nomask)
            premask = sequence[0:]
            maskseq = mask * len(sequence)
            for region in disobj.list['RegionDisorder']: maskseq = maskseq[:region[0]-1] + premask[region[0]-1:region[1]] + maskseq[region[1]:]
            sequence = maskseq
            maskx = sequence.count(mask) - premask.count(mask)
            if maskx: masked.append('DisMask (%d regions): %d %s' % (len(disobj.list['RegionFold']),maskx,mask))

        ### ~ [3] ~ Relative conservation masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not self.obj['SLiMCalc']:
            slimcalc = self.obj['SLiMCalc'] = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc=','usealn=T'])
            slimcalc.loadBLOSUM()
        else: slimcalc = self.obj['SLiMCalc']
        if self.getBool('ConsMask'):
            premask = sequence[0:]
            if not disobj:
                if self.getStrLC('MegaSLiM'): disobj = self.megaSLiMSeqDisorderObj(sname,nomask)
                else: disobj = self.seqDisorderObj(sname,nomask)
            relcon = []
            if self.getStrLC('MegaSLiM'):
                if not self.getStr('rlc',default=''):
                    self.setStr({'rlc':'%s.rlc.txt' % self.getStrBase('MegaSLiM')})
                dp = self.getInt('MegaSLiMdp')
                if rje.exists(self.getStr('rlc')): fline = self.findline('rlc',sname)
                else: fline = ''
                if fline:
                    for dstr in string.split(fline)[1:]: relcon.append(string.atof(dstr))
                    if len(relcon) != len(sequence):
                        self.errorLog('%s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),len(sequence)),printerror=False)
                        fline = ''
            else: fline = ''
            if not fline:
                seq = self.getSeqObj(sname,nomask)
                seq.obj['Disorder'] = disobj
                relcon = slimcalc.relConListFromSeq(seq,window=30)
                if relcon and self.getStrLC('MegaSLiM'):
                    dlist = []
                    for x in relcon:
                        if dp > 0: dlist.append('%s' % rje.dp(x,dp))
                        else: dlist.append('%f' % x)
                    open(self.getStr('rlc'),'a').write('%s\t%s\n' % (sname,string.join(dlist)))
            if len(relcon) != len(nomask):
                self.errorLog('%s RLC score length mismatch (%d score vs %d pos)' % (sname,len(relcon),len(nomask)),printerror=False)
                relcon = []
            ## ~ [1b] Mask sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not relcon:
                warntype = 'RLC Failure'
                self.warnLog('RLC scoring for %s failed!' % sname,warntype,quitchoice=warntype not in self.log.no_suppression,suppress=False)
                if warntype not in self.log.no_suppression: self.log.no_suppression.append(warntype)
            else:
                sequence = ''
                for i in range(len(premask)):
                    if relcon[i] < 0: sequence += mask
                    else: sequence += premask[i]
                maskx = sequence.count(mask) - premask.count(mask)
                if maskx: masked.append('ConsMask (RLC): %d %s' % (maskx,mask))
        ### ~ [4] ~ UniProt Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.list['IMask'] or self.list['FTMask']:
            self.warnLog('Uniprot masking (imask/ftmask) not permitted for single sequences. (Switched off.)')
            self.list['IMask'] = []; self.list['FTMask'] = []
        ### ~ [5] ~ Low Complexity Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if rje.matchExp('(\d+),(\d+)',self.getStr('CompMask')):
            (x,y) = rje.matchExp('(\d+),(\d+)',self.getStr('CompMask'))
            premask = sequence[0:]
            sequence = rje_sequence.maskLowComplexity(sequence,lowfreq=int(x),winsize=int(y),mask=mask)
            maskx = sequence.count(mask) - premask.count(mask)
            if maskx: masked.append('CompMask (%s,%s): %d %s' % (x,y,maskx,mask))
        elif self.getStr('CompMask').lower()[:1] not in ['f','n','']:
            self.errorLog('CompMask "%s" wrong format for complexity masking. (Switching off.)' % self.getStr('CompMask'),printerror=False)
            self.setStr({'CompMask':'None'})
        ### ~ [6] ~ N-terminal Met ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('MaskM'):
            if sequence[:1] == 'M': sequence = mask + sequence[1:]
            masked.append('MetMask: 1 %s' % mask)
        ### ~ [7] ~ Position-specific masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.dict['MaskPos']:
            premask = sequence[0:]
            for pos in rje.sortKeys(self.dict['MaskPos']):
                try: (x,y) = (int(pos),rje.strList(self.dict['MaskPos'][pos]))
                except: self.errorLog('Problem with MaskPos entry "%s:%s"' % (pos,self.dict['MaskPos'].pop(pos)))
            sequence = rje_sequence.maskPosAA(sequence,self.dict['MaskPos'],mask='X')
            maskx = sequence.count(mask) - premask.count(mask)
            if maskx: masked.append('MaskPos: %d %s' % (maskx,mask))
        ### ~ [8] ~ Motif masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getStr('MotifMask'):
            premask = sequence[0:]
            if 'MotifMask' not in self.obj or not self.obj['MotifMask']:
                self.obj['MotifMask'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
                self.obj['MotifMask'].loadMotifs(self.getStr('MotifMask'))
            maskslims = self.obj['MotifMask']
            for slim in maskslims.slims():
                for hit in slim.searchSequence(sequence=nomask):    # {Pos,Variant,Match,ID,MisMatch}
                    try:
                        r = hit['Pos'] - 1
                        for i in range(len(hit['Variant'])):
                            if hit['Variant'][i] != '.': sequence = rje.strSub(sequence,r+i,r+i,'X')
                    except: self.errorLog('Problem with motifmask %s' % hit)
            maskx = sequence.count(mask) - premask.count(mask)
            if maskx: masked.append('MotifMask (%d SLiMs): %d %s' % (maskslims.motifNum(),maskx,mask))
        ### ~ [9] ~ AA Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.list['AAMask']:
            premask = sequence[0:]
            sequence = rje_sequence.maskAA(sequence,self.list['AAMask'])
            maskx = sequence.count(mask) - premask.count(mask)
            if maskx: masked.append('AAMask (%s): %d %s' % (string.join(self.list['AAMask'],'|'),maskx,mask))
        ### ~ [10] ~ QRegion ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
            self.warnLog('No QRegion masking for single sequences. Switching off.')
            self.list['QRegion'] = []
        ### ~ [11] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if log:
            if masked: self.printLog('#MASK','%s masked: %s.' % (sname, string.join(masked,'; ')),screen=screen)
            else:  self.printLog('#MASK','%s: No positions masked.' % sname,screen=screen)
        return sequence
#########################################################################################################################
    def getSeqObj(self,name,sequence):  ### Returns an rje_sequence.Sequence object for given (name,sequence)
        '''Returns an rje_sequence.Sequence object for given (name,sequence).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            newseq = rje_sequence.Sequence(log=self.log,cmd_list=self.cmd_list,parent=self)
            newseq.setStr({'Name':name,'Type':self.getStr('SeqType')})
            newseq.addSequence(sequence)
            newseq.extractDetails(gnspacc=True)
            return newseq
        except: self.errorLog('SLiMCore.getSeqObj() error'); return None
#########################################################################################################################
    def seqDisorderObj(self,sname,sequence):    ### Returns a disorder object for given (name,sequence).
        '''
        Returns a disorder object for given (sname,sequence).
        '''
        try:### ~ [0] Setup Disorder Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disobj = rje_disorder.Disorder(self.log,self.cmd_list)
            dkey = disobj.getStrLC('Disorder')
            ### ~ [1] Make disorder scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            disobj.disorder(sequence=sequence,name=sname)
            if disobj.list['ResidueDisorder']:  # Should have scores
                if len(disobj.list['ResidueDisorder']) != len(sequence): raise ValueError('%s disorder score length mismatch! (%d score vs %d pos)!' % (sname,len(disobj.list['ResidueDisorder']),len(sequence)))
            else: raise ValueError('Disorder scoring for %s failed!' % sname)
            ## ~ [1a] Make RegionDisorder ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            disobj.makeRegions(sequence)
            return disobj
        except: self.errorLog('%s.seqDisorderObj() failure.' % self.prog()); raise
#########################################################################################################################
    ### <12> ### Extra Functions                                                                                        #
#########################################################################################################################
    def motifSeq(self):     ### Outputs sequence files for given motifs
        '''Outputs sequence files for given motifs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['MotifSeq']: return False  # Dictionary of {pattern:output filename}
            filedict = {}                               # Motif:FileName - Motif made from pattern
            slimlist = self.obj['SlimList']             # Should already contain loaded motifs
            patternlist = rje.sortKeys(self.dict['MotifSeq'])
            for slim in slimlist.slims():
                pattern = slim.pattern()
                if pattern in patternlist:
                    filedict[slim] = self.dict['MotifSeq'][pattern]
                    patternlist.remove(pattern)
            self.printLog('#DEV','%s patterns mapped to SLiMs -> %s remain' % (len(self.dict['MotifSeq']),len(patternlist)))

            for pattern in patternlist:
                slim = slimlist._addMotif(pattern,pattern,reverse=False,check=False,logrem=True)
                if slim: filedict[slim] = self.dict['MotifSeq'][pattern]
                else: self.errorLog('Cannot make SLiM object for "%s"' % pattern,printerror=False)

            if not filedict: return True

            ### Sequences ###
            for slim in slimlist.slims():
                if slim not in filedict: continue
                ## Make SeqList ##
                mseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                if slim.slim() not in self.dict['Slim']:
                    try:
                        for seq in self.seqs():
                            if slim.searchSequence(sequence=seq.info['MaskSeq']): mseq.seq.append(seq)
                    except:
                        self.errorLog('SLiM Error (%s)' % slim,quitchoice=False)
                        continue
                else:
                    for (seq,pos) in self.dict['Slim'][slim.slim()]['Occ']:
                        if seq not in mseq.seq: mseq.seq.append(seq)
                ## Save ##
                OUT = open(filedict[slim],'w')
                for seq in mseq.seq:
                    if self.getBool('Masked'): OUT.write(seq.fasta('PreMask'))
                    else: OUT.write(seq.fasta())
                OUT.close()
                self.printLog('#FAS','%s (unmasked) sequences output to %s' % (rje.iLen(mseq.seq),filedict[slim]))
                del mseq
            return True
        except: self.errorLog('Major problem with %s.motifSeq()' % self.prog())
        return False
#########################################################################################################################
    def randomise(self):    ### Makes random datasets using batch file UPCs
        '''Makes random datasets using batch file UPCs.'''
        try:### ~ [1] ~ Setup Random sequence source ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.checkForFile(self.getStr('RandSource')):
                randsource = rje_seq.SeqList(self.log,['accnr=F','seqnr=F']+self.cmd_list+['seqin=%s' % self.getStr('RandSource'),'autoload=T'])
                if randsource.seqNum() < 1: return self.errorLog('Problem with %s' % self.getStr('RandSource'))
            else:
                randsource = None
                self.printLog('#RAND','No Random Sequence source - will use actual UPC')
            ### ~ [2] ~ Setup UP Lists and make *.upc if missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fullupc = []    # Full list of UPCs (containing Sequence objects!)
            datsize = []    # List of dataset sizes (UPC number) to make
            if self.list['Batch']:
                batchfiles = rje.getFileList(self,filelist=self.list['Batch'],subfolders=False,summary=True,filecount=0)
            if not self.list['Batch'] or not batchfiles:
                self.errorLog('No input file(s) found!',printerror=False)
                return False
            ## ~ [2a] ~ Read *.upc and check/create *.slimdb ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for infile in batchfiles:
                ## Basefile ##
                #self.info['Basefile'] = rje.baseFile(infile)
                self.setStr({'Input':infile})
                if self.getStrLC('ResDir'):
                    rje.mkDir(self,self.getStr('ResDir'))
                ## SeqNum ##
                seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['seqin=%s' % infile,'autoload=T']
                self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
                if self.seqNum() > self.getInt('MaxSeq') > 0:
                    self.printLog('#SKIP','Skipping %s: %s seq > %d MaxSeq' % (self.dataset(),self.seqNum(),self.getInt('MaxSeq')))
                    continue
                ## RandSource/UPC ##
                if randsource:
                    datsize.append(self.seqNum())
                    for i in range(self.seqNum()):
                        r = random.randint(0,randsource.seqNum()-1)
                        fullupc.append((randsource.seq[r],))                        
                else:
                    self.makeUPC()
                    datsize.append(self.UPNum())
                    for upc in self.list['UP']: fullupc.append(upc)
            ### ~ [3] ~ Make Random Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fullupc = rje.randomList(fullupc)
            rje.mkDir(self,self.getStr('RanDir'))
            rbase = self.getStr('RanDir') + self.getStr('Randbase')
            r = 0
            for d in datsize:
                r += 1
                dseq = rje_seq.SeqList(self.log,self.cmd_list+['seqin=None'])
                for u in range(d):
                    upc = fullupc.pop(0)
                    for seq in upc: dseq.seq.append(seq)
                dseq._checkForDup(remdup=True)
                dseq.setStr({'Name': '%s_%s_%d-%d.fas' % (rbase,rje.preZero(r,len(datsize)),dseq.seqNum(),d)})
                dseq.saveFasta()
            ### ~ [4] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('Win32') and len(sys.argv): self.verbose(0,0,'Finished!',1)
        except:
            self.errorLog('Major problem with %s.randomise()' % self.prog())
            raise
#########################################################################################################################
    def BLOSUM2equiv(self,blosumfile,cutoff=1,aalist=[]):   ### Makes equivalence lists from BLOSUM using score cutoff
        '''
        Makes equivalence lists from BLOSUM using score cutoff.
        >> blosumfile:str = File containing BLOSUM matrix
        >> cutoff:int [1] = Minimum BLOSUM score allowed for AA grouping
        >> aalist:list [] = Restrict to listed AAs. Will use default_aa if empty.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not aalist: aalist = rje_slim.default_aas
            blosum = rje_dismatrix.DisMatrix(self.log,self.cmd_list+['Symmetric=T'])
            blosum.loadMatrix(blosumfile,checksym=True,default=0)
            self.printLog('#BLO','Loaded BLOSUM matrix from %s' % blosumfile)
            ### ~ [1] ~ Make lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            equiv = []
            ## ~ [1a] ~ Make BLOSUM neighbours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            neighbours = {}
            for a1 in aalist:
                neighbours[a1] = []
                for a2 in aalist:
                    if a1 == a2: continue
                    if blosum.getDis(a1,a2) >= cutoff: neighbours[a1].append(a2)
            ## ~ [1b] ~ Make equivalence groups ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            growth = aalist[0:]
            while growth:
                grown = []
                for group in growth:
                    stop = True
                    for a1 in aalist:
                        if a1 in group: continue
                        add = True
                        for a2 in group:
                            if a2 not in neighbours[a1]: add = False
                        if add: grown.append(rje.strSort('%s%s' % (group,a1))); stop = False
                    if stop and len(group) > 1: equiv.append(group)
                growth = rje.sortUnique(grown)
            return equiv
        except: self.errorLog('BLOSUM2equiv error!')
#########################################################################################################################
    ### <14> ### REST server output methods. (Should be replaced in mature Classes)                                     #
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT for:

        disorder = List of predicted disorder scores for proteins (if consmask=T or using special disorder rest call)
        rlc = List of RLC scores for proteins (if consmask=T or using special rlc rest call)
        upc = Groupings of unrelated proteins (if efilter=T)

        Note that SLiMCore can either be run as http://rest.slimsuite.unsw.edu.au/slimcore or special runs can be used to
        try and directly access RLC or Disorder scores for an individual protein:

        http://rest.slimsuite.unsw.edu.au/disorder&acc=X
        http://rest.slimsuite.unsw.edu.au/rlc&acc=X&spcode=X

        If these data already exist, they will be returned directly as plain text. If not, a jobid will be returned,
        which will have the desired output once run.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self): return rje.sortKeys(self.dict['Output'])
#########################################################################################################################
### End of SECTION II: SLiMCore Class                                                                                   #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def patternFromCode(slim):  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
    '''Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG).'''
    pattern = ''
    slist = string.split(slim,'-')
    i = 0
    while i < (len(slist)):
        a = slist[i]
        i += 2
        x = '0'
        if i < len(slist): x = slist[i-1]
        if len(a) == 1: pattern += a 
        else: pattern += '[%s]' % a
        if len(x) == 1: pattern += '.' * int(x)
        else: pattern += '.{%d,%d}' % (int(x[0]),int(x[-1]))
    return pattern
#########################################################################################################################
def slimPos(slim):  ### Returns the number of positions in a slim
    return (string.count(slim,'-') / 2) + 1
#########################################################################################################################
def slimLen(slim):  ### Returns length of slim
    return len(patternFromCode(slim))
#########################################################################################################################
def slimDif(slim1,slim2):   ### Returns number of different positins between slim1 and slim2
    '''Returns number of different positions between slim1 and slim2.'''
    if slim1 == slim2: return 0
    (split1,split2) = (string.split(slim1,'-'),string.split(slim2,'-'))
    if len(split1) != len(split2): return rje.modulus(len(split1) - len(split2))
    d = 0
    for i in range(len(split1)):
        if split1[i] != split2[i]: d += 1
    return d
#########################################################################################################################
### END OF SECTION III                                                                                                  #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION IV: MAIN PROGRAM                                                                                            #
#########################################################################################################################
def runMain():
    ### Basic Setup of Program ###
    try: [info,out,mainlog,cmd_list] = setupProgram()
    except SystemExit: return  
    except:
        print 'Unexpected error during program setup:', sys.exc_info()[0]
        return
        
    ### Rest of Functionality... ###
    try: SLiMCore(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: pass    #!#return  # Fork exit etc.
    except KeyboardInterrupt: mainlog.errorLog('User terminated.')
    except: mainlog.errorLog('Fatal error in main %s run.' % info.program)
    mainlog.printLog('#LOG', '%s V:%s End: %s\n' % (info.program, info.version, time.asctime(time.localtime(time.time()))))
#########################################################################################################################
if __name__ == "__main__":      ### Call runMain 
    try: runMain()
    except: print 'Cataclysmic run error:', sys.exc_info()[0]
    sys.exit()
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
