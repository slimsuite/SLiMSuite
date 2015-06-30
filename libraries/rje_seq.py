#!/usr/bin/python

# rje_seq - DNA/Protein sequence module
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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 31 Shanagarry, Milltown Road, Milltown, Dublin 6, Ireland.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      RJE_SEQ
Description:  DNA/Protein sequence list module
Version:      3.23.0
Last Edit:    12/04/15
Copyright (C) 2005  Richard J. Edwards - See source code for GNU License Notice

Function:
    Contains Classes and methods for sets of DNA and protein sequences. 

Sequence Input/Output Options: 
    seqin=FILE      : Loads sequences from FILE (fasta,phylip,aln,uniprot or fastacmd names from fasdb) [None]
    query=X         : Selects query sequence by name [None]
    acclist=LIST    : Extract only AccNums in list. LIST can be FILE or list of AccNums X,Y,.. [None]
    fasdb=FILE      : Fasta format database to extract sequences from [None]
    mapseq=FILE     : Maps sequences from FILE to sequences of same name [None]
    mapdna=FILE     : Map DNA sequences from FILE onto sequences of same name in protein alignment [None]
    seqout=FILE     : Saves 'tidied' sequences to FILE after loading and manipulations [None]
    reformat=X      : Outputs sequence in a particular format, where X is:
                     - fasta/fas/phylip/scanseq/acclist/speclist/acc/idlist/fastacmd/teiresias/mysql/nexus/3rf/6rf/est6rf     [None]
                     - if no seqout=FILE given, will use input file name as base and add appropriate exension.
    #!# reformat=X may not be fully implemented. Report erroneous behaviour! #!#
    logrem=T/F      : Whether to log removed sequences [True] - suggest False with filtering of large files!

Sequence Loading/Formatting Options:
    alphabet=LIST   : Alphabet allowed in sequences [standard 1 letter AA codes]
    replacechar=T/F : Whether to remove numbers and replace characters not found in the given alphabet with 'X' [True]
    autofilter=T/F  : Whether to automatically apply sequence filters etc. upon loading sequence [True]
    autoload=T/F    : Whether to automatically load sequences upon initiating object [True]
    memsaver=T/F    : Minimise memory usage. Input sequences must be fasta. [False]
    degap=T/F       : Degaps each sequence [False]
    tidygap=T/F     : Removes any columns from alignments that are 100% gap [True]
    ntrim=X         : Trims of regions >= X proportion N bases (X residues for protein) [0.0]
    gnspacc=T/F     : Convert sequences into gene_SPECIES__AccNum format wherever possible. [False] 
    seqtype=X       : Force program to read as DNA, RNA, Protein or Mixed (case insensitive; read=Will work it out) [None]
    dna=T/F         : Alternative identification of sequences as DNA [False]
    mixed=T/F       : Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    align=T/F       : Whether the sequences should be aligned. Will align if unaligned. [False]
    rna2dna=T/F     : Converts RNA to DNA [False]
    trunc=X         : Truncates each sequence to the first X aa. (Last X aa if -ve) (Useful for webservers like SingalP.) [0]
    usecase=T/F     : Whether to output sequences in mixed case rather than converting all to upper case [False]
    case=LIST       : List of positions to switch case, starting with first lower case (e.g case=20,-20 will have ends UC) []
    countspec=T/F   : Generate counts of different species and output in log [False]

Sequence Filtering Options:
    filterout=FILE  : Saves filtered sequences (as fasta) into FILE. *NOTE: File is appended if append=T* [None]
    minlen=X        : Minimum length of sequences [0]
    maxlen=X        : Maximum length of sequences (<=0 = No maximum) [0]
    maxgap=X        : Maximum proportion of sequence that may be gaps (<=0 = No maximum) [0]
    maxx=X          : Maximum proportion of sequence that may be Xs (<=0 = No maximum; >=1 = Absolute no.) [0]
    maxglob=X       : Maximum proportion of sequence predicted to be ordered (<=0 = None; >=1 = Absolute) [0]
    minorf=X        : Minimum ORF length for a DNA/EST translation (reformatting only) [0]
    minpoly=X       : Minimum length of poly-A tail for 3rf / 6rf EST translation (reformatting only) [20]
    gapfilter=T/F   : Whether to filter gappy sequences upon loading [True]
    nosplice=T/F    : If nosplice=T, UniProt splice variants will be filtered out [False]
    dblist=LIST     : List of databases in order of preference (good to bad)
                      [sprot,ipi,uniprot,trembl,ens_known,ens_novel,ens_scan]
    dbonly=T/F      : Whether to only allow sequences from listed databases [False]
    unkspec=T/F     : Whether sequences of unknown species are allowed [True]
    accnr=T/F       : Check for redundant Accession Numbers/Names on loading sequences. [True]
    seqnr=T/F       : Make sequences Non-Redundant [False]
    nrid=X          : %Identity cut-off for Non-Redundancy (GABLAMO) [100.0]
    nrsim=X         : %Similarity cut-off for Non-Redundancy (GABLAMO) [None]      
    nralign=T/F     : Use ALIGN for non-redundancy calculations rather than GABLAMO [False]
    specnr=T/F      : Non-Redundancy within same species only [False]
    querynr=T/F     : Perform Non-Redundancy on Query species (True) or limit to non-Query species (False) [True]
    nrkeepann=T/F   : Append annotation of redundant sequences onto NR sequences [False]
    goodX=LIST      : Filters where only sequences meeting the requirement of LIST are kept.
                      LIST may be a list X,Y,..,Z or a FILE which contains a list [None]
                        - goodacc  = list of accession numbers
                        - goodseq  = list of sequence names
                        - goodspec = list of species codes
                        - gooddb   = list of source databases
                        - gooddesc = list of terms that, at least one of which must be in description line
    badX=LIST       : As goodX but excludes rather than retains filtered sequences

System Info Options:
    * Use forward slashes for paths (/)
    blastpath=PATH  : Path to BLAST programs ['']           
    blast+path=PATH : Path to BLAST+ programs ['']           
    fastapath=PATH  : Path to FASTA programs ['']          
    clustalw=PATH   : Path to CLUSTALW program ['clustalw']
    clustalo=PATH   : Path to CLUSTAL Omega alignment program ['clustalo']
    mafft=PATH      : Path to MAFFT alignment program ['mafft']
    muscle=PATH     : Path to MUSCLE alignment program ['muscle']
    fsa=PATH        : Path to FSA alignment program ['fsa']            
    pagan=PATH      : Path to PAGAN alignment program ['pagan']            
    win32=T/F       : Run in Win32 Mode [False]
    alnprog=X       : Choice of alignment program to use (clustalw/clustalo/muscle/mafft/fsa/pagan) [clustalo]
    
Sequence Manipulation/Function Options:
    pamdis      : Makes an all by all PAM distance matrix
    split=X     : Splits file into numbered files of X sequences. (Useful for webservers like TMHMM.)
    relcons=FILE: Returns a file containing Pos AbsCons RelCons [None]
    relconwin=X : Window size for relative conservation scoring [30]
    makepng=T/F : Whether to make RelCons PNG files [False]
    seqname=X   : Output sequence names for PNG files etc. (short/Name/Number/AccNum/ID) [short]

DisMatrix Options:
    outmatrix=X : Type for output matrix - text / mysql / phylip

Special Options:
    blast2fas=FILE1,FILE2,...,FILEn : Will blast sequences against list of databases and compile a fasta file of results per query
        - use options from rje_blast.py for each individual blast (blastd=FILE will be over-ridden)
        - saves results in AccNum.blast.fas and will append existing files!
    keepblast=T/F   : Whether to keep BLAST results files for blast2fas searches [True]
    haqbat=FILE     : Generate a batch file (FILE) to run HAQESAC on generated BLAST files, with seqin as queries [None]

Classes:
    SeqList(rje.RJE_Object):     
        - Sequence List Class. Holds a list of Sequence Objects and has methods for manipulation etc.
    Sequence(rje_sequence.Sequence):     
        - Individual Sequence Class.    
    DisMatrix(rje_dismatrix.DisMatrix):     
        - Sequence Distance Matrix Class.
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                              #
#########################################################################################################################
import copy, glob, math, os, re, shutil, string, sys, time
#########################################################################################################################
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below
import rje
import rje_blast_V1
import rje_blast_V2 as rje_blast
import rje_dismatrix_V2 as rje_dismatrix    #!# Check that this is OK! #!# #import rje_dismatrix_V1 as rje_dismatrix
import rje_pam
import rje_sequence, rje_uniprot, rje_slimcalc
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 0.1 - Renamed major attributes
    # 0.2 - New implementation on more generic OO approach. Non-Redundancy Output
    # 0.3 - No Out Object in Objects
    # 1.0 - Better Documentation to go with GASP V:1.2
    # 1.1 - Better DNA stuff
    # 1.2 - Added ClustalW align
    # 1.3 - Separated Sequence object into rje_sequence.py
    # 1.4 - Add rudimentary gnspacc=T/F
    # 1.5 - Changed pwAln to use popen()
    # 1.6 - Fixed nrdic problem in Redundancy check and added user-definition of database list
    # 1.8 - Added UniProt input and acclist reading
    # 1.9 - Added 'reformat=scanseq' option but not properly implemented. Added align=T/F.
    # 2.0 - Major reworking of commandline options and introduction of self.list dictionary (rje v3.0)
    # 2.1 - Added reformat of UniProt with memsaver=T.
    # 2.2 - Added GABLAM non-redundancy
    # 2.3 - Added NR in memsaver mode
    # 2.4 - Changed some of the log output (REM and redundancy) to look better.
    # 2.5 - Added nr_qry to makeNR()
    # 2.6 - Added mysql reformat output: fastacmd, protein_id, acc_num, spec_code, description (delimited)
    # 2.7 - Added SeqCount() method. Incorporated reading of sequence case.
    # 2.8 - Added NEXUS output for MrBayes compatibility
    # 2.9 - Added setupSubDict(masking=True) for use in probabilistic conservation scores
    # 3.0 - Start of improvements for DNA sequences with dna=T.
    # 3.1 - Added relative conservation calculations for a whole alignment.
    # 3.2 - Added output of sequences for making alignments in R.
    # 3.3 - Added 6RF reformatting for DNA sequences.
    # 3.4 - Added HAQBAT option
    # 3.5 - Added extra alignment program, MAFFT
    # 3.6 - Added stripGap() method. Replaced self.seq with self.seqs() for reading. (Replace with list at some point.)
    # 3.7 - Added raw option for single sequence load.
    # 3.8 - Added maxGlob setting for screening out globular proteins.
    # 3.9 - Added reading of mafft format when not producing standard fasta.
    # 3.10- Added ntrim=X : Trims of regions >= X proportion N bases (X residues for protein) [0.5]
    # 3.11- Added mapdna=FILE option to map DNA sequences onto protein alignment
    # 3.12- Added countspec=T/F   : Generate counts of different species and output in log [False]
    # 3.13- Updated sequence type checking for use with GABLAM 2.10.
    # 3.14- Added CLUSTAL Omega alignment program ['clustalo']
    # 3.15- Added PAGAN alignment program ['pagan'] and (hopefully) fixed minor Windows fastacmd bug.
    # 3.16- Added BLAST+ path and seqFromBlastDBCmd()
    # 3.17- Updated to use BLAST+ and rje_blast_V2
    # 3.18- Minor BLAST+ bug fixes. Added exceptions to readBLAST failure.
    # 3.19- Fixed BLAST+ sequence extraction name truncation error.
    # 3.20- Added run() method for SeqSuite.
    # 3.21.0 - Added extraction of uniprot IDs for seqin.
    # 3.22.0 - Added loading sequences from provided sequence files contents directly, bypassing file reading.
    # 3.22.1 - Fixed problem if seqin is blank, triggering odd Uniprot download.
    # 3.23.0 - Add speclist to reformat options.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [y] Read in sequences from Fasta/Phylip/ClustalW format
    # [y] Output to Fasta Format
    # [y] Will check (for) alignment [checkAln()]
    # [y] Chop sequence up into fragments [winChop()]
    # [x] Interactive menus for module and major activities
    # [ ] Tidy and __doc__ all methods
    # [ ] Profiles
    # [y] Distance Matrices
    # [y] SwissProt Download Format input
    # [y] Phylip Format 'phy' output
    # [ ] Clustalw Align format 'aln' output
    # [ ] MapNames methods - either from order or beginnings of names
    # [y] Redundant Accession number filter
    # [ ] Species from speclist.txt (grep option)
    # [ ] SplitByList method - splits sequences according to those present/absent in another SeqList
    # [ ] SplitByDBase method - like above but divides according to source database (combine in splitBy(split=None)
    # [y] Making of PAM all by all distance matrix
    # [y] Make a split=X command - splits file up into chunks of X [4000 for TMHMM] sequences
    # [ ] Add winchop command and output
    # [y] Add a seqout command - saves results of manipulations to file when called directly from command line
    # [y] Add a extra method to remove 100% gapped columns from alignment 'tidygap'
    # [ ] Allow different cases in sequences?
    # [Y] Remove filterout
    # [Y] Allow MemSaver filtering and sequence manipulation
    # [Y] Move dismatrix into different module
    # [Y] reformat=X option
    # [Y] Add good- and badacc, badspec & badseq instead of filter...
    # [ ] Replace ALIGN with GABLAM for most major applications, such as redundancy filtering etc.
    # [ ] Think about moving methods to an rje_seq_filter methods module?
    # [ ] Incorporate rje_hmm
    # [ ] Make a user menu
    # [ ] Read in alignment etc but extract sequence details from UniProt file?
    # [Y] Extend reformat=X => fasta/phylip/scanseq/acclist/fastacmd/teiresias
    # [Y] Add nextUniProtSeq, like nextFasSeq, for memsaver=T: check format in memSaverReFormat method
    # [ ] Add more details of functionality in docstring
    # [ ] Add masking as in SLiMCore
    # [ ] Tidy this ToDo list!
    # [ ] Add in extra alignment commands as commandline option
    # [ ] Replace self.seq with self.list['Sequence']
    # [ ] Update memsaveNR to use BLAST+. (If not superseded by rje_seqlist.)
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, cyear) = ('RJE_SEQ', '3.23.0', 'April 2015', '2005')
    description = 'RJE Sequence Dataset Manipulation Module'
    author = 'Dr Richard J. Edwards.'
    comments = ['Please report bugs to r.edwards@soton.ac.uk']
    return rje.Info(program,version,last_edit,description,author,time.time(),cyear,comments)
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
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
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
### Regular Expressions for Sequence Details Extraction
re_acconly = re.compile('^(\S+)\s*$')
re_plain = re.compile('^(\S+)\s+(\S.*)$')
re_gnspacc = re.compile('^(\S+)_([A-Z0-9]+)\/(\S+)')
re_gn_sp__acc = re.compile('^(\S+)_([A-Z0-9]+)__(\S+)')
re_uniprot = re.compile('^(sp|tr)\|(\S+)\|(\S+)_(\S+)')
re_uniprot2 = re.compile('^(\S+)_(\S+)\s+\((\S+)\)')
re_unirefprot = re.compile('^([A-Za-z0-9\-]+)\s+([A-Za-z0-9]+)_([A-Za-z0-9]+)\s+')
re_genbank = re.compile('^(\S+)\|(\d+)\|')
re_unigene = re.compile('^\S+\|UG\|(\S+)')
### Standard Alphabets
alph_protx = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X']
alph_dna = ['A','C','G','T','N']
alph_rna = ['A','C','G','U','N']
### Standard Databases (in order of preference, good to bad)
default_db = 'sprot,ipi,uniprot,trembl,ens_known,ens_novel,ens_scan'
### Sequence Extensions ###
seq_ext = {'fasta':'fas','fas':'fas','clustalo':'aln','clustalw':'aln','phylip':'phy','scanseq':'scanseq','teiresias':'teiresias.fas',
           'acclist':'acc','acc':'acc','fastacmd':'fastacmd.txt','uniprot':'dat','idlist':'id','nexus':'nex',
           '6rf':'trans.fas','est6rf':'trans.fas','3rf':'trans.fas','speclist':'spcode'}
#########################################################################################################################
### END OF SECTION I
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: CLASSES                                                                                                     #
#########################################################################################################################

#########################################################################################################################
##  SeqList Class: list of sequences
#########################################################################################################################
### Upgrades to make:
#   - MakeProfile option (HMMER?)
#   - MakeDistanceMarix (and read)
#########################################################################################################################
class SeqList(rje.RJE_Object):     
    '''
    Sequence List Class. Author: Rich Edwards (2005).

    Info:str
    - AccList = Extract only AccNums in list. X can be FILE or list of AccNums X,Y,.. [None]
    - AlnProg = Choice of alignment program to use (clustalw/clustalo/muscle/mafft/fsa) [clustalo]
    - Basefile = Base for filenames (e.g. extension removed)    !New addition so may have quirky functionality!
    - BLAST Path = path to blast programs
    - BLAST+ Path = path to blast+ programs
    - ClustalO = Path to CLUSTAL Omega program ['clustalo']
    - ClustalW = Path to ClustalW
    - DBList = List of databases in order of prefernce X,Y,Z (lower case)
    - Description = Description of sequence list (if desired)
    - FasDB = Fasta format database to extract sequences from [None]
    - FilterOut = Saves filtered sequences (as fasta) into FILE. *NOTE: This file is appended* [None]
    - FASTA Path = path to fasta programs
    - FSA = Path to FSA alignment program ['fsa']            
    - HAQBAT = Generate a batch file (FILE) to run HAQESAC on generated BLAST files, with seqin as queries [None]
    - MapDNA = Map DNA sequences from FILE onto sequences of same name in protein alignment [None]
    - MapSeq = Maps sequences from FILE to sequences of same name [None]
    - MAFFT = path to MAFFT
    - MUSCLE = path to MUSCLE
    - Name = Name of sequence group (Generally filename)
    - PAGAN = path to PAGAN
    - Query = Query sequence name for selection [None]
    - ReFormat = Outputs sequence in a particular format, where X is fasta/phylip/scanseq/acclist/fastacmd/teiresias [None]
    - RelCons = Returns a file containing Pos AbsCons RelCons [None]
    - SeqName = Output sequence names for PNG files etc. (short/Name/Number/AccNum/ID) [short]
    - SeqOut = Saves 'tidied' sequence file after loading and manipulations
    - Type = Type of sequences
    
    Opt:boolean
    - AccNR = Whether to check for redundant Accession Numbers/Names on loading sequences. [False]
    - Align = Whether the sequences should be aligned. Will align if unaligned. [False]
    - Aligned = Whether sequences are aligned (i.e. same length)
    - AutoFilter = Whether to automatically apply sequence filters etc. upon loading sequence [True]
    - AutoLoad = Whether to automatically load sequences upon initiating object [True]
    - CountSpec = Generate counts of different species and output in log [False]
    - DBOnly = Whether to only allow sequences from listed databases
    - Degap = Degaps each sequence. [False]
    - DNA = Alternative identification of sequences as DNA [False]
    - GapFilter = Whether to filter gappy sequences upon loading [True]
    - Gapped = Whether sequences have gaps
    - GeneSpAcc = Converts sequence names into gene_SPEC/Acc format [False]
    - KeepBlast = Whether to keep BLAST results files for blast2fas searches [True]
    - LogRem = Whether to log removed sequences [True] - suggest use of this with filtering of large files!
    - MakePNG = Whether to make RelCons PNG files [False]
    - Mixed = Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    - NoSplice = If nosplice=T, UniProt splice variants will be filtered out [False]
    - NRKeepAnn = Append annotation of redundant sequences onto NR sequences [False]
    - QueryNR = Perform Non-Redundancy on Query species (True) or limit to non-Query species (False) [True]
    - RNAtoDNA = Converts RNA to DNA [False]
    - SeqNR = Make sequences Non-Redundant [False]
    - SpecNR = Non-Redundancy within same species only
    - TidyGap = Removes any columns from alignments that are 100% gap [True]
    - UnkSpec = Whether sequences of unknown species are allowed [True]    
    - NR Align = Use ALIGN for non-redundancy calculations rather than GABLAMO [False]
    - UseCase = Whether to output sequences in mixed case rather than converting all to upper case [False]
    - ReplaceChar = Whether to replace characters not found in the given alphabet with 'X' [True]
    - StripNum = Whether to remove numbers from the input sequence [True]
    
    Stat:numeric 
    - MaxGap = Maximum proportion of gaps in sequences (<=0 = No maximum) [0]
    - MaxGlob = Maximum proportion of sequence predicted to be ordered (<=0 = None; >=1 = Absolute) [0]
    - MaxLen = Maximum length of sequences (<=0 = No maximum) [0]
    - MinLen = Minimum length of sequences [0]
    - MinORF = Minimum ORF length for a DNA/EST translation (reformatting only) [0]
    - MinPoly = Minimum length of poly-A tail for 3rf / 6rf EST translation (reformatting only) [20]
    - MaxX = Maximum proportion of sequence that may be Xs (<=0 = No maximum) [0]
    - NR ID = %Identity cut-off for Non-Redundancy (GABLAMO) [100.0]
    - NR Sim = %Similarity cut-off for Non-Redundancy (GABLAMO) [0.0 (not used)]
    - NTrim = Trims of regions >= X proportion N bases [0.0]
    - RelConWin = Determines window size for RelCons calculation [30]
    - Split = Splits file into numbered files of X sequences. (Useful for webservers like TMHMM.)
    - Trunc = Truncates each sequence to the first X aa. (Last X aa if -ve) (Useful for webservers like SingalP.) [0]

    List:list
    - Alphabet = Alphabet allowed in sequences [standard 1 letter AA codes]
    #!# GoodX and BadX only retained in memory if MemSaver=False. Otherwise, will load, filter, and clear afterwards #!#
    - GoodAcc = list of good accession numbers
    - GoodSeq = list of good sequence names
    - GoodSpec = list of good species codes
    - GoodDB   = list of good source databases
    - GoodDesc = list of terms that, at least one of which must be in description line
    - BadAcc = list of Bad accession numbers
    - BadSeq = list of Bad sequence names
    - BadSpec = list of Bad species codes
    - BadDB   = list of Bad source databases
    - BadDesc = list of terms that, at least one of which must be in description line
    - Case = List of positions to switch case, starting with first lower case (e.g case=20,-20 will have ends UC) []
    - Blast2Fas = FILE1,FILE2,...,FILEn : Will blast sequences against list of databases and compile a fasta file of results per query

    Obj:RJE_Objects
    - QuerySeq = Sequence Object to be used as query
    - PWAln ID = DisMatrix of %ID from Pairwise global alignments (ALIGN)
    - PWAln Sim = DisMatrix of %Sim from Pairwise global alignments (ALIGN)
    - MSA ID = DisMatrix of %ID from Multiple Sequence Alignment
    - MSA Gaps = DisMatrix of %Gaps from Multiple Sequence Alignment
    - MSA Extra = DisMatrix of %Extra aa from Multiple Sequence Alignment
    - PAM Dis = DisMatrix of PAM Distances
    - UniProt = rje_uniprot.UniProt object if UniProt input and memsaver=F

    Other:
    seq:list of Sequence objects
    '''
    ### Old attributes
    ## Info
    # Type : seq_type = None   # Type of sequences (DNA/Protein/Mix etc.)
    ## Opt
    # Aligned & Gapped : aligned = None    # Whether sequences are aligned (True/False/Maybe)
    ## Obj
    # QuerySeq : queryseq = None     # Query Sequence Object
#########################################################################################################################
    ### Attributes
    seq = []          # List of sequences
#########################################################################################################################
    def seqs(self,x=None):
        try: return self.seq[x]
        except: return self.seq
    def seqNum(self): return len(self.seqs())
    def units(self): return {True:'NT',False:'AA'}[self.opt['DNA']]
    def dna(self): return self.opt['DNA']
#########################################################################################################################
    def seqLen(self,seqkey='Sequence'):   ### Returns alignment length - if aligned!
        '''Returns alignment length - if aligned!.'''
        if self.seqNum() == 0: return 0
        if not self.opt['Aligned'] and not self.checkAln(tidygaps=False,seqkey=seqkey):
            self.errorLog('SeqList.seqLen() called but sequences not aligned!',printerror=False)
            raise ValueError
        return self.seqs()[0].seqLen()
#########################################################################################################################
    def aaTotal(self,nonx=False):  ### Returns total number of AA in dataset
        '''Returns total number of AA in dataset.'''
        ax = 0
        for seq in self.seqs():
            if nonx and self.opt['DNA']: ax += seq.nonN()
            elif nonx: ax += seq.nonX()
            else: ax += seq.aaLen()
        return ax
#########################################################################################################################
    def accList(self,seqs=None):
        acclist = []
        if seqs == None: seqs = self.seqs()
        for seq in seqs: acclist.append(seq.info['AccNum'])
        return acclist
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes and adds objects
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### <a> ### Basics
        self.infolist = ['Name','Type','Description','BLAST Path','FASTA Path','MUSCLE','SeqOut','ClustalO','ClustalW',
                         'Basefile','DBList','FasDB','AccList','FilterSpec','Query','MapSeq','ReFormat','FilterOut',
                         'RelCons','SeqName','HAQBAT','AlnProg','MAFFT','FSA','MapDNA','PAGAN','BLAST+ Path']
        self.statlist = ['NR ID','NR Sim','MinLen','MaxLen','MaxGap','RelConWin','Split','Trunc','MaxX','MinORF',
                         'MinPoly','MaxGlob','NTrim']
        self.optlist = ['AutoFilter','AutoLoad','Align','Aligned','GapFilter','Gapped','LogRem','SpecNR','AccNR',
                        'GeneSpAcc','TidyGap','RNAtoDNA','SeqNR','DBOnly','UnkSpec','Degap','NR Align','ReplaceChar',
                        'UseCase','StripNum','DNA','MakePNG','QueryNR','Mixed','NoSplice','NRKeepAnn','Special',
                        'CountSpec','KeepBlast']
        self.objlist = ['QuerySeq','PWAln ID','PWAln Sim','MSA ID','MSA Gaps','MSA Extra','PAM Dis','UniProt']
        self.listlist = ['Alphabet','GoodAcc','GoodSeq','GoodSpec','GoodDB','GoodDesc','BadAcc','BadSeq','BadSpec',
                         'BadDB','BadDesc','Case','Blast2Fas']
        ### <b> ### Defaults
        self._setDefaults(info='None',opt=False,stat=0.0,obj=None,setlist=True)
        if 'win32' in self.cmd_list or 'win32=T' in self.cmd_list:
            self.setInfo({'BLAST Path':rje.makePath('c:/bioware/blast/'),'FASTA Path':rje.makePath('c:/bioware/fasta/'),
                          'MUSCLE':rje.makePath('c:/bioware/muscle.exe',wholepath=True),
                          'ClustalW':rje.makePath('c:/bioware/clustalw.exe',wholepath=True),
                          'PAGAN':rje.makePath('pagan',wholepath=True),
                          'ClustalO':rje.makePath('c:/bioware/clustalo.exe',wholepath=True)})
        else:
            self.setInfo({'BLAST Path':rje.makePath(''),'FASTA Path':rje.makePath(''),'MAFFT':'mafft','FSA':'fsa',
                          'MUSCLE':rje.makePath('/home/bioware/muscle3.52/muscle',wholepath=True),
                          'ClustalW':rje.makePath('clustalw',wholepath=True),
                          'PAGAN':rje.makePath('pagan',wholepath=True),
                          'ClustalO':rje.makePath('clustalo',wholepath=True)})
        self.setOpt({'AccNR':True,'TidyGap':True,'UnkSpec':True,'AutoFilter':True,'AutoLoad':True,'GapFilter':False,
                     'LogRem':True,'ReplaceChar':True,'StripNum':True,'QueryNR':True,'NRKeepAnn':False})
        self.setStat({'NR ID':100.0,'NR Sim':0.0,'RelConWin':30,'MinPoly':20,'NTrim':0.0})
        self.setInfo({'DBList':default_db,'AlnProg':'clustalo','SeqName':'short'})    #,'Type':'Protein'})
        ### <c> ### Other Attributes
        self.seq = []
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        ### Reading of commandline options ###
        for cmd in self.cmd_list:
            try:
                ### General Options ###
                self._generalCmd(cmd)
                
                ### Input/Output Options ###
                self._cmdRead(cmd,type='info',att='Name',arg='seqin')
                self._cmdReadList(cmd,'info',['Query','AccList','MapSeq','ReFormat','AlnProg','SeqName'])
                self._cmdRead(cmd,type='info',att='AlnProg',arg='alignprog')
                self._cmdReadList(cmd,'file',['FasDB','SeqOut','RelCons','HAQBAT','MapDNA'])
                self._cmdReadList(cmd,'opt',['LogRem','MakePNG','CountSpec'])
                self._cmdRead(cmd,type='opt',att='LogRem',arg='remlog')
                
                ### Loading/Formatting Options ###
                self._cmdReadList(cmd,'opt',['AutoFilter','AutoLoad','GapFilter','Align','ReplaceChar','UseCase',
                                             'StripNum','DNA','Degap','Mixed','NoSplice','NRKeepAnn','Special',
                                             'KeepBlast'])
                self._cmdRead(cmd,type='opt',att='GeneSpAcc',arg='gnspacc')
                self._cmdRead(cmd,type='info',att='Type',arg='seqtype')
                self._cmdRead(cmd,type='info',att='Type')
                self._cmdRead(cmd,type='opt',att='RNAtoDNA',arg='rna2dna')
                self._cmdReadList(cmd,'int',['Trunc','RelConWin','MinORF','MinPoly'])
                self._cmdRead(cmd,type='int',att='RelConWin',arg='relwin')
                self._cmdReadList(cmd,'list',['Alphabet','Case'])
                self._cmdReadList(cmd,'stat',['NTrim'])
                self._cmdReadList(cmd,'glist',['Blast2Fas'])

                ### System Info Commands ###
                self._cmdRead(cmd,type='path',att='BLAST Path',arg='blastpath')
                self._cmdRead(cmd,type='path',att='BLAST+ Path',arg='blast+path')
                self._cmdRead(cmd,type='path',att='FASTA Path',arg='fastapath')
                self._cmdReadList(cmd,'file',['MUSCLE','ClustalO','ClustalW','MAFFT','FSA','PAGAN'])
                
                ### Sequence Manipulation/Function Commands ###
                self._cmdRead(cmd,type='int',att='Split')
            except: self.errorLog('Problem with cmd:%s' % cmd)
        self._filterCmd()

        ### Additional processing of commandline options ###
        try:
            if self.stat['NTrim'] > 100.0: self.printLog('#NTRIM','NTrim (%s) cannot exceed 100%%! Switched off.' % self.stat['NTrim']); self.stat['NTrim'] = 0.0
            elif self.stat['NTrim'] >= 1.0: self.stat['NTrim'] /= 100.0
            ## Sequence Type ##
            if self.info['Type'] != 'None':
                self.info['Type'] = self.info['Type'].upper()
                if self.info['Type'] == 'PROTEIN': self.info['Type'] = 'Protein'
            if self.opt['DNA']: self.info['Type'] = 'DNA'
            if len(self.list['Alphabet']) == 1: self.list['Alphabet'] = string.split(string.join(self.list['Alphabet'][0]).upper())
            elif self.list['Alphabet']: self.list['Alphabet'] = string.split(string.join(self.list['Alphabet']).upper())
            #if not self.list['Alphabet']:
            #    if self.info['Type'] == 'DNA': self.list['Alphabet'] = alph_dna + ['-']
            #    elif self.info['Type'] == 'RNA': self.list['Alphabet'] = alph_rna + ['-']
            #    else: self.list['Alphabet'] = alph_protx + ['-']
            ## Non-Redundancy ##                
            if self.opt['SeqNR']: self.opt['AccNR'] = True
            ## AutoLoad ##
            self.autoLoad()
        except:
            self.errorLog('Major problem during rje_seq additional processing of commandline options',quitchoice=True)
#########################################################################################################################
    def getAlphabet(self):  ### Returns appropriate alphabet for use
        '''Returns appropriate alphabet for use.'''
        if self.list['Alphabet']: return self.list['Alphabet']
        if self.info['Type'] == 'DNA': return alph_dna + ['-']
        elif self.info['Type'] == 'RNA': return alph_rna + ['-']
        else: return alph_protx + ['-']        
#########################################################################################################################
    def printLog(self, id='#ERR', text='Log Text Missing!', timeout=True, screen=True, log=True, newline=True, clear=0):   ### Modified from rje
        '''
        Prints text to log with or without run time.
        >> id:str = identifier for type of information
        >> text:str = log text
        >> timeout:boolean = whether to print run time
        >> screen:boolean = whether to print to screen (v>=0)
        >> log:boolean = whether to print to log file
        '''
        try:
            if id in ['\r#REM','#REM'] and not self.opt['LogRem']: return
            try: return self.log.printLog(id, text, timeout, screen and not self.opt['Silent'], log and not self.opt['Silent'], newline)
            except: return self.log.printLog(id, text, timeout, screen, log, newline)
        except: self.errorLog('printLog() problem')
#########################################################################################################################
    def _filterCmd(self,cmd_list=None):   ### Reads filter commands into attributes
        '''
        Reads filter commands into attributes.
        >> cmd_list:list of commands from which to get filter options [None = self.cmd_list]
        '''
        if cmd_list == None: cmd_list = self.cmd_list
        for cmd in cmd_list:
            try:
                self._cmdRead(cmd,type='info',att='FilterOut')
                self._cmdReadList(cmd,'int',['MinLen','MaxLen'])
                self._cmdReadList(cmd,'stat',['MaxGap','MaxX','MaxGlob'])
                self._cmdRead(cmd,type='stat',att='NR ID',arg='nrid')
                self._cmdRead(cmd,type='stat',att='NR Sim',arg='nrsim')
                self._cmdRead(cmd,type='opt',att='NR Align',arg='nralign')
                self._cmdRead(cmd,type='clist',att='DBList')
                self._cmdReadList(cmd,'opt',['LogRem','DBOnly','UnkSpec','AccNR','SeqNR','SpecNR','QueryNR'])
                for filt in ['Acc','Seq','Spec','DB','Desc']:
                    self._cmdRead(cmd,type='list',att='Good%s' % filt)
                    self._cmdRead(cmd,type='list',att='Bad%s' % filt)
            except: self.errorLog('Problem with cmd:%s' % cmd)        
#########################################################################################################################
    def addMatrix(self,idkey,sym=False):    ### Adds distance matrix object if missing
        '''
        Adds distance matrix object if missing.
        >> idkey:str = Key for DisMatrix
        >> sym:boolean = Whether symmetrical
        '''
        try:
            if idkey not in self.obj.keys() or not self.obj[idkey]:
                newmatrix = rje_dismatrix.DisMatrix(log=self.log)
                self.obj[idkey] = copy.deepcopy(newmatrix)
                self.obj[idkey].dict['Matrix'] = {}     # V2
                self.obj[idkey].info['Type'] = idkey
                self.obj[idkey].opt['Symmetric'] = sym
        except:
            self.errorLog('Disaster during addMatrix().')
            raise
#########################################################################################################################
    def makeBaseFile(self): ### Makes self.info['Basefile'] based on self.info['Name']
        '''Makes self.info['Basefile'] based on self.info['Name'].'''
        self.info['Basefile'] = rje.baseFile(self.info['Name'],strip_path=False,extlist=['.fasta','.fa','.fas','.aln','.phy'])
#########################################################################################################################
    def run(self): return self.reFormat()   # Limited function for SeqSuite.
#########################################################################################################################
    ### <2> ### AutoLoad and AutoFilter                                                                                 #
#########################################################################################################################
    def autoLoad(self,autoload=None,autofilter=None,memsaver=None):     ### Controls automatic loading and filtering
        '''
        Controls automatic loading and filtering. Arguments will use self.opt values if None.
        >> autoload:boolean = Whether to load sequences [None]
        >> autofilter:boolean = Whether to filter sequences [None]
        >> memsaver:boolean = Whether to perform filtering using MemSaver option [None]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.opt['MemSaver']: self._filterCmd()
            if self.info['Basefile'] == 'None' and self.getStrLC('Name'): self.makeBaseFile()
            ### ~ [1] ~ AutoLoad ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['Name'].lower() not in ['','none'] and self.opt['AutoLoad'] and not self.opt['MemSaver']:
                self.loadSeqs(seqtype=self.info['Type'])
            ### ~ [2] ~ Alternative AutoFilter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
            if self.info['Name'].lower() not in ['','none'] and self.opt['AutoLoad'] and self.opt['AutoFilter'] and self.opt['MemSaver']:
                self.memSaveAutoFilter()
            #!# Add Reformat #!#
            ### ~ [3] ~ Special CountSpec method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['CountSpec']: self._countSpec()
            ### ~ [4] ~ Special relative conservation method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['RelCons'].lower() not in ['none','t','f','']: self.relCons()
        except: self.errorLog('Pretty bad whoopsie in rje_seq autoLoad.',quitchoice=True)
#########################################################################################################################
    def autoFilter(self,cmd_list=None):    ### AutoFilters sequences based on cmd_list.
        '''
        Performs automatic sequence filtering based on cmd_list.
        >> cmd_list:list of commandline options [self.cmd_list]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cmd_list == None: cmd_list = self.cmd_list
            if cmd_list != self.cmd_list: self._filterCmd(cmd_list)
            ### ~ [1] ~ FilterSeqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
            for filt in ['Acc','Seq','Spec','DB','Desc']:
                if self.list['Bad%s' % filt] and not self.opt['LogRem']: self.printLog('#FILT','Filtering on Bad%s.' % filt)
                if self.list['Good%s' % filt] and not self.opt['LogRem']: self.printLog('#FILT','Filtering on Good%s.' % filt)
            if self.opt['NoSplice']: self.printLog('#FILT','Filtering out UniProt-style splice variants (AccNum-X).')
            self._filterSeqs()
            ### ~ [2] ~ SeqNR ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['SeqNR']: self.makeNR()
            ### ~ [3] ~ SeqOut/SplitSeq/Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['ReFormat'] == 'None' and (int(self.stat['Split']) > 0 or self.info['SeqOut'] != 'None'):
                self.info['ReFormat'] = 'fasta'
            #!# Reformat and splitseq taken care of by module main run #!#
        except: self.errorLog('Minor disaster in rje_seq autoFilter!',quitchoice=True)
#########################################################################################################################
    def memSaveAutoFilter(self,cmd_list=None):    ### AutoFilter sequence using MemSaver mode
        '''
        Performs automatic sequence filtering based on cmd_list.
        >> cmd_list:list of commandline options [self.cmd_list]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Setup Commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if cmd_list == None: cmd_list = self.cmd_list
            self._filterCmd(cmd_list)
            ## ~ [0b] ~ Check file format (Fasta or UniProt) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            re_header = {}
            re_header['fas'] = re.compile('^>\S+')
            re_header['phy'] = re.compile('^\s+(\d+)\s+(\d+)')
            re_header['aln'] = re.compile('^CLUSTAL')
            re_header['uniprot'] = re.compile('^ID')
            filetype = self._readFileType(self.info['Name'],None,re_header)
            if filetype not in ['fas','uniprot']:
                self.errorLog('Cannot Filter sequences of %s format in MemSaver mode. Use memsaver=F.' % filetype,printerror=False)
                raise ValueError
            ## ~ [0c] ~ SetUp MemSaver SeqOut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while self.info['SeqOut'] == 'None' and self.stat['Interactive'] >= 0:
                if rje.yesNo('Overwrite %s following MemSaver filtering?' % self.info['Name'],default='N'):
                    self.info['SeqOut'] = self.info['Name']
                    break
                else:
                    self.info['SeqOut'] = rje.choice('New file name?',default='%s.fas' % self.info['Basefile'])
                    if self.info['SeqOut'] == self.info['Name']: self.info['SeqOut'] = 'None'
            if self.info['SeqOut'] == self.info['Name']: memout = '%s.%s.tmp.fas' % (self.info['Basefile'],rje.randomString(6))
            else: memout = self.info['SeqOut']

            ### ~ [1] ~ Filter Seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for filt in ['Acc','Seq','Spec','DB','Desc']:
                if self.list['Bad%s' % filt]: self.printLog('#FILT','Filtering on Bad%s.' % filt)
                if self.list['Good%s' % filt]: self.printLog('#FILT','Filtering on Good%s.' % filt)
            if self.opt['NoSplice']: self.printLog('#FILT','Filtering out UniProt-style splice variants (AccNum-X).')
            SEQIN = open(self.info['Name'],'r')
            lastline = 'Firstline'
            filter_txt = 'Filtering sequences from %s' % self.info['Name']
            if self.opt['Append']: SEQOUT = open(memout, 'a')
            else: SEQOUT = open(memout, 'w')
            if filetype == 'fas': (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            else: (seq,lastline) = self.nextUniProtSeq(SEQIN,lastline)
            sx = 0  # No. sequences read
            rx = 0  # No. retained
            while seq:
                self._filterSeqs()
                if self.seqNum() > 0:
                    SEQOUT.write('>%s\n%s\n' % (self.seqs()[-1].info['Name'],self.seqs()[-1].getSequence(case=self.opt['UseCase'])))
                    rx += 1
                sx += 1
                self.printLog('\r#FILT','%s: %s seqs, %s retained.' % (filter_txt,rje.integerString(sx),rje.integerString(rx)),log=False,newline=False)
                self.seq = []
                if filetype == 'fas': (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
                else: (seq,lastline) = self.nextUniProtSeq(SEQIN,lastline)
            SEQIN.close()
            SEQOUT.close()
            filter_txt = 'Filtered sequences from %s to %s' % (self.info['Name'],memout)
            self.printLog('\r#FILT','%s: %s seqs, %s retained.' % (filter_txt,rje.integerString(sx),rje.integerString(rx)))
            
            ### ~ [2] ~ Post-filtering cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [2a] ~ Clear FilterCmd Attributes to save memory ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for filt in ['Acc','Seq','Spec','DB','Desc']:
                self.list['Good%s' % filt] = []
                self.list['Bad%s' % filt] = []
            ## ~ [2b] ~ SeqOut ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if memout != self.info['SeqOut']: os.rename(memout,self.info['SeqOut'])
            self.info['Name'] = self.info['SeqOut']
            ## ~ [2c] ~ SeqNR: Not available in MemSaver mode ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['AccNR']: self.printLog('#MEM','AccNR Redundancy filter not available in MemSaver mode!')
            if self.opt['SeqNR']: self.memSaveNR()
            ## ~ [2d] ~ SplitSeq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.info['ReFormat'] == 'None' and int(self.stat['Split']) > 0: self.info['ReFormat'] = 'fasta'
            #!# Reformat and splitseq taken care of by module main run #!#
        except: self.errorLog('Minor disaster in rje_seq memSaveAutoFilter!',quitchoice=True)
#########################################################################################################################
    def _filterSeqs(self):    ### AutoFilters sequences based on self.attributes.
        '''
        Performs automatic sequence filtering based on self.attributes.
        '''
        try:
            ### Filter ###
            for seq in self.seqs()[0:]:
                remseq = None

                ## Bad Filters ##
                if self.list['BadDB'] and seq.info['DBase'] in self.list['BadDB']:
                    remseq = 'BadDB'
                elif self.list['BadSpec'] and (seq.info['SpecCode'] in self.list['BadSpec'] or seq.info['Species'] in self.list['BadSpec']):
                    remseq = 'BadSpec'
                elif self.list['BadAcc'] and seq.info['AccNum'] in self.list['BadAcc']:
                    remseq = 'BadAcc'
                elif self.list['BadSeq'] and (seq.shortName() in self.list['BadSeq'] or seq.info['Name'] in self.list['BadSeq'] or seq.info['ID'] in self.list['BadSeq'] or seq.info['AccNum'] in self.list['BadSeq'] or '%s_%s' % (seq.info['AccNum'],seq.info['SpecCode']) in self.list['BadSeq']):   
                    remseq = 'BadSeq'   #!# Add rem from logs? #!#
                elif self.list['BadDesc']:
                    for desc in self.list['BadDesc']:
                        if seq.info['Name'].find(desc) >= 0:
                            remseq = 'BadDesc'
                            break
                if not remseq and self.opt['NoSplice'] and rje.matchExp('(\S+)-(\d+)',seq.info['AccNum']): remseq = 'Splice'
                if remseq:
                    if self.info['FilterOut'].lower() != 'none':
                        self.saveFasta(seqs=[seq],seqfile=self.info['FilterOut'],append=True,id=False,log=False)                        
                    self.removeSeq(text='Excluded by %s filter' % remseq,seq=seq)
                    continue

                ## Good Filters ##
                remseq = 'Good'
                goodseq = False
                if self.list['GoodDB']:
                    remseq += '-DB'
                    if seq.info['DBase'] in self.list['GoodDB']:
                        goodseq = True
                if self.list['GoodSpec'] and not goodseq:
                    remseq += '-Spec'
                    if seq.info['SpecCode'] in self.list['GoodSpec'] or seq.info['Species'] in self.list['GoodSpec']:
                        goodseq = True
                if self.list['GoodSeq'] and not goodseq:
                    remseq += '-Seq'
                    if seq.shortName() in self.list['GoodSeq'] or seq.info['Name'] in self.list['GoodSeq'] or seq.info['ID'] in self.list['GoodSeq'] or seq.info['AccNum'] in self.list['GoodSeq'] or '%s_%s' % (seq.info['AccNum'],seq.info['SpecCode']) in self.list['GoodSeq']:   
                        goodseq = True
                if self.list['GoodAcc'] and not goodseq:
                    remseq += '-Acc'
                    if seq.info['AccNum'] in self.list['GoodAcc']:
                        goodseq = True
                if self.list['GoodDesc'] and not goodseq:
                    remseq += '-Desc'
                    gooddesc = False
                    for desc in self.list['GoodDesc']:
                        if seq.info['Name'].find(desc) >= 0:
                            gooddesc = True
                    if gooddesc:
                        goodseq = True
                if remseq != 'Good' and not goodseq:
                    if self.info['FilterOut'].lower() != 'none':
                        self.saveFasta(seqs=[seq],seqfile=self.info['FilterOut'],append=True,id=False,log=False)                        
                    self.removeSeq(text='Not found in %s filter' % remseq,seq=seq)
                    continue

                ### Length and Gap Filters ###
                if self.stat['MinLen'] > 0 or self.stat['MaxLen'] > 0:
                    if seq.aaLen() < self.stat['MinLen']:
                        if self.info['FilterOut'].lower() != 'none':
                            self.saveFasta(seqs=[seq],seqfile=self.info['FilterOut'],append=True,id=False,log=False)                        
                        self.removeSeq(text='Too short! (<%daa)' % self.stat['MinLen'],seq=seq)
                        continue
                    elif self.stat['MaxLen'] > 0 and seq.aaLen() > self.stat['MaxLen']:
                        if self.info['FilterOut'].lower() != 'none':
                            self.saveFasta(seqs=[seq],seqfile=self.info['FilterOut'],append=True,id=False,log=False)                        
                        self.removeSeq(text='Too long! (>%daa)' % self.stat['MaxLen'],seq=seq)
                        continue
                if self.opt['GapFilter'] and self.stat['MaxGap'] > 0:
                    self.gapSeqFilter(relative='self')
                self.maxX()
                self.maxGlob()

            ### Additional Sequence Manipulations ###
            if self.opt['RNAtoDNA']:
                rna2dna(seqs=self.seqs())
            if int(self.stat['Trunc']) != 0:
                self.truncSeq(int(self.stat['Trunc']))
            if self.opt['Degap']:
                self.degapSeq()
        except:
            self.errorLog('Minor disaster in rje_seq _filterSeqs!',quitchoice=True)
#########################################################################################################################
    ### <3> ### Sequence Loading Methods                                                                                #
#########################################################################################################################
    def loadSeqs(self,seqfile=None,filetype=None,seqtype=None,aln=None,nodup=None,clearseq=True,seqlines=[]):     ### Loads sequences from file
        '''
        Loads sequences from file.
        >> seqfile:str = file name
        >> filetype:str = format of sequence file
        - 'fas' = fasta, 'phy' = phylip, 'aln' = clustal alignment
        >> seqtype:str = type of sequence in file
        >> aln:Boolean = whether sequences should be aligned
        >> nodup:Boolean = whether to check for (and remove) duplicate sequences.
        >> clearseq:Boolean = whether to clear existing sequences prior to loading [True]
        >> seqlines:list [] = List of lines from sequence "file" to over-ride reading from seqfile.
        '''
        try:### ~ [0] ~ SetUp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if clearseq: self.seq = []
            if aln == None: aln = self.opt['Align']
            if seqtype: 
                if seqtype.lower()[:4] == 'prot': self.info['Type'] = 'Protein'; self.opt['DNA'] = False
                if seqtype.lower()[1:3] == 'na': self.opt['DNA'] = True
            ## ~ [0a] ~ Check File ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not seqlines:
                if not seqfile: seqfile = self.info['Name']
                if seqfile.lower() in ['','none']:
                    if self.stat['Interactive'] >= 0: #?# and len(sys.argv) < 2:
                        seqfile = rje.choice(text='\nInput filename? (Blank to exit.)')
                        if seqfile == '': sys.exit()
                    else:
                        self.errorLog('No file name given: cannot load sequences!',printerror=False,quitchoice=True)
                        return False
                while not os.path.exists(seqfile):
                    if ',' in seqfile or '.' not in seqfile: # Interpret as uniprot extraction list
                        self.obj['UniProt'] = uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
                        uniprot._extractProteinsFromURL(string.split(seqfile,','))
                        sx = 0
                        for uentry in uniprot.entries(): self.seq.append(uentry.obj['Sequence'])
                        if sx:
                            self.printLog('\r#SEQ','%s sequences loaded from Uniprot download (%s accnum).' % (rje.iStr(self.seqNum()),rje.iLen(string.split(seqfile,','))))
                            if self.getStrLC('Basefile'): self.setStr({'Name':self.basefile()})
                            else: self.setStr({'Name':'seqin'})
                            return True
                        else: self.printLog('\r#FAIL','%s not recognised as Uniprot accnum.' % seqfile)
                    if self.i() >= 0: #?# and len(sys.argv) < 2:
                        seqfile = rje.choice(text='Input file "%s" not found. Input filename? (Blank to exit.)' % seqfile)
                        if seqfile == '': sys.exit()
                    else:
                        self.errorLog('File %s not found. Cannot load sequences!' % seqfile,printerror=False,quitchoice=True)
                        return False
                self.info['Name'] = seqfile
            ## ~ [0b] ~ File Type ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            re_header = {}
            re_header['fas'] = re.compile('^>\S+')
            re_header['phy'] = re.compile('^\s+(\d+)\s+(\d+)')
            re_header['aln'] = re.compile('^CLUSTAL')
            re_header['uniprot'] = re.compile('^ID\s+\S')
            re_header['mafft'] = re.compile('^##### atgcfreq')
            if seqlines: filetype = self._readFileType(seqfile,filetype,re_header,seqlines[0])
            else: filetype = self._readFileType(seqfile,filetype,re_header)
            if not filetype: return False
            re_header['fastacmd'] = re.compile('^(\S+)')
            self.makeBaseFile()

            ### ~ [1] ~ Read file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            acclist = rje.listFromCommand(self.info['AccList'])
            if filetype != 'uniprot' and not seqlines:
                seqlines = self.loadFromFile(seqfile,v=2,checkpath=True,chomplines=True)
                if not seqlines: self.errorLog('No lines in file "%s".' % seqfile,printerror=False); return False
                if not re_header[filetype].search(seqlines[0]):
                    self.errorLog("First line of %s is not %s format" % (seqfile,filetype),printerror=False)
                    return False
            elif filetype == 'uniprot' and seqlines:
                self.errorLog('Cannot load Uniprot from seqlines',printerror=False)
                return False

            ### ~ [2] ~ Read sequences according to format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ### 
            p = 0
            if seqfile:
                logseqfile = seqfile
                if len(logseqfile) > 64: logseqfile = '...%s' % seqfile[-64:]
            else: logseqfile = '%d input lines' % len(seqlines)
            ## ~ [2a] ~ Fasta Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if filetype == 'fas':
                sequence = None; name = None
                for line in seqlines:
                    if line.find('>') == 0:   # New Sequence
                        p += 1
                        self.printLog('\r#SEQ','Loading sequences from %s ... %s' % (logseqfile,rje.integerString(p)),log=False,newline=False)
                        if sequence and name: self._addSeq(name,sequence)   # Previous Sequence to Add
                        name = line[1:]
                        sequence = ''
                    else: sequence += line[0:]
                if sequence and name: self._addSeq(name,sequence)   # Previous Sequence to Add
                self.printLog('\r#SEQ','%s sequences loaded from %s (Format: %s).' % (rje.integerString(self.seqNum()),seqfile,filetype))
            ##  ~ [2b] ~ Phylip Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'phy':
                aln = True
                phynum = rje.regExp(re_header['phy'],seqlines[0])
                re_phy = re.compile('^(\S+)\s+(\S+)')
                seqlines.pop(0)
                for line in seqlines:
                    self.printLog('\r#SEQ','Loading sequences from %s ... %s' % (logseqfile,rje.integerString(p)),log=False,newline=False)
                    if re_phy.search(line):   # New Sequence
                        p += 1
                        physeq = rje.regExp(re_phy,line)
                        if len(physeq[1]) != int(phynum[1]):   # Wrong length!
                            raise ValueError, "Phylip Format error: %s is not declared length of %s aa." % (physeq[1],phynum[1])
                        self._addSeq(name=physeq[0], sequence=physeq[1])
                if self.seqNum() != int(phynum[0]):   # Wrong sequence number!
                    raise ValueError, "Phylip Format error: %s seqs declared, %d seqs read." % (physeq[0],self.seqNum())
                self.printLog('\r#SEQ','%s sequences loaded from %s (Format: %s).' % (rje.integerString(self.seqNum()),seqfile,filetype))
            ##  ~ [2c] ~ Aln Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'aln':
                aln = True
                seqlines.pop(0)
                re_aln = re.compile('^(\S+)\s+(\S+)')
                alnseq = {}
                for line in seqlines:
                    self.printLog('\r#SEQ','Loading sequences from %s ... %s' % (logseqfile,rje.integerString(p)),log=False,newline=False)
                    if re_aln.search(line):   # Sequence Line
                        p += 1
                        seq_aln = rje.regExp(re_aln,line)
                        if alnseq.has_key(seq_aln[0]): alnseq[seq_aln[0]] = alnseq[seq_aln[0]] + seq_aln[1]
                        else: alnseq[seq_aln[0]] = seq_aln[1]
                for key in alnseq.keys(): self._addSeq(name=key, sequence=alnseq[key])
                self.printLog('\r#SEQ','%s sequences loaded from %s (Format: %s).' % (rje.integerString(self.seqNum()),seqfile,filetype))
            ## ~ [2d] ~ Fasta Format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'mafft':
                sequence = None; name = None; mheader = True
                for line in seqlines:
                    if line.find('>') == 0:   # New Sequence
                        p += 1; mheader = False
                        self.printLog('\r#SEQ','Loading sequences from %s ... %s' % (logseqfile,rje.integerString(p)),log=False,newline=False)
                        if sequence and name: self._addSeq(name,sequence)   # Previous Sequence to Add
                        name = line[1:]
                        sequence = ''
                    elif mheader: continue
                    elif not line: self._addSeq(name,sequence); sequence = None; name = None
                    else: sequence += line[0:]
                if sequence and name: self._addSeq(name,sequence)   # Previous Sequence to Add
                self.printLog('\r#SEQ','%s sequences loaded from %s (Format: %s).' % (rje.integerString(self.seqNum()),seqfile,filetype))
            ## ~ [2e] ~ UniProt DAT file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'uniprot':
                uniprot = rje_uniprot.UniProt(log=self.log,cmd_list=self.cmd_list)
                uniprot.readUniProt(filename=seqfile,acclist=acclist,logft=False,use_index=False)
                if not self.opt['MemSaver']: self.obj['UniProt'] = uniprot
                for entry in uniprot.list['Entry']:
                    self.seq.append(entry.obj['Sequence'])
                    if not self.opt['MemSaver']: self.seqs()[-1].obj['Entry'] = entry
                    p += 1
                self.printLog('\r#SEQ','%s sequences loaded from %s (Format: %s).' % (rje.integerString(self.seqNum()),seqfile,filetype))
            ## ~ [2f] ~ UniProt DAT file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif filetype == 'fastacmd':
                ## Check Database ##
                if not os.path.exists(self.info['FasDB']):
                    self.errorLog('Fasta database file %s missing! Cannot load sequences.' % self.info['FasDB'],printerror=False)
                    raise IOError
                protdb = True
                if self.info['Type'].lower() in ['rna','dna']: protdb = False
                if not rje_blast.checkForDB(dbfile=self.info['FasDB'],log=self.log,protein=protdb,oldblast=True):
                    rje_blast.formatDB(self.info['FasDB'],self.info['BLAST Path'],protein=protdb,log=self.log,oldblast=True)
                ## Try to extract sequences ##
                failures = 0
                for line in seqlines:
                    fastacmd = rje.regExp(re_header['fastacmd'],line)
                    if fastacmd:
                        if self.seqFromFastaCmd(fastacmd[0],self.info['FasDB']): p += 1
                        else:
                            self.errorLog('Could not extract %s from %s using fastacmd. Check win32=T/F!' % (fastacmd[0],self.info['FasDB']),printerror=False)
                            failures += 1
                        self.printLog('\r#SEQ','Loading sequences from %s ... %s' % (logseqfile,rje.integerString(p)),log=False,newline=False)
                self.printLog('\r#SEQ','%s lines read from %s. %s sequences extracted from %s using fastacmd. %s failures.' %
                                  (rje.integerString(len(seqlines)),self.info['Name'],rje.integerString(self.seqNum()),self.info['FasDB'],rje.integerString(failures)))

            ### ~ [3] ~ AccList Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if acclist:
                ## ~ [3a] ~ Perform Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                missing_list = acclist[0:]
                for nextseq in self.seqs()[0:]:
                    seqfound = False
                    for acc in acclist[0:]:
                        if nextseq.hasID(acc):
                            seqfound = True
                            if acc in missing_list: missing_list.remove(acc)
                    if not seqfound and filetype != 'uniprot': self.seq.remove(nextseq) # Not in list
                if filetype != 'uniprot': self.printLog('#SEQ','%d kept after acclist check.' % (self.seqNum()))
                self.printLog('#SEQ','%d missing accnum not found in %s.' % (len(missing_list),self.info['Name']))
                ## ~ [3b] ~ Report missing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for missing in missing_list:
                    if filetype == 'uniprot': self.printLog('#ERR','%s not found in %s: check for UniProt secondary AccNum.' % (missing,self.info['Name']))
                    else: self.printLog('#ERR','%s not found in %s!' % (missing,self.info['Name']))
                if self.seqNum() == 0: self.errorLog('No sequences retained!',printerror=False); return False

            ### ~ [4] ~ Sequence Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqtype and seqtype.lower() not in ['none','read']: self.info['Type'] = seqtype
            elif seqtype and seqtype.lower() == 'read': self.readSeqType()  # Try to work out type

            ### ~ [5] ~ AccNR & CheckAln ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.querySeq()     # Must make query before checking duplicates
            if nodup == None: nodup = self.opt['AccNR']
            if nodup: self._checkForDup(True)
            if self.stat['Verbose'] >= 2:
                for s in range(self.seqNum()):
                    self.verbose(2,4,"%d: %s (Len=%s)" % ((s+1),self.seqs()[s].info['Name'],rje.integerString(self.seqs()[s].aaLen())),1)
                    self.verbose(3,3,self.seqs()[s].details(),1)
                self.verbose(2,4,self.details(),1)

            ### ~ [6] ~ Map Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['MapSeq'] != 'None':
                tmpseq = SeqList(cmd_list=['v=-1','i=-1','autofilter=F'])
                tmpseq.loadSeqs(seqfile=self.info['MapSeq'])
                self.mapSeq(tmpseq)    
                self.querySeq()

            ### ~ [7] ~ Filters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['AutoFilter']: self.autoFilter()
            self._checkAln(aln=aln,realign=True)

            ### ~ [8] ~ Map Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.info['MapDNA'].lower() not in ['','none']:
                tmpseq = SeqList(cmd_list=['v=-1','i=-1','autofilter=F'])
                tmpseq.loadSeqs(seqfile=self.info['MapDNA'])
                self.mapSeq(tmpseq,mapaln=True)    
                self.querySeq()

            ### ~ [9] ~ Sequence Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            _stage = 'Sequence Details'
            if self.stat['Verbose'] >= 2:
                for s in range(self.seqNum()):
                    self.verbose(2,4,"%d: %s" % ((s+1),self.seqs()[s].info['Name']),1)
                    self.verbose(3,3,self.seqs()[s].details(),1)
                self.verbose(2,4,self.details(),1)

            return True
        except: self.errorLog('Error in "%s" loadSeqs' % self.info['Name']); raise
#########################################################################################################################
    def readSeqType(self,log=True): ### Calculates sequence format from sequences
        '''Calculates sequence format from sequences.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prot = dna = rna = 0
            if self.opt['Mixed']: tseq = self.seqs()
            else: tseq = self.seqs()[:1]
            ### ~ [1] ~ Read Sequence Types ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in tseq:
                if seq.seqType() == 'DNA': dna = 1
                elif seq.seqType() == 'Protein': prot = 2
                elif seq.seqType() == 'RNA': rna = 4
            total = dna + rna + prot
            if total == 1: self.info['Type'] = 'DNA'
            elif total == 2: self.info['Type'] = 'Protein'
            elif total == 4: self.info['Type'] = 'RNA'
            elif total == 5: self.info['Type'] = 'NA'
            else: self.info['Type'] = 'Mixed'
            if log: self.printLog('#TYPE','Sequence type: %s' % self.info['Type'])
            return self.info['Type']
        except: self.errorLog('Error in "%s" readSeqType' % self.info['Name']); raise
#########################################################################################################################
    def querySeq(self,query=None):     ### Sets Query Sequence if appropriate
        '''
        Sets Query Sequence if appropriate.
        >> qrycmd:str [None] = Name to be used for query search (equiv. cmd query=X)
        << qry if query selected. None if not.
        '''
        old_qry = self.obj['QuerySeq']
        self.obj['QuerySeq'] = None
        if not query:
            for cmd in self.cmd_list:
                if cmd.find('query=') == 0: query = cmd[len('query='):]     # Command found
        if query:
            rawq = query
            for seq in self.seqs():
                prenum = rje.matchExp('^(\d+_)(.+)$',seq.info['ID'])
                if prenum and prenum[1] == query: query = prenum[0] + query
                if query in [seq.info['ID'],seq.info['AccNum'],seq.shortName()]:
                    if self.obj['QuerySeq'] != seq and old_qry != seq:     # Not Already query
                        self.printLog('#QRY','Query Sequence = %s (%s %s).' % (seq.shortName(),seq.aaNum(),self.units()))
                    self.obj['QuerySeq'] = seq
                    self.seq.remove(seq)
                    self.seq = [seq] + self.seqs()
            if not self.obj['QuerySeq'] and rje.matchExp('^(\d+)$',rawq):
                qx = string.atoi(rje.matchExp('^(\d+)$',rawq)[0]) - 1
                try: seq = self.seqs()[qx]
                except: seq = None
                if seq:
                    if old_qry != seq:     # Not Already query
                        self.printLog('#QRY','Query Sequence = Seq %d = %s (%s %s).' % (qx+1,seq.shortName(),seq.aaNum(),self.units()))
                    self.obj['QuerySeq'] = seq
                    self.seq.remove(seq)
                    self.seq = [seq] + self.seqs()
                else: self.printLog('#QRY','Cannot find query Sequence "%s".' % (rawq))
        else: self.obj['QuerySeq'] = old_qry
        return self.obj['QuerySeq']
#########################################################################################################################
    def _readFileType(self,seqfile,filetype,re_header,firstline=None):       ### Determines format of file
        '''
        Determines format of input file for loading sequences.
        >> seqfile:str = file name
        >> filetype:str = format of sequence file
        - 'fas' = fasta, 'phy' = phylip, 'aln' = clustal alignment, 'uniprot' = UniProt, 'fastacmd' = Fastacommand names
        >> re_header:dict = Regex recognition for filetypes
        >> firstline:str [None] = Optional file first line to over-ride use of seqfile.
        '''
        ### <a> ### Look for Recognised format given       
        for ftype in re_header.keys():
            if filetype == ftype:
                return filetype
        ### <b> ### Discern from file   
        if not firstline:
            try:
                SEQFILE = open(seqfile, 'r')
                firstline = SEQFILE.readline()
                SEQFILE.close()
            except(IOError):
                self.errorLog("File %s not found" % seqfile)
                raise
        for ftype in re_header.keys():
            if re_header[ftype].search(firstline):
                return ftype
        if self.info['FasDB'].lower() != 'none':    ### Fastacmd extraction
            return 'fastacmd'
        self.errorLog('%s not a recognied format (Fasta/Phylip/Clustal) and no FasDB given for fastacmd extraction.' % seqfile,printerror=False)
        return None
#########################################################################################################################
    def _addSeq(self,name,sequence):    ### Adds a new Sequence Object to list
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
            #self.debug(newseq.info)
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
    def nextFasSeq(self,fileobject=None,lastline=None,raw=False):  ### Returns sequence object of next sequence from passed File Object
        '''
        Returns sequence object of next sequence from passed File Object. Returns None if end of file.
        If lastline=None, file MUST be FASTA format and ONE LINE PER SEQUENCE, else Fasta only.
        If lastline given, will return a tuple of (sequence object, nextline)
        Sequence object is also placed in self.seq.
        >> fileobject: File Object from which sequence to be read
        >> lastline:str = last line read from FILE object, typically the next description line
        >> raw:bool = Whether to return (name,sequence) instead of sequence object
        << sequence object or tuple of (sequence object, nextline)
        '''
        try:### ~ [1] ~ Without lastline: Fasta files, one line per sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if lastline == None:
                line = fileobject.readline().strip('\r\n')
                if not line: return None
                if line.find('>') != 0:   # Wrong format
                    self.errorLog("Format Problem! '>' Expected!",printerror=False)
                    raise ValueError
                name = line[1:]
                sequence = fileobject.readline().strip('\r\n')
                if raw: return (name,sequence)
                self._addSeq(name=name, sequence=sequence)
                return self.seqs()[-1]
            ### ~ [2] ~ With lastline: fasta file but sequence can be multi-line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                line = lastline
                while line.find('>') != 0:  # Find next description line
                    line = fileobject.readline()
                    if not line: return (None,'')
                name = line.strip('\r\n')[1:]
                sequence = fileobject.readline().strip('\r\n')
                line = fileobject.readline()
                while line.find('>') != 0:
                    sequence = '%s%s' % (sequence,line.strip('\r\n'))
                    line = fileobject.readline()
                    if not line: break
                if raw: return ((name,sequence),line)
                self._addSeq(name=name, sequence=sequence)
                return (self.seqs()[-1],line)
        except:
            self.errorLog('Major error during nextFasSeq()')
            raise
#########################################################################################################################
    def nextUniProtSeq(self,fileobject=None,lastline=None):  ### Returns sequence object of next sequence from passed File Object
        '''
        Returns sequence object of next sequence from passed File Object. Returns None if end of file.
        Will return a tuple of (sequence object, nextline)
        Sequence object is also placed in self.seq.
        >> fileobject: File Object from which sequence to be read
        >> lastline:str = last line read from FILE object, typically the next description line
        << tuple of (sequence object, nextline)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            uniprot = rje_uniprot.UniProt(log=self.log,cmd_list=self.cmd_list)
            line = lastline
            if line == '': return (None,'')
            ## ~ [0a] ~ Find next ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            while line.find('ID') != 0:
                line = fileobject.readline()
                if not line: return (None,'')
            ### ~ [1] ~ Read entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            unilines = [line]
            while line.find('//') != 0:
                line = fileobject.readline()
                if not line:
                    retline = ''
                    break
                unilines.append(line)
                retline = line
            ### ~ [2] ~ Process Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###                
            i = -1
            _ex = 0
            _entry = None
            _reading = False
            while i < (len(unilines)-1):
                ## ~ [2a] ~ New Line ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                i += 1
                line = string.strip(unilines[i],'\r\n')
                type = line[0:2]
                if type == '': continue
                rest = line[5:]
                ## ~ [2b] ~ New Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if type == 'ID':
                    _reading = True
                    _entry = rje_uniprot.UniProtEntry(log=self.log,cmd_list=self.cmd_list)
                ## ~ [1c] ~ End of Entry ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if type == '//':
                    uniprot.list['Entry'].append(_entry)
                    _entry.process(logft=False)    # Extracts details from uniprot
                    _ex += 1
                    _reading = False
                    continue
                ## ~ [2d] ~ Entry of Details ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if _entry.dict['Data'].has_key(type):   # Append list
                    if rest[0] != ' ': _entry.dict['Data'][type].append(rest)   # New entry
                    else: _entry.dict['Data'][type][-1] = '%s %s' % (_entry.dict['Data'][type][-1], rest)
                elif type == '  ':
                    if _entry.dict['Data'].has_key('SEQ'): _entry.dict['Data']['SEQ'][0] = '%s %s' % (_entry.dict['Data']['SEQ'][0],rest)
                    elif _entry.dict['Data'].has_key('SQ'): _entry.dict['Data']['SEQ'] = [rest]
                else: _entry.dict['Data'][type] = [rest]
            ### ~ [3] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.seq.append(uniprot.list['Entry'][0].obj['Sequence'])
            return (self.seqs()[-1],retline)
        except:
            self.errorLog('Major error during nextUniProtSeq()')
            raise
#########################################################################################################################
    def seqFromFastaCmd(self,id,dbase=None):  ### Returns sequence object of sequence as obtained with fastacmd
        '''
        Returns sequence object of sequence as obtained with fastacmd.
        >> id:str = id of sequence to pass to fastacmd
        >> dbase:str = formatted database
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not dbase: dbase = self.info['Name']
            blastpath = rje.makePath(self.info['BLAST Path']) + 'fastacmd'
            if self.opt['Win32']: FASTACMD = os.popen("%s -s \"%s\" -d %s" % (blastpath,id,dbase))
            else: FASTACMD = os.popen("%s -s '%s' -d %s" % (blastpath,id,dbase))
            ## ~ [1a] Check success ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            line = FASTACMD.readline().strip('\r\n')
            if not line or line.find('>') != 0:
                FASTACMD.close()
                self.errorLog('No sequence returned by "%s -s %s -d %s"' % (blastpath,id,dbase),printerror=False)
                #self.deBug('BLAST Path = %s' % self.info['BLAST Path'])
                #self.deBug('CmdList = %s' % self.cmd_list)
                return None
            ### ~ [2] Extract details and add sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            if line.startswith('>lcl|'): name = line[5:]
            elif line.startswith('>'): name = line[1:]
            else: name = line[0:]
            sequence = ''
            while 1:
                line = FASTACMD.readline().strip('\r\n')
                if not line: break
                sequence += line[0:]
            FASTACMD.close()
            self._addSeq(name=name, sequence=sequence)
            if self.seqNum() > 0: return self.seqs()[-1]
            else:
                self.errorLog('No sequences in %s?! (Tried to add %s, %d aa).' % (self.info['Name'],name,len(sequence)),printerror=False)
                return None
        except:
            self.errorLog('Major error during seqFromFastaCmd()')
            raise
#########################################################################################################################
    def seqFromBlastDBCmd(self,id,dbase=None,expect=True):  ### Returns sequence object(s) of sequence as obtained with blastdbcmd
        '''
        Returns sequence object of sequence as obtained with blastdbcmd.
        >> id:str = id (or list of ids) of sequence to pass to blastdbcmd
        >> dbase:str = formatted database
        >> expect:bool [True] = whether sequences are expected to be returned.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cmdseq = []
            try: id = id.split(); singleseq = True
            except: singleseq = False    # Already a list of IDs
            id = string.join(id,',')
            if not dbase: dbase = self.info['Name']
            blastpath = rje.makePath(self.info['BLAST+ Path']) + 'blastdbcmd'
            if self.opt['Win32']: BLASTDBCMD = os.popen("%s -entry \"%s\" -db %s" % (blastpath,id,dbase))
            else: BLASTDBCMD = os.popen('%s -entry "%s" -db %s' % (blastpath,id,dbase))
            ### ~ [2] Extract details and add sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###            
            lastline = 'BLASTDBCMD'; sx = 0
            (blastseq,lastline) = self.nextFasSeq(BLASTDBCMD,lastline,raw=True)
            while blastseq:
                (name,sequence) = blastseq; sx += 1
                if name.startswith('lcl|'): name = name[4:]
                cmdseq.append(self._addSeq(name=name, sequence=sequence))
                (blastseq,lastline) = self.nextFasSeq(BLASTDBCMD,lastline,raw=True)
            BLASTDBCMD.close()
            ### ~ [3] Return sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if cmdseq and singleseq: return cmdseq[0]
            elif cmdseq: return cmdseq
            else:
                if expect: self.errorLog('No sequences extracted from %s. (%d returned by blastdbcmd)' % (dbase,sx),printerror=False)
                return []
        except:
            self.errorLog('Major error during seqFromBlastDBCmd()')
            raise
#########################################################################################################################
    def mapSeq(self,seqlist,mapaln=False):  ### Maps sequences (info['Sequence']) from seqlist onto self Sequence Objects
        '''
        Maps sequences (info['Sequence']) from seqlist onto self Sequence Objects.
        - looks for matching shortName() or matching AccNum
        >> seqlist:SeqList Object that has sequences to map
        >> mapaln:bool [False] = Map the new sequences onto an existing alignment.
        << returns True if all sequences mapped, else False
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            map = {}
            mapx = 0
            mapped = {}
            for seq in self.seqs(): mapped[seq] = False
            ### ~ [1] ~ Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for mapseq in seqlist.seq:
                map[mapseq] = True      # Sequence needs to be mapped
                for seq in self.seqs():
                    if mapped[seq]: continue
                    if seq.shortName() == mapseq.shortName():
                        if mapaln: seq.info['Unmapped'] = seq.info['Sequence']
                        seq.info['Sequence'] = mapseq.info['Sequence']
                        mapped[seq] = True
                        map[mapseq] = False
                        mapx += 1
                        break
                if map[mapseq]: # Try AccNum
                    for seq in self.seqs():
                        if mapped[seq]: continue
                        if seq.info['AccNum'] == mapseq.info['AccNum']:
                            if mapaln: seq.info['Unmapped'] = seq.info['Sequence']
                            seq.info['Sequence'] = mapseq.info['Sequence']
                            mapped[seq] = True
                            map[mapseq] = False
                            mapx += 1
                            break
            for mapseq in seqlist.seq:
                if map[mapseq]: # Try partial match
                    for seq in self.seqs():
                        if mapped[seq]: continue
                        if seq.shortName().find(mapseq.shortName()) == 0:
                            if mapaln: seq.info['Unmapped'] = seq.info['Sequence']
                            seq.info['Sequence'] = mapseq.info['Sequence']
                            mapped[seq] = True
                            map[mapseq] = False
                            mapx += 1
                            break
            ### ~ [2] ~ Check Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            allmapped = False
            if seqlist.seqNum() == self.seqNum():   # Assume 1 to 1 mapping
                if mapx == self.seqNum(): allmapped = True
                else:
                    print 'No match for:'
                    for mapseq in seqlist.seq:
                        if map[mapseq]:  print mapseq.shortName(),
                    print '\n\nNo match to:'
                    for seq in self.seqs():
                        if mapped[seq] == False: print seq.shortName(),
                    print '\n'                            
                    self.errorLog('Only %d of %d sequences mapped from %s.' % (mapx,seqlist.seqNum(),seqlist.info['Name']),printerror=False,quitchoice=True)
            else:
                if mapx == self.seqNum(): allmapped = True
                else: self.printLog('#SEQ','Only %d of %d sequences mapped from %s.' % (mapx,self.seqNum(),seqlist.info['Name']))
            if not (allmapped and mapaln): return allmapped
            ### ~ [3] ~ Attempt to map sequence onto alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                seqlen = string.replace(seq.info['Unmapped'],'-','')
                seq.info['Sequence'] = string.replace(seq.info['Sequence'],'-','')
                maplen = seq.aaLen()
                if seqlen == maplen:    # Simple 1to1 mapping
                    newseq = ''; i = 0
                    for r in seq.info['Unmapped']:
                        if r == '-': newseq += '-'
                        else: newseq += seq.info['Sequence'][i]; i += 1
                    if i == seq.aaLen(): seq.info['Sequence'] = newseq
                    else: self.errorLog('Wrong number of %s in %s' % (seq.unit(),seq.shortName()),printerror=False); allmapped = False
                elif seqlen == (maplen * 3):    # Map protein onto DNA
                    newseq = ''; i = 0; n = 0
                    for r in seq.info['Unmapped']:
                        if r == '-': newseq += '-'
                        elif n == 1: n += 1
                        elif n == 2: n = 0
                        else: newseq += seq.info['Sequence'][i]; i += 1; n += 1
                    if i == seq.aaLen(): seq.info['Sequence'] = newseq
                    else: self.errorLog('Wrong number of %s in %s' % (seq.unit(),seq.shortName()),printerror=False); allmapped = False
                elif (seqlen * 3) == maplen:        # Map DNA onto protein
                    newseq = ''; i = 0
                    for r in seq.info['Unmapped']:
                        if r == '-': newseq += '-'
                        else: newseq += seq.info['Sequence'][i:i+3]; i += 3
                    if i == seq.aaLen(): seq.info['Sequence'] = newseq
                    else: self.errorLog('Wrong number of %s in %s' % (seq.unit(),seq.shortName()),printerror=False); allmapped = False
                else: self.errorLog('Wrong number of %s in %s' % (seq.unit(),seq.shortName()),printerror=False); allmapped = False
            return allmapped
        except: self.errorLog('Major error during mapSeq()'); raise
#########################################################################################################################
    def _countSpec(self):   ### Counts sequences for each species and outputs numbers in log file
        '''Counts sequences for each species and outputs numbers in log file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['MemSaver']: return self.printLog('#SPEC','Cannot count species in MemSaver mode')
            specx = {}
            ### ~ [1] ~ Count species ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                spec = seq.info['SpecCode']
                if spec in specx: specx[spec] += 1
                else: specx[spec] = 1
            ### ~ [2] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for spec in rje.sortKeys(specx): self.printLog('#SPEC','%s: %s sequences' % (spec,rje.iStr(specx[spec])))                               
        except: self.errorLog('Major error during countSpec()'); raise
#########################################################################################################################
    ### <4> ### Sequence Output                                                                                         #
#########################################################################################################################
    def saveFasta(self,seqs=[],seqfile=None,linelen=0,name='Name',namelen=0,append=False,id=False,log=True,case=None,screen=None):  ### Saves sequences in fasta format
        '''
        Saves sequences in SeqList object in fasta format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> linelen:int [0] = max seqline length [0 = all on one line]
        >> name:str = Type of name to use as sequence name: 'short'=shortName(), 'AccNum'=AccNum,
                    'Teiresias'=Teiresias format, 'Number' = Number only
        >> namelen:int [0] = max length of sequence name [0 = no max]
        >> append:boolean [False] = append, do not overwrite, file
        >> id:boolean [False] = Appends sequence number to start of name.
        >> log:boolean [True] = Whether to log output
        >> case:boolen [False] = Whether to use self.dict['Case'] to set output case
        >> screen:bool [None] = Whether to print log output to screen (None will use log setting)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqs: seqs = self.seqs()
            if not seqfile: seqfile = self.info['Name']
            if case == None: case = self.opt['UseCase']
            if screen == None: screen = log
            outlist = []    # Store lines to be printed to file
            if name.lower() == 'teiresias': linelen = 0; id = False
            ### ~ [1] ~ Build output list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in seqs:
                ## ~ [1a] ~ Sequence name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if name.lower() == 'short': outname = seq.shortName()
                elif name.lower() == 'teiresias': outname = '%s 1' % seq.shortName()
                elif name.lower() in ['num','number']: outname = '%d' % (self.seqs().index(seq)+1); id = False
                else: outname = seq.info[name]
                if id: outname = '%d %s' % (self.seqs().index(seq)+1,outname)
                if (namelen > 0) & (len(outname) >= namelen): outname = outname[0:(namelen-3)] + '...'
                outlist.append('>%s\n' % outname)
                ## ~ [1b] ~ Sequence data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                sequence = seq.getSequence(case)
                if (linelen > 0) and (seq.seqLen() > linelen):
                    r = linelen
                    while r < seq.seqLen():
                        outlist.append('%s\n' % sequence[(r-linelen):r])
                        r += linelen
                    if seq.seq[(r-linelen):]: outlist.append('%s\n' % sequence[(r-linelen):])
                else: outlist.append('%s\n' % sequence)
            ### ~ [2] ~ Open file & write lines ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if append: open(seqfile, 'a').writelines(outlist)
            else: open(seqfile, 'w').writelines(outlist)
            self.printLog('\r#FAS',"%d Sequences output to %s.     " % (len(seqs),seqfile),log=log,screen=screen)
        except(IOError): self.errorLog("Cannot create %s" % seqfile); raise
        except: self.errorLog("Problem saving sequences to %s" % seqfile); raise
#########################################################################################################################
    def fasta(self,seqs=[]):    ### Returns text of sequences in fasta format
        '''Returns text of sequences in fasta format.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqs: seqs = self.seqs()
            fastxt = ''
            for seq in seqs: fastxt += '>%s\n%s\n' % (seq.info['Name'],seq.getSequence(case=self.getBool('UseCase')))
            return fastxt
        except: self.errorLog("Problem with fasta()"); raise
#########################################################################################################################
    def savePhylip(self,seqs=[],seqfile=None,name='num',id=False,log=True):  ### Saves sequences in phylip format
        '''
        Saves sequences in SeqList object in fasta format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].phy] = filename
        >> name:str = Type of name to use as sequence name: 'short'=shortName(), 'AccNum'=AccNum, 'num'=Number
            - Note that if any names are >10 characters long, Numbers will be used instead
        >> id:boolean [False] = Appends sequence number to start of name.
        >> log:boolean [True] = Whether to log output
        '''
        ### <0> ### Setup
        if seqs == []:
            seqs = self.seqs()
        baselen = seqs[0].seqLen()
        for seq in seqs:
            if seq.seqLen() != baselen:
                self.errorLog('Phylip output selected but sequences are not equal lengths.',True,False)
                raise ValueError

        ### <1> ### Sort output names
        output_names = {}   # seq:name
        numbered = False
        ## Numbered or names? ##
        for seq in seqs:
            if name == 'short':
                outname = seq.shortName()
            elif name == 'AccNum':
                outname = seq.info[name]
            else:
                numbered = True
                break
            if id:
                outname = '%d_%s' % (self.seqs().index(seq)+1,outname)
            if len(outname) > 10:
                numbered = True
                break
            output_names[seq] = outname
        ## Replace with numbered if appropriate ##
        if numbered:
            for seq in seqs:
                output_names[seq] = '%d' % (self.seqs().index(seq)+1)
        ## Make correct length ##
        for seq in seqs:
            while len(output_names[seq]) < 10:
                output_names[seq] += ' '

        ### <2> ### List of lines to output
        outlist = ['  %d  %d\n' % (len(seqs),seqs[0].seqLen())]
        for seq in seqs:
            outlist.append('%s%s\n' % (output_names[seq],seq.info['Sequence']))

        ## <b> ## Open file & write lines
        if seqfile == None:
            seqfile = '%s.phy' % rje.baseFile(self.info['Name'],True)
        try:
            SEQOUT = open(seqfile, 'w')
            SEQOUT.writelines(outlist)
            SEQOUT.close()
        except(IOError):
            self.errorLog("Cannot create %s" % seqfile)
            raise
        except:
            self.errorLog("Problem creating %s" % seqfile)
            raise

        if log:
            self.printLog('#PHY',"%d Sequences output to %s\n" % (len(seqs),seqfile),1)
#########################################################################################################################
    def saveNexus(self,seqs=[],seqfile=None,name='short',id=False,log=True):  ### Saves sequences in Nexus format
        '''
        Saves sequences in SeqList object in nexus format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].nex] = filename
        >> name:str = Type of name to use as sequence name: 'short'=shortName(), 'AccNum'=AccNum, 'num'=Number
        >> id:boolean [False] = Appends sequence number to start of name.
        >> log:boolean [True] = Whether to log output
        '''
        try:
            ### Setup ###
            if seqs == []:
                seqs = self.seqs()
            if not self._checkAln(aln=True,tidygaps=False):
                self.printLog('#ALN','Nexus reformatting not supported for unaligned files.')
                raise ValueError
                
            ### Output File ###
            if not seqfile:
                seqfile = '%s.nex' % rje.baseFile(self.info['Name'],True)
            NEX = open(seqfile,'w')
            NEX.write('#NEXUS\nbegin data;\ndimensions ntax=%d nchar=%d;\n' % (self.seqNum(),self.seqLen()))
            NEX.write('format datatype=%s interleave=no gap=-;\nmatrix\n' % self.info['Type'])

            ### Output Sequences ###
            for seq in seqs:
                snum = '%d' % (self.seqs().index(seq)+1)
                if name == 'short':
                    outname = seq.shortName()
                else:
                    outname = seq.getData(name,str=True,default=snum)
                if id and outname != snum:
                    outname = '%s_%s' % (snum,outname)
                outname = string.join(string.split(outname),'_')
                NEX.write('%s %s\n' % (outname,seq.getSequence(case=self.opt['UseCase'])))
            NEX.write(';\nend;\n')
            NEX.close()

            ### Log ###
            if log:
                self.printLog('#NEX',"%d Sequences output to %s\n" % (len(seqs),seqfile),1)
            return True
        except:
            self.errorLog('Problem during rje_seq.saveNexus(%s)' % seqfile,quitchoice=True)
            return False
#########################################################################################################################
    def saveScanSeq(self,seqs=[],seqfile=None): ### Saves sequences in format for Scansite Parallel upload
        '''
        Saves sequences in format for Scansite Parallel upload.
        >> seqs:list of Sequence Objects. [self.seq]
        >> seqfile:str = Name of file ['%s.scanseq' % self.info['Basefile']]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqs == []: seqs = self.seqs()
            if seqfile == None: seqfile = '%s.scanseq' % self.info['Basefile']
            ### ~ [1] ~ Output sequences to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SEQOUT = open(seqfile, 'w')
            for seq in seqs: SEQOUT.write('%s %s\n' % (seq.info['AccNum'],re.sub('-','',seq.info['Sequence'])))
            SEQOUT.close()
            self.printLog('#OUT',"%d Sequences output to %s in scansite parallel format\n" % (len(seqs),seqfile),1)
        except:
            self.errorLog("Problem with saveScanSeq (%d seqs, file: %s)." % (len(seqs),seqfile))
            raise
#########################################################################################################################
    def saveAcc(self,seqs=[],accfile=None,scansite=False,uniprot=False,log=True): ### Saves accession numbers for UniProt retrieval etc
        '''
        Saves accession numbers for UniProt retrieval Scansite Parallel upload
        >> seqs:list of Sequence Objects. [self.seq]
        >> accfile:str = Name of file ['%s.acc' % self.info['Basefile']]
        >> scansite:boolean = whether to append a database identifier [False]
        >> uniprot:boolean = whether to output UniProt AccNum only [False]
        >> log:boolean = whether to output report to log [False]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqs == []: seqs = self.seqs()[0:]
            if accfile == None: accfile = '%s.acc' % self.info['Basefile']
            ### ~ [1] ~ Output accnum to file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ACCOUT = open(accfile, 'w')
            for seq in seqs[0:]:
                if uniprot and seq.info['DBase'][0:3] == 'ens': seqs.remove(seq)
                else: ACCOUT.write(seq.info['AccNum'])
                if scansite:
                    if seq.info['DBase'][0:3] == 'ens': ACCOUT.write(' EN')
                    elif seq.info['DBase'] in ['sprot', 'trembl']: ACCOUT.write(' ST')
                ACCOUT.write('\n')
            ACCOUT.close()
            if log: self.printLog('#OUT',"%d Sequence AccNums output to %s.\n" % (len(seqs),accfile),1)
        except:
            self.errorLog("Problem with saveScanSeq (%d seqs, file: %s)." % (len(seqs),seqfile))
            raise
#########################################################################################################################
    def saveR(self,seqs=[],seqfile=None,name='Name',namelen=0,id=False,log=True,case=False):  ### Saves sequences in R-tdt format
        '''
        Saves sequences in SeqList object in TDT format for R conversion to PNG.
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> name:str = Type of name to use as sequence name: 'short'=shortName(), 'AccNum'=AccNum, 'Number' = Number only
        >> namelen:int [0] = max length of sequence name [0 = no max]
        >> id:boolean [False] = Appends sequence number to start of name.
        >> log:boolean [True] = Whether to log output
        >> case:boolen [False] = Whether to use self.dict['Case'] to set output case
        '''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if seqs == []: seqs = self.seqs()
        if not seqfile: seqfile = self.info['Name'] + '.tdt'
        seqdat = {}
        maxlen = 0
        headers = []
        for seq in seqs:
            if name == 'short': outname = seq.shortName()
            elif name == 'Number': outname = '%d' % (self.seqs().index(seq)+1); id = False
            else: outname = seq.info[name]
            if id: outname = '%d %s' % (self.seqs().index(seq)+1,outname)
            if (namelen > 0) & (len(outname) >= namelen): outname = outname[0:(namelen-3)] + '...'
            headers.append(outname)
            seqdat[seq] = seq.getSequence(case)
            maxlen = max(maxlen,len(seqdat[seq]))
        rje.delimitedFileOutput(self,seqfile,headers,'\t',rje_backup=True)
        ### ~ [2] ~ Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for r in range(maxlen):
            seqtdt = []
            for seq in seqs:
                try: seqtdt.append(seqdat[seq][r])
                except: seqtdt.append('-')
            open(seqfile,'a').write('%s\n' % string.join(seqtdt,'\t'))
            #rje.delimitedFileOutput(self,seqfile,headers,'\t',datadict)
        if log: self.printLog('#FAS',"%d Sequences output to %s\n" % (len(seqs),seqfile),1)
#########################################################################################################################
    def seqNameDic(self,key='short',proglog=True):   ### Returns a dictionary of seqName:Sequence object for self.seq
        '''
        Returns a dictionary of seqName:Sequence object for self.seq.
        >> key:str = type of name to use as key:
            'short' = seq.shortName(), else uses seq.info[key]
            'NumName' = trim off leading 'X '
            'UniProt' = Original UniProt IDs.
            'Max' = return a dictionary that has shortNames, IDs and AccNums as keys!
        >> proglog:bool [True] = whether to print output to log
        '''
        try:
            seqdic = {}
            (sx,seqnum) = (0.0,self.seqNum())
            for seq in self.seqs():
                if proglog: self.printLog('\r#DICT','Making "%s" seqName dictionary: %.1f%%' % (key,sx/seqnum),newline=False,log=False)
                sx += 100.0
                if key == 'short': seqdic[seq.shortName()] = seq
                elif key == 'NumName':
                    name = seq.info['Name']
                    if re.search('^\d+\s(\S.+)$',name):
                        name = rje.matchExp('^\d+\s(\S.+)$',name)[0]
                    seqdic[name] = seq
                elif key == 'UniProt':
                    seqdic[seq.info['ID']] = seq
                    seqdic['%s_%s' % (seq.info['AccNum'],seq.info['SpecCode'])] = seq
                elif key == 'Max':
                    seqdic[seq.shortName()] = seq
                    seqdic[seq.info['ID']] = seq
                    seqdic[seq.info['AccNum']] = seq
                    try:
                        if seq.info['NCBI']: seqdic[seq.info['NCBI']] = seq
                    except: pass
                elif seq.info.has_key(key): seqdic[seq.info[key]] = seq
                #!# else will return an empty dictionary #!#
            if seqnum and proglog: self.printLog('\r#DICT','Made "%s" seqName dictionary: %s seq; %s keys.' % (key,rje.iStr(seqnum),rje.iStr(len(seqdic))))
            return seqdic
        except:
            self.errorLog('Problem during rje_seq.seqNameDic()')
            raise
#########################################################################################################################
    def reFormat(self,outfile=None,reformat=None,split=True): ### Saves sequences in appropriate format using self.attributes
        '''
        Saves sequences in appropriate format using self.attributes.
        >> outfile:str = Name of file  ['%s.*' % self.info['Basefile']]
        >> reformat:str = New format of file (fasta/phylip/scanseq/acclist/idlist/fastacmd/teiresias/6rf/3rf/est6rf)
        >> split:boolean = whether to use self.stat['Split'] [True]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            delimit = rje.getDelimit(self.cmd_list,',')
            seq_ext['mysql'] = rje.delimitExt(delimit)
            if not outfile: outfile = self.info['SeqOut']       # Set filename
            if not reformat: reformat = self.info['ReFormat']   # Changed during autofilter if split > 0 or seqout
            if reformat == 'None': return                       # Quit if no need to reformat
            elif reformat not in seq_ext.keys():
                return self.errorLog('Cannot reformat: "%s" format not recognised!' % reformat,printerror=False)
            if outfile == 'None': outfile = self.info['Name']
            if outfile == self.info['Name']: outfile = '%s.%s' % (self.info['Basefile'],seq_ext[reformat])
            ## ~ [0a] ~ Splitseq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            splitfile = 0
            splitx = 0
            self.stat['Split'] = int(self.stat['Split'])
            if self.stat['Split'] > 0: splitfile = 1
            if splitfile > 0:
                (outbase,outext) = os.path.splitext(outfile)
                outfile = '%s.%d%s' % (outbase,splitfile,outext)

            ### ~ [1] ~ Special Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] ~ Phylip format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if reformat == 'phylip':
                if self.opt['MemSaver']: return self.printLog('#MEM','Phylip reformatting not supported for MemSaver=T.')
                if splitfile > 0: self.printLog('#SPLIT','Phylip reformatting not supported for file splitting.')
                return self.savePhylip(seqfile=outfile)
            ## ~ [1b] ~ UniProt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if reformat == 'uniprot':
                if self.opt['MemSaver']: return self.printLog('#MEM','UniProt reformatting not supported for MemSaver=T.')
                uniprot = rje_uniprot.UniProt(self.log,self.cmd_list)
                for seq in self.seqs(): uniprot.addFromSeq(seq)
                return uniprot.saveUniProt(outfile,append=self.opt['Append'])
            ## ~ [1c] ~ Nexus ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if reformat == 'nexus':
                if self.opt['MemSaver']: self.printLog('#MEM','Nexus reformatting not supported for MemSaver=T.')
                elif splitfile > 0: self.printLog('#SPLIT','Nexus reformatting not supported for file splitting.')
                else: self.saveNexus(seqfile=outfile)
                return
                    
            ### ~ [2] ~ General Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.verbose(0,3,'Reformatting %s to %s...' % (self.info['Name'],outfile),0)
            ## ~ [2a] ~ Prepare ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0
            tmpout = False  # Whether had to make a temp output file to save reading and writing from same file
            seq = None
            if self.opt['MemSaver']:
                if outfile == self.info['Name']:
                    outfile = '%s.%s.tmp' % (outfile,rje.randomString(6))
                    tmpout = True
                self.seq = []
                SEQIN = open(self.info['Name'],'r')
                lastline = 'Firstline'
                (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            elif self.seqNum() > 0: seq = self.seqs()[0]
            ## ~ [2b] ~ Loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.opt['Append'] and splitfile < 1: SEQOUT = open(outfile, 'a')
            else:
                SEQOUT = open(outfile, 'w')
                if reformat == 'mysql':
                    rje.writeDelimit(SEQOUT,['fastacmd','protein_id','acc_num','spec_code','description'],delimit)
            while seq:
                ## Sequence Output ##
                self.formatOut(SEQOUT,seq,reformat,delimit)
                sx += 1
                splitx += 1
                rje.progressPrint(self,sx)
                if self.opt['MemSaver']:
                    self.seq = []
                    (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
                elif sx < self.seqNum(): seq = self.seqs()[sx]
                else: seq = None
                ## Check splitfile
                if splitfile > 0 and splitx >= self.stat['Split']:
                    self.printLog('#OUT','%d Sequences output to %s in %s format' % (splitx,outfile,reformat))
                    splitfile += 1
                    splitx = 0
                    outfile = '%s.%d%s' % (outbase,splitfile,outext)
                    self.verbose(0,3,'... %s...' % (outfile),0)
                    SEQOUT.close()
                    SEQOUT = open(outfile, 'w')

            ### ~ [3] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.verbose(0,1,'Done!',1)
            SEQOUT.close()
            if self.opt['MemSaver']:
                SEQIN.close()
                if tmpout:
                    os.rename(outfile,self.info['Name'])
                    outfile = self.info['Name']
            self.printLog('#OUT','%d Sequences output to %s in %s format' % (splitx,outfile,reformat)) 
        except:
            self.errorLog('Calamity during rje_seq.reFormat(%s)' % reformat)
#########################################################################################################################
    def formatOut(self,SEQOUT,seq,format,delimit=','):    ### Saves seq to open file SEQOUT in format
        '''
        Saves seq to open file SEQOUT in format.
        >> SEQOUT:file handle = output file, open for writing
        >> seq:Sequence Object to be saved
        >> format:str = format of file (fasta/phylip/scanseq/acclist/idlist/fastacmd/teiresias)
        >> delimit:str = text delimiter for MySQL format
        '''
        try:
            ### fasta ###
            if format in ['fas','fasta']:
                SEQOUT.write('>%s\n%s\n' % (seq.info['Name'],seq.getSequence(case=self.opt['UseCase'])))
            ### 6RF Fasta ###
            elif format in ['6rf','est6rf','3rf']:
                dnaseq = seq.getSequence()
                if format == '6rf': rfdic = rje_sequence.sixFrameTranslation(dnaseq)
                elif format == '6rf': rfdic = threeFrameTranslation(dnaseq,minpoly=self.stat['MinPoly'])
                else: rfdic = rje_sequence.estTranslation(dnaseq,minpoly=self.stat['MinPoly'])
                for rf in [1,2,3,-1,-2,-3]:
                    if rf in rfdic and rfdic[rf] and rje_sequence.bestORF(rfdic[rf])[1] >= max(1,self.stat['MinORF']):
                        SEQOUT.write('>%s-RF%s %s\n%s\n' % (seq.shortName(),rf,seq.info['Description'],rfdic[rf]))
            ### phylip ###
            elif format == 'phylip':
                #!# Not supported #!#           
                self.errorLog('Phylip output format not yet supported!',printerror=False,quitchoice=True)
                raise ValueError
            ### scanseq ###
            elif format == 'scanseq':
                SEQOUT.write('%s %s\n' % (seq.info['AccNum'],re.sub('-','',seq.info['Sequence'])))
            ### acclist ###
            elif format in ['acc','acclist']:
                SEQOUT.write('%s\n' % seq.info['AccNum'])
            ### idlist ###
            elif format == 'idlist':
                SEQOUT.write('%s\n' % seq.info['ID'])
            ### speclist ###
            elif format == 'speclist':
                SEQOUT.write('%s\n' % seq.info['SpecCode'])
            ### fastacmd ###
            elif format == 'fastacmd':
                SEQOUT.write('%s\n' % seq.shortName())
            ### teiresias ###
            elif format == 'teiresias':
                SEQOUT.write('>%s 1\n%s\n' % (seq.shortName(),re.sub('-','',seq.info['Sequence'])))
            ### MySQL ###
            elif format == 'mysql':
                formatlist = [seq.shortName(),seq.info['ID'],seq.info['AccNum'],seq.info['SpecCode'],seq.info['Description']]
                rje.writeDelimit(SEQOUT,formatlist,delimit)
        except: self.errorLog('Problem during rje_seq.formatOut(%s)' % format)
#########################################################################################################################
    ### <5> ### Sequence Information
#########################################################################################################################
    def numbersForNames(self,check1toN=False):  ### Returns True if all names are purely numbers, else False
        '''
        Returns True if all names are purely numbers, else False.
        >> check1toN:boolean = whether to also check that the numbers are 1 to N [False]
        '''
        numlist = []
        for seq in self.seqs():
            if re.search('\D',seq.info['Name']): return False
            if check1toN: numlist.append(string.atoi(seq.info['Name']))
        numlist.sort()
        if check1toN and numlist != range(1,self.seqNum()+1): return False
        return True
#########################################################################################################################
    def fullDetails(self):  ### Displays details of SeqList and all Sequences
        '''Displays details of SeqList and all Sequences.'''
        self.verbose(0,2,self.details(),0)
        for seq in self.seqs(): self.verbose(2,2,seq.details(),0)
        self.verbose(0,1,'\n%s. %d Sequences.' % (self.info['Name'],self.seqNum()),2)
#########################################################################################################################
    def _checkAln(self,aln=False,realign=False,tidygaps=True,seqkey='Sequence'):    ### Checks whether the sequences are aligned.
        '''
        Checks whether the sequences are aligned using:
        (1) Presence of Gaps
        (2) Equal lengths of sequences.
        >> aln:boolean [False] = whether sequences are 'meant' to be aligned
        >> realign:boolean [False] = whether to realign if aln=True and sequences aren't aligned
        >> tidygaps:boolean [True] = whether to tidy 100% gapped columns if already aligned and meant be
        >> seqkey:str ['Sequence'] = seq.info key to use to check alignment
        '''
        try:### ~ [1] Setup Alignment/Gap attributes. Check sequences present. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tidygaps = tidygaps and not self.opt['UseCase']
            self.opt['Aligned'] = False
            self.opt['Gapped'] = False
            if self.seqNum() < 1: return not aln      # No sequences!
            has_gaps = False   # Whether sequences contains gaps
            equal_len = True  # Whether all sequences of same length
            baselen = self.seqs()[0].seqLen()
            ### ~ [2] Check for equal lengths and gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                if seq.info[seqkey].find('-') >= 0: has_gaps = True
                if seq.seqLen() != baselen: equal_len = False
                if not aln and not equal_len and not has_gaps: return True      # Quick exit if alignment not desired
            ### ~ [3] Process state of alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.opt['Aligned'] = equal_len
            self.opt['Gapped'] = has_gaps
            if aln and equal_len and tidygaps: self.tidyGaps()
            if aln and not equal_len:
                if realign:
                    self.align(outfile=None,mapseq=True)
                    return self._checkAln(aln=True)
                else:
                    self.errorLog("Sequences should be aligned but are of different lengths!",printerror=False)
                    return False
            return True
        except:
            self.errorLog('Problem checking Alignment Status.')
            if aln: raise
#########################################################################################################################
    def aaFreq(self,alphabet=None,fromfile=None,loadfile=None,total=False):     ### Returns dictionary of AA (& gap etc.) frequencies
        '''
        Returns dictionary of AA (& gap etc.) frequencies.
        >> alphabet:list [None] = list of characters of interest
        - if alphabet == None, will return all characters found in seq
        >> fromfile:str = File from which to read sequences. If None, will use self.seq. [None]
        >> loadfile:str = File of aa frequencies from which to read frequencies. If None, will generate from self.seq [None]
        >> total:boolean = Whether to return additional element {'Total':aax} [False]
        << aafreq:dic = dictionary of frequency of characters
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aax = 0
            aafreq = {}
            ## ~ [1a] Basefile for AA Frequency output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if fromfile and not os.path.exists(fromfile):
                self.errorLog('Sequence file %s for AA Frequency calculation is missing!' % fromfile,printerror=False)
                raise IOError
            if fromfile: basefile = rje.baseFile(fromfile)      # Use given sequence file as source
            else:                                               # Else, use self
                fromfile = self.info['Name']
                basefile = self.info['Basefile']
            ## ~ [1b] Check existence of file from which to load frequencies ~~~~~~~~~~~~~~~~~~~~~~ ##
            if loadfile:
                if not os.path.exists(loadfile):
                    self.errorLog('AA Frequency file %s missing!' % loadfile,printerror=False)
                    raise IOError
                aafile = loadfile
            else: aafile = '%s.aafreq.txt' % basefile
            ## ~ [1c] Setup alphabet for aafreq dictionary keys ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if alphabet == None: newkeys = True
            else:
                newkeys = False
                for letter in alphabet: aafreq[letter] = 0.0

            ### ~ [2] Load AA Frequencies if present ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if loadfile or rje.isYounger(aafile,fromfile):      # Freq file exists (and is newer than sequence file)
                ## ~ [2a] Load data into dictionary and establish Frequency column ~~~~~~~~~~~~~~~~ ##
                loadfreq = rje.dataDict(self,aafile,getheaders=True)
                if len(loadfreq['Headers']) == 1: loadkey = loadfreq['Headers'][0:]
                else:
                    for trykey in loadfreq['Headers']:
                        if trykey.lower()[:4] == 'freq':
                            loadkey = trykey
                            break
                loadfreq.pop('Headers')
                ## ~ [2b] Convert to numerical dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for aa in loadfreq:     #!# At some point, read whether AA or NT and use "NT" #!#
                    try:
                        if aa != 'Total' and (aafreq.has_key(aa) or newkeys): aafreq[aa] = string.atof(loadfreq[aa][loadkey])
                        elif total: aafreq[aa] = string.atoi(loadfreq[aa][loadkey])
                    except: self.errorLog('Problem with AA key "%s" (="%s")' % (aa,loadfreq[aa][loadkey]))
                    #X#self.log.printLog('\r#AA','Loading AA Frequencies from %s: %.1f%%' % (aafile,float(aalines.index(line))/len(aalines)),log=False,newline=False)
                if aafreq.has_key('Total'): self.printLog('\r#AA','Loaded %d AA Frequencies from %s. %s aa total.' % (len(aafreq)-1,aafile,rje.integerString(aafreq['Total'])))
                else: self.printLog('\r#AA','Loaded %d AA Frequencies from %s.' % (len(aafreq),aafile))
                ## ~ [2c] Special addition of Total if desired ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if total and not aafreq.has_key('Total'):
                    if fromfile == self.info['Name']: aafreq['Total'] = self.aaTotal()
                    else: aafreq['Total'] = DBSize(self,fromfile)
                ## ~ [2d] Normalise and reduce alphabet if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                try: atot = aafreq.pop('Total')
                except: atot = 0.0
                aafreq = rje.dictFreq(aafreq,total=False)
                if total: aafreq['Total'] = atot
                return aafreq

            #!# Tidy from here onwards... #!#
            ### <2> ### Count AAs
            if fromfile != self.info['Name']:
                _stage = '<1a> %s AAs' % fromfile
                aaseq = SeqList(cmd_list=['i=-1','v=-1'])
                AASEQ = open(fromfile,'r')
                (seq,lastline) = aaseq.nextFasSeq(fileobject=AASEQ,lastline='')
                while seq:
                    self.printLog('\r#AA','Calculating AA Frequencies from %s: %s aa' % (fromfile,rje.integerString(aax)),log=False,newline=False)
                    for r in range(seq.seqLen()):
                        aa = seq.info['Sequence'][r]
                        if aafreq.has_key(aa):
                            aafreq[aa] += 1
                            aax += 1
                        elif newkeys:
                            aafreq[aa] = 1.0
                            aax += 1
                    aaseq.seq = []
                    (seq,lastline) = aaseq.nextFasSeq(fileobject=AASEQ,lastline=lastline)
                self.printLog('\r#AA','Calculated AA Frequencies from %s: %s aa.' % (fromfile,rje.integerString(aax)))
                AASEQ.close()
            else:                
                _stage = '<1b> self.seq AAs'
                if self.seqNum() == 0:
                    self.errorLog('No sequences from which to calculate AAFreq!',printerror=False)
                    return {}
                for seq in self.seqs():
                    self.printLog('\r#AA','Calculating AA Frequencies from %s seq: %s aa' % (rje.integerString(self.seqNum()),rje.integerString(aax)),log=False,newline=False)
                    for r in range(seq.seqLen()):
                        aa = seq.info['Sequence'][r]
                        if aafreq.has_key(aa):
                            aafreq[aa] += 1
                            aax += 1
                        elif newkeys:
                            aafreq[aa] = 1.0
                            aax += 1
                self.printLog('\r#AA','Calculated AA Frequencies from %s seq: %s aa.' % (rje.integerString(self.seqNum()),rje.integerString(aax)))

            ### <3> ### Convert to Freqs
            for k in aafreq.keys(): aafreq[k] /= aax
            aafreq['Total'] = aax

            ### <4> ### Save Frequencies
            try:
                AAFREQ = open(aafile,'w')
                AAFREQ.write('AA\tFREQ\n')
                for aa in rje.sortKeys(aafreq):
                    if aa != 'Total': AAFREQ.write('%s\t%f\n' % (aa,aafreq[aa]))
                AAFREQ.write('Total\t%d\n' % aafreq['Total'])
                AAFREQ.close()
                self.printLog('\r#AA','AA Frequencies saved to %s.' % aafile)
            except:
                pass
                        
            if not total:
                aafreq.pop('Total')
            return aafreq
        except:
            self.errorLog('Major error during SeqList.aaFreq()')
            raise
#########################################################################################################################
    def seqAlnPos(self,seq,alnpos,next=True):   ### Returns corresponding position of sequence in alignment
        '''
        Returns corresponding position of sequence in alignment.
        >> seq:Sequence object
        >> alnpos:int = Position (0->L) in alignment
        >> next:bool [True] = Whether to return next position if gap (else returns -1)
        << returns position in sequence (0->L)
        '''
        try:### ~ [1] Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seq.info['Sequence'][alnpos] == '-' and not next: return -1
            p = -1
            for r in range(alnpos+1):
                if seq.info['Sequence'][r] != '-': p += 1
            if seq.info['Sequence'][alnpos] == '-': p += 1
            return p
        except:
            self.errorLog('Problem during SeqList.seqAlnPos()')
            return -1
#########################################################################################################################
    def aaCount(self):     ### Returns total count of AA in seqlist
        '''Returns total count of AA in seqlist.'''
        try:
            aacount = 0
            for seq in self.seqs(): aacount += seq.aaLen()
            return aacount
        except:
            self.errorLog("Disaster during aaCount! Continuing with count of zero.")
            return 0
#########################################################################################################################
    def setupSubDict(self,masking=True,alphabet=[]):    ### Sets up and returns query subsitution frequency dictionary
        '''
        Sets up and returns query subsitution frequency dictionary.
        >> masking:bool [True] = whether to use qry.info['MaskSeq'] if found.
        >> alphabet:list [] = Alphabet used for dictionary keys. 
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self._checkAln(True,tidygaps=False): raise ValueError    # Must be aligned!
            if self.obj['QuerySeq']: qry = self.obj['QuerySeq']
            else: qry = self.querySeq()
            if not qry:
                self.errorLog('Cannot setup SubDict without query!',printerror=False)
                raise ValueError
            subdict = {}
            if masking and 'MaskSeq' in qry.info:
                qryseq = qry.info['MaskSeq'].upper()
                if 'X' in alphabet: alphabet.remove('X')
            else: qryseq = qry.info['Sequence'].upper()
            for seq in self.seqs():
                if seq == qry: continue
                subdict[seq] = {}
                for aa in string.split(string.join(alphabet).upper()): subdict[seq][aa] = {}

            ### ~ [2] Make raw number dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for r in range(len(qryseq)):
                aa = qryseq[r]
                if aa in ['-','X']: continue    # Gap or masked residue
                ## ~ [2a] Look at residue in each aligned sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##                           
                for seq in self.seqs():
                    if seq == qry: continue
                    sa = seq.info['Sequence'].upper()[r]
                    if aa not in subdict[seq]: subdict[seq][aa] = {sa:1}
                    elif sa in subdict[seq][aa]: subdict[seq][aa][sa] += 1
                    else: subdict[seq][aa][sa] = 1
                
            ### ~ [3] Normalise substitution numbers to frequencies ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                if seq == qry: continue
                for aa in subdict[seq]: subdict[seq][aa] = rje.dictFreq(subdict[seq][aa],total=False)
            return subdict

        except:
            self.errorLog('Problem during SeqList.setupSubDict()')
            raise
#########################################################################################################################
    ### <6> ### Sequence Filters
#########################################################################################################################
    def maxX(self):     ### Removes sequences with too many Xs
        '''Removes sequences with too many Xs.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MaxX'] <= 0: return
            ### ~ [1] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs()[0:]:
                xx = string.count(seq.info['Sequence'].upper(),'X')
                if (float(xx) / seq.aaNum()) > self.stat['MaxX']:
                    self.removeSeq('Proportion of Xs too great for MaxX (%.2f%%)' % (100.0*self.stat['MaxX']),seq=seq)
                elif self.stat['MaxX'] >= 1 and (xx) > self.stat['MaxX']:
                    self.removeSeq('Too many Xs for MaxX (%d)' % self.stat['MaxX'],seq=seq)
        except: self.errorLog('Problem during maxX()')
#########################################################################################################################
    def maxGlob(self,inverse=False):  ### Removes sequences with too much predicted order
        '''Removes sequences with too much predicted order.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MaxGlob'] <= 0: return
            ### ~ [1] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs()[0:]:
                g = seq.globProportion(absolute=True)
                if (float(g) / seq.aaNum()) > self.stat['MaxGlob']:
                    self.removeSeq('Proportion of Order too great for MaxGlob (%.2f%%)' % (100.0*self.stat['MaxGlob']),seq=seq)
                elif self.stat['MaxGlob'] >= 1 and (xx) > self.stat['MaxGlob']:
                    self.removeSeq('Too many Ordered residues for MaxGlob (%d)' % self.stat['MaxGlob'],seq=seq)
        except: self.errorLog('Problem during maxGlob()')
#########################################################################################################################
    def gapSeqFilter(self,relative='query',keepqry=True):     ### Removes gappy sequences, relative to query if given
        '''
        Removes gappy sequences, relative to query if given.
        >> relative:str = gap measure relative to:
            - self = own sequence alone (as part of alignment)
            - query = Query sequence
            - neighbour = closest neighbour in alignment (%ID)
        >> keepqry:bool [True] = whether to keep query no matter what
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.stat['MaxGap'] <= 0 or not self.opt['Aligned']: return
            gcut = self.stat['MaxGap']
            if gcut < 1: gcut *= 100.0
            gtxt = 'Too many gaps (> %.1f%%) relative to %s' % (gcut,relative)
            ltxt = 'Removing gappy (> %.1f%%) sequences relative to %s' % (gcut,relative)
            if relative == 'query':
                stxt = ' (%s)' % self.obj['QuerySeq'].shortName()
                ltxt += ' (%s)' % self.obj['QuerySeq'].shortName()
            else: stxt = ''
            if not self.obj['MSA Gaps']: self.addMatrix('MSA Gaps',sym=False)
            if not self.obj['MSA ID']: self.addMatrix('MSA ID',sym=False)
            self.verbose(1,3,'%s...' % ltxt,1)
            ### ~ [1] ~ Gappy removal ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs()[0:]:
                if seq == self.obj['QuerySeq'] and keepqry: continue
                gcount = string.count(seq.info['Sequence'],'-')
                if relative == 'query':
                    if seq == self.obj['QuerySeq'] or (100.0 * float(gcount) / self.obj['QuerySeq'].aaLen()) <= gcut: continue
                    gaps = self.getDis(seq,self.obj['QuerySeq'],key='MSA Gaps')
                elif relative == 'self': gaps = 100.0 * gcount / seq.seqLen()
                else:
                    ## Quick check ##
                    quick_ok = True
                    for otherseq in self.seqs()[0:]:
                        if seq != otherseq and (100.0 * float(gcount) / otherseq.aaLen()) > gcut:
                            quick_ok = False
                            break
                    if quick_ok: continue
                    bestid = 0
                    gaps = 100.0
                    for otherseq in self.seqs()[0:]:
                        if seq == otherseq: continue
                        compid = self.getDis(seq,otherseq,key='MSA ID')
                        compgap = self.getDis(seq,otherseq,key='MSA Gaps')
                        if compid >= bestid:
                            bestid = compid
                            if compgap < gaps: gaps = compgap
                        stxt = ' (%s)' % otherseq.shortName()
                if gaps > gcut: self.removeSeq(text='%s%s! (%d%%)' % (gtxt,stxt,gaps),seq=seq)
            ### ~ [2] ~ Finish and check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.verbose(1,3,'Gappy Sequence removal complete! %d sequences remain.' % self.seqNum(),1)
            self._checkAln(aln=True)
        except:
            self.errorLog('Major Problem in rje_seq.GapFilter()')
            raise
#########################################################################################################################
    ### <7> ### Sequence Redundancy
#########################################################################################################################
    def memSaveNR(self):    ### Unaligned sequence redundancy method for memsaver mode (add more options with time)
        '''
        Unaligned sequence redundancy method for memsaver mode.
        '''
        try:
            ### Setup ###
            _stage = 'Setup'
            ## Removal/Redundancy Output ##
            remtext = 'Redundancy Check'    #!# No option for extra text yet
            #!# Removed possibility to save filtered (redundant) sequences: add again at end #!#
            ## NR ID ##
            nrid = float(self.stat['NR ID'])
            if nrid <= 0:
                nrid = 101.0
            elif nrid < 1:
                nrid = 100.0 * nrid
            nrsim = float(self.stat['NR Sim'])
            if nrsim <= 0:
                nrsim = 101.0
            elif nrsim < 1:
                nrsim = 100.0 * nrsim
            ## Species Filter ##
            spectxt = ''
            samespec = self.opt['SpecNR']
            if samespec:
                spectxt = ' (same species only)'

            ### AccNum NR ### First pass is to remove redundant AccNums ###
            accnr = False   #!# Too annoying for now. Add later? #!#
            _stage = 'AccNum Redundancy'
            fullacclist = []
            redacc = []
            SEQIN = open(self.info['Name'],'r')
            lastline = 'Firstline'
            self.verbose(0,3,'Checking %s for AccNum redundancy...' % (self.info['Name']),0)
            (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            sx = 0
            while seq:
                #if seq.info['AccNum'] in fullacclist:   # Redundant AccNum! #
                #    redacc.append(seq.info['AccNum'])
                #else:
                #    fullacclist.append(seq.info['AccNum'])
                sx += 1
                rje.progressPrint(self,sx)
                self.seq = []
                (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            SEQIN.close()
            self.verbose(0,3,'%d of %d AccNums redundant.' % (len(redacc),sx),1)
            if len(redacc):
                self.verbose(0,0,'Redundant AccNums cannot yet be handled in MemSaver mode! Sorry!!',1)
                sys.exit
                #!# NB. These will be filtered later! #

            ### Sequence Redundancy ###            
            _stage = 'Sequence Redundancy'
            if nrid > 100.0 and nrsim > 100.0:
                self.verbose(0,3,'NR ID and NR Sim both <= 0%!!',1)
                return
            redseq = []     # List of redundant Sequences
            blastformat = True  # Possibly add a counter as before to reduce DB Size once certain no. of redundancies removed?
            SEQIN = open(self.info['Name'],'r')
            lastline = 'Firstline'
            self.verbose(0,3,'Checking %s for Sequence redundancy...' % (self.info['Name']),0)
            (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            blastx = sx - len(redacc)
            sx = 0            
            while seq:
                ## Housekeep sequence ##
                sx += 1
                rje.progressPrint(self,sx,dotx=10,numx=100,v=0)
                self.seq = []
                ## Check if already removed ##
                if seq.shortName() in redseq:   # Redundant and identified! #
                    (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
                    continue
                ## GABLAM against rest of seqlist ##
                ## <i> ## BLAST Format Database for searching ##
                if blastformat:   # Need to make new BLAST Database (of smaller size)
                    blastdb = re.sub('\.fas$','.blastdb',self.info['Name'])
                    shutil.copyfile(self.info['Name'],blastdb)
                    if self.info['Type'] == 'DNA':
                        rje_blast.formatDB(blastdb,self.info['BLAST Path'],protein=False,log=None,oldblast=True)
                    else:
                        rje_blast.formatDB(blastdb,self.info['BLAST Path'],protein=True,log=None,oldblast=True)
                ## <ii> ## Perform BLAST & GABLAMO ##
                seqblast = rje_blast.blastObj(log=self.log,cmd_list=self.cmd_list+['v=%d' % (self.stat['Verbose']-2)],type='Old')
                acc = seq.info['AccNum']
                tempseq = SeqList(cmd_list=['i=-1','v=-1'])
                tempseq.seq = [seq]
                tempseq.saveFasta(seqfile='%s.tmp.fas' % acc)
                seqblast.setInfo({'Name': '%s.tmp.txt' % acc, 'DBase': blastdb, 'Type': 'blastp', 'InFile': '%s.tmp.fas' % acc})
                if self.info['Type'] in ['DNA','RNA']:
                    seqblast.info['Type'] = 'blastn'
                seqblast.setStat({'OneLine': blastx, 'HitAln': blastx})   # Need alignments for GABLAMO
                seqblast.blast()
                if not seqblast.readBLAST(gablam=True): raise ValueError    ## Seqs now in order of seqblast.search[0].hit
                os.unlink('%s.tmp.txt' % acc)
                os.unlink('%s.tmp.fas' % acc)
                ## <iii> ## Work through hits
                hits = seqblast.search[0].hit
                for hit in hits:
                    hitname = hit.info['Name']
                    if hitname == seq.shortName():     # Ignore as query
                        continue
                    elif hitname in redseq:    # Already Redundant!
                        continue
                    elif samespec and (seq.info['SpecCode'] == 'UNK' or seq.info['SpecCode'] != rje_sequence.specCodeFromName(hitname)):
                        self.verbose(2,4,'%s Hit %d (%s): Wrong Species' % (seq.shortName(),hits.index(hit),hitname),1)
                        continue    # Not worth checking as different species
                    remredseq = False
                    remredtext = ''
                    gdict = hit.dict['GABLAM']
                    idq = 100 * float(gdict['Query']['GABLAMO ID']) / seq.aaLen()
                    simq = 100 * float(gdict['Query']['GABLAMO Sim']) / seq.aaLen()
                    idh = 100 * float(gdict['Hit']['GABLAMO ID']) / hit.stat['Length']
                    simh = 100 * float(gdict['Hit']['GABLAMO Sim']) / hit.stat['Length']
                    if idq >= nrid or idh >= nrid:
                        remredtext = '>=%.2f%% (%.2f%%/%.2f%%) ID' % (nrid,idq,idh)
                        remredseq = True
                    elif simq >= nrsim or simh >= nrsim:
                        remredtext = '>=%.2f%% (%.2f%%/%.2f%%) Sim' % (nrsim,simq,simh)
                        remredseq = True
                        
                ## <iv> ## Remove Seq ##
                    if remredseq:
                        otherseq = self.seqFromFastaCmd(hitname,blastdb)
                        rem = self._worstSeq(seq, otherseq)
                        redseq.append(rem[0].shortName())
                        self.printLog('\r#REM','Deleted redundant %s %s vs %s: %s\n' % (rem[0].shortName(),remredtext,rem[2],rem[1]),1)
                        if seq in rem:  ## Removed BLAST query! ##
                            break

                ## <vi> ## Next sequence
                (seq,lastline) = self.nextFasSeq(SEQIN,lastline)

            ## <vi> ## Tidy Up
            self.verbose(0,2,'%d Redundant Sequences.' % len(redseq),1)
            try:    #!# Use RJE_BLAST code #!#
                os.unlink('%s' % blastdb)
                if self.info['Type'] in ['DNA','RNA']:
                    os.unlink('%s.nsd' % blastdb)
                    os.unlink('%s.nhr' % blastdb)
                    os.unlink('%s.nin' % blastdb)
                    os.unlink('%s.nsi' % blastdb)
                    os.unlink('%s.nsq' % blastdb)
                else:
                    os.unlink('%s.psd' % blastdb)
                    os.unlink('%s.phr' % blastdb)
                    os.unlink('%s.pin' % blastdb)
                    os.unlink('%s.psi' % blastdb)
                    os.unlink('%s.psq' % blastdb)
            except(OSError):
                self.errorLog('Problem unlinking (deleting) blast (%s) files.' % blastdb)

            ### Finish & filter out bad sequences ###
            SEQIN.close()

            memout = '%s.%s.tmp.fas' % (self.info['Basefile'],rje.randomString(6))
            SEQIN = open(self.info['Name'],'r')
            lastline = 'Firstline'
            self.verbose(0,3,'Filtering redundant sequences from %s...' % (self.info['Name'],memout),0)
            SEQOUT = open(memout, 'w')
            (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            sx = 0
            okx = 0
            while seq:
                if seq.shortName() in redseq:
                    redseq.remove(seq.shortName())
                else:
                    SEQOUT.write('>%s\n%s\n' % (self.seqs()[-1].info['Name'],self.seqs()[-1].info['Sequence']))
                    okx += 1
                sx += 1
                rje.progressPrint(self,sx)
                self.seq = []
                (seq,lastline) = self.nextFasSeq(SEQIN,lastline)
            SEQIN.close()
            SEQOUT.close()
            os.rename(memout,self.info['Name'])
            self.verbose(0,3,'Done! %d NR seqs remain.' % okx,1)
                                    
        except:
            self.errorLog('Major problem during rje_seq.memSaveNR(%s).' % _stage,True)
            return       
#########################################################################################################################
    def makeNR(self,best_ann=True,text='',nrid=None,samespec=None,blast=100,pw_aln=True,save_red=False,check100=True,skip_to_seq=None,nrsim=None,nr_qry=None):   ### Makes list non-redundant
        '''
        Makes sequence list non-redundant. Works on %ID, calculating with align if necessary.
        If best_ann is selected, will choose the sequence with the best annotation:
        - sprot > trEMBL > ens_known > ens_novel > ens_scan
        Otherwise (or if even), will choose the longer sequence else the first sequence.
        >> best_ann:boolean [True] = whether to preferentially keep shorter but better annotated sequences.
        >> text:str = Extra text for removeSeq()
        >> nrid:float [100.0] = percentage identity cut-off for redundancy decision (self.stat['NR ID'])
        >> nrsim:float [100.0] = percentage similarity cut-off for redundancy decision (self.stat['NR Sim'])
        >> check100:boolean [True] = Whether to check for 100% matches by simple sequence comparison
        >> samespec:boolean [None] = Whether to apply cut-off to same species only. [Not check100] (self.opt['SpecNR'])
        >> blast:integer [100] = number of hits is BLAST screen used to identify most similar proteins before applying ID cut-off
            - 0 = no BLAST screen
        >> pw_aln:boolean [True] = will use ALIGN to calculate %ID (<100). If False will use current sequence alignment.
        >> save_red:boolean [False] = if True, will return Sequence Object containing redundant sequences
        >> skip_to_seq:Sequence Object [None] = if given, will not make any tests until Sequence reached.
        >> nr_qry:Sequence object [None] = if given, will only compare sequences to this one sequence (not each other)
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            silence = self.log.opt['Silent']
            ## ~ [0a] Removal/Redundancy Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            remtext = 'Redundancy Check'
            if text != '': remtext += ' (%s)' % text
            redseq = None
            if save_red:
                redseq = SeqList(log=self.log)
                redseq.info['Name'] = '%s.redseq.fas' % rje.baseFile(self.info['Name'],extlist=['fas'])
            ## ~ [0b] NR ID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if nrid == None: nrid = float(self.stat['NR ID'])
            if nrid <= 0: nrid = 101.0
            if nrid < 1: nrid = 100.0 * nrid
            if nrsim == None: nrsim = float(self.stat['NR Sim'])
            if nrsim <= 0: nrsim = 101.0
            if nrsim < 1: nrsim = 100.0 * nrsim
            ## ~ [0c] Species Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            spectxt = ''
            if samespec == None: samespec = self.opt['SpecNR']
            if samespec: spectxt = ' (same species only)'
            ## ~ [0d] DisMatrices ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if pw_aln:
                if not self.obj['PWAln ID']: self.addMatrix('PWAln ID')
            elif self.opt['Aligned']:
                if not self.obj['MSA ID']: self.addMatrix('MSA ID')
            ### ~ [1] 100% Identity first ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if check100 and not nr_qry: #!# nr_query not implemented! #!#
                ## ~ [1a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nrdic = {'All':{}}  # Dictionary per species code of sequence:seq object                
                red = 0
                ## ~ [1b] "Forward" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                nridseq = self.seqs()[0:]  #!# nr_query not implemented! #!#
                (sx,stot) = (0.0,len(nridseq))
                for seq in nridseq:
                    progress = (rje.integerString(self.seqNum()),spectxt,sx/stot,rje.integerString(red)); sx += 100.0
                    self.progLog('\r#NR','NR @ 100%% identity, %s Sequences%s ... %.1f%% (Phase I): %s removed' % progress)
                    ## Get appropriate species code for nrdic key
                    if samespec:
                        _seqspec = seq.info['SpecCode']
                        if _seqspec in ['UNK','Unknown']: # Skip unknown species
                            nrdic[seq.info['Name']] = {seq.info['Sequence']:seq}
                            continue
                        elif self.obj['QuerySeq'] and not self.opt['QueryNR'] and _seqspec == self.obj['QuerySeq'].info['SpecCode']:
                            nrdic[seq.info['Name']] = {seq.info['Sequence']:seq}
                            continue
                        elif _seqspec not in nrdic.keys():  # Make blank entry
                            nrdic[seq.info['SpecCode']] = {}
                    else: _seqspec = 'All'
                    
                    ## Skip to sequence if desired
                    if seq == skip_to_seq:
                        skip_to_seq = None
                        nrnum = 0
                        for code in nrdic.keys():
                            nrnum += len(nrdic[code])
                        self.verbose(0,4,'Reached %s' % seq.shortName(),0)
                        self.verbose(0,2,'(Index:%d. %d seqs in nrdic.' % (nridseq.index(seq),nrnum),0)
                    elif skip_to_seq:
                        nrdic[_seqspec][seq.info['Sequence']] = seq
                        continue
                        
                    ## Look for 100% Redundancy
                    nrstr = string.join(['']+nrdic[_seqspec].keys()+[''],' ')
                    while seq in self.seqs() and rje.matchExp('\s(\S*%s\S*)\s' % seq.info['Sequence'],nrstr):  # Sequence is as Fragment or has total overlap with a previous sequence
                        otherseq = nrdic[_seqspec][rje.matchExp('\s(\S*%s\S*)\s' % seq.info['Sequence'],nrstr)[0]]
                        rem = self._worstSeq(otherseq, seq, best_ann)
                        #self.verbose(1,5,'%s of %d : ' % (rje.preZero(self.seq.index(seq)+1,self.seqNum()),self.seqNum()),0)
                        red += 1
                        if seq in rem:
                            if save_red: redseq.seq.append(seq)
                            if self.opt['NRKeepAnn']: otherseq.info['Name'] += ' >%s' % seq.info['Name']
                            self.removeSeq(text='100%% ID vs %s: %s' % (otherseq.shortName(),rem[1]),seq=seq)
                        else:
                            if save_red: redseq.seq.append(otherseq)
                            if self.opt['NRKeepAnn']: seq.info['Name'] += ' >%s' % otherseq.info['Name']
                            self.removeSeq(text='100%% ID vs %s: %s' % (seq.shortName(),rem[1]),seq=otherseq)
                            nrdic[_seqspec].pop(otherseq.info['Sequence'])
                            nrstr = string.join(['']+nrdic[_seqspec].keys()+[''],' ')
                    if seq in self.seqs():
                        #print nrdic, seq.info['SpecCode'], nrdic[seq.info['SpecCode']], seq.info['Sequence']
                        nrdic[_seqspec][seq.info['Sequence']] = seq

                fwdseqn = 0
                for code in nrdic.keys(): fwdseqn += len(nrdic[code])
                if self.seqNum() != fwdseqn:
                    self.errorLog('NR WARNING: Should have %d seqs in nrdic but have %d. Safe to continue but some redundant fragments may be missed.' % (self.seqNum(),fwdseqn),quitchoice=False,printerror=False)
                progress = (rje.integerString(self.seqNum()),spectxt,100.0,rje.integerString(red))
                self.printLog('\r#NR','Searched for 100%% identity among %s Sequences%s ... %.1f%% (Phase I): %s removed.  ' % progress,log=False,newline=True)

                ## <c> ## "Backward"
                _stage = '<1c> 100% Identity - Bwd'
                #X#self.verbose(0,3,'...Phase 2 (%d seqs)...' % fwdseqn,0)
                perc = rje.setPerc(fwdseqn,100,1000,self.stat['Verbose'])
                for seq in nridseq:
                    progress = (rje.integerString(self.seqNum()),spectxt,(100.0 * (1 + nridseq.index(seq)) / len(nridseq)),rje.integerString(red))
                    self.printLog('\r#NR','NR @ 100%% identity, %s Sequences%s ... %.1f%% (Phase II): %s removed' % progress,log=False,newline=False)
                    if seq not in self.seqs(): # Already removed
                        continue
                    ## Get appropriate species code for nrdic key
                    if samespec:
                        _seqspec = seq.info['SpecCode']
                        if _seqspec in ['UNK','Unknown']: # Skip unknown species
                            nrdic[seq.info['Name']].pop(seq.info['Sequence'])
                            continue
                        elif self.obj['QuerySeq'] and not self.opt['QueryNR'] and _seqspec == self.obj['QuerySeq'].info['SpecCode']:
                            nrdic[seq.info['Name']].pop(seq.info['Sequence'])
                            continue
                    else: _seqspec = 'All'
                    ## Remove sequence to prevent self-hit
                    if nrdic.has_key(_seqspec) and nrdic[_seqspec].has_key(seq.info['Sequence']):
                        nrdic[_seqspec].pop(seq.info['Sequence'])
                    else:
                        self.errorLog('Sequence %s missing from nrdic[%s]!' % (seq.shortName(),_seqspec),quitchoice=False,printerror=False)
                        continue
                    ## Look for sequence in later members
                    nrstr = string.join(['']+nrdic[_seqspec].keys()+[''],' ')
                    while rje.matchExp('\s(\S*%s\S*)\s' % seq.info['Sequence'],nrstr):  # Fragment or total overlap with sequence
                        otherseq = nrdic[_seqspec][rje.matchExp('\s(\S*%s\S*)\s' % seq.info['Sequence'],nrstr)[0]]
                        rem = self._worstSeq(otherseq, seq, best_ann)
                        self.verbose(1,5,'%s of %d : ' % (rje.preZero(self.seqs().index(seq)+1,self.seqNum()),self.seqNum()),0)
                        red += 1
                        if seq in rem:
                            if save_red: redseq.seq.append(seq)
                            if self.opt['NRKeepAnn']: otherseq.info['Name'] += ' >%s' % seq.info['Name']
                            self.removeSeq(text='100%% ID vs %s: %s' % (otherseq.shortName(),rem[1]),seq=seq)
                            break
                        else:
                            if save_red: redseq.seq.append(otherseq)
                            if self.opt['NRKeepAnn']: seq.info['Name'] += ' >%s' % otherseq.info['Name']
                            self.removeSeq(text='100%% ID vs %s: %s' % (seq.shortName(),rem[1]),seq=otherseq)
                            nrdic[_seqspec].pop(otherseq.info['Sequence'])
                            nrstr = string.join(['']+nrdic[_seqspec].keys()+[''],' ')
                ## <d> ## Finish
                bwdseqn = 0
                for code in nrdic.keys(): bwdseqn += len(nrdic[code])
                if bwdseqn: self.errorLog('NR WARNING: Should have 0 seqs left in nrdic but have %d. Safe to continue but some redundant fragments may be missed.' % (bwdseqn),quitchoice=False,printerror=False)
                progress = (rje.integerString(self.seqNum()),spectxt,rje.integerString(red))
                self.printLog('\r#NR','Searched for 100%% identity among %s Sequences%s (Both Phases): %s removed.          ' % progress,log=True,newline=True)

        ### <2> ### Lesser Identity Cut-off
            blastdb = None
            _stage = '<2> < 100% Identity'
            if (nrid < 100.0 or nrsim < 100.0 or check100==False) and (pw_aln or self.opt['Aligned']):
                (red,blastred) = (0,blast)    # (Total,since_formatDB)
                if nr_qry: nridseq = [nr_qry]
                else: nridseq = self.seqs()[0:]
                for seq in nridseq:
                    progress = (nrid,nrsim,rje.integerString(self.seqNum()),spectxt,(100.0 * (1 + nridseq.index(seq)) / len(nridseq)),rje.integerString(red))
                    self.log.opt['Silent'] = silence
                    self.printLog('\r#NR','Searching for %.2f%% identity / %.2f%% similarity among %s Sequences%s ... %.1f%% : %s removed.' % progress,log=False,newline=False)
                    fullhits = True     # Whether full number of BLAST hits has been reached and still removing redundancy
                    if seq == skip_to_seq:
                        skip_to_seq = None
                    if (seq != nr_qry and seq not in self.seqs()) or skip_to_seq != None:  # Not reached or been removed!
                        continue
                    
                    while fullhits:
                        self.log.opt['Silent'] = True
                        fullhits = False

                    ## <i> ## BLAST Format Database for searching
                        if blastred >= blast and blast > 0:   # Need to make new BLAST Database (of smaller size)
                            blastdb = re.sub('\.fas$','.blastdb',self.info['Name'])
                            blastseq = SeqList(cmd_list=['i=-1','v=-1'])
                            blastseq.seq = self.seqs()[0:]
                            blastseq.degapSeq()
                            blastseq.saveFasta(seqfile=blastdb)
                            if self.info['Type'] == 'DNA':
                                rje_blast.formatDB(blastdb,self.info['BLAST Path'],protein=False,log=None,oldblast=True)
                            else:
                                rje_blast.formatDB(blastdb,self.info['BLAST Path'],protein=True,log=None,oldblast=True)
                            blastred = 0

                    ## <ii> ## Work through sequence comparisons by BLAST or ALL
                        if blast > 0:
                            seqblast = rje_blast.blastObj(log=self.log,cmd_list=self.cmd_list+['v=%d' % (self.stat['Verbose']-2)],type='Old')
                            acc = seq.info['AccNum']
                            tempseq = SeqList(cmd_list=['i=-1','v=-1'])
                            tempseq.seq = [seq]
                            tempseq.degapSeq()
                            tempseq.saveFasta(seqfile='%s.tmp.fas' % acc)
                            seqblast.setInfo({'Name': '%s.tmp.txt' % acc, 'DBase': blastdb, 'Type': 'blastp', 'InFile': '%s.tmp.fas' % acc})
                            if self.info['Type'] in ['DNA','RNA']:
                                seqblast.info['Type'] = 'blastn'
                            if self.opt['NR Align']:
                                seqblast.setStat({'OneLine': blast, 'HitAln': 0})
                                seqblast.blast()
                                if not seqblast.readBLAST(): raise ValueError    ## Seqs now in order of seqblast.search[0].hit
                            else:
                                seqblast.setStat({'OneLine': blast, 'HitAln': blast})   # Need alignments for GABLAMO
                                seqblast.blast()
                                if not seqblast.readBLAST(gablam=True): raise ValueError    ## Seqs now in order of seqblast.search[0].hit
                            os.unlink('%s.tmp.txt' % acc)
                            os.unlink('%s.tmp.fas' % acc)
                            hits = seqblast.search[0].hit
                            hitdic = seqblast.search[0].hitSeq(self)   # Returns dictionary of hits as sequences from seqlist (or None if missing)
                            self.log.opt['Silent'] = silence
                            for hit in hits:
                                otherseq = hitdic[hit]
                                if otherseq == seq:
                                    continue
                                elif otherseq == None:  # Same sequence or deleted sequence
                                    self.verbose(2,4,'%s Hit %d: Deleted' % (seq.shortName(),hits.index(hit)),1)
                                    continue
                                elif samespec and not seq.sameSpec(otherseq):
                                    self.verbose(2,4,'%s Hit %d (%s): Wrong Species' % (seq.shortName(),hits.index(hit),otherseq.shortName()),1)
                                    continue    # Not worth checking as different species
                                elif samespec and self.obj['QuerySeq'] and seq.sameSpec(otherseq) and seq.sameSpec(self.obj['QuerySeq']) and not self.opt['QueryNR']:
                                    self.verbose(2,4,'%s Hit %d (%s): Query Species (querynr=F)' % (seq.shortName(),hits.index(hit),otherseq.shortName()),1)
                                    continue    # Not worth checking as query species

                            ## <iiia> ## Use ALIGN to align sequences -> self.obj['PWAln ID']
                                remredseq = False
                                remredtext = ''
                                if pw_aln and self.opt['NR Align']:
                                    dis1 = self.getDis(seq,otherseq,'PWAln ID')
                                    dis2 = self.getDis(otherseq,seq,'PWAln ID')
                                    self.verbose(1,3,'%s of %d %s Hit %d (%s): PWAln = %.2f%%ID & %.2f%%ID' % (rje.preZero(self.seqs().index(seq)+1,self.seqNum()),self.seqNum(),seq.shortName(),hits.index(hit),otherseq.shortName(),dis1,dis2),1)
                                    if dis1 >= nrid or dis2 >= nrid: ## Remove Sequence
                                        remredtext = '>=%.2f%% (%.2f%%/%.2f%%) ID' % (nrid,dis1,dis2)
                                        remredseq = True
                                    else:
                                        break
                            ## <iiib> ## Use GABLAMO Statistics
                                elif pw_aln:
                                    gdict = hit.dict['GABLAM']
                                    idq = 100.0 * float(gdict['Query']['GABLAMO ID']) / seq.aaLen() #X#gdict['Query']['GABLAMO Len']
                                    simq = 100.0 * float(gdict['Query']['GABLAMO Sim']) / seq.aaLen() #X#gdict['Query']['GABLAMO Len']
                                    idh = 100.0 * float(gdict['Hit']['GABLAMO ID']) / hit.stat['Length']    #X# gdict['Hit']['GABLAMO Len']
                                    simh = 100.0 * float(gdict['Hit']['GABLAMO Sim']) / hit.stat['Length']    #X# gdict['Hit']['GABLAMO Len']
                                    if idq >= nrid or idh >= nrid:
                                        remredtext = '>=%.2f%% (%.2f%%/%.2f%%) ID' % (nrid,idq,idh)
                                        remredseq = True
                                    elif simq >= nrsim or simh >= nrsim:
                                        remredtext = '>=%.2f%% (%.2f%%/%.2f%%) Sim' % (nrsim,simq,simh)
                                        remredseq = True
                                    #self.deBug(gdict)
                                    #self.deBug(idq)
                                    #self.deBug(idh)
                                    
                            ## <iv> ## Based on current alignment
                                else:   ### Use MSA -> self.obj['MSA ID']
                                    continue    #!# Add

                            ## <v> ## Remove Seq ##
                                if remredseq:
                                    rem = self._worstSeq(seq, otherseq, best_ann)
                                    red += 1
                                    blastred += 1
                                    if seq in rem:
                                        if save_red: redseq.seq.append(seq)
                                        if self.opt['NRKeepAnn']: otherseq.info['Name'] += ' >%s' % seq.info['Name']
                                        self.removeSeq(text='%s vs %s: %s' % (remredtext,otherseq.shortName(),rem[1]),seq=seq)
                                        break
                                    else:
                                        if save_red: redseq.seq.append(otherseq)
                                        if self.opt['NRKeepAnn']: seq.info['Name'] += ' >%s' % otherseq.info['Name']
                                        self.removeSeq(text='%s vs %s: %s' % (remredtext,seq.shortName(),rem[1]),seq=otherseq)
                                        if hit == hits[-1] and len(hits) == blast:  # Read all hits - generate more!
                                            fullhits = True     
                                            blast += 10
                progress = (nrid,nrsim,rje.integerString(self.seqNum()),spectxt,rje.integerString(red))
                self.log.opt['Silent'] = silence
                self.printLog('\r#NR','Searched for %.2f%% identity / %.2f%% similarity among %s Sequences%s : %s removed.          ' % progress,log=True,newline=True)

                ## <vi> ## Tidy Up
                if blast > 0 and blastdb: rje_blast.cleanupDB(self,blastdb,True)

        ### <4> ### Finish        
            self.log.opt['Silent'] = silence
            return redseq
        except:
            self.opt['Silent'] = silence
            self.errorLog('Major problem during SeqList.makeNR()',quitchoice=True)
            return       
#########################################################################################################################
    def _worstSeq(self,seq1,seq2,best_ann=False):   ### Returns the sequence that should be removed and why
        '''
        Returns the sequence that should be removed and reason, based on annotation etc.
        >> seq1:Seq Object
        >> seq2:Seq Object
        >> best_ann:boolean [False] = whether to preferentially keep shorter but better annotated sequences.
        << tuple:(Worst Sequence, reason (string), good sequence shortname)
        '''
        try:
            ### <a> ### Keep Master Sequence
            if seq1 == self.obj['QuerySeq']:
                return (seq2,'Not Query Sequence',seq1.shortName())
            elif seq2 == self.obj['QuerySeq']:
                return (seq1,'Not Query Sequence',seq2.shortName())
            ### <b> ### Keep Longer Sequence or better annotated
            len1 = seq1.aaLen()
            len2 = seq2.aaLen()                
            bestdb = self._bestDB(seq1,seq2)
            ## <i> ## Fragments First
            if seq1.info['Name'].lower().find('(fragment)') > 0 and len2 > len1:
                return (seq1,'Shorter Fragment sequence (%daa vs %daa)' % (len1,len2),seq2.shortName())
            elif seq2.info['Name'].lower().find('(fragment)') > 0 and len1 > len2:
                return (seq2,'Shorter Fragment sequence (%daa vs %daa)' % (len2,len1),seq1.shortName())
            ## <ii> ## Annotation Option if 'full length'
            if best_ann == False or bestdb == None:
                if len1 > len2:
                    return (seq2,'Shorter sequence (%daa vs %daa)' % (len2,len1),seq1.shortName())
                elif len2 > len1:
                    return (seq1,'Shorter sequence (%daa vs %daa)' % (len1,len2),seq2.shortName())
            ### <c> ### Keep Best Sequence from Description/Database
            ## <i> ## DBase
            if bestdb == seq1:
                return (seq2,'Kept %s sequence over %s.' % (seq1.info['DBase'],seq2.info['DBase']),seq1.shortName())
            elif bestdb == seq2:
                return (seq1,'Kept %s sequence over %s.' % (seq2.info['DBase'],seq1.info['DBase']),seq2.shortName())
            ## <ii> ## Type of sequence from description
            types = ['genscan','status:novel','Status:known','hypothetical','putative','probable']   # Order of 'badness'
            for type in types:
                if seq1.info['Name'].find(type) > 0 and seq2.info['Name'].find(type) < 0:
                    return (seq1,type,seq2.shortName())
                elif seq2.info['Name'].find(type) > 0 and seq1.info['Name'].find(type) < 0:
                    return (seq2,type,seq1.shortName())
            ### <d> ### Manual Choice
            while self.stat['Interactive'] >= 1:
                print '\n%s and %s seem equally good.' % (seq1.info['Name'],seq2.info['Name'])
                rem = rje.choice('<1> Remove %s, <2> Remove %s, <A>utomatic' % (seq1.shortName(),seq2.shortName()),default='A')
                if rem == '1':
                    return (seq1, 'Manual Choice',seq2.shortName())
                elif rem == '2':
                    return (seq2, 'Manual Choice',seq1.shortName())
                elif rem.lower() == 'a':
                    break
            ### <e> ### #!# Add option to reject less common species
            ### <f> ### Keep earlier sequence in list
            return (seq2,'Arbitrary decision',seq1.shortName())           
        except:
            self.errorLog('Major problem during _worstSeq().',True)
            return (seq2,'Problem with _worstSeq() - arbitrary choice.',seq1.shortName())
#########################################################################################################################
    def _bestDB(self,seq1,seq2,dblist=None):    ### Returns sequence from 'best' database or None if the same
        '''
        Returns sequence from 'best' database or None if the same.
        Used for determining which sequence to remove in the case of redundancy.
        >> seq1:Sequence Object
        >> seq2:Sequence Object
        >> dblist:list of strings = list of DBase annotations, best to worst (lower case)
        << Sequence Object or None
        '''
        try:
            if dblist == None:
                dblist = self.info['DBList']
            dblist = string.split(dblist.lower(),',')
            ### <a> ### Same DBase
            db1 = seq1.info['DBase'].lower()
            db2 = seq2.info['DBase'].lower()
            if db1 == db2:
                return None
            ### <b> ### DBases in list
            if db1 in dblist and db2 in dblist:
                if dblist.index(db1) > dblist.index(db2):   # DB List is good to bad
                    return seq2
                else:
                    return seq1
            elif db1 in dblist:
                return seq1
            elif db2 in dblist:
                return seq2
            else:
                return None
        except:
            self.errorLog('Major problem during _bestDB().',True)
            return None
#########################################################################################################################
    def _checkForDup(self,remdup=True):     ### Checks for and removes Duplicate Sequences
        '''
        Checks for and removes Duplicate Sequences. Checks Name, AccNum and Sequence.
        - renames if not removed or sequence different
        >> remdup:Boolean = whether to remove duplicate sequences
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            dup = 0     # No. duplicate sequences
            ren = 0     # No. renamed sequences
            preseqx = self.seqNum()
            checkseq = self.seqs()[0:]
            seqdict = self.seqNameDic(proglog=self.v()>0)
            accdict = self.seqNameDic('AccNum',proglog=self.v()>0)
            ### ~ [1] ~ Check for Unique names and accnum ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(accdict) == self.seqNum():
                self.printLog('\r#NR','No duplicate AccNum detected.',screen=self.v()>0)
                return False
            ### ~ [2] ~ Duplicate sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            sx = 0.0; stot = self.seqNum()
            for seq in self.seq:
                self.progLog('\r#NR','Removing/renaming duplicate sequences: %.2f%%' % (sx/stot)); sx += 100.0
                if seq in accdict.values(): continue
                acc = seq.info['AccNum']
                otherseq = accdict[acc]
                if seq == otherseq: continue
                remove = seq.info['Sequence'] == otherseq.info['Sequence']  # Same!
                if remove:
                    self.removeSeq('Duplicate Sequence Found!',otherseq)
                    dup += 1
                else:
                    self.printLog('\r#SEQ','WARNING! %s and %s have same Acc#!' % (seq.shortName(),otherseq.shortName()))
                    i = 1
                    while '%s-%d' % (otherseq.shortName(),i) in seqdict or '%s-%d' % (otherseq.info['AccNum'],i) in accdict: i += 1
                    otherseq.info['Name'] = '%s-%d %s' % (otherseq.shortName(),i,otherseq.info['Description'])
                    otherseq.info['AccNum'] = '%s-%d' % (otherseq.info['AccNum'],i)
                    self.printLog('\r#SEQ','Duplicate sequence renamed %s!' % otherseq.shortName())
                    accdict['%s-%d' % (acc,i)] = otherseq
                    seqdict[otherseq.shortName()] = otherseq
                    ren += 1
                accdict[acc] = seq
                seqdict[seq.shortName()] = seq
            self.printLog('\r#NR','Removed/renamed duplicate sequences: %d Removed, %d Renamed.' % (dup,ren))            
            ### ~ [3] ~ Check for duplicate sequence objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(accdict) != self.seqNum():
                sx = 0.0; stot = self.seqNum(); checkseq = self.seqs()[0:]; nrseq = []
                while checkseq:
                    self.progLog('\r#NR','Checking for duplicate sequence objects: %.2f%%' % (sx/stot)); sx += 100.0
                    seq = checkseq.pop(0); nrseq.append(seq)
                    while checkseq.count(seq): checkseq.remove(seq); sx += 100.0
                self.printLog('\r#NR','Checked for duplicate sequence objects: %d removed' % (stot-len(nrseq)))
                self.seq = nrseq[0:]
            self.printLog('\r#NR','%s unique sequence AccNum.' % self.seqNum()); return True
        except:
            self.errorLog('Major problem during _checkForDup().')
            return None
#########################################################################################################################
    ### <8> ### Sequence Manipulation                                                                                   #
#########################################################################################################################
    def MWt(self,sequence): return rje_sequence.MWt(sequence)
#########################################################################################################################
    def winChop(self,seq,win,slide):    ### Chops a sequence up using sliding window and returns SeqList
        """
        Chops a sequence into small sequence windows and returns a new SeqList object
        >> seq:str = sequence
        >> win:int = window size
        >> slide:int = steps to slide window by
        << chopped:SeqList = sequence list of fragments
        """
        ### <1> ### Create new sequence list
        self.verbose(0,2,"Chopping %s (%s)" % (seq.info['Name'], seq.info['Sequence']),1)
        chopped = SeqList(log=self.log)
        chopped.info['Type'] = seq.info['Type']
        ### <2> ### Chop and populate
        r = 0
        while (r+win) <= seq.seqLen():
            start = rje.preZero(r+1,seq.seqLen())
            end = rje.preZero(r+win,seq.seqLen())
            chopname = "%s_%s-%s" % (seq.shortName(), start, end)
            chop = seq.info['Sequence'][r:(r+win)]
            chopped._addSeq(name=chopname,sequence=chop)
            r += slide
        ### <3> ### Hangover
        if r > seq.seqLen():
            r = seq.seqLen() - win
            start = rje.preZero(r+1,seq.seqLen())
            end = rje.preZero(r+win,seq.seqLen())
            chopname = "%s_%s-%s" % (seq.shortName(), start, end)
            chop = seq.info['Sequence'][r:(r+win)]
            chopped._addSeq(name=chopname,sequence=chop)
        ### <4> ### Finish
        self.printLog('#SEQ',"Chopped %s (win=%daa, slide=%daa) => %d Sequences\n" % (seq.shortName(), win, slide, chopped.seqNum()),1)
        return chopped
#########################################################################################################################
    def splitSeq(self,split=4000):  ### Splits seqList into numbered files of X sequences.
        '''
        Splits seqList into numbered files of X sequences.
        >> split:int = Number of sequences per files
        << returns list of output files
        '''
        try:
            ### <a> ### Setup
            self.makeBaseFile()
            basefile = self.info['Basefile']
            split = int(split)
            seqnum = self.seqNum()
            ### <b> ### Split
            if seqnum <= split:
                self.verbose(0,2,'No more than %d sequences - no need to split!' % split,1)
                return
            self.verbose(0,2,'Splitting %s sequences into files of %s sequences...' % (rje.integerString(seqnum),rje.integerString(split)),1)
            s = [0,split]
            filelist = []
            splitseq = copy.deepcopy(self)
            while s[0] < seqnum:
                if s[1] > seqnum:
                    s[1] = seqnum
                splitseq.info['Name'] = '%s.%s-%s.fas' % (basefile,rje.preZero(s[0]+1,seqnum),rje.preZero(s[1],seqnum))
                splitseq.seq = self.seqs()[s[0]:s[1]]
                splitseq.saveFasta()
                s[0] = s[1]
                s[1] += split
                filelist.append(splitseq.info['Name'])
            self.printLog('#SEQ','%s split into %s files of <=%s sequences.' % (self.info['Name'],rje.integerString(len(filelist)),rje.integerString(split)))
            return filelist
        except:
            self.errorLog('Major Problem with splitSeq().',True)
#########################################################################################################################
    def truncSeq(self,trunc=0):  ### Truncates sequences in seqList to first X AAs.
        '''
        Splits seqList into numbered files of X sequences.
        >> trunc:int = Length of truncated sequences. (-ve values with return last X AAs)
        '''
        try:
            for seq in self.seqs():    
                if trunc > 0: seq.info['Sequence'] = seq.info['Sequence'][:trunc]
                else: seq.info['Sequence'] = seq.info['Sequence'][trunc:]
        except: self.errorLog('Major Problem with truncSeq().',True)
#########################################################################################################################
    def tidyQueryGaps(self,key='Sequence',backup=''):   ### Removes columns that are gapped in the query.
        if not self.obj['QuerySeq'] and not self.querySeq(): return self.errorLog('Cannot tidyQueryGaps without Query!',printerror=False)
        self.printLog('#DEGAP','Removing alignment columns with gaps in %s' % self.obj['QuerySeq'].shortName())
        return self.tidyGaps(key,[self.obj['QuerySeq']],backup)
#########################################################################################################################
    def tidyGaps(self,key='Sequence',seqs=[],backup=''):  ### Removes 100% gap columns from the alignment
        '''
        Removes 100% gap columns from the alignment.
        >> key:str ['Sequence'] = seq.info Key to tidy. (May be MaskSeq etc.)
        >> seqs:list [] = List of seqs used to judge gappiness. e.g. Could be just Query. If empty, will use all.
        >> backup:str [''] = If given, will copy full length info to seq.info[backup]
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.list['GapBorder'] = []
            if self.seqNum() < 1 or not self._checkAln(aln=True,realign=False,tidygaps=False): return
            if not seqs: seqs = self.seqs()[0:]
            if backup:
                for seq in self.seqs():
                    seq.info[backup] = seq.info[key][0:]
                    #self.deBug(rje.sortKeys(seq.info))
            ### ~ [1] ~ Select 100% gap columns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gapcol = []
            for r in range(self.seqs()[0].seqLen()):
                allgap = True
                for seq in seqs:
                    if seq.info[key][r] != '-':
                        allgap = False
                        break
                if allgap: gapcol.append(r)
            ### ~ [2] ~ Remove gapcol ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gapcol.reverse(); gapborder = []
            for col in gapcol:
                for seq in self.seqs():
                    seq.info[key] = seq.info[key][:col] + seq.info[key][(col+1):]
                seq = self.seqs()[0]
                g = len(seq.info[key][(col+1):])
                if g not in gapborder: gapborder.append(g)
                if (g+1) not in gapborder: gapborder.append(g+1)
            for g in gapborder:
                if g <= self.seqLen(): self.list['GapBorder'].append(self.seqLen()-g)
            self.list['GapBorder'].sort()
            if gapcol and (self.v() > 0 or backup): self.printLog('#DEGAP','%d gapped columns removed.' % len(gapcol))
        except:
            self.errorLog('Major Problem with tidyGaps(%s)' % key)
            raise
#########################################################################################################################
    def tidyXGaps(self,key='Sequence'):  ### Removes 100% gap/"X" columns from the alignment
        '''Removes 100% gap/"X" columns from the alignment.'''
        try:
            gapcolx = self.stripGap(self.seqNum(),gaps=['-','X'],seqkey=key)
            return self.verbose(1,3,'%d 100%% X-gapped columns removed.' % gapcolx,1)
        except:
            self.errorLog('Major Problem with tidyGaps(%s)' % key)
            raise
#########################################################################################################################
    def stripGap(self,stripgap=0,codons=False,backup='',gaps=['-'],seqkey='Sequence'): ### Removes columns containing gaps from the alignment
        '''
        Removes columns containing gaps from the alignment.
        >> stripgap:num [0] = Number of sequences with gaps before stripping from alignment. Proportion if < 1.
        >> codons:bool [False] = Whether to treat alignment using a codon model (i.e. strip sets of three bases)
        >> backup:str [''] = Backup full-length sequences in seq.info[backup]
        >> gaps:list ['-'] = Characters to recognise as gaps (e.g. could add 'X' for tidyXGap equivalent)
        >> seqkey:str ['Sequence'] = seq.info key to use to check alignment
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if stripgap < 0: return    # No gap stripping
            if not self.opt['Aligned'] and not self.checkAln(aln=True,tidygaps=False,seqkey=seqkey):
                self.errorLog('SeqList.stripGap() called but sequences not aligned!',printerror=False)
                raise ValueError
            if stripgap < 1: stripgap = max(1,int(stripgap*self.seqNum()+0.5))  # Must need at least one seq with gap!
            if backup:
                for seq in self.seqs()(): seq.info[backup] = seq.info[seqkey]
            gaplist = [0] * self.seqLen(seqkey=seqkey)
            ### ~ [2] Count gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                for r in range(seq.seqLen()):
                    if seq.info[seqkey][r] in gaps: gaplist[r] += 1
            ## ~ [2a] Adjust codon counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if codons:
                cnum = len(gaplist) / 3
                if len(gaplist) / 3 != len(gaplist) / 3.0: self.errorLog('Warning! Codons stripGap model invoked but sequence length not multiple of 3!', printerror=False); cnum += 1
                for cx in range(cnum):
                    ci = cx * 3
                    cgap = max(gaplist[ci:ci+3])
                    for i in range(3):
                        try: gaplist[ci+i] = cgap
                        except: pass
            ### ~ [3] Remake sequences without gappy columns ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                oldseq = seq.info[seqkey]
                seq.info[seqkey] = ''
                for r in range(len(gaplist)):
                    if gaplist[r] < stripgap: seq.info[seqkey] += oldseq[r]; gaplist[r] = 0
            stripx = len(gaplist) - gaplist.count(0)
            if codons and stripx: self.printLog('#GAPS','Stripped %s gapped codons (>%s gapped sequences)' % (stripx/3,stripgap))
            elif stripx: self.printLog('#GAPS','Stripped %s gapped positions (>%s gapped sequences)' % (stripx,stripgap))
            return stripx
        except ValueError: raise
        except: self.errorLog('Major problem with SeqList.stripGap()'); return -1
#########################################################################################################################
    def degapSeq(self,log=True):  ### Removes gaps from sequences
        '''
        Removes gaps from sequences.
        >> log:boolean = Whether to report degapping in log.
        '''
        try:### ~ [1] ~ Remove all gap characters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            gx = 0
            for seq in self.seqs():
                gx += seq.info['Sequence'].count('-')
                seq.info['Sequence'] = string.join(string.split(seq.info['Sequence'],'-'),'')
            self.opt['Gapped'] = self.opt['Aligned'] = False
            if log and gx: self.printLog('#GAP','%s gaps removed from %s sequences.' % (rje.iStr(gx),rje.iStr(self.seqNum())))
        except: self.errorLog('Major Problem with degapSeq().',True)
#########################################################################################################################
    def removeSeq(self,text='',seq=None,checkAln=False):     ### Removes a sequence and logs reason
        '''
        Removes a sequence and logs the reason.
        >> text:str = Reason for removal
        >> seq:Sequence object
        >> checkAln:Boolean = whether to CheckAln (includes TidyGap) after seq removal
        '''
        try:
            if text != '':
                text = ': ' + text
                if text[-1] != '.': text = '%s.' % text
            if seq == None:
                self.errorLog('Non-existent sequence object called for removal%s' % text)
                if seq in self.seqs():     # Sequence to remove
                    self.printLog('\r#EEK','Deleted Null sequence object')
                    self.seq.remove(seq)
            elif seq in self.seqs():     # Sequence to remove
                self.printLog('\r#REM','Deleted %s%s   \n' % (seq.shortName(),text),1)
                self.seq.remove(seq)
            else: self.errorLog('Sequence %s called for removal but not in SeqList%s' % (seq.shortName(),text))
            if checkAln: self._checkAln(aln=self.opt['Aligned'])
        except: self.errorLog('Major Problem with removeSeq. No sequence removed.')
#########################################################################################################################
    def mapX(self,query=None,seqlist=[],skipgaps=True,qtrim=False,focus=[0,0]):  ### Maps Xs from query sequence onto other sequences
        '''
        Maps Xs from query sequence onto other sequences.
        >> query:Sequence Object from which to take Xs (Returns if None) [None]
        >> seqlist:list of Sequence objects to map Xs onto (uses self.seq if []) []
        >> skipgaps:boolean = Whether to skip gaps (True) or replace gaps with Xs (False) [True]
        >> qtrim:boolean= Whether to replace residues outside of Query with Xs [False]
        >> focus:list of range positions [X:Y] to look at. If Y=0 then [X:]. Outside will be Xd.
        '''
        try:
            _stage = '<0> Setup'
            if seqlist == []:
                seqlist = self.seqs()
            if focus[1] < 1:
                focus[1] = seqlist[0].seqLen() + focus[1]

            _stage = '<1> Query'
            if query and qtrim:
                [qstart, qend] = [0, query.seqLen()]
                r = 0
                while r < query.seqLen():
                    if query.info['Sequence'][r] != '-':
                        qstart = r
                        break
                    r += 1
                r = 0
                while r < query.seqLen():
                    r += 1                    
                    if query.info['Sequence'][-r] != '-':
                        qend = query.seqLen() - r 
                        break
                if focus[0] < qstart:
                    focus[0] = qstart
                if focus[1] > qend:
                    focus[1] = qend 

            _stage = '<2> Map Xs'
            for r in range(seqlist[0].seqLen()):
                if (query and query.info['Sequence'][r] == 'X') or r < focus[0] or r > focus[1]:
                    for seq in seqlist:
                        a = seq.info['Sequence'][r]
                        if (a != 'X') and (a != '-' or skipgaps == False):
                            seq.info['Sequence'] = rje.strSub(seq.info['Sequence'],r,r,'X')
        except:
            self.errorLog('Major problem with mapX(%s)' % _stage)
            raise
#########################################################################################################################
    def sortByLen(self,longfirst=True,proglog=False):    ### Sort sequences according to length
        '''Sort sequences according to length.'''
        ### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        lendict = {}
        neworder = []
        (sx,stot) = (0.0,self.seqNum())
        ### ~ [2] ~ Sort ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for seq in self.seqs():
            if proglog: self.progLog('\r#SORT','Sorting %s sequences by length: %.2f%%' % (rje.integerString(stot),sx/stot)); sx += 50.0
            if seq.nonN() not in lendict: lendict[seq.nonN()] = []
            lendict[seq.nonN()].append(seq)
        for nlen in rje.sortKeys(lendict,revsort=longfirst):
            if proglog: self.progLog('\r#SORT','Sorting %s sequences by length: %.2f%%' % (rje.integerString(stot),sx/stot)); sx += len(lendict[nlen])*50.0
            neworder += lendict[nlen]
        self.seq = neworder
        if proglog: self.printLog('\r#SORT','Sorting %s sequences by length complete' % rje.integerString(stot),log=False)
#########################################################################################################################
    ### <9> ### Sequence Identity and Distance Matrices etc.                                                            #
#########################################################################################################################
    def getDis(self,seq1,seq2,key=None,unlink=True):  ### Returns/generates appropriate distance
        '''
        Returns/generates appropriate distance.
        >> seq1 & seq2:Sequence Objects
        >> key:str = key to self.obj list (DisMatrix object)
        >> unlink:Boolean = whether to delete tmp files.
        << float = Distance
        '''
        try:### ~ Setup: Check sequences and matrices ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seq1 or not seq2: return 0.0
            if key not in self.obj.keys():
                self.errorLog('Key (%s) called for getDis() but no DisMatrix Object.' % key,True)
                raise ValueError
            if key.find('MSA') == 0 and seq1.seqLen() != seq2.seqLen():
                self.errorLog('Key (%s) called for getDis() but sequences and are different lengths! %s:%daa and %s:%daa. (Not aligned?!)' % (key,seq1.shortName(),seq1.seqLen(),seq2.shortName(),seq2.seqLen()),True)
                raise ValueError
            ### ~ Look for existing distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj[key]:
                dis = self.obj[key].getDis(seq1,seq2)
                if dis != None: return dis
            ### ~ Else make/calculate distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ PWAln Identity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if key == 'PWAln ID':
                pwaln = self.pwAln(seq1,seq2,unlink)
                self.obj[key].addDis(seq1,seq2,pwaln.stat['Identity'] * pwaln.stat['Length'] / pwaln.stat['QryEnd'])
                self.obj[key].addDis(seq2,seq1,pwaln.stat['Identity'] * pwaln.stat['Length'] / pwaln.stat['SbjEnd'])
            ## ~ %Identity, Gaps and Extra of seq1 vs seq2 based on MSA ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            elif key.find('MSA') == 0:
                for msa in ['MSA ID','MSA Gaps','MSA Extra']:
                    if self.obj[msa] == None:
                        self.addMatrix(msa)
                seq = [seq1.info['Sequence'],seq2.info['Sequence']]
                id = [0,0]
                gaps = [0,0]
                extra = [0,0]
                seqlen = (seq1.aaLen(),seq2.aaLen())
                if seq[0] == seq[1]: id[0:] = seqlen[0:]
                else:
                    for r in range(len(seq[0])):
                        aa = (seq[0][r],seq[1][r])
                        if aa[0] == aa[1]:  # Matching residue
                            if aa[0] != '-' and aa[0] != 'X':    # Not a gap or X
                                id[0] += 1
                                id[1] += 1
                        elif aa[0] == '-':  # Gap in seq1 but not seq2
                            gaps[0] += 1
                            extra[1] += 1
                        elif aa[1] == '-':  # Gap in seq2 but not seq1
                            gaps[1] += 1
                            extra[0] += 1
                self.obj['MSA ID'].addDis(seq1,seq2,100.0 * float(id[0]) / seqlen[0])
                self.obj['MSA ID'].addDis(seq2,seq1,100.0 * float(id[1]) / seqlen[1])
                self.obj['MSA Gaps'].addDis(seq1,seq2,100.0 * float(gaps[0]) / seqlen[1])
                self.obj['MSA Gaps'].addDis(seq2,seq1,100.0 * float(gaps[1]) / seqlen[0])
                self.obj['MSA Extra'].addDis(seq1,seq2,100.0 * float(extra[0]) / seqlen[0])
                self.obj['MSA Extra'].addDis(seq2,seq1,100.0 * float(extra[1]) / seqlen[1])
            ## ~ No other methods currently supported ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else: return 0.0
            return self.obj[key].getDis(seq1,seq2,default=0.0)      
        except SystemExit: raise
        except KeyboardInterrupt: raise
        except:
            self.errorLog('Major problem during getDis(%s vs %s. %s).' % (seq1.shortName(),seq2.shortName(),key))
            return 0.0
#########################################################################################################################
    def bestMeanID(self,queries,seqlist,key=None):     ### Returns sequence with best mean %ID to other sequences
        '''
        >> queries:list of Sequence Objects to compare
        >> seqlist:list of Sequence Objects used for comparison
        >> key:str = key to self.obj list (DisMatrix object)
        << bestseq = Sequence Object
        '''
        try:
            bestseq = queries[0]
            bestid = 0.0
            for seq in queries:
                meanid = 0.0
                for otherseq in seqlist:
                    meanid += self.getDis(seq,otherseq,key=key)
                if meanid > bestid:
                    bestid = meanid
                    bestseq = seq
            return bestseq                    
        except:
            self.errorLog('Major Problem with bestMeanID(%s). Possibly arbitrary best sequence returned.' % key)
            return bestseq
#########################################################################################################################
    ### <10> ### Alignments                                                                                             #
#########################################################################################################################
    def pwAln(self,seq1,seq2,unlink=True,retry=5):  ### Uses align to make a pairwise alignment of seq1 and seq2
        '''
        Uses align to make a pairwise alignment of seq1 and seq2.
        >> seq1 & seq2:Sequence Objects
        >> unlink:Boolean = whether to delete tmp files.
        >> retry:int = Descrease with each retry: if retry = 0, give up!
        << PWAln object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Sequence check ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not (seq1 and seq2):
                self.errorLog('Sequence Missing!',printerror=False)
                pwaln.stat['Identity'] = 0.0
                pwaln.stat['Length'] = 100.0
                pwaln.stat['QryEnd'] = pwaln.stat['Length']
                pwaln.stat['SbjEnd'] = pwaln.stat['Length']
                return pwaln
            elif seq1 == seq2 or seq1.info['Sequence'] == seq2.info['Sequence']: # Same Sequence - no need for align!
                pwaln = rje_blast_V1.PWAln(log=self.log)
                pwaln.info['Name'] = '%s Self' % seq1.shortName()
                pwaln.stat['Identity'] = 100.0
                pwaln.stat['Length'] = seq1.aaLen()
                pwaln.stat['QryEnd'] = pwaln.stat['Length']
                pwaln.stat['SbjEnd'] = pwaln.stat['Length']
                return pwaln
            ## ~ [0b] ~ ALIGN setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pwaln = rje_blast_V1.PWAln(log=self.log)
            acc1 = '%s' % seq1.info['AccNum']
            acc2 = '%s' % seq2.info['AccNum']
            if unlink:  # Random numbers are desirable to help avoid accidental deletion of files by parallel forks
                acc1 = '%s.%s' % (acc1,rje.randomString(6))
                acc2 = '%s.%s' % (acc2,rje.randomString(6))
            acc1 += '.tmp.fas'
            acc2 += '.tmp.fas'
            command = '%align %s %s' % (self.info['FASTA Path'],acc1,acc2)

            ### ~ [1] ~ Align Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###        
            ## ~ [1a] ~ Save sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            SeqList(log=self.log,cmd_list=['v=-1']).saveFasta(seqs=[seq1],seqfile=acc1,log=False)
            SeqList(log=self.log,cmd_list=['v=-1']).saveFasta(seqs=[seq2],seqfile=acc2,log=False)
            ## ~ [1b] ~ Call Align ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.verbose(2,3,command,1)
            if rje.checkForFile(acc1) and rje.checkForFile(acc2):
                lines = os.popen(command).readlines()
                if unlink:
                    os.unlink(acc1)
                    os.unlink(acc2)
            else: self.errorLog('File %s or %s missing!' % (acc1,acc2),printerror=False); raise IOError

            ### ~ [3] ~ Read Results into PWAln Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            resline = []
            for line in lines: resline.append(rje.chomp(line))
            ## ~ [3a] ~ Process Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            pwaln.setStat({'QryStart': 1, 'QryEnd': seq1.aaLen(), 'SbjStart': 1, 'SbjEnd': seq2.aaLen()})
            i = 0
            while i < len(resline):
                line = resline[i]
                if line.find('Global alignment score:') >= 0:
                    stats = rje.matchExp('(\S+)\% identity;\s+Global alignment score:\s+(\-?\d+)',line)
                    pwaln.setStat({'Identity': string.atof(stats[0]), 'BitScore': string.atof(stats[1])})
                    break
                i += 1
            ## ~ [3b] ~ Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            leader = 0
            portion = 0
            i += 1
            while i < len(resline):
                line = resline[i]
                if re.search('^(\S+\s+)(\S+)',line):    # Qry
                    qaln = rje.matchExp('^(\S+\s+)(\S+)',line)
                    leader = len(qaln[0])
                    portion = len(qaln[1])
                    pwaln.info['QrySeq'] += qaln[1]
                    i += 1
                    line = resline[i]
                    pwaln.info['AlnSeq'] += line[leader:(leader+portion)]
                    i += 1
                    line = resline[i]
                    pwaln.info['SbjSeq'] += line[leader:(leader+portion)]
                i += 1
            pwaln.stat['Length'] = len(pwaln.info['QrySeq'])
            if pwaln.stat['Length'] == 0:
                self.errorLog('Major problem with Align - no length!.')
                raise ValueError
            ## ~ [3c] ~ Return Object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            return pwaln
        except KeyboardInterrupt: raise
        except:
            if seq1: errseq = seq1.shortName()
            else: errseq = 'None'
            if seq2: errseq += ' v %s' % seq2.shortName()
            else: errseq = ' v None'
            if command: print 'Align command: %s' % command
            if retry > 0:
                self.errorLog('Problem during pwAln(%s). Trying again... (%d)' % (errseq,retry))
                return self.pwAln(seq1,seq2,unlink=unlink,retry=(retry-1))
            else: self.errorLog('Major problem during pwAln(%s)' % (errseq)); raise
#########################################################################################################################
    def Align2Seq(self,sequence1,sequence2,unlink=True,retry=5):  ### Uses align to make a pairwise alignment of seq1 and seq2
        '''
        Uses align to make a pairwise alignment of seq1 and seq2.
        >> sequence1 & sequence2:string Sequences to align
        >> unlink:Boolean = whether to delete tmp files.
        >> retry:int = Descrease with each retry: if retry = 0, give up!
        << PWAln object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            tempseq = SeqList(log=self.log,cmd_list=['v=-1'])
            seq1 = tempseq._addSeq('tmp%s' % rje.randomString(8),sequence1)
            seq2 = tempseq._addSeq('tmp%s' % rje.randomString(8),sequence2)
            return self.pwAln(seq1,seq2,unlink,retry)
        except: self.errorLog('Problem during Align2Seq(). Should be using pwAln()?'); raise
#########################################################################################################################
    def align(self,alnprog=None,outfile=None,mapseq=True,alncommands='',log=True,clustalbackup=True):  ### Uses program of choice to align self and return object
        '''
        Uses program of choicle to align self and returns aligned sequence object.
        >> alnprog:str = Alignment program to use (self.info['AlnProg'] if None)
        >> outfile:str = name for output file. [self.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        >> log:bool [True] = Whether to report alignment activity to Log
        >> clustalbackup [True] = Whether to use ClustalO and ClustalW as backup in case of failure
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [0a] ~ Check sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.seqNum() < 2:   # Cannot align <2 sequences!
                self.printLog('#ALN','Less than 2 sequences - cannot align!')
                self.degapSeq(log=False)
                self.opt['Gapped'] = self.opt['Aligned'] = True
                self.obj['MSA Gaps'] = self.obj['MSA ID'] = self.obj['MSA Extra'] = None
                if outfile: self.saveFasta(seqfile=outfile)                
                return self
            ## ~ [0b] ~ Select alignment program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not alnprog: alnprog = self.info['AlnProg']
            if alnprog.lower() not in ['muscle','maaft','mafft','fsa','clustalw','clustalo','pagan']:
                self.errorLog('Alignment program "%s" not recognised - using ClustalW' % alnprog,printerror=False)
                alnprog = 'clustalw'
            #self.deBug('%s -> %s' % (self.info['AlnProg'],alnprog))
            ## ~ [0c] ~ Save sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.degapSeq(log=False)
            infile = self.info['Name']
            if infile.lower() in ['none','']: infile = 'tmp%s.fas' % rje.randomString(5)
            if mapseq and alnprog.lower() in ['clustalw']:
                infile = 'tmp%s.fas' % rje.randomString(5)
                self.saveFasta(seqfile=infile,name='AccNum')
            else: self.saveFasta(seqfile=infile)

            ### ~ [1] ~ Call alignment program ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alnprog.lower() == 'muscle': alnseq = self.alignWithMUSCLE(infile,outfile,mapseq,alncommands)
            elif alnprog.lower() in ['maaft','mafft']: alnseq = self.alignWithMAFFT(infile,outfile,mapseq,alncommands)
            elif alnprog.lower() == 'fsa': alnseq = self.alignWithFSA(infile,outfile,mapseq,alncommands)
            elif alnprog.lower() == 'pagan': alnseq = self.alignWithPAGAN(infile,outfile,mapseq,alncommands)
            else: alnseq = self.alignWithClustal(infile,outfile,mapseq,alncommands,alnprog=alnprog)
            ## ~ [1a] ~ Cleanup and check alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##            
            if infile not in [self.info['Name'],outfile] and os.path.exists(infile): os.unlink(infile)
            if not (alnseq and alnseq._checkAln(aln=True,realign=False)):
                if clustalbackup and alnprog.lower() not in ['clustalo','clustalw']:
                    self.printLog('#ALN','Align with %s failed. Trying ClustalO.' % alnprog)
                    return self.align('clustalo',outfile,mapseq,'',log,clustalbackup=False)
                elif clustalbackup and alnprog.lower() != 'clustalw':
                    self.printLog('#ALN','Align with %s failed. Trying ClustalW.' % alnprog)
                    return self.align('clustalw',outfile,mapseq,'',log,clustalbackup=False)
                else: self.printLog('#ALN','Align with %s failed.' % alnprog); raise ValueError
            ### ~ [2] ~ Map sequences if appropriate and return ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if mapseq:
                if self.mapSeq(alnseq): self.printLog('#MAP','Aligned sequences successfully mapped onto originals')
                self.opt['Gapped'] = self.opt['Aligned'] = True
                self.obj['MSA Gaps'] = None
                self.obj['MSA ID'] = None
                self.obj['MSA Extra'] = None
                return self
            else:
                alnseq.setOpt(self.opt)
                alnseq.opt['Gapped'] = alnseq.opt['Aligned'] = True
                return alnseq
        except: self.errorLog('Problem during align(%s).' % (alnprog)); self.degapSeq(log=False); return None
#########################################################################################################################
    def alignWithMUSCLE(self,infile,outfile=None,mapseq=True,alncommands='',log=True):   ### Uses muscle to align infile and returns aligned sequence object
        '''
        Uses muscle to align self and returns aligned sequence object.
        >> infile:str = name of input file. 
        >> outfile:str = name for output file. [input.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aln_out = '%s.muscle' % rje.baseFile(infile)
            ### ~ [1] ~ Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mpath = self.info['MUSCLE']
            if self.opt['Win32']:
                alncommands = '-maxmb 128 %s' % alncommands
                command = mpath + ' -in "%s" -out "%s" %s' % (infile,aln_out,alncommands)
            else: command = mpath + ' -in "%s" -out "%s" %s' % (infile,aln_out,alncommands)
            if self.stat['Verbose'] < 0:
                command = command + ' -quiet'
                if log: self.printLog('#ALN','MUSCLE alignment: %s' % command)
                os.popen(command).readlines()
            else:
                if log: self.printLog('#ALN','MUSCLE alignment: %s' % command)
                os.system(command)
            ### ~ [2] ~ New SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cmd = self.cmd_list + ['seqin=%s' % aln_out,'autoload=T','degap=F','autofilter=F','seqnr=F','unkspec=T']
            if mapseq: cmd += ['v=-1']
            alnseq = SeqList(log=self.log,cmd_list=cmd)
            if os.path.exists(aln_out):
                if outfile: os.rename(aln_out, outfile)
                elif aln_out != outfile and os.path.exists(aln_out): os.unlink(aln_out)
            return alnseq
        except: self.errorLog('Problem during alignWithMuscle(%s).' % (self.info['Name'])); return None
#########################################################################################################################
    def alignWithPAGAN(self,infile,outfile=None,mapseq=True,alncommands='',log=True):   ### Uses pagan to align infile and returns aligned sequence object
        '''
        Uses PAGAN to align self and returns aligned sequence object.
        >> infile:str = name of input file. 
        >> outfile:str = name for output file. [input.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aln_out = '%s.pagan' % rje.baseFile(infile)
            ### ~ [1] ~ Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mpath = self.info['PAGAN']
            command = mpath + ' --seqfile "%s" --outfile "%s" %s' % (infile,aln_out,alncommands)
            if self.stat['Verbose'] < 0:
                command = command + ' --silent'
                if log: self.printLog('#ALN','PAGAN alignment: %s' % command)
                os.popen(command).readlines()
            else:
                if log: self.printLog('#ALN','PAGAN alignment: %s' % command)
                os.system(command)
            ### ~ [2] ~ New SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aln_fas = '%s.fas' % aln_out
            aln_tre = '%s.tre' % aln_out
            cmd = self.cmd_list + ['seqin=%s' % aln_fas,'autoload=T','degap=F','autofilter=F','seqnr=F','unkspec=T']
            if mapseq: cmd += ['v=-1']
            alnseq = SeqList(log=self.log,cmd_list=cmd)
            if os.path.exists(aln_fas):
                if outfile: os.rename(aln_fas, outfile)
                elif aln_fas != outfile and os.path.exists(aln_fas): os.unlink(aln_fas)
            if os.path.exists(aln_tre): os.unlink(aln_tre)
            return alnseq
        except: self.errorLog('Problem during alignWithMuscle(%s).' % (self.info['Name'])); return None
#########################################################################################################################
    def alignWithMAFFT(self,infile,outfile=None,mapseq=True,alncommands='',log=True):   ### Uses muscle to align infile and returns aligned sequence object
        '''
        Uses MAFFT to align self and returns aligned sequence object.
        >> infile:str = name of input file. 
        >> outfile:str = name for output file. [input.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aln_out = '%s.mafft' % rje.baseFile(infile)
            ### ~ [1] ~ Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mpath = self.info['MAFFT']
            command = mpath + ' --auto "%s" > "%s"' % (infile,aln_out)
            if log: self.printLog('#ALN','MAFFT alignment: %s' % command)
            if self.stat['Verbose'] < 0: os.popen(command).readlines()
            else: os.system(command)
            ### ~ [2] ~ New SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cmd = self.cmd_list + ['seqin=%s' % aln_out,'autoload=T','degap=F','autofilter=F','seqnr=F','unkspec=T']
            if mapseq: cmd += ['v=-1']
            alnseq = SeqList(log=self.log,cmd_list=cmd)
            if os.path.exists(aln_out):
                if outfile: os.rename(aln_out, outfile)
                elif aln_out != outfile and os.path.exists(aln_out): os.unlink(aln_out)
            return alnseq
        except: self.errorLog('Problem during alignWithMAFFT(%s).' % (self.info['Name'])); return None
#########################################################################################################################
    def alignWithFSA(self,infile,outfile=None,mapseq=True,alncommands='',log=True):   ### Uses muscle to align infile and returns aligned sequence object
        '''
        Uses FSA to align self and returns aligned sequence object.
        >> infile:str = name of input file. 
        >> outfile:str = name for output file. [input.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            aln_out = '%s.fsa_aln' % rje.baseFile(infile)
            ### ~ [1] ~ Alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            mpath = self.info['FSA']
            command = mpath + ' "%s" > "%s"' % (infile,aln_out)
            if log: self.printLog('#ALN','FSA alignment: %s' % command)
            if self.stat['Verbose'] < 0: os.popen(command).readlines()
            else: os.system(command)
            ### ~ [2] ~ New SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cmd = self.cmd_list + ['seqin=%s' % aln_out,'autoload=T','degap=F','autofilter=F','seqnr=F','unkspec=T']
            if mapseq: cmd += ['v=-1']
            alnseq = SeqList(log=self.log,cmd_list=cmd)
            if os.path.exists(aln_out):
                if outfile: os.rename(aln_out, outfile)
                elif aln_out != outfile and os.path.exists(aln_out): os.unlink(aln_out)
            return alnseq
        except: self.errorLog('Problem during alignWithFSA(%s).' % (self.info['Name'])); return None
#########################################################################################################################
    def alignWithClustal(self,infile,outfile=None,mapseq=True,alncommands='',log=True,alnprog='clustalw'):    ### Uses clustalW/O to align self and returns aligned sequence object
        '''
        Uses clustalw or clustalo to align self and returns aligned sequence object.
        >> infile:str = name of input file. 
        >> outfile:str = name for output file. [base.aln] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        >> log:bool = Whether to report progress to log
        >> alngprog:str = Whether to use ClustalW or ClustalO
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            base = infile
            if base[-4:] == '.fas': base = base[:-4]
            aln_out = '%s.aln' % base   #rje.baseFile(infile)
            ### ~ [1] ~ Alignment with Clustalw ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if alnprog == 'clustalo':
                mpath = self.info['ClustalO']
                command = mpath + ' -i %s -o %s --outfmt=clu %s' % (infile,aln_out,alncommands)
            else:
                mpath = self.info['ClustalW']
                command = mpath + ' %s %s' % (infile,alncommands)
            if log and alnprog == 'clustalo': self.printLog('#ALN','Clustal Omega alignment: %s' % command)
            elif log: self.printLog('#ALN','ClustalW alignment: %s' % command)
            if self.stat['Verbose'] < 0: os.popen(command).readlines()
            else: os.system(command)
            ### ~ [2] ~ New SeqList object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cmd = self.cmd_list + ['seqin=%s' % aln_out,'autoload=T','degap=F','autofilter=F','seqnr=F','unkspec=T']
            if mapseq: cmd += ['v=-1']
            alnseq = SeqList(log=self.log,cmd_list=cmd)
            if os.path.exists('%s.dnd' % base): os.unlink('%s.dnd' % base)
            if os.path.exists(aln_out):
                if outfile: os.rename(aln_out, outfile)
                elif aln_out != outfile and os.path.exists(aln_out): os.unlink(aln_out)
            return alnseq
        except: self.errorLog('Problem during alignWithClustal(%s).' % (self.info['Name'])); return None
#########################################################################################################################
    def clustalAln(self,outfile=None,mapseq=True,alncommands=''):    ### Uses clustalW to align self and returns aligned sequence object
        '''
        Uses clustalw to align self and returns aligned sequence object.
        >> outfile:str = name for output file. [base.aln] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.align(alnprog='clustalw',outfile=outfile,mapseq=mapseq,alncommands=alncommands)
        except:
            self.errorLog('Problem during clustalAln(%s). (May be name problem if not *.fas?)' % (self.info['Name']))
            return None
#########################################################################################################################
    def muscleAln(self,outfile=None,mapseq=True,alncommands=''):    ### Uses muscle to align self and returns aligned sequence object
        '''
        Uses muscle to align self and returns aligned sequence object.
        >> outfile:str = name for output file. [self.muscle] if None
        >> mapseq:boolean = whether to map sequences back on to (and return) self [True]
        >> alncommands:str = additional commandline commands
        << SeqList object
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.align(outfile=outfile,mapseq=mapseq,alncommands=alncommands)    # Will use self.info['AlnProg']
        except:
            self.errorLog('Problem during muscleAln(%s).' % (self.info['Name']))
            return None
#########################################################################################################################
    def relCons(self,outfile=None):  ### Calculates and outputs alignment relative conservation based on Shannon Entropy
        '''
        Calculates and outputs alignment relative conservation based on Shannon Entropy.
        >> outfile:str [None] = Name of file to use for output in place of self.info['RelCons']
        '''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            slimcalc = rje_slimcalc.SLiMCalc(self.log,self.cmd_list+['conspec=','slimcalc='])
            slimcalc.loadBLOSUM()
            if not outfile: outfile = self.info['RelCons']
            headers = ['Pos','Abs','Rel','AbsVNE','RelVNE']
            if self.obj['QuerySeq']:
                self.tidyQueryGaps(backup='PreStrip')
                qryseq = self.obj['QuerySeq'].info['Sequence']
                self.obj['QuerySeq'].disorder()
                self.obj['QuerySeq'].info['Sequence'] = qryseq
                headers.insert(1,self.obj['QuerySeq'].shortName())
                headers += ['IUPred','Disorder']
            else: qryseq = ''
            rje.delimitedFileOutput(self,outfile,headers,rje_backup=True)
            ### ~ [2] ~ Calculate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            abs = self.relConList(0)
            rel = self.relWinCon(abs,self.stat['RelConWin'])
            abv = self.relConList(0,slimcalc)
            vne = self.relWinCon(abv,self.stat['RelConWin'],self.obj['QuerySeq'])
            q = 0
            for i in range(len(abs)):
                if not qryseq: rje.delimitedFileOutput(self,outfile,headers,datadict={'Pos':i+1,'Abs':'%.5f' % abs[i],'Rel':'%.5f' % rel[i],'AbsVNE':'%.5f' % abv[i],'RelVNE':'%.5f' % vne[i]})
                elif qryseq[i] != '-':
                    q += 1
                    datadict={'Pos':q,'Abs':'%.5f' % abs[i],'Rel':'%.5f' % rel[i],'AbsVNE':'%.5f' % abv[i],'RelVNE':'%.5f' % vne[i]}
                    try:
                        datadict['Disorder'] = {True:'Disorder',False:'Order'}[self.obj['QuerySeq'].isDisordered(q-1)]
                        datadict['IUPred'] = self.obj['QuerySeq'].obj['Disorder'].list['ResidueDisorder'][q-1]
                    except: pass
                    try: datadict[self.obj['QuerySeq'].shortName()] = self.obj['QuerySeq'].info['Sequence'][i]
                    except: pass
                    rje.delimitedFileOutput(self,outfile,headers,datadict=datadict)
            if qryseq: self.printLog('#CONS','Relative conservation (window=%d) for Query "%s" output to %s' % (self.stat['RelConWin'],self.obj['QuerySeq'].shortName(),outfile))
            else: self.printLog('#CONS','Relative conservation (window=%d) output to %s' % (self.stat['RelConWin'],outfile))
            ### ~ [3] ~ Generate PNG files with R (optional) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.opt['MakePNG']:
                ## ~ [3a] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                basefile = string.replace(outfile,'.rel.tdt','')
                relfile = '%s.rel.tdt' % basefile
                alnfile = '%s.aln.tdt' % basefile
                ## ~ [3b] ~ Sort out relfile for R ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if outfile != relfile: open(relfile,'w').write(open(outfile,'r').read())
                ## ~ [3c] ~ Sort out alnfile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.saveR(seqfile=alnfile,name=self.info['SeqName'])
                ## ~ [3d] ~ Call R ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                rcmd = '%s --no-restore --no-save --args "rel2png" "%s"' % (self.info['RPath'],basefile)
                rslimjim = '%srje.r' % self.info['Path']
                if not os.path.exists(rslimjim): rslimjim = '%s/rje.r' % self.info['Path']
                rcmd += ' < "%s" > "%s.r.tmp.txt" 2>&1' % (rslimjim,basefile)
                self.printLog('#RPNG',rcmd)
                problems = os.popen(rcmd).read()
                if problems: self.errorLog(problems,printerror=False)
                pngx = len(glob.glob('%s*png' % basefile))
                self.printLog('#PNG','%d RelCons PNG files made for %s' % (pngx,basefile))
                if pngx and os.path.exists('%s.r.tmp.txt' % basefile): os.unlink('%s.r.tmp.txt' % basefile)
            ### ~ [4] ~ Finish Off ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.obj['QuerySeq']:    # Return original sequence
                for seq in self.seqs():
                    try: seq.info['Sequence'] = seq.info.pop('PreStrip')
                    except: pass
        except: self.errorLog('Problem with rje_seq.relCons()')
#########################################################################################################################
    def relConList(self,window=35,slimcalc=None):   ### Returns list of relative conservation per residue.
        '''
        Returns list of relative conservation per residue, incorporating sequence weighting, using qry and seqs.
        >> window:int = number of residues either side to consider. If window=0, will return absolute values. 
        '''
        try:### ~ [1] Calculate Conservation List ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self._checkAln(aln=True,realign=False):
                self.printLog('#ERR','Cannot generate relConList data from unaligned sequences!')
                raise ValueError
            xse = -1                                # Entropy for use with 100% X columns
            conslist = []
            qry = self.seqs()[0]
            for i in range(qry.seqLen()):   # Each residue in turn
                ## ~ [1a] Shannon Entropy score ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                cfreq = {}          # Dictionary of {aa:column frequency}
                for seq in self.seqs():
                    aa = seq.info['Sequence'][i]
                    if aa not in cfreq: cfreq[aa] = 0.0
                    cfreq[aa] += 1.0
                ## ~ [1b] ~ Gap Penalty and Masked residues ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                gap_pen = 1
                if '-' in cfreq:
                    gap_pen = 1 - (cfreq['-'] / sum(cfreq.values()))    # Proportion of ungapped sequences
                    cfreq.pop('-')
                if not slimcalc or not slimcalc.opt['RelGapPen']: gap_pen = 1
                if 'X' in cfreq: cfreq.pop('X') # X's count towards ungapped positions but contribute nothing to entropy
                if qry.dna() and 'N' in cfreq: cfreq.pop('N')
                ## ~ [1c] Calculate entropy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                #x#self.opt['DeBug'] = True
                base = {True:4,False:20}[qry.dna()]
                cfreq = rje.dictFreq(cfreq,total=False)
                if slimcalc: cfreq = slimcalc.vneEigen(cfreq)
                se = 1.0    
                for aa in cfreq:
                    if cfreq[aa] > 1e-10: se += cfreq[aa] * math.log(cfreq[aa],base)     # Shannon entropy is -ve
                #x#self.deBug('%d: %s = %.3f' % (base,cfreq,se))
                if not cfreq:
                    if xse < 0: xse = self.aaFreqEntropy()
                    se = xse
                if se < 0.0: raise ValueError
                conslist.append(gap_pen * se)   # Final entropy, weighted by non-gap proportion
            ### ~ [2] Convert into relative conservation across window ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            return self.relWinCon(conslist,window,qry)
        except: self.errorLog('Error in rje_seq.relConList()',quitchoice=True)
        return [0.0] * seq.aaLen()              # Return flat nothing if failure.
#########################################################################################################################
    def relWinCon(self,conslist,window=30,qry=None):    ### Returns relative conservation (x-mean)/sd for windows
        '''Returns relative conservation (x-mean)/sd for +/- windows.'''
        try:### ~ [1] ~ Calculate relative conservation (x-mean)/sd for +/- windows ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if window == 0: return conslist
            relcon = []
            for i in range(len(conslist)):      # Each residue in turn
                if qry:
                    win = []
                    for w in range(max(0,i-window),min(len(conslist),i+window+1)):
                        if qry.isDisordered(w) == qry.isDisordered(i): win.append(conslist[w])
                else: win = conslist[max(0,i-window):i+window+1]
                (wmean,wsd) = rje.meansd(win)
                if wsd >= 1e-6: relcon.append((conslist[i] - wmean) / wsd)
                else: relcon.append(0.0)
            return relcon
        except: self.errorLog('Error in rje_seq.relCon()',quitchoice=True)
        return [0.0] * len(conslist)            # Return flat nothing if failure.
#########################################################################################################################
    def aaFreqEntropy(self):   ### Returns entropy based on whole alignment AAFreqs
        '''Returns entropy based on whole alignment AAFreqs.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.seqs()[0].dna(): alph = alph_protx[:-1]
            else: alph = alph_dna[1:]
            freq = {}
            for x in alph: freq[x] = 0.0
            ### ~ [2] Add Freqs from Seqs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs(): freq = seq.aaFreq(freq,newkeys=False)  # Adds to aafreq dictionary
            freq = rje.dictFreq(freq,total=False)                           # Convert to frequencies
            ### ~ [3] Calculate entropy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            base = len(alph)
            se = 1.0    
            for aa in freq: se += freq[aa] * math.log(freq[aa],base)        # Shannon entropy is -ve
            if se < 0.0: raise ValueError
            return se
        except: self.errorLog('Problem with rje_seq.aaFreqEntropy()')
        return 0.0
#########################################################################################################################
### End of SeqList Class                                                                                                #
#########################################################################################################################
    
#########################################################################################################################
##  Sequence Class: Individual sequences
#########################################################################################################################
class Sequence(rje_sequence.Sequence):     
    '''
    Individual Sequence Class. Author: Rich Edwards (2005).
    See rje_sequence.py for details.
    '''
#########################################################################################################################
### End of Sequence Class
#########################################################################################################################

#########################################################################################################################
### DisMatrix Class: 
#########################################################################################################################
class DisMatrix(rje_dismatrix.DisMatrix):     
    '''
    Sequence Distance Matrix Class. Author: Rich Edwards (2005).
    See rje_dismatrix.py for details.
    '''
#########################################################################################################################
### End of DisMatrix Class
#########################################################################################################################
    
#########################################################################################################################
## End of SECTION II
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def pamDis(seqlist,pam):   ### Makes an all-by-all PAM Distance Matrix
    '''
    Makes an all-by-all PAM Distance Matrix.
    - seqin=FILE    : Sequence file (aligned)
    >> seqlist:SeqList Object
    >> pam:rje_pam.PamCtrl Object
    '''
    ### <a> ### Check
    print 'Calculating',
    if seqlist.seqNum() == 0:
        print 'Must load (aligned) sequences using \'seqin=FILE\' command.'
        sys.exit()
    elif seqlist.opt['Aligned'] == False:
        print 'Input sequences must be aligned! Will attempt alignment.'
        seqlist = seqlist.align()
    ### <b> ### Make matrix
    matrix = rje_dismatrix.DisMatrix(log=seqlist.log)
    matrix.info['Type'] = 'PAM'
    matrix.opt['Symmetric'] = True
    perc = rje.setPerc(seqlist.seqNum(),1,seqlist.seqNum()/4,0)
    for s1 in range(seqlist.seqNum()):
        perc = rje.perCounter(perc)
        seq1 = seqlist.seq[s1]
        matrix.addDis(seq1,seq1,0)  # Self
        for s2 in range(s1+1,seqlist.seqNum()):
            seq2 = seqlist.seq[s2]
            p = pam.pamML(ancseq=seq1.info['Sequence'],descseq=seq2.info['Sequence'])
            matrix.addDis(seq1,seq2,p)
    seqlist.obj['PAM Dis'] = matrix
    print '...Saving as %s.pamdis' % seqlist.info['Name']
    matrix.saveMatrix(seqlist.seq,filename='%s.pamdis' % seqlist.info['Name'],delimit=',')
#########################################################################################################################
def MWt(sequence=''):   ### Returns Molecular Weight of Sequence
    '''Returns Molecular Weight of Sequence.'''
    mwtdic = {'A':89, 'V':117, 'L':131, 'I':131, 'P':115,
              'F':165, 'W':204, 'M':149, 'G':75, 'S':105,
              'T':119, 'C':121, 'Y':181, 'N':132, 'Q':146,
              'D':133, 'E':147, 'K':146, 'R':174, 'H':155,
              'X':136.75}
    _mwt = 0.0
    for aa in sequence:
        if aa in mwtdic.keys():
            _mwt += mwtdic[aa]
            _mwt += 18  # H2O
    if _mwt > 18:
        _mwt -= 18
    return _mwt
#########################################################################################################################
def surfaceAccessibility(sequence='',returnlist=True):  ### Returns a list of Janin et al SA values for each residue
    '''
    Returns a list of Janin et al SA values for each residue. Based on method of Norman Davey.
    >> sequence:str = input sequence
    >> returnlist:boolean [True] = Returns list of hydropathies rather than total
    << sa_prob:list of SA values for each residue
    '''
    ### Setup Dictionary ###
    janin_sa = {'A':0.49,'R':0.95,'N':0.81,'D':0.78,'C':0.26,'Q':0.84,'E':0.84,'G':0.48,'H':0.66,'I':0.34,'L':0.40,
                'K':0.97,'M':0.48,'F':0.42,'P':0.75,'S':0.65,'T':0.70,'W':0.51,'Y':0.76,'V':0.36}
    janin_sa['X'] = sum(janin_sa.values())/len(janin_sa)
    ### Make List ###
    sa_prob = []
    for i in range(len(sequence)):
        for j in range(1,7):
            if j == 1:
                try:
                    productL = janin_sa[sequence[i + 3 - j]]    # -3 to +2
                except:
                    productL = janin_sa['X']
                try:
                    productR = janin_sa[sequence[i + 4 - j]]    # -2 to +3
                except:
                    productR = janin_sa['X']
            else:
                try:
                    productL *= janin_sa[sequence[i + 3 - j]]
                except:
                    productL *= janin_sa['X']
                try:
                    productR *= janin_sa[sequence[i + 4 - j]]
                except:
                    productR *= janin_sa['X']
        sa_prob.append((productL+productR)*17.6/2)    # Gets a 1+/- scale
    if returnlist:
        return sa_prob
    return sum(sa_prob)
#########################################################################################################################
def eisenbergHydropathy(sequence='',returnlist=True):  ### Returns the Eisenberg Hydropathy for the sequence
    '''
    Returns the Eisenberg Hydropathy for the sequence.
    >> sequence:str = AA sequence
    >> returnlist:boolean [True] = Returns list of hydropathies rather than total
    << hyd:float = Eisenberg Hydropathy for the sequence
    '''
    ### Setup ###
    eishyd = {'A':0.62, 'R':-2.53, 'N':-0.78, 'D':-0.9, 'C':0.29, 'Q':-0.85, 'E':-0.74, 'G':0.48, 'H':-0.4, 'I':1.38,
              'L':1.06, 'K':-1.5, 'M':0.64, 'F':1.19, 'P':0.12, 'S':-0.18, 'T':-0.05, 'W':0.81, 'Y':0.26, 'V':1.08,
              '-':0, 'X':0}
    eislist = []
    ### Calculate ###
    for aa in sequence: eislist.append(rje.getFromDict(eishyd,aa,False,False,0.0))
    if returnlist: return eislist
    return sum(eislist)
#########################################################################################################################
def rna2dna(seqs=[]):    ### Converts sequences in list from RNA to DNA (U -> T)
    '''
    Converts sequences in list from RNA to DNA (U -> T).
    >> seq:list of Sequence Objects
    '''
    for seq in seqs:
        seq.info['Sequence'] = string.replace(seq.info['Sequence'],'U','T')
#########################################################################################################################
def deGap(seqstring):    ### Degaps sequence
    '''Degaps sequence.'''
    return string.join(string.split(seqstring,'-'),'')
#########################################################################################################################
def pwIDGapExtra(seq1,seq2,nomatch=['X']):    ### Returns a dictionary of stats: 'ID','Gaps','Extra'
    '''
    Returns a dictionary of stats: 'ID','Gaps','Extra', 'Len'. Each stat is a two-element list of [seq1Vseq2,seq2Vseq1].
    Numbers returned are absolute counts and should be divided by 'Len' stat for %.
    >> seq1 & seq2 = strings to be compared.
    >> nomatch = list of characters not to count as ID
    '''
    seq = [seq1,seq2]
    id = [0,0]
    gaps = [0,0]
    extra = [0,0]
    seqlen = [len(seq1) - string.count(seq1,'-'),len(seq2) - string.count(seq2,'-')]
    if seq[0] == seq[1]:
        id[0:] = seqlen[0:]
    else:
        for r in range(len(seq[0])):
            aa = (seq[0][r],seq[1][r])
            if aa[0] == aa[1]:  # Matching residue
                if aa[0] != '-' and aa[0] not in nomatch:    # Not a gap or X
                    id[0] += 1
                    id[1] += 1
            elif aa[0] == '-':  # Gap in seq1 but not seq2
                gaps[0] += 1
                extra[1] += 1
            elif aa[1] == '-':  # Gap in seq2 but not seq1
                gaps[1] += 1
                extra[0] += 1
    return {'ID':id,'Gaps':gaps,'Extra':extra,'Len':seqlen}
#########################################################################################################################
def Blast2Fas(seqlist,outdir='./'): ### Will blast sequences against list of databases and compile a fasta file of results per query
    '''Will blast sequences against list of databases and compile a fasta file of results per query.'''
    try:### ~ [0] ~ General Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        blastdbs = seqlist.list['Blast2Fas'][0:]
        blast2fas = len(seqlist.list['Blast2Fas'])
        for db in blastdbs[0:]:
            if db.lower() in ['','none']: blastdbs.remove(db)
            elif not os.path.exists(db):
                seqlist.errorLog('Blast2Fas Database "%s" not found.' % db,printerror=False)
                blastdbs.remove(db)
        if blast2fas and len(blastdbs) == 0:
            seqlist.errorLog('Blast2Fas option called but no Database files given/found.',printerror=False)
            sys.exit()
        ## ~ [0a] Pre-existing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        fx = 0; rx = 0; bx = 0
        for seq in seqlist.seqs():
            filename = '%s%s.blast.fas' % (outdir,seq.info['AccNum'])
            if os.path.exists(filename):
                fx += 1
                if rje.backup(seqlist,filename): bx += 1
                if not os.path.exists(filename): rx += 1
        seqlist.printLog('#INFO','%s BLAST2Fas outputs found; %s backed up; %s deleted. (Append=%s)' % (rje.iStr(fx),rje.iStr(bx),rje.iStr(rx),seqlist.getBool('Append')))
        ### ~ [1] ~ BLAST Runs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for dbase in blastdbs:
            dbnum = SeqCount(seqlist,dbase)
            bcmd = ['blastv=%d' % dbnum,'blastb=%d' % dbnum] + seqlist.cmd_list + ['blastd=%s' % dbase,'blasti=%s' % seqlist.getStr('Name')]
            blast = rje_blast.blastObj(log=seqlist.log,cmd_list=bcmd,type='Dev')
            blast.setStr({'Name':'%s.%s.blast' % (rje.baseFile(seqlist.info['Name'],True),rje.baseFile(dbase,True))})
            if blast.getStr('Type') in ['blastn','tblastx','tblastn']: blast.formatDB(protein=False,force=False)  
            else: blast.formatDB(force=False)
            seqlist.printLog('#BLAST','BLASTing %s vs %s (%s; e=%s)' % (seqlist.info['Name'],dbase,blast.getStr('Type'),blast.getNum('E-Value')),log=True)
            blast.blast()
            if not blast.readBLAST(clear=True,unlink=not seqlist.getBool('KeepBlast')): raise ValueError
            tmpseq = SeqList(log=seqlist.log,cmd_list=['append=T']+seqlist.cmd_list+['seqin=None'])
            ## ~ [1a] Old BLAST method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if blast.getBool('OldBLAST'):
                for search in blast.search:
                    seq = seqlist.seq[blast.search.index(search)]
                    tmpseq.seq = []
                    if seq.shortName().find(search.info['Name']) == 0:
                        blast.hitToSeq(seqlist=tmpseq,searchlist=[search],filename='%s%s.blast.fas' % (outdir,seq.info['AccNum']),appendfile=True)
                    else:
                        seqlist.errorLog('Search %s does not match sequence %s.' % (search.info['Name'],seq.shortName()))
                    if not tmpseq.seq:
                        if seqlist.info['Name'] == dbase: seqlist.printLog('#ERR','%s does not find itself! Check sequence.' % (seq.info['AccNum']))
                        else: seqlist.printLog('#NULL','No hits for %s vs %s' % (seq.info['AccNum'],dbase))
            ## ~ [1b] New BLAST method ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            else:
                for seq in seqlist.seqs():
                    tmpseq.seq = []
                    filename = '%s%s.blast.fas' % (outdir,seq.info['AccNum'])
                    hitlist = blast.queryHits(seq.shortName())
                    if hitlist: blast.hitToSeq(tmpseq,hitlist,filename=filename,appendfile=True)
                    if tmpseq.seqNum() != len(hitlist): tmpseq.warnLog('Only %d sequences extracted for %d %s hits!' % (tmpseq.seqNum(),len(hitlist),seq.info['AccNum']),'hit_num_mismatch',quitchoice=True,suppress=True)
                    if not tmpseq.seq:
                        if seqlist.info['Name'] == dbase: seqlist.printLog('#ERR','%s does not find itself! Check sequence.' % (seq.info['AccNum']))
                        else: seqlist.printLog('#NULL','No hits for %s vs %s' % (seq.info['AccNum'],dbase))
        ### ~ [2] ~ Generate HAQESAC Batch files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if seqlist.info['HAQBAT'].lower() not in ['none','']:
            rje.backup(seqlist,seqlist.info['HAQBAT'])
            for seq in seqlist.seqs():
                acc = seq.info['AccNum']
                haqcmd = '%shaqesac.py seqin=%s.blast.fas query=%s basefile=%s' % (seqlist.info['Path'],acc,acc,acc)
                open(seqlist.info['HAQBAT'],'a').write('python %s\n' % haqcmd)
            seqlist.printLog('#HAQ','HAQESAC batch commands output to %s' % seqlist.info['HAQBAT'])
    except: seqlist.errorLog('Catastrophic Error in rje_seq.Blast2Fas()')
#########################################################################################################################
def seqInfoList(seqs=[],key='short'): ### Returns list of sequence info with key
    '''Returns list of sequence info with key or 'short' for shortName()'''
    ilist = []
    for seq in seqs:
        if key == 'short':
            ilist.append(seq.shortName())
        elif key in seq.info.keys():
            ilist.append(seq.info[key])
        else:
            ilist.append(None)
    return ilist
#########################################################################################################################
def SeqInfoListFromFile(callobj,filename,key='short',startfrom=None):  ### Zips through file and generates list
    '''
    Zips through file and counts number of sequences.
    >> callobj=Calling Object
    >> filename=Name of Sequence File
    << returns number of sequences
    '''
    try:
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(filename): return 0
        seqlist = SeqList(callobj.log,['i=-1','v=-1','autoload=F','seqin=%s' % filename,'replacechar=F','accnr=F'])
        SEQFILE = open(filename,'r')
        lastline = ''
        sx = 0
        ilist = []
        ### ~ [2] Count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while 1:
            (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
            if startfrom:
                if nextseq.shortName().find('%s ' % startfrom) >= 0: startfrom = None
                else: continue
            seqlist.seq = []
            if nextseq: ilist += seqInfoList([nextseq],key)
            else: break
        SEQFILE.close()
        return ilist
    except:
        callobj.errorLog('rje_seq.SeqInfoListFromFile(%s) buggered.' % filename)
        return 0
#########################################################################################################################
def SeqCount(callobj,filename):     ### Zips through file and counts number of sequences
    '''
    Zips through file and counts number of sequences.
    >> callobj=Calling Object
    >> filename=Name of Sequence File
    << returns number of sequences
    '''
    try:
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(filename): return 0
        seqlist = SeqList(callobj.log,['i=-1','v=-1','autoload=F','seqin=%s' % filename,'replacechar=F','accnr=F'])
        SEQFILE = open(filename,'r')
        lastline = SEQFILE.readline()
        sx = 0
        ### ~ [2] Count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while lastline:
            if lastline[:1] == '>': sx += 1
            lastline = SEQFILE.readline()
        SEQFILE.close()
        return sx
    except:
        callobj.errorLog('rje_seq.SeqCount(%s) buggered.' % filename)
        return 0
#########################################################################################################################
def DBSize(callobj,filename,seqsize=False,nonx=False):   ### Zips through file and counts number of amino acids/nucleotides
    '''
    Zips through file and counts number of amino acids/nucleotides.
    >> callobj=Calling Object
    >> filename=Name of Sequence File
    << returns number of amino acids/nucleotides.
    '''
    try:
        ### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if not os.path.exists(filename): return 0
        seqlist = SeqList(callobj.log,['i=-1','v=-1','autoload=F','seqin=%s' % filename])
        SEQFILE = open(filename,'r')
        lastline = ''
        sx = 0; sizedict = {}
        ### ~ [2] Count ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        while 1:
            (nextseq,lastline) = seqlist.nextFasSeq(SEQFILE,lastline)
            if nextseq:
                if nonx: alen = nextseq.nonX()
                else: alen = nextseq.aaLen()
                sx += alen
                if seqsize:
                    if alen not in sizedict: sizedict[alen] = 0
                    sizedict[alen] += 1
            else: break
            seqlist.seq = []
        SEQFILE.close()
        if seqsize: return sizedict
        return sx
    except:
        callobj.errorLog('rje_seq.SeqCount(%s) buggered.' % filename)
        return 0
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
    try:
        ## Create SeqList ##
        cmd_list = ['i=1'] + cmd_list
        seqlist = SeqList(log=mainlog,cmd_list=cmd_list)

        ## BLAST to Fasta file, special activity ##        
        Blast2Fas(seqlist)

        ## PAM Distance Matrix ###        
        if 'pamdis' in cmd_list:
            mainlog.printLog('#EXE', "Making all-by-all PAM Distance Matrix.", 1)
            pam = rje_pam.PamCtrl(log=mainlog,cmd_list=['pammax=50']+cmd_list)
            pamDis(seqlist,pam)

        ### SeqOut, Reformat and SplitSeq ###
        seqlist.reFormat()

        ### SeqSize ###
        if 'seqsize' in cmd_list:
            sizedb = DBSize(seqlist,seqlist.info['Name'],seqsize=True,nonx='nonx' in cmd_list)
            dfile = seqlist.info['Basefile'] + '.seqsize.tdt'; dhead = ['SeqLen','SeqN']
            rje.delimitedFileOutput(seqlist,dfile,dhead)
            for size in rje.sortKeys(sizedb): rje.delimitedFileOutput(seqlist,dfile,dhead,datadict={'SeqLen':size,'SeqN':sizedb[size]})
                
        # self.saveFasta(seqfile=self.info['SeqOut'])
        # self.splitSeq(self,split=int(self.stat['Split']))

        ### Temp Special Feature ###
        if None and 'reformat=scanseq' in cmd_list: #!# Can this be deleted now?! #!#
            if seqlist.opt['AutoLoad'] and not seqlist.opt['MemSaver'] and seqlist.stat['Split'] < 1:
                seqlist.saveScanSeq()
            elif os.path.exists(seqlist.info['Name']):
                ## Setup ##
                good_acc = rje.listFromCommand(seqlist.info['AccList'])
                bad_spec = rje.listFromCommand(seqlist.info['FilterSpec'])
                splitfile = 0
                splitx = 0
                if seqlist.stat['Split'] > 0:
                    splitfile = 1
                ## Process ##
                seqlist.seq = []
                SEQIN = open(seqlist.info['Name'],'r')
                lastline = 'Firstline'
                if splitfile > 0:
                    seqfile = '%s.%d.scanseq' % (rje.baseFile(seqlist.info['Name'],True),splitfile)
                else:
                    seqfile = '%s.scanseq' % rje.baseFile(seqlist.info['Name'],True)
                out.verbose(0,3,'Converting %s to %s...' % (seqlist.info['Name'],seqfile),0)
                SEQOUT = open(seqfile, 'w')
                (seq,lastline) = seqlist.nextFasSeq(SEQIN,lastline)
                sx = 0
                acc_out = []
                while seq:
                    process = True
                    if good_acc and seq.info['AccNum'] not in good_acc:
                        process = False
                    if process and seq.info['SpecCode'] not in bad_spec and seq.info['AccNum'] not in acc_out:
                        SEQOUT.write('%s %s\n' % (seq.info['AccNum'],re.sub('-','',seq.info['Sequence'])))
                        acc_out.append(seq.info['AccNum'])
                        sx += 1
                        splitx += 1
                        rje.progressPrint(out,sx)
                    seqlist.seq = []
                    (seq,lastline) = seqlist.nextFasSeq(SEQIN,lastline)
                    if splitfile > 0 and splitx >= seqlist.stat['Split']:
                        mainlog.printLog('#OUT',"%d Sequences output to %s in scansite parallel format\n" % (splitx,seqfile),1)
                        splitfile += 1
                        splitx = 0
                        seqfile = '%s.%d.scanseq' % (rje.baseFile(seqlist.info['Name'],True),splitfile)
                        out.verbose(0,3,'... %s...' % (seqfile),0)
                        SEQOUT.close()
                        SEQOUT = open(seqfile, 'w')
                SEQOUT.close()
                SEQIN.close()
                out.verbose(0,1,'Done!',1)
                mainlog.printLog('#OUT',"%d Sequences output to %s in scansite parallel format\n" % (splitx,seqfile),1)
        
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
### END OF SECTION IV
#########################################################################################################################
