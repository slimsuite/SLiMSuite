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
Version:      1.45.2
Last Edit:    14/12/20
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

Edit:
    Version 1.16 introduced an interactive edit mode (`edit=T`) that gives users the options to rearrange, copy, delete,
    split, truncate, rename, join, merge (as consensus) etc. Please contact the author for more details.

    From Version 1.20, a delimited text file can also be given
    as `edit=FILE`, which should contain: Locus, Pos, Edit, Details. Edit is the type of change (INS/DEL/SUB) and Details
    contains the nature of the change (ins/sub sequence or del length). Edits are made in reverse order per locus to
    avoid position conflicts and overlapping edits should be avoided. WARNING: These will not be checked for! An optional
    Notes field will be used if present for annotating changes in the log file.

    From Version 1.22, a delimited file can be given in place of Start,End for region=X. This file should contain Locus,
    Start, End and NewAcc fields. If No NewAcc field is present, the new accession number will be the previous accnum
    (extracted from the sequence name) with '.X-Y' appended.

Commandline:
    ### ~ INPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FILE      : Sequence input file name. [None]
    seqmode=X       : Sequence mode, determining method of sequence storage (full/list/file/index/db/filedb). [file]
    seqdb=FILE      : Sequence file from which to extract sequences (fastacmd/index formats) [None]
    seqindex=T/F    : Whether to save (and load) sequence index file in file mode. [True]
    seqformat=X     : Expected format of sequence file [None]
    seqtype=X       : Sequence type (prot(ein)/dna/rna/mix(ed)) [None]
    mixed=T/F       : Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    dna=T/F         : Alternative option to indicate dealing with nucleotide sequences [False]
    autoload=T/F    : Whether to automatically load sequences upon initialisation. [True]
    autofilter=T/F  : Whether to automatically apply sequence filtering. [True]
    duperr=T/F      : Whether identification of duplicate sequence names should raise an error [True]

    ### ~ SEQUENCE FORMATTING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    reformat=X      : Output format for sequence files (fasta/short/acc/acclist/accdesc/speclist/index/dna2prot/peptides/(q)region/revcomp/reverse/descaffold) [fasta]
    rename=T/F      : Whether to rename sequences [False]
    spcode=X        : Species code for non-gnspacc format sequences [None]
    newacc=X        : New base for sequence accession numbers - will rename sequences [None]
    newgene=X       : New gene for renamed sequences (if blank will use newacc or 'seq' if none read) [None]
    genecounter=T/F : Whether new gene have a numbered suffix (will match newacc numbering) [False]
    newdesc=FILE    : File of new names for sequences (over-rules other naming). First word should match input [None]
    keepname=T/F    : Whether to keep the original name (first word) when mapping with newdesc=FILE [True]
    concatenate=T   : Concatenate sequences into single output sequence named after file [False]
    split=X         : String to be inserted between each concatenated sequence [''].
    seqshuffle=T/F  : Randomly shuffle each sequence without replacement (maintains monomer composition) [False]
    region=X,Y      : Alignment/Query region to use for reformat=peptides/(q)region reformatting of fasta alignment (1-L) [1,-1]
    edit=T/F/FILE   : Enter sequence edit mode upon loading (will switch seqmode=list) (see above) [False]
    gnspacc=T/F     : Whether to automatically try to enforce SLiMSuite gene_SPCODE__AccNum format [True]

    ### ~ DNA TRANSLATIONS (reformat=dna2prot) ~~~~~~~~~~~~ ###
    minorf=X        # Min. ORF length for translated sequences output. -1 for single translation inc stop codons [-1]
    terminorf=X     # Min. length for terminal ORFs, only if no minorf=X ORFs found (good for short sequences) [-1]
    orfmet=T/F      # Whether ORFs must start with a methionine (before minorf cutoff) [True]
    rftran=X        # No. reading frames (RF) into which to translate (1,3,6) [1]

    ### ~ FILTERING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqnr=T/F       : Whether to check for redundancy on loading. (Will remove, save and reload if found) [False]
    grepnr=T/F      : Whether to use grep based forking NR mode (needs sized-sorted one-line-per-sequence fasta) [True]
    twopass=T/F     : Whether to perform second pass looking for redundancy of earlier sequences within later ones [True]
    revcompnr=T/F   : Whether to check reverse complement for redundancy too [True]
    goodX=LIST      : Inclusive filtering, only retaining sequences matching list []
    badX=LIST       : Exclusive filtering, removing sequences matching list []
    - where X is 'Acc', Accession number; 'Seq', Sequence name; 'Spec', Species code; 'Desc', part of name;
    minlen=X        : Minimum sequence length [0]
    maxlen=X	    : Maximum length of sequences (<=0 = No maximum) [0]

    ### ~ EXTRACT/MASK OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    maskseq=TDTFILE : File of Locus, Start, End positions for masking [None]
    grabseq=TDTFILE : File of Locus, Start, End positions for region extraction [None]
    posfields=LIST  : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    addflanks=INT   : Length of flanking sequence to also extract/mask [0]

    ### ~ SEQUENCE TILING OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    tile=INT        : Tile sequences into INT bp chunks (<1 for no tiling) [0]
    mintile=X       : Min. length for tile, else appended to previous (<1 for proportion of tile=INT) [0.1]
    tilestep=0      : Gap between end of one tile and start of next. Can be negative [0]
    tilename=STR    : Tile naming strategy (pos, start, num, purepos, purestart, purenum) [pos]

    ### ~ OUTPUT OPTIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqout=FILE     : Whether to output sequences to new file after loading and filtering [None]
    usecase=T/F     : Whether to return sequences in the same case as input (True), or convert to Upper (False) [False]
    sortseq=X       : Whether to sort sequences prior to output (size/invsize/accnum/name/seq/species/desc) [None]
    sampler=N(,X)   : Generate (X) file(s) sampling a random N sequences from input into seqout.N.X.fas [0]
    summarise=T/F   : Generate some summary statistics in log file for sequence data after loading [False]
    genomesize=X    : Genome size for NG50 and LG50 summary output (if >0) [0]
    raw=T/F         : Adjust summary statistics for raw sequencing data [False]
    fracstats=T/F   : Output a table of N(G)XX and L(G)XX statistics for a range of XX [False]
    fracstep=INT    : Step size for NXX and LXX fractions (1/2/5/10/25) [5]
    lenstats=LIST   : List of min sequence lengths to output stats for (raw=T) []
    gapstats=T/F    : Output a summary of assembly gap sizes and positions [False]
    mingap=INT      : Minimum length of a stretch of N bases to count as a gap (0=None unless gapstats=T) [10]
    gapfix=X:Y(,X:Y): List of gap lengths X to convert to different lengths Y []
    maker=T/F       : Whether to extract MAKER2 statistics (AED, eAED, QI) from sequence names [False]
    splitseq=X      : Split output sequence file according to X (gene/species) [None]
    tmpdir=PATH     : Directory used for temporary files ['./tmp/']

See also rje.py generic commandline options.

Uses general modules: copy, glob, os, string, sys, time
Uses RJE modules: rje, rje_obj, rje_zen
Other modules needed: None
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, random, re, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_db, rje_menu, rje_obj, rje_sequence, rje_uniprot, rje_zen
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
    # 1.15.4 - Fixed REST server output bug.
    # 1.15.5 - Fixed reformat=fasta default issue introduced from fixing REST output bug.
    # 1.16.0 - Added edit=T sequence edit mode upon loading (will switch seqmode=list).
    # 1.17.0 - Added additional summarise=T output for seqmode=db.
    # 1.18.0 - Added revcomp to reformat options.
    # 1.19.0 - Added option log description for deleting sequence during edit.
    # 1.20.0 - Added option to give a file of changes for edit mode.
    # 1.20.1 - Fixed edit=FILE deletion bug.
    # 1.21.0 - Added capacity to add/update database object from self.summarise() even if not seqmode=db. Added filedb mode.
    # 1.22.0 - Added geneDic() method.
    # 1.23.0 - Added seqSequence() method.
    # 1.24.0 - Add NNN gaps option and "delete rest of sequences" to edit().
    # 1.24.1 - Minor edit bug fix and DNA toggle option.
    # 1.25.0 - Added loading of FASTQ files in seqmode=file mode.
    # 1.26.0 - Updated sequence statistics and fixed N50 underestimation bug.
    # 1.26.1 - Fixed median length overestimation bug.
    # 1.26.2 - Fixed sizesort bug. (Now big to small as advertised.)
    # 1.27.0 - Added grepNR() method (dev only). Switched default to RevCompNR=T.
    # 1.28.0 - Fixed second pass NR naming bug and added option to switch off altogether.
    # 1.29.0 - Added maker=T/F : Whether to extract MAKER2 statistics (AED, eAED, QI) from sequence names [False]
    # 1.30.0 - Updated and improved DNA2Protein.
    # 1.31.0 - Added genecounter to rename option for use with other programs, e.g. PAGSAT.
    # 1.31.1 - Fixed edit bug when not in DNA mode.
    # 1.32.0 - Added genomesize and NG50/LG50 to DNA summarise.
    # 1.32.1 - Fixed LG50/L50 bug.
    # 1.32.2 - Added reformat=accdesc to generate output without gene and species code.
    # 1.32.3 - Added checkNames() to check for duplicate sequence names and/or lack of gnspacc format.
    # 1.32.3 - Added duperr=T/F : Whether identification of duplicate sequence names should raise an error [True]
    # 1.33.0 - Added newdesc=FILE : File of new names for sequences (over-rules other naming). First word should match input [None]
    # 1.33.1 - Fixed bug with appending sequences with gap insertion.
    # 1.34.0 - Added genecounter=T/F : Whether new gene have a numbered suffix (will match newacc numbering) [False]
    # 1.35.0 - Added initial extraction of sequences from BLASTDB from rje_seq.
    # 1.36.0 - Added bpFromStr(seqlen)
    # 1.36.1 - Changed default duplicate suffix to X2.
    # 1.37.0 - Added masking and extraction from loaded table of positions.
    # 1.38.0 - Added assembly gap summary and manipulation fundtions.
    # 1.39.0 - Added descaffolding, tiling output and gnspacc=T/F to control edit renaming.
    # 1.40.0 - Added keepname=T/F : Whether to keep the original name (first word) when mapping with newdesc=FILE [True]
    # 1.41.0 - Added contig N50 and L50 output. Tweaked tiling output to leave off name suffix when full length sequence.
    # 1.41.1 - Fixed contig N50 and L50 output. (Previously not sorted!)
    # 1.42.0 - Added tabular summary output for different L/N(G) values.
    # 1.42.1 - Switched mingap=INT to 0=None unless gapstats=T.
    # 1.43.0 - Added raw=T/F and lenstats=LIST to adjust summary statistics for raw sequencing data
    # 1.43.1 - Added sequence reversal (not complemented) to reformat and edit
    # 1.44.0 - Added some additional parsing of common sequence formats from rje_sequence: need to expand.
    # 1.45.0 - Modified the newDesc() method for updating descriptions.
    # 1.45.1 - Added CtgNum to output stats.
    # 1.45.2 - Slight increase of gap extraction speed.
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
    # [Y] : Add warning if region set but no reformat? (Check whether other programs use it.) Add seqin to REST?
    # [ ] : Add seqbatch=FILES capacity for reading in sequences from multiple files.
    # [Y] : Add handling of FASTQ files in addition to FASTA.
    # [ ] : Add stripgap=X method, based on rje_seq stripGap() method.
    # [ ] : Add align=T [alnprog=X] methods directly using rje_seq.
    # [Y] : Add raw=T/F setting to slightly adjust the summarise stats. (X coverage and no gaps.)
    '''
#########################################################################################################################
def makeInfo(): ### Makes Info object which stores program details, mainly for initial print to screen.
    '''Makes Info object which stores program details, mainly for initial print to screen.'''
    (program, version, last_edit, copy_right) = ('SeqList', '1.45.2', 'December 2020', '2011')
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
            print('\n\nHelp for %s %s: %s\n' % (info.program, info.version, time.asctime(time.localtime(info.start_time))))
            out.verbose(-1,4,text=__doc__)
            if rje.yesNo('Show general commandline options?'): out.verbose(-1,4,text=rje.__doc__)
            if rje.yesNo('Quit?'): sys.exit()           # Option to quit after help
            cmd_list += rje.inputCmds(out,cmd_list)     # Add extra commands interactively.
        elif out.stat['Interactive'] > 1: cmd_list += rje.inputCmds(out,cmd_list)    # Ask for more commands
        ### ~ [3] ~ Return commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        return cmd_list
    except SystemExit: sys.exit()
    except KeyboardInterrupt: sys.exit()
    except: print('Major Problem with cmdHelp()')
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
    except: print('Problem during initial setup.'); raise
#########################################################################################################################
file_ext = {'fasta':'fas','short':'fas','acc':'fas','acclist':'acc','speclist':'txt',
            'index':'fas',  # This does not seem to work? (At least for REST output.)
            'dna2prot':'fas','peptides':'txt','qregion':'fas','region':'fas'}
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
    - Edit = Option edit file for sequence editing [None]
    - GrabSeq=TDTFILE : File of Locus, Start, End positions for region extraction [None]
    - MaskSeq=TDTFILE    : File of Locus, Start, End positions for masking [None]
    - Name = Sequence file name - specifies output. [None]
    - NewAcc = New base for sequence accession numbers - will rename sequences [None]
    - NewDesc=FILE    : File of new names for sequences (over-rules other naming). First word should match input [None]
    - NewGene = New gene for renamed sequences (if blank will use newacc) [None]
    - Region = Query region to use for peptides/qregion reformatting of fasta alignment [1,-1]
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
    - TileName = Tile naming strategy (pos, start, num, purepos, purestart, purenum) [pos]
    - TmpDir = Directory used for temporary files ['./tmp/']

    Bool:boolean
    - AutoFilter = Whether to automatically apply sequence filtering. [True]
    - AutoLoad = Whether to automatically load sequences upon initialisation. [True]
    - Concatenate = Concatenate sequences into single output sequence named after file [False]
    - DNA = Alternative option to indicate dealing with nucleotide sequences [False]
    - DupErr = Whether identification of duplicate sequence names should raise an error [True]
    - Edit = Enter sequence edit mode upon loading (will switch seqmode=list) [False]
    - FracStats=T/F      : Output a table of N(G)XX and L(G)XX statistics for a range of XX [False]
    - GapStats=T/F    : Output a summary of assembly gap sizes and positions [False]
    - GeneCounter=T/F : Whether new gene have a numbered suffix (will match newacc numbering) [False]
    - GeneSpAcc=T/F     : Whether to automatically try to enforce SLiMSuite gene_SPCODE__AccNum format [True]
    - GrepNR = Whether to use grep based forking NR mode (needs sized-sorted one-line-per-sequence fasta) [True]
    - KeepName=T/F    : Whether to keep the original name (first word) when mapping with newdesc=FILE [True]
    - Maker = Whether to extract MAKER2 statistics (AED, eAED, QI) from sequence names [False]
    - Mixed = Whether to allow auto-identification of mixed sequences types (else uses first seq only) [False]
    - ORFMet = Whether ORFs must start with a methionine (before minorf cutoff) [True]
    - Raw=T/F         : Adjust summary statistics for raw sequencing data [False]
    - ReName = Whether to rename sequences (will need newacc and spcode) [False]
    - RevCompNR = Whether to check reverse complement for redundancy too [True]
    - SeqIndex = Whether to save (and load) sequence index file in file mode. [True]
    - SeqNR = Whether to check for redundancy on loading. (Will remove, save and reload if found) [False]
    - SeqShuffle = Randomly shuffle each sequence (cannot use file or index mode) [False]
    - SizeSort = Sort sequences by size big -> small re-output prior to loading/filtering [False]
    - Summarise = Generate some summary statistics in log file for sequence data after loading [False]
    - TwoPass=T/F     : Whether to perform second pass looking for redundancy of earlier sequences within later ones [True]
    - UseCase = Whether to return sequences in the same case as input (True), or convert to Upper (False) [False]

    Int:integer
    - AddFlanks=INT   : Length of flanking sequence to also extract/mask [0]
    - FracStep=INT         : Step size for NXX and LXX fractions (1/2/5/10/25) [5]
    - LenStats=LIST : List of min sequence lengths to output stats for (raw=T) []
    - MinGap=INT      : Minimum length of a stretch of N bases to count as a gap [10]
    - MinLen = Minimum sequence length [0]
    - MaxLen = Maximum length of sequences (<=0 = No maximum) [0]
    - MinORF = Min. ORF length for translated sequences output. -1 for single translation inc stop codons [-1]
    - RFTran = No. reading frames (RF) into which to translate (1,3,6) [1]
    - TerMinORF = Min. length for terminal ORFs, only if no minorf=X ORFs found (good for short sequences) [-1]
    - Tile=INT        : Tile sequences into INT bp chunks (<1 for no tiling) [0]
    - TileStep=0      : Gap between end of one tile and start of next. Can be negative [0]

    Num:float
    - GenomeSize=X    : Genome size for NG50 and LG50 summary output (if >0) [0]
    - MinTile=X       : Min. length for tile, else appended to previous (<1 for proportion of tile=INT) [0.1]

    List:list
    - Edit = Stores edits made during self.edit() as (Type,AccNum,Details) tuples. (Enables extraction of edits made.)
    - PosFields=LIST  : Fields in checkpos file to give Locus, Start and End for checking [Locus,Start,End]
    - Sampler = N(,X) = Generate (X) file(s) sampling a random N sequences from input into seqout.N.X.fas [0]
    - Seq = List of sequences - nature depends on SeqMode []

    Dict:dictionary
    - GapFix=X:Y(,X:Y): List of gap lengths X to convert to different lengths Y []
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
        self.strlist = ['Edit','Name','NameFormat','NewAcc','NewDesc','Region','GrabSeq','MaskSeq',
                        'SeqDB','SeqDictType','SeqFormat','SeqIn','SeqMode','SeqType','SeqOut',
                        'Reformat','SpCode','SeqNR','NewGene','Split','SortSeq','SplitSeq','TileName','TmpDir']
        self.boollist = ['AutoFilter','AutoLoad','Concatenate','DNA','DupErr','Edit','GapStats','GeneCounter','GrepNR',
                         'GeneSpAcc','Maker','Mixed','FracStats','ORFMet','ReName','RevCompNR','SizeSort','TwoPass','KeepName',
                         'Raw','SeqIndex','SeqShuffle','Summarise','UseCase']
        self.intlist = ['AddFlanks','FracStep','MinGap','MinLen','MaxLen','MinORF','RFTran','TerMinORF','Tile','TileStep']
        self.numlist = ['GenomeSize','MinTile']
        self.listlist = ['Edit','LenStats','PosFields','Sampler','Seq']
        self.dictlist = ['Filter','GapFix','SeqDict']
        self.objlist = ['Current','CurrSeq','DB','SEQFILE','INDEX']
        ### ~ Defaults ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self._setDefaults(str='None',bool=False,int=0,num=0.0,obj=None,setlist=True,setdict=True)
        self.setStr({'SeqMode':'file','ReFormat':'None','Region':'1,-1','TileName':'pos','TmpDir':rje.makePath('./tmp/')})
        self.setBool({'AutoFilter':True,'AutoLoad':True,'DupErr':True,'GapStats':False,'GeneSpAcc':True,'KeepName':True,
                      'ORFMet':True,'Raw':False,'SeqIndex':True,'GrepNR':True,'RevCompNR':True,'TwoPass':True})
        self.setInt({'FracStep':5,'MinGap':10,'MinORF':-1,'RFTran':1,'TerMinORF':-1,'Tile':0,'TileStep':0})
        self.setNum({'MinTile':0.1})
        self.list['PosFields'] = ['Locus','Start','End']
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
                self._cmdReadList(cmd,'str',['NewAcc','NewGene','Region','SeqFormat','SeqMode','ReFormat','SpCode','SeqType','Split','SortSeq','SplitSeq','TileName'])
                self._cmdReadList(cmd,'file',['Edit','NewDesc','SeqDB','SeqIn','SeqOut','GrabSeq','MaskSeq'])
                self._cmdReadList(cmd,'path',['TmpDir'])
                self._cmdReadList(cmd,'int',['AddFlanks','FracStep','MinGap','MinLen','MaxLen','MinORF','RFTran','TerMinORF','Tile','TileStep'])
                self._cmdReadList(cmd,'num',['GenomeSize','MinTile'])
                self._cmdReadList(cmd,'list',['PosFields'])
                self._cmdReadList(cmd,'ilist',['LenStats'])
                self._cmdReadList(cmd,'nlist',['Sampler'])
                self._cmdReadList(cmd,'cdict',['GapFix'])
                self._cmdRead(cmd,type='bool',att='GeneSpAcc',arg='gnspacc')
                self._cmdReadList(cmd,'bool',['Align','AutoFilter','AutoLoad','Concatenate','DNA','DupErr','Edit','GapStats','GeneCounter','GeneSpAcc','GrepNR','KeepName','Maker','Mixed','FracStats','ORFMet','Raw','ReName','RevCompNR','SizeSort','SeqIndex','SeqNR','SeqShuffle','Summarise','TwoPass','UseCase'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ## ~ [1a] ~ Tidy Commands ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.getStrLC('SeqMode') == 'tuple': self.setStr({'SeqMode':'list'})
        if self.getStrLC('Edit') and not rje.exists(self.getStr('Edit')) and self.getStrLC('Edit') not in ['t','f','true','false']:
            raise IOError('Cannot find edit file (edit=%s)' % rje.exists(self.getStr('Edit')))
        elif self.getStrLC('Edit') and rje.exists(self.getStr('Edit')): self.setBool({'Edit':True})
        if self.getBool('Edit') and self.i() < 0:
            self.warnLog('Edit mode requested but i<0. Need interactivity for edit=T. Edit:False.')
            self.setBool({'Edit':False})
        elif self.getBool('Edit') and not self.getBool('AutoLoad'):
            self.warnLog('Edit mode requested but autoload=F. Edit:False.')
            self.setBool({'Edit':False})
        elif self.getBool('Edit') and self.getStr('SeqMode') != 'list':
            self.setStr({'SeqMode':'list'})
            self.printLog('#EDIT','Edit mode detected: seqmode=list')
        if not self.list['Sampler']: self.list['Sampler'] = [0,1]
        elif len(self.list['Sampler']) < 2: self.list['Sampler'].append(1)
        elif len(self.list['Sampler']) > 2:
            self.warnLog('Too many values given to sampler=N,X. Will use sampler=%s,%s' % self.list['Sampler'][:2])
            self.list['Sampler'] = self.list['Sampler'][:2]
        if self.getBool('Edit') and int(self.list['Sampler'][0]) > 0:
            self.list['Sampler'] = [0,1]
            self.printLog('#EDIT','Edit mode not compatible with Sampler: sampler=0[,1]')
        if not self.getStrLC('SortSeq') and self.getBool('SizeSort'): self.setStr({'SortSeq':'invsize'})
        else: self.setStr({'SortSeq':self.getStrLC('SortSeq')})
        if self.getInt('RFTran') not in [1,3,6]:
            self.warnLog('rftran=%d not recognised: will use rftran=1' % self.getInt('RFTran'))
            self.setInt({'RFTran':1})
        ## ~ [1b] ~ REST Command setup/adjustment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        if self.getStrLC('Rest') in string.split('fasta/short/acc/acclist/accdesc/speclist/index/dna2prot/rna2prot/translate/nt2prot/peptides/qregion/region/descaffold','/'):
            self.setStr({'ReFormat':self.getStrLC('Rest')})
            self.dict['Output'][self.getStrLC('Rest')] = 'SeqOut'
        elif self.getStrLC('Rest') and self.getStrLC('ReFormat'):
            if self.getStrLC('ReFormat') == 'rest': self.setStr({'ReFormat':'fasta'})
            self.dict['Output']['seqout'] = 'SeqOut'
        elif self.getStrLC('ReFormat') == 'rest': self.setStr({'ReFormat':'None'})
        if self.getStrLC('Region') and self.getStr('Region') != '1,-1' and self.getStrLC('ReFormat') not in ['qregion','region']:
            self.warnLog('Region=%s but no reformat=(q)region: no region trimming.' % self.getStr('Region'))
        if self.getStrLC('ReFormat') == 'dna2prot' and self.getStrLC('Rest'): self.setBool({'DNA':True})
        ### ~ [2] ~ AutoLoad option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.getBool('DNA') and not self.nt(): self.str['SeqType'] = 'dna'
        self._filterCmd()
        if self.getBool('AutoLoad'):
            if self.getBool('SeqShuffle'): self.shuffleSeq()
            #if self.getBool('SizeSort'): self.sizeSort()
            if self.getStr('SortSeq'): self.seqSort()
            if self.getBool('SeqNR'):
                if self.getBool('GrepNR'): self.grepNR()
                else: self.seqNR(twopass=not self.getStrLC('SortSeq'),grepnr=self.getBool('GrepNR'))
            if self.getBool('Concatenate'): self.concatenate()
            if self.getBool('ReName'): self.rename(genecounter=self.getBool('GeneCounter'))
            self.loadSeq()
            if self.getBool('AutoFilter'): self.filterSeqs(screen=self.v()>0)
            if self.getStrLC('MaskSeq'): self.grabSeq(mask=True)
            if self.getStrLC('GrabSeq'): self.grabSeq(mask=False)
            if self.getBool('Edit'):
                if not self.edit(): return False
            if self.dict['GapFix']: self.gapFix()
            if self.getBool('Summarise'): self.summarise()
            if self.list['Sampler'][0] > 0: self.sampler()
            elif self.getStrLC('SplitSeq'): self.splitSeq()
            elif self.getInt('Tile') > 0:
                self.tileSeq()
                #,'MinTile'
            elif self.getStrLC('SeqOut'): self.saveSeq()        # Will handle any reformatting
            elif self.getStrLC('ReFormat'):
                if self.getStrLC('ReFormat') == 'fasta': linkext = 'out'
                else: linkext = self.getStrLC('ReFormat')
                if self.getStrLC('Basefile'):  seqout = '%s.%s.%s' % (self.baseFile(),linkext,file_ext[self.getStrLC('ReFormat')])
                else: seqout = '%s.%s.%s' % (rje.baseFile(self.getStr('SeqIn')),linkext,file_ext[self.getStrLC('ReFormat')])
                if self.i() < 0 or rje.yesNo('Reformat (%s) and save to %s?' % (self.getStrLC('ReFormat'),seqout)):
                    self.setStr({'SeqOut':seqout})
                    self.saveSeq(seqfile=seqout)
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        The `seqlist` server is primarily for simple reformatting and sequence manipulation tasks:
        - fasta = standard `gene_SPECIES__AccNum Description` fasta format
        - short = fasta format without any description
        - acc = fasta format with accession numbers (only) as sequence names
        - acclist = plain text list of accession numbers of sequences
        - speclist = plain text list of Uniprot species codes for sequences
        - dna2protrna2prot/translate/nt2prot = translation of DNA (or RNA) sequence into protein
        - peptides = plain list (without names) of protein sequences
        - qregion = fasta alignment restricted to the columns incorporating the given sequence region of the query (sequence 1)
        - region = fasta alignment restricted to the given columns of the alignment

        If &rest=X, where `X` is in the above list, the relevant reformatting will be triggered and the resulting text
        output returned. Otherwise, output is &rest=seqout.

        Additional sequence filtering (degapping etc.) can be performed with the related `seq` server, which has
        &rest=seqout output only. (See: http://rest.slimsuite.unsw.edu.au/seq for commandline options.)

        Run with &rest=help for general options. Run with &rest=full to get full server output as text or &rest=format
        for more user-friendly formatted output. Individual outputs can be identified/parsed using &rest=OUTFMT.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Output']['seqin'] = 'SeqIn'
            #for outfmt in self.restOutputOrder(): self.dict['Output'][outfmt] = 'No output generated.'
            #!# Add specific program output here. Point self.dict['Output'][&rest=X] to self.str key.
            # NOTE: Currently, this method is not called (it is handled by _cmdList above) but is present for docs.
            return
        except: self.errorLog('RestSetup error')
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
    def units(self): return {True:'nt',False:'aa'}[self.nt()]
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
    def geneDic(self):  ### Returns a dictionary of {gene:[seqs]} (rather than one-to-one mapping of seqNameDic())
        '''Returns a dictionary of {gene:[seqs]} (rather than one-to-one mapping of seqNameDic()).'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            genedic = {}
            ### ~ [1] Make dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                gene = self.seqGene(seq)
                if gene not in genedic: genedic[gene] = []
                genedic[gene].append(seq)
            return genedic
        except:
            self.errorLog('geneDic problem')
            raise ValueError
#########################################################################################################################
    def makeSeqNameDic(self,keytype=None,clear=True,warnings=True):  ### Make SeqDict sequence name dictionary.
        '''
        >> keytype:str [None] = Type of data to use as key for dictionary (accnum/short/name/max/loci)
        >> clear:bool [True] = whether to clear self.dict['SeqDict'] before filling
        >> warnings:bool [True] = whether to warn when keys are getting over-written
        << returns self.dict['SeqDict']
        '''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if keytype: self.setStr({'SeqDictType':keytype.lower()})
            if not self.getStrLC('SeqDictType'): self.setStr({'SeqDictType':'short'})
            keytype = self.getStrLC('SeqDictType')
            if keytype not in ['name','full','max','short','acc','accnum','id','loci']: raise ValueError('SeqNameDic keytype "%s" not recognised!' % keytype)
            if clear: self.dict['SeqDict'] = {}
            ### ~ [1] Populate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                skeys = []
                name = self.getSeq(seq,'tuple')[0]
                if not self.gnSpAcc(name) and self.getBool('GnSpAcc'):
                    data = self.extractNameDetails(name)
                    name = '{0}_{1}__{2} {3}'.format(data['Gene'],data['SpecCode'],data['AccNum'],data['Description'])
                if keytype in ['name','full','max']: skeys.append(name)
                if keytype in ['short','full','max','loci']: skeys.append(string.split(name)[0])
                if keytype in ['acc','accnum','full','max','loci']: skeys.append(string.split(string.split(name)[0],'__')[-1])
                if keytype in ['id','full','max']: skeys.append(string.split(string.split(name)[0],'__')[0])
                skeys = rje.sortUnique(skeys)
                #self.bugPrint('%s >> %s' % (name,skeys))
                for skey in skeys:
                    if warnings and skey in self.dict['SeqDict']:
                        self.warnLog('Sequence "%s" already in SeqDict.' % skey,suppress=True)
                        self.debug('???')
                    self.dict['SeqDict'][skey] = seq
            #self.debug('%s: %d seq -> %d seqdict' % (keytype,self.seqNum(),len(self.dict['SeqDict'])))
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
    def seqSequence(self,seq=None):
        if seq == None: seq = self.obj['Current']
        return self.getSeq(seq,'tuple')[1]
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
    def names(self,short=True,seqs=None): ### Returns list of sequence (short) names
        '''Returns list of sequence (short) names.'''
        if not seqs: seqs = self.seqs()
        names = []
        for seq in seqs:
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
            if type(seq) == type(()): mode = 'list'
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
                elif format in ['pos','index']: return self.list['Seq'].index(seq)
                else: return seq
            elif mode == 'db':       # Store sequence data in database object.
                (name,sequence) = (seq['Name'],seq['Sequence'])
                if format == 'obj': return self.getSeqObj(name,sequence)
                elif format == 'tuple':
                    if case: return (name,sequence)
                    else: return (name,sequence.upper())
                elif format == 'short': return string.split(name)[0]
                else: return seq
            elif mode == 'filedb':       # Store sequence data in database object.
                return self.getSeq(seq['FPos'],format,'file',case)
            elif mode == 'file':     # List of file positions.
                SEQFILE = self.SEQFILE()
                SEQFILE.seek(seq)
                name = rje.chomp(SEQFILE.readline())
                if name[:1] == '@': # FASTQ
                    name = name[1:]
                    sequence = ''; line = rje.chomp(SEQFILE.readline())
                    while line and line[:1] != '+':
                        sequence += line
                        line = rje.chomp(SEQFILE.readline())
                elif name[:1] == '>':
                    name = name[1:]
                    sequence = ''; line = rje.chomp(SEQFILE.readline())
                    while line and line[:1] != '>':
                        sequence += line
                        line = rje.chomp(SEQFILE.readline())
                else:
                    self.deBug(name)
                    self.errorLog('Given file position that is not name line. May have SeqList objects mixed up?',printerror=False)
                    raise ValueError    # Must be Fasta or Fastq!
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
            #!# Find out if any programs use this code, then change to .X-Y not -X.Y
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
    def gnSpAcc(self,name): ### Returns whether the sequence name is in GeneSpAcc format
        '''Returns whether the sequence name is in GeneSpAcc format.'''
        if rje.matchExp('^(\S+)_(\S+)__(\S+)',name): return True
        else: return False
#########################################################################################################################
    def extractNameDetails(self,name): return rje_sequence.extractNameDetails(name,self)
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
    def summarise(self,seqs=None,basename=None,sumdb=False,save=True):    ### Generates summary statistics for sequences
        '''
        Generates summary statistics for sequences.
        >> seqs:list None = list of seqs if not using all
        >> basename:str = prefix for output files
        >> sumdb:bool = Whether to store individual sequence stats in 'sequences' table
        >> save:bool = Save stat tables to files
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            raw = self.getBool('Raw') and self.dna()
            gapstats = self.dna() and self.getBool('GapStats') and not raw     # Might need to move this to function argument?
            fracstats = self.dna() and self.getBool('FracStats')
            fracstep = self.getInt('FracStep')
            if fracstats and fracstep not in [1,2,5,10,25] and not raw:
                self.warnLog('fracstep=%d not a recognised value: defaulting to 5' % fracstep)
                fracstep = 5
            if fracstats and raw and self.getNum('GenomeSize') <= 0:
                self.warnLog('fracstats=T raw=T needs genomesize=INT set')
                fracstats = False
            lenstats = self.list['LenStats']    # List of min sequence lengths to output stats for
            #!# Add lendb if raw=F
            if not raw: lenstats = []
            mingap = max(0,self.getInt('MinGap'))
            if gapstats: mingap = max(1,mingap)
            if raw: mingap = 0
            gapre = re.compile('N{%d,}' % mingap)
            seqbase = rje.baseFile(self.getStr('SeqIn'),strip_path=True)
            if self.getStrLC('SeqOut'):
                seqbase = rje.baseFile(self.getStr('SeqOut'),strip_path=True)
            if not self.getStrLC('SeqIn'): seqbase = self.getStr('Basefile')
            if not basename: basename = seqbase
            if sumdb or gapstats or fracstats or lenstats or self.mode().endswith('db'): # {'Name':name,'Sequence':sequence,'FPos':fpos}
                if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
                db = self.db()
            seqdata = {'File':self.getStr('SeqIn'),'GC':0,'GapLength':0,'NCount':0,'GapCount':0}    # SeqNum, TotLength, MinLength, MaxLength, MeanLength, MedLength, N50Length, L50Count, NCount, GapLength, GapCount, GC
            if sumdb or self.mode().endswith('db'): # {'Name':name,'Sequence':sequence,'FPos':fpos}
                sdb = db.addEmptyTable('sequences',['name','desc','gene','spec','accnum','length'],['name'])
                if self.dna():
                    sdb.addField('gc')      # GC content
                    sdb.addField('basen')    # Number of Ns (ignoring mingap)
                    if mingap:
                        sdb.addField('gapn')    # Number of Ns (applying mingap)
                        sdb.addField('gapnum')  # Number of gaps (applying mingap)
                if self.getBool('Maker'): sdb.addFields(['AED','eAED','QIUTR5','QISpliceEST','QIExonEST','QIExonAlign','QISpliceSNAP','QIExonSNAP','QIExonmRNA','QIUTR3','QIProtLen'])
            gapdb = None
            if gapstats: gapdb = db.addEmptyTable('gaps',['seqname','start','end','seqlen','gaplen'],['seqname','start','end'])
            fracdb = None
            if fracstats or lenstats:
                if raw: fracdb = db.addEmptyTable('fracstats',['file','X','seqlen','seqnum','sumbp'],['file','X','seqlen'])
                else:
                    fracdb = db.addEmptyTable('fracstats',['file','fraction','LXX','NXX','CtgLXX','CtgNXX'],['file','fraction'])
                    if self.getNum('GenomeSize') > 0: fracdb.addFields(['LGXX','NGXX'])
            if self.dna() and mingap: self.printLog('#GAPS','Identifying all runs of %d+ Ns as gaps. (Warning: can be slow for large, gappy genomes.)' % mingap)
            #?# Should we add an option to trim/ignore terminal gaps? #?#
            #self.printLog('#~~#','# ~~~~~~~~~~~~~~~~~~~~~~~ SEQUENCE SUMMARY FOR %s ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #' % basename)
            self.headLog('Sequence Summary for %s' % basename)
            seqlen = []
            ctglen = []
            self.progLog('\r#SUM','Total number of sequences:')
            if not seqs: seqs = self.seqs()
            # Total number of sequences
            for seq in seqs:
                if self.seqNum() > 1000:
                    self.progLog('\r#SUM','Total number of sequences: %s' % rje.iLen(seqlen),rand=0.01)
                else:
                    self.progLog('\r#SUM','Total number of sequences: %s' % rje.iLen(seqlen))
                sname = self.shortName(seq)
                seqlen.append(self.seqLen(seq))
                gapn = 0
                if self.dna():
                    sequence = self.seqSequence(seq).upper()
                    seqdata['GC'] += sequence.count('G') + sequence.count('C')
                    seqdata['NCount'] += sequence.count('N')
                    gaps = []
                    if mingap:
                        gaps = re.findall('N{%d,}' % mingap,sequence)
                    gapn = len(''.join(gaps))
                    seqdata['GapLength'] += gapn
                    seqdata['GapCount'] += len(gaps)
                    gapy = 0
                    if gapn:
                        for m in gapre.finditer(sequence):
                            ctg = m.start() - gapy
                            if ctg: ctglen.append(ctg)
                            gapx = m.start() + 1
                            gapy = m.end()
                            gapl = len(m.group())
                            if gapstats:
                                gapdb.dict['Data'][(sname,gapx,gapy)] = {'seqname':sname,'start':gapx,'end':gapy,'seqlen':len(sequence),'gaplen':gapl}
                                #gapdb.addEntry({'seqname':self.shortName(seq),'start':gapx,'end':gapy,'seqlen':len(sequence),'gaplen':gapl})
                    if gapy < len(sequence):
                        ctglen.append(len(sequence) - gapy)
                if sumdb or self.mode().endswith('db'):
                    entry = {'name':sname,'desc':self.seqDesc(seq),'gene':self.seqGene(seq),
                             'spec':self.seqSpec(seq),'accnum':self.seqAcc(seq),'length':self.seqLen(seq)}
                    if self.dna():
                        nonN = (len(sequence) - sequence.count('N'))
                        if nonN:
                            entry['gc'] = float(sequence.count('G') + sequence.count('C')) / nonN
                        else:
                            entry['gc'] = -1
                            if len(sequence) > 0:
                                self.warnLog('%s is 100%% Ns!' % sname,warntype="alln")
                            else:
                                self.warnLog('%s is zero length!' % sname,warntype="zerolen")
                        entry['basen'] = float(sequence.count('N')) / len(sequence)
                        entry['gapn'] = gapn
                        entry['gapnum'] = len(gaps)
                    if self.getBool('Maker'):
                        name = self.seqName(seq)
                        if rje.matchExp('\sAED:(\S+)',name): entry['AED'] = rje.matchExp('\sAED:(\S+)',name)[0]
                        if rje.matchExp('eAED:(\S+)',name): entry['eAED'] = rje.matchExp('eAED:(\S+)',name)[0]
                        if rje.matchExp('QI:(\S+)',name):
                            qi = string.split(rje.matchExp('QI:(\S+)',name)[0],'|')
                            for field in ['QIUTR5','QISpliceEST','QIExonEST','QIExonAlign','QISpliceSNAP','QIExonSNAP','QIExonmRNA','QIUTR3','QIProtLen']:
                                entry[field] = qi.pop(0)
                    sdb.addEntry(entry)
            self.printLog('#SUM','Total number of sequences: %s' % rje.iLen(seqlen))
            if not seqs: return {}
            if (sumdb or self.mode().endswith('db')) and save: sdb.saveToFile('%s.sequences.tdt' % seqbase)
            if gapstats:
                if save: gapdb.saveToFile('%s.gaps.tdt' % seqbase)
                gapdb.addField('N',evalue=1)
                gapdb.addField('term',evalue=0)
                for gapentry in gapdb.entries():
                    if gapentry['start'] == 1 or gapentry['end'] == gapentry['seqlen']: gapentry['term'] = 1
                gapdb.compress(['gaplen'],{'N':'sum','term':'sum'})
                gapdb.keepFields(['gaplen','N','term'])
                if save: gapdb.saveToFile('%s.gaplen.tdt' % seqbase)

            # Total sequence length
            sumlen = sum(seqlen)
            self.printLog('#SUM','Total length of sequences: %s' % rje.iStr(sumlen))
            seqdata['SeqNum'] = len(seqlen)
            seqdata['CtgNum'] = seqdata['SeqNum'] + seqdata['GapCount']
            seqdata['TotLength'] = sumlen
            # Min, Max
            seqlen.sort()
            ctglen.sort()
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
            else: median = sum(seqlen[(len(seqlen)/2)-1:][:2]) / 2.0
            self.printLog('#SUM','Median length of sequences: %s' % (rje.iStr(median)))
            seqdata['MedLength'] = median
            ## N50 calculation
            n50len = sumlen / 2.0
            n50 = seqlen[0:]
            l50 = 0
            while n50 and n50len > n50[-1]: n50len -= n50.pop(-1); l50 += 1
            if n50:
                self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(n50[-1]))
                seqdata['N50Length'] = n50[-1]
                l50 += 1
                self.printLog('#SUM','L50 count of sequences: %s' % rje.iStr(l50))
                seqdata['L50Count'] = l50
            else:
                raise ValueError('Half sum of sequences not reached!')
                #self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(seqlen[-1]))
                #seqdata['N50Length'] = seqlen[-1]
            ## Contig N50 calculation
            if self.dna():
                self.printLog('#SUM','Total number of contigs: %s' % rje.iStr(seqdata['CtgNum']))
            if self.dna() and seqdata['GapCount'] > 0:
                n50len = sumlen / 2.0
                n50 = ctglen[0:]
                l50 = 0
                while n50 and n50len > n50[-1]: n50len -= n50.pop(-1); l50 += 1
                if n50:
                    self.printLog('#SUM','Contig N50 length of sequences: %s' % rje.iStr(n50[-1]))
                    seqdata['N50Ctg'] = n50[-1]
                    l50 += 1
                    self.printLog('#SUM','Contig L50 count of sequences: %s' % rje.iStr(l50))
                    seqdata['L50Ctg'] = l50
                else:
                    raise ValueError('Half sum of sequences not reached!')
                    #self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(seqlen[-1]))
                    #seqdata['N50Length'] = seqlen[-1]
            elif self.dna():
                seqdata['N50Ctg'] = seqdata['N50Length']
                seqdata['L50Ctg'] = seqdata['L50Count']
            ## NG50 calculation
            ng50len = 0.0
            lg50 = 0
            if self.getNum('GenomeSize') > 0: ng50len = self.getNum('GenomeSize')/2.0
            ng50 = seqlen[0:]
            while ng50 and ng50len > ng50[-1]: ng50len -= ng50.pop(-1); lg50 += 1
            if self.getNum('GenomeSize') > 0 and ng50:
                self.printLog('#SUM','NG50 length of sequences (%s): %s' % (dnaLen(self.getNum('GenomeSize')),rje.iStr(ng50[-1])))
                seqdata['NG50Length'] = ng50[-1]
                lg50 += 1
                self.printLog('#SUM','LG50 count of sequences (%s): %s' % (dnaLen(self.getNum('GenomeSize')),rje.iStr(lg50)))
                seqdata['LG50Count'] = lg50
            elif ng50len:
                self.printLog('#SUM','Half genome size (%s) not reached for NG50!' % (dnaLen(self.getNum('GenomeSize'))))
                #self.printLog('#SUM','N50 length of sequences: %s' % rje.iStr(seqlen[-1]))
                seqdata['NG50Length'] = 0
                seqdata['LG50Count'] = -1
            ## FracStats calculations
            if fracdb and not raw:
                slen = seqlen[0:]    # Scaffolds Sorted small to big
                stot = 0             # Summed scaffold length counter
                snum = 0
                clen = ctglen[0:]    # Contigs Sorted small to big
                ctot = 0            # Summed contig length counter
                cnum = 0
                glen = seqlen[0:]    # Scaffolds Sorted small to big
                gtot = 0            # Summed contig length counter
                gnum = 0
                fXX = 0             # Current fraction
                tXX = 0             # Current fraction target length
                gXX = 0             # Current genome fraction target length
                while fXX <= 100:
                    fentry = {'file':self.getStr('SeqIn'),'fraction':fXX}
                    tXX = sumlen * fXX / 100.0
                    gXX = self.getNum('GenomeSize') * fXX / 100.0
                    # Strip long scaffolds until target size met
                    while slen and stot + slen[-1] < tXX: stot += slen.pop(-1); snum += 1
                    fentry['LXX'] = snum
                    if slen:
                        fentry['NXX'] = slen[-1]
                        fentry['LXX'] += 1
                    elif fXX == 100:
                        fentry['NXX'] = seqlen[0]
                    # Strip long contigs until target size met
                    while clen and ctot + clen[-1] < tXX: ctot += clen.pop(-1); cnum += 1
                    fentry['CtgLXX'] = cnum
                    if clen:
                        fentry['CtgNXX'] = clen[-1]
                        fentry['CtgLXX'] += 1
                    else: fentry['CtgNXX'] = ctglen[0]
                    # Strip long scaffolds until target size met
                    if self.getNum('GenomeSize') > 0:
                        while glen and gtot + glen[-1] < gXX: gtot += glen.pop(-1); gnum += 1
                        if glen:
                            fentry['LGXX'] = gnum + 1
                            fentry['NGXX'] = glen[-1]
                    # Increment
                    fracdb.addEntry(fentry)
                    fXX += fracstep
                if save: fracdb.saveToFile('%s.fracstats.tdt' % seqbase)
            ## Raw FracStats and LenStats calculations
            elif fracdb:
                slen = seqlen[0:]   # Scaffolds Sorted small to big
                stot = 0            # Summed sequence length counter
                snum = 0            # Number of sequences
                lenstats.sort(reverse=True)
                while lenstats:
                    while lenstats and slen and slen[-1] > lenstats[0]: stot += slen.pop(-1); snum += 1
                    #i# 'seqlen','seqnum','sumbp'
                    fentry = {'file':self.getStr('SeqIn'),'X':0,'seqlen':lenstats.pop(0),'sumbp':stot,'seqnum':snum}
                    if self.getNum('GenomeSize') > 0: fentry['X'] = float(stot) / self.getNum('GenomeSize')
                    fracdb.addEntry(fentry)
                xlen = seqlen[0:]   # Scaffolds Sorted small to big
                xcov = fracstep     # Target X coverage
                xsum = 0            # Target summed bases
                xtot = 0            # Summed sequence length counter
                xnum = 0            # Number of sequences
                prevlen = 0         # Previous seq length
                while fracstep and xlen:
                    xsum = self.getNum('GenomeSize') * xcov
                    while xlen and xtot < xsum: prevlen = xlen.pop(-1); xtot += prevlen; xnum += 1
                    fentry = {'file':self.getStr('SeqIn'),'X':xcov,'seqlen':prevlen,'sumbp':xtot,'seqnum':xnum}
                    fracdb.addEntry(fentry)
                    xcov += fracstep
                if save: fracdb.saveToFile('%s.lenstats.tdt' % seqbase,sfdict={'X':4})
            # GC content
            if self.dna():
                if (float(sumlen) - seqdata['GapLength']) > 0:
                    seqdata['GCPC'] = 100.0 * seqdata['GC'] / (float(sumlen) - seqdata['GapLength'])
                else: seqdata['GCPC'] = -1
                self.printLog('#SUM','GC content: %.2f%%' % seqdata['GCPC'])
                self.printLog('#SUM','N bases: %s (%.2f%%)' % (rje.iStr(seqdata['NCount']),100.0 * seqdata['NCount'] / float(sumlen)))
                if mingap:
                    self.printLog('#SUM','Gap (%d+ N) length: %s (%.2f%%)' % (mingap,rje.iStr(seqdata['GapLength']),100.0 * seqdata['GapLength'] / float(sumlen)))
                    self.printLog('#SUM','Gap (%d+ N) count: %s' % (mingap,rje.iStr(seqdata['GapCount'])))
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
    def rename(self,keepsprotgene=False,genecounter=False):   ### Renames sequences and saves them to file before reloading as normal
        '''
        Renames sequences and saves them to file.
        >> keepsprotgene:bool[False] = Whether to keep SwissProt gene if recognised.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('NewDesc'): return self.newDesc()
            #x#if self.getStr('NewAcc').lower() in ['','none']: return False
            newgene = self.getStrLC('NewGene')
            if not self.getStrLC('NewGene'):
                self.str['NewGene'] = self.getStr('NewAcc').lower()
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
                        if newgene and (gene.upper() != gene or not keepsprotgene):
                            gene = self.getStr('NewGene')
                            if gene in ['acc','accnum']: gene = accnum
                        if genecounter: gene = '%s%s' % (gene,rje.preZero(sx,seqx))
                        if acc in ['','none']: name = '%s_%s__%s %s' % (gene,spcode,accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                        if acc in ['','none']: name = '%s_%s__%s %s' % (gene,spcode,accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                        else: name = '%s_%s__%s%s %s %s' % (gene,spcode,self.getStr('NewAcc'),rje.preZero(sx,seqx),accnum,string.join(string.split(rje.chomp(iline[1:]))[1:]))
                    except:
                        accnum = string.split(rje.chomp(iline)[1:])[0]
                        spcode = self.getStr('SpCode').upper()
                        gene = self.getStr('NewGene')
                        if gene in ['acc','accnum']: gene = accnum
                        if genecounter: gene = '%s%s' % (gene,rje.preZero(sx,seqx))
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
    def newDesc(self,newdesc=None,keepname=None):   ### Renames sequences and saves them to file before reloading as normal
        '''
        Renames sequences and saves them to file.
        >> newdesc:dict[None] = Pre-populated dictionary of {shortname:description}
        >> keepname:bool[None] = Whether to keep original sequence name (uses keepname=T/F if None)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not newdesc and not rje.exists(self.getStr('NewDesc')):
                raise IOError('New description file newdesc="%s" not found!' % self.getStr('NewDesc'))
            if keepname == None: keepname = self.getBool('KeepName')
            ## ~ [0a] ~ Make description dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not newdesc:
                newdesc = {}
                for dline in open(self.getStr('NewDesc'),'r').readlines():
                    if not dline: continue
                    descdata = string.split(rje.chomp(dline),maxsplit=1)
                    if not descdata: continue
                    if descdata[0]:
                        if len(descdata) > 1: newdesc[descdata[0]] = descdata[1]
                        else: newdesc[descdata[0]] = ''
            ## ~ [0b] ~ Setup output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']: outfile = '%s.new.fas' % rje.baseFile(self.getStr('SeqIn'))
            else: self.setStr({'SeqOut':'None'})
            rje.backup(self,outfile)
            ### ~ [1] Rename Sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            IN = open(self.getStr('SeqIn'),'r')
            OUT = open(outfile,'a')
            seqx = 0; descx = 0; blankx = 0
            ## ~ [1a] Process sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            iline = IN.readline(); name = None; seq = None
            while iline:
                if iline[:1] == '>':
                    name = rje.chomp(iline)[1:]
                    seqx += 1
                    sname = string.split(name)[0]
                    if sname in newdesc:
                        if newdesc[sname]:
                            if keepname: name = '%s %s' % (sname,newdesc[sname])
                            else: name = newdesc[sname]
                        else: name = sname; blankx += 1
                        descx += 1
                    OUT.write('>%s\n' % name)
                    self.progLog('\r#DESC','Updating %s descriptions' % rje.iStr(seqx))
                else: OUT.write(iline)
                iline = IN.readline()
            self.printLog('\r#DESC','Updated %s descriptions for %s sequences (%s blank) -> %s' % (rje.iStr(descx),rje.iStr(seqx),rje.iStr(blankx),outfile))
            IN.close(); OUT.close()
            self.setStr({'SeqIn':outfile})
        except: self.errorLog('Problem during %s newDesc.' % self); return False   # Rename failed
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
                        if nodup and name in dupcheck:
                            self.dict['Filter']['Duplicate'] += 1
                            if self.getBool('DupErr'): raise ValueError('Duplicate sequence name: "%s". Switch duperr=F to try to ignore.' % name)
                        else:
                            slen = len(sequence)
                            if slen not in sizedict: sizedict[slen] = []
                            sizedict[slen].append(fprev)
                            if nodup: dupcheck.append(name)
                    name = string.split(line[1:])[0]; sequence = ''; fprev = fpos; sx += 1
                else: sequence += line[0:]
                fpos = SEQ.tell()
            if sequence and name:
                if nodup and name in dupcheck:
                    self.dict['Filter']['Duplicate'] += 1
                    if self.getBool('DupErr'): raise ValueError('Duplicate sequence name: "%s". Switch duperr=F to try to ignore.' % name)
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
                        if nodup and name in dupcheck:
                            self.dict['Filter']['Duplicate'] += 1
                            if self.getBool('DupErr'): raise ValueError('Duplicate sequence name: "%s". Switch duperr=F to try to ignore.' % name)
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
            iline = IN.readline(); name = None; seq = None; fullname = None
            while seqx:
                if iline[:1] == '>' or not iline:
                    if grepnr and name and name not in grepnr:
                        goodrev = goodname = ''
                        foundme = False
                        #self.deBug("grep %s %s -B 1 | grep '>'" % (seq,self.getStr('SeqIn')))
                        #!# Possibly too long!
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
            if twopass and not self.getBool('TwoPass'): self.printLog('#SEQNR','Redundancy removal Pass II deactivated (twopass=F).')
            if twopass:
                goodseq = []
                while seqlist:
                    check = seqlist.pop(0)
                    sequences = string.join(seqlist)
                    if self.getBool('RevCompNR') and self.nt(): revseq = rje_sequence.reverseComplement(check,rna=self.rna())
                    else: revseq = ''
                    if check in sequences:
                        matched = rje.matchExp('\s(\S*%s\S*)\s' % check,' %s ' % sequences)[0]
                        if self.getBool('TwoPass'):
                            name = seqdict.pop(check)
                            self.printLog('\r#SEQNR','%s removed: contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                            self.dict['Filter']['NR'] += 1
                        else:
                            name = seqdict[check]
                            self.warnLog('#SEQNR: %s contained within %s; not removed (twopass=F)' % (string.split(name)[0],string.split(seqdict[matched])[0]),warntype='twopass',suppress=True)
                            goodseq.append(check)
                    elif revseq and revseq in sequences:
                        matched = rje.matchExp('\s(\S*%s\S*)\s' % revseq,' %s ' % sequences)[0]
                        if self.getBool('TwoPass'):
                            name = seqdict.pop(check)
                            self.printLog('\r#SEQNR','%s removed: reverse complement contained within %s' % (string.split(name)[0],string.split(seqdict[matched])[0]))
                            self.dict['Filter']['NR'] += 1
                        else:
                            name = seqdict[check]
                            self.warnLog('#SEQNR: %s reverse complement contained within %s; not removed (twopass=F)' % (string.split(name)[0],string.split(seqdict[matched])[0]),warntype='twopass',suppress=True)
                            goodseq.append(check)
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
        except: self.errorLog('Problem during %s seqNR.' % self); return False   # Tidy failed
#########################################################################################################################
    def grepNR(self,seqnr=True):  ### Removes redundancy from sequences and saves them to file before reloading as normal
        '''
        Removes redundancy from sequences and saves them to file before reloading as normal.
        >> seqnr:bool [True] = Whether to use old seqNR() process if fails. Else with raise an error.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fastgrep = True #!# Replace with commandline
            tmpdir = self.getStr('TmpDir')  #'./tmp/'
            maxgrep = 1000
            redseq = {}     # Dictionary of redundant sequences {bad sequence:(good sequence, description)}
            revnr = self.getBool('RevCompNR') and self.nt()     # Whether to also look at reverse complements
            try:
                self.progLog('#GREP','Testing grep...')
                seqx = string.atoi(rje.chomp(os.popen("grep -c '>' %s" % self.getStr('SeqIn')).readlines()[0]))
                self.printLog('#GREP','Identified %s sequences using grep' % rje.iStr(seqx))
            except:
                raise ValueError('grep failure: cannot use grepNR mode')
            if not self.force(): self.printLog('#TMP','Using tmpdir=%s - will reuse *.nr files (force=F).' % self.getStr('TmpDir'))
            else: self.printLog('#TMP','Using tmpdir=%s - will remake *.nr files (force=T).' % self.getStr('TmpDir'))
            if self.nt(): self.printLog('#REVNR','Reverse complement NR checks: %s' % revnr)
            else: self.printLog('#REVNR','Protein NR (dna=F): no reverse complement NR checks.')
            if not fastgrep: self.warnLog('Note: grep and/or head may report some write or broken pipe errors. These are harmless and can be ignored.')
            ## ~ [0a] ~ Check single line sorting and make seqname list ~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            seqname = []    # List of sequence names in order
            prevseq = ''    # Previous sequence
            SEQIN = open(self.getStr('SeqIn'),'r')
            sx = 0.0
            while SEQIN:
                self.progLog('\r#CHECK','Checking sorted fasta: %.2f%%' % (sx/seqx)); sx += 100.0
                sline = SEQIN.readline()
                if not sline: break
                if not sline.startswith('>'):
                    raise ValueError('Expected ">" line: fasta files needs to be one line per sequence')
                sname = string.split(sline)[0][1:]
                seqname.append(sname)
                seq = string.split(SEQIN.readline())[0]
                if prevseq:
                    if len(seq) > len(prevseq):
                        raise ValueError('Expected ">" line: fasta files needs to be one line per sequence')
                    elif seq == prevseq:
                        redseq[sname] = (seqname[-2],'100% identical to')
                        continue
                    #i# Could also check revcomp but I think it would be too slow.
                prevseq = seq
            SEQIN.close()
            self.printLog('\r#CHECK','Checking sorted fasta complete: %s sequences; %s redundant' % (rje.iLen(seqname),rje.iLen(redseq)))
            self.debug('End of a Check.')
            ## ~ [0b] SeqFile Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            rje.mkDir(self,tmpdir,log=True)
            outfile = self.getStr('SeqOut')
            if outfile == self.getStr('SeqIn'): outfile = ''
            if self.getBool('AutoFilter') or outfile.lower() in ['','none']: outfile = '%s.nr.fas' % rje.baseFile(self.getStr('SeqIn'))
            else: self.setStr({'SeqOut':'None'})
            self.dict['Filter']['NR'] = 0   # Counter for NR removal
            ## ~ [0c] Forking Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            forkx = max(1,self.getInt('Forks'))      # Number of forks to have running at one time
            forks = []      # List of active fork PIDs
            forked = {}     # Dictionary of {pid:sname}
            killforks = self.getInt('KillForks')     # Time in seconds to wait after main thread has apparently finished
            killtime = time.time()

            ### ~ [1] Forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SEQIN = open(self.getStr('SeqIn'),'r')
            sx = 0.0    # Sequence counter
            hx = 0      # Head length counter
            fx = 0      # Forking cleanup loop counter
            readingseq = True   # Still reading in sequences from file
            while readingseq or len(forks):
                ## ~ [1a] Read sequence and start new fork ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                while readingseq and (len(forks) < forkx):     # Add more forks
                    if self.debugging(): self.progLog('\r#GREPNR','Checking for redundancy: %s|%d' % (rje.iStr(int(sx/100)),fx))
                    else: self.progLog('\r#GREPNR','Checking for redundancy: %.2f%%' % (sx/seqx))
                    killtime = time.time()  # Reset killtime - still doing stuff

                    sline = SEQIN.readline()
                    if not sline: readingseq = False; break
                    sname = string.split(sline)[0][1:]
                    if not sname: readingseq = False; break
                    seq = string.split(SEQIN.readline())[0]
                    sx += 100.0; hx += 2
                    if sname in redseq: continue    # Sequence already marked as redundant
                    if sname == seqname[0]: continue    # First Sequence - cannot be redundant

                    #self.debug("head -n %d %s | grep -B 1 %s" % (hx,self.getStr('SeqIn'),seq[:maxgrep]))

                    newpid = os.fork()
                    if newpid == 0: # child
                        self.setBool({'Child':True})
                        self.log.stat['Interactive'] = -1
                        self.log.stat['Verbose'] = -1
                        sfile = '%s%s.nr' % (tmpdir,sname)
                        #i# If file exists, use (unless force=T)
                        if rje.exists(sfile) and not self.force(): raise KeyboardInterrupt    # Exit process

                        #!# Add grepnr code
                        goodname = None
                        desc = None
                        #GREP = os.popen("head -n %d %s | grep -B 1 %s 2>&1" % (hx,self.getStr('SeqIn'),seq[:maxgrep]))
                        GREP = os.popen("head -n %d %s | grep -B 1 %s" % (hx,self.getStr('SeqIn'),seq[:maxgrep]))
                        while GREP:
                            gname = string.split(GREP.readline())[0]
                            if not gname: break
                            elif gname[:1] == '>': gname = gname[1:]
                            else: continue
                            gseq = rje.chomp(GREP.readline())
                            if seq in gseq:
                                goodname = gname
                                if gseq == seq: desc = '100% identical to'
                                else: desc = 'a subsequence of'
                                break
                        #?# Add a quick method that gets more pipefail errors without this while?
                        while GREP.readline() and not fastgrep: pass
                        GREP.close()
                        if not goodname:
                            goodname = sname
                            desc = 'Grep error! Self not found'

                        if goodname == sname and revnr:
                            revseq = rje_sequence.reverseComplement(seq,rna=self.rna())
                            #GREP = os.popen("head -n %d %s | grep -B 1 %s 2>&1" % (hx,self.getStr('SeqIn'),revseq[:maxgrep]))
                            GREP = os.popen("head -n %d %s | grep -B 1 %s" % (hx,self.getStr('SeqIn'),revseq[:maxgrep]))
                            while GREP:
                                gline = GREP.readline()
                                if not gline: break
                                gname = string.split(gline)[0]
                                if not gname: break
                                elif gname[:1] == '>': gname = gname[1:]
                                else: continue
                                gseq = rje.chomp(GREP.readline())
                                if revseq in gseq:
                                    goodname = gname
                                    if gseq == revseq: desc = 'RevComp 100% identical to'
                                    else: desc = 'RevComp a subsequence of'
                                    break
                            #?# Add a quick method that gets more pipefail errors without this while?
                            while GREP.readline() and not fastgrep: pass
                            GREP.close()

                        open(sfile,'w').write('%s\n%s\n' % (goodname,desc))

                        raise KeyboardInterrupt    # Exit process
                    elif newpid == -1: self.errorLog('Problem forking %s.' % sname,printerror=False)  # Error!
                    else:
                        forked[newpid] = sname
                        forks.append(newpid)    # Add fork to list

                ## ~ [1b] Monitor and remove finished forks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                time.sleep(0.01)       # Sleep for 1s
                forklist = self._activeForks(forks,nowarn=[0])
                if len(forklist) != len(forks):
                    #self.verbose(1,2,' => %d of %d forks finished!' % (len(forks) - len(forklist),len(forks)),1)
                    forks = forklist[0:]
                    for pid in forked.keys():   # Go through current forks
                        if pid not in forks:
                            sname = forked.pop(pid)
                            sfile = '%s%s.nr' % (tmpdir,sname)
                            fx += 1
                            #i# Process sname.nr file: first line is goodseq; second line (if redundant) is description
                            NR = open(sfile,'r')
                            goodname = rje.chomp(NR.readline())
                            if goodname !=  sname:  # Redundant
                                desc = rje.chomp(NR.readline())
                                redseq[sname] = (goodname,desc)
                            NR.close()
                            killtime = time.time()  # Reset killtime - still doing stuff
                #self.bugPrint(redseq)
                #self.debug('End of a Cycle.')

                ## ~ [1c] Look for eternal hanging of threads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if (time.time() - killtime) > killforks:
                    self.warnLog('%d seconds of main thread inactivity. %d forks still active!' % (killforks,len(forks)))
                    for fork in forks:
                        self.warnLog('GrepNRFork %s, PID %d still Active!' % (forked[fork],fork),1)
                    if self.i() < 0 or rje.yesNo('Kill Main Thread?'): break   #!# killing options
                    elif rje.yesNo('Kill hanging forks?'):
                        for fork in forks:
                            self.printLog('#KILL','Killing GopherFork %s, PID %d.' % (forked[fork],fork))
                            os.system('kill %d' % fork)
                    else: killtime = time.time()

            ### ~ [2] Finish Forking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            SEQIN.close()
            self.printLog('\r#GREPNR','Checking for redundancy complete: %s sequences; %s redundant' % (rje.iLen(seqname),rje.iLen(redseq)))
            if len(forks) > 0: self.errorLog('%d GrepNR Forks still active after %d seconds of mainthread inactivity' % (len(forks),killforks),quitchoice=True,printerror=False)

            ### ~ [3] Save if filtering done ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if redseq:
                self.printLog('\r#SEQNR','%s redundant sequences flagged for removal' % (rje.iLen(redseq)))
                self._filterCmd(clear=True)
                self.list['BadSeq'] = rje.sortKeys(redseq)
                RED = open('%s.redundant.txt' % rje.baseFile(outfile),'w')
                for gname in self.list['BadSeq']:
                    gnr = redseq[gname]
                    goodname = gnr[0]
                    while goodname in redseq: goodname = redseq[goodname][0]
                    RED.write('%s is %s %s\n' % (gname,gnr[1],goodname))
                    self.printLog('#NR','%s is %s %s' % (gname,gnr[1],goodname),screen=False)
                RED.close()
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
            ## ~ [3a] Tidied tmpdir ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            sx = 0.0; stot = len(seqname)
            for sname in seqname:
                self.progLog('\r#TMP','Cleaning up tmpdir: %.2f%%' % (sx/stot)); sx += 100.0
                sfile = '%s%s.nr' % (tmpdir,sname)
                if os.path.exists(sfile): os.unlink(sfile)
            self.printLog('\r#TMP','Cleaning up tmpdir complete!')
            try: os.rmdir(tmpdir)
            except: self.printLog('\r#TMP','Failed to delete tmpdir.')

        except KeyboardInterrupt:
            if self.getBool('Child'): os._exit(0)
            else: raise
        except ValueError:
            if self.getBool('Child'): os._exit(1)
            self.errorLog('Problem during %s grepNR.' % self)
            if seqnr: return self.seqNR(twopass=True,grepnr=False)
            else: return False   # Tidy failed
        except:
            if self.getBool('Child'): os._exit(1)
            self.errorLog('Cannot use grepNR mode.'); return False   # Tidy failed
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
                    seqfile = rje.choice(text='\nNo file name given: Input filename? (Blank to exit.)')
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
                if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
                self.obj['DB'].addEmptyTable(self.getStr('Name'),['Name','Sequence','FPos'],['Name'])
            elif mode == 'filedb':
                if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
                self.obj['DB'].addEmptyTable(self.getStr('Name'),['Name','FPos'],['Name'])

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
            ## ~ [1c] ~ Fastq format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if filetype == 'fastq':
                sequence = None; name = None; fprev = 0
                if mode == 'index' and not makeindex: fpos = fend
                while fpos < fend:
                    self.progLog('\r#SEQ','Loading seq from %s: %.2f%%' % (logseqfile,(100.0*fpos/fend)),rand=0.2,screen=screen)
                    line = rje.chomp(SEQ.readline())
                    if line.find('@') == 0:                                     # New Sequence
                        if sequence and name: self._addSeq(name,sequence,fprev,makeindex,nodup=nodup)  # Previous Sequence to Add
                        name = line[1:]
                        sequence = rje.chomp(SEQ.readline())
                        SEQ.readline()  # +
                        SEQ.readline()  # Quality scores
                        fprev = fpos; sx += 1
                    elif line: raise ValueError('Unexpected fastq line where shoud be sequence name: "%s"' % line)
                    fpos = SEQ.tell()
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
            if nodup and self.dict['Filter']['Duplicate']:
                if self.getBool('DupErr'): raise ValueError('%s+ duplicate sequences. Switch duperr=F to try to ignore.' % rje.iStr(self.dict['Filter']['Duplicate']))
                self.printLog('\r#DUP','%s duplicate sequences ignored' % (rje.iStr(self.dict['Filter']['Duplicate'])),screen=screen)
            if 'Sequence unavailable' in self.dict['Filter']: self.printLog('\r#NOSEQ','%s unavailable sequences ignored' % (rje.iStr(self.dict['Filter']['Sequence unavailable'])),screen=screen)
            #!# Add rest of code #!#
            return True
        except: self.errorLog('%s.loadSeq error' % self.prog()); return False
#########################################################################################################################
    def checkNames(self,gnspacc=False,dupnames=True,warnonly=False):    ### Check sequence names
        '''
        Check sequence names.
        :param gnspacc:bool [False] = Whether to check for Gene_SPEC__AccNum format
        :param dupnames:bool [True] = Whether to check for duplicate short names
        :param warnonly:bool [False] = Whether to simply warn about issues, or raise ValueError
        :return:
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add option to convert using self.extractNameDetails(name)
            names = []
            for seq in self.seqs():
                sname = self.shortName(seq)
                if dupnames and sname in names:
                    if not warnonly: raise ValueError('Duplicated sequence name "%s" detected!' % sname)
                    self.warnLog('Duplicated sequence name "%s" detected!' % sname)
                if gnspacc and not rje.matchExp('^(\S+_\S+__\S+)',sname):
                    if not warnonly: raise ValueError('Sequence "%s" not in SLiMSuite gene_SPEC__AccNum format!' % sname)
                    self.warnLog('Sequence "%s" not in SLiMSuite gene_SPEC__AccNum format!' % sname,warntype='gnspacc',suppress=True)
                names.append(sname)
            return True
        except: self.errorLog('Problem encountered during checkNames()'); raise
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
                       ('fastq','^@(\S+)'),             # Fastq format
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
                        if nodup:
                            self.dict['Filter']['Duplicate'] += 1
                            if self.getBool('DupErr'): raise ValueError('Duplicate sequence name: "%s". Switch duperr=F to try to ignore.' % fkey)
                            return
                        else: self.warnLog('Sequence name "%s" occurs 2+ times!' % fkey)
                    else: self.dict['SeqDict'][fkey] = fpos
                    self.list['Seq'].append(fpos)
                else: self.list['Seq'].append(fpos)     #!# No duplicate checking if index not used #!#
                return
            elif self.mode() == 'index':    # No loading of sequences. Use index file to find sequences on the fly.
                fkey = string.split(name)[0]
                if fkey in self.dict['SeqDict']:
                    if nodup:
                        self.dict['Filter']['Duplicate'] += 1
                        if self.getBool('DupErr'): raise ValueError('Duplicate sequence name: "%s". Switch duperr=F to try to ignore.' % fkey)
                        return
                    else: self.warnLog('Sequence name "%s" occurs 2+ times!' % fkey)
                else: self.dict['SeqDict'][fkey] = fpos
                return
            ### ~ [3] ~ Sequence Processing and filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# At some point, add processing and filtering based on name and sequence #!#
            ### ~ [4] ~ Full ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.mode() == 'full':   # Full loading into Sequence Objects.
                return self.OLD_addSeq(name,sequence)     # Temporary method
            ### ~ [5] ~ Reduced ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            elif self.mode() in ['list','tuple']:     # Lists of (name,sequence) tuples only.
                if nodup and (name,sequence) in self.list['Seq']:
                    self.dict['Filter']['Duplicate'] += 1
                    if self.getBool('DupErr'): raise ValueError('Duplicate sequence: "%s". Switch duperr=F to try to ignore.' % name)
                else: self.list['Seq'].append((name,sequence))
            elif self.mode() == 'db':       # Store sequence data in database object.
                db = self.db(self.getStr('Name'))
                entry = {'Name':name,'Sequence':sequence,'FPos':fpos}
                db.addEntry(entry)  #!# Later, add method to expand/extract data
                self.list['Seq'].append(entry)
            elif self.mode() == 'filedb':       # Store sequence data in database object.
                db = self.db(self.getStr('Name'))
                entry = {'Name':name,'FPos':fpos}
                db.addEntry(entry)  #!# Later, add method to expand/extract data
                self.list['Seq'].append(entry)
        except: self.errorLog('%s._addSeq error (%s)' % (self.prog(),name))
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
            #!# Add extractNameDetails(name) filtering for other formats
            minlen = max(0,self.getInt('MinLen')); maxlen = max(0,self.getInt('MaxLen'))
            #self.debug('Min: %s; Max: %s' % (minlen,maxlen))
            seqfilter = minlen + maxlen
            ### ~ [1] ~ Simple Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            simplefilter = self.mode() == 'file' and self.dict['SeqDict'] and not seqfilter and self.getStr('SeqDictType') in ['short','name']
            if simplefilter:
                sx = 0.0; stot = len(self.dict['SeqDict'])
                for name in self.dict['SeqDict'].keys()[0:]:
                    if screen: self.progLog('\r#FILT','Filtering sequences %.2f%%' % (sx/stot),rand=0.01); sx += 100
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
                    if screen: self.progLog('\r#FILT','Filtering sequences %.2f%%' % (sx/stot),rand=0.01); sx += 100
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
    def grabSeq(self,mask=False):  ### Grabs and/or masks certain regions of sequences. Converts to seqmode=list.
        '''
        Grabs and/or masks certain regions of sequences. Converts to seqmode=list.
        >> mask:bool [False] = whether to mask [True] or extract [False] sequences.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ltype = '#Grab'
            posfile = self.getStr('GrabSeq')
            if mask:
                posfile = self.getStr('MaskSeq')
                ltype = '#Mask'
            if not rje.exists(posfile):
                raise IOError('%sSeq file "%s" not found!' % (ltype[1:],posfile))
            ## ~ [1a] ~ Setup position table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            db = self.db()
            covflanks = self.getInt('AddFlanks')    #i# Negative value will shrink regions, not grow them
            if not len(self.list['PosFields']) == 3:
                raise ValueError('posfields=LIST must have exactly 3 elements: Locus, Start, End. %d found!' % len(self.list['PosFields']))
            [locusfield,startfield,endfield] = self.list['PosFields']
            #cdb = db.addTable(self.getStr('CheckPos'),mainkeys=self.list['CheckFields'],name='check',expect=True)
            cdb = db.addTable(posfile,mainkeys=self.list['PosFields'],name='pos',expect=True,ignore=[])
            cdb.dataFormat({startfield:'int',endfield:'int'})
            if not cdb.entryNum():
                self.warnLog('No %sSeq positions identified to %s' % (ltype[1:],ltype[1:].lower()))
                return False
            ## ~ [1b] Setup new sequence list ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            newseq = []     # List of (seqname, sequence) to replace self.list['Seq']
            if mask and self.getStr('SeqMode') != 'list':
                for seq in self.seqs(): newseq.append(self.getSeq(seq))
                self.list['Seq'] = newseq
                self.setStr({'SeqMode':'list'})
                self.printLog(ltype.upper(),'%s mode detected: seqmode=list' % ltype[1:])
            seqdict = self.seqNameDic()     # Will map Locus onto shortname onto sequence using seqdict

            ### ~ [2] Grab/Mask ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cx = 0.0; ctot = cdb.entryNum(); masked = []; nullx = 0
            for centry in cdb.entries():
                self.progLog('\r%s' % ltype,'Processing %sSeq positions: %.2f%%' % (ltype[1:],cx/ctot)); cx += 100.0
                locus = centry[locusfield]
                if locus not in seqdict:
                    nullx += 1
                    continue
                (sname, sequence) = self.getSeq(seqdict[locus])
                cstart = min(centry[startfield],centry[endfield])   # Can accept either direction
                cend = max(centry[startfield],centry[endfield])   # Can accept either direction
                cstart = max(1,cstart-covflanks)
                seqlen = len(sequence)
                cend = min(seqlen,cend+covflanks)
                if cend < cstart:
                    self.warnLog('Region skipped: flanked end position before flanked start position')
                    continue
                seqi = self.seqNum()    # This will be replaced will masked sequence
                if mask:
                    seqi = self.list['Seq'].index((sname, sequence))
                    sequence = sequence[:cstart-1] + 'N' * (cend-cstart+1) + sequence[cend:]
                    if len(sequence) != seqlen: raise ValueError('Masking problem!')
                    sname = '%s (masked %s-%s)' % (sname,rje.iStr(cstart),rje.iStr(cend))
                    self.list['Seq'][i] = (seqname, sequence)
                    masked.append(i)
                else:
                    sequence = sequence[cstart-1:cend]
                    sname = '%s (region %s-%s)' % (sname,rje.iStr(cstart),rje.iStr(cend))
                    newseq.append((sname, sequence))
            self.printLog('\r%s' % ltype.upper(),'Processing %sSeq positions complete!' % (ltype[1:]))
            if nullx: self.warnLog('%s regions did not match known sequence. (Checking naming and filtering.)' % rje.iStr(nullx))
            if mask:
                regx = len(masked); maskx = len(rje.sortUnique(masked))
                self.printLog('\r%s' % ltype.upper(),'%s regions in %s sequences masked.' % (rje.iStr(regx),rje.iStr(maskx)))
            ## Update for grab ##
            if not mask:
                self.list['Seq'] = newseq
                if self.getStr('SeqMode') != 'list':
                    self.setStr({'SeqMode':'list'})
                self.printLog(ltype.upper(),'%s regions extracted: seqmode=list' % (rje.iLen(newseq)))
            return True

        except: self.errorLog('Major error during grabSeq()')
        return False
#########################################################################################################################
    def seqFromBlastDBCmd(self,id,dbase=None,expect=True,add=False):  ### Returns sequence object(s) of sequence as obtained with blastdbcmd
        '''
        Returns sequence object of sequence as obtained with blastdbcmd.
        >> id:str = id (or list of ids) of sequence to pass to blastdbcmd
        >> dbase:str = formatted database
        >> expect:bool [True] = whether sequences are expected to be returned.
        '''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dev(): raise ValueError('seqFromBlastDBCmd() not yet implemented. Contact author for alternative processing options.')
            cmdseq = []
            try: id = id.split(); singleseq = True
            except: singleseq = False    # Already a list of IDs
            id = string.join(id,',')
            if not dbase: dbase = self.getStr('SeqIn')
            blastpath = rje.makePath(self.getStr('BLAST+ Path')) + 'blastdbcmd'
            if self.getBool('Win32'): BLASTDBCMD = os.popen("%s -entry \"%s\" -db %s -long_seqids" % (blastpath,id,dbase))
            else:
                self.debug('%s -entry "%s" -db %s -long_seqids' % (blastpath,id,dbase))
                BLASTDBCMD = os.popen('%s -entry "%s" -db %s -long_seqids' % (blastpath,id,dbase))
            ### ~ [2] Extract details and add sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lastline = 'BLASTDBCMD'; sx = 0
            (blastseq,lastline) = self.nextFasSeq(BLASTDBCMD,lastline,raw=True)
            while blastseq:
                (name,sequence) = blastseq; sx += 1
                self.bugPrint(name)
                if name.startswith('lcl|'): name = name[4:]
                if add: cmdseq.append(self._addSeq(name=name, sequence=sequence))
                else: cmdseq.append((name,sequence))
                (blastseq,lastline) = self.nextFasSeq(BLASTDBCMD,lastline,raw=True)
            BLASTDBCMD.close()
            ### ~ [3] Return sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.debug('%d seqs' % len(cmdseq))
            if cmdseq and singleseq: return cmdseq[0]
            elif cmdseq: return cmdseq
            else:
                if expect: self.errorLog('No sequences extracted from %s. (%d returned by blastdbcmd)' % (dbase,sx),printerror=False)
                return []
        except:
            self.errorLog('Major error during seqFromBlastDBCmd()')
            raise
#########################################################################################################################
    def nextFasSeq(self,fileobject=None,lastline=None,raw=True):  ### Returns sequence object of next sequence from passed File Object
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
     ### <5> ### Class Saving Methods                                                                                   #
#########################################################################################################################
    def saveRegions(self,seqs=None,seqfile=None,append=None,log=True,screen=None,backup=True):
        '''
        Saves sequences in SeqList object in fasta format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> reformat:str = Type of output format (fasta/short/acc/acclist/speclist/peptide/qregion/region)
        >> append:boolean [None] = append, do not overwrite, file. Use self.getBool('Append') if None
        >> log:boolean [True] = Whether to log output
        >> screen:bool [None] = Whether to print log output to screen (None will use log setting)
        >> backup:bool [True] = Whether to ask about backing up existing file.
        >> seqtuples:bool [False] = Whether given sequences are tuples to be used directly (otherwise will fetch from self.seq)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if seqs is None: seqs = self.seqs()
            seqdict = self.makeSeqNameDic('max')
            if not seqfile: seqfile = self.getStr('SeqOut')
            if screen == None: screen = log
            if append == None: append = self.getBool('Append')
            # Make region database Table if not existing
            if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            db = self.db()
            regdb = self.db('region')
            if not regdb: regdb = db.addTable(filename=self.getStr('Region'),mainkeys='auto',name='region',expect=True)
            if 'NewAcc' in regdb.fields(): regdb.newKey(['NewAcc'])
            else:
                self.warnLog('NewAcc not in Region file fields: will use unique Locus,Start.End combos.')
                regdb.newKey(['Locus','Start','End'])
            regdb.dataFormat({'Start':'int','End':'int'})
            # add findField() method to look for field and ask if missing and i>=0
            ### ~ [1] ~ Generate Sequence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            # Use regiondb to make new sequence and output
            if backup and not append: rje.backup(self,seqfile)
            if append: SEQOUT = open(seqfile,'a')
            else: SEQOUT = open(seqfile,'w')
            fasx = 0
            loci = rje.listIntersect(regdb.indexKeys('Locus'),seqdict.keys())
            self.printLog('#REGLOC','%s of %s region loci mapped to sequences' % (rje.iLen(loci),rje.iLen(regdb.index('Locus'))),log=log,screen=screen)
            for rkey in regdb.dataKeys():
                entry = regdb.data(rkey)
                if entry['Locus'] not in loci: continue
                seq = seqdict[entry['Locus']]
                if seq not in seqs: continue
                (name,sequence) = self.getSeq(seq)
                seqlen = len(sequence)
                begx = entry['Start']
                if begx < 0: begx = max(1,seqlen-begx+1)
                endx = entry['End']
                if endx < 0: endx = min(seqlen,endx)
                if 'NewAcc' in regdb.fields(): newacc = entry['NewAcc']
                else: newacc = '%s.%s-%s' % (self.seqAcc(seq),rje.preZero(begx,seqlen),rje.preZero(endx,seqlen))
                newname = '%s_%s__%s %s Region %s-%s; %s' % (self.seqGene(seq),self.seqSpec(seq),newacc,self.seqAcc(seq),rje.iStr(begx),rje.iStr(endx),self.seqDesc(seq))
                SEQOUT.write('>%s\n%s\n' % (newname,sequence[begx-1:endx]))
                fasx += 1
            SEQOUT.close()
            self.printLog('\r#OUT','%s regions (%s loci) output to %s' % (rje.iStr(fasx),rje.iLen(loci),seqfile),log=log,screen=screen)
        except: self.errorLog('Major error during saveRegions()')
#########################################################################################################################
    def saveSeq(self,seqs=None,seqfile=None,reformat=None,append=None,log=True,screen=None,backup=True,seqtuples=False):  ### Saves sequences in fasta format
        '''
        Saves sequences in SeqList object in fasta format
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> reformat:str = Type of output format (fasta/short/acc/acclist/speclist/peptide/qregion/region)
        >> append:boolean [None] = append, do not overwrite, file. Use self.getBool('Append') if None
        >> log:boolean [True] = Whether to log output
        >> screen:bool [None] = Whether to print log output to screen (None will use log setting)
        >> backup:bool [True] = Whether to ask about backing up existing file.
        >> seqtuples:bool [False] = Whether given sequences are tuples to be used directly (otherwise will fetch from self.seq)
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqtuples and seqs is None: seqs = self.seqs()
            if not seqfile: seqfile = self.getStr('SeqOut')
            if not seqs:
                if self.getStrLC('Rest'):
                    if reformat in self.dict['Output']: self.dict['Output'][reformat] = 'ERROR: No input sequences.'
                    else: self.dict['Output']['seqout'] = 'ERROR: No input sequences.'
                return self.errorLog('Cannot output to "%s": No sequences.' % seqfile,printerror=False)
            if screen == None: screen = log
            if append == None: append = self.getBool('Append')
            if reformat == None:
                if self.getStrLC('ReFormat'): reformat = self.getStr('ReFormat')
                else: reformat = 'fasta'
            if reformat == 'speclist': speclist = []
            if reformat in ['dna2prot','rna2prot','translate','nt2prot'] and not self.nt():
                self.readSeqType()
                if not self.nt(): self.warnLog('Trying to translate possible non-nucleotide sequence (dna=F).')
            ifile = '%s.index' % rje.baseFile(seqfile)
            if rje.exists(ifile):
                self.warnLog('%s deleted.' % ifile)
                os.unlink(ifile)
            if self.getStrLC('Region') and rje.exists(rje.makePath(self.getStr('Region'),wholepath=True)):
                return self.saveRegions(seqs,seqfile,append,log,screen,backup)
            mingap = max(1,self.getInt('MinGap'))
            gapre = re.compile('N{%d,}' % mingap)
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
            sx = 0.0; stot = len(seqs); fasx = 0
            for seq in seqs:
                if screen: self.progLog('\r#OUT','Sequence output: %.2f%%' % (sx/stot)); sx += 100.0
                if seqtuples: (name,sequence) = seq
                else: (name,sequence) = self.getSeq(seq,format='tuple')
                slen = len(sequence); startchop = 0; endchop = 0
                if qregion:
                    startchop = qstart - sequence[:qstart].count('-')   # Positions lost
                    if qend != -1:
                        endchop = slen - len(sequence[:qend]) - sequence[qend:].count('-')   # Positions lost
                        sequence = sequence[qstart:qend]
                    else: sequence = sequence[qstart:]
                if reformat[:3] in ['fas']: SEQOUT.write('>%s\n%s\n' % (name,sequence)); fasx += 1
                elif reformat == 'peptides':
                    sequence = string.replace(sequence,'-','')
                    if sequence: SEQOUT.write('%s\n' % (sequence)); fasx += 1
                    else: self.warnLog('No peptide for %s (%d -> %d): 100%% gaps' % (string.split(name)[0],qstart+1,qend))
                elif reformat == 'short': SEQOUT.write('>%s\n%s\n' % (string.split(name)[0],sequence)); fasx += 1
                elif reformat.endswith('region'):
                    sstart = startchop + 1
                    send = self.aaLen(seq) - endchop
                    SEQOUT.write('>%s.%s-%s\n%s\n' % (string.split(name)[0],rje.preZero(sstart,slen),rje.preZero(send,slen),sequence)); fasx += 1
                elif reformat[:3] in ['acc','spe']:
                    try: (gene,spec,acc) = rje.matchExp('^(\S+)_(\S+)__(\S+)',name)
                    except: acc = string.split(name)[0]; gene = 'seq'; spec = 'UNKSP'
                    if reformat in ['acc','accfas']: SEQOUT.write('>%s\n%s\n' % (acc,sequence)); fasx += 1
                    elif reformat == 'accdesc':
                        try:
                            desc = string.split(name,maxsplit=1)[1]
                            SEQOUT.write('>%s %s\n' % (acc,desc)); fasx += 1
                        except:
                            SEQOUT.write('>%s\n' % (acc)); fasx += 1
                        SEQOUT.write('%s\n' % (sequence)); fasx += 1
                    elif reformat == 'acclist': SEQOUT.write('%s\n' % (acc)); fasx += 1
                    elif reformat == 'speclist':
                        if spec not in speclist: speclist.append(spec)
                elif reformat in ['dna2prot','rna2prot','translate','nt2prot']:
                    #SEQOUT.write('>%s\n%s\n' % (name,rje_sequence.dna2prot(sequence)))
                    pfasta = self.dna2protFasta(name,sequence)
                    fasx = fasx + (pfasta.count('\n')/2) - 1
                    SEQOUT.write(pfasta); fasx += 1
                elif reformat in ['revcomp']:
                    #SEQOUT.write('>%s\n%s\n' % (name,rje_sequence.dna2prot(sequence)))
                    name = string.split(name)
                    name.insert(1,'RevComp')
                    name = string.join(name)
                    sequence = rje_sequence.reverseComplement(sequence,rna=self.rna())
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fasx += 1
                elif reformat in ['reverse']:
                    #SEQOUT.write('>%s\n%s\n' % (name,rje_sequence.dna2prot(sequence)))
                    name = string.split(name)
                    name.insert(1,'Reversed')
                    name = string.join(name)
                    sequence = rje.strReverse(sequence)
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fasx += 1
                elif reformat == 'descaffold':
                    sequence = sequence.upper()
                    seqi = 0; seqj = 0
                    for m in gapre.finditer(sequence):
                        seqj = m.start()
                        contig = sequence[seqi:seqj]
                        cname = string.split(name) + ['(Contig %s..%s)' % (seqi+1,seqj)]
                        cname[0] = '{0}.{1}-{2}'.format(cname[0],seqi+1,seqj)
                        cname = ' '.join(cname)
                        if contig:
                            SEQOUT.write('>%s\n%s\n' % (cname,contig)); fasx += 1
                        seqi = m.end()
                    if seqi:
                        seqj = len(sequence)
                        contig = sequence[seqi:seqj]
                        cname = string.split(name)
                        cname[0] = '{0}.{1}-{2}'.format(cname[0],seqi+1,seqj)
                        cname = ' '.join(cname)
                        SEQOUT.write('>%s\n%s\n' % (cname,contig)); fasx += 1
                    elif not seqj:
                        SEQOUT.write('>%s\n%s\n' % (name,sequence)); fasx += 1
                else:
                    self.errorLog('Cannot recognise reformat type "%s": changing to Fasta' % reformat)
                    reformat = 'fasta'
                    SEQOUT.write('>%s\n%s\n' % (name,sequence)); fasx += 1
            ### ~ [3] ~ Close files and tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if reformat == 'speclist': SEQOUT.write(string.join(speclist+[''],'\n')); fasx = stot
            SEQOUT.close()
            self.printLog('\r#OUT','%s Sequences %s %s' % (rje.iStr(fasx),outlog,seqfile),log=log,screen=screen)
            if reformat == 'descaffold':
                self.printLog('\r#CONTIG','%s contigs generated from %s scaffolds' % (rje.iStr(fasx),rje.iStr(stot)))
            
        except(IOError): self.errorLog("Cannot create %s" % seqfile); raise
        except: self.errorLog("Problem saving sequences to %s" % seqfile); raise
#########################################################################################################################
    def tileSeq(self,seqs=None,seqfile=None,append=None,log=True,screen=None,backup=True,seqtuples=False):  ### Saves sequences in fasta format
        '''
        Tiles sequences in SeqList object and saves in fasta format (unless seqfile=False).
        >> seqs:list of Sequence Objects (if none, use self.seq)
        >> seqfile:str [self.info['Name'].fas] = filename
        >> append:boolean [None] = append, do not overwrite, file. Use self.getBool('Append') if None
        >> log:boolean [True] = Whether to log output
        >> screen:bool [None] = Whether to print log output to screen (None will use log setting)
        >> backup:bool [True] = Whether to ask about backing up existing file.
        >> seqtuples:bool [False] = Whether given sequences are tuples to be used directly (otherwise will fetch from self.seq)
        << Returns True/False or updated tuples if seqtuples=True or seqfile=False
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not seqtuples and seqs is None: seqs = self.seqs()
            if not seqfile and seqfile != False:
                if self.getStrLC('SeqOut'): seqfile = self.getStr('SeqOut')
                else: seqfile = '%s.tiled.fasta' % self.baseFile()
            if not seqs:
                if self.getStrLC('Rest'):
                    if reformat in self.dict['Output']: self.dict['Output'][reformat] = 'ERROR: No input sequences.'
                    else: self.dict['Output']['seqout'] = 'ERROR: No input sequences.'
                return self.errorLog('Cannot output to "%s": No sequences.' % seqfile,printerror=False)
            if screen == None: screen = log
            if append == None: append = self.getBool('Append')
            ifile = '%s.index' % rje.baseFile(seqfile)
            if rje.exists(ifile):
                self.warnLog('%s deleted.' % ifile)
                os.unlink(ifile)
            ## ~ [0a] ~ Setup tilinng ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            tilelen = self.getInt('Tile')
            if tilelen < 1: raise ValueError('Cannot tile sequences with tile length < 1!')
            tilestep = self.getInt('TileStep')
            if tilelen + tilestep < 1: raise ValueError('Cannot backtrack tilestep=INT further than tile length (tile=INT)!')
            mintile = self.getNum('MinFile')
            if mintile < 1: mintile = max(1,int(0.5+mintile*tilelen))
            tilename = self.getStrLC('TileName')
            if tilename not in ['start','num','count','pos','purepos','purestart','purenum']:
                self.warnLog('\r#TILE','Tiling name setting tilename={0} not recognised: setting tilename=pos.' % (tilename))
                tilename = 'pos'
                self.setStr({'TileName':tilename})
            ## ~ [0b] ~ Open output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            outlog = 'output to'
            if seqfile == False:
                outlog = 'generated (no output)'
                #self.printLog('#OUT','Append sequences to "%s"' % seqfile,log=log,screen=screen)
                SEQOUT = None
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
            sx = 0.0; stot = len(seqs); fasx = 0
            tiled = []  # Filled in with tuples if seqtuples or seqfile=False
            self.printLog('\r#TILE','Tiling %s sequences (%s tiles, %s shifted, min size=%s)' % (rje.iStr(stot),dnaLen(tilelen),dnaLen(tilestep),dnaLen(mintile)),log=log,screen=screen)
            fullx = 0   # Number of full-length sequences
            for seq in seqs:
                if screen: self.progLog('\r#OUT','Tiling sequences: %.2f%%' % (sx/stot)); sx += 100.0
                if seqtuples: (name,sequence) = seq
                else: (name,sequence) = self.getSeq(seq,format='tuple')
                seqlen = len(sequence)
                maxtile = int(seqlen / float(tilelen+tilestep)) + 1
                tilei = 0; tilex = 0
                while tilei < seqlen:
                    #tilej = min(tilei + tilelen,seqlen)
                    tilej = tilei + tilelen
                    remains = sequence[tilej+tilestep:]
                    if len(remains) < mintile: tilej = seqlen     # Extend tile to end of sequence
                    tileseq = sequence[tilei:tilej]
                    tilex += 1
                    # Add postion to description
                    tname = string.split(name) + ['(Tile %d: %s..%s of %s)' % (tilex,tilei+1,min(tilej,seqlen),seqlen)]
                    # Add accnum suffix
                    if tilei == 0 and tilej >= seqlen:
                        fullx += 1    # No need to modify name as this is the entire sequence
                        tname[-1] = '(Full-length: %s bp)' % seqlen
                    elif tilename == 'purepos':
                        tname[0] = '{0}.{1}-{2}'.format(tname[0],tilei+1,min(tilej,seqlen))
                    elif tilename == 'pos':
                        tname[0] = '{0}.{1}-{2}'.format(tname[0],rje.preZero(tilei+1,seqlen),rje.preZero(min(tilej,seqlen),seqlen))
                    elif tilename == 'purestart':
                        tname[0] = '{0}.{1}'.format(tname[0],tilei+1)
                    elif tilename == 'start':
                        tname[0] = '{0}.{1}'.format(tname[0],rje.preZero(tilei+1,seqlen))
                    elif tilename in ['num','count']:
                        tname[0] = '{0}.{1}'.format(tname[0],rje.preZero(tilex,maxtile))
                    elif tilename == 'purenum':
                        tname[0] = '{0}.{1}'.format(tname[0],tilex)
                    # Regenerate name
                    tname = ' '.join(tname)
                    # Save/Output sequence
                    if SEQOUT:
                        SEQOUT.write('>%s\n%s\n' % (tname,tileseq)); fasx += 1
                        self.deBug('>%s\n%s\n' % (tname,tileseq))
                    if seqtuples or not SEQOUT:
                        tiled.append((tname,tileseq))
                    if tilei >= tilej + tilestep: break
                    tilei = tilej + tilestep
                    if tilei < 1: break
            ### ~ [3] ~ Close files and tidy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if SEQOUT:
                SEQOUT.close()
                self.printLog('\r#OUT','%s Sequences %s %s (%s full length)' % (rje.iStr(fasx),outlog,seqfile,rje.iStr(fullx)),log=log,screen=screen)
            else:
                self.printLog('\r#OUT','%s Sequences %s (%s full length)' % (rje.iStr(fasx),outlog,rje.iStr(fullx)),log=log,screen=screen)

        except(IOError): self.errorLog("Cannot create %s" % seqfile); raise
        except: self.errorLog("Problem tiling sequences (seqfile=%s)" % seqfile); raise
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
            return splits
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
            #X# Not any more: if terminorf < 0: terminorf = minorf
            rforfs = {}; longorf = False; fullorfs = {}
            for frame in rje.sortKeys(rfseq):
                self.progLog('\r#ORF','Sequence RF%s ORFs...   ' % frame)
                osplit = string.split(string.replace(rfseq[frame].upper(),'*','*|'),'|')
                orfs = []
                for orf in osplit[:-1]:
                    if self.getBool('ORFMet'):  # Trim to Nterminal Met
                        if 'M' in orf: orf = orf[orf.find('M'):]
                        else: continue
                    if len(orf) < (minorf + 1): continue
                    orfs.append(orf); longorf = True    # More than termini - internal ORFs meet length requirement
                fullorfs[frame] = orfs[0:]
                if terminorf >= 0 and not longorf:
                    if len(osplit) == 1:    # Full read-through
                        if len(osplit[0]) > terminorf: orfs.append(osplit[0])
                    else:
                        if len(osplit[0]) >= terminorf: orfs.append(osplit[0])
                        if len(osplit[-1]) > terminorf: orfs.append(osplit[-1])
                rforfs[frame] = orfs

            for frame in rje.sortKeys(rfseq):
                self.progLog('\r#ORF','Sequence ORF Fasta...   ')
                if longorf: orfs = fullorfs[frame]
                else: orfs = rforfs[frame]
                # Remaining ORFs should meet length requirements
                for i in range(len(orfs)):
                    tranfas += '>%s\n%s\n' % (string.join([namesplit[0]+'.RF%d.ORF%d' % (frame,i+1)]+namesplit[1:]+['Length=%d' % len(orfs[i])]),orfs[i])
            return tranfas
        except: self.errorLog("Problem with dna2protFasta()"); raise
#########################################################################################################################
    def stripGap(self,stripgap=0,codons=False,gaps='-'): ### Removes columns containing gaps from the alignment
        '''
        Removes columns containing gaps from the alignment.
        >> stripgap:num [0] = Number of sequences with gaps before stripping from alignment. Proportion if < 1.
        >> codons:bool [False] = Whether to treat alignment using a codon model (i.e. strip sets of three bases)
        >> gaps:list ['-'] = Characters to recognise as gaps (e.g. could add 'X' for tidyXGap equivalent)
        << tuplist:list of (name,sequence) tuples with gapped columns stripped.
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
     ### <6> ### Menu-based Sequence Editing                                                                            #
#########################################################################################################################
    def edit(self,save=True): ### Gives menu-based sequence editing options. Saves to *.edit.fas by default.
        '''Gives menu-based sequence editing options. Saves to *.edit.fas by default.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Edit') and rje.exists(self.getStr('Edit')):
                edited = self.fileEdit()
                if not rje.yesNo('Continue to regular interactive edit menu?'): return edited
            headtext = 'Sequence edit menu\n(NOTE: in development)'
            self.dict['SeqDict'] = {}
            self.setStr({'EditFile':'%s.edit.fas' % rje.baseFile(self.getStr('SeqIn'))})
            # List of menu item tuples (edit code,description,optiontype,optionkey)
            if self.dna(): dna_menu = [('R','<R>everse complement','return','R'),
                                       ('T','<T>ranslate','return','T'),
                                       ('V','Re<V>erse (not complemented)','return','V')]
            else: dna_menu = [('DNA','Switch to <DNA> mode','DNA')]
            menulist = [('','# ~ NAME/DESCRIPTION ~ #','',''),
                        ('G','Edit <G>ene','return','G'),
                        ('S','Edit <S>pecies','return','S'),
                        ('A','Edit <A>ccession','return','A'),
                        ('D','Edit <D>escription','return','D'),
                        ('','# ~ SEQUENCE ~ #','',''),
                        ('U','<U>ngap','return','U'),
                        ('M','<M>utate','return','M'),
                        ('/','Divide sequence','return','/'),
                        ('X','E<x>tract subsequence/trim','return','X')] + dna_menu
            menulist += [('','# ~ ORGANISATION ~ #','',''),
                         ('<','Move up sequence list','return','<'),
                         ('>','Move down sequence list','return','>'),
                         ('<<','Move to start of sequence list','return','<<'),
                         ('>>','Move to end of sequence list','return','>>'),
                         ('-','Remove sequence from list','return','-'),
                         ('--','Remove this and following sequences from list','return','--'),
                         ('+','Add sequences from file','return','+'),
                         ('++','Duplicate sequence','return','++'),
                         ('&','Append (and remove) following sequence(s)','return','&'),
                         ('C','Combine with following sequence(s) to create <C>onsensus','return','C'),
                         ('','# ~ MENU ~ #','',''),
                         ('L','Add note to <L>og file','return','L'),
                         ('P','<P>revious sequence','return','P'),
                         ('N','<N>ext sequence','return','N'),
                         ('J','<J>ump to sequence','return','J'),
                         ('?','Find sequence','return','?'),
                         ('F','Output <F>ilename','str','EditFile'),
                         ('O','<O>utput sequences to file','return','O'),
                         ('Q','<Q>uit edit menu','return','Q')]
            ### ~ [1] ~ Cycle through menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ei = 0      # Sequence index position
            while True and self.list['Seq']:
                ## ~ [1a] Setup sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.obj['Current'] = self.list['Seq'][ei]
                (seqname,sequence) = self.current()
                seqacc = self.seqAcc()
                gnspacc = self.gnSpAcc(seqname) or self.getBool('GeneSpAcc')
                self.deBug('%s -> %s' % (seqname,self.gnSpAcc(seqname)))
                #!# Add wrapping of seqDesc
                edesc = '\n\n%d of %d: %s (%s %s)\n%s\n' % (ei+1,self.seqNum(),self.shortName(),self.seqLen(),self.units(),self.screenWrap(self.seqDesc(),prefix='    '))
                choice = rje_menu.menu(self,headtext+edesc,menulist,choicetext='Please select:',changecase=True,default='N')
                ## ~ [1b] Name/Description ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if choice == 'G':
                    newgene = rje.choice('New Gene (no spaces) for %s?' % self.shortName(),default=self.seqGene())
                    newgene = string.join(string.split(newgene),'')
                    newshort = '%s_%s__%s' % (newgene,self.seqSpec(),self.seqAcc())
                    self.printLog('#EDIT','%s -> %s' % (self.shortName(),newshort))
                    self.list['Edit'].append(('Gene',seqacc,'%s -> %s' % (self.shortName(),newshort)))
                    self.list['Seq'][ei] = ('%s %s' % (newshort,self.seqDesc()),sequence)
                elif choice == 'S':
                    spcode = rje.choice('New species code (no spaces) for %s?' % self.shortName(),default=self.seqSpec())
                    spcode = string.join(string.split(spcode),'').upper()
                    newshort = '%s_%s__%s' % (self.seqGene(),spcode,self.seqAcc())
                    self.printLog('#EDIT','%s -> %s' % (self.shortName(),newshort))
                    self.list['Edit'].append(('Species',seqacc,'%s -> %s' % (self.shortName(),newshort)))
                    self.list['Seq'][ei] = ('%s %s' % (newshort,self.seqDesc()),sequence)
                elif choice == 'A':
                    newacc = rje.choice('New accession (no spaces) for %s?' % self.shortName(),default=self.seqAcc())
                    newacc = string.join(string.split(newacc),'')
                    if gnspacc: newshort = '%s_%s__%s' % (self.seqGene(),self.seqSpec(),newacc)
                    else: newshort = newacc
                    self.printLog('#EDIT','%s -> %s' % (self.shortName(),newshort))
                    self.list['Edit'].append(('AccNum',seqacc,'%s -> %s' % (self.shortName(),newshort)))
                    self.list['Seq'][ei] = ('%s %s' % (newshort,self.seqDesc()),sequence)
                elif choice == 'D':
                    newdesc = rje.choice('New description for %s?' % self.shortName(),default=self.seqDesc())
                    self.printLog('#EDIT','"%s" -> "%s"' % (self.seqDesc(),newdesc))
                    self.list['Edit'].append(('Desc',seqacc,'"%s" -> "%s"' % (self.seqDesc(),newdesc)))
                    self.list['Seq'][ei] = ('%s %s' % (self.shortName(),newdesc),sequence)
                ## ~ [1c] Sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif choice == 'U':
                    gapx = sequence.count('-')
                    if not gapx: self.vPrint('No gaps to remove!')
                    else:
                        self.list['Seq'][ei] = (seqname,sequence.replace('-',''))
                        self.printLog('#EDIT','%s gaps removed from %s' % (gapx,self.shortName()))
                        self.list['Edit'].append(('Ungap',seqacc,'%s gaps removed from %s' % (gapx,self.shortName())))
                elif choice == 'X':
                    x = rje.getInt('New start position (1-L)?',default=1)
                    y = rje.getInt('New end position (1-L)?',default=len(sequence))
                    newseq = sequence[x-1:y]
                    if not rje.yesNo('Extract subsequence %d to %d (=%d %s)?)' % (x,y,len(newseq),self.units())): continue
                    seqname = '%s (Region %d to %d)' % (seqname,x,y)
                    self.printLog('#EDIT','%s -> Region %d to %d.' % (seqname,x,y))
                    self.list['Edit'].append(('Subsequence',seqacc,'%s -> Region %d to %d.' % (seqname,x,y)))
                    self.list['Seq'][ei] = (seqname,newseq)
                elif choice == '/':
                    divc = rje.getInt('Start position for  new sequence (1-L)?',default=1)
                    if divc > 0: divi = divc - 1
                    else: divi = divc
                    newseq = sequence[:divi]
                    addseq = sequence[divi:]
                    if newseq and addseq:
                        divi = len(newseq)
                        if rje.yesNo('Split %s into %d and %d %s sequences?' % (self.shortName(),len(newseq),len(addseq),self.units())):
                            newacc = rje.choice('New accession (no spaces) for %s %d+?' % (self.shortName(),divi),default='%s.%d' % (self.seqAcc(),divi+1))
                            newacc = string.join(string.split(newacc),'')
                            if gnspacc: newshort = '%s_%s__%s' % (self.seqGene(),self.seqSpec(),newacc)
                            else: newshort = newacc
                            self.printLog('#EDIT','Split %s at %d -> %s' % (self.shortName(),divi,newshort))
                            self.list['Edit'].append(('Divide',seqacc,'Split %s at %d -> %s' % (self.shortName(),divi,newshort)))
                            self.list['Seq'] = self.list['Seq'][:ei] + [('%s (Region 1 to %d)' % (seqname,divi),newseq), ('%s %s (Region %d to %d)' % (newshort,self.seqDesc(),divi+1,len(sequence)),addseq)] + self.list['Seq'][ei+1:]
                        else: self.vPrint('Split aborted.')
                    else: self.verbose(0,0,'Cannot split %s at %d: produces empty sequence.' % (self.shortName(),divc))
                elif choice == 'R':
                    self.printLog('#EDIT','RevComp %s' % self.shortName())
                    self.list['Edit'].append(('RevComp',seqacc,'RevComp %s' % self.shortName()))
                    sequence = rje_sequence.reverseComplement(sequence,rna=self.rna())
                    if self.seqDesc() and string.split(self.seqDesc())[0] == 'RevComp':
                        if len(string.split(self.seqDesc())) > 1: newdesc = string.join(string.split(self.seqDesc())[1:])
                        else: newdesc = ''
                    else: newdesc = 'RevComp %s' % self.seqDesc()
                    self.list['Seq'][ei] = ('%s %s' % (self.shortName(),newdesc),sequence)
                elif choice == 'V':
                    self.printLog('#EDIT','Reversed %s' % self.shortName())
                    self.list['Edit'].append(('Reversed',seqacc,'Reversed %s' % self.shortName()))
                    sequence = rje.strReverse(sequence)
                    if self.seqDesc() and string.split(self.seqDesc())[0] == 'Reversed':
                        if len(string.split(self.seqDesc())) > 1: newdesc = string.join(string.split(self.seqDesc())[1:])
                        else: newdesc = ''
                    else: newdesc = 'Reversed %s' % self.seqDesc()
                    self.list['Seq'][ei] = ('%s %s' % (self.shortName(),newdesc),sequence)
                elif choice == 'T': print('Not yet implemented!')
                elif choice == 'DNA':
                    self.setStr({'SeqType':'dna'})
                    self.setBool({'DNA':True})
                elif choice == 'M':
                    seqx = len(sequence)
                    mutx = rje.getFloat('Number of mutations? (<1 for proportion)',default=0.01,confirm=True)
                    multi = rje.yesNo('Allow multiple hits?',default='Y')
                    if not multi and mutx > len(sequence): self.warnLog('Cannot have more mutations than positions without multiple hits!'); mutx = len(sequence)
                    if mutx < 1: mutx = max(0,int(mutx*len(sequence)+0.5))
                    mutpos = rje.randomList(range(len(sequence)),listlen=mutx,replace=multi)
                    for m in mutpos:
                        new = sequence[m]
                        while new == sequence[m]:
                            if self.dna(): new = random.choice('GATC')
                            else: new = random.choice('ACDEFGHKLMNPQRSTVWY')
                        sequence = sequence[:m] + new + sequence[m+1:]
                    if seqx != len(sequence): raise ValueError
                    self.printLog('#EDIT','%s mutations added to %s' % (mutx,self.shortName()))
                    self.list['Edit'].append(('Mutate',seqacc,'%s mutations added to %s' % (mutx,self.shortName())))
                ## ~ [1d] Sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif choice == '<':
                    if ei:
                        self.list['Seq'] = self.list['Seq'][:ei-1] + [self.current(),self.list['Seq'][ei-1]] + self.list['Seq'][ei+1:]
                        ei -= 1
                    else: self.verbose(0,0,'Already first sequence in list!')
                elif choice == '>':
                    if ei < self.seqNum() - 1:
                        self.list['Seq'] = self.list['Seq'][:ei] + [self.list['Seq'][ei+1],self.current()] + self.list['Seq'][ei+2:]
                        ei += 1
                    else: self.verbose(0,0,'Already last sequence in list!')
                elif choice == '<<':
                    if ei:
                        self.list['Seq'] = [self.current()] + self.list['Seq'][:ei] + self.list['Seq'][ei+1:]
                        ei = 0
                    else: self.verbose(0,0,'Already first sequence in list!')
                elif choice == '>>':
                    if ei < self.seqNum() - 1:
                        self.list['Seq'] = self.list['Seq'][:ei] + self.list['Seq'][ei+1:] + [self.current()]
                        ei = self.seqNum() - 1
                    else: self.verbose(0,0,'Already last sequence in list!')
                elif choice == '-':
                    if not rje.yesNo('Remove %s?' % self.shortName()): continue
                    logtxt = rje.choice('Note for log file (optional)?:',default='',confirm=True)
                    if logtxt: self.printLog('#EDIT','Removed %s - %s' % (seqname,logtxt))
                    else: self.printLog('#EDIT','Removed %s' % seqname)
                    self.list['Edit'].append(('Remove',seqacc,'Removed %s' % seqname))
                    self.list['Seq'].pop(ei)
                    ei = min(ei,self.seqNum()-1)
                elif choice == '--':
                    # Select sequences
                    cseqx = min(rje.getInt('Number of sequences to delete?',default=self.seqNum()-ei),self.seqNum()-ei)
                    cseq = self.list['Seq'][ei:ei+cseqx]
                    cnames = self.names(seqs=cseq)
                    if not rje.yesNo('Delete %s sequences: %s?' % (cseqx,string.join(cnames,', '))):
                        self.vPrint('Deletion aborted.')
                        continue
                    logtxt = rje.choice('Note for log file (optional)?:',default='',confirm=True)
                    for (seqname,sseq) in cseq:
                        if logtxt: self.printLog('#EDIT','Removed %s - %s' % (seqname,logtxt))
                        else: self.printLog('#EDIT','Removed %s' % seqname)
                        self.list['Edit'].append(('Remove',seqacc,'Removed %s' % seqname))
                    self.list['Seq'] = self.list['Seq'][:ei] + self.list['Seq'][ei+cseqx:]
                    ei -= 1
                elif choice == '++' and rje.yesNo('Duplicate %s?' % self.shortName()):
                    newacc = rje.choice('New accession (no spaces) for duplicate %s?' % self.shortName(),default='%sX2' % self.seqAcc())
                    newacc = string.join(string.split(newacc),'')
                    if gnspacc: newshort = '%s_%s__%s' % (self.seqGene(),self.seqSpec(),newacc)
                    else: newshort = newacc
                    self.printLog('#EDIT','%s +> %s' % (self.shortName(),newshort))
                    self.list['Edit'].append(('Duplicate',seqacc,'%s +> %s' % (self.shortName(),newshort)))
                    self.list['Seq'].insert(ei+1,('%s %s' % (newshort,self.seqDesc()),sequence))
                    ei += 1
                elif choice == '+':
                    filename = None
                    try:
                        while not filename:
                            filename = rje.choice('Sequence file to add? [? to list directory content]')
                            if filename == '?':
                                self.vPrint(os.popen('ls').read())
                                filename = None
                                continue
                            if os.path.exists(filename): break
                            else:
                                self.errorLog('File "%s" not found!' % filename,printerror=False)
                                if not rje.yesNo('Try again?'): raise KeyboardInterrupt
                            filename = None
                        addcmd = ['seqin=%s' % filename,'seqmode=list','autoload=F']
                        newseq = SeqList(self.log,addcmd)
                        newseq.loadSeq()
                        if not newseq.seqNum(): self.errorLog('No Sequences loaded!',printerror=False); continue
                        if not rje.yesNo('Add %s sequences (%s...)?' % (newseq.seqNum(),string.join(newseq.names()[:3],','))):
                            self.printLog('#EDIT','Addition of %s sequences aborted.' % newseq.seqNum())
                            continue
                        self.list['Seq'] = self.list['Seq'][:ei] + newseq.list['Seq'] + self.list['Seq'][ei:]
                        self.printLog('#EDIT','Addition of %s sequences.' % newseq.seqNum())
                    except KeyboardInterrupt: continue
                    except: raise
                    print('Not yet implemented!')
                elif choice == '&':
                    # Select sequences
                    cseqx = min(rje.getInt('Number of sequences to append?',default=self.seqNum()-ei-1),self.seqNum()-ei-1)
                    cseq = self.list['Seq'][ei+1:ei+1+cseqx]
                    cnames = self.names(seqs=cseq)
                    if not rje.yesNo('Concatenate %s with %s sequences: %s?' % (self.shortName(),cseqx,string.join(cnames,', '))):
                        self.vPrint('Concatenation aborted.')
                        continue
                    gapx = {True:'N',False:'X'}[self.nt()]
                    gapn = rje.getInt('Length of "%s" linker' % gapx,blank0=True,default='0',confirm=True)
                    for (sname,sseq) in cseq:
                        if gapn > 0: sequence += gapx * gapn
                        sequence += sseq
                        seqname += ' & %s' % sname
                    self.list['Seq'][ei] = (seqname,sequence)
                    self.printLog('#EDIT','Joined %s' % seqname)
                    self.list['Edit'].append(('Join',seqacc,'Joined %s' % seqname))
                    self.list['Seq'] = self.list['Seq'][:ei+1] + self.list['Seq'][ei+1+cseqx:]
                elif choice == 'C':
                    # Select sequences
                    cseqx = min(rje.getInt('Number of sequences to combine?',default=self.seqNum()-ei-1),self.seqNum()-ei-1)
                    cseq = self.list['Seq'][ei+1:ei+1+cseqx]
                    cnames = self.names(seqs=cseq)
                    if not rje.yesNo('Combine %s with %s sequences: %s?' % (self.shortName(),cseqx,string.join(cnames,', '))):
                        self.vPrint('Consensus generation aborted.')
                        continue
                    # Run self.makeConsensus(refseq,subseq)
                    subseq = [sequence]     # The main sequence should be in the consensus sequences too!
                    aligned = True
                    reflen = len(sequence)
                    maxlen = reflen
                    for (sname,sseq) in cseq:
                        aligned = aligned and len(sseq) == len(sequence)
                        maxlen = max(maxlen,len(sseq))
                        subseq.append(sseq)
                    try:
                        if not aligned and rje.yesNo('Subsequences of different %s lengths! Reference = %s; Max = %s. Generate %s %s consensus?' % (self.units(),reflen,maxlen,maxlen,self.units()),default='N'):
                            aligned = True
                            sequence = sequence + '-' * (maxlen-reflen)
                        if not aligned and reflen < maxlen and rje.yesNo('Subsequences of different %s lengths! Reference = %s; Max = %s. Trim to %s %s consensus?' % (self.units(),reflen,maxlen,reflen,self.units()),default='N'):
                            aligned = True
                        if not aligned: raise KeyboardInterrupt
                        newseq = self.makeConsensus(sequence,subseq,mindepth=rje.getInt('Min. no. sequences for non-gap position?',default=1))
                        newdesc = 'Consensus (%s + %s) %s' % (self.shortName(),string.join(cnames,','),self.seqDesc())
                        self.list['Seq'][ei] = ('%s %s' % (self.shortName(),newdesc),newseq)
                        self.printLog('#EDIT',newdesc)
                        self.list['Edit'].append(('Combine',seqacc,newdesc))
                        self.list['Seq'] = self.list['Seq'][:ei+1] + self.list['Seq'][ei+1+cseqx:]
                    except KeyboardInterrupt:
                        self.printLog('#EDIT','Consensus generation from unaligned sequences aborted!')
                    except: raise
                ## ~ [1e] Menu ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                elif choice == 'L':
                    logtxt = rje.choice('Note for log file:',default='',confirm=True)
                    if logtxt: self.printLog('#NOTE',logtxt)
                elif choice == 'P':
                    ei -= 1
                    if ei < 0: ei = self.seqNum() - 1
                elif choice == 'N':
                    ei += 1
                    if ei >= self.seqNum(): ei = 0
                elif choice == 'J':
                    ei = max(1,min(self.seqNum(),rje.getInt('Jump to sequence (1-%d)?' % self.seqNum(),default=1))) - 1
                elif choice == '?': print('Not yet implemented!')
                elif choice == 'O': # Output
                    self.saveSeq(seqfile=self.getStr('EditFile'))
                    self.setStr({'SeqIn':self.getStr('EditFile')})
                elif choice == 'Q' or not choice:    # Cancel/quit
                    if save and rje.yesNo('Save to %s?' % self.getStr('EditFile')):
                        self.saveSeq(seqfile=self.getStr('EditFile'))
                        self.setStr({'SeqIn':self.getStr('EditFile')})
                    break
            if not self.list['Seq']: self.printLog('#EDIT','No sequences for edit!'); return False
            return True
        except KeyboardInterrupt:
            return False
        except SystemExit: raise
        except: self.errorLog('Problem during %s edit.' % self); return False   # Edit failed
#########################################################################################################################
    def makeConsensus(self,refseq,subseq,mindepth=1):  ### Makes a consensus sequence from reference and subsequences
        '''
        Makes a consensus sequence from reference and subsequences.
        >> refseq:str = Full length Reference sequence
        >> subseq:list of strings making consensus (may be missing Cterm/3')
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if len(subseq) == 1: return refseq
            conseq = ''
            wild = {True:'N',False:'X'}[self.nt()]
            nocons = ['*',wild,'']
            if mindepth < 2: nocons.append('-')
            if mindepth > len(subseq):
                self.warnLog('Reduced consensus mindepth=%d to sequence number (mindepth=%d)' % (mindepth,len(subseq)))
                mindepth = len(subseq)
            ### ~ [1] ~ Compile ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for i in range(len(refseq)):
                r = refseq[i]
                ilist = []
                nongapx = 0
                for seq in subseq:
                    if seq[i:i+1] not in nocons: ilist.append(seq[i])
                    if seq[i] != '-': nongapx += 1
                for x in ilist:
                    if x == '-' and nongapx >= mindepth: continue
                    if ilist.count(x) > ilist.count(r): r = x
                    elif nongapx >= mindepth and r == '-': r = x
                conseq = conseq + r
            return conseq
        except: self.errorLog(rje_zen.Zen().wisdom()); return refseq
#########################################################################################################################
    def fileEdit(self): ### Makes changes from edit file.
        '''Makes changes from edit file.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.obj['DB']: self.obj['DB'] = rje_db.Database(self.log,self.cmd_list+['tuplekeys=T'])
            db = self.db()
            edb = db.addTable(self.getStr('Edit'),mainkeys=['Locus','Pos'],name='edit')
            edb.dataFormat({'Pos':'int'})
            #ekeys = edb.dataKeys()
            seqdict = self.makeSeqNameDic('max')
            if self.getStrLC('SeqOut'):
                self.setStr({'EditFile':self.getStr('SeqOut')})
            else:
                self.setStr({'EditFile':'%s.edit.fas' % rje.baseFile(self.getStr('SeqIn'))})

            ### ~ [2] ~ Make changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for locus in edb.index('Locus'):
                ekeys = edb.index('Locus')[locus]
                ekeys.reverse()     # Later positions edited first
                (seqname,seq) = self.getSeq(seqdict[locus])
                seqi = self.getSeq((seqname,seq),format='index')
                for ekey in ekeys:
                    entry = edb.data(ekey)
                    try:
                        slen = len(seq)
                        etype = entry['Edit'][:3].upper()
                        if etype == 'DEL':
                            newseq = ''
                            elen = len(entry['Details'])
                            oldseq = seq[entry['Pos']-1:entry['Pos']+elen-1]
                            if oldseq != entry['Details']: self.warnLog('Deletion sequence mismatch for %s %s (%s)' % (entry['Locus'],entry['Pos'],entry['Details']))
                        elif etype == 'INS':
                            newseq = seq[entry['Pos']-1] + entry['Details']
                            elen = 1    # Needs to be 1 because the original position is part of newseq and must be skipped.
                            oldseq = seq[entry['Pos']-1]
                        elif etype == 'SUB':
                            newseq = entry['Details']
                            elen = 1
                            oldseq = seq[entry['Pos']-1:entry['Pos']]
                        else: raise ValueError('Edit Type "%s" not recognised' % etype)
                        seq = seq[:entry['Pos']-1] + newseq + seq[entry['Pos']+elen-1:]
                        if etype == 'DEL': newseq = '-'
                        if 'Notes' in edb.fields():
                            self.printLog('#%s' % etype,'%s: %s -> %s. %s (%s -> %s %s)' % (entry['Pos'],oldseq,newseq,entry['Notes'],rje.iStr(slen),rje.iLen(seq),self.units()))
                        else:
                            self.printLog('#%s' % etype,'%s: %s -> %s. (%s -> %s %s)' % (entry['Pos'],oldseq,newseq,rje.iStr(slen),rje.iLen(seq),self.units()))
                        if not seqname.endswith(' [edited]'): seqname += ' [edited]'
                        self.list['Seq'][seqi] = (seqname,seq)
                    except:
                        self.errorLog('Problem with edit: %s' % entry,quitchoice=True)
                        self.debug(self.list['Seq'])
                seqdict[locus] = (seqname,seq)
            ### ~ [3] ~ Save changes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if rje.yesNo('Save to %s?' % self.getStr('EditFile')):
                self.saveSeq(seqfile=self.getStr('EditFile'))
                self.setStr({'SeqIn':self.getStr('EditFile')})
            else: self.warnLog('Edits not saved!')
            return True
        except: self.errorLog(rje_zen.Zen().wisdom()); return False
#########################################################################################################################
    def gapFix(self): ### Resizes gaps where appropriate.
        '''Resizes gaps where appropriate.'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            seqlist = []
            gapfix = self.dict['GapFix']
            if not gapfix: return True
            if not self.dna(): self.printLog('#GAPFIX','Gap fixing only available for dna=T'); return True
            # Convert to integers
            for gaplen in rje.sortKeys(gapfix):
                try: gapfix[int(gaplen)] = int(gapfix.pop(gaplen))
                except: raise ValueError('gapfix=X:Y dictionary must contain integers')
            for gaplen in rje.sortKeys(gapfix):
                self.printLog('#GAPFIX','Converting gaps: {0}bp -> {1}bp'.format(gaplen,gapfix[gaplen]))
            mingap = min(gapfix.keys())
            gapre = re.compile('N{%d,}' % mingap)
            ### ~ [1] ~ Convert gaps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            fixed = 0; seqx = 0.0; seqtot = self.seqNum()
            for seq in self.seqs():
                self.progLog('\r#GAPFIX','Fixing gap lengths: %.1f%% (%s fixed)' % (seqx/seqtot,rje.iStr(fixed))); seqx += 100.0
                (sname,sequence) = self.getSeq(seq,'tuple')
                sequence = sequence.upper()
                tofix = []
                for m in gapre.finditer(sequence):
                    if len(m.group()) in gapfix: tofix.append((m.start(),m.end(),len(m.group())))
                if tofix:
                    tofix.reverse()
                    for (gstart,gend,glen) in tofix:
                        fixed += 1
                        sequence = sequence[:gstart] + 'N' * gapfix[glen] + sequence[gend:]
                    if len(tofix) > 1:
                        sname = sname + ' (%d gaps resized)' % len(tofix)
                    else:
                        sname = sname + ' (%d gap resized)' % len(tofix)
                seqlist.append((sname,sequence))
                if self.mode() == 'db':
                    seq['Name'] = sname
                    seq['Sequence'] = sequence
            if self.mode() != 'db':
                self.setStr({'SeqMode':'list'})
                self.list['Seq'] = seqlist
            self.printLog('\r#GAPFIX','Fixed %s gap lengths (%d alt gap size)' % (rje.iStr(fixed),len(gapfix)))
            return True
        except: self.errorLog(rje_zen.Zen().wisdom()); return False
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
def phredScore(qchar,qscale=33): ### Returns the Phred quality score for a given character.
    '''Returns the Phred quality score for a given character.'''
    return ord(qchar) - qscale
#########################################################################################################################
def bpFromStr(seqlen):  ### Returns the number of basepairs as an integer from a str with possible [TGMk][b/bp] suffix
    '''
    Returns the number of basepairs as an integer from a str with possible [TGMk][b/bp] suffix
    :param seqlen: string version of seqlen
    :return: integer version of seqlen
    '''
    seqlen = ''.join(seqlen.lower().split())
    if seqlen.endswith('p'): seqlen = seqlen[:-1]
    if seqlen.endswith('b'): seqlen = seqlen[:-1]
    mult = 1
    if seqlen.endswith('t'): seqlen = seqlen[:-1]; mult = 1e12
    elif seqlen.endswith('g'): seqlen = seqlen[:-1]; mult = 1e9
    elif seqlen.endswith('m'): seqlen = seqlen[:-1]; mult = 1e6
    elif seqlen.endswith('k'): seqlen = seqlen[:-1]; mult = 1e3
    try:
        seqlen = int(float(seqlen * mult)+0.5)
    except:
        raise ValueError('Expected number followed by T/G/M/k')
    return seqlen
#########################################################################################################################
def dnaLen(seqlen,dp=2,sf=3):
    units = ['bp','kb','Mb','Gb','Tb']
    while seqlen > 1000 and len(units) > 1:
        units.pop(0); seqlen /= 1000.0
    units = units[0]
    if units == 'bp': return '%d bp' % seqlen
    if dp == 0: return '%s %s' % (rje.sf(seqlen,sf),units)
    if dp == 1: return '%.1f %s' % (rje.dp(seqlen,dp),units)
    if dp == 2: return '%.2f %s' % (rje.dp(seqlen,dp),units)
    if dp == 3: return '%.3f %s' % (rje.dp(seqlen,dp),units)
    if dp == 4: return '%.4f %s' % (rje.dp(seqlen,dp),units)
    return '%s %s' % (rje.dp(seqlen,dp),units)
#########################################################################################################################
def batchSummarise(callobj,seqfiles,save=True,overwrite=False):   ### Batch run seqlist summarise on seqfiles and output table of results
    '''
    Batch run seqlist summarise on seqfiles and output table of results. General options (dna=T etc.) should be given
    as part of the callobj.cmd_list.
    >> callobj:object controlling log output and providing database object. (Will be added if None.)
    >> seqfiles:list of sequence files to summarise
    >> save:bool [True] = Whether to save the summarise table to a file.
    >> overwrite:bool [False] = Whether to overwrite (True) or update (False) existing summary table.
    << summarise db Table
    '''
    try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if 'DB' in callobj.obj and callobj.obj['DB']: db = callobj.obj['DB']
        else:
            callobj.obj['DB'] = db = rje_db.Database(callobj.log,callobj.cmd_list)
        sdb = None
        if not overwrite:
            sdb = callobj.db('summarise')
            if not sdb:
                sdb = db.addTable(mainkeys=['File'],name='summarise',expect=False)
        if not sdb: sdb = db.addEmptyTable('summarise',['File'],['File'])
        ### ~ [2] Run Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        callobj.printLog('#BATCH','Batch summarising %s input files' % rje.iLen(seqfiles))
        for file in seqfiles:
            seqdata = SeqList(callobj.log,callobj.cmd_list+['seqin=%s' % file,'autoload=T','summarise=F','sortseq=None']).summarise()
            if seqdata:
                if 'GC' in seqdata:
                    seqdata.pop('GC')
                    seqdata['GCPC'] = '%.2f' % seqdata['GCPC']
                if 'GapLength' in seqdata: seqdata['GapPC'] = '%.2f' % (100.0*seqdata['GapLength']/seqdata['TotLength'])
                seqdata['MeanLength'] = '%.1f' % seqdata['MeanLength']
                for field in string.split('SeqNum, TotLength, MinLength, MaxLength, MeanLength, MedLength, N50Length, L50Count, CtgNum, N50Ctg, L50Ctg, NG50Length, LG50Count, GapLength, GapPC, GCPC',', '):
                    if field in seqdata and field not in sdb.fields(): sdb.addField(field)
                for field in seqdata.keys():
                    if field not in sdb.fields(): sdb.addField(field)
                sdb.addEntry(seqdata)
            else: callobj.errorLog('Summarise failed for %s' % file,printerror=False)
        ### ~ [3] Output Summarise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if save: sdb.saveToFile()
        return sdb
    except: callobj.errorLog('%s.batchSummarise error' % callobj.prog()); return False
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
    except: print('Unexpected error during program setup:', sys.exc_info()[0]); return
    
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
    except: print('Cataclysmic run error:', sys.exc_info()[0])
    os._exit(0)
#########################################################################################################################
### END OF SECTION IV                                                                                                   #
#########################################################################################################################
