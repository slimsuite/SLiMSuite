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
# Author contact: <redwards@cabbagesofdoom.co.uk> / 99 Wilton Road, Southampton SO15 5JH.
#
# To incorporate this module into your own programs, please see GNU Lesser General Public License disclaimer in rje.py

"""
Program:      SLiMFinder
Description:  Short Linear Motif Finder
Version:      5.2.3
Last Edit:    07/12/15
Citation:     Edwards RJ, Davey NE & Shields DC (2007), PLoS ONE 2(10): e967. [PMID: 17912346]
ConsMask Citation: Davey NE, Shields DC & Edwards RJ (2009), Bioinformatics 25(4): 443-50. [PMID: 19136552]
SigV/SigPrime Citation: Davey NE, Edwards RJ & Shields DC (2010), BMC Bioinformatics 11: 14. [PMID: 20055997]
SLiMScape/REST Citation: Olorin E, O'Brien KT, Palopoli N, Perez-Bercoff A & Shields DC, Edwards RJ (2015), F1000Research 4:477.
SLiMMaker Citation: Palopoli N, Lythgow KT & Edwards RJ (2015), Bioinformatics 31(14): 2284-2293. [PMID: 25792551]
Webserver:    http://bioware.ucd.ie/
Copyright (C) 2007  Richard J. Edwards - See source code for GNU License Notice

Function:
    Short linear motifs (SLiMs) in proteins are functional microdomains of fundamental importance in many biological
    systems. SLiMs typically consist of a 3 to 10 amino acid stretch of the primary protein sequence, of which as few
    as two sites may be important for activity, making identification of novel SLiMs extremely difficult. In particular,
    it can be very difficult to distinguish a randomly recurring "motif" from a truly over-represented one. Incorporating
    ambiguous amino acid positions and/or variable-length wildcard spacers between defined residues further complicates
    the matter.

    SLiMFinder is an integrated SLiM discovery program building on the principles of the SLiMDisc software for accounting
    for evolutionary relationships [Davey NE, Shields DC & Edwards RJ (2006): Nucleic Acids Res. 34(12):3546-54].
    SLiMFinder is comprised of two algorithms:

    1. `SLiMBuild` identifies convergently evolved, short motifs in a dataset. Motifs with fixed amino acid positions are
    identified and then combined to incorporate amino acid ambiguity and variable-length wildcard spacers. Unlike
    programs such as TEIRESIAS, which return all shared patterns, SLiMBuild accelerates the process and reduces returned
    motifs by explicitly screening out motifs that do not occur in enough unrelated proteins. For this, SLiMBuild uses
    the "Unrelated Proteins" (UP) algorithm of SLiMDisc in which BLAST is used to identify pairwise relationships.
    Proteins are then clustered according to these relationships into "Unrelated Protein Clusters" (UPC), which are
    defined such that no protein in a UPC has a BLAST-detectable relationship with a protein in another UPC.  If desired,
    `SLiMBuild` can be used as a replacement for TEIRESIAS in other software (teiresias=T slimchance=F).

    2. `SLiMChance` estimates the probability of these motifs arising by chance, correcting for the size and composition
    of the dataset, and assigns a significance value to each motif. Motif occurrence probabilities are calculated
    independently for each UPC, adjusted for the size of a UPC using the Minimum Spanning Tree algorithm from SLiMDisc.
    These individual occurrence probabilities are then converted into the total probability of the seeing the observed
    motifs the observed number of (unrelated) times. These probabilities assume that the motif is known before the
    search. In reality, only over-represented motifs from the dataset are looked at, so these probabilities are adjusted
    for the size of motif-space searched to give a significance value. The returned corrected probability is an estimate
    of the probability of seeing ANY motif with that significance (or greater) from the dataset (i.e. an estimate of the
    probability of seeing that motif, *or another one like it*). These values are calculated separately for each length
    of motif.

    SLiMFinder version 4.0 introduced a more precise (but more computationally intensive) statistical model, which can
    be switched on using sigprime=T. Likewise, the more precise (but more computationally intensive) correction to the
    mean UPC probability heuristic can be switched on using sigv=T. (Note that the other `SLiMChance` options may not
    work with either of these options.) The allsig=T option will output all four scores. In this case, SigPrimeV will be
    used for ranking etc. unless probscore=X is used.

    ### Clouds and Statistics:
    Where significant motifs are returned, SLiMFinder will group them into Motif "Clouds", which consist of physically
    overlapping motifs (2+ non-wildcard positions are the same in the same sequence). This provides an easy indication
    of which motifs may actually be variants of a larger SLiM and should therefore be considered together. From version
    V4.7, `*.cloud.txt` output includes a `SLiMMaker` summary Regex for the whole cloud. NOTE: This may not necessarily
    match all occurrences in the cloud.

    Additional Motif Occurrence Statistics, such as motif conservation, are handled by the `rje_slimlist` module and
    `rje_slimcalc` modules. Please see the documentation for these module for a full list of commandline options. These
    options have not been fully tested in SLiMFinder, so please report issues and/or request desired functions. Note that
    occfilter=LIST *does* affect the motifs returned by SLiMBuild and thus the TEIRESIAS output (as does min. IC and min.
    Support) but the overall Motif slimfilter=LIST *only* affects SLiMFinder output following SLiMChance calculations.

    ### Secondary Functions:
    The "MotifSeq" option will output fasta files for a list of X:Y, where X is a motif pattern and Y is the output file.

    The "Randomise" function will take a set of input datasets (as in Batch Mode) and regenerate a set of new datasets
    by shuffling the UPC among datasets. Note that, at this stage, this is quite crude and may result in the final
    datasets having fewer UPC due to common sequences and/or relationships between UPC clusters in different datasets.

    Where pre-known motifs are also of interest, these can be given with the slimcheck=MOTIFS option and will be added to
    the output. In general, it is better to use `SLiMProb` to look for enrichment (or depletion) of pre-defined motifs.

Input/Output:
    ### SLiMFinder Input:
    The main input for SLiMFinder is the `seqin=SEQFILE` file of protein sequences, which can be Uniprot plain text
    (`DATFILE`) or fasta (`FASFILE`) format. A batch of files (incorporating wildcards) can be given using
    batch=FILELIST. Alternative primary input is uniprotid=LIST. This requires an active internet connection to retrieve
    the corresponding Uniprot entries.

    ### SLiMFinder Output:
    Please see Manual for details.

Commandline: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Basic Input/Output Options: ###
    seqin=SEQFILE   : Sequence file to search. Over-rules batch=FILE and uniprotid=LIST [None]
    batch=FILELIST  : List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    uniprotid=LIST  : Extract IDs/AccNums in list from Uniprot into BASEFILE.dat and use as seqin=FILE. []
    maxseq=X        : Maximum number of sequences to process [500]
    maxupc=X        : Maximum UPC size of dataset to process [0]
    sizesort=X      : Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    walltime=X      : Time in hours before program will abort search and exit [1.0]
    resfile=FILE    : Main SLiMFinder results table [slimfinder.csv]
    resdir=PATH     : Redirect individual output files to specified directory (and look for intermediates) [SLiMFinder/]
    buildpath=PATH  : Alternative path to look for existing intermediate files [SLiMFinder/]
    force=T/F       : Force re-running of BLAST, UPC generation and SLiMBuild [False]
    pickup=T/F      : Pick-up from aborted batch run by identifying datasets in resfile [False]
    pickid=T/F      : Whether to use RunID to identify run datasets when using pickup [True]
    pickall=T/F     : Whether to skip aborted runs (True) or only those datasets that ran to completion (False) [True]
    dna=T/F         : Whether the sequences files are DNA rather than protein [False]
    alphabet=LIST   : List of characters to include in search (e.g. AAs or NTs) [default AA or NT codes]
    megaslim=FILE   : Make/use precomputed results for a proteome (FILE) in fasta format [None]
    megablam=T/F    : Whether to create and use all-by-all GABLAM results for (gablamdis) UPC generation [False]
    ptmlist=LIST    : List of PTM letters to add to alphabet for analysis and restrict PTM data []
    ptmdata=DSVFILE : File containing PTM data, including AccNum, ModType, ModPos, ModAA, ModCode

SLiMBuild: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### SLiMBuild Options I (Evolutionary Filtering): ###
    efilter=T/F     : Whether to use evolutionary filter [True]
    blastf=T/F      : Use BLAST Complexity filter when determining relationships [True]
    blaste=X        : BLAST e-value threshold for determining relationships [1e=4]
    altdis=DSVFILE  : Alternative all by all distance matrix for relationships [None]
    gablamdis=FILE  : Alternative GABLAM results file [None] (!!!Experimental feature!!!)
    homcut=X        : Max number of homologues to allow (to reduce large multi-domain families) [0]
    newupc=PATH     : Look for alternative UPC file and calculate Significance using new clusters [None]

    ### SLiMBuild Options II (Input Masking): ###
    masking=T/F     : Master control switch to turn off all masking if False [True]
    dismask=T/F     : Whether to mask ordered regions (see rje_disorder for options) [False]
    consmask=T/F    : Whether to use relative conservation masking [False]
    ftmask=LIST     : UniProt features to mask out (True=EM,DOMAIN,TRANSMEM) []
    imask=LIST      : UniProt features to inversely ("inclusively") mask. (Seqs MUST have 1+ features) []
    compmask=X,Y    : Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    casemask=X      : Mask Upper or Lower case [None]
    motifmask=X     : List (or file) of motifs to mask from input sequences []
    metmask=T/F     : Masks the N-terminal M (can be useful if termini=T) [True]
    posmask=LIST    : Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    aamask=LIST     : Masks list of AAs from all sequences (reduces alphabet) []
    qregion=X,Y     : Mask all but the region of the query from (and including) residue X to residue Y [0,-1]

    ### SLiMBuild Options III (Basic Motif Construction): ###
    termini=T/F     : Whether to add termini characters (^ & $) to search sequences [True]
    minwild=X       : Minimum number of consecutive wildcard positions to allow [0]
    maxwild=X       : Maximum number of consecutive wildcard positions to allow [2]
    slimlen=X       : Maximum length of SLiMs to return (no. non-wildcard positions) [5]
    minocc=X        : Minimum number of unrelated occurrences for returned SLiMs. (Proportion of UP if < 1) [0.05]
    absmin=X        : Used if minocc<1 to define absolute min. UP occ [3]
    alphahelix=T/F  : Special i, i+3/4, i+7 motif discovery [False]
    fixlen=T/F      : If true, will use maxwild and slimlen to define a fixed total motif length [False]
    palindrome=T/F  : Special DNA mode that will search for palindromic sequences only [False]

    ### SLiMBuild Options IV (Ambiguity): ###
    ambiguity=T/F   : (preamb=T/F) Whether to search for ambiguous motifs during motif discovery [True]
    ambocc=X        : Min. UP occurrence for subvariants of ambiguous motifs (minocc if 0 or > minocc) [0.05]
    absminamb=X     : Used if ambocc<1 to define absolute min. UP occ [2]
    equiv=LIST      : List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMVF,FYW,FYH,KRH,DE,ST]
    wildvar=T/F     : Whether to allow variable length wildcards [True]
    combamb=T/F     : Whether to search for combined amino acid degeneracy and variable wildcards [False]

    ### SLiMBuild Options V (Advanced Motif Filtering): ###
    altupc=PATH     : Look for alternative UPC file and filter based on minocc [None]
    musthave=LIST   : Returned motifs must contain one or more of the AAs in LIST (reduces search space) []
    query=LIST      : Return only SLiMs that occur in 1+ Query sequences (Name/AccNum) []
    focus=FILE      : FILE containing focal groups for SLiM return (see Manual for details) [None]
    focusocc=X      : Motif must appear in X+ focus groups (0 = all) [0]
    * See also rje_slimcalc options for occurrence-based calculations and filtering *
    
SLiMChance: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    cloudfix=T/F    : Restrict output to clouds with 1+ fixed motif (recommended) [False]
    slimchance=T/F  : Execute main SLiMFinder probability method and outputs [True]
    sigprime=T/F    : Calculate more precise (but more computationally intensive) statistical model [False]
    sigv=T/F        : Use the more precise (but more computationally intensive) fix to mean UPC probability [False]
    dimfreq=T/F     : Whether to use dimer masking pattern to adjust number of possible sites for motif [True]
    probcut=X       : Probability cut-off for returned motifs (sigcut=X also recognised) [0.1]
    maskfreq=T/F    : Whether to use masked AA Frequencies (True), or (False) mask after frequency calculations [True]
    aafreq=AAFILE   : Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    aadimerfreq=FILE: Use empirical dimer frequencies from FILE (fasta or *.aadimer.tdt) (!!!Experimental!!!) [None]
    negatives=SEQFILE : Multiply raw probabilities by under-representation in FILE (!!!Experimental!!!) [None]
    smearfreq=T/F   : Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    seqocc=T/F      : Whether to upweight for multiple occurrences in same sequence (heuristic) [False]
    probscore=X     : Score to be used for probability cut-off and ranking (Prob/Sig/S/R) [Sig]

Advanced: ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    ### Advanced Output Options I (Output data): ###
    clouds=X        : Identifies motif "clouds" which overlap at 2+ positions in X+ sequences (0=minocc / -1=off) [2]
    runid=X         : Run ID for resfile (allows multiple runs on same data) [DATE]
    logmask=T/F     : Whether to log the masking of individual sequences [True]
    slimcheck=MOTIFS : Motif file/list to add to resfile output []

    ### Advanced Output Options II (Output formats): ###
    teiresias=T/F   : Replace TEIRESIAS, making *.out and *.mask.fasta files [False]
    slimdisc=T/F    : Emulate SLiMDisc output format (*.rank & *.dat.rank + TEIRESIAS *.out & *.fasta) [False]
    extras=X        : Whether to generate additional output files (alignments etc.) [1]
                        --1 = No output beyond main results file
                        - 0 = Generate occurrence file
                        - 1 = Generate occurrence file, alignments and cloud file
                        - 2 = Generate all additional SLiMFinder outputs
                        - 3 = Generate SLiMDisc emulation too (equiv extras=2 slimdisc=T)
    targz=T/F       : Whether to tar and zip dataset result files (UNIX only) [False]
    savespace=0     : Delete "unneccessary" files following run (best used with targz): [0]
                        - 0 = Delete no files
                        - 1 = Delete all bar *.upc and *.pickle
                        - 2 = Delete all bar *.upc (pickle added to tar)
                        - 3 = Delete all dataset-specific files including *.upc and *.pickle (not *.tar.gz)

    ### Advanced Output Options III (Additional Motif Filtering): ###
    topranks=X      : Will only output top X motifs meeting probcut [1000]
    oldscores=T/F   : Whether to also output old SLiMDisc score (S) and SLiMPickings score (R) [False]
    allsig=T/F      : Whether to also output all SLiMChance combinations (Sig/SigV/SigPrime/SigPrimeV) [False]
    minic=X         : Minimum information content for returned motifs [2.1]
    * See also rje_slimcalc options for occurrence-based calculations and filtering *

    ### Additional Functions I (MotifSeq): ###
    motifseq=LIST   : Outputs fasta files for a list of X:Y, where X is the pattern and Y is the output file []
    slimbuild=T/F   : Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]

    ### Additional Functions II (Randomised datasets): ###
    randomise=T/F   : Randomise UPC within batch files and output new datasets [False]
    randir=PATH     : Output path for creation of randomised datasets [Random/]
    randbase=X      : Base for random dataset name [rand]
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
"""
#########################################################################################################################
### SECTION I: GENERAL SETUP & PROGRAM DETAILS                                                                          #
#########################################################################################################################
import os, pickle, string, sys, time
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../libraries/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),'../tools/'))
### User modules - remember to add *.__doc__ to cmdHelp() below ###
import rje, rje_seq, rje_sequence, rje_scoring, rje_xgmml
import rje_slim, rje_slimcalc, rje_slimcore, rje_slimlist, slimmaker
import rje_motif_V3 as rje_motif            # Used for expect method only 
#import rje_dismatrix_V2 as rje_dismatrix
import comparimotif_V3 as comparimotif
import ned_rankbydistribution
#########################################################################################################################
def history():  ### Program History - only a method for PythonWin collapsing! ###
    '''
    # 0.0 - Initial Compilation.
    # 1.0 - Preliminary working version with Poisson probabilities
    # 1.1 - Binomial probabilities, bonferroni corrections and complexity masking
    # 1.2 - Added musthave=LIST option and denferroni correction.
    # 1.3 - Added resfile=FILE output
    # 1.4 - Added option for termini
    # 1.5 - Reworked slim mechanics to be ai-x-aj strings for future ambiguity (split on '-' to make list)
    # 1.6 - Added basic ambiguity and flexible wildcards plus MST weighting for UP clusters
    # 1.7 - Added counting of generic dimer frequencies for improved Bonferroni and probability calculation (No blockmask.)
    #     - Added topranks=X and query=X
    # 1.8 - Added *.upc rather than *.self.blast. Added basic randomiser function.
    # 1.9 - Added MotifList object to handle extra calculations and occurrence filtering.
    # 2.0 - Tidied up and standardised output. Implemented extra filtering and scoring options.
    # 2.1 - Changed defaults. Removed poisson as option and other obseleted functions.
    # 2.2 - Tidied and reorganised code using SLiMBuild/SLiMChance subdivision of labour. Removed rerun=T/F (just Force.)
    # 2.3 - Added AAFreq "smear" and "better" p1+ calculation. Added extra cloud summary output.
    # 2.4 - Minor bug fixes and tidying. Removed power output. (Rubbish anyway!) Can read UPC from distance matrix.
    # 3.0 - Dumped useless stats and calculations. Simplified output. Improved ambiguity & clouds.
    # 3.1 - Added minwild and alphahelix options. (Partial aadimerfreq & negatives)
    # 3.2 - Tidied up with SLiMCore, replaced old Motif objects with SLiM objects and SLiMCalc.
    # 3.3 - Added XGMML output. Added webserver option with additional output.
    # 3.4 - Added consmask relative conservation masking.
    # 3.5 - Standardised masking options. Add motifmask and motifcull.
    # 3.6 - Added aamasking and alphabet.
    # 3.7 - Added option to switch off dimfreq and better handling of given aafreq
    # 3.8 - Added SLiMDisc & SLiMPickings scores and options to rank on them.
    # 3.9 - Added clouding consensus information. [Aborted due to technical challenges.]
    # 3.10- Added differentiation of methods for pickling and tarring.
    # 4.0 - Added SigPrime and SigV calculation from Norman. Added graded extras output.
    # 4.1 - Added SizeSort, AltUPC and NewUPC options. Added #END output for webserver.
    # 4.2 - Added fixlen option and improved Alphahelix option
    # 4.3 - Updated the output for Max/Min filtering and the pickup options. Removed TempMaxSetting.
    # 4.4 - Modified to work with GOPHER V3.0.
    # 4.5 - Minor modifications to fix sigV and sigPrime bugs. Modified extras setting. Added palindrome setting for DNA motifs.
    # 4.6 - Minor modification to seqocc=T function. !Experimental! Added main occurrence output and modified savespace.
    # 4.7 - Added SLiMMaker generation to motif clouds. Added Q and Occ to Chance column.
    # 4.8 - Modified cloud generation to avoid issues with flexible-length wildcards.
    # 4.9 - Preparation for SLiMFinder V5.0 & SLiMCore V2.0 using newer RJE_Object.
    # 5.0 - Converted to use rje_obj.RJE_Object as base. Version 4.9 moved to legacy/.
    # 5.1 - Modified SLiMChance slightly to catch missing aafreq.
    # 5.1.1 - Modified alphabet handling and fixed musthave bug.
    # 5.2.0 - Added PTMList and PTMData modes (dev only).
    # 5.2.1 - Fixed ambocc<1 and minocc<1 issue. (Using integers rather than floats.) Fixed OccRes Sig output format.
    # 5.2.2 - Added warnings for ambocc and minocc that exceed the absolute minima. Updated docstring.
    # 5.2.3 - Switched feature masking OFF by default to give consistent Uniprot versus FASTA behaviour. Fixed FTMask=T/F bug.
    '''
#########################################################################################################################
def todo():     ### Major Functionality to Add - only a method for PythonWin collapsing! ###
    '''
    # [ ] : Interactive user Menu. (SLiMSuite Wrapper Module?)
    # [ ] : Add iterative search mode that re-searches with remaining sequences
    # [Y] : Improve the use of the pickling - check any options that might change output etc.
    # [?] : Add generation of data summary from Log file (from rje_misc) and tidy resfile=FILE output
    # [Y] : Incorporate MotifList object for ease of use of methods (e.g. IC) and future conservation etc.?
    # [Y] : Check and define when SeqOcc is used - used for calculating "Support" column
    # [ ] : Consider adding SlimCheck SLiMs to the alignment outputs - will need to find occurrences with PRESTO.
    # [ ] : Consider adding an MST file for the Focus Group that can be read like the UPC
    # [x] : Add gablam=FILE option to read gablam results from a file (and generate this file first if missing).
    # [Y] : Replace the rje_motif Expectation calculation for p1+ with an internal calculation.
    # [Y] : Add the number of sequences (absolute support) to output
    # [?] : Add positional variation from Nterm and Cterm to output? (Can calculate from Occ File)
    # [Y] : Add an option to "smear" aa frequencies over all UPC & output masked sequences in masked.fas (extras=T).
    # [Y] : Add gzipping of pickles in UNIX.
    # [Y] : Remove Normferroni and Occferroni measures. (Crap!)
    # [Y] : Make improved use of Query and Focus groups - binomial-based ajustment rather than reduced dataset size.
    # [ ] : Add line to output if walltime reached.
    # [ ] : Code-up negatives option - use MotifSeq as base?
    # [~] : Move stats and filtering on stats out into slimcalc and correctly use SLiMList.
    # [ ] : Add occscale=LIST : Rescale probabilities according to occstats []
    # [ ] : Add output compatible with Norman's motif drawing thingamie.
    # [ ] : Check what motifcull is and whether it is implemented.
    # [?] : Fix the bug that outputs some occurrences twice. (Variable wildcards only.)
    # [Y] : Sort out output and better pickup options.
    # [Y] : Add palindrome setting for DNA motifs.
    # [ ] : Add revcomp mode for DNA motifs
    # [ ] : Add recode=LIST option.
    # [Y] : Add SLiMMaker to cloud output.
    # [Y] : Add cloudfix=T for restricting output to ambiguous motifs ONLY if the cloud has a significant fixed variant.
    # [Y] : Modify basefile=X to NOT affect the individual results outputs.
    # [ ] : Consolidate and sort out use of basefile and dataset naming across SLiMSuite.
    '''
#########################################################################################################################
def makeInfo():     ### Makes Info object
    '''Makes rje.Info object for program.'''
    (program, version, last_edit, copyyear) = ('SLiMFinder', '5.2.3', 'December 2015', '2007')
    description = 'Short Linear Motif Finder'
    author = 'Richard J. Edwards, Norman E. Davey & Denis C. Shields'
    comments = ['Cite: Edwards, Davey & Shields (2007), PLoS ONE 2(10): e967. [PMID: 17912346]',
                'ConsMask: Davey NE, Shields DC & Edwards RJ (2009), Bioinformatics 25(4): 443-50. [PMID: 19136552]',
                'SigV/SigPrime: Davey NE, Edwards RJ & Shields DC (2010), BMC Bioinformatics 11: 14. [PMID: 20055997]',
                'SLiMMaker Citation: Palopoli N, Lythgow KT & Edwards RJ (2015), Bioinformatics 31(14): 2284-2293. [PMID: 25792551]',
                'Please report bugs to R.Edwards@UNSW.edu.au']
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
            if rje.yesNo('Show SLiMCalc commandline options?',default='N'): out.verbose(-1,4,text=rje_slimcalc.__doc__)
            if rje.yesNo('Show RJE_SEQ commandline options?',default='N'): out.verbose(-1,4,text=rje_seq.__doc__)
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
default_equiv = 'AGS,ILMVF,FYW,FYH,KRH,DE,ST'
probscores = {'none':'Sig','uncorrected':'Prob','raw':'Prob','bonferroni':'Sig','sig':'Sig','':'Sig','s':'S','r':'R',
              'slimdisc':'S','slimpicks':'R','slimpickings':'R','sigprime':'SigPrime','sigv':'SigV',
              'sigprimev':'SigPrimeV','sigvprime':'SigPrimeV','1':'Sig','2':'SigV','3':'SigPrime','4':'SigPrimeV'}
basic_headers = ['Rank','Sig','Pattern','IC','Occ','Support','UP','ExpUP','Prob','S','R','Cloud','CloudSeq','CloudUP']
#########################################################################################################################
### END OF SECTION I                                                                                                    #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION II: SLiMFinder Class                                                                                        #
#########################################################################################################################
class SLiMFinder(rje_slimcore.SLiMCore):     
    '''
    SLiMFinder Class. Author: Rich Edwards (2007).

    Str:str
    - AAFreq = Use FILE to replace individual sequence AAFreqs (FILE can be sequences or aafreq) [None]
    - AddQuery = Adds query sequence(s) to batch jobs from FILE [None]  (QSLiMFinder only)
    - AltDis = Alternative all by all distance matrix for relationships [None]
    - AltUPC = Look for alternative UPC file and filter based on minocc [None]
    - Build = String giving summary of key SLiMBuild options
    - BuildPath = Alternative path to look for existing intermediate files [SLiMFinder/]
    - CaseMask = Mask Upper or Lower case [None]
    - Chance = String giving summary of key SLiMChance options
    - CompMask = Mask low complexity regions (same AA in X+ of Y consecutive aas) [5,8]
    - Focus = FILE containing focal groups for SLiM return (see Manual for details) [None]
    - GablamDis = Alternative GABLAM results file [None]
    - Input = Original name (and path) of input file
    - AADimerFreq = Use empirical dimer frequencies from FILE (fasta or *.dimer.tdt) [None]
    - Negatives = Multiply raw probabilities by under-representation in FILE [None]
    - NewUPC = Look for alternative UPC file and calculate Significance using new clusters [None]
    - ProbScore = Score to be used for probability cut-off and ranking (Prob/Sig) [Sig]
    - RanDir = Output path for creation of randomised datasets [./]
    - Randbase = Base for random dataset name [rand_]
    - ResDir = Redirect individual output files to specified directory [SLiMFinder/]
    - ResFile = If FILE is given, will also produce a table of results in resfile [slimfinder.csv]
    - RunID = Run ID for resfile (allows multiple runs on same data) [DATE]
    - SlimCheck = Motif file/list to add to resfile output []
    
    Bool:boolean
    - AllSig = Whether to also output all SLiMChance combinations (Sig/SigV/SigPrime/SigPrimeV) [False]
    - AlphaHelix = Special i, i+3/4, i+7 motif discovery [False]
    - CloudFix = Restrict output to clouds with 1+ fixed motif (recommended) [False]
    - CombAmb = Whether to search for combined amino acid degeneracy and variable wildcards [True]
    - DimFreq = Whether to use dimer masking pattern to adjust number of possible sites for motif [True]
    - DisMask = Whether to mask ordered regions (see rje_disorder for options) [False]
    - DNA = Whether the sequences files are DNA rather than protein [False]
    - Force = whether to force recreation of key files [False]
    - EFilter = Whether to use evolutionary filter [True]
    - Extras = Whether to generate additional output files (alignments etc.) [False]
    - FixLen = If true, will use maxwild and slimlen to define a fixed total motif length  [False]
    - LogMask = Whether to log the masking of individual sequences [True]
    - OldScores = Whether to also output old SLiMDisc score (S) and SLiMPickings score (R) [False]
    - Palindrome = Special DNA mode that will search for palindromic sequences only [False]
    - Masked = Whether dataset has been masked [False]
    - Masking = Master control switch to turn off all masking if False [True]
    - MaskM = Masks the N-terminal M (can be useful if termini=T) [False]
    - PickAll = Whether to skip aborted runs (True) or only those datasets that ran to completion (False) [True]
    - PickID = Whether to use RunID to identify run datasets when using pickup [True]
    - PreAmb = Whether to search for ambiguous motifs during motif discovery [False]
    - QExact = Calculate exact Query motif space (True) or over-estimate from dimers (False) (quicker)[True]
    - OccStatsCalculated = Whether OccStats have been calculated for all occurrence [False]
    - MaskFreq = Whether to mask input before any analysis, or after frequency calculations [True]
    - SigV = Use the more precise (but more computationally intensive) fix to mean UPC probability [False]
    - Randomise = Randomise UPC within batch files [False]
    - SeqOcc = Whether to upweight for multiple occurrences in same sequence [True]
    - SigPrime = Calculate more precise (but more computationally intensive) statistical model [False]
    - SlimBuild = Whether to build motifs with SLiMBuild. (For combination with motifseq only.) [True]
    - SlimDisc = Output in SLiMDisc format (*.rank & *.dat.rank) [False]
    - SlimChance = Output in SLiMFinder Format (*.rank & *.occ.txt) [True]
    - SmearFreq = Whether to "smear" AA frequencies across UPC rather than keep separate AAFreqs [False]
    - TarGZ = Whether to tar and zip dataset result files (UNIX only) [False]
    - Teiresias = Replace TEIRESIAS only, making *.out and *.mask.fas files [False]
    - Termini = Whether to add termini characters (^ & $) to search sequences [True]
    - Test = Special Test parameter for experimentation with code [False]
    - WildVar = Whether to allow variable length wildcards [False]
    - Pickup = Pick-up from aborted batch run by identifying last dataset output in resfile [False]
    - Webserver = Generate additional webserver-specific output [False]

    Int:numeric
    - AbsMin = Used if minocc<1 to define absolute min. UP occ [3]
    - AbsMinAmb = Used if ambocc<1 to define absolute min. UP occ [2]    
    - Clouds = Identifies motif "clouds" which overlap at 2+ positions in X+ sequences (0=minocc) [2]
    - Extras = Whether to generate additional output files (alignments etc.) [1]
    - FocusOcc = Motif must appear in X+ focus groups (0 = all) [0]
    - Minwild = Minimum number of consecutive wildcard positions to allow [0]
    - MaxWild = Maximum number of consecutive wildcard positions to allow [3]
    - SlimLen = Maximum length of SLiMs to return (no. non-wildcard positions) [5]
    - MaxSeq = Maximum number of sequences to process [500]
    - MaxUPC = Maximum UPC size of dataset to process [0]
    - SaveSpace = Delete "unneccessary" files following run (see Manual for details) [0]
    - SizeSort = Sorts batch files by size prior to running (+1 small->big; -1 big->small; 0 none) [0]
    - TopRanks = Will only output top X motifs meeting probcut [0]

    Num:numeric
    - AmbOcc = Min. UP occurrence for subvariants of ambiguous motifs (minocc if 0 or > minocc) [0]
    - MinOcc = Minimum number of unrelated occurrences for returned SLiMs [2]
    - MinIC = Minimum information content for returned motifs [1.1]
    - MST = MST corrected size for whole dataset
    - ProbCut = Probability cut-off for returned motifs [0.01]
    - StartTime = Starting time in seconds (for output when using shared log file)
    - WallTime = Time in hours before program will terminate [1.0]

    List:list
    - AAMask = Masks list of AAs from all sequences (reduces alphabet) []
    - Alphabet = List of characters to include in search (e.g. AAs or NTs)
    - Batch = List of files to search, wildcards allowed. (Over-ruled by seqin=FILE.) [*.dat,*.fas]
    - Equiv = List (or file) of TEIRESIAS-style ambiguities to use [AGS,ILMV,FYW,KRH,DE,ST]
    - FTMask = UniProt features to mask out []
    - Headers = Headers for main SLiMFinder output table
    - HomCut = Max number of homologues to allow (to reduce large multi-domain families) [0]
    - IMask = UniProt features to inversely ("inclusively") mask [IM]
    - MustHave = Returned motifs must contain one or more of the AAs in LIST (reduces search space) []
    - OccScale = Rescale probabilities according to occstats (see manual for details) []
    - Query = Return only SLiMs that occur in 1+ Query sequences (Name/AccNum) []
    - SigSlim = List of significant SLiMs - matches keys to self.dict['Slim(Freq)'] - *in rank order*
    - SlimCheckExtra = List of extra SLiMs from SLiMCheck - added to extra outputs
    - UP = List of UP cluster tuples
    - Warning = List of text (log) warnings to reproduce at end of run
     
    Dict:dictionary
    - AADimerFreq = Empirical dimer frequencies {i:{x:{j:freq}}}
    - AAFreq = AA frequency dictionary for each seq / UPC
    - Clouds = Dictionary of Motif Clouds made by the makeClouds() method {ID':{stats}}
    - DimFreq = Frequency of dimers of each X length per upc {upc:[freqs]}
    - Dimers = main nested dictionary for SLiMFinder {Ai:{X:{Aj:{'UP':[UPC],'Occ':[(Seq,Pos)]}}}}
    - ElementIC = dictionary of {Motif Element:Information Content}
    - Focus = Dictionary of {focal group:list of seqs}
    - Extremf. = Dictionary of {length:extremferroni correction}
    - MaskPos = Masks list of position-specific aas, where list = pos1:aas,pos2:aas  [2:A]
    - MotifSeq = Dictionary of {pattern:output file for sequences}
    - MST = MST corrected size for UPC {UPC:MST}
    - Slim = main dictionary containing SLiMs with enough support {Slim:{'UPC':[UPC],'Occ':[Seq,Pos]}}
    - SeqOcc = dictionary of {Slim:{Seq:Count}} 

    Obj:RJE_Objects
    - SlimCheck = MotifList object handling motifs to check
    - SeqList = main SeqList object containing dataset to be searched
    - SlimList = SLiMList object handling motif stats and filtering options
    '''
#########################################################################################################################
    ### <1> ### Class Initiation etc.: sets attributes                                                                  #
#########################################################################################################################
    def _setAttributes(self):   ### Sets Attributes of Object
        '''Sets Attributes of Object.'''
        ### Basics ###
        self.strlist = ['AAFreq','Input','CaseMask','CompMask','ResFile','RunID','ResDir','RanDir','Randbase','Build',
                         'SlimCheck','Focus','BuildPath','ProbScore','AltDis','GablamDis','AADimerFreq','Negatives',
                         'Chance','SeqIn','AltUPC','NewUPC','AddQuery']
        self.boollist = ['BuildProb','DisMask','Force','MaskFreq','SeqOcc','SlimDisc','Teiresias','Test','AlphaHelix',
                        'LogMask','Termini','MaskM','PreAmb','WildVar','Masked','FocusOcc','Masking','Pickup',
                        'EFilter','SlimChance','SlimBuild','Randomise','Extras','TarGZ','CombAmb','SmearFreq',
                        'OccStatsCalculated','Webserver','DNA','ConsMask','DimFreq','OldScores','SigV','SigPrime',
                        'AllSig','PickAll','PickID','FixLen','QExact','CloudFix']
        self.numlist = ['AmbOcc','MinOcc','ProbCut','MinIC','StartTime','MST','WallTime']
        self.intlist = ['AbsMin','MaxWild','SlimLen','AbsMinAmb','TopRanks','MaxSeq','SaveSpace',
                        'FocusOcc','Clouds','MaxUPC','MinWild','HomCut','Extras']
        self.listlist = ['AAMask','Alphabet','Batch','Equiv','FTMask','IMask','MustHave','SigSlim','UP','NewScore',
                         'OccScale','Headers','Query','Warning','SlimCheckExtra']
        self.dictlist = ['AAFreq','DimFreq','Dimers','Slim','SeqOcc','Extremf.','ElementIC','MST','NewScore',
                         'Focus','MotifSeq','MaskPos','AADimerFreq']
        self.objlist = ['SeqList','SlimCheck','SlimList']
        ### Defaults ###
        self._setDefaults(str='None',bool=True,num=0.0,int=0,obj=None,setlist=True,setdict=True)
        self.coreDefaults()
        self.setStr({'CompMask':'5,8','ResFile':'slimfinder.csv','ResDir':rje.makePath('SLiMFinder/'),'Input':'',
                      'RanDir':rje.makePath('Random/'),'Randbase':'rand','Basefile':'','SlimCheck':'',
                      'BuildPath':rje.makePath('SLiMFinder/'),'ProbScore':'None','AltUPC':'None','NewUPC':'None'})
        self.setInt({'AbsMin':3,'MaxWild':2,'SlimLen':5,'AbsMinAmb':2,'TopRanks':1000,'MaxSeq':500,
                     'SaveSpace':0,'Clouds':2,'MaxUPC':0,'MinWild':0,'HomCut':0,'Extras':1})
        self.setNum({'MinOcc':0.05,'ProbCut':0.1,'MinIC':2.1,'StartTime':time.time(),'WallTime':1.0,'AmbOcc':0.05})
        self.setBool({'MaskFreq':True,'SlimDisc':False,'Teiresias':False,'SeqOcc':False,'Extras':True,'Masked':False,
                     'MaskM':True,'PreAmb':True,'WildVar':True,'Randomise':False,'TarGZ':False,'Force':False,
                     'CombAmb':False,'SmearFreq':False,'DisMask':False,'Test':False,'Pickup':False,'AlphaHelix':False,
                     'OccStatsCalculated':False,'Webserver':False,'DNA':False,'ConsMask':False,'DimFreq':True,
                     'OldScores':False,'SigV':False,'SigPrime':False,'AllSig':False,'FixLen':False,
                     'PickAll':True,'PickID':True,'QExact':True,'Palindrome':False,'CloudFix':False})
        self.list['Equiv'] = string.split(default_equiv,',')
        self.list['Batch'] = ['*.dat','*.fas']
        self.dict['MaskPos'] = {2:'A'}
        if self.log.info['Name'] == 'QSLiMFinder':
            self.setStr({'BuildPath':rje.makePath('QSLiMFinder/')})
            self.setStr({'ResDir': self.getStr('BuildPath'),'ResFile':'qslimfinder.csv'})
        ### Other Attributes ###
        self.obj['SlimList'] = rje_slimlist.SLiMList(self.log,self.cmd_list)
        self.obj['SlimList'].obj['SLiMCalc'].setupFilters(slimheaders=[],occheaders=[])    ### Sets up SLiM/Occ Filters
#########################################################################################################################
    def _cmdList(self):     ### Sets Attributes from commandline
        '''
        Sets attributes according to commandline parameters:
        - see .__doc__ or run with 'help' option
        '''
        self.coreCmd()  #!# Remove these core commands from repetition below #!#
        for cmd in self.cmd_list:
            try:
                ### General Options ### 
                self._generalCmd(cmd)
                ### Class Options ###
                self._cmdReadList(cmd,'file',['AAFreq','ResFile','Focus','AltDis','GablamDis','AADimerFreq','Negatives','AddQuery'])
                self._cmdReadList(cmd,'str',['CaseMask','CompMask','RunID','Randbase','SlimCheck','ProbScore'])
                self._cmdReadList(cmd,'path',['ResDir','RanDir','BuildPath','AltUPC','NewUPC'])
                self._cmdReadList(cmd,'int',['AbsMin','MaxWild','SaveSpace','SlimLen','AbsMinAmb','TopRanks','MaxSeq',
                                             'FocusOcc','Clouds','MaxUPC','MinWild','HomCut','Extras'])
                self._cmdReadList(cmd,'num',['AmbOcc','MinOcc','ProbCut','WallTime','MinIC'])
                self._cmdReadList(cmd,'bool',['DisMask','Force','MaskFreq','SeqOcc','Extras','Test','AlphaHelix',
                                             'SlimDisc','Teiresias','LogMask','Randomise','Pickup','Webserver','SigPrime',
                                             'Termini','MaskM','PreAmb','WildVar','EFilter','TarGZ','Masking','SigV',
                                             'CombAmb','SlimChance','SlimBuild','SmearFreq','DimFreq','OldScores',
                                             'AllSig','PickAll','PickID','FixLen','QExact','Palindrome','CloudFix'])
                self._cmdRead(cmd,'bool','PreAmb','ambiguity')
                self._cmdRead(cmd,'num','ProbCut','sigcut')
                self._cmdReadList(cmd,'list',['MustHave','Equiv','Batch','OccScale','Query'])
                self._cmdReadList(cmd,'cdict',['MotifSeq','MaskPos'])
            except: self.errorLog('Problem with cmd:%s' % cmd)
        ### Special Conversion ###
        self.list['MustHave'] = string.split(string.join(self.list['MustHave']).upper())
        if self.getStrLC('Basefile'): self.setStr({'ResFile':'%s.csv' % self.baseFile()})
        if not self.getBool('CloudFix'): self.warnLog('NOTE: cloudfix=F. Be wary of ambiguity over-predictions.')
        if self.getBool('Pickup'): self.setBool({'Append':True})
        if not self.getStrLC('ProbScore') and self.getBool('AllSig'): self.setStr({'ProbScore':'SigPrimeV'})
        if self.getBool('Palindrome'):
            if not self.getBool('DNA'):
                self.errorLog('Palindrome mode only available for DNA sequences. Switched off.', printerror=False)
                self.setBool({'Palindrome': False})
            elif self.getBool('Termini'):
                self.printLog('#TERM','Cannot search for termini in Palindrome mode. Termini=False.')
                self.setBool({'Termini':False})
        if self.getStr('ProbScore') not in basic_headers:
            if self.getStrLC('ProbScore') in probscores: self.setStr({'ProbScore': probscores[self.getStrLC('ProbScore')]})
            else:
                self.errorLog('ProbScore "%s" not recognised: will use SLiMChance significance.' % self.getStr('ProbScore'),printerror=False)
                self.setStr({'ProbScore':'Sig'})
        if self.getInt('MinWild') > self.getInt('MaxWild'):
            self.errorLog('MinWild (%d) > MaxWild (%d): MinWild reduced to %d' % (self.getInt('MinWild'),self.getInt('MaxWild'),self.getInt('MaxWild')),printerror=False)
            self.setInt({'MinWild':self.getInt('MaxWild')})
        if self.getBool('AlphaHelix'):
            self.setInt({'MinWild':0,'MaxWild':3,'SlimLen':4})
        elif self.getBool('FixLen'):
            self.setInt({'MinWild':0}); self.setBool({'WildVar': False})
            self.printLog('#FIX','Will fix motifs at %daa (%d defined; %d wildcard)' % (self.getInt('MaxWild')+self.getInt('SlimLen'),self.getInt('SlimLen'),self.getInt('MaxWild')))
        if self.getBool('Webserver') and self.getNum('ProbCut') > 0.05 and self.getInt('TopRanks') > 100:
            self.printLog('#TOP','Webserver output limited to top 100 SLiMs. (Reduced from %d)' % self.getInt('TopRanks'))
            self.setInt({'TopRanks':100})
        if not self.getBool('DisMask'):
            self.warnLog('No disorder masking. Recommended setting: dismask=T')
            if self.i() >= 0 and rje.yesNo('Switch on disorder masking?'): self.setBool({'DisMask':True})
#########################################################################################################################
    ### <2> ### Simple Stats Methods                                                                                    #
#########################################################################################################################
    def slimNum(self): return len(self.dict['Slim'])
#########################################################################################################################
    def slimUP(self,slim):  ### Returns UP Num for SLiM if in dictionary, else 0.
        if self.dict['Slim'].has_key(slim): return len(self.dict['Slim'][slim]['UP'])
        return 0
#########################################################################################################################
    def slimOccNum(self,slim,upc=None,mstx=True):    ### Returns number of occ of Slim in given UPC
        '''
        Returns number of occ of Slim in given UPC.
        >> slim:str = SLiM to return stats for
        >> upc:list = UPC to be queried. If no UPC, return total Occ count across all sequences.
        >> mstx:bool = Whether to return an MST-adjusted mean (rounded down) or the Max support in the UPC.
        '''
        ## No UPC - return total! ##
        if not upc: return len(self.dict['Slim'][slim]['Occ'])
        ## See if SLiM occurs at all! ##
        if not self.dict['SeqOcc'].has_key(slim): return 0
        ## MST-adjusted mean over UPC ##
        sx = 0; mx = 0.0
        for seq in upc:
            if self.dict['SeqOcc'][slim].has_key(seq):
                mx += self.dict['SeqOcc'][slim][seq]
                sx = max(sx,self.dict['SeqOcc'][slim][seq])     # Old (flawed) way.
        mx = int(mx * self.dict['MST'][upc] / len(upc))
        if mstx: return mx
        return sx
#########################################################################################################################
    def slimSeqNum(self,slim):  ### Returns the number of sequences SLiM occurs in.
        '''Returns the number of sequences SLiM occurs in.'''
        if not self.dict['SeqOcc'].has_key(slim): return 0      ## See if SLiM occurs at all! ##
        return len(self.dict['SeqOcc'][slim])
#########################################################################################################################
    ### <3> ### General Run Methods                                                                                     #
#########################################################################################################################
    def run(self,batch=False):  ### Main SLiMFinder Run Method
        '''
        Main SLiMFinder Run Method:
        0. PreCheck:
            - Check for randomise function and execute if appropriate
        1. Input:
            - Read sequences into SeqList
            - or - Identify appropriate Batch datasets and rerun each with batch=True
        2. SLiMBuild:
            - Check for existing Pickle and load if found. Check appropriate parameter settings and re-run if desired.
            - or - Save sequences as fasta and mask sequences in SeqList
            -  Perform BLAST and generate UPC based on saved fasta.
            - Calculate AAFreq for each sequence and UPC.
            - Find all dimer motifs in dataset using MinWild/MaxWild parameters.
            - Extend to SLiMs and add ambiguity
        5. Identify significant SLiMs.
        6. Output results and tidy files.
        >> batch:bool [False] = whether this run is already an individual batch mode run.
        '''
        try:###~PRECHECK~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ### Randomise Function ###
            if self.getBool('Randomise') and not batch: return self.randomise()

            ###~INPUT~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            #seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None','autofilter=F']
            #self.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            self.setupSeqIn()
            self.setupBasefile()
            self.loadAADimerFreq()
            ## Batch Mode ##
            if not batch and not self.getStrLC('SeqIn'):   # No sequences loaded - use batch mode
                pickup = self.pickup()
                #self.debug(pickup)
                batchfiles = self.batchFiles(pickup)
                if not batchfiles: self.errorLog('No input files found!',printerror=False)
                else:
                    self.setupResults()                 ## Sets up OccStats filter etc. - check against Pickle ##
                    self.backupOrCreateResFile()
                    self.setBool({'Append': True})
                    self.list['Batch'] = []
                    bx = 0
                    for infile in batchfiles:
                        bx += 1
                        if pickup:
                            next = os.path.split(rje.baseFile(infile))[1]
                            if next in pickup:
                                self.printLog('#PICKUP','Skipping batch file %s %s' % (rje.integerString(bx),infile),log=False)
                                continue
                        self.printLog('#BATCH','Batch running %s' % infile)
                        bsf = self.newBatchRun(infile)
                        bsf.dict['AADimerFreq'] = self.dict['AADimerFreq']
                        bsf.run(batch=True)
                        self.printLog('#BATCH','Batch file %s run. Cleaning up for next file.' % infile)
                        del bsf.obj
                        del bsf.list
                        del bsf.dict
                        del bsf
                        self.printLog('#BATCH','|---------- %s run <<<|>>> %s to go -----------|' % (rje.integerString(bx),rje.integerString(len(batchfiles)-bx)),log=False)
                if self.getBool('Win32') and len(sys.argv) < 2: self.verbose(0,0,'Finished!',1) # Optional pause for win32
                return not not batchfiles
            else:
                self.setupResults()                 ## Sets up OccStats filter etc. - check against Pickle ##
                self.backupOrCreateResFile()
            ## Check whether to bother running dataset at all - Check Input versus Min and Max Seq ##
            if 0 < self.getInt('MaxSeq') < self.obj['SeqList'].seqNum():
                self.printLog('#SEQ','%s = %s seqs > Max %s seq. Analysis terminated.' % (self.dataset(),rje.iStr(self.seqNum()),rje.iStr(self.getInt('MaxSeq'))))
                self.serverEnd('MaxSeq',exit=False)
                self.results(aborted='>');
                return False
            if self.seqNum() < self.getInt('MinOcc'):
                self.printLog('#SEQ','Insufficient Sequences (%d) for MinOcc setting (%d). Run aborted.' % (self.seqNum(),self.getInt('MinOcc')))
                self.serverEnd('FewSeq',exit=False); self.results(aborted='<');
                return False

            ###~SLiMBuild~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.setNum({'StartTime':time.time()})
            ## UPC and MinOcc settings: needed to identify correct pickle so must be done first ##
            if not self.makeUPC():
                self.errorLog('Error during makeUPC(). Abandoning %s run' % self.dataset(),printerror=False)
                self.serverEnd('Crash','makeUPC()')
            if not self.setupFocus(): self.serverEnd('Crash','setupFocus()') ## Setup Focus after UPC & MST - also check against Pickle ##
            if not self.setupMinOcc():
                self.printLog('#UPC','Insufficient UPC (%d) for MinOcc setting (%d). Run aborted.' % (self.UPNum(),self.getInt('MinOcc')))
                self.serverEnd('FewUPC',exit=False); self.results(aborted='<'); return False
            if self.getInt('MaxUPC') >= self.getInt('MinOcc') and self.getInt('MaxUPC') < self.UPNum():
                self.printLog('#UPC','Too many UPC (%d) for MaxUPC setting (%d). Run aborted.' % (self.UPNum(),self.getInt('MaxUPC')))
                self.serverEnd('MaxUPC',exit=False)
                self.results(aborted='>');
                return False

            ## Check for existing pickle to replace SLiMBuild portion ##
            pickled = self.searchPickleMe(load=self.getBool('Pickle'))  # Returns appropriate pickled SLiMFinder Object, else None
            if pickled: self = pickled  ## Replace me with my pickle!
            #self.deBug('>>%s' % self.list['Headers'])

            ## Setup Main Results File early in case of user intervention ##
            #x#self.backupOrCreateResFile()

            ## AA Frequency Calculations made early as needed superficially in SLiMBuild ##
            if not pickled: self.maskInput()      ## Mask Input Data - makes info['PreMask'] and info['MaskSeq']
            if self.getBool('MaskFreq'): self.makeAAFreq()
            else:
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                self.makeAAFreq()
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
            self.adjustAATotals()

            ## Execute SLiMBuild if pickle not loaded, else recalculate Bonferroni ##
            if self.getBool('SlimBuild') or self.getBool('SlimChance') or self.getBool('SlimDisc'):
                self.makeBonferroni()   # Estimates and reports total no. motifs in dataset
                if not pickled:     ## Find all dimer motifs in dataset using MaxWild parameters. ##
                    self.makeDimers()       # Makes all ai.{0,x}aj dimers
                    self.reduceDimers()     # Reduces to interesting subset
                    self.makeSLiMs()        # Makes all SLiMs with sufficient support
                    self.searchPickleMe()         # Generates pickling for speedy re-running

            ### Special MotifSeq Output ###
            # This must occur after Input masking but needs no AA Frequencies or SLiMBuild #
            if self.motifSeq() and not self.getBool('SlimBuild'): return

            ###~Post-SLiMBuild Processing/Filtering before SLiMChance and Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
            ## Non-SLiMChance filtering of motifs ##
            self.dict['ElementIC'] = {}
            self.filterSLiMs()
            ## TEIRESIAS Output ##
            if self.getBool('Teiresias') or self.getBool('SlimDisc'): self.teiresias()
            if not self.getBool('SlimChance') and not self.getBool('SlimDisc'):
                self.printLog('#PROB','SlimChance=F and SlimDisc=F : No SLiM probability calculations')
                return

            ###~SLiMChance Probability and Significance Calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.slimChance()

            ###~SLiMFinder Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ### Output results and tidy files. ###
            if not (self.occFilter() or self.statFilter()): self.calculateSLiMOccStats()
            self.obj['SlimList'].combMotifOccStats()
            self.tidyMotifObjects()     # Temporary solution to problem with unknown cause
            self.makeClouds()           # Identifies "clouds" of motifs - similar patterns that overlap
            self.cloudConsensi()        # Generates SLiMMaker consensus motifs
            self.rankScore()            # Converts rankings into Numeric
            self.results()              # Controls SLiMFinder results output
            self.slimCheck()            # Additional SlimCheck Motifs
            if self.extras(0): self.extraOutput()   # MotifList Outputs
            self.tarZipSaveSpace()      # Tarring, Zipping and Saving Space

            ### End ###
            self.printLog('#RES','%s results output to %s and %s.' % (self.prog(),self.getStr('ResFile'),self.getStr('OccResFile')))
            self.printLog('#RES','Additional dataset results output to %s.*' % (self.seqBaseFile()))
            targz = '%s.tgz' % self.runBase()
            if self.getBool('TarGZ') and not self.getBool('Win32') and rje.exists(targz):
                self.printLog('#TGZ','Additional dataset results tarred to %s' % targz)
            if self.getBool('Win32') and len(sys.argv) < 2 and not batch: self.verbose(0,0,'Finished!',1)
            if self.list['SigSlim'] or self.list['SlimCheckExtra']: self.serverEnd('Finished!',exit=False); return True
            else: self.serverEnd('NoSLiM',exit=False); return False
        except KeyboardInterrupt: raise  # Killed
        except SystemExit:
            if self.getNum('WallTime') <= 0 or (time.time() - self.getNum('StartTime')) < (self.getNum('WallTime')*3600): self.serverEnd('Crash'); raise
            if self.list['Headers']: self.results(aborted='!')
            return False # Walltime reached
        except:
            self.errorLog('Error in %s.run()' % self.prog(),printerror=True,quitchoice=False)
            self.serverEnd(endcause='Crash',details='main run')
            return False
#########################################################################################################################
    def newBatchRun(self,infile):   ### Returns SLiMFinder object for new batch run
        '''Returns SLiMFinder object for new batch run.'''
        return SLiMFinder(self.log,self.cmd_list[0:] + ['seqin=%s' % infile,'append=%s' % self.getBool('Append')])
#########################################################################################################################
    def pickup(self):   ### Returns a lists of datasets to skip due to pickup settings                              #V4.3
        '''Returns a lists of datasets to skip due to pickup settings.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not (self.getBool('Pickup') and os.path.exists(self.getStr('ResFile'))): return []     # No pickup needed
            ### ~ [1] Easy Pickups First ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('PickAll'):
                ## ~ [1a] With RunID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getBool('PickID'):
                    pickup = rje.dataDict(self,self.getStr('ResFile'),['RunID'],['Dataset'],lists=True)
                    if self.getStr('RunID') in pickup: pickup = pickup[self.getStr('RunID')]['Dataset']
                    else: pickup = []
                    self.printLog('#PICKUP','Pickup %s %s datasets' % (rje.integerString(len(pickup)),self.getStr('RunID')))
                    return pickup
                ## ~ [1b] All Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:
                    pickup = rje.dataDict(self,self.getStr('ResFile'),['Dataset'],['Pattern'],lists=True).keys()
                    self.printLog('#PICKUP','Pickup %s datasets (No RunID)' % (rje.integerString(len(pickup))))
                    return pickup
            ### ~ [2] Complex Pickup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            else:
                ## ~ [1a] With RunID ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                pickup = []
                if self.getBool('PickID'):
                    pickdat = rje.dataDict(self,self.getStr('ResFile'),['Dataset','RunID'],['Dataset','RunID','Pattern'],lists=True)
                    for k in pickdat.keys()[0:]:
                        pickrun = pickdat.pop(k)
                        if self.getStr('RunID') in pickrun['RunID']:
                            for x in ['<','>','!']:
                                while x in pickrun['Pattern']: pickrun['Pattern'].remove(x)
                            if pickrun['Pattern']: pickup += pickrun['Dataset']
                    self.printLog('#PICKUP','Pickup %s %s datasets (clean runs only)' % (rje.integerString(len(pickup)),self.getStr('RunID')))
                    return pickup
                ## ~ [1b] All Datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                else:
                    pickdat = rje.dataDict(self,self.getStr('ResFile'),['Dataset'],['Dataset','Pattern'],lists=True)
                    for k in pickdat.keys()[0:]:
                        pickrun = pickdat.pop(k)
                        for x in ['<','>','!']:
                            while x in pickrun['Pattern']: pickrun['Pattern'].remove(x)
                        if pickrun['Pattern']: pickup += pickrun['Dataset']
                    self.printLog('#PICKUP','Pickup %s datasets (No RunID, clean runs only)' % (rje.integerString(len(pickup))))
                    return pickup            
        except: self.errorLog('Problem with %s pickup()' % self.prog()); return []
#########################################################################################################################
    ### <4> ### Setup/Input Methods                                                                                     #
#########################################################################################################################
    def setupMinOcc(self):  ### Adjusts MinOcc and AmbOcc settings according to UPNum etc.
        '''Adjusts MinOcc and AmbOcc settings according to UPNum etc.'''
        if self.getNum('AmbOcc') <= 0: self.setNum({'AmbOcc':self.getInt('MinOcc')})
        for occ in ['AmbOcc','MinOcc']:
            if self.getNum(occ) < 1:
                newmin = self.getNum(occ) * self.UPNum()
                if int(newmin) == newmin: newmin = int(newmin)
                else: newmin = int(newmin) + 1
                if occ == 'MinOcc':
                    if newmin <= self.getInt('AbsMin'): newmin = self.getInt('AbsMin')
                    else: self.warnLog('Adjusted MinOcc exceeds absolute minimum (%d) - may be sacrificing power for speed.' % self.getInt('AbsMin'))
                if occ == 'AmbOcc':
                    if newmin <= self.getInt('AbsMinAmb'): newmin = self.getInt('AbsMinAmb')
                    else: self.warnLog('Adjusted AmbOcc exceeds absolute minimum (%d) - may be sacrificing power for speed.' % self.getInt('AbsMinAmb'))
                self.printLog('#MIN','%s %.2f & %d UP => %s %d UPC.' % (occ,self.getNum(occ),self.UPNum(),occ,newmin))
                self.setNum({occ:newmin})
        if self.getInt('AmbOcc') > self.getInt('MinOcc'):
            self.printLog('#MIN','AmbOcc > MinOcc => reduced to %d UPC.' % self.getInt('MinOcc'))
            self.setInt({'AmbOcc':self.getInt('MinOcc')})
        if not self.getBool('PreAmb') and self.getInt('AmbOcc') != self.getInt('MinOcc'):
            self.printLog('#MIN','Preamb=False: AmbOcc -> MinOcc (%d UPC).' % self.getInt('MinOcc'))
            self.setInt({'AmbOcc':self.getInt('MinOcc')})
        if self.UPNum() < self.getInt('MinOcc'): return False
        return True
#########################################################################################################################
    ### <5> ### SLiMBuild Generation Methods                                                                            #
#########################################################################################################################
    def makeDimers(self):   ### Finds all possible dimers with wildcards, using MaxWild stat
        '''Finds all possible dimers with wildcards, using MinWild/MaxWild stat.'''
        try:
            ### Setup ###
            self.dict['Dimers'] = {}
            self.dict['DimFreq'] = {}
            uplist = self.list['UP'][0:]
            nonx = {}   # Total count of non-X positions in UPC
            for upc in uplist + self.seqs():
                nonx[upc] = 0.0
                self.dict['DimFreq'][upc] = [0] * (self.getInt('MaxWild') + 1)
            dx = 0
            sx = 0

            ### Read ###
            for seq in self.seqs():
                ## Setup Sequence ##
                sx += 1
                sequence = seq.info['Sequence'].upper()
                #self.debug(sequence)
                if self.getBool('DNA'): sequence = string.replace(sequence,'N','X')
                if self.getBool('Termini'): sequence = '^%s$' % sequence
                ## Setup UPC and DimFreq ##
                try:
                    upc = self.getUP(seq)
                    if not upc: raise ValueError
                except: self.errorLog('ERR','Something has gone wrong with UPC generation. Redundancy in input data?'); raise
                ## Find Dimers ##
                for i in range(len(sequence)):
                    ## Choose first position and check for wildcard ##
                    r = i
                    if self.getBool('Termini'): r = i - 1
                    ai = sequence[i]
                    if ai in wildcards: continue
                    nonx[upc] += 1
                    nonx[seq] += 1
                    ## Examine each wildcard length in turn ##
                    for x in range(self.getInt('MinWild'),self.getInt('MaxWild')+1):
                        if self.getBool('AlphaHelix') and x == 1: continue
                        j = i + x + 1
                        if len(sequence) <= j: continue
                        aj = sequence[j]
                        if aj in wildcards: continue
                        ## Add Dimer ##
                        self.dict['DimFreq'][seq][x] += 1
                        self.dict['DimFreq'][upc][x] += 1
                        if not self.dict['Dimers'].has_key(ai): self.dict['Dimers'][ai] = {}
                        if not self.dict['Dimers'][ai].has_key(x): self.dict['Dimers'][ai][x] = {}
                        if not self.dict['Dimers'][ai][x].has_key(aj):
                            self.dict['Dimers'][ai][x][aj] = {'UP':[],'Occ':[]}
                            dx += 1
                        if upc not in self.dict['Dimers'][ai][x][aj]['UP']: self.dict['Dimers'][ai][x][aj]['UP'].append(upc)
                        self.dict['Dimers'][ai][x][aj]['Occ'].append((seq,r))
                        newslim = '%s-%s-%s' % (ai,x,aj)   #!# Str->List mod 1.5 #!#
                        if not self.dict['SeqOcc'].has_key(newslim): self.dict['SeqOcc'][newslim] = {seq:1}
                        elif not self.dict['SeqOcc'][newslim].has_key(seq): self.dict['SeqOcc'][newslim][seq] = 1
                        else: self.dict['SeqOcc'][newslim][seq] += 1
                self.progLog('\r#DIM','Reading dimers (%d seq) %s dimers' % (sx,rje.integerString(dx)))
            self.printLog('\r#DIM','Read dimers from %d seq: %s dimers' % (sx,rje.integerString(dx)))
            self.setInt({'Dimers': dx})

            ### Adjust DimFreq ###
            for upc in uplist + self.seqs():
                for x in range(self.getInt('MinWild'),(self.getInt('MaxWild') + 1)):
                    if nonx[upc]: self.dict['DimFreq'][upc][x] = self.dict['DimFreq'][upc][x] / nonx[upc]
                    else:
                        if upc in self.list['UP']: self.printLog('#NONX','WARNING! UPC cluster %d has zero unmasked residues!' % uplist.index(upc))
                        else: self.printLog('#NONX','WARNING! Sequence %s has zero unmasked residues!' % upc.shortName())
                        self.dict['DimFreq'][upc][x] = 0.0

        except:
            self.errorLog('Major problem during makeDimers()')
            raise
#########################################################################################################################
    def reduceDimers(self):     ### Reduces Dimers to those with enough Support
        '''Reduces Dimers to those with enough Support.'''
        try:### Select Dimers ###
            dx = 0
            self.dict['FullDimers'] = {}    ### Store a full copy of dimers for finding later ambiguity!
            for ai in self.dict['Dimers'].keys()[0:]:
                self.dict['FullDimers'][ai] = {}
                for x in self.dict['Dimers'][ai].keys()[0:]:
                    self.dict['FullDimers'][ai][x] = {}
                    for aj in self.dict['Dimers'][ai][x].keys()[0:]:
                        self.dict['FullDimers'][ai][x][aj] = self.dict['Dimers'][ai][x][aj]
                        ox = len(self.dict['Dimers'][ai][x][aj]['UP'])
                        if ox < self.getInt('AmbOcc'): self.dict['Dimers'][ai][x].pop(aj)
                        else: dx += 1
                    self.progLog('\r#DIM','Reducing dimers: %s >= %d of %d UPC ' % (rje.integerString(dx),self.getInt('AmbOcc'),self.UPNum()))
                    if not self.dict['Dimers'][ai][x]: self.dict['Dimers'][ai].pop(x)
                if not self.dict['Dimers'][ai]: self.dict['Dimers'].pop(ai)
            self.printLog('\r#DIM','Reducing dimers: %s >= %d of %d UPC ' % (rje.integerString(dx),self.getInt('AmbOcc'),self.UPNum()))
        except: self.errorLog('Problem reducing Dimers to AmbOcc+')
#########################################################################################################################
    def makeSLiMs(self):    ### Makes SLiMs with enough support from Dimers
        '''Makes SLiMs with enough support from Dimers.'''
        try:
            ### Setup ###
            self.dict['Slim'] = {}
            prevslim = []

            ### Select Dimers ###
            for ai in self.dict['Dimers']:
                for x in self.dict['Dimers'][ai]:
                    if self.getBool('AlphaHelix') and x < 2: continue   # Cannot start with x=0
                    for aj in self.dict['Dimers'][ai][x].keys()[0:]:
                        slim = '%s-%s-%s' % (ai,x,aj)
                        prevslim.append(slim)
                        self.dict['Slim'][slim] = self.dict['Dimers'][ai][x][aj]
                    self.progLog('\r#SLIM','Selecting 2aa SLiMs: %s >= %d of %d UPC ' % (rje.integerString(self.slimNum()),self.getInt('AmbOcc'),self.UPNum()))
            self.printLog('\r#SLIM','Selecting 2aa SLiMs: %s >= %d of %d UPC' % (rje.integerString(self.slimNum()),self.getInt('AmbOcc'),self.UPNum()))

            ### Ambiguity ###
            self.ambSLiM(prevslim)

            ### Extend ###
            for f in range(3,self.getInt('SlimLen')+1):
                if not prevslim: break
                newslim = []
                ex = 0.0
                for slim in prevslim:
                    self.progLog('\r#SLIM','Extending %daa SLiMs >= %d of %d UPC: %.1f%%' % (f,self.getInt('AmbOcc'),self.UPNum(),ex/len(prevslim)))
                    ex += 100.0
                    newslim += self.extendSLiM(slim)    # Why was there "if self.opt['PreAmb']" ???
                    self.wallTime()
                prevslim = newslim
                self.printLog('\r#SLIM','Extending %daa SLiMs >= %d of %d UPC: %s SLiMs' % (f,self.getInt('AmbOcc'),self.UPNum(),rje.integerString(len(prevslim))))
                ## Add ambiguity ##
                self.ambSLiM(prevslim)
            self.printLog('#SLIM','%s SLiMs >= %d of %d UPC' % (rje.integerString(self.slimNum()),self.getInt('AmbOcc'),self.UPNum()))
        except SystemExit: raise
        except:
            self.errorLog('Fatal error making SLiMs from Dimers')
            raise
#########################################################################################################################
    def extendSLiM(self,slim):  ### Finds and returns extensions of SLiM with sufficient support
        '''
        Finds and returns extensions of SLiM with sufficient support.
        >> slim:str = SLiM to extend (using dimers)
        '''
        try:### Setup ###
            ai = slim[-1]
            if not self.dict['Dimers'].has_key(ai): return []
            extend = []
            slimocc = self.dict['Slim'][slim]['Occ']
            if self.getBool('AlphaHelix'): prevx = string.atoi(string.split(slim,'-')[-2])
            elif self.getBool('FixLen'):
                prevx = 0
                for x in string.split(slim,'-'):
                    try: prevx += string.atoi(x)
                    except: pass

            ### Try dimers ###
            for x in self.dict['Dimers'][ai]:
                if self.getBool('AlphaHelix') and (x == prevx or (prevx + x) == 3): continue                 # and slim.find('-%d-' % x) > 0: continue   # Not i,i+3/4,i+7
                elif self.getBool('FixLen') and (prevx + x) > self.getInt('MaxWild'): continue
                for aj in self.dict['Dimers'][ai][x]:   # No longer check minocc as removed in prev makeSLiMs
                    newslim = slim + '-%s-%s' % (x,aj)   #!# Str->List mod 1.5 #!#
                    newocc = []
                    newup = []
                    for (seq,pos) in slimocc:
                        if (seq,pos+slimLen(slim)-1) in self.dict['Dimers'][ai][x][aj]['Occ']:
                            newocc.append((seq,pos))
                            upc = self.getUP(seq)
                            if upc not in newup: newup.append(upc)
                    if len(newup) >= self.getInt('AmbOcc'):
                        extend.append(newslim)
                        self.dict['Slim'][newslim] = {'Occ':newocc,'UP':newup}
                        for (seq,pos) in newocc:
                            if not self.dict['SeqOcc'].has_key(newslim): self.dict['SeqOcc'][newslim] = {seq:1}
                            elif not self.dict['SeqOcc'][newslim].has_key(seq): self.dict['SeqOcc'][newslim][seq] = 1
                            else: self.dict['SeqOcc'][newslim][seq] += 1
            ### Return ###
            return extend
        except:
            try: self.errorLog('Problem extending SLiM "%s"' % slim,quitchoice=True)
            except: raise KeyboardInterrupt
        return []
#########################################################################################################################
    def ambSLiM(self,prevslim): ### Combines SLiMs from prevslim list into ambiguous SLiMs
        '''Combines SLiMs from prevslim list into ambiguous SLiMs.'''
        try:
            ### Setup ###
            if not self.getBool('PreAmb') or (not self.getBool('WildVar') and not self.list['Equiv']): return
            ### Wildcard ###
            if self.getBool('WildVar'):
                wildamb = ''
                for i in range(self.getInt('MinWild'),self.getInt('MaxWild')+1): wildamb += str(i)
                self.addAmb(prevslim,[wildamb],'Wildcard')
            ### Degeneracy ###
            self.addAmb(prevslim,self.list['Equiv'],'Degeneracy')
            if self.getBool('CombAmb') and self.list['Equiv'] and self.getBool('WildVar'):
                self.addAmb(prevslim,self.list['Equiv']+[wildamb])
        except: self.errorLog('Major problem with %s.ambSLiM()' % self.prog())
#########################################################################################################################
    def addAmb(self,prevslim,equivlist,type='Combined'): ### Combines SLiMs from prevslim list into ambiguous SLiMs
        '''Combines SLiMs from prevslim list into ambiguous SLiMs.'''
        try:
            ### Setup ###
            if not equivlist or not prevslim: return

            ### Take each SLiM in turn ###
            sx = 0.0
            ax = 0
            for slim in prevslim:
                self.wallTime()
                self.progLog('\r#AMB','Adding %s ambiguity: %.1f%% (%s amb motifs)' % (type,sx/len(prevslim),rje.integerString(ax)))
                sx += 100.0
                if not self.slimFocus(slim): continue  # Don't look for ambiguity in non-focal SLiMs
                ## Check IC - if not enough then degenerate SLiM will be worse! ##
                if self.slimIC(slim) < self.getNum('MinIC'): continue
                ## Make options for each position ##
                ambopt = []
                for pos in string.split(slim,'-'):
                    newamb = []
                    for equiv in equivlist:
                        if equiv.find(pos) >= 0: newamb.append(equiv)
                    if not newamb: newamb = [pos]
                    ambopt.append(newamb)
                ## Make List of option combos ##
                combos = rje.listCombos(ambopt)
                ## Assess each combo for possible ambiguous Motifs ##
                for combo in combos:
                    varlist = rje.listCombos(combo)
                    for v in range(len(varlist)): varlist[v] = string.join(varlist[v],'-')
                    ## Order so that least expected comes first ##
                    varexp = {}
                    varup = {}
                    varlist.remove(slim)
                    for var in varlist[0:]:
                        if self.slimUP(var) < self.getInt('AmbOcc'):    # Insufficient support
                            varlist.remove(var)
                            continue
                        varexp[var] = rje_motif.expect(patternFromCode(var),self.dict['AAFreq']['Dataset'],self.dict['AAFreq']['Dataset']['Total'],1)
                        varup[var] = self.dict['Slim'][var]['UP'][0:]
                        for upc in self.dict['Slim'][slim]['UP'][0:]:
                            if upc in varup[var]: varup[var].remove(upc)
                        if not varup[var]: varlist.remove(var)   # No new UP
                    for i in range(len(varlist)):
                        for j in range(i):
                            if varexp[varlist[i]] <= varexp[varlist[j]]:
                                (varlist[i],varlist[j]) = (varlist[j],varlist[i])
                    ## Order by Total UP / Total Occ ##
                    for i in range(len(varlist)):
                        for j in range(i):
                            if self.slimOccNum(varlist[i]) >= self.slimOccNum(varlist[j]):
                                (varlist[i],varlist[j]) = (varlist[j],varlist[i])
                    for i in range(len(varlist)):
                        for j in range(i):
                            if self.slimUP(varlist[i]) >= self.slimUP(varlist[j]):
                                (varlist[i],varlist[j]) = (varlist[j],varlist[i])
                    ## Order by no. of extra UP added and only add if adding UP ##
                    occvar = [slim]     # List of variants with sufficient support #
                    while varlist:
                        ## Order by extra UP ##
                        for i in range(len(varlist)):
                            for j in range(i):
                                if (len(varup[varlist[i]])/float(slimDif(slim,varlist[i]))) >= (len(varup[varlist[j]])/float(slimDif(slim,varlist[j]))):
                                    (varlist[i],varlist[j]) = (varlist[j],varlist[i])
                        ## Add best ##
                        occvar.append(varlist.pop(0))
                        ## Reduce UP accordingly ##
                        for var in varlist[0:]:
                            for upc in self.dict['Slim'][occvar[-1]]['UP'][0:]:
                                if upc in varup[var]: varup[var].remove(upc)
                            if not varup[var]: varlist.remove(var)     # No new UP
                    ## Combine Variants into Ambiguous Motif ##
                    ambslim = self.combineAmb(occvar)
                    if ambslim: ax += 1
            self.printLog('\r#AMB','Adding %s ambiguity: %.1f%% (%s amb motifs)' % (type,sx/len(prevslim),rje.integerString(ax)))
        except SystemExit: raise
        except: self.errorLog('Major problem with %s.addAmb()' % self.prog())
#########################################################################################################################
    def combineAmb(self,occvar):    ### Combines motif variants in occvar into an ambiguous motif
        '''Combines motif variants in occvar into an ambiguous motif, incorporating "missing" variants as appropriate.'''
        ### ~ [1] Combine slim variants into a new SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if len(occvar) < 2: return None     # Nothing to combine
        ## ~ [1a] Make new SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        newslim = string.split(occvar[0],'-')
        for var in occvar[1:]:      ## Add variants to each position ##
            varsplit = string.split(var,'-')
            for i in range(len(varsplit)):
                pos = varsplit[i]
                if newslim[i].find(pos) < 0: newslim[i] = newslim[i] + pos
        ## ~ [1b] Expand wildcards and Reformat ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        slim = ''
        for pos in newslim:
            newpos = rje.strSort(pos)
            if rje.matchExp('(\d)\d*(\d)',newpos):
                (imin,imax) = rje.matchExp('(\d)\d*(\d)',newpos)
                newpos = ''
                for x in range(int(imin),int(imax)+1): newpos += '%d' % x
            slim += newpos + '-'
        slim = slim[:-1]
        ### ~ [2] Check for existence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        if self.dict['Slim'].has_key(slim): return None
        ### ~ [3] Add New SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        self.dict['Slim'][slim] = {'Occ':[],'UP':[]}
        ## ~ [3a] Occurrences and UPC for all variants ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        for var in rje.listCombos(string.split(slim,'-')):
            varslim = string.join(var,'-')
            if varslim not in self.dict['Slim']: continue
            for occ in self.dict['Slim'][varslim]['Occ']:
                if occ not in self.dict['Slim'][slim]['Occ']: self.dict['Slim'][slim]['Occ'].append(occ)
            for upc in self.dict['Slim'][varslim]['UP']:
                if upc not in self.dict['Slim'][slim]['UP']: self.dict['Slim'][slim]['UP'].append(upc)
        ## ~[3b] SeqOcc dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        self.dict['SeqOcc'][slim] = {}
        for (seq,pos) in self.dict['Slim'][slim]['Occ']:
            if not self.dict['SeqOcc'][slim].has_key(seq): self.dict['SeqOcc'][slim][seq] = 1
            else: self.dict['SeqOcc'][slim][seq] += 1
        ## ~[3c] Return new SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
        return slim
#########################################################################################################################
    ### <6> ### SLiM Filtering Methods                                                                                  #
#########################################################################################################################
    def motifOccStats(self,slim):    ### Calculates Motif OccStats and Filters Occ if appropriate
        '''Calculates Motif OccStats and Filters Occ if appropriate. Returns Motif/None if slim still OK/Not.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            Motif = self.addSLiMToList(slim)
            if not self.getBool('OccStatsCalculated'):
                slimocc = Motif.dict['Occ']
                for seq in slimocc: self.obj['SlimList'].obj['SLiMCalc'].occStats(slimocc[seq],xpad=0,progress=False,silent=True)
            ### ~ [2] Filter Occurrences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            prefiltx = Motif.occNum()
            Motif.occFilter(self.occFilter())
            ## ~ [2a] Redefine Occ and slimUP if necessary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if Motif.occNum() < prefiltx:  ### Occurrences have been lost
                self.dict['Slim'][slim]['UP'] = []
                self.dict['Slim'][slim]['Occ'] = []
                for Seq in Motif.dict['Occ']:
                    for occ in Motif.dict['Occ'][Seq]:
                        pos = occ['Pos'] - 1
                        self.dict['Slim'][slim]['Occ'].append((Seq,pos))
                        upc = self.getUP(Seq)
                        if upc not in self.dict['Slim'][slim]['UP']: self.dict['Slim'][slim]['UP'].append(upc)
            ## ~ [2b] Check for loss of Motif ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.slimUP(slim) < self.getInt('MinOcc'):
                self.dict['Slim'].pop(slim)
                self.obj['SlimList'].removeMotif(Motif)
                return None
            return Motif
        except:
            self.errorLog('Problem with SLiMCore.motifOccStats(%s)' % slim)
            return False
#########################################################################################################################
    def _remakeSLiMDictUP(self):   ### Remakes SLiM UP dictionaries from SLiM Occ
        '''Remakes SLiM UP dictionaries from SLiM Occ.'''
        sx = 0.0; stot = self.slimNum()
        for slim in self.dict['Slim'].keys():
            self.progLog('#UPX','Redefining SLiM UP Support: %.2f%%' % (sx/stot)); sx += 100.0
            self.dict['Slim'][slim]['UP'] = []
            for (Seq,pos) in self.dict['Slim'][slim]['Occ']:
                upc = self.getUP(Seq)
                if upc not in self.dict['Slim'][slim]['UP']: self.dict['Slim'][slim]['UP'].append(upc)
        self.printLog('#UPX','Redefined SLiM UP Support from SLiM Occ.')
#########################################################################################################################
    def filterSLiMs(self):  ### Filters SLiMs on non-SLiMChance parameters
        '''Filters SLiMs on non-SLiMChance parameters.'''
        try:### ~ [1] Reduce to MinOcc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('NewUPC'):
                ufile = '%s%s.upc' % (self.getStr('NewUPC'),self.dataset())
                try:
                    self.loadUPCFromFile(ufile)
                    if self.getBool('MaskFreq'): self.makeAAFreq()
                    else:
                        for seq in self.seqs(): seq.info['Sequence'] = seq.info['PreMask'][0:]
                        self.makeAAFreq()
                        for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                    self.adjustAATotals()
                    self._remakeSLiMDictUP()
                except: self.errorLog('Problem during NewUPC loading')
            (sx,stot) = (0.0,self.slimNum())
            self.dict['AllSlims'] = {}
            for slim in self.dict['Slim'].keys()[0:]:
                self.dict['AllSlims'][slim] = self.dict['Slim'][slim]
                self.progLog('\r#SLIM','Filtering %s SLiMs: %.1f%%' % (rje.integerString(stot),sx/stot)); sx += 100.0
                ### Full Length ###
                if self.getBool('FixLen') or self.getBool('AlphaHelix'):
                    totlen = 0
                    for x in string.split(slim,'-'):
                        try: totlen += string.atoi(x)
                        except: totlen += 1
                    #self.deBug('%s = %d' % (slim,totlen))
                    if self.getBool('AlphaHelix') and totlen != 8: self.dict['Slim'].pop(slim); continue
                    elif self.getBool('FixLen'):
                        if totlen != (self.getInt('SlimLen') + self.getInt('MaxWild')): self.dict['Slim'].pop(slim); continue
                        #i# Should only keep SLiMs that have gone through all Slimlen cycles and added maxwild Xs.
                ### Palindrome Mode ###
                if self.getBool('DNA') and self.getBool('Palindrome'):
                    slimpos = string.split(slim,'-')
                    for i in range(len(slimpos)):
                        try: x = '%d' % int(slimpos[i])
                        except: x = rje.strSort(rje_sequence.reverseComplement(slimpos[i]))
                        if x != slimpos[-(i+1)]: self.dict['Slim'].pop(slim); break
                    if slim not in self.dict['Slim']: continue
                    self.deBug(slim)
                ### Check Query ###
                if not self.slimFocus(slim): self.dict['Slim'].pop(slim)#; self.deBug('slimFocus')
                ### Check IC ###
                elif self.slimIC(slim) < self.getNum('MinIC'): self.dict['Slim'].pop(slim)#; self.deBug('slimIC')
                ### Check MustHave ###
                elif not self.mustHave(slim): self.dict['Slim'].pop(slim)#; self.deBug('slimMustHave')
                ### Check Occ ##
                elif self.slimUP(slim) < self.getInt('MinOcc'): self.dict['Slim'].pop(slim)#; self.deBug('slimOcc')
                #self.deBug('%s: %s' % (slim,slim in self.dict['Slim']))
            ### ~ [2] Additional Motif Stats if used for OccFilter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.occFilter():
                self.calculateSLiMOccStats()
                for slim in self.dict['Slim'].keys()[0:]:
                    if not self.motifOccStats(slim): self.dict['Slim'].pop(slim)
            self.printLog('\r#SLIM','Filtering %s SLiMs: %s retained.' % (rje.integerString(stot),rje.integerString(self.slimNum())))
            ### ~ [3] AltUPC Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('AltUPC'):
                ## ~ [3a] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                my_uplist = self.list['UP']
                my_mstdict = self.dict['MST']
                my_mst = self.getNum('MST')
                ufile = '%s.upc' % self.seqBaseFile('AltUPC')
                try: ## ~ [3b] Filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                    self.loadUPCFromFile(ufile)
                    (sx,stot) = (0.0,self.slimNum())
                    for slim in self.dict['Slim'].keys()[0:]:
                        self.progLog('\r#ALT','AltUPC filtering %s SLiMs: %.1f%%' % (rje.integerString(stot),sx/stot)); sx += 100.0
                        slimocc = self.dict['Slim'][slim]['Occ']
                        newup = []
                        for (seq,pos) in slimocc:
                            upc = self.getUP(seq)
                            if upc not in newup: newup.append(upc)
                        if len(newup) < self.getInt('MinOcc'): self.dict['Slim'].pop(slim)
                    self.printLog('\r#ALT','AltUPC filtering %s SLiMs: %s retained.' % (rje.integerString(stot),rje.integerString(self.slimNum())))
                except: self.errorLog('Problem during AltUPC filtering')
                ## ~ [3c] Restore data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.list['UP'] = my_uplist
                self.dict['MST'] = my_mstdict
                self.setNum({'MST': my_mst})
        except:
            self.errorLog('Fatal error filtering SLiMs before SLiMChance')
            raise
#########################################################################################################################
    def setupFocus(self):   ### Sets up Focus dictionary
        '''Sets up Focus dictionary. Returns True if OK, else False (which cancels run).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Focus'] = {}
            if not self.list['Query'] and not self.getStrLC('Focus'): return True   # No focus
            if self.getStrLC('Focus') and not os.path.exists(self.getStr('Focus')):
                self.errorLog('Cannot open Focus file "%s". Quitting run.' % self.getStr('Focus'),printerror=False)
                return False
            myseq = self.obj['SeqList'].seq[0:]
            fx = 0
            ### ~ [1] ~ Focus Dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('Focus'):
                self.obj['SeqList'].opt['LogRem'] = False
                flines = self.loadFromFile(self.getStr('Focus'),chomplines=True)
                (group,type,queries) = ('','Seq',[])
                for fline in flines:
                    if rje.matchExp('^#+(\S+):(\S+)',fline):    #Group:Type (Seq/Spec/Desc/DB/Acc)
                        (group,type) = rje.matchExp('^#+(\S+):(\S+)',fline)
                    elif fline[:2] == '//' and group: # Update
                        self.obj['SeqList']._filterCmd(['good%s=%s' % (type.lower(),string.join(queries,','))])
                        self.obj['SeqList']._filterSeqs()
                        self.printLog('#FOCUS','Group "%s" = %s: %d sequences identified' % (group,type,self.obj['SeqList'].seqNum()))
                        if self.obj['SeqList'].seqNum() < 1:
                            self.errorLog('Focus group "%s" (%s) mapped to no sequences.' % (group,type),printerror=False)
                            continue
                        fx += self.obj['SeqList'].seqNum()
                        self.dict['Focus'][group] = self.obj['SeqList'].seq[0:]
                        self.obj['SeqList'].seq = myseq[0:]
                        (group,type,queries) = ('','Seq',[])
                    elif fline[:1] not in ['#','']: queries.append(fline)
            ## ~ [1a] ~ Queries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if self.list['Query']:
                self.dict['Focus']['Query'] = []
                for qry in self.list['Query']:
                    if self.obj['SeqList'].querySeq(qry):
                        qseq = self.obj['SeqList'].obj['QuerySeq']
                        if qseq not in self.dict['Focus']['Query']:
                            fx += 1
                            self.dict['Focus']['Query'].append(qseq)
                if self.dict['Focus']['Query']: self.printLog('#QRY','%d sequences mapped to focal group "Query"' % len(self.dict['Focus']['Query']))
                else:
                    self.printLog('#ERR','No sequences mapped for focal group "Query"')
                    self.dict['Focus'].pop('Query')
            ## ~ [1b] ~ FocusUPC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.dict['FocusUPC'] = {}
            for grp in self.dict['Focus']:
                self.dict['FocusUPC'][grp] = []
                for seq in self.dict['Focus'][grp]:
                    u = self.getUP(seq)
                    if u and u not in self.dict['FocusUPC'][grp]: self.dict['FocusUPC'][grp].append(u)

            ### ~ [2] ~ Finish ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SeqList'].seq = myseq
            self.printLog('#FOCUS','Mapped %d focal sequences in %d group(s)' % (fx,len(self.dict['Focus'])))
            return True
        except:
            self.errorLog('Problem setting up Focus for %s' % self.dataset())
            return False
#########################################################################################################################
    def slimFocus(self,slim):   ### Returns True if slim if Focal sequence groups, else False
        '''Returns True if slim if Focal sequence groups, else False.'''
        ### Setup ###
        if not self.dict['Focus']: return True
        maxfail = 0
        if self.getInt('FocusOcc') > 0: maxfail = len(self.dict['Focus']) - self.getInt('FocusOcc')
        slimgrp = self.dict['Focus'].keys()     # Groups not accounted for
        for (seq,occ) in self.dict['Slim'][slim]['Occ']:
            for grp in slimgrp[0:]:
                if seq in self.dict['Focus'][grp]: slimgrp.remove(grp)
        if len(slimgrp) > maxfail: return False    # Too many group(s) not accounted for by occs
        return True
#########################################################################################################################
    def mustHave(self,slim):    ### Looks at SLiM w.r.t. MustHave list and returns True/False if OK or not
        '''Looks at SLiM w.r.t. MustHave list and returns True/False if OK or not.'''
        if not self.list['MustHave']: return True
        for a in self.list['MustHave']:
            if slim.count(a) > 0: return True
        return False
#########################################################################################################################
    ### <7> ### SLiMChance Probability Methods                                                                          #
#########################################################################################################################
    def loadAADimerFreq(self):  ### Calculates/loads AA Dimer Frequencies
        '''Calculates/loads AA Dimer Frequencies.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['AADimerFreq'] = {}   # Empirical dimer frequencies {i:{x:{j:freq}}}
            if not self.getStrLC('AADimerFreq'): return
            if not os.path.exists(self.getStr('AADimerFreq')): raise IOError
            ## ~ [1a] New SLiMBuild (SLiMFinder) object ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slimbuild = self.newBatchRun(self.getStr('AADimerFreq'))
            seqcmd = ['gnspacc=T','usecase=T'] + self.cmd_list + ['autoload=T','query=None']
            slimbuild.obj['SeqList'] = rje_seq.SeqList(self.log,seqcmd)
            slimbuild.setupBasefile()
            slimbuild.setInt({'MinOcc':1,'MinWild':0,'MaxWild':max(self.getInt('MaxWild'),4)})
            slimbuild.setBool({'EFilter':False})
            aadfile = '%s.w%d.aadimer.tdt' % (rje.baseFile(self.getStr('AADimerFreq'),self.getInt('MaxWild')))

            ### ~ [2] Check for and load *.*.aadimer.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if os.path.exists(aadfile):     ### Load dimer frequencies
                loadfreq = rje.dataDict(self,aadfile,['dimer'])     # {dimer:{wild:freq}}
                ## ~ [2a] Convert to dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                for d in loadfreq:
                    (i,j) = (d[0],d[1])
                    if i not in self.dict['AADimerFreq']: self.dict['AADimerFreq'][i] = {}
                    for x in loadfreq[d]:
                        if x not in self.dict['AADimerFreq'][i]: self.dict['AADimerFreq'][i][x] = {}
                        self.dict['AADimerFreq'][i][x][j] = string.atof(loadfreq[d][x])
                return 

            ### ~ [3] Make dimer dictionary and save to *.*.aadimer.tdt ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not slimbuild.makeUPC(): raise ValueError     #!# Check this
            slimbuild.maskInput()    # Dataset is masked - run through first without masking if that's what is wanted
            slimbuild.makeDimers()       # Makes all ai.{0,x}aj dimers
            ## ~ [3a] Convert SLiMBuild dimer frequencies to counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            dimercount = {}
            for i in slimbuild.dict['Dimers']:    # Main dimer dictionary {Ai:{X:{Aj:{'UP':[UPC],'Occ':[(Seq,Pos)]}}}}
                dimercount[i] = {}
                for x in slimbuild.dict['Dimers'][i]:
                    dimercount[i][x] = {}
                    for j in slimbuild.dict['Dimers'][i][x]:
                        dimercount[i][x][j] = 0.0
                        for (seq,pos) in slimbuild.dict['Dimers'][i][x][j]['Occ']:
                            upc = self.getUP(seq)
                            dimercount[i][x][j] += self.dict['MST'][upc]       # Add MST-weighted occurrence
            ## ~ [3b] Convert counts to AADimerFreq dictionary ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for i in dimercount:
                self.dict['AADimerFreq'][i] = {}
                for x in dimercount[i]:
                    self.dict['AADimerFreq'][i][str(x)] = {}
                    dx = sum(dimercount[i][x].values())
                    for j in slimbuild.dict['Dimers'][i][x]:
                        self.dict['AADimerFreq'][i][str(x)][j] = dimercount[i][x][j] / dx
            ## ~ [3c] Output AADimerFreq file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            aalist = rje_slim.default_aas
            aahead = ['dimer']
            for x in range(slimbuild.stat['MaxWild']): aahead.append(str(x))
            rje.delimitedFileOutput(self,aadfile,aahead,'\t',{})
            for i in aalist:
                for j in aalist:
                    datadict = {'dimer':'%s%s' % (i,j)}
                    for x in aahead[1:]: datadict[x] = '%.5f' % self.dict['AADimerFreq'][i][x][j]
                    rje.delimitedFileOutput(self,aadfile,aahead,'\t',datadict)
            
        except:
            self.errorLog('Error in %s.makeAADimerFreq(%s).' % (self.prog(),self.getStr('AADimerFreq')))
            if self.obj['Interactive'] < 0 or not rje.yesNo('Continue without using AADimerFreq?'): raise 
#########################################################################################################################
    def adjustAATotals(self):   ### Adjusts AA Totals following masking
        '''Adjusts AA Totals following masking and makes appropriate frequencies. If maskfreq=T then the amino acid counts
        will be exactly as they were. If maskfreq=F, however, frequencies will need to be adjust for the new number of
        masked and non-masked amino acids.'''
        try:
            ### Simple PreMasking Procedure ###
            if self.getBool('MaskFreq'):
                for upc in self.list['UP']:
                    if self.dict['AAFreq'][upc].has_key('X'): self.dict['AAFreq'][upc].pop('X') # Ignore Xs
                    if self.getBool('DNA') and self.dict['AAFreq'][upc].has_key('N'): self.dict['AAFreq'][upc].pop('N') # Ignore Ns from DNA sequences
                    self.dict['AAFreq'][upc].pop('Total')  #!# Total remade by dictFreq #!#
                    if self.getBool('Termini'):
                        self.dict['AAFreq'][upc]['^'] = len(upc)
                        self.dict['AAFreq'][upc]['$'] = len(upc)
                    rje.dictFreq(self.dict['AAFreq'][upc])
                    ## Make MST adjustments for UPC ##
                    self.dict['AAFreq'][upc]['Total'] = int(0.5+(self.dict['AAFreq'][upc]['Total']*self.dict['MST'][upc]))
                if self.getBool('SmearFreq'): self.smearAAFreq()
                return

            ### More complicated PostMasking Procedure ###
            (prex,postx) = (0,0)
            for upc in self.list['UP']:
                x = 0
                if self.dict['AAFreq'][upc].has_key('X'): x = self.dict['AAFreq'][upc].pop('X') # Ignore Xs
                if self.getBool('DNA') and self.dict['AAFreq'][upc].has_key('N'): x += self.dict['AAFreq'][upc].pop('N') # Ignore Ns from DNA sequences
                totalaa = self.dict['AAFreq'][upc].pop('Total')
                preaa = totalaa - x   # Want to calculate new total
                prex += preaa
                nonx = 0.0
                for seq in upc: nonx += seq.aaLen() - seq.info['Sequence'].count('X')
                rje.dictFreq(self.dict['AAFreq'][upc])
                ## Termini ##
                if self.getBool('Termini'):
                    nonx += 2 * len(upc)
                    self.dict['AAFreq'][upc]['^'] = len(upc) / nonx
                    self.dict['AAFreq'][upc]['$'] = len(upc) / nonx
                ## Make MST adjustments for UPC ##
                self.dict['AAFreq'][upc]['Total'] = int(0.5+(nonx*self.dict['MST'][upc])) 
                postx += self.dict['AAFreq'][upc]['Total']
                    
            ### Finish ###
            if self.getStrLC('AAFreq'): prex = self.dict['AAFreq']['Dataset']['Total']
            if self.getBool('DNA'): self.printLog('#ADJ','Effective dataset size reduced from %s nt to %s nt' % (rje.integerString(prex),rje.integerString(postx)))
            else: self.printLog('#ADJ','Effective dataset size reduced from %s AA to %s AA' % (rje.integerString(prex),rje.integerString(postx)))
            if self.getBool('SmearFreq'): self.smearAAFreq()
        except:
            self.errorLog('Problem during %s.adjustAATotals()' % self.prog())
            raise
#########################################################################################################################
    def slimProb(self,slim): ### Calculate Probabilities for given SLiM
        '''Calculate Probabilities for given SLiM.'''
        try:
            ###~Calculate prob of 1+ occ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            p1 = {}         # Dictionary of {upc:chance of 1+ occ in upc}
            ##~~Setup pattern and variable-lenght multiplier~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            poslist = []    # List of AA positions in SLiM
            wildlist = []   # List of wildcard lengths in SLiM
            wild = False    # Whether next part is a wildcard length
            mult = 1        # Variable-length multiplier
            minslimlen = 0  # Minimum SLiM length
            for part in string.split(slim,'-'):      # Split SLiM code in components
                ## Update lists ##
                if wild: wildlist.append(part)
                else: poslist.append(part); minslimlen += 1
                ## Calculate multiplier ##
                if wild:
                    (minx,maxx) = (self.getInt('MaxWild'),0)
                    for x in part:
                        minx = min(minx,int(x))
                        maxx = max(maxx,int(x))
                    mult *= (int(maxx) - int(minx) + 1)
                    minslimlen += minx
                wild = not wild
            ##~~Calculate p1+ for each UPC~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            for upc in self.list['UP']:
                if self.dict['AADimerFreq']: (k,N,p) = self.aaDp1(slim,upc)
                else:
                    ## Setup  parameters for binomial ##
                    N = self.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
                    p = 1.0                                 # Probability of SLiM at each position
                    k = 1                                   # Number of successful trials (occurrences)
                    if self.getBool('SeqOcc'): k = max(1,self.slimOccNum(slim,upc))
                    ## Calculate p and N from AAFreq and DimFreq ##
                    for pos in poslist:     # AA position
                        posfreq = 0.0
                        for aa in pos:
                            if aa not in self.dict['AAFreq'][upc]:
                                self.warnLog('AAFrequency for %s not found in %s UPC. Will be given arbitrary freq 0.001.' % (aa,upc[0].shortName()),warntype='aamissing',quitchoice=True,suppress=True)
                                self.dict['AAFreq'][upc][aa] = 0.001
                            posfreq += rje.getFromDict(self.dict['AAFreq'][upc],aa,returnkey=False,default=0.0)  # Options for ambiguity
                        p *= posfreq
                    if self.getBool('DimFreq'):
                        for dim in wildlist:    # DimerFreq
                            dimfreq = 0.0
                            for x in dim:
                                try: dimfreq += self.dict['DimFreq'][upc][int(x)]   # Options for wildcard length
                                except: pass
                            N *= (dimfreq / len(dim))       # Mutliply by mean dimer frequency
                    else: N -= ((minslimlen-1) * self.dict['MST'][upc])
                    N = max(0,N)
                    N *= mult       # Each length variant is effectively another position the SLiM could occur
                    if p > 1: p = 1.0   # Cannot in reality have p > 1!
                    ## Calculate binomial ##
                    p1[upc] = rje.binomial(k,N,p,usepoisson=False,callobj=self)
                    #if slim == 'K-0-L-0-Y': open('kly.tmp','a').write('%s::SF| k = %d; p = %s; N = %d; p1+ = %s\n' % (upc[0].shortName(),k,p,N,p1[upc]))                        
            ## Extra verbosity. Remove at some point? ##
            self.verbose(2,3,'%s: %s' % (patternFromCode(slim),p1.values()),1)
            self.verbose(2,3,'%s: %s vs %s\n' % (patternFromCode(slim),self.slimUP(slim),sum(p1.values())),2)

            ###~Calculate overall probability of observed support~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## All observed occurrences ##
            self.dict['Slim'][slim]['ExpUP'] = sum(p1.values())    # Expected number of observed UPCs
            (k,n,p) = (self.slimUP(slim), self.UPNum(), self.dict['Slim'][slim]['ExpUP']/self.UPNum()) # Use mean p1+
            if k <= 0: self.dict['Slim'][slim]['Prob'] = 1.0
            else: self.dict['Slim'][slim]['Prob'] = rje.binomial(k,n,p,usepoisson=False,callobj=self)
            ## Catch for binomial problems. Should no longer happen. ##
            if self.dict['Slim'][slim]['Prob'] <= 0:    # Shouldn't happen now! #
                self.errorLog('Probability for %s <= 0.0 due to numerical limitations: Given arbitrary 1e-16.!' % (patternFromCode(slim)),printerror=False)
                self.dict['Slim'][slim]['Prob'] = 1e-16
            ## Correction for restricted focal sequences ##
            self.focusAdjustment(slim)

            ###~Old Score calculations~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.dict['Slim'][slim]['S'] = self.slimUP(slim) * self.slimIC(slim)
            self.dict['Slim'][slim]['R'] = self.dict['Slim'][slim]['S'] * self.slimUP(slim) / self.dict['Slim'][slim]['ExpUP']
        except:
            self.errorLog('Error with slimProb(%s)' % slim)
            self.dict['Slim'][slim]['Prob'] = 1.0
#########################################################################################################################
    def aaDp1(self,slim,upc):  ### Setup  parameters for p1+ binomial using AADimerFreq
        '''Setup  parameters for p1+ binomial.'''
        ### ~ [1] Setup Parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        N = self.dict['AAFreq'][upc]['Total']   # Number of possible sites for SLiM to occur
        p = 0.0                                 # Probability of SLiM at each position
        k = 1                                   # Number of successful trials (occurrences)
        if self.getBool('SeqOcc'): k = max(1,self.slimOccNum(slim,upc))
        ### ~ [2] Special AADimerFreq option ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
        for slimvar in rje.listCombos(string.split(slim,'-')):
            v = rje.getFromDict(self.dict['AAFreq'][upc],slimvar[0],returnkey=False,default=0.0)
            slist = slimvar[0:]
            while slist:
                [i,x,j] = slist[:3]
                slist = slist[2:]
                try: v *= self.dict['AADimerFreq'][i][x][j]
                except: v = 0.0
            p += v
        return (k,N,p)
#########################################################################################################################
    def focusAdjustment(self,slim): ### Adjust raw probabilities according to focus dictionary
        '''
        Adjust raw probabilities according to focus dictionary.
        >> slim:str = SLiM for probability adjustment
        '''
        try:### ~ [1] Calculate probabilities for each focus group ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Focus']: return
            pgroup = {}
            for grp in self.dict['Focus']:
                a = len(self.dict['FocusUPC'][grp])     # No. of UPC in focal group
                b = self.slimUP(slim)                   # No. of UPC that the SLiM occurs in
                N = self.UPNum()                        # Total number of UPC
                pgroup[grp] = self.abNprob(a,b,N,1,'more')
            ### ~ [2] Adjust probability ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            meanp = rje.meanse(pgroup.values())[0]
            self.dict['Slim'][slim]['Prob'] *= rje.binomial(self.getInt('FocusOcc'),len(pgroup),meanp,callobj=self)
        except:
            self.errorLog('Major problem with %s.focusAdjustment()' % self.prog())
            raise
#########################################################################################################################
    def makeBonferroni(self):   ### Calculates Bonferroni & Lenferroni corrections number for dataset using DimFreq.
        '''Calculates Bonferroni corrections number for dataset.'''
        try:### ~ [1] Setup MustHave restriction adjustment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            alphx = musthave = len(self.list['Alphabet'])
            for aa in rje.sortUnique(self.list['MustHave']):
                while aa not in self.list['Alphabet'] and aa in self.list['MustHave']:
                    self.list['MustHave'].remove(aa)
                    self.warnLog('MustHave aa "%s" removed: not in alphabet' % aa)
            if self.list['MustHave']:
                self.list['MustHave'] = rje.sortUnique(self.list['MustHave'])
                musthave = len(self.list['MustHave'])
                self.printLog('#MUST','MustHave: %s' % string.join(self.list['MustHave']))
            ### ~ [2] Calculate motif space ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AlphaHelix'):
                self.dict['Extremf.'][2] = 1
                i = 3; iwild = 2
                maxi = pow(alphx,i)
                if musthave: maxi -= pow(alphx-musthave,i)   # 1+ position must have musthave aa
                maxi *= iwild
                #maxi = musthave * pow(len(self.list['Alphabet']),i-1) * iwild       # 1+ position must have musthave aa
                self.dict['Extremf.'][i] = maxi
                self.printLog('#SPACE','Motif Space, %d positions: %s motifs' % (i,rje.integerString(maxi)))
                i = 4; iwild = 1
                maxi = pow(alphx,i)
                if musthave: maxi -= pow(alphx-musthave,i)   # 1+ position must have musthave aa
                maxi *= iwild
                #maxi = musthave * pow(len(self.list['Alphabet']),i-1) * iwild       # 1+ position must have musthave aa
                self.dict['Extremf.'][i] = maxi
                self.printLog('#SPACE','Motif Space, %d positions: %s motifs' % (i,rje.integerString(maxi)))                
            elif self.getBool('DNA') and self.getBool('Palindrome'):
                for i in range(2,self.getInt('SlimLen')+1):
                    iwild = pow(self.getInt('MaxWild')-self.getInt('MinWild')+1,(i-1))
                    if rje.isOdd(i): maxi = (pow(alphx,((i+1)/2)) - pow(alphx-musthave,((i+1)/2))) * iwild
                    else: maxi = (pow(alphx,(i/2)) - pow(alphx-musthave,(i/2))) * iwild
                    #if rje.isOdd(i): maxi = musthave * pow(len(self.list['Alphabet']),((i+1)/2)-1) * iwild       # 1+ position must have musthave aa
                    #else: maxi = musthave * pow(len(self.list['Alphabet']),(i/2)-1) * iwild       # 1+ position must have musthave aa
                    self.dict['Extremf.'][i] = maxi
                    self.printLog('#SPACE','Motif Space, %d positions: %s palindrome motifs' % (i,rje.integerString(maxi)))
            elif self.getBool('FixLen'):
                iwild = (self.getInt('SlimLen') + self.getInt('MaxWild')) - 2
                for i in range(1,self.getInt('MaxWild')): iwild *= i
                maxi = pow(alphx,self.getInt('SlimLen'))
                if musthave: maxi -= pow(alphx-musthave,self.getInt('SlimLen')) # 1+ position must have musthave aa
                maxi *= iwild
                #maxi = musthave * pow(len(self.list['Alphabet']),self.getInt('SlimLen')) * iwild       # 1+ position must have musthave aa
                self.dict['Extremf.'][self.getInt('SlimLen')] = maxi
                self.printLog('#SPACE','Motif Space, %d positions: %s motifs' % (self.getInt('SlimLen'),rje.integerString(maxi)))
            else:
                for i in range(2,self.getInt('SlimLen')+1):
                    iwild = pow(self.getInt('MaxWild')-self.getInt('MinWild')+1,(i-1))
                    maxi = pow(alphx,i)
                    if musthave: maxi -= pow(alphx-musthave,i)   # 1+ position must have musthave aa
                    maxi *= iwild
                    #maxi = musthave * pow(len(self.list['Alphabet']),i-1) * iwild       # 1+ position must have musthave aa
                    self.dict['Extremf.'][i] = maxi
                    self.printLog('#SPACE','Motif Space, %d positions: %s motifs' % (i,rje.integerString(maxi)))
        except:
            self.errorLog('Problem with Motif Space calculation. Will not use.')
            self.setStr({'ProbScore':'Prob'})
            for i in range(1,self.getInt('SlimLen')+1): self.dict['Extremf.'][i] = 1.0
#########################################################################################################################
    def slimChanceNorman(self,fixes,replacesig=False):  ### Calculates special corrected significance of Norman Davey
        '''Calculates special corrected significance of Norman Davey.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            probscore = self.getStr('ProbScore')
            models = ['Error!','Sig','SigV','SigPrime','SigPrimeV']
            ## ~ [0a] Model Object options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            options = {}
            options["fixes"] = '%d' % fixes
            options["walltime"] = self.getNum('WallTime')
            options["maxLength"] = self.getInt('SlimLen')
            options["minLength"] = max(2,int(self.getNum('MinIC')))
            options["full"] = "F"
            options["precision"] = 100000
            options["wildcard"] = self.getInt('MaxWild')
            if fixes == 2: options["sample_percentage"] = {3:1,4:0.1,5:0.005,6:0.0005}
            ## ~ [0b] Score to be used temporarily ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            self.setNum({'NormSig': 0})    # Count of motifs <= probcut
            for slim in self.dict['Slim']:
                self.dict['Slim'][slim]["Normf."] = 1.0
            self.setStr({'ProbScore':'Sig'})
            ### ~ [1] Run Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            oldfreq = rje.combineDict({},self.dict['AAFreq'])
            rankerTrue = ned_rankbydistribution.rankByDistribution(options)     # Creates a rankByDistribution object
            if fixes == 2: self = rankerTrue.calculateRankingsMeanFix(self)     # SigV correction only
            else: self = rankerTrue.calculateRankings(self)                     # SigPrime correction
            for slim in self.dict['Slim']:
                self.dict['Slim'][slim][models[fixes]] = self.dict['Slim'][slim]['Normf.']
            if self.getInt('TopRanks') > 0: self.setNum({'NormSig':min(self.getInt('TopRanks'),self.getInt('NormSig'))})
            self.printLog('\r#PROB','Calculation of %s Probability complete: %s Sig SLiMs.' % (models[fixes],rje.iStr(self.getInt('NormSig'))))
            self.setStr({'ProbScore':probscore})
            ### ~ [2] Update Sig? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if replacesig:
                self.printLog('\r#SIG','Replacing basic "Sig" with "%s"' % models[fixes])
                self.list['SigSlim'] = []; self.dict['PCut'] = {}
                for slim in self.dict['Slim']:
                    self.dict['Slim'][slim]['Sig'] = self.dict['Slim'][slim]['Normf.']
                    self.sigSlim(slim,calculate=False)
            else: self.dict['AAFreq'] = oldfreq
        except: self.errorLog('Error in %s.slimChanceNorman(%s)' % (self.prog(),fixes),printerror=True)            
        self.setStr({'ProbScore':probscore})
#########################################################################################################################
    def slimChance(self):   ### SLiMChance Probability and Significance Calculations
        '''SLiMChance Probability and Significance Calculations.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['LogFactorial'] = {}
            ptxt = 'Calculating [%s] SLiM Probabilities (%s' % (self.chanceText(),self.getStr('ProbScore'))
            if not self.getNum('ProbCut'): ptxt = '%s ranking)' % ptxt
            elif self.getStr('ProbScore')[:3] in ['Sig','Pro']: ptxt = '%s<=%s)' % (ptxt,rje_slim.expectString(self.getNum('ProbCut')))
            else: ptxt = '%s>=%s)' % (ptxt,rje_slim.expectString(self.getNum('ProbCut')))
            sx = 0.0
            self.dict['PCut'] = {}
            ### ~ [1] Special Statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AllSig'):
                self.slimChanceNorman(2)
                self.slimChanceNorman(3)
                self.slimChanceNorman(4)
            elif self.getBool('SigV') and self.getBool('SigPrime'): self.slimChanceNorman(4,True)
            elif self.getBool('SigV'): self.slimChanceNorman(2,True)
            elif self.getBool('SigPrime'): self.slimChanceNorman(3,True)
            ### ~ [2] Original statistics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('AllSig') or (not self.getBool('SigV') and not self.getBool('SigPrime')):
                self.list['SigSlim'] = []; self.dict['PCut'] = {}
                for slim in self.dict['Slim']:
                    self.progLog('\r#PROB','%s %.2f%% (%s ranked SLiMs)' % (ptxt,sx/self.slimNum(),rje.integerString(len(self.list['SigSlim'])))); sx += 100.0
                    self.sigSlim(slim)
                    self.wallTime()
            self.printLog('\r#PROB','%s complete: %s ranked SLiMs.' % (ptxt,rje.integerString(len(self.list['SigSlim']))))
            ### ~ [3] Reduce to Significant SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.obj['SlimList'].list['Motif'] = []
            for slim in self.dict['Slim'].keys()[0:]:
                if slim not in self.list['SigSlim']: self.dict['Slim'].pop(slim)
            for slim in self.list['SigSlim']: self.obj['SlimList'].list['Motif'].append(self.addSLiMToList(slim))
        except SystemExit:
            self.list['SigSlim'] = []
            raise
        except: self.errorLog('Error in %s.slimChance()' % self.prog(),printerror=True,quitchoice=True)
#########################################################################################################################
    def chanceText(self):    ### Makes/Returns self.info['Chance']
        '''Makes/Returns self.info['Chance'].'''
        if self.getStrLC('Chance'): return self.getStr('Chance')
        self.setStr({'Chance':'Sig'})
        if self.getBool('SigPrime'): self.setStr({'Chance':self.getStr('Chance') + 'Prime'})
        if self.getBool('SigV'): self.setStr({'Chance':self.getStr('Chance') + 'V'})
        if self.getBool('AllSig'): self.setStr({'Chance':'All'})
        if self.prog()[:1] == 'Q': self.setStr({'Chance':'Q' + self.getStr('Chance')})
        if self.getBool('SeqOcc'): self.setStr({'Chance':'Occ' + self.getStr('Chance')})
        return self.getStr('Chance')
#########################################################################################################################
    def buildText(self):  ### Returns pickle identifier, also used for Outputs "Build" column (self.info['Build'])
        '''Returns pickle identifier, also used for Outputs "Build" column.'''
        ## Pickle name ##
        if not self.getStrLC('Build'):
            amb = 0
            if self.getBool('PreAmb') and self.list['Equiv']: amb += 1
            if self.getBool('PreAmb') and self.getBool('WildVar'): amb += 2
            #x#self.deBug('%s: %d' % (self.list['Equiv'],amb))
            if amb == 3 and self.getBool('CombAmb'): amb = 4
            if self.getBool('AlphaHelix'): self.setStr({'Build':'alpha-o%da%d' % (self.getInt('AmbOcc'),amb)})
            elif self.getBool('FixLen'): self.setStr({'Build':'fix-o%da%d' % (self.getInt('AmbOcc'),amb)})
            elif self.getInt('MinWild'): self.setStr({'Build':'l%dw%d-%do%da%d' % (self.getInt('SlimLen'),self.getInt('MinWild'),self.getInt('MaxWild'),self.getInt('AmbOcc'),amb)})
            else: self.setStr({'Build':'l%dw%do%da%d' % (self.getInt('SlimLen'),self.getInt('MaxWild'),self.getInt('AmbOcc'),amb)})
            if self.prog()[:1] == 'Q': self.setStr({'Build': 'Q' + self.getStr('Build')})
        return self.getStr('Build')
#########################################################################################################################
    def sigSlim(self,slim,calculate=True,rankfilter=True):  ### Adds SLiM to self.list['SigSlim'] providing it meets requirements.
        '''Adds SLiM to self.list['SigSlim'] providing it meets requirements.'''
        try:### ~ [1] Calculate SLiMChance Probability ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if calculate:
                self.slimProb(slim)
                prob = self.dict['Slim'][slim]['Prob']
                mlen = slimPos(slim)
                if mlen not in self.dict['PCut']: self.dict['PCut'][mlen] = 1.0
                if prob > self.dict['PCut'][mlen]: return False     # Cannot possibly be significant!
                ## ~ [1a] Binomial Bonferroni correct ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                e = self.dict['Extremf.'][mlen] * prob
                self.dict['Slim'][slim]['E'] = e
                self.dict['Slim'][slim]['Sig'] = rje.poisson(1,e,callobj=self)
                #x#if self.opt['AllSig']: self.dict['Slim'][slim]['Sigu'] = self.dict['Slim'][slim]['Sig']
            if not rankfilter: return
            #self.deBug(self.dict['Slim'][slim])
            ### ~ [2] Select ProbScore ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            try: prob = self.dict['Slim'][slim][self.getStr('ProbScore')]
            except:
                self.errorLog('ProbScore "%s" not recognised for %s. Will use SLiMChance Sig' % (self.getStr('ProbScore'),patternFromCode(slim)))
                self.setStr({'ProbScore':'Sig'})
                prob = self.dict['Slim'][slim][self.getStr('ProbScore')]
            if self.getNum('ProbCut') != 0:
                if self.getStr('ProbScore')[:3] in ['Sig','Pro'] and prob > self.getNum('ProbCut'): return False
                if self.getStr('ProbScore')[:3] not in ['Sig','Pro'] and prob < self.getNum('ProbCut'): return False
            if self.getStr('ProbScore')[:3] in ['Sig','Pro']: self.dict['Slim'][slim]['Ranker'] = prob
            else: self.dict['Slim'][slim]['Ranker'] = -prob
            #self.deBug(self.dict['Slim'][slim])
            ### ~ [3] Motif Occurrence Statistics & Filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.statFilter():
                if not self.occFilter(): Motif = self.motifOccStats(slim)
                else: Motif = self.addSLiMToList(slim)
                if not Motif:
                    self.errorLog('SLiM "%s" disappeared making Motif Object in sigSlim()' % slim,printerror=False)
                    return False
                Motif.stat['OccNum'] = self.slimOccNum(slim)
                Motif.stat['OccSeq'] = self.slimUP(slim)
                #!# What is this...?
                Motif.dict['Expect'] = {self.getStr('Basefile'):self.dict['Slim'][slim]['ExpUP']}
                try: self.obj['SlimList'].obj['SLiMCalc'].combMotifOccStats(Motif.occList())
                except: self.errorLog('Cannot combine OccStats!')
                if not rje_scoring.statFilterObj(self,[Motif],self.statFilter()):    
                    self.obj['SlimList'].removeMotif(Motif)
                    return False                    
            #self.deBug(slim)
            ### ~ [4] Score and Ranking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.list['SigSlim']:
                for i in range(len(self.list['SigSlim'])-1,-1,-1):  # Counts down through indices
                    sig = self.list['SigSlim'][i]
                    if self.dict['Slim'][slim]['Ranker'] > self.dict['Slim'][sig]['Ranker']:
                        self.list['SigSlim'].insert(i+1,slim)
                        break
                    elif self.dict['Slim'][slim]['Ranker'] == self.dict['Slim'][sig]['Ranker'] and self.dict['Slim'][slim]['Prob'] >= self.dict['Slim'][sig]['Prob']:
                        self.list['SigSlim'].insert(i+1,slim)
                        break
                    if i == 0: self.list['SigSlim'].insert(0,slim)
                if self.getInt('TopRanks') > 0: self.list['SigSlim'] = self.list['SigSlim'][:self.getInt('TopRanks')]
            else: self.list['SigSlim'] = [slim]
            #self.deBug(self.list['SigSlim'])
            ## ~ [4a] Update SlimList ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            if slim in self.list['SigSlim']: return True
            elif self.getStr('ProbScore') in ['Prob','Sig']: self.dict['PCut'][mlen] = self.dict['Slim'][slim]['Prob']
        except ValueError: raise
        except: self.errorLog('Error with sigSlim(%s)' % slim)
        return False
#########################################################################################################################
    ### <8> ### SLiMBuild Pickle Methods                                                                                #
#########################################################################################################################
    def OLDpickleMe(self,load=False):  ### Loads existing pickle, or saves pickle for later!
        '''Saves pickle for later!.'''
        try:
            ### Setup ###
            if load and self.force(): return None      # Re-run SLiMBuild!
            elif not load and self.getBool('MemSaver') and not self.extras(): return None  # Do not save pickle
            ## Pickle name ##
            mypickle = self.myPickle()

            ### Load ###
            if load:
                newme = None    ## New SLiMFinder object, loaded from Pickle
                ## Check for file and load ##
                for pfile in ['%s%s' % (self.getStr('ResDir'),self.dataset()),'%s%s' % (self.getStr('BuildPath'),self.dataset())]:
                    if not self.getBool('Win32') and os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)):
                        if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            if rje.isYounger('%s.%s.pickle.gz' % (pfile,mypickle),'%s.%s.pickle' % (pfile,mypickle)) == '%s.%s.pickle.gz' % (pfile,mypickle):
                                os.unlink('%s.%s.pickle' % (pfile,mypickle))
                        if not os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                            try: os.system('gunzip %s.%s.pickle.gz' % (pfile,mypickle))
                            except: self.errorLog('Cannot unzip %s.%s.pickle.gz' % (pfile,mypickle))
                    if os.path.exists('%s.%s.pickle' % (pfile,mypickle)):
                        self.printLog('#LOAD','Attempting to load %s pickle.' % self.prog(),log=False)
                        newme = pickle.load(open('%s.%s.pickle' % (pfile,mypickle),'r'))
                        self.printLog('#LOAD','%s intermediate loaded: %s.%s.pickle.' % (self.prog(),pfile,mypickle))
                        if not self.getBool('Win32'):
                            try:
                                if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                                os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                                self.printLog('#GZIP','%s %s.%s.pickle zipped.' % (self.prog(),pfile,mypickle))
                            except: self.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
                        break
                    elif self.dev(): self.printLog('#PICK','Cannot find %s.%s.pickle' % (pfile,mypickle))
                if not newme: return None
                ## Check other pertinent attributes - masking and additional filtering ##
                ## Note that MustHave and OccFilter filtering currently occur *after* SLiMBuild only ##
                changes = []
                for var in ['CompMask','CaseMask','MotifMask']:         # Info
                    if self.getStr(var) != newme.getStr(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['Masking','DisMask','ConsMask','MaskM']:   # Opt
                    if self.getBool(var) != newme.getBool(var): changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                for var in ['FTMask','IMask','Equiv']:      # List
                    slist = self.list[var][0:]
                    nlist = newme.list[var][0:]
                    slist.sort()
                    nlist.sort()
                    if slist != nlist:
                        changes.append(self.errorLog('Warning: "%s" parameter mismatch' % var, printerror=False, nextline=False))
                if newme.dict['Focus']:     # Post-Focus filtering is OK! Only worry if pickle has focus filtering
                    if rje.sortKeys(newme.dict['Focus']) != rje.sortKeys(self.dict['Focus']): # Assume all else is the same! #
                        changes.append(self.errorLog('Warning: "Focus" parameter mismatch', printerror=False, nextline=False))
                if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
                    if newme.list['QRegion'] and 'Focus' in newme.dict and 'Query' in newme.dict['Focus']:
                        if newme.list['QRegion'] != self.list['QRegion']: changes.append(self.errorLog('Warning: "QRegion" parameter mismatch', printerror=False, nextline=False))
                    else: changes.append(self.errorLog('Warning: "QRegion" parameter mismatch', printerror=False, nextline=False))
                ## Recreate or use pickle but add new commands ##
                #x#self.deBug(changes)
                if changes and (self.i() < 0 or rje.yesNo('%d SLiMBuild parameter mismatches with pickle. Create new pickle?' % len(changes))):
                    self.printLog('#PICKLE','Parameters changed. Making new pickle.')
                    return None
                self.list['Warning'] += changes
                newme.cmd_list = self.cmd_list
                newme.setStr(self.str)
                newme.setInt(self.int)
                newme.setNum(self.num)
                self.setBool({'Masked':newme.getBool('Masked'),'DNA':newme.getBool('DNA')})
                newme.setBool(self.bool)
                newme.setStr({'ResFile': self.getStr('ResFile'),'OccResFile': self.getStr('OccResFile'),'ResDir': self.getStr('ResDir'),'BuildPath': self.getStr('BuildPath')})
                newme.setNum({'StartTime':self.getNum('StartTime')})
                newme.obj['SlimList'] = self.obj['SlimList']    # Should take SLiMCalc with it
                newme.setLog(self.log)
                self.setLog(self.log)
                for mylist in ['Headers','MustHave','NewScore']: newme.list[mylist] = self.list[mylist]
                for mydict in ['Focus','NewScore']: newme.dict[mydict] = self.dict[mydict]
                newme.setupFocus()  #!# Need to convert to new Seq Objects - clean up at some point! #!#
                return newme

            ### Save ###
            if not self.getBool('Pickle'):
                self.printLog('#PICKLE','%s pickling disabled with pickle=F.' % self.prog(),log=True)
                return None
            self.printLog('#SAVE','Attempting to save %s with pickle.' % self.prog(),log=False)
            pickle.dump(self,open('%s.%s.pickle' % (self.getStr('Basefile'),mypickle),'w'))
            self.printLog('#SAVE','%s intermediate saved as %s.%s.pickle (Python pickle).' % (self.prog(),self.getStr('Basefile'),mypickle))
            if not self.getBool('Win32'):
                try:
                    pfile = self.getStr('Basefile')
                    if os.path.exists('%s.%s.pickle.gz' % (pfile,mypickle)): os.unlink('%s.%s.pickle.gz' % (pfile,mypickle))
                    os.system('gzip %s.%s.pickle' % (pfile,mypickle))
                    self.printLog('#GZIP','%s %s.%s.pickle zipped.' % (self.prog(),pfile,mypickle))
                except: self.errorLog('Cannot gzip %s.%s.pickle' % (pfile,mypickle))
            return None

        except:
            self.errorLog('Major problem with %s pickling!' % self.prog())
            return None
#########################################################################################################################
    def processPickle(self,newme):  ### Changes attributes accordingly
        '''Changes attributes accordingly. Replace this method in subclasses.'''
        try:### ~ [1] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## Check other pertinent attributes - masking and additional filtering ##
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
            if newme.dict['Focus']:     # Post-Focus filtering is OK! Only worry if pickle has focus filtering
                if rje.sortKeys(newme.dict['Focus']) != rje.sortKeys(self.dict['Focus']): # Assume all else is the same! #
                    changes.append(self.errorLog('Warning: "Focus" parameter mismatch', printerror=False, nextline=False))
            if self.list['QRegion'] and 'Focus' in self.dict and 'Query' in self.dict['Focus']:
                if newme.list['QRegion'] and 'Focus' in newme.dict and 'Query' in newme.dict['Focus']:
                    if newme.list['QRegion'] != self.list['QRegion']: changes.append(self.errorLog('Warning: "QRegion" parameter mismatch', printerror=False, nextline=False))
                else: changes.append(self.errorLog('Warning: "QRegion" parameter mismatch', printerror=False, nextline=False))
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
            for mylist in ['Headers','MustHave','NewScore']: newme.list[mylist] = self.list[mylist]
            for mydict in ['Focus','NewScore']: newme.dict[mydict] = self.dict[mydict]
            newme.setupFocus()  #!# Need to convert to new Seq Objects - clean up at some point! #!#
            return newme
        except: self.errorLog('Problem during %s.processPickle()' % self); return None
#########################################################################################################################
    ### <9> ### Results Output Methods                                                                                  #
#########################################################################################################################
    def getSlimProb(self,slim):   ### Returns appropriate SLiM Probability given settings
        '''Returns appropriate SLiM Score given settings.'''
        return self.dict['Slim'][slim][self.getStr('ProbScore')]
#########################################################################################################################
    def OLDmakeClouds(self):   ### Identifies "clouds" of motifs - similar patterns that overlap  (from SigSlim)
        '''Identifies "clouds" of motifs - similar patterns that overlap (from SigSlim).'''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Clouds'] = {}
            # 'Best': pattern for most significant cloud SLiM
            # 'Seq' : List of sequences in cloud
            # 'Slim': List of patterns for SLiMs in code
            # 'Code': List of SLiM codes to match self.dict['Slim'][slim]['Occ']
            # 'Sig' : Sig value for most significant cloud SLiM
            # 'UPC' : Total UPC coverage of Cloud
            # 'Fix' : Boolean value whether cloud includes a fixed (non-ambiguous) SLiM
            ctxt = 'Making "Motif Clouds" for %d Sig Motifs' % len(self.list['SigSlim'])
            if self.getInt('Clouds') < 0:
                self.printLog('#CLOUD','Not %s (clouds=X < 0)' % ctxt)
                ctxt = 'Preparing Cloud dictionary for SLiMMaker only'
            elif self.getInt('Clouds') < 1: self.setInt({'Clouds': self.getInt('MinOcc')})
            clouds = {}                             # Dictionary of slim:[slim cloud list] used in construction
            ctot = len(self.list['SigSlim']) * 3    # Attributes for Clouding progress

            ### ~ [1] ~ Rework Occ in searchable dictionary with every defined position ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lx = 0.0
            occ = {}        # New dictionary of {slim:{seq:[pos]}}
            for slim in self.list['SigSlim']:
                self.progLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot)); lx += 100.0
                occ[slim] = {}
                for (seq,pos) in self.dict['Slim'][slim]['Occ']:
                    if seq not in occ[slim]: occ[slim][seq] = []     # List of defined positions in SLiM-seq pair
                    occ[slim][seq] = rje.listUnion(occ[slim][seq],rje_slim.slimDefPos(slim,pos,seq.getSequence()))

            ### ~ [2] ~ Make cloud partnerships ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for slim1 in self.list['SigSlim']:
                clouds[slim1] = []  # List of slims that cloud with slim1
                self.progLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot)); lx += 100.0
                if self.getInt('Clouds') < 0: continue
                for slim2 in self.list['SigSlim']:
                    if slim1 == slim2: continue
                    cx = 0  # Shared sequence count
                    for seq in occ[slim1]:
                        if seq not in occ[slim2]: continue    # No occ for this seq
                        px = 0
                        for pos in occ[slim1][seq]:
                            if pos in occ[slim2][seq]: px += 1
                        if px >= 2: cx += 1
                        #X#print slim1, 'v', slim2, seq, px, cx, '>=', self.stat['Clouds']
                        if cx >= self.getInt('Clouds'):    # We have a match!
                            clouds[slim1].append(slim2)
                            break   # Check no more sequences
            cbug = ''
            for slim1 in self.list['SigSlim']:
                cbug += patternFromCode(slim1) + ': '
                for slim2 in clouds[slim1]: cbug += patternFromCode(slim2) + '; '
                cbug += '\n'
            #self.deBug(cbug)

            ###~Make Actual Clouds~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            cnum = 0
            for slim in self.list['SigSlim'][0:]:
                ## Make a cloud for each Sig SLiM not already in one ##
                self.progLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot)); lx += 100.0
                if slim not in clouds: continue     # Absorbed into another cloud
                ## Make a new cloud with this SLiM as the "best" variant ## 
                cnum += 1   # New cloud
                self.dict['Clouds'][cnum] = {'Best':patternFromCode(slim),'Seq':[],'Slim':[patternFromCode(slim)],
                                             'Sig':rje_slim.expectString(self.dict['Slim'][slim]['Sig']),
                                             'Code':[slim],'Fix':rje_slim.slimFix(slim)==rje_slim.slimPos(slim) and not rje_slim.varWild(slim)}
                for seq in occ[slim]:
                    if occ[slim][seq]: self.dict['Clouds'][cnum]['Seq'].append(seq)
                self.dict['Slim'][slim]['Cloud'] = cnum
                cx = 0
                while cx != len(clouds[slim]):      # Keep absorbing and reducing clouds till no more match
                    cx = len(clouds[slim])
                    for slim2 in clouds[slim][0:]:
                        if slim2 == slim: continue
                        if slim2 not in clouds: continue    # Already taken
                        self.dict['Slim'][slim2]['Cloud'] = cnum
                        self.dict['Clouds'][cnum]['Slim'].append(patternFromCode(slim2))
                        self.dict['Clouds'][cnum]['Code'].append(slim2)
                        for seq in occ[slim2]:
                            if occ[slim2][seq] and seq not in self.dict['Clouds'][cnum]['Seq']: self.dict['Clouds'][cnum]['Seq'].append(seq)
                        for newslim in clouds.pop(slim2):   # Combine clouds
                            if newslim not in clouds[slim]: clouds[slim].append(newslim)
                            #X# self.dict['Slim'][newslim]['Cloud'] = cnum
                        self.dict['Clouds'][cnum]['Fix'] = self.dict['Clouds'][cnum]['Fix'] or (rje_slim.slimFix(slim2) == rje_slim.slimPos(slim2) and not rje_slim.varWild(slim2))
                ## Convert Cloud seqs into UPC ##
                self.dict['Clouds'][cnum]['UPC'] = []
                for seq in self.dict['Clouds'][cnum]['Seq']:
                    upc = self.getUP(seq)
                    if upc not in self.dict['Clouds'][cnum]['UPC']: self.dict['Clouds'][cnum]['UPC'].append(upc)
            self.printLog('\r#CLOUD','%s: %d Clouds' % (ctxt,cnum))

            ### ~ CloudFix reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('CloudFix'):
                srx = 0; crx = 0    # Number of SLiMs/Clouds removed
                for cnum in range(1,len(self.dict['Clouds'])+1):
                    if self.dict['Clouds'][cnum]['Fix']:
                        if crx:
                            if (cnum - crx) in self.dict['Clouds']: raise ValueError
                            self.dict['Clouds'][cnum-crx] = self.dict['Clouds'].pop(cnum)
                            for slim in self.dict['Clouds'][cnum-crx]['Code']: self.dict['Slim'][slim]['Cloud'] = cnum-crx
                        continue
                    crx += 1
                    for slim in self.dict['Clouds'].pop(cnum)['Code']:
                        self.dict['Slim'].pop(slim); srx += 1
                        self.list['SigSlim'].remove(slim)
                        pattern = rje_slim.patternFromCode(slim)
                        for Motif in self.obj['SlimList'].slims()[0:]:
                            if pattern == Motif.info['Sequence']: self.obj['SlimList'].removeMotif(Motif)
                self.printLog('#FIX','CloudFix filter: %d SLiM (%d cloud) removed.' % (srx,crx))

            ### Finish ###
            #x#self.deBug(self.dict['Clouds'])
        except: self.errorLog('Problem %s' % ctxt)
#########################################################################################################################
    def makeClouds(self):   ### Identifies "clouds" of motifs - similar patterns that overlap  (from SigSlim)
        '''
        Identifies "clouds" of motifs - similar patterns that overlap (from SigSlim) - and makes self.dict['Clouds'].
        # 'Best': pattern for most significant cloud SLiM
        # 'Seq' : List of sequences in cloud
        # 'Slim': List of patterns for SLiMs in code
        # 'Code': List of SLiM codes to match self.dict['Slim'][slim]['Occ']
        # 'Sig' : Sig value for most significant cloud SLiM
        # 'UPC' : Total UPC coverage of Cloud
        # 'Fix' : Boolean value whether cloud includes a fixed (non-ambiguous) SLiM
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            self.dict['Clouds'] = {}
            ctxt = 'Making "Motif Clouds" for %d Sig Motifs' % len(self.list['SigSlim'])
            if self.getInt('Clouds') < 0:
                self.printLog('#CLOUD','Not %s (clouds=X < 0)' % ctxt)
                ctxt = 'Preparing Cloud dictionary for SLiMMaker only'
            elif self.getInt('Clouds') < 1: self.setInt({'Clouds': self.getInt('MinOcc')})
            clouds = {}                             # Dictionary of slim:[slim cloud list] used in construction
            ctot = len(self.list['SigSlim']) * 2    # Attributes for Clouding progress

            ### ~ [1] ~ Make cloud partnerships ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            lx = 0.0
            occ = {}        # New dictionary of {slim:{seq:[pos]}}
            for slim1 in self.list['SigSlim']:
                self.progLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot)); lx += 100.0
                clouds[slim1] = []  # List of slims that cloud with slim1
                if self.getInt('Clouds') < 0: continue      # Single motif in cloud (needed for SLiMMaker)
                for slim2 in self.list['SigSlim']:
                    if slim1 == slim2: continue
                    sharedseq = []      # Shared sequences
                    for (seq1,pos1) in self.dict['Slim'][slim1]['Occ']:
                        if seq1 in sharedseq: continue
                        for (seq2,pos2) in self.dict['Slim'][slim2]['Occ']:
                            if seq1 != seq2 or seq2 in sharedseq: continue
                            if len(rje.listIntersect(rje_slim.slimDefPos(slim1,pos1,seq1.getSequence()),rje_slim.slimDefPos(slim2,pos2,seq2.getSequence()))) >= 2: sharedseq.append(seq1)
                    if len(sharedseq) >= self.getInt('Clouds'): clouds[slim1].append(slim2)

            ### ~ [2] ~ Make Actual Clouds ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cnum = 0
            for slim in self.list['SigSlim'][0:]:
                ## Make a cloud for each Sig SLiM not already in one ##
                self.progLog('\r#CLOUD','%s: %.1f%%' % (ctxt,lx/ctot)); lx += 100.0
                if slim not in clouds: continue     # Absorbed into another cloud
                ## Make a new cloud with this SLiM as the "best" variant ##
                cnum += 1   # New cloud
                self.dict['Clouds'][cnum] = {'Best':patternFromCode(slim),'Seq':[],'Slim':[patternFromCode(slim)],
                                             'Sig':rje_slim.expectString(self.dict['Slim'][slim]['Sig']),
                                             'Code':[slim],'Fix':rje_slim.slimFix(slim)==rje_slim.slimPos(slim) and not rje_slim.varWild(slim)}
                for (seq,pos) in self.dict['Slim'][slim]['Occ']:
                    if seq not in self.dict['Clouds'][cnum]['Seq']: self.dict['Clouds'][cnum]['Seq'].append(seq)
                self.dict['Slim'][slim]['Cloud'] = cnum
                cx = 0
                while cx != len(clouds[slim]):      # Keep absorbing and reducing clouds till no more match
                    cx = len(clouds[slim])
                    for slim2 in clouds[slim][0:]:
                        if slim2 == slim: continue
                        if slim2 not in clouds: continue    # Already taken
                        self.dict['Slim'][slim2]['Cloud'] = cnum
                        self.dict['Clouds'][cnum]['Slim'].append(patternFromCode(slim2))
                        self.dict['Clouds'][cnum]['Code'].append(slim2)
                        for (seq2,pos2) in self.dict['Slim'][slim2]['Occ']:
                            if seq2 not in self.dict['Clouds'][cnum]['Seq']: self.dict['Clouds'][cnum]['Seq'].append(seq2)
                        for newslim in clouds.pop(slim2):   # Combine clouds
                            if newslim not in clouds[slim]: clouds[slim].append(newslim)
                            #X# self.dict['Slim'][newslim]['Cloud'] = cnum
                        self.dict['Clouds'][cnum]['Fix'] = self.dict['Clouds'][cnum]['Fix'] or (rje_slim.slimFix(slim2) == rje_slim.slimPos(slim2) and not rje_slim.varWild(slim2))
                ## Convert Cloud seqs into UPC ##
                self.dict['Clouds'][cnum]['UPC'] = []
                for seq in self.dict['Clouds'][cnum]['Seq']:
                    upc = self.getUP(seq)
                    if upc not in self.dict['Clouds'][cnum]['UPC']: self.dict['Clouds'][cnum]['UPC'].append(upc)
            self.printLog('\r#CLOUD','%s: %d Clouds' % (ctxt,cnum))

            ### ~ [3] CloudFix reduction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getBool('CloudFix'):
                srx = 0; crx = 0    # Number of SLiMs/Clouds removed
                for cnum in range(1,len(self.dict['Clouds'])+1):
                    if self.dict['Clouds'][cnum]['Fix']:
                        if crx:
                            if (cnum - crx) in self.dict['Clouds']: raise ValueError
                            self.dict['Clouds'][cnum-crx] = self.dict['Clouds'].pop(cnum)
                            for slim in self.dict['Clouds'][cnum-crx]['Code']: self.dict['Slim'][slim]['Cloud'] = cnum-crx
                        continue
                    crx += 1
                    for slim in self.dict['Clouds'].pop(cnum)['Code']:
                        self.dict['Slim'].pop(slim); srx += 1
                        self.list['SigSlim'].remove(slim)
                        pattern = rje_slim.patternFromCode(slim)
                        for Motif in self.obj['SlimList'].slims()[0:]:
                            if pattern == Motif.info['Sequence']: self.obj['SlimList'].removeMotif(Motif)
                self.printLog('#FIX','CloudFix filter: %d SLiM (%d cloud) removed.' % (srx,crx))

            ### Finish ###
            #self.deBug(self.dict['Clouds'])
        except: self.errorLog('Problem %s' % ctxt)
#########################################################################################################################
    def cloudConsensi(self):    ### Generate consensus motifs from clouds using SLiMMaker.
        '''Generate consensus motifs from clouds using SLiMMaker.'''
        try:### ~ [0] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            cdict = self.dict['Clouds'] # {cnum:{'Code': List of SLiM codes to match self.dict['Slim'][slim]['Occ'] = (seq,occ)}}

            ### ~ [1] Sequence positions & Offsets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for c in cdict:
                offset = {}     # {slim1:{slim2:[offset list]}}
                for i in range(len(cdict[c]['Code'])):
                    slimi = cdict[c]['Code'][i]
                    if slimi not in offset: offset[slimi] = {}
                    for slimj in cdict[c]['Code'][i+1:]:
                        if slimj not in offset: offset[slimj] = {}
                        if slimj not in offset[slimi]: offset[slimi][slimj] = []
                        if slimi not in offset[slimj]: offset[slimj][slimi] = []
                        leni = rje_slim.slimLenFromCode(slimi)
                        lenj = rje_slim.slimLenFromCode(slimj)
                        for (seqi,posi) in self.dict['Slim'][slimi]['Occ']:
                            for (seqj,posj) in self.dict['Slim'][slimj]['Occ']:
                                if seqi != seqj: continue
                                #if posj > posi + leni - 3: continue   # Need 2 position overlap
                                #if posi > posj + lenj - 3: continue   # Need 2 position overlap
                                #if len(rje.listIntersect(rje_slim.slimDefPos(slimi,posi,seqi.getSequence(),self),rje_slim.slimDefPos(slimj,posj,seqj.getSequence(),self))) < 2: continue
                                if len(rje.listIntersect(rje_slim.slimDefPos(slimi,posi,seqi.getSequence()),rje_slim.slimDefPos(slimj,posj,seqj.getSequence()))) < 2: continue
                                offset[slimi][slimj].append(posj-posi)
                                offset[slimj][slimi].append(posi-posj)
                        if not offset[slimi][slimj]:
                            offset[slimj].pop(slimi)
                            offset[slimi].pop(slimj)
                            continue
                        ## Restrict offset to most frequent offset, else smallest offset
                        if len(rje.sortUnique(offset[slimi][slimj])) > 1: self.warnLog('Multiple motif offsets for %s vs %s' % (rje_slim.patternFromCode(slimi),rje_slim.patternFromCode(slimj)))
                        olist = rje.listMax(offset[slimi][slimj])
                        while len(olist) > 1:
                            if rje.modulus(olist[0]) > rje.modulus(olist[1]): olist.pop(1)
                            else: olist.pop(0)
                        offset[slimi][slimj] = olist[0]
                        olist = rje.listMax(offset[slimj][slimi])
                        while len(olist) > 1:
                            if rje.modulus(olist[0]) > rje.modulus(olist[1]): olist.pop(1)
                            else: olist.pop(0)
                        offset[slimj][slimi] = olist[0]

            ### ~ [2] Convert offsets to first SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                slim1 = cdict[c]['Code'][0]
                offset[slim1][slim1] = 0
                while len(offset) > 1:  # Want to convert to offset[slim1] dictionary
                    #  offset[slimi][slimj] = posj - posi
                    for slimi in cdict[c]['Code'][1:]:
                        if slimi not in offset: continue
                        for slimj in rje.sortKeys(offset[slimi])[0:]:
                            if slimj in offset[slim1] and slimi in offset[slim1]: pass
                            # Either slimi or slimj is missing from slim1
                            elif slimi in offset[slim1]:
                                offset[slim1][slimj] = offset[slim1][slimi] + offset[slimi][slimj]
                            else:
                                offset[slim1][slimi] = offset[slim1][slimj] + offset[slimj][slimi]
                            offset[slimi].pop(slimj)
                        if not offset[slimi]: offset.pop(slimi)

            ### ~ [3] Assemble peptide lists ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
                offset = offset[slim1]
                cdict[c]['PepSeq'] = []
                ## ~ [3a] ~ Adjust offsets to be positive ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                offmin = min(offset.values())
                for slimi in offset: offset[slimi] -= offmin
                ## ~ [3b] ~ Peptide length ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                peplen = 0
                for slimi in offset: peplen = max(peplen,offset[slimi] + rje_slim.slimLenFromCode(slimi))
                ## ~ [3c] ~ Generate peptides ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                seqpos = []     # Combined (seq,offset pos)
                for slim in offset:
                    for (seq,pos) in self.dict['Slim'][slim]['Occ']:
                        pos -= offset[slim]
                        if (seq,pos) not in seqpos:
                            seqpos.append((seq,pos))
                            if pos >=0 : cdict[c]['PepSeq'].append(seq.info['Sequence'][pos:pos+peplen])     #!# Check masking!
                            else: cdict[c]['PepSeq'].append('-' * -pos + seq.info['Sequence'][:pos+peplen])
                            if not cdict[c]['PepSeq'][-1]: self.bugPrint('%s+%s' % (pos,peplen)); self.debug(seq.info)
                scmd = self.cmd_list+['peptides=%s' % string.join(cdict[c]['PepSeq'],',')]
                (cdict[c]['SLiMMaker'],cdict[c]['SLiMMakerStr']) = slimmaker.SLiMMaker(self.log,scmd).run(log=False)
                self.printLog('#SLIM','Cloud %d (%s) SLiMMaker: %s' % (c,string.join(cdict[c]['Slim'],'|'),cdict[c]['SLiMMaker']))
                self.printLog('#FREQ','SLiMMaker %s' % (cdict[c]['SLiMMakerStr']))
        except: self.errorLog('Error in %s.cloudConsensi()' % self.prog())
#########################################################################################################################       
    def rankScore(self):  ### Scores and Ranks Sig Motifs
        '''Scores and Ranks Sig Motifs.'''
        try:
            ### Assign numerical rankings ###
            if not self.list['SigSlim']: return
            (prev,rank,prevp) = (1,1,1)
            self.printLog('\r#RANK','Rank calculations...',newline=False,log=False)
            for slim in self.list['SigSlim']:
                ## Assign Rank ##
                if self.getBool('SlimDisc') or self.dict['Slim'][slim][self.getStr('ProbScore')] > prev or self.dict['Slim'][slim]['Prob'] != prevp:
                    self.dict['Slim'][slim]['Rank'] = self.list['SigSlim'].index(slim) + 1
                else: self.dict['Slim'][slim]['Rank'] = rank
                (prev,rank,prevp) = (self.dict['Slim'][slim][self.getStr('ProbScore')],self.dict['Slim'][slim]['Rank'],self.dict['Slim'][slim]['Prob'])
            self.printLog('\r#RANK','Rank calculations complete')
        except:
            self.errorLog('Major disaster during SLiMFinder.rankScore()')
            raise
#########################################################################################################################
    def setupResults(self):     ### Sets up Main Results File as well as StatFilters etc.
        '''Sets up Main Results File as well as StatFilters etc.'''
        try:### ~ [1] Setup Initial Headers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [1a] Main Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            slist = self.obj['SlimList']
            self.setupBasefile()
            self.list['Headers'] = basic_headers[0:]
            if self.getBool('AllSig'):
                si = self.list['Headers'].index('Pattern')
                self.list['Headers'] = self.list['Headers'][:si] + ['SigV','SigPrime','SigPrimeV'] + self.list['Headers'][si:]
            if not self.getBool('OldScores'):
                self.list['Headers'].remove('S')
                self.list['Headers'].remove('R')
            ## ~ [1b] Occurrence Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            occhead = self.list['OccHeaders'] = ['Dataset','RunID','Rank','Pattern','Sig','Seq','Start_Pos','End_Pos','Prot_Len','Match','Variant','MisMatch','Desc']
            #!# Consider adding Basefile too but Dataset and RunID should be enough #!#

            ### ~ [2] Special Stats ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for h in slist.obj['SLiMCalc'].list['Headers']:
                if h not in self.list['Headers']: self.list['Headers'].append(h)
            for h in slist.obj['SLiMCalc'].list['OccHeaders']:
                if h not in occhead: occhead.append(h)
            occhead += ['PepSeq','PepDesign']

            ### ~ [3] Custom Scores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Modify slimcalc to deal with this - both occurrence and combined new scores? #!#
            #!# >> Make a method that takes the existing headers as input along with NewScore cmd? #!#
            #!# For now, Custom scores have been retired as redundant (original a SLiMPickings heuristic feature #!#
            #!#(newheads,self.list['NewScore'],self.dict['NewScore']) = rje_scoring.setupCustomScores(self,self.list['Headers'],self.list['NewScore'],self.dict['NewScore'])
            #!#for new in newheads:
            #!#    if new not in self.list['Headers'] and new not in self.list['OccStats']:
            #!#        if rje.formula(self,self.dict['NewScore'][new],varlist=self.list['Headers'],check=False,calculate=False): self.list['Headers'].append(new)
            #X#self.deBug(self.list['OccStats'])   # This is used to determine which OccStats to combine (without NewScore)
            #X#self.deBug(self.list['Headers'])    # These are the headers for the main output - including combinations
            #X#self.deBug(self.list['NewScore'])   # These are new scores to calculate - compare with OccStats/Headers for when calculating

            ### ~ [4] ResFile(s) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.getStrLC('ResFile') and self.getStr('ResFile').find('.') < 0: self.setStr({'ResFile': self.getStr('ResFile') + '.csv'})
            resfilesplit = string.split(self.getStr('ResFile'),'.')
            resfilesplit = resfilesplit[:-1] + ['occ',resfilesplit[-1]]
            self.setStr({'OccResFile':string.join(resfilesplit,'.')})
            #self.deBug(self.list['Headers'])
        except:
            self.errorLog('Problem with SLiMFinder.setupResults()')
            raise
#########################################################################################################################
    def backupOrCreateResFile(self):    ### Backups up and/or creates main results file
        '''Backups up and/or creates main results files.'''
        if self.getStrLC('ResFile') and self.getBool('SlimChance'):
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.getStr('ResFile'),write=True))
            rje.delimitedFileOutput(self,self.getStr('ResFile'),self.resHead(),delimit,rje_backup=True)
            rje.delimitedFileOutput(self,self.getStr('OccResFile'),self.list['OccHeaders'],delimit,rje_backup=True)
        self.dict['Output']['main'] = 'ResFile'
        self.dict['Output']['occ'] = 'OccResFile'
#########################################################################################################################
    def resHead(self):  ### Returns main Output headers
        '''Returns main Output headers.'''
        heads = ['Dataset','RunID','Masking','Build','Chance','RunTime','SeqNum','UPNum','AANum','MotNum'] + self.list['Headers']
        if self.getBool('Test'): heads.insert(heads.index('Sig'),'E')
        if self.getBool('AllSig'): heads.remove('Chance')
        return heads
#########################################################################################################################
    def results(self,aborted=None):  ### Main SLiMFinder Results Output
        '''Main SLiMFinder Results Output.'''
        try:
            ###~Setup General Dataset Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            self.printLog('#OUT','Generating outputs...',log=False)
            self.setupBasefile()    # Just in case this is called from outside main SLiMFinder run
            extrabase = self.runBase()
            ## General Output Stats ##
            sec = int(time.time() - self.getNum('StartTime') + 0.5)
            (hour,min) = (0,0)
            while sec >= 60: (min,sec) = (min+1,sec-60)
            while min >= 60: (hour,min) = (hour+1,min-60)
            self.setStr({'RunTime': '%s:%s:%s' % (rje.preZero(hour,24),rje.preZero(min,60),rje.preZero(sec,60))})
            totalaa = 0     # Total number of AA in dataset
            for seq in self.seqs():
                if self.getBool('Masked'): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                totalaa += seq.nonX()
            t = time.localtime(self.log.stat['StartTime'])
            self.setStr({'Date': '%s-%s-%s' % (str(t[0])[-4:],rje.preZero(t[1],12),rje.preZero(t[2],31))})
            if not self.getStr('RunID'): self.setStr({'RunID': self.getStr('Date')})
            masking = self.maskText()   # Summary of Masking Options

            ###~Setup SLiMDisc Output (and make dat.rank)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            if (self.getBool('SlimDisc') or self.extras(3)) and not aborted:
                ## Rank File ##
                rankfile = extrabase + '.rank'; self.dict['Output']['rank'] = rankfile
                open(rankfile,'w')
                headers = basic_headers[0:]
                prehead = ['#---------------------------Input stats----------------------------',
                           '#\tInput dataset:   %s' % os.path.abspath(self.getStr('Input')),
                           '#\tNo. of proteins: %d' % self.seqNum(),
                           '#\tOverall MST:     %f' % self.getNum('MST'),
                           '#']
                for seq in self.seqs(): prehead.append('#\t-%s\t%d a.a' % (seq.shortName(),seq.aaLen()))
                prehead += ['#                               -----', '#\tNo. of residues:\t%d' % totalaa]
                ## UPC ##
                prehead += ['#','#\tNo. of UP Clusters:\t%d' %  self.UPNum()]
                upout = []
                for seq in self.seqs():
                    upc = self.getUP(seq)
                    if upc in upout: continue
                    upout.append(upc)
                    upstr = []
                    for useq in upc: upstr.append(useq.shortName())
                    prehead.append('#\t-[%s]' % (string.join(upstr,',')))
                ## No. Motifs ##
                if self.getStr('ProbScore')[:3] in ['Sig','Pro']: prehead += ['#','#\tNo. of motifs %s <= %s: %d' % (self.getStr('ProbScore'),self.getNum('ProbCut'),len(self.list['SigSlim']))]
                else: prehead += ['#','#\tNo. of motifs %s >= %s: %d' % (self.getStr('ProbScore'),self.getNum('ProbCut'),len(self.list['SigSlim']))]
                prehead += ['#------------------------------------------------------------------']
            
                ## Dat Rank File ##
                datfile = extrabase + '.dat.rank'; self.dict['Output']['dat.rank'] = rankfile
                dathead = ['#Proteins']
                for i in range(self.seqNum()): dathead.append('%d\t%s' % (i,self.seqs()[i].shortName()))
                dathead.append('#Motifs')
                DAT = open(datfile,'w')
                DAT.write('%s\n' % string.join(dathead,'\n'))
                for slim in self.list['SigSlim']:
                    DAT.write('%d\t%s\t' % (self.dict['Slim'][slim]['Rank'],patternFromCode(slim)))
                    for (seq,pos) in self.dict['Slim'][slim]['Occ']: DAT.write('%d:%d ' % (self.seqs().index(seq),pos))
                    DAT.write('\n')
                DAT.close()
            
            ###~Main SLiMFinder Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            ## Setup ##
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.getStr('ResFile'),write=True))
            reshead = self.resHead()    #!# Combine with general setup of masking & totalaa etc. setup #!#
            #self.deBug(self.list['Headers'])
            #self.deBug(reshead)
            ## Main Output and SLiMDisc Rank File ##
            self.setBool({'Append': False})
            if (self.getBool('SlimDisc') or self.extras(3)) and not aborted:
                rje.delimitedFileOutput(self,rankfile,headers,'\t',rje_backup=False)
                rlines = open(rankfile,'r').read()
                open(rankfile,'w').write('%s\n%s' % (string.join(prehead,'\n'),rlines))
            ## SLiM Data ##
            for slim in self.list['SigSlim']:   #!# Add SigNum? #!#
                pattern = patternFromCode(slim)
                Motif = self.addSLiMToList(slim)    #x# self.obj['MotifList'].mapPattern(pattern,update=False)
                datadict = {'Dataset':self.dataset(),'RunID':self.getStr('RunID'),
                            'Masking':masking,'Build':self.getStr('Build'),'RunTime':self.getStr('RunTime'),'Chance':self.chanceText(),
                            'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa,
                            'MotNum':len(self.dict['Slim']),'Rank':self.dict['Slim'][slim]['Rank'],
                            'Pattern':pattern,'Occ':len(self.dict['Slim'][slim]['Occ']),'Support':self.slimSeqNum(slim),
                            'IC':'%.2f' % self.slimIC(slim),'UP':self.slimUP(slim), 'Norm':self.slimUP(slim),
                            }
                if self.dict['Clouds']:
                    datadict['Cloud'] = self.dict['Slim'][slim]['Cloud']
                    datadict['CloudSeq'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['Seq'])
                    datadict['CloudUP'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['UPC'])
                for p in ['ExpUP','Prob','Sig','E','R','S','SigV','SigPrime','SigPrimeV']:
                    if self.dict['Slim'][slim].has_key(p): datadict[p] = rje_slim.expectString(self.dict['Slim'][slim][p])
                for h in self.list['Headers']:
                    if Motif and h not in datadict and h in Motif.stat:
                        if h.endswith('_mean'): datadict[h] = rje_slim.expectString(Motif.stat[h])
                        else: datadict[h] = Motif.stat[h]
                #X#print Motif.stat
                ## Main Output ##
                if self.getStrLC('ResFile') and self.getBool('SlimChance'):
                    rje.delimitedFileOutput(self,self.getStr('ResFile'),reshead,delimit,datadict)
                ## Rank File ##
                if (self.getBool('SlimDisc') or self.extras(3)) and not aborted:
                    datadict['Rank'] = '(%d)' % self.dict['Slim'][slim]['Rank']
                    rje.delimitedFileOutput(self,rankfile,headers,'\t',datadict)
            ## Clouds ##
            if not aborted: self.cloudOutput()
            ## No SLiMs? ##
            if not self.list['SigSlim'] and self.getStrLC('ResFile') and self.getBool('SlimChance'):
                datadict = {'Dataset':self.dataset(),'RunID':self.getStr('RunID'),'Sig':self.getNum('ProbCut'),
                            'Masking':masking,'Build':self.getStr('Build'),'RunTime':self.getStr('RunTime'),'Chance':self.chanceText(),
                            'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa,
                            'MotNum':len(self.dict['Slim']),'Rank':0,'Pattern':'-'}
                if aborted: datadict['Pattern'] = aborted
                #self.deBug(datadict)
                rje.delimitedFileOutput(self,self.getStr('ResFile'),reshead,delimit,datadict)
                    
            ###~Additional Screen Output~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
            if (self.getBool('SlimDisc') or self.extras(3)) and not aborted:
                ox = self.seqNum() + 11
                self.verbose(1,1,string.join(['\n\n']+open(datfile,'r').readlines()[:ox+2],''),2)
                self.verbose(1,1,string.join(open(rankfile,'r').readlines()[:ox+16+self.UPNum()],''),2)
            if self.list['Warning']: self.verbose(0,0,string.join(['']+self.list['Warning'],'\n'),1)
            return True            
        except: self.errorLog('Error with SLiMFinder results output.')
        return False
#########################################################################################################################
    def getSigSlim(self,pattern):   ### Returns slimcode for given pattern (should be in SigSlim)
        '''Returns slimcode for given pattern (should be in SigSlim).'''
        for slim in self.list['SigSlim']:
            if patternFromCode(slim) == pattern: return slim
        return None
#########################################################################################################################
    def cloudOutput(self):  ### Generates motif cloud output
        '''Generates motif cloud output.'''
        try:
            ### Motif Cloud Summary Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.dict['Clouds'] or not self.extras(): return
            extrabase = self.runBase()
            cloudfile = '%s.cloud.txt' % extrabase; self.dict['Output']['cloud'] = cloudfile
            ## General Cloud Data ##
            CLOUD = open(cloudfile,'w')
            CLOUD.write('%s: %d Seq, %d UPC, %d Sig SLiMs in %d Clouds\n\n' % (self.dataset(),self.seqNum(),self.UPNum(),len(self.list['SigSlim']),len(self.dict['Clouds'])))
            for c in rje.sortKeys(self.dict['Clouds']):
                CLOUD.write('Cloud %d (%d SLiMs):\t%s' % (c,len(self.dict['Clouds'][c]['Slim']),string.join(self.dict['Clouds'][c]['Slim'],'\t')))
                CLOUD.write('\t(p = %s)\n' % self.dict['Clouds'][c]['Sig'])
                seqnames = []
                for seq in self.dict['Clouds'][c]['Seq']: seqnames.append(seq.shortName())
                seqnames.sort()
                CLOUD.write('%d UPC\t%d Seq:\t%s\n\n' % (len(self.dict['Clouds'][c]['UPC']),len(self.dict['Clouds'][c]['Seq']),string.join(seqnames,'\t')))
                CLOUD.write('SLiMMaker Consensus: %s\n' % self.dict['Clouds'][c]['SLiMMaker'])
                CLOUD.write('- %s\n\n' % self.dict['Clouds'][c]['SLiMMakerStr'])
                CLOUD.write('%s\n\n' % string.join(self.dict['Clouds'][c]['PepSeq'],'\n'))
            CLOUD.close()
            ## Cloud Coverage ##
            for type in ['Seq','UPC']:
                if type == 'Seq': dnum = self.seqNum()
                else: dnum = self.UPNum()
                open(cloudfile,'a').write('\n#Cloud %s coverage/overlap\n' % type)
                chead = ['Cloud','Dataset'] + rje.sortKeys(self.dict['Clouds'])
                rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict={})
                cloudcoverage = {}
                for c in rje.sortKeys(self.dict['Clouds']):
                    cloudcoverage[c] = {'Cloud':'%d:%s' % (c,self.dict['Clouds'][c]['Best']),'Dataset':'%.3f' % (len(self.dict['Clouds'][c][type]) / float(dnum))}
                    for j in rje.sortKeys(self.dict['Clouds']):
                        common = [] 
                        for seq in self.dict['Clouds'][j][type]:
                            if seq in self.dict['Clouds'][c][type]: common.append(seq)
                        cloudcoverage[c][j] = '%.3f' % (len(common) / float(len(self.dict['Clouds'][j][type])))
                    rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict=cloudcoverage[c])
                ## Cloud coverage probabilites - overlaps ##
                if len(self.dict['Clouds']) > 1:
                    open(cloudfile,'a').write('\n#Probabilities of motif clouds overlapping %s this much or more\n' % type)
                    rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict={})
                    for c in rje.sortKeys(self.dict['Clouds']):
                        cloudcoverage[c] = {'Cloud':'%d:%s' % (c,self.dict['Clouds'][c]['Best']),'Dataset':'1.000'}
                        for j in rje.sortKeys(self.dict['Clouds']):
                            common = [] 
                            for seq in self.dict['Clouds'][j][type]:
                                if seq in self.dict['Clouds'][c][type]: common.append(seq)
                            a = len(self.dict['Clouds'][c][type])
                            b = len(self.dict['Clouds'][j][type])
                            cloudcoverage[c][j] = '%.3f' % self.abNprob(a,b,dnum,overlap=len(common),dirn='more')
                        rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict=cloudcoverage[c])
                    ## Cloud coverage probabilites - overlaps ##
                    open(cloudfile,'a').write('\n#Probabilities of motif clouds overlapping %s this little or less\n' % type)
                    rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict={})
                    for c in rje.sortKeys(self.dict['Clouds']):
                        cloudcoverage[c] = {'Cloud':'%d:%s' % (c,self.dict['Clouds'][c]['Best']),'Dataset':'1.000'}
                        for j in rje.sortKeys(self.dict['Clouds']):
                            common = [] 
                            for seq in self.dict['Clouds'][j][type]:
                                if seq in self.dict['Clouds'][c][type]: common.append(seq)
                            a = len(self.dict['Clouds'][c][type])
                            b = len(self.dict['Clouds'][j][type])
                            cloudcoverage[c][j] = '%.3f' % self.abNprob(a,b,dnum,overlap=len(common),dirn='less')
                        rje.delimitedFileOutput(self,cloudfile,chead,delimit='\t',datadict=cloudcoverage[c])
                open(cloudfile,'a').write('\n')
            self.printLog('#CLOUD','Motif cloud data output to %s for %d clouds.' % (cloudfile,len(self.dict['Clouds'])))
        except: self.errorLog('Problem with SLiMFinder cloud output for %s' % self.dataset())
#########################################################################################################################
    def extraOutput(self):  ### Method controlling additional outputs (primarily MotifList alignments)
        '''
        Method controlling additional outputs (primarily MotifList alignments):
        - Full protein alignments (in subdirectory) with orthologues (no masking)
        - Protein Alignments with Motifs and masking marked
        - Motif alignments with and without masking
        - Occurrence data table
        - CompariMotif analysis?
        '''
        try:### ~ [1] Masked Sequence output (fasta) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            extrabase = self.runBase()
            if self.extras(1) and self.getBool('Masked'):
                for seq in self.seqs(): seq.info['Sequence'] = seq.info['MaskSeq']
                self.obj['SeqList'].saveFasta(seqfile='%s.masked.fas' % extrabase)
                self.dict['Output']['masked'] = '%s.masked.fas' % extrabase

            ### ~ [2] Setup SLiM List = SigSLiM + SLiMCheck~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            extraslim = self.list['SigSlim'][0:] + self.list['SlimCheckExtra'][0:]
            if extraslim: self.printLog('#EXTRA','Additional outputs for %d SLiMs...' % len(extraslim),log=False)
            else: return self.printLog('#EXTRA','No SLiMs for additional outputs')
            
            ### ~ [3] CompariMotif self comparison ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.extras(2):
                sigpat = []
                self.dict['Output']['motifs'] = '%s.motifs' % extrabase
                SIGSLIM = open('%s.motifs' % extrabase,'w')
                for slim in extraslim:
                    sigpat.append(patternFromCode(slim))
                    if not self.dict['Slim'][slim].has_key('Rank'): self.dict['Slim'][slim]['Rank'] = 0
                    if slim in self.list['SigSlim'][0:]:
                        SIGSLIM.write('%s # %s SLiMFinder (%s) rank %s' % (patternFromCode(slim),self.dataset(),self.getStr('RunID'),self.dict['Slim'][slim]['Rank']))
                        try: SIGSLIM.write(' [p=%s]\n' % (rje_slim.expectString(self.dict['Slim'][slim]['LSig'])))
                        except: SIGSLIM.write('\n')
                    else: SIGSLIM.write('%s # SLiMCheck Motif\n' % (patternFromCode(slim)))
                SIGSLIM.close()
                compcmd = ['motifs=%s.motifs' % extrabase,'searchdb=','resfile=%s' % extrabase,
                           'reverse=F','outstyle=reduced','dna=F','aafreq=None']
                xgmml = None
                if len(extraslim) > 1000: self.printLog('#COMP','CompariMotif skipped due to high motif numbers. (No %s.compare.tdt)' % extrabase)
                elif extraslim:
                    xgmml = comparimotif.CompariMotif(self.log,self.cmd_list+compcmd).run()
                    self.dict['Output']['compare'] = '%s.compare.tdt' % extrabase
            sigpat = []
            #!# NB. This method and these outputs have not yet been subject to general tidying and complete debugging #!#
            #!# Incorporate SlimCheck extraslims into self.obj['MotifList'] with occurrence data at some point? (PRESTO!) 
            if not self.list['SigSlim']: return self.printLog('#EXTRA','No Significant SLiMs for alignment outputs')

            ### ~ [4] Protein Alignments with Motifs and Masking marked ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.extras(1) and not self.getBool('PTMod'):  #!# Add boolean recognition of sequence modification
                (usealn,self.obj['SlimList'].opt['UseAln']) = (self.obj['SlimList'].opt['UseAln'],False)  # No alignments used
                alncmd = ['seqin=None','accnr=F','seqnr=F','autofilter=F','align=F','gnspacc=F','unkspec=T'] 
                paln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd)
                paln.info['Name'] = '%s.mapping.fas' % extrabase
                self.dict['Output']['mapping'] = paln.info['Name']
                if self.getBool('Webserver'): waln = rje_seq.SeqList(log=self.log,cmd_list=self.cmd_list+alncmd+['replacechar=F'])
                ## ~ [4a] Generate Alignment Objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                motif_occ = self.obj['SlimList'].motifOcc(byseq=True)
                for seq in self.seqs():
                    if self.getBool('Masked'): seq.info['Sequence'] = seq.info['PreMask']   #x# seq.info['MaskSeq']
                    occlist = []
                    if motif_occ.has_key(seq):
                        for Motif in motif_occ[seq].keys():
                            occlist += motif_occ[seq][Motif]
                    saln = self.obj['SlimList'].obj['SLiMCalc'].singleProteinAlignment(seq,occlist,usegopher=False,savefasta=False,wintuple=30)
                    saln.seq[0].info['Name'] = '%s-%s Motifs' % (seq.shortName(),self.dataset())
                    #self.deBug('%d:%s' % (saln.seqNum(),saln.accList()))
                    if self.getBool('Masked'):
                        alnseq = saln.seq[1].info['Sequence']
                        #x# saln.seq[1].info['Sequence'] = rje_sequence.mapGaps(seq.info['PreMask'],alnseq,self)
                        saln._addSeq('%s-masked' % seq.shortName(),rje_sequence.mapGaps(seq.info['MaskSeq'],alnseq,self))
                        saln.seq = saln.seq[:1] + saln.seq[-1:] + saln.seq[1:-1]
                    paln.seq += saln.seq[0:]
                    ## ~ [4b] Special webserver output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
                ## ~ [4c] Mapping Alignment Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                paln.saveFasta()
                self.obj['SlimList'].opt['UseAln'] = usealn

            ### ~ [5] Motif alignments with and without masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if self.extras(1) and not self.getBool('PTMod'):
                ## ~ [5a] No Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                self.obj['SlimList'].info['Name'] = '%s %s SLiMFinder Motifs' % (self.dataset(), self.getStr('RunID'))
                neworder = [None] * len(self.list['SigSlim'])
                for Motif in self.obj['SlimList'].slims()[0:]:
                    slim = self.getSigSlim(Motif.info['Sequence'])
                    if not slim or slim not in self.list['SigSlim']:
                        self.errorLog('Trying to make extra output for "%s" but not in SigSlim!' % slim,printerror=False)
                        print len(self.list['SigSlim']), self.obj['SlimList'].motifNum()
                        self.obj['SlimList'].removeMotif(Motif)
                        continue    #!# Problem?! #!#
                    elif not slim in self.dict['Slim']:
                        self.errorLog('Trying to make extra output for "%s" but missing from dictionary!' % slim,printerror=False)
                        self.obj['SlimList'].removeMotif(Motif)
                        continue
                    neworder[self.list['SigSlim'].index(slim)] = Motif
                    Motif.info['Name'] = '%s' % (rje.preZero(self.dict['Slim'][slim]['Rank'],len(self.list['SigSlim'])))
                self.obj['SlimList'].list['Motif'] = neworder[0:]
                self.obj['SlimList'].motifAlignments('%s.motifaln.fas' % extrabase)
                self.dict['Output']['motifaln'] = '%s.motifaln.fas' % extrabase
                ## ~ [5b] With Masking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
                if self.getBool('Masked'):
                    for seq in self.seqs():
                        seq.info['Sequence'] = seq.info['MaskSeq']
                    self.obj['SlimList'].motifAlignments('%s.maskaln.fas' % extrabase)
                    self.dict['Output']['maskaln'] = '%s.maskaln.fas' % extrabase
                    for seq in self.seqs():
                        seq.info['Sequence'] = seq.info['PreMask']


            ### ~ [6] Occurrence data table ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [6a] Setup XGMML - details added during Occurrence data table ~~~~~~~~~~~~~~~~~~~~ ##
            if self.extras(2): xgmml = self.setupXGMML(xgmml)
            ## ~ [6b] Setup occurrence file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            occfile = '%s.occ.csv' % extrabase  #!# Tidy this up to use delimit option?
            open(occfile,'w')
            occhead = ['Dataset','Rank','Pattern','Sig','Seq','Start_Pos','End_Pos','Prot_Len','Match','Variant','MisMatch','Desc']
            slist = self.obj['SlimList']
            for h in slist.obj['SLiMCalc'].list['OccHeaders']:
                if h not in occhead: occhead.append(h)
            occhead += ['PepSeq','PepDesign']
            rje.delimitedFileOutput(self,occfile,occhead,',')
            ## ~ [6c] Occurrence Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            motocc = self.obj['SlimList'].motifOcc()
            ox = 0
            if self.extras(2): edges = xgmml.dict['Edge']
            for slim in self.list['SigSlim']:   #!# Add SigNum? #!#
                pattern = patternFromCode(slim)
                Motif = self.obj['SlimList'].mapPattern(pattern,update=False)
                if not Motif in motocc:
                    self.errorLog('Motif "%s" missing from Motif Occurrence dictionary.' % pattern,printerror=False)
                    continue
                for Seq in motocc[Motif]:
                    for Occ in motocc[Motif][Seq]:
                        ## Add Occurrence to Data Files ##
                        occdata = {'Dataset':self.dataset(),'RunID':self.getStr('RunID'),'Rank':self.dict['Slim'][slim]['Rank'],'Desc':Seq.info['Description'],
                                   'Pattern':pattern,'PepDesign':rje_sequence.peptideDetails(Occ['PepSeq'],self),
                                   'Sig':rje_slim.expectString(self.dict['Slim'][slim]['Sig'])}
                        rje.combineDict(occdata,Occ)
                        for occstat in ['Cons','SA','Hyd','Fold','IUP','Comp']:
                            if occstat in occdata: occdata[occstat] = rje_slim.expectString(occdata[occstat])
                        occdata['Seq'] = Seq.shortName()
                        rje.delimitedFileOutput(self,occfile,occhead,',',occdata)
                        rje.delimitedFileOutput(self,self.getStr('OccResFile'),self.list['OccHeaders'],datadict=occdata)
                        ox += 1
                        ## Add Edge to XGMML ##
                        if self.extras(2): 
                            edge = (Seq.info['ID'],occdata['Pattern'])
                            if edge in edges['occ']: continue   # Only keep first occurrence
                            edges['occ'][edge] = {}
                            for stat in ['Mean_Cons','Start_Pos','Match','Variant']:
                                if stat in occdata: edges['occ'][edge][stat] = occdata[stat]
                            edges['occ'][edge]['edge.lineStyle'] = 'SOLID'
                            edges['occ'][edge]['edge.color'] = '0,0,255'
                            edges['occ'][edge]['edge.lineWidth'] = '1.0'
                            if 'Mean_Cons' in occdata: edges['occ'][edge]['edge.lineWidth'] = '%.2f' % (0.5 + string.atof(occdata['Mean_Cons']))
                            edges['occ'][edge]['edge.sourceArrowShape'] = 'DELTA'
                            edges['occ'][edge]['edge.sourceArrowColor'] = '0,0,255'
                self.printLog('\r#OCC','Occurrence data output for %s motif occurrences.' % rje.integerString(ox),newline=False,log=False)
            self.printLog('\r#OCC','Occurrence data output for %s motif occurrences.' % rje.integerString(ox))

            ### ~ [7] Output XGMML ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #!# Add option to map motif nodes onto clouds and use instead?! #!#
            if self.extras(2):
                xfile = '%s.xgmml' % xgmml.info['Name']
                self.dict['Output']['xgmml'] = xfile
                xgmml.saveXGMML(xfile)
                if self.dict['Clouds']: self.cloudXGMML(xgmml)
        except: self.errorLog('Problem with additional SLiMFinder output for %s' % self.dataset())
#########################################################################################################################
    def setupXGMML(self,xgmml=None):    ### Adds proteins and UPC to XGMML (new one if no object given)
        '''Adds proteins and UPC to XGMML (new one if no object given).'''
        try:
            ### ~ [1] Create new XGMML object if none given & update Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not xgmml: xgmml = rje_xgmml.XGMML(self.log,self.cmd_list)
            nodeatt = rje.combineDict(xgmml.dict['NodeAtt'],{'Motif':'string','Pattern':'string','Description':'string',
                                                             'PosLength':'real','MinLength':'real','MaxLength':'real',
                                                             'FixLength':'real','IC':'real','Dataset':'string',
                                                             'TYPE':'string','Sig':'string',
                                                             'Protein':'string','MaskLen':'real','UnmaskLen':'real',
                                                             'AccNum':'string','Species':'string','Rank':'real'},
                                      overwrite=False)
            for cyt in ['shape','fillColor','size','fontSize','borderColor','font']: xgmml.dict['NodeAtt']['node.%s' % cyt] = 'string' 
            edgeatt = rje.combineDict(xgmml.dict['EdgeAtt'],{'Mean_Cons':'real','Start_Pos':'real','Match':'string','Variant':'string'},overwrite=False)
            for cyt in ['lineStyle','color','sourceArrowShape','sourceArrowColor','targetArrowShape','targetArrowColor','lineWidth']: xgmml.dict['EdgeAtt']['edge.%s' % cyt] = 'string' 
            ## ~ [1b] XGMML File name ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            xgmml.info['Name'] = self.runBase()
            ## ~ [1c] CompariMotif Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            nodes = xgmml.dict['Node']
            for name in nodes:  # These will be node patterns
                nodes[name]['TYPE'] = 'Motif'
                nodes[name]['node.shape'] = 'diamond'
                nodes[name]['node.fillColor'] = '255,127,127'
                nodes[name]['Dataset'] = self.dataset()
                nsize = 25
                try:
                    try: slim = self.getSigSlim(name)   # Sig SLiM
                    except: slim = rje_slim.slimFromPattern(name)   # MotifCheck?
                    nodes[name]['Rank'] = self.dict['Slim'][slim]['Rank']
                    nodes[name]['Sig'] = rje_slim.expectString(self.dict['Slim'][slim]['Sig'])
                    sig = self.dict['Slim'][slim]['Sig']
                    if sig <= 0.05: nsize = 30
                    if sig <= 0.01: nsize = 35
                    if sig <= 0.001: nsize = 40
                    if sig <= 0.0001: nsize = 50
                except: pass
                nodes[name]['node.size'] = '%d.0' % nsize
            ## ~ [1d] CompariMotif Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edges = xgmml.dict['Edge']
            for type in edges:
                for edge in edges[type]:
                    edges[type][edge]['edge.lineStyle'] = 'SOLID'
                    edges[type][edge]['edge.color'] = '255,0,0'
                    edges[type][edge]['edge.lineWidth'] = '1.0'
                    edges[type][edge]['edge.targetArrowColor'] = edges[type][edge]['edge.color']
                    edges[type][edge]['edge.sourceArrowColor'] = edges[type][edge]['edge.color']
            edges['occ'] = {}
                
            ### ~ [2] Add Proteins and Protein Attributes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for seq in self.seqs():
                name = seq.info['ID']
                nodes[name] = {'Protein':seq.shortName(),'Description':seq.info['Description'],
                               'Species':seq.info['Species'],'AccNum':seq.info['AccNum'],
                               'MaskLen':seq.aaLen() - string.count(seq.info['MaskSeq'].upper(),'X'),'TYPE':'Protein',
                               'UnmaskLen':seq.aaLen()}
                nodes[name]['node.shape'] = 'ellipse'
                nodes[name]['node.fillColor'] = '127,127,255'
                nodes[name]['Dataset'] = self.dataset()
                nodes[name]['node.font'] = 'Impact,plain,10'
                nodes[name]['node.fontSize'] = '10.0'
                nodes[name]['node.size'] = '40.0'
            ## ~ [2a] Add UPC links between proteins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            edges['upc'] = {}
            for upc in self.list['UP']:
                for s1 in range(len(upc)):
                    for s2 in range(s1+1,len(upc)):
                        #(seq1,seq2) = (upc[s1].shortName(),upc[s2].shortName())
                        (seq1,seq2) = (upc[s1].info['ID'],upc[s2].info['ID'])
                        edges['upc'][(seq1,seq2)] = {}
                        edges['upc'][(seq1,seq2)]['edge.lineStyle'] = 'SOLID'
                        edges['upc'][(seq1,seq2)]['edge.color'] = '0,0,0'
                        edges['upc'][(seq1,seq2)]['edge.lineWidth'] = '1.0'
            
        except: self.errorLog('Problem with SLiMFinder.setupXGMML()')
        return xgmml
#########################################################################################################################
    def cloudXGMML(self,xgmml):     ### Converts full XGMML into Cloud XGMMML
        '''Converts full XGMML into Cloud XGMMML.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xfile = '%s.cloud.xgmml' % xgmml.info['Name']
            nodes = xgmml.dict['Node']
            edges = xgmml.dict['Edge']
            ## ~ [1a] Remove motif nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            for Motif in self.obj['SlimList'].motifs():
                if Motif.pattern() in nodes: nodes.pop(Motif.pattern())

            ### ~ [2] Make Cloud Nodes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for cnum in self.dict['Clouds']:
                plus = len(self.dict['Clouds'][cnum]['Slim']) - 1
                name = self.dict['Clouds'][cnum]['Best']
                if plus: name = '%s +%d' % (name,plus)
                self.dict['Clouds'][cnum]['Node'] = name
                nodes[name] = {'Motif':name,'Description':'%s SLiMFinder Cloud %s' % (self.dataset(),cnum),
                               'TYPE':'Cloud','node.shape':'hexagon','node.fillColor':'255,127,127',
                               'Dataset':self.dataset(),'node.font':'Impact,plain,10','node.fontSize':'10.0',
                               'Sig':self.dict['Clouds'][cnum]['Sig']}
                nsize = 25
                try:
                    try: slim = self.getSigSlim(self.dict['Clouds'][cnum]['Best'])   # Sig SLiM
                    except: slim = rje_slim.slimFromPattern(self.dict['Clouds'][cnum]['Best'])   # MotifCheck?
                    nodes[name]['Rank'] = self.dict['Slim'][slim]['Rank']
                    sig = self.dict['Slim'][slim]['Sig']
                    if sig <= 0.05: nsize = 30
                    if sig <= 0.01: nsize = 35
                    if sig <= 0.001: nsize = 40
                    if sig <= 0.0001: nsize = 50
                except: pass
                nodes[name]['node.size'] = '%d.0' % nsize

            ### ~ [3] Add Cloud Edges ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            for type in edges.keys()[0:]:
                if type not in ['upc']: edges.pop(type)      # Remove comparimotif edges
            edges['occ'] = {}    # Delete motif occurrences
            motocc = self.obj['SlimList'].motifOcc()
            for slim in self.list['SigSlim']:   #!# Add SigNum? #!#
                pattern = patternFromCode(slim)
                Motif = self.obj['SlimList'].mapPattern(pattern,update=False)
                cnum = self.dict['Slim'][slim]['Cloud']
                node = self.dict['Clouds'][cnum]['Node']
                if not Motif in motocc:
                    self.errorLog('Motif "%s" missing from Motif Occurrence dictionary.' % pattern,printerror=False)
                    continue
                for Seq in motocc[Motif]:
                    for Occ in motocc[Motif][Seq]:
                        edge = (Seq.info['ID'],name)
                        if edge in edges['occ']: continue   # Only keep first occurrence of cloud in sequence
                        edges['occ'][edge] = {}
                        for stat in ['Mean_Cons','Start_Pos','Match','Variant']:
                            if stat in Occ: edges['occ'][edge][stat] = Occ[stat]
                        edges['occ'][edge]['edge.lineStyle'] = 'SOLID'
                        edges['occ'][edge]['edge.color'] = '0,0,255'
                        edges['occ'][edge]['edge.lineWidth'] = '1.0'
                        if 'Mean_Cons' in Occ: edges['occ'][edge]['edge.lineWidth'] = '%.2f' % (0.5 + string.atof(Occ['Mean_Cons']))
                        edges['occ'][edge]['edge.sourceArrowShape'] = 'DELTA'
                        edges['occ'][edge]['edge.sourceArrowColor'] = '0,0,255'

            ### ~ [4] Output ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            xgmml.saveXGMML(xfile)
        except: self.errorLog('Problem during SLiMFinder.cloudXGMML()')
#########################################################################################################################
    def teiresias(self):     ### Output in TEIRESIAS format along with masked sequence (input) fasta file.
        '''
        Output in TEIRESIAS format along with masked sequence (input) fasta file.
        '''
        try:
            ### Setup ###
            if self.getStrLC('Basefile'):
                self.baseFile(rje.baseFile(self.obj['SeqList'].info['Name']))
                if self.getStrLC('ResDir'):
                    rje.mkDir(self,self.getStr('ResDir'))
                    self.baseFile(self.getStr('ResDir') + self.dataset())

            ### Output ###
            extrabase = self.runBase()
            if self.getBool('SlimDisc'): self.obj['SeqList'].info['Name'] = '%s.fasta' % extrabase
            else: self.obj['SeqList'].info['Name'] = '%s.masked.fasta' % extrabase
            for seq in self.obj['SeqList'].seq: seq.info['Sequence'] = seq.info['MaskSeq']
            self.obj['SeqList'].saveFasta()

            ## Rank File ##
            outfile = extrabase + '.out'
            OUT = open(outfile,'w')
            ## Header ##
            OUT.write(string.join(['##########################################################',
                                   '#                                                        #',
                                   '#                       FINAL RESULTS                    #',
                                   '#                                                        #',
                                   '##########################################################',''],'\n'))
            ## Motifs ##
            for slim in rje.sortKeys(self.dict['Slim']):    ### These have been pre-filtered using self.filterSLiMs() 
                OUT.write('%d\t%d\t%s' % (len(self.dict['Slim'][slim]['Occ']),self.slimUP(slim),patternFromCode(slim)))
                for (seq,pos) in self.dict['Slim'][slim]['Occ']: OUT.write(' %d %d' % (self.seqs().index(seq),pos))
                OUT.write('\n')
            OUT.close()
            self.printLog('#OUT','%s shared patterns output to %s' % (rje.integerString(len(self.dict['Slim'])),outfile))
                        
        except:
            self.errorLog('Major disaster with SLiMFinder.teiresias()')
            raise
#########################################################################################################################
    def slimCheck(self):    ### Checks given list of Motifs and adds to output
        '''Checks given list of Motifs and adds to output.'''
        try:### ~ [1] Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            if not self.getStrLC('SlimCheck'): return
            self.obj['SlimCheck'] = rje_slimlist.SLiMList(self.log,self.cmd_list+['motifs=%s' % self.getStr('SlimCheck')])
            self.obj['SlimCheck'].loadMotifs()
            if not self.obj['SlimCheck'].motifs(): return
            if self.getStrLC('ResFile'):
                self.errorLog('Cannot check SLiMs without resfile=FILE',printerror=False)
                return

            ### ~ [2] Make into SLiMFinder SLiMs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            #X#aa = string.split('A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y,^,$',',')
            checklist = []
            for Motif in self.obj['SlimCheck'].slims()[0:]:
                if Motif.slim(): checklist.append(Motif.slim())
                else: self.obj['SlimCheck'].removeMotif(Motif)
            if not self.obj['SlimCheck'].slims(): return

            ### ~ [3] Calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            X = self.obj['SlimCheck'].motifNum() * self.UPNum()
            x = 0.0
            self.list['SlimCheckExtra'] = []   # Extra SLiMs
            for Motif in self.obj['SlimCheck'].motifs()[0:]:
                self.progLog('\r#CHECK','SLiM Check (%d motifs) %.2f%%' % (self.obj['SlimCheck'].motifNum(),x/X))
                slim = Motif.slim()
                if slim not in self.dict['Slim']:
                    try:
                        ux = 0
                        sx = 0
                        ox = 0
                        for upc in self.list['UP']:
                            uph = False
                            for Seq in upc:
                                hits = Motif.searchSequence(sequence=Seq.info['MaskSeq'])
                                if hits:
                                    uph = True
                                    sx += 1
                                    ox += len(hits)
                            if uph: ux += 1
                        self.dict['Slim'][slim] = {'Occ':[1] * ox,'UP':[1] * ux,'Support':sx}
                        self.slimProb(slim)
                        e = self.dict['Extremf.'][slimPos(slim)] * self.dict['Slim'][slim]['Prob']
                        self.dict['Slim'][slim]['Sig'] = rje.poisson(1,e,callobj=self)
                        self.list['SlimCheckExtra'].append(slim)
                    except:
                        self.errorLog('SLiM Error (%s)' % slim,quitchoice=False)
                        self.obj['SlimCheck'].removeMotif(Motif)
                        continue
                else:
                    self.dict['Slim'][slim]['Support'] = self.slimSeqNum(slim)
                    self.dict['Slim'][slim]['UP'] = [1] * self.slimUP(slim)
                    if slim not in self.list['SigSlim']: self.list['SlimCheckExtra'].append(slim)
            self.printLog('\r#CHECK','SLiM Check (%d motifs) complete: %d extra SLiMs.' % (self.obj['SlimCheck'].motifNum(),len(self.list['SlimCheckExtra'])))
                    
            ### ~ [4] Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            ## ~ [4a] Setup Results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            delimit = rje.getDelimit(self.cmd_list,rje.delimitFromExt(filename=self.getStr('ResFile')))
            reshead = self.resHead()
            totalaa = 0     # Total number of AA in dataset
            for seq in self.seqs():
                if self.getBool('Masked'): seq.info['Sequence'] = seq.info['MaskSeq'][0:]
                totalaa += seq.nonX()
            masking = self.maskText()   # Summary of Masking Options
            ## ~ [4b] Output for each checks SLiM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
            #!# Add SliMCalc, OccFilter and StatFilter #!#
            for Motif in self.obj['SlimCheck'].motifs():
                slim = Motif.slim()
                pattern = patternFromCode(slim)
                datadict = {'Dataset':self.dataset(),'RunID':self.getStr('RunID'),
                            'Masking':masking,'Build':self.buildText(),'RunTime':self.getStr('RunTime'),
                            'SeqNum':self.seqNum(),'UPNum':self.UPNum(),'AANum':totalaa,
                            'MotNum':len(self.dict['Slim']),'Rank':'*',
                            'Pattern':pattern,'Occ':len(self.dict['Slim'][slim]['Occ']),'Support':self.dict['Slim'][slim]['Support'],
                            'IC':self.slimIC(slim),'UP':self.slimUP(slim), 'Norm':self.slimUP(slim),
                            }
                if self.dict['Clouds'] and self.dict['Slim'][slim].has_key('Cloud'):
                    try:
                        datadict['Cloud'] = self.dict['Slim'][slim]['Cloud']
                        datadict['CloudSeq'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['Seq'])
                        datadict['CloudUP'] = len(self.dict['Clouds'][self.dict['Slim'][slim]['Cloud']]['UPC'])
                    except: pass
                for p in ['ExpUP','Prob','Sig','E','SigV','SigPrime','SigPrimeV','S','R']:
                    if self.dict['Slim'][slim].has_key(p): datadict[p] = rje_slim.expectString(self.dict['Slim'][slim][p])
                for h in self.list['Headers']:
                    if Motif and h not in datadict and h in Motif.stat: datadict[h] = Motif.stat[h]
                rje.delimitedFileOutput(self,self.getStr('ResFile'),reshead,delimit,datadict)

            self.printLog('\r#OUT','SLiMCheck Output for %d motifs into %s complete.' % (self.obj['SlimCheck'].motifNum(),self.getStr('ResFile')))

        except: self.errorLog('Major problem with SLiMFinder.slimCheck(%s)' % self.getStr('SlimCheck'))
#########################################################################################################################
    def restSetup(self):    ### Sets up self.dict['Output'] and associated output options if appropriate.
        '''
        The standared REST server call for SLiMFinder is in the form:
        `slimfinder`&uniprotid=LIST&dismask=T/F&consmask=T/F

        Different sources of input can also be given with:
        `slimfinder`&seqin=LIST&dismask=T/F&consmask=T/F

        Run with &rest=help for general options. Run with &rest=full to get full server output as plain text. Otherwise,
        individual outputs are parsed and presented in different tabs:

        ### Outputs available:
            main = main results file (extras=-1)
            seqin = Input file (extras=-1)
            occ = occurrence file (extras=0)
            upc = UPC file (extras=0)
            slimdb = Fasta file used for UPC generation etc. (extras=0)
            cloud = cloud.txt (extras=1)
            masked = masked.fas (extras=1)
            mapping = mapping.fas file (extras=1)
            motifaln = motif alignments (extras=1)
            maskaln = masked motif alignments (extras=1)
            motifs = motifs file for CompariMotif (extras=2)
            compare = CM compare.tdt file (extras=2)
            xgmml = XGMML file (extras=2)
            dismatrix = *.dis.tdt file (extras=3)
            rank = optional SLiMDisc output (extras=3)
            dat.rank = optional SLiMDisc output (extras=3)

        &rest=OUTFMT can then be used to retrieve individual parts of the output in future.
        '''
        try:### ~ [0] ~ Setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
            outfmt = self.getStrLC('Rest')
            if not self.extras(0) and outfmt in ['main','occ','upc','slimdb','seqin']: self.setInt({'Extras':0})
            if not self.extras(1) and outfmt in ['cloud','masked','mapping','motifaln','maskaln']: self.setInt({'Extras':1})
            if not self.extras(2) and outfmt in ['motifs','compare','xgmml']: self.setInt({'Extras':2})
            if not self.extras(3) and outfmt in ['dismatrix','rank','dat.rank']: self.setInt({'Extras':3})
            if outfmt == 'default' and not self.extras(2): self.setInt({'Extras':2})
            return
        except: self.errorLog('RestSetup error')
#########################################################################################################################
    def restOutputOrder(self):
        output = ['csv','seqin']
        if self.extras(0): output = ['main','occ','upc']
        if self.extras(1): output += ['cloud','masked','mapping','motifaln','maskaln']
        if self.extras(2): output += ['motifs','compare','xgmml']
        if self.extras(3): output += ['dismatrix','rank','dat.rank']
        if self.extras(0): output += ['slimdb','seqin']
        return output
#########################################################################################################################
    def restOutputDefault(self): return 'full'
#########################################################################################################################
### End of SECTION II: SLiMFinder Class                                                                                 #
#########################################################################################################################

                                                    ### ~ ### ~ ###

#########################################################################################################################
### SECTION III: SPECIFIC METHODS                                                                                       #
#########################################################################################################################
def patternFromCode(slim): return rje_slim.patternFromCode(slim)  ### Returns pattern with wildcard for iXj formatted SLiM (e.g. A-3-T-0-G becomes A...TG)
#########################################################################################################################
def slimPos(slim): return (string.count(slim,'-') / 2) + 1  ### Returns the number of positions in a slim
#########################################################################################################################
def slimLen(slim): return len(patternFromCode(slim))    ### Returns length of slim
#########################################################################################################################
def slimDif(slim1,slim2): return rje_slimcore.slimDif(slim1,slim2)  ### Returns no. of dif. pos. between slim1 and slim2
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
    try: SLiMFinder(mainlog,cmd_list).run()
        
    ### End ###
    except SystemExit: pass    #!#return  # Fork exit etc.
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
